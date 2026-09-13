"""
Leakage-controlled evaluation of the screening pipeline (manuscript Methods
2.5 / Results 3.1-3.2). Rewritten after external review.

Three things the previous version got wrong:

  1. No separation between the seed chemistry used to build the derivative
     library and the "held-out" actives used to evaluate retrieval. Derivatives
     generated from seed A being used to validate retrieval of seed-A-like
     actives inflates every metric reported.
  2. A single ROC curve, with no baseline. "Our final compounds score well"
     is not evidence that a multi-stage pipeline adds anything over plain
     docking.
  3. ROC-AUC alone under heavy class imbalance (1 active : 50 decoys) is
     optimistic; PR-AUC and early-enrichment metrics are the informative ones.

This module therefore provides: Bemis-Murcko scaffold splitting, an explicit
seed<->evaluation leakage audit, baseline-vs-pipeline comparison, and
bootstrap confidence intervals.

PRE-SPECIFICATION: splits, metrics, and comparisons are fixed here before
results are seen, and the outcome is reported in whichever direction it falls.
A finding that the full pipeline does NOT outperform baseline docking is a
result to report, not a failure to iterate away (see docs/04_review_response.md
section B).
"""
from __future__ import annotations

import random
from dataclasses import dataclass

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.metrics import average_precision_score, roc_auc_score

RDLogger.DisableLog("rdApp.*")

_FPGEN = AllChem.GetMorganGenerator(radius=2, fpSize=2048)


def _fingerprint(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    return _FPGEN.GetFingerprint(mol) if mol else None


def murcko_scaffold(smiles: str) -> str | None:
    """Bemis-Murcko scaffold SMILES, used as the grouping key for splitting."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return MurckoScaffold.MurckoScaffoldSmiles(mol=mol)


# --------------------------------------------------------------------------
# 1. Scaffold splitting
# --------------------------------------------------------------------------

def scaffold_split(actives: list[str], test_fraction: float = 0.3,
                   seed: int = 20260913) -> tuple[list[str], list[str]]:
    """Split actives so that no Murcko scaffold appears in both partitions.

    A random split would place close analogs on both sides and make retrieval
    look far better than it is -- this is the single most important control in
    the whole evaluation.
    """
    by_scaffold: dict[str, list[str]] = {}
    for smi in actives:
        scaf = murcko_scaffold(smi)
        if scaf is None:
            continue
        by_scaffold.setdefault(scaf, []).append(smi)

    scaffolds = sorted(by_scaffold, key=lambda s: (-len(by_scaffold[s]), s))
    rng = random.Random(seed)
    rng.shuffle(scaffolds)

    n_target_test = int(round(len(actives) * test_fraction))
    test, train = [], []
    for scaf in scaffolds:
        if len(test) < n_target_test:
            test.extend(by_scaffold[scaf])
        else:
            train.extend(by_scaffold[scaf])
    return train, test


# --------------------------------------------------------------------------
# 2. Leakage audit
# --------------------------------------------------------------------------

@dataclass
class LeakageReport:
    n_evaluated: int
    n_scaffold_collisions: int
    max_tanimoto_to_seed: float
    median_tanimoto_to_seed: float
    n_above_threshold: int
    threshold: float

    def as_dict(self) -> dict:
        return self.__dict__


def audit_leakage(seed_smiles: list[str], evaluation_smiles: list[str],
                  tanimoto_threshold: float = 0.7) -> LeakageReport:
    """Quantify how close the evaluation set is to the seed/derivative chemistry.

    Reports scaffold collisions and nearest-neighbour Tanimoto similarity of
    every evaluation compound to its closest seed compound. These numbers
    belong in the manuscript Results, not in Supporting Information -- they
    are what tells a reader whether the retrieval metrics mean anything.
    """
    seed_scaffolds = {murcko_scaffold(s) for s in seed_smiles} - {None}
    seed_fps = [fp for fp in (_fingerprint(s) for s in seed_smiles) if fp]

    collisions = 0
    nn_similarities = []
    for smi in evaluation_smiles:
        if murcko_scaffold(smi) in seed_scaffolds:
            collisions += 1
        fp = _fingerprint(smi)
        if fp and seed_fps:
            nn_similarities.append(max(DataStructs.BulkTanimotoSimilarity(fp, seed_fps)))

    sims = np.array(nn_similarities) if nn_similarities else np.array([np.nan])
    return LeakageReport(
        n_evaluated=len(evaluation_smiles),
        n_scaffold_collisions=collisions,
        max_tanimoto_to_seed=float(np.nanmax(sims)),
        median_tanimoto_to_seed=float(np.nanmedian(sims)),
        n_above_threshold=int(np.nansum(sims >= tanimoto_threshold)),
        threshold=tanimoto_threshold,
    )


# --------------------------------------------------------------------------
# 3. Decoy construction
# --------------------------------------------------------------------------

def generate_property_matched_decoys(actives_smiles: list[str], decoy_pool_smiles: list[str],
                                      n_decoys_per_active: int = 50,
                                      mw_tol: float = 25.0, logp_tol: float = 1.0,
                                      tanimoto_max: float = 0.35, seed: int = 42) -> list[str]:
    """DUD-E-style decoys: matched on MW/logP, topologically dissimilar.

    `decoy_pool_smiles` should be a large drug-like background set (e.g. a
    random ZINC or ChEMBL subset) supplied by the caller.
    """
    rng = random.Random(seed)

    def props(smi):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return None
        return (Descriptors.MolWt(mol), Descriptors.MolLogP(mol), _FPGEN.GetFingerprint(mol))

    active_props = [p for p in (props(s) for s in actives_smiles) if p]
    pool_props = [(s, p) for s, p in ((s, props(s)) for s in decoy_pool_smiles) if p]

    decoys: set[str] = set()
    for mw_a, logp_a, fp_a in active_props:
        candidates = [
            s for s, (mw_d, logp_d, fp_d) in pool_props
            if abs(mw_d - mw_a) <= mw_tol and abs(logp_d - logp_a) <= logp_tol
            and DataStructs.TanimotoSimilarity(fp_a, fp_d) <= tanimoto_max
        ]
        rng.shuffle(candidates)
        decoys.update(candidates[:n_decoys_per_active])
    return sorted(decoys)


# --------------------------------------------------------------------------
# 4. Metrics
# --------------------------------------------------------------------------

def enrichment_factor(labels_ranked: np.ndarray, top_fraction: float) -> float:
    """labels_ranked: 1=active, 0=decoy, ordered best-scoring first."""
    n = len(labels_ranked)
    top_n = max(1, int(round(n * top_fraction)))
    total_actives = labels_ranked.sum()
    expected = total_actives * top_fraction
    return float(labels_ranked[:top_n].sum() / expected) if expected > 0 else float("nan")


def compute_metrics(scores: np.ndarray, labels: np.ndarray) -> dict:
    """`scores`: higher = better ranked. Convert Vina affinities with -affinity."""
    order = np.argsort(-scores)
    ranked_labels = labels[order]
    return {
        "roc_auc": float(roc_auc_score(labels, scores)),
        "pr_auc": float(average_precision_score(labels, scores)),
        "ef_1pct": enrichment_factor(ranked_labels, 0.01),
        "ef_5pct": enrichment_factor(ranked_labels, 0.05),
    }


def bootstrap_metrics(scores: np.ndarray, labels: np.ndarray, n_boot: int = 2000,
                       seed: int = 20260913) -> pd.DataFrame:
    """Bootstrap CIs -- required to claim one method beats another."""
    rng = np.random.default_rng(seed)
    rows = []
    n = len(scores)
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        if labels[idx].sum() == 0 or labels[idx].sum() == len(idx):
            continue
        rows.append(compute_metrics(scores[idx], labels[idx]))
    boot = pd.DataFrame(rows)
    return boot.agg(["mean", lambda s: s.quantile(0.025), lambda s: s.quantile(0.975)]).T.set_axis(
        ["mean", "ci_low", "ci_high"], axis=1)


# --------------------------------------------------------------------------
# 5. Baseline vs. pipeline comparison
# --------------------------------------------------------------------------

def compare_methods(results: pd.DataFrame, label_col: str = "is_active",
                     score_cols: tuple[str, ...] = ("baseline_vina_score",
                                                     "consensus_score",
                                                     "full_pipeline_score")) -> pd.DataFrame:
    """Evaluate each stage of the pipeline on identical compounds and splits.

    Without this comparison the paper cannot distinguish 'the pipeline works'
    from 'the pipeline is an elaborate way of running Vina'.
    """
    labels = results[label_col].to_numpy()
    rows = []
    for col in score_cols:
        if col not in results.columns:
            continue
        scores = results[col].to_numpy()
        metrics = compute_metrics(scores, labels)
        ci = bootstrap_metrics(scores, labels)
        metrics.update({f"{m}_ci": (round(ci.loc[m, "ci_low"], 3),
                                     round(ci.loc[m, "ci_high"], 3))
                        for m in ("roc_auc", "pr_auc") if m in ci.index})
        metrics["method"] = col
        rows.append(metrics)
    return pd.DataFrame(rows).set_index("method")


def paired_bootstrap_difference(results: pd.DataFrame, method_a: str, method_b: str,
                                 metric: str = "roc_auc", label_col: str = "is_active",
                                 n_boot: int = 2000, seed: int = 20260913) -> dict:
    """Paired bootstrap on the difference (method_a - method_b).

    If the 95% CI of the difference includes zero, the pipeline does not
    demonstrably outperform the baseline -- report that plainly.
    """
    rng = np.random.default_rng(seed)
    labels = results[label_col].to_numpy()
    a, b = results[method_a].to_numpy(), results[method_b].to_numpy()
    diffs = []
    n = len(labels)
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        if labels[idx].sum() in (0, len(idx)):
            continue
        diffs.append(compute_metrics(a[idx], labels[idx])[metric]
                      - compute_metrics(b[idx], labels[idx])[metric])
    diffs = np.array(diffs)
    ci_low, ci_high = np.quantile(diffs, [0.025, 0.975])
    return {
        "metric": metric,
        "comparison": f"{method_a} - {method_b}",
        "mean_difference": float(diffs.mean()),
        "ci_low": float(ci_low),
        "ci_high": float(ci_high),
        "significant_at_95pct": bool(ci_low > 0 or ci_high < 0),
    }


if __name__ == "__main__":
    # Expects docking already run on actives + decoys under the SAME protocol
    # and box as the real screen, with one score column per pipeline stage.
    results = pd.read_csv("data/validation_docking_results_summary.csv")
    seed_compounds = pd.read_csv("data/seed_named_compounds.csv")["standardized_smiles"].dropna().tolist()
    actives = pd.read_csv("data/chembl_actives.csv")["standardized_smiles"].dropna().tolist()

    train_actives, test_actives = scaffold_split(actives)
    leakage = audit_leakage(seed_compounds, test_actives)
    print("Leakage audit:", leakage.as_dict())

    comparison = compare_methods(results)
    print(comparison)

    delta = paired_bootstrap_difference(results, "full_pipeline_score", "baseline_vina_score")
    print("Pipeline vs. baseline:", delta)
    if not delta["significant_at_95pct"]:
        print("NOTE: no demonstrable advantage over baseline docking. "
              "Report this as the result (docs/04_review_response.md section B).")

    comparison.to_csv("data/pipeline_validation_metrics.csv")
    pd.DataFrame([leakage.as_dict()]).to_csv("data/leakage_audit.csv", index=False)
    pd.DataFrame([delta]).to_csv("data/baseline_comparison.csv", index=False)
