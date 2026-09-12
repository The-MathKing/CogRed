"""
Pipeline validation via actives/decoys enrichment (manuscript Section
2.5/3.2) -- "how do you know your docking pipeline works?" is close to a
guaranteed reviewer question, and this is the standard, citable answer.

Uses known ChEMBL actives for the target (from `02_seed_compound_library.
fetch_chembl_actives`) plus property-matched decoys generated in the DUD-E
style (same MW/logP/HBD/HBA bins, dissimilar topology), docks both sets with
the same protocol as the real screen, and reports ROC-AUC and enrichment
factor at 1%/5% -- report these numbers in the manuscript before presenting
any results on the novel derivative library.
"""
from __future__ import annotations

import random

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors
from sklearn.metrics import roc_auc_score


def generate_property_matched_decoys(actives_smiles: list[str], decoy_pool_smiles: list[str],
                                      n_decoys_per_active: int = 50,
                                      mw_tol: float = 25.0, logp_tol: float = 1.0,
                                      tanimoto_max: float = 0.35, seed: int = 42) -> list[str]:
    """Simplified DUD-E-style decoy selection: for each active, pull decoys
    from a large background pool (e.g., a random ZINC/ChEMBL "drug-like"
    subset supplied by the caller) matched on MW/logP but topologically
    dissimilar (low Tanimoto similarity on Morgan fingerprints)."""
    rng = random.Random(seed)
    fpgen = AllChem.GetMorganGenerator(radius=2, fpSize=2048)

    def props(smi):
        mol = Chem.MolFromSmiles(smi)
        return (Descriptors.MolWt(mol), Descriptors.MolLogP(mol),
                fpgen.GetFingerprint(mol)) if mol else None

    active_props = [p for p in (props(s) for s in actives_smiles) if p]
    pool_props = [(s, props(s)) for s in decoy_pool_smiles]
    pool_props = [(s, p) for s, p in pool_props if p]

    decoys = set()
    for mw_a, logp_a, fp_a in active_props:
        candidates = [
            s for s, (mw_d, logp_d, fp_d) in pool_props
            if abs(mw_d - mw_a) <= mw_tol and abs(logp_d - logp_a) <= logp_tol
            and DataStructs.TanimotoSimilarity(fp_a, fp_d) <= tanimoto_max
        ]
        rng.shuffle(candidates)
        decoys.update(candidates[:n_decoys_per_active])
    return sorted(decoys)


def compute_enrichment_factor(labels_sorted_by_score: np.ndarray, top_fraction: float) -> float:
    """labels_sorted_by_score: 1=active, 0=decoy, sorted best-affinity first."""
    n = len(labels_sorted_by_score)
    top_n = max(1, int(round(n * top_fraction)))
    hits_in_top = labels_sorted_by_score[:top_n].sum()
    total_actives = labels_sorted_by_score.sum()
    expected_hits = total_actives * top_fraction
    return (hits_in_top / expected_hits) if expected_hits > 0 else float("nan")


def validate_pipeline(docking_results: pd.DataFrame, actives_smiles: set[str]) -> dict:
    """`docking_results` must have columns: smiles, affinity_best (kcal/mol,
    more negative = better), for the combined actives+decoys docking run."""
    df = docking_results.copy()
    df["is_active"] = df["smiles"].isin(actives_smiles).astype(int)
    df = df.sort_values("affinity_best")  # best (most negative) first

    auc = roc_auc_score(df["is_active"], -df["affinity_best"])
    ef1 = compute_enrichment_factor(df["is_active"].to_numpy(), 0.01)
    ef5 = compute_enrichment_factor(df["is_active"].to_numpy(), 0.05)

    return {"roc_auc": auc, "ef_1pct": ef1, "ef_5pct": ef5, "n_actives": int(df["is_active"].sum()),
            "n_decoys": int((1 - df["is_active"]).sum())}


if __name__ == "__main__":
    # Expects docking already run on the combined actives+decoys set using
    # the same protocol/box as 05_docking_pipeline.py, with results at:
    docking_results = pd.read_csv("data/validation_docking_results_summary.csv")
    actives = set(pd.read_csv("data/chembl_actives.csv")["standardized_smiles"])

    metrics = validate_pipeline(docking_results, actives)
    print(metrics)
    pd.DataFrame([metrics]).to_csv("data/pipeline_validation_metrics.csv", index=False)
