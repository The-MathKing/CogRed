"""
Predicted cross-kinase selectivity analysis (manuscript Methods 2.9).
Rewritten after external review.

Two corrections from the previous version:

  1. TERMINOLOGY. This is not a kinome selectivity panel -- no kinase assay is
     performed. Everything here is a prediction from cross-docking.

  2. SCORE COMPARABILITY. The previous implementation computed
     `selectivity = offtarget_affinity - target_affinity` on raw Vina scores.
     That is invalid: Vina scores are not calibrated across different
     receptors. Binding-site volume, buriedness, and atom count all shift a
     receptor's score distribution, so -9.1 against kinase A and -8.3 against
     kinase B does not mean 0.8 kcal/mol of selectivity.

     The fix is to normalize WITHIN each receptor: dock a shared
     property-matched background set into every kinase, then express each
     candidate's score as a z-score (or percentile) against that receptor's own
     background distribution. Selectivity is then a comparison of standardized
     ranks, which is a defensible quantity.

  3. METRIC VALIDATION. A selectivity metric that has never been checked
     against known answers is decoration. `benchmark_against_reference_profiles`
     scores reference inhibitors whose experimental kinome profiles are
     published (e.g. Davis et al., Nat Biotechnol 2011; Klaeger et al.,
     Science 2017) and reports whether the predicted ordering recovers the
     experimental one. If it does not, the heatmap in Figure 6 should not be
     interpreted quantitatively -- say so in the manuscript.

Two further corrections (external review round 2, manuscript §2.9):

  4. PKA MISCLASSIFICATION. An earlier version of this module called PKA a
     "structurally distant negative control." That is wrong: CHK1 belongs to
     the CAMK group of the human kinome and PKA is the AGC-group archetype of
     the SAME bilobal eukaryotic protein kinase fold (Manning et al., Science
     2002). PKA is a near neighbour of CHK1 by fold, not a distant one. It is
     kept in the panel as a FOLD COMPARATOR (does a candidate's selectivity
     liability track fold similarity rather than target identity?), and
     `DISTANT_CONTROL` below adds a genuinely fold-unrelated ATP-binding
     protein so the panel actually has the negative control it claimed to
     have.

  5. REFERENCE-SET SIZE. Three reference compounds cannot produce an
     interpretable Spearman correlation -- its confidence interval spans most
     of [-1, 1]. `REFERENCE_COMPOUNDS` is now sized to require the caller to
     supply 20-50 compounds drawn from Davis et al. (2011) / Klaeger et al.
     (2017) before `benchmark_against_reference_profiles` is treated as
     validation rather than illustration; the three-compound dict below is
     kept only as a placeholder/smoke-test and `MIN_REFERENCE_COMPOUNDS`
     enforces this at call time.
"""
from __future__ import annotations

from importlib import import_module

import numpy as np
import pandas as pd
import pypdb
from scipy.stats import spearmanr

docking = import_module("05_docking_pipeline")

# Representative structures -- verify current best-resolution, ligand-bound
# depositions on RCSB before running; each needs its own ATP-pocket box.
OFFTARGET_PANEL = {
    "mTOR": "4JSX",
    "DNA-PKcs": "7NI9",
    "ATM": "7SID",
    "CDK1": "4YC3",
    "CDK2": "1HCK",
    "WEE1": "5V5V",
    "PLK1": "2OWB",
    "PKA": "1ATP",  # AGC-group FOLD COMPARATOR to CHK1 (CAMK group), not a
                    # distant negative control -- see module docstring item 4.
}

# A genuinely fold-distant ATP-binding protein, to supply the negative
# control PKA does not provide. [PENDING: pick and verify a specific PDB ID
# and confirm the fold classification before running -- not yet selected.]
DISTANT_CONTROL: dict[str, str] = {
    # e.g. "Hsp90": "<PDB ID>"  -- an ATP-binding chaperone with no
    # relationship to the protein kinase fold, unlike every entry above.
}

# Reference compounds with published experimental kinome selectivity profiles,
# used to test whether the predicted metric tracks reality at all. This
# three-compound dict is a placeholder/smoke-test only -- see item 5 above.
# Before treating `benchmark_against_reference_profiles` as validation
# (rather than illustration), expand to 20-50 compounds drawn from Davis et
# al. (2011) / Klaeger et al. (2017).
MIN_REFERENCE_COMPOUNDS = 20

REFERENCE_COMPOUNDS = {
    "staurosporine": "broadly promiscuous",
    "berzosertib": "ATR-selective",
    "prexasertib": "CHK1/CHK2-selective",
}


def fetch_pdb_structure(pdb_id: str, out_dir: str = "data/receptors") -> str:
    pdb_text = pypdb.get_pdb_file(pdb_id, filetype="pdb", compression=False)
    out_path = f"{out_dir}/{pdb_id}.pdb"
    with open(out_path, "w") as fh:
        fh.write(pdb_text)
    return out_path


def dock_background_set(receptor_pdb: str, kinase: str, box: docking.DockingBox,
                         background_csv: str = "data/background_property_matched.csv") -> pd.DataFrame:
    """Dock a shared background compound set into one receptor.

    The SAME background set must be used for every kinase -- it is the ruler
    against which each receptor's score distribution is measured.
    """
    _, summary = docking.dock_library(
        compound_csv=background_csv, receptor_pdb=receptor_pdb,
        box=box, target_name=f"{kinase}_background", n_replicates=1,
    )
    return summary


def normalize_within_receptor(candidate_scores: pd.Series,
                               background_scores: pd.Series) -> pd.DataFrame:
    """Express candidate affinities as z-scores and percentiles against the
    background distribution docked into that same receptor.

    More negative Vina affinity = better, so the sign is flipped such that
    higher z = stronger predicted engagement.
    """
    mu, sigma = background_scores.mean(), background_scores.std(ddof=1)
    z = -(candidate_scores - mu) / sigma if sigma > 0 else pd.Series(np.nan, index=candidate_scores.index)
    pct = candidate_scores.apply(lambda s: (background_scores > s).mean())
    return pd.DataFrame({"affinity": candidate_scores, "z_score": z, "percentile": pct})


def run_selectivity_analysis(shortlist_csv: str, primary_target: str,
                              boxes: dict[str, docking.DockingBox],
                              panel: dict[str, str] = OFFTARGET_PANEL) -> pd.DataFrame:
    """Cross-dock shortlisted candidates across the panel, normalizing per receptor.

    `boxes` must supply a per-kinase ATP-pocket box (derived from each
    structure's co-crystal ligand centroid). Reusing one box across different
    receptors is not valid.
    """
    frames = []
    for kinase, pdb_id in panel.items():
        if kinase not in boxes:
            raise ValueError(f"No docking box supplied for {kinase}; derive it from {pdb_id}.")
        receptor_path = fetch_pdb_structure(pdb_id)
        box = boxes[kinase]

        _, candidates = docking.dock_library(
            compound_csv=shortlist_csv, receptor_pdb=receptor_path,
            box=box, target_name=kinase, n_replicates=3,
        )
        background = dock_background_set(receptor_path, kinase, box)

        norm = normalize_within_receptor(
            candidates.set_index("smiles")["affinity_best"],
            background["affinity_best"],
        ).reset_index()
        norm["kinase"] = kinase
        frames.append(norm)

    return pd.concat(frames, ignore_index=True)


def compute_selectivity_scores(normalized: pd.DataFrame, primary_target: str) -> pd.DataFrame:
    """Selectivity = primary-target z-score minus the best off-target z-score.

    Because both terms are standardized within their own receptor, the
    difference is interpretable as a relative-rank margin -- NOT as a free
    energy difference, and it must never be reported in kcal/mol.
    """
    primary = normalized[normalized["kinase"] == primary_target].set_index("smiles")["z_score"]
    offtarget = normalized[normalized["kinase"] != primary_target]
    best_off = offtarget.groupby("smiles")["z_score"].max()

    df = pd.DataFrame({
        "primary_z": primary,
        "best_offtarget_z": best_off,
    }).dropna()
    df["selectivity_margin_z"] = df["primary_z"] - df["best_offtarget_z"]
    df["most_liable_offtarget"] = (
        offtarget.loc[offtarget.groupby("smiles")["z_score"].idxmax()]
        .set_index("smiles")["kinase"]
    )
    return df.sort_values("selectivity_margin_z", ascending=False)


def benchmark_against_reference_profiles(normalized: pd.DataFrame,
                                          experimental_ranks: dict[str, dict[str, float]],
                                          require_min_compounds: bool = True) -> dict:
    """Does the predicted selectivity ordering recover experimental kinome data?

    `experimental_ranks`: {compound: {kinase: experimental_affinity_or_rank}}
    taken from published profiling (Davis 2011 / Klaeger 2017). Returns Spearman
    correlation per reference compound. Weak or negative correlation means the
    predicted heatmap must not be read quantitatively -- report that outcome
    rather than dropping the benchmark.

    `require_min_compounds`: with fewer than MIN_REFERENCE_COMPOUNDS reference
    compounds, any correlation computed is not interpretable as validation
    (its CI spans most of [-1, 1]) -- by default this raises rather than
    silently returning a number that looks like validation but isn't. Pass
    False only to run an explicitly-labelled illustrative smoke test.
    """
    if require_min_compounds and len(experimental_ranks) < MIN_REFERENCE_COMPOUNDS:
        raise ValueError(
            f"Only {len(experimental_ranks)} reference compounds supplied; "
            f"need >= {MIN_REFERENCE_COMPOUNDS} for an interpretable "
            f"correlation (manuscript §2.9). Pass require_min_compounds=False "
            f"to run anyway as an explicitly illustrative-only smoke test."
        )
    out = {}
    for compound, exp_profile in experimental_ranks.items():
        pred = normalized[normalized["smiles"] == compound].set_index("kinase")["z_score"]
        shared = [k for k in exp_profile if k in pred.index]
        if len(shared) < 3:
            out[compound] = {"n_shared": len(shared), "spearman": None}
            continue
        rho, p = spearmanr([pred[k] for k in shared], [exp_profile[k] for k in shared])
        out[compound] = {"n_shared": len(shared), "spearman": float(rho), "p_value": float(p)}
    return out


if __name__ == "__main__":
    # Each kinase needs its own box derived from its co-crystal ligand centroid.
    boxes = {k: docking.DockingBox(0.0, 0.0, 0.0) for k in OFFTARGET_PANEL}
    boxes["CHK1"] = docking.DockingBox(0.0, 0.0, 0.0)

    normalized = run_selectivity_analysis(
        shortlist_csv="data/shortlist_top_hits.csv",
        primary_target="CHK1", boxes=boxes,
    )
    normalized.to_csv("data/selectivity_normalized_scores.csv", index=False)

    selectivity = compute_selectivity_scores(normalized, primary_target="CHK1")
    selectivity.to_csv("data/selectivity_scores.csv")
    print(selectivity.head(10))
