"""
Off-target kinome selectivity panel (manuscript Section 2.9).

Re-uses the docking wrapper from `05_docking_pipeline.py` to cross-dock the
shortlisted hits against a curated panel of structurally related kinases,
then computes a selectivity score. Panel choice rationale:

  - PIKK family (shares ATR's fold/ATP pocket architecture): mTOR, DNA-PKcs,
    ATM. Cross-reactivity here is the most mechanistically likely liability.
  - Cell-cycle kinases downstream/parallel to CHK1 signaling: CDK1, CDK2,
    WEE1, PLK1. Off-target hits here would confound phenotypic interpretation
    in any follow-up cellular assay.
  - A distant-fold negative control kinase (e.g., PKA/PRKACA) to sanity-check
    that the panel isn't trivially promiscuous.

PDB IDs below are illustrative -- verify current best-resolution, ligand-bound
structures on RCSB PDB for each kinase before running (structures are
periodically superseded by higher-resolution depositions).
"""
from __future__ import annotations

import pandas as pd
import pypdb

from importlib import import_module

docking = import_module("05_docking_pipeline")

OFFTARGET_PANEL = {
    # kinase: representative PDB ID (verify/replace with current best structure)
    "mTOR": "4JSX",
    "DNA-PKcs": "7NI9",
    "ATM": "7SID",
    "CDK1": "4YC3",
    "CDK2": "1HCK",
    "WEE1": "5V5V",
    "PLK1": "2OWB",
    "PKA": "1ATP",  # distant-fold negative control
}


def fetch_pdb_structure(pdb_id: str, out_dir: str = "data/receptors") -> str:
    """Download a PDB structure file via pypdb/RCSB for receptor prep."""
    pdb_text = pypdb.get_pdb_file(pdb_id, filetype="pdb", compression=False)
    out_path = f"{out_dir}/{pdb_id}.pdb"
    with open(out_path, "w") as fh:
        fh.write(pdb_text)
    return out_path


def run_selectivity_panel(shortlist_csv: str, target_name: str,
                           target_box: docking.DockingBox,
                           panel: dict[str, str] = OFFTARGET_PANEL) -> pd.DataFrame:
    all_results = []
    for kinase, pdb_id in panel.items():
        receptor_path = fetch_pdb_structure(pdb_id)
        # NOTE: each off-target needs its own docking box, ideally centered on
        # its ATP pocket (derive from co-crystal ligand centroid, as with the
        # primary target) -- `target_box` here is a placeholder reused across
        # the panel and MUST be replaced per-kinase before real runs.
        _, summary = docking.dock_library(
            compound_csv=shortlist_csv, receptor_pdb=receptor_path,
            box=target_box, target_name=kinase, n_replicates=3,
        )
        summary["kinase"] = kinase
        all_results.append(summary)

    panel_df = pd.concat(all_results, ignore_index=True)
    return panel_df


def compute_selectivity_scores(target_summary: pd.DataFrame,
                                panel_df: pd.DataFrame) -> pd.DataFrame:
    """Selectivity score = best off-target affinity - primary target affinity
    (positive = more selective; the compound binds its intended target more
    tightly than any off-target in the panel)."""
    best_offtarget = (panel_df.groupby(["compound_index", "smiles"])
                      ["affinity_best"].min()
                      .reset_index().rename(columns={"affinity_best": "best_offtarget_affinity"}))
    merged = target_summary.merge(best_offtarget, on=["compound_index", "smiles"])
    merged["selectivity_score"] = (
        merged["best_offtarget_affinity"] - merged["affinity_best"])
    return merged.sort_values("selectivity_score", ascending=False)


if __name__ == "__main__":
    # Placeholder box -- replace with the real primary-target box used in 05.
    box = docking.DockingBox(center_x=0.0, center_y=0.0, center_z=0.0)

    panel_results = run_selectivity_panel(
        shortlist_csv="data/shortlist_top_hits.csv",
        target_name="CHK1", target_box=box,
    )
    panel_results.to_csv("data/selectivity_panel_results.csv", index=False)

    target_summary = pd.read_csv("data/docking_results_summary.csv")
    selectivity = compute_selectivity_scores(target_summary, panel_results)
    selectivity.to_csv("data/selectivity_scores.csv", index=False)
    print(selectivity.head(10))
