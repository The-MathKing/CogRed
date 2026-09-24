# Pilot Run v1 — Executed Results

All files in this directory are actual outputs from the cognate redocking
validation pipeline. See the repository README for reproduction instructions.

## Subdirectories

- `ligands/` — Co-crystallized ligands extracted from PDB structures
- `validation_analysis/` — Symmetry-corrected RMSD, wwPDB reports, redocked poses
- `cak_survey/` — 14-structure CDK7/cyclin H/MAT1 survey (manuscript Table 2)
- `flexible_redocking/` — Flexible-sidechain redocking test (manuscript §3.4)
- `rank_concordance/` — 9L40 vs. 9L4B rank comparison (manuscript §3.5)
- `casf_2016_control/` — 20-complex CASF-2016 positive control (manuscript §3.3)

## Top-level files

- `seed_compounds.csv` — 5 seed compounds with SMILES
- `derivative_library.csv` — 786 BRICS derivatives
- `developability_predictions.csv` — Drug-likeness scores for all compounds
- `docking_candidate_set.csv` — 40-compound pilot set
- `redocking_validation_results.json` — Primary pass/fail results
- `atr_docking_results_raw.csv` — Library docking raw output
- `atr_docking_results_summary.csv` — Library docking summary

Receptor structures are not included (regenerable from RCSB PDB AWS mirror,
snapshot 2026-01-01). See repository README for retrieval instructions.
