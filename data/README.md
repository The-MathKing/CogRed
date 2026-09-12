# Data directory

This directory is populated by running the scripts in `scripts/`, in order.
Nothing here is committed as tracked binary data except small CSVs generated
by the pipeline itself — large trajectories (`*.dcd`), checkpoints (`*.chk`),
and downloaded receptor structures should stay out of version control (see
`.gitignore`) and instead be archived on Zenodo/OSF alongside the paper, with
the DOI referenced in Supporting Information.

Expected contents once the pipeline has been run:

| File | Produced by |
|---|---|
| `fragile_site_features.csv` | `scripts/01_fragile_site_mapping.py` |
| `seed_named_compounds.csv`, `chembl_actives.csv` | `scripts/02_seed_compound_library.py` |
| `derivative_library.csv` | `scripts/03_derivative_generation.py` |
| `admet_triage_results.csv`, `admet_ready_compounds.csv` | `scripts/04_admet_filtering.py` |
| `docking_results_raw.csv`, `docking_results_summary.csv` | `scripts/05_docking_pipeline.py` |
| `md/<compound>/production.dcd` (untracked) | `scripts/06_md_simulation_setup.py` |
| `mmpbsa_results.csv` | `scripts/07_mmpbsa_binding_energy.py` |
| `selectivity_panel_results.csv`, `selectivity_scores.csv` | `scripts/08_kinase_selectivity_profiling.py` |
| `pipeline_validation_metrics.csv` | `scripts/09_pipeline_validation.py` |
| `receptors/*.pdb` (untracked) | downloaded from RCSB PDB |
