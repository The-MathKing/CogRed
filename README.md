# CogRed — Cognate Redocking Validation Pipeline

Pipeline code and analysis scripts supporting:

> **Cognate redocking as a quality-control step for cryo-EM kinase receptor structures: a study of ATR/CHK1 and CDK-activating kinase**
>
> Ye E, Padarthi A, Kim S, Elie-Dit-Cosaque A, Jo E (2026)
>
> *Journal of Computer-Aided Molecular Design* (submitted)

Raw data and results are available in the companion repository:
[ethanyyy123/research-molecule-editing](https://github.com/ethanyyy123/research-molecule-editing)

---

## Installation

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

### External dependencies (not pip-installable)

| Tool | Version | Purpose |
|------|---------|---------|
| [AutoDock Vina](https://github.com/ccsb-scripps/AutoDock-Vina) | 1.2.7 (Python bindings) | Molecular docking |
| [Open Babel](https://openbabel.org/) | 3.1.1 | Receptor/ligand preparation, format conversion |
| [Meeko](https://github.com/forlilab/Meeko) | 2025 | Flexible-receptor PDBQT preparation |
| [spyrmsd](https://github.com/RMeli/spyrmsd) | 0.6.2 | Symmetry-corrected RMSD |

### Receptor structures

Receptor structures are **not** included in this repository (large, regenerable). Retrieve from the RCSB PDB AWS Open Data mirror:

```
s3://pdbsnapshots (snapshot 2026-01-01)
```

PDB IDs used: `9L40`, `9L4B`, `2YM8`, and 14 CAK structures (`8P6V`–`8P7L`; see Table 2 in the manuscript).

---

## Pipeline Scripts

The `scripts/` directory contains the full computational pipeline. Scripts are numbered and designed to run in order:

| # | Script | Manuscript § | Description |
|---|--------|-------------|-------------|
| 01 | `01_fragile_site_mapping.py` | SI | CFS genomic features (FRA3B/FRA16D/FRA7H) |
| 02 | `02_seed_compound_library.py` | 2.2 | Seed compound SMILES retrieval from PubChem/ChEMBL |
| 03 | `03_derivative_generation.py` | 2.3 | BRICS fragmentation + recombination (capped at 5,000 candidates) |
| 04 | `04_developability_prediction.py` | 2.6 | Lipinski/Veber/Ghose/Egan drug-likeness + PAINS/Brenk filtering |
| 05 | `05_docking_pipeline.py` | 2.4 | Vina batch docking (exhaustiveness 32); 3 seeded replicates; redocking validation required before scoring |
| 06 | `06_md_simulation_setup.py` | 2.7 | OpenMM MD setup (OpenFF ligand parametrization) |
| 07 | `07_mmpbsa_binding_energy.py` | 2.8 | MM-GBSA relative binding energy estimation |
| 08 | `08_kinase_selectivity_profiling.py` | 2.9 | Cross-kinase selectivity profiling |
| 09 | `09_pipeline_validation.py` | 2.5 | Scaffold-split retrieval benchmark with leakage exclusion |
| 10 | `10_integrated_ranking.py` | 2.10 | Multi-objective ranking with Dirichlet weight-sensitivity analysis |

> **Note:** Scripts 06–10 require additional compute resources (GPU for MD, AmberTools for MM-GBSA) and were not executed in the pilot study described in the manuscript. They represent the full-scope pipeline design.

---

## Data Directory

The `data/pilot_run_v1/` directory contains all executed analyses. Each subdirectory includes the Python script that generated its results:

```
data/pilot_run_v1/
├── seed_compounds.csv              # 5 seed compounds (3 co-crystallized + 2 known inhibitors)
├── derivative_library.csv          # 786 BRICS derivatives passing pre-filters
├── developability_predictions.csv  # Drug-likeness scores for all 791 compounds
├── docking_candidate_set.csv       # 40-compound pilot set selected for docking
├── redocking_validation_results.json  # Pass/fail for 9L40, 9L4B, 2YM8
├── atr_docking_results_raw.csv     # 40 compounds × 3 seeds vs. validated 9L4B
├── atr_docking_results_summary.csv # Mean ± SD Vina scores
├── box_centers.json                # Docking box centers (ligand centroids)
│
├── ligands/                        # Extracted co-crystallized ligands (.sdf)
├── validation_analysis/            # Symmetry-corrected RMSD, wwPDB reports, redocked poses
├── cak_survey/                     # 14-structure CDK7/cyclin H/MAT1 survey (Table 2)
├── flexible_redocking/             # Flexible-sidechain redocking test (§3.4)
├── rank_concordance/               # 9L40 vs. 9L4B rank comparison (§3.5)
└── casf_2016_control/              # 20-complex CASF-2016 positive control (§3.3)
```

---

## Reproducing Key Results

### 1. Redocking validation (Table 1)
```bash
cd data/pilot_run_v1
python run_redocking_validation.py    # Initial validation
python run_redocking_replicates.py    # 3-seed replicates
cd validation_analysis
python compute_symmrmsd.py            # Symmetry-corrected RMSD
```

### 2. CAK multi-structure survey (Table 2)
```bash
cd data/pilot_run_v1/cak_survey
python run_batch.py                   # All 14 CAK structures
```

### 3. CASF-2016 positive control (Table 3)
```bash
cd data/pilot_run_v1/casf_2016_control
python run_casf_batch.py              # 20 CASF-2016 complexes
```

### 4. Flexible-sidechain analysis (§3.4)
```bash
cd data/pilot_run_v1/flexible_redocking
python run_flexible_redocking_one_seed.py   # 3 seeds, 5 flexible residues
cd ../casf_2016_control
python run_rigid_control_matched.py         # Matched rigid control (30 Å box)
```

### 5. Rank concordance / error propagation (§3.5)
```bash
cd data/pilot_run_v1/rank_concordance
python run_library_docking_9l40_3seed.py    # 40 compounds vs. excluded 9L40
```

---

## Key Results

| Analysis | Result | Manuscript |
|----------|--------|-----------|
| 9L40 redocking (2.87 Å resolution) | **Fail** — RMSD 2.53 ± 0.01 Å | Table 1 |
| 9L4B redocking (3.20 Å resolution) | **Pass** — RMSD 0.61 ± 0.02 Å | Table 1 |
| 2YM8 redocking (2.07 Å, X-ray) | **Pass** — RMSD 0.76 ± 0.01 Å | Table 1 |
| CAK survey (14 cryo-EM structures) | 13/14 fail (7.1% pass rate) | Table 2 |
| CASF-2016 control (20 complexes) | 8/20 pass (40%; 95% CI 19–64%) | Table 3 |
| Flexible redocking of 9L40 | 3.00 ± 0.92 Å (improved vs. matched rigid control 4.50 Å, still fails) | §3.4 |
| Rank concordance 9L40 vs. 9L4B | Spearman ρ = 0.21; 5/10 shared top-10 | §3.5, Table 4 |

---

## Citation

If you use this pipeline or data, please cite the manuscript (reference to be updated upon publication).

## License

This project is provided for academic and research use. See the manuscript for full methodological details.
