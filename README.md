# In Silico Optimization of ATR/CHK1 Modulators for Fragile-Site Stability

Open-science computational pipeline supporting the manuscript *"Structure-
Based Optimization and Selectivity Profiling of ATR/CHK1 Pathway Modulators
to Suppress Replication Stress at Common Fragile Sites."*

Start here:

- **[`docs/01_critique_and_rigor_elevation.md`](docs/01_critique_and_rigor_elevation.md)**
  — PI-level critique of the project's scientific framing, the specific gaps
  that would draw reviewer pushback (especially the Cri du Chat mechanistic
  claim and the inhibitor-vs-modulator pharmacology question), and the
  additions that elevate this from a screening exercise to a publishable
  study. **Read this first** — it explains the design choices in every script
  below.
- **[`docs/02_manuscript_outline.md`](docs/02_manuscript_outline.md)** —
  granular IMRaD outline, section-by-section, mapped to the scripts that
  produce each result.
- **[`docs/03_journal_targeting_strategy.md`](docs/03_journal_targeting_strategy.md)**
  — five target journals (JCIM, PLOS Comp Biol, Scientific Reports, Frontiers
  in Pharmacology, ACS Omega) with audience, formatting, and the specific
  technical benchmarks each expects.

## Pipeline (`scripts/`, run in order)

| # | Script | Manuscript section | What it does |
|---|---|---|---|
| 01 | `01_fragile_site_mapping.py` | 2.10 / 3.7 | Genomic coordinates + AT-content/flexibility features for FRA3B/FRA16D/FRA7H |
| 02 | `02_seed_compound_library.py` | 2.2 | Reproducible SMILES pull for known ATR/CHK1 modulators (PubChem) + broader ChEMBL actives |
| 03 | `03_derivative_generation.py` | 2.3 | BRICS-based derivative generation + Lipinski/Veber/PAINS/Brenk pre-filter |
| 04 | `04_admet_filtering.py` | 2.6 | Descriptor-based ADMET triage (SwissADME-equivalent rules) + hooks for merging SwissADME/ProTox-II exports |
| 05 | `05_docking_pipeline.py` | 2.4 | AutoDock Vina batch docking wrapper (ligand/receptor prep via Open Babel, triplicate seeded runs) |
| 06 | `06_md_simulation_setup.py` | 2.7 | OpenMM MD (OpenFF ligand parametrization) for shortlisted hits |
| 07 | `07_mmpbsa_binding_energy.py` | 2.8 | MM-GBSA/MM-PBSA via ParmEd → MMPBSA.py on MD trajectories |
| 08 | `08_kinase_selectivity_profiling.py` | 2.9 | Off-target kinome panel docking + selectivity scoring |
| 09 | `09_pipeline_validation.py` | 2.5 / 3.2 | Actives/decoys enrichment (ROC-AUC, EF1%/EF5%) to validate the docking protocol itself |

## Setup

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

External CLI tools (not pip-installable, install separately): Open Babel
(`obabel`), AutoDock Vina (`vina`), AmberTools (`MMPBSA.py`, for script 07
only — used purely as an MM-PBSA post-processing tool, the MD engine
throughout is OpenMM).

## Status / open parameters

See `docs/01_critique_and_rigor_elevation.md` §6 for the assumptions currently
baked into the scripts (seed compound strategy, OpenMM/MMPBSA.py choice,
compute-budget scoping for MD) — flag if any should change before running the
pipeline for real.
