# Candidate ATR/CHK1 Ligand Prioritization for Replication-Stress Research

Open-science computational pipeline supporting the manuscript *"Leakage-
Controlled Benchmarking and Structure-Based Prioritization of Candidate ATR
and CHK1 Ligands for Replication-Stress Research."*

> **Evidentiary boundary — read before using or citing anything here.**
> This pipeline supports claims about **predicted** binding, relative ranking,
> retrieval performance under scaffold-split evaluation, conformational
> stability, estimated relative energetics, predicted selectivity, and
> predicted developability. It **cannot** establish whether a ligand activates,
> inhibits, or partially modulates ATR/CHK1, nor any statement about pathway
> output, fork stability, or cellular or genomic phenotype. ATR activation is
> protein-mediated (ATRIP/TOPBP1/ETAA1) and is not determined by ATP-pocket
> occupancy. Outputs are **candidates and hypotheses for experimental
> testing**, not demonstrated modulators.

Start here:

- **[`manuscript/submission/paper.pdf`](manuscript/submission/paper.pdf)** —
  a real, fully executed, submittable pilot paper (LaTeX source in the same
  folder). This "sandboxed" environment actually has pip access to PyPI
  (real `vina`, `openmm`, `openbabel-wheel`, `spyrmsd` install from there),
  read access to the RCSB PDB's AWS Open Data mirror (`s3://pdbsnapshots`,
  reachable even though `rcsb.org` is not — including its wwPDB validation
  reports), and apt access (a full LaTeX + PyMOL toolchain installs
  cleanly) — so this paper's central finding is real, not a placeholder:
  **a redocking validation control reproducibly passed one 2025 ATR
  structure (9L4B, 3.20 Å) and failed another (9L40, 2.87 Å — the
  *higher*-resolution of the two) solved in the same study**, a finding
  that survives symmetry-corrected RMSD, ruling out the adjacent allosteric
  site as a confound, and re-validation at the exhaustiveness actually used
  for screening. A 40-compound BRICS-derived pipeline sanity-check screen
  was then docked against the validated structure. As of round 4, this
  n=1 finding was extended to a real 17-structure survey (14 additional
  cryo-EM structures of an unrelated kinase, CDK7/cyclin H/MAT1): **82%
  fail redocking, and wwPDB Q-score/residue inclusion do not predict which
  ones will** (Q-score actually trends the wrong way) — a negative result
  for the metric, reported as found. Flexible-sidechain docking (5
  active-site residues) does not rescue 9L40's failing pose either, and
  docking the same pilot library against 9L40 instead of the validated
  9L4B produces a substantially different, less trustworthy ranking (only
  5/10 shared top compounds) — a concrete demonstration of what the
  validation step prevents. Scope is deliberately narrower than the full
  pipeline below (no MD — no GPU available — no kinome panel, and the
  17-structure survey is a convenience sample, not a systematic one); see
  `data/pilot_run_v1/README.md` for the real data behind every number in
  it, `docs/06_review_response_round3.md` for what a third external
  review tried to break and what survived, and
  `docs/07_review_response_round4.md` for the round-4 extension (survey,
  flexible docking, error-propagation demo) and what did *not* come out
  the way it was hoped to.
- **[`manuscript/draft_v1.md`](manuscript/draft_v1.md)** — the earlier,
  fully-`[PENDING]` draft against the *original*, larger-scope pipeline
  below, written when PubChem/ChEMBL/UCSC/RCSB and Vina/OpenMM/AmberTools
  really were unreachable in that session. Superseded as the "submittable"
  artifact by the real pilot paper above, but still the right reference for
  the full multi-target pipeline's intended design once more compute/access
  is available.
- **[`docs/05_review_response_round2.md`](docs/05_review_response_round2.md)**
  — response to a second external review, of the manuscript draft itself.
  Several structural-biology and benchmarking-literature claims in that review
  were independently verified before acting on them (one — a specific analog
  count in a cited leakage audit — matched exactly). Real gaps closed as a
  result: a docking-accuracy validation control, an actual (not just reported)
  leakage-exclusion filter, a corrected MM-GBSA/MM-PBSA naming inconsistency,
  and a corrected PKA fold classification.
- **[`docs/04_review_response.md`](docs/04_review_response.md)** — response to
  external review (Round 1): what was accepted, what was contested, and every
  change it forced. **Read this first** — it explains why the project's central
  claim was lowered and why several scripts were rewritten.
- **[`docs/02_manuscript_outline.md`](docs/02_manuscript_outline.md)** —
  current IMRaD outline (Revision 2), restructured so the pipeline must earn
  trust via benchmarking before its candidates are presented.
- **[`docs/03_journal_targeting_strategy.md`](docs/03_journal_targeting_strategy.md)**
  — journal strategy (Revision 2). Scientific Reports is now the primary
  target; JCIM was demoted for a verified scope reason.
- **[`docs/01_critique_and_rigor_elevation.md`](docs/01_critique_and_rigor_elevation.md)**
  — original rigor audit. Sections 2 and 4 are **superseded**; retained for the
  revision record.

## Pipeline (`scripts/`, run in order)

| # | Script | Methods § | What it does |
|---|---|---|---|
| 01 | `01_fragile_site_mapping.py` | SI | CFS genomic features (FRA3B/FRA16D/FRA7H). Supporting Information — defines the terminal experimental prediction, does **not** drive compound selection |
| 02 | `02_seed_compound_library.py` | 2.2 | Reproducible SMILES retrieval (PubChem/ChEMBL); provenance recorded for the leakage audit |
| 03 | `03_derivative_generation.py` | 2.3 | BRICS derivative generation + property pre-filters. Standard cheminformatics, **not claimed as novelty** |
| 04 | `04_developability_prediction.py` | 2.6 | In silico developability/toxicity **prediction** (renamed from "ADMET triage" — nothing here is measured) |
| 05 | `05_docking_pipeline.py` | 2.4 | Vina batch docking (exhaustiveness 32 default); ≥3 seeded runs per ligand; **redocking validation control** (RMSD ≤2.0 Å pre-specified pass threshold) required before any candidate is scored |
| 06 | `06_md_simulation_setup.py` | 2.7 | OpenMM MD (OpenFF ligand parametrization) |
| 07 | `07_mmpbsa_binding_energy.py` | 2.8 | **MM-GBSA** relative energetic estimation (fixed naming — the method run here is GB, not PB; `MMPBSA.py` remains the correct tool name) + trajectory-window sensitivity. Never reported as ΔG |
| 08 | `08_kinase_selectivity_profiling.py` | 2.9 | **Predicted** cross-kinase selectivity, z-normalized *within each receptor*; PKA correctly framed as a fold comparator (CAMK/AGC), not a distant control; reference-set benchmark now requires ≥20 compounds before being treated as validation |
| 09 | `09_pipeline_validation.py` | 2.5 | Scaffold-split retrieval with a pre-specified Tanimoto ≥0.4 exclusion that actually removes leaked evaluation compounds (not just reports them); baseline-vs-pipeline comparison with bootstrap CIs; dual decoy-standard comparison |
| 10 | `10_integrated_ranking.py` | 2.10 | Pre-specified desirability weights **plus** Dirichlet Monte-Carlo weight-sensitivity analysis and objective-correlation reporting |

### Pre-specified analysis plan

Splits, metrics, and desirability weights are fixed in code (`09_*`, `10_*`)
before results are seen. **If the benchmark shows the multi-stage pipeline does
not outperform baseline docking, that is the reported result** — see
`docs/04_review_response.md` §B. Both scripts print an explicit notice when
their outputs fail to support the more favorable interpretation.

## Setup

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

External CLI tools (not pip-installable): Open Babel (`obabel`), AutoDock Vina
(`vina`), AmberTools (`MMPBSA.py`, script 07 only — post-processing only; the
MD engine throughout is OpenMM).

## Open items

- **Highest-value next step:** a wet-lab collaborator for a biochemical
  ATR/CHK1 kinase assay on the top candidates. It is the only thing that
  restores a functional-pharmacology claim and reopens JCIM/PLOS Comp Biol.
- Assumptions currently baked into the scripts (seed strategy, OpenMM +
  MMPBSA.py choice, MD compute scoping) are listed in
  `docs/01_critique_and_rigor_elevation.md` §6.
