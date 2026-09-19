<!--
DRAFT STATUS — READ BEFORE USING OR CIRCULATING THIS FILE

What is real in this draft:
  - Introduction: literature-grounded biology, written to the tiered
    disease-relevance framing established in docs/04_review_response.md.
  - Methods: describes the pipeline exactly as implemented in scripts/
    01-10 of this repository. It is a protocol description of code that
    exists and compiles, not a report of a completed run.
  - Discussion 4.2, 4.4, 4.5: these do not depend on numerical results and
    are written in full.
  - Reference list: real papers, cited from memory. Author/title/journal/
    year are given in good faith; volume, page, and DOI are deliberately
    omitted where not independently verified in this session and must be
    confirmed against the primary source before submission. Two DOIs
    (Vina, JCIM MD guidelines) were verified via WebSearch in this session
    and are marked as such.

What is NOT real and must not be mistaken for it:
  - No compound has been docked, no MD trajectory has been run, no MM-PBSA
    value, developability prediction, or genomic feature number in this
    file is a real computed output. This sandboxed session has no network
    egress to PubChem, ChEMBL, UCSC, or RCSB PDB, and no AutoDock Vina,
    OpenMM, or AmberTools installation -- so no stage of the pipeline
    could actually be executed here.
  - Every place a number, table, or figure would go, this draft uses an
    explicit [PENDING: ...] marker naming the exact script that produces
    it. Do not fill those with invented numbers; run the script.

Do not submit, cite, or share this file's quantitative content as if it
were data. It is a manuscript skeleton with the prose that does not
depend on data already written.
-->

# Leakage-Controlled Benchmarking and Structure-Based Prioritization of Candidate ATR and CHK1 Ligands for Replication-Stress Research

**Authors:** [PENDING — author list and affiliations]
**Corresponding author:** [PENDING]
**Target journal:** *Scientific Reports* (title ≤20 words: 14 words here;
abstract ≤200 words; main text ≤4,500 words excluding Abstract, Methods,
References, and figure legends — see `docs/03_journal_targeting_strategy.md`)

---

## Abstract

Common fragile sites (CFSs) such as FRA3B (*FHIT*) and FRA16D (*WWOX*) are
chromosomal regions prone to breakage under replication stress, and the
ATR-CHK1 checkpoint axis contributes to their stability by restraining origin
firing and protecting stalled replication forks. Almost all existing chemical
probes of this pathway are ATP-competitive kinase inhibitors developed for
oncology synthetic-lethality strategies, and no openly available,
leakage-controlled computational framework exists for prioritizing candidate
ATR/CHK1 ligands independent of that inhibitor-centric chemical space. Here we
present an open-source pipeline that assembles a seed library from
PubChem/ChEMBL, generates derivatives by BRICS recombination, and evaluates
retrieval performance under Bemis–Murcko scaffold-split cross-validation with
an explicit seed-to-evaluation leakage audit, benchmarking a multi-stage
screen (docking, in silico developability prediction, molecular dynamics,
MM-PBSA relative energetics, and predicted cross-kinase selectivity) against
baseline docking alone. **[PENDING: retrieval performance (ROC-AUC, PR-AUC,
enrichment factors) and the outcome of the pipeline-vs-baseline comparison —
`scripts/09_pipeline_validation.py`].** We report **[PENDING: N]** prioritized
candidates with predicted binding, selectivity, and developability profiles.
These are offered as testable hypotheses for functional evaluation, not as
demonstrated modulators: computational docking cannot establish whether a
ligand activates, inhibits, or partially modulates ATR/CHK1 signaling. Code,
data, and the pre-specified analysis plan are openly available.

*(Word count target: ≤200 for Scientific Reports. Current draft, including
placeholder bracket text: 194 words — re-check after replacing placeholders,
since real numbers are typically shorter than the bracket text describing
them.)*

---

## 1. Introduction

### 1.1 Replication stress and common fragile sites

Common fragile sites are chromosomal loci that show a reproducible,
sequence-associated tendency to form gaps and breaks when cells are exposed to
mild replication stress, classically induced experimentally with the DNA
polymerase-α inhibitor aphidicolin¹. FRA3B, spanning the *FHIT* gene on
3p14.2, and FRA16D, spanning *WWOX* on 16q23.1, are among the most extensively
characterized CFSs and are recurrently deleted in a wide range of cancers².
The determinants of fragility that are well established in the literature
include late or delayed replication timing relative to the rest of the
genome, a paucity of licensed replication origins across large genomic
distances, transcription–replication conflicts (many CFS-associated genes,
including *FHIT* and *WWOX*, are unusually large), and the resulting
requirement that a small number of forks complete an unusually long stretch of
synthesis before mitosis³. A further class of *candidate* sequence-level
contributors — AT-rich flexibility peaks and regions predicted to form
non-B DNA secondary structure — has also been described at CFS loci, but
should be treated as compound-specific and contributing rather than as a
universal defining property of all fragile sites, since more recent work has
emphasized that fragility emerges from the interaction of several of these
features rather than from any single sequence code³.

### 1.2 The ATR–CHK1 axis

ATR (ataxia telangiectasia and Rad3-related kinase) is the principal sensor of
replication stress and is recruited to stalled forks and resected DNA ends
coated with RPA, where it is activated through protein-mediated mechanisms —
principally the ATRIP–TOPBP1 axis and, independently, ETAA1 — rather than by
occupation of its own ATP-binding pocket by a small molecule⁴,⁵. Once active,
ATR phosphorylates and activates CHK1, which in turn restrains further origin
firing, stabilizes replication forks against nucleolytic collapse, and
enforces intra-S and G2/M checkpoint arrest until replication is complete⁶.
Genetic and pharmacological attenuation of ATR or CHK1 activity increases CFS
expression under replication stress, consistent with this axis being
protective at these loci⁷.

It is essential to keep two distinct pharmacologies separate. Physiological
ATR activation is a regulated, protein-complex-dependent process. Nearly all
existing ATR/CHK1-directed chemical matter — berzosertib, ceralasertib,
elimusertib, prexasertib, and related compounds — instead works by
**ATP-competitive inhibition** of kinase catalytic activity, developed to
remove checkpoint function in cancer cells that already carry high endogenous
replication stress (a synthetic-lethality strategy)⁸. There is no established
continuum in which higher or lower binding affinity of an ATP-pocket ligand
maps onto a graded range from "inhibition" through "normal signaling" to
"activation." This distinction is the central constraint on what a
docking-based computational study of this target class can and cannot claim,
and it is stated here rather than left implicit.

### 1.3 A chemical-matter gap, not a modulator-discovery claim

Because the existing chemical matter for ATR/CHK1 is inhibitor-centric and
oncology-directed, there is no systematically characterized, openly
reproducible library of candidate ligands assembled and evaluated
independently of that inhibitor screening context, and no published
leakage-controlled benchmark of standard virtual-screening practice
(docking, molecular dynamics, MM-PBSA, cross-target selectivity scoring) for
this specific target pair. This is the gap addressed here: a reproducible
**ligand prioritization** framework, not a modulator-discovery claim.

### 1.4 What computation can and cannot establish here

ATP-pocket occupancy, as estimated by docking, molecular dynamics, and
MM-PBSA, does not by itself determine whether a ligand activates or inhibits
ATR or CHK1, because — as noted in §1.2 — physiological activation of these
kinases is governed by upstream protein interactions largely absent from a
static or short-timescale structural model. Accordingly, this study
prioritizes candidate ligands and defines the experiment that would be needed
to resolve their functional pharmacology; it does not claim to have resolved
it computationally.

### 1.5 Disease relevance, stated at the tier it is supported

**Tier 1 (established):** replication stress destabilizes common fragile
sites, and ATR/CHK1 signaling contributes to their protection³,⁷. **Tier 2
(hypothesis motivated by, but not tested in, this work):** replication fork
instability more broadly is mechanistically implicated in the formation of
non-recurrent genomic rearrangements through fork stalling and template
switching (FoSTeS) and microhomology-mediated break-induced replication
(MMBIR)⁹,¹⁰; by extension, it is plausible — but not established here — that
ATR/CHK1 pathway status could influence the probability of such rearrangements
at vulnerable loci more generally. **Tier 3 (distant motivation):** this
could eventually bear on constitutional genomic disorders arising from
non-recurrent deletions¹¹. This manuscript does not claim to establish a
causal link to any specific constitutional syndrome, and no such syndrome is
required to motivate the computational work that follows.

### 1.6 Objectives

We built and benchmarked an open pipeline (Figure 1) that: (i) assembles a
seed compound set for ATR and CHK1 from public databases; (ii) generates
structural derivatives; (iii) evaluates retrieval performance of the docking
protocol under scaffold-split cross-validation with an explicit leakage audit,
comparing baseline docking against a multi-stage screen; (iv) predicts
developability and toxicity liabilities in silico; (v) assesses pose stability
by molecular dynamics and estimates relative binding energetics by MM-PBSA for
a shortlist; (vi) predicts cross-kinase selectivity against a panel of
structurally related kinases, normalized within each receptor and benchmarked
against reference compounds with published experimental kinome profiles; and
(vii) produces an integrated ranking with a pre-specified weighting scheme and
a Monte Carlo weight-sensitivity analysis. All code is openly available (see
Data and Code Availability).

---

## 2. Methods

*Word count for this section is excluded from the Scientific Reports 4,500-word
main-text limit.*

### 2.1 Structural preparation

ATR and CHK1 receptor structures are obtained from RCSB PDB. For each
structure the following are recorded prior to use, per the reproducibility
checklist in `docs/01_critique_and_rigor_elevation.md` §5: PDB identifier,
resolution, species and construct boundaries, kinase activation state,
identity of any co-crystallized ligand, treatment of cofactors and metal ions,
retention or removal of crystallographic waters, resolution of alternate
conformations, modeling of missing residues, and the protonation/tautomer
assignment protocol and target pH. **[PENDING: specific PDB IDs selected for
ATR and CHK1, with the justification for why each represents a biologically
relevant state — to be finalized before docking is run; see
`data/receptors/` once populated.]**

### 2.2 Seed compound library (`scripts/02_seed_compound_library.py`)

A named seed set of known ATR/CHK1-pathway ligands (berzosertib,
ceralasertib, elimusertib, gartisertib for ATR; prexasertib, SRA737,
rabusertib, MK-8776 for CHK1) is resolved to canonical SMILES via the PubChem
PUG-REST API and standardized in RDKit (salt stripping, charge
neutralization, canonicalization). A broader bioactivity-defined seed set is
retrieved from the ChEMBL REST API for the same two targets at a stated
potency threshold (IC50 ≤ 1 µM in this draft; **[PENDING: confirm final
threshold and exact ChEMBL target IDs before running]**). Provenance
(source database, retrieval date, original and standardized SMILES) is
recorded for every seed compound, which is required for the leakage audit in
§2.5.

### 2.3 Derivative generation (`scripts/03_derivative_generation.py`)

The seed set is fragmented with the BRICS algorithm¹² and recombined
(`BRICS.BRICSBuild`) to propose new structures. Candidates are pre-filtered
before any downstream compute: a relaxed Lipinski rule (≤1 violation of
molecular weight ≤500, cLogP ≤5, H-bond donors ≤5, H-bond acceptors ≤10)¹³,
Veber's rule (rotatable bonds ≤10, TPSA ≤140 Å²)¹⁴, and removal of any
structure matching the RDKit PAINS¹⁵ or Brenk structural-alert catalogs. This
is standard cheminformatics practice and is not presented as a methodological
contribution in its own right (see §2.5 for what is).

### 2.4 Molecular docking (`scripts/05_docking_pipeline.py`)

Ligands are 3D-embedded and protonated at pH 7.4 with Open Babel and converted
to PDBQT; receptors are stripped of waters/heteroatoms and converted
similarly. Docking is performed with AutoDock Vina¹⁶. Each ligand is docked in
three independent replicates with different random seeds (seeds and the exact
Vina version are logged per run into the results table, not retyped by hand)
at a fixed exhaustiveness setting. The search box is derived from the
co-crystallized ligand centroid of the chosen receptor structure for each
target. **[PENDING: final exhaustiveness value, box dimensions per receptor,
and Vina version string — populated automatically by the script's
`get_vina_version()` call at run time.]**

### 2.5 Evaluation protocol: scaffold-split benchmarking with a leakage audit
(`scripts/09_pipeline_validation.py`)

This is the methodological core of the study. Known actives are split into
train/evaluation partitions by Bemis–Murcko scaffold¹⁷ (`scaffold_split`),
so that no scaffold present in the seed/derivative-generation set can also
appear in the evaluation set — a standard random split would place close
analogs of seed compounds on both sides and inflate every retrieval metric
that follows. An explicit leakage audit (`audit_leakage`) reports, for every
evaluation compound, its nearest-neighbor Tanimoto similarity (radius-2
Morgan fingerprints) to the seed/derivative set and flags scaffold collisions
directly, rather than assuming the split succeeded. Property-matched decoys
are constructed in the DUD-E style (matched on molecular weight and cLogP,
Tanimoto similarity ≤0.35 to the corresponding active)¹⁸. Retrieval
performance — ROC-AUC, PR-AUC, and enrichment factor at the top 1% and 5% of
ranked compounds, each with bootstrap 95% confidence intervals — is computed
for baseline Vina docking alone, for consensus scoring, and for the complete
multi-stage pipeline, and the pipeline is compared against the baseline with a
paired bootstrap test on the difference in ROC-AUC (`paired_bootstrap_difference`).
All splits, decoy-generation parameters, and comparisons were fixed in code
before any result was generated. **A result in which the confidence interval
of the pipeline-minus-baseline difference includes zero — i.e., no
demonstrable advantage over baseline docking — is reported as the finding,
not treated as a reason to keep adjusting the pipeline** (see
`docs/04_review_response.md` §B for the rationale).

### 2.6 In silico developability and toxicity prediction
(`scripts/04_developability_prediction.py`)

Every retained candidate is scored against four descriptor-based drug-likeness
rules computed directly in RDKit — Lipinski¹³, Veber¹⁴, Ghose¹⁹, and Egan²⁰ —
yielding a continuous count (0–4 rules satisfied) alongside the PAINS/Brenk
structural-alert flag from §2.3. This reproduces what SwissADME reports as
rule-based filters without depending on an unofficial scraping workflow.
Compounds retained after this stage are additionally submitted to SwissADME²¹
and ProTox-II²² for endpoints RDKit cannot approximate (gastrointestinal
absorption, CYP inhibition profile, predicted LD50 and hepatotoxicity
probability); results are merged by SMILES (`merge_external_admet_results`)
with the access date and tool version recorded. All values from this section
are predictions, not measurements, and are reported as such throughout.

### 2.7 Molecular dynamics (`scripts/06_md_simulation_setup.py`)

For the shortlist surviving docking, consensus scoring, and developability
prediction, receptor–ligand complexes are solvated (TIP3P, 1.0 nm padding,
0.15 M NaCl, neutralized) and parametrized with Amber ff14SB (protein) and
OpenFF Sage 2.1.0 (ligand, via `openmmforcefields`). Each system is
energy-minimized, equilibrated under NVT (100 ps) and NPT (100 ps, Monte
Carlo barostat, 1 atm, 300 K), and carried into production with OpenMM²³. Per
the JCIM guidelines for reporting molecular dynamics simulations²⁴, **each
system is run in at least three replicates started from different initial
velocities**, random seeds are disclosed, and replicate-to-replicate
statistical variance is reported alongside every summary statistic rather
than a single trajectory's value. **[PENDING: final production length per
replicate — 100 ns assumed as the default in the script; state explicitly if
reduced for compute reasons, per `docs/01_critique_and_rigor_elevation.md`
§6.4.]** Analyses include protein and ligand RMSD, RMSF, protein–ligand
contact persistence, and hydrogen-bond occupancy; an RMSD plateau is reported
as conformational stability and is not interpreted as evidence of binding
strength.

### 2.8 MM-PBSA relative energetic estimation
(`scripts/07_mmpbsa_binding_energy.py`)

OpenMM systems and trajectories are converted to Amber prmtop/inpcrd format
with ParmEd and post-processed with `MMPBSA.py`²⁵ (single-trajectory MM-GBSA,
igb=5, 0.15 M salt) purely as an end-state free-energy estimator; the
molecular dynamics engine throughout remains OpenMM. Entropy is neglected by
default, which is stated explicitly rather than left implicit. Estimates are
recomputed over three disjoint windows of each trajectory
(`assess_window_sensitivity`) to assess sensitivity to the choice of
analysis window; a result that varies substantially across windows indicates
the trajectory is not converged for this purpose and is reported as such. All
values are reported as **MM-PBSA estimated relative energetics** for ranking
purposes and never as a binding free energy (ΔG) or an experimental affinity.

### 2.9 Predicted cross-kinase selectivity analysis
(`scripts/08_kinase_selectivity_profiling.py`)

Shortlisted candidates are cross-docked, under the identical protocol
described in §2.4, into a panel of structurally related kinases chosen for
mechanistic relevance: the PIKK family (mTOR, DNA-PKcs, ATM, which share
ATR's ATP-pocket fold), cell-cycle kinases signaling in parallel to or
downstream of CHK1 (CDK1, CDK2, WEE1, PLK1), and a structurally distant
negative control (PKA). Because Vina scores are not calibrated across
different receptors, raw affinities are never subtracted directly between
targets. Instead, a shared property-matched background compound set is docked
into every receptor in the panel, and each candidate's score is expressed as
a z-score and percentile against that receptor's own background distribution
(`normalize_within_receptor`). A selectivity margin is then computed as the
candidate's standardized score against its intended target minus its best
standardized score against any off-target kinase. The metric itself is
validated, not assumed to be meaningful: reference compounds with published
experimental kinome-wide selectivity profiles (staurosporine as a broadly
promiscuous control; berzosertib as ATR-selective; prexasertib as
CHK1/CHK2-selective)²⁶,²⁷ are run through the identical protocol, and the
Spearman correlation between predicted and experimental rank order is reported
(`benchmark_against_reference_profiles`). **A weak or absent correlation is
reported as a limitation on how the selectivity heatmap may be interpreted,
not omitted.**

### 2.10 Integrated ranking (`scripts/10_integrated_ranking.py`)

Five objectives — predicted target engagement, predicted selectivity margin,
conformational stability (inverse RMSD), predicted developability, and
synthetic accessibility — are min–max normalized and combined into a
desirability score using pre-specified weights (engagement 0.35, selectivity
0.30, stability 0.15, developability 0.12, synthetic accessibility 0.08; full
rationale in the script docstring), fixed before any ranking was computed. To
guard against the ranking being an artifact of this particular weighting, a
Monte Carlo sensitivity analysis samples 5,000 alternative weight vectors from
a Dirichlet distribution concentrated around the pre-specified weights and
reports, for every candidate, the fraction of alternative weightings under
which it remains in the top five (`weight_sensitivity_analysis`). Candidates
reaching the top five under a small minority of sampled weightings are
reported as weight-sensitive rather than presented as robust hits. The
Spearman correlation matrix between the five objectives is also reported,
since docking score, MM-PBSA estimate, and MD stability share force-field and
pose assumptions and should not be read as independent confirmations of one
result.

### 2.11 Genomic feature analysis of common fragile sites
(`scripts/01_fragile_site_mapping.py`; Supporting Information)

Gene coordinates for *FHIT* (FRA3B), *WWOX* (FRA16D), and *CNTNAP2* (FRA7H,
comparison) are taken from the hg38 RefSeq track, and genomic sequence is
retrieved via the UCSC REST API to compute gene span, AT content, and GC skew.
This analysis is reported in Supporting Information, not in the main Results,
because these genomic features do not inform compound selection or docking in
this pipeline (see `docs/04_review_response.md` §A for why that separation is
deliberate). Its role is to specify the loci and expected direction of effect
for the translational experiment proposed in §4.5 of the Discussion.

### 2.12 Reproducibility

Software versions, hardware, all random seeds, and the pre-specified analysis
plan (splits, decoy parameters, weighting scheme) are recorded in the
repository accompanying this manuscript and archived at a versioned release
tag (Data and Code Availability).

---

## 3. Results

*Ordering follows `docs/02_manuscript_outline.md` Revision 2: the pipeline
must earn trust via benchmarking (§3.1–3.2) before its candidates are
presented (§3.3 onward).*

### 3.1 Scaffold-split retrieval of known actives

**[PENDING — `scripts/09_pipeline_validation.py`]** ROC-AUC, PR-AUC, and
enrichment factors at 1%/5%, with bootstrap 95% CIs, for the full pipeline on
scaffold-split held-out actives vs. property-matched decoys. Leakage audit
results (scaffold collisions, nearest-neighbor Tanimoto distribution) reported
alongside, not deferred to Supporting Information.

*Table 1 [PENDING]: Retrieval performance metrics with 95% CIs.*
*Figure 3 (panel A–B) [PENDING]: ROC and PR curves.*

### 3.2 Pipeline vs. baseline docking

**[PENDING — `scripts/09_pipeline_validation.py`, `compare_methods` /
`paired_bootstrap_difference`]** Per-stage contribution (baseline Vina →
consensus scoring → full pipeline) and the paired-bootstrap comparison against
baseline. *This result is reported whichever direction it falls, per the
pre-specified analysis plan (§2.5); a null result here is expected to be
plausible and will be presented as the finding, not treated as pipeline
failure.*

*Figure 3 (panel C–E) [PENDING]: baseline vs. consensus vs. full-pipeline
ROC/PR, and the difference distribution from the paired bootstrap.*

### 3.3 Chemical space of the derivative library

**[PENDING — `scripts/03_derivative_generation.py`]** Whether BRICS-based
derivative generation expanded accessible chemical space relative to the seed
set while retaining ATR/CHK1-like physicochemical character (molecular
weight, cLogP, TPSA distributions).

*Figure 2 [PENDING]: chemical-space projection (e.g., PCA over descriptors),
seed vs. derivative compounds.*

### 3.4 Docking and consensus scoring of the candidate set

**[PENDING — `scripts/05_docking_pipeline.py`]** Ranked docking results for
the derivative library against ATR and CHK1, with replicate variance across
the three seeded runs per ligand.

### 3.5 Predicted developability outcomes

**[PENDING — `scripts/04_developability_prediction.py`]** Number of
candidates retained after descriptor-based triage and (where available)
SwissADME/ProTox-II merge; property distributions before and after; reasons
for exclusion.

### 3.6 Conformational robustness (molecular dynamics)

**[PENDING — `scripts/06_md_simulation_setup.py`]** RMSD/RMSF and replicate
variance across ≥3 independent MD runs per shortlisted compound. Reported as
stability, not as evidence of binding strength.

### 3.7 MM-PBSA relative energetics

**[PENDING — `scripts/07_mmpbsa_binding_energy.py`]** Estimated relative
energetics with uncertainty and trajectory-window sensitivity
(`assess_window_sensitivity`). Axis/column labels retain "estimated relative
energetics," never ΔG.

*Figure 5 [PENDING]: MD RMSD/RMSF + MM-PBSA estimated relative energetics with
replicate variance shown explicitly.*

### 3.8 Predicted cross-kinase selectivity

**[PENDING — `scripts/08_kinase_selectivity_profiling.py`]** Selectivity
margins (within-receptor z-scores) across the panel, and — reported whether
favorable or not — the Spearman correlation between predicted and
experimental rank order for the reference-compound benchmark
(`benchmark_against_reference_profiles`).

*Figure 6 [PENDING]: predicted selectivity heatmap, with the reference-compound
benchmark panel shown alongside rather than presented separately.*

### 3.9 Integrated ranking and its stability

**[PENDING — `scripts/10_integrated_ranking.py`]** Final candidate ranking
under the pre-specified weights, with each candidate's top-5 stability
frequency under the Monte Carlo weight-sensitivity analysis, and the
Spearman correlation matrix between the five ranking objectives.

*Table 2 [PENDING]: top candidates with all objective scores, desirability
score, and top-5 stability frequency.*

---

## 4. Discussion

### 4.1 Comparison with known ATR/CHK1 pharmacology

**[PENDING — depends on §3.4 results.]** Quantitative comparison of
prioritized candidates against berzosertib/ceralasertib/prexasertib/SRA737 as
reference points, not as structures the pipeline is expected to rediscover.

### 4.2 The engagement-range hypothesis (does not depend on results)

The candidate set, once ranked, is expected to span a range of predicted
target engagement rather than clustering at a single extreme. This provides a
framework for a future functional experiment testing whether pathway output
(CHK1 phosphorylation, checkpoint signaling amplitude) varies non-monotonically
with target engagement — i.e., whether intermediate engagement is associated
with sustained, lower-amplitude signaling relative to both no engagement and
maximal engagement. This is stated as an explicitly untested hypothesis
motivating the translational path in §4.5, not as a conclusion this
computational work has reached. Binding affinity, by itself, is not equivalent
to kinase activation, pathway output, fork stability, or any cellular
phenotype.

### 4.3 Predicted selectivity liabilities

**[PENDING — depends on §3.8 results.]** Whether observed selectivity
liabilities plausibly arise from structural similarity of ATP pockets across
the PIKK family, which would itself be a mechanistic finding rather than only
a practical caveat.

### 4.4 Limitations (does not depend on results)

- **Chemical-space bias.** The seed library is drawn from ChEMBL/PubChem
  bioactivity data for ATR/CHK1, which is itself dominated by ATP-competitive
  inhibitor chemotypes; derivative generation inherits this bias.
- **Residual leakage risk.** Scaffold splitting and the leakage audit (§2.5)
  reduce but cannot fully eliminate similarity between seed and evaluation
  compounds, particularly for privileged kinase-inhibitor scaffolds (e.g.,
  aminopyrazines, pyrimidines) that recur across many kinase programs.
- **Non-equivalence of docking scores across targets.** Addressed by
  within-receptor normalization (§2.9), but the normalization itself depends
  on the chosen background compound set and is not a substitute for
  experimental selectivity data.
- **Predicted-ADMET uncertainty.** Descriptor-rule and QSAR-based
  developability/toxicity predictions (§2.6) are screening heuristics; a
  predicted absence of liability is not evidence of safety.
- **MM-PBSA model dependence.** Entropy is neglected by default; dielectric
  and dielectric-boundary assumptions, and sensitivity to sampling window, are
  reported explicitly (§2.8) but remain approximations.
- **Protein-state and regulatory-context uncertainty.** ATR/CHK1 function
  depends on upstream protein complexes (ATRIP, TOPBP1, ETAA1) absent from the
  structural models used here.
- **No functional directionality.** This is the limitation that governs every
  other claim in this manuscript: computational binding prediction does not
  establish activation versus inhibition, and no result in this study should
  be read as having done so.

### 4.5 Translational path

The path from this computational work to a disease-relevant claim is staged
and specific, not a single jump from ranking to phenotype: (1) a biochemical
ATR or CHK1 kinase assay on the top-ranked candidates from §3.9; (2) a
cellular checkpoint-signaling readout (phospho-CHK1 Ser345, phospho-ATR
Thr1989) for any compound active in (1); (3) a replication-fork phenotype
assay (DNA fiber analysis) for any compound active in (2); and (4) a
CFS-specific readout — gap-and-break frequency at FRA3B and FRA16D under
aphidicolin challenge, compared against replication-timing-matched non-fragile
control loci identified from the Supporting Information genomic analysis
(§2.11) — for any compound active in (3). This is the point at which the
genomic feature analysis does real work in this manuscript: it specifies which
loci and controls the terminal experiment should use, rather than serving as
decorative context for the chemistry.

---

## 5. Conclusion

**[PENDING — final wording depends on results, but the claim boundary below
must not be exceeded regardless of outcome.]**

Drafted floor for this section, to be filled in rather than expanded beyond:
*"We report a reproducible, leakage-controlled computational framework that
prioritizes candidate ATR and CHK1 ligands, benchmarks its own retrieval
performance against baseline docking, and yields [a candidate set / a
benchmarking result showing no demonstrable advantage over baseline docking —
delete as applicable] together with a specified experimental test of the
engagement-range hypothesis."* Do not write "we developed modulators that
suppress replication stress," or any equivalent claim of demonstrated
function — see §1.4 and §4.4.

---

## Data and Code Availability

All pipeline code (`scripts/01`–`10`), the pre-specified analysis plan, and
this manuscript source are available at the repository accompanying this
submission: [PENDING — repository URL / archived release DOI, e.g. via
Zenodo]. Large derived data (MD trajectories, downloaded receptor structures)
are excluded from version control per `.gitignore` and will be deposited at
[PENDING] with the archived release.

## Author Contributions

[PENDING]

## Competing Interests

[PENDING]

## Funding

[PENDING]

## Acknowledgements

[PENDING]

---

## References

*Cited from memory; author, approximate year, and venue are given in good
faith. Volume, page numbers, and DOIs are omitted except where noted as
verified in this session, and every entry must be checked against the primary
source before submission.*

1. Glover TW, Berger C, Coyle J, Echo B. DNA polymerase α inhibition by
   aphidicolin induces gaps and breaks at common fragile sites in human
   chromosomes. *Human Genetics*, 1984.
2. Bignell GR, et al. Signatures of mutation and selection in the cancer
   genome / fragile-site tumor suppressor studies on *FHIT* and *WWOX*
   (representative citation; verify specific primary reference).
3. Debatisse M, Le Tallec B, Letessier A, Dutrillaux B, Brison O. Common
   fragile sites: mechanisms of instability revisited. *Trends in Genetics*,
   2012.
4. Kumagai A, Lee J, Yoo HY, Dunphy WG. TopBP1 activates the ATR-ATRIP
   complex. *Cell*, 2006.
5. Bass TE, Luzwick JW, Kavanaugh G, et al. ETAA1 acts at stalled replication
   forks to maintain genome integrity. *Nature Cell Biology*, 2016.
6. Saldivar JC, Cortez D, Cimprich KA. The essential kinase ATR: ensuring
   faithful duplication of a challenging genome. *Nature Reviews Molecular
   Cell Biology*, 2017.
7. Casper AM, Nghiem P, Arlt MF, Glover TW. ATR regulates fragile site
   stability. *Cell*, 2002.
8. Fokas E, et al. Targeting ATR in vivo using the novel inhibitor
   berzosertib / representative ATR-inhibitor pharmacology review (verify
   specific primary reference for the pharmacology claim in §1.2).
9. Lee JA, Carvalho CM, Lupski JR. A DNA replication mechanism for generating
   nonrecurrent rearrangements associated with genomic disorders. *Cell*,
   2007.
10. Hastings PJ, Ira G, Lupski JR. A microhomology-mediated break-induced
    replication model for the origin of human copy number variation. *PLoS
    Genetics*, 2009.
11. Niebuhr E. The Cri du Chat syndrome: epidemiology, cytogenetics, and
    clinical features. *Human Genetics*, 1978. *(Cited only as background on
    the Tier-3 disease class in §1.5, per the scoping in
    docs/04_review_response.md; not used as a claim of causal mechanism.)*
12. Degen J, Wegscheid-Gerlach C, Zaliani A, Rarey M. On the art of compiling
    and using 'drug-like' chemical fragment spaces. *ChemMedChem*, 2008.
    (BRICS.)
13. Lipinski CA, Lombardo F, Dominy BW, Feeney PJ. Experimental and
    computational approaches to estimate solubility and permeability in drug
    discovery and development settings. *Advanced Drug Delivery Reviews*,
    1997/2001 update.
14. Veber DF, et al. Molecular properties that influence the oral
    bioavailability of drug candidates. *Journal of Medicinal Chemistry*,
    2002.
15. Baell JB, Holloway GA. New substructure filters for removal of pan assay
    interference compounds (PAINS). *Journal of Medicinal Chemistry*, 2010.
16. Trott O, Olson AJ. AutoDock Vina: improving the speed and accuracy of
    docking with a new scoring function, efficient optimization, and
    multithreading. *Journal of Computational Chemistry*, 2010. DOI:
    10.1002/jcc.21334 *(verified)*.
17. Bemis GW, Murcko MA. The properties of known drugs. 1. Molecular
    frameworks. *Journal of Medicinal Chemistry*, 1996.
18. Mysinger MM, Carchia M, Irwin JJ, Shoichet BK. Directory of Useful
    Decoys, Enhanced (DUD-E): better ligands and decoys for better
    benchmarking. *Journal of Medicinal Chemistry*, 2012.
19. Ghose AK, Viswanadhan VN, Wendoloski JJ. A knowledge-based approach in
    designing combinatorial or medicinal chemistry libraries for drug
    discovery. *Journal of Combinatorial Chemistry*, 1999.
20. Egan WJ, Merz KM, Baldwin JJ. Prediction of drug absorption using
    multivariate statistics. *Journal of Medicinal Chemistry*, 2000.
21. Daina A, Michielin O, Zoete V. SwissADME: a free web tool to evaluate
    pharmacokinetics, druglikeness and medicinal chemistry friendliness of
    small molecules. *Scientific Reports*, 2017.
22. Banerjee P, Eckert AO, Schrey AK, Preissner R. ProTox-II: a webserver for
    the prediction of toxicity of chemicals. *Nucleic Acids Research*, 2018.
23. Eastman P, et al. OpenMM 7: Rapid development of high performance
    algorithms for molecular dynamics. *PLoS Computational Biology*, 2017.
24. Guidelines for Reporting Molecular Dynamics Simulations in JCIM
    Publications (editorial). *Journal of Chemical Information and
    Modeling*, 2023. DOI: 10.1021/acs.jcim.3c00599 *(verified)*.
25. Miller BR 3rd, McGee TD Jr, Swails JM, Homeyer N, Gohlke H, Roitberg AE.
    MMPBSA.py: an efficient program for end-state free energy calculations.
    *Journal of Chemical Theory and Computation*, 2012.
26. Davis MI, Hunt JP, Herrgard S, et al. Comprehensive analysis of kinase
    inhibitor selectivity. *Nature Biotechnology*, 2011.
27. Klaeger S, Heinzlmeir S, Wilhelm M, et al. The target landscape of
    clinical kinase drugs. *Science*, 2017.
