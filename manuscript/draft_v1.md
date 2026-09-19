<!--
DRAFT STATUS — READ BEFORE USING OR CIRCULATING THIS FILE (Revision 2,
post external review round 2 — see docs/05_review_response_round2.md)

What is real in this draft:
  - Introduction, Methods, and the results-independent parts of Discussion:
    literature-grounded prose describing the pipeline exactly as implemented
    in scripts/01-10 of this repository. It is a protocol description of
    code that exists and compiles, not a report of a completed run.
  - Reference list: real papers. Several were independently verified via
    WebSearch this session (marked "verified" — DOI/resolution/finding
    checked against RCSB PDB, PubMed, or the publisher directly) after a
    second external review flagged specific, checkable claims (structural
    biology, benchmarking literature) rather than accepting them on say-so.
    Entries not marked "verified" are cited from memory in good faith and
    must still be checked against the primary source before submission.
  - Citation style: author-year (e.g., "Glover et al., 1984") rather than
    numbered superscripts, deliberately, so that references can be added,
    removed, or corrected during revision without renumbering every in-text
    mark. Convert to Scientific Reports' numbered format at final
    formatting, after the reference list is frozen.

What is NOT real and must not be mistaken for it:
  - No compound has been docked, no MD trajectory has been run, no MM-GBSA
    value, developability prediction, or genomic feature number in this
    file is a real computed output. This sandboxed session has no network
    egress to PubChem, ChEMBL, UCSC, or RCSB PDB (only a small allowlist,
    e.g. PyPI, is reachable), and no AutoDock Vina, OpenMM, or AmberTools
    installation -- so no stage of the pipeline could actually be executed
    here. "Code that compiles" is not evidence the pipeline runs end to end;
    treat every stage as undemonstrated until it has actually been run once.
  - Every place a number, table, or figure would go, this draft uses an
    explicit [PENDING: ...] marker naming the exact script that produces
    it. Do not fill those with invented numbers; run the script.
  - Target venue is, and has been throughout this repository's history,
    Scientific Reports (see docs/03_journal_targeting_strategy.md). No
    instruction in this project's history has ever named a different venue.

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
ATR–CHK1 checkpoint axis contributes to their stability by restraining origin
firing and protecting stalled replication forks. Nearly all existing chemical
probes of this pathway are ATP-competitive kinase inhibitors — though
berzosertib was recently shown to also engage an allosteric site at the
ATR–ATR dimer interface (Wang et al., 2025) — and, to our knowledge, no
leakage-controlled benchmark of virtual-screening practice has been reported
for this target pair. We present an open-source pipeline, building on
established scaffold-split and decoy-bias-control methods, that assembles a
seed library from PubChem/ChEMBL, generates derivatives by BRICS
recombination, and evaluates retrieval under scaffold-split cross-validation
with a pre-specified similarity-exclusion threshold, benchmarking a
multi-stage screen (redocking-validated docking, in silico developability
prediction, molecular dynamics, MM-GBSA relative energetics, predicted
cross-kinase selectivity) against baseline docking alone. **[PENDING:
retrieval performance and outcome of the pipeline-vs-baseline comparison —
`scripts/09_pipeline_validation.py`].** We report **[PENDING: N]** prioritized
candidates with predicted binding, selectivity, and developability profiles,
offered as testable hypotheses, not demonstrated modulators: computational
docking cannot establish whether a ligand activates, inhibits, or modulates
ATR/CHK1 signaling. Code, data, and a timestamped analysis plan are openly
available.

*(Word count target: ≤200 for Scientific Reports; recheck once placeholders
are replaced with real numbers, which are typically shorter than the
bracket text describing them.)*

---

## 1. Introduction

### 1.1 Replication stress and common fragile sites

Common fragile sites are chromosomal loci that show a reproducible,
sequence-associated tendency to form gaps and breaks when cells are exposed to
mild replication stress, classically induced experimentally with the DNA
polymerase-α inhibitor aphidicolin (Glover et al., 1984). FRA3B, spanning the
*FHIT* gene on 3p14.2, and FRA16D, spanning *WWOX* on 16q23.1, are among the
most extensively characterized CFSs and are recurrently deleted in a wide
range of cancers (Huebner & Croce, 2001). The determinants of fragility that
are well established in the literature include late or delayed replication
timing relative to the rest of the genome, a paucity of licensed replication
origins across large genomic distances, transcription–replication conflicts
(many CFS-associated genes, including *FHIT* and *WWOX*, are unusually
large), and the resulting requirement that a small number of forks complete
an unusually long stretch of synthesis before mitosis (Debatisse et al.,
2012). A further class of *candidate* sequence-level contributors — AT-rich
flexibility peaks and regions predicted to form non-B DNA secondary structure
— has also been described at CFS loci, but should be treated as
compound-specific and contributing rather than as a universal defining
property of all fragile sites, since more recent work has emphasized that
fragility emerges from the interaction of several of these features rather
than from any single sequence code (Debatisse et al., 2012).

### 1.2 The ATR–CHK1 axis

ATR (ataxia telangiectasia and Rad3-related kinase) is the principal sensor of
replication stress and is recruited to stalled forks and resected DNA ends
coated with RPA, where it is activated through protein-mediated mechanisms —
principally the ATRIP–TOPBP1 axis (Kumagai et al., 2006) and, independently,
ETAA1 (Bass et al., 2016) — rather than by occupation of its own ATP-binding
pocket by a small molecule. Once active, ATR phosphorylates and activates
CHK1, which in turn restrains further origin firing, stabilizes replication
forks against nucleolytic collapse, and enforces intra-S and G2/M checkpoint
arrest until replication is complete (Saldivar, Cortez & Cimprich, 2017).
Genetic and pharmacological attenuation of ATR or CHK1 activity increases CFS
expression under replication stress, consistent with this axis being
protective at these loci (Casper et al., 2002).

It is essential to keep two distinct pharmacologies separate. Physiological
ATR activation is a regulated, protein-complex-dependent process. Existing
ATR/CHK1-directed chemical matter — berzosertib, ceralasertib, elimusertib,
prexasertib, and related compounds — was developed to remove checkpoint
function in cancer cells that already carry high endogenous replication
stress (a synthetic-lethality strategy) via **ATP-competitive inhibition** of
kinase catalytic activity (Fokas et al., 2012). This picture is, however,
more structurally complex than a single ATP-pocket story: near-atomic-
resolution cryo-EM of the full human ATR–ATRIP complex bound to berzosertib
shows one complex engaging **four** berzosertib molecules — two in the
canonical ATP-competitive active sites, and two at a distinct allosteric
pocket formed at the ATR–ATR dimer interface (Wang et al., 2025, *verified*
— DOI 10.1016/j.scib.2025.05.009). This is directly relevant to the framing
of this study: it demonstrates that at least one clinical-stage ATR ligand
already engages a non-ATP-competitive site, which is a more plausible
structural entry point for a genuinely modulatory (as opposed to purely
inhibitory) pharmacology than the ATP pocket this pipeline currently docks
against. We do not pursue that allosteric site computationally here — the
receptor structures, box definitions, and scoring protocol below are
ATP-pocket-focused throughout — but we flag it explicitly as the most
promising specific direction for extending this work (§4.3, Future Work),
rather than only in general terms. There remains no established continuum in
which higher or lower ATP-pocket binding affinity alone maps onto a graded
range from "inhibition" through "normal signaling" to "activation," and nothing
in the paragraph above should be read as claiming otherwise.

### 1.3 A chemical-matter gap, scoped to what a documented search supports

Because the existing ATP-pocket-targeted chemical matter for ATR/CHK1 is
inhibitor-centric and oncology-directed, there is, to our knowledge,
no openly reproducible, leakage-controlled benchmark of standard
virtual-screening practice (docking, molecular dynamics, MM-GBSA,
cross-target selectivity scoring) specifically for the ATR/CHK1 target pair.
**This is a narrower claim than an earlier draft of this manuscript made**,
and it needs to stay narrow: leakage-controlled virtual-screening
benchmarking is an active, well-developed subfield in its own right, and this
work is an instantiation of that methodology for one target pair, not a claim
to have invented the methodology. In particular, this pipeline's evaluation
protocol (§2.5) builds directly on: benchmark-bias diagnostics originally
raised for the DUD-E decoy standard (Sieg, Flachsenberg & Rarey, 2019; Chen et
al., 2019), bias-corrected benchmark construction such as LIT-PCBA
(Tran-Nguyen, Jacquemard & Rognan, 2020) and DeepCoy-generated decoys (Imrie et
al., 2021), a recent audit documenting that scaffold splitting alone does not
prevent analog-level leakage in LIT-PCBA itself — over 350 train/validation
analog pairs at ECFP4 Tanimoto ≥ 0.6, 323 in the ALDH1 target alone
(anonymous, 2025, *verified* — arXiv:2507.21404) — and fixed-budget,
multi-strategy virtual-screening benchmarking as recently instantiated for
other target classes (anonymous, 2026, *verified* — DOI
10.3390/ijms27167497). **[PENDING: the negative claim above — that no
leakage-controlled ATR/CHK1-specific benchmark has been published — is stated
without a documented systematic search; before submission, record the
databases, date, and query strings used to check this, or soften the claim
further.]** This is the gap addressed here: a reproducible **ligand
prioritization** framework for this specific target pair, not a claim to
have invented leakage-controlled benchmarking, and not a modulator-discovery
claim.

### 1.4 What computation can and cannot establish here

ATP-pocket occupancy, as estimated by docking, molecular dynamics, and
MM-GBSA, does not by itself determine whether a ligand activates or inhibits
ATR or CHK1, because — as noted in §1.2 — physiological activation of these
kinases is governed by upstream protein interactions largely absent from a
static or short-timescale structural model, and because (per §1.2) even the
existing structural picture of ATP-pocket engagement is more complex than a
single-site model. Accordingly, this study prioritizes candidate ligands and
defines the experiment that would be needed to resolve their functional
pharmacology; it does not claim to have resolved it computationally.

### 1.5 Disease relevance, stated at the tier it is supported

**Tier 1 (established):** replication stress destabilizes common fragile
sites, and ATR/CHK1 signaling contributes to their protection (Debatisse et
al., 2012; Casper et al., 2002). **Tier 2 (hypothesis motivated by, but not
tested in, this work):** replication fork instability more broadly is
mechanistically implicated in the formation of non-recurrent genomic
rearrangements through fork stalling and template switching (FoSTeS) and
microhomology-mediated break-induced replication (MMBIR) (Lee, Carvalho &
Lupski, 2007; Hastings, Ira & Lupski, 2009); by extension, it is plausible —
but not established here — that ATR/CHK1 pathway status could influence the
probability of such rearrangements at vulnerable loci more generally. This
manuscript does not extend the disease-relevance framing beyond this tier: no
specific constitutional syndrome is named or required to motivate the
computational work that follows, and no such claim should be added in later
revision without independent support for it.

### 1.6 Objectives

We built and benchmarked an open pipeline (Figure 1) that: (i) assembles a
seed compound set for ATR and CHK1 from public databases; (ii) generates
structural derivatives; (iii) evaluates retrieval performance of the docking
protocol under scaffold-split cross-validation with a pre-specified
similarity-exclusion threshold, comparing baseline docking against a
multi-stage screen; (iv) predicts developability and toxicity liabilities in
silico; (v) validates the docking protocol by redocking co-crystallized
ligands, then assesses pose stability by molecular dynamics and estimates
relative binding energetics by MM-GBSA for a shortlist; (vi) predicts
cross-kinase selectivity against a panel of structurally related kinases,
normalized within each receptor and benchmarked against reference compounds
with published experimental kinome profiles; and (vii) produces an integrated
ranking with a pre-specified weighting scheme, a uniform-weight and
Pareto-front comparison, and a Monte Carlo weight-sensitivity analysis. All
code is openly available (see Data and Code Availability).

---

## 2. Methods

*Word count for this section is excluded from the Scientific Reports 4,500-word
main-text limit.*

### 2.1 Structural preparation

ATR and CHK1 receptor structures are obtained from RCSB PDB, chosen for
resolution and construct suitability rather than by default to the
best-known structure. For ATR specifically, the original 2018 cryo-EM
structure of the full ATR–ATRIP complex (PDB 5YZ0, 4.7 Å overall, *verified*)
is **excluded by design** as a docking receptor: at that resolution, side-chain
positions in the ATP pocket are not reliably resolved, which would make any
docking result against it uninterpretable regardless of everything else in
this pipeline. Instead, the 2025 near-atomic-resolution structures of the
ATR kinase domain bound to clinical-stage inhibitors — PDB 9L40 (bound to
VE-822/berzosertib) and 9L4B (bound to RP-3500), both at approximately 3 Å,
*verified* (Wang et al., 2025) — are the working default for the ATP-pocket
docking described in this manuscript. CHK1 has numerous high-resolution
(sub-2.5 Å) ligand-bound crystal structures publicly available;
**[PENDING: specific CHK1 PDB ID, selected by the same resolution and
construct-completeness criterion applied to ATR above, and recorded here
before docking is run]**. For each structure used, the following are
recorded per the reproducibility checklist in
`docs/01_critique_and_rigor_elevation.md` §5: PDB identifier, resolution,
species and construct boundaries, kinase activation state, identity of any
co-crystallized ligand, treatment of cofactors and metal ions, retention or
removal of crystallographic waters, resolution of alternate conformations,
modeling of missing residues, and the protonation/tautomer assignment
protocol and target pH.

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
§2.5. The resulting library size at this and every subsequent stage is
**[PENDING]** and will be reported in a compound-flow diagram (Figure 4),
since retrieval metrics such as enrichment factor at 1% are only
interpretable once the evaluated library size is known.

### 2.3 Derivative generation (`scripts/03_derivative_generation.py`)

The seed set is fragmented with the BRICS algorithm (Degen et al., 2008) and
recombined (`BRICS.BRICSBuild`) to propose new structures. Candidates are
pre-filtered before any downstream compute: a relaxed Lipinski rule (≤1
violation of molecular weight ≤500, cLogP ≤5, H-bond donors ≤5, H-bond
acceptors ≤10; Lipinski et al., 1997), Veber's rule (rotatable bonds ≤10,
TPSA ≤140 Å²; Veber et al., 2002), and removal of any structure matching the
RDKit PAINS (Baell & Holloway, 2010) or Brenk structural-alert catalogs. This
is standard cheminformatics practice and is not presented as a methodological
contribution in its own right (see §2.5 for what is). Note that BRICS
recombination can produce products whose Bemis–Murcko scaffold differs from
their parent seed compound while remaining close chemical neighbors by
fingerprint similarity; this is a known limitation of pairing BRICS-based
generation with scaffold-only leakage control and is addressed procedurally,
not eliminated, by the similarity-based exclusion threshold in §2.5.

### 2.4 Molecular docking (`scripts/05_docking_pipeline.py`)

Ligands are 3D-embedded and protonated at pH 7.4 with Open Babel and converted
to PDBQT; receptors are stripped of waters/heteroatoms and converted
similarly. Docking is performed with AutoDock Vina (Trott & Olson, 2010,
*verified* — DOI 10.1002/jcc.21334) at exhaustiveness 32 by default (raised
from Vina's own default of 8, which is not considered adequate for a result
intended to be reported as a reproducible benchmark score). The search box for
each target is derived from its co-crystallized ligand centroid.

**Docking validation.** Before any candidate is scored, the protocol is
validated by redocking: the co-crystallized ligand of each receptor structure
is extracted, re-docked into its own receptor under the identical protocol,
and the RMSD between the redocked and native pose is reported against a
pre-specified pass threshold of ≤2.0 Å (heavy atoms). A receptor structure or
box definition that fails this control is not used for candidate scoring
without an explicit stated reason. This is reported separately from, and
should not be conflated with, the three independent seeded replicate runs per
ligand described next: **replicate agreement measures search reproducibility
(does Vina return the same answer on repeated runs), not pose accuracy (is
the answer correct)** — redocking RMSD measures the latter, and only the
latter licenses any claim about docking accuracy. Seeds and the exact Vina
version are logged per run into the results table, not retyped by hand.
**[PENDING: final box dimensions per receptor and Vina version string —
populated automatically by the script's `get_vina_version()` call at run
time.]**

### 2.5 Evaluation protocol: scaffold-split benchmarking with a leakage audit
and exclusion (`scripts/09_pipeline_validation.py`)

This is the methodological core of the study. Known actives are split into
train/evaluation partitions by Bemis–Murcko scaffold (Bemis & Murcko, 1996;
`scaffold_split`), so that no scaffold present in the seed/derivative-generation
set can also appear in the evaluation set — a standard random split would
place close analogs of seed compounds on both sides and inflate every
retrieval metric that follows. Because scaffold splitting alone has been shown
to leave substantial analog-level leakage in comparable benchmarks (over 350
train/validation analog pairs at ECFP4 Tanimoto ≥ 0.6 in LIT-PCBA; anonymous,
2025, *verified* — arXiv:2507.21404), the leakage audit here does not stop at
reporting: any evaluation compound whose nearest-neighbor Tanimoto similarity
(radius-2 Morgan fingerprints) to the seed/derivative set is **≥ 0.4** —
pre-specified before any result was generated — is **excluded from the
evaluation set**, not merely flagged, and the number excluded is reported.
Property-matched decoys are constructed in the DUD-E style (matched on
molecular weight and cLogP, Tanimoto similarity ≤0.35 to the corresponding
active; Mysinger et al., 2012). **This decoy standard has documented,
systematic biases** — property-matching alone can be learned and exploited
independent of true activity (Sieg, Flachsenberg & Rarey, 2019; Chen et al.,
2019) — so results are also reported, where feasible, against an
experimentally confirmed inactive set and/or DeepCoy-generated decoys (Imrie
et al., 2021) as a secondary check, and any large discrepancy between the two
decoy standards is reported rather than the more favorable number alone.
Retrieval performance — ROC-AUC, PR-AUC, and enrichment factor at the top 1%
and 5% of ranked compounds, each with bootstrap 95% confidence intervals — is
computed for baseline Vina docking alone, for consensus scoring, and for the
complete multi-stage pipeline, and the pipeline is compared against the
baseline with a paired bootstrap test on the difference in ROC-AUC
(`paired_bootstrap_difference`).

**Scope of what this benchmark can show.** Because the evaluation actives are
drawn from the same ChEMBL/PubChem bioactivity pull used to build the seed and
derivative library (§2.2), a strong retrieval result here demonstrates that
the pipeline recovers **ATP-competitive kinase-inhibitor chemotypes** — the
same chemical space the field already has — and does not by itself
demonstrate an ability to find novel chemotypes outside it. This is stated
explicitly rather than left for a reader to infer. A temporally held-out
active set (compounds first published after a stated cutoff date) would
address this and is noted as a planned extension, not yet implemented
**[PENDING]**.

All splits, exclusion thresholds, decoy-generation parameters, and comparisons
were fixed in code before any result was generated, and this analysis plan
will be timestamped (a Zenodo/OSF registration or a signed, dated git tag
preceding the first pipeline execution — **[PENDING]**, see §2.12) so that
the pre-specification claim is independently checkable rather than only
asserted in the text. **A result in which the confidence interval of the
pipeline-minus-baseline difference includes zero — i.e., no demonstrable
advantage over baseline docking — is reported as the finding, not treated as
a reason to keep adjusting the pipeline** (see `docs/04_review_response.md`
§B for the original rationale).

### 2.6 In silico developability and toxicity prediction
(`scripts/04_developability_prediction.py`)

Every retained candidate is scored against four descriptor-based drug-likeness
rules computed directly in RDKit — Lipinski (Lipinski et al., 1997), Veber
(Veber et al., 2002), Ghose (Ghose, Viswanadhan & Wendoloski, 1999), and Egan
(Egan, Merz & Baldwin, 2000) — yielding a discrete count (0–4 rules satisfied)
alongside the PAINS/Brenk structural-alert flag from §2.3. This reproduces
what SwissADME reports as rule-based filters without depending on an
unofficial scraping workflow. Compounds retained after this stage are
additionally submitted to SwissADME (Daina, Michielin & Zoete, 2017) and
ProTox-II (Banerjee et al., 2018) for endpoints RDKit cannot approximate
(gastrointestinal absorption, CYP inhibition profile, predicted LD50 and
hepatotoxicity probability); results are merged by SMILES
(`merge_external_admet_results`) with the access date and tool version
recorded. All values from this section are predictions, not measurements, and
are reported as such throughout.

### 2.7 Molecular dynamics (`scripts/06_md_simulation_setup.py`)

For the shortlist surviving docking (including the redocking validation
control), consensus scoring, and developability prediction, receptor–ligand
complexes are solvated (TIP3P, 1.0 nm padding, 0.15 M NaCl, neutralized) and
parametrized with Amber ff14SB (protein) and OpenFF Sage 2.1.0 (ligand, via
`openmmforcefields`). Each system is energy-minimized, equilibrated under NVT
(100 ps) and NPT (100 ps, Monte Carlo barostat, 1 atm, 300 K), and carried
into production with OpenMM (Eastman et al., 2017). Per the JCIM guidelines
for reporting molecular dynamics simulations (JCIM Editorial, 2023, *verified*
— DOI 10.1021/acs.jcim.3c00599), **each system is run in at least three
replicates started from different initial velocities**, random seeds are
disclosed, and replicate-to-replicate statistical variance is reported
alongside every summary statistic rather than a single trajectory's value.
**[PENDING: final production length per replicate — 100 ns assumed as the
default in the script; state explicitly if reduced for compute reasons, per
`docs/01_critique_and_rigor_elevation.md` §6.4.]** Analyses include protein
and ligand RMSD, RMSF, protein–ligand contact persistence, and hydrogen-bond
occupancy; an RMSD plateau is reported as conformational stability and is not
interpreted as evidence of binding strength.

### 2.8 MM-GBSA relative energetic estimation
(`scripts/07_mmpbsa_binding_energy.py`)

**Naming note:** the method implemented and reported here is MM-GBSA
(Generalized Born), not MM-PBSA (Poisson–Boltzmann) — an earlier draft of
this manuscript used the two terms inconsistently. `MMPBSA.py` (Miller et al.,
2012) is the correct name of the AmberTools *program* used for this step
regardless of which implicit-solvent model it is configured to run; the
*method* it is configured to run here, via the `&gb` input block with
`igb=5`, is MM-GBSA, and every other reference to this analysis in this
manuscript — abstract, this heading, §3.7, Figure 5 — uses that name. OpenMM
systems and trajectories are converted to Amber prmtop/inpcrd format with
ParmEd and post-processed with `MMPBSA.py` (0.15 M salt) purely as an
end-state free-energy estimator; the molecular dynamics engine throughout
remains OpenMM. If the ParmEd conversion from an OpenMM/OpenFF system proves
unreliable in practice — a known friction point for this specific toolchain
combination — `gmx_MMPBSA` is the fallback path and will be noted explicitly
if used instead. Entropy is neglected by default, which is stated explicitly
rather than left implicit. Estimates are recomputed over three disjoint
windows of each trajectory (`assess_window_sensitivity`) to assess sensitivity
to the choice of analysis window; a result that varies substantially across
windows indicates the trajectory is not converged for this purpose and is
reported as such. All values are reported as **MM-GBSA estimated relative
energetics** for ranking purposes and never as a binding free energy (ΔG) or
an experimental affinity.

### 2.9 Predicted cross-kinase selectivity analysis
(`scripts/08_kinase_selectivity_profiling.py`)

Shortlisted candidates are cross-docked, under the identical protocol
described in §2.4, into a panel of structurally related kinases chosen for
mechanistic relevance: the PIKK family (mTOR, DNA-PKcs, ATM, which share
ATR's ATP-pocket fold) and cell-cycle kinases signaling in parallel to or
downstream of CHK1 (CDK1, CDK2, WEE1, PLK1). PKA is included as a **canonical
AGC-group fold comparator**, not, as an earlier draft mischaracterized it, a
structurally distant negative control: CHK1 belongs to the CAMK group of the
human kinome and PKA is the AGC-group archetype of the same bilobal
eukaryotic protein kinase fold (Manning et al., 2002), so PKA is a near
neighbor of CHK1 by fold, not a distant one. It is retained specifically to
test whether selectivity liabilities are fold-driven rather than
target-specific, which is a different and useful question, but a genuinely
fold-distant ATP-binding-protein control is **not yet included in the panel**
and is flagged as a needed addition **[PENDING]**. Because Vina scores are not
calibrated across different receptors, raw affinities are never subtracted
directly between targets. Instead, a shared property-matched background
compound set is docked into every receptor in the panel, and each candidate's
score is expressed as a z-score and percentile against that receptor's own
background distribution (`normalize_within_receptor`). A selectivity margin
is then computed as the candidate's standardized score against its intended
target minus its best standardized score against any off-target kinase. The
metric's correspondence to reality is checked, not assumed: reference
compounds with published experimental kinome-wide selectivity profiles
(Davis et al., 2011; Klaeger et al., 2017) are run through the identical
protocol, and the Spearman correlation between predicted and experimental
rank order is reported (`benchmark_against_reference_profiles`). **A
reference set of three compounds (staurosporine, berzosertib, prexasertib), as
used in the initial pipeline design, is acknowledged as too small to produce
an interpretable correlation** — its confidence interval would span most of
the possible range — and is being expanded to 20–50 compounds drawn from
Davis et al. (2011) and Klaeger et al. (2017) before this validation is run
for real **[PENDING]**. Until that expansion, any correlation computed from
the smaller set is reported as illustrative only, not as validation, and the
selectivity heatmap is not presented as quantitatively interpretable on that
basis alone.

### 2.10 Integrated ranking (`scripts/10_integrated_ranking.py`)

Five objectives — predicted target engagement, predicted selectivity margin,
conformational stability (inverse RMSD), predicted developability, and
synthetic accessibility — are min–max normalized and combined into a
desirability score using pre-specified weights (engagement 0.35, selectivity
0.30, stability 0.15, developability 0.12, synthetic accessibility 0.08; full
rationale in the script docstring), fixed before any ranking was computed.
Docking score, MM-GBSA estimate, and MD stability share force-field and pose
assumptions and are not independent evidence; a Spearman correlation matrix
between all five objectives is computed, and where any pair correlates at ρ >
0.7, both the original weighted score and a redundancy-corrected variant
(the correlated pair's combined weight halved and redistributed) are reported
side by side, so a reader can see whether the ranking depends on treating
correlated objectives as independent. To guard against the ranking being an
artifact of the pre-specified weighting itself, three checks are reported
together rather than the first alone: (i) a Monte Carlo sensitivity analysis
sampling 5,000 alternative weight vectors from a Dirichlet distribution
concentrated around the pre-specified weights, reporting the fraction of
alternative weightings under which each candidate remains in the top five
(`weight_sensitivity_analysis`); (ii) a uniform-weight ranking (all five
objectives weighted equally) reported alongside for comparison; and (iii) the
Pareto-optimal (non-dominated) candidate set, which requires no weighting
scheme at all. Rank agreement across all three is reported; a candidate
robust only under the pre-specified weights and not under (ii) or (iii) is
reported as weight-dependent rather than presented as a robust hit. Synthetic
accessibility is currently a soft-weighted heuristic score rather than a hard
filter; adding a retrosynthesis-based hard filter (e.g., AiZynthFinder;
Genheden et al., 2020) so that unmakeable BRICS products cannot reach the top
ranks regardless of their other scores is a planned addition, not yet
implemented **[PENDING]**.

### 2.11 Genomic feature analysis of common fragile sites
(`scripts/01_fragile_site_mapping.py`; Supporting Information)

Gene coordinates for *FHIT* (FRA3B), *WWOX* (FRA16D), and *CNTNAP2* (FRA7H,
comparison) are taken from the hg38 RefSeq track, and genomic sequence is
retrieved via the UCSC REST API to compute gene span, AT content, and GC skew.
This analysis is reported in Supporting Information, not in the main Results,
because these genomic features do not inform compound selection or docking in
this pipeline (see `docs/04_review_response.md` §A for why that separation is
deliberate). **Its intended role — narrowed from an earlier draft — is to
provide descriptive genomic context for §1.1, not to identify
replication-timing-matched control loci by itself**: gene span, AT content,
and GC skew cannot substitute for actual replication-timing data (e.g.,
Repli-seq from ENCODE or the 4D Nucleome Consortium), which this analysis does
not currently incorporate. Identifying appropriate non-fragile control loci
for the translational assay proposed in §4.4 requires that data and is listed
there as a dependency still to be added, not as something this Supporting
Information analysis already delivers.

### 2.12 Reproducibility

Software versions, hardware, all random seeds, and the pre-specified analysis
plan (splits, similarity-exclusion threshold, decoy parameters, weighting
scheme) are recorded in the repository accompanying this manuscript.
**Pre-specification is claimed only where it is independently checkable**: the
analysis plan will be registered with a timestamp preceding the first
pipeline execution (a Zenodo or OSF deposit, or at minimum a signed, dated git
tag) **[PENDING]** — an in-text assertion of pre-specification without such a
record is not evidence of it, and this manuscript should not rely on the
assertion alone once that registration exists.

---

## 3. Results

*Ordering follows `docs/02_manuscript_outline.md` Revision 2: the pipeline
must earn trust via benchmarking (§3.1–3.2) before its candidates are
presented (§3.3 onward).*

### 3.1 Docking validation and scaffold-split retrieval of known actives

**[PENDING — `scripts/05_docking_pipeline.py` redocking control;
`scripts/09_pipeline_validation.py` for retrieval.]** Redocking RMSD for each
receptor structure against the ≤2.0 Å pass threshold (§2.4), reported before
any retrieval metric, since a receptor that fails this control invalidates
everything downstream. ROC-AUC, PR-AUC, and enrichment factors at 1%/5%, with
bootstrap 95% CIs, for the full pipeline on scaffold-split held-out actives
(post similarity-exclusion) vs. both DUD-E-style and DeepCoy-style decoys.
Leakage audit results (number of evaluation compounds excluded at the ≥0.4
Tanimoto threshold) reported alongside, not deferred to Supporting
Information.

*Table 1 [PENDING]: Retrieval performance metrics with 95% CIs, both decoy
standards.*
*Figure 3 (panel A–B) [PENDING]: ROC and PR curves.*
*Figure 4 [PENDING]: compound-flow diagram (seed → derivatives → docking
survivors → developability survivors → MD shortlist → final candidates),
with N reported at every stage.*

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
the three seeded runs per ligand, reported separately from the redocking
accuracy control (§3.1).

### 3.5 Predicted developability outcomes

**[PENDING — `scripts/04_developability_prediction.py`]** Number of
candidates retained after descriptor-based triage and (where available)
SwissADME/ProTox-II merge; property distributions before and after; reasons
for exclusion. See Figure 4 for stage-by-stage N.

### 3.6 Conformational robustness (molecular dynamics)

**[PENDING — `scripts/06_md_simulation_setup.py`]** RMSD/RMSF and replicate
variance across ≥3 independent MD runs per shortlisted compound. Reported as
stability, not as evidence of binding strength.

### 3.7 MM-GBSA relative energetics

**[PENDING — `scripts/07_mmpbsa_binding_energy.py`]** Estimated relative
energetics with uncertainty and trajectory-window sensitivity
(`assess_window_sensitivity`). Axis/column labels retain "estimated relative
energetics" and the method is labelled MM-GBSA throughout, never ΔG and
never MM-PBSA (see §2.8).

*Figure 5 [PENDING]: MD RMSD/RMSF + MM-GBSA estimated relative energetics with
replicate variance shown explicitly.*

### 3.8 Predicted cross-kinase selectivity

**[PENDING — `scripts/08_kinase_selectivity_profiling.py`]** Selectivity
margins (within-receptor z-scores) across the panel, with PKA's results
interpreted as a fold-similarity comparator rather than a negative control
(§2.9), and — reported whether favorable or not — the reference-compound
correlation, explicitly labelled illustrative rather than validating until
the reference set is expanded (§2.9).

*Figure 6 [PENDING]: predicted selectivity heatmap, with the reference-compound
comparison shown alongside rather than presented separately.*

### 3.9 Integrated ranking and its stability

**[PENDING — `scripts/10_integrated_ranking.py`]** Final candidate ranking
under the pre-specified weights, the uniform-weight ranking, and the
Pareto-optimal set, with rank agreement across all three reported explicitly;
each candidate's top-5 stability frequency under the Monte Carlo
weight-sensitivity analysis; and the Spearman correlation matrix between the
five ranking objectives, with the redundancy-corrected score reported
wherever any pair correlates at ρ > 0.7.

*Table 2 [PENDING]: top candidates with all objective scores, desirability
score (original and redundancy-corrected), uniform-weight rank, Pareto
membership, and top-5 stability frequency.*

---

## 4. Discussion

### 4.1 Comparison with known ATR/CHK1 pharmacology

**[PENDING — depends on §3.4 results.]** Quantitative comparison of
prioritized candidates against berzosertib/ceralasertib/prexasertib/SRA737 as
reference points, not as structures the pipeline is expected to rediscover.

### 4.2 Predicted selectivity liabilities

**[PENDING — depends on §3.8 results.]** Whether observed selectivity
liabilities plausibly arise from structural similarity of ATP pockets across
the PIKK family and, given §2.9, whether PKA-level scores track fold
similarity (CAMK/AGC) rather than target-specific promiscuity — itself a
mechanistic finding rather than only a practical caveat.

### 4.3 Limitations

- **Chemical-space bias.** The seed library is drawn from ChEMBL/PubChem
  bioactivity data for ATR/CHK1, which is itself dominated by ATP-competitive
  inhibitor chemotypes; derivative generation inherits this bias, and the
  retrieval benchmark in §3.1–3.2 measures recovery of that same chemical
  space (§2.5), not generalization beyond it.
- **Residual leakage risk.** The ≥0.4 Tanimoto exclusion (§2.5) reduces but
  does not eliminate similarity between seed and evaluation compounds,
  particularly for privileged kinase-inhibitor scaffolds (aminopyrazines,
  pyrimidines) that recur across many kinase programs, and BRICS-generated
  derivatives can defeat scaffold-only splitting by design (§2.3).
- **Docking-score accuracy vs. reproducibility.** Independent Vina seeds
  measure whether the search is reproducible, not whether the returned pose is
  correct; only the redocking control (§2.4) speaks to accuracy, and a
  receptor/box combination that fails it should not be trusted for anything
  downstream.
- **Decoy-standard bias.** DUD-E-style property-matched decoys have
  documented, exploitable biases (Sieg, Flachsenberg & Rarey, 2019; Chen et
  al., 2019); results are cross-checked against a second decoy standard
  (§2.5), but neither standard is bias-free.
- **Non-equivalence of docking scores across targets.** Addressed by
  within-receptor normalization (§2.9), but the normalization itself depends
  on the chosen background compound set and is not a substitute for
  experimental selectivity data; the current fold-comparator panel lacks a
  truly fold-distant negative control.
- **Predicted-ADMET uncertainty.** Descriptor-rule and QSAR-based
  developability/toxicity predictions (§2.6) are screening heuristics; a
  predicted absence of liability is not evidence of safety.
- **MM-GBSA model dependence.** Entropy is neglected by default; dielectric
  assumptions and sensitivity to sampling window are reported explicitly
  (§2.8) but remain approximations, and single-trajectory MM-GBSA is known to
  have weak, system-dependent ranking power outside congeneric series.
- **Protein-state and regulatory-context uncertainty.** ATR/CHK1 function
  depends on upstream protein complexes (ATRIP, TOPBP1, ETAA1) absent from the
  structural models used here, and the ATR structures used dock against the
  ATP-competitive site only, not the allosteric dimer-interface site
  identified by Wang et al. (2025) (§1.2).
- **Reference list under active verification.** Several citations in this
  draft are given from memory pending confirmation against the primary
  source (see header note); none should be treated as confirmed until that
  check is complete.
- **No functional directionality.** This is the limitation that governs every
  other claim in this manuscript: computational binding prediction does not
  establish activation versus inhibition, and no result in this study should
  be read as having done so.

### 4.4 Translational path

The path from this computational work to a disease-relevant claim is staged
and specific, not a single jump from ranking to phenotype: (1) a biochemical
ATR or CHK1 kinase assay on the top-ranked candidates from §3.9; (2) a
cellular checkpoint-signaling readout (phospho-CHK1 Ser345, phospho-ATR
Thr1989) for any compound active in (1); (3) a replication-fork phenotype
assay (DNA fiber analysis) for any compound active in (2); and (4) a
CFS-specific readout — gap-and-break frequency at FRA3B and FRA16D under
aphidicolin challenge — for any compound active in (3). Stage (4) requires
non-fragile control loci matched on replication timing, which — per the
narrowed scope of §2.11 — this manuscript's genomic analysis does not yet
supply; incorporating public Repli-seq data (ENCODE or 4D Nucleome) to
identify such loci is listed here as a concrete dependency for stage (4), not
assumed to already be in hand.

### 4.5 Future work

Two extensions are flagged here rather than folded into the main claims
above, because neither is currently supported by data or citation in this
manuscript: (i) computationally investigating the allosteric ATR dimer-
interface pocket identified by Wang et al. (2025) (§1.2) is, in our
assessment, the most promising specific route to a genuinely modulatory (as
opposed to purely inhibitory) pharmacology for this target, and is a more
concrete direction than a general appeal to "partial modulation"; (ii) whether
pathway output (CHK1 phosphorylation, checkpoint signaling amplitude) varies
non-monotonically with predicted target engagement is a plausible-sounding
but currently uncited hypothesis — no evidence or precedent for it is offered
in this manuscript — and it is recorded here as a question for a future
functional experiment, not as a framework this computational work has
established.

---

## 5. Conclusion

**[PENDING — final wording depends on results, but the claim boundary below
must not be exceeded regardless of outcome.]**

Drafted floor for this section, to be filled in rather than expanded beyond:
*"We report a reproducible, leakage-controlled computational framework that
prioritizes candidate ATR and CHK1 ligands for a target pair, benchmarks its
own retrieval performance against baseline docking and against known
decoy-bias concerns, and yields [a candidate set / a benchmarking result
showing no demonstrable advantage over baseline docking — delete as
applicable] together with specific, staged experimental next steps."* Do not
write "we developed modulators that suppress replication stress," or any
equivalent claim of demonstrated function — see §1.4 and §4.3.

---

## Data and Code Availability

All pipeline code (`scripts/01`–`10`), the pre-specified analysis plan, and
this manuscript source are available at the repository accompanying this
submission: [PENDING — repository URL / archived release DOI, e.g. via
Zenodo]. The analysis plan (splits, exclusion thresholds, decoy parameters,
weighting scheme) will be timestamped via that same archived release, dated
before the first pipeline execution, so pre-specification is independently
checkable rather than only asserted in the text (§2.5, §2.12) **[PENDING]**.
Large derived data (MD trajectories, downloaded receptor structures) are
excluded from version control per `.gitignore` and will be deposited at
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

*Author-year style for drafting (see header note); convert to Scientific
Reports' numbered format once the list is frozen. Entries marked "verified"
were checked via WebSearch this session (DOI, resolution, or reported finding
confirmed against RCSB PDB, PubMed, the publisher, or arXiv); all other
entries are cited from memory in good faith and must be checked against the
primary source before submission.*

Anonymous [authors not verified this session]. Data Leakage and Redundancy in
the LIT-PCBA Benchmark. *arXiv:2507.21404* (2025). **[verified: exists;
finding of 323 ALDH1 analog pairs at ECFP4 Tc≥0.6 confirmed via WebSearch;
author names not confirmed]**

Anonymous [authors not verified this session]. Beyond the Score: Fixed-Budget
Benchmarking of Virtual Screening Integration Strategies for Decision-Centric
Drug Discovery. *International Journal of Molecular Sciences* 27(16):7497
(2026). **[verified DOI: 10.3390/ijms27167497; author names not confirmed]**

Baell JB, Holloway GA. New substructure filters for removal of pan assay
interference compounds (PAINS). *Journal of Medicinal Chemistry* (2010).

Banerjee P, Eckert AO, Schrey AK, Preissner R. ProTox-II: a webserver for the
prediction of toxicity of chemicals. *Nucleic Acids Research* (2018).

Bass TE, Luzwick JW, Kavanaugh G, et al. ETAA1 acts at stalled replication
forks to maintain genome integrity. *Nature Cell Biology* (2016).

Bemis GW, Murcko MA. The properties of known drugs. 1. Molecular frameworks.
*Journal of Medicinal Chemistry* (1996).

Casper AM, Nghiem P, Arlt MF, Glover TW. ATR regulates fragile site stability.
*Cell* (2002).

Chen L, Cruz A, Ramsey S, et al. Hidden bias in the DUD-E dataset leads to
misleading performance of deep learning in structure-based virtual screening.
*PLOS ONE* 14(8):e0220113 (2019).

Daina A, Michielin O, Zoete V. SwissADME: a free web tool to evaluate
pharmacokinetics, druglikeness and medicinal chemistry friendliness of small
molecules. *Scientific Reports* (2017).

Davis MI, Hunt JP, Herrgard S, et al. Comprehensive analysis of kinase
inhibitor selectivity. *Nature Biotechnology* (2011).

Debatisse M, Le Tallec B, Letessier A, Dutrillaux B, Brison O. Common fragile
sites: mechanisms of instability revisited. *Trends in Genetics* (2012).

Degen J, Wegscheid-Gerlach C, Zaliani A, Rarey M. On the art of compiling and
using 'drug-like' chemical fragment spaces. *ChemMedChem* (2008).

Eastman P, Swails J, Chodera JD, et al. OpenMM 7: Rapid development of high
performance algorithms for molecular dynamics. *PLoS Computational Biology*
(2017).

Egan WJ, Merz KM, Baldwin JJ. Prediction of drug absorption using
multivariate statistics. *Journal of Medicinal Chemistry* (2000).

Fokas E, Prevo R, Pollard JR, et al. Targeting ATR in vivo using the novel
inhibitor VE-822 results in selective sensitization of pancreatic tumors to
radiation. *Cell Death & Disease* (2012).

Genheden S, Thakkar A, Chadimová V, et al. AiZynthFinder: a fast, robust and
flexible open-source software for retrosynthetic planning. *Journal of
Cheminformatics* (2020).

Ghose AK, Viswanadhan VN, Wendoloski JJ. A knowledge-based approach in
designing combinatorial or medicinal chemistry libraries for drug discovery.
*Journal of Combinatorial Chemistry* (1999).

Glover TW, Berger C, Coyle J, Echo B. DNA polymerase α inhibition by
aphidicolin induces gaps and breaks at common fragile sites in human
chromosomes. *Human Genetics* (1984).

Hastings PJ, Ira G, Lupski JR. A microhomology-mediated break-induced
replication model for the origin of human copy number variation. *PLoS
Genetics* (2009).

Huebner K, Croce CM. FRA3B and other common fragile sites: the weakest links.
*Nature Reviews Cancer* (2001).

Imrie F, Bradley AR, Deane CM. Generating property-matched decoy molecules
using deep learning (DeepCoy). *Bioinformatics* 37(15):2134–2141 (2021).

JCIM Editorial. Guidelines for Reporting Molecular Dynamics Simulations in
JCIM Publications. *Journal of Chemical Information and Modeling* (2023).
**[verified DOI: 10.1021/acs.jcim.3c00599]**

Klaeger S, Heinzlmeir S, Wilhelm M, et al. The target landscape of clinical
kinase drugs. *Science* (2017).

Kumagai A, Lee J, Yoo HY, Dunphy WG. TopBP1 activates the ATR-ATRIP complex.
*Cell* (2006).

Lee JA, Carvalho CM, Lupski JR. A DNA replication mechanism for generating
nonrecurrent rearrangements associated with genomic disorders. *Cell* (2007).

Lipinski CA, Lombardo F, Dominy BW, Feeney PJ. Experimental and computational
approaches to estimate solubility and permeability in drug discovery and
development settings. *Advanced Drug Delivery Reviews* 23(1-3):3-25 (1997).

Manning G, Whyte DB, Martinez R, Hunter T, Sudarsanam S. The protein kinase
complement of the human genome. *Science* 298(5600):1912-1934 (2002).

Miller BR 3rd, McGee TD Jr, Swails JM, Homeyer N, Gohlke H, Roitberg AE.
MMPBSA.py: an efficient program for end-state free energy calculations.
*Journal of Chemical Theory and Computation* (2012).

Mysinger MM, Carchia M, Irwin JJ, Shoichet BK. Directory of Useful Decoys,
Enhanced (DUD-E): better ligands and decoys for better benchmarking. *Journal
of Medicinal Chemistry* (2012).

Saldivar JC, Cortez D, Cimprich KA. The essential kinase ATR: ensuring
faithful duplication of a challenging genome. *Nature Reviews Molecular Cell
Biology* (2017).

Sieg J, Flachsenberg F, Rarey M. In need of bias control: evaluating chemical
data for machine learning in structure-based virtual screening. *Journal of
Chemical Information and Modeling* 59(6):947-961 (2019).

Tran-Nguyen VK, Jacquemard C, Rognan D. LIT-PCBA: An unbiased data set for
machine learning and virtual screening. *Journal of Chemical Information and
Modeling* 60(9):4263-4273 (2020).

Trott O, Olson AJ. AutoDock Vina: improving the speed and accuracy of docking
with a new scoring function, efficient optimization, and multithreading.
*Journal of Computational Chemistry* (2010). **[verified DOI:
10.1002/jcc.21334]**

Veber DF, Johnson SR, Cheng HY, et al. Molecular properties that influence
the oral bioavailability of drug candidates. *Journal of Medicinal Chemistry*
(2002).

Wang X, et al. Molecular architecture and inhibition mechanism of human
ATR-ATRIP. *Science Bulletin* 70(13):2137-2146 (2025). **[verified DOI:
10.1016/j.scib.2025.05.009; finding of four VE-822 molecules per complex (two
active-site, two dimer-interface) confirmed via WebSearch; full author list
and residue-level detail (e.g., specific TRD-domain boundaries) not
independently confirmed in this session — check the primary text]**
