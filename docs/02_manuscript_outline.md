# Manuscript Outline (Revision 2)

Revised after external review — see `04_review_response.md` for what changed
and why. The central claim is now **candidate ligand prioritization**, not
functional modulation.

**The evidentiary boundary, stated once and enforced everywhere below:**
this study can support claims about predicted binding, relative ranking,
retrieval performance under scaffold-split evaluation, conformational
stability, estimated relative energetics, predicted selectivity, and predicted
developability. It cannot support claims about kinase activation, inhibition,
partial modulation, pathway output, fork stability, or any cellular or
genomic phenotype. Every sentence in the manuscript must fall on the
supportable side of that line or be explicitly marked as hypothesis.

Format sized for Scientific Reports (primary target): **title ≤20 words,
abstract ≤200 words, main text ≤4,500 words** excluding Abstract, Methods,
References, and figure legends. See `03_journal_targeting_strategy.md` for
per-journal adjustments.

## Title (working, 14 words)
*"Leakage-Controlled Benchmarking and Structure-Based Prioritization of
Candidate ATR and CHK1 Ligands for Replication-Stress Research"*

Alternatives if the benchmark returns a negative result (likely — see
`04_review_response.md` §B), in which case the benchmark becomes the headline:
*"Scaffold-Split Benchmarking Limits the Apparent Advantage of Multi-Stage
Virtual Screening for ATR and CHK1"*

Note what the title no longer claims: no "modulator", no "suppress replication
stress", no disease name. Fragile-site biology is context in the Introduction
and the source of the terminal experimental prediction, not a result.

## Abstract (≤200 words for Scientific Reports)
Sequence: problem → gap → approach → quantitative result → scoped significance.
Every predicted quantity must carry the word *predicted* or *estimated*.

- Problem (1 sentence): replication stress at common fragile sites; ATR-CHK1
  as the checkpoint axis maintaining replication integrity.
- Gap (1 sentence): existing chemical matter is overwhelmingly inhibitory and
  oncology-directed; candidate ligands have not been systematically
  prioritized under leakage-controlled evaluation.
- Approach (2-3 sentences): ChEMBL/PubChem seed set → derivative generation →
  scaffold-split benchmarking against baseline docking → consensus scoring →
  predicted developability → MD/MM-PBSA on shortlist → predicted cross-kinase
  selectivity.
- Result (2 sentences): benchmark numbers (ROC-AUC, PR-AUC, EF1% vs. baseline,
  with CIs) **and** the candidate set. Report the benchmark result whichever
  direction it goes.
- Significance (1 sentence): candidates and the open, reusable evaluation
  protocol are offered as **testable hypotheses for functional evaluation** —
  not as demonstrated modulators.

## 1. Introduction
1.1 **Replication stress and common fragile sites.** Distinguish *established
    determinants* of fragility — late/delayed replication timing, origin
    paucity, transcription-replication conflict, large-gene architecture —
    from *candidate sequence-level contributors* (AT-rich flexibility peaks,
    non-B DNA structure). Do not present the latter as universal defining
    properties; CFS fragility is now understood as emergent from multiple
    interacting features.
1.2 **The ATR-CHK1 axis.** Fork stabilization, origin firing suppression,
    S/G2-M checkpoint. Draw a hard distinction between (a) physiological ATR
    activation, which proceeds through ATRIP/TOPBP1/ETAA1 protein-mediated
    mechanisms, and (b) pharmacological inhibition of kinase activity by
    ATP-competitive ligands. Do **not** imply a single continuum
    (inhibition ← normal ← activation) traversable by tuning binding affinity.
1.3 **Chemical-matter gap.** Known ATR/CHK1 ligands are inhibitors developed
    for oncology synthetic-lethality strategies. State the gap as one of
    *systematic, leakage-controlled prioritization of candidate ligands*, not
    as a gap in modulator discovery — we are not demonstrating modulators.
1.4 **Scope of what computation can establish here.** Short, explicit
    paragraph: ATP-pocket occupancy does not determine functional direction.
    The study therefore prioritizes candidates and defines the experiment that
    would resolve their pharmacology. Stating this in the Introduction (rather
    than burying it in Limitations) is what earns reviewer trust for
    everything that follows.
1.5 **Disease relevance, tiered and terminating early.** Tier 1 (literature-
    supported): replication stress → CFS instability; ATR/CHK1 contributes to
    replication integrity. Tier 2 (hypothesis generated here): altered
    ATR/CHK1 signaling *might* influence the probability of
    replication-associated rearrangement at vulnerable loci. Tier 3 (distant
    motivation, one sentence maximum): eventual relevance to constitutional
    genomic disorders. **Cri du Chat is not required to motivate this work and
    should not appear in the Abstract or Results.**
1.6 Objectives and pipeline overview (Figure 1).

## 2. Methods
2.1 **Structural preparation.** Per structure: PDB ID, resolution, species/
    isoform, biological construct, kinase activation state, co-crystallized
    ligand, cofactor/metal treatment, retained vs. removed waters, alternate
    conformations, missing residues and how modeled, protonation/tautomer
    assignment protocol and pH. Justify why each structure represents a
    biologically relevant state — this is a scope question, not a formality.
2.2 **Seed library.** PubChem/ChEMBL programmatic retrieval, bioactivity
    threshold, standardization. Record provenance per compound; the leakage
    audit in 2.5 depends on knowing exactly what entered here.
2.3 **Derivative generation.** BRICS/MMP enumeration with parameters, plus
    property pre-filters. Explicitly **not claimed as methodological novelty**
    — these are standard cheminformatics operations and are described as such.
2.4 **Docking.** Vina version, exhaustiveness, ≥3 independent seeded runs per
    ligand (seeds disclosed), box coordinates and their derivation, ensemble
    docking across conformers if used.
2.5 **Evaluation protocol (the methodological core).** Bemis-Murcko
    scaffold-split of actives; property-matched decoy construction; explicit
    seed↔evaluation-set leakage audit with nearest-neighbour Tanimoto
    reported for every retained candidate; baseline (Vina alone) vs. consensus
    vs. full-pipeline comparison; ROC-AUC, PR-AUC, EF1%/EF5% with bootstrap
    confidence intervals. Weights, splits, and metrics **pre-specified before
    running** — state this, and report the result in whichever direction it
    falls.
2.6 **In silico developability and toxicity prediction.** (Renamed from
    "ADMET triage.") Descriptor rules computed reproducibly in RDKit;
    SwissADME/ProTox-II outputs reported as *predictions* with access dates
    and versions. Thresholds justified or presented as continuous scores.
    State plainly: predicted absence of toxicity is not evidence of safety.
2.7 **Molecular dynamics.** OpenMM version, force fields (protein + ligand
    parametrization), solvation/ion model, system size, minimization/
    equilibration/production protocol. **≥3 replicas from different initial
    velocities/coordinates, RNG seeds disclosed, convergence assessment and
    replicate-to-replicate variance reported** (per JCIM MD reporting
    guidelines, DOI 10.1021/acs.jcim.3c00599 — adopted regardless of target
    journal). Analyses: protein and ligand RMSD, RMSF, contact persistence,
    H-bond occupancy, clustering. An RMSD plateau is reported as conformational
    stability, never as evidence of binding strength.
2.8 **MM-PBSA relative energetic estimation.** Tool, settings, frames used,
    frame correlation, dielectric assumptions, entropy treatment (neglected by
    default — stated, not hidden), error estimation, sensitivity to trajectory
    window. Reported as relative ranking energetics, never as ΔG or affinity.
2.9 **Predicted cross-kinase selectivity analysis.** (Renamed from "kinome
    panel.") Off-target kinases with PDB IDs and rationale (PIKK family first).
    Scores z-normalized **within each receptor** against a property-matched
    background set docked into that same structure, because raw scores are not
    comparable across proteins. Metric benchmarked against reference
    inhibitors with published experimental kinome profiles.
2.10 **Integrated ranking.** Pre-specified desirability weights with stated
    rationale, **plus** Monte-Carlo sensitivity analysis over the weight
    simplex reporting how often each candidate remains top-ranked.
2.11 **Reproducibility.** Software version table, hardware, all seeds, code
    repository and archived release DOI.
2.12 *(Supporting Information)* CFS genomic feature analysis — data sources
    and features computed. Moved out of the main Methods; its role is to
    specify the terminal experimental prediction (§4.5), not to inform
    compound selection.

## 3. Results
Ordered so the pipeline must earn trust before it is used. This ordering is
the paper's main structural argument.

3.1 **Does the pipeline recover known chemistry?** Scaffold-split retrieval on
    held-out actives: ROC-AUC, PR-AUC, EF1%/EF5% with CIs. Leakage audit
    reported here, not in SI.
3.2 **Does it beat simpler approaches?** Baseline Vina vs. consensus vs. full
    pipeline, per-stage contribution. Report faithfully — a null result here
    is the finding, not a failure (Figure 3).
3.3 **Chemical space of the derivative set.** Framed as a question — did
    derivative generation expand accessible chemical space while retaining
    ATR/CHK1-like character? — not as a decorative PCA plot.
3.4 **Candidate ranking under docking and consensus scoring.**
3.5 **Predicted developability outcomes.** How many candidates survived, why
    others were excluded, property distributions before/after.
3.6 **Conformational robustness (MD).** Replicate variance shown, not hidden.
3.7 **MM-PBSA relative energetics** with uncertainty.
3.8 **Predicted cross-kinase selectivity**, including the reference-compound
    benchmark that establishes whether the metric is trustworthy at all.
3.9 **Integrated ranking and its stability** under weight perturbation.

Presented as five distinct evidence dimensions — pose plausibility, retrieval
performance, conformational robustness, predicted developability, predicted
selectivity — **not** as four independent confirmations of one truth. Their
assumptions overlap and their uncertainties are correlated; say so.

## 4. Discussion
4.1 Comparison with known ATR/CHK1 pharmacology, quantitatively — avoid
    anecdotal "our docking recovered known inhibitors" claims.
4.2 The engagement-range hypothesis, explicitly labelled untested: the
    candidate set spans a range of predicted target engagement, providing a
    framework for determining experimentally whether functional effects vary
    non-monotonically with engagement. Never stated as a demonstrated
    mechanism.
4.3 Predicted selectivity liabilities — and the more interesting mechanistic
    question of whether they arise because kinase ATP pockets are intrinsically
    hard to discriminate structurally. Turns a negative result into insight.
4.4 **Limitations.** Chemical-space bias from ChEMBL-derived seeds; residual
    benchmark leakage risk; non-equivalence of docking scores across targets;
    predicted-ADMET uncertainty; MM-PBSA model dependence; protein-state and
    regulatory-context uncertainty (ATR function depends on ATRIP/TOPBP1/ETAA1
    context absent from these models); and — stated so it cannot be missed —
    **no functional directionality: computational binding does not establish
    activation versus inhibition.**
4.5 **Translational path.** Staged and specific: biochemical ATR/CHK1 kinase
    assay → cellular checkpoint signaling readout → replication-fork phenotype
    → CFS-specific readout (FRA3B/FRA16D gap-and-break frequency under
    aphidicolin versus replication-timing-matched non-fragile control loci,
    the loci specified by the SI genomic analysis). This is where the CFS
    biology does real work: it makes the hypothesis falsifiable at a named
    locus with a named assay.

## 5. Conclusion
Restrained. Something close to: *"We report a reproducible, leakage-controlled
computational framework that prioritizes candidate ATR and CHK1 ligands and
benchmarks its own retrieval performance against baseline docking, yielding a
candidate set and a specified experimental test."* Never: *"we developed
modulators that suppress replication stress."*

## Figures (6 main text)
1. Pipeline schematic, with the evidentiary boundary drawn on it.
2. Chemical space, seed vs. derivative set.
3. **Benchmarking (the methodological heart):** (A) actives/decoy
   construction, (B) baseline docking ROC/PR, (C) consensus ROC/PR, (D)
   held-out performance, (E) scaffold-split performance. Note: "held-out",
   never "prospective" — there is no prospective validation without experiments.
4. Binding poses — illustrative of predicted interactions, explicitly not
   offered as mechanistic proof.
5. MD replicate variance + MM-PBSA **estimated relative energetics** with
   uncertainty (axis label must not read ΔG).
6. Predicted cross-kinase selectivity heatmap, with the reference-compound
   benchmark panel alongside it.

## Supporting Information
CFS genomic feature analysis; full compound table (SMILES + all scores);
software/version table; per-replicate MD analyses for all shortlisted
compounds; pre-specified analysis plan; code repository and release DOI.
