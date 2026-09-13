# Response to External Review (Round 1)

Log of an external methodological review of the project outline, what was
accepted, what was contested, and what changed in this repository as a result.
Kept under version control so the revision history of the *argument* is as
traceable as the revision history of the code.

**Headline outcome:** the review's central objection is correct and forced a
change of the paper's central claim, not merely its wording. Docking, MD, and
MM-PBSA can support statements about *predicted binding, ranking, pose
stability, estimated energetics, and predicted selectivity*. They cannot
establish whether a ligand **activates, inhibits, or partially modulates**
ATR/CHK1. The prior outline was written as a functional-modulation paper while
the methods could only support a ligand-prioritization paper. That gap is now
closed by moving the claim down, not by adding more computation.

---

## Accepted in full

| # | Review point | Change made |
|---|---|---|
| 1 | "Modulator / stabilize rather than abolish signaling" is unproven pharmacology | Central claim changed to **candidate ligand prioritization**; "modulator" removed from title, abstract, and Results throughout `02_manuscript_outline.md` |
| 2 | ATR activation is protein-mediated (ATRIP/TOPBP1/ETAA1), not a function of ATP-pocket occupancy | New Introduction §1.4 states this explicitly and drops the implied inhibition↔activation affinity continuum; see "Contested" §C below for what survives of the original hypothesis |
| 3 | Benchmark leakage between ChEMBL seeds, BRICS derivatives, and validation actives | `scripts/09_pipeline_validation.py` rewritten: Bemis-Murcko scaffold-split, explicit seed↔validation leakage audit, nearest-neighbour Tanimoto reporting for every held-out hit |
| 4 | Single ROC curve is insufficient; need baseline comparison | Same script now runs **baseline Vina vs. consensus vs. full pipeline**, reports ROC-AUC + PR-AUC + EF at 1/5%, with bootstrap CIs |
| 5 | "ADMET/toxicity triage" implies measurement | Renamed throughout to **in silico developability and toxicity prediction**; script renamed `04_developability_prediction.py`, output column `predicted_developability_pass` |
| 6 | "Kinome selectivity panel" implies an assay | Renamed **predicted cross-kinase selectivity analysis** in docs and `scripts/08_*` |
| 7 | Docking scores across different proteins are not directly comparable | `scripts/08_*` no longer subtracts raw affinities; it z-normalizes each kinase's score against a property-matched background set docked into that same receptor, and benchmarks the metric against reference inhibitors with published kinome profiles |
| 8 | MM-PBSA values must not be presented as ΔG / experimental affinity | `scripts/07_*` output column renamed `mmpbsa_estimated_energy_kcal_mol`; docs specify "relative energetic estimation", entropy neglect stated, sensitivity to trajectory window required |
| 9 | Desirability weights can determine the result | New `scripts/10_integrated_ranking.py`: pre-specified weights **plus** Dirichlet Monte-Carlo sensitivity analysis reporting rank stability of each candidate across weighting schemes |
| 10 | MD spec too thin (≥3 replicas, convergence, seeds) | Methods §2.7 expanded; verified against JCIM's MD reporting guidelines (Guidelines for Reporting Molecular Dynamics Simulations in JCIM Publications, *JCIM* 2023, DOI 10.1021/acs.jcim.3c00599) — ≥3 replicas from different velocities/coordinates, disclosed RNG seeds, statistical variance discussed |
| 11 | "Mid-affinity window" does not follow from MD/MM-PBSA | Demoted from proposed in silico proxy to an explicitly labelled untested hypothesis in Discussion §4.2, with the reviewer's suggested non-monotonic phrasing adopted |
| 12 | CFS sequence features are not universal determinants | Introduction §1.1 now separates *established determinants* (replication timing, origin paucity, transcription-replication conflict, large-gene architecture) from *candidate sequence-level contributors* (AT-rich flexibility, non-B structure) |
| 13 | Scientific Reports limits (20-word title, 200-word abstract) | Verified independently — also a **4,500-word main-text limit** the review did not mention. All three now constrain the outline |
| 14 | JCIM scope excludes straightforward docking without experimental validation | Verified against JCIM's stated scope. Journal strategy reordered — see below |

## Verified independently (not taken on the reviewer's word)

- **Scientific Reports:** title ≤20 words, abstract ≤200 words, main text ≤4,500
  words (excluding Abstract, Methods, References, figure legends).
- **JCIM scope:** explicitly does not consider straightforward applications of
  docking to a single target system without adequate experimental validation.
- **JCIM MD reporting:** ≥3 replica simulations, ideally from different
  coordinates/velocities, with RNG seeds disclosed and statistical variance
  discussed.

All three check out, and they change the journal ranking materially (see
`03_journal_targeting_strategy.md`, now reordered).

---

## Contested / modified rather than adopted

### A. "Make the CFS analysis drive the computation, or move it to SI"

The review prefers Option 1 (make it causal). **We take Option 2, plus a
specific job for the CFS analysis** — and the reasoning matters more than the
choice.

Genomic features of FRA3B/FRA16D cannot legitimately determine which kinase
ligand ranks highest. Replication timing and AT-content do not inform ATP-pocket
complementarity. Constructing a pipeline in which the genomic analysis
*appears* to drive compound prioritization would manufacture exactly the
decorative integration the review objects to, only harder to detect. That is a
worse outcome than an honest structural separation.

Instead the CFS analysis is given a real, non-decorative function: **it defines
the falsifiable experimental prediction the computational work terminates in.**
The genomic features identify which loci should show a differential response
(FRA3B/FRA16D gap-and-break frequency under aphidicolin, versus non-fragile
control loci matched for replication timing) if any prioritized compound turns
out to engage the pathway functionally. The analysis moves to Supporting
Information as a descriptive characterization, and its result is *used* in
Discussion §4.5 to specify the assay, the loci, and the expected direction of
effect. It earns its place by making the paper's hypothesis testable rather
than by pretending to have shaped the docking.

### B. "Demonstrate the pipeline outperforms conventional baselines" → novelty

Accepted as the right experiment, with one addition the review did not make:
**we pre-commit to reporting the comparison honestly if it comes out negative.**

The realistic outcome is that a Vina + RDKit + MD + MM-PBSA stack does *not*
substantially outperform well-tuned baseline docking on scaffold-split
retrieval. If we only report the benchmark when it favors the pipeline, we
reproduce the selective-reporting problem under a rigor label. So the analysis
plan is fixed in advance (`09_pipeline_validation.py`, weights and splits
pre-specified), and a null or negative result is reported as the finding. A
leakage-controlled negative benchmark of a widely used open-source CADD stack
is a genuinely useful contribution — arguably more useful than another
pipeline that claims to win.

This also settles the novelty question honestly: **assembling standard tools is
not a methodological contribution.** What can be one is the leakage-controlled
benchmarking protocol and the cross-target score normalization for selectivity.
If those do not survive their own tests, the paper is a well-scoped
ligand-prioritization study and should be submitted as such, not dressed up.

### C. What survives of the original modulator hypothesis

Not "mid-affinity binding produces intermediate signaling." What survives is
narrower and testable: prioritized compounds span a range of predicted target
engagement, which makes them a **usable input set** for a future functional
experiment that asks whether pathway output varies non-monotonically with
engagement. The computation supplies the compounds and the rationale for the
experiment; it does not supply the answer. This is stated as such in Discussion
§4.2 and the Conclusion.

### D. Terminology: "prospective"

The review uses "prospective/holdout performance." We avoid "prospective"
entirely. Without experimental testing there is no prospective validation —
only held-out retrospective evaluation. Calling it prospective is the same
category of overclaim as calling predictions triage.

---

## Net effect on the project

Nothing computational was removed. The pipeline scripts are largely unchanged
in what they compute; what changed is the claim they support, the validation
that surrounds them, and the language used to report them. Publication-readiness
by the review's own scoring should move primarily through scope correction, not
new experiments — with the significant caveat that the benchmarking in §B may
return a negative result, which we will report.

**Remaining unresolved:** whether the project should acquire a wet-lab
collaborator. If functional validation of even one or two top compounds
(biochemical ATR/CHK1 kinase assay) becomes available, the modulation claim
becomes defensible and JCIM/PLOS CB re-enter scope. Worth pursuing in parallel
with drafting rather than after it.
