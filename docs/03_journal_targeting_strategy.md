# Journal Targeting Strategy (Revision 2)

Reordered after external review. The previous version named JCIM as the primary
target; **that was wrong** and is corrected below. The correction follows from
a verified scope statement, not from a change of ambition.

## Verified constraints (checked against journal sources, September 2026)

- **Scientific Reports:** title ≤20 words; abstract ≤200 words; main text
  ≤4,500 words excluding Abstract, Methods, References, figure legends.
- **JCIM scope:** explicitly does **not** consider straightforward applications
  of molecular docking to a single target system without adequate experimental
  validation.
- **JCIM MD reporting guidelines** (*JCIM* 2023, DOI 10.1021/acs.jcim.3c00599):
  ≥3 replica simulations, ideally from different coordinates/velocities, RNG
  seeds disclosed, statistical variance discussed. **Adopted for this project
  regardless of target journal** — it is simply correct practice.

---

## 1. Scientific Reports — most realistic target for the current study
- **Audience:** Broad multidisciplinary; covers computational biology, drug
  discovery, and molecular modeling.
- **OA status:** Fully open access; waiver program available.
- **Why it fits:** the study is technically sound computational work with
  honestly scoped claims. Scientific Reports evaluates technical soundness
  rather than demanding novelty or experimental validation, which matches what
  this project can support today.
- **Must hit:** the hard format limits above (the previous 250-word abstract
  and long title both violated them); code availability; honest scoping of
  every predicted quantity; the scaffold-split benchmark with baselines.
- **Risk:** low, provided the modulator language is gone.

## 2. PLOS Computational Biology — possible, but only with a unified biological story
- **Audience:** Computational biologists; weighs biological insight and
  significance heavily, not just technical execution.
- **OA status:** Fully open access; fee waivers available.
- **Why it might fit:** if the paper delivers a genuine biological insight
  rather than a candidate list. The journal also expects computational
  discovery to be validated or enriched by experiment or real-world data
  *where possible*.
- **Blocker:** as outlined, the ligand-prioritization work and the CFS biology
  are structurally separate (see `04_review_response.md` §A — we judge that
  separation to be honest rather than fixable by reorganization). Without a
  single integrated biological argument, this is a moderate-to-low fit.
- **Would become a strong fit if:** even minimal functional data (a biochemical
  ATR/CHK1 kinase assay on 2-3 top candidates) were added, which would also
  restore the mechanistic claim the review correctly stripped out.

## 3. Frontiers in Pharmacology — viable, with the selectivity work as the draw
- **Audience:** Pharmacology-focused; receptive to target-mechanism-centered
  computational work.
- **OA status:** Fully open access.
- **Why it fits:** the predicted cross-kinase selectivity analysis — especially
  benchmarked against reference compounds with published kinome profiles — is
  well matched to this readership.
- **Must hit:** this audience will scrutinize mechanism claims hardest, so the
  Introduction §1.4 scope paragraph (ATP-pocket occupancy ≠ functional
  direction) is load-bearing here. Interactive review can help a
  methodologically careful but mechanistically modest paper.

## 4. ACS Omega — honest fallback
- **Audience:** Broad chemistry; accepts exploratory computational work.
- **OA status:** Fully open access.
- **Why it fits:** if MD/MM-PBSA must be scoped down to a small shortlist for
  compute reasons, this venue accommodates that — provided the limitation is
  stated rather than concealed.

## 5. JCIM — **not a current target**; re-enters scope only under specific conditions
- **Why it was demoted:** its scope statement excludes straightforward docking
  applications to a single target without experimental validation. A pipeline
  assembled from Vina + RDKit + SwissADME + OpenMM + MM-PBSA is not, by itself,
  a new methodology. The previous ranking of JCIM as primary target reflected
  ambition rather than the journal's stated criteria.
- **Re-enters scope if either:**
  - **(a) Experimental validation** is obtained for top candidates (even a
    single biochemical kinase assay), **or**
  - **(b) Demonstrated methodological contribution** — the leakage-controlled
    scaffold-split benchmarking protocol and the within-receptor score
    normalization for cross-kinase selectivity, shown quantitatively to change
    conclusions relative to conventional practice. Note this is a real
    empirical question: if the benchmark shows the multi-stage pipeline does
    *not* outperform baseline docking (a likely outcome), that is a publishable
    negative result but not a methodological novelty claim, and the paper
    should go to Scientific Reports instead.

---

## Recommended sequencing

1. **Draft to Scientific Reports constraints** (20-word title, 200-word
   abstract, 4,500-word main text) as the default target. These limits are the
   tightest of the five, so drafting to them keeps every other venue reachable
   without restructuring.
2. **Run the benchmark before choosing the final venue.** The comparison in
   Results §3.2 determines what kind of paper this is. Strong, honest
   outperformance → consider JCIM route (b) or Frontiers. Null or negative →
   Scientific Reports, with the benchmark as the headline contribution.
3. **Pursue a wet-lab collaborator in parallel with drafting, not after.** One
   biochemical kinase assay on the top candidates is the single highest-value
   addition available to this project: it unlocks JCIM route (a), materially
   strengthens PLOS Computational Biology, and restores the functional claim
   that the review correctly required us to drop.
