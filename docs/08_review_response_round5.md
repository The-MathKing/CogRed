# Response to External Review — Round 5 (Re-review of the Survey)

This round attacked the round-4 survey directly, computed a real
statistical baseline from a cited paper, and independently verified two
specific claims from primary sources. Two of its findings were confirmed
and required withdrawing or correcting a claim; one raised a decisive
question that could only be answered by actually running the control it
demanded, which we did.

## §1: Is the 82%/93% CAK failure rate real, or a pipeline bug? — Resolved by running the control

The review computed that a 3/17 (or 1/14 CAK-only) pass rate is
~3.5σ below AutoDock Vina's documented 58% self-docking success rate on
the 285-complex PDBbind Core Set (Röhrig et al. 2023, *J. Chem. Inf.
Model.* 63(12):3925–3940 — verified: real paper, DOI
10.1021/acs.jcim.3c00054 confirmed, and its reported Vina/GOLD/AC-2.0
figures (58.0%/63.9%/73.3%) confirmed by web search), and demanded a
positive control: run the unmodified pipeline on the Astex Diverse Set or
an equivalent external benchmark before trusting the 82% claim.

We ran it. The exact Astex Diverse Set list could not be retrieved (most
academic domains — arxiv, PMC, Zenodo, deepwiki, Springer — are blocked
in this sandbox), but a real, redistributed CASF-2016 core-set listing
was found on GitHub (`OptiMaL-PSE-Lab/DeepDock`) — the exact benchmark
the review's own 58% figure is measured on, which if anything makes the
comparison more direct, not less. The first 20 PDB IDs in that file's
order (an objective rule, fixed before any of the 20 was inspected) were
run through the identical, unmodified pipeline, with the ligand of
interest auto-detected (largest non-excluded heteroatom group) rather
than looked up, since these targets were not chosen with any specific
ligand in mind.

**Result: 8/20 pass (40%).** Two comparisons follow:

- Against the literature's 58% baseline: not a significant departure
  (one-sided binomial p=0.08) — this pipeline is not badly broken; a
  modest gap from a simpler, generic protocol (fixed 22 Å box, Open
  Babel/Gasteiger charges, no local refinement) is expected.
- Against this pipeline's *own* measured 40% baseline: the CAK survey's
  1/14 remains a significant outlier (p=0.008), while the ATR pair's 1/2
  is fully consistent (p=0.84).

The CAK finding survives the review's own proposed decisive test. One
real bug was caught while reviewing this batch: 3DX1's auto-detected
ligand was a glycosylation sugar (NAG) instead of the real inhibitor
(YHO), because glycan codes were missing from the exclusion list. Fixed;
the corrected result (3.29 Å, still failing) did not change the count.
This is reported in full, including the bug, rather than only the
favorable headline number.

## §2: The 8P74/8P75 "flagship example" — withdrawn, verified wrong

The review checked both structures' deposition titles directly and found
they are *ring-up* and *ring-down* conformational classes of the same
ligand (ICEC0880), computationally separated by the depositors from one
dataset because it populates both conformations — not independent
depositions, as our text claimed. We verified this ourselves from the
raw mmCIF `_struct.title` records (exact match to the review's account)
before acting on it, per this project's standing rule of checking
external claims against primary sources rather than deferring to the
reviewer's account alone.

We also tested the review's own suggested alternative reading — that
Vina's "failing" pose for 8P74 is actually a correct match to the *other*
conformer (8P75's), unfairly scored against the wrong reference — by
computing symmetry-corrected RMSD between 8P74's docked pose and 8P75's
native ligand directly. It does not hold: 6.76 Å, no closer than 8P74's
own native (6.67 Å). The docked pose matches neither deposited conformer.
Both the withdrawal and this additional check are written into Results
and Discussion in full, rather than only removing the wrong claim.

## §3: The 11/14 vs. 13/14 arithmetic error — confirmed and fixed

The review recounted Table 4 and found 13/14 (93%) CAK structures fail,
not 11/14 (79%) as the abstract and discussion stated. We recounted
directly from `survey_combined.json` and confirmed: the review is right,
this was a real arithmetic error, not a matter of interpretation. Fixed
everywhere it appeared (abstract, results, discussion, conclusions).

## Major findings (M1–M5): accepted and acted on

| # | Finding | Verified | Action |
|---|---|---|---|
| M1 | 2YM8 resolution given as 2.00 Å in Table 4 vs. 2.07 Å elsewhere; Table 2's caption wrongly implies all resolutions come from the EM-only mmCIF field | Confirmed from 2YM8's own `_reflns.d_resolution_high` (2.07); the 2.00 was a transcription error introduced when hardcoding the original three structures into the round-4 stats script | Fixed in Table 4 and the resolution-correlation stats (ρ=−0.23→−0.28, still non-significant); Table 2 caption corrected to name both mmCIF fields |
| M2 | The flexible-sidechain test changed box size, receptor preparation, and flexibility all at once, then attributed the RMSD difference to flexibility alone | Confirmed by inspection of the two protocols | Ran the matched control (same 30 Å box, same Meeko prep, no flexible residues): 4.50±0.05 Å — worse than both the original rigid (2.53 Å) and the flexible (3.00 Å) result. Against this correct baseline, flexibility genuinely *improves* the pose (4.50→3.00 Å) — the opposite conclusion from what was first reported. Rewritten in Methods/Results/Discussion/Conclusions as a correction |
| M3 | A point-biserial r/p reported on a 2-vs-10 Q-score split is not statistically interpretable, especially in the abstract | Confirmed — a 2-item group cannot support an interpretable p-value | Dropped the statistic; now reported as a descriptive group-mean comparison only, in the abstract and throughout |
| M4 | The error-propagation comparison used 1 seed for 9L40 against a 3-seed-averaged 9L4B | Confirmed as a real design asymmetry | Re-ran 9L40 with the same 3 seeds. Result: materially unchanged (ρ=0.21 vs. 0.19, same 5/10 top-10 overlap) — the asymmetry was a valid concern to check, but correcting it did not change the conclusion |
| M5 | Title/scope mismatch — the title still described the v2 pilot, not this version's 17-structure survey | Agreed | Retitled to reflect the actual current scope, now that the survey is externally validated |

## Minor findings

- m1 (mixing X-ray/cryo-EM resolution conventions in one correlation): noted explicitly in the text now.
- m2 (Limitations vs. Conclusions scale consistency): kept the CAK-series/convenience-sample caveat at the same strength in both places.
- m3 (more figures would help): not actioned this round — noted as a real gap, deprioritized against the correctness issues above given session time.
- m4/m5 (author placeholders, Zenodo DOI): unchanged, as before — these require the actual authors and cannot be filled in on their behalf.

## Net effect

Two claims were withdrawn or reversed outright (the CAK fail-count
arithmetic; the flexible-docking "does not rescue, if anything worse"
reading, which is now "helps, but not enough"), one flagship example was
withdrawn as a mischaracterization (8P74/8P75), and the paper's single
most consequential claim (that the 82%/93% failure rate is real and not
a pipeline artifact) was subjected to the exact decisive test the review
demanded and survived it. Nothing in this round was defended without
being checked; two things that looked defensible on first read (the
flexible-docking comparison, the 8P74/8P75 pair) turned out not to be,
and are reported as corrections rather than quietly re-worded.
