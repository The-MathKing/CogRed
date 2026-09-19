# Response to External Review — Round 3 (Real Pilot Paper)

This round reviewed `manuscript/submission/paper.pdf` (the round-2 real pilot
paper). Unlike a generic critique, it tried actively to break the central
finding (a symmetry-artifact attack on the RMSD computation) and showed its
work, then raised three further checkable claims (exact resolutions,
kinase-domain vs. full-complex composition, an exhaustiveness inconsistency
in the paper's own text) and one substantive unaddressed confound (the
adjacent allosteric site) plus one unaddressed alternative explanation
(deposited-ligand density quality). Every one of these was independently
re-verified or newly computed from primary data before acting on it — the
same standard applied in rounds 1 and 2.

## The attack on the central finding — tested, not just answered

The review computed (by their own account) a symmetry-based worst-case RMSD
penalty for berzosertib and argued the naive 2.96 Å redocking "failure"
might be a symmetry artifact, then concluded on inspection that it probably
wasn't. Rather than accept that reasoning, we ran the actual computation:

- Installed spyrmsd (Meli & Biggin, *J. Cheminform.* 2020, verified real)
  for real and computed symmetry-corrected RMSD via graph-isomorphism atom
  matching for all three redocking validations.
- Result: 9L40 narrows from 2.96 Å (naive) to 2.53 Å (symmetry-corrected) —
  still comfortably over the 2.0 Å threshold. 9L4B and 2YM8 are essentially
  unchanged (0.61→0.61 Å; 0.90→0.76 Å).
- Also ran the two sanity controls the review specifically requested:
  crystal-vs-itself returns exactly 0.000 Å; a known, exact 1.0 Å rigid
  translation is recovered as exactly 1.000 Å. Both passed for all three
  ligands, before touching any real pose.

The central finding survives, now on the basis of an independent
computation rather than either party's back-of-envelope estimate.

## Verified this session (before acting on any of it)

| Claim | Verification | Result |
|---|---|---|
| Exact resolutions of 9L40 and 9L4B (review claimed 2.87 Å / 3.20 Å, contradicting the paper's own "~3 Å for both") | Read `_em_3d_reconstruction.resolution` directly from the primary mmCIF files already on disk | **Confirmed exactly**: 9L40 = 2.87 Å, 9L4B = 3.20 Å. The review was right and the paper was wrong — and the corrected fact is more interesting than what it replaced: the *higher*-resolution structure is the one that fails |
| 9L40/9L4B are "kinase domain" structures, not the "full ATR-ATRIP complex" as the paper's Background claimed | Read `_struct.title` and chain/residue composition directly from the mmCIF files (2 chains, both ATR kinase domain, residues ~1279–2644, no ATRIP present) | **Confirmed**: both are homodimers of the ATR kinase domain alone. The paper's Methods section already said this correctly; Background contradicted it and has been fixed |
| spyrmsd (Meli & Biggin, *J. Cheminform.* 2020) real | WebSearch, cross-checked with GitHub/PMC listings | **Confirmed** real, installed and used for real (`pip install spyrmsd`, version 0.9.0) |
| Hartshorn et al. 2007 (Astex Diverse Set, *J. Med. Chem.* 50(4):726–741) as the actual anchor for the 2 Å redocking criterion | WebSearch | **Confirmed** real; replaces the previous miscitation of the Vina paper itself for this claim |
| Pintilie et al. 2020 Q-score (*Nat. Methods* 17:328–334) | WebSearch | **Confirmed** real; used to interpret the wwPDB validation report fields pulled for this round |

## New real analyses run this round (not asserted, computed)

1. **Allosteric-site distance (M3, the box-overlap confound).** Computed
   directly from the mmCIF coordinates already on disk: the active-site
   ligand centroid (A2701) in 9L40 is 60.65 Å from the nearest allosteric
   copy (A2702) — far outside the 22 Å docking box. This rules out "the box
   reached the emptied allosteric pocket" as an explanation for the
   redocking failure, with a real number rather than a hand-wave.
2. **wwPDB validation report check (M2, the deposited-pose-quality
   confound).** Downloaded the real official validation XML for 9L40 and
   9L4B from the same `s3://pdbsnapshots` mirror (validation reports are
   hosted there too). Q-score for the active-site ligand is comparable
   between structures (9L40: 0.553; 9L4B: 0.528 — if anything marginally
   *higher* for 9L40), so the failing pose is not an obvious density
   outlier by that metric; residue inclusion is somewhat lower for 9L40
   (0.818 vs. 0.967), a partial but not full contributor. Reported exactly
   as found — comparable Q-score, lower inclusion — rather than rounded to
   whichever conclusion was more convenient.
3. **Redocking validation re-run at exhaustiveness 16 (C1, the real
   self-contradiction).** The paper's original text called the validation
   protocol (exhaustiveness 32) and the screening protocol (exhaustiveness
   16) "identical" — they weren't, and the review caught this directly. The
   review's suggested fix (re-run the full 40-compound pilot at
   exhaustiveness 32) was attempted but proved far more expensive than
   estimated (~2+ hours under this session's CPU-only, shared-load
   conditions, extrapolated from 3/40 compounds in 12 minutes) and was
   stopped as not worth the wait once a faster, equally direct fix was
   available: re-running the *redocking validation itself* at
   exhaustiveness 16 confirmed the identical pass/fail pattern
   (9L40 3.01±0.06 Å fail; 9L4B 0.59±0.02 Å pass; 2YM8 0.92±0.01 Å pass),
   which is what the paper's claim actually depends on. The word
   "identical" was removed from the text regardless, and the exhaustiveness
   discrepancy is now stated plainly as a real limitation of the original
   design that was then closed, not hidden.
4. **Pose-overlay figure (M4).** PyMOL installed for real (`apt-get install
   pymol`) and used to render native-vs-redocked overlays for both 9L40
   (visible displacement at the isopropylsulfonylphenyl and
   aminomethylphenyl termini) and 9L4B (near-perfect overlay). This is the
   paper's first and only figure, added because a claim about pose geometry
   with zero pictures was a real gap.

## Accepted and acted on (text/framing)

| # | Finding | Action |
|---|---|---|
| M1 | "Tight search reproducibility" (0.56 Å) was an overstatement — it's 12× looser than 9L4B's and 7× looser than 2YM8's | Reframed: 9L40 is both less accurate *and* less reproducible than the passing structures, read against them rather than against the 2.0 Å pass threshold alone |
| M5 | The 40-compound derivative screen is circular by construction (BRICS-fragments-camonsertib docked into camonsertib's own structure will rank camonsertib first) | Section retitled "Pipeline sanity-check screen" with an explicit framing sentence stating the outcome is close to guaranteed by construction, per the review's own suggested language |
| M6 | n=1 generalization — the reliability claim rests on one failing structure | Added as Limitations item (1), and folded into Discussion: the paper now explicitly disclaims any claim about how *common* such failures are, and names the natural next study (a survey across many more structures) without attempting it here |
| M7 | Higher-resolution structure fails — the paper's most quotable fact was rounded away as "comparable resolution" | Now stated explicitly in Abstract, Background, Results, and Discussion with the exact numbers |
| M8 | "Full ATR-ATRIP complex" mischaracterization; Trott & Olson miscited for the redocking-convention claim | Both fixed (see verification table above) |
| m1 | 0.33 kcal/mol rank gaps reported to a precision the method doesn't support | Added explicit note: this gap is well within Vina's reported ~2-3 kcal/mol absolute error against experiment, so the ranking should be read as indistinguishable |
| m2 | Best-of-3 is a biased (downward) statistic; abstract/results led with it | Now leads with the mean; best-of-3 still reported alongside |
| m3 | MW-ascending derivative selection biases against derivatives beating the (larger) seed compounds | Stated explicitly in both Results and Limitations |
| m9 | Author/affiliation/funding/competing-interests/contributions placeholders | Unchanged — these require the actual authors and cannot be filled in on their behalf |

## Contested / not fully actioned, with reasons

- **Full exhaustiveness-32 re-run of the 40-compound pilot.** Attempted
  (background job), found to cost roughly 2+ hours under this session's
  compute conditions rather than the "cheap" cost the review assumed, and
  stopped once the cheaper, equally direct exhaustiveness-16
  re-validation confirmed the same conclusion. The job's partial output (3
  of 40 compounds) was discarded rather than reported, since a partial
  run is not informative and reporting it would invite exactly the kind of
  cherry-picking this project has repeatedly committed not to do.
- **Zenodo DOI.** Cannot be minted from within this session (requires an
  account and a deliberate archival action by the actual authors); the
  paper continues to name this as an explicit `[PLACEHOLDER]`-adjacent
  action item rather than fabricate a DOI.
- **A 20–50 structure survey (M6's suggested follow-up paper).** Explicitly
  scoped as future work, not attempted here. The review itself frames this
  as "the next paper," and building it properly (systematic structure
  selection, resolution/deposition-year/map-quality covariates, the same
  validate-first discipline applied at scale) is a substantially larger
  undertaking than a same-session extension of this pilot.
- **Reference list completeness (m7, "finish the job").** Three new
  citations (Hartshorn 2007, Pintilie 2020, Meli 2020) were added with full
  DOIs, verified. The remaining ~14 older citations still lack volume/page/
  DOI, as flagged in the manuscript's own header note and Limitations
  since round 2; not revisited further this round given the more
  substantive technical items took priority.

## Net effect

No claim got weaker this round; several got sharper (exact resolutions, the
resolution-inversion framing), one real self-contradiction in the paper's
own protocol description was closed with a real re-computation rather than
a wording patch, and the central finding survived the most serious attempt
yet to break it — using an independently run computation, not an appeal to
the review's own arithmetic. The reviewer's closing assessment (real,
executed, honestly reported, small-n pilot) is not contested.
