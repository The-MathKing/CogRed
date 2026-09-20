# Response to External Review — Round 4 (JCAMD Elevation Strategy)

This round was not an external review of the manuscript, but a strategy
document proposing four concrete ways to elevate the pilot to a Journal of
Computer-Aided Molecular Design (JCAMD) standard: (1) actually execute a
larger-n survey and correlate redocking pass/fail with wwPDB Q-score/
residue inclusion, (2) test flexible-sidechain docking as a low-compute
alternative to MD, (3) reframe the BRICS pilot screen as a demonstration of
error propagation, (4) submit as a Short Communication. All four were
attempted for real; the results below are as found, not as hoped for.

## 1. The larger-n survey: run, and it does not confirm the hoped-for hypothesis

The redocking-validation pipeline was generalized into a single script
(`data/pilot_run_v1/survey_round4/survey_structure.py`) and applied to 14
additional, recently deposited cryo-EM structures — not more ATR/CHK1
structures (none were readily identifiable within this session's network
access beyond the two already used), but 14 structures of an unrelated
kinase complex, CDK7/cyclin H/MAT1 (the CDK-activating kinase, CAK), each
bound to a distinct ATP-competitive inhibitor from a fragment-elaboration
series (Cushing et al. 2024, *Nat. Commun.* 15:2265, DOI
10.1038/s41467-024-46375-9 — real, verified directly from the mmCIF
`_citation` records of the deposited structures, then cross-checked that
the DOI resolves to that title). Before trusting any result from the
generalized script, it was validated by re-running it, unmodified, on
9L4B — a structure with an independently known outcome from the original
hand-curated pipeline — and it reproduced that result almost exactly
(symmetry-corrected RMSD 0.60±0.03 Å vs. the original 0.61±0.02 Å;
Q-score 0.528, identical; residue inclusion 0.967, identical).

Combined with the original 3 structures, this gives n=17: **14/17 (82%)
fail redocking**, only 9L4B, 2YM8, and one CAK structure (8P75) pass.
Two things did *not* come out the way the strategy document hoped:

- **Resolution does not predict outcome** (Spearman ρ=−0.23, p=0.38).
- **wwPDB Q-score does not predict outcome, and trends the wrong way**:
  the 2 passing structures with a reported Q-score average 0.569; the 10
  failing ones average 0.688 (point-biserial r=−0.62, p=0.03). Residue
  inclusion shows no relationship at all (r=−0.14, p=0.61).

This is reported exactly as computed, including the fact that it argues
*against* the strategy document's implicit hope that these metrics would
"predict rigid-docking failure statistically." They don't, at least not
in this sample. The one qualifier that matters: 14 of the 17 structures
come from a single target's fragment-elaboration series (one research
group, one deposition window), so this is a sample of convenience, not a
systematic survey — the true failure rate and the true relationship (if
any) between validation-report metrics and redocking reliability could
differ in a properly stratified sample across many independent targets.
This is now stated explicitly in the paper's Limitations, not glossed
over.

One incidental finding sharpens the point without needing any statistic:
**8P74 and 8P75 are independent depositions of the same CAK complex bound
to the same inhibitor** (comp_id X3Z) — one fails redocking at 6.67 Å, the
other passes at 1.03 Å. Whatever drives redocking reliability here is a
property of the individual deposited coordinate set, not of the target or
the ligand.

A real bug was found and fixed en route: the wwPDB validation-XML
attribute regex assumed `resname` appears before `resnum` in the tag —
it doesn't (the real attribute order is `.../chain/resnum/.../resname/...`)
— so the original extraction silently returned no match for several
structures. Fixed by matching each `<ModelledSubgroup>` tag whole and
filtering attributes as an order-independent dict, keyed on
chain+resnum+resname together (resnum alone is not unique across chains).
Four CAK structures still genuinely lack a Q-score (verified by reading
the raw XML directly): they have dual-occupancy (split-altloc) ligands,
and the wwPDB pipeline does not report Q-score for those — a real absence
in the primary data, not an extraction failure.

## 2. Flexible-sidechain docking: does not rescue the pose

Vina's flexible-sidechain feature was used to let 5 active-site residues
(K2327, Y2365, E2378, W2379, H2477) move during 9L40 redocking, prepared
with Meeko (`mk_prepare_receptor.py`/`mk_prepare_ligand.py`). Result:
mean symmetry-corrected RMSD across 3 seeds = 3.00±0.92 Å (individual
seeds: 4.29, 2.23, 2.49 Å) — **worse, not better, than the fully rigid
result (2.53 Å)**, and no seed passed 2.0 Å. This is the "still fails"
outcome the strategy document itself anticipated as the alternative to
rescue, and it is what was found: it weighs against a simple
rigid-conformer artifact and toward either a genuinely less-supported
deposited pose or a conformational change larger than 5 side chains can
capture.

Two real bugs were hit and fixed while building this: (a) mixing an Open
Babel-typed ligand with a Meeko-typed flexible receptor crashed Vina's
native layer outright (`swig::stop_iteration`) rather than erroring
cleanly — fixed by re-preparing the ligand with Meeko too; (b) a
box-containment error ("ligand is outside the grid box") appeared only
once side chains were allowed to flex, even with unchanged starting
coordinates — fixed by enlarging the box from 22 Å to 30 Å to give
flexible rotamers room. A separate, unrelated crash (same
`swig::stop_iteration` message, initially misattributed to reusing
multiple `Vina()` objects in one process) was root-caused by systematic
step-by-step tracing to **import order**: `import vina` followed by
`from spyrmsd import ...` aborts at the spyrmsd import itself, before any
user code runs; importing spyrmsd first is a clean, reliable fix, now
applied everywhere both packages are used together.

## 3. Error propagation: demonstrated, not just asserted

The same 40-compound pilot library was docked against 9L40 (the excluded
structure) instead of 9L4B (the validated one), one seed per compound.
Result: Spearman rank correlation between the two receptors' scores across
all 40 compounds is ρ=0.19 (p=0.24, not significant); only 5 of the top-10
compounds by 9L4B score also appear in 9L40's top-10; mean absolute rank
shift is 12.3 positions (of 40), maximum 29. This is a real, computed
demonstration that skipping the validation step and using 9L40 by default
would have prioritized a substantially different, and on the evidence
above less trustworthy, compound set — not a hypothetical.

## 4. Short Communication reformatting: not done, and said so rather than guessed

`link.springer.com` (JCAMD's actual submission-guidelines page) is blocked
in this session, and a web search for JCAMD's specific Short Communication
word limit/format did not surface the actual figure. Rather than invent a
plausible-sounding word count or silently relabel the manuscript, the
paper keeps its current title, structure, and target formatting (Journal
of Cheminformatics-style) unchanged. If the actual authors want to pursue
JCAMD specifically, the real instructions-for-authors page needs to be
read directly (or the editorial office contacted) to get the Short
Communication requirements right rather than approximated.

## Net effect

Two of the four strategy items came back negative or partially negative
relative to what they were hoped to show (Q-score does not predict
redocking failure; flexible docking does not rescue the pose), and are
reported that way rather than reframed as confirmations. The other two
(the survey's scale, and the error-propagation demonstration) are real,
computed strengthenings of the paper's evidentiary base. The paper's
central claim is now broader and, on the evidence, *more* important than
before this round (structure-dependent redocking failure recurs across an
unrelated target at a similar rate, and no cheap validation-report metric
reliably flags it), while its Limitations are correspondingly more
explicit about the survey's convenience-sample nature. The fourth item
(venue-specific reformatting) was not actioned because it could not be
verified against a primary source, and this is stated rather than
papered over.
