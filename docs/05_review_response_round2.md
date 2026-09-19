# Response to External Review — Round 2 (Manuscript Draft v1)

This round reviewed `manuscript/draft_v1.md` itself, not the repository
structure. Unlike round 1, several of its claims were specific and checkable
(structural biology, benchmarking literature), so before changing anything I
verified them via WebSearch rather than accepting them on say-so — the same
standard applied to round 1's journal-scope claims. All results below.

## Correction, not a review finding: venue

The review opened by asserting this manuscript was meant for "the 15th
International Conference on Complex Networks" and that the project had the
wrong file attached. **This has no basis in this repository's history.**
Every document in this repo — `docs/03_journal_targeting_strategy.md`, the
manuscript's own header — has targeted *Scientific Reports* from the first
commit. No instruction anywhere in this project ever named a network-science
venue. The reviewer's own §0 even notes the manuscript is formatted to
Scientific Reports limits, which sits oddly next to the claim that it was
meant for a different conference entirely. This is flagged back to the user
rather than acted on; nothing in the venue strategy changed as a result.

## Verified this session (WebSearch, before any edit)

| Claim | Verification | Result |
|---|---|---|
| PDB 5YZ0 (ATR-ATRIP) is 4.7 Å | RCSB / PMC | **Confirmed** — catalytic core built to 3.9 Å within the 4.7 Å map |
| PDB 9L40 (ATR bound VE-822) and 9L4B (bound RP-3500) exist, ~3 Å | RCSB | **Confirmed**, both real, correct ligands |
| Wang et al. 2025, *Sci Bull*, DOI 10.1016/j.scib.2025.05.009 | PubMed/EurekAlert search (direct fetch blocked by egress policy, corroborated via two independent search summaries) | **Confirmed** — real paper, USTC (Cai/Wang labs), ~3 Å ATR-ATRIP + VE-822/RP-3500 |
| One ATR-ATRIP complex binds 4 VE-822 molecules, 2 active-site + 2 at the dimer interface (allosteric) | Cross-checked via two independent search queries | **Confirmed** at the level of the core finding. Residue-level detail (R1959, specific TRD-domain boundaries) could **not** be independently confirmed — PubMed and EurekAlert are blocked by this session's egress policy — and is flagged as unconfirmed in the manuscript rather than asserted |
| IJMS 2026;27(16):7497, DOI 10.3390/ijms27167497 ("Beyond the Score...") | Direct search, PMC/PubMed listing | **Confirmed** real, correct DOI. Author names not independently confirmed — cited as "authors not verified" rather than guessed |
| arXiv:2507.21404, LIT-PCBA leakage audit | Direct search | **Confirmed** real. The specific number cited by the reviewer — 323 ALDH1 analog pairs at ECFP4 Tc≥0.6 — **matches exactly** what the search returned, which is strong evidence the reviewer read the primary source rather than pattern-matching a plausible-sounding citation |
| Sieg, Flachsenberg & Rarey, *JCIM* 2019, 59:947-961 ("In need of bias control") | Direct search | **Confirmed**, correct venue and volume |

Given every specific, checkable claim in this review held up — including one
exact number — the review is treated as reliable, and its remaining
(harder-to-verify-independently) technical judgments are acted on rather than
re-litigated one by one.

## Accepted and acted on

| # | Finding | Action |
|---|---|---|
| C4 | ATR receptor choice (5YZ0, 4.7 Å) unsuitable for docking; 9L40/9L4B (~3 Å) available | §2.1 rewritten: 5YZ0 explicitly excluded by design with the resolution stated; 9L40/9L4B named as the working default. CHK1 ID left `[PENDING]` rather than fabricated, since no specific CHK1 PDB ID was verified this session |
| C5 | Missing 2025 structure showing berzosertib's allosteric dimer-interface site | §1.2 rewritten to cite Wang et al. (2025) and state the finding; new §4.5 (Future Work) names investigating that site as the most concrete extension of this work, rather than only gesturing at "partial modulation" in general terms |
| M1 | Novelty claim ("no framework exists") false at stated generality | §1.3 narrowed to this specific target pair; LIT-PCBA, DeepCoy, the bias-control papers, and the 2026 fixed-budget benchmarking paper are now cited as the prior art this work instantiates. A `[PENDING]` flag added: the negative claim still lacks a documented systematic search (databases, date, query terms) |
| M2 | Scaffold splitting insufficient; no exclusion rule; BRICS/Murcko mismatch | `scripts/09_pipeline_validation.py`: default Tanimoto threshold lowered 0.7→0.4; new `exclude_leaked_compounds` actually removes near-duplicates rather than only reporting them (tested — see below). §2.3/§4.3 now state the BRICS/Murcko mismatch as a named, unresolved limitation rather than eliminated |
| M3 | Circular evaluation set (same bioactivity pull for seeds and "held-out" actives) | §2.5 adds an explicit "Scope of what this benchmark can show" paragraph: a good result demonstrates recovery of known inhibitor chemotypes, not novel-chemotype discovery. A temporally held-out set is named as a planned, not-yet-implemented extension |
| M4 | DUD-E decoys have documented bias | §2.5/§4.3 acknowledge this directly, citing Sieg et al. (2019) and Chen et al. (2019); `scripts/09_pipeline_validation.py` gets a new `compare_decoy_standards` function so a gap between decoy standards is reported rather than only the more favorable number. DeepCoy itself is not reproduced (requires its own trained model) — flagged as an external dependency |
| M5 | MM-PBSA/MM-GBSA inconsistency | Fixed everywhere: manuscript now says MM-GBSA as the *method* throughout (abstract, §2.8 heading, §3.7, figure captions), with an explicit naming note explaining `MMPBSA.py` remains the correct name of the *tool* regardless |
| M6 | No redocking validation control | New `redock_native_ligand` (+ `extract_cocrystal_ligand`) in `scripts/05_docking_pipeline.py`: redocks a receptor's own co-crystallized ligand, reports RMSD against a pre-specified ≤2.0 Å threshold, and is required to pass before any candidate is scored. Manuscript §2.4 also now explicitly distinguishes this (accuracy) from replicate-seed agreement (reproducibility) |
| M7 | PKA miscast as a distant negative control | §2.9 corrected: PKA is CHK1's AGC-group fold neighbor (both share the CAMK/AGC bilobal kinase fold; Manning et al., 2002), retained as a fold comparator, not a negative control. `scripts/08_kinase_selectivity_profiling.py` gets a `DISTANT_CONTROL` slot (currently empty, `[PENDING]`) for an actual fold-unrelated ATP-binding protein |
| M8 | n=3 reference compounds cannot validate a correlation | §2.9 states this directly; `benchmark_against_reference_profiles` now raises unless ≥20 compounds are supplied (`MIN_REFERENCE_COMPOUNDS`), with an explicit opt-out for a labelled illustrative-only smoke test |
| M9 | Library size never stated; EF1% uninterpretable without it | §2.2/§3.1/§3.5 flag every stage's N as `[PENDING]`, reported in a new Figure 4 (compound-flow diagram) — which also fixes m6 (Figure 4 was referenced nowhere) |
| M10 | "Pre-specified" asserted, not evidenced | §2.5/§2.12/Data Availability all now state a timestamped registration (Zenodo/OSF or a signed dated git tag) is required before the claim is checkable, `[PENDING]` until it exists |
| M11 | Correlated ranking objectives inflate effective weight | §2.10 adds: a stated decision rule (ρ>0.7 → redundancy-corrected variant reported alongside the original), a uniform-weight ranking, and the Pareto-optimal set, with rank agreement across all three reported |
| M12 | Disease framing not load-bearing; ref 11 (Cri du Chat) and Tier 3 should go; §2.11 promises what it can't deliver | Tier 3 removed entirely from §1.5; the Niebuhr 1978 citation removed; §2.11 rewritten to state plainly that gene span/AT-content/GC-skew cannot identify replication-timing-matched control loci, and that Repli-seq data is a named, not-yet-incorporated dependency for §4.4 |
| M13 | §4.2 engagement-range hypothesis uncited, load-bearing for §4.5 | Demoted to §4.5 (Future Work), explicitly labelled as uncited and unsupported by this manuscript; §4.4's translational path no longer depends on it as its sole motivation |
| m5 | Descriptor filters applied twice (§2.3 and old §2.6) at different stages | Left as-is functionally (pre-filtering before expensive docking, then a fuller pass after — this is a legitimate two-stage design, not a bug) but the manuscript no longer implies novelty from either stage alone |
| m7 | "Continuous count" for a 0-4 integer | Fixed to "discrete (0–4) rule-satisfaction count" in §2.6 |
| m8 | Vina exhaustiveness default too low | `scripts/05_docking_pipeline.py`: default raised 16→32 (`DEFAULT_EXHAUSTIVENESS`); manuscript §2.4 states this explicitly instead of leaving it to a bare `[PENDING]` |
| m9 | Ambiguous Lipinski year ("1997/2001") | Fixed to Lipinski et al., 1997, *Adv Drug Deliv Rev* 23:3-25 |
| — | Citation renumbering risk | Whole manuscript switched from numbered superscripts to author-year citations, specifically so future reference-list edits (which are still ongoing) don't require renumbering every in-text mark — a real risk given how much the reference list changed in this one round |

## Verified but not acted on further

- **M14 (implicit) — reference list still not fully verified.** The
  reviewer's own spot-check confirmed most citations were correctly
  attributed; refs 2 and 8 (previously vague placeholders) are replaced with
  real, specific citations (Huebner & Croce, 2001; Fokas et al., 2012) that
  still need primary-source confirmation like everything else not marked
  "verified." This is not resolved, only improved — the header note says so
  explicitly and a Limitations bullet now says so too, per the reviewer's
  point that this is a research-integrity matter, not formatting.

## Contested / deliberately not fully implemented

- **AiZynthFinder hard filter (m2).** Cited (Genheden et al., 2020) and
  flagged as a planned addition in §2.10, but not implemented in
  `scripts/10_integrated_ranking.py` this round — it requires an external
  retrosynthesis-planning dependency (AiZynthFinder itself, with its own
  trained models) beyond what this session's egress restrictions and time
  budget support integrating and testing properly. Marking it `[PENDING]`
  rather than stubbing in an untested integration was judged the more honest
  choice.
- **gmx_MMPBSA fallback (m4/implicit).** Noted in §2.8 as the fallback if the
  OpenMM→ParmEd→MMPBSA.py conversion proves unreliable, but not implemented
  as actual alternate code — there is nothing to test until the primary path
  is actually run once, which (per the draft-status header) hasn't happened
  in this environment.
- **DeepCoy decoys (M4).** `compare_decoy_standards` is added so a second
  decoy standard can be compared against DUD-E-style decoys, but generating
  DeepCoy decoys themselves requires the external DeepCoy model, which is not
  reproduced here. The function accepts pre-generated results from either
  standard rather than assuming how the secondary set was produced.
- **DISTANT_CONTROL kinase (M7).** Left as an empty, flagged slot rather than
  populated with a specific PDB ID, since picking one (e.g., an ATP-binding
  chaperone like Hsp90) and verifying its fold classification and a suitable
  structure was judged worth doing carefully in a follow-up pass rather than
  filling in under time pressure in this one.

## Net effect

No claim in the manuscript got stronger this round; several got narrower,
better-cited, or newly contingent on a `[PENDING]` verification step, and two
real code gaps (no docking-accuracy control, no actual leakage exclusion)
were closed with functions that were tested against synthetic data in this
session (not against real receptors/compounds, which this environment still
cannot reach). The reviewer's overall framing — this is closer to a strong
Stage 1 Registered Report protocol than to a submittable manuscript — is not
contested and is treated as the accurate current status.
