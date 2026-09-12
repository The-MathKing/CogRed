# Journal Targeting Strategy

All five are peer-reviewed and either fully open-access or open-access-optioned,
and all regularly publish pure-computational CADD studies (no wet-lab data
required). Ranked roughly by fit for this specific paper.

## 1. Journal of Chemical Information and Modeling (JCIM, ACS) — primary target
- **Audience:** Computational/medicinal chemists, cheminformaticians; the
  disciplinary home for docking/MD/ADMET methods papers.
- **OA status:** Hybrid — OA optional via ACS AuthorChoice (fee applies); free
  to publish subscription route if OA isn't required.
- **Formatting:** ACS format, ~6000-8000 words + SI; strict on methods
  reproducibility (they will desk-reject or require revision if version
  numbers, seeds, and grid box parameters are missing).
- **Benchmarks this paper must hit:** enrichment/validation of the docking
  pipeline on known actives vs. decoys (Section 2.5/3.2 above) is close to
  mandatory here; MM-PBSA/MM-GBSA with proper error bars; full SI compound
  table. This is the most technically demanding reviewer pool of the five —
  strongest fit if the MD+selectivity additions are actually completed.

## 2. PLOS Computational Biology — strong secondary target
- **Audience:** Broad computational biology, mechanistic/systems framing valued
  over pure method novelty; reviewers weigh biological significance heavily.
- **OA status:** Fully open access (APC ~$3000, fee waivers available).
- **Formatting:** No strict word limit; structured abstract optional; strong
  preference for a clear "Author Summary" (plain-language paragraph) — this
  suits the disease-relevance narrative (fragile sites, genomic instability)
  well if framed carefully per the Tier 1/Tier 2 distinction in the critique.
- **Benchmarks:** the genomic feature analysis (§2.10/3.7) and the mechanistic
  narrative matter more here than at JCIM; still expect rigorous methods
  reporting and code availability (they require a public repository — this
  one qualifies).

## 3. Scientific Reports (Nature Portfolio)
- **Audience:** Broad, multidisciplinary; higher volume/faster turnaround than
  the above two; good fallback or simultaneous-tier option.
- **OA status:** Fully open access (APC ~$2390, waiver program).
- **Formatting:** Flexible IMRaD, generally more lenient on scope/length;
  editorial bar for "significance" is lower than PLOS Comp Biol, but technical
  soundness review is still real (statistics, reproducibility).
- **Benchmarks:** Full pipeline is not strictly required to pass review, but
  omitting MD/selectivity entirely would likely draw a "computational
  validation is limited to docking" reviewer comment — same core additions
  still recommended.

## 4. Frontiers in Pharmacology (Experimental Pharmacology and Drug Discovery
   section)
- **Audience:** Pharmacology-focused, receptive to target-mechanism-centered
  CADD papers, including kinase modulator discovery.
- **OA status:** Fully open access (APC ~$2950).
- **Formatting:** Structured, section-based Frontiers template; interactive
  review process (reviewers and authors correspond directly), which can be an
  advantage for a methodologically ambitious but first-time-in-this-subfield
  submission.
- **Benchmarks:** Selectivity profiling against kinase off-targets is
  particularly well-received in this venue given its pharmacology audience;
  make sure the "modulator vs. inhibitor" pharmacological rationale (critique
  §2) is stated clearly since this audience will scrutinize mechanism claims.

## 5. ACS Omega
- **Audience:** Broad chemistry, higher acceptance of preliminary/exploratory
  computational work; good option if compute constraints limit the study to
  docking + ADMET + a smaller MD subset.
- **OA status:** Fully open access (APC ~$2500, lower than most ACS journals).
- **Formatting:** ACS format, shorter/faster review cycle than JCIM; less
  demanding on enrichment/validation benchmarks.
- **Benchmarks:** This is the realistic fallback if MD/MM-PBSA and kinome
  selectivity can't be completed at full scale before submission — but the
  paper should still be explicit about that scope limitation in Discussion
  rather than silently omitting the analyses.

## Recommended sequencing
1. Build the full pipeline (docking + ADMET + MD/MM-PBSA + selectivity +
   genomic feature analysis) targeting **JCIM** as the primary submission.
2. If compute/time constraints force cutting MD/MM-PBSA down to a small
   shortlist only, retarget to **PLOS Computational Biology** or **Frontiers
   in Pharmacology**, which weigh mechanistic/biological narrative alongside
   (not solely) computational exhaustiveness.
3. Keep **Scientific Reports** and **ACS Omega** as fallback venues if the
   primary submission is rejected on scope grounds rather than technical
   grounds — both have faster cycles and would not require re-scoping the
   science.
