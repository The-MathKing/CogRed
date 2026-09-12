# Manuscript Outline

Target format: standard IMRaD, sized for JCIM / PLOS Computational Biology /
Scientific Reports (see `03_journal_targeting_strategy.md` for per-journal
adjustments). Section numbers map to files/functions in `scripts/` where
applicable so results are traceable to code.

## Title (working)
*"Structure-Based Optimization and Selectivity Profiling of ATR/CHK1
Pathway Modulators to Suppress Replication Stress at Common Fragile Sites:
An Open-Source Computational Pipeline"*

(Note: title drops the direct "Cri du Chat" claim per the rigor critique —
keep the disease-relevance framing for the Introduction/Discussion, scoped as
in §1 of the critique doc.)

## Abstract (250 words)
- 1 sentence: replication stress at common fragile sites (CFS) as a driver of
  genomic instability; ATR-CHK1 as the checkpoint axis.
- 1-2 sentences: gap — most existing ATR/CHK1 chemical matter is inhibitory
  (oncology-oriented); need for modulators that stabilize rather than abolish
  checkpoint signaling.
- 2-3 sentences: pipeline summary (ChEMBL seed set → derivative generation →
  ensemble docking → consensus scoring → ADMET/tox triage → MD/MM-PBSA on
  shortlist → kinome selectivity panel).
- 1-2 sentences: headline quantitative results (placeholder until data exists).
- 1 sentence: significance / open-science contribution (fully reproducible,
  code released).

## 1. Introduction
1.1 Replication stress and common fragile sites — biology of FRA3B/FHIT,
    FRA16D/WWOX, FRA7H; late replication timing, paucity of origins, AT-rich
    flexibility/secondary-structure motifs (cite Glover, Durkin, Casper,
    Debatisse & Rosselli reviews).
1.2 ATR-CHK1 axis — fork stabilization, origin firing suppression, S/G2-M
    checkpoint; distinguish physiological activation from chronic
    hyperactivation/senescence.
1.3 Disease relevance, precisely scoped — Tier 1 (CFS expression, cancer
    genomic instability, well-cited) and Tier 2 (hypothesis: fork-stability
    modulation as a contributor to non-recurrent germline microdeletion
    formation via FoSTeS/MMBIR, explicitly flagged as motivating future
    experimental work, Cri du Chat used as one illustrative example of this
    disease class).
1.4 Rationale for a modulator (not inhibitor) chemical strategy — the
    therapeutic-window argument from the critique doc.
1.5 Study objectives and pipeline overview (one summary figure).

## 2. Computational Methods
2.1 Target selection and structural preparation
   - PDB IDs for ATR and CHK1 (state resolution, apo/holo, co-crystal ligand,
     protonation protocol, missing-loop handling).
   - Binding site definition (co-crystal ligand-derived box vs. cavity
     detection tool).
2.2 Seed compound library construction
   - ChEMBL/PubChem programmatic pull of known ATR/CHK1 actives (bioactivity
     threshold stated, e.g., IC50 < 1 µM).
   - Deduplication, structure standardization (RDKit `MolStandardize`).
2.3 Derivative / analog generation
   - Method used (matched molecular pairs, BRICS recombination, R-group
     enumeration) with parameters.
   - Pre-filter: Lipinski/Veber/PAINS/Brenk (RDKit FilterCatalog) before
     docking to cut compute.
2.4 Molecular docking
   - AutoDock Vina/Webina version, exhaustiveness, number of independent runs,
     seed handling, box coordinates.
   - Ensemble docking against multiple conformers/structures (if used).
2.5 Consensus scoring and enrichment validation
   - Second scoring function used; DUD-E-style actives/decoys AUC/EF1% for
     the pipeline itself.
2.6 ADMET and toxicity triage
   - SwissADME parameters reported, ProTox-II endpoints reported, and/or
     RDKit-computed descriptor filters as a reproducible fallback; explicit
     pass/fail thresholds.
2.7 Molecular dynamics
   - OpenMM version, force fields (protein + ligand parametrization method),
     system size, solvation/ion model, minimization/equilibration/production
     protocol, number of replicates, simulation length, analysis metrics
     (RMSD/RMSF, ligand residence proxy, H-bond occupancy).
2.8 Binding free energy calculation
   - MM-PBSA/MM-GBSA tool and settings, number of frames used, entropy
     treatment (or explicit statement that entropy was neglected and why).
2.9 Off-target kinome selectivity panel
   - List of off-target kinases and PDB IDs chosen (rationale: PIKK family +
     cell-cycle kinases), same docking/scoring protocol applied, selectivity
     metric defined.
2.10 Fragile-site genomic feature analysis
   - Data sources (UCSC, HumCFS, public Repli-seq), features computed
     (replication timing, AT content, gene size, flexibility index), and how
     this connects back to the biological narrative (not a docking result —
     a supporting genomics analysis).
2.11 Statistical analysis and reproducibility statement
   - Software version table, hardware, random seeds, code/data availability
     (this repo + release DOI).

## 3. Results
3.1 Seed library and derivative set characterization (chemical space plot —
    PCA/UMAP over descriptors, scaffold diversity).
3.2 Docking enrichment validation (pipeline AUC/EF on known actives/decoys) —
    establishes the pipeline works before reporting novel hits.
3.3 Docking + consensus scoring results for derivative library — ranked
    table, top N hits.
3.4 ADMET/toxicity triage outcomes — how many survived, why compounds were
    excluded, property distributions before/after.
3.5 MD stability and MM-PBSA binding free energies for shortlisted hits —
    ΔG with error bars, per-residue energy decomposition for top 2-3 hits.
3.6 Selectivity profiling — heatmap of predicted affinity across target +
    off-target panel; selectivity score ranking.
3.7 Fragile-site genomic feature results — descriptive statistics tying
    FRA3B/FRA16D features to the mechanistic narrative (supports Introduction/
    Discussion, doesn't need to "prove" anything about the compounds).
3.8 Integrated multi-objective ranking (desirability function combining
    affinity, selectivity, ADMET, synthetic accessibility) → final candidate
    list (e.g., top 3-5 "lead-like modulators").

## 4. Discussion
4.1 Interpretation of top candidates relative to known ATR/CHK1 pharmacology
    (compare to berzosertib/ceralasertib/prexasertib/SRA737 as reference
    points, not templates being re-discovered).
4.2 Revisit the activation-vs-hyperactivation framing with the MM-PBSA/MD
    data — does anything in the results support the "mid-affinity window"
    hypothesis?
4.3 Selectivity liabilities and structural rationale (e.g., hinge-region
    similarity across PIKK family).
4.4 Explicit limitations: single/limited-conformer docking bias, MM-PBSA
    entropy approximation, no explicit cellular/phenotypic validation, the
    Tier 2 disease-relevance hypothesis remains untested in vivo.
4.5 Translational path — what wet-lab experiments would be the immediate next
    step (kinase assay, cellular replication-stress reporter, CFS expression
    assay under aphidicolin) — framing this as hypothesis-generating in
    silico work, consistent with journal expectations for pure computational
    studies.

## 5. Conclusion
- 3-4 sentences: what was found, what the pipeline contributes as an
  open-source resource, and the explicit next experimental step.

## Supporting Information
- Full compound table (SMILES, all scores) as CSV.
- Full software/version table.
- MD trajectory analysis plots for all shortlisted compounds (not just top
  hits).
- Code repository link + release tag/DOI.

## Figures (suggested, 5-6 main text)
1. Pipeline schematic (target → seed library → derivatives → docking →
   ADMET → MD/MM-PBSA → selectivity → ranked candidates).
2. Chemical space plot of seed vs. derivative library.
3. Enrichment (ROC) curve for pipeline validation.
4. Top-hit binding pose figure(s) with key interactions annotated.
5. MD RMSD/RMSF + MM-PBSA ΔG bar chart with error bars.
6. Selectivity heatmap across kinase panel.
