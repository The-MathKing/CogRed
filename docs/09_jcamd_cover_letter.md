# Cover Letter for JCAMD Submission

**Journal:** *Journal of Computer-Aided Molecular Design* (JCAMD)  
**Article Type:** Short Communication  
**Title:** Redocking Failure Is Common Across Recent Kinase Structures and Is Not Predicted by Map-Quality Metrics: A Validated Survey of ATR/CHK1 and CDK-Activating Kinase  

---

**To:**  
Editor-in-Chief  
*Journal of Computer-Aided Molecular Design*  

Dear Editor-in-Chief,

We are pleased to submit our manuscript titled **"Redocking Failure Is Common Across Recent Kinase Structures and Is Not Predicted by Map-Quality Metrics: A Validated Survey of ATR/CHK1 and CDK-Activating Kinase"** for consideration as a **Short Communication** in the *Journal of Computer-Aided Molecular Design*.

### Scientific Context and Key Findings

With the rapid influx of near-atomic cryo-EM models into the Protein Data Bank (PDB), structure-based drug discovery and virtual screening campaigns increasingly adopt these structures as docking receptors, often assuming that high global resolution ($\leq 3.2$ Å) and deposition recency ensure structural reliability for docking.

In this work, we conduct a rigorous, pre-registered redocking benchmarking survey across 17 recently deposited kinase structures (the ATR/CHK1 checkpoint kinases and a CDK-activating kinase [CAK] fragment-elaboration series):

1. **High Empirical Failure Rate:** Across 17 structures, over 80% (14/17, and 13/14 CAK structures) reproducibly fail a standard 2.0 Å symmetry-corrected RMSD redocking threshold, despite identical, objective receptor/ligand preparation.
2. **Independence from wwPDB Validation Metrics:** Neither nominal resolution nor wwPDB map-model validation metrics (Q-scores, residue inclusion) predict redocking success or failure ($p > 0.25$). Counterintuitively, lower-resolution structures can pass where higher-resolution counterparts in the same series fail.
3. **External Positive Control Validation:** To ensure this failure rate reflects structural properties rather than a defective preparation script, we subjected the identical, unmodified pipeline to an objective positive control from the CASF-2016 core set ($n=20$). The pipeline recovered a 40% pass rate ($p = 0.08$ vs. literature baseline), while the CAK failure rate remained a statistically significant outlier ($p = 0.008$), ruling out pipeline artifacts.
4. **Consequences for Virtual Screening:** In a matched pilot screen, scoring against an unvalidated failing structure versus a validated passing structure scrambled compound prioritization (Spearman $\rho = 0.21$, $p = 0.19$, with rank shifts up to 29 positions out of 40).

### Journal Scope & Format Alignment

This manuscript is submitted as a **Short Communication**, presenting a timely, methodologically disciplined cautionary benchmark directly relevant to the computational chemistry and molecular modeling community. 

- **Open Science & Reproducibility:** All code, notebooks, raw PDBQTs, and validation tables are archived in a public GitHub repository with a permanent release DOI.
- **Exclusivity:** This manuscript is original, has not been published previously, and is not currently under consideration elsewhere.

Thank you for your consideration of our work.

Sincerely,

Ethan Ye, Aryan Padarthi, Sungju Kim, Alexandre Elie-Dit-Cosaque  
Allen High School, Allen ISD, Allen, TX 75013, USA  
Corresponding author email: ethan.ye@student.allenisd.org, aryan.padarthi@student.allenisd.org  
