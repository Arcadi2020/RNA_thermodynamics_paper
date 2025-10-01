RNA Biological Functions and Secondary Structure Stability
Overview

This section of the repository accompanies the manuscript analyses showing how RNA biological functions correlate with the stability of secondary structures. The work integrates thermodynamic modeling of RNA folding with functional and evolutionary biology, revealing how structural stability encodes a conserved regulatory dimension across transcript biotypes, tissues, signaling pathways, and evolutionary lineages.

Stable RNAs (more negative ΔG_fold): compact, energetically favorable, and resistant to degradation or thermal denaturation; common in tRNAs, rRNAs, and metabolically demanding processes such as extracellular matrix remodeling.
Less stable RNAs (less negative ΔG_fold): conformationally flexible, higher entropy, and thermodynamically reversible; enriched in stress-response and error-prone repair pathways.

These properties link structural stability to cellular states:
Proliferation: favors stable transcripts to sustain heat dissipation and maintain irreversibility.
Stress/adaptation: favors less stable transcripts, enhancing buffering capacity and reversibility.

Computational Workflow
The script in this folder performs Biotype-level comparisons:
Protein-coding mRNAs vs. tRNAs vs. rRNAs (Kolmogorov–Smirnov tests).
Protein-coding transcripts display less negative ΔG_fold compared to tRNAs and rRNAs.
Extremes of stability
5th vs. 95th percentile transcripts analyzed by GSEA.

Tissue specificity
ΔG_fold values integrated with Human Protein Atlas expression data.

Pathway-level analyses
Curated pathways and MSigDB gene sets (KEGG, Reactome, GO).
Stability signatures distinguish DNA repair processes:

Stress-responsive pathways

Evolutionary analysis
Phylostratigraphy (Trigos et al., PNAS 2017) integrated with ΔG_fold values.


Figures Generated
Fig. 1A: Biotype-level ΔG_fold distributions.
Fig. 1B: Structural exemplars from extreme tails of the distribution.
Fig. 1C: Tissue-specific ΔG_fold profiles.
Fig. 1D: Pathway-level ΔG_fold values (curated sets).
Fig. 1E: Error-prone vs. high-fidelity DNA repair pathways.
Fig. 1F: Evolutionary trajectory of RNA structural stability.
Suppl. Fig. 1–2: GSEA validation and extended pathway analyses.
