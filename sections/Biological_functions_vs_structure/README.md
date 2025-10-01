RNA Biological Functions and Secondary Structure Stability
Overview

This section of the repository accompanies the manuscript analyses showing how RNA biological functions correlate with the stability of secondary structures. The work integrates thermodynamic modeling of RNA folding with functional and evolutionary biology, revealing how structural stability encodes a conserved regulatory dimension across transcript biotypes, tissues, signaling pathways, and evolutionary lineages.

Biological Rationale

RNA molecules spontaneously fold into ensembles of alternative conformations, each defined by a folding free energy (ΔG_fold). The probability of a given conformation follows Boltzmann statistics, with lower ΔG_fold values indicating more stable structures. While ensemble properties can in principle be captured by partition functions, transcriptome-scale studies typically approximate stability using the Minimum Free Energy (MFE) structure.

Stable RNAs (more negative ΔG_fold): compact, energetically favorable, and resistant to degradation or thermal denaturation; common in tRNAs, rRNAs, and metabolically demanding processes such as extracellular matrix remodeling.

Less stable RNAs (less negative ΔG_fold): conformationally flexible, higher entropy, and thermodynamically reversible; enriched in stress-response and error-prone repair pathways.

These properties link structural stability to cellular states:

Proliferation: favors stable transcripts to sustain heat dissipation and maintain irreversibility.

Stress/adaptation: favors less stable transcripts, enhancing buffering capacity and reversibility.

Computational Workflow

The scripts in this folder perform genome-wide secondary structure predictions and downstream functional analyses:

Transcriptome-wide MFE calculation

Predicted with ViennaRNA (RNAfold).

Normalized per nucleotide length (ΔG_fold / nt).

Transcript annotations from GENCODE.

Biotype-level comparisons

Protein-coding mRNAs vs. tRNAs vs. rRNAs (Kolmogorov–Smirnov tests).

Protein-coding transcripts display less negative ΔG_fold compared to tRNAs and rRNAs.

Extremes of stability

5th vs. 95th percentile transcripts analyzed by GSEA.

Less stable RNAs enriched for DNA repair and chromosome segregation.

Highly stable RNAs enriched for migration and extracellular matrix remodeling.

Tissue specificity

ΔG_fold values integrated with Human Protein Atlas expression data.

More stable RNAs enriched in high-metabolic, higher-temperature tissues (brain, liver, muscle).

Less stable RNAs enriched in lower-temperature tissues (testis, retina).

Pathway-level analyses

Curated pathways and MSigDB gene sets (KEGG, Reactome, GO).

Stability signatures distinguish DNA repair processes:

BER transcripts: highly stable (–0.391 ± 0.053 kcal/mol).

HR transcripts: less stable (–0.304 ± 0.069 kcal/mol).

Provides mechanistic rationale for synergy of hyperthermia + PARP inhibition.

Stress-responsive pathways

Error-prone DNA polymerases and NHEJ transcripts show less negative ΔG_fold, consistent with adaptive flexibility.

Evolutionary analysis

Phylostratigraphy (Trigos et al., PNAS 2017) integrated with ΔG_fold values.

Primitive unicellular phyla: highly negative ΔG_fold (stable).

Multicellular/metazoan lineages: progressively less negative ΔG_fold, reflecting adaptive plasticity.

Figures Generated

Fig. 1A: Biotype-level ΔG_fold distributions.

Fig. 1B: Structural exemplars from extreme tails of the distribution.

Fig. 1C: Tissue-specific ΔG_fold profiles.

Fig. 1D: Pathway-level ΔG_fold values (curated sets).

Fig. 1E: Error-prone vs. high-fidelity DNA repair pathways.

Fig. 1F: Evolutionary trajectory of RNA structural stability.

Suppl. Fig. 1–2: GSEA validation and extended pathway analyses.
