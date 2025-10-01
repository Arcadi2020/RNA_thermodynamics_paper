# Cell Cycle Transcriptome–Thermodynamics Analysis

## Overview
This repository contains R scripts designed to integrate **RNA secondary structure energetics** with **cell cycle transcriptomics** in human cancer cell lines.  
The overall goal is to understand how transcript stability (ΔG<sub>fold</sub>) and thermal buffering (heat capacity) relate to gene expression programs during different phases of the cell cycle.

## Scripts

### 1. `CellCycle_MCF7_HeLa_Thermodynamics.R`
- **Aim:** Analyze transcriptome-wide coupling of RNA folding energetics with gene expression changes across cell cycle phases.  
- **Datasets:**  
  - **MCF-7** (breast cancer, GSE94479, Liu et al. PNAS 2017)  
  - **HeLa** (cervical cancer, GSE104736, Boström et al. PLOS ONE 2017)  
- **Main steps:**  
  - Load differential expression data and join with RNAfold-derived metrics (MFE, heat capacity).  
  - Compare up- vs down-regulated transcripts across G0/G1, S, and G2/M phases.  
  - Compute **expression-weighted mean MFE** per phase.  
  - Perform **KEGG GSEA** (clusterProfiler) to identify enriched signaling pathways.  
- **Outputs:** Plots linking RNA structure with expression, weighted thermodynamic summaries per phase, and pathway enrichment results.  

### 2. `CellCycle_ExpressedGenes_GSE94479.R`
- **Aim:** Identify genes expressed in each phase of the cell cycle (G0/G1, S, M) using BigWig RNA-seq coverage files from GSE94479.  
- **Datasets:**  
  - BigWig coverage tracks from GEO (G0/G1, S, M phases).  
  - Gene annotations (refFlat, hg19).  
  - Precomputed DESeq2 contrasts and normalized GRO-seq counts.  
- **Main steps:**  
  - Import BigWig files and overlap with gene annotations to extract phase-specific gene lists.  
  - Integrate precomputed differential expression results and GRO-seq counts.  
  - Map RefSeq IDs to Ensembl gene/transcript IDs (biomaRt).  
- **Outputs:** Gene lists per phase, annotated DE tables with gene IDs, and GRO-seq–integrated datasets.  

## Purpose
Together, these scripts provide a **framework to connect transcriptome remodeling with RNA structural stability** across the human cell cycle, enabling comparisons across cell lines (MCF-7 vs HeLa) and integration of expressed gene sets with thermodynamic and pathway-level insights.

