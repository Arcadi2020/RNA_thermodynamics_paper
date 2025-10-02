# This script implements the ssGSEA function to assign a proliferation score to the transcriptonal profile of a specific sample.
# The proliferation signature has been obtained from:
# Jiang, Z. et al. 'A conserved expression signature predicts growth rate and reveals cell & lineage-specific differences', PLOS COMPUTATIONAL BIOLOGY 2021


# Aim
# Compute single-sample GSEA (ssGSEA)–based proliferation scores for one or more
# RNA-seq expression profiles, using conserved proliferation signatures
# (Carey 2021; Science paper set from All_gene_sets_AC.R).
#
# What this script provides
# - A self-contained ssGSEA implementation (rank-based enrichment per sample).
# - Utilities to build gene sets from publications:
#   * Carey 2021 list (via Entrez → Ensembl → Gene Symbol mapping with biomaRt).
#   * Proliferation_Science set loaded from All_gene_sets_AC.R.
# - A unified gene-set object you can pass to ssGSEA or GSVA.
#
# Inputs (set paths in RStudio: Session → Set Working Directory)
# - Expression matrix X: genes (rows) × samples (cols), rownames = Gene Symbols.
# - "Proliferation_signature_carey_2021.csv" (Entrez IDs) and
#   "Proliferation_Carey_2021_gene_names.csv" (gene symbols) OR use All_gene_sets_AC.R.
# - "All_gene_sets_AC.R" to access Proliferation_Science_paper.
#
# Core workflow
# 1) Define ssGSEA():
#    - Ranks genes per sample, computes running-sum enrichment against each set.
#    - Key args: alpha (weighting), scale (normalize by gene count), norm (range),
#      single (classic ssGSEA sum vs max deviation).
# 2) Build proliferation gene sets:
#    - Read Carey 2021; map IDs with biomaRt; deduplicate symbols.
#    - Load Proliferation_Science from All_gene_sets_AC.R.
#    - Optionally combine into a single data.frame/list "gene_sets".
# 3) Score samples:
#    - Ensure rownames(X) are Gene Symbols matching the gene sets.
#    - Run: es <- ssgsea(X, gene_sets, alpha = 0.25, scale = TRUE, norm = FALSE).
#    - (Alternative) GSVA: gsva(X, gene_sets, method = "ssgsea", …).
# 4) Visualize / export (optional):
#    - Heatmap of enrichment scores (ComplexHeatmap / circlize).
#    - Write matrix to CSV for downstream stats or plotting.
#
# Outputs
# - Matrix of enrichment scores: rows = gene sets; cols = samples.
# - Optional heatmap objects and CSV exports.
#
# Practical notes & pitfalls
# - ID consistency: rownames(X) must be the same ID type as the gene sets
#   (this script standardizes to Gene Symbols). Re-map if needed.
# - NA handling: remove genes with all-NA or zero variance before ranking.
# - Duplicates: if multiple rows share the same symbol, aggregate (e.g., max/mean).
# - Alpha & scaling: alpha = 0.25 is standard; set scale = TRUE for sample-size
#   comparability; set norm = TRUE only if you need min–max range normalization.
#
# Minimal example (after preparing X and gene_sets)
#   es <- ssgsea(X, gene_sets, alpha = 0.25, scale = TRUE, norm = FALSE)
#   ComplexHeatmap::Heatmap(es, name = "ssGSEA_ES")
#
# Reuse
# - Swap in any expression matrix (bulk or single-sample), keep symbols consistent.
# - Drop/replace signatures as needed; the ssGSEA() is general-purpose.



library(matrixStats)
library(circlize)
library(ComplexHeatmap)
library(data.table)
library(dplyr)
library(tidyverse)

# SSGEA function
ssgsea = function(X, gene_sets, alpha = 0.25, scale = T, norm = F, single = T) {
  row_names = rownames(X)
  num_genes = nrow(X)
  gene_sets = lapply(gene_sets, function(genes) {which(row_names %in% genes)})
  
  # Ranks for genes
  R = matrixStats::colRanks(X, preserveShape = T, ties.method = 'average')
  
  # Calculate enrichment score (es) for each sample (column)
  es = apply(R, 2, function(R_col) {
    gene_ranks = order(R_col, decreasing = TRUE)
    
    # Calc es for each gene set
    es_sample = sapply(gene_sets, function(gene_set_idx) {
      # pos: match (within the gene set)
      # neg: non-match (outside the gene set)
      indicator_pos = gene_ranks %in% gene_set_idx
      indicator_neg = !indicator_pos
      
      rank_alpha  = (R_col[gene_ranks] * indicator_pos) ^ alpha
      
      step_cdf_pos = cumsum(rank_alpha)    / sum(rank_alpha)
      step_cdf_neg = cumsum(indicator_neg) / sum(indicator_neg)
      
      step_cdf_diff = step_cdf_pos - step_cdf_neg
      
      # Normalize by gene number
      if (scale) step_cdf_diff = step_cdf_diff / num_genes
      
      # Use ssGSEA or not
      if (single) {
        sum(step_cdf_diff)
      } else {
        step_cdf_diff[which.max(abs(step_cdf_diff))]
      }
    })
    unlist(es_sample)
  })
  
  if (length(gene_sets) == 1) es = matrix(es, nrow = 1)
  
  # Normalize by absolute diff between max and min
  if (norm) es = es / diff(range(es))
  
  # Prepare output
  rownames(es) = names(gene_sets)
  colnames(es) = colnames(X)
  return(es)
}




# Open the proliferation signature in Carey 2021
Carey_2021 <- read.delim("path_to/Proliferation_signature_carey_2021.csv", sep = ",", header = T, stringsAsFactors = F)
# For each gene symbol identify the Ensembl canonical transcript ID
library(biomaRt)
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host="https://www.ensembl.org")
transcript_ids <- (Carey_2021$entrezgene_id) # This reads the gene symbols in the Cell_proliferation set

res_Carey_2021 <- getBM(attributes = c('ensembl_transcript_id_version', 
                                       'ensembl_gene_id',
                                       'ensembl_gene_id_version',
                                       'external_gene_name',
                                       'ensembl_transcript_id',
                                       'entrezgene_id',
                                       'transcript_is_canonical',
                                       'external_gene_name'),
                        filters = 'entrezgene_id', 
                        values = transcript_ids,
                        mart = mart)

Proliferation_Carey_2021 <- dplyr::select(res_Carey_2021, external_gene_name, ensembl_gene_id_version, ensembl_gene_id, entrezgene_id) %>%
  dplyr::distinct(external_gene_name, .keep_all = FALSE) %>%
  dplyr::rename(Proliferation_Carey = external_gene_name)
Proliferation_Carey_2021 <- as.data.frame(Proliferation_Carey_2021)


# Get the vector with the gene list
Proliferation_Carey_2021 <- read.delim(file="path_to/Proliferation_Carey_2021_gene_names.csv", sep = ",",
                                  header = T, stringsAsFactors = F)

Proliferation_Carey_2021 <- dplyr::pull(Proliferation_Carey_2021, Proliferation_Carey)

# Add the quotes and the coma to each name of the vector
Proliferation_Carey_2021 <- paste0('"', gsub('"', '""', Proliferation_Carey_2021), '"', collapse = ',')
cat(Proliferation_Carey_2021, '\n')




# Open the gene_set file with the proliferation signature in the Science paper
source("path_to/All_gene_sets_AC.R")
Proliferation_Science <- Proliferation_Science_paper
Proliferation_Science <- as.data.frame(Proliferation_Science)


# I joined the two signatures in a single .csv file
Prolif_signatures <- read.delim("path_to/Proliferation_signatures.csv", sep = ",", header = T, stringsAsFactors = F)
gene_sets <- Prolif_signatures


# Assign proliferation score

# Use GSVA to assign the proliferation score
# res_GSE182018_GSVA <- gsva(dds_gene_rl_matrix, Prolif_signatures, verbose=FALSE)

# Use ssGSEA to assign the proliferation score
# res_GSE182018_ssGSEA <- (ssgsea(dds_gene_rl_matrix, Prolif_signatures, scale = TRUE, norm = FALSE))









