
# Aim
# Build the figure panels for the manuscript by (i) assembling gene-level log2FC
# across many drug/cell-line datasets, (ii) running GSEA (KEGG + custom
# proliferation sets), and (iii) generating heatmaps, dot plots, violins, and
# bar charts summarizing MFE and pathway activity (including “Science models”).
#
# Inputs (set paths in RStudio: Session → Set Working Directory)
# - Precomputed DESeq2 objects/data frames loaded via:
#     load("path_to/Cross_validation_MFE_AC_environment.RData")
#   (provides res_* data frames, MFE_* summaries, *_counts_fraction, etc.)
# - Per-dataset tab-delimited log2FC tables for “Science models” (72h & EP).
# - KEGG/GSEA dependencies: clusterProfiler, enrichplot; ID mapping via biomaRt.
#
# Core workflow
# 1) Collect log2FC per dataset:
#    - Select (Gene_symbol, log2FoldChange) from each res_* table.
#    - Read additional 72h/EP Nanopore exports; average duplicates by gene.
#    - Full-join into Expression_changes_targeted_agents (wide matrix).
#
# 2) Prepare IDs for GSEA:
#    - Map Gene Symbol → Entrez (biomaRt).
#    - Merge Entrez IDs into the log2FC matrix.
#
# 3) KEGG GSEA per dataset:
#    - For each log2FC vector, sort decreasing, run gseKEGG(hsa).
#    - Save per-sample dotplots (PDF) and return combined table (NES, p.adjust).
#    - Filter to Selected_KEGG_pathways; rescale NES to ±4; draw heatmaps.
#
# 4) Select “cell-cycle arrest with no mitogenic signaling” samples:
#    - Start from pathways with NES<0 in hsa04110 (Cell cycle).
#    - Exclude samples with NES≥0 in any pro-mitogenic set.
#    - Heatmap of filtered samples (pathways × samples).
#
# 5) Proliferation signatures (custom GSEA):
#    - Build TERM2GENE for 3 signatures (Signt_1/2/3).
#    - Run GSEA (clusterProfiler::GSEA) on the same ranked vectors.
#    - Combine with KEGG results; bubble (dot) plot of NES vs –log10(p.adjust).
#
# 6) Violin plots of Normalized_MFE vs sign(log2FC):
#    - For selected GREIN datasets and Science models (72h, EP).
#    - Add companion stats: KS and Wilcoxon tests; export small summary tables.
#
# 7) Fraction plots (high/low MFE as % of total counts):
#    - Use precomputed fractions per condition; two-proportion Z-tests.
#    - Bar charts (control vs treated) for GREIN and Science models.
#
# 8) Science models comparisons (72h vs EP):
#    - Combine signature-level log2FC; boxplots with means.
#    - Paired/unpaired Wilcoxon tests across time points; save PDFs.
#
# Outputs
# - Per-sample KEGG dotplots (PDFs), combined GSEA CSVs, filtered heatmaps,
#   bubble plot of pathways, violin panels, and bar charts of high/low MFE
#   fractions; summary CSVs of NES/p.adjust and test p-values.
#
# Usage notes
# - Object names are assumed from the loaded .RData (res_*, MFE_res_*, *_fraction,
#   *_dexseq_ANNOTATED, etc.). If you rename upstream objects, update selectors.
# - Ensure Gene Symbol consistency before joins and GSEA; remove duplicates.
# - Long path strings include escaped spaces; adjust to your filesystem.




library(DESeq2)
library(GEOquery)
library(canvasXpress)
library(ggplot2)
library(clinfun)
library(GGally)
library(factoextra)
library(dplyr)
library(tidyverse)
library(RColorBrewer)
library(gplots)
library(GSVA)
library(RMySQL)
library(ggrepel)
library(ggpubr)
library(decoupleR)
library(OmnipathR)
library(tibble)
library(tidyr)
library(pheatmap)
library(enrichR)
library(biomaRt)
library(gridExtra)


# It is also possible to load the environment of the script Cross_validation_MFE.R
load("path_to/Cross_validation_MFE_AC_environment.RData")



# Create a data frame with the gene-level DESeq2 log_FC values of all the cell lines exposed to targeted agents.
# Take all the log_FC values (not only the transcripts with a significant log_FC)

# Generate the data frame with the gene-level DESeq2 log_FC values of all the cell lines exposed to single targeted agents:
# "GSE137535_SH454", "GSE143462_KYSE150_omipalisib", "GSE148096_TNBC_erlotinib", "GSE178755_HCC827_erlotinib",
# "GSE103021_LU65_ARS_1620",
# "GSE118548_T84_Trametinib", "GSE177054_MDAMB231_Ribociclib",
# "GSE93124_SW480_crizotinib",
# "GSE135514_MDAMB468_Das_24h", "GSE135514_MDAMB468_Das_72h",
# "GSE121105_BT474_trastuzumab", 
# "GSE105083_MIA_PACA2_CM03_6h", "GSE105083_MIA_PACA2_CM03_24h", "GSE105083_PANC1_CM03_24h", 
# "X778_nutlin3a_72h", "X778_palbo_72h", "X778_tunic_72h", "SKBR3_LAP_72h", "SKMEL28_VEM_72h", "X778_palbo_EP", "SKBR3_LAP_72h", "SKMEL28_VEM_EP"



# Load the DESeq2 data at the gene level
Expression_changes_GSE137535_SH454 <- dplyr::select(res_GSE137535_gene_SH454_df, Gene_symbol, log2FoldChange) %>%
  dplyr::rename(GSE137535_SH454_Log_FC = log2FoldChange)

Expression_changes_GSE143462_KYSE150_omipalisib <- dplyr::select(res_GSE143462_gene_df, Gene_symbol, log2FoldChange) %>%
  dplyr::rename(GSE143462_KYSE150_omipalisib_Log_FC = log2FoldChange)

Expression_changes_GSE148096_TNBC_erlotinib <- dplyr::select(res_GSE148096_gene_df, Gene_symbol, log2FoldChange) %>%
  dplyr::rename(GSE148096_TNBC_erlotinib_Log_FC = log2FoldChange)

Expression_changes_GSE178755_HCC827_erlotinib <- dplyr::select(res_GSE178755_gene_df, Gene_symbol, log2FoldChange) %>%
  dplyr::rename(GSE178755_HCC827_erlotinib_Log_FC = log2FoldChange)

Expression_changes_GSE103021_LU65_ARS_1620 <- dplyr::select(res_GSE103021_ARS_1620_LU65_gene_df, Gene_symbol, log2FoldChange) %>%
  dplyr::rename(GSE103021_LU65_ARS_1620_Log_FC = log2FoldChange)

Expression_changes_GSE118548_T84_Trametinib <- dplyr::select(res_GSE118548_Trametinib_gene_df, Gene_symbol, log2FoldChange) %>%
  dplyr::rename(GSE118548_T84_Trametinib_Log_FC = log2FoldChange)

Expression_changes_GSE93124_SW480_crizotinib <- dplyr::select(res_GSE93124_gene_df, Gene_symbol, log2FoldChange) %>%
  dplyr::rename(GSE93124_SW480_crizotinib_Log_FC = log2FoldChange)

Expression_changes_GSE135514_MDAMB468_Das_24h <- dplyr::select(res_GSE135514_Das_24h_gene_df, Gene_symbol, log2FoldChange) %>%
  dplyr::rename(GSE135514_MDAMB468_Das_24h_Log_FC = log2FoldChange)

Expression_changes_GSE135514_MDAMB468_Das_72h <- dplyr::select(res_GSE135514_Das_72h_gene_df, Gene_symbol, log2FoldChange) %>%
  dplyr::rename(GSE135514_MDAMB468_Das_72h_Log_FC = log2FoldChange)

Expression_changes_GSE121105_BT474_trastuzumab <- dplyr::select(res_GSE121105_trastuzumab_gene_df, Gene_symbol, log2FoldChange) %>%
  dplyr::rename(GSE121105_BT474_trastuzumab_Log_FC = log2FoldChange)

Expression_changes_GSE105083_MIA_PACA2_CM03_6h <- dplyr::select(res_GSE105083_MIA_PACA2_gene_6h_df, Gene_symbol, log2FoldChange) %>%
  dplyr::rename(GSE105083_MIA_PACA2_CM03_6h_Log_FC = log2FoldChange)

Expression_changes_GSE105083_MIA_PACA2_CM03_24h <- dplyr::select(res_GSE105083_MIA_PACA2_gene_24h_df, Gene_symbol, log2FoldChange) %>%
  dplyr::rename(GSE105083_MIA_PACA2_CM03_24h_Log_FC = log2FoldChange)

Expression_changes_GSE105083_PANC1_CM03_24h <- dplyr::select(res_GSE105083_PANC1_gene_24h_df, Gene_symbol, log2FoldChange) %>%
  dplyr::rename(GSE105083_PANC1_CM03_24h_Log_FC = log2FoldChange)

Expression_changes_GSE129221_PC9GR <- dplyr::select(res_GSE129221_PC9GR_gene_df, Gene_symbol, log2FoldChange)  %>%
  dplyr::rename(GSE129221_PC9GR_apatinib_Log_FC = log2FoldChange)




# Science models 72h time point
Expression_changes_X778_PALBO <- read.delim(file="path_to/72h_778_palbociclib.tx_counts.both_conditions.treated_untreated_deseq2.tab", header = TRUE, sep = "\t", quote = "") %>%
  dplyr::rename(Gene_symbol = gene_name) %>%
  dplyr::rename(X778_PALBO_Log_FC = log2FoldChange) %>%
  dplyr::select(Gene_symbol, X778_PALBO_Log_FC) %>%
  group_by(Gene_symbol) %>%
  summarize(X778_PALBO_Log_FC = mean(X778_PALBO_Log_FC)) %>%
  dplyr::filter(!is.na(X778_PALBO_Log_FC))



Expression_changes_X778_TUNIC <- read.delim(file="path_to/72h_778_tunicamycin.tx_counts.both_conditions.treated_untreated_deseq2.tab", header = TRUE, sep = "\t", quote = "") %>%
  dplyr::rename(Gene_symbol = gene_name) %>%
  dplyr::rename(X778_TUNIC_Log_FC = log2FoldChange) %>%
  dplyr::select(Gene_symbol, X778_TUNIC_Log_FC) %>%
  group_by(Gene_symbol) %>%
  summarize(X778_TUNIC_Log_FC = mean(X778_TUNIC_Log_FC)) %>%
  dplyr::filter(!is.na(X778_TUNIC_Log_FC))


Expression_changes_X778_nutlin <- read.delim(file="path_to/72h_778_nutlin3a.tx_counts.both_conditions.treated_untreated_deseq2.tab", header = TRUE, sep = "\t", quote = "") %>%
  dplyr::rename(Gene_symbol = gene_name) %>%
  dplyr::rename(X778_nutlin_Log_FC = log2FoldChange) %>%
  dplyr::select(Gene_symbol, X778_nutlin_Log_FC) %>%
  group_by(Gene_symbol) %>%
  summarize(X778_nutlin_Log_FC = mean(X778_nutlin_Log_FC)) %>%
  dplyr::filter(!is.na(X778_nutlin_Log_FC))


Expression_changes_SKBR3_LAP <- read.delim(file="path_to/72h_SKBR3_lapatinib.tx_counts.both_conditions.treated_untreated_deseq2.tab", header = TRUE, sep = "\t", quote = "") %>%
  dplyr::rename(Gene_symbol = gene_name) %>%
  dplyr::rename(SKBR3_Log_FC = log2FoldChange) %>%
  dplyr::select(Gene_symbol, SKBR3_Log_FC) %>%
  group_by(Gene_symbol) %>%
  summarize(SKBR3_Log_FC = mean(SKBR3_Log_FC)) %>%
  dplyr::filter(!is.na(SKBR3_Log_FC))


Expression_changes_SKMEL28_VEM <- read.delim(file="path_to/72h_SKMEL28_vemurafenib.tx_counts.both_conditions.treated_untreated_deseq2.tab", header = TRUE, sep = "\t", quote = "") %>%
  dplyr::rename(Gene_symbol = gene_name) %>%
  dplyr::rename(SKMEL28_Log_FC = log2FoldChange) %>%
  dplyr::select(Gene_symbol, SKMEL28_Log_FC) %>%
  group_by(Gene_symbol) %>%
  summarize(SKMEL28_Log_FC = mean(SKMEL28_Log_FC)) %>%
  dplyr::filter(!is.na(SKMEL28_Log_FC))







# Science models early-phase time point (the nutlin3a and the tunicamycin conditions have not been sequenced during the early phase response)
Expression_changes_X778_PALBO_EP <- read.delim(file="path_to/early.phase_778_palbociclib.tx_counts.both_conditions.treated_untreated_deseq2.tab", header = TRUE, sep = "\t", quote = "") %>%
  dplyr::rename(Gene_symbol = gene_name) %>%
  dplyr::rename(X778_PALBO_EP_Log_FC = log2FoldChange) %>%
  dplyr::select(Gene_symbol, X778_PALBO_EP_Log_FC) %>%
  group_by(Gene_symbol) %>%
  summarize(X778_PALBO_EP_Log_FC = mean(X778_PALBO_EP_Log_FC)) %>%
  dplyr::filter(!is.na(X778_PALBO_EP_Log_FC))


Expression_changes_SKBR3_LAP_EP <- read.delim(file="path_to/early.phase_SKBR3_lapatinib.tx_counts.both_conditions.treated_untreated_deseq2.tab", header = TRUE, sep = "\t", quote = "") %>%
  dplyr::rename(Gene_symbol = gene_name) %>%
  dplyr::rename(SKBR3_EP_Log_FC = log2FoldChange) %>%
  dplyr::select(Gene_symbol, SKBR3_EP_Log_FC) %>%
  group_by(Gene_symbol) %>%
  summarize(SKBR3_EP_Log_FC = mean(SKBR3_EP_Log_FC)) %>%
  dplyr::filter(!is.na(SKBR3_EP_Log_FC))


Expression_changes_SKMEL28_VEM_EP <- read.delim(file="path_to/early.phase_SKMEL28_vemurafenib.tx_counts.both_conditions.treated_untreated_deseq2.tab", header = TRUE, sep = "\t", quote = "") %>%
  dplyr::rename(Gene_symbol = gene_name) %>%
  dplyr::rename(SKMEL28_EP_Log_FC = log2FoldChange) %>%
  dplyr::select(Gene_symbol, SKMEL28_EP_Log_FC) %>%
  group_by(Gene_symbol) %>%
  summarize(SKMEL28_EP_Log_FC = mean(SKMEL28_EP_Log_FC)) %>%
  dplyr::filter(!is.na(SKMEL28_EP_Log_FC))






# Complete expression changes data frame

Expression_changes_targeted_agents <- dplyr::full_join(Expression_changes_GSE137535_SH454, Expression_changes_GSE143462_KYSE150_omipalisib, by = "Gene_symbol") %>%
  dplyr::full_join(Expression_changes_GSE148096_TNBC_erlotinib, by = "Gene_symbol") %>%
  dplyr::full_join(Expression_changes_GSE178755_HCC827_erlotinib, by = "Gene_symbol") %>%
  dplyr::full_join(Expression_changes_GSE103021_LU65_ARS_1620, by = "Gene_symbol") %>%
  dplyr::full_join(Expression_changes_GSE118548_T84_Trametinib, by = "Gene_symbol") %>%
  dplyr::full_join(Expression_changes_GSE93124_SW480_crizotinib, by = "Gene_symbol") %>%
  dplyr::full_join(Expression_changes_GSE135514_MDAMB468_Das_24h, by = "Gene_symbol") %>%
  dplyr::full_join(Expression_changes_GSE135514_MDAMB468_Das_72h, by = "Gene_symbol") %>%
  dplyr::full_join(Expression_changes_GSE121105_BT474_trastuzumab, by = "Gene_symbol") %>%
  dplyr::full_join(Expression_changes_GSE105083_MIA_PACA2_CM03_6h, by = "Gene_symbol") %>%
  dplyr::full_join(Expression_changes_GSE105083_MIA_PACA2_CM03_24h, by = "Gene_symbol") %>%
  dplyr::full_join(Expression_changes_GSE105083_PANC1_CM03_24h, by = "Gene_symbol") %>%
  dplyr::full_join(Expression_changes_GSE129221_PC9GR, by = "Gene_symbol") %>%
  dplyr::full_join(Expression_changes_X778_PALBO, by = "Gene_symbol") %>%
  dplyr::full_join(Expression_changes_X778_TUNIC, by = "Gene_symbol") %>%
  dplyr::full_join(Expression_changes_X778_nutlin, by = "Gene_symbol") %>%
  dplyr::full_join(Expression_changes_SKBR3_LAP, by = "Gene_symbol") %>%
  dplyr::full_join(Expression_changes_SKMEL28_VEM, by = "Gene_symbol") %>%
  dplyr::full_join(Expression_changes_X778_PALBO_EP, by = "Gene_symbol") %>%
  dplyr::full_join(Expression_changes_SKBR3_LAP_EP, by = "Gene_symbol") %>%
  dplyr::full_join(Expression_changes_SKMEL28_VEM_EP, by = "Gene_symbol")









# Fig.2B GSEA cell proliferation and mitogenic signaling ####
# Supplementary Fig. 4 ####

# Analysis cell proliferation with clusterProfiler GSE KEGG
# Perform GSEA on the complete Log_FC values of each cell line
# Also in this case use all the Log_FC values as clusterProfiler requires all of them ranked
# The input for clusterProfiler must be a vector of entrez gene names


library(biomaRt)
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host="https://useast.ensembl.org")
gene_ids <- (Expression_changes_targeted_agents$Gene_symbol) # This reads the Ensembl transcripts IDs in the specific column of the 'Table_ANNOTATED file.

res_gene_entrez_id <- getBM(attributes = c('entrezgene_id',
                                           'external_gene_name'),
                            filters = 'external_gene_name',
                            values = gene_ids,
                            mart = mart)
res_gene_entrez_id <- dplyr::filter(res_gene_entrez_id, !is.na(entrezgene_id)) %>%
  dplyr::rename(Gene_symbol = external_gene_name)

Expression_changes_targeted_agents_genes <- dplyr::left_join(Expression_changes_targeted_agents, res_gene_entrez_id, by = "Gene_symbol")



# Generate a function which performs GSEA across all cell lines and saves the plots

library(clusterProfiler)
library(pathview)
library(enrichplot)


process_and_save_plots <- function(data_frame, log_fc_cols, plot_titles, output_dir) {
  # Create the output directory if it does not exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  plot_list <- list()  # Initialize a list to store the plots
  kk2_list <- list()  # Initialize a list to store the GSEA results
  
  for (i in seq_along(log_fc_cols)) {
    log_fc_col <- log_fc_cols[i]
    plot_title <- plot_titles[i]
    file_name <- paste0(output_dir, "/", plot_title, "_GSEA_KEGG.pdf")
    
    # Process data and generate plot
    all_values_entrez <- data_frame %>%
      dplyr::select(entrezgene_id, !!sym(log_fc_col)) %>%
      dplyr::filter(!is.na(!!sym(log_fc_col))) %>%
      dplyr::filter(!is.na(entrezgene_id))
    
    all_values_Log2_FC_entrez <- all_values_entrez[[log_fc_col]]
    names(all_values_Log2_FC_entrez) <- all_values_entrez$entrezgene_id
    
    # Remove duplicates
    all_values_Log2_FC_entrez <- all_values_Log2_FC_entrez[!duplicated(names(all_values_Log2_FC_entrez))]
    
    all_values_Log2_FC_entrez <- na.omit(all_values_Log2_FC_entrez)
    all_values_Log2_FC_entrez <- sort(all_values_Log2_FC_entrez, decreasing = TRUE)
    
    # Debugging: Print the number of genes before and after filtering
    cat("Sample:", plot_titles[i], "\n")
    cat("Number of genes:", length(all_values_Log2_FC_entrez), "\n")
    
    if (length(all_values_Log2_FC_entrez) == 0) {
      cat("No genes available for", plot_title, "\n")
      next
    }
    
    kk2 <- gseKEGG(all_values_Log2_FC_entrez,
                   organism = "hsa",
                   keyType = "kegg",
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "BH",
                   minGSSize = 10,
                   maxGSSize = 500)
    
    
    if (is.null(kk2) || nrow(kk2) == 0) {
      cat("No significant pathways found for", plot_title, "\n")
      next
    }
    
    # Add kk2 results to the list
    kk2_list[[plot_title]] <- kk2
    
    p <- dotplot(kk2, 
                 x = "GeneRatio",
                 color = "p.adjust",
                 showCategory = 50,
                 size = NULL,
                 split = ".sign",
                 font.size = 8,
                 title = plot_title,
                 orderBy = "x",
                 label_format = 100) + 
      facet_grid(~.sign) +
      theme(panel.spacing = unit(0.5, "cm"))
    
    # Save plot
    ggsave(file_name, plot = p, width = 8, height = 12)
  }
  
  # Combine all kk2 results into a single data frame
  combined_kk2 <- do.call(rbind, lapply(names(kk2_list), function(name) {
    df <- as.data.frame(kk2_list[[name]])
    df$Sample <- name
    return(df)
  }))
  
  return(combined_kk2)
}
# Define your data frame, log FC columns, and plot titles
data_frame <- Expression_changes_targeted_agents_genes
log_fc_cols <- c("GSE93124_SW480_crizotinib_Log_FC", "GSE178755_HCC827_erlotinib_Log_FC", 
                 "GSE137535_SH454_Log_FC", "GSE105083_PANC1_CM03_24h_Log_FC", "GSE118548_T84_Trametinib_Log_FC", "GSE148096_TNBC_erlotinib_Log_FC", 
                 "GSE143462_KYSE150_omipalisib_Log_FC", "GSE121105_BT474_trastuzumab_Log_FC",
                 "GSE103021_LU65_ARS_1620_Log_FC", "GSE135514_MDAMB468_Das_24h_Log_FC", "GSE135514_MDAMB468_Das_72h_Log_FC",
                 "GSE105083_MIA_PACA2_CM03_6h_Log_FC", "GSE105083_MIA_PACA2_CM03_24h_Log_FC", 
                 "GSE129221_PC9GR_apatinib_Log_FC",
                 "X778_nutlin_Log_FC", "X778_TUNIC_Log_FC", "X778_PALBO_Log_FC", "SKBR3_Log_FC", "SKMEL28_Log_FC", 
                 "X778_PALBO_EP_Log_FC", "SKBR3_EP_Log_FC", "SKMEL28_EP_Log_FC")

plot_titles <- c("GSE93124", "GSE178755", 
                 "GSE137535_SH454", "GSE105083_PANC1_CM03_24h", "GSE118548_tra", "GSE148096", 
                 "GSE143462", "GSE121105",
                 "GSE103021", "GSE135514_24h", "GSE135514_72h",
                 "GSE105083_MIA_PACA2_6h", "GSE105083_MIA_PACA2_24h",
                 "GSE129221_PC9GR_apatinib",
                 "X778_nutlin_72h", "X778_TUNIC_72h", "X778_PALBO_72h", "SKBR3_72h", "SKMEL28_72h", 
                 "X778_PALBO_EP", "SKBR3_EP", "SKMEL28_EP")

# Define the output directory
output_dir <- "path_to/GSEA_KEGG_clusterProfiler"

# Call the function to process and save the plots
# process_and_save_plots(data_frame, log_fc_cols, plot_titles, output_dir)

# Combine all the saved plots into a single supplementary figure
# library(gridExtra)
# library(grid)
# 
# # Arrange plots in 3 columns and 6 rows (modify as necessary based on the number of plots)
# combined_plot <- grid.arrange(
#   grobs = plot_list,
#   ncol = 3,
#   nrow = 6
# )
# 
# # Save the combined plot as a PDF
# ggsave(
#   filename = paste0(output_dir, "/Supplementary_Fig_2_combined_plots.pdf"), 
#   plot = combined_plot, 
#   width = 25, 
#   height = 40
# )

# Call the function to get also the combined kk2 results
combined_kk2_df <- process_and_save_plots(data_frame, log_fc_cols, plot_titles, output_dir)
# Totally 22 samples of which 14 are GREIN datasets






# Generate a heat map with the kk2 results
library(pheatmap)

# Additional code to generate the heatmap
# Filter only the pathways of interest

# Selected KEGG signaling pathways
Selected_KEGG_pathways <- c("hsa04150", "hsa04668", "hsa04657", "hsa03030", "hsa03410", "hsa03430", "hsa03420",
                            "hsa03320", "hsa03460", "hsa03440", "hsa04014", "hsa04066", "hsa04070",
                            "hsa04310", "hsa03082", "hsa04015", "hsa04625", "hsa04371", "hsa04390", "hsa04010", "hsa04622", "hsa04350",
                            "hsa04152", "hsa04926", "hsa04630", "hsa04115", "hsa04151", "hsa04910", "hsa04330", "hsa04012",
                            "hsa04110") # hsa04110 is the cell cycle gene set

combined_kk2_df_filtered <- dplyr::filter(combined_kk2_df, ID %in% Selected_KEGG_pathways) %>%
  dplyr::mutate(NES = scales::rescale(NES, to = c(-4, 4))) # set a specific range of values to have the 0 centered to white


generate_heatmap <- function(kk2_df) {
  # Prepare data for the heatmap
  kk2_df <- kk2_df %>%
    dplyr::select(Description, Sample, NES) %>%
    tidyr::spread(Sample, NES)
  
  # Replace NA values with a large number for better visualization in the heatmap
  kk2_df[is.na(kk2_df)] <- 0
  
  rownames(kk2_df) <- kk2_df$Description
  kk2_df <- kk2_df[,-1]  # Remove the Description column
  
  # Create a matrix for display numbers
  display_matrix <- as.matrix(kk2_df)
  display_numbers <- ifelse(display_matrix == 0, "", sprintf("%.2g", display_matrix)) # do not show values if they are 0
  
  # Generate heatmap
  pheatmap(as.matrix(kk2_df), 
           color = colorRampPalette(c("blue", "white", "red"))(50), 
           cluster_rows = TRUE, 
           cluster_cols = TRUE, 
           display_numbers = display_numbers, 
           number_format = "%.2g", 
           main = "Heatmap of KEGG signaling pathways")
}

# Generate the heatmap using the combined_kk2_df
Heatmap_cell_cycle <- generate_heatmap(combined_kk2_df_filtered)
print(Heatmap_cell_cycle)

# Save the plot
# pdf("/Users/z3514253/Desktop/Heatmap_cell_cycle.pdf", width = 9, height = 12)
# generate_heatmap(combined_kk2_df_filtered)
# dev.off()








# Select only samples with a negative NES in the cell cycle gene set
Cell_cycle_arrest_samples <- dplyr::filter(combined_kk2_df_filtered, ID == "hsa04110" & NES < 0) %>%
  dplyr::select(Sample)

Seleted_non_proliferating_samples <- dplyr::filter(combined_kk2_df_filtered, Sample %in% Cell_cycle_arrest_samples$Sample) %>%
  dplyr::mutate(NES = scales::rescale(NES, to = c(-4, 4))) # set a specific range of values to have the 0 centered to white

generate_heatmap <- function(kk2_df) {
  # Prepare data for the heatmap
  kk2_df <- kk2_df %>%
    dplyr::select(Description, Sample, NES) %>%
    tidyr::spread(Sample, NES)
  
  # Replace NA values with a large number for better visualization in the heatmap
  kk2_df[is.na(kk2_df)] <- 0
  
  rownames(kk2_df) <- kk2_df$Description
  kk2_df <- kk2_df[,-1]  # Remove the Description column
  
  # Create a matrix for display numbers
  display_matrix <- as.matrix(kk2_df)
  display_numbers <- ifelse(display_matrix == 0, "", sprintf("%.2g", display_matrix)) # do not show values if they are 0
  
  # Generate heatmap
  pheatmap(as.matrix(kk2_df), 
           color = colorRampPalette(c("blue", "white", "red"))(50), 
           cluster_rows = TRUE, 
           cluster_cols = TRUE, 
           display_numbers = display_numbers, 
           number_format = "%.2g", 
           main = "Heatmap of KEGG signaling pathways")
}

# Generate the heatmap using the combined_kk2_df
Heatmap_cell_cycle <- generate_heatmap(Seleted_non_proliferating_samples)
print(Heatmap_cell_cycle)





# Identify the samples with a significant down-regulation of cell cycle and no mitogenic signaling
# Check all the pathways that may foster cell proliferation
ids_to_check <- c(
  "hsa04150", "hsa04668", "hsa04657", "hsa03030", "hsa03410", "hsa03430", "hsa03420",
  "hsa03320", "hsa03460", "hsa03440", "hsa04014", "hsa04066", "hsa04070",
  "hsa04310", "hsa03082", "hsa04015", "hsa04625", "hsa04371", "hsa04390", "hsa04010", "hsa04622", "hsa04350",
  "hsa04152", "hsa04926", "hsa04630", "hsa04115", "hsa04151", "hsa04910", "hsa04330", "hsa04012"
)

# Identify samples with NES > 0 for the specified IDs
samples_to_exclude <- Seleted_non_proliferating_samples %>%
  dplyr::filter(ID %in% ids_to_check, NES >= 0) %>%
  dplyr::select(Sample) %>%
  dplyr::distinct()

# Remove samples with actime mitogenic signaling
Cell_cycle_arrest_mitogenic_signaling <- dplyr::filter(Seleted_non_proliferating_samples, !Sample %in% samples_to_exclude$Sample)


generate_heatmap <- function(kk2_df) {
  # Prepare data for the heatmap
  kk2_df <- kk2_df %>%
    dplyr::select(Description, Sample, NES) %>%
    tidyr::spread(Sample, NES)

  # Replace NA values with a large number for better visualization in the heatmap
  kk2_df[is.na(kk2_df)] <- 0


  rownames(kk2_df) <- kk2_df$Description
  kk2_df <- kk2_df[,-1]  # Remove the Description column

  # Create a matrix for display numbers
  display_matrix <- as.matrix(kk2_df)
  # display_numbers <- ifelse(display_matrix == 0, "", sprintf("%.2g", display_matrix)) # do not show values if they are 0

  # Generate heatmap
  pheatmap(as.matrix(kk2_df),
           color = colorRampPalette(c("blue", "white"))(50), 
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           display_numbers = TRUE,
           number_format = "%.2g",
           main = "Heatmap of KEGG pathways cell cycle arrest and no mitogenic signaling",
           fontsize_row = 12,
           fontsize_col = 12,
           fontsize_number = 0)
}

# Generate the heatmap using the combined_kk2_df
Heatmap_cell_cycle <- generate_heatmap(Cell_cycle_arrest_mitogenic_signaling)
print(Heatmap_cell_cycle)






# Filter the samples with a significant down-regulation of additional proliferation gene sets
# Load the complete table containing the p values of additional proliferation gene sets of the seleted samples

# Prolif_p_val_gene_sets_seleted_samples <- GREIN_datasets_complete_samples %>%
#   dplyr::mutate(Samples = recode(Samples, "GSE135514_MDAMB468_Das_24h" = "GSE135514_24h", 
#                           "GSE135514_MDAMB468_Das_72h" = "GSE135514_72h",
#                           "GSE103021_LU65_ARS_1620" = "GSE103021",
#                           "GSE143462_KYSE150_omipalisib" = "GSE143462",
#                           "GSE177054_MDAMB231_Ribociclib" = "GSE177054",
#                           "GSE121105_BT474_trastuzumab" = "GSE121105")) %>%
#   dplyr::filter(Samples %in% c("GSE103021", "GSE121105", "GSE143462",
#                                "GSE177054", "GSE135514_24h", "GSE135514_72h")) %>%
#   dplyr::filter(Wilcox_test_Proliferation_Science_adjp < 0.05 & Wilcox_test_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019_adjp < 0.05 & Wilcox_test_Meta_PCNA_adjp < 0.05)

# The samples have a significant down-regulation of the additional proliferation gene sets.










# Generate a plot showing the NES values of all the signaling pathways for the selected samples
# Prepare the data frame for the plot by selecting all the signaling pathways of interest
Plot_GREIN_selected_samples <- dplyr::select(Cell_cycle_arrest_mitogenic_signaling, ID, Description, NES, p.adjust, Sample)

all_ids <- c("hsa04110", "hsa03030", "hsa04010", "hsa03082", "hsa04625", 
             "hsa03460", "hsa03440", "hsa03420", "hsa03430", "hsa03410", 
             "hsa04115", "hsa04014", "hsa03320", "hsa04150", "hsa04390", 
             "hsa04310", "hsa04657", "hsa04015", "hsa04668", "hsa04371", 
             "hsa04151", "hsa04910", "hsa04330", "hsa04350",
             "hsa04630", "hsa04012", "hsa04152")

# Create a complete data frame with all the pathways and samples
df_complete <- Plot_GREIN_selected_samples %>%
  dplyr::select(Sample, ID) %>%
  dplyr::distinct() %>%
  complete(Sample, ID = all_ids) %>%
  dplyr::left_join(Plot_GREIN_selected_samples, by = c("Sample", "ID")) %>%
  replace_na(list(Description = "", NES = 0, p.adjust = 0))


# Fill in the missing description names
Description_names <- dplyr::select(combined_kk2_df_filtered, ID, Description) %>%
  dplyr::distinct()

Plot_GREIN_selected_samples_complete <- dplyr::left_join(df_complete, Description_names, by = "ID") %>%
  dplyr::select(-Description.x) %>%
  dplyr::rename(Description = Description.y) %>%
  dplyr::mutate(Description = gsub("signaling ", "", Description)) %>% # Remove "signaling pathway" from the Description column
  filter(!Sample %in% c("X778_PALBO_72h", "X778_PALBO_EP", "X778_nutlin_72h", "X778_TUNIC_72h"))



# Generate the heat map
ggplot(Plot_GREIN_selected_samples_complete, aes(x = Sample, y = Description, fill = NES)) +
  geom_tile(color = "black") +  # Add black lines around the tiles
  scale_fill_gradient(low = "blue", high = "white") +  # Adjust colors as needed
  #geom_text(aes(label = ifelse(NES < 0, sprintf("%.2g", NES), "")), color = "black") +  # Display NES values only if < 0
  theme_minimal() +
  labs(title = "",
       x = "",
       y = "",
       fill = "NES") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12),  # Increase size of x-axis text
        axis.text.y = element_text(size = 16),  # Increase size of y-axis text
        axis.title.x = element_text(size = 16),  # Increase size of x-axis title
        axis.title.y = element_text(size = 14)) 

















# GSEA with specific proliferation gene sets ####
# I use clusterProfiler to determine the enrichment of specific proliferation signatures
# Load the proliferation signatures


# Creating the TERM2GENE data frame
# Define the gene sets
proliferation_science_genes <- c(
  "TPX2", "PRC1", "FOXM1", "CDK1", "C20orf24", "MCM2", "H2AFZ", "TOP2A", "PCNA", "UBE2C",
  "MELK", "TRIP13", "NCAPD2", "MCM7", "RNASEH2A", "RAD51AP1", "KIF20A", "CDC45", "MAD2L1", "ESPL1",
  "CCNB2", "FEN1", "TTK", "CCT5", "RFC4", "ATAD2", "CKAP5", "NUP205", "CDC20", "CKS2",
  "RRM2", "ELAVL1", "CCNB1", "RRM1", "AURKB", "MSH6", "EZH2", "CTPS", "DKC1", "OIP5",
  "CDCA8", "PTTG1", "CEP55", "H2AFX", "CMAS", "NCAPH", "MCM10", "LSM4", "NCAPG2", "ASF1B",
  "ZWINT", "PBK", "ZWILCH", "CDCA3", "ECT2", "CDC6", "UNG", "MTCH2", "RAD21", "ACTL6A",
  "GPI", "SRSF2", "HDGF", "NXT1", "NEK2", "DHCR7", "AURKA", "NDUFAB1", "TMEM194A", "KIF4A"
)

meta_pcna_genes <- c(
  "PCNA", "CKS2", "NUSAP1", "RRM2", "ZWINT", "PRC1", "TFDP1", "CCNA2", "CCNB1", "MELK",
  "TPX2", "BIRC5", "NCAPG2", "RFWD3", "TACC3", "CDK1", "PCLAF", "MCM2", "MCM7", "MCM5",
  "ASF1B", "CCNB2", "GTPBP2", "KIF20A", "PTTG1", "AURKA", "CDC20", "FOXM1", "GINS2", "MAD2L1",
  "UBE2C", "GINS1", "MCM6", "NCAPD2", "NUP37", "CKS1B", "MCM7", "ARID3A", "AURKB", "CDCA8",
  "SNRPD1", "SNRPD2", "TROAP", "MIS18A", "DDX39A", "AHSP", "ESPL1", "HMBS", "LSM6", "MCM4",
  "NCAPD3", "PPIH", "CDC45", "DTL", "EPB42", "HMGN2", "MCM3", "PF4", "PPBP", "RFC4",
  "TOP2A", "BUB1B", "CDT1", "FEN1", "GATA1", "GYPB", "HMGB2", "KIF22", "KLF1", "CENPU",
  "RFC3", "RHCE", "TAL1", "TCF3", "TCF7L1", "ALAS2", "TSPO2", "CDCA3", "CENPA", "CKLF",
  "GTSE1", "NFE2", "OIP5", "SHCBP1", "SNF8", "SNRPB", "SPTA1", "KEL", "KIF2C", "LYL1",
  "PSMD9", "RACGAP1", "RPIA", "TIMELESS", "TRMT5", "TYMS", "VRK1", "FBXO7", "H3-3A", "NUDT1",
  "PLEK", "POLE2", "RHAG", "NSD2", "APOBEC3B", "BPGM", "CHAF1A", "DNAJC9", "FBXO5", "KIF4A",
  "KLF15", "MKI67", "PGD", "RPA3", "SRSF2", "TRIM58", "ADAMTS13", "CDCA4", "CDKN3", "FECH",
  "LBR", "LIG1", "LMNB1", "MICB", "NUP210", "ORC6", "RAD51AP1", "RHD", "SMC4", "GYPA",
  "RPP30", "TRIM10"
)

essential_cell_cycle_genes <- c(
  "ANAPC2", "ATP6V1B2", "ATR", "AURKA", "BIRC5", "BUB3", "CCDC101", "CDC16", "CDC23", "CDC45",
  "CDC5L", "CDC6", "CDK1", "CDT1", "CHAF1A", "CHEK1", "DAD1", "DCTN2", "DDB1", "DDX54",
  "DGCR14", "DHX30", "DONSON", "DTL", "EEF2", "EIF3C", "EP300", "ESPL1", "FAU", "GINS2",
  "HINFP", "INTS3", "KIF11", "KIF23", "NDC80", "NELFB", "NOL9", "NOP9", "NSF", "ORAOV1",
  "ORC3", "PCID2", "PLK1", "POLA1", "POLE2", "PPP1CB", "PPP2CA", "PRMT1", "PRPF19", "PSMA3",
  "PSMB7", "RACGAP1", "RFC2", "RPA1", "RPS13", "RPS2", "RPS6", "RRM1", "RRM2", "SEC61A1",
  "SLC7A6OS", "SMARCA4", "SRSF3", "TARS", "TICRR", "TIMELESS", "TONSL", "TOP2A", "TPX2", "U2AF2",
  "UPF1", "WDR74", "WEE1", "YME1L1", "ZNF830"
)

# Create the TERM2GENE data frame
TERM2GENE <- data.frame(
  term = c(rep("Proliferation_Science_paper", length(proliferation_science_genes)),
           rep("Meta_PCNA", length(meta_pcna_genes)),
           rep("Essential_cell_cycle_Viner_Breuer_Cell_Division_2019", length(essential_cell_cycle_genes))),
  Gene_symbol = c(proliferation_science_genes, meta_pcna_genes, essential_cell_cycle_genes)
)


# Rename the values in the 'term' column
TERM2GENE <- TERM2GENE %>%
  mutate(term = recode(term,
                       "Proliferation_Science_paper" = "Signt_1",
                       "Meta_PCNA" = "Signt_2",
                       "Essential_cell_cycle_Viner_Breuer_Cell_Division_2019" = "Signt_3"))
# View the TERM2GENE data frame
print(TERM2GENE)



library(biomaRt)
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host="https://useast.ensembl.org")
gene_ids <- (TERM2GENE$Gene_symbol) # This reads the Ensembl transcripts IDs in the specific column of the 'Table_ANNOTATED file.

res_gene_entrez_id <- getBM(attributes = c('entrezgene_id',
                                           'external_gene_name'),
                            filters = 'external_gene_name',
                            values = gene_ids,
                            mart = mart)
res_gene_entrez_id <- dplyr::filter(res_gene_entrez_id, !is.na(entrezgene_id)) %>%
  dplyr::rename(Gene_symbol = external_gene_name) %>%
  distinct(Gene_symbol, .keep_all = TRUE)

TERM2GENE_entrez <- dplyr::left_join(TERM2GENE, res_gene_entrez_id, by = "Gene_symbol") %>%
  dplyr::filter(!is.na(entrezgene_id)) %>%
  dplyr::select(-Gene_symbol)
  TERM2GENE <- TERM2GENE_entrez
  print(TERM2GENE)
  



  library(clusterProfiler)
  library(dplyr)
  
  process_and_get_gsea_results <- function(data_frame, log_fc_cols, TERM2GENE) {
    gsea_results_list <- list()  # Initialize a list to store the GSEA results
    
    for (i in seq_along(log_fc_cols)) {
      log_fc_col <- log_fc_cols[i]
      plot_title <- plot_titles[i]
      
      # Process data
      all_values_entrez <- data_frame %>%
        dplyr::select(entrezgene_id, !!sym(log_fc_col)) %>%
        dplyr::filter(!is.na(!!sym(log_fc_col))) %>%
        dplyr::filter(!is.na(entrezgene_id))
      
      all_values_Log2_FC_entrez <- all_values_entrez[[log_fc_col]]
      names(all_values_Log2_FC_entrez) <- all_values_entrez$entrezgene_id
      
      # Remove duplicates
      all_values_Log2_FC_entrez <- all_values_Log2_FC_entrez[!duplicated(names(all_values_Log2_FC_entrez))]
      
      all_values_Log2_FC_entrez <- na.omit(all_values_Log2_FC_entrez)
      all_values_Log2_FC_entrez <- sort(all_values_Log2_FC_entrez, decreasing = TRUE)
      
      # Debugging: Print the number of genes before and after filtering
      cat("Sample:", plot_titles[i], "\n")
      cat("Number of genes:", length(all_values_Log2_FC_entrez), "\n")
      
      if (length(all_values_Log2_FC_entrez) == 0) {
        cat("No genes available for", plot_title, "\n")
        next
      }
      
      # Run GSEA
      gsea_res <- GSEA(geneList = all_values_Log2_FC_entrez,
                       TERM2GENE = TERM2GENE,
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "BH",
                       minGSSize = 10,
                       maxGSSize = 1000)
      
      if (is.null(gsea_res) || nrow(gsea_res) == 0) {
        cat("No significant pathways found for", plot_title, "\n")
        next
      }
      
      # Add GSEA results to the list
      gsea_results_list[[plot_title]] <- gsea_res
    }
    
    # Combine all GSEA results into a single data frame
    combined_gsea <- do.call(rbind, lapply(names(gsea_results_list), function(name) {
      df <- as.data.frame(gsea_results_list[[name]])
      df$Sample <- name
      return(df)
    }))
    
    return(combined_gsea)
  }
  
  # Define your data frame, log FC columns, and plot titles
  data_frame <- Expression_changes_targeted_agents_genes
  log_fc_cols <- c("GSE93124_SW480_crizotinib_Log_FC", "GSE178755_HCC827_erlotinib_Log_FC", 
                   "GSE137535_SH454_Log_FC", "GSE105083_PANC1_CM03_24h_Log_FC", "GSE118548_T84_Trametinib_Log_FC", "GSE148096_TNBC_erlotinib_Log_FC", 
                   "GSE143462_KYSE150_omipalisib_Log_FC", "GSE121105_BT474_trastuzumab_Log_FC",
                   "GSE103021_LU65_ARS_1620_Log_FC", "GSE135514_MDAMB468_Das_24h_Log_FC", "GSE135514_MDAMB468_Das_72h_Log_FC",
                   "GSE105083_MIA_PACA2_CM03_6h_Log_FC", "GSE105083_MIA_PACA2_CM03_24h_Log_FC",
                   "X778_nutlin_Log_FC", "X778_TUNIC_Log_FC", "X778_PALBO_Log_FC", "SKBR3_Log_FC", "SKMEL28_Log_FC", 
                   "X778_PALBO_EP_Log_FC", "SKBR3_EP_Log_FC", "SKMEL28_EP_Log_FC")
  
  plot_titles <- c("GSE93124", "GSE178755", 
                   "GSE137535_SH454", "GSE105083_PANC1_CM03_24h", "GSE118548_tra", "GSE148096", 
                   "GSE143462", "GSE121105",
                   "GSE103021", "GSE135514_24h", "GSE135514_72h",
                   "GSE105083_MIA_PACA2_6h", "GSE105083_MIA_PACA2_24h",
                   "X778_nutlin_72h", "X778_TUNIC_72h", "X778_PALBO_72h", "SKBR3_72h", "SKMEL28_72h", 
                   "X778_PALBO_EP", "SKBR3_EP", "SKMEL28_EP")
  
  # Call the function to get the results
  combined_gsea <- process_and_get_gsea_results(data_frame, log_fc_cols, TERM2GENE)

  # # Replace the names of the proliferation signatures with shorter ones
  #   combined_gsea <- combined_gsea %>%
  #   mutate(Description = case_when(
  #     Description == "Proliferation_Science_paper" ~ "Signt_1",
  #     Description == "Meta_PCNA" ~ "Signt_2",
  #     Description == "Essential_cell_cycle_Viner_Breuer_Cell_Division_2019" ~ "Signt_3",
  #     TRUE ~ Description  # Retain other values as they are
  #   ))
  # 



# Add these data to those obtained using the KEGG pathways
# Filter the combined_gsea object for the data and the samples of interest
combined_gsea_filtered <- dplyr::select(combined_gsea, Sample, ID, NES, p.adjust, Description)
# Delete rownames
rownames(combined_gsea_filtered) <- NULL
# Select the samples of interest
combined_gsea_seleted <- dplyr::filter(combined_gsea_filtered, Sample %in% Plot_GREIN_selected_samples_complete$Sample)
  
Plot_GREIN_selected_samples_complete_plot <- dplyr::bind_rows(combined_gsea_seleted, Plot_GREIN_selected_samples_complete)



# Plot the data
desired_order <- c("Signt_1", "Signt_2", "Signt_3", "Cell cycle", "DNA replication",
                   "Homologous recombination", "Mismatch repair", "Fanconi anemia pathway", "Base excision repair", "Nucleotide excision repair", "ATP-dependent chromatin remodeling", "p53 pathway",
                   "PPAR pathway", "Apelin pathway", "C-type lectin receptor pathway", "ErbB pathway", "Hippo pathway",
                   "IL-17 pathway", "Insulin pathway", "JAK-STAT pathway", "MAPK pathway",  "AMPK pathway", "mTOR pathway", "PI3K-Akt pathway", "Ras pathway",
                   "Rap1 pathway", "TGF-beta pathway", "TNF pathway", "Wnt pathway", "Hedgehog pathway", "Phospholipase D pathway", "Notch pathway")


# Reorder the 'Description' column based on the desired order
Plot_GREIN_selected_samples_complete_plot <- Plot_GREIN_selected_samples_complete_plot %>%
  mutate(Description = factor(Description, levels = desired_order),
         neg_log10_p_adjust = -log10(p.adjust))

# Reorder the 'Sample' column based on the desired order
Plot_GREIN_selected_samples_complete_plot <- Plot_GREIN_selected_samples_complete_plot %>%
  mutate(Sample = factor(Sample, levels = c("GSE103021", "GSE121105", "GSE143462", "GSE135514_24h", "GSE135514_72h")))

# Generate the plot with circles
Plot_cell_proliferation <- ggplot(Plot_GREIN_selected_samples_complete_plot, aes(x = Sample, y = Description, fill = NES, size = neg_log10_p_adjust)) +
  geom_point(shape = 21, color = "black", stroke = 0.5) +  # Use shape = 21 for circles with border color
  scale_fill_gradient2(low = "dodgerblue4", mid = "white", high = "firebrick3", midpoint = 0, limits = c(-3, 3)) +  # Set colors and midpoint for the fill gradient
  scale_size_continuous(range = c(2, 10)) +  # Adjust the range of circle sizes as needed
  theme_bw() +
  labs(title = "",
       x = "",
       y = "",
       fill = "NES",
       size = "-log10(p.adjust)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16),  # Increase size of x-axis text
        axis.text.y = element_text(size = 14),  # Increase size of y-axis text
        axis.title.x = element_text(size = 16),  # Increase size of x-axis title
        axis.title.y = element_text(size = 10))  # Increase size of y-axis title

print(Plot_cell_proliferation)




























# Fig.2C violin plots MFE vs Log2_FC ####
# Load required packages
library(gridExtra)
library(grid)

# Load the original violin plots without the jitter function:
Plot_MFE_GSE143462 <- ggplot(MFE_res_GSE143462, aes(x = factor(log2FoldChange > 0), y = Normalized_MFE)) +
  geom_violin(fill = "white", size = 1.5) +
  #geom_jitter(aes(color = above_quantile_075), width = 0.1, size = 1, alpha = 1) +
  #scale_color_manual(values = c("FALSE" = "black", "TRUE" = "dodgerblue3")) +
  theme_classic() +
  ylim(-0.6, 0) +
  theme(
    legend.position = "none", 
    axis.title = element_text(size = 20), 
    plot.title = element_text(size = 10),
    axis.text.x = element_text(size = 16, color = "black")
  ) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  labs(
    y = "Normalized_MFE", 
    x = NULL, 
    title = "GSE143462_Log2_FC_adj.P<0.05"
  ) +
  scale_x_discrete(labels = c("log2_FC < 0", "log2_FC > 0"))




Plot_MFE_GSE121105 <- ggplot(MFE_res_GSE121105_trastuzumab, aes(x = factor(log2FoldChange > 0), y = Normalized_MFE)) +
  geom_violin(fill = "white", size = 1.5) +
  #geom_jitter(aes(color = above_quantile_075), width = 0.1, size = 1, alpha = 1) +
  #scale_color_manual(values = c("FALSE" = "black", "TRUE" = "dodgerblue3")) +
  theme_classic() +
  ylim(-0.6, 0) +
  theme(
    legend.position = "none", 
    axis.title = element_text(size = 20), 
    plot.title = element_text(size = 10),
    axis.text.x = element_text(size = 16, color = "black")
  ) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  labs(
    y = "Normalized_MFE", 
    x = NULL, 
    title = "GSE121105_trastuzumab_Log2_FC_adj.P<0.05"
  ) +
  scale_x_discrete(labels = c("log2_FC < 0", "log2_FC > 0"))





Plot_MFE_GSE103021_ARS_1620_LU65 <- ggplot(MFE_res_GSE103021_ARS_1620_LU65, aes(x = factor(log2FoldChange > 0), y = Normalized_MFE)) +
  geom_violin(fill = "white", size = 1.5) +
  #geom_jitter(aes(color = above_quantile_075), width = 0.1, size = 1, alpha = 1) +
  #scale_color_manual(values = c("FALSE" = "black", "TRUE" = "dodgerblue3")) +
  theme_classic() +
  ylim(-0.6, 0) +
  theme(
    legend.position = "none", 
    axis.title = element_text(size = 20), 
    plot.title = element_text(size = 10),
    axis.text.x = element_text(size = 16, color = "black")
  ) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  labs(
    y = "Normalized_MFE", 
    x = NULL, 
    title = "GSE103021_ARS_1620_LU65_Log2_FC_adj.P<0.05"
  ) +
  scale_x_discrete(labels = c("log2_FC < 0", "log2_FC > 0"))







Plot_MFE_GSE135514_24h <- ggplot(MFE_res_GSE135514_Das_24h, aes(x = factor(log2FoldChange > 0), y = Normalized_MFE)) +
  geom_violin(fill = "white", size = 1.5) +
  #geom_jitter(aes(color = above_quantile_075), width = 0.1, size = 1, alpha = 1) +
  #scale_color_manual(values = c("FALSE" = "black", "TRUE" = "dodgerblue3")) +
  theme_classic() +
  ylim(-0.6, 0) +
  theme(
    legend.position = "none", 
    axis.title = element_text(size = 20), 
    plot.title = element_text(size = 10),
    axis.text.x = element_text(size = 16, color = "black")
  ) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  labs(
    y = "Normalized_MFE", 
    x = NULL, 
    title = "GSE135514_Das_24h_Log2_FC_adj.P<0.05"
  ) +
  scale_x_discrete(labels = c("log2_FC < 0", "log2_FC > 0"))



Plot_MFE_GSE135514_72h <- ggplot(MFE_res_GSE135514_Das_72h, aes(x = factor(log2FoldChange > 0), y = Normalized_MFE)) +
  geom_violin(fill = "white", size = 1.5) +
  #geom_jitter(aes(color = above_quantile_075), width = 0.1, size = 1, alpha = 1) +
  #scale_color_manual(values = c("FALSE" = "black", "TRUE" = "dodgerblue3")) +
  theme_classic() +
  ylim(-0.6, 0) +
  theme(
    legend.position = "none", 
    axis.title = element_text(size = 20), 
    plot.title = element_text(size = 10),
    axis.text.x = element_text(size = 16, color = "black")
  ) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  labs(
    y = "Normalized_MFE", 
    x = NULL, 
    title = "GSE135514_Das_72h_Log2_FC_adj.P<0.05"
  ) +
  scale_x_discrete(labels = c("log2_FC < 0", "log2_FC > 0"))






# Combine all plots
plots <- list(
  Plot_MFE_GSE143462,
  Plot_MFE_GSE121105,
  Plot_MFE_GSE103021_ARS_1620_LU65,
  Plot_MFE_GSE135514_24h,
  Plot_MFE_GSE135514_72h
)

# Add margin to each plot
plots_with_margin <- lapply(plots, function(plot) {
  plot + 
    geom_violin(fill = "lightgrey", color = "black", size = 0.1) +
    #geom_jitter(position = position_jitter(width = 0.1), size = 0.05, color = "black") +
    geom_boxplot(width = 0.1, color = "black", outlier.size = 0) +
    #scale_color_manual(values = c("FALSE" = "black", "TRUE" = "black")) +  # Update colors
    theme(
      plot.margin = unit(c(4, 0.2, 0.2, 0.2), "cm"),
      axis.title = element_blank(),  # Remove y-axis text)  # Adjust the margin values as needed
      axis.text.x = element_text(size = 14), # Set x-axis text size to 12
      plot.title = element_blank() # Remove titles
    )
})

# Arrange plots in 5 columns and 2 rows with spacing
combined_plots <- grid.arrange(
  grobs = plots_with_margin,
  ncol = 5,
  nrow = 1
)



# These plots have been excluded because will be analysed in separately
# Plot_MFE_GSE105083_MIA_PACA2_6h
# Plot_MFE_GSE105083_MIA_PACA2_24h


# Mean MFE values and KS test p values on the violin plots
Mean_MFE_res_GSE143462_UP <- mean(MFE_res_GSE143462_UP$Normalized_MFE)
Mean_MFE_res_GSE143462_DOWN <- mean(MFE_res_GSE143462_DOWN$Normalized_MFE)
ks_test_res_GSE143462_all_UP_vs_DOWN$p.value

Mean_MFE_res_GSE121105_trastuzumab_UP <- mean(MFE_res_GSE121105_trastuzumab_UP$Normalized_MFE)
Mean_MFE_res_GSE121105_trastuzumab_DOWN <- mean(MFE_res_GSE121105_trastuzumab_DOWN$Normalized_MFE)
ks_test_res_GSE121105_trastuzumab_all_UP_vs_DOWN$p.value

Mean_MFE_res_GSE103021_ARS_1620_LU65_UP <- mean(MFE_res_GSE103021_ARS_1620_LU65_UP$Normalized_MFE)
Mean_MFE_res_GSE103021_ARS_1620_LU65_DOWN <- mean(MFE_res_GSE103021_ARS_1620_LU65_DOWN$Normalized_MFE)
ks_test_res_GSE103021_ARS_1620_LU65_all_UP_vs_DOWN$p.value

Mean_MFE_res_GSE135514_Das_24h_UP <- mean(MFE_res_GSE135514_Das_24h_UP$Normalized_MFE)
Mean_MFE_res_GSE135514_Das_24h_DOWN <- mean(MFE_res_GSE135514_Das_24h_DOWN$Normalized_MFE)
ks_test_res_GSE135514_Das_24h_all_UP_vs_DOWN$p.value

Mean_MFE_res_GSE135514_Das_24h_UP <- mean(MFE_res_GSE135514_Das_24h_UP$Normalized_MFE)
Mean_MFE_res_GSE135514_Das_24h_DOWN <- mean(MFE_res_GSE135514_Das_24h_DOWN$Normalized_MFE)
ks_test_res_GSE135514_Das_24h_all_UP_vs_DOWN$p.value

Mean_MFE_res_GSE135514_Das_72h_UP <- mean(MFE_res_GSE135514_Das_72h_UP$Normalized_MFE)
Mean_MFE_res_GSE135514_Das_72h_DOWN <- mean(MFE_res_GSE135514_Das_72h_DOWN$Normalized_MFE)
ks_test_res_GSE135514_Das_72h_all_UP_vs_DOWN$p.value



# Generate the data frame
results <- data.frame(
  Dataset = c(
    "GSE143462",
    "GSE121105_trastuzumab",
    "GSE103021_ARS_1620_LU65",
    "GSE135514_Das_24h",
    "GSE135514_Das_72h"
  ),
  Mean_MFE_UP = c(
    mean(MFE_res_GSE143462_UP$Normalized_MFE),
    mean(MFE_res_GSE121105_trastuzumab_UP$Normalized_MFE),
    mean(MFE_res_GSE103021_ARS_1620_LU65_UP$Normalized_MFE),
    mean(MFE_res_GSE135514_Das_24h_UP$Normalized_MFE),
    mean(MFE_res_GSE135514_Das_72h_UP$Normalized_MFE)
  ),
  Mean_MFE_DOWN = c(
    mean(MFE_res_GSE143462_DOWN$Normalized_MFE),
    mean(MFE_res_GSE121105_trastuzumab_DOWN$Normalized_MFE),
    mean(MFE_res_GSE103021_ARS_1620_LU65_DOWN$Normalized_MFE),
    mean(MFE_res_GSE135514_Das_24h_DOWN$Normalized_MFE),
    mean(MFE_res_GSE135514_Das_72h_DOWN$Normalized_MFE)
  ),
  KS_Test_p_value = c(
    ks_test_res_GSE143462_all_UP_vs_DOWN$p.value,
    ks_test_res_GSE121105_trastuzumab_all_UP_vs_DOWN$p.value,
    ks_test_res_GSE103021_ARS_1620_LU65_all_UP_vs_DOWN$p.value,
    ks_test_res_GSE135514_Das_24h_all_UP_vs_DOWN$p.value,
    ks_test_res_GSE135514_Das_72h_all_UP_vs_DOWN$p.value
  )
)

# Display the data frame
print(results)






# Create a data frame with only the Wilcoxon rank-sum (U test) results
wilcoxon_results <- data.frame(
  Dataset = c(
    "GSE143462",
    "GSE121105_trastuzumab",
    "GSE103021_ARS_1620_LU65",
    "GSE135514_Das_24h",
    "GSE135514_Das_72h"
  ),
  Wilcoxon_p_value = c(
    wilcox.test(MFE_res_GSE143462_DOWN$Normalized_MFE,
                MFE_res_GSE143462_UP$Normalized_MFE,
                paired = FALSE)$p.value,
    
    wilcox.test(MFE_res_GSE121105_trastuzumab_DOWN$Normalized_MFE,
                MFE_res_GSE121105_trastuzumab_UP$Normalized_MFE,
                paired = FALSE)$p.value,
    
    wilcox.test(MFE_res_GSE103021_ARS_1620_LU65_DOWN$Normalized_MFE,
                MFE_res_GSE103021_ARS_1620_LU65_UP$Normalized_MFE,
                paired = FALSE)$p.value,
    
    wilcox.test(MFE_res_GSE135514_Das_24h_DOWN$Normalized_MFE,
                MFE_res_GSE135514_Das_24h_UP$Normalized_MFE,
                paired = FALSE)$p.value,
    
    wilcox.test(MFE_res_GSE135514_Das_72h_DOWN$Normalized_MFE,
                MFE_res_GSE135514_Das_72h_UP$Normalized_MFE,
                paired = FALSE)$p.value
  )
)

# Display the data frame
print(wilcoxon_results)





















# Fig.2D high- and low-MFE transcripts as fraction of total counts ####

# Generate a data frame with the following values:
# High MFE transcripts in control samples as fraction of total counts
# High MFE transcripts in drug-treated samples as fraction of total counts
# Fisher's exact test to compare the proportions of high MFE transcripts counts between the two groups

Data_high_MFE_fraction_total_counts <- data.frame(
Samples = c("GSE143462",
            "GSE121105", 
            "GSE103021",
            "GSE135514_24h", 
            "GSE135514_72h"),
High_MFE_of_total_control = c(GSE143462_control_high_MFE_TX_counts_fraction,
                              GSE121105_trastuzumab_control_high_MFE_TX_counts_fraction,
                              GSE103021_ARS_1620_LU65_control_high_MFE_TX_counts_fraction,
                              GSE135514_Das_24h_control_high_MFE_TX_counts_fraction,
                              GSE135514_Das_72h_control_high_MFE_TX_counts_fraction),
High_MFE_of_total_drug_treated = c(GSE143462_Omipalisib_high_MFE_TX_counts_fraction,
                                   GSE121105_trastuzumab_trastuzumab_high_MFE_TX_counts_fraction,
                                   GSE103021_ARS_1620_LU65_ARS_1620_LU65_high_MFE_TX_counts_fraction,
                                   GSE135514_Das_24h_Das_24h_high_MFE_TX_counts_fraction,
                                   GSE135514_Das_72h_Das_72h_high_MFE_TX_counts_fraction),
Z_test_result_75th = c(GSE143462_z_test_counts_75_p_value,
                 GSE121105_trastuzumab_z_test_counts_75_p_value,
                 GSE103021_ARS_1620_LU65_z_test_counts_75_p_value,
                 GSE135514_Das_24h_z_test_counts_75_p_value,
                 GSE135514_Das_72h_z_test_counts_75_p_value))

print(Data_high_MFE_fraction_total_counts)


# Generate a bar plot with samples in the x-axis and the fraction of high MFE transcripts as y-axis

# Reshape the data frame for ggplot2
library(tidyr)
df_long <- pivot_longer(Data_high_MFE_fraction_total_counts,
                        cols = c(High_MFE_of_total_control, High_MFE_of_total_drug_treated),
                        names_to = "Condition",
                        values_to = "Value")


# Generate the bar plot


# Load necessary libraries
library(ggplot2)
library(tidyr)

# Assuming your data frame is named Data_high_MFE_fraction_total_counts
# Here's a sample structure of your data frame:
# Data_high_MFE_fraction_total_counts <- data.frame(
#     Samples = c("GSE105083", "GSE118548_tra", "GSE148096", "GSE143462", "GSE118548_JQ1", "GSE121105"),
#     High_MFE_of_total_control = c(0.1, 0.2, 0.15, 0.2, 0.25, 0.3),
#     High_MFE_of_total_drug_treated = c(0.05, 0.25, 0.2, 0.3, 0.35, 0.4),
#     Fisher_result = c("ns", "significant", "ns", "significant", "ns", "significant")
# )

# Define the desired order of samples
sample_order <- c("GSE103021", "GSE121105", "GSE143462", "GSE135514_24h", "GSE135514_72h")

# Convert Samples column to a factor with the specified order
Data_high_MFE_fraction_total_counts$Samples <- factor(Data_high_MFE_fraction_total_counts$Samples, levels = sample_order)

# Reshape the data frame for ggplot2
df_long <- pivot_longer(Data_high_MFE_fraction_total_counts,
                        cols = c(High_MFE_of_total_control, High_MFE_of_total_drug_treated),
                        names_to = "Condition",
                        values_to = "Value")

# Plot the data
Plot_high_fraction_tot_counts <- ggplot(df_long, aes(x = Samples, y = Value, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = round(Value, 3)), position = position_dodge(width = 0.8), hjust = -0.3, vjust = 0.4, angle = 90, size = 4) +
  labs(title = "",
       x = "",
       y = "Fraction of total counts",
       fill = "") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.title.y = element_text(size = 14),
        legend.position = "none") +
  ylim(0, 0.4) +
  scale_fill_manual(values = c("High_MFE_of_total_control" = "#8cc5e3",
                               "High_MFE_of_total_drug_treated" = "#2066a8"),
                    labels = c("High_MFE_of_total_control" = "Untreated",
                               "High_MFE_of_total_drug_treated" = "Drug_treated"))

print(Plot_high_fraction_tot_counts)














Data_low_MFE_fraction_total_counts <- data.frame(
  Samples = c("GSE143462",
              "GSE121105", 
              "GSE103021",
              "GSE135514_24h", 
              "GSE135514_72h"),
  low_MFE_of_total_control = c(GSE143462_control_low_MFE_TX_counts_fraction,
                                GSE121105_trastuzumab_control_low_MFE_TX_counts_fraction,
                                GSE103021_ARS_1620_LU65_control_low_MFE_TX_counts_fraction,
                               GSE135514_Das_24h_control_low_MFE_TX_counts_fraction,
                               GSE135514_Das_72h_control_low_MFE_TX_counts_fraction),
  low_MFE_of_total_drug_treated = c(GSE143462_Omipalisib_low_MFE_TX_counts_fraction,
                                    GSE121105_trastuzumab_trastuzumab_low_MFE_TX_counts_fraction,
                                     GSE103021_ARS_1620_LU65_low_MFE_TX_counts_fraction,
                                    GSE135514_Das_24h_low_MFE_TX_counts_fraction,
                                    GSE135514_Das_72h_low_MFE_TX_counts_fraction),
  Z_test_result_25th = c(GSE143462_z_test_counts_25_p_value,
                    GSE121105_trastuzumab_z_test_counts_25_p_value,
                    GSE103021_ARS_1620_LU65_z_test_counts_25_p_value,
                    GSE135514_Das_24h_z_test_counts_25_p_value,
                    GSE135514_Das_72h_z_test_counts_25_p_value))

print(Data_low_MFE_fraction_total_counts)



# Generate a bar plot with samples in the x-axis and the fraction of high MFE transcripts as y-axis

# Reshape the data frame for ggplot2
library(tidyr)
df_long <- pivot_longer(Data_low_MFE_fraction_total_counts,
                        cols = c(low_MFE_of_total_control, low_MFE_of_total_drug_treated),
                        names_to = "Condition",
                        values_to = "Value")


# Generate the bar plot


# Load necessary libraries
library(ggplot2)
library(tidyr)

# Assuming your data frame is named Data_low_MFE_fraction_total_counts
# Here's a sample structure of your data frame:
# Data_low_MFE_fraction_total_counts <- data.frame(
#     Samples = c("GSE105083", "GSE118548_tra", "GSE148096", "GSE143462", "GSE118548_JQ1", "GSE121105"),
#     low_MFE_of_total_drug_treated = c(0.1, 0.2, 0.15, 0.2, 0.25, 0.3),
#     low_MFE_of_total_drug_treated = c(0.05, 0.25, 0.2, 0.3, 0.35, 0.4),
#     Fisher_result = c("ns", "significant", "ns", "significant", "ns", "significant")
# )

# Define the desired order of samples
sample_order <- c("GSE103021", "GSE121105", "GSE143462", "GSE135514_24h", "GSE135514_72h")

# Convert Samples column to a factor with the specified order
Data_low_MFE_fraction_total_counts$Samples <- factor(Data_low_MFE_fraction_total_counts$Samples, levels = sample_order)

# Reshape the data frame for ggplot2
df_long <- pivot_longer(Data_low_MFE_fraction_total_counts,
                        cols = c(low_MFE_of_total_control, low_MFE_of_total_drug_treated),
                        names_to = "Condition",
                        values_to = "Value")

# Plot the data
Plot_low_fraction_tot_counts <- ggplot(df_long, aes(x = Samples, y = Value, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = round(Value, 3)), position = position_dodge(width = 0.8), hjust = -0.3, vjust = 0.4, angle = 90, size = 4) +
  labs(title = "",
       x = "",
       y = "Fraction of total counts",
       fill = "") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.title.y = element_text(size = 14),
        legend.position = "none") +
  ylim(0, 0.4) +
  scale_fill_manual(values = c("low_MFE_of_total_control" = "#8cc5e3",
                               "low_MFE_of_total_drug_treated" = "#2066a8"),
                    labels = c("low_MFE_of_total_control" = "Untreated",
                               "low_MFE_of_total_drug_treated" = "Drug_treated"))

print(Plot_low_fraction_tot_counts)




























# Fig.2E plot mitogenic signaling Science models ####

# Generate a plot showing the NES values of all the signaling pathways for the selected samples
# Prepare the data frame for the plot by selecting all the signaling pathways of interest

# IMPORTANT
# Cell_cycle_arrest_mitogenic_signaling contains aonly samamples with a down-regulation of all proliferation pathways and absence of any mitogenic signaling
# The samples that are not here have mitogenic signaling and have been excluded
# These are the samples of the Science models with mitogenic signaling:
# SKBR3_72h
# SKMEL28_72h
# SKMEL28_EP
# X778_PALBO_EP
# X778_nutlin_72h are cells treated with the MDM2-TP53 inhibitor Nutlin-3a which is not directly inhibiting mitogenic signaling but rather activating the p53 pathway
Plot_GREIN_selected_samples_Science_models <- dplyr::select (Cell_cycle_arrest_mitogenic_signaling, ID, Description, NES, p.adjust, Sample)

all_ids <- c("hsa04110", "hsa03030", "hsa04010", "hsa03082", "hsa04625", 
             "hsa03460", "hsa03440", "hsa03420", "hsa03430", "hsa03410", 
             "hsa04115", "hsa04014", "hsa03320", "hsa04150", "hsa04390", 
             "hsa04310", "hsa04657", "hsa04015", "hsa04668", "hsa04371", 
             "hsa04151", "hsa04910", "hsa04330", "hsa04350",
             "hsa04630", "hsa04012", "hsa04152")

# Create a complete data frame with all the pathways and samples
df_complete_Science_models <- Plot_GREIN_selected_samples_Science_models %>%
  dplyr::select(Sample, ID) %>%
  dplyr::distinct() %>%
  complete(Sample, ID = all_ids) %>%
  dplyr::left_join(Plot_GREIN_selected_samples_Science_models, by = c("Sample", "ID")) %>%
  replace_na(list(Description = "", NES = 0, p.adjust = 0))


# Fill in the missing description names
Description_names <- dplyr::select(combined_kk2_df_filtered, ID, Description) %>%
  dplyr::distinct()

Plot_GREIN_selected_samples_complete_Science_models <- dplyr::left_join(df_complete_Science_models, Description_names, by = "ID") %>%
  dplyr::select(-Description.x) %>%
  dplyr::rename(Description = Description.y) %>%
  dplyr::mutate(Description = gsub("signaling ", "", Description)) %>% # Remove "signaling pathway" from the Description column
  filter(Sample %in% "X778_PALBO_72h")



# Generate the heat map
ggplot(Plot_GREIN_selected_samples_complete_Science_models, aes(x = Sample, y = Description, fill = NES)) +
  geom_tile(color = "black") +  # Add black lines around the tiles
  scale_fill_gradient(low = "blue", high = "white") +  # Adjust colors as needed
  #geom_text(aes(label = ifelse(NES < 0, sprintf("%.2g", NES), "")), color = "black") +  # Display NES values only if < 0
  theme_minimal() +
  labs(title = "",
       x = "",
       y = "",
       fill = "NES") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12),  # Increase size of x-axis text
        axis.text.y = element_text(size = 16),  # Increase size of y-axis text
        axis.title.x = element_text(size = 16),  # Increase size of x-axis title
        axis.title.y = element_text(size = 14)) 




# Add these data to those obtained using the KEGG pathways
# Filter the combined_gsea object for the data and the samples of interest
combined_gsea_filtered_Science_models <- dplyr::select(combined_gsea, Sample, ID, NES, p.adjust, Description)
# Delete rownames
rownames(combined_gsea_filtered_Science_models) <- NULL
# Select the samples of interest
combined_gsea_seleted_Science_models <- dplyr::filter(combined_gsea_filtered_Science_models, Sample %in% Plot_GREIN_selected_samples_complete_Science_models$Sample)

Plot_GREIN_selected_samples_complete_plot_Science_models <- dplyr::bind_rows(combined_gsea_seleted_Science_models, Plot_GREIN_selected_samples_complete_Science_models)



# Plot the data
desired_order <- c("Signt_1", "Signt_2", "Signt_3", "Cell cycle", "DNA replication",
                   "Homologous recombination", "Mismatch repair", "Fanconi anemia pathway", "Base excision repair", "Nucleotide excision repair", "ATP-dependent chromatin remodeling", "p53 pathway",
                   "PPAR pathway", "Apelin pathway", "C-type lectin receptor pathway", "ErbB pathway", "Hippo pathway",
                   "IL-17 pathway", "Insulin pathway", "JAK-STAT pathway", "MAPK pathway",  "AMPK pathway", "mTOR pathway", "PI3K-Akt pathway", "Ras pathway",
                   "Rap1 pathway", "TGF-beta pathway", "TNF pathway", "Wnt pathway", "Hedgehog pathway", "Phospholipase D pathway", "Notch pathway")


# Reorder the 'Description' column based on the desired order
Plot_GREIN_selected_samples_complete_plot_Science_models <- Plot_GREIN_selected_samples_complete_plot_Science_models %>%
  mutate(Description = factor(Description, levels = desired_order),
         neg_log10_p_adjust = -log10(p.adjust))

# Reorder the 'Sample' column based on the desired order
Plot_GREIN_selected_samples_complete_plot_Science_models <- Plot_GREIN_selected_samples_complete_plot_Science_models %>%
  mutate(Sample = factor(Sample, levels = "X778_PALBO_72h"))

# Generate the plot with circles
Plot_cell_proliferation_Science_models <- ggplot(Plot_GREIN_selected_samples_complete_plot_Science_models, aes(x = Sample, y = Description, fill = NES, size = neg_log10_p_adjust)) +
  geom_point(shape = 21, color = "black", stroke = 0.5) +  # Use shape = 21 for circles with border color
  scale_fill_gradient2(low = "dodgerblue4", mid = "white", high = "firebrick3", midpoint = 0, limits = c(-3, 3)) +  # Set colors and midpoint for the fill gradient
  scale_size_continuous(range = c(2, 10)) +  # Adjust the range of circle sizes as needed
  theme_bw() +
  labs(title = "",
       x = "",
       y = "",
       fill = "NES",
       size = "-log10(p.adjust)") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 16),  # Increase size of x-axis text
        axis.text.y = element_text(size = 14),  # Increase size of y-axis text
        axis.title.x = element_text(size = 16),  # Increase size of x-axis title
        axis.title.y = element_text(size = 10))  # Increase size of y-axis title

print(Plot_cell_proliferation_Science_models)


























# Fig.2F violin plots Science models 72h ####

# Load the data of the Science models generated at 72h:
source("path_to/All_model_systems_72h_summary.R")

# Load the data of the Science models generated during the early phase response
source("path_to/All_model_systems_early_phase_summary.R")




# Load required packages
library(gridExtra)
library(grid)


# Plot_X778_NUTLIN_72h <- ggplot(X778_nutlin_dexseq_ANNOTATED, aes(x = factor(log2fold_treated_untreated > 0), y = Normalized_MFE)) +
#   geom_violin(fill = "lightgrey", size = 1) +
#   geom_boxplot(width = 0.1, color = "black", outlier.size = 0, size = 0.3) +
#   #geom_jitter(width = 0.1, size = 1, alpha = 0.3) +
#   theme_classic() +
#   theme(legend.position = "none", axis.title = element_text(size = 20), plot.title = element_text(size = 10),
#         axis.text.x = element_text(size = 16, color = "black")) +
#   geom_vline(xintercept = 0, size = 0.3) +
#   ylim(-0.6, 0) +
#   labs(y = "Normalized_MFE", x = NULL, title = "778_NUTLIN_72h_Log2_FC_adj.P<0.05") +
#   scale_x_discrete(labels = c("Log2_FC < 0", "Log2_FC > 0"))
# print(Plot_X778_NUTLIN_72h)
# 
# # Save the combined plot
# ggsave("Plot_X778_NUTLIN_72h.pdf", plot = Plot_X778_NUTLIN_72h, width = 8, height = 8, path = "/Users/z3514253/Desktop")





# Plot_SKBR3_LAP_72h <- ggplot(SKBR3_dexseq_ANNOTATED, aes(x = factor(log2fold_treated_untreated > 0), y = Normalized_MFE)) +
#   geom_violin(fill = "lightgrey", size = 1) +
#   geom_boxplot(width = 0.1, color = "black", outlier.size = 0, size = 0.3) +
#   #geom_jitter(width = 0.1, size = 1, alpha = 0.3) +
#   theme_classic() +
#   theme(legend.position = "none", axis.title = element_text(size = 20), plot.title = element_text(size = 10),
#         axis.text.x = element_text(size = 16, color = "black")) +
#   geom_vline(xintercept = 0, size = 0.3) +
#   ylim(-0.6, 0) +
#   labs(y = "Normalized_MFE", x = NULL, title = "SKBR3_LAP_72h_Log2_FC_adj.P<0.05") +
#   scale_x_discrete(labels = c("Log2_FC < 0", "Log2_FC > 0"))
# print(Plot_SKBR3_LAP_72h)
# 
# # Save the combined plot
# ggsave("Plot_SKBR3_LAP_72h", plot = Plot_SKBR3_LAP_72h, width = 8, height = 8, path = "/Users/z3514253/Desktop")






# Plot_SKMEL28_VEM_72h <- ggplot(SKMEL28_dexseq_ANNOTATED, aes(x = factor(log2fold_treated_untreated > 0), y = Normalized_MFE)) +
#   geom_violin(fill = "lightgrey", size = 1) +
#   geom_boxplot(width = 0.1, color = "black", outlier.size = 0, size = 0.3) +
#   #geom_jitter(width = 0.1, size = 1, alpha = 0.3) +
#   theme_classic() +
#   theme(legend.position = "none", axis.title = element_text(size = 20), plot.title = element_text(size = 10),
#         axis.text.x = element_text(size = 16, color = "black")) +
#   geom_vline(xintercept = 0, size = 0.3) +
#   ylim(-0.6, 0) +
#   labs(y = "Normalized_MFE", x = NULL, title = "SKMEL28_VEM_72h_Log2_FC_adj.P<0.05") +
#   scale_x_discrete(labels = c("Log2_FC < 0", "Log2_FC > 0"))
# print(Plot_SKMEL28_VEM_72h)
# 
# # Save the combined plot
# ggsave("Plot_SKMEL28_VEM_72h", plot = Plot_SKMEL28_VEM_72h, width = 8, height = 8, path = "/Users/z3514253/Desktop")






Plot_X778_PALBO_72h <- ggplot(X778_PALBO_dexseq_ANNOTATED, aes(x = factor(log2fold_treated_untreated > 0), y = Normalized_MFE)) +
  geom_violin(fill = "lightgrey", size = 1) +
  geom_boxplot(width = 0.1, color = "black", outlier.size = 0, size = 0.3) +
  #geom_jitter(width = 0.1, size = 1, alpha = 0.3) +
  theme_classic() +
  theme(legend.position = "none", axis.title = element_text(size = 20), plot.title = element_text(size = 10),
        axis.text.x = element_text(size = 16, color = "black")) +
  geom_vline(xintercept = 0, size = 0.3) +
  ylim(-0.6, 0) +
  labs(y = "", x = NULL, title = "") +
  scale_x_discrete(labels = c("Log2_FC < 0", "Log2_FC > 0"))
print(Plot_X778_PALBO_72h)

# Save the combined plot
ggsave("Plot_X778_PALBO72h.pdf", plot = Plot_X778_PALBO_72h, width = 6, height = 8, path = "/Users/z3514253/Desktop")



ks_test_MFE_X778_PALBO_LAP_72h_all_UP_vs_DOWN <- ks.test(X778_PALBO_dexseq_ANNOTATED_DOWN$Normalized_MFE, X778_PALBO_dexseq_ANNOTATED_UP$Normalized_MFE, alternative = "greater")
print(ks_test_MFE_X778_PALBO_LAP_72h_all_UP_vs_DOWN)




# Generate the data frame
results <- data.frame(
  Dataset = c(
    "X778_PALBO_72h"
  ),
  Mean_MFE_UP = c(
    Mean_X778_PALBO_dexseq_ANNOTATED_UP
    ),
  Mean_MFE_DOWN = c(
    Mean_X778_PALBO_dexseq_ANNOTATED_DOWN
    ),
  KS_Test_p_value = c(
    ks_test_MFE_X778_PALBO_LAP_72h_all_UP_vs_DOWN$p.value
    )
)

# Display the data frame
print(results)




# Wilcoxon rank-sum test (Mann–Whitney U) one-sided:
# Tests if UP group has greater (less negative) values than DOWN group
wilcox_test_X778_PALBO_72h <- wilcox.test(
  X778_PALBO_dexseq_ANNOTATED_UP$Normalized_MFE,
  X778_PALBO_dexseq_ANNOTATED_DOWN$Normalized_MFE,
  paired = FALSE,
  alternative = "greater"
)

# Display test result
print(wilcox_test_X778_PALBO_72h)


wilcoxon_results <- data.frame(
  Dataset = "X778_PALBO_72h",
  Wilcoxon_p_value = wilcox_test_X778_PALBO_72h$p.value
)

print(wilcoxon_results)


















# Fig.2G high- and low-MFE transcripts as fraction of total counts ####

# Create data frame with control and drug-treated fractions and Fisher’s test result
Data_high_MFE_fraction_total_counts <- data.frame(
  Samples = "X778_palbo_72h",
  High_MFE_of_total_control = X778_PALBO_dexseq_ANNOTATED_control_high_MFE_TX_counts_fraction,
  High_MFE_of_total_drug_treated = X778_PALBO_dexseq_ANNOTATED_palbociclib_high_MFE_TX_counts_fraction,
  Fisher_result = X778_PALBO_z_test_counts_75_p_value
)

print(Data_high_MFE_fraction_total_counts)

# Load required libraries
library(ggplot2)
library(tidyr)

# Set factor level for sample ordering
Data_high_MFE_fraction_total_counts$Samples <- factor(Data_high_MFE_fraction_total_counts$Samples,
                                                      levels = "X778_palbo_72h")

# Reshape data to long format for plotting
df_long_high_MFE <- pivot_longer(Data_high_MFE_fraction_total_counts,
                                 cols = starts_with("High_MFE_of_total"),
                                 names_to = "Condition",
                                 values_to = "Value")

# Define fill color and label mapping
fill_colors <- c("High_MFE_of_total_control" = "#8cc5e3",
                 "High_MFE_of_total_drug_treated" = "#2066a8")

condition_labels <- c("High_MFE_of_total_control" = "Untreated",
                      "High_MFE_of_total_drug_treated" = "Drug-treated")

# Create bar plot
Plot_high_fraction_tot_counts <- ggplot(df_long_high_MFE, aes(x = Samples, y = Value, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = round(Value, 3)), 
            position = position_dodge(width = 0.8), 
            vjust = 0.4, hjust = -0.3, angle = 90, size = 4) +
  labs(x = NULL,
       y = "",
       fill = NULL) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),
        axis.title.y = element_text(size = 14),
        legend.position = "none") +
  ylim(0, 0.6) +
  scale_fill_manual(values = fill_colors, labels = condition_labels)

# Display plot
print(Plot_high_fraction_tot_counts)











Data_low_MFE_fraction_total_counts <- data.frame(
  Samples = "X778_palbo_72h",
  low_MFE_of_total_control = X778_PALBO_dexseq_ANNOTATED_control_low_MFE_TX_counts_fraction,
  low_MFE_of_total_drug_treated = X778_PALBO_dexseq_ANNOTATED_palbociclib_low_MFE_TX_counts_fraction,
  Fisher_result = X778_PALBO_dexseq_ANNOTATED_fisher_result_counts_low$p.value
)

print(Data_low_MFE_fraction_total_counts)

# Load required libraries
library(ggplot2)
library(tidyr)

# Set factor level for sample ordering
Data_low_MFE_fraction_total_counts$Samples <- factor(Data_low_MFE_fraction_total_counts$Samples,
                                                     levels = "X778_palbo_72h")

# Reshape data to long format for plotting
df_long_low_MFE <- pivot_longer(Data_low_MFE_fraction_total_counts,
                                cols = starts_with("low_MFE_of_total"),
                                names_to = "Condition",
                                values_to = "Value")

# Define fill color and label mapping
fill_colors <- c("low_MFE_of_total_control" = "#8cc5e3",
                 "low_MFE_of_total_drug_treated" = "#2066a8")

condition_labels <- c("low_MFE_of_total_control" = "Untreated",
                      "low_MFE_of_total_drug_treated" = "Drug-treated")

# Create bar plot
Plot_low_fraction_tot_counts <- ggplot(df_long_low_MFE, aes(x = Samples, y = Value, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = round(Value, 3)), 
            position = position_dodge(width = 0.8), 
            vjust = 0.4, hjust = -0.3, angle = 90, size = 4) +
  labs(x = NULL,
       y = "",
       fill = NULL) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),
        axis.title.y = element_text(size = 14),
        legend.position = "none") +
  ylim(0, 0.6) +
  scale_fill_manual(values = fill_colors, labels = condition_labels)

# Display plot
print(Plot_low_fraction_tot_counts)


























# Fig.2H Science models proliferation ####

library(gridExtra)
library(grid)


# Combine all the data into a single data frame
combine_log2FoldChange <- function(data_frames) {
  # Find the maximum length of all log2FoldChange vectors
  max_length <- max(sapply(data_frames, function(df) length(df$log2FoldChange)))
  
  # Pad each log2FoldChange vector to the maximum length with NA
  padded_log2FC <- lapply(data_frames, function(df) {
    len <- length(df$log2FoldChange)
    c(df$log2FoldChange, rep(NA, max_length - len))
  })
  
  # Combine the padded vectors into a data frame and name the columns
  combined_df <- as.data.frame(do.call(cbind, padded_log2FC))
  colnames(combined_df) <- names(data_frames)
  
  return(combined_df)
}

# Generate the data frame
data_frames_Proliferation_Science_mode_systems <- list(
  X778_palbo_72h_Signt_1 = df_X778_PALBO_Proliferation_Science,
  X778_palbo_EP_Signt_1 = df_X778_PALBO_Proliferation_Science_EP,
  X778_palbo_72h_Signt_2 = df_X778_PALBO_Meta_PCNA,
  X778_palbo_EP_Signt_2 = df_X778_PALBO_Meta_PCNA_EP,
  X778_palbo_72h_Signt_3 = df_X778_PALBO_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019,
  X778_palbo_EP_Signt_3 = df_X778_PALBO_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019_EP
  )

# Rename the list elements
names(data_frames_Proliferation_Science_mode_systems) <- c(
  "WDLPS_72h_Signt_1", 
  "WDLPS_EP_Signt_1", 
  "WDLPS_72h_Signt_2", 
  "WDLPS_EP_Signt_2", 
  "WDLPS_72h_Signt_3", 
  "WDLPS_EP_Signt_3"
)


combined_df_Proliferation_Science_model_systems <- combine_log2FoldChange(data_frames_Proliferation_Science_mode_systems)



# Calculate the mean for each sample
means <- colMeans(combined_df_Proliferation_Science_model_systems, na.rm = TRUE)

# Generate the box plot with the extended left margin
Proliferation_Science_model_systems <- ggplot(stack(combined_df_Proliferation_Science_model_systems), aes(x = ind, y = values)) +
  geom_boxplot(aes(fill = ifelse(ind %in% c("WDLPS_72h_Signt_1", "WDLPS_72h_Signt_2", "WDLPS_72h_Signt_3"), "lightgrey", "gray45"))) + 
  stat_summary(fun = mean, geom = "point", shape = 20, size = 2, color = "black") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1)) +
  theme(axis.text.y = element_text(size = 10)) +
  xlab("") + ylab("Log2_FC Cell Proliferation") +
  ylim(-6, 4) +
  geom_hline(yintercept = 0, linetype = 1, color = "black") +
  #ggtitle("Proliferation_Science Wilcox p < 0.05 targeted agents") +
  theme(plot.title = element_text(color = "black", size = 14, face = "bold")) +
  scale_x_discrete(limits = c("WDLPS_72h_Signt_1", "WDLPS_EP_Signt_1", "WDLPS_72h_Signt_2","WDLPS_EP_Signt_2", "WDLPS_72h_Signt_3", "WDLPS_EP_Signt_3")) +
  
  scale_fill_identity() +  # Use identity scale for fill
  theme(plot.margin = unit(c(0, 0, 0, 1), "cm"))  # Adjust the margin (top, right, bottom, left)

# Print the plot
print(Proliferation_Science_model_systems)







# Wilcoxon test proliferation gene sets Science models - individual samples
# IMPORTANT: some of the samples may not have a significant repression of proliferation gene sets, therefore, they may be absent in the final plots
Wilcox_test_Proliferation_Science_models <- data.frame(
  Samples = c("WDLPS_palbo_72h_Signt_1", 
              "WDLPS_palbo_EP_Signt_1", 
              "WDLPS_palbo_72h_Signt_2", 
              "WDLPS_palbo_EP_Signt_2", 
              "WDLPS_palbo_72h_Signt_3", 
              "WDLPS_palbo_EP_Signt_3"),
  wilcox_test_Proliferation_Science = c(wilcox_test_Proliferation_Science_X778_PALBO_72h$p.value, 
                                        wilcox_test_Proliferation_Science_X778_PALBO_EP$p.value), 
  wilcox_test_Meta_PCNA = c(wilcox_test_Meta_PCNA_X778_PALBO_72h$p.value, 
                            wilcox_test_Meta_PCNA_X778_PALBO_EP$p.value),
  wilcox_test_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019 = c(wilcox_test_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019_X778_PALBO_72h$p.value, 
                                                                       wilcox_test_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019_X778_PALBO_EP$p.value))
print(Wilcox_test_Proliferation_Science_models)






# Check whether the difference between the proliferation signatures in aggreate is significant between the 72h and the EP time points
wilcox_test_results <- list(
  wilcox_test_Proliferation_Science_Signt_1 = wilcox.test(combined_df_Proliferation_Science_model_systems$WDLPS_72h_Signt_1, combined_df_Proliferation_Science_model_systems$WDLPS_EP_Signt_1),
  wilcox_test_Proliferation_Science_Signt_2 = wilcox.test(combined_df_Proliferation_Science_model_systems$WDLPS_72h_Signt_2, combined_df_Proliferation_Science_model_systems$WDLPS_EP_Signt_2),
  wilcox_test_Proliferation_Science_Signt_3 = wilcox.test(combined_df_Proliferation_Science_model_systems$WDLPS_72h_Signt_3, combined_df_Proliferation_Science_model_systems$WDLPS_EP_Signt_3)
)

# Extract p-values from each Wilcoxon test result
wilcox_p_values <- sapply(wilcox_test_results, function(test) test$p.value)

# Combine the p-values into a data frame for easier viewing
p_value_df <- data.frame(
  Signature = c("Signt_1", "Signt_2", "Signt_3"),
  P_Value = wilcox_p_values
)

# Print the p-value results
print(p_value_df)























# Fig.2I Science models 72h vs early phase violin plots ####

# IMPORTANT:
# Only the palbociclib condition has been profiled during the early phase response with the Nanopore sequencing technology
# The SKBR3 and the SKMEL28 do not show a significant repression of the proliferative signatures, for this reason they have been not shown

# Plot_X778_PALBO_PALBO_72h <- ggplot(X778_PALBO_dexseq_ANNOTATED, aes(x = factor(log2fold_treated_untreated > 0), y = Normalized_MFE)) +
#   geom_violin(fill = "lightgrey", size = 0.6) +
#   geom_boxplot(width = 0.1, color = "black", outlier.size = 0, size = 0.3) +
#   #eom_jitter(width = 0.1, size = 1, alpha = 0.3) +
#   theme_classic() +
#   theme(legend.position = "none", axis.title = element_text(size = 20), plot.title = element_text(size = 10),
#         axis.text.x = element_text(size = 16, color = "black")) +
#   geom_vline(xintercept = 0, size = 0.3) +
#   ylim(-0.6, 0) +
#   labs(y = "Normalized_MFE", x = NULL, title = "778_PALBO_PALBO_72h_Log2_FC_adj.P<0.05") +
#   scale_x_discrete(labels = c("Log2_FC < 0", "Log2_FC > 0"))
# print(Plot_X778_PALBO_PALBO_72h)
# # Save the combined plot
# ggsave("Plot_X778_PALBO_PALBO_72h.pdf", plot = Plot_X778_PALBO_PALBO_72h, width = 8, height = 8, path = "/Users/z3514253/Desktop")
# 
# 
# ks_test_MFE_X778_PALBO_LAP_72h_all_UP_vs_DOWN







Plot_X778_PALBO_EP <- ggplot(X778_counts_EP, aes(x = factor(log2fold_treated_untreated > 0), y = Normalized_MFE)) +
  geom_violin(fill = "lightgrey", size = 0.6) +
  geom_boxplot(width = 0.1, color = "black", outlier.size = 0, size = 0.3) +
  #geom_jitter(width = 0.1, size = 1, alpha = 0.3) +
  theme_classic() +
  theme(legend.position = "none", axis.title = element_text(size = 20), plot.title = element_text(size = 10),
        axis.text.x = element_text(size = 16, color = "black")) +
  geom_vline(xintercept = 0, size = 0.3) +
  ylim(-0.6, 0) +
  labs(y = "", x = NULL, title = "") +
  scale_x_discrete(labels = c("Log2_FC < 0", "Log2_FC > 0"))
print(Plot_X778_PALBO_EP)

# Save the combined plot
ggsave("Plot_X778_PALBO_EP.pdf", plot = Plot_X778_PALBO_EP, width = 6, height = 8, path = "/Users/z3514253/Desktop")


ttest_MFE_X778_PALBO_PALBO_72h_all_UP_EP_vs_DOWN_EP
ks_test_MFE_X778_PALBO_PALBO_72h_all_UP_EP_vs_DOWN_EP








# Combine all plots
# plots <- list(
#   Plot_X778_PALBO_PALBO_72h,
#   Plot_X778_PALBO_EP
# )
# 
# # Add margin to each plot
# plots_with_margin <- lapply(plots, function(plot) {
#   plot + 
#     # geom_violin(fill = "lightgrey", color = "black", size = 0.1) +
#     #geom_jitter(position = position_jitter(width = 0.1), size = 0.05, color = "black") +
#     # geom_boxplot(width = 0.1, color = "black", outlier.size = 0) +
#     #scale_color_manual(values = c("FALSE" = "black", "TRUE" = "black")) +  # Update colors
#     theme(
#       plot.margin = unit(c(5, 0.2, 0.2, 0.2), "cm"),
#       axis.title = element_blank(),  # Remove y-axis text)  # Adjust the margin values as needed
#       axis.text.x = element_text(size = 14), # Set x-axis text size to 12
#       plot.title = element_blank() # Remove titles
#     )
# })
# 
# # Arrange plots in 5 columns and 2 rows with spacing
# combined_plot <- grid.arrange(
#   grobs = plots_with_margin,
#   ncol = 2,
#   nrow = 1
# )
# 
# # Save the combined plot
# ggsave("combined_plot.pdf", plot = combined_plot, width = 8, height = 7, path = "/Users/z3514253/Desktop")







# Generate the data frame
results <- data.frame(
  Dataset = c(
    "778_PALBO_72h",
    "778_PALBO_EP"
  ),
  Mean_MFE_UP = c(
    Mean_X778_PALBO_dexseq_ANNOTATED_UP,
    mean(X778_PALBO_EP_dexseq_ANNOTATED_UP_EP$Normalized_MFE)
  ),
  Mean_MFE_DOWN = c(
    Mean_X778_PALBO_dexseq_ANNOTATED_DOWN,
    mean(X778_PALBO_EP_dexseq_ANNOTATED_DOWN_EP$Normalized_MFE)
  ),
  KS_Test_p_value = c(
    ks_test_MFE_X778_PALBO_LAP_72h_all_UP_vs_DOWN$p.value,
    ks_test_MFE_X778_PALBO_UP_vs_DOWN_EP$p.value
  )
)

# Display the data frame
print(results)





# Wilcoxon rank-sum test for X778 PALBO 72h (one-sided: up > down)
wilcox_test_X778_PALBO_72h <- wilcox.test(
  X778_PALBO_dexseq_ANNOTATED_UP$Normalized_MFE,
  X778_PALBO_dexseq_ANNOTATED_DOWN$Normalized_MFE,
  paired = FALSE,
  alternative = "greater"
)

# Wilcoxon rank-sum test for X778 PALBO EP (one-sided: up > down)
wilcox_test_X778_PALBO_EP <- wilcox.test(
  X778_PALBO_EP_dexseq_ANNOTATED_UP_EP$Normalized_MFE,
  X778_PALBO_EP_dexseq_ANNOTATED_DOWN_EP$Normalized_MFE,
  paired = FALSE,
  alternative = "greater"
)



wilcoxon_results <- data.frame(
  Dataset = c("778_PALBO_72h", "778_PALBO_EP"),
  Wilcoxon_p_value = c(
    wilcox_test_X778_PALBO_72h$p.value,
    wilcox_test_X778_PALBO_EP$p.value
  )
)

# Display results
print(wilcoxon_results)




















