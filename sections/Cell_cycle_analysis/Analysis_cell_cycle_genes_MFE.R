# CELL CYCLE TRANSCRIPTOME–THERMODYNAMICS PIPELINE
# Aim:
#   To investigate how RNA folding energetics (Minimum Free Energy [MFE] and heat capacity) couple with 
#   differential gene expression across cell cycle phases in two human cancer cell lines: MCF-7 (breast cancer) 
#   and HeLa (cervical cancer).
#
# Datasets:
#   - MCF-7: GSE94479 (Liu et al., PNAS 2017), transcriptome profiling across G0/G1, S, and G2/M phases.
#   - HeLa: GSE104736 (Boström et al., PLOS ONE 2017), FUCCI-sorted cells in G1, S, and G2 phases with three replicates/condition.
#   - RNA structural annotations (Table_ANNOTATED) from RNAfold, providing transcript-level Normalized_MFE (ΔGfold/nt) 
#     and Specific_heat.
#
# Workflow:
#   1) Data preparation
#        - Load raw or processed counts from GREIN/GEO.
#        - Harmonize sample metadata (cell line, condition, phase).
#        - Construct DESeq2 datasets (design ~ condition) and perform QC (sample clustering, PCA, heatmaps).
#   2) Differential expression
#        - Perform pairwise contrasts between phases (G0/G1 vs G2/M, S vs G0/G1, G2/M vs S for MCF-7; G1 vs G2, S vs G1, G2 vs S for HeLa).
#        - Prefilter low counts, normalize, run DESeq().
#        - Extract significant transcripts (padj < 0.05) and join with RNA structural metrics.
#   3) Thermodynamic integration
#        - Scatterplots: log2FC vs Normalized_MFE and vs Specific_heat.
#        - Violin + boxplots: distribution of MFE/heat capacity in up- vs down-regulated transcripts.
#        - Statistical testing: group means, two-sample t-tests, Kolmogorov–Smirnov tests, and quantile-tail comparisons (0.25, 0.75).
#   4) Phase-level weighted profiles
#        - Load TPMs and average replicates by phase.
#        - Compute expression-weighted mean Normalized_MFE ± weighted SD for each phase.
#        - Plot phase-wise dot plots with error bars and transcript-level weighted distributions.
#   5) Gene set enrichment (KEGG GSEA)
#        - Collapse to gene-level log2FC (with baseMean weights where relevant).
#        - Map Ensembl IDs → Entrez (biomaRt).
#        - Run clusterProfiler::gseKEGG per contrast.
#        - Visualize pathway-level results with dot plots, NES heatmaps, and bubble plots.



library(dplyr)
library(readr)
library(ggplot2)
library(scales)
library(gtools)
library(tibble)
library(reshape2)
library(tidyverse)
library(ggpubr)
library(tidyr)
library(biomaRt)

Gencode_TX_norm_MFE_HC_table <- read.delim("path_to/Gencode_TX_norm_MFE_HC_table.csv",
                                           sep = ",", header = T, stringsAsFactors = F)
Table_ANNOTATED <- Gencode_TX_norm_MFE_HC_table




# Cell cycle MCF-7 cell line ####
# Data have been obtained from: Yin Liu et al. PNAS 2017 'Transcriptional landscape of the human cell cycle'. In this paper the authors identify the genes expressed in each phase of the cell cycle.


# IMPORTANT
# These data have been analysed in the GREIN database GSE94479, but there are only single replicates for the different phases of the cell cycle (G0/G1, S and M phases), therefore it is not possible to perform differential expression analysis using the DESeq2 package.
# I have accessed the original data in Gene Expression Omnibus (GEO) and I have downloaded the differential gene expression data that are present in the GEO database.
# In GEO it is also possible to download the gene-level normalized counts of the MCF7 cell line profiled in G0/G1, S and M phases. There are two replicates for each cell cycle phase. However this is useless as DESeq2 requires the raw counts for differential gene expression analysis.


# Load the script which gives the differential gene expression data of the cell cycle phases
source("path_to/Analysis_RNAseq_PNAS_2017.R")


## Differential gene expression G0/G1 vs G2/M ####
# IMPORTANT:
# The DeSeq2 analysis has been done by contrasting G2M vs G0G1, we need the G0G1 vs G2M. 
# By multiplying by -1 will invert the direction of the fold changes, effectively switching the comparison.
# and filter the significant expression changes and join with table containing the MFE values
GSE94479_G0G1vsG2M <- dplyr::filter(GSE94479_Deseq_G0G1vsG2M_IDs, padj < 0.05) %>%
  dplyr::left_join(Table_ANNOTATED, by="ensembl_transcript_id") %>%
  dplyr::filter(!is.na(Normalized_MFE))


# Generate a scatter plot Log2_FC vs Normalized_MFE
ggplot(GSE94479_G0G1vsG2M, aes(y=Normalized_MFE, x=log2FoldChange)) +
  geom_point(size=1.5, alpha=0.3) +
  theme_classic() +
  theme(legend.position="none") +
  geom_vline(xintercept = 0, linewidth=0.3) +
  #xlim(-5.5,5.5) +
  theme(axis.title = element_text(size = 20)) +
  theme(plot.title = element_text(size = 15)) +
  labs(y="Normalized_MFE",x="Log2_FC_expression",title="G0G1vsG2M_MCF7")


# Generate a violin plot of the Normalized_MFE values of the genes up-regulated in G0G1 vs G2M
MCF7_G0G1_vs_G2M_plot <- ggplot(GSE94479_G0G1vsG2M, aes(x = factor(log2FoldChange > 0), y = Normalized_MFE)) +
  geom_violin(fill = "lightgrey", color = "black", size = 0.5) +  # Set violin fill to grey
  geom_boxplot(width = 0.1, color = "black", fill = "white", outlier.shape = NA) +  # Add white boxplot inside violins
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
    title = "MCF7_G0G1_vs_G2M_Log2_FC_adj.P<0.05"
  ) +
  scale_x_discrete(labels = c("log2_FC < 0", "log2_FC > 0"))

print(MCF7_G0G1_vs_G2M_plot)

# Save the combined plot
ggsave("MCF7_G0G1_vs_G2M_plot.pdf", plot = MCF7_G0G1_vs_G2M_plot, width = 5, height = 6, path = "path_to/Desktop")




# Statistical analysis of the MFE values of the genes up-regulated in G0G1 vs G2M
MFE_MCF7_G0G1vsG2M_UP <- dplyr::filter(GSE94479_G0G1vsG2M, log2FoldChange > 0)
Mean_MFE_MCF7_G0G1vsG2M_UP <- mean(MFE_MCF7_G0G1vsG2M_UP$Normalized_MFE, na.rm = TRUE)

MFE_MCF7_G0G1vsG2M_DOWN <- dplyr::filter(GSE94479_G0G1vsG2M, log2FoldChange < 0)
Mean_MFE_MCF7_G0G1vsG2M_DOWN <- mean(MFE_MCF7_G0G1vsG2M_DOWN$Normalized_MFE, na.rm = TRUE)

TTEST_MFE_MCF7_G0G1vsG2M_all_UP_vs_DOWN <- t.test(MFE_MCF7_G0G1vsG2M_DOWN$Normalized_MFE, MFE_MCF7_G0G1vsG2M_UP$Normalized_MFE, paired=FALSE)
print(TTEST_MFE_MCF7_G0G1vsG2M_all_UP_vs_DOWN)



# Add K-S test to look at the differences between the up-regukated and down-regulated genes
ks_test_MCF7_G0G1vsG2M_all_UP_vs_DOWN <- ks.test(MFE_MCF7_G0G1vsG2M_DOWN$Normalized_MFE, MFE_MCF7_G0G1vsG2M_UP$Normalized_MFE)
print(ks_test_MCF7_G0G1vsG2M_all_UP_vs_DOWN)

quantile_075_MFE_MCF7_G0G1vsG2M_UP <- quantile(MFE_MCF7_G0G1vsG2M_UP$Normalized_MFE, probs = c(0.75)) 
quantile_075_MFE_MCF7_G0G1vsG2M_DOWN <- quantile(MFE_MCF7_G0G1vsG2M_DOWN$Normalized_MFE, probs = c(0.75)) 
# Values above 0.75 quantile
ks_test_MFE_MCF7_G0G1vsG2M_075_above <- ks.test(MFE_MCF7_G0G1vsG2M_UP$Normalized_MFE[MFE_MCF7_G0G1vsG2M_UP$Normalized_MFE >= quantile_075_MFE_MCF7_G0G1vsG2M_UP],
                                           MFE_MCF7_G0G1vsG2M_DOWN$Normalized_MFE[MFE_MCF7_G0G1vsG2M_DOWN$Normalized_MFE >= quantile_075_MFE_MCF7_G0G1vsG2M_DOWN])
print(ks_test_MFE_MCF7_G0G1vsG2M_075_above)

# Look at the differences between the 0.25 quantile 
quantile_025_MFE_MCF7_G0G1vsG2M_UP <- quantile(MFE_MCF7_G0G1vsG2M_UP$Normalized_MFE, probs = c(0.25)) 
quantile_025_MFE_MCF7_G0G1vsG2M_DOWN <- quantile(MFE_MCF7_G0G1vsG2M_DOWN$Normalized_MFE, probs = c(0.25)) 
# Values below 0.25 quantile
ks_test_MF_MCF7_G0G1vsG2M_025_below <- ks.test(MFE_MCF7_G0G1vsG2M_UP$Normalized_MFE[MFE_MCF7_G0G1vsG2M_UP$Normalized_MFE <= quantile_025_MFE_MCF7_G0G1vsG2M_UP],
                                          MFE_MCF7_G0G1vsG2M_DOWN$Normalized_MFE[MFE_MCF7_G0G1vsG2M_DOWN$Normalized_MFE <= quantile_025_MFE_MCF7_G0G1vsG2M_DOWN])
print(ks_test_MF_MCF7_G0G1vsG2M_025_below)







# Differential gene expression G0G1 vs G2M heat capacity changes

# This is calculated on the total n. of transcripts with a significant log2 fold change
# Generate a scatter plot Log2_FC vs heat capacity
ggplot(GSE94479_G0G1vsG2M, aes(y=Cp_37C, x=log2FoldChange)) +
  geom_point(size=1.5, alpha=0.1)+
  theme_classic() +
  theme(legend.position="none") +
  geom_vline(xintercept = 0, size=0.3) +
  xlim(-12, 12) +
  ylim(0,400) +
  theme(axis.title = element_text(size = 20)) +
  theme(plot.title = element_text(size = 15)) +
  labs(y="Heat capacity (J/C)",x="Log2_FC_expression",title="G0G1vsG2M_MCF7")


# Show the data with a box plot
Plot_GSE94479_G0G1vsG2M_heat_capacity <- ggplot(GSE94479_G0G1vsG2M, aes(x = factor(log2FoldChange > 0), y = Cp_37C)) +
  geom_boxplot(fill = "lightgrey", alpha = 0.7) +
  #geom_jitter(width = 0.1, size = 1, alpha = 0.3) +
  theme_classic() +
  ylim(0, 400) +
  theme(legend.position = "none", 
        axis.title = element_text(size = 15), 
        plot.title = element_text(size = 10),
        axis.text.x = element_text(size = 10, color = "black")) +
  theme(plot.title = element_text(size = 20)) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  labs(y="Heat capacity (J/C)",x="",title="G0G1vsG2M_MCF7") +
  scale_x_discrete(labels = c("log2FoldChange > 0", "log2FoldChange < 0"))

print(Plot_GSE94479_G0G1vsG2M_heat_capacity)
ggsave("Plot_GSE94479_G0G1vsG2M_heat_capacity.pdf", plot = Plot_GSE94479_G0G1vsG2M_heat_capacity, width = 5, height = 6, path = "path_to/Desktop")




Mean_heat_capacity_MFE_MCF7_G0G1vsG2M_UP <- mean(MFE_MCF7_G0G1vsG2M_UP$Cp_37C, na.rm = TRUE)
Mean_heat_capacity_MFE_MCF7_G0G1vsG2M_DOWN <- mean(MFE_MCF7_G0G1vsG2M_DOWN$Cp_37C, na.rm = TRUE)
TTEST_heat_capacity_MFE_MCF7_G0G1vsG2M_all_UP_vs_DOWN <- t.test(MFE_MCF7_G0G1vsG2M_DOWN$Cp_37C, MFE_MCF7_G0G1vsG2M_UP$Cp_37C, paired=FALSE)
print(TTEST_heat_capacity_MFE_MCF7_G0G1vsG2M_all_UP_vs_DOWN)
# Add K-S test to look at the differences between the up-regukated and down-regulated genes
ks_test_MFE_MCF7_G0G1vsG2M_all_UP_vs_DOWN_heat_capacity <- ks.test(MFE_MCF7_G0G1vsG2M_DOWN$Cp_37C, MFE_MCF7_G0G1vsG2M_UP$Cp_37C)
print(ks_test_MFE_MCF7_G0G1vsG2M_all_UP_vs_DOWN_heat_capacity)

















## Differential gene expression S vs G0/G1 ####

# Open the differential gene expression G2M vs S and filter the significant expression changes and join with table containing the MFE values
GSE94479_SvsG0G1 <- dplyr::filter(GSE94479_Deseq_SvsG0G1_IDs, padj < 0.05) %>%
  dplyr::left_join(Table_ANNOTATED, by="ensembl_transcript_id") %>%
  dplyr::filter(!is.na(Normalized_MFE))

# Generate a scatter plot Log2_FC vs Normalized_MFE
library(ggplot2)
ggplot(GSE94479_SvsG0G1, aes(y=Normalized_MFE, x=log2FoldChange)) +
  geom_point(size=1.5, alpha=0.3) +
  theme_classic() +
  theme(legend.position="none") +
  geom_vline(xintercept = 0, size=0.3) +
  #xlim(-5.5,5.5) +
  theme(axis.title = element_text(size = 20)) +
  theme(plot.title = element_text(size = 15)) +
  labs(y="Normalized_MFE",x="Log2_FC_expression",title="SvsG0G1_MCF-7")



# Generate a violin plot of the Normalized_MFE values of the genes up-regulated in G0G1 vs G2M
MCF7_SvsG0G1_plot <- ggplot(GSE94479_SvsG0G1, aes(x = factor(log2FoldChange > 0), y = Normalized_MFE)) +
  geom_violin(fill = "lightgrey", color = "black", size = 0.5) +  # Set violin fill to grey
  geom_boxplot(width = 0.1, color = "black", fill = "white", outlier.shape = NA) +  # Add white boxplot inside violins
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
    title = "MCF7_S_vsG0G1_Log2_FC_adj.P<0.05"
  ) +
  scale_x_discrete(labels = c("log2_FC < 0", "log2_FC > 0"))

print(MCF7_SvsG0G1_plot)

# Save the combined plot
ggsave("MCF7_SvsG0G1_plot.pdf", plot = MCF7_SvsG0G1_plot, width = 5, height = 6, path = "path_to/Desktop")




# Statistical analysis of the MFE values of the genes up-regulated in G0G1 vs G2M
MFE_MCF7_SvsG0G1_UP <- dplyr::filter(GSE94479_SvsG0G1, log2FoldChange > 0)
Mean_MFE_MCF7_SvsG0G1_UP <- mean(MFE_MCF7_SvsG0G1_UP$Normalized_MFE)

MFE_MCF7_SvsG0G1_DOWN <- dplyr::filter(GSE94479_SvsG0G1, log2FoldChange < 0)
Mean_MFE_MCF7_SvsG0G1_DOWN <- mean(MFE_MCF7_SvsG0G1_DOWN$Normalized_MFE)

TTEST_MFE_MCF7_SvsG0G1_all_UP_vs_DOWN <- t.test(MFE_MCF7_SvsG0G1_DOWN$Normalized_MFE, MFE_MCF7_SvsG0G1_UP$Normalized_MFE, paired=FALSE)
print(TTEST_MFE_MCF7_SvsG0G1_all_UP_vs_DOWN)



# Add K-S test to look at the differences between the up-regukated and down-regulated genes
ks_test_MCF7_SvsG0G1_all_UP_vs_DOWN <- ks.test(MFE_MCF7_SvsG0G1_DOWN$Normalized_MFE, MFE_MCF7_SvsG0G1_UP$Normalized_MFE)
print(ks_test_MCF7_SvsG0G1_all_UP_vs_DOWN)

quantile_075_MFE_MCF7_SvsG0G1_UP <- quantile(MFE_MCF7_SvsG0G1_UP$Normalized_MFE, probs = c(0.75)) 
quantile_075_MFE_MCF7_SvsG0G1_DOWN <- quantile(MFE_MCF7_SvsG0G1_DOWN$Normalized_MFE, probs = c(0.75)) 
# Values above 0.75 quantile
ks_test_MFE_MCF7_SvsG0G1_075_above <- ks.test(MFE_MCF7_SvsG0G1_UP$Normalized_MFE[MFE_MCF7_SvsG0G1_UP$Normalized_MFE >= quantile_075_MFE_MCF7_SvsG0G1_UP],
                                             MFE_MCF7_SvsG0G1_DOWN$Normalized_MFE[MFE_MCF7_SvsG0G1_DOWN$Normalized_MFE >= quantile_075_MFE_MCF7_SvsG0G1_DOWN])
print(ks_test_MFE_MCF7_SvsG0G1_075_above)

# Look at the differences between the 0.25 quantile 
quantile_025_MFE_MCF7_SvsG0G1_UP <- quantile(MFE_MCF7_SvsG0G1_UP$Normalized_MFE, probs = c(0.25)) 
quantile_025_MFE_MCF7_SvsG0G1_DOWN <- quantile(MFE_MCF7_SvsG0G1_DOWN$Normalized_MFE, probs = c(0.25)) 
# Values below 0.25 quantile
ks_test_MF_MCF7_SvsG0G1_025_below <- ks.test(MFE_MCF7_SvsG0G1_UP$Normalized_MFE[MFE_MCF7_SvsG0G1_UP$Normalized_MFE <= quantile_025_MFE_MCF7_SvsG0G1_UP],
                                            MFE_MCF7_SvsG0G1_DOWN$Normalized_MFE[MFE_MCF7_SvsG0G1_DOWN$Normalized_MFE <= quantile_025_MFE_MCF7_SvsG0G1_DOWN])
print(ks_test_MF_MCF7_SvsG0G1_025_below)












# Differential gene expression S vs G0G1 heat capacity changes

# This is calculated on the total n. of transcripts with a significant log2 fold change
# Generate a scatter plot Log2_FC vs heat capacity
ggplot(GSE94479_SvsG0G1, aes(y=Cp_37C, x=log2FoldChange)) +
  geom_point(size=1.5, alpha=0.1)+
  theme_classic() +
  theme(legend.position="none") +
  geom_vline(xintercept = 0, size=0.3) +
  xlim(-12, 12) +
  ylim(0,400) +
  theme(axis.title = element_text(size = 20)) +
  theme(plot.title = element_text(size = 15)) +
  labs(y="Heat capacity (J/C)",x="Log2_FC_expression",title="SvsG0G1_MCF7")


# Show the data with a box plot
Plot_GSE94479_SvsG0G1_heat_capacity <- ggplot(GSE94479_SvsG0G1, aes(x = factor(log2FoldChange > 0), y = Cp_37C)) +
  geom_boxplot(fill = "lightgrey", alpha = 0.7) +
  #geom_jitter(width = 0.1, size = 1, alpha = 0.3) +
  theme_classic() +
  ylim(0, 400) +
  theme(legend.position = "none", 
        axis.title = element_text(size = 20), 
        plot.title = element_text(size = 10),
        axis.text.x = element_text(size = 16, color = "black")) +
  theme(plot.title = element_text(size = 20)) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  labs(y="Heat capacity (J/C)",x="",title="SvsG0G1_MCF7") +
  scale_x_discrete(labels = c("log2FoldChange > 0", "log2FoldChange < 0"))


print(Plot_GSE94479_SvsG0G1_heat_capacity)

# Save the combined plot
ggsave("Plot_GSE94479_SvsG0G1_heat_capacity.pdf", plot = Plot_GSE94479_SvsG0G1_heat_capacity, width = 5, height = 6, path = "path_to/Desktop")



Mean_heat_capacity_MFE_MCF7_SvsG0G1_UP <- mean(MFE_MCF7_SvsG0G1_UP$Cp_37C, na.rm = TRUE)
Mean_heat_capacity_MFE_MCF7_SvsG0G1_DOWN <- mean(MFE_MCF7_SvsG0G1_DOWN$Cp_37C, na.rm = TRUE)
TTEST_heat_capacity_MFE_MCF7_SvsG0G1_all_UP_vs_DOWN <- t.test(MFE_MCF7_SvsG0G1_DOWN$Cp_37C, MFE_MCF7_SvsG0G1_UP$Cp_37C, paired=FALSE)
print(TTEST_heat_capacity_MFE_MCF7_SvsG0G1_all_UP_vs_DOWN)
# Add K-S test to look at the differences between the up-regukated and down-regulated genes
ks_test_MFE_MCF7_SvsG0G1_all_UP_vs_DOWN_heat_capacity <- ks.test(MFE_MCF7_SvsG0G1_DOWN$Cp_37C, MFE_MCF7_SvsG0G1_UP$Cp_37C)
print(ks_test_MFE_MCF7_SvsG0G1_all_UP_vs_DOWN_heat_capacity)













## Differential gene expression G2/M vs S ####

# Open the differential gene expression G2M vs S and filter the significant expression changes and join with table containing the MFE values
GSE94479_G2MvsS <- dplyr::filter(GSE94479_Deseq_G2MvsS_IDs, padj < 0.05) %>%
  dplyr::left_join(Table_ANNOTATED, by="ensembl_transcript_id") %>%
  dplyr::filter(!is.na(Normalized_MFE))


# Generate a scatter plot Log2_FC vs Normalized_MFE
library(ggplot2)
ggplot(GSE94479_G2MvsS, aes(y=Normalized_MFE, x=log2FoldChange)) +
  geom_point(size=1.5, alpha=0.3) +
  theme_classic() +
  theme(legend.position="none") +
  geom_vline(xintercept = 0, size=0.3) +
  #xlim(-5.5,5.5) +
  theme(axis.title = element_text(size = 20)) +
  theme(plot.title = element_text(size = 15)) +
  labs(y="Normalized_MFE",x="Log2_FC_expression",title="G2MvsS_MCF7")



# Generate a violin plot of the Normalized_MFE values of the genes up-regulated in G0G1 vs G2M
MCF7_G2MvsS_plot <- ggplot(GSE94479_G2MvsS, aes(x = factor(log2FoldChange > 0), y = Normalized_MFE)) +
  geom_violin(fill = "lightgrey", color = "black", size = 0.5) +  # Set violin fill to grey
  geom_boxplot(width = 0.1, color = "black", fill = "white", outlier.shape = NA) +  # Add white boxplot inside violins
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
    title = "MCF7_G2M_vs_S_Log2_FC_adj.P<0.05"
  ) +
  scale_x_discrete(labels = c("log2_FC < 0", "log2_FC > 0"))

print(MCF7_G2MvsS_plot)

# Save the combined plot
ggsave("MCF7_G2MvsS_plot.pdf", plot = MCF7_G2MvsS_plot, width = 5, height = 6, path = "path_to/Desktop")





# Statistical analysis of the MFE values of the genes up-regulated in G0G1 vs G2M
MFE_MCF7_G2MvsS_UP <- dplyr::filter(GSE94479_G2MvsS, log2FoldChange > 0)
Mean_MFE_MCF7_G2MvsS_UP <- mean(MFE_MCF7_G2MvsS_UP$Normalized_MFE)

MFE_MCF7_G2MvsS_DOWN <- dplyr::filter(GSE94479_G2MvsS, log2FoldChange < 0)
Mean_MFE_MCF7_G2MvsS_DOWN <- mean(MFE_MCF7_G2MvsS_DOWN$Normalized_MFE)

TTEST_MFE_MCF7_G2MvsS_all_UP_vs_DOWN <- t.test(MFE_MCF7_G2MvsS_DOWN$Normalized_MFE, MFE_MCF7_G2MvsS_UP$Normalized_MFE, paired=FALSE)
print(TTEST_MFE_MCF7_G2MvsS_all_UP_vs_DOWN)



# Add K-S test to look at the differences between the up-regukated and down-regulated genes
ks_test_MCF7_G2MvsS_all_UP_vs_DOWN <- ks.test(MFE_MCF7_G2MvsS_DOWN$Normalized_MFE, MFE_MCF7_G2MvsS_UP$Normalized_MFE)
print(ks_test_MCF7_G2MvsS_all_UP_vs_DOWN)

quantile_075_MFE_MCF7_G2MvsS_UP <- quantile(MFE_MCF7_G2MvsS_UP$Normalized_MFE, probs = c(0.75)) 
quantile_075_MFE_MCF7_G2MvsS_DOWN <- quantile(MFE_MCF7_G2MvsS_DOWN$Normalized_MFE, probs = c(0.75)) 
# Values above 0.75 quantile
ks_test_MFE_MCF7_G2MvsS_075_above <- ks.test(MFE_MCF7_G2MvsS_UP$Normalized_MFE[MFE_MCF7_G2MvsS_UP$Normalized_MFE >= quantile_075_MFE_MCF7_G2MvsS_UP],
                                             MFE_MCF7_G2MvsS_DOWN$Normalized_MFE[MFE_MCF7_G2MvsS_DOWN$Normalized_MFE >= quantile_075_MFE_MCF7_G2MvsS_DOWN])
print(ks_test_MFE_MCF7_G2MvsS_075_above)

# Look at the differences between the 0.25 quantile 
quantile_025_MFE_MCF7_G2MvsS_UP <- quantile(MFE_MCF7_G2MvsS_UP$Normalized_MFE, probs = c(0.25)) 
quantile_025_MFE_MCF7_G2MvsS_DOWN <- quantile(MFE_MCF7_G2MvsS_DOWN$Normalized_MFE, probs = c(0.25)) 
# Values below 0.25 quantile
ks_test_MF_MCF7_G2MvsS_025_below <- ks.test(MFE_MCF7_G2MvsS_UP$Normalized_MFE[MFE_MCF7_G2MvsS_UP$Normalized_MFE <= quantile_025_MFE_MCF7_G2MvsS_UP],
                                            MFE_MCF7_G2MvsS_DOWN$Normalized_MFE[MFE_MCF7_G2MvsS_DOWN$Normalized_MFE <= quantile_025_MFE_MCF7_G2MvsS_DOWN])
print(ks_test_MF_MCF7_G2MvsS_025_below)










# Differential gene expression G2MvsS heat capacity changes

# This is calculated on the total n. of transcripts with a significant log2 fold change
# Generate a scatter plot Log2_FC vs heat capacity
ggplot(GSE94479_G2MvsS, aes(y=Cp_37C, x=log2FoldChange)) +
  geom_point(size=1.5, alpha=0.1)+
  theme_classic() +
  theme(legend.position="none") +
  geom_vline(xintercept = 0, size=0.3) +
  xlim(-12, 12) +
  ylim(0,400) +
  theme(axis.title = element_text(size = 20)) +
  theme(plot.title = element_text(size = 15)) +
  labs(y="Heat capacity (J/C)",x="Log2_FC_expression",title="G2MvsS_MCF7")


# Show the data with a box plot
Plot_GSE94479_G2MvsS_heat_capacity <- ggplot(GSE94479_G2MvsS, aes(x = factor(log2FoldChange > 0), y = Cp_37C)) +
  geom_boxplot(fill = "lightgrey", alpha = 0.7) +
  #geom_jitter(width = 0.1, size = 1, alpha = 0.3) +
  theme_classic() +
  ylim(0, 400) +
  theme(legend.position = "none", 
        axis.title = element_text(size = 20), 
        plot.title = element_text(size = 10),
        axis.text.x = element_text(size = 16, color = "black")) +
  theme(plot.title = element_text(size = 20)) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  labs(y="Heat capacity (J/C)",x="",title="G2MvsS_MCF7") +
  scale_x_discrete(labels = c("log2FoldChange > 0", "log2FoldChange < 0"))

print(Plot_GSE94479_G2MvsS_heat_capacity)

# Save the combined plot
ggsave("Plot_GSE94479_G2MvsS_heat_capacity.pdf", plot = Plot_GSE94479_G2MvsS_heat_capacity, width = 5, height = 6, path = "path_to/Desktop")

Mean_heat_capacity_MFE_MCF7_G2MvsS_UP <- mean(MFE_MCF7_G2MvsS_UP$Cp_37C, na.rm = TRUE)
Mean_heat_capacity_MFE_MCF7_G2MvsS_DOWN <- mean(MFE_MCF7_G2MvsS_DOWN$Cp_37C, na.rm = TRUE)
TTEST_heat_capacity_MFE_MCF7_G2MvsS_all_UP_vs_DOWN <- t.test(MFE_MCF7_G2MvsS_DOWN$Cp_37C, MFE_MCF7_G2MvsS_UP$Cp_37C, paired=FALSE)
print(TTEST_heat_capacity_MFE_MCF7_G2MvsS_all_UP_vs_DOWN)
# Add K-S test to look at the differences between the up-regukated and down-regulated genes
ks_test_MFE_MCF7_G2MvsS_all_UP_vs_DOWN_heat_capacity <- ks.test(MFE_MCF7_G2MvsS_DOWN$Cp_37C, MFE_MCF7_G2MvsS_UP$Cp_37C)
print(ks_test_MFE_MCF7_G2MvsS_all_UP_vs_DOWN_heat_capacity)




















# GSEA cell proliferation and mitogenic signaling cell cycle phases ####

# G2M vs G0G1
GSE94479_Deseq_G0G1vsG2M_IDs_gene_level_df <- GSE94479_Deseq_G0G1vsG2M_IDs %>%
  dplyr::filter(!is.na(log2FoldChange), !is.na(baseMean), baseMean > 0) %>%
  dplyr::group_by(ensembl_gene_id) %>%
  dplyr::summarise(
    weighted_log2FC = sum(log2FoldChange * baseMean) / sum(baseMean),
    min_pvalue = min(pvalue, na.rm = TRUE),
    mean_baseMean = mean(baseMean, na.rm = TRUE),
    n_transcripts = n()
  ) %>%
  ungroup() %>%
  dplyr::filter(!is.na(ensembl_gene_id))


# Check uniqueness of ensembl_gene_id
ensembl_unique <- sum(duplicated(GSE94479_Deseq_G0G1vsG2M_IDs_gene_level_df$ensembl_gene_id))
cat("Number of duplicated ensembl_gene_id:", ensembl_unique, "\n")



# S vs G0G1
GSE94479_Deseq_SvsG0G1_IDs_gene_level_df <- GSE94479_Deseq_SvsG0G1_IDs %>%
  dplyr::filter(!is.na(log2FoldChange), !is.na(baseMean), baseMean > 0) %>%
  dplyr::group_by(ensembl_gene_id) %>%
  dplyr::summarise(
    weighted_log2FC = sum(log2FoldChange * baseMean) / sum(baseMean),
    min_pvalue = min(pvalue, na.rm = TRUE),
    mean_baseMean = mean(baseMean, na.rm = TRUE),
    n_transcripts = n()
  ) %>%
  ungroup() %>%
  dplyr::filter(!is.na(ensembl_gene_id))


# Check uniqueness of ensembl_gene_id
ensembl_unique <- sum(duplicated(GSE94479_Deseq_SvsG0G1_IDs_gene_level_df$ensembl_gene_id))
cat("Number of duplicated ensembl_gene_id:", ensembl_unique, "\n")




# G2M vs S
GSE94479_Deseq_G2MvsS_IDs_gene_level_df <- GSE94479_Deseq_G2MvsS_IDs %>%
  dplyr::filter(!is.na(log2FoldChange), !is.na(baseMean), baseMean > 0) %>%
  dplyr::group_by(ensembl_gene_id) %>%
  dplyr::summarise(
    weighted_log2FC = sum(log2FoldChange * baseMean) / sum(baseMean),
    min_pvalue = min(pvalue, na.rm = TRUE),
    mean_baseMean = mean(baseMean, na.rm = TRUE),
    n_transcripts = n()
  ) %>%
  ungroup() %>%
  dplyr::filter(!is.na(ensembl_gene_id))


# Check uniqueness of ensembl_gene_id
ensembl_unique <- sum(duplicated(GSE94479_Deseq_G2MvsS_IDs_gene_level_df$ensembl_gene_id))
cat("Number of duplicated ensembl_gene_id:", ensembl_unique, "\n")









# Add a column in each data frame with the entrez gene IDs required for the GSEA analysis
library(biomaRt)
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host="https://useast.ensembl.org")

gene_ids <- (GSE94479_Deseq_G0G1vsG2M_IDs_gene_level_df$ensembl_gene_id) 
res_gene_entrez_id_G0G1vsG2M <- getBM(attributes = c('entrezgene_id',
                                           'ensembl_gene_id'),
                            filters = 'ensembl_gene_id',
                            values = gene_ids,
                            mart = mart)
res_gene_entrez_id_G0G1vsG2M <- dplyr::distinct(res_gene_entrez_id_G0G1vsG2M, ensembl_gene_id, .keep_all = TRUE) %>%
  dplyr::filter(!is.na(entrezgene_id)) 
GSE94479_Deseq_G0G1vsG2M_entrez_IDs <- dplyr::left_join(GSE94479_Deseq_G0G1vsG2M_IDs_gene_level_df, res_gene_entrez_id_G0G1vsG2M, by = "ensembl_gene_id")
# Some entrez_id map on multiple ensembl_gene_id
# Collapse to one row per entrezgene_id using the highest absolute log2FC
GSE94479_Deseq_G0G1vsG2M_entrez_IDs <- GSE94479_Deseq_G0G1vsG2M_entrez_IDs %>%
  group_by(entrezgene_id) %>%
  slice_max(order_by = abs(weighted_log2FC), n = 1, with_ties = FALSE) %>%
  ungroup()





gene_ids <- (GSE94479_Deseq_SvsG0G1_IDs_gene_level_df$ensembl_gene_id)
res_gene_entrez_id_SvsG0G1 <- getBM(attributes = c('entrezgene_id',
                                           'ensembl_gene_id'),
                            filters = 'ensembl_gene_id',
                            values = gene_ids,
                            mart = mart)
res_gene_entrez_id_SvsG0G1 <- dplyr::distinct(res_gene_entrez_id_SvsG0G1, ensembl_gene_id, .keep_all = TRUE) %>%
  dplyr::filter(!is.na(entrezgene_id))
GSE94479_Deseq_SvsG0G1_entrez_IDs <- dplyr::left_join(GSE94479_Deseq_SvsG0G1_IDs_gene_level_df, res_gene_entrez_id_SvsG0G1, by = "ensembl_gene_id")
# Some entrez_id map on multiple ensembl_gene_id
# Collapse to one row per entrezgene_id using the highest absolute log2FC
GSE94479_Deseq_SvsG0G1_entrez_IDs <- GSE94479_Deseq_SvsG0G1_entrez_IDs %>%
  group_by(entrezgene_id) %>%
  slice_max(order_by = abs(weighted_log2FC), n = 1, with_ties = FALSE) %>%
  ungroup()







gene_ids <- (GSE94479_Deseq_G2MvsS_IDs_gene_level_df$ensembl_gene_id)
res_gene_entrez_id_G2MvsS <- getBM(attributes = c('entrezgene_id',
                                                   'ensembl_gene_id'),
                                    filters = 'ensembl_gene_id',
                                    values = gene_ids,
                                    mart = mart)
res_gene_entrez_id_G2MvsS <- dplyr::distinct(res_gene_entrez_id_G2MvsS, ensembl_gene_id, .keep_all = TRUE) %>%
  dplyr::filter(!is.na(entrezgene_id))
GSE94479_Deseq_G2MvsS_entrez_IDs <- dplyr::left_join(GSE94479_Deseq_G2MvsS_IDs_gene_level_df, res_gene_entrez_id_G2MvsS, by = "ensembl_gene_id")
# Some entrez_id map on multiple ensembl_gene_id
# Collapse to one row per entrezgene_id using the highest absolute log2FC
GSE94479_Deseq_G2MvsS_entrez_IDs <- GSE94479_Deseq_G2MvsS_entrez_IDs %>%
  group_by(entrezgene_id) %>%
  slice_max(order_by = abs(weighted_log2FC), n = 1, with_ties = FALSE) %>%
  ungroup()











# Generate a data frame with the gene-level Log2 FC values for each comparison
GSE94479_gene_level_Log2_FC_df <- GSE94479_Deseq_G0G1vsG2M_entrez_IDs %>%
  dplyr::select(ensembl_gene_id, entrezgene_id, weighted_log2FC) %>%
  dplyr::rename(G0G1vsG2M_Log_FC = weighted_log2FC) %>%
  dplyr::left_join(
    GSE94479_Deseq_SvsG0G1_entrez_IDs %>%
      dplyr::select(ensembl_gene_id, entrezgene_id, weighted_log2FC) %>%
      dplyr::rename(SvsG0G1_Log_FC = weighted_log2FC),
    by = c("ensembl_gene_id", "entrezgene_id")
  ) %>%
  dplyr::left_join(
    GSE94479_Deseq_G2MvsS_entrez_IDs %>%
      dplyr::select(ensembl_gene_id, entrezgene_id, weighted_log2FC) %>%
      dplyr::rename(G2MvsS_Log_FC = weighted_log2FC),
    by = c("ensembl_gene_id", "entrezgene_id")
  )










# Generate a function which performs GSEA of the samples

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
data_frame <- GSE94479_gene_level_Log2_FC_df
log_fc_cols <- c("G0G1vsG2M_Log_FC", "SvsG0G1_Log_FC", "G2MvsS_Log_FC")

plot_titles <- c("GSE94479_G0G1vsG2M", "GSE94479_SvsG0G1", "GSE94479_G2MvsS")

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

# Write the file in a csv format
write.csv(combined_kk2_df, file = "path_to/combined_kk2_df.csv")









# Generate a heat map with the kk2 results
library(pheatmap)

# Selected KEGG signaling pathways
Selected_KEGG_pathways <- c("hsa04150", "hsa04668", "hsa04657", "hsa03030", "hsa03410", "hsa03430", "hsa03420", "hsa04370", "hsa4072",
                            "hsa03320", "hsa03460", "hsa03440", "hsa04014", "hsa04070", "hsa04020", "hsa01230",
                            "hsa04310", "hsa03082", "hsa04015", "hsa04625", "hsa04371", "hsa04390", "hsa04010", "hsa04622", "hsa04350",
                            "hsa04152", "hsa04926", "hsa04630", "hsa04115", "hsa04151", "hsa04910", "hsa04330", "hsa04012",
                            "hsa03440", "hsa03460") # hsa04110 is the cell cycle gene set

combined_kk2_df_filtered <- dplyr::filter(combined_kk2_df, ID %in% Selected_KEGG_pathways) %>%
  dplyr::mutate(NES = scales::rescale(NES, to = c(-4, 4))) # set a specific range of values to have the 0 centered to white



generate_heatmap <- function(kk2_df) {
  kk2_df <- kk2_df %>%
    dplyr::select(Description, Sample, NES) %>%
    tidyr::spread(Sample, NES)
  
  kk2_df[is.na(kk2_df)] <- 0
  rownames(kk2_df) <- kk2_df$Description
  kk2_df <- kk2_df[, -1]
  
  desired_order <- c("GSE94479_G0G1vsG2M", "GSE94479_SvsG0G1", "GSE94479_G2MvsS")
  kk2_df <- kk2_df[, desired_order]
  
  pathway_order <- c(
    "Wnt signaling pathway",
    "Ras signaling pathway",
    "PI3K-Akt signaling pathway",
    "MAPK signaling pathway",
    "Rap1 signaling pathway",
    "VEGF signaling pathway",
    "Calcium signaling pathway",
    "p53 signaling pathway",
    "DNA replication",
    "Homologous recombination",
    "Fanconi anemia pathway"
  )
  
  kk2_df <- kk2_df[pathway_order, , drop = FALSE]
  
  display_matrix <- as.matrix(kk2_df)
  display_numbers <- ifelse(display_matrix == 0, "", sprintf("%.2g", display_matrix))
  
  heatmap_plot <- pheatmap(as.matrix(kk2_df), 
                           color = colorRampPalette(c("blue", "white", "red"))(50), 
                           cluster_rows = FALSE, 
                           cluster_cols = FALSE, 
                           display_numbers = display_numbers, 
                           number_format = "%.2g", 
                           main = "Heatmap of KEGG signaling pathways",
                           silent = TRUE)
  
  return(list(plot = heatmap_plot, data = kk2_df))
}

# Now you can access both:
heatmap_result <- generate_heatmap(combined_kk2_df_filtered)

# Access the plot
Heatmap_cell_cycle <- heatmap_result$plot

# Access the matrix
heatmap_df <- heatmap_result$data
View(heatmap_df)













# Plot the data
desired_order <- c( "Wnt signaling pathway",
                    "Ras signaling pathway",
                    "PI3K-Akt signaling pathway",
                    "MAPK signaling pathway",
                    "Rap1 signaling pathway",
                    "VEGF signaling pathway",
                    "Calcium signaling pathway",
                    "p53 signaling pathway",
                    "DNA replication",
                    "Homologous recombination",
                    "Fanconi anemia pathway")


# Reorder the 'Description' column based on the desired order
GSE94479_combined_kk2_df_filtered_plot <- combined_kk2_df_filtered %>%
  dplyr::mutate(Description = factor(Description, levels = desired_order),
                neg_log10_p_adjust = -log10(p.adjust))

GSE94479_combined_kk2_df_filtered_plot <- GSE94479_combined_kk2_df_filtered_plot %>%
  dplyr::mutate(Sample = factor(Sample, levels = c("GSE94479_G0G1vsG2M", "GSE94479_SvsG0G1", "GSE94479_G2MvsS")))


# Generate the plot with circles
Plot_cell_proliferation <- ggplot(GSE94479_combined_kk2_df_filtered_plot, aes(x = Sample, y = Description, fill = NES, size = neg_log10_p_adjust)) +
  geom_point(shape = 21, color = "black", stroke = 0.5) +  # Use shape = 21 for circles with border color
  scale_fill_gradient2(low = "dodgerblue4", mid = "white", high = "firebrick3", midpoint = 0, limits = c(-4, 4)) +  # Set colors and midpoint for the fill gradient
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

# Save the plot
ggsave("Plot_cell_proliferation.pdf", plot = Plot_cell_proliferation, width = 6, height = 10, path = "path_to/Desktop")
















## Weighted mean MFE of each cell cycle phase ####
# Compute the weighted mean MFE of the full transcriptome at each phase of the cell cycle, where the transcript-specific MFE values are weighted by the normalized counts for each transcript

# Load the transcript-level normalized counts
GSE94479_MCF7_TX_norm_counts <- read.delim("/Volumes/AC_HD_8TB/Exp.\ 277\ Analysis\ of\ Gibbs\ free\ energy\ in\ response\ to\ drug\ selection/Analysis_cell_cycle/Analysis_RNAseq_PNAS_2017/GSE94479_raw_and_normalized_counts/GSE94479_TranscriptLevel_TPM_data.csv",
                                           sep = ",", header = T, stringsAsFactors = F)

# Access the Series Matrix
library(GEOquery)
Data_GSE94479 <- getGEO(GEO = "GSE94479")

# Extract the phenotypic/clinical data matrix from the series matrix
Expdata_GSE94479 <- Data_GSE94479[["GSE94479_series_matrix.txt.gz"]]@phenoData@data
colnames(Expdata_GSE94479)[colnames(Expdata_GSE94479) == "description"] <- "condition" # Rename the column to have an easier identifier for the conditions

# Select the samples profiled with RNAseq
Expdata_GSE94479 <- dplyr::filter(Expdata_GSE94479, title %in% c("MCF7 RNA-seq G0/G1", "MCF7 RNA-seq M phase", "MCF7 RNA-seq S phase"))

Expdata_GSE94479 <- Expdata_GSE94479 %>% dplyr::mutate(condition = ifelse(condition =="RNA-seq G0/G1","G0_G1_phase", condition)) %>%
  dplyr::mutate(condition = ifelse(condition =="RNA-seq S phase","S_phase", condition)) %>%
  dplyr::mutate(condition = ifelse(condition =="RNA-seq M phase","M_phase", condition))


# Rename the columns of the normalised counts with the cell cycle phases in the Expdata_GSE94479 data frame
GSE94479_MCF7_TX_norm_counts <- dplyr::rename(GSE94479_MCF7_TX_norm_counts, G0_G1_phase = GSM2476255, S_phase = GSM2476256, M_phase = GSM2476257, ensembl_transcript_id = X)

# Join the data with the table containing the MFE values
GSE94479_MCF7_TX_norm_counts <- dplyr::left_join(GSE94479_MCF7_TX_norm_counts, Table_ANNOTATED, by="ensembl_transcript_id") %>%
  dplyr::filter(!is.na(Normalized_MFE))



# Generate a dotplot with the weighted mean MFE of each phase of the cell cycle
# Function to calculate weighted standard deviation
weighted_sd <- function(x, w) {
  sum_w <- sum(w)
  mean_w <- sum(x * w) / sum_w
  sqrt(sum(w * (x - mean_w)^2) / sum_w)
}

# Calculate weighted mean and weighted SD for each phase
G0_G1_phase_weighted_mean_MFE <- sum(GSE94479_MCF7_TX_norm_counts$Normalized_MFE * GSE94479_MCF7_TX_norm_counts$G0_G1_phase) / sum(GSE94479_MCF7_TX_norm_counts$G0_G1_phase)
G0_G1_phase_weighted_sd_MFE <- weighted_sd(GSE94479_MCF7_TX_norm_counts$Normalized_MFE, GSE94479_MCF7_TX_norm_counts$G0_G1_phase)

S_phase_weighted_mean_MFE <- sum(GSE94479_MCF7_TX_norm_counts$Normalized_MFE * GSE94479_MCF7_TX_norm_counts$S_phase) / sum(GSE94479_MCF7_TX_norm_counts$S_phase)
S_phase_weighted_sd_MFE <- weighted_sd(GSE94479_MCF7_TX_norm_counts$Normalized_MFE, GSE94479_MCF7_TX_norm_counts$S_phase)

M_phase_weighted_mean_MFE <- sum(GSE94479_MCF7_TX_norm_counts$Normalized_MFE * GSE94479_MCF7_TX_norm_counts$M_phase) / sum(GSE94479_MCF7_TX_norm_counts$M_phase)
M_phase_weighted_sd_MFE <- weighted_sd(GSE94479_MCF7_TX_norm_counts$Normalized_MFE, GSE94479_MCF7_TX_norm_counts$M_phase)

# Create a data frame including standard deviations
weighted_mean_MFE_df <- data.frame(
  Phase = c("G0/G1 Phase", "S Phase", "M Phase"),
  Weighted_MFE = c(G0_G1_phase_weighted_mean_MFE, S_phase_weighted_mean_MFE, M_phase_weighted_mean_MFE),
  SD = c(G0_G1_phase_weighted_sd_MFE, S_phase_weighted_sd_MFE, M_phase_weighted_sd_MFE)
)

# Plot with error bars
Plot_weighted_MFE <- ggplot(weighted_mean_MFE_df, aes(x = Phase, y = Weighted_MFE)) +
  geom_point(size = 6, color = "black") +  # Dots
  geom_errorbar(aes(ymin = Weighted_MFE - SD, ymax = Weighted_MFE + SD), width = 0.2, color = "black") +  # Error bars
  theme_classic() +
  ylim(-0.45, -0.22) +
  theme(
    legend.position = "none", 
    axis.title = element_text(size = 20), 
    axis.text.x = element_text(size = 16, color = "black"),
    plot.title = element_text(size = 20)
  ) +
  scale_x_discrete(limits = c("M Phase", "G0/G1 Phase", "S Phase")) +  # Specify Phase order
  labs(y = "Weighted Mean DGfold / nt (kcal/mol)", x = "", title = "Weighted Mean MFE Across Cell Cycle Phases of MCF7") 

print(Plot_weighted_MFE)











# Generate a boxplot with the transcript-level weighted MFE values for each phase
# Create a long-format data frame with transcript MFE values weighted by expression
long_df <- GSE94479_MCF7_TX_norm_counts %>%
  dplyr::mutate(G0_G1_phase = Normalized_MFE * G0_G1_phase) %>%
  dplyr::mutate(S_phase = Normalized_MFE * S_phase) %>%
  dplyr::mutate(M_phase = Normalized_MFE * M_phase) %>%
  dplyr::select(ensembl_transcript_id, Normalized_MFE, G0_G1_phase, S_phase, M_phase) %>%
  pivot_longer(cols = c(G0_G1_phase, S_phase, M_phase),
               names_to = "Phase",
               values_to = "Weighted_MFE")
  # dplyr::filter(Weighted_MFE < 0)  # Filter out negative values
  
# Plot the distribution as a box plot
library(ggbreak)

boxplot_weighted_MFE <- ggplot(long_df, aes(x = Phase, y = Weighted_MFE)) +
  geom_boxplot(fill = "lightgrey", outlier.color = "black", outlier.shape = 1) +
  theme_classic() +
  ylim(-50000, 0.25) +
  theme(
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 8, color = "black", angle = 45, hjust = 1),
    axis.text.x.top = element_blank(),   # <== removes the top x-axis labels
    axis.ticks.x.top = element_blank(),  # <== removes the top x-axis ticks
    axis.line.x.top = element_blank(),        # Remove top x-axis line
    axis.line.y.right = element_blank(),      # Remove right-side y-axis line
    axis.text.y.right = element_blank(),        # Remove right-side y-axis labels
    axis.ticks.y.right = element_blank(),    # emove right-side y-axis ticks
    plot.title = element_text(size = 20)
  ) +
  scale_x_discrete(limits = c("M_phase", "G0_G1_phase", "S_phase")) +
  labs(y = expression(Delta*G[fold]~"/ nt (kcal/mol)"), x = "", title = "") +
  scale_y_break(c(-10000, -0.75), scales = 1.5) 

#print(boxplot_weighted_MFE)

# Save the combined plot
ggsave("boxplot_weighted_MFE.pdf", plot = boxplot_weighted_MFE, width = 3, height = 5, path = "path_to/Desktop")















# Cell cycle HeLa cell line ####
# These data have been obtained from:
# 'Comparative cell cycle transcriptomics reveals synchronization of developmental transcription factor networks in cancer cells' Johan Bostrom et al. PLOS ONE 2017
# The raw and normalized counts have been downloaded from: http://www.ilincs.org/apps/grein/?gse=GSE104736

# The transcript-level normalized counts of the HeLa cell line profiled in G1, S and G2 phases are in GREIN. There are three replicates / condition

library(DESeq2)
library(GEOquery)
library(canvasXpress)
library(clinfun)
library(GGally)
library(factoextra)
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
library(data.table)




# Load the raw counts downloaded from GREIN
Counts_GSE104736_GREIN <- read.delim("/Volumes/AC_HD_8TB/Exp.\ 277\ Analysis\ of\ Gibbs\ free\ energy\ in\ response\ to\ drug\ selection/Analysis_cell_cycle/Analysis_HeLa_and_U2OS_cell_cycle/GSE104736_raw_and\ normalized_counts/GSE104736_TranscriptLevel_Raw_data.csv",
                                     sep = ",", header = T, stringsAsFactors = F)
Counts_GSE104736_GREIN <- Counts_GSE104736_GREIN %>% remove_rownames %>% column_to_rownames(var="X") # Convert first column into row names in data frame




# Access the Series Matrix
Data_GSE104736 <- getGEO(GEO = "GSE104736")
# Extract the phenotypic/clinical data matrix from the series matrix
Expdata_GSE104736 <- Data_GSE104736[["GSE104736_series_matrix.txt.gz"]]@phenoData@data

# Select the HeLa cells
Expdata_GSE104736 = Expdata_GSE104736[c(1:9),] # Remove the samples that are not necessary for the analysis
colnames(Expdata_GSE104736)[colnames(Expdata_GSE104736) == "phase:ch1"] <- "Cell_cycle_phase"
colnames(Expdata_GSE104736)[colnames(Expdata_GSE104736) == "cell line:ch1"] <- "Cell_line"
colnames(Expdata_GSE104736)[colnames(Expdata_GSE104736) == "title"] <- "condition"

# Expdata_GSE104736$condition<-gsub(", 8hr","",as.character(Expdata_GSE104736$condition)) # Remove the last characters from the description of the conditions
Expdata_GSE104736$condition <- gsub('HeLa_Fucci_G1_1', 'HeLa_G1', Expdata_GSE104736$condition) # Edit the name of the description in the conditions
Expdata_GSE104736$condition <- gsub('HeLa_Fucci_G1_2', 'HeLa_G1', Expdata_GSE104736$condition) # Edit the name of the description in the conditions
Expdata_GSE104736$condition <- gsub('HeLa_Fucci_G1_3', 'HeLa_G1', Expdata_GSE104736$condition) # Edit the name of the description in the conditions
Expdata_GSE104736$condition <- gsub('HeLa_Fucci_S_1', 'HeLa_S', Expdata_GSE104736$condition) # Edit the name of the description in the conditions
Expdata_GSE104736$condition <- gsub('HeLa_Fucci_S_2', 'HeLa_S', Expdata_GSE104736$condition) # Edit the name of the description in the conditions
Expdata_GSE104736$condition <- gsub('HeLa_Fucci_S_3', 'HeLa_S', Expdata_GSE104736$condition) # Edit the name of the description in the conditions
Expdata_GSE104736$condition <- gsub('HeLa_Fucci_G2_1', 'HeLa_G2', Expdata_GSE104736$condition) # Edit the name of the description in the conditions
Expdata_GSE104736$condition <- gsub('HeLa_Fucci_G2_2', 'HeLa_G2', Expdata_GSE104736$condition) # Edit the name of the description in the conditions
Expdata_GSE104736$condition <- gsub('HeLa_Fucci_G2_3', 'HeLa_G2', Expdata_GSE104736$condition) # Edit the name of the description in the conditions



# Analysis HeLa cells
# Data are converted and organized to the format required
Counts_GSE104736_GREIN <- as.matrix(Counts_GSE104736_GREIN)
Counts_GSE104736_GREIN <- Counts_GSE104736_GREIN[,c(1:9)] # Remove the samples that are not necessary for the analysis

# Round counts down to the nearest integer
Counts_GSE104736_GREIN <- floor(Counts_GSE104736_GREIN)

# IMPORTANT:
# It is absolutely critical that the columns of the count matrix and the rows of the column data (information about samples) are in the same order.

# Replace the row names of Expdata_GSE104736 (sampleID) with the same sample name [title] of [Counts_GSE104736]. 
# This will help to match sample names in both matrices.
# The outcome should be TRUE
all(rownames(Expdata_GSE104736) %in% colnames(Counts_GSE104736_GREIN))
all(colnames(Counts_GSE104736_GREIN) %in% rownames(Expdata_GSE104736))

# With the count matrix, Counts_GSE104736_GREIN, and the sample information, Expdata_GSE104736, we can construct a DESeqDataSet:
dds_GSE104736 <- DESeqDataSetFromMatrix(countData = Counts_GSE104736_GREIN,
                                        colData = Expdata_GSE104736,
                                        design = ~ condition)
dds_GSE104736



# DATA QUALITY ASSESSMENT
# Before performing differential expression analysis, we perform data quality assessment by sample clustering and visualization
library("pheatmap")
# We log transform the counts
ntd <- normTransform(dds_GSE104736)
# Heatmap of the count matrix
select <- order(rowMeans(counts(dds_GSE104736,normalized=FALSE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds_GSE104736)[,c("condition","type")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)



# Heatmap of the sample-to-sample distances
# A heatmap of this distance matrix gives us an overview over similarities and dissimilarities between samples
# Transform the values
vsd <- vst(dds_GSE104736, blind=FALSE)
# rld <- rlog(dds_GSE104736, blind=FALSE)
head(assay(vsd), 3)

# Apply the dist function to the transpose of the transformed count matrix to get sample-to-sample distances
sampleDists <- dist(t(assay(vsd)))

library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

# Principal component plot of the samples
plotPCA(vsd, intgroup="condition")


# DIFFERENTIAL EXPRESSION ANALYSIS
# We perform a minimal pre-filtering to keep only rows that have at least 3 reads total
keep <- rowSums(counts(dds_GSE104736)) >= 3
dds_GSE104736 <- dds_GSE104736[keep,]
table(keep) # This shows how many transcripts have been removed after the filtering
nrow(dds_GSE104736) # to confirm the filtered transcripts

# DIFFERENTIAL EXPRESSION ANALYSIS
dds_GSE104736 <- DESeq(dds_GSE104736) # The DESeq() function normalizes the read counts, estimates dispersions, and fits the linear model, all in one go
# dds_GSE104736 <- DESeq(dds_GSE104736, betaPrior = T)  The betaPrior = T opton tells DESeq2 to squeeze the log Fold Changes of lowly expressed genes toward zero.








## Differential gene expression G1 vs G2 ####

# Specify the conditions to contrast for differential expression analysis
res_GSE104736_G1_vs_G2 <- results(dds_GSE104736, contrast=c("condition","HeLa_G1","HeLa_G2")) # to specify the samples to contrast
res_GSE104736_G1_vs_G2
# Get a summary of the results
summary(res_GSE104736_G1_vs_G2)



# Export results to CSV files
# write.csv(as.data.frame(res), file="Treated_vs_Controls.csv")

# Plot the log2 Fold Changes as function of expresson:
plotMA(res_GSE104736_G1_vs_G2,alpha = 0.05)
abline(h=c(-2,2),col = "red")

# Transform the row names into a column
res_GSE104736_df <- data.frame(res_GSE104736_G1_vs_G2)
res_GSE104736_G1_vs_G2 <- tibble::rownames_to_column(res_GSE104736_df, "ensembl_transcript_id")



# Join the expression change values with the Table_ANNOTATED containing MFE values
MFE_GSE104736_G1_vs_G2 <- dplyr::inner_join(Table_ANNOTATED, res_GSE104736_G1_vs_G2, by="ensembl_transcript_id") %>%
  dplyr::filter(padj < 0.05) 


# Generate a scatter plot Log2_FC vs Normalized_MFE
library(ggplot2)
Plot_GSE104736 <- ggplot(MFE_GSE104736_G1_vs_G2, aes(y=Normalized_MFE, x=log2FoldChange)) +
  geom_point(size=1.5, alpha=0.2)+
  theme_classic() +
  theme(legend.position="none") +
  geom_vline(xintercept = 0, size=0.3) +
  #xlim(-5.5,5.5) +
  theme(axis.title = element_text(size = 20)) +
  theme(plot.title = element_text(size = 15)) +
  labs(y="DGfold / nt (kcal/mol)",x="Log2_FC_expression",
       title="GSE104736_HeLa_Log2_FC_adj.P<0.05_G1_vs_G2")
print(Plot_GSE104736)



# Look at the mean log2_FC
MFE_GSE104736_G1_vs_G2_UP <- dplyr::filter(MFE_GSE104736_G1_vs_G2, log2FoldChange > 0)
Mean_MFE_GSE104736_G1_vs_G2_UP <- mean(MFE_GSE104736_G1_vs_G2_UP$Normalized_MFE)

MFE_GSE104736_G1_vs_G2_DOWN <- dplyr::filter(MFE_GSE104736_G1_vs_G2, log2FoldChange < 0)
Mean_MFE_GSE104736_G1_vs_G2_DOWN <- mean(MFE_GSE104736_G1_vs_G2_DOWN$Normalized_MFE)

TTEST_MFE_GSE104736_G1_vs_G2_all_UP_vs_DOWN <- t.test(MFE_GSE104736_G1_vs_G2_DOWN$Normalized_MFE, MFE_GSE104736_G1_vs_G2_UP$Normalized_MFE, paired=FALSE)
print(TTEST_MFE_GSE104736_G1_vs_G2_all_UP_vs_DOWN)



# Add K-S test to look at the differences between the up-regukated and down-regulated genes
ks_test_GSE104736_all_UP_vs_DOWN <- ks.test(MFE_GSE104736_G1_vs_G2_DOWN$Normalized_MFE, MFE_GSE104736_G1_vs_G2_UP$Normalized_MFE)
print(ks_test_GSE104736_all_UP_vs_DOWN)


quantile_075_MFE_GSE104736_G1_vs_G2_UP <- quantile(MFE_GSE104736_G1_vs_G2_UP$Normalized_MFE, probs = c(0.75)) 
quantile_075_MFE_GSE104736_G1_vs_G2_DOWN <- quantile(MFE_GSE104736_G1_vs_G2_DOWN$Normalized_MFE, probs = c(0.75)) 
# Values above 0.75 quantile
ks_test_MFE_GSE104736_G1_vs_G2_075_above <- ks.test(MFE_GSE104736_G1_vs_G2_UP$Normalized_MFE[MFE_GSE104736_G1_vs_G2_UP$Normalized_MFE >= quantile_075_MFE_GSE104736_G1_vs_G2_UP],
                                           MFE_GSE104736_G1_vs_G2_DOWN$Normalized_MFE[MFE_GSE104736_G1_vs_G2_DOWN$Normalized_MFE >= quantile_075_MFE_GSE104736_G1_vs_G2_DOWN])
print(ks_test_MFE_GSE104736_G1_vs_G2_075_above)

# Look at the differences between the 0.25 quantile 
quantile_025_MFE_GSE104736_G1_vs_G2_UP <- quantile(MFE_GSE104736_G1_vs_G2_UP$Normalized_MFE, probs = c(0.25)) 
quantile_025_MFE_GSE104736_G1_vs_G2_DOWN <- quantile(MFE_GSE104736_G1_vs_G2_DOWN$Normalized_MFE, probs = c(0.25)) 
# Values below 0.25 quantile
ks_test_MF_GSE104736_025_below <- ks.test(MFE_GSE104736_G1_vs_G2_UP$Normalized_MFE[MFE_GSE104736_G1_vs_G2_UP$Normalized_MFE <= quantile_025_MFE_GSE104736_G1_vs_G2_UP],
                                          MFE_GSE104736_G1_vs_G2_DOWN$Normalized_MFE[MFE_GSE104736_G1_vs_G2_DOWN$Normalized_MFE <= quantile_025_MFE_GSE104736_G1_vs_G2_DOWN])
print(ks_test_MF_GSE104736_025_below)




# Violin plots with box plots inside
GSE104736_HeLa_G1_vs_G2_plot <- ggplot(MFE_GSE104736_G1_vs_G2, aes(x = factor(log2FoldChange > 0), y = Normalized_MFE)) +
  geom_violin(fill = "lightgrey") +  # Set violin plot fill to light grey
  geom_boxplot(width = 0.1, position = position_dodge(0.9), color = "black", fill = "white", outlier.shape = NA) +  # Set box plot fill to white
  theme_classic() +
  theme(
    legend.position = "none", 
    axis.title = element_text(size = 20), 
    plot.title = element_text(size = 10),
    axis.text.x = element_text(size = 16, color = "black")
  ) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  ylim(-0.6, 0) +
  labs(
    y = "Normalized_MFE", 
    x = NULL, 
    title = "GSE104736_Log2_FC_adj.P<0.05_G1_vs_G2"
  ) +
  scale_x_discrete(labels = c("log2FoldChange < 0", "log2FoldChange > 0"))

print(GSE104736_HeLa_G1_vs_G2_plot)

# Save the combined plot
ggsave("GSE104736_HeLa_G1_vs_G2_plot.pdf", plot = GSE104736_HeLa_G1_vs_G2_plot, width = 5, height = 6, path = "path_to/Desktop")






# Differential gene expression G1 vs G2 heat capacity changes

# This is calculated on the total n. of transcripts with a significant log2 fold change
# Generate a scatter plot Log2_FC vs heat capacity
ggplot(MFE_GSE104736_G1_vs_G2, aes(y=Cp_37C, x=log2FoldChange)) +
  geom_point(size=1.5, alpha=0.1)+
  theme_classic() +
  theme(legend.position="none") +
  geom_vline(xintercept = 0, size=0.3) +
  xlim(-12, 12) +
  ylim(0,400) +
  theme(axis.title = element_text(size = 20)) +
  theme(plot.title = element_text(size = 15)) +
  labs(y="Heat capacity (J/C)",x="Log2_FC_expression",title="MFE_GSE104736_G1_vs_G2")


# Show the data with a box plot
Plot_MFE_GSE104736_G1_vs_G2_heat_capacity <- ggplot(MFE_GSE104736_G1_vs_G2, aes(x = factor(log2FoldChange > 0), y = Cp_37C)) +
  geom_boxplot(fill = "lightgrey", alpha = 0.7) +
  #geom_jitter(width = 0.1, size = 1, alpha = 0.3) +
  theme_classic() +
  ylim(0, 400) +
  theme(legend.position = "none", 
        axis.title = element_text(size = 20), 
        plot.title = element_text(size = 10),
        axis.text.x = element_text(size = 16, color = "black")) +
  theme(plot.title = element_text(size = 20)) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  labs(y="Heat capacity (J/C)",x="",title="MFE_GSE104736_G1_vs_G2") +
  scale_x_discrete(labels = c("log2FoldChange > 0", "log2FoldChange < 0"))

print(Plot_MFE_GSE104736_G1_vs_G2_heat_capacity)

ggsave("Plot_MFE_GSE104736_G1_vs_G2_heat_capacity.pdf", plot = Plot_MFE_GSE104736_G1_vs_G2_heat_capacity, width = 5, height = 6, path = "path_to/Desktop")


Mean_heat_capacity_MFE_GSE104736_G1_vs_G2_UP <- mean(MFE_GSE104736_G1_vs_G2_UP$Cp_37C, na.rm = TRUE)
Mean_heat_capacity_MFE_GSE104736_G1_vs_G2_DOWN <- mean(MFE_GSE104736_G1_vs_G2_DOWN$Cp_37C, na.rm = TRUE)
TTEST_heat_capacity_MFE_GSE104736_G1_vs_G2_all_UP_vs_DOWN <- t.test(MFE_GSE104736_G1_vs_G2_DOWN$Cp_37C, MFE_GSE104736_G1_vs_G2_UP$Cp_37C, paired=FALSE)
print(TTEST_heat_capacity_MFE_GSE104736_G1_vs_G2_all_UP_vs_DOWN)
# Add K-S test to look at the differences between the up-regukated and down-regulated genes
ks_test_MFE_GSE104736_G1_vs_G2_all_UP_vs_DOWN_heat_capacity <- ks.test(MFE_GSE104736_G1_vs_G2_DOWN$Cp_37C, MFE_GSE104736_G1_vs_G2_UP$Cp_37C)
print(ks_test_MFE_GSE104736_G1_vs_G2_all_UP_vs_DOWN_heat_capacity)



















## Differential gene expression S vs G1 ####

# Specify the conditions to contrast for differential expression analysis
res_GSE104736_S_vs_G1 <- results(dds_GSE104736, contrast=c("condition","HeLa_S","HeLa_G1")) # to specify the samples to contrast
res_GSE104736_S_vs_G1
# Get a summary of the results
summary(res_GSE104736_S_vs_G1)



# Export results to CSV files
# write.csv(as.data.frame(res), file="Treated_vs_Controls.csv")

# Plot the log2 Fold Changes as function of expresson:
plotMA(res_GSE104736_S_vs_G1,alpha = 0.05)
abline(h=c(-2,2),col = "red")

# Transform the row names into a column
res_GSE104736_df <- data.frame(res_GSE104736_S_vs_G1)
res_GSE104736_S_vs_G1 <- tibble::rownames_to_column(res_GSE104736_df, "ensembl_transcript_id")



# Join the expression change values with the Table_ANNOTATED containing MFE values
MFE_GSE104736_S_vs_G1 <- dplyr::inner_join(Table_ANNOTATED, res_GSE104736_S_vs_G1, by="ensembl_transcript_id") %>%
  dplyr::filter(padj < 0.05) %>%
  dplyr::filter(log2FoldChange < 15 & log2FoldChange > -15) # Remove the few transcripts with extreme log2FoldChange values


# Generate a scatter plot Log2_FC vs Normalized_MFE
library(ggplot2)
Plot_GSE104736 <- ggplot(MFE_GSE104736_S_vs_G1, aes(y=Normalized_MFE, x=log2FoldChange)) +
  geom_point(size=1.5, alpha=0.2)+
  theme_classic() +
  theme(legend.position="none") +
  geom_vline(xintercept = 0, size=0.3) +
  #xlim(-5.5,5.5) +
  theme(axis.title = element_text(size = 20)) +
  theme(plot.title = element_text(size = 15)) +
  labs(y="DGfold / nt (kcal/mol)",x="Log2_FC_expression",
       title="GSE104736_Log2_FC_adj.P<0.05_S_vs_G1")
print(Plot_GSE104736)

# Save plot
# ggsave("Plot_GSE104736.pdf", plot = Plot_GSE104736, path = "path_to/Desktop")


# Look at the mean log2_FC
MFE_GSE104736_S_vs_G1_UP <- dplyr::filter(MFE_GSE104736_S_vs_G1, log2FoldChange > 0)
Mean_MFE_GSE104736_S_vs_G1_UP <- mean(MFE_GSE104736_S_vs_G1_UP$Normalized_MFE)

MFE_GSE104736_S_vs_G1_DOWN <- dplyr::filter(MFE_GSE104736_S_vs_G1, log2FoldChange < 0)
Mean_MFE_GSE104736_S_vs_G1_DOWN <- mean(MFE_GSE104736_S_vs_G1_DOWN$Normalized_MFE)

TTEST_MFE_GSE104736_S_vs_G1_all_UP_vs_DOWN <- t.test(MFE_GSE104736_S_vs_G1_DOWN$Normalized_MFE, MFE_GSE104736_S_vs_G1_UP$Normalized_MFE, paired=FALSE)
print(TTEST_MFE_GSE104736_S_vs_G1_all_UP_vs_DOWN)



# Add K-S test to look at the differences between the up-regukated and down-regulated genes
ks_test_GSE104736_all_UP_vs_DOWN <- ks.test(MFE_GSE104736_S_vs_G1_DOWN$Normalized_MFE, MFE_GSE104736_S_vs_G1_UP$Normalized_MFE)
print(ks_test_GSE104736_all_UP_vs_DOWN)


quantile_075_MFE_GSE104736_S_vs_G1_UP <- quantile(MFE_GSE104736_S_vs_G1_UP$Normalized_MFE, probs = c(0.75)) 
quantile_075_MFE_GSE104736_S_vs_G1_DOWN <- quantile(MFE_GSE104736_S_vs_G1_DOWN$Normalized_MFE, probs = c(0.75)) 
# Values above 0.75 quantile
ks_test_MFE_GSE104736_S_vs_G1_075_above <- ks.test(MFE_GSE104736_S_vs_G1_UP$Normalized_MFE[MFE_GSE104736_S_vs_G1_UP$Normalized_MFE >= quantile_075_MFE_GSE104736_S_vs_G1_UP],
                                           MFE_GSE104736_S_vs_G1_DOWN$Normalized_MFE[MFE_GSE104736_S_vs_G1_DOWN$Normalized_MFE >= quantile_075_MFE_GSE104736_S_vs_G1_DOWN])
print(ks_test_MFE_GSE104736_S_vs_G1_075_above)

# Look at the differences between the 0.25 quantile 
quantile_025_MFE_GSE104736_S_vs_G1_UP <- quantile(MFE_GSE104736_S_vs_G1_UP$Normalized_MFE, probs = c(0.25)) 
quantile_025_MFE_GSE104736_S_vs_G1_DOWN <- quantile(MFE_GSE104736_S_vs_G1_DOWN$Normalized_MFE, probs = c(0.25)) 
# Values below 0.25 quantile
ks_test_MF_GSE104736_025_below <- ks.test(MFE_GSE104736_S_vs_G1_UP$Normalized_MFE[MFE_GSE104736_S_vs_G1_UP$Normalized_MFE <= quantile_025_MFE_GSE104736_S_vs_G1_UP],
                                          MFE_GSE104736_S_vs_G1_DOWN$Normalized_MFE[MFE_GSE104736_S_vs_G1_DOWN$Normalized_MFE <= quantile_025_MFE_GSE104736_S_vs_G1_DOWN])
print(ks_test_MF_GSE104736_025_below)




# Violin plots with box plots inside
GSE104736_HeLa_S_vs_G1_plot <- ggplot(MFE_GSE104736_S_vs_G1, aes(x = factor(log2FoldChange > 0), y = Normalized_MFE)) +
  geom_violin(fill = "lightgrey") +  # Set violin plot fill to light grey
  geom_boxplot(width = 0.1, position = position_dodge(0.9), color = "black", fill = "white", outlier.shape = NA) +  # Set box plot fill to white
  theme_classic() +
  theme(
    legend.position = "none", 
    axis.title = element_text(size = 20), 
    plot.title = element_text(size = 10),
    axis.text.x = element_text(size = 16, color = "black")
  ) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  ylim(-0.6, 0) +
  labs(
    y = "Normalized_MFE", 
    x = NULL, 
    title = "GSE104736_Log2_FC_adj.P<0.05_S_vs_G1"
  ) +
  scale_x_discrete(labels = c("log2FoldChange < 0", "log2FoldChange > 0"))

print(GSE104736_HeLa_S_vs_G1_plot)

# Save the combined plot
ggsave("GSE104736_HeLa_S_vs_G1_plot.pdf", plot = GSE104736_HeLa_S_vs_G1_plot, width = 5, height = 6, path = "path_to/Desktop")






# Differential gene expression S vs G1 heat capacity changes

# This is calculated on the total n. of transcripts with a significant log2 fold change
# Generate a scatter plot Log2_FC vs heat capacity
ggplot(MFE_GSE104736_S_vs_G1, aes(y=Cp_37C, x=log2FoldChange)) +
  geom_point(size=1.5, alpha=0.1)+
  theme_classic() +
  theme(legend.position="none") +
  geom_vline(xintercept = 0, size=0.3) +
  xlim(-12, 12) +
  ylim(0,400) +
  theme(axis.title = element_text(size = 20)) +
  theme(plot.title = element_text(size = 15)) +
  labs(y="Heat capacity (J/C)",x="Log2_FC_expression",title="MFE_GSE104736_S_vs_G1")


# Show the data with a box plot
Plot_MFE_GSE104736_S_vs_G1_heat_capacity <- ggplot(MFE_GSE104736_S_vs_G1, aes(x = factor(log2FoldChange > 0), y = Cp_37C)) +
  geom_boxplot(fill = "lightgrey", alpha = 0.7) +
  #geom_jitter(width = 0.1, size = 1, alpha = 0.3) +
  theme_classic() +
  ylim(0, 400) +
  theme(legend.position = "none", 
        axis.title = element_text(size = 20), 
        plot.title = element_text(size = 10),
        axis.text.x = element_text(size = 16, color = "black")) +
  theme(plot.title = element_text(size = 20)) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  labs(y="Heat capacity (J/C)",x="",title="MFE_GSE104736_S_vs_G1") +
  scale_x_discrete(labels = c("log2FoldChange > 0", "log2FoldChange < 0"))

print(Plot_MFE_GSE104736_S_vs_G1_heat_capacity)

ggsave("Plot_MFE_GSE104736_S_vs_G1_heat_capacity.pdf", plot = Plot_MFE_GSE104736_S_vs_G1_heat_capacity, width = 5, height = 6, path = "path_to/Desktop")


Mean_heat_capacity_MFE_GSE104736_S_vs_G1_UP <- mean(MFE_GSE104736_S_vs_G1_UP$Cp_37C, na.rm = TRUE)
Mean_heat_capacity_MFE_GSE104736_S_vs_G1_DOWN <- mean(MFE_GSE104736_S_vs_G1_DOWN$Cp_37C, na.rm = TRUE)
TTEST_heat_capacity_MFE_GSE104736_S_vs_G1_all_UP_vs_DOWN <- t.test(MFE_GSE104736_S_vs_G1_DOWN$Cp_37C, MFE_GSE104736_S_vs_G1_UP$Cp_37C, paired=FALSE)
print(TTEST_heat_capacity_MFE_GSE104736_S_vs_G1_all_UP_vs_DOWN)
# Add K-S test to look at the differences between the up-regukated and down-regulated genes
ks_test_MFE_GSE104736_S_vs_G1_all_UP_vs_DOWN_heat_capacity <- ks.test(MFE_GSE104736_S_vs_G1_DOWN$Cp_37C, MFE_GSE104736_S_vs_G1_UP$Cp_37C)
print(ks_test_MFE_GSE104736_S_vs_G1_all_UP_vs_DOWN_heat_capacity)























## Differential gene expression G2 vs S ####

# Specify the conditions to contrast for differential expression analysis
res_GSE104736_G2_vs_S <- results(dds_GSE104736, contrast=c("condition","HeLa_G2","HeLa_S")) # to specify the samples to contrast
res_GSE104736_G2_vs_S
# Get a summary of the results
summary(res_GSE104736_G2_vs_S)

# Export results to CSV files
# write.csv(as.data.frame(res), file="Treated_vs_Controls.csv")

# Plot the log2 Fold Changes as function of expresson:
plotMA(res_GSE104736_G2_vs_S,alpha = 0.05)
abline(h=c(-2,2),col = "red")

# Transform the row names into a column
res_GSE104736_df <- data.frame(res_GSE104736_G2_vs_S)
res_GSE104736_G2_vs_S <- tibble::rownames_to_column(res_GSE104736_df, "ensembl_transcript_id")



# Join the expression change values with the Table_ANNOTATED containing MFE values
MFE_GSE104736_G2_vs_S <- dplyr::inner_join(Table_ANNOTATED, res_GSE104736_G2_vs_S, by="ensembl_transcript_id") %>%
  dplyr::filter(padj < 0.05) %>%
  dplyr::filter(log2FoldChange < 15 & log2FoldChange > -15) # Remove the few transcripts with extreme log2FoldChange values



# Generate a scatter plot Log2_FC vs Normalized_MFE
library(ggplot2)
Plot_GSE104736 <- ggplot(MFE_GSE104736_G2_vs_S, aes(y=Normalized_MFE, x=log2FoldChange)) +
  geom_point(size=1.5, alpha=0.2)+
  theme_classic() +
  theme(legend.position="none") +
  geom_vline(xintercept = 0, size=0.3) +
  #xlim(-5.5,5.5) +
  theme(axis.title = element_text(size = 20)) +
  theme(plot.title = element_text(size = 15)) +
  labs(y="DGfold / nt (kcal/mol)",x="Log2_FC_expression",
       title="GSE104736_Log2_FC_adj.P<0.05_G2_vs_S")
print(Plot_GSE104736)

# Save plot
# ggsave("Plot_GSE104736.pdf", plot = Plot_GSE104736, path = "path_to/Desktop")


# Look at the mean log2_FC
MFE_GSE104736_G2_vs_S_UP <- dplyr::filter(MFE_GSE104736_G2_vs_S, log2FoldChange > 0)
Mean_MFE_GSE104736_G2_vs_S_UP <- mean(MFE_GSE104736_G2_vs_S_UP$Normalized_MFE)

MFE_GSE104736_G2_vs_S_DOWN <- dplyr::filter(MFE_GSE104736_G2_vs_S, log2FoldChange < 0)
Mean_MFE_GSE104736_G2_vs_S_DOWN <- mean(MFE_GSE104736_G2_vs_S_DOWN$Normalized_MFE)

TTEST_MFE_GSE104736_G2_vs_S_all_UP_vs_DOWN <- t.test(MFE_GSE104736_G2_vs_S_DOWN$Normalized_MFE, MFE_GSE104736_G2_vs_S_UP$Normalized_MFE, paired=FALSE)
print(TTEST_MFE_GSE104736_G2_vs_S_all_UP_vs_DOWN)



# Add K-S test to look at the differences between the up-regukated and down-regulated genes
ks_test_GSE104736_all_UP_vs_DOWN <- ks.test(MFE_GSE104736_G2_vs_S_DOWN$Normalized_MFE, MFE_GSE104736_G2_vs_S_UP$Normalized_MFE)
print(ks_test_GSE104736_all_UP_vs_DOWN)


quantile_075_MFE_GSE104736_G2_vs_S_UP <- quantile(MFE_GSE104736_G2_vs_S_UP$Normalized_MFE, probs = c(0.75)) 
quantile_075_MFE_GSE104736_G2_vs_S_DOWN <- quantile(MFE_GSE104736_G2_vs_S_DOWN$Normalized_MFE, probs = c(0.75)) 
# Values above 0.75 quantile
ks_test_MFE_GSE104736_G2_vs_S_075_above <- ks.test(MFE_GSE104736_G2_vs_S_UP$Normalized_MFE[MFE_GSE104736_G2_vs_S_UP$Normalized_MFE >= quantile_075_MFE_GSE104736_G2_vs_S_UP],
                                           MFE_GSE104736_G2_vs_S_DOWN$Normalized_MFE[MFE_GSE104736_G2_vs_S_DOWN$Normalized_MFE >= quantile_075_MFE_GSE104736_G2_vs_S_DOWN])
print(ks_test_MFE_GSE104736_G2_vs_S_075_above)

# Look at the differences between the 0.25 quantile 
quantile_025_MFE_GSE104736_G2_vs_S_UP <- quantile(MFE_GSE104736_G2_vs_S_UP$Normalized_MFE, probs = c(0.25)) 
quantile_025_MFE_GSE104736_G2_vs_S_DOWN <- quantile(MFE_GSE104736_G2_vs_S_DOWN$Normalized_MFE, probs = c(0.25)) 
# Values below 0.25 quantile
ks_test_MF_GSE104736_025_below <- ks.test(MFE_GSE104736_G2_vs_S_UP$Normalized_MFE[MFE_GSE104736_G2_vs_S_UP$Normalized_MFE <= quantile_025_MFE_GSE104736_G2_vs_S_UP],
                                          MFE_GSE104736_G2_vs_S_DOWN$Normalized_MFE[MFE_GSE104736_G2_vs_S_DOWN$Normalized_MFE <= quantile_025_MFE_GSE104736_G2_vs_S_DOWN])
print(ks_test_MF_GSE104736_025_below)




# Violin plots with box plots inside
GSE104736_HeLa_G2_vs_S_plot <- ggplot(MFE_GSE104736_G2_vs_S, aes(x = factor(log2FoldChange > 0), y = Normalized_MFE)) +
  geom_violin(fill = "lightgrey") +  # Set violin plot fill to light grey
  geom_boxplot(width = 0.1, position = position_dodge(0.9), color = "black", fill = "white", outlier.shape = NA) +  # Set box plot fill to white
  theme_classic() +
  theme(
    legend.position = "none", 
    axis.title = element_text(size = 20), 
    plot.title = element_text(size = 10),
    axis.text.x = element_text(size = 16, color = "black")
  ) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  ylim(-0.6, 0) +
  labs(
    y = "Normalized_MFE", 
    x = NULL, 
    title = "GSE104736_Log2_FC_adj.P<0.05_G2_vs_S"
  ) +
  scale_x_discrete(labels = c("log2FoldChange < 0", "log2FoldChange > 0"))

print(GSE104736_HeLa_G2_vs_S_plot)

# Save the combined plot
ggsave("GSE104736_HeLa_G2_vs_S_plot.pdf", plot = GSE104736_HeLa_G2_vs_S_plot, width = 5, height = 6, path = "path_to/Desktop")








# Differential gene expression G2 vs S heat capacity changes

# This is calculated on the total n. of transcripts with a significant log2 fold change
# Generate a scatter plot Log2_FC vs heat capacity
ggplot(MFE_GSE104736_G2_vs_S, aes(y=Cp_37C, x=log2FoldChange)) +
  geom_point(size=1.5, alpha=0.1)+
  theme_classic() +
  theme(legend.position="none") +
  geom_vline(xintercept = 0, size=0.3) +
  xlim(-12, 12) +
  ylim(0,400) +
  theme(axis.title = element_text(size = 20)) +
  theme(plot.title = element_text(size = 15)) +
  labs(y="Heat capacity (J/C)",x="Log2_FC_expression",title="MFE_GSE104736_G2_vs_S")


# Show the data with a box plot
Plot_MFE_GSE104736_G2_vs_S_heat_capacity <- ggplot(MFE_GSE104736_G2_vs_S, aes(x = factor(log2FoldChange > 0), y = Cp_37C)) +
  geom_boxplot(fill = "lightgrey", alpha = 0.7) +
  #geom_jitter(width = 0.1, size = 1, alpha = 0.3) +
  theme_classic() +
  ylim(0, 400) +
  theme(legend.position = "none", 
        axis.title = element_text(size = 20), 
        plot.title = element_text(size = 10),
        axis.text.x = element_text(size = 16, color = "black")) +
  theme(plot.title = element_text(size = 20)) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  labs(y="Heat capacity (J/C)",x="",title="MFE_GSE104736_G2_vs_S") +
  scale_x_discrete(labels = c("log2FoldChange > 0", "log2FoldChange < 0"))

print(Plot_MFE_GSE104736_G2_vs_S_heat_capacity)

ggsave("Plot_MFE_GSE104736_G2_vs_S_heat_capacity.pdf", plot = Plot_MFE_GSE104736_G2_vs_S_heat_capacity, width = 5, height = 6, path = "path_to/Desktop")


Mean_heat_capacity_MFE_GSE104736_G2_vs_S_UP <- mean(MFE_GSE104736_G2_vs_S_UP$Cp_37C, na.rm = TRUE)
Mean_heat_capacity_MFE_GSE104736_G2_vs_S_DOWN <- mean(MFE_GSE104736_G2_vs_S_DOWN$Cp_37C, na.rm = TRUE)
TTEST_heat_capacity_MFE_GSE104736_G2_vs_S_all_UP_vs_DOWN <- t.test(MFE_GSE104736_G2_vs_S_DOWN$Cp_37C, MFE_GSE104736_G2_vs_S_UP$Cp_37C, paired=FALSE)
print(TTEST_heat_capacity_MFE_GSE104736_G2_vs_S_all_UP_vs_DOWN)
# Add K-S test to look at the differences between the up-regukated and down-regulated genes
ks_test_MFE_GSE104736_G2_vs_S_all_UP_vs_DOWN_heat_capacity <- ks.test(MFE_GSE104736_G2_vs_S_DOWN$Cp_37C, MFE_GSE104736_G2_vs_S_UP$Cp_37C)
print(ks_test_MFE_GSE104736_G2_vs_S_all_UP_vs_DOWN_heat_capacity)
















## Weighted mean MFE of each cell cycle phase ####
# Compute the weighted mean MFE of the full transcriptome at each phase of the cell cycle, where the transcript-specific MFE values are weighted by the normalized counts for each transcript

# Load the transcript-level normalized counts
GSE104736_HeLa_TX_norm_counts <- read.delim("/Volumes/AC_HD_8TB/Exp.\ 277\ Analysis\ of\ Gibbs\ free\ energy\ in\ response\ to\ drug\ selection/Analysis_cell_cycle/Analysis_HeLa_and_U2OS_cell_cycle/GSE104736_raw_and\ normalized_counts/GSE104736_TranscriptLevel_TPM_data.csv",
                                           sep = ",", header = T, stringsAsFactors = F)


# Based on the Expdata_GSE104736 metadata file, select the samples corresponding to the different phases of the cell cycle
GSE104736_HeLa_TX_norm_counts <- dplyr::select(GSE104736_HeLa_TX_norm_counts,
                                               X,
                                               GSM2806897,
                                               GSM2806898,
                                               GSM2806899,
                                               GSM2806900,
                                               GSM2806901,
                                               GSM2806902,
                                               GSM2806903,
                                               GSM2806904,
                                               GSM2806905) %>%
  dplyr::rename(ensembl_transcript_id = X)

# Samples GSM2806897, GSM2806898, GSM2806899 are the triplicates of the G1 phase
# Samples GSM2806900, GSM2806901, GSM2806902 are the triplicates of the S phase
# Samples GSM2806903, GSM2806904, GSM2806905 are the triplicates of the G2 phase

# Average the counts across the triplicates for each phase of the cell cycle
GSE104736_HeLa_TX_norm_counts_averaged <- GSE104736_HeLa_TX_norm_counts %>%
  dplyr::mutate(
    G1_phase = rowMeans(dplyr::select(., GSM2806897, GSM2806898, GSM2806899), na.rm = TRUE),
    S_phase  = rowMeans(dplyr::select(., GSM2806900, GSM2806901, GSM2806902), na.rm = TRUE),
    G2_phase = rowMeans(dplyr::select(., GSM2806903, GSM2806904, GSM2806905), na.rm = TRUE)
  ) %>%
  dplyr::select(ensembl_transcript_id, G1_phase, S_phase, G2_phase)  # Keep only the new mean columns


# Join the data with the table containing the MFE values
GSE104736_HeLa_TX_norm_counts <- dplyr::left_join(GSE104736_HeLa_TX_norm_counts_averaged, Table_ANNOTATED, by="ensembl_transcript_id") %>%
  dplyr::filter(!is.na(Normalized_MFE))










# Generate a dotplot with the weighted mean MFE of each phase of the cell cycle
# Function to calculate weighted standard deviation
weighted_sd <- function(x, w) {
  sum_w <- sum(w)
  mean_w <- sum(x * w) / sum_w
  sqrt(sum(w * (x - mean_w)^2) / sum_w)
}

# Calculate weighted mean and weighted SD for each phase
G1_phase_weighted_mean_MFE <- sum(GSE104736_HeLa_TX_norm_counts$Normalized_MFE * GSE104736_HeLa_TX_norm_counts$G1_phase) / sum(GSE104736_HeLa_TX_norm_counts$G1_phase)
G1_phase_weighted_sd_MFE <- weighted_sd(GSE104736_HeLa_TX_norm_counts$Normalized_MFE, GSE104736_HeLa_TX_norm_counts$G1_phase)

S_phase_weighted_mean_MFE <- sum(GSE104736_HeLa_TX_norm_counts$Normalized_MFE * GSE104736_HeLa_TX_norm_counts$S_phase) / sum(GSE104736_HeLa_TX_norm_counts$S_phase)
S_phase_weighted_sd_MFE <- weighted_sd(GSE104736_HeLa_TX_norm_counts$Normalized_MFE, GSE104736_HeLa_TX_norm_counts$S_phase)

G2_phase_weighted_mean_MFE <- sum(GSE104736_HeLa_TX_norm_counts$Normalized_MFE * GSE104736_HeLa_TX_norm_counts$G2_phase) / sum(GSE104736_HeLa_TX_norm_counts$G2_phase)
G2_phase_weighted_sd_MFE <- weighted_sd(GSE104736_HeLa_TX_norm_counts$Normalized_MFE, GSE104736_HeLa_TX_norm_counts$G2_phase)

# Create a data frame including standard deviations
weighted_mean_MFE_df <- data.frame(
  Phase = c("G1 Phase", "S Phase", "G2 Phase"),
  Weighted_MFE = c(G1_phase_weighted_mean_MFE, S_phase_weighted_mean_MFE, G2_phase_weighted_mean_MFE),
  SD = c(G1_phase_weighted_sd_MFE, S_phase_weighted_sd_MFE, G2_phase_weighted_sd_MFE)
)

# Plot with error bars
Plot_weighted_MFE <- ggplot(weighted_mean_MFE_df, aes(x = Phase, y = Weighted_MFE)) +
  geom_point(size = 6, color = "black") +  # Dots
  geom_errorbar(aes(ymin = Weighted_MFE - SD, ymax = Weighted_MFE + SD), width = 0.2, color = "black") +  # Error bars
  theme_classic() +
  ylim(-0.45, -0.22) +
  theme(
    legend.position = "none", 
    axis.title = element_text(size = 20), 
    axis.text.x = element_text(size = 16, color = "black"),
    plot.title = element_text(size = 12)
  ) +
  scale_x_discrete(limits = c("G2 Phase", "G1 Phase", "S Phase")) +  # Specify Phase order
  labs(y = "Weighted Mean DGfold / nt (kcal/mol)", x = "", title = "Weighted Mean MFE Across Cell Cycle Phases of HeLa") 

print(Plot_weighted_MFE)











# Generate a boxplot with the transcript-level weighted MFE values for each phase
# Create a long-format data frame with transcript MFE values weighted by expression
long_df <- GSE104736_HeLa_TX_norm_counts %>%
  dplyr::mutate(G1_phase = Normalized_MFE * G1_phase) %>%
  dplyr::mutate(S_phase = Normalized_MFE * S_phase) %>%
  dplyr::mutate(G2_phase = Normalized_MFE * G2_phase) %>%
  dplyr::select(ensembl_transcript_id, Normalized_MFE, G1_phase, S_phase, G2_phase) %>%
  pivot_longer(cols = c(G1_phase, S_phase, G2_phase),
               names_to = "Phase",
               values_to = "Weighted_MFE")
# dplyr::filter(Weighted_MFE < 0)  # Filter out negative values

# Plot the distribution as a box plot
library(ggbreak)

boxplot_weighted_MFE <- ggplot(long_df, aes(x = Phase, y = Weighted_MFE)) +
  geom_boxplot(fill = "lightgrey", outlier.color = "black", outlier.shape = 1) +
  theme_classic() +
  ylim(-50000, 0.25) +
  theme(
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 8, color = "black", angle = 45, hjust = 1),
    axis.text.x.top = element_blank(),   # <== removes the top x-axis labels
    axis.ticks.x.top = element_blank(),  # <== removes the top x-axis ticks
    axis.line.x.top = element_blank(),        # Remove top x-axis line
    axis.line.y.right = element_blank(),      # Remove right-side y-axis line
    axis.text.y.right = element_blank(),        # Remove right-side y-axis labels
    axis.ticks.y.right = element_blank(),    # emove right-side y-axis ticks
    plot.title = element_text(size = 20)
  ) +
  scale_x_discrete(limits = c("G2_phase", "G1_phase", "S_phase")) +
  labs(y = expression(Delta*G[fold]~"/ nt (kcal/mol)"), x = "", title = "") +
  scale_y_break(c(-10000, -0.75), scales = 1.5) 

#print(boxplot_weighted_MFE)

# Save the combined plot
ggsave("boxplot_weighted_MFE.pdf", plot = boxplot_weighted_MFE, width = 3, height = 5, path = "path_to/Desktop")















## GSEA HeLa cell line ####

# Load the gene-level raw counts downloaded from GREIN
Counts_GSE104736_GREIN_gene_level <- read.delim("/Volumes/AC_HD_8TB/Exp.\ 277\ Analysis\ of\ Gibbs\ free\ energy\ in\ response\ to\ drug\ selection/Analysis_cell_cycle/Analysis_HeLa_and_U2OS_cell_cycle/GSE104736_raw_and\ normalized_counts/GSE104736_GeneLevel_Raw_data.csv",
                                     sep = ",", header = T, stringsAsFactors = F)
Counts_GSE104736_GREIN_gene_level <- Counts_GSE104736_GREIN_gene_level %>% remove_rownames %>% column_to_rownames(var="X") %>% # Convert first column into row names in data frame
  dplyr::select(-c(gene_symbol, GSM2806906, GSM2806907, GSM2806908, GSM2806909, GSM2806910, GSM2806911)) # Remove the columns that are not necessary for the analysis

# Analysis HeLa cells
# Data are converted and organized to the format required
Counts_GSE104736_GREIN_gene_level <- as.matrix(Counts_GSE104736_GREIN_gene_level)

# Round counts down to the nearest integer
Counts_GSE104736_GREIN_gene_level <- floor(Counts_GSE104736_GREIN_gene_level)

# IMPORTANT:
# It is absolutely critical that the columns of the count matrix and the rows of the column data (information about samples) are in the same order.

# Replace the row names of Expdata_GSE104736 (sampleID) with the same sample name [title] of [Counts_GSE104736]. 
# This will help to match sample names in both matrices.
# The outcome should be TRUE

all(rownames(Expdata_GSE104736) %in% colnames(Counts_GSE104736_GREIN_gene_level))
all(colnames(Counts_GSE104736_GREIN_gene_level) %in% rownames(Expdata_GSE104736))

# With the count matrix, Counts_GSE104736_GREIN_gene_level, and the sample information, Expdata_GSE104736, we can construct a DESeqDataSet:
dds_GSE104736_gene <- DESeqDataSetFromMatrix(countData = Counts_GSE104736_GREIN_gene_level,
                                        colData = Expdata_GSE104736,
                                        design = ~ condition)
dds_GSE104736_gene



# DATA QUALITY ASSESSMENT
# Before performing differential expression analysis, we perform data quality assessment by sample clustering and visualization
library("pheatmap")
# We log transform the counts
ntd <- normTransform(dds_GSE104736_gene)
# Heatmap of the count matrix
select <- order(rowMeans(counts(dds_GSE104736_gene,normalized=FALSE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds_GSE104736_gene)[,c("condition","type")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)



# Heatmap of the sample-to-sample distances
# A heatmap of this distance matrix gives us an overview over similarities and dissimilarities between samples
# Transform the values
vsd <- vst(dds_GSE104736_gene, blind=FALSE)
# rld <- rlog(dds_GSE104736_gene, blind=FALSE)
head(assay(vsd), 3)

# Apply the dist function to the transpose of the transformed count matrix to get sample-to-sample distances
sampleDists <- dist(t(assay(vsd)))

library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

# Principal component plot of the samples
plotPCA(vsd, intgroup="condition")

# Add the names of the samples to the PCA plot
# Extract PCA data
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
# Add sample names from colnames
pcaData$sample <- rownames(pcaData)
# Plot with sample labels
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color = condition, label = sample)) +
  geom_point(size = 3) +
  geom_text_repel(size = 3) +  # Use geom_text() instead if you don't want repel
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal()





# DIFFERENTIAL EXPRESSION ANALYSIS
# We perform a minimal pre-filtering to keep only rows that have at least 3 reads total
keep <- rowSums(counts(dds_GSE104736_gene)) >= 3
dds_GSE104736_gene <- dds_GSE104736_gene[keep,]
table(keep) # This shows how many transcripts have been removed after the filtering
nrow(dds_GSE104736_gene) # to confirm the filtered transcripts

# DIFFERENTIAL EXPRESSION ANALYSIS
dds_GSE104736_gene <- DESeq(dds_GSE104736_gene) # The DESeq() function normalizes the read counts, estimates dispersions, and fits the linear model, all in one go
# dds_GSE104736_gene <- DESeq(dds_GSE104736_gene, betaPrior = T)  The betaPrior = T opton tells DESeq2 to squeeze the log Fold Changes of lowly expressed genes toward zero.





# Differential gene expression G1 vs G2
# Specify the conditions to contrast for differential expression analysis
res_GSE104736_G1_vs_G2_gene <- results(dds_GSE104736_gene, contrast=c("condition","HeLa_G1","HeLa_G2")) # to specify the samples to contrast
res_GSE104736_G1_vs_G2_gene
# Get a summary of the results
summary(res_GSE104736_G1_vs_G2_gene)

# Plot the log2 Fold Changes as function of expresson:
plotMA(res_GSE104736_G1_vs_G2_gene,alpha = 0.05)
abline(h=c(-2,2),col = "red")

# Transform the row names into a column
res_GSE104736_df <- data.frame(res_GSE104736_G1_vs_G2_gene)
res_GSE104736_G1_vs_G2_gene_df <- tibble::rownames_to_column(res_GSE104736_df, "ensembl_gene_id") %>%
  dplyr::select(ensembl_gene_id, log2FoldChange)

# Add a column in each data frame with the entrez gene IDs required for the GSEA analysis
library(biomaRt)
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host="https://useast.ensembl.org")

gene_ids <- (res_GSE104736_G1_vs_G2_gene_df$ensembl_gene_id) 
res_GSE104736_gene_entrez_id_G1vsG2 <- getBM(attributes = c('entrezgene_id',
                                                     'ensembl_gene_id'),
                                      filters = 'ensembl_gene_id',
                                      values = gene_ids,
                                      mart = mart)
res_GSE104736_gene_entrez_id_G1vsG2 <- dplyr::distinct(res_GSE104736_gene_entrez_id_G1vsG2, ensembl_gene_id, .keep_all = TRUE) %>%
  dplyr::filter(!is.na(entrezgene_id)) 
GSE104736_gene_entrez_id_G1vsG2_entrez_IDs <- dplyr::left_join(res_GSE104736_G1_vs_G2_gene_df, res_GSE104736_gene_entrez_id_G1vsG2, by = "ensembl_gene_id") %>%
  dplyr::rename(G1vsG2_log2FoldChange = log2FoldChange)









# Differential gene expression S vs G1
# Specify the conditions to contrast for differential expression analysis
res_GSE104736_S_vs_G1_gene <- results(dds_GSE104736_gene, contrast=c("condition","HeLa_S","HeLa_G1")) # to specify the samples to contrast
res_GSE104736_S_vs_G1_gene
# Get a summary of the results
summary(res_GSE104736_S_vs_G1_gene)

# Plot the log2 Fold Changes as function of expresson:
plotMA(res_GSE104736_S_vs_G1_gene,alpha = 0.05)
abline(h=c(-2,2),col = "red")

# Transform the row names into a column
res_GSE104736_df <- data.frame(res_GSE104736_S_vs_G1_gene)
res_GSE104736_S_vs_G1_gene_df <- tibble::rownames_to_column(res_GSE104736_df, "ensembl_gene_id") %>%
  dplyr::select(ensembl_gene_id, log2FoldChange)

# Add a column in each data frame with the entrez gene IDs required for the GSEA analysis
library(biomaRt)
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host="https://useast.ensembl.org")

gene_ids <- (res_GSE104736_S_vs_G1_gene_df$ensembl_gene_id) 
res_GSE104736_gene_entrez_id_SvsG1 <- getBM(attributes = c('entrezgene_id',
                                                            'ensembl_gene_id'),
                                             filters = 'ensembl_gene_id',
                                             values = gene_ids,
                                             mart = mart)
res_GSE104736_gene_entrez_id_SvsG1 <- dplyr::distinct(res_GSE104736_gene_entrez_id_SvsG1, ensembl_gene_id, .keep_all = TRUE) %>%
  dplyr::filter(!is.na(entrezgene_id)) 
GSE104736_gene_entrez_id_SvsG1_entrez_IDs <- dplyr::left_join(res_GSE104736_S_vs_G1_gene_df, res_GSE104736_gene_entrez_id_SvsG1, by = "ensembl_gene_id") %>%
  dplyr::rename(SvsG1_log2FoldChange = log2FoldChange)










# Differential gene expression G2 vs S
# Specify the conditions to contrast for differential expression analysis
res_GSE104736_G2_vs_S_gene <- results(dds_GSE104736_gene, contrast=c("condition","HeLa_G2","HeLa_S")) # to specify the samples to contrast
res_GSE104736_G2_vs_S_gene
# Get a summary of the results
summary(res_GSE104736_G2_vs_S_gene)

# Plot the log2 Fold Changes as function of expresson:
plotMA(res_GSE104736_G2_vs_S_gene,alpha = 0.05)
abline(h=c(-2,2),col = "red")

# Transform the row names into a column
res_GSE104736_df <- data.frame(res_GSE104736_G2_vs_S_gene)
res_GSE104736_G2_vs_S_gene_df <- tibble::rownames_to_column(res_GSE104736_df, "ensembl_gene_id") %>%
  dplyr::select(ensembl_gene_id, log2FoldChange)

# Add a column in each data frame with the entrez gene IDs required for the GSEA analysis
library(biomaRt)
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host="https://useast.ensembl.org")

gene_ids <- (res_GSE104736_G2_vs_S_gene_df$ensembl_gene_id) 
res_GSE104736_gene_entrez_id_G2vsS <- getBM(attributes = c('entrezgene_id',
                                                           'ensembl_gene_id'),
                                            filters = 'ensembl_gene_id',
                                            values = gene_ids,
                                            mart = mart)
res_GSE104736_gene_entrez_id_G2vsS <- dplyr::distinct(res_GSE104736_gene_entrez_id_G2vsS, ensembl_gene_id, .keep_all = TRUE) %>%
  dplyr::filter(!is.na(entrezgene_id)) 
GSE104736_gene_entrez_id_G2vsS_entrez_IDs <- dplyr::left_join(res_GSE104736_G2_vs_S_gene_df, res_GSE104736_gene_entrez_id_G2vsS, by = "ensembl_gene_id") %>%
  dplyr::rename(G2vsS_log2FoldChange = log2FoldChange)









# Generate a data frame with the gene-level Log2 FC values for each comparison
GSE104736_gene_level_Log2_FC_df <- GSE104736_gene_entrez_id_G1vsG2_entrez_IDs %>%
  dplyr::select(ensembl_gene_id, entrezgene_id, G1vsG2_log2FoldChange) %>%
  dplyr::left_join(
    GSE104736_gene_entrez_id_SvsG1_entrez_IDs %>%
      dplyr::select(ensembl_gene_id, entrezgene_id, SvsG1_log2FoldChange),
    by = c("ensembl_gene_id", "entrezgene_id")
  ) %>%
  dplyr::left_join(
    GSE104736_gene_entrez_id_G2vsS_entrez_IDs %>%
      dplyr::select(ensembl_gene_id, entrezgene_id, G2vsS_log2FoldChange),
    by = c("ensembl_gene_id", "entrezgene_id")
  )














# Generate a function which performs GSEA of the samples

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
data_frame <- GSE104736_gene_level_Log2_FC_df
log_fc_cols <- c("G1vsG2_log2FoldChange", "SvsG1_log2FoldChange", "G2vsS_log2FoldChange")

plot_titles <- c("GSE104736_G1vsG2", "GSE104736_SvsG1", "GSE104736_G2vsS")

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

# Write the file in a csv format
write.csv(combined_kk2_df, file = "path_to/combined_kk2_df.csv")







# Generate a heat map with the kk2 results
library(pheatmap)

# Additional code to generate the heatmap
# Filter only the pathways of interest

# Selected KEGG signaling pathways
Selected_KEGG_pathways <- c("hsa04150", "hsa04668", "hsa04657", "hsa03030", "hsa03410", "hsa03430", "hsa03420", "hsa04370", "hsa4072",
                            "hsa03320", "hsa03460", "hsa03440", "hsa04014", "hsa04070", "hsa04020", "hsa01230",
                            "hsa04310", "hsa03082", "hsa04015", "hsa04625", "hsa04371", "hsa04390", "hsa04010", "hsa04622", "hsa04350",
                            "hsa04152", "hsa04926", "hsa04630", "hsa04115", "hsa04151", "hsa04910", "hsa04330", "hsa04012",
                            "hsa03440", "hsa03460") # hsa04110 is the cell cycle gene set

combined_kk2_df_filtered <- dplyr::filter(combined_kk2_df, ID %in% Selected_KEGG_pathways) %>%
  dplyr::mutate(NES = scales::rescale(NES, to = c(-4, 4))) # set a specific range of values to have the 0 centered to white



generate_heatmap <- function(kk2_df) {
  kk2_df <- kk2_df %>%
    dplyr::select(Description, Sample, NES) %>%
    tidyr::spread(Sample, NES)
  
  kk2_df[is.na(kk2_df)] <- 0
  rownames(kk2_df) <- kk2_df$Description
  kk2_df <- kk2_df[, -1]
  
  desired_order <- c("GSE104736_SvsG1", "GSE104736_G2vsS")
  kk2_df <- kk2_df[, desired_order]
  
  pathway_order <- c(
    "Biosynthesis of amino acids",
    "DNA replication",
    "Homologous recombination",
    "Fanconi anemia pathway",
    "Base excision repair"
  )
  
  kk2_df <- kk2_df[pathway_order, , drop = FALSE]
  
  display_matrix <- as.matrix(kk2_df)
  display_numbers <- ifelse(display_matrix == 0, "", sprintf("%.2g", display_matrix))
  
  heatmap_plot <- pheatmap(as.matrix(kk2_df), 
                           color = colorRampPalette(c("blue", "white", "red"))(50), 
                           cluster_rows = FALSE, 
                           cluster_cols = FALSE, 
                           display_numbers = display_numbers, 
                           number_format = "%.2g", 
                           main = "Heatmap of KEGG signaling pathways",
                           silent = TRUE)
  
  return(list(plot = heatmap_plot, data = kk2_df))
}

# Now you can access both:
heatmap_result <- generate_heatmap(combined_kk2_df_filtered)

# Access the plot
Heatmap_cell_cycle <- heatmap_result$plot

# Access the matrix
heatmap_df <- heatmap_result$data
View(heatmap_df)






# Plot the data
desired_order <- c(     "Biosynthesis of amino acids",
                        "DNA replication",
                        "Homologous recombination",
                        "Fanconi anemia pathway",
                        "Base excision repair"
                    )


# Reorder the 'Description' column based on the desired order
GSE104736_combined_kk2_df_filtered_plot <- combined_kk2_df_filtered %>%
  dplyr::mutate(Description = factor(Description, levels = desired_order),
                neg_log10_p_adjust = -log10(p.adjust))

GSE104736_combined_kk2_df_filtered_plot <- GSE104736_combined_kk2_df_filtered_plot %>%
  dplyr::mutate(Sample = factor(Sample, levels = c("GSE104736_SvsG1", "GSE104736_G2vsS")))


# Generate the plot with circles
Plot_cell_proliferation <- ggplot(GSE104736_combined_kk2_df_filtered_plot, aes(x = Sample, y = Description, fill = NES, size = neg_log10_p_adjust)) +
  geom_point(shape = 21, color = "black", stroke = 0.5) +  # Use shape = 21 for circles with border color
  scale_fill_gradient2(low = "dodgerblue4", mid = "white", high = "firebrick3", midpoint = 0, limits = c(-4, 4)) +  # Set colors and midpoint for the fill gradient
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

# Save the plot
ggsave("Plot_cell_proliferation.pdf", plot = Plot_cell_proliferation, width = 6, height = 10, path = "path_to/Desktop")








# To save the environment
save.image(file = "path_to/Analysis_cell_cycle_genes_MFE_environment.RData")















