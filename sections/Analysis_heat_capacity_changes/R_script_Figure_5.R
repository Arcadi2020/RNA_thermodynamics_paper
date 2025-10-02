# This script generates the plots of Fig.5


# Overview: Fig. 5 – GREIN/Science datasets: logFC aggregations, heat-capacity contrasts, and summary plots
#
# What this script does (high level)
# • Loads precomputed environments (GREIN, iLINCS, Science models, cell-cycle sets).
# • Builds a wide gene × sample table of DESeq2 log2FC across many GREIN studies + Science models (72h & early phase).
# • For selected datasets, contrasts “heat capacity” (Specific_heat) between UP (logFC>0) vs DOWN (logFC<0):
#   - All significant transcripts and “high-MFE” subsets (where provided).
#   - Statistics: two-sample t-test, Kolmogorov–Smirnov, and Wilcoxon.
# • Assembles multi-panel box/violin plots (GREIN & Science models) and compact summary tables of p-values/means.
# • Generates iLINCS “heat capacity volcano” (Log10_Specific_heat_change vs −log10 p) + top-aligned boxplot and one-sample t-test.
# • Benchmarks cell-cycle contrasts (MCF7 & HeLa phases) with the same Specific_heat UP vs DOWN tests and plots.
#
# Required inputs (update 'path_to/...'):
# • Cross_validation_MFE_AC_environment.RData (GREIN/Science preprocessed objects)
# • iLINCS environment + CSVs: Results_all_cancer_cell_lines_10uM.csv, cellinfo_beta.csv,
#   iLINCs_complete_data_selected_samples.csv
# • Science model DESeq/DEXSeq tables (72h & early phase) for X778, SKBR3, SKMEL28
# • Table_ANNOTATED_iLINCs.csv (ensembl_id, Normalized_MFE, Specific_heat)
# • Cell-cycle environments: Analysis_cell_cycle_genes_MFE_environment.RData
#
# Key steps / flow
# 1) GREIN logFC matrix
#    - Select {Gene_symbol, log2FoldChange} per dataset; rename columns; full_join by Gene_symbol ➜ Expression_changes_targeted_agents.
# 2) Heat capacity (GREIN)
#    - For selected GSEs: boxplots of Specific_heat ~ I(log2FC>0) (all sig TX); optional “high-MFE TX” subset.
#    - Compute t-test, KS, Wilcoxon; bind into TTEST_KS_results (+ high-MFE table).
#    - Arrange 5 GREIN panels in a single row (gridExtra / patchwork).
# 3) Science models (72h & early phase)
#    - Plot Specific_heat UP vs DOWN for X778 palbociclib (72h & EP) + stats table.
#    - Load summaries; compare proliferation signatures (Signt_1/2/3) 72h vs EP (paired Wilcoxon per signature).
# 4) iLINCS heat-capacity volcano (Fig. 5E)
#    - Compute Log10_Specific_heat_change = log10(MeanSpecific_heat_DOWN/MeanSpecific_heat_UP).
#    - Volcano: Log10_Specific_heat_change vs −log10(Wilcoxon p); threshold at 0.05 (≈1.301).
#    - Add top boxplot; one-sample t-test for mean(Log10_Specific_heat_change) > 0.
# 5) iLINCS per-sample boxplots (Fig. 5F)
#    - Select exemplars (e.g., PC3_torin-1, A375_tozasertib, ES2_lapatinib) with ks p<0.05.
#    - Join to Table_ANNOTATED_iLINCs; boxplots of Specific_heat by UP/DOWN; save one PDF per sample.
# 6) Cell-cycle contrasts (Fig. 5G)
#    - MCF7 (G0/G1 vs G2/M; S vs G0/G1) and HeLa (G1 vs G2; S vs G1): boxplots of Specific_heat; t-test, KS, Wilcoxon.
#
# Outputs
# • Expression_changes_targeted_agents (wide logFC table).
# • GREIN/Science boxplot panels (combined and per-dataset PDFs).
# • TTEST_KS_results (+ optional high-MFE table) and per-model “Statistics” data frames.
# • iLINCS heat-capacity volcano + boxplot composite; t-test readout.
# • Per-sample Specific_heat PDFs for selected iLINCS signatures.




# Analysis GREIN data

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



# Load the environment:
load("path_to/Cross_validation_MFE_AC_environment.RData")




# Here below I perform the analysis of all the datasets treated with agents targeting specific signaling pathways regulating cell proliferation.



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



















# Fig. 5B Heat_capacity GREIN dataset ####
# Heat capacity plots





# Heat capacity changes GSE143462_omipalisib 12h
# GSE143462_omipalisib box plot heat capacity all transcripts UP- and DOWN-regulated
Plot_Heat_capacity_GSE143462_omipalisib_12h <- ggplot(MFE_res_GSE143462, aes(x = factor(log2FoldChange > 0), y = Specific_heat)) +
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
  labs(y="Heat capacity (J/C)",x="",title="GSE143462_Log2_FC_adj.P<0.05 heat_capacity all transcripts") +
  scale_x_discrete(labels = c("log2FoldChange < 0", "log2FoldChange > 0"))

print(Plot_Heat_capacity_GSE143462_omipalisib_12h)

# TTEST and KS test of the heat capacity values of all transcripts significantly UP- and DOWN-regulated
TTEST_heat_capacity_GSE143462_all_UP_vs_DOWN
ks_test_GSE143462_all_UP_vs_DOWN_heat_capacity



# GSE143462_omipalisib 12h box plot heat capacity high MFE transcripts
Plot_Heat_capacity_high_GSE143462_omipalisib_12h <- ggplot(GSE143462_Specific_heat_data, aes(x = regulation, y = GSE143462_Specific_heat, fill = regulation)) +
  geom_boxplot(fill = "lightgrey", alpha = 0.7) +
  labs(title = "GSE143462 specific Heat of High MFE Transcripts",
       x = "",
       y = "") +
  scale_x_discrete(labels = c("log2_FC < 0", "log2_FC > 0")) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 16),
        axis.title = element_text(size = 25))

print(Plot_Heat_capacity_high_GSE143462_omipalisib_12h)

# TTEST and KS test of the heat capacity values of high MFE transcripts
TTEST_heat_capacity_GSE143462_all_high_MFE_TX_UP_vs_DOWN
ks_test_heat_capacity_GSE143462_all_high_MFE_TX_UP_vs_DOWN










# Heat capacity changes GSE121105_trastuzumab 72h
# GSE121105_trastuzumab box plot heat capacity all transcripts UP- and DOWN-regulated
Plot_Heat_capacity_GSE121105_trastuzumab_72h <- ggplot(MFE_res_GSE121105_trastuzumab, aes(x = factor(log2FoldChange > 0), y = Specific_heat)) +
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
  labs(y="Heat capacity (J/C)",x="",title="GSE121105_trastuzumab_Log2_FC_adj.P<0.05 heat_capacity all transcripts") +
  scale_x_discrete(labels = c("log2FoldChange < 0", "log2FoldChange > 0"))

print(Plot_Heat_capacity_GSE121105_trastuzumab_72h)

# TTEST and KS test of the heat capacity values of all transcripts significantly UP- and DOWN-regulated
TTEST_heat_capacity_GSE121105_trastuzumab_all_UP_vs_DOWN
ks_test_GSE121105_trastuzumab_all_UP_vs_DOWN_heat_capacity



# GSE121105_trastuzumab 72h box plot heat capacity high MFE transcripts
Plot_Heat_capacity_high_GSE121105_trastuzumab_72h <- ggplot(GSE121105_trastuzumab_Specific_heat_data, aes(x = regulation, y = GSE121105_trastuzumab_Specific_heat, fill = regulation)) +
  geom_boxplot(fill = "lightgrey", alpha = 0.7) +
  labs(title = "GSE121105_trastuzumab specific Heat of High MFE Transcripts",
       x = "",
       y = "") +
  scale_x_discrete(labels = c("log2_FC < 0", "log2_FC > 0")) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 16),
        axis.title = element_text(size = 25))

print(Plot_Heat_capacity_high_GSE121105_trastuzumab_72h)

# TTEST and KS test of the heat capacity values of high MFE transcripts
TTEST_heat_capacity_GSE121105_trastuzumab_all_high_MFE_TX_UP_vs_DOWN
ks_test_heat_capacity_GSE121105_trastuzumab_all_high_MFE_TX_UP_vs_DOWN









# Heat capacity changes GSE103021_ARS_1620_LU65_24h
# GSE103021_ARS_1620_LU65_24h box plot heat capacity all transcripts UP- and DOWN-regulated
Plot_Heat_capacity_GSE103021_ARS_1620_LU65_24h <- ggplot(MFE_res_GSE103021_ARS_1620_LU65, aes(x = factor(log2FoldChange > 0), y = Specific_heat)) +
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
  labs(y="Heat capacity (J/C)",x="",title="GSE103021_ARS_1620_LU65_Log2_FC_adj.P<0.05 heat_capacity all transcripts") +
  scale_x_discrete(labels = c("log2FoldChange < 0", "log2FoldChange > 0"))

print(Plot_Heat_capacity_GSE103021_ARS_1620_LU65_24h)

# TTEST and KS test of the heat capacity values of all transcripts significantly UP- and DOWN-regulated
TTEST_heat_capacity_GSE103021_ARS_1620_LU65_all_UP_vs_DOWN
ks_test_GSE103021_ARS_1620_LU65_all_UP_vs_DOWN_heat_capacity



# GSE103021_ARS_1620_LU65_24h box plot heat capacity high MFE transcripts
Plot_Heat_capacity_high_GSE103021_ARS_1620_LU65_24h <- ggplot(GSE103021_ARS_1620_LU65_Specific_heat_data, aes(x = regulation, y = GSE103021_ARS_1620_LU65_Specific_heat, fill = regulation)) +
  geom_boxplot(fill = "lightgrey", alpha = 0.7) +
  labs(title = "GSE103021_ARS_1620_LU65 specific Heat of High MFE Transcripts",
       x = "",
       y = "") +
  scale_x_discrete(labels = c("log2_FC < 0", "log2_FC > 0")) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 16),
        axis.title = element_text(size = 25))

print(Plot_Heat_capacity_high_GSE103021_ARS_1620_LU65_24h)

# TTEST and KS test of the heat capacity values of high MFE transcripts
TTEST_heat_capacity_GSE103021_ARS_1620_LU65_all_high_MFE_TX_UP_vs_DOWN
ks_test_heat_capacity_GSE103021_ARS_1620_LU65_all_high_MFE_TX_UP_vs_DOWN












# GSE135514_Das_24h heat capacity changes
# GSE135514_Das_24h box plot heat capacity all transcripts UP- and DOWN-regulated
Plot_Heat_capacity_GSE135514_Das_24h <- ggplot(MFE_res_GSE135514_Das_24h, aes(x = factor(log2FoldChange > 0), y = Specific_heat)) +
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
  labs(y="Heat capacity (J/C)",x="",title="GSE135514_Das_24h_Log2_FC_adj.P<0.05 heat_capacity all transcripts") +
  scale_x_discrete(labels = c("log2FoldChange < 0", "log2FoldChange > 0"))

print(Plot_Heat_capacity_GSE135514_Das_24h)

# TTEST and KS test of the heat capacity values of all transcripts significantly UP- and DOWN-regulated
TTEST_heat_capacity_GSE135514_Das_24h_all_UP_vs_DOWN
ks_test_GSE135514_Das_24h_all_UP_vs_DOWN_heat_capacity



# GSE135514_Das_24h 24h box plot heat capacity high MFE transcripts
Plot_Heat_capacity_high_GSE135514_Das_24h <- ggplot(GSE135514_Das_24h_Specific_heat_data, aes(x = regulation, y = GSE135514_Das_24h_Specific_heat, fill = regulation)) +
  geom_boxplot(fill = "lightgrey", alpha = 0.7) +
  labs(title = "GSE135514_Das_24h specific Heat of High MFE Transcripts",
       x = "",
       y = "") +
  scale_x_discrete(labels = c("log2_FC < 0", "log2_FC > 0")) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 16),
        axis.title = element_text(size = 25))

print(Plot_Heat_capacity_high_GSE135514_Das_24h)

# TTEST and KS test of the heat capacity values of high MFE transcripts
TTEST_heat_capacity_GSE135514_Das_24h_all_high_MFE_TX_UP_vs_DOWN
ks_test_heat_capacity_GSE135514_Das_24h_all_high_MFE_TX_UP_vs_DOWN









# GSE135514_Das_72h heat capacity changes
# GSE135514_Das_72h box plot heat capacity all transcripts UP- and DOWN-regulated
Plot_Heat_capacity_GSE135514_Das_72h <- ggplot(MFE_res_GSE135514_Das_72h, aes(x = factor(log2FoldChange > 0), y = Specific_heat)) +
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
  labs(y="Heat capacity (J/C)",x="",title="GSE135514_Das_72h_Log2_FC_adj.P<0.05 heat_capacity all transcripts") +
  scale_x_discrete(labels = c("log2FoldChange < 0", "log2FoldChange > 0"))

print(Plot_Heat_capacity_GSE135514_Das_72h)

# TTEST and KS test of the heat capacity values of all transcripts significantly UP- and DOWN-regulated
TTEST_heat_capacity_GSE135514_Das_72h_all_UP_vs_DOWN
ks_test_GSE135514_Das_72h_all_UP_vs_DOWN_heat_capacity



# GSE135514_Das_72h 72h box plot heat capacity high MFE transcripts
Plot_Heat_capacity_high_GSE135514_Das_72h <- ggplot(GSE135514_Das_72h_Specific_heat_data, aes(x = regulation, y = GSE135514_Das_72h_Specific_heat, fill = regulation)) +
  geom_boxplot(fill = "lightgrey", alpha = 0.7) +
  labs(title = "GSE135514_Das_72h specific Heat of High MFE Transcripts",
       x = "",
       y = "") +
  scale_x_discrete(labels = c("log2_FC < 0", "log2_FC > 0")) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 16),
        axis.title = element_text(size = 25))

print(Plot_Heat_capacity_high_GSE135514_Das_72h)

# TTEST and KS test of the heat capacity values of high MFE transcripts
TTEST_heat_capacity_GSE135514_Das_72h_all_high_MFE_TX_UP_vs_DOWN
ks_test_heat_capacity_GSE135514_Das_72h_all_high_MFE_TX_UP_vs_DOWN







# Generate a data frame with the ttest and KS test results
TTEST_KS_results <- data.frame(
  Sample = c("GSE143462_omopalisib_12h", "GSE121105_trastuzumab_72h", 
             "GSE103021_LU65_ARS_1620_24h", "GSE135514_dasatinib_24h", 
             "GSE135514_dasatinib_72h"),
  Mean_heat_capacity_TX_UP = c(Mean_heat_capacity_GSE143462_UP,
                               Mean_heat_capacity_GSE121105_trastuzumab_UP,
                               Mean_heat_capacity_GSE103021_ARS_1620_LU65_UP,
                               Mean_heat_capacity_GSE135514_Das_24h_UP,
                               Mean_heat_capacity_GSE135514_Das_72h_UP
  ),
  Mean_heat_capacity_TX_DOWN = c(Mean_heat_capacity_GSE143462_DOWN,
                                 Mean_heat_capacity_GSE121105_trastuzumab_DOWN,
                                 Mean_heat_capacity_GSE103021_ARS_1620_LU65_DOWN,
                                 Mean_heat_capacity_GSE135514_Das_24h_DOWN,
                                 Mean_heat_capacity_GSE135514_Das_72h_DOWN
  ),
  TTEST_p_value = c(TTEST_heat_capacity_GSE143462_all_UP_vs_DOWN$p.value,
                    TTEST_heat_capacity_GSE121105_trastuzumab_all_UP_vs_DOWN$p.value,
                    TTEST_heat_capacity_GSE103021_ARS_1620_LU65_all_UP_vs_DOWN$p.value,
                    TTEST_heat_capacity_GSE135514_Das_24h_all_UP_vs_DOWN$p.value,
                    TTEST_heat_capacity_GSE135514_Das_72h_all_UP_vs_DOWN$p.value
  ),
  KS_p_value = c(ks_test_GSE143462_all_UP_vs_DOWN_heat_capacity$p.value,
                 ks_test_GSE121105_trastuzumab_all_UP_vs_DOWN_heat_capacity$p.value,
                 ks_test_GSE103021_ARS_1620_LU65_all_UP_vs_DOWN_heat_capacity$p.value,
                 ks_test_GSE135514_Das_24h_all_UP_vs_DOWN_heat_capacity$p.value,
                 ks_test_GSE135514_Das_72h_all_UP_vs_DOWN_heat_capacity$p.value
  ))

print(TTEST_KS_results)






# Combine all plots heat capacity all transcripts UP- and DOWN- regulated
plots <- list(
  Plot_Heat_capacity_GSE143462_omipalisib_12h,
  Plot_Heat_capacity_GSE121105_trastuzumab_72h,
  Plot_Heat_capacity_GSE103021_ARS_1620_LU65_24h,
  Plot_Heat_capacity_GSE135514_Das_24h,
  Plot_Heat_capacity_GSE135514_Das_72h
)

# Add margin and box plots to each plot
plots_with_margin <- lapply(plots, function(plot) {
  plot + 
    geom_boxplot(width = 0.1, color = NA, outlier.size = 0, alpha = 0) +
    theme(
      plot.margin = unit(c(4, 0.2, 0.2, 0.2), "cm"),
      axis.title = element_blank(),  # Remove axis titles
      axis.text.x = element_text(size = 6),
      plot.title = element_blank() # Remove titles
    )
})

# Arrange plots in 5 columns and 1 row with spacing
combined_plots <- grid.arrange(
  grobs = plots_with_margin,
  ncol = 5,
  nrow = 1
)










# Plot the data relative to the high MFE transcripts

# Generate a data frame with the ttest and KS test results
TTEST_KS_results_high_MFE_TX <- data.frame(
  Sample = c("GSE143462_omopalisib_12h", 
             "GSE121105_trastuzumab_72h", 
             "GSE103021_LU65_ARS_1620_24h",
             "GSE135514_dasatinib_24h", 
             "GSE135514_dasatinib_72h"),
  Mean_heat_capacity_TX_UP = c(Mean_heat_capacity_GSE143462_high_MFE_TX_UP,
                               Mean_heat_capacity_GSE121105_trastuzumab_high_MFE_TX_UP,
                               Mean_heat_capacity_GSE103021_ARS_1620_LU65_high_MFE_TX_UP,
                               Mean_heat_capacity_GSE135514_Das_24h_high_MFE_TX_UP,
                               Mean_heat_capacity_GSE135514_Das_72h_high_MFE_TX_UP
  ),
  Mean_heat_capacity_TX_DOWN = c(Mean_heat_capacity_GSE143462_high_MFE_TX_DOWN,
                                 Mean_heat_capacity_GSE121105_trastuzumab_high_MFE_TX_DOWN,
                                 Mean_heat_capacity_GSE103021_ARS_1620_LU65_high_MFE_TX_DOWN,
                                 Mean_heat_capacity_GSE135514_Das_24h_high_MFE_TX_DOWN,
                                 Mean_heat_capacity_GSE135514_Das_72h_high_MFE_TX_DOWN
  ),
  TTEST_p_value = c(TTEST_heat_capacity_GSE143462_all_high_MFE_TX_UP_vs_DOWN$p.value,
                    TTEST_heat_capacity_GSE121105_trastuzumab_all_high_MFE_TX_UP_vs_DOWN$p.value,
                    TTEST_heat_capacity_GSE103021_ARS_1620_LU65_all_high_MFE_TX_UP_vs_DOWN$p.value,
                    TTEST_heat_capacity_GSE135514_Das_24h_all_high_MFE_TX_UP_vs_DOWN$p.value,
                    TTEST_heat_capacity_GSE135514_Das_72h_all_high_MFE_TX_UP_vs_DOWN$p.value
  ),
  KS_p_value = c(ks_test_GSE143462_all_UP_vs_DOWN_heat_capacity$p.value,
                 ks_test_heat_capacity_GSE121105_trastuzumab_all_high_MFE_TX_UP_vs_DOWN$p.value,
                 ks_test_heat_capacity_GSE103021_ARS_1620_LU65_all_high_MFE_TX_UP_vs_DOWN$p.value,
                 ks_test_heat_capacity_GSE135514_Das_24h_all_high_MFE_TX_UP_vs_DOWN$p.value,
                 ks_test_heat_capacity_GSE135514_Das_72h_all_high_MFE_TX_UP_vs_DOWN$p.value
  ))

print(TTEST_KS_results_high_MFE_TX)







# Repeat the statistics using the Wilcoxon test

# Wilcoxon test for all transcripts (UP vs DOWN)
wilcox_heat_capacity_GSE143462_all_UP_vs_DOWN <- wilcox.test(
  Specific_heat ~ factor(log2FoldChange > 0), 
  data = MFE_res_GSE143462
)

wilcox_heat_capacity_GSE121105_trastuzumab_all_UP_vs_DOWN <- wilcox.test(
  Specific_heat ~ factor(log2FoldChange > 0), 
  data = MFE_res_GSE121105_trastuzumab
)

wilcox_heat_capacity_GSE103021_ARS_1620_LU65_all_UP_vs_DOWN <- wilcox.test(
  Specific_heat ~ factor(log2FoldChange > 0), 
  data = MFE_res_GSE103021_ARS_1620_LU65
)

wilcox_heat_capacity_GSE135514_Das_24h_all_UP_vs_DOWN <- wilcox.test(
  Specific_heat ~ factor(log2FoldChange > 0), 
  data = MFE_res_GSE135514_Das_24h
)

wilcox_heat_capacity_GSE135514_Das_72h_all_UP_vs_DOWN <- wilcox.test(
  Specific_heat ~ factor(log2FoldChange > 0), 
  data = MFE_res_GSE135514_Das_72h
)


# Wilcoxon test for high MFE transcripts (UP vs DOWN)
wilcox_heat_capacity_GSE143462_high_MFE_TX_UP_vs_DOWN <- wilcox.test(
  GSE143462_Specific_heat ~ regulation, 
  data = GSE143462_Specific_heat_data
)

wilcox_heat_capacity_GSE121105_trastuzumab_high_MFE_TX_UP_vs_DOWN <- wilcox.test(
  GSE121105_trastuzumab_Specific_heat ~ regulation, 
  data = GSE121105_trastuzumab_Specific_heat_data
)

wilcox_heat_capacity_GSE103021_ARS_1620_LU65_high_MFE_TX_UP_vs_DOWN <- wilcox.test(
  GSE103021_ARS_1620_LU65_Specific_heat ~ regulation, 
  data = GSE103021_ARS_1620_LU65_Specific_heat_data
)

wilcox_heat_capacity_GSE135514_Das_24h_high_MFE_TX_UP_vs_DOWN <- wilcox.test(
  GSE135514_Das_24h_Specific_heat ~ regulation, 
  data = GSE135514_Das_24h_Specific_heat_data
)

wilcox_heat_capacity_GSE135514_Das_72h_high_MFE_TX_UP_vs_DOWN <- wilcox.test(
  GSE135514_Das_72h_Specific_heat ~ regulation, 
  data = GSE135514_Das_72h_Specific_heat_data
)


# Add Wilcoxon p-values to the summary table
TTEST_KS_results$Wilcoxon_p_value <- c(
  wilcox_heat_capacity_GSE143462_all_UP_vs_DOWN$p.value,
  wilcox_heat_capacity_GSE121105_trastuzumab_all_UP_vs_DOWN$p.value,
  wilcox_heat_capacity_GSE103021_ARS_1620_LU65_all_UP_vs_DOWN$p.value,
  wilcox_heat_capacity_GSE135514_Das_24h_all_UP_vs_DOWN$p.value,
  wilcox_heat_capacity_GSE135514_Das_72h_all_UP_vs_DOWN$p.value
)

# # For high MFE transcripts
# TTEST_KS_results_high_MFE_TX$Wilcoxon_p_value <- c(
#   wilcox_heat_capacity_GSE143462_high_MFE_TX_UP_vs_DOWN$p.value,
#   wilcox_heat_capacity_GSE121105_trastuzumab_high_MFE_TX_UP_vs_DOWN$p.value,
#   wilcox_heat_capacity_GSE103021_ARS_1620_LU65_high_MFE_TX_UP_vs_DOWN$p.value,
#   wilcox_heat_capacity_GSE135514_Das_24h_high_MFE_TX_UP_vs_DOWN$p.value,
#   wilcox_heat_capacity_GSE135514_Das_72h_high_MFE_TX_UP_vs_DOWN$p.value
# )






# Combine all plots heat capacity high MFE transcripts

library(patchwork)

combined_plot <- (Plot_Heat_capacity_GSE143462_omipalisib_12h +
                    Plot_Heat_capacity_GSE121105_trastuzumab_72h +
                    Plot_Heat_capacity_GSE103021_ARS_1620_LU65_24h +
                    Plot_Heat_capacity_GSE135514_Das_24h +
                    Plot_Heat_capacity_GSE135514_Das_72h) +
  plot_layout(ncol = 5, guides = "collect") &  # 5 columns, collect legends
  theme(
    axis.text.x = element_blank(),    # remove x-axis text
    axis.ticks.x = element_blank(),   # remove x-axis ticks
    axis.title.x = element_blank(),   # remove x-axis title
    axis.title.y = element_blank(),   # remove y-axis title
    plot.title = element_blank()      # remove plot titles
    # leave y-axis text and ticks intact
  )























# Fig. 5C Heat_capacity ####

# Analysis heat capacity Science model systems

# Analysis 778 palbociclib
# Analysis heat capacity 778 palbociclib 72h all UP-regulated vs DOWN-regulated transcripts
Plot_heat_capacity_778_PALBO <- ggplot(X778_PALBO_dexseq_ANNOTATED, aes(x = factor(log2fold_treated_untreated > 0), y = Specific_heat)) +
  geom_boxplot(fill = "lightgrey", alpha = 0.7) +
  #geom_jitter(width = 0.1, size = 1, alpha = 0.3) +
  theme_classic() +
  ylim(0, 400) +
  theme(legend.position = "none", 
        axis.title = element_text(size = 20), 
        plot.title = element_text(size = 10),
        axis.text.x = element_text(size = 12, color = "black")) +
  theme(plot.title = element_text(size = 20)) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  labs(y="",x="",title="") +
  scale_x_discrete(labels = c("log2FoldChange < 0", "log2FoldChange > 0"))
print(Plot_heat_capacity_778_PALBO)

TTEST_heat_capacity_X778_PALBO_dexseq_ANNOTATED_all_UP_vs_DOWN
ks_test_X778_PALBO_dexseq_ANNOTATED_all_UP_vs_DOWN_heat_capacity





# Wilcoxon test of the heat capacity values of all transcripts UP- vs DOWN-regulated
wilcox_heat_capacity_X778_PALBO_dexseq_ANNOTATED_all_UP_vs_DOWN <- wilcox.test(
  Specific_heat ~ factor(log2fold_treated_untreated > 0),
  data = X778_PALBO_dexseq_ANNOTATED
)



# Generate a data frame with the ttest and KS test results
# Generate a data frame with the TTEST, KS, and Wilcoxon test results
Statistics <- data.frame(
  Sample = "WDLPS_PALBO_72h",
  Mean_heat_capacity_TX_UP = Mean_heat_capacity_X778_PALBO_dexseq_ANNOTATED_UP,
  Mean_heat_capacity_TX_DOWN = Mean_heat_capacity_X778_PALBO_dexseq_ANNOTATED_DOWN,
  TTEST_p_value = TTEST_heat_capacity_X778_PALBO_dexseq_ANNOTATED_all_UP_vs_DOWN$p.value,
  KS_p_value = ks_test_X778_PALBO_dexseq_ANNOTATED_all_UP_vs_DOWN_heat_capacity$p.value,
  Wilcoxon_p_value = wilcox_heat_capacity_X778_PALBO_dexseq_ANNOTATED_all_UP_vs_DOWN$p.value
)

print(Statistics)












































# Fig. 5D Heat capacity Science models EP ####


# Load the data of the Science models generated at 72h:
source("path_to/All_model_systems_72h_summary.R")

# Load the data of the Science models generated during the early phase response
source("path_to/All_model_systems_early_phase_summary.R")




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











# Analysis 778 palbociclib early-phase response

# Analysis 778 palbociclib early-phase response
Plot_heat_capacity_778_PALBO_EP <- ggplot(X778_counts_EP, aes(x = factor(log2fold_treated_untreated > 0), y = Specific_heat)) +
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
  labs(y="",x="",title="") +
  scale_x_discrete(labels = c("log2FoldChange < 0", "log2FoldChange > 0"))
print(Plot_heat_capacity_778_PALBO_EP)

TTEST_heat_capacity_X778_PALBO_EP_dexseq_ANNOTATED_all_UP_vs_DOWN
ks_test_X778_PALBO_EP_dexseq_ANNOTATED_all_UP_vs_DOWN_heat_capacity

# Wilcoxon test of heat capacity values (UP vs DOWN) for EP dataset
wilcox_heat_capacity_X778_PALBO_EP_dexseq_ANNOTATED_all_UP_vs_DOWN <- wilcox.test(
  Specific_heat ~ factor(log2fold_treated_untreated > 0),
  data = X778_counts_EP
)



# Generate a data frame with the ttest and KS test results
Statistics_results <- data.frame(
  Sample = "WDLPS_PALBO_EP",
  Mean_heat_capacity_TX_UP = c(Mean_heat_capacity_X778_PALBO_EP_dexseq_ANNOTATED_UP
  ),
  Mean_heat_capacity_TX_DOWN = c(Mean_heat_capacity_X778_PALBO_EP_dexseq_ANNOTATED_DOWN
  ),
  TTEST_p_value = c(TTEST_heat_capacity_X778_PALBO_EP_dexseq_ANNOTATED_all_UP_vs_DOWN$p.value
  ),
  Wilcoxon_p_value = c(wilcox_heat_capacity_X778_PALBO_EP_dexseq_ANNOTATED_all_UP_vs_DOWN$p.value
  ),
  KS_p_value = c(ks_test_X778_PALBO_EP_dexseq_ANNOTATED_all_UP_vs_DOWN_heat_capacity$p.value
  ))

print(Statistics_results)











# Combine all plots heat capacity all transcripts UP- and DOWN- regulated
library(patchwork)

combined_plot <- (Plot_heat_capacity_778_PALBO + 
                    Plot_heat_capacity_778_PALBO_EP) +
  plot_layout(ncol = 2, guides = "collect") & 
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.title = element_blank(),
    legend.position = "none"
  )


















# Fig. 5E volcano plot iLINCs Wilcoxon test ####

# Load the environment of the script output
load("path_to/Analysis_MFE_iLINCS_all_datasets_10uM_cancer_lines_environment_selected.RData")


library(readr)
library(dplyr)
library(biomaRt)
library(ggplot2)
library(ggpmisc)
library(ggrepel)
library(forcats)  
library(stringr)
library(pheatmap)
library(patchwork)
library(tidyr)
library(gridExtra)


# Load the complete iLINCs data
all_results_2_complete <- read.csv("path_to/Results_all_cancer_cell_lines_10uM.csv",
                                   header = TRUE, sep = ",", stringsAsFactors = FALSE)


# Add a new column called 'Log10_MFE_change' containing the log10 of the ratio between the mean MFE of the DOWN-regulated genes and the mean MFE of the UP-regulated genes
iLINCS_analysis <- dplyr::mutate(all_results_2_complete, Log10_MFE_change = log10(MeanMFE_DOWN/MeanMFE_UP)) %>%
  dplyr::relocate(MeanMFE_UP, .after = MeanMFE_DOWN) %>%
  dplyr::relocate(Log10_MFE_change, .after = MeanMFE_UP)






# Add information on cell lines and drugs and group the samples accordingly

# Open the iLINCs cell_info file
iLINCs_cell_info <- read.delim("path_to/cellinfo_beta.csv",
                               sep = ",", header = T, stringsAsFactors = F) %>%
  dplyr::rename(iLINCs_cell_line = cell_iname)


iLINCS_analysis <- dplyr::rename(iLINCS_analysis, iLINCs_cell_line = cell_line) %>%
  dplyr::select(-primary_disease)

# Merge the iLINCS_analysis data frame with iLINCs_cell_info on the iLINCs_cell_line and cell_iname columns
iLINCs_cell_info_complete <- dplyr::left_join(iLINCS_analysis, iLINCs_cell_info, by = "iLINCs_cell_line") %>%
  dplyr::select(iLINCs_cell_line, Drug, primary_disease) %>%
  dplyr::mutate(primary_disease = ifelse(iLINCs_cell_line == "P1A82", "teratocarcinoma", primary_disease))




# Add the drug targets
# In the original data frame, some of the drug targets, were missing. They have been added together with the relative signaling pathways. The file has been save as .csv file and loaded here below.
iLINCs_drug_info_complete <- read.delim("path_to/iLINCs_complete_data_selected_samples.csv",
                                        sep = ",", header = T, stringsAsFactors = F) %>%
  dplyr::select("Anti_cancer_agents", "TARGET", "TARGET_PATHWAY") %>%
  dplyr::rename(Drug = Anti_cancer_agents)



# Join the data frame to generate the complete plot
IlINCs_complete_df <- iLINCs_cell_info_complete %>%
  dplyr::left_join(iLINCs_drug_info_complete, by = c("Drug"), relationship = "many-to-many") %>%
  dplyr::mutate(TARGET_PATHWAY = str_replace(TARGET_PATHWAY, 'ATR-Chk1 DNA damage response', 'DNA damage response')) %>%
  dplyr::distinct(iLINCs_cell_line, Drug, primary_disease, TARGET, TARGET_PATHWAY)




# Plot the data
# Create a new column with primary_disease for grouping iLINCs_cell_lines
IlINCs_complete_df <- IlINCs_complete_df %>%
  mutate(primary_disease_label = factor(primary_disease, levels = unique(primary_disease)),
         TARGET_PATHWAY_label = factor(TARGET_PATHWAY, levels = unique(TARGET_PATHWAY)))


# Total n. of cell lines:
Total_cell_lines <- unique(IlINCs_complete_df$iLINCs_cell_line) # 98
# Total n. of anti-cancer agents:
Total_anti_cancer_agents <- unique(IlINCs_complete_df$Drug) # 42











iLINCS_analysis_plot_membrane_receptors <- ggplot(iLINCS_analysis, aes(x = Log10_Specific_heat_change, y = -log10(iLINCS_analysis$wilcox_test_all_UP_vs_DOWN_Specific_heat))) +
  theme_classic() +
  geom_point(aes(
    fill = ifelse(Log10_Specific_heat_change > 0 & -log10(iLINCS_analysis$wilcox_test_all_UP_vs_DOWN_Specific_heat) > 1.301, "blue", 
                  ifelse(Log10_Specific_heat_change > 0 & -log10(iLINCS_analysis$wilcox_test_all_UP_vs_DOWN_Specific_heat) < 1.301, "lightgrey",
                         ifelse(Log10_Specific_heat_change < 0 & -log10(iLINCS_analysis$wilcox_test_all_UP_vs_DOWN_Specific_heat) > 1.301, "lightgrey", 
                                ifelse(Log10_Specific_heat_change < 0 & -log10(iLINCS_analysis$wilcox_test_all_UP_vs_DOWN_Specific_heat) < 1.301, "lightgrey", "lightgrey"))))),
    color = "black",  # Border color for all points
    size = 4,  # Set a larger fixed size for all points
    shape = 21,  # Set a consistent shape for all points
    alpha = ifelse(-log10(iLINCS_analysis$wilcox_test_all_UP_vs_DOWN_Specific_heat) < 1.301, 0.6, 1),
    stroke = 0.6) +  # Stroke for empty circles' border
  scale_fill_identity() +
  scale_color_identity() +
  scale_alpha_identity() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 1.301, linetype = "dashed", color = "black") +
  labs(x = "Log10_Specific_heat_change", y = "-Log10_ks_test_all_UP_vs_DOWN_pValue", title = "") +
  theme(plot.title = element_text(color = "black", size = 14, face = "bold"),
        legend.position = "none") +
  xlim(-0.2, 0.2) +
  geom_text_repel(aes(label = ifelse(-log10(iLINCS_analysis$wilcox_test_all_UP_vs_DOWN_Specific_heat) > 1.301 & 
                                       paste(iLINCs_cell_line, Drug, sep = "_") %in% c("HT29_vemurafenib", "LNCAP_sorafenib", "HCC515_trametinib"),
                                     paste(iLINCs_cell_line, Drug, sep = "_"), "")), 
                  size = 3, box.padding = 0.5, point.padding = 0.5, nudge_x = 0.004, nudge_y = 0.001, max.overlaps = Inf, 
                  segment.color = 'lightgrey', segment.size = 0.5)

print(iLINCS_analysis_plot_membrane_receptors)
















# Add the name of the specific cell lines and pharmacologic agents to the plot
iLINCS_analysis_plot_membrane_receptors <- ggplot(iLINCS_analysis, aes(x = Log10_Specific_heat_change, y = -log10(wilcox_test_all_UP_vs_DOWN_Specific_heat))) +
  theme_classic() +
  geom_point(
    aes(
      fill = ifelse(
        Log10_Specific_heat_change > 0 & -log10(wilcox_test_all_UP_vs_DOWN_Specific_heat) > 1.301, "blue", 
        "lightgrey"
      ),
      #alpha = ifelse(-log10(wilcox_test_all_UP_vs_DOWN_Specific_heat) < 1.301, 0.6, 1)
    ),
    shape = 21,     # fixed value, placed outside aes()
    size = 4,
    stroke = 0.6,
    color = "black"
  ) +
  scale_fill_identity() +
  scale_alpha_identity() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 1.301, linetype = "dashed", color = "black") +
  labs(x = "Log10_Specific_heat_change", y = "-Log10_ks_test_all_UP_vs_DOWN_pValue", title = "") +
  theme(
    plot.title = element_text(color = "black", size = 14, face = "bold"),
    legend.position = "none"
  ) +
  xlim(-0.2, 0.2) +
  geom_text_repel(
    aes(label = ifelse(
      Log10_Specific_heat_change > 0 &
        -log10(wilcox_test_all_UP_vs_DOWN_Specific_heat) > 1.301 &
        paste(iLINCs_cell_line, Drug, sep = "_") %in% c("ES2_lapatinib", "A375_tozasertib", "PC3_torin-1"),
      paste(iLINCs_cell_line, Drug, sep = "_"),
      ""
    )),
    size = 3,
    box.padding = 0.5,
    point.padding = 0.5,
    nudge_x = 0.004,
    nudge_y = 0.001,
    max.overlaps = Inf,
    segment.color = 'lightgrey',
    segment.size = 0.5
  )

# Print the plot
print(iLINCS_analysis_plot_membrane_receptors)












# Classify samples into heat capacity quadrants
iLINCS_analysis <- iLINCS_analysis %>%
  mutate(HeatCapacity_Quadrant = case_when(
    Log10_Specific_heat_change > 0 & -log10(wilcox_test_all_UP_vs_DOWN_Specific_heat) > 1.301 ~ "UP_regulated_heat_capacity",
    Log10_Specific_heat_change > 0 & -log10(wilcox_test_all_UP_vs_DOWN_Specific_heat) < 1.301 ~ "Q2",
    Log10_Specific_heat_change < 0 & -log10(wilcox_test_all_UP_vs_DOWN_Specific_heat) > 1.301 ~ "DOWN_regulated_heat_capacity",
    Log10_Specific_heat_change < 0 & -log10(wilcox_test_all_UP_vs_DOWN_Specific_heat) < 1.301 ~ "Q4",
    TRUE ~ "Other"
  ))

# Percentages for UP-regulated heat capacity
N_total_UP_heat <- iLINCS_analysis %>%
  filter(HeatCapacity_Quadrant %in% c("UP_regulated_heat_capacity", "Q2")) %>%
  nrow()
print(N_total_UP_heat)
# Totally 76 cell lines UP-regulate transcripts with high heat capacity

Perc_total_UP_heat <- 100 * N_total_UP_heat / nrow(iLINCS_analysis)
print(Perc_total_UP_heat)
# 54.2% of samples showed upregulation of transcripts with high predicted heat capacity


N_sig_UP_heat <- iLINCS_analysis %>%
  filter(HeatCapacity_Quadrant == "UP_regulated_heat_capacity") %>%
  nrow()
print(N_sig_UP_heat)
# 22 cell lines show significant upregulation of transcripts with high predicted heat capacity

Perc_sig_UP_heat <- 100 * N_sig_UP_heat / nrow(iLINCS_analysis)
print(Perc_sig_UP_heat)

Perc_sig_within_UP_heat <- 100 * N_sig_UP_heat / N_total_UP_heat
print(Perc_sig_within_UP_heat)
# 28.9% of the cell lines that upregulate transcripts with high predicted heat capacity show significant upregulation of these transcripts

# Percentages for DOWN-regulated heat capacity
N_total_DOWN_heat <- iLINCS_analysis %>%
  filter(HeatCapacity_Quadrant %in% c("DOWN_regulated_heat_capacity", "Q4")) %>%
  nrow()

Perc_total_DOWN_heat <- 100 * N_total_DOWN_heat / nrow(iLINCS_analysis)
print(Perc_total_DOWN_heat)

N_sig_DOWN_heat <- iLINCS_analysis %>%
  filter(HeatCapacity_Quadrant == "DOWN_regulated_heat_capacity") %>%
  nrow()
print(N_sig_DOWN_heat)

Perc_sig_within_DOWN_heat <- 100 * N_sig_DOWN_heat / N_total_DOWN_heat
print(Perc_sig_within_DOWN_heat)

Perc_sig_DOWN_heat <- 100 * N_sig_DOWN_heat / nrow(iLINCS_analysis)
print(Perc_sig_DOWN_heat)

# ---- Volcano plot ----
volcano_plot <- ggplot(iLINCS_analysis, aes(x = Log10_Specific_heat_change, y = -log10(wilcox_test_all_UP_vs_DOWN_Specific_heat))) +
  theme_classic() +
  geom_point(aes(
    fill = case_when(
      Log10_Specific_heat_change > 0 & -log10(wilcox_test_all_UP_vs_DOWN_Specific_heat) > 1.301 ~ "blue",
      TRUE ~ "lightgrey"
    ),
    #alpha = ifelse(-log10(wilcox_test_all_UP_vs_DOWN_Specific_heat) < 1.301, 0.6, 1)
  ),
  shape = 21, size = 4, stroke = 0.6, color = "black") +
  scale_fill_identity() +
  scale_alpha_identity() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 1.301, linetype = "dashed", color = "black") +
  
  # Uncomment the geom_text_repel to remove the labels
  # # Label only selected samples
  # geom_text_repel(
  #   aes(label = ifelse(
  #     -log10(wilcox_test_all_UP_vs_DOWN_Specific_heat) > 1.301 &
  #       paste(iLINCs_cell_line, Drug, sep = "_") %in%
  #       c("ES2_lapatinib", "A375_tozasertib", "PC3_torin-1"),
  #     paste(iLINCs_cell_line, Drug, sep = "_"), "")
  #   ),
  #   size = 3,
  #   box.padding = 0.5,
  #   point.padding = 0.5,
  #   nudge_x = 0.004,
  #   nudge_y = 0.001,
  #   max.overlaps = Inf,
  #   segment.color = 'lightgrey',
  #   segment.size = 0.5
  # ) +
  
  labs(
    x = "Log10_Specific_heat_change",
    y = "-Log10_ks_test_all_UP_vs_DOWN_pValue",
    title = ""
  ) +
  xlim(-0.2, 0.2) +
  theme(
    plot.title = element_text(color = "black", size = 14, face = "bold"),
    legend.position = "none"
  )


# Horizontal boxplot for Log10_Specific_heat_change
boxplot_heat_capacity <- ggplot(iLINCS_analysis, aes(y = "", x = Log10_Specific_heat_change)) +
  theme_classic() +
  geom_boxplot(fill = "lightgrey", color = "black", width = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.8, color = "black") +
  labs(x = "", y = "") +
  xlim(-0.2, 0.2) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank()
  )

# Combine boxplot and volcano plot using patchwork
combined_plot <- boxplot_heat_capacity / volcano_plot + plot_layout(heights = c(0.07, 0.93))


print(combined_plot)



# Null hypothesis (H0): mean(Log10_Specific_heat_change) = 0
# Alternative hypothesis (HA): mean(Log10_Specific_heat_change) > 0

test_result <- t.test(iLINCS_analysis$Log10_Specific_heat_change, mu = 0, alternative = "greater")
print(test_result)





















# Fig. 5F Box plots heat capacity ####


# Edit the names of the samples to match the iLINCs format
all_results_2_complete$iLINCs_cell_line_drug <- sub("^(([^_]*_[^_]*)_).*", "\\1", all_results_2_complete$SampleName)

# Remove the trailing underscore
all_results_2_complete$iLINCs_cell_line_drug <- gsub("_$", "", all_results_2_complete$iLINCs_cell_line_drug)

iLINCs_all_results_2_complete_filtered_Specific_heat <- dplyr::filter(all_results_2_complete,  iLINCs_cell_line_drug %in% c("PC3_torin-1", "A375_tozasertib", "ES2_lapatinib")) %>%
  #dplyr::filter(Log10_Specific_heat_change > 0) %>%
  dplyr::filter(ks_test_all_UP_vs_DOWN_Specific_heat < 0.05)  # Filter for significant samples


# Filter the result_list_2 list object to only include the selected samples
iLINCs_selected_samples_Specific_heat <- result_list_2[grep(paste(iLINCs_all_results_2_complete_filtered_Specific_heat$SampleName, collapse = "|"), names(result_list_2))]



# Open the Table_ANNOTATED_iLINCs
Table_ANNOTATED_iLINCs <- read.delim("path_to/Table_ANNOTATED_iLINCs.csv",
                                     sep = ",", header = T, stringsAsFactors = F)
Table_ANNOTATED_iLINCs <- dplyr::rename(Table_ANNOTATED_iLINCs, ensembl_id = ensembl_gene_id)




# Generate the data of the Specific_heat values
# Initialize plot_list as an empty list
plot_list <- list()

# Initialize all_results as an empty data frame with appropriate column names
all_results <- data.frame(SampleName = character(),
                          MeanSpecific_heat_UP = numeric(),
                          MeanSpecific_heat_DOWN = numeric(),
                          ks_test_all_UP_vs_DOWN_Specific_heat = numeric(),
                          stringsAsFactors = FALSE)

# Function to analyze the data and store plots
process_result_list <- function(Sample_name) { 
  message("Processing ", Sample_name)
  df <- iLINCs_selected_samples_Specific_heat[[Sample_name]] # Access individual data frames
  
  # Modify title_plot_1 to remove text after the second underscore
  split_title <- strsplit(Sample_name, "_")[[1]]
  if (length(split_title) > 2) {
    title_plot_1 <- paste(split_title[1:2], collapse = "_")
  } else {
    title_plot_1 <- Sample_name  # If there are less than 2 underscores, keep the original
  }
  
  title_plot_1 <- paste0(title_plot_1)  # Optionally add the rest of the title
  
  # Join the data with Table_ANNOTATED_iLINCs
  Sample_name_TX <- dplyr::right_join(df, Table_ANNOTATED_iLINCs, by = "ensembl_id")
  
  # Filter only statistically significant z-scores
  Sample_name_TX_filtered <- dplyr::filter(Sample_name_TX, p_value < 0.05)
  
  # Calculate means for up-regulated and down-regulated genes
  Specific_heat_Sample_name_UP <- dplyr::filter(Sample_name_TX_filtered, z_score > 0)
  Mean_Specific_heat_Sample_name_UP <- mean(Specific_heat_Sample_name_UP$Specific_heat, na.rm = TRUE)
  
  Specific_heat_Sample_name_DOWN <- dplyr::filter(Sample_name_TX_filtered, z_score < 0)
  Mean_Specific_heat_Sample_name_DOWN <- mean(Specific_heat_Sample_name_DOWN$Specific_heat, na.rm = TRUE)
  
  # Perform statistical tests
  TTEST <- t.test(Specific_heat_Sample_name_DOWN$Specific_heat, Specific_heat_Sample_name_UP$Specific_heat, paired = FALSE)
  
  ks_test_all_UP_vs_DOWN <- ks.test(Specific_heat_Sample_name_DOWN$Specific_heat, Specific_heat_Sample_name_UP$Specific_heat)
  
  wilcox_test_all_UP_vs_DOWN <- wilcox.test(Specific_heat_Sample_name_DOWN$Specific_heat, 
                                            Specific_heat_Sample_name_UP$Specific_heat)
  
  # Calculate the net Specific_heat change
  Sample_name_TX_net_Specific_heat_change <- Sample_name_TX_filtered %>% 
    dplyr::mutate(Net_Specific_heat_change = (Specific_heat * z_score) * (-1))
  Net_Specific_heat_change <- sum(Sample_name_TX_net_Specific_heat_change$Net_Specific_heat_change, na.rm = TRUE)
  
  ## Generate box plot 
  Specific_heat_Sample_name <- ggplot(Sample_name_TX_filtered, aes(x = factor(z_score > 0), y = Specific_heat)) +
    geom_boxplot(fill = "lightgrey", alpha = 0.7) +
    #geom_jitter(width = 0.1, size = 1, alpha = 0.3) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.title = element_text(size = 20),
      axis.text.x = element_blank(),  # Removes x-axis text
      plot.title = element_text(size = 10)
    ) +
    geom_vline(xintercept = 0, linewidth = 0.3) +
    labs(y = "", x = NULL, title = title_plot_1) +  # Use modified title
    scale_x_discrete(labels = c("z_score < 0", "z_score > 0"))
  
  
  # Save the plot individually for each sample
  ggsave(
    filename = paste0(Sample_name, "_Specific_heat_plot.pdf"),  # Create a unique filename
    plot = Specific_heat_Sample_name,
    width = 4,
    height = 5,
    path = "/Users/z3514253/Desktop"  # Specify the path to save the plot
  )
  
  # Store the plot in the list
  plot_list[[Sample_name]] <<- Specific_heat_Sample_name
  
  # Create a data frame with the results for the current sample
  # Create a data frame with the results for the current sample
  sample_results <- data.frame(
    SampleName = Sample_name,
    MeanSpecific_heat_UP = Mean_Specific_heat_Sample_name_UP,
    MeanSpecific_heat_DOWN = Mean_Specific_heat_Sample_name_DOWN,
    ks_test_all_UP_vs_DOWN_Specific_heat = ks_test_all_UP_vs_DOWN$p.value,
    wilcox_test_all_UP_vs_DOWN_Specific_heat = wilcox_test_all_UP_vs_DOWN$p.value  # Add Wilcoxon p-value
  )
  
  
  # Append the sample results to the all_results data frame
  all_results <<- rbind(all_results, sample_results)
}

# Run the function for each iLINCS data frame
Data_result_list <- lapply(names(iLINCs_selected_samples_Specific_heat), process_result_list)

# Convert to data frame to clearly see the data
# Bind rows and include a source column
Data_combined <- bind_rows(Data_result_list, .id = "Source") %>%
  dplyr::distinct(SampleName, .keep_all = TRUE)
print(Data_combined)















# Fig. 5G box plots heat capacity cell cycle ####

# To load the environment
load("path_to/Analysis_cell_cycle_genes_MFE_environment.RData")


# Heat capacity changes MCF7 G0G1 vs G2M
# MCF7 box plot heat capacity all transcripts UP- and DOWN-regulated


Plot_Heat_capacity_MCF7_G0G1_vs_G2M <- ggplot(GSE94479_G0G1vsG2M, aes(x = factor(log2FoldChange > 0), y = Specific_heat)) +
  geom_boxplot(fill = "lightgrey", alpha = 0.7) +
  # geom_jitter(width = 0.1, size = 1, alpha = 0.3) +  # optional
  theme_classic() +
  ylim(0, 400) +  # adjust if needed
  theme(
    legend.position = "none", 
    axis.title = element_text(size = 20), 
    axis.text.x = element_text(size = 10, color = "black"),
    plot.title = element_text(size = 20)
  ) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  labs(
    y = "", 
    x = "", 
    title = ""
  ) +
  scale_x_discrete(labels = c("log2FoldChange < 0", "log2FoldChange > 0"))

print(Plot_Heat_capacity_MCF7_G0G1_vs_G2M)



# TTEST and KS test of the heat capacity values of all transcripts significantly UP- and DOWN-regulated
HeatCap_MCF7_G0G1vsG2M_UP <- dplyr::filter(GSE94479_G0G1vsG2M, log2FoldChange > 0)
HeatCap_MCF7_G0G1vsG2M_DOWN <- dplyr::filter(GSE94479_G0G1vsG2M, log2FoldChange < 0)

TTEST_heat_capacity_MCF7_all_UP_vs_DOWN <- t.test(
  HeatCap_MCF7_G0G1vsG2M_DOWN$Specific_heat,
  HeatCap_MCF7_G0G1vsG2M_UP$Specific_heat,
  paired = FALSE
)
print(TTEST_heat_capacity_MCF7_all_UP_vs_DOWN)

ks_test_MCF7_all_UP_vs_DOWN_heat_capacity <- ks.test(
  HeatCap_MCF7_G0G1vsG2M_DOWN$Specific_heat,
  HeatCap_MCF7_G0G1vsG2M_UP$Specific_heat
)
print(ks_test_MCF7_all_UP_vs_DOWN_heat_capacity)

wilcox_test_MCF7_all_UP_vs_DOWN_heat_capacity <- wilcox.test(
  HeatCap_MCF7_G0G1vsG2M_DOWN$Specific_heat,
  HeatCap_MCF7_G0G1vsG2M_UP$Specific_heat
)
print(wilcox_test_MCF7_all_UP_vs_DOWN_heat_capacity)










# Heat capacity changes GSE104736 HeLa G1 vs G2
# Box plot heat capacity for all transcripts UP- and DOWN-regulated

Plot_Heat_capacity_HeLa_G1_vs_G2 <- ggplot(MFE_GSE104736_G1_vs_G2, aes(x = factor(log2FoldChange > 0), y = Specific_heat)) +
  geom_boxplot(fill = "lightgrey", alpha = 0.7) +
  # geom_jitter(width = 0.1, size = 1, alpha = 0.3) +  # optional
  theme_classic() +
  ylim(0, 400) +  # adjust based on data range
  theme(
    legend.position = "none", 
    axis.title = element_text(size = 20), 
    axis.text.x = element_text(size = 10, color = "black"),
    plot.title = element_text(size = 20)
  ) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  labs(
    y = "", 
    x = "", 
    title = ""
  ) +
  scale_x_discrete(labels = c("log2FoldChange < 0", "log2FoldChange > 0"))

print(Plot_Heat_capacity_HeLa_G1_vs_G2)



# TTEST and KS test of the heat capacity values of all transcripts significantly UP- and DOWN-regulated
HeatCap_HeLa_G1_vs_G2_UP <- dplyr::filter(MFE_GSE104736_G1_vs_G2, log2FoldChange > 0)
HeatCap_HeLa_G1_vs_G2_DOWN <- dplyr::filter(MFE_GSE104736_G1_vs_G2, log2FoldChange < 0)

TTEST_heat_capacity_HeLa_G1_vs_G2_UP_vs_DOWN <- t.test(
  HeatCap_HeLa_G1_vs_G2_DOWN$Specific_heat,
  HeatCap_HeLa_G1_vs_G2_UP$Specific_heat,
  paired = FALSE
)
print(TTEST_heat_capacity_HeLa_G1_vs_G2_UP_vs_DOWN)

ks_test_HeLa_G1_vs_G2_UP_vs_DOWN_heat_capacity <- ks.test(
  HeatCap_HeLa_G1_vs_G2_DOWN$Specific_heat,
  HeatCap_HeLa_G1_vs_G2_UP$Specific_heat
)
print(ks_test_HeLa_G1_vs_G2_UP_vs_DOWN_heat_capacity)

wilcox_test_HeLa_G1_vs_G2_UP_vs_DOWN_heat_capacity <- wilcox.test(
  HeatCap_HeLa_G1_vs_G2_DOWN$Specific_heat,
  HeatCap_HeLa_G1_vs_G2_UP$Specific_heat
)
print(wilcox_test_HeLa_G1_vs_G2_UP_vs_DOWN_heat_capacity)
















# Heat capacity changes MCF7 S vs G0G1
# MCF7 box plot heat capacity all transcripts UP- and DOWN-regulated

Plot_Heat_capacity_MCF7_SvsG0G1 <- ggplot(GSE94479_SvsG0G1, aes(x = factor(log2FoldChange > 0), y = Specific_heat)) +
  geom_boxplot(fill = "lightgrey", alpha = 0.7) +
  # geom_jitter(width = 0.1, size = 1, alpha = 0.3) +  # optional
  theme_classic() +
  ylim(0, 400) +  # adjust based on your data distribution
  theme(
    legend.position = "none", 
    axis.title = element_text(size = 20), 
    axis.text.x = element_text(size = 10, color = "black"),
    plot.title = element_text(size = 20)
  ) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  labs(
    y = "", 
    x = "", 
    title = ""
  ) +
  scale_x_discrete(labels = c("log2FoldChange < 0", "log2FoldChange > 0"))

print(Plot_Heat_capacity_MCF7_SvsG0G1)




# TTEST and KS test of the heat capacity values of all transcripts significantly UP- and DOWN-regulated
HeatCap_MCF7_SvsG0G1_UP <- dplyr::filter(GSE94479_SvsG0G1, log2FoldChange > 0)
HeatCap_MCF7_SvsG0G1_DOWN <- dplyr::filter(GSE94479_SvsG0G1, log2FoldChange < 0)

TTEST_heat_capacity_MCF7_SvsG0G1_UP_vs_DOWN <- t.test(
  HeatCap_MCF7_SvsG0G1_DOWN$Specific_heat,
  HeatCap_MCF7_SvsG0G1_UP$Specific_heat,
  paired = FALSE
)
print(TTEST_heat_capacity_MCF7_SvsG0G1_UP_vs_DOWN)

ks_test_MCF7_SvsG0G1_UP_vs_DOWN_heat_capacity <- ks.test(
  HeatCap_MCF7_SvsG0G1_DOWN$Specific_heat,
  HeatCap_MCF7_SvsG0G1_UP$Specific_heat
)
print(ks_test_MCF7_SvsG0G1_UP_vs_DOWN_heat_capacity)

wilcox_test_MCF7_SvsG0G1_UP_vs_DOWN_heat_capacity <- wilcox.test(
  HeatCap_MCF7_SvsG0G1_DOWN$Specific_heat,
  HeatCap_MCF7_SvsG0G1_UP$Specific_heat
)
print(wilcox_test_MCF7_SvsG0G1_UP_vs_DOWN_heat_capacity)












# Heat capacity changes GSE104736 HeLa S vs G1
# Box plot heat capacity for all transcripts UP- and DOWN-regulated

Plot_Heat_capacity_HeLa_S_vs_G1 <- ggplot(MFE_GSE104736_S_vs_G1, aes(x = factor(log2FoldChange > 0), y = Specific_heat)) +
  geom_boxplot(fill = "lightgrey", alpha = 0.7) +
  # geom_jitter(width = 0.1, size = 1, alpha = 0.3) +  # optional
  theme_classic() +
  ylim(0, 400) +  # adjust based on your data range
  theme(
    legend.position = "none", 
    axis.title = element_text(size = 20), 
    axis.text.x = element_text(size = 10, color = "black"),
    plot.title = element_text(size = 20)
  ) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  labs(
    y = "", 
    x = "", 
    title = ""
  ) +
  scale_x_discrete(labels = c("log2FoldChange < 0", "log2FoldChange > 0"))

print(Plot_Heat_capacity_HeLa_S_vs_G1)





# TTEST and KS test of the heat capacity values of all transcripts significantly UP- and DOWN-regulated
HeatCap_HeLa_S_vs_G1_UP <- dplyr::filter(MFE_GSE104736_S_vs_G1, log2FoldChange > 0)
HeatCap_HeLa_S_vs_G1_DOWN <- dplyr::filter(MFE_GSE104736_S_vs_G1, log2FoldChange < 0)

TTEST_heat_capacity_HeLa_S_vs_G1_UP_vs_DOWN <- t.test(
  HeatCap_HeLa_S_vs_G1_DOWN$Specific_heat,
  HeatCap_HeLa_S_vs_G1_UP$Specific_heat,
  paired = FALSE
)
print(TTEST_heat_capacity_HeLa_S_vs_G1_UP_vs_DOWN)

ks_test_HeLa_S_vs_G1_UP_vs_DOWN_heat_capacity <- ks.test(
  HeatCap_HeLa_S_vs_G1_DOWN$Specific_heat,
  HeatCap_HeLa_S_vs_G1_UP$Specific_heat
)
print(ks_test_HeLa_S_vs_G1_UP_vs_DOWN_heat_capacity)

wilcox_test_HeLa_S_vs_G1_UP_vs_DOWN_heat_capacity <- wilcox.test(
  HeatCap_HeLa_S_vs_G1_DOWN$Specific_heat,
  HeatCap_HeLa_S_vs_G1_UP$Specific_heat
)
print(wilcox_test_HeLa_S_vs_G1_UP_vs_DOWN_heat_capacity)






