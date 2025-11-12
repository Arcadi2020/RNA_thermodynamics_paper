# This script generates the plots for Figure 5
#
# Overview: Figure 5 – Cross-validation of transcriptome thermodynamic shifts using GREIN datasets,
#           Science model systems, iLINCS drug-response profiles, and cell-cycle contrasts
#
# High-level workflow
# • Load precomputed RData environments and external tables (GREIN, iLINCS, Science models, cell-cycle gene sets).
# • Construct a wide log2FC matrix: DESeq2 gene-level fold changes from multiple GREIN studies
#   (erlotinib, trametinib, dasatinib, trastuzumab, etc.) and from Science model systems (X778 palbociclib,
#   nutlin-3a, tunicamycin, SKBR3 lapatinib, SKMEL28 vemurafenib; both 72h and early-phase).
# • For selected datasets, compute RNA “heat capacity” (Specific_heat) contrasts between UP-regulated (log2FC > 0)
#   and DOWN-regulated (log2FC < 0) transcripts:
#   – Analyze both all significant transcripts and high-MFE subsets.
#   – Apply multiple tests: two-sample t-test, Kolmogorov–Smirnov, Wilcoxon rank-sum.
# • Assemble per-dataset and combined boxplots (GREIN and Science models) and summarize mean values and p-values.
# • Generate iLINCS analyses:
#   – Volcano plots of Log10_Specific_heat_change vs −log10(Wilcoxon p) at 3.33 µM and 10 µM,
#     with quadrant classification of samples.
#   – Overlay horizontal boxplots and test for global shifts in Specific_heat distribution (Wilcoxon test).
#   – Produce per-sample boxplots for selected significant drug–cell line combinations (e.g. PC3_torin-1, A375_tozasertib, ES2_lapatinib).
# • Benchmark cell-cycle transitions in MCF7 and HeLa (G0/G1 vs G2/M, S vs G0/G1, G1 vs G2, S vs G1)
#   using the same UP vs DOWN heat-capacity contrasts.


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




# Load Table_ANNOTATED_MFE_heat_capacity ####
# The table has been generated with the script:
# path_to/R_scripts/MFE_Human_transcriptome.R
Gencode_TX_norm_MFE_HC_table <- read.delim("path_to/Gencode_TX_norm_MFE_HC_table.csv",
                                           sep = ",", header = T, stringsAsFactors = F)
Table_ANNOTATED_gene_names <- Gencode_TX_norm_MFE_HC_table



load("path_to/Cross_validation_MFE_AC_environment.RData")


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
Expression_changes_X778_PALBO <- read.delim(file="/Volumes/AC_HD_8TB/Exp.\ 271\ Profiling\ of\ alternative\ splicing\ landscape_Nanopore\ sequencing\ drug\ treatment\ 72h/Data_Andre_modified_analysis_26_Sep_2022/72h.files/72h_778_palbociclib.tx_counts.both_conditions.treated_untreated_deseq2.tab", header = TRUE, sep = "\t", quote = "") %>%
  dplyr::rename(Gene_symbol = gene_name) %>%
  dplyr::rename(X778_PALBO_Log_FC = log2FoldChange) %>%
  dplyr::select(Gene_symbol, X778_PALBO_Log_FC) %>%
  group_by(Gene_symbol) %>%
  summarize(X778_PALBO_Log_FC = mean(X778_PALBO_Log_FC)) %>%
  dplyr::filter(!is.na(X778_PALBO_Log_FC))


Expression_changes_SKBR3_LAP <- read.delim(file="/Volumes/AC_HD_8TB/Exp.\ 271\ Profiling\ of\ alternative\ splicing\ landscape_Nanopore\ sequencing\ drug\ treatment\ 72h/Data_Andre_modified_analysis_26_Sep_2022/72h.files/72h_SKBR3_lapatinib.tx_counts.both_conditions.treated_untreated_deseq2.tab", header = TRUE, sep = "\t", quote = "") %>%
  dplyr::rename(Gene_symbol = gene_name) %>%
  dplyr::rename(SKBR3_Log_FC = log2FoldChange) %>%
  dplyr::select(Gene_symbol, SKBR3_Log_FC) %>%
  group_by(Gene_symbol) %>%
  summarize(SKBR3_Log_FC = mean(SKBR3_Log_FC)) %>%
  dplyr::filter(!is.na(SKBR3_Log_FC))


Expression_changes_SKMEL28_VEM <- read.delim(file="/Volumes/AC_HD_8TB/Exp.\ 271\ Profiling\ of\ alternative\ splicing\ landscape_Nanopore\ sequencing\ drug\ treatment\ 72h/Data_Andre_modified_analysis_26_Sep_2022/72h.files/72h_SKMEL28_vemurafenib.tx_counts.both_conditions.treated_untreated_deseq2.tab", header = TRUE, sep = "\t", quote = "") %>%
  dplyr::rename(Gene_symbol = gene_name) %>%
  dplyr::rename(SKMEL28_Log_FC = log2FoldChange) %>%
  dplyr::select(Gene_symbol, SKMEL28_Log_FC) %>%
  group_by(Gene_symbol) %>%
  summarize(SKMEL28_Log_FC = mean(SKMEL28_Log_FC)) %>%
  dplyr::filter(!is.na(SKMEL28_Log_FC))





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
  dplyr::full_join(Expression_changes_SKBR3_LAP, by = "Gene_symbol") %>%
  dplyr::full_join(Expression_changes_SKMEL28_VEM, by = "Gene_symbol")
















# Fig. 5B Heat_capacity GREIN dataset ####
# Heat capacity plots


# Heat capacity changes GSE143462_omipalisib 12h
# GSE143462_omipalisib box plot heat capacity all transcripts UP- and DOWN-regulated
Plot_Heat_capacity_GSE143462_omipalisib_12h <- ggplot(MFE_res_GSE143462, aes(x = factor(log2FoldChange > 0), y = Cp_37C)) +
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
Plot_Heat_capacity_high_GSE143462_omipalisib_12h <- ggplot(GSE143462_Cp_37C_data, aes(x = regulation, y = GSE143462_Cp_37C, fill = regulation)) +
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
Plot_Heat_capacity_GSE121105_trastuzumab_72h <- ggplot(MFE_res_GSE121105_trastuzumab, aes(x = factor(log2FoldChange > 0), y = Cp_37C)) +
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
Plot_Heat_capacity_high_GSE121105_trastuzumab_72h <- ggplot(GSE121105_trastuzumab_Cp_37C_data, aes(x = regulation, y = GSE121105_trastuzumab_Cp_37C, fill = regulation)) +
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
Plot_Heat_capacity_GSE103021_ARS_1620_LU65_24h <- ggplot(MFE_res_GSE103021_ARS_1620_LU65, aes(x = factor(log2FoldChange > 0), y = Cp_37C)) +
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
Plot_Heat_capacity_high_GSE103021_ARS_1620_LU65_24h <- ggplot(GSE103021_ARS_1620_LU65_Cp_37C_data, aes(x = regulation, y = GSE103021_ARS_1620_LU65_Cp_37C, fill = regulation)) +
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
Plot_Heat_capacity_GSE135514_Das_24h <- ggplot(MFE_res_GSE135514_Das_24h, aes(x = factor(log2FoldChange > 0), y = Cp_37C)) +
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
Plot_Heat_capacity_high_GSE135514_Das_24h <- ggplot(GSE135514_Das_24h_Cp_37C_data, aes(x = regulation, y = GSE135514_Das_24h_Cp_37C, fill = regulation)) +
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
Plot_Heat_capacity_GSE135514_Das_72h <- ggplot(MFE_res_GSE135514_Das_72h, aes(x = factor(log2FoldChange > 0), y = Cp_37C)) +
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
Plot_Heat_capacity_high_GSE135514_Das_72h <- ggplot(GSE135514_Das_72h_Cp_37C_data, aes(x = regulation, y = GSE135514_Das_72h_Cp_37C, fill = regulation)) +
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

# Save the combined plot
ggsave("combined_plots.pdf", plot = combined_plots, width = 12, height = 6, path = "path_to/Desktop")












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
  Cp_37C ~ factor(log2FoldChange > 0), 
  data = MFE_res_GSE143462
)

wilcox_heat_capacity_GSE121105_trastuzumab_all_UP_vs_DOWN <- wilcox.test(
  Cp_37C ~ factor(log2FoldChange > 0), 
  data = MFE_res_GSE121105_trastuzumab
)

wilcox_heat_capacity_GSE103021_ARS_1620_LU65_all_UP_vs_DOWN <- wilcox.test(
  Cp_37C ~ factor(log2FoldChange > 0), 
  data = MFE_res_GSE103021_ARS_1620_LU65
)

wilcox_heat_capacity_GSE135514_Das_24h_all_UP_vs_DOWN <- wilcox.test(
  Cp_37C ~ factor(log2FoldChange > 0), 
  data = MFE_res_GSE135514_Das_24h
)

wilcox_heat_capacity_GSE135514_Das_72h_all_UP_vs_DOWN <- wilcox.test(
  Cp_37C ~ factor(log2FoldChange > 0), 
  data = MFE_res_GSE135514_Das_72h
)


# Wilcoxon test for high MFE transcripts (UP vs DOWN)
wilcox_heat_capacity_GSE143462_high_MFE_TX_UP_vs_DOWN <- wilcox.test(
  GSE143462_Cp_37C ~ regulation, 
  data = GSE143462_Cp_37C_data
)

wilcox_heat_capacity_GSE121105_trastuzumab_high_MFE_TX_UP_vs_DOWN <- wilcox.test(
  GSE121105_trastuzumab_Cp_37C ~ regulation, 
  data = GSE121105_trastuzumab_Cp_37C_data
)

wilcox_heat_capacity_GSE103021_ARS_1620_LU65_high_MFE_TX_UP_vs_DOWN <- wilcox.test(
  GSE103021_ARS_1620_LU65_Cp_37C ~ regulation, 
  data = GSE103021_ARS_1620_LU65_Cp_37C_data
)

wilcox_heat_capacity_GSE135514_Das_24h_high_MFE_TX_UP_vs_DOWN <- wilcox.test(
  GSE135514_Das_24h_Cp_37C ~ regulation, 
  data = GSE135514_Das_24h_Cp_37C_data
)

wilcox_heat_capacity_GSE135514_Das_72h_high_MFE_TX_UP_vs_DOWN <- wilcox.test(
  GSE135514_Das_72h_Cp_37C ~ regulation, 
  data = GSE135514_Das_72h_Cp_37C_data
)


# Add Wilcoxon p-values to the summary table
TTEST_KS_results$Wilcoxon_p_value <- c(
  wilcox_heat_capacity_GSE143462_all_UP_vs_DOWN$p.value,
  wilcox_heat_capacity_GSE121105_trastuzumab_all_UP_vs_DOWN$p.value,
  wilcox_heat_capacity_GSE103021_ARS_1620_LU65_all_UP_vs_DOWN$p.value,
  wilcox_heat_capacity_GSE135514_Das_24h_all_UP_vs_DOWN$p.value,
  wilcox_heat_capacity_GSE135514_Das_72h_all_UP_vs_DOWN$p.value
)







# Combine all plots heat capacity high MFE transcripts

library(patchwork)

combined_plot <- (Plot_Heat_capacity_GSE143462_omipalisib_12h +
                    Plot_Heat_capacity_GSE121105_trastuzumab_72h +
                    Plot_Heat_capacity_GSE103021_ARS_1620_LU65_24h +
                    Plot_Heat_capacity_GSE135514_Das_24h +
                    Plot_Heat_capacity_GSE135514_Das_72h) +
  plot_layout(ncol = 5, guides = "collect") &  # 5 columns, collect legends
  theme(
    axis.text.x = element_blank(),    
    axis.ticks.x = element_blank(),   
    axis.title.x = element_blank(),   
    axis.title.y = element_blank(),   
    plot.title = element_blank()      
  )

ggsave("GREIN_combined_plots.pdf", combined_plot, width = 12, height = 6, path = "path_to/Desktop")























# Fig. 5C Heat_capacity ####

# Analysis heat capacity Science model systems

# Analysis 778 palbociclib
# Analysis heat capacity 778 palbociclib 72h all UP-regulated vs DOWN-regulated transcripts
Plot_heat_capacity_778_PALBO <- ggplot(X778_PALBO_dexseq_ANNOTATED, aes(x = factor(log2fold_treated_untreated > 0), y = Cp_37C)) +
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


# Save the plot
ggsave("Plot_heat_capacity_778_PALBO.pdf", plot = Plot_heat_capacity_778_PALBO, width = 4, height = 6, path = "path_to/Desktop")




# Wilcoxon test of the heat capacity values of all transcripts UP- vs DOWN-regulated
wilcox_heat_capacity_X778_PALBO_dexseq_ANNOTATED_all_UP_vs_DOWN <- wilcox.test(
  Cp_37C ~ factor(log2fold_treated_untreated > 0),
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















# Fig. 5D volcano plot iLINCS 3.33 µM ####

# Load the environment of the script output
load("path_to/Analysis_MFE_iLINCS_all_datasets_3.33uM_cancer_lines_environment.RData")

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
all_results_2_complete_3.33uM <- read.csv("path_to/Results_all_cancer_cell_lines_3.33uM.csv",
                                   header = TRUE, sep = ",", stringsAsFactors = FALSE) %>%
  dplyr::rename(iLINCs_cell_line = cell_line)

iLINCS_analysis_3.33uM <- all_results_2_complete_3.33uM



library(ggrepel)

iLINCS_analysis_plot_membrane_receptors <- ggplot(
  iLINCS_analysis_3.33uM,
  aes(x = Log10_Cp_37C_change,
      y = -log10(wilcox_test_all_UP_vs_DOWN_Cp_37C))
) +
  theme_classic() +
  geom_point(
    aes(
      fill = ifelse(
        Log10_Cp_37C_change > 0 &
          -log10(wilcox_test_all_UP_vs_DOWN_Cp_37C) > 1.301,
        "blue", "lightgrey"
      )
    ),
    color = "black",
    size = 4,
    shape = 21,
    stroke = 0.6,
    alpha = 1   # <- constant alpha so all greys look the same
  ) +
  scale_fill_identity() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 1.301, linetype = "dashed", color = "black") +
  labs(
    x = "Log10_Cp_37C_change",
    y = "-Log10_wilcox_test_all_UP_vs_DOWN_pValue",
    title = ""
  ) +
  theme(
    plot.title = element_text(color = "black", size = 14, face = "bold"),
    legend.position = "none"
  ) +
  xlim(-0.2, 0.2) +
  geom_text_repel(
    aes(label = ifelse(
      -log10(wilcox_test_all_UP_vs_DOWN_Cp_37C) > 1.301 &
        paste(iLINCs_cell_line, Drug, sep = "_") %in%
        c("A375_MLN-2480", "PC3_PI-103", "PC3_torin-1"),
      paste(iLINCs_cell_line, Drug, sep = "_"),
      ""
    )),
    size = 3,
    box.padding = 0.5,
    point.padding = 0.5,
    nudge_x = 0.004,
    nudge_y = 0.001,
    max.overlaps = Inf,
    segment.color = "lightgrey",
    segment.size = 0.5
  )

print(iLINCS_analysis_plot_membrane_receptors)














# Classify samples into heat capacity quadrants
library(dplyr)

# --- Parameters
alpha_thr <- 0.05
y_thr <- -log10(alpha_thr)   # = 1.30103

# --- Clean data (drop NAs just in case)
iLINCS_analysis_3.33uM_clean <- iLINCS_analysis_3.33uM %>%
  filter(is.finite(Log10_Cp_37C_change),
         is.finite(wilcox_test_all_UP_vs_DOWN_Cp_37C))

# --- Classify quadrants (include boundary lines)
iLINCS_analysis_3.33uM_clean <- iLINCS_analysis_3.33uM_clean %>%
  mutate(
    HeatCapacity_Quadrant = case_when(
      Log10_Cp_37C_change >= 0 & -log10(wilcox_test_all_UP_vs_DOWN_Cp_37C) >= y_thr ~ "UP_regulated_heat_capacity",
      Log10_Cp_37C_change >= 0 & -log10(wilcox_test_all_UP_vs_DOWN_Cp_37C) <  y_thr ~ "Q2",
      Log10_Cp_37C_change  < 0 & -log10(wilcox_test_all_UP_vs_DOWN_Cp_37C) >= y_thr ~ "DOWN_regulated_heat_capacity",
      Log10_Cp_37C_change  < 0 & -log10(wilcox_test_all_UP_vs_DOWN_Cp_37C) <  y_thr ~ "Q4",
      TRUE ~ "Other"
    )
  )

# --- Overall counts and percentages by quadrant
quad_table <- iLINCS_analysis_3.33uM_clean %>%
  count(HeatCapacity_Quadrant, name = "N") %>%
  mutate(Percent = round(100 * N / sum(N), 1)) %>%
  arrange(factor(HeatCapacity_Quadrant,
                 levels = c("UP_regulated_heat_capacity","Q2",
                            "DOWN_regulated_heat_capacity","Q4","Other")))

# --- Side totals (RIGHT = UP; LEFT = DOWN)
up_total   <- iLINCS_analysis_3.33uM_clean %>% filter(Log10_Cp_37C_change >= 0) %>% nrow()
down_total <- iLINCS_analysis_3.33uM_clean %>% filter(Log10_Cp_37C_change  < 0) %>% nrow()
N_all      <- nrow(iLINCS_analysis_3.33uM_clean)

# --- Significant subsets (Wilcoxon)
sig_up   <- iLINCS_analysis_3.33uM_clean %>% filter(HeatCapacity_Quadrant == "UP_regulated_heat_capacity") %>% nrow()
sig_down <- iLINCS_analysis_3.33uM_clean %>% filter(HeatCapacity_Quadrant == "DOWN_regulated_heat_capacity") %>% nrow()

side_summary <- tibble::tibble(
  Side = c("UP", "DOWN"),
  Total_N = c(up_total, down_total),
  Total_Percent = round(100 * c(up_total, down_total) / N_all, 1),
  Sig_N = c(sig_up, sig_down),
  Sig_Percent_Overall = round(100 * c(sig_up, sig_down) / N_all, 1),
  Sig_Percent_WithinSide = round(100 * c(sig_up / up_total, sig_down / down_total), 1)
)

# --- Outputs (data frames)
quad_table
side_summary












# Volcano plot
volcano_plot_3.33uM <- ggplot(iLINCS_analysis_3.33uM, aes(x = Log10_Cp_37C_change, y = -log10(wilcox_test_all_UP_vs_DOWN_Cp_37C))) +
  theme_classic() +
  geom_point(aes(
    fill = case_when(
      Log10_Cp_37C_change > 0 & -log10(wilcox_test_all_UP_vs_DOWN_Cp_37C) > 1.301 ~ "blue",
      TRUE ~ "lightgrey"
    ),
    #alpha = ifelse(-log10(wilcox_test_all_UP_vs_DOWN_Cp_37C) < 1.301, 0.6, 1)
  ),
  shape = 21, size = 4, stroke = 0.6, color = "black") +
  scale_fill_identity() +
  scale_alpha_identity() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 1.301, linetype = "dashed", color = "black") +
  
  # Uncomment the geom_text_repel to remove the labels
  # # Label only selected samples
  geom_text_repel(
    aes(label = ifelse(
      -log10(wilcox_test_all_UP_vs_DOWN_Cp_37C) > 1.301 &
        paste(iLINCs_cell_line, Drug, sep = "_") %in%
        c("A375_MLN-2480", "PC3_PI-103", "PC3_torin-1"),
      paste(iLINCs_cell_line, Drug, sep = "_"), "")
    ),
    size = 3,
    box.padding = 0.5,
    point.padding = 0.5,
    nudge_x = 0.004,
    nudge_y = 0.001,
    max.overlaps = Inf,
    segment.color = 'lightgrey',
    segment.size = 0.5
  ) +
  labs(
    x = "Log10_Cp_37C_change",
    y = "-Log10_wilcox_test_all_UP_vs_DOWN_pValue",
    title = ""
  ) +
  xlim(-0.2, 0.2) +
  ylim(0, 10) +
  theme(
    plot.title = element_text(color = "black", size = 14, face = "bold"),
    legend.position = "none"
  )

print(volcano_plot_3.33uM)



# Horizontal boxplot for Log10_Cp_37C_change
boxplot_heat_capacity <- ggplot(iLINCS_analysis_3.33uM, aes(y = "", x = Log10_Cp_37C_change)) +
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
combined_plot <- boxplot_heat_capacity / volcano_plot_3.33uM + plot_layout(heights = c(0.07, 0.93))


print(combined_plot)

ggsave("combined_heat_capacity_volcano_plot.tiff",
       plot = combined_plot, width = 5, height = 8.5,
       path = "path_to/Desktop", dpi = 300)


Heat_capacity_shift_significance_3.33uM <- wilcox.test(iLINCS_analysis_3.33uM$Log10_Cp_37C_change,
            mu = 0, alternative = "greater")
print(Heat_capacity_shift_significance_3.33uM)





















# Fig. 5E Box plots iLINCS 3.33 µM heat capacity ####

library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(patchwork)



# Load inputs
directory <- "path_to/Output_script_cancer_cell_lines_3.33uM"
result_list_2 <- readRDS(file.path(directory, "iLINCs_filtered_samples.rds"))

Table_ANNOTATED_iLINCs <- read.delim(
  "path_to/Table_ANNOTATED_iLINCs.csv",
  sep = ",", header = TRUE, stringsAsFactors = FALSE
) %>% dplyr::rename(ensembl_id = ensembl_gene_id)

# Mapping (REP-aware)

# Helper: extract canonical key "CELL_DRUG" from a list element name
canon_from_listname <- function(x) {
  parts <- strsplit(x, "_", fixed = TRUE)[[1]]
  paste(parts[1:min(2, length(parts))], collapse = "_")
}

# Prefer a candidate that starts with "<key>_REP"; fallback to shortest matching name
prefer_rep_pick <- function(key, candidates) {
  rep_cand <- grep(paste0("^", key, "_REP"), candidates, value = TRUE)
  pick_set <- if (length(rep_cand) > 0) rep_cand else candidates
  pick_set[order(nchar(pick_set))][1]
}

pool_names <- names(result_list_2)
all_keys <- unique(vapply(pool_names, canon_from_listname, character(1)))

# Build lookup: canonical_key -> chosen listname (REP preferred)
lookup <- tibble::tibble(canonical_key = all_keys) %>%
  rowwise() %>%
  mutate(listname = {
    cand <- grep(paste0("^", canonical_key, "($|_)"), pool_names, value = TRUE)
    if (length(cand) == 0) NA_character_ else prefer_rep_pick(canonical_key, cand)
  }) %>%
  ungroup()

# Requested exemplars for 3.33 µM
samples_of_interest_3.33uM <- c("A375_MLN-2480", "PC3_PI-103", "PC3_torin-1")

# Normalize requested keys to canonical form "CELL_DRUG"
norm_key <- function(key) {
  parts <- strsplit(key, "_", fixed = TRUE)[[1]]
  paste(parts[1:min(2, length(parts))], collapse = "_")
}

requested_canon <- tibble::tibble(
  requested = samples_of_interest_3.33uM,
  canonical_key = vapply(samples_of_interest_3.33uM, norm_key, character(1))
)

# Resolve requested to list names via REP-aware lookup (and fallback if needed)
resolved <- requested_canon %>% left_join(lookup, by = "canonical_key")

if (any(is.na(resolved$listname))) {
  for (i in which(is.na(resolved$listname))) {
    key <- resolved$canonical_key[i]
    cand <- grep(paste0("^", key, "($|_)"), pool_names, value = TRUE)
    if (length(cand) > 0) resolved$listname[i] <- prefer_rep_pick(key, cand)
  }
}

if (all(is.na(resolved$listname))) {
  stop("No result_list_2 entries matched. Inspect names(result_list_2)[1:30] and refine sample keys.")
}

# Per-sample boxplots and summary table
plot_list_3.33uM <- list()
all_results_3.33uM <- data.frame(
  SampleName = character(),
  MeanCp_37C_UP = numeric(),
  MeanCp_37C_DOWN = numeric(),
  ks_test_all_UP_vs_DOWN_Cp_37C = numeric(),
  wilcox_test_all_UP_vs_DOWN_Cp_37C = numeric(),
  stringsAsFactors = FALSE
)

# Build a canonical panel title "CELL_DRUG" from the list element name
title_from_listname <- function(x) {
  parts <- strsplit(x, "_", fixed = TRUE)[[1]]
  paste(parts[1:min(2, length(parts))], collapse = "_")
}

plotted <- unique(stats::na.omit(resolved$listname))

for (ln in plotted) {
  df <- result_list_2[[ln]]
  if (is.null(df)) next
  
  title_plot <- title_from_listname(ln)
  if (title_plot %in% all_results_3.33uM$SampleName) next
  
  # Join with annotation to get Cp_37C
  Sample_TX <- dplyr::right_join(df, Table_ANNOTATED_iLINCs, by = "ensembl_id")
  
  # Keep only significant z-scores
  Sample_TX_sig <- dplyr::filter(Sample_TX, p_value < 0.05)
  
  # Split by direction and compute means
  df_up   <- dplyr::filter(Sample_TX_sig, z_score > 0)
  df_down <- dplyr::filter(Sample_TX_sig, z_score < 0)
  
  Mean_UP   <- mean(df_up$Cp_37C,   na.rm = TRUE)
  Mean_DOWN <- mean(df_down$Cp_37C, na.rm = TRUE)
  
  # Stats (KS + Wilcoxon)
  KS    <- tryCatch(ks.test(df_down$Cp_37C, df_up$Cp_37C), error = function(e) NA)
  WILCX <- tryCatch(wilcox.test(df_down$Cp_37C, df_up$Cp_37C), error = function(e) NA)
  
  # Plot with y-axis fixed to 0..300
  p <- ggplot(
    Sample_TX_sig,
    aes(
      x = factor(z_score > 0,
                 levels = c(FALSE, TRUE),
                 labels = c("z_score < 0", "z_score > 0")),
      y = Cp_37C
    )
  ) +
    geom_boxplot(fill = "lightgrey", alpha = 0.7, outlier.shape = 16, outlier.size = 0.8) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.title = element_text(size = 20),
      axis.text.x = element_text(size = 12, color = "black"),
      plot.title = element_text(size = 12, face = "bold")
    ) +
    geom_vline(xintercept = 1.5, linewidth = 0.3) +
    labs(y = "", x = NULL, title = title_plot) +
    scale_y_continuous(limits = c(0, 300))
  
  # Save each panel (optional)
  ggsave(filename = paste0(title_plot, "_Cp_37C_plot_3.33uM.pdf"),
         plot = p, width = 4, height = 5, path = "path_to/Desktop")
  
  # Store and append summary
  plot_list_3.33uM[[title_plot]] <- p
  
  all_results_3.33uM <- rbind(
    all_results_3.33uM,
    data.frame(
      SampleName = title_plot,
      MeanCp_37C_UP = Mean_UP,
      MeanCp_37C_DOWN = Mean_DOWN,
      ks_test_all_UP_vs_DOWN_Cp_37C = if (is.list(KS)) KS$p.value else NA_real_,
      wilcox_test_all_UP_vs_DOWN_Cp_37C = if (is.list(WILCX)) WILCX$p.value else NA_real_,
      stringsAsFactors = FALSE
    )
  )
}

# Ensure unique rows
all_results_3.33uM <- dplyr::distinct(all_results_3.33uM, SampleName, .keep_all = TRUE)
print(all_results_3.33uM)





















# Fig. 5F volcano plot iLINCS 10 µM ####

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

iLINCS_analysis <- all_results_2_complete



# Add information on cell lines and drugs and group the samples accordingly

# Open the iLINCs cell_info file
iLINCs_cell_info <- read.delim("path_to/cellinfo_beta.csv",
                               sep = ",", header = T, stringsAsFactors = F) %>%
  dplyr::rename(cell_line = cell_iname)



# Merge the iLINCS_analysis data frame with iLINCs_cell_info on the cell_line and cell_iname columns
iLINCs_cell_info_complete <- dplyr::left_join(iLINCS_analysis, iLINCs_cell_info, by = c("cell_line", "primary_disease")) %>%
  dplyr::select(cell_line, Drug, primary_disease) %>%
  dplyr::mutate(primary_disease = ifelse(cell_line == "P1A82", "teratocarcinoma", primary_disease))




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
  dplyr::distinct(cell_line, Drug, primary_disease, TARGET, TARGET_PATHWAY)




# Plot the data
# Create a new column with primary_disease for grouping iLINCs_cell_lines
IlINCs_complete_df <- IlINCs_complete_df %>%
  mutate(primary_disease_label = factor(primary_disease, levels = unique(primary_disease)),
         TARGET_PATHWAY_label = factor(TARGET_PATHWAY, levels = unique(TARGET_PATHWAY)))


# Total n. of cell lines:
Total_cell_lines <- unique(IlINCs_complete_df$cell_line)
print(length(Total_cell_lines))
# Total n. of anti-cancer agents:
Total_anti_cancer_agents <- unique(IlINCs_complete_df$Drug)
print(length(Total_anti_cancer_agents))











iLINCS_analysis_plot_membrane_receptors <- ggplot(iLINCS_analysis, aes(x = Log10_Cp_37C_change, y = -log10(iLINCS_analysis$wilcox_test_all_UP_vs_DOWN_Cp_37C))) +
  theme_classic() +
  geom_point(aes(
    fill = ifelse(Log10_Cp_37C_change > 0 & -log10(iLINCS_analysis$wilcox_test_all_UP_vs_DOWN_Cp_37C) > 1.301, "blue", 
                  ifelse(Log10_Cp_37C_change > 0 & -log10(iLINCS_analysis$wilcox_test_all_UP_vs_DOWN_Cp_37C) < 1.301, "lightgrey",
                         ifelse(Log10_Cp_37C_change < 0 & -log10(iLINCS_analysis$wilcox_test_all_UP_vs_DOWN_Cp_37C) > 1.301, "lightgrey", 
                                ifelse(Log10_Cp_37C_change < 0 & -log10(iLINCS_analysis$wilcox_test_all_UP_vs_DOWN_Cp_37C) < 1.301, "lightgrey", "lightgrey"))))),
    color = "black",  # Border color for all points
    size = 4,  # Set a larger fixed size for all points
    shape = 21,  # Set a consistent shape for all points
    alpha = ifelse(-log10(iLINCS_analysis$wilcox_test_all_UP_vs_DOWN_Cp_37C) < 1.301, 0.6, 1),
    stroke = 0.6) +  # Stroke for empty circles' border
  scale_fill_identity() +
  scale_color_identity() +
  scale_alpha_identity() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 1.301, linetype = "dashed", color = "black") +
  labs(x = "Log10_Cp_37C_change", y = "-Log10_ks_test_all_UP_vs_DOWN_pValue", title = "") +
  theme(plot.title = element_text(color = "black", size = 14, face = "bold"),
        legend.position = "none") +
  xlim(-0.2, 0.2) +
  geom_text_repel(
    aes(label = ifelse(
      Log10_Cp_37C_change > 0 &
        -log10(wilcox_test_all_UP_vs_DOWN_Cp_37C) > 1.301,
      paste(cell_line, Drug, sep = "_"),
      NA_character_  # use NA so unlabeled points are skipped cleanly
    )),
    size = 3, box.padding = 0.5, point.padding = 0.5,
    nudge_x = 0.004, nudge_y = 0.001, max.overlaps = Inf,
    segment.color = "lightgrey", segment.size = 0.5, na.rm = TRUE
  )


print(iLINCS_analysis_plot_membrane_receptors)












# Add the name of the specific cell lines and pharmacologic agents to the plot
iLINCS_analysis_plot_membrane_receptors <- ggplot(iLINCS_analysis, aes(x = Log10_Cp_37C_change, y = -log10(wilcox_test_all_UP_vs_DOWN_Cp_37C))) +
  theme_classic() +
  geom_point(
    aes(
      fill = ifelse(
        Log10_Cp_37C_change > 0 & -log10(wilcox_test_all_UP_vs_DOWN_Cp_37C) > 1.301, "blue", 
        "lightgrey"
      ),
      #alpha = ifelse(-log10(wilcox_test_all_UP_vs_DOWN_Cp_37C) < 1.301, 0.6, 1)
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
  labs(x = "Log10_Cp_37C_change", y = "-Log10_ks_test_all_UP_vs_DOWN_pValue", title = "") +
  theme(
    plot.title = element_text(color = "black", size = 14, face = "bold"),
    legend.position = "none"
  ) +
  xlim(-0.2, 0.2) +
  geom_text_repel(
    aes(label = ifelse(
      Log10_Cp_37C_change > 0 &
        -log10(wilcox_test_all_UP_vs_DOWN_Cp_37C) > 1.301 &
        paste(cell_line, Drug, sep = "_") %in% c("ES2_lapatinib", "A375_tozasertib", "PC3_torin-1"),
      paste(cell_line, Drug, sep = "_"),
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
library(dplyr)

alpha_thr <- 0.05
y_thr <- -log10(alpha_thr)   # 1.30103

# Make sure it's a tibble and drop NAs IN-PLACE
iLINCS_analysis <- iLINCS_analysis %>%
  as_tibble() %>%
  filter(is.finite(Log10_Cp_37C_change),
         is.finite(wilcox_test_all_UP_vs_DOWN_Cp_37C))

# Classify
iLINCS_analysis <- iLINCS_analysis %>%
  mutate(
    HeatCapacity_Quadrant = case_when(
      Log10_Cp_37C_change >= 0 & -log10(wilcox_test_all_UP_vs_DOWN_Cp_37C) >= y_thr ~ "UP_regulated_heat_capacity",
      Log10_Cp_37C_change >= 0 & -log10(wilcox_test_all_UP_vs_DOWN_Cp_37C) <  y_thr ~ "Q2",
      Log10_Cp_37C_change  < 0 & -log10(wilcox_test_all_UP_vs_DOWN_Cp_37C) >= y_thr ~ "DOWN_regulated_heat_capacity",
      Log10_Cp_37C_change  < 0 & -log10(wilcox_test_all_UP_vs_DOWN_Cp_37C) <  y_thr ~ "Q4",
      TRUE ~ "Other"
    )
  )

# Counts/percentages
quad_table <- iLINCS_analysis %>%
  dplyr::count(HeatCapacity_Quadrant, name = "N") %>%
  mutate(Percent = round(100 * N / sum(N), 1)) %>%
  arrange(factor(HeatCapacity_Quadrant,
                 levels = c("UP_regulated_heat_capacity","Q2",
                            "DOWN_regulated_heat_capacity","Q4","Other")))

# Side totals
up_total   <- sum(iLINCS_analysis$Log10_Cp_37C_change >= 0, na.rm = TRUE)
down_total <- sum(iLINCS_analysis$Log10_Cp_37C_change  < 0, na.rm = TRUE)
N_all      <- nrow(iLINCS_analysis)

sig_up   <- sum(iLINCS_analysis$HeatCapacity_Quadrant == "UP_regulated_heat_capacity", na.rm = TRUE)
sig_down <- sum(iLINCS_analysis$HeatCapacity_Quadrant == "DOWN_regulated_heat_capacity", na.rm = TRUE)

side_summary <- tibble::tibble(
  Side = c("UP", "DOWN"),
  Total_N = c(up_total, down_total),
  Total_Percent = round(100 * c(up_total, down_total) / N_all, 1),
  Sig_N = c(sig_up, sig_down),
  Sig_Percent_Overall = round(100 * c(sig_up, sig_down) / N_all, 1),
  Sig_Percent_WithinSide = round(100 * c(ifelse(up_total==0, NA, sig_up / up_total),
                                         ifelse(down_total==0, NA, sig_down / down_total)), 1)
)

quad_table
side_summary










# ---- Volcano plot ----
volcano_plot <- ggplot(iLINCS_analysis, aes(x = Log10_Cp_37C_change, y = -log10(wilcox_test_all_UP_vs_DOWN_Cp_37C))) +
  theme_classic() +
  geom_point(aes(
    fill = case_when(
      Log10_Cp_37C_change > 0 & -log10(wilcox_test_all_UP_vs_DOWN_Cp_37C) > 1.301 ~ "blue",
      TRUE ~ "lightgrey"
    ),
    #alpha = ifelse(-log10(wilcox_test_all_UP_vs_DOWN_Cp_37C) < 1.301, 0.6, 1)
  ),
  shape = 21, size = 4, stroke = 0.6, color = "black") +
  scale_fill_identity() +
  scale_alpha_identity() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 1.301, linetype = "dashed", color = "black") +
  
  
  labs(
    x = "Log10_Cp_37C_change",
    y = "-Log10_ks_test_all_UP_vs_DOWN_pValue",
    title = ""
  ) +
  xlim(-0.2, 0.2) +
  theme(
    plot.title = element_text(color = "black", size = 14, face = "bold"),
    legend.position = "none"
  )


# Horizontal boxplot for Log10_Cp_37C_change
boxplot_heat_capacity <- ggplot(iLINCS_analysis, aes(y = "", x = Log10_Cp_37C_change)) +
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

ggsave("combined_heat_capacity_volcano_plot.tiff",
       plot = combined_plot, width = 5, height = 8.5,
       path = "path_to/Desktop", dpi = 300)



Heat_capacity_shift_significance_10uM <- wilcox.test(iLINCS_analysis$Log10_Cp_37C_change,
                                                mu = 0, alternative = "greater")
print(Heat_capacity_shift_significance_10uM)


















# Fig. 5G Box plots iLINCS 10 µM heat capacity ####


# Box plots heat capacity

# Edit the names of the samples to match the iLINCs format
all_results_2_complete$iLINCs_cell_line_drug <- sub("^(([^_]*_[^_]*)_).*", "\\1", all_results_2_complete$SampleName)

# Remove the trailing underscore
all_results_2_complete$iLINCs_cell_line_drug <- gsub("_$", "", all_results_2_complete$iLINCs_cell_line_drug)

iLINCs_all_results_2_complete_filtered_Cp_37C <- dplyr::filter(all_results_2_complete,  iLINCs_cell_line_drug %in% c("PC3_torin-1", "A375_tozasertib", "ES2_lapatinib")) %>%
  #dplyr::filter(Log10_Cp_37C_change > 0) %>%
  dplyr::filter(ks_test_all_UP_vs_DOWN_Cp_37C < 0.05)  # Filter for significant samples


# Filter the result_list_2 list object to only include the selected samples
iLINCs_selected_samples_Cp_37C <- result_list_2[grep(paste(iLINCs_all_results_2_complete_filtered_Cp_37C$SampleName, collapse = "|"), names(result_list_2))]



# Open the Table_ANNOTATED_iLINCs
Table_ANNOTATED_iLINCs <- read.delim("path_to/Table_ANNOTATED_iLINCs.csv",
                                     sep = ",", header = T, stringsAsFactors = F)
Table_ANNOTATED_iLINCs <- dplyr::rename(Table_ANNOTATED_iLINCs, ensembl_id = ensembl_gene_id)




# Generate the data of the Cp_37C values
# Initialize plot_list as an empty list
plot_list <- list()

# Initialize all_results as an empty data frame with appropriate column names
all_results <- data.frame(SampleName = character(),
                          MeanCp_37C_UP = numeric(),
                          MeanCp_37C_DOWN = numeric(),
                          ks_test_all_UP_vs_DOWN_Cp_37C = numeric(),
                          stringsAsFactors = FALSE)

# Function to analyze the data and store plots
process_result_list <- function(Sample_name) { 
  message("Processing ", Sample_name)
  df <- iLINCs_selected_samples_Cp_37C[[Sample_name]] # Access individual data frames
  
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
  Cp_37C_Sample_name_UP <- dplyr::filter(Sample_name_TX_filtered, z_score > 0)
  Mean_Cp_37C_Sample_name_UP <- mean(Cp_37C_Sample_name_UP$Cp_37C, na.rm = TRUE)
  
  Cp_37C_Sample_name_DOWN <- dplyr::filter(Sample_name_TX_filtered, z_score < 0)
  Mean_Cp_37C_Sample_name_DOWN <- mean(Cp_37C_Sample_name_DOWN$Cp_37C, na.rm = TRUE)
  
  # Perform statistical tests
  TTEST <- t.test(Cp_37C_Sample_name_DOWN$Cp_37C, Cp_37C_Sample_name_UP$Cp_37C, paired = FALSE)
  
  ks_test_all_UP_vs_DOWN <- ks.test(Cp_37C_Sample_name_DOWN$Cp_37C, Cp_37C_Sample_name_UP$Cp_37C)
  
  wilcox_test_all_UP_vs_DOWN <- wilcox.test(Cp_37C_Sample_name_DOWN$Cp_37C, 
                                            Cp_37C_Sample_name_UP$Cp_37C)
  
  # Calculate the net Cp_37C change
  Sample_name_TX_net_Cp_37C_change <- Sample_name_TX_filtered %>% 
    dplyr::mutate(Net_Cp_37C_change = (Cp_37C * z_score) * (-1))
  Net_Cp_37C_change <- sum(Sample_name_TX_net_Cp_37C_change$Net_Cp_37C_change, na.rm = TRUE)
  
  ## Generate box plot 
  Cp_37C_Sample_name <- ggplot(Sample_name_TX_filtered, aes(x = factor(z_score > 0), y = Cp_37C)) +
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
    filename = paste0(Sample_name, "_Cp_37C_plot.pdf"),  # Create a unique filename
    plot = Cp_37C_Sample_name,
    width = 4,
    height = 5,
    path = "path_to/Desktop"  # Specify the path to save the plot
  )
  
  # Store the plot in the list
  plot_list[[Sample_name]] <<- Cp_37C_Sample_name
  
  # Create a data frame with the results for the current sample
  sample_results <- data.frame(
    SampleName = Sample_name,
    MeanCp_37C_UP = Mean_Cp_37C_Sample_name_UP,
    MeanCp_37C_DOWN = Mean_Cp_37C_Sample_name_DOWN,
    ks_test_all_UP_vs_DOWN_Cp_37C = ks_test_all_UP_vs_DOWN$p.value,
    wilcox_test_all_UP_vs_DOWN_Cp_37C = wilcox_test_all_UP_vs_DOWN$p.value  # Add Wilcoxon p-value
  )
  
  
  # Append the sample results to the all_results data frame
  all_results <<- rbind(all_results, sample_results)
}

# Run the function for each iLINCS data frame
Data_result_list <- lapply(names(iLINCs_selected_samples_Cp_37C), process_result_list)

# Convert to data frame to clearly see the data
# Bind rows and include a source column
Data_combined <- bind_rows(Data_result_list, .id = "Source") %>%
  dplyr::distinct(SampleName, .keep_all = TRUE)
print(Data_combined)
























# Fig. 5H box plots heat capacity cell cycle ####

# To load the environment
load("path_to/Analysis_cell_cycle_genes_MFE_environment.RData")


# Heat capacity changes MCF7 G0G1 vs G2M
# MCF7 box plot heat capacity all transcripts UP- and DOWN-regulated


Plot_Heat_capacity_MCF7_G0G1_vs_G2M <- ggplot(GSE94479_G0G1vsG2M, aes(x = factor(log2FoldChange > 0), y = Cp_37C)) +
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

# Save plot
ggsave("Plot_Heat_capacity_MCF7_G0G1_vs_G2M.pdf", plot = Plot_Heat_capacity_MCF7_G0G1_vs_G2M, width = 4, height = 6, path = "path_to/Desktop")



# TTEST and KS test of the heat capacity values of all transcripts significantly UP- and DOWN-regulated
HeatCap_MCF7_G0G1vsG2M_UP <- dplyr::filter(GSE94479_G0G1vsG2M, log2FoldChange > 0)
HeatCap_MCF7_G0G1vsG2M_DOWN <- dplyr::filter(GSE94479_G0G1vsG2M, log2FoldChange < 0)

TTEST_heat_capacity_MCF7_all_UP_vs_DOWN <- t.test(
  HeatCap_MCF7_G0G1vsG2M_DOWN$Cp_37C,
  HeatCap_MCF7_G0G1vsG2M_UP$Cp_37C,
  paired = FALSE
)
print(TTEST_heat_capacity_MCF7_all_UP_vs_DOWN)

ks_test_MCF7_all_UP_vs_DOWN_heat_capacity <- ks.test(
  HeatCap_MCF7_G0G1vsG2M_DOWN$Cp_37C,
  HeatCap_MCF7_G0G1vsG2M_UP$Cp_37C
)
print(ks_test_MCF7_all_UP_vs_DOWN_heat_capacity)

wilcox_test_MCF7_all_UP_vs_DOWN_heat_capacity <- wilcox.test(
  HeatCap_MCF7_G0G1vsG2M_DOWN$Cp_37C,
  HeatCap_MCF7_G0G1vsG2M_UP$Cp_37C
)
print(wilcox_test_MCF7_all_UP_vs_DOWN_heat_capacity)










# Heat capacity changes GSE104736 HeLa G1 vs G2
# Box plot heat capacity for all transcripts UP- and DOWN-regulated

Plot_Heat_capacity_HeLa_G1_vs_G2 <- ggplot(MFE_GSE104736_G1_vs_G2, aes(x = factor(log2FoldChange > 0), y = Cp_37C)) +
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

# Save the plot
ggsave("Plot_Heat_capacity_HeLa_G1_vs_G2.pdf", plot = Plot_Heat_capacity_HeLa_G1_vs_G2, width = 4, height = 6, path = "path_to/Desktop")



# TTEST and KS test of the heat capacity values of all transcripts significantly UP- and DOWN-regulated
HeatCap_HeLa_G1_vs_G2_UP <- dplyr::filter(MFE_GSE104736_G1_vs_G2, log2FoldChange > 0)
HeatCap_HeLa_G1_vs_G2_DOWN <- dplyr::filter(MFE_GSE104736_G1_vs_G2, log2FoldChange < 0)

TTEST_heat_capacity_HeLa_G1_vs_G2_UP_vs_DOWN <- t.test(
  HeatCap_HeLa_G1_vs_G2_DOWN$Cp_37C,
  HeatCap_HeLa_G1_vs_G2_UP$Cp_37C,
  paired = FALSE
)
print(TTEST_heat_capacity_HeLa_G1_vs_G2_UP_vs_DOWN)

ks_test_HeLa_G1_vs_G2_UP_vs_DOWN_heat_capacity <- ks.test(
  HeatCap_HeLa_G1_vs_G2_DOWN$Cp_37C,
  HeatCap_HeLa_G1_vs_G2_UP$Cp_37C
)
print(ks_test_HeLa_G1_vs_G2_UP_vs_DOWN_heat_capacity)

wilcox_test_HeLa_G1_vs_G2_UP_vs_DOWN_heat_capacity <- wilcox.test(
  HeatCap_HeLa_G1_vs_G2_DOWN$Cp_37C,
  HeatCap_HeLa_G1_vs_G2_UP$Cp_37C
)
print(wilcox_test_HeLa_G1_vs_G2_UP_vs_DOWN_heat_capacity)
















# Heat capacity changes MCF7 S vs G0G1
# MCF7 box plot heat capacity all transcripts UP- and DOWN-regulated

Plot_Heat_capacity_MCF7_SvsG0G1 <- ggplot(GSE94479_SvsG0G1, aes(x = factor(log2FoldChange > 0), y = Cp_37C)) +
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

# Save plot
ggsave("Plot_Heat_capacity_MCF7_SvsG0G1.pdf", plot = Plot_Heat_capacity_MCF7_SvsG0G1, width = 4, height = 6, path = "path_to/Desktop")



# TTEST and KS test of the heat capacity values of all transcripts significantly UP- and DOWN-regulated
HeatCap_MCF7_SvsG0G1_UP <- dplyr::filter(GSE94479_SvsG0G1, log2FoldChange > 0)
HeatCap_MCF7_SvsG0G1_DOWN <- dplyr::filter(GSE94479_SvsG0G1, log2FoldChange < 0)

TTEST_heat_capacity_MCF7_SvsG0G1_UP_vs_DOWN <- t.test(
  HeatCap_MCF7_SvsG0G1_DOWN$Cp_37C,
  HeatCap_MCF7_SvsG0G1_UP$Cp_37C,
  paired = FALSE
)
print(TTEST_heat_capacity_MCF7_SvsG0G1_UP_vs_DOWN)

ks_test_MCF7_SvsG0G1_UP_vs_DOWN_heat_capacity <- ks.test(
  HeatCap_MCF7_SvsG0G1_DOWN$Cp_37C,
  HeatCap_MCF7_SvsG0G1_UP$Cp_37C
)
print(ks_test_MCF7_SvsG0G1_UP_vs_DOWN_heat_capacity)

wilcox_test_MCF7_SvsG0G1_UP_vs_DOWN_heat_capacity <- wilcox.test(
  HeatCap_MCF7_SvsG0G1_DOWN$Cp_37C,
  HeatCap_MCF7_SvsG0G1_UP$Cp_37C
)
print(wilcox_test_MCF7_SvsG0G1_UP_vs_DOWN_heat_capacity)












# Heat capacity changes GSE104736 HeLa S vs G1
# Box plot heat capacity for all transcripts UP- and DOWN-regulated

Plot_Heat_capacity_HeLa_S_vs_G1 <- ggplot(MFE_GSE104736_S_vs_G1, aes(x = factor(log2FoldChange > 0), y = Cp_37C)) +
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

# Save the plot
ggsave("Plot_Heat_capacity_HeLa_S_vs_G1.pdf", plot = Plot_Heat_capacity_HeLa_S_vs_G1, width = 4, height = 6, path = "path_to/Desktop")




# TTEST and KS test of the heat capacity values of all transcripts significantly UP- and DOWN-regulated
HeatCap_HeLa_S_vs_G1_UP <- dplyr::filter(MFE_GSE104736_S_vs_G1, log2FoldChange > 0)
HeatCap_HeLa_S_vs_G1_DOWN <- dplyr::filter(MFE_GSE104736_S_vs_G1, log2FoldChange < 0)

TTEST_heat_capacity_HeLa_S_vs_G1_UP_vs_DOWN <- t.test(
  HeatCap_HeLa_S_vs_G1_DOWN$Cp_37C,
  HeatCap_HeLa_S_vs_G1_UP$Cp_37C,
  paired = FALSE
)
print(TTEST_heat_capacity_HeLa_S_vs_G1_UP_vs_DOWN)

ks_test_HeLa_S_vs_G1_UP_vs_DOWN_heat_capacity <- ks.test(
  HeatCap_HeLa_S_vs_G1_DOWN$Cp_37C,
  HeatCap_HeLa_S_vs_G1_UP$Cp_37C
)
print(ks_test_HeLa_S_vs_G1_UP_vs_DOWN_heat_capacity)

wilcox_test_HeLa_S_vs_G1_UP_vs_DOWN_heat_capacity <- wilcox.test(
  HeatCap_HeLa_S_vs_G1_DOWN$Cp_37C,
  HeatCap_HeLa_S_vs_G1_UP$Cp_37C
)
print(wilcox_test_HeLa_S_vs_G1_UP_vs_DOWN_heat_capacity)






