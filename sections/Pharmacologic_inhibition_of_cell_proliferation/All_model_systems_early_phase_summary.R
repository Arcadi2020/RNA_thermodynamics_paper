# This script generates the analysis of the immediate response of all model systems sequenced during the early phase of the adaptive response with Nanopore

# Early-Phase Nanopore Response Analysis (example shown: X778 + palbociclib)
# Aim
#   Integrate transcript-level thermodynamics (Normalized MFE, Specific_heat)
#   with Nanopore DEXSeq counts to profile the EARLY adaptive response.
#
# Samples in scope
#   Sequenced at 72h: 778+palbo, 778+nutlin, 778+tunic, SKBR3+lap, SKMEL28+vem
#   Early-phase subset analyzed here: 778+palbo, SKBR3+lap, SKMEL28+vem
#
# Workflow (key steps)
#   1) Load/rename DEXSeq transcript-count tables per model; harmonize replicate
#      names and order (untreated vs treated).
#   2) Load annotated thermodynamics table (MFE, Specific_heat, length) and
#      QC correlation: Transcript_length vs Normalized_MFE (Pearson).
#   3) Join DEXSeq (padj < 0.05) to thermodynamics; split UP vs DOWN by
#      log2fold_treated_untreated; compare MFE via t.test + KS; plot scatter,
#      violin/box, and density.
#   4) High-/low-MFE composition:
#        • Define thresholds: 75th (high) and 25th (low) percentiles.
#        • Count enrichment among UP vs DOWN (Fisher, two-proportion tests).
#        • Convert counts to long format; aggregate mean counts by condition;
#          compute fractions of total counts attributable to high/low MFE
#          (Control vs Drug) and test (Fisher / prop.test).
#        • Plot MFE vs average counts for Control and Drug panels.
#   5) Heat capacity:
#        • Compare Specific_heat between UP and DOWN (all DE) and within
#          high-MFE subsets; report mean differences; t.test + KS; box plots.
#   6) Biology checks (gene-level, DESeq2):
#        • Load curated gene sets (proliferation, G1/S, HR, Fanconi, etc.).
#        • For each set: filter, compute mean ± SD log2FC; visualize (dot+error,
#          merged box plots); Wilcoxon tests vs 0 to verify expected early arrest.
#
# Outputs
#   - Plots: length–MFE QC; MFE vs log2FC; violin/box/density; MFE vs counts;
#     pathway dot+error and box plots.
#   - Stats: t.test/KS for MFE & heat capacity; Fisher & proportion tests for
#     high/low-MFE enrichment and counts fractions; Wilcoxon per pathway.
#
# Note
#   This file shows the early-phase pipeline for selected models. The same logic,
#   thresholds, and statistics are applied consistently across other cell lines
#   and conditions in the full script.


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




# Function to calculate the mode
calculate_mode <- function(x) {
  unique_x <- unique(x)
  unique_x[which.max(tabulate(match(x, unique_x)))]
}



# We have sequenced 5 samples at 72h (778 palbo, 778 nutlin, 778 tunic, SKBR3 lap and SKMEL28 vem)
# We have only 3 samples sequenced during the early phase response (778 palbo, SKBR3 lap and SKMEL28 vem)

# 778 PALBO cell line early phase response
X778_counts_EP_dexseq <-read.delim("path_to/early.phase_778_palbociclib.tx_counts.both_conditions.treated_untreated_dexseq.tab",
                            header = TRUE, sep = "\t",quote = "") %>%
  dplyr::rename(GENEID = groupID) %>%
  dplyr::rename(ensembl_transcript_id_version = featureID) %>%
  dplyr::relocate(c(countData.early.phase_778_palbociclib_R1, countData.early.phase_778_palbociclib_R2, countData.early.phase_778_palbociclib_R3), .after =  c(countData.early.phase_778_untreated_R1, countData.early.phase_778_untreated_R2, countData.early.phase_778_untreated_R3)) %>%
  dplyr::rename(X778_untreated_R1 = countData.early.phase_778_untreated_R1) %>%
  dplyr::rename(X778_untreated_R2 = countData.early.phase_778_untreated_R2) %>%
  dplyr::rename(X778_untreated_R3 = countData.early.phase_778_untreated_R3) %>%
  dplyr::rename(X778_palbociclib_R1 = countData.early.phase_778_palbociclib_R1) %>%
  dplyr::rename(X778_palbociclib_R2 = countData.early.phase_778_palbociclib_R2) %>%
  dplyr::rename(X778_palbociclib_R3 = countData.early.phase_778_palbociclib_R3) # I rename the columns to be consistent with the analysis here below which was previously generated with the old counts files



# SKBR3 LAP cell line early phase response
SKBR3_counts_EP_dexseq<-read.delim("path_to/early.phase_SKBR3_lapatinib.tx_counts.both_conditions.treated_untreated_dexseq.tab",
                             header = TRUE, sep = "\t",quote = "") %>%
  dplyr::rename(ensembl_transcript_id_version = featureID) %>%
  dplyr::relocate(c(countData.early.phase_SKBR3_lapatinib_R1, countData.early.phase_SKBR3_lapatinib_R2, countData.early.phase_SKBR3_lapatinib_R3), .after =  c(countData.early.phase_SKBR3_untreated_R1, countData.early.phase_SKBR3_untreated_R2, countData.early.phase_SKBR3_untreated_R3)) %>%
  dplyr::rename(SKBR3_untreated_R1 = countData.early.phase_SKBR3_untreated_R1) %>%
  dplyr::rename(SKBR3_untreated_R2 = countData.early.phase_SKBR3_untreated_R2) %>%
  dplyr::rename(SKBR3_untreated_R3 = countData.early.phase_SKBR3_untreated_R3) %>%
  dplyr::rename(SKBR3_lapatinib_R1 = countData.early.phase_SKBR3_lapatinib_R1) %>%
  dplyr::rename(SKBR3_lapatinib_R2 = countData.early.phase_SKBR3_lapatinib_R2) %>%
  dplyr::rename(SKBR3_lapatinib_R3 = countData.early.phase_SKBR3_lapatinib_R3) # I rename the columns to be consistent with the analysis here below which was previously generated with the old counts files



# SKMEL28 VEM cell line early phase response
SKMEL28_counts_EP_dexseq <- read.delim("path_to/early.phase_SKMEL28_vemurafenib.tx_counts.both_conditions.treated_untreated_dexseq.tab",
                             header = TRUE, sep = "\t",quote = "") %>%
  dplyr::rename(ensembl_transcript_id_version = featureID) %>%
  dplyr::relocate(c(countData.early.phase_SKMEL28_vemurafenib_R1, countData.early.phase_SKMEL28_vemurafenib_R2, countData.early.phase_SKMEL28_vemurafenib_R3), .after =  c(countData.early.phase_SKMEL28_untreated_R1, countData.early.phase_SKMEL28_untreated_R2, countData.early.phase_SKMEL28_untreated_R3)) %>%
  dplyr::rename(SKMEL28_untreated_R1 = countData.early.phase_SKMEL28_untreated_R1) %>%
  dplyr::rename(SKMEL28_untreated_R2 = countData.early.phase_SKMEL28_untreated_R2) %>%
  dplyr::rename(SKMEL28_untreated_R3 = countData.early.phase_SKMEL28_untreated_R3) %>%
  dplyr::rename(SKMEL28_vemurafenib_R1 = countData.early.phase_SKMEL28_vemurafenib_R1) %>%
  dplyr::rename(SKMEL28_vemurafenib_R2 = countData.early.phase_SKMEL28_vemurafenib_R2) %>%
  dplyr::rename(SKMEL28_vemurafenib_R3 = countData.early.phase_SKMEL28_vemurafenib_R3) # I rename the columns to be consistent with the analysis here below which was previously generated with the old counts files














# Load Table_ANNOTATED_MFE_heat_capacity ####
Gencode_TX_norm_MFE_HC_table <- read.delim("path_to/Gencode_TX_norm_MFE_HC_table.csv",
                                           sep = ",", header = T, stringsAsFactors = F)
Table_ANNOTATED <- Gencode_TX_norm_MFE_HC_table

# This is the complete table with MFE values of the entire Gencode transcriptome and the heat capacity values limited to the Nanopore annotated transcripts




# Scatter plot sequence length vs MFE values
library("ggpubr")

Plot <- ggscatter(Table_ANNOTATED, x="Transcript_length", y="Normalized_MFE",
                  color = "black", shape = 21, size = 1.5, # Points color, shape and size
                  add = "reg.line", conf.int = TRUE, 
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  xlab = "Sequence_length (nt)", ylab ="LncFinder_MFE",title="Correlation transcript length and MFE values")
Plot

# Calculate correlation coefficient between sequence length and MFE values.
Correlation <- cor.test(Table_ANNOTATED$Normalized_MFE, Table_ANNOTATED$Transcript_length, 
                        method = "pearson")
Correlation













# 778 palbociclib ANNOTATED transcripts ####
# 778 PALBO cell line early phase response
X778_counts_EP<- dplyr::inner_join(X778_counts_EP_dexseq, Table_ANNOTATED, by="ensembl_transcript_id_version") %>%
  dplyr::filter(padj < 0.05) # The file 'Table' contains the Normalized_MFE values of all ANNOTATED transcripts


# Generate a scatter plot Log2_FC vs Normalized_MFE
library(ggplot2)
ggplot(X778_counts_EP, aes(y=Normalized_MFE, x=log2fold_treated_untreated)) +
  geom_point(size=1.5, alpha=0.6)+
  theme_classic() +
  theme(legend.position="none") +
  geom_vline(xintercept = 0, size=0.3) +
  #xlim(-5.5,5.5) +
  theme(axis.title = element_text(size = 20)) +
  theme(plot.title = element_text(size = 15)) +
  labs(y="Normalized_MFE",x="Log2_FC_expression",title="X778_PALBO_shared_ANNOTATED_TX_early_phase_Log2_FC_adj.P<0.05")


X778_PALBO_EP_dexseq_ANNOTATED_UP_EP <- dplyr::filter(X778_counts_EP, log2fold_treated_untreated > 0)
Mean_X778_PALBO_EP_dexseq_ANNOTATED_UP_EP <- mean(X778_PALBO_EP_dexseq_ANNOTATED_UP_EP$Normalized_MFE, na.rm = TRUE)

X778_PALBO_EP_dexseq_ANNOTATED_DOWN_EP <- dplyr::filter(X778_counts_EP, log2fold_treated_untreated < 0)
Mean_X778_PALBO_EP_dexseq_ANNOTATED_DOWN_EP <- mean(X778_PALBO_EP_dexseq_ANNOTATED_DOWN_EP$Normalized_MFE, na.rm = TRUE)

ttest_MFE_X778_PALBO_PALBO_72h_all_UP_EP_vs_DOWN_EP <- t.test(X778_PALBO_EP_dexseq_ANNOTATED_DOWN_EP$Normalized_MFE, X778_PALBO_EP_dexseq_ANNOTATED_UP_EP$Normalized_MFE, paired=FALSE)
print(ttest_MFE_X778_PALBO_PALBO_72h_all_UP_EP_vs_DOWN_EP)

# Perform K-S test to compare distributions
ks_test_MFE_X778_PALBO_UP_vs_DOWN_EP <- ks.test(
  na.omit(X778_PALBO_EP_dexseq_ANNOTATED_DOWN_EP$Normalized_MFE),
  na.omit(X778_PALBO_EP_dexseq_ANNOTATED_UP_EP$Normalized_MFE)
)
print(ks_test_MFE_X778_PALBO_UP_vs_DOWN_EP)


# Show the data with a violin plot
ggplot(X778_counts_EP, aes(x = factor(log2fold_treated_untreated > 0), y = Normalized_MFE)) +
  geom_violin(fill = "lightgrey", size = 1.5) +
  geom_boxplot(width = 0.1, outlier.shape = NA) + 
  #geom_jitter(width = 0.1, size = 1, alpha = 0.3) +
  theme_classic() +
  theme(legend.position = "none", axis.title = element_text(size = 20), plot.title = element_text(size = 10),
        axis.text.x = element_text(size = 16, color = "black")) +
  geom_vline(xintercept = 0, size = 0.3) +
  labs(y = "Normalized_MFE", x = NULL, title = "778_shared_ANNOTATED_TX_early_phase_Log2_FC_adj.P<0.05") +
  scale_x_discrete(labels = c("Log2_FC < 0", "Log2_FC > 0"))
# In this modified code, the x aesthetic is set as a factor variable factor(log2fold_treated_untreated > 0). 
# This creates two separate groups based on the condition Log2_FC > 0 or Log2_FC < 0. Each group is represented by a different violin shape.

ks_test_MFE_X778_PALBO_PALBO_72h_all_UP_EP_vs_DOWN_EP <- ks.test(X778_PALBO_EP_dexseq_ANNOTATED_DOWN_EP$Normalized_MFE, X778_PALBO_EP_dexseq_ANNOTATED_UP_EP$Normalized_MFE)
print(ks_test_MFE_X778_PALBO_PALBO_72h_all_UP_EP_vs_DOWN_EP)

# Look at the differences between the 0.75 quantile of UP- and DOWN-regulated transcripts
quantile_075_778_UP_EP <- quantile(X778_PALBO_EP_dexseq_ANNOTATED_UP_EP$Normalized_MFE, probs = c(0.75)) 
quantile_075_778_DOWN_EP <- quantile(X778_PALBO_EP_dexseq_ANNOTATED_DOWN_EP$Normalized_MFE, probs = c(0.75)) 
ks_test_778_075 <- ks.test(X778_PALBO_EP_dexseq_ANNOTATED_UP_EP$Normalized_MFE[X778_PALBO_EP_dexseq_ANNOTATED_UP_EP$Normalized_MFE >= quantile_075_778_UP_EP], X778_PALBO_EP_dexseq_ANNOTATED_DOWN_EP$Normalized_MFE[X778_PALBO_EP_dexseq_ANNOTATED_DOWN_EP$Normalized_MFE >= quantile_075_778_DOWN_EP])
print(ks_test_778_075)

# Look at the differences between the 0.25 quantile of UP- and DOWN-regulated transcripts
quantile_025_778_UP_EP <- quantile(X778_PALBO_EP_dexseq_ANNOTATED_UP_EP$Normalized_MFE, probs = c(0.25)) 
quantile_025_778_DOWN_EP <- quantile(X778_PALBO_EP_dexseq_ANNOTATED_DOWN_EP$Normalized_MFE, probs = c(0.25)) 
ks_test_778_025 <- ks.test(X778_PALBO_EP_dexseq_ANNOTATED_UP_EP$Normalized_MFE[X778_PALBO_EP_dexseq_ANNOTATED_UP_EP$Normalized_MFE <= quantile_025_778_UP_EP], X778_PALBO_EP_dexseq_ANNOTATED_DOWN_EP$Normalized_MFE[X778_PALBO_EP_dexseq_ANNOTATED_DOWN_EP$Normalized_MFE <= quantile_025_778_DOWN_EP])
print(ks_test_778_025)




# Show the data in a density plot
max_length_778_EP_MFE <- max(c(length(X778_PALBO_EP_dexseq_ANNOTATED_UP_EP$Normalized_MFE), length(X778_PALBO_EP_dexseq_ANNOTATED_DOWN_EP$Normalized_MFE)))

data_778_EP_MFE_density <- data.frame(MFE_TX_UP_EP_778 = c(X778_PALBO_EP_dexseq_ANNOTATED_UP_EP$Normalized_MFE, rep(NA, max_length_778_EP_MFE - length(X778_PALBO_EP_dexseq_ANNOTATED_UP_EP$Normalized_MFE))), # This is useful when we need to create a data frame and the number of values between columns are different
                                   MFE_TX_DOWN_EP_778 = c(X778_PALBO_EP_dexseq_ANNOTATED_DOWN_EP$Normalized_MFE, rep(NA, max_length_778_EP_MFE- length(X778_PALBO_EP_dexseq_ANNOTATED_DOWN_EP$Normalized_MFE))))
library(reshape)
data_778_EP_MFE_density <- melt(data_778_EP_MFE_density)

library(ggplot2)
ggplot(data= data_778_EP_MFE_density, aes(x=value, group=variable, color=variable)) +
  geom_density (size=1) +
  theme_classic() +
  #xlim(0,1.5) +
  labs(title = "MFE_distribution_ANNOTATED_TX_778_early_phase")+
  labs(y="Density")+
  labs(x="DGfold / nt (kcal/mol)")






## Analysis net MFE and heat capacity ####
# Determine whether within the differentially expressed genes there is a significant up-regulation of high MFE transcripts:
# Define high MFE transcripts by taking the 75th percentile of the combined MFE distributions of control and drug-treated samples.
# Count the number of transcripts in the significantly up- and down-regulated genes that have MFE values above this threshold.
# Perform a Fisher’s exact test to compare the proportions of high MFE transcripts between the two groups.


X778_PALBO_EP_dexseq_ANNOTATED_norm_TX_counts_MFE_df <- X778_counts_EP %>%
  dplyr::select(ensembl_transcript_id_version, ensembl_transcript_id, ensembl_gene_id, external_gene_name, transcript_is_canonical, Transcript_length,
                X778_untreated_R1, X778_untreated_R2, X778_untreated_R3,
                X778_palbociclib_R1, X778_palbociclib_R2, X778_palbociclib_R3,
                Normalized_MFE, Temperature, Specific_heat)


# Identify the 75th percentile of the combined MFE distributions of control and drug-treated samples
quantile_075_X778_PALBO_EP_dexseq_ANNOTATED_norm_TX_counts_MFE_df <- quantile(X778_PALBO_EP_dexseq_ANNOTATED_norm_TX_counts_MFE_df$Normalized_MFE, probs = c(0.75))

# Select all the transcripts with Normalized_MFE values above the 75th percentile
X778_PALBO_EP_dexseq_ANNOTATED_high_MFE_TX <- dplyr::filter(X778_PALBO_EP_dexseq_ANNOTATED_norm_TX_counts_MFE_df, Normalized_MFE > quantile_075_X778_PALBO_EP_dexseq_ANNOTATED_norm_TX_counts_MFE_df)

# Count the number of transcripts in the significantly up- and down-regulated pools that have MFE values above this threshold.
X778_PALBO_EP_dexseq_ANNOTATED_high_MFE_TX_in_UP <- dplyr::filter(X778_PALBO_EP_dexseq_ANNOTATED_UP_EP, ensembl_transcript_id_version %in% X778_PALBO_EP_dexseq_ANNOTATED_high_MFE_TX$ensembl_transcript_id_version)
X778_PALBO_EP_dexseq_ANNOTATED_high_MFE_TX_in_DOWN <- dplyr::filter(X778_PALBO_EP_dexseq_ANNOTATED_DOWN_EP, ensembl_transcript_id_version %in% X778_PALBO_EP_dexseq_ANNOTATED_high_MFE_TX$ensembl_transcript_id_version)


# Perform a Fisher’s exact test to compare the proportions of high MFE transcripts between the two groups.
X778_PALBO_EP_dexseq_ANNOTATED_total_UP <- nrow(X778_PALBO_EP_dexseq_ANNOTATED_UP_EP)
X778_PALBO_EP_dexseq_ANNOTATED_total_DOWN <- nrow(X778_PALBO_EP_dexseq_ANNOTATED_DOWN_EP)
X778_PALBO_EP_dexseq_ANNOTATED_count_high_MFE_UP <- nrow(X778_PALBO_EP_dexseq_ANNOTATED_high_MFE_TX_in_UP)
X778_PALBO_EP_dexseq_ANNOTATED_count_high_MFE_DOWN <- nrow(X778_PALBO_EP_dexseq_ANNOTATED_high_MFE_TX_in_DOWN)

X778_PALBO_EP_dexseq_ANNOTATED_contingency_table <- matrix(c(X778_PALBO_EP_dexseq_ANNOTATED_count_high_MFE_UP, X778_PALBO_EP_dexseq_ANNOTATED_total_UP - X778_PALBO_EP_dexseq_ANNOTATED_count_high_MFE_UP,
                                        X778_PALBO_EP_dexseq_ANNOTATED_count_high_MFE_DOWN, X778_PALBO_EP_dexseq_ANNOTATED_total_DOWN - X778_PALBO_EP_dexseq_ANNOTATED_count_high_MFE_DOWN),
                                      nrow = 2, byrow = TRUE,
                                      dimnames = list(UP = c("UP-regulated", "DOWN-regulated"),
                                                      DOWN = c("High MFE TX", "Not High MFE TX")))

X778_PALBO_EP_dexseq_ANNOTATED_contingency_table_df <- as.data.frame(X778_PALBO_EP_dexseq_ANNOTATED_contingency_table)

X778_PALBO_EP_dexseq_ANNOTATED_fisher_result <- fisher.test(X778_PALBO_EP_dexseq_ANNOTATED_contingency_table)

cat("X778_PALBO_EP_dexseq_ANNOTATED P-value from Fisher test:", X778_PALBO_EP_dexseq_ANNOTATED_fisher_result$p.value, "\n")


# Look at the fractions
X778_PALBO_EP_dexseq_ANNOTATED_high_MFE_in_UP_fraction <- X778_PALBO_EP_dexseq_ANNOTATED_count_high_MFE_UP / X778_PALBO_EP_dexseq_ANNOTATED_total_UP
X778_PALBO_EP_dexseq_ANNOTATED_high_MFE_in_DOWN_fraction <- X778_PALBO_EP_dexseq_ANNOTATED_count_high_MFE_DOWN / X778_PALBO_EP_dexseq_ANNOTATED_total_DOWN

# Output the results
cat("Fraction of X778_PALBO_EP_dexseq_ANNOTATED high MFE TX in UP-regulated genes:", X778_PALBO_EP_dexseq_ANNOTATED_high_MFE_in_UP_fraction, "\n") # This value goes in the summary table
cat("Fraction of X778_PALBO_EP_dexseq_ANNOTATED high MFE TX in DOWN-regulated genes:", X778_PALBO_EP_dexseq_ANNOTATED_high_MFE_in_DOWN_fraction, "\n") # This value goes in the summary table
cat("Mean of X778_PALBO_EP_dexseq_ANNOTATED high MFE TX in UP-regulated genes:", mean(X778_PALBO_EP_dexseq_ANNOTATED_high_MFE_TX_in_UP$Normalized_MFE), "\n")
cat("Mean of X778_PALBO_EP_dexseq_ANNOTATED high MFE TX in DOWN-regulated genes:", mean(X778_PALBO_EP_dexseq_ANNOTATED_high_MFE_TX_in_DOWN$Normalized_MFE), "\n")
cat("Mode of X778_PALBO_EP_dexseq_ANNOTATED high MFE TX in UP-regulated genes:", calculate_mode(X778_PALBO_EP_dexseq_ANNOTATED_high_MFE_TX_in_UP$Normalized_MFE), "\n")
cat("Mode of X778_PALBO_EP_dexseq_ANNOTATED high MFE TX in DOWN-regulated genes:", calculate_mode(X778_PALBO_EP_dexseq_ANNOTATED_high_MFE_TX_in_DOWN$Normalized_MFE), "\n")





# Perform a two-proportion z-test to compare proportions
# Calculate the proportions and total counts
X778_PALBO_EP_dexseq_ANNOTATED_high_MFE_in_UP_proportion <- X778_PALBO_EP_dexseq_ANNOTATED_high_MFE_in_UP_fraction
X778_PALBO_EP_dexseq_ANNOTATED_high_MFE_in_DOWN_proportion <- X778_PALBO_EP_dexseq_ANNOTATED_high_MFE_in_DOWN_fraction

# Combined proportion for pooled standard error
combined_proportion <- (X778_PALBO_EP_dexseq_ANNOTATED_count_high_MFE_UP + X778_PALBO_EP_dexseq_ANNOTATED_count_high_MFE_DOWN) /
  (X778_PALBO_EP_dexseq_ANNOTATED_total_UP + X778_PALBO_EP_dexseq_ANNOTATED_total_DOWN)

# Standard error of the difference in proportions
standard_error <- sqrt(
  combined_proportion * (1 - combined_proportion) * 
    (1 / X778_PALBO_EP_dexseq_ANNOTATED_total_UP + 1 / X778_PALBO_EP_dexseq_ANNOTATED_total_DOWN)
)

# Z-statistic
z_statistic <- (X778_PALBO_EP_dexseq_ANNOTATED_high_MFE_in_UP_proportion - X778_PALBO_EP_dexseq_ANNOTATED_high_MFE_in_DOWN_proportion) /
  standard_error

# Two-tailed p-value
X778_PALBO_EP_dexseq_ANNOTATED_z_test_MFE_p_value <- 2 * (1 - pnorm(abs(z_statistic)))






# Show data as violin plot
# Add a new column to indicate if the value is above the 75th percentile
X778_counts_EP$above_quantile_075 <- X778_PALBO_EP_dexseq_ANNOTATED_norm_TX_counts_MFE_df$Normalized_MFE > quantile_075_X778_PALBO_EP_dexseq_ANNOTATED_norm_TX_counts_MFE_df

# Generate violin plot with conditional jitter coloring
ggplot(X778_counts_EP, aes(x = factor(log2fold_treated_untreated > 0), y = Normalized_MFE)) +
  geom_violin(fill = "white") +
  geom_jitter(aes(color = above_quantile_075), width = 0.1, size = 1, alpha = 1) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "blue")) +
  theme_classic() +
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
    title = "X778_PALBO_EP_dexseq_ANNOTATED_early_phase_Log2_FC_adj.P<0.05"
  ) +
  scale_x_discrete(labels = c("log2FoldChange < 0", "log2FoldChange > 0"))














# Determine the impact of these changes in terms of fraction of total counts:
# Calculate the fraction of total counts of high MFE transcripts in control and drug-treated groups. 
# This is done by summing the counts of the high MFE transcripts and dividing by the total counts in each group.
# Fisher’s exact test to compare the proportions of high MFE transcripts counts between the two groups.


# Convert from wide to long format
X778_PALBO_EP_dexseq_ANNOTATED_long_data <- X778_PALBO_EP_dexseq_ANNOTATED_norm_TX_counts_MFE_df %>%
  pivot_longer(
    cols = starts_with("X778"),      # Select all columns starting with "X778"
    names_to = "replicate",         # Name of the new column for replicate identifiers
    values_to = "counts"            # Name of the new column for the counts
  )

# Add a condition column based on the replicate
X778_PALBO_EP_dexseq_ANNOTATED_long_data <- X778_PALBO_EP_dexseq_ANNOTATED_long_data %>%
  mutate(condition = ifelse(replicate %in% c("X778_untreated_R1", "X778_untreated_R2", "X778_untreated_R3"), "Control", "palbociclib"))


# Aggregate the replicates for each condition
X778_PALBO_EP_dexseq_ANNOTATED_norm_TX_aggregated_data_avg_counts <- X778_PALBO_EP_dexseq_ANNOTATED_long_data %>%
  group_by(ensembl_transcript_id_version, condition, Normalized_MFE, Specific_heat) %>%
  summarise(
    avg_counts = mean(counts),
    .groups = 'drop'
  )


# Define the 75th percentile of the entire distribution of MFE values for the X778_PALBO_EP_dexseq_ANNOTATED data
X778_PALBO_EP_dexseq_ANNOTATED_MFE_all_TX_75_percentile <- quantile(X778_PALBO_EP_dexseq_ANNOTATED_norm_TX_aggregated_data_avg_counts$Normalized_MFE, 0.75)

# Identify all the transcripts with MFE values above the 75th percentile
X778_PALBO_EP_dexseq_ANNOTATED_high_MFE_TX_total <- dplyr::filter(X778_PALBO_EP_dexseq_ANNOTATED_norm_TX_counts_MFE_df, Normalized_MFE > X778_PALBO_EP_dexseq_ANNOTATED_MFE_all_TX_75_percentile)





# Fraction of total high MFE transcripts counts in control and drug-treated groups
X778_PALBO_EP_dexseq_ANNOTATED_control_total_counts <- dplyr::filter(X778_PALBO_EP_dexseq_ANNOTATED_norm_TX_aggregated_data_avg_counts, condition == "Control") %>%
  dplyr::summarise(total_counts = sum(avg_counts))

X778_PALBO_EP_dexseq_ANNOTATED_control_high_MFE_TX_counts <- dplyr::filter(X778_PALBO_EP_dexseq_ANNOTATED_norm_TX_aggregated_data_avg_counts, condition == "Control") %>%
  dplyr::filter(ensembl_transcript_id_version %in% X778_PALBO_EP_dexseq_ANNOTATED_high_MFE_TX_total$ensembl_transcript_id_version) %>%
  dplyr::summarise(total_counts_high_MFE_TX = sum(avg_counts))

X778_PALBO_EP_dexseq_ANNOTATED_control_high_MFE_TX_counts_fraction_EP <- X778_PALBO_EP_dexseq_ANNOTATED_control_high_MFE_TX_counts$total_counts_high_MFE_TX / X778_PALBO_EP_dexseq_ANNOTATED_control_total_counts$total_counts
print(X778_PALBO_EP_dexseq_ANNOTATED_control_high_MFE_TX_counts_fraction_EP) # This value goes in the summary table





X778_PALBO_EP_dexseq_ANNOTATED_palbociclib_total_counts <- dplyr::filter(X778_PALBO_EP_dexseq_ANNOTATED_norm_TX_aggregated_data_avg_counts, condition == "palbociclib") %>%
  dplyr::summarise(total_counts = sum(avg_counts))

X778_PALBO_EP_dexseq_ANNOTATED_palbociclib_high_MFE_TX_counts <- dplyr::filter(X778_PALBO_EP_dexseq_ANNOTATED_norm_TX_aggregated_data_avg_counts, condition == "palbociclib") %>%
  dplyr::filter(ensembl_transcript_id_version %in% X778_PALBO_EP_dexseq_ANNOTATED_high_MFE_TX_total$ensembl_transcript_id_version) %>%
  dplyr::summarise(total_counts_high_MFE_TX = sum(avg_counts))

X778_PALBO_EP_dexseq_ANNOTATED_palbociclib_high_MFE_TX_counts_fraction_EP <- X778_PALBO_EP_dexseq_ANNOTATED_palbociclib_high_MFE_TX_counts$total_counts_high_MFE_TX / X778_PALBO_EP_dexseq_ANNOTATED_palbociclib_total_counts$total_counts
print(X778_PALBO_EP_dexseq_ANNOTATED_palbociclib_high_MFE_TX_counts_fraction_EP) # This value goes in the summary table





# Fisher's exact test to compare the proportions of high MFE transcripts counts between the two groups
# Create the groups
X778_PALBO_EP_dexseq_ANNOTATED_high_MFE_TX_counts_control_CT <- X778_PALBO_EP_dexseq_ANNOTATED_control_high_MFE_TX_counts_fraction_EP * X778_PALBO_EP_dexseq_ANNOTATED_control_total_counts$total_counts
X778_PALBO_EP_dexseq_ANNOTATED_low_MFE_TX_counts_control_CT <- X778_PALBO_EP_dexseq_ANNOTATED_control_total_counts$total_counts - X778_PALBO_EP_dexseq_ANNOTATED_control_high_MFE_TX_counts$total_counts_high_MFE_TX

X778_PALBO_EP_dexseq_ANNOTATED_high_MFE_TX_counts_palbociclib_CT <- X778_PALBO_EP_dexseq_ANNOTATED_palbociclib_high_MFE_TX_counts_fraction_EP * X778_PALBO_EP_dexseq_ANNOTATED_palbociclib_total_counts$total_counts
X778_PALBO_EP_dexseq_ANNOTATED_low_MFE_TX_counts_palbociclib_CT <- X778_PALBO_EP_dexseq_ANNOTATED_palbociclib_total_counts$total_counts - X778_PALBO_EP_dexseq_ANNOTATED_palbociclib_high_MFE_TX_counts$total_counts_high_MFE_TX

# Construct the contingency table with raw counts
X778_PALBO_EP_dexseq_ANNOTATED_counts_contingency_table <- matrix(
  c(
    X778_PALBO_EP_dexseq_ANNOTATED_control_high_MFE_TX_counts$total_counts_high_MFE_TX, 
    X778_PALBO_EP_dexseq_ANNOTATED_palbociclib_high_MFE_TX_counts$total_counts_high_MFE_TX,
    X778_PALBO_EP_dexseq_ANNOTATED_control_total_counts$total_counts - X778_PALBO_EP_dexseq_ANNOTATED_control_high_MFE_TX_counts$total_counts_high_MFE_TX, 
    X778_PALBO_EP_dexseq_ANNOTATED_palbociclib_total_counts$total_counts - X778_PALBO_EP_dexseq_ANNOTATED_palbociclib_high_MFE_TX_counts$total_counts_high_MFE_TX
  ),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(Category = c("High_MFE_TX_counts", "Low_MFE_TX_counts"),
                  Group = c("Control", "palbociclib"))
)

# Ensure the counts are rounded to integers (they should already be integers if the input data is correct)
# X778_PALBO_EP_dexseq_ANNOTATED_counts_contingency_table <- round(X778_PALBO_EP_dexseq_ANNOTATED_counts_contingency_table)

X778_PALBO_EP_dexseq_ANNOTATED_counts_contingency_table_df <- as.data.frame(X778_PALBO_EP_dexseq_ANNOTATED_counts_contingency_table)

X778_PALBO_EP_dexseq_ANNOTATED_fisher_result_counts <- fisher.test(X778_PALBO_EP_dexseq_ANNOTATED_counts_contingency_table_df)
print(X778_PALBO_EP_dexseq_ANNOTATED_fisher_result_counts)




cat("X778_PALBO_EP_dexseq_ANNOTATED Fisher P-value Fisher proportions of high MFE transcripts counts in drug-treated and control :", X778_PALBO_EP_dexseq_ANNOTATED_fisher_result_counts$p.value, "\n")  # This value goes in the summary table


# Output the results
cat("X778_PALBO_EP_dexseq_ANNOTATED control high MFE TX as fraction of total counts:", X778_PALBO_EP_dexseq_ANNOTATED_control_high_MFE_TX_counts_fraction_EP, "\n")
cat("X778_PALBO_EP_dexseq_ANNOTATED palbociclib high MFE TX as fraction of total counts:", X778_PALBO_EP_dexseq_ANNOTATED_palbociclib_high_MFE_TX_counts_fraction_EP, "\n")
cat("X778_PALBO_EP_dexseq_ANNOTATED net change counts high MFE TX drug vs control", X778_PALBO_EP_dexseq_ANNOTATED_palbociclib_high_MFE_TX_counts_fraction_EP-X778_PALBO_EP_dexseq_ANNOTATED_control_high_MFE_TX_counts_fraction_EP, "\n")




# Perform a two-proportion Z-test to compare the proportions of high MFE transcripts counts between the two groups
# The arguments are the counts of successes and the total number of trials for each group

# Number of high MFE TX counts in control and palbociclib groups
high_MFE_control <- X778_PALBO_EP_dexseq_ANNOTATED_control_high_MFE_TX_counts$total_counts_high_MFE_TX
high_MFE_palbociclib <- X778_PALBO_EP_dexseq_ANNOTATED_palbociclib_high_MFE_TX_counts$total_counts_high_MFE_TX

# Total counts in control and palbociclib groups
total_control <- X778_PALBO_EP_dexseq_ANNOTATED_control_total_counts$total_counts
total_palbociclib <- X778_PALBO_EP_dexseq_ANNOTATED_palbociclib_total_counts$total_counts

# Perform the two-proportion Z-test
X778_PALBO_EP_dexseq_ANNOTATED_z_test_result <- prop.test(
  x = c(high_MFE_control, high_MFE_palbociclib),    # successes: high MFE counts in both groups
  n = c(total_control, total_palbociclib)           # total counts in both groups
)

# Output the result of the two-proportion Z-test
X778_PALBO_EP_dexseq_ANNOTATED_z_test_counts_75_p_value <- X778_PALBO_EP_dexseq_ANNOTATED_z_test_result$p.value








# Plot Normalized_MFE vs transcript counts of both control and drug-treated conditions
# Create separate data frames for Control and palbociclib conditions
X778_PALBO_EP_dexseq_ANNOTATED_control_data <- X778_PALBO_EP_dexseq_ANNOTATED_norm_TX_aggregated_data_avg_counts %>%
  dplyr::filter(condition == "Control")

X778_PALBO_EP_dexseq_ANNOTATED_palbociclib_data <- X778_PALBO_EP_dexseq_ANNOTATED_norm_TX_aggregated_data_avg_counts %>%
  dplyr::filter(condition == "palbociclib") 


# Generate the plot
Plot_X778_PALBO_EP_dexseq_ANNOTATED_control <- ggplot(data = X778_PALBO_EP_dexseq_ANNOTATED_control_data, aes(x = avg_counts, y = Normalized_MFE, color = condition)) +
  geom_point(alpha = 0.6) +
  theme_classic() +
  labs(
    title = "X778_PALBO_EP_dexseq_ANNOTATED MFE vs. total avg counts Control",
    x = "Transcript Counts",
    y = "Normalized MFE (kcal/mol)"
  ) +
  #xlim(0, 150000) +
  ylim(-0.8, 0) +
  theme(
    axis.title = element_text(size = 15),
    plot.title = element_text(size = 10)
  ) +
  guides(color = FALSE)


Plot_X778_PALBO_EP_dexseq_ANNOTATED_palbociclib <- ggplot(data = X778_PALBO_EP_dexseq_ANNOTATED_palbociclib_data, aes(x = avg_counts, y = Normalized_MFE, color = condition)) +
  geom_point(alpha = 0.6) +
  theme_classic() +
  labs(
    title = "X778_PALBO_EP_dexseq_ANNOTATED MFE vs. total avg counts palbociclib",
    x = "Transcript Counts",
    y = "Normalized MFE (kcal/mol)"
  ) +
  #xlim(0, 150000) +
  ylim(-0.8, 0) +
  theme(
    axis.title = element_text(size = 15),
    plot.title = element_text(size = 10)
  ) +
  guides(color = FALSE)

library(patchwork)
X778_PALBO_EP_dexseq_ANNOTATED_combined_plots <- Plot_X778_PALBO_EP_dexseq_ANNOTATED_control + Plot_X778_PALBO_EP_dexseq_ANNOTATED_palbociclib
X778_PALBO_EP_dexseq_ANNOTATED_combined_plots








# Define the 25th percentile of the entire distribution of MFE values for the X778_PALBO_EP_dexseq_ANNOTATED data
X778_PALBO_EP_dexseq_ANNOTATED_MFE_all_TX_25_percentile <- quantile(X778_PALBO_EP_dexseq_ANNOTATED_norm_TX_aggregated_data_avg_counts$Normalized_MFE, 0.25)

# Identify all the transcripts with MFE values below the 25th percentile
X778_PALBO_EP_dexseq_ANNOTATED_low_MFE_TX_total <- dplyr::filter(X778_PALBO_EP_dexseq_ANNOTATED_norm_TX_counts_MFE_df, Normalized_MFE < X778_PALBO_EP_dexseq_ANNOTATED_MFE_all_TX_25_percentile)





# Fraction of total low MFE transcripts counts in control and drug-treated groups
X778_PALBO_EP_dexseq_ANNOTATED_control_total_counts <- dplyr::filter(X778_PALBO_EP_dexseq_ANNOTATED_norm_TX_aggregated_data_avg_counts, condition == "Control") %>%
  dplyr::summarise(total_counts = sum(avg_counts))

X778_PALBO_EP_dexseq_ANNOTATED_control_low_MFE_TX_counts <- dplyr::filter(X778_PALBO_EP_dexseq_ANNOTATED_norm_TX_aggregated_data_avg_counts, condition == "Control") %>%
  dplyr::filter(ensembl_transcript_id_version %in% X778_PALBO_EP_dexseq_ANNOTATED_low_MFE_TX_total$ensembl_transcript_id_version) %>%
  dplyr::summarise(total_counts_low_MFE_TX = sum(avg_counts))

X778_PALBO_EP_dexseq_ANNOTATED_control_low_MFE_TX_counts_fraction_EP <- X778_PALBO_EP_dexseq_ANNOTATED_control_low_MFE_TX_counts$total_counts_low_MFE_TX / X778_PALBO_EP_dexseq_ANNOTATED_control_total_counts$total_counts
print(X778_PALBO_EP_dexseq_ANNOTATED_control_low_MFE_TX_counts_fraction_EP) # This value goes in the summary table





X778_PALBO_EP_dexseq_ANNOTATED_palbociclib_total_counts <- dplyr::filter(X778_PALBO_EP_dexseq_ANNOTATED_norm_TX_aggregated_data_avg_counts, condition == "palbociclib") %>%
  dplyr::summarise(total_counts = sum(avg_counts))

X778_PALBO_EP_dexseq_ANNOTATED_palbociclib_low_MFE_TX_counts <- dplyr::filter(X778_PALBO_EP_dexseq_ANNOTATED_norm_TX_aggregated_data_avg_counts, condition == "palbociclib") %>%
  dplyr::filter(ensembl_transcript_id_version %in% X778_PALBO_EP_dexseq_ANNOTATED_low_MFE_TX_total$ensembl_transcript_id_version) %>%
  dplyr::summarise(total_counts_low_MFE_TX = sum(avg_counts))

X778_PALBO_EP_dexseq_ANNOTATED_palbociclib_low_MFE_TX_counts_fraction_EP <- X778_PALBO_EP_dexseq_ANNOTATED_palbociclib_low_MFE_TX_counts$total_counts_low_MFE_TX / X778_PALBO_EP_dexseq_ANNOTATED_palbociclib_total_counts$total_counts
print(X778_PALBO_EP_dexseq_ANNOTATED_palbociclib_low_MFE_TX_counts_fraction_EP) # This value goes in the summary table





# Fisher's exact test to compare the proportions of low MFE transcripts counts between the two groups
# Create the groups
X778_PALBO_EP_dexseq_ANNOTATED_low_MFE_TX_counts_control_CT <- X778_PALBO_EP_dexseq_ANNOTATED_control_low_MFE_TX_counts_fraction_EP * X778_PALBO_EP_dexseq_ANNOTATED_control_total_counts$total_counts
X778_PALBO_EP_dexseq_ANNOTATED_low_MFE_TX_counts_control_CT <- X778_PALBO_EP_dexseq_ANNOTATED_control_total_counts$total_counts - X778_PALBO_EP_dexseq_ANNOTATED_control_low_MFE_TX_counts$total_counts_low_MFE_TX

X778_PALBO_EP_dexseq_ANNOTATED_low_MFE_TX_counts_palbociclib_CT <- X778_PALBO_EP_dexseq_ANNOTATED_palbociclib_low_MFE_TX_counts_fraction_EP * X778_PALBO_EP_dexseq_ANNOTATED_palbociclib_total_counts$total_counts
X778_PALBO_EP_dexseq_ANNOTATED_low_MFE_TX_counts_palbociclib_CT <- X778_PALBO_EP_dexseq_ANNOTATED_palbociclib_total_counts$total_counts - X778_PALBO_EP_dexseq_ANNOTATED_palbociclib_low_MFE_TX_counts$total_counts_low_MFE_TX

# Construct the contingency table with raw counts
X778_PALBO_EP_dexseq_ANNOTATED_counts_contingency_table_low <- matrix(
  c(
    X778_PALBO_EP_dexseq_ANNOTATED_control_low_MFE_TX_counts$total_counts_low_MFE_TX, 
    X778_PALBO_EP_dexseq_ANNOTATED_palbociclib_low_MFE_TX_counts$total_counts_low_MFE_TX,
    X778_PALBO_EP_dexseq_ANNOTATED_control_total_counts$total_counts - X778_PALBO_EP_dexseq_ANNOTATED_control_low_MFE_TX_counts$total_counts_low_MFE_TX, 
    X778_PALBO_EP_dexseq_ANNOTATED_palbociclib_total_counts$total_counts - X778_PALBO_EP_dexseq_ANNOTATED_palbociclib_low_MFE_TX_counts$total_counts_low_MFE_TX
  ),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(Category = c("Low_MFE_TX_counts", "Other_TX_counts"),
                  Group = c("Control", "palbociclib"))
)

X778_PALBO_EP_dexseq_ANNOTATED_counts_contingency_table_df_low <- as.data.frame(X778_PALBO_EP_dexseq_ANNOTATED_counts_contingency_table_low)

X778_PALBO_EP_dexseq_ANNOTATED_fisher_result_counts_low <- fisher.test(X778_PALBO_EP_dexseq_ANNOTATED_counts_contingency_table_low)
print(X778_PALBO_EP_dexseq_ANNOTATED_fisher_result_counts_low)


cat("X778_PALBO_EP_dexseq_ANNOTATED Fisher P-value Fisher proportions of low MFE transcripts counts in drug-treated and control :", X778_PALBO_EP_dexseq_ANNOTATED_fisher_result_counts$p.value, "\n")  # This value goes in the summary table


# Output the results
cat("X778_PALBO_EP_dexseq_ANNOTATED control low MFE TX as fraction of total counts:", X778_PALBO_EP_dexseq_ANNOTATED_control_low_MFE_TX_counts_fraction_EP, "\n")
cat("X778_PALBO_EP_dexseq_ANNOTATED palbociclib low MFE TX as fraction of total counts:", X778_PALBO_EP_dexseq_ANNOTATED_palbociclib_low_MFE_TX_counts_fraction_EP, "\n")
cat("X778_PALBO_EP_dexseq_ANNOTATED net change counts low MFE TX drug vs control", X778_PALBO_EP_dexseq_ANNOTATED_palbociclib_low_MFE_TX_counts_fraction_EP-X778_PALBO_EP_dexseq_ANNOTATED_control_low_MFE_TX_counts_fraction_EP, "\n")





# Perform a two-proportion Z-test to compare the proportions of low MFE transcripts counts between the two groups
# The arguments are the counts of successes and the total number of trials for each group

# Number of low MFE TX counts in control and palbociclib groups
low_MFE_control <- X778_PALBO_EP_dexseq_ANNOTATED_control_low_MFE_TX_counts$total_counts_low_MFE_TX
low_MFE_palbociclib <- X778_PALBO_EP_dexseq_ANNOTATED_palbociclib_low_MFE_TX_counts$total_counts_low_MFE_TX

# Total counts in control and palbociclib groups
total_control <- X778_PALBO_EP_dexseq_ANNOTATED_control_total_counts$total_counts
total_palbociclib <- X778_PALBO_EP_dexseq_ANNOTATED_palbociclib_total_counts$total_counts

# Perform the two-proportion Z-test
X778_PALBO_EP_dexseq_ANNOTATED_z_test_result_low <- prop.test(
  x = c(low_MFE_control, low_MFE_palbociclib),    # successes: low MFE counts in both groups
  n = c(total_control, total_palbociclib)         # total counts in both groups
)

# Output the result of the two-proportion Z-test
X778_PALBO_EP_dexseq_ANNOTATED_z_test_count_25_p_value <- X778_PALBO_EP_dexseq_ANNOTATED_z_test_result_low$p.value













## X778_counts_EP heat capacity changes ####

# This is calculated on the total n. of transcripts with a significant log2 fold change
# Generate a scatter plot Log2_FC vs heat capacity
ggplot(X778_counts_EP, aes(y=Specific_heat, x=log2fold_treated_untreated)) +
  geom_point(size=1.5, alpha=0.1)+
  theme_classic() +
  theme(legend.position="none") +
  geom_vline(xintercept = 0, size=0.3) +
  xlim(-12, 12) +
  ylim(0,400) +
  theme(axis.title = element_text(size = 20)) +
  theme(plot.title = element_text(size = 15)) +
  labs(y="Heat capacity (J/C)",x="Log2_FC_expression",title="X778_PALBO_EP_dexseq_ANNOTATED_Log2_FC_adj.P<0.05")


# Show the data with a box plot
ggplot(X778_counts_EP, aes(x = factor(log2fold_treated_untreated > 0), y = Specific_heat)) +
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
  labs(y="Heat capacity (J/C)",x="",title="X778_PALBO_EP_dexseq_ANNOTATED_Log2_FC_adj.P<0.05 heat_capacity all transcripts") +
  scale_x_discrete(labels = c("log2FoldChange < 0", "log2FoldChange > 0"))

Mean_heat_capacity_X778_PALBO_EP_dexseq_ANNOTATED_UP <- mean(X778_PALBO_EP_dexseq_ANNOTATED_UP_EP$Specific_heat, na.rm = TRUE)
Mean_heat_capacity_X778_PALBO_EP_dexseq_ANNOTATED_DOWN <- mean(X778_PALBO_EP_dexseq_ANNOTATED_DOWN_EP$Specific_heat, na.rm = TRUE)
TTEST_heat_capacity_X778_PALBO_EP_dexseq_ANNOTATED_all_UP_vs_DOWN <- t.test(X778_PALBO_EP_dexseq_ANNOTATED_DOWN_EP$Specific_heat, X778_PALBO_EP_dexseq_ANNOTATED_UP_EP$Specific_heat, paired=FALSE)
print(TTEST_heat_capacity_X778_PALBO_EP_dexseq_ANNOTATED_all_UP_vs_DOWN)
# Add K-S test to look at the differences between the up-regukated and down-regulated genes
ks_test_X778_PALBO_EP_dexseq_ANNOTATED_all_UP_vs_DOWN_heat_capacity <- ks.test(X778_PALBO_EP_dexseq_ANNOTATED_DOWN_EP$Specific_heat, X778_PALBO_EP_dexseq_ANNOTATED_UP_EP$Specific_heat)
print(ks_test_X778_PALBO_EP_dexseq_ANNOTATED_all_UP_vs_DOWN_heat_capacity)




# Calculate the changes of heat capacity only in the high MFE transcripts

# Contribution of the UP or DOWN-regulated high MFE transcripts to the heat capacity
# Look at the heat capacity of the UP-regulated high MFE transcripts
X778_PALBO_EP_dexseq_ANNOTATED_high_MFE_TX_in_UP_heat_capacity <- mean(X778_PALBO_EP_dexseq_ANNOTATED_high_MFE_TX_in_UP$Specific_heat)

# Look at the heat capacity of the DOWN-regulated high MFE transcripts
X778_PALBO_EP_dexseq_ANNOTATED_high_MFE_TX_in_DOWN_heat_capacity <- mean(X778_PALBO_EP_dexseq_ANNOTATED_high_MFE_TX_in_DOWN$Specific_heat)

# Heat capacity gain in UP-regulated high MFE transcripts vs DOWN-regulated high MFE transcripts
X778_PALBO_EP_dexseq_ANNOTATED_heat_capacity_gain_UP_vs_DOWN <- X778_PALBO_EP_dexseq_ANNOTATED_high_MFE_TX_in_UP_heat_capacity - X778_PALBO_EP_dexseq_ANNOTATED_high_MFE_TX_in_DOWN_heat_capacity
cat("X778_PALBO_EP_dexseq_ANNOTATED net heat capacity UP vs DOWN high MFE TX:", X778_PALBO_EP_dexseq_ANNOTATED_heat_capacity_gain_UP_vs_DOWN, "\n")


# Extract specific heat data
X778_PALBO_EP_dexseq_ANNOTATED_Specific_heat_UP <- X778_PALBO_EP_dexseq_ANNOTATED_high_MFE_TX_in_UP$Specific_heat
X778_PALBO_EP_dexseq_ANNOTATED_Specific_heat_DOWN <- X778_PALBO_EP_dexseq_ANNOTATED_high_MFE_TX_in_DOWN$Specific_heat

# Combine data into a single data frame for easier plotting
X778_PALBO_EP_dexseq_ANNOTATED_Specific_heat_data <- data.frame(
  X778_PALBO_EP_dexseq_ANNOTATED_Specific_heat = c(X778_PALBO_EP_dexseq_ANNOTATED_Specific_heat_UP, X778_PALBO_EP_dexseq_ANNOTATED_Specific_heat_DOWN),
  regulation = factor(c(rep("UP", length(X778_PALBO_EP_dexseq_ANNOTATED_Specific_heat_UP)), rep("DOWN", length(X778_PALBO_EP_dexseq_ANNOTATED_Specific_heat_DOWN))))
)

# Generate box plot
ggplot(X778_PALBO_EP_dexseq_ANNOTATED_Specific_heat_data, aes(x = regulation, y = X778_PALBO_EP_dexseq_ANNOTATED_Specific_heat, fill = regulation)) +
  geom_boxplot(fill = "lightgrey", alpha = 0.7) +
  labs(title = "X778_PALBO_EP_dexseq_ANNOTATED specific Heat of High MFE Transcripts",
       x = "",
       y = "Specific Heat (J/C)") +
  scale_x_discrete(labels = c("log2_FC < 0", "log2_FC > 0")) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 16),
        axis.title = element_text(size = 25))

Mean_heat_capacity_X778_PALBO_EP_dexseq_ANNOTATED_high_MFE_TX_UP <- mean(X778_PALBO_EP_dexseq_ANNOTATED_high_MFE_TX_in_UP$Specific_heat)
Mean_heat_capacity_X778_PALBO_EP_dexseq_ANNOTATED_high_MFE_TX_DOWN <- mean(X778_PALBO_EP_dexseq_ANNOTATED_high_MFE_TX_in_DOWN$Specific_heat)
TTEST_heat_capacity_X778_PALBO_EP_dexseq_ANNOTATED_all_high_MFE_TX_UP_vs_DOWN <- t.test(X778_PALBO_EP_dexseq_ANNOTATED_high_MFE_TX_in_DOWN$Specific_heat, X778_PALBO_EP_dexseq_ANNOTATED_high_MFE_TX_in_UP$Specific_heat, paired=FALSE)
print(TTEST_heat_capacity_X778_PALBO_EP_dexseq_ANNOTATED_all_high_MFE_TX_UP_vs_DOWN)

# Add K-S test to look at the differences between the up-regukated and down-regulated genes
ks_test_heat_capacity_X778_PALBO_EP_dexseq_ANNOTATED_all_high_MFE_TX_UP_vs_DOWN <- ks.test(X778_PALBO_EP_dexseq_ANNOTATED_high_MFE_TX_in_DOWN$Specific_heat, X778_PALBO_EP_dexseq_ANNOTATED_high_MFE_TX_in_UP$Specific_heat)
print(ks_test_heat_capacity_X778_PALBO_EP_dexseq_ANNOTATED_all_high_MFE_TX_UP_vs_DOWN)










# Analysis cell proliferation Nanopore data 778 PALBO ####
# It is important to demonstrate that cells are not proliferating 72h after drug treatment
# I also add the DNA repair pathways as positive control. We know that some of them are down-regulated when cells undergo cell cycle arrest.
# Open the file with the gene sets
source("path_to/All_gene_sets_AC.R")

# Open the DESeq2 file with the log2 fold changes at the gene level
X778_PALBO_DESeq2 <- read.delim(file="path_to/early.phase_778_palbociclib.tx_counts.both_conditions.treated_untreated_deseq2.tab", header = TRUE, sep = "\t", quote = "")

# Filter the proliferation gene sets by using a function
# Use a function also to determine the mean and standard deviation of each set and the plot all of them together

# Generate the function that filters each gene set in the target data frame
# This function will be used in the following commands
filtering <- function(x) {
  Filtering <- dplyr::filter(X778_PALBO_DESeq2, gene_name %in% x)
  return(Filtering)
}

# Generate the list containing all the gene sets to be used for the filtering.
CURATED_PATHWAYS.input <- list(MTOR_positive_regulators, MTOR_negative_regulators, Proliferation_Science_paper, Meta_PCNA, Proliferation_GO_genes_Magill_Sack_2018, Proliferation_STOP_genes_Magill_Sack_2018, Reactome_G1_phase, Proliferation_Carey_2021, Essential_cell_cycle_Viner_Breuer_Cell_Division_2019, Essential_S_Phase_Viner_Breuer_Cell_Division_2019, Homologous_recombination, Fanconi_anemia, High_fidelity_polymerases)

# This creates the names that will be used for the final data frames of each file resulting from the filtering of the target data frame with each gene set
CURATED_PATHWAYS.output <- c("X778_PALBO_MTOR_POS_EP", "X778_PALBO_MTOR_NEG_EP", "X778_PALBO_Proliferation_Science_EP", "X778_PALBO_Meta_PCNA_EP", "X778_PALBO_Proliferation_GO_EP", "X778_PALBO_Proliferation_STOP_EP", "X778_PALBO_G1_phase_EP", "X778_PALBO_Proliferation_Carey_EP", "X778_PALBO_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019_EP", "X778_PALBO_Essential_S_Phase_Viner_Breuer_Cell_Division_2019_EP", "X778_PALBO_HR_EP", "X778_PALBO_Fanconi_anemia_EP", "X778_PALBO_HF_POL_EP")

# Assign the names of the gene sets to each entry in the list
names(CURATED_PATHWAYS.input) <- CURATED_PATHWAYS.output

# Apply the filtering function to create a list of data frames
Curated_test_df <- lapply(CURATED_PATHWAYS.input, function(x) filtering(x))

# creating a vector with desired names of filtered data frames
Curated_dfs.output <- c("df_X778_PALBO_MTOR_POS_EP", "df_X778_PALBO_MTOR_NEG_EP", "df_X778_PALBO_Proliferation_Science_EP", "df_X778_PALBO_Meta_PCNA_EP", "df_X778_PALBO_Proliferation_GO_EP", "df_X778_PALBO_Proliferation_STOP_EP", "df_X778_PALBO_G1_phase_EP", "df_X778_PALBO_Proliferation_Carey_EP", "df_X778_PALBO_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019_EP", "df_X778_PALBO_Essential_S_Phase_Viner_Breuer_Cell_Division_2019_EP", "df_X778_PALBO_HR_EP", "df_X778_PALBO_Fanconi_anemia_EP", "df_X778_PALBO_HF_POL_EP")

# loop to assign data frames to the desired variable name
for (i in seq_along(Curated_dfs.output)) {
  assign(Curated_dfs.output[i], Curated_test_df[[i]])
}





# Generate the function to calculate the MEAN of the log2FoldChange of each gene set resulting from the filtering above.
Mean_curated <- function(x) {
  mean_curated <- mean(x$log2FoldChange, na.rm = TRUE)
  return(mean_curated)
}

# Create an input of the data frames generated above to calculate the mean of log2FoldChange
Curated_input_dfs <- list(df_X778_PALBO_MTOR_POS_EP, df_X778_PALBO_MTOR_NEG_EP, df_X778_PALBO_Proliferation_Science_EP, df_X778_PALBO_Meta_PCNA_EP, df_X778_PALBO_Proliferation_GO_EP, df_X778_PALBO_Proliferation_STOP_EP, df_X778_PALBO_G1_phase_EP, df_X778_PALBO_Proliferation_Carey_EP, df_X778_PALBO_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019_EP, df_X778_PALBO_Essential_S_Phase_Viner_Breuer_Cell_Division_2019_EP, df_X778_PALBO_HR_EP, df_X778_PALBO_Fanconi_anemia_EP, df_X778_PALBO_HF_POL_EP)

# Apply the mean function to create a list of means
Curated_mean_test_df <- lapply(Curated_input_dfs, function(x) Mean_curated(x))

# creating a vector with desired names of the files containing the mean of each gene set log2FoldChange
Curated_Mean_X778_PALBO_dfs.output <- c("Mean_X778_PALBO_MTOR_POS_EP", "Mean_X778_PALBO_MTOR_NEG_EP", "Mean_X778_PALBO_Proliferation_Science_EP", "Mean_X778_PALBO_Meta_PCNA_EP", "Mean_X778_PALBO_Proliferation_GO_EP", "Mean_X778_PALBO_Proliferation_STOP_EP", "Mean_X778_PALBO_G1_phase_EP", "Mean_X778_PALBO_Proliferation_Carey_EP", "Mean_X778_PALBO_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019_EP", "Mean_X778_PALBO_Essential_S_Phase_Viner_Breuer_Cell_Division_2019_EP", "Mean_X778_PALBO_HR_EP", "Mean_X778_PALBO_Fanconi_anemia_EP", "Mean_X778_PALBO_HF_POL_EP")

# loop to assign means to the desired variable name
for (i in seq_along(Curated_Mean_X778_PALBO_dfs.output)) {
  assign(Curated_Mean_X778_PALBO_dfs.output[i], Curated_mean_test_df[[i]])
}




# Generate the function to calculate the STANDARD DEVIATION of the log2FoldChange of each gene set resulting from the filtering above.
STDEV_curated <- function(x) {
  stdev_curated <- sd(x$log2FoldChange, na.rm = TRUE)
  return(stdev_curated)
}

# Apply the standard deviation function to create a list of standard deviations
Curated_stdev_test_df <- lapply(Curated_input_dfs, function(x) STDEV_curated(x))

# creating a vector with desired names of the files containing the standard deviation of each gene set log2FoldChange
Curated_STDEV_X778_PALBO_dfs.output <- c("STDEV_X778_PALBO_MTOR_POS_EP", "STDEV_X778_PALBO_MTOR_NEG_EP", "STDEV_X778_PALBO_Proliferation_Science_EP", "STDEV_X778_PALBO_Meta_PCNA_EP", "STDEV_X778_PALBO_Proliferation_GO_EP", "STDEV_X778_PALBO_Proliferation_STOP_EP", "STDEV_X778_PALBO_G1_phase_EP", "STDEV_X778_PALBO_Proliferation_Carey_EP", "STDEV_X778_PALBO_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019_EP", "STDEV_X778_PALBO_Essential_S_Phase_Viner_Breuer_Cell_Division_2019_EP", "STDEV_X778_PALBO_HR_EP", "STDEV_X778_PALBO_Fanconi_anemia_EP", "STDEV_X778_PALBO_HF_POL_EP")
# loop to assign data.frames to the desired variable name
for (i in (seq_along(Curated_STDEV_X778_PALBO_dfs.output))) {
  do.call("<-",list(Curated_STDEV_X778_PALBO_dfs.output[i], Curated_stdev_test_df[[i]]))
}




# Generate the data frame for the plot

# Generate the data frame for the plot
df_778_PALBO_proliferation <- data.frame(mean=c(Mean_X778_PALBO_MTOR_POS_EP, Mean_X778_PALBO_MTOR_NEG_EP, Mean_X778_PALBO_Proliferation_Science_EP, Mean_X778_PALBO_Meta_PCNA_EP, Mean_X778_PALBO_Proliferation_GO_EP, Mean_X778_PALBO_Proliferation_STOP_EP, Mean_X778_PALBO_G1_phase_EP, Mean_X778_PALBO_Proliferation_Carey_EP, Mean_X778_PALBO_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019_EP, Mean_X778_PALBO_Essential_S_Phase_Viner_Breuer_Cell_Division_2019_EP, Mean_X778_PALBO_HR_EP, Mean_X778_PALBO_Fanconi_anemia_EP, Mean_X778_PALBO_HF_POL_EP),
                                   sd=c(STDEV_X778_PALBO_MTOR_POS_EP, STDEV_X778_PALBO_MTOR_NEG_EP, STDEV_X778_PALBO_Proliferation_Science_EP, STDEV_X778_PALBO_Meta_PCNA_EP, STDEV_X778_PALBO_Proliferation_GO_EP, STDEV_X778_PALBO_Proliferation_STOP_EP, STDEV_X778_PALBO_G1_phase_EP, STDEV_X778_PALBO_Proliferation_Carey_EP, STDEV_X778_PALBO_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019_EP, STDEV_X778_PALBO_Essential_S_Phase_Viner_Breuer_Cell_Division_2019_EP, STDEV_X778_PALBO_HR_EP, STDEV_X778_PALBO_Fanconi_anemia_EP, STDEV_X778_PALBO_HF_POL_EP),
                                   Category=as.factor(c("X778_PALBO_MTOR_POS_EP", "X778_PALBO_MTOR_NEG_EP", "X778_PALBO_Proliferation_Science_EP", "X778_PALBO_Meta_PCNA_EP", "X778_PALBO_Proliferation_GO_EP", "X778_PALBO_Proliferation_STOP_EP", "X778_PALBO_G1_phase_EP", "X778_PALBO_Proliferation_Carey_EP", "X778_PALBO_Essential_cell_cycle_EP", "X778_PALBO_Essential_S_Phase_EP", "X778_PALBO_HR_EP", "X778_PALBO_FA_EP", "X778_PALBO_HF_POL_EP")))
df_778_PALBO_proliferation




library(ggplot2)
ggplot(df_778_PALBO_proliferation, aes(x=reorder(Category, -mean), y=mean)) +
  theme_bw() +
  geom_point(position=position_dodge(width = 1), stat="identity",colour='black') +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.4) +
  theme(axis.text.x = element_text(size=12, angle=45,hjust=0.95,vjust=0.98)) +
  theme(axis.text.y = element_text(size=14)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 10))) +
  geom_point(colour = "black", size = 3) +
  xlab("") + ylab("Relative log2 fold change") +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  ggtitle("Analysis cell proliferation 778 PALBO NANOPORE 72h vs T0") +
  scale_x_discrete(limits = c("X778_PALBO_Proliferation_Science_EP", "X778_PALBO_Meta_PCNA_EP", "X778_PALBO_Proliferation_GO_EP", "X778_PALBO_Proliferation_STOP_EP", "X778_PALBO_Proliferation_Carey_EP", "X778_PALBO_G1_phase_EP", "X778_PALBO_MTOR_POS_EP", "X778_PALBO_MTOR_NEG_EP", "X778_PALBO_Essential_cell_cycle_EP", "X778_PALBO_Essential_S_Phase_EP", "X778_PALBO_HR_EP", "X778_PALBO_FA_EP", "X778_PALBO_HF_POL_EP")) +
  theme(plot.title = element_text(color="black", size=14, face="bold"))


# Merge all the data frames to generate box plots
# Add a column with the specific name of the gene set
df_X778_PALBO_MTOR_POS_EP$Gene_set <- "X778_PALBO_MTOR_POS_EP"
df_X778_PALBO_MTOR_NEG_EP$Gene_set <- "X778_PALBO_MTOR_NEG_EP"
df_X778_PALBO_Proliferation_Science_EP$Gene_set <- "X778_PALBO_Proliferation_Science_EP"
df_X778_PALBO_Meta_PCNA_EP$Gene_set <- "X778_PALBO_Meta_PCNA_EP"
df_X778_PALBO_Proliferation_GO_EP$Gene_set <- "X778_PALBO_Proliferation_GO_EP"
df_X778_PALBO_Proliferation_STOP_EP$Gene_set <- "X778_PALBO_Proliferation_STOP_EP"
df_X778_PALBO_G1_phase_EP$Gene_set <- "X778_PALBO_G1_phase_EP"
df_X778_PALBO_Proliferation_Carey_EP$Gene_set <- "X778_PALBO_Proliferation_Carey_EP"
df_X778_PALBO_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019_EP$Gene_set <- "X778_PALBO_Essential_cell_cycle_EP"
df_X778_PALBO_Essential_S_Phase_Viner_Breuer_Cell_Division_2019_EP$Gene_set <- "X778_PALBO_Essential_S_Phase_EP"
df_X778_PALBO_HR_EP$Gene_set <- "X778_PALBO_HR_EP"
df_X778_PALBO_Fanconi_anemia_EP$Gene_set <- "X778_PALBO_Fanconi_anemia_EP"
df_X778_PALBO_HF_POL_EP$Gene_set <- "X778_PALBO_HF_POL_EP"



# Combine the data frames
combined_df <- rbind(
  df_X778_PALBO_MTOR_POS_EP,
  df_X778_PALBO_MTOR_NEG_EP,
  df_X778_PALBO_Proliferation_Science_EP,
  df_X778_PALBO_Meta_PCNA_EP,
  df_X778_PALBO_Proliferation_GO_EP,
  df_X778_PALBO_Proliferation_STOP_EP,
  df_X778_PALBO_G1_phase_EP,
  df_X778_PALBO_Proliferation_Carey_EP,
  df_X778_PALBO_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019_EP,
  df_X778_PALBO_Essential_S_Phase_Viner_Breuer_Cell_Division_2019_EP,
  df_X778_PALBO_HR_EP,
  df_X778_PALBO_Fanconi_anemia_EP,
  df_X778_PALBO_HF_POL_EP
)

# Ensure the Gene_set column is a factor with the correct order
combined_df$Gene_set <- factor(combined_df$Gene_set, levels = c("X778_PALBO_Proliferation_Science_EP", "X778_PALBO_Meta_PCNA_EP", "X778_PALBO_Proliferation_GO_EP", "X778_PALBO_Proliferation_STOP_EP", "X778_PALBO_Proliferation_Carey_EP", "X778_PALBO_G1_phase_EP", "X778_PALBO_MTOR_POS_EP", "X778_PALBO_MTOR_NEG_EP", "X778_PALBO_Essential_cell_cycle_EP", "X778_PALBO_Essential_S_Phase_EP", "X778_PALBO_HR_EP", "X778_PALBO_Fanconi_anemia_EP", "X778_PALBO_HF_POL_EP"))


# Generate the box plot
library(ggplot2)
ggplot(combined_df, aes(x=Gene_set, y=log2FoldChange)) +
  geom_boxplot(fill="lightgrey") + 
  stat_summary(fun.y=mean, geom="point", shape=20, size=4, color="blue") +
  # geom_jitter(width=0.2, alpha=0.5) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle=45, hjust=0.95, vjust=0.98)) +
  theme(axis.text.y = element_text(size=14)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 10))) +
  xlab("") + ylab("Relative log2 fold change") +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  ggtitle("Analysis cell proliferation 778 palbociclib early phase") +
  scale_x_discrete(limits = levels(combined_df$Gene_set)) +
  theme(plot.title = element_text(color="black", size=14, face="bold")) 






Mean_df_MTOR_POS_X778_PALBO_EP <- mean(df_X778_PALBO_MTOR_POS_EP$log2FoldChange)
Mean_df_MTOR_NEG_X778_PALBO_EP <- mean(df_X778_PALBO_MTOR_NEG_EP$log2FoldChange)
Mean_df_Proliferation_Science_X778_PALBO_EP <- mean(df_X778_PALBO_Proliferation_Science_EP$log2FoldChange)
Mean_df_Meta_PCNA_X778_PALBO_EP <- mean(df_X778_PALBO_Meta_PCNA_EP$log2FoldChange)
Mean_df_Proliferation_GO_X778_PALBO_EP <- mean(df_X778_PALBO_Proliferation_GO_EP$log2FoldChange)
Mean_df_Proliferation_STOP_X778_PALBO_EP <- mean(df_X778_PALBO_Proliferation_STOP_EP$log2FoldChange)
Mean_df_G1_phase_X778_PALBO_EP <- mean(df_X778_PALBO_G1_phase_EP$log2FoldChange)
Mean_df_Proliferation_Carey_X778_PALBO_EP <- mean(df_X778_PALBO_Proliferation_Carey_EP$log2FoldChange)
Mean_df_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019_X778_PALBO_EP <- mean(df_X778_PALBO_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019_EP$log2FoldChange)
Mean_df_Essential_S_Phase_Viner_Breuer_Cell_Division_2019_X778_PALBO_EP <- mean(df_X778_PALBO_Essential_S_Phase_Viner_Breuer_Cell_Division_2019_EP$log2FoldChange)
Mean_df_HR_X778_PALBO_EP <- mean(df_X778_PALBO_HR_EP$log2FoldChange)
Mean_df_FA_X778_PALBO_EP <- mean(df_X778_PALBO_Fanconi_anemia_EP$log2FoldChange)
Mean_df_HF_POL_X778_PALBO_EP <- mean(df_X778_PALBO_HF_POL_EP$log2FoldChange)

# Perform the Wilcoxon signed-rank test
wilcox_test_MTOR_POS_X778_PALBO_EP <- wilcox.test(df_X778_PALBO_MTOR_POS_EP$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_MTOR_NEG_X778_PALBO_EP <- wilcox.test(df_X778_PALBO_MTOR_NEG_EP$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_Proliferation_Science_X778_PALBO_EP <- wilcox.test(df_X778_PALBO_Proliferation_Science_EP$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_Meta_PCNA_X778_PALBO_EP <- wilcox.test(df_X778_PALBO_Meta_PCNA_EP$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_Proliferation_GO_X778_PALBO_EP <- wilcox.test(df_X778_PALBO_Proliferation_GO_EP$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_Proliferation_STOP_X778_PALBO_EP <- wilcox.test(df_X778_PALBO_Proliferation_STOP_EP$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_G1_phase_X778_PALBO_EP <- wilcox.test(df_X778_PALBO_G1_phase_EP$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_Proliferation_Carey_X778_PALBO_EP <- wilcox.test(df_X778_PALBO_Proliferation_Carey_EP$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019_X778_PALBO_EP <- wilcox.test(df_X778_PALBO_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019_EP$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_Essential_S_Phase_Viner_Breuer_Cell_Division_2019_X778_PALBO_EP <- wilcox.test(df_X778_PALBO_Essential_S_Phase_Viner_Breuer_Cell_Division_2019_EP$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_HR_X778_PALBO_EP <- wilcox.test(df_X778_PALBO_HR_EP$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_FA_X778_PALBO_EP <- wilcox.test(df_X778_PALBO_Fanconi_anemia_EP$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_HF_POL_X778_PALBO_EP <- wilcox.test(df_X778_PALBO_HF_POL_EP$log2FoldChange, mu = 0, alternative = "less")





































# SKBR3 cell line ANNOTATED transcripts ####
SKBR3_counts_EP <- dplyr::inner_join(SKBR3_counts_EP_dexseq, Table_ANNOTATED, by="ensembl_transcript_id_version") %>%
  dplyr::filter(padj < 0.05) # The file 'Table' contains the Normalized_MFE values of all ANNOTATED transcripts


# Generate a scatter plot Log2_FC vs Normalized_MFE
library(ggplot2)
ggplot(SKBR3_counts_EP, aes(y=Normalized_MFE, x=log2fold_treated_untreated)) +
  geom_point(size=1.5, alpha=0.3)+
  theme_classic() +
  theme(legend.position="none") +
  geom_vline(xintercept = 0, size=0.3) +
  #xlim(-5.5,5.5) +
  theme(axis.title = element_text(size = 20)) +
  theme(plot.title = element_text(size = 15)) +
  labs(y="Normalized_MFE",x="Log2_FC_expression",title="SKBR3_shared_ANNOTATED_transcripts_72h_Log2_FC_adj.P<0.05")


# It is well known that ANNOTATED transcripts, in general have lower MFEs than the UNANNOTATED transcripts

SKBR3_EP_dexseq_ANNOTATED_UP_EP <- dplyr::filter(SKBR3_counts_EP, log2fold_treated_untreated > 0)
Mean_SKBR3_EP_dexseq_ANNOTATED_UP_EP <- mean(SKBR3_EP_dexseq_ANNOTATED_UP_EP$Normalized_MFE)

SKBR3_EP_dexseq_ANNOTATED_DOWN_EP <- dplyr::filter(SKBR3_counts_EP, log2fold_treated_untreated < 0)
Mean_SKBR3_EP_dexseq_ANNOTATED_DOWN_EP <- mean(SKBR3_EP_dexseq_ANNOTATED_DOWN_EP$Normalized_MFE)

t.test(SKBR3_EP_dexseq_ANNOTATED_DOWN_EP$Normalized_MFE, SKBR3_EP_dexseq_ANNOTATED_UP_EP$Normalized_MFE, paired=FALSE)

# There is a significant increase of the mean MFE of the UP-regulated vs the DOWN-regulated ANNOTATED transcripts
# p-value = 0.002695



# Show the data with a violin plot
ggplot(SKBR3_counts_EP, aes(x = factor(log2fold_treated_untreated > 0), y = Normalized_MFE)) +
  geom_violin(fill = "lightgrey", alpha = 0.7) +
  geom_jitter(width = 0.1, size = 1, alpha = 0.3) +
  theme_classic() +
  theme(legend.position = "none", axis.title = element_text(size = 20), plot.title = element_text(size = 10),
        axis.text.x = element_text(size = 16, color = "black")) +  # Increase font size and make x-axis labels black
  geom_vline(xintercept = 0, size = 0.3) +
  labs(y = "Normalized_MFE", x = NULL, title = "SKBR3_shared_ANNOTATED_transcripts_72h_Log2_FC_adj.P<0.05") +
  scale_x_discrete(labels = c("Log2_FC < 0", "Log2_FC > 0"))
# In this modified code, the x aesthetic is set as a factor variable factor(log2fold_treated_untreated > 0). 
# This creates two separate groups based on the condition Log2_FC > 0 or Log2_FC < 0. Each group is represented by a different violin shape.



# Look at the differences between the 0.75 quantile of UP- and DOWN-regulated transcripts
quantile_075_SKBR3_UP_EP <- quantile(SKBR3_EP_dexseq_ANNOTATED_UP_EP$Normalized_MFE, probs = c(0.75)) 
quantile_075_SKBR3_DOWN_EP <- quantile(SKBR3_EP_dexseq_ANNOTATED_DOWN_EP$Normalized_MFE, probs = c(0.75)) 
ks_test_SKBR3_075 <- ks.test(SKBR3_EP_dexseq_ANNOTATED_UP_EP$Normalized_MFE[SKBR3_EP_dexseq_ANNOTATED_UP_EP$Normalized_MFE >= quantile_075_SKBR3_UP_EP], SKBR3_EP_dexseq_ANNOTATED_DOWN_EP$Normalized_MFE[SKBR3_EP_dexseq_ANNOTATED_DOWN_EP$Normalized_MFE >= quantile_075_SKBR3_DOWN_EP])
print(ks_test_SKBR3_075)

# Look at the differences between the 0.25 quantile of UP- and DOWN-regulated transcripts
quantile_025_SKBR3_UP_EP <- quantile(SKBR3_EP_dexseq_ANNOTATED_UP_EP$Normalized_MFE, probs = c(0.25)) 
quantile_025_SKBR3_DOWN_EP <- quantile(SKBR3_EP_dexseq_ANNOTATED_DOWN_EP$Normalized_MFE, probs = c(0.25)) 
ks_test_SKBR3_025 <- ks.test(SKBR3_EP_dexseq_ANNOTATED_UP_EP$Normalized_MFE[SKBR3_EP_dexseq_ANNOTATED_UP_EP$Normalized_MFE <= quantile_025_SKBR3_UP_EP], SKBR3_EP_dexseq_ANNOTATED_DOWN_EP$Normalized_MFE[SKBR3_EP_dexseq_ANNOTATED_DOWN_EP$Normalized_MFE <= quantile_025_SKBR3_DOWN_EP])
print(ks_test_SKBR3_025)



# Generate a data frame with all the CUB values and plot the data in a density plot
max_length_SKBR3_EP_MFE <- max(c(length(SKBR3_EP_dexseq_ANNOTATED_UP_EP$Normalized_MFE), length(SKBR3_EP_dexseq_ANNOTATED_DOWN_EP$Normalized_MFE)))

data_SKBR3_EP_MFE_density <- data.frame(MFE_TX_UP_EP_SKBR3 = c(SKBR3_EP_dexseq_ANNOTATED_UP_EP$Normalized_MFE, rep(NA, max_length_SKBR3_EP_MFE - length(SKBR3_EP_dexseq_ANNOTATED_UP_EP$Normalized_MFE))), # This is useful when we need to create a data frame and the number of values between columns are different
                                     MFE_TX_DOWN_EP_SKBR3 = c(SKBR3_EP_dexseq_ANNOTATED_DOWN_EP$Normalized_MFE, rep(NA, max_length_SKBR3_EP_MFE- length(SKBR3_EP_dexseq_ANNOTATED_DOWN_EP$Normalized_MFE))))
library(reshape)
data_SKBR3_EP_MFE_density <- melt(data_SKBR3_EP_MFE_density)

library(ggplot2)
ggplot(data= data_SKBR3_EP_MFE_density, aes(x=value, group=variable, color=variable)) +
  geom_density (size=1) +
  theme_classic() +
  #xlim(0,1.5) +
  labs(title = "MFE_distribution_ANNOTATED_Transcripts_SKBR3_72h")+
  labs(y="Density")+
  labs(x="DGfold / nt (kcal/mol)")










## Analysis net MFE and heat capacity ####
# Determine whether within the differentially expressed genes there is a significant up-regulation of high MFE transcripts:
# Define high MFE transcripts by taking the 75th percentile of the combined MFE distributions of control and drug-treated samples.
# Count the number of transcripts in the significantly up- and down-regulated genes that have MFE values above this threshold.
# Perform a Fisher’s exact test to compare the proportions of high MFE transcripts between the two groups.


SKBR3_EP_dexseq_ANNOTATED_norm_TX_counts_MFE_df <- SKBR3_counts_EP %>%
  dplyr::select(ensembl_transcript_id_version, ensembl_transcript_id, ensembl_gene_id, external_gene_name, transcript_is_canonical, Transcript_length,
                SKBR3_untreated_R1, SKBR3_untreated_R2, SKBR3_untreated_R3,
                SKBR3_lapatinib_R1, SKBR3_lapatinib_R2, SKBR3_lapatinib_R3,
                Normalized_MFE, Temperature, Specific_heat) 



# Identify the 75th percentile of the combined MFE distributions of control and drug-treated samples
quantile_075_SKBR3_EP_dexseq_ANNOTATED_norm_TX_counts_MFE_df <- quantile(SKBR3_EP_dexseq_ANNOTATED_norm_TX_counts_MFE_df$Normalized_MFE, probs = c(0.75))

# Select all the transcripts with Normalized_MFE values above the 75th percentile
SKBR3_EP_dexseq_ANNOTATED_high_MFE_TX <- dplyr::filter(SKBR3_EP_dexseq_ANNOTATED_norm_TX_counts_MFE_df, Normalized_MFE > quantile_075_SKBR3_EP_dexseq_ANNOTATED_norm_TX_counts_MFE_df)

# Count the number of transcripts in the significantly up- and down-regulated pools that have MFE values above this threshold.
SKBR3_EP_dexseq_ANNOTATED_high_MFE_TX_in_UP <- dplyr::filter(SKBR3_EP_dexseq_ANNOTATED_UP_EP, ensembl_transcript_id_version %in% SKBR3_EP_dexseq_ANNOTATED_high_MFE_TX$ensembl_transcript_id_version)
SKBR3_EP_dexseq_ANNOTATED_high_MFE_TX_in_DOWN <- dplyr::filter(SKBR3_EP_dexseq_ANNOTATED_DOWN_EP, ensembl_transcript_id_version %in% SKBR3_EP_dexseq_ANNOTATED_high_MFE_TX$ensembl_transcript_id_version)


# Perform a Fisher’s exact test to compare the proportions of high MFE transcripts between the two groups.
SKBR3_EP_dexseq_ANNOTATED_total_UP <- nrow(SKBR3_EP_dexseq_ANNOTATED_UP_EP)
SKBR3_EP_dexseq_ANNOTATED_total_DOWN <- nrow(SKBR3_EP_dexseq_ANNOTATED_DOWN_EP)
SKBR3_EP_dexseq_ANNOTATED_count_high_MFE_UP <- nrow(SKBR3_EP_dexseq_ANNOTATED_high_MFE_TX_in_UP)
SKBR3_EP_dexseq_ANNOTATED_count_high_MFE_DOWN <- nrow(SKBR3_EP_dexseq_ANNOTATED_high_MFE_TX_in_DOWN)

SKBR3_EP_dexseq_ANNOTATED_contingency_table <- matrix(c(SKBR3_EP_dexseq_ANNOTATED_count_high_MFE_UP, SKBR3_EP_dexseq_ANNOTATED_total_UP - SKBR3_EP_dexseq_ANNOTATED_count_high_MFE_UP,
                                                           SKBR3_EP_dexseq_ANNOTATED_count_high_MFE_DOWN, SKBR3_EP_dexseq_ANNOTATED_total_DOWN - SKBR3_EP_dexseq_ANNOTATED_count_high_MFE_DOWN),
                                                         nrow = 2, byrow = TRUE,
                                                         dimnames = list(UP = c("UP-regulated", "DOWN-regulated"),
                                                                         DOWN = c("High MFE TX", "Not High MFE TX")))

SKBR3_EP_dexseq_ANNOTATED_contingency_table_df <- as.data.frame(SKBR3_EP_dexseq_ANNOTATED_contingency_table)

SKBR3_EP_dexseq_ANNOTATED_fisher_result <- fisher.test(SKBR3_EP_dexseq_ANNOTATED_contingency_table)

cat("SKBR3_counts_EP P-value from Fisher test:", SKBR3_EP_dexseq_ANNOTATED_fisher_result$p.value, "\n")


# Look at the fractions
SKBR3_EP_dexseq_ANNOTATED_high_MFE_in_UP_fraction <- SKBR3_EP_dexseq_ANNOTATED_count_high_MFE_UP / SKBR3_EP_dexseq_ANNOTATED_total_UP
SKBR3_EP_dexseq_ANNOTATED_high_MFE_in_DOWN_fraction <- SKBR3_EP_dexseq_ANNOTATED_count_high_MFE_DOWN / SKBR3_EP_dexseq_ANNOTATED_total_DOWN

# Output the results
cat("Fraction of SKBR3_counts_EP high MFE TX in UP-regulated genes:", SKBR3_EP_dexseq_ANNOTATED_high_MFE_in_UP_fraction, "\n") # This value goes in the summary table
cat("Fraction of SKBR3_counts_EP high MFE TX in DOWN-regulated genes:", SKBR3_EP_dexseq_ANNOTATED_high_MFE_in_DOWN_fraction, "\n") # This value goes in the summary table
cat("Mean of SKBR3_counts_EP high MFE TX in UP-regulated genes:", mean(SKBR3_EP_dexseq_ANNOTATED_high_MFE_TX_in_UP$Normalized_MFE), "\n")
cat("Mean of SKBR3_counts_EP high MFE TX in DOWN-regulated genes:", mean(SKBR3_EP_dexseq_ANNOTATED_high_MFE_TX_in_DOWN$Normalized_MFE), "\n")
cat("Mode of SKBR3_counts_EP high MFE TX in UP-regulated genes:", calculate_mode(SKBR3_EP_dexseq_ANNOTATED_high_MFE_TX_in_UP$Normalized_MFE), "\n")
cat("Mode of SKBR3_counts_EP high MFE TX in DOWN-regulated genes:", calculate_mode(SKBR3_EP_dexseq_ANNOTATED_high_MFE_TX_in_DOWN$Normalized_MFE), "\n")




# Perform a two-proportion z-test to compare proportions
# Calculate the proportions and total counts
SKBR3_EP_dexseq_ANNOTATED_high_MFE_in_UP_proportion <- SKBR3_EP_dexseq_ANNOTATED_high_MFE_in_UP_fraction
SKBR3_EP_dexseq_ANNOTATED_high_MFE_in_DOWN_proportion <- SKBR3_EP_dexseq_ANNOTATED_high_MFE_in_DOWN_fraction

# Combined proportion for pooled standard error
combined_proportion <- (SKBR3_EP_dexseq_ANNOTATED_count_high_MFE_UP + SKBR3_EP_dexseq_ANNOTATED_count_high_MFE_DOWN) /
  (SKBR3_EP_dexseq_ANNOTATED_total_UP + SKBR3_EP_dexseq_ANNOTATED_total_DOWN)

# Standard error of the difference in proportions
standard_error <- sqrt(
  combined_proportion * (1 - combined_proportion) * 
    (1 / SKBR3_EP_dexseq_ANNOTATED_total_UP + 1 / SKBR3_EP_dexseq_ANNOTATED_total_DOWN)
)

# Z-statistic
z_statistic <- (SKBR3_EP_dexseq_ANNOTATED_high_MFE_in_UP_proportion - SKBR3_EP_dexseq_ANNOTATED_high_MFE_in_DOWN_proportion) /
  standard_error

# Two-tailed p-value
SKBR3_EP_dexseq_ANNOTATED_z_test_MFE_p_value <- 2 * (1 - pnorm(abs(z_statistic)))








# Show data as violin plot
# Add a new column to indicate if the value is above the 75th percentile
SKBR3_counts_EP$above_quantile_075 <- SKBR3_counts_EP$Normalized_MFE > quantile_075_SKBR3_EP_dexseq_ANNOTATED_norm_TX_counts_MFE_df

# Generate violin plot with conditional jitter coloring
ggplot(SKBR3_counts_EP, aes(x = factor(log2fold_treated_untreated > 0), y = Normalized_MFE)) +
  geom_violin(fill = "white") +
  geom_jitter(aes(color = above_quantile_075), width = 0.1, size = 1, alpha = 0.3) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "blue")) +
  theme_classic() +
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
    title = "SKBR3_EP_dexseq_ANNOTATED_Log2_FC_adj.P<0.05"
  ) +
  scale_x_discrete(labels = c("log2FoldChange < 0", "log2FoldChange > 0"))













# Determine the impact of these changes in terms of fraction of total counts:
# Calculate the fraction of total counts of high MFE transcripts in control and drug-treated groups. 
# This is done by summing the counts of the high MFE transcripts and dividing by the total counts in each group.
# Fisher’s exact test to compare the proportions of high MFE transcripts counts between the two groups.


# Convert from wide to long format
SKBR3_EP_dexseq_ANNOTATED_long_data <- SKBR3_EP_dexseq_ANNOTATED_norm_TX_counts_MFE_df %>%
  pivot_longer(
    cols = starts_with("SKBR3"),      # Select all columns starting with "SKBR3"
    names_to = "replicate",         # Name of the new column for replicate identifiers
    values_to = "counts"            # Name of the new column for the counts
  )

# Add a condition column based on the replicate
SKBR3_EP_dexseq_ANNOTATED_long_data <- SKBR3_EP_dexseq_ANNOTATED_long_data %>%
  mutate(condition = ifelse(replicate %in% c("SKBR3_untreated_R1", "SKBR3_untreated_R2", "SKBR3_untreated_R3"), "Control", "lapatinib"))


# Aggregate the replicates for each condition
SKBR3_EP_dexseq_ANNOTATED_norm_TX_aggregated_data_avg_counts <- SKBR3_EP_dexseq_ANNOTATED_long_data %>%
  group_by(ensembl_transcript_id_version, condition, Normalized_MFE, Specific_heat) %>%
  summarise(
    avg_counts = mean(counts),
    .groups = 'drop'
  )


# Define the 75th percentile of the entire distribution of MFE values for the SKBR3_EP_dexseq_ANNOTATED data
SKBR3_EP_dexseq_ANNOTATED_MFE_all_TX_75_percentile <- quantile(SKBR3_EP_dexseq_ANNOTATED_norm_TX_aggregated_data_avg_counts$Normalized_MFE, 0.75)

# Identify all the transcripts with MFE values above the 75th percentile
SKBR3_EP_dexseq_ANNOTATED_high_MFE_TX_total <- dplyr::filter(SKBR3_EP_dexseq_ANNOTATED_norm_TX_counts_MFE_df, Normalized_MFE > SKBR3_EP_dexseq_ANNOTATED_MFE_all_TX_75_percentile)





# Fraction of total high MFE transcripts counts in control and drug-treated groups
SKBR3_EP_dexseq_ANNOTATED_control_total_counts <- dplyr::filter(SKBR3_EP_dexseq_ANNOTATED_norm_TX_aggregated_data_avg_counts, condition == "Control") %>%
  dplyr::summarise(total_counts = sum(avg_counts))

SKBR3_EP_dexseq_ANNOTATED_control_high_MFE_TX_counts <- dplyr::filter(SKBR3_EP_dexseq_ANNOTATED_norm_TX_aggregated_data_avg_counts, condition == "Control") %>%
  dplyr::filter(ensembl_transcript_id_version %in% SKBR3_EP_dexseq_ANNOTATED_high_MFE_TX_total$ensembl_transcript_id_version) %>%
  dplyr::summarise(total_counts_high_MFE_TX = sum(avg_counts))

SKBR3_EP_dexseq_ANNOTATED_control_high_MFE_TX_counts_fraction_EP <- SKBR3_EP_dexseq_ANNOTATED_control_high_MFE_TX_counts$total_counts_high_MFE_TX / SKBR3_EP_dexseq_ANNOTATED_control_total_counts$total_counts
print(SKBR3_EP_dexseq_ANNOTATED_control_high_MFE_TX_counts_fraction_EP) # This value goes in the summary table





SKBR3_EP_dexseq_ANNOTATED_lapatinib_total_counts <- dplyr::filter(SKBR3_EP_dexseq_ANNOTATED_norm_TX_aggregated_data_avg_counts, condition == "lapatinib") %>%
  dplyr::summarise(total_counts = sum(avg_counts))

SKBR3_EP_dexseq_ANNOTATED_lapatinib_high_MFE_TX_counts <- dplyr::filter(SKBR3_EP_dexseq_ANNOTATED_norm_TX_aggregated_data_avg_counts, condition == "lapatinib") %>%
  dplyr::filter(ensembl_transcript_id_version %in% SKBR3_EP_dexseq_ANNOTATED_high_MFE_TX_total$ensembl_transcript_id_version) %>%
  dplyr::summarise(total_counts_high_MFE_TX = sum(avg_counts))

SKBR3_EP_dexseq_ANNOTATED_lapatinib_high_MFE_TX_counts_fraction_EP <- SKBR3_EP_dexseq_ANNOTATED_lapatinib_high_MFE_TX_counts$total_counts_high_MFE_TX / SKBR3_EP_dexseq_ANNOTATED_lapatinib_total_counts$total_counts
print(SKBR3_EP_dexseq_ANNOTATED_lapatinib_high_MFE_TX_counts_fraction_EP) # This value goes in the summary table





# Fisher's exact test to compare the proportions of high MFE transcripts counts between the two groups
# Create the groups
SKBR3_EP_dexseq_ANNOTATED_high_MFE_TX_counts_control_CT <- SKBR3_EP_dexseq_ANNOTATED_control_high_MFE_TX_counts_fraction_EP * SKBR3_EP_dexseq_ANNOTATED_control_total_counts$total_counts
SKBR3_EP_dexseq_ANNOTATED_low_MFE_TX_counts_control_CT <- SKBR3_EP_dexseq_ANNOTATED_control_total_counts$total_counts - SKBR3_EP_dexseq_ANNOTATED_control_high_MFE_TX_counts$total_counts_high_MFE_TX

SKBR3_EP_dexseq_ANNOTATED_high_MFE_TX_counts_lapatinib_CT <- SKBR3_EP_dexseq_ANNOTATED_lapatinib_high_MFE_TX_counts_fraction_EP * SKBR3_EP_dexseq_ANNOTATED_lapatinib_total_counts$total_counts
SKBR3_EP_dexseq_ANNOTATED_low_MFE_TX_counts_lapatinib_CT <- SKBR3_EP_dexseq_ANNOTATED_lapatinib_total_counts$total_counts - SKBR3_EP_dexseq_ANNOTATED_lapatinib_high_MFE_TX_counts$total_counts_high_MFE_TX

# Construct the contingency table with raw counts
SKBR3_EP_dexseq_ANNOTATED_counts_contingency_table <- matrix(
  c(
    SKBR3_EP_dexseq_ANNOTATED_control_high_MFE_TX_counts$total_counts_high_MFE_TX, 
    SKBR3_EP_dexseq_ANNOTATED_lapatinib_high_MFE_TX_counts$total_counts_high_MFE_TX,
    SKBR3_EP_dexseq_ANNOTATED_control_total_counts$total_counts - SKBR3_EP_dexseq_ANNOTATED_control_high_MFE_TX_counts$total_counts_high_MFE_TX, 
    SKBR3_EP_dexseq_ANNOTATED_lapatinib_total_counts$total_counts - SKBR3_EP_dexseq_ANNOTATED_lapatinib_high_MFE_TX_counts$total_counts_high_MFE_TX
  ),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(Category = c("High_MFE_TX_counts", "Low_MFE_TX_counts"),
                  Group = c("Control", "lapatinib"))
)

# Ensure the counts are rounded to integers (they should already be integers if the input data is correct)
# SKBR3_EP_dexseq_ANNOTATED_counts_contingency_table <- round(SKBR3_EP_dexseq_ANNOTATED_counts_contingency_table)

SKBR3_EP_dexseq_ANNOTATED_counts_contingency_table_df <- as.data.frame(SKBR3_EP_dexseq_ANNOTATED_counts_contingency_table)

SKBR3_EP_dexseq_ANNOTATED_fisher_result_counts <- fisher.test(SKBR3_EP_dexseq_ANNOTATED_counts_contingency_table_df)
print(SKBR3_EP_dexseq_ANNOTATED_fisher_result_counts)


cat("SKBR3_EP_dexseq_ANNOTATED Fisher P-value Fisher proportions of high MFE transcripts counts in drug-treated and control :", SKBR3_EP_dexseq_ANNOTATED_fisher_result_counts$p.value, "\n")  # This value goes in the summary table


# Output the results
cat("SKBR3_EP_dexseq_ANNOTATED control high MFE TX as fraction of total counts:", SKBR3_EP_dexseq_ANNOTATED_control_high_MFE_TX_counts_fraction_EP, "\n")
cat("SKBR3_EP_dexseq_ANNOTATED lapatinib high MFE TX as fraction of total counts:", SKBR3_EP_dexseq_ANNOTATED_lapatinib_high_MFE_TX_counts_fraction_EP, "\n")
cat("SKBR3_EP_dexseq_ANNOTATED net change counts high MFE TX drug vs control", SKBR3_EP_dexseq_ANNOTATED_lapatinib_high_MFE_TX_counts_fraction_EP-SKBR3_EP_dexseq_ANNOTATED_control_high_MFE_TX_counts_fraction_EP, "\n")





# Perform a two-proportion z-test instead of Fisher's exact test
# Extract the raw counts for high MFE transcripts in both groups
high_MFE_control <- SKBR3_EP_dexseq_ANNOTATED_control_high_MFE_TX_counts$total_counts_high_MFE_TX
high_MFE_lapatinib <- SKBR3_EP_dexseq_ANNOTATED_lapatinib_high_MFE_TX_counts$total_counts_high_MFE_TX

# Total counts in both groups
total_control <- SKBR3_EP_dexseq_ANNOTATED_control_total_counts$total_counts
total_lapatinib <- SKBR3_EP_dexseq_ANNOTATED_lapatinib_total_counts$total_counts

# Perform the two-proportion z-test
z_test_result <- prop.test(
  c(high_MFE_control, high_MFE_lapatinib), 
  c(total_control, total_lapatinib)
)

# Output the result
SKBR3_EP_dexseq_ANNOTATED_z_test_counts_75_p_value <- z_test_result$p.value










# Plot Normalized_MFE vs transcript counts of both control and drug-treated conditions
# Create separate data frames for Control and lapatinib conditions
SKBR3_EP_dexseq_ANNOTATED_control_data <- SKBR3_EP_dexseq_ANNOTATED_norm_TX_aggregated_data_avg_counts %>%
  dplyr::filter(condition == "Control")

SKBR3_EP_dexseq_ANNOTATED_lapatinib_data <- SKBR3_EP_dexseq_ANNOTATED_norm_TX_aggregated_data_avg_counts %>%
  dplyr::filter(condition == "lapatinib") 


# Generate the plot
Plot_SKBR3_EP_dexseq_ANNOTATED_control <- ggplot(data = SKBR3_EP_dexseq_ANNOTATED_control_data, aes(x = avg_counts, y = Normalized_MFE, color = condition)) +
  geom_point(alpha = 0.6) +
  theme_classic() +
  labs(
    title = "SKBR3_EP_dexseq_ANNOTATED MFE vs. total avg counts Control",
    x = "Transcript Counts",
    y = "Normalized MFE (kcal/mol)"
  ) +
  #xlim(0, 150000) +
  ylim(-0.8, 0) +
  theme(
    axis.title = element_text(size = 15),
    plot.title = element_text(size = 10)
  ) +
  guides(color = FALSE)


Plot_SKBR3_EP_dexseq_ANNOTATED_lapatinib <- ggplot(data = SKBR3_EP_dexseq_ANNOTATED_lapatinib_data, aes(x = avg_counts, y = Normalized_MFE, color = condition)) +
  geom_point(alpha = 0.6) +
  theme_classic() +
  labs(
    title = "SKBR3_EP_dexseq_ANNOTATED MFE vs. total avg counts lapatinib",
    x = "Transcript Counts",
    y = "Normalized MFE (kcal/mol)"
  ) +
  #xlim(0, 150000) +
  ylim(-0.8, 0) +
  theme(
    axis.title = element_text(size = 15),
    plot.title = element_text(size = 10)
  ) +
  guides(color = FALSE)

library(patchwork)
SKBR3_EP_dexseq_ANNOTATED_combined_plots <- Plot_SKBR3_EP_dexseq_ANNOTATED_control + Plot_SKBR3_EP_dexseq_ANNOTATED_lapatinib
SKBR3_EP_dexseq_ANNOTATED_combined_plots







# Define the 25th percentile of the entire distribution of MFE values for the SKBR3_EP_dexseq_ANNOTATED data
SKBR3_EP_dexseq_ANNOTATED_MFE_all_TX_25_percentile <- quantile(SKBR3_EP_dexseq_ANNOTATED_norm_TX_aggregated_data_avg_counts$Normalized_MFE, 0.25)

# Identify all the transcripts with MFE values below the 25th percentile
SKBR3_EP_dexseq_ANNOTATED_low_MFE_TX_total <- dplyr::filter(SKBR3_EP_dexseq_ANNOTATED_norm_TX_counts_MFE_df, Normalized_MFE < SKBR3_EP_dexseq_ANNOTATED_MFE_all_TX_25_percentile)





# Fraction of total low MFE transcripts counts in control and drug-treated groups
SKBR3_EP_dexseq_ANNOTATED_control_total_counts <- dplyr::filter(SKBR3_EP_dexseq_ANNOTATED_norm_TX_aggregated_data_avg_counts, condition == "Control") %>%
  dplyr::summarise(total_counts = sum(avg_counts))

SKBR3_EP_dexseq_ANNOTATED_control_low_MFE_TX_counts <- dplyr::filter(SKBR3_EP_dexseq_ANNOTATED_norm_TX_aggregated_data_avg_counts, condition == "Control") %>%
  dplyr::filter(ensembl_transcript_id_version %in% SKBR3_EP_dexseq_ANNOTATED_low_MFE_TX_total$ensembl_transcript_id_version) %>%
  dplyr::summarise(total_counts_low_MFE_TX = sum(avg_counts))

SKBR3_EP_dexseq_ANNOTATED_control_low_MFE_TX_counts_fraction_EP <- SKBR3_EP_dexseq_ANNOTATED_control_low_MFE_TX_counts$total_counts_low_MFE_TX / SKBR3_EP_dexseq_ANNOTATED_control_total_counts$total_counts
print(SKBR3_EP_dexseq_ANNOTATED_control_low_MFE_TX_counts_fraction_EP) # This value goes in the summary table





SKBR3_EP_dexseq_ANNOTATED_lapatinib_total_counts <- dplyr::filter(SKBR3_EP_dexseq_ANNOTATED_norm_TX_aggregated_data_avg_counts, condition == "lapatinib") %>%
  dplyr::summarise(total_counts = sum(avg_counts))

SKBR3_EP_dexseq_ANNOTATED_lapatinib_low_MFE_TX_counts <- dplyr::filter(SKBR3_EP_dexseq_ANNOTATED_norm_TX_aggregated_data_avg_counts, condition == "lapatinib") %>%
  dplyr::filter(ensembl_transcript_id_version %in% SKBR3_EP_dexseq_ANNOTATED_low_MFE_TX_total$ensembl_transcript_id_version) %>%
  dplyr::summarise(total_counts_low_MFE_TX = sum(avg_counts))

SKBR3_EP_dexseq_ANNOTATED_lapatinib_low_MFE_TX_counts_fraction_EP <- SKBR3_EP_dexseq_ANNOTATED_lapatinib_low_MFE_TX_counts$total_counts_low_MFE_TX / SKBR3_EP_dexseq_ANNOTATED_lapatinib_total_counts$total_counts
print(SKBR3_EP_dexseq_ANNOTATED_lapatinib_low_MFE_TX_counts_fraction_EP) # This value goes in the summary table





# Fisher's exact test to compare the proportions of low MFE transcripts counts between the two groups
# Create the groups
SKBR3_EP_dexseq_ANNOTATED_low_MFE_TX_counts_control_CT <- SKBR3_EP_dexseq_ANNOTATED_control_low_MFE_TX_counts_fraction_EP * SKBR3_EP_dexseq_ANNOTATED_control_total_counts$total_counts
SKBR3_EP_dexseq_ANNOTATED_low_MFE_TX_counts_control_CT <- SKBR3_EP_dexseq_ANNOTATED_control_total_counts$total_counts - SKBR3_EP_dexseq_ANNOTATED_control_low_MFE_TX_counts$total_counts_low_MFE_TX

SKBR3_EP_dexseq_ANNOTATED_low_MFE_TX_counts_lapatinib_CT <- SKBR3_EP_dexseq_ANNOTATED_lapatinib_low_MFE_TX_counts_fraction_EP * SKBR3_EP_dexseq_ANNOTATED_lapatinib_total_counts$total_counts
SKBR3_EP_dexseq_ANNOTATED_low_MFE_TX_counts_lapatinib_CT <- SKBR3_EP_dexseq_ANNOTATED_lapatinib_total_counts$total_counts - SKBR3_EP_dexseq_ANNOTATED_lapatinib_low_MFE_TX_counts$total_counts_low_MFE_TX

# Construct the contingency table with raw counts
SKBR3_EP_dexseq_ANNOTATED_counts_contingency_table_low <- matrix(
  c(
    SKBR3_EP_dexseq_ANNOTATED_control_low_MFE_TX_counts$total_counts_low_MFE_TX, 
    SKBR3_EP_dexseq_ANNOTATED_lapatinib_low_MFE_TX_counts$total_counts_low_MFE_TX,
    SKBR3_EP_dexseq_ANNOTATED_control_total_counts$total_counts - SKBR3_EP_dexseq_ANNOTATED_control_low_MFE_TX_counts$total_counts_low_MFE_TX, 
    SKBR3_EP_dexseq_ANNOTATED_lapatinib_total_counts$total_counts - SKBR3_EP_dexseq_ANNOTATED_lapatinib_low_MFE_TX_counts$total_counts_low_MFE_TX
  ),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(Category = c("Low_MFE_TX_counts", "Other_TX_counts"),
                  Group = c("Control", "lapatinib"))
)

SKBR3_EP_dexseq_ANNOTATED_counts_contingency_table_df_low <- as.data.frame(SKBR3_EP_dexseq_ANNOTATED_counts_contingency_table_low)

SKBR3_EP_dexseq_ANNOTATED_fisher_result_counts_low <- fisher.test(SKBR3_EP_dexseq_ANNOTATED_counts_contingency_table_low)
print(SKBR3_EP_dexseq_ANNOTATED_fisher_result_counts_low)

cat("SKBR3_EP_dexseq_ANNOTATED Fisher P-value Fisher proportions of low MFE transcripts counts in drug-treated and control :", SKBR3_EP_dexseq_ANNOTATED_fisher_result_counts$p.value, "\n")  # This value goes in the summary table


# Output the results
cat("SKBR3_EP_dexseq_ANNOTATED control low MFE TX as fraction of total counts:", SKBR3_EP_dexseq_ANNOTATED_control_low_MFE_TX_counts_fraction_EP, "\n")
cat("SKBR3_EP_dexseq_ANNOTATED palbociclib low MFE TX as fraction of total counts:", SKBR3_EP_dexseq_ANNOTATED_lapatinib_low_MFE_TX_counts_fraction_EP, "\n")
cat("SKBR3_EP_dexseq_ANNOTATED net change counts low MFE TX drug vs control", SKBR3_EP_dexseq_ANNOTATED_lapatinib_low_MFE_TX_counts_fraction_EP-SKBR3_EP_dexseq_ANNOTATED_control_low_MFE_TX_counts_fraction_EP, "\n")




# Perform a two-proportion z-test instead of Fisher's exact test
# Extract the raw counts for low MFE transcripts in both groups
low_MFE_control <- SKBR3_EP_dexseq_ANNOTATED_control_low_MFE_TX_counts$total_counts_low_MFE_TX
low_MFE_lapatinib <- SKBR3_EP_dexseq_ANNOTATED_lapatinib_low_MFE_TX_counts$total_counts_low_MFE_TX

# Total counts in both groups
total_control <- SKBR3_EP_dexseq_ANNOTATED_control_total_counts$total_counts
total_lapatinib <- SKBR3_EP_dexseq_ANNOTATED_lapatinib_total_counts$total_counts

# Perform the two-proportion z-test
z_test_result_low <- prop.test(
  c(low_MFE_control, low_MFE_lapatinib), 
  c(total_control, total_lapatinib)
)

# Output the result
SKBR3_EP_dexseq_ANNOTATED_z_test_counts_25_p_value <- z_test_result_low$p.value









## SKBR3_counts_EP heat capacity changes ####

# This is calculated on the total n. of transcripts with a significant log2 fold change
# Generate a scatter plot Log2_FC vs heat capacity
ggplot(SKBR3_counts_EP, aes(y=Specific_heat, x=log2fold_treated_untreated)) +
  geom_point(size=1.5, alpha=0.1)+
  theme_classic() +
  theme(legend.position="none") +
  geom_vline(xintercept = 0, size=0.3) +
  xlim(-12, 12) +
  ylim(0,400) +
  theme(axis.title = element_text(size = 20)) +
  theme(plot.title = element_text(size = 15)) +
  labs(y="Heat capacity (J/C)",x="Log2_FC_expression",title="SKBR3_EP_dexseq_ANNOTATED_Log2_FC_adj.P<0.05")


# Show the data with a box plot
ggplot(SKBR3_counts_EP, aes(x = factor(log2fold_treated_untreated > 0), y = Specific_heat)) +
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
  labs(y="Heat capacity (J/C)",x="",title="SKBR3_EP_dexseq_ANNOTATED_Log2_FC_adj.P<0.05 heat_capacity all transcripts") +
  scale_x_discrete(labels = c("log2FoldChange < 0", "log2FoldChange > 0"))

Mean_heat_capacity_SKBR3_EP_dexseq_ANNOTATED_UP_EP <- mean(SKBR3_EP_dexseq_ANNOTATED_UP_EP$Specific_heat)
Mean_heat_capacity_SKBR3_EP_dexseq_ANNOTATED_DOWN_EP <- mean(SKBR3_EP_dexseq_ANNOTATED_DOWN_EP$Specific_heat)
TTEST_heat_capacity_SKBR3_EP_dexseq_ANNOTATED_all_UP_vs_DOWN <- t.test(SKBR3_EP_dexseq_ANNOTATED_DOWN_EP$Specific_heat, SKBR3_EP_dexseq_ANNOTATED_UP_EP$Specific_heat, paired=FALSE)
print(TTEST_heat_capacity_SKBR3_EP_dexseq_ANNOTATED_all_UP_vs_DOWN)
# Add K-S test to look at the differences between the up-regukated and down-regulated genes
ks_test_SKBR3_EP_dexseq_ANNOTATED_all_UP_vs_DOWN_heat_capacity <- ks.test(SKBR3_EP_dexseq_ANNOTATED_DOWN_EP$Specific_heat, SKBR3_EP_dexseq_ANNOTATED_UP_EP$Specific_heat)
print(ks_test_SKBR3_EP_dexseq_ANNOTATED_all_UP_vs_DOWN_heat_capacity)




# Calculate the changes of heat capacity only in the high MFE transcripts

# Contribution of the UP or DOWN-regulated high MFE transcripts to the heat capacity
# Look at the heat capacity of the UP-regulated high MFE transcripts
SKBR3_EP_dexseq_ANNOTATED_high_MFE_TX_in_UP_heat_capacity <- mean(SKBR3_EP_dexseq_ANNOTATED_high_MFE_TX_in_UP$Specific_heat)

# Look at the heat capacity of the DOWN-regulated high MFE transcripts
SKBR3_EP_dexseq_ANNOTATED_high_MFE_TX_in_DOWN_heat_capacity <- mean(SKBR3_EP_dexseq_ANNOTATED_high_MFE_TX_in_DOWN$Specific_heat)

# Heat capacity gain in UP-regulated high MFE transcripts vs DOWN-regulated high MFE transcripts
SKBR3_EP_dexseq_ANNOTATED_heat_capacity_gain_UP_vs_DOWN <- SKBR3_EP_dexseq_ANNOTATED_high_MFE_TX_in_UP_heat_capacity - SKBR3_EP_dexseq_ANNOTATED_high_MFE_TX_in_DOWN_heat_capacity
cat("SKBR3_EP_dexseq_ANNOTATED net heat capacity UP vs DOWN high MFE TX:", SKBR3_EP_dexseq_ANNOTATED_heat_capacity_gain_UP_vs_DOWN, "\n")


# Extract specific heat data
SKBR3_EP_dexseq_ANNOTATED_Specific_heat_UP <- SKBR3_EP_dexseq_ANNOTATED_high_MFE_TX_in_UP$Specific_heat
SKBR3_EP_dexseq_ANNOTATED_Specific_heat_DOWN <- SKBR3_EP_dexseq_ANNOTATED_high_MFE_TX_in_DOWN$Specific_heat

# Combine data into a single data frame for easier plotting
SKBR3_EP_dexseq_ANNOTATED_Specific_heat_data <- data.frame(
  SKBR3_EP_dexseq_ANNOTATED_Specific_heat = c(SKBR3_EP_dexseq_ANNOTATED_Specific_heat_UP, SKBR3_EP_dexseq_ANNOTATED_Specific_heat_DOWN),
  regulation = factor(c(rep("UP", length(SKBR3_EP_dexseq_ANNOTATED_Specific_heat_UP)), rep("DOWN", length(SKBR3_EP_dexseq_ANNOTATED_Specific_heat_DOWN))))
)

# Generate box plot
ggplot(SKBR3_EP_dexseq_ANNOTATED_Specific_heat_data, aes(x = regulation, y = SKBR3_EP_dexseq_ANNOTATED_Specific_heat, fill = regulation)) +
  geom_boxplot(fill = "lightgrey", alpha = 0.7) +
  labs(title = "SKBR3_EP_dexseq_ANNOTATED specific Heat of High MFE Transcripts",
       x = "",
       y = "Specific Heat (J/C)") +
  scale_x_discrete(labels = c("log2_FC < 0", "log2_FC > 0")) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 16),
        axis.title = element_text(size = 25))

Mean_heat_capacity_SKBR3_EP_dexseq_ANNOTATED_high_MFE_TX_UP <- mean(SKBR3_EP_dexseq_ANNOTATED_high_MFE_TX_in_UP$Specific_heat)
Mean_heat_capacity_SKBR3_EP_dexseq_ANNOTATED_high_MFE_TX_DOWN <- mean(SKBR3_EP_dexseq_ANNOTATED_high_MFE_TX_in_DOWN$Specific_heat)
TTEST_heat_capacity_SKBR3_EP_dexseq_ANNOTATED_all_high_MFE_TX_UP_vs_DOWN <- t.test(SKBR3_EP_dexseq_ANNOTATED_high_MFE_TX_in_DOWN$Specific_heat, SKBR3_EP_dexseq_ANNOTATED_high_MFE_TX_in_UP$Specific_heat, paired=FALSE)
print(TTEST_heat_capacity_SKBR3_EP_dexseq_ANNOTATED_all_high_MFE_TX_UP_vs_DOWN)

# Add K-S test to look at the differences between the up-regukated and down-regulated genes
ks_test_heat_capacity_SKBR3_EP_dexseq_ANNOTATED_all_high_MFE_TX_UP_vs_DOWN <- ks.test(SKBR3_EP_dexseq_ANNOTATED_high_MFE_TX_in_DOWN$Specific_heat, SKBR3_EP_dexseq_ANNOTATED_high_MFE_TX_in_UP$Specific_heat)
print(ks_test_heat_capacity_SKBR3_EP_dexseq_ANNOTATED_all_high_MFE_TX_UP_vs_DOWN)



















# Analysis cell proliferation Nanopore data SKBR3 LAP ####

SKBR3_LAP_DESeq2 <- read.delim(file="path_to/early.phase_SKBR3_lapatinib.tx_counts.both_conditions.treated_untreated_deseq2.tab", header = TRUE, sep = "\t", quote = "")


# Filter the proliferation gene sets by using a function
# Use a function also to determine the mean and standard deviation of each set and the plot all of them together

# Generate the function that filters each gene set in the target data frame
# This function will be used in the following commands
filtering <- function(x) {
  Filtering <- dplyr::filter(SKBR3_LAP_DESeq2, gene_name %in% x)
  return(Filtering)
}

# Generate the list containing all the gene sets to be used for the filtering.
CURATED_PATHWAYS.input <- list(MTOR_positive_regulators, MTOR_negative_regulators, Proliferation_Science_paper, Meta_PCNA, Proliferation_GO_genes_Magill_Sack_2018, Proliferation_STOP_genes_Magill_Sack_2018, Reactome_G1_phase, Proliferation_Carey_2021, Essential_cell_cycle_Viner_Breuer_Cell_Division_2019, Essential_S_Phase_Viner_Breuer_Cell_Division_2019, Homologous_recombination, Fanconi_anemia, High_fidelity_polymerases)

# This creates the names that will be used for the final data frames of each file resulting from the filtering of the target data frame with each gene set
CURATED_PATHWAYS.output <- c("SKBR3_LAP_MTOR_POS_EP", "SKBR3_LAP_MTOR_NEG_EP", "SKBR3_LAP_Proliferation_Science_EP", "SKBR3_LAP_Meta_PCNA_EP", "SKBR3_LAP_Proliferation_GO_EP", "SKBR3_LAP_Proliferation_STOP_EP", "SKBR3_LAP_G1_phase_EP", "SKBR3_LAP_Proliferation_Carey_EP", "SKBR3_LAP_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019_EP", "SKBR3_LAP_Essential_S_Phase_Viner_Breuer_Cell_Division_2019_EP", "SKBR3_LAP_HR_EP", "SKBR3_LAP_Fanconi_anemia_EP", "SKBR3_LAP_HF_POL_EP")

# Assign the names of the gene sets to each entry in the list
names(CURATED_PATHWAYS.input) <- CURATED_PATHWAYS.output

# Apply the filtering function to create a list of data frames
Curated_test_df <- lapply(CURATED_PATHWAYS.input, function(x) filtering(x))

# creating a vector with desired names of filtered data frames
Curated_dfs.output <- c("df_SKBR3_LAP_MTOR_POS_EP", "df_SKBR3_LAP_MTOR_NEG_EP", "df_SKBR3_LAP_Proliferation_Science_EP", "df_SKBR3_LAP_Meta_PCNA_EP", "df_SKBR3_LAP_Proliferation_GO_EP", "df_SKBR3_LAP_Proliferation_STOP_EP", "df_SKBR3_LAP_G1_phase_EP", "df_SKBR3_LAP_Proliferation_Carey_EP", "df_SKBR3_LAP_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019_EP", "df_SKBR3_LAP_Essential_S_Phase_Viner_Breuer_Cell_Division_2019_EP", "df_SKBR3_LAP_HR_EP", "df_SKBR3_LAP_Fanconi_anemia_EP", "df_SKBR3_LAP_HF_POL_EP")

# loop to assign data frames to the desired variable name
for (i in seq_along(Curated_dfs.output)) {
  assign(Curated_dfs.output[i], Curated_test_df[[i]])
}





# Generate the function to calculate the MEAN of the log2FoldChange of each gene set resulting from the filtering above.
Mean_curated <- function(x) {
  mean_curated <- mean(x$log2FoldChange, na.rm = TRUE)
  return(mean_curated)
}

# Create an input of the data frames generated above to calculate the mean of log2FoldChange
Curated_input_dfs <- list(df_SKBR3_LAP_MTOR_POS_EP, df_SKBR3_LAP_MTOR_NEG_EP, df_SKBR3_LAP_Proliferation_Science_EP, df_SKBR3_LAP_Meta_PCNA_EP, df_SKBR3_LAP_Proliferation_GO_EP, df_SKBR3_LAP_Proliferation_STOP_EP, df_SKBR3_LAP_G1_phase_EP, df_SKBR3_LAP_Proliferation_Carey_EP, df_SKBR3_LAP_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019_EP, df_SKBR3_LAP_Essential_S_Phase_Viner_Breuer_Cell_Division_2019_EP, df_SKBR3_LAP_HR_EP, df_SKBR3_LAP_Fanconi_anemia_EP, df_SKBR3_LAP_HF_POL_EP)

# Apply the mean function to create a list of means
Curated_mean_test_df <- lapply(Curated_input_dfs, function(x) Mean_curated(x))

# creating a vector with desired names of the files containing the mean of each gene set log2FoldChange
Curated_Mean_SKBR3_LAP_dfs.output <- c("Mean_SKBR3_LAP_MTOR_POS_EP", "Mean_SKBR3_LAP_MTOR_NEG_EP", "Mean_SKBR3_LAP_Proliferation_Science_EP", "Mean_SKBR3_LAP_Meta_PCNA_EP", "Mean_SKBR3_LAP_Proliferation_GO_EP", "Mean_SKBR3_LAP_Proliferation_STOP_EP", "Mean_SKBR3_LAP_G1_phase_EP", "Mean_SKBR3_LAP_Proliferation_Carey_EP", "Mean_SKBR3_LAP_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019_EP", "Mean_SKBR3_LAP_Essential_S_Phase_Viner_Breuer_Cell_Division_2019_EP", "Mean_SKBR3_LAP_HR_EP", "Mean_SKBR3_LAP_Fanconi_anemia_EP", "Mean_SKBR3_LAP_HF_POL_EP")

# loop to assign means to the desired variable name
for (i in seq_along(Curated_Mean_SKBR3_LAP_dfs.output)) {
  assign(Curated_Mean_SKBR3_LAP_dfs.output[i], Curated_mean_test_df[[i]])
}




# Generate the function to calculate the STANDARD DEVIATION of the log2FoldChange of each gene set resulting from the filtering above.
STDEV_curated <- function(x) {
  stdev_curated <- sd(x$log2FoldChange, na.rm = TRUE)
  return(stdev_curated)
}

# Apply the standard deviation function to create a list of standard deviations
Curated_stdev_test_df <- lapply(Curated_input_dfs, function(x) STDEV_curated(x))

# creating a vector with desired names of the files containing the standard deviation of each gene set log2FoldChange
Curated_STDEV_SKBR3_LAP_dfs.output <- c("STDEV_SKBR3_LAP_MTOR_POS_EP", "STDEV_SKBR3_LAP_MTOR_NEG_EP", "STDEV_SKBR3_LAP_Proliferation_Science_EP", "STDEV_SKBR3_LAP_Meta_PCNA_EP", "STDEV_SKBR3_LAP_Proliferation_GO_EP", "STDEV_SKBR3_LAP_Proliferation_STOP_EP", "STDEV_SKBR3_LAP_G1_phase_EP", "STDEV_SKBR3_LAP_Proliferation_Carey_EP", "STDEV_SKBR3_LAP_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019_EP", "STDEV_SKBR3_LAP_Essential_S_Phase_Viner_Breuer_Cell_Division_2019_EP", "STDEV_SKBR3_LAP_HR_EP", "STDEV_SKBR3_LAP_Fanconi_anemia_EP", "STDEV_SKBR3_LAP_HF_POL_EP")
# loop to assign data.frames to the desired variable name
for (i in (seq_along(Curated_STDEV_SKBR3_LAP_dfs.output))) {
  do.call("<-",list(Curated_STDEV_SKBR3_LAP_dfs.output[i], Curated_stdev_test_df[[i]]))
}




# Generate the data frame for the plot

# Generate the data frame for the plot
df_778_PALBO_proliferation <- data.frame(mean=c(Mean_SKBR3_LAP_MTOR_POS_EP, Mean_SKBR3_LAP_MTOR_NEG_EP, Mean_SKBR3_LAP_Proliferation_Science_EP, Mean_SKBR3_LAP_Meta_PCNA_EP, Mean_SKBR3_LAP_Proliferation_GO_EP, Mean_SKBR3_LAP_Proliferation_STOP_EP, Mean_SKBR3_LAP_G1_phase_EP, Mean_SKBR3_LAP_Proliferation_Carey_EP, Mean_SKBR3_LAP_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019_EP, Mean_SKBR3_LAP_Essential_S_Phase_Viner_Breuer_Cell_Division_2019_EP, Mean_SKBR3_LAP_HR_EP, Mean_SKBR3_LAP_Fanconi_anemia_EP, Mean_SKBR3_LAP_HF_POL_EP),
                                         sd=c(STDEV_SKBR3_LAP_MTOR_POS_EP, STDEV_SKBR3_LAP_MTOR_NEG_EP, STDEV_SKBR3_LAP_Proliferation_Science_EP, STDEV_SKBR3_LAP_Meta_PCNA_EP, STDEV_SKBR3_LAP_Proliferation_GO_EP, STDEV_SKBR3_LAP_Proliferation_STOP_EP, STDEV_SKBR3_LAP_G1_phase_EP, STDEV_SKBR3_LAP_Proliferation_Carey_EP, STDEV_SKBR3_LAP_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019_EP, STDEV_SKBR3_LAP_Essential_S_Phase_Viner_Breuer_Cell_Division_2019_EP, STDEV_SKBR3_LAP_HR_EP, STDEV_SKBR3_LAP_Fanconi_anemia_EP, STDEV_SKBR3_LAP_HF_POL_EP),
                                         Category=as.factor(c("SKBR3_LAP_MTOR_POS_EP", "SKBR3_LAP_MTOR_NEG_EP", "SKBR3_LAP_Proliferation_Science_EP", "SKBR3_LAP_Meta_PCNA_EP", "SKBR3_LAP_Proliferation_GO_EP", "SKBR3_LAP_Proliferation_STOP_EP", "SKBR3_LAP_G1_phase_EP", "SKBR3_LAP_Proliferation_Carey_EP", "SKBR3_LAP_Essential_cell_cycle_EP", "SKBR3_LAP_Essential_S_Phase_EP", "SKBR3_LAP_HR_EP", "SKBR3_LAP_FA_EP", "SKBR3_LAP_HF_POL_EP")))
df_778_PALBO_proliferation




library(ggplot2)
ggplot(df_778_PALBO_proliferation, aes(x=reorder(Category, -mean), y=mean)) +
  theme_bw() +
  geom_point(position=position_dodge(width = 1), stat="identity",colour='black') +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.4) +
  theme(axis.text.x = element_text(size=12, angle=45,hjust=0.95,vjust=0.98)) +
  theme(axis.text.y = element_text(size=14)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 10))) +
  geom_point(colour = "black", size = 3) +
  xlab("") + ylab("Relative log2 fold change") +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  ggtitle("Analysis cell proliferation 778 PALBO NANOPORE 72h vs T0") +
  scale_x_discrete(limits = c("SKBR3_LAP_Proliferation_Science_EP", "SKBR3_LAP_Meta_PCNA_EP", "SKBR3_LAP_Proliferation_GO_EP", "SKBR3_LAP_Proliferation_STOP_EP", "SKBR3_LAP_Proliferation_Carey_EP", "SKBR3_LAP_G1_phase_EP", "SKBR3_LAP_MTOR_POS_EP", "SKBR3_LAP_MTOR_NEG_EP", "SKBR3_LAP_Essential_cell_cycle_EP", "SKBR3_LAP_Essential_S_Phase_EP", "SKBR3_LAP_HR_EP", "SKBR3_LAP_FA_EP", "SKBR3_LAP_HF_POL_EP")) +
  theme(plot.title = element_text(color="black", size=14, face="bold"))


# Merge all the data frames to generate box plots
# Add a column with the specific name of the gene set
df_SKBR3_LAP_MTOR_POS_EP$Gene_set <- "SKBR3_LAP_MTOR_POS_EP"
df_SKBR3_LAP_MTOR_NEG_EP$Gene_set <- "SKBR3_LAP_MTOR_NEG_EP"
df_SKBR3_LAP_Proliferation_Science_EP$Gene_set <- "SKBR3_LAP_Proliferation_Science_EP"
df_SKBR3_LAP_Meta_PCNA_EP$Gene_set <- "SKBR3_LAP_Meta_PCNA_EP"
df_SKBR3_LAP_Proliferation_GO_EP$Gene_set <- "SKBR3_LAP_Proliferation_GO_EP"
df_SKBR3_LAP_Proliferation_STOP_EP$Gene_set <- "SKBR3_LAP_Proliferation_STOP_EP"
df_SKBR3_LAP_G1_phase_EP$Gene_set <- "SKBR3_LAP_G1_phase_EP"
df_SKBR3_LAP_Proliferation_Carey_EP$Gene_set <- "SKBR3_LAP_Proliferation_Carey_EP"
df_SKBR3_LAP_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019_EP$Gene_set <- "SKBR3_LAP_Essential_cell_cycle_EP"
df_SKBR3_LAP_Essential_S_Phase_Viner_Breuer_Cell_Division_2019_EP$Gene_set <- "SKBR3_LAP_Essential_S_Phase_EP"
df_SKBR3_LAP_HR_EP$Gene_set <- "SKBR3_LAP_HR_EP"
df_SKBR3_LAP_Fanconi_anemia_EP$Gene_set <- "SKBR3_LAP_Fanconi_anemia_EP"
df_SKBR3_LAP_HF_POL_EP$Gene_set <- "SKBR3_LAP_HF_POL_EP"



# Combine the data frames
combined_df <- rbind(
  df_SKBR3_LAP_MTOR_POS_EP,
  df_SKBR3_LAP_MTOR_NEG_EP,
  df_SKBR3_LAP_Proliferation_Science_EP,
  df_SKBR3_LAP_Meta_PCNA_EP,
  df_SKBR3_LAP_Proliferation_GO_EP,
  df_SKBR3_LAP_Proliferation_STOP_EP,
  df_SKBR3_LAP_G1_phase_EP,
  df_SKBR3_LAP_Proliferation_Carey_EP,
  df_SKBR3_LAP_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019_EP,
  df_SKBR3_LAP_Essential_S_Phase_Viner_Breuer_Cell_Division_2019_EP,
  df_SKBR3_LAP_HR_EP,
  df_SKBR3_LAP_Fanconi_anemia_EP,
  df_SKBR3_LAP_HF_POL_EP
)

# Ensure the Gene_set column is a factor with the correct order
combined_df$Gene_set <- factor(combined_df$Gene_set, levels = c("SKBR3_LAP_Proliferation_Science_EP", "SKBR3_LAP_Meta_PCNA_EP", "SKBR3_LAP_Proliferation_GO_EP", "SKBR3_LAP_Proliferation_STOP_EP", "SKBR3_LAP_Proliferation_Carey_EP", "SKBR3_LAP_G1_phase_EP", "SKBR3_LAP_MTOR_POS_EP", "SKBR3_LAP_MTOR_NEG_EP", "SKBR3_LAP_Essential_cell_cycle_EP", "SKBR3_LAP_Essential_S_Phase_EP", "SKBR3_LAP_HR_EP", "SKBR3_LAP_Fanconi_anemia_EP", "SKBR3_LAP_HF_POL_EP"))


# Generate the box plot
library(ggplot2)
ggplot(combined_df, aes(x=Gene_set, y=log2FoldChange)) +
  geom_boxplot(fill="lightgrey") + 
  stat_summary(fun.y=mean, geom="point", shape=20, size=4, color="blue") +
  # geom_jitter(width=0.2, alpha=0.5) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle=45, hjust=0.95, vjust=0.98)) +
  theme(axis.text.y = element_text(size=14)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 10))) +
  xlab("") + ylab("Relative log2 fold change") +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  ggtitle("Analysis cell proliferation 778 palbociclib early phase") +
  scale_x_discrete(limits = levels(combined_df$Gene_set)) +
  theme(plot.title = element_text(color="black", size=14, face="bold")) 






Mean_df_MTOR_POS_SKBR3_LAP_EP <- mean(df_SKBR3_LAP_MTOR_POS_EP$log2FoldChange)
Mean_df_MTOR_NEG_SKBR3_LAP_EP <- mean(df_SKBR3_LAP_MTOR_NEG_EP$log2FoldChange)
Mean_df_Proliferation_Science_SKBR3_LAP_EP <- mean(df_SKBR3_LAP_Proliferation_Science_EP$log2FoldChange)
Mean_df_Meta_PCNA_SKBR3_LAP_EP <- mean(df_SKBR3_LAP_Meta_PCNA_EP$log2FoldChange)
Mean_df_Proliferation_GO_SKBR3_LAP_EP <- mean(df_SKBR3_LAP_Proliferation_GO_EP$log2FoldChange)
Mean_df_Proliferation_STOP_SKBR3_LAP_EP <- mean(df_SKBR3_LAP_Proliferation_STOP_EP$log2FoldChange)
Mean_df_G1_phase_SKBR3_LAP_EP <- mean(df_SKBR3_LAP_G1_phase_EP$log2FoldChange)
Mean_df_Proliferation_Carey_SKBR3_LAP_EP <- mean(df_SKBR3_LAP_Proliferation_Carey_EP$log2FoldChange)
Mean_df_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019_SKBR3_LAP_EP <- mean(df_SKBR3_LAP_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019_EP$log2FoldChange)
Mean_df_Essential_S_Phase_Viner_Breuer_Cell_Division_2019_SKBR3_LAP_EP <- mean(df_SKBR3_LAP_Essential_S_Phase_Viner_Breuer_Cell_Division_2019_EP$log2FoldChange)
Mean_df_HR_SKBR3_LAP_EP <- mean(df_SKBR3_LAP_HR_EP$log2FoldChange)
Mean_df_FA_SKBR3_LAP_EP <- mean(df_SKBR3_LAP_Fanconi_anemia_EP$log2FoldChange)
Mean_df_HF_POL_SKBR3_LAP_EP <- mean(df_SKBR3_LAP_HF_POL_EP$log2FoldChange)

# Perform the Wilcoxon signed-rank test
wilcox_test_MTOR_POS_SKBR3_LAP_EP <- wilcox.test(df_SKBR3_LAP_MTOR_POS_EP$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_MTOR_NEG_SKBR3_LAP_EP <- wilcox.test(df_SKBR3_LAP_MTOR_NEG_EP$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_Proliferation_Science_SKBR3_LAP_EP <- wilcox.test(df_SKBR3_LAP_Proliferation_Science_EP$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_Meta_PCNA_SKBR3_LAP_EP <- wilcox.test(df_SKBR3_LAP_Meta_PCNA_EP$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_Proliferation_GO_SKBR3_LAP_EP <- wilcox.test(df_SKBR3_LAP_Proliferation_GO_EP$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_Proliferation_STOP_SKBR3_LAP_EP <- wilcox.test(df_SKBR3_LAP_Proliferation_STOP_EP$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_G1_phase_SKBR3_LAP_EP <- wilcox.test(df_SKBR3_LAP_G1_phase_EP$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_Proliferation_Carey_SKBR3_LAP_EP <- wilcox.test(df_SKBR3_LAP_Proliferation_Carey_EP$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019_SKBR3_LAP_EP <- wilcox.test(df_SKBR3_LAP_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019_EP$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_Essential_S_Phase_Viner_Breuer_Cell_Division_2019_SKBR3_LAP_EP <- wilcox.test(df_SKBR3_LAP_Essential_S_Phase_Viner_Breuer_Cell_Division_2019_EP$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_HR_SKBR3_LAP_EP <- wilcox.test(df_SKBR3_LAP_HR_EP$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_FA_SKBR3_LAP_EP <- wilcox.test(df_SKBR3_LAP_Fanconi_anemia_EP$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_HF_POL_SKBR3_LAP_EP <- wilcox.test(df_SKBR3_LAP_HF_POL_EP$log2FoldChange, mu = 0, alternative = "less")



































# SKMEL28 cell line ANNOTATED transcripts ####
SKMEL28_counts_EP <- dplyr::inner_join(SKMEL28_counts_EP_dexseq, Table_ANNOTATED, by="ensembl_transcript_id_version") %>%
  dplyr::filter(padj < 0.05) # The file 'Table' contains the Normalized_MFE values of all ANNOTATED transcripts


# Generate a scatter plot Log2_FC vs Normalized_MFE
library(ggplot2)
ggplot(SKMEL28_counts_EP, aes(y=Normalized_MFE, x=log2fold_treated_untreated)) +
  geom_point(size=1.5, alpha=0.3)+
  theme_classic() +
  theme(legend.position="none") +
  geom_vline(xintercept = 0, size=0.3) +
  #xlim(-5.5,5.5) +
  theme(axis.title = element_text(size = 20)) +
  theme(plot.title = element_text(size = 15)) +
  labs(y="Normalized_MFE",x="Log2_FC_expression",title="SKMEL28_shared_ANNOTATED_TX_early_phase_Log2_FC_adj.P<0.05")


# It is well known that ANNOTATED transcripts, in general have lower MFEs than the UNANNOTATED transcripts

SKMEL28_EP_dexseq_ANNOTATED_UP_EP <- dplyr::filter(SKMEL28_counts_EP, log2fold_treated_untreated > 0)
Mean_SKMEL28_EP_dexseq_ANNOTATED_UP_EP <- mean(SKMEL28_EP_dexseq_ANNOTATED_UP_EP$Normalized_MFE)

SKMEL28_EP_dexseq_ANNOTATED_DOWN_EP <- dplyr::filter(SKMEL28_counts_EP, log2fold_treated_untreated < 0)
Mean_SKMEL28_EP_dexseq_ANNOTATED_DOWN_EP <- mean(SKMEL28_EP_dexseq_ANNOTATED_DOWN_EP$Normalized_MFE)

t.test(SKMEL28_EP_dexseq_ANNOTATED_DOWN_EP$Normalized_MFE, SKMEL28_EP_dexseq_ANNOTATED_UP_EP$Normalized_MFE, paired=FALSE)

# There is a significant increase of the mean MFE of the UP-regulated vs the DOWN-regulated ANNOTATED transcripts
# p-value = 0.0002245

# Show the data with a violin plot
ggplot(SKMEL28_counts_EP, aes(x = factor(log2fold_treated_untreated > 0), y = Normalized_MFE)) +
  geom_violin(fill = "lightgrey", alpha = 0.7) +
  geom_jitter(width = 0.1, size = 1, alpha = 0.3) +
  theme_classic() +
  theme(legend.position = "none", axis.title = element_text(size = 20), plot.title = element_text(size = 10),
        axis.text.x = element_text(size = 16, color = "black")) +  # Increase font size and make x-axis labels black
  geom_vline(xintercept = 0, size = 0.3) +
  labs(y = "Normalized_MFE", x = NULL, title = "SKMEL28_shared_ANNOTATED_TX_early_phase_Log2_FC_adj.P<0.05") +
  scale_x_discrete(labels = c("Log2_FC < 0", "Log2_FC > 0"))
# In this modified code, the x aesthetic is set as a factor variable factor(log2fold_treated_untreated > 0). 
# This creates two separate groups based on the condition Log2_FC > 0 or Log2_FC < 0. Each group is represented by a different violin shape.


# Look at the differences between the 0.75 quantile of UP- and DOWN-regulated transcripts
quantile_075_SKMEL28_UP_EP <- quantile(SKMEL28_EP_dexseq_ANNOTATED_UP_EP$Normalized_MFE, probs = c(0.75)) 
quantile_075_SKMEL28_DOWN_EP <- quantile(SKMEL28_EP_dexseq_ANNOTATED_DOWN_EP$Normalized_MFE, probs = c(0.75)) 
ks_test_SKMEL28_075 <- ks.test(SKMEL28_EP_dexseq_ANNOTATED_UP_EP$Normalized_MFE[SKMEL28_EP_dexseq_ANNOTATED_UP_EP$Normalized_MFE >= quantile_075_SKMEL28_UP_EP], SKMEL28_EP_dexseq_ANNOTATED_DOWN_EP$Normalized_MFE[SKMEL28_EP_dexseq_ANNOTATED_DOWN_EP$Normalized_MFE >= quantile_075_SKMEL28_DOWN_EP])
print(ks_test_SKMEL28_075)

# Look at the differences between the 0.25 quantile of UP- and DOWN-regulated transcripts
quantile_025_SKMEL28_UP_EP <- quantile(SKMEL28_EP_dexseq_ANNOTATED_UP_EP$Normalized_MFE, probs = c(0.25)) 
quantile_025_SKMEL28_DOWN_EP <- quantile(SKMEL28_EP_dexseq_ANNOTATED_DOWN_EP$Normalized_MFE, probs = c(0.25)) 
ks_test_SKMEL28_025 <- ks.test(SKMEL28_EP_dexseq_ANNOTATED_UP_EP$Normalized_MFE[SKMEL28_EP_dexseq_ANNOTATED_UP_EP$Normalized_MFE <= quantile_025_SKMEL28_UP_EP], SKMEL28_EP_dexseq_ANNOTATED_DOWN_EP$Normalized_MFE[SKMEL28_EP_dexseq_ANNOTATED_DOWN_EP$Normalized_MFE <= quantile_025_SKMEL28_DOWN_EP])
print(ks_test_SKMEL28_025)




# Generate a data frame with all the CUB values and plot the data in a density plot
max_length_SKMEL28_MFE <- max(c(length(SKMEL28_EP_dexseq_ANNOTATED_UP_EP$Normalized_MFE), length(SKMEL28_EP_dexseq_ANNOTATED_DOWN_EP$Normalized_MFE)))

data_SKMEL28_MFE_density <- data.frame(MFE_TX_UP_EP_SKMEL28 = c(SKMEL28_EP_dexseq_ANNOTATED_UP_EP$Normalized_MFE, rep(NA, max_length_SKMEL28_MFE - length(SKMEL28_EP_dexseq_ANNOTATED_UP_EP$Normalized_MFE))), # This is useful when we need to create a data frame and the number of values between columns are different
                                       MFE_TX_DOWN_EP_SKMEL28 = c(SKMEL28_EP_dexseq_ANNOTATED_DOWN_EP$Normalized_MFE, rep(NA, max_length_SKMEL28_MFE- length(SKMEL28_EP_dexseq_ANNOTATED_DOWN_EP$Normalized_MFE))))
library(reshape)
data_SKMEL28_MFE_density <- melt(data_SKMEL28_MFE_density)

library(ggplot2)
ggplot(data= data_SKMEL28_MFE_density, aes(x=value, group=variable, color=variable)) +
  geom_density (size=1) +
  theme_classic() +
  #xlim(0,1.5) +
  labs(title = "MFE_distribution_ANNOTATED_Transcripts_SKMEL28_72h")+
  labs(y="Density")+
  labs(x="DGfold / nt (kcal/mol)")









## Analysis net MFE and heat capacity ####
# Determine whether within the differentially expressed genes there is a significant up-regulation of high MFE transcripts:
# Define high MFE transcripts by taking the 75th percentile of the combined MFE distributions of control and drug-treated samples.
# Count the number of transcripts in the significantly up- and down-regulated genes that have MFE values above this threshold.
# Perform a Fisher’s exact test to compare the proportions of high MFE transcripts between the two groups.



SKMEL28_EP_dexseq_ANNOTATED_norm_TX_counts_MFE_df <- SKMEL28_counts_EP %>%
  dplyr::select(ensembl_transcript_id_version, ensembl_transcript_id, ensembl_gene_id, external_gene_name, transcript_is_canonical, Transcript_length,
                SKMEL28_untreated_R1, SKMEL28_untreated_R2, SKMEL28_untreated_R3,
                SKMEL28_vemurafenib_R1, SKMEL28_vemurafenib_R2, SKMEL28_vemurafenib_R3,
                Normalized_MFE, Temperature, Specific_heat) 


# Identify the 75th percentile of the combined MFE distributions of control and drug-treated samples
quantile_075_SKMEL28_EP_dexseq_ANNOTATED_norm_TX_counts_MFE_df <- quantile(SKMEL28_EP_dexseq_ANNOTATED_norm_TX_counts_MFE_df$Normalized_MFE, probs = c(0.75))

# Select all the transcripts with Normalized_MFE values above the 75th percentile
SKMEL28_EP_dexseq_ANNOTATED_high_MFE_TX <- dplyr::filter(SKMEL28_EP_dexseq_ANNOTATED_norm_TX_counts_MFE_df, Normalized_MFE > quantile_075_SKMEL28_EP_dexseq_ANNOTATED_norm_TX_counts_MFE_df)

# Count the number of transcripts in the significantly up- and down-regulated pools that have MFE values above this threshold.
SKMEL28_EP_dexseq_ANNOTATED_high_MFE_TX_in_UP <- dplyr::filter(SKMEL28_EP_dexseq_ANNOTATED_UP_EP, ensembl_transcript_id_version %in% SKMEL28_EP_dexseq_ANNOTATED_high_MFE_TX$ensembl_transcript_id_version)
SKMEL28_EP_dexseq_ANNOTATED_high_MFE_TX_in_DOWN <- dplyr::filter(SKMEL28_EP_dexseq_ANNOTATED_DOWN_EP, ensembl_transcript_id_version %in% SKMEL28_EP_dexseq_ANNOTATED_high_MFE_TX$ensembl_transcript_id_version)


# Perform a Fisher’s exact test to compare the proportions of high MFE transcripts between the two groups.
SKMEL28_EP_dexseq_ANNOTATED_total_UP <- nrow(SKMEL28_EP_dexseq_ANNOTATED_UP_EP)
SKMEL28_EP_dexseq_ANNOTATED_total_DOWN <- nrow(SKMEL28_EP_dexseq_ANNOTATED_DOWN_EP)
SKMEL28_EP_dexseq_ANNOTATED_count_high_MFE_UP <- nrow(SKMEL28_EP_dexseq_ANNOTATED_high_MFE_TX_in_UP)
SKMEL28_EP_dexseq_ANNOTATED_count_high_MFE_DOWN <- nrow(SKMEL28_EP_dexseq_ANNOTATED_high_MFE_TX_in_DOWN)

SKMEL28_EP_dexseq_ANNOTATED_contingency_table <- matrix(c(SKMEL28_EP_dexseq_ANNOTATED_count_high_MFE_UP, SKMEL28_EP_dexseq_ANNOTATED_total_UP - SKMEL28_EP_dexseq_ANNOTATED_count_high_MFE_UP,
                                                     SKMEL28_EP_dexseq_ANNOTATED_count_high_MFE_DOWN, SKMEL28_EP_dexseq_ANNOTATED_total_DOWN - SKMEL28_EP_dexseq_ANNOTATED_count_high_MFE_DOWN),
                                                   nrow = 2, byrow = TRUE,
                                                   dimnames = list(UP = c("UP-regulated", "DOWN-regulated"),
                                                                   DOWN = c("High MFE TX", "Not High MFE TX")))

SKMEL28_EP_dexseq_ANNOTATED_contingency_table_df <- as.data.frame(SKMEL28_EP_dexseq_ANNOTATED_contingency_table)

SKMEL28_EP_dexseq_ANNOTATED_fisher_result <- fisher.test(SKMEL28_EP_dexseq_ANNOTATED_contingency_table)

cat("SKMEL28_counts_EP P-value from Fisher test:", SKMEL28_EP_dexseq_ANNOTATED_fisher_result$p.value, "\n")


# Look at the fractions
SKMEL28_EP_dexseq_ANNOTATED_high_MFE_in_UP_fraction <- SKMEL28_EP_dexseq_ANNOTATED_count_high_MFE_UP / SKMEL28_EP_dexseq_ANNOTATED_total_UP
SKMEL28_EP_dexseq_ANNOTATED_high_MFE_in_DOWN_fraction <- SKMEL28_EP_dexseq_ANNOTATED_count_high_MFE_DOWN / SKMEL28_EP_dexseq_ANNOTATED_total_DOWN

# Output the results
cat("Fraction of SKMEL28_counts_EP high MFE TX in UP-regulated genes:", SKMEL28_EP_dexseq_ANNOTATED_high_MFE_in_UP_fraction, "\n") # This value goes in the summary table
cat("Fraction of SKMEL28_counts_EP high MFE TX in DOWN-regulated genes:", SKMEL28_EP_dexseq_ANNOTATED_high_MFE_in_DOWN_fraction, "\n") # This value goes in the summary table
cat("Mean of SKMEL28_counts_EP high MFE TX in UP-regulated genes:", mean(SKMEL28_EP_dexseq_ANNOTATED_high_MFE_TX_in_UP$Normalized_MFE), "\n")
cat("Mean of SKMEL28_counts_EP high MFE TX in DOWN-regulated genes:", mean(SKMEL28_EP_dexseq_ANNOTATED_high_MFE_TX_in_DOWN$Normalized_MFE), "\n")
cat("Mode of SKMEL28_counts_EP high MFE TX in UP-regulated genes:", calculate_mode(SKMEL28_EP_dexseq_ANNOTATED_high_MFE_TX_in_UP$Normalized_MFE), "\n")
cat("Mode of SKMEL28_counts_EP high MFE TX in DOWN-regulated genes:", calculate_mode(SKMEL28_EP_dexseq_ANNOTATED_high_MFE_TX_in_DOWN$Normalized_MFE), "\n")





# Perform a two-proportion z-test to compare proportions
# Calculate the proportions and total counts
SKMEL28_EP_dexseq_ANNOTATED_high_MFE_in_UP_proportion <- SKMEL28_EP_dexseq_ANNOTATED_high_MFE_in_UP_fraction
SKMEL28_EP_dexseq_ANNOTATED_high_MFE_in_DOWN_proportion <- SKMEL28_EP_dexseq_ANNOTATED_high_MFE_in_DOWN_fraction

# Combined proportion for pooled standard error
combined_proportion <- (SKMEL28_EP_dexseq_ANNOTATED_count_high_MFE_UP + SKMEL28_EP_dexseq_ANNOTATED_count_high_MFE_DOWN) /
  (SKMEL28_EP_dexseq_ANNOTATED_total_UP + SKMEL28_EP_dexseq_ANNOTATED_total_DOWN)

# Standard error of the difference in proportions
standard_error <- sqrt(
  combined_proportion * (1 - combined_proportion) * 
    (1 / SKMEL28_EP_dexseq_ANNOTATED_total_UP + 1 / SKMEL28_EP_dexseq_ANNOTATED_total_DOWN)
)

# Z-statistic
z_statistic <- (SKMEL28_EP_dexseq_ANNOTATED_high_MFE_in_UP_proportion - SKMEL28_EP_dexseq_ANNOTATED_high_MFE_in_DOWN_proportion) /
  standard_error

# Two-tailed p-value
SKMEL28_EP_dexseq_ANNOTATED_z_test_MFE_p_value <- 2 * (1 - pnorm(abs(z_statistic)))








# Show data as violin plot
# Add a new column to indicate if the value is above the 75th percentile
SKMEL28_counts_EP$above_quantile_075 <- SKMEL28_counts_EP$Normalized_MFE > quantile_075_SKMEL28_EP_dexseq_ANNOTATED_norm_TX_counts_MFE_df

# Generate violin plot with conditional jitter coloring
ggplot(SKMEL28_counts_EP, aes(x = factor(log2fold_treated_untreated > 0), y = Normalized_MFE)) +
  geom_violin(fill = "white") +
  geom_jitter(aes(color = above_quantile_075), width = 0.1, size = 1, alpha = 0.3) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "blue")) +
  theme_classic() +
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
    title = "SKMEL28_EP_dexseq_ANNOTATED_Log2_FC_adj.P<0.05"
  ) +
  scale_x_discrete(labels = c("log2FoldChange < 0", "log2FoldChange > 0"))













# Determine the impact of these changes in terms of fraction of total counts:
# Calculate the fraction of total counts of high MFE transcripts in control and drug-treated groups. 
# This is done by summing the counts of the high MFE transcripts and dividing by the total counts in each group.
# Fisher’s exact test to compare the proportions of high MFE transcripts counts between the two groups.


# Convert from wide to long format
SKMEL28_EP_dexseq_ANNOTATED_long_data <- SKMEL28_EP_dexseq_ANNOTATED_norm_TX_counts_MFE_df %>%
  pivot_longer(
    cols = starts_with("SKMEL28"),      # Select all columns starting with "countData.72h"
    names_to = "replicate",         # Name of the new column for replicate identifiers
    values_to = "counts"            # Name of the new column for the counts
  )

# Add a condition column based on the replicate
SKMEL28_EP_dexseq_ANNOTATED_long_data <- SKMEL28_EP_dexseq_ANNOTATED_long_data %>%
  mutate(condition = ifelse(replicate %in% c("SKMEL28_untreated_R1", "SKMEL28_untreated_R2", "SKMEL28_untreated_R3"), "Control", "vemurafenib"))


# Aggregate the replicates for each condition
SKMEL28_EP_dexseq_ANNOTATED_norm_TX_aggregated_data_avg_counts <- SKMEL28_EP_dexseq_ANNOTATED_long_data %>%
  group_by(ensembl_transcript_id_version, condition, Normalized_MFE, Specific_heat) %>%
  summarise(
    avg_counts = mean(counts),
    .groups = 'drop'
  )


# Define the 75th percentile of the entire distribution of MFE values for the SKMEL28_EP_dexseq_ANNOTATED data
SKMEL28_EP_dexseq_ANNOTATED_MFE_all_TX_75_percentile <- quantile(SKMEL28_EP_dexseq_ANNOTATED_norm_TX_aggregated_data_avg_counts$Normalized_MFE, 0.75)

# Identify all the transcripts with MFE values above the 75th percentile
SKMEL28_EP_dexseq_ANNOTATED_high_MFE_TX_total <- dplyr::filter(SKMEL28_EP_dexseq_ANNOTATED_norm_TX_counts_MFE_df, Normalized_MFE > SKMEL28_EP_dexseq_ANNOTATED_MFE_all_TX_75_percentile)





# Fraction of total high MFE transcripts counts in control and drug-treated groups
SKMEL28_EP_dexseq_ANNOTATED_control_total_counts <- dplyr::filter(SKMEL28_EP_dexseq_ANNOTATED_norm_TX_aggregated_data_avg_counts, condition == "Control") %>%
  dplyr::summarise(total_counts = sum(avg_counts))

SKMEL28_EP_dexseq_ANNOTATED_control_high_MFE_TX_counts <- dplyr::filter(SKMEL28_EP_dexseq_ANNOTATED_norm_TX_aggregated_data_avg_counts, condition == "Control") %>%
  dplyr::filter(ensembl_transcript_id_version %in% SKMEL28_EP_dexseq_ANNOTATED_high_MFE_TX_total$ensembl_transcript_id_version) %>%
  dplyr::summarise(total_counts_high_MFE_TX = sum(avg_counts))

SKMEL28_EP_dexseq_ANNOTATED_control_high_MFE_TX_counts_fraction_EP <- SKMEL28_EP_dexseq_ANNOTATED_control_high_MFE_TX_counts$total_counts_high_MFE_TX / SKMEL28_EP_dexseq_ANNOTATED_control_total_counts$total_counts
print(SKMEL28_EP_dexseq_ANNOTATED_control_high_MFE_TX_counts_fraction_EP) # This value goes in the summary table





SKMEL28_EP_dexseq_ANNOTATED_vemurafenib_total_counts <- dplyr::filter(SKMEL28_EP_dexseq_ANNOTATED_norm_TX_aggregated_data_avg_counts, condition == "vemurafenib") %>%
  dplyr::summarise(total_counts = sum(avg_counts))

SKMEL28_EP_dexseq_ANNOTATED_vemurafenib_high_MFE_TX_counts <- dplyr::filter(SKMEL28_EP_dexseq_ANNOTATED_norm_TX_aggregated_data_avg_counts, condition == "vemurafenib") %>%
  dplyr::filter(ensembl_transcript_id_version %in% SKMEL28_EP_dexseq_ANNOTATED_high_MFE_TX_total$ensembl_transcript_id_version) %>%
  dplyr::summarise(total_counts_high_MFE_TX = sum(avg_counts))

SKMEL28_EP_dexseq_ANNOTATED_vemurafenib_high_MFE_TX_counts_fraction_EP <- SKMEL28_EP_dexseq_ANNOTATED_vemurafenib_high_MFE_TX_counts$total_counts_high_MFE_TX / SKMEL28_EP_dexseq_ANNOTATED_vemurafenib_total_counts$total_counts
print(SKMEL28_EP_dexseq_ANNOTATED_vemurafenib_high_MFE_TX_counts_fraction_EP) # This value goes in the summary table





# Fisher's exact test to compare the proportions of high MFE transcripts counts between the two groups
# Create the groups
SKMEL28_EP_dexseq_ANNOTATED_high_MFE_TX_counts_control_CT <- SKMEL28_EP_dexseq_ANNOTATED_control_high_MFE_TX_counts_fraction_EP * SKMEL28_EP_dexseq_ANNOTATED_control_total_counts$total_counts
SKMEL28_EP_dexseq_ANNOTATED_low_MFE_TX_counts_control_CT <- SKMEL28_EP_dexseq_ANNOTATED_control_total_counts$total_counts - SKMEL28_EP_dexseq_ANNOTATED_control_high_MFE_TX_counts$total_counts_high_MFE_TX

SKMEL28_EP_dexseq_ANNOTATED_high_MFE_TX_counts_vemurafenib_CT <- SKMEL28_EP_dexseq_ANNOTATED_vemurafenib_high_MFE_TX_counts_fraction_EP * SKMEL28_EP_dexseq_ANNOTATED_vemurafenib_total_counts$total_counts
SKMEL28_EP_dexseq_ANNOTATED_low_MFE_TX_counts_vemurafenib_CT <- SKMEL28_EP_dexseq_ANNOTATED_vemurafenib_total_counts$total_counts - SKMEL28_EP_dexseq_ANNOTATED_vemurafenib_high_MFE_TX_counts$total_counts_high_MFE_TX

# Construct the contingency table with raw counts
SKMEL28_EP_dexseq_ANNOTATED_counts_contingency_table <- matrix(
  c(
    SKMEL28_EP_dexseq_ANNOTATED_control_high_MFE_TX_counts$total_counts_high_MFE_TX, 
    SKMEL28_EP_dexseq_ANNOTATED_vemurafenib_high_MFE_TX_counts$total_counts_high_MFE_TX,
    SKMEL28_EP_dexseq_ANNOTATED_control_total_counts$total_counts - SKMEL28_EP_dexseq_ANNOTATED_control_high_MFE_TX_counts$total_counts_high_MFE_TX, 
    SKMEL28_EP_dexseq_ANNOTATED_vemurafenib_total_counts$total_counts - SKMEL28_EP_dexseq_ANNOTATED_vemurafenib_high_MFE_TX_counts$total_counts_high_MFE_TX
  ),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(Category = c("High_MFE_TX_counts", "Low_MFE_TX_counts"),
                  Group = c("Control", "vemurafenib"))
)

# Ensure the counts are rounded to integers (they should already be integers if the input data is correct)
# SKMEL28_EP_dexseq_ANNOTATED_counts_contingency_table <- round(SKMEL28_EP_dexseq_ANNOTATED_counts_contingency_table)

SKMEL28_EP_dexseq_ANNOTATED_counts_contingency_table_df <- as.data.frame(SKMEL28_EP_dexseq_ANNOTATED_counts_contingency_table)

SKMEL28_EP_dexseq_ANNOTATED_fisher_result_counts <- fisher.test(SKMEL28_EP_dexseq_ANNOTATED_counts_contingency_table_df)
print(SKMEL28_EP_dexseq_ANNOTATED_fisher_result_counts)


cat("SKMEL28_EP_dexseq_ANNOTATED Fisher P-value Fisher proportions of high MFE transcripts counts in drug-treated and control :", SKMEL28_EP_dexseq_ANNOTATED_fisher_result_counts$p.value, "\n")  # This value goes in the summary table


# Output the results
cat("SKMEL28_EP_dexseq_ANNOTATED control high MFE TX as fraction of total counts:", SKMEL28_EP_dexseq_ANNOTATED_control_high_MFE_TX_counts_fraction_EP, "\n")
cat("SKMEL28_EP_dexseq_ANNOTATED vemurafenib high MFE TX as fraction of total counts:", SKMEL28_EP_dexseq_ANNOTATED_vemurafenib_high_MFE_TX_counts_fraction_EP, "\n")
cat("SKMEL28_EP_dexseq_ANNOTATED net change counts high MFE TX drug vs control", SKMEL28_EP_dexseq_ANNOTATED_vemurafenib_high_MFE_TX_counts_fraction_EP-SKMEL28_EP_dexseq_ANNOTATED_control_high_MFE_TX_counts_fraction_EP, "\n")




# Perform a two-proportion z-test instead of Fisher's exact test
# Extract the raw counts for high MFE transcripts in both groups
high_MFE_control <- SKMEL28_EP_dexseq_ANNOTATED_control_high_MFE_TX_counts$total_counts_high_MFE_TX
high_MFE_vemurafenib <- SKMEL28_EP_dexseq_ANNOTATED_vemurafenib_high_MFE_TX_counts$total_counts_high_MFE_TX

# Total counts in both groups
total_control <- SKMEL28_EP_dexseq_ANNOTATED_control_total_counts$total_counts
total_vemurafenib <- SKMEL28_EP_dexseq_ANNOTATED_vemurafenib_total_counts$total_counts

# Perform the two-proportion z-test
z_test_result_high <- prop.test(
  c(high_MFE_control, high_MFE_vemurafenib), 
  c(total_control, total_vemurafenib)
)

# Output the result
SKMEL28_EP_dexseq_ANNOTATED_z_test_counts_75_p_value <- z_test_result_high$p.value






# Plot Normalized_MFE vs transcript counts of both control and drug-treated conditions
# Create separate data frames for Control and vemurafenib conditions
SKMEL28_EP_dexseq_ANNOTATED_control_data <- SKMEL28_EP_dexseq_ANNOTATED_norm_TX_aggregated_data_avg_counts %>%
  dplyr::filter(condition == "Control")

SKMEL28_EP_dexseq_ANNOTATED_vemurafenib_data <- SKMEL28_EP_dexseq_ANNOTATED_norm_TX_aggregated_data_avg_counts %>%
  dplyr::filter(condition == "vemurafenib") 


# Generate the plot
Plot_SKMEL28_EP_dexseq_ANNOTATED_control <- ggplot(data = SKMEL28_EP_dexseq_ANNOTATED_control_data, aes(x = avg_counts, y = Normalized_MFE, color = condition)) +
  geom_point(alpha = 0.6) +
  theme_classic() +
  labs(
    title = "SKMEL28_EP_dexseq_ANNOTATED MFE vs. total avg counts Control",
    x = "Transcript Counts",
    y = "Normalized MFE (kcal/mol)"
  ) +
  #xlim(0, 150000) +
  ylim(-0.8, 0) +
  theme(
    axis.title = element_text(size = 15),
    plot.title = element_text(size = 10)
  ) +
  guides(color = FALSE)


Plot_SKMEL28_EP_dexseq_ANNOTATED_vemurafenib <- ggplot(data = SKMEL28_EP_dexseq_ANNOTATED_vemurafenib_data, aes(x = avg_counts, y = Normalized_MFE, color = condition)) +
  geom_point(alpha = 0.6) +
  theme_classic() +
  labs(
    title = "SKMEL28_EP_dexseq_ANNOTATED MFE vs. total avg counts vemurafenib",
    x = "Transcript Counts",
    y = "Normalized MFE (kcal/mol)"
  ) +
  #xlim(0, 150000) +
  ylim(-0.8, 0) +
  theme(
    axis.title = element_text(size = 15),
    plot.title = element_text(size = 10)
  ) +
  guides(color = FALSE)

library(patchwork)
SKMEL28_EP_dexseq_ANNOTATED_combined_plots <- Plot_SKMEL28_EP_dexseq_ANNOTATED_control + Plot_SKMEL28_EP_dexseq_ANNOTATED_vemurafenib
SKMEL28_EP_dexseq_ANNOTATED_combined_plots








# Define the 25th percentile of the entire distribution of MFE values for the SKMEL28_EP_dexseq_ANNOTATED data
SKMEL28_EP_dexseq_ANNOTATED_MFE_all_TX_25_percentile <- quantile(SKMEL28_EP_dexseq_ANNOTATED_norm_TX_aggregated_data_avg_counts$Normalized_MFE, 0.25)

# Identify all the transcripts with MFE values below the 25th percentile
SKMEL28_EP_dexseq_ANNOTATED_low_MFE_TX_total <- dplyr::filter(SKMEL28_EP_dexseq_ANNOTATED_norm_TX_counts_MFE_df, Normalized_MFE < SKMEL28_EP_dexseq_ANNOTATED_MFE_all_TX_25_percentile)





# Fraction of total low MFE transcripts counts in control and drug-treated groups
SKMEL28_EP_dexseq_ANNOTATED_control_total_counts <- dplyr::filter(SKMEL28_EP_dexseq_ANNOTATED_norm_TX_aggregated_data_avg_counts, condition == "Control") %>%
  dplyr::summarise(total_counts = sum(avg_counts))

SKMEL28_EP_dexseq_ANNOTATED_control_low_MFE_TX_counts <- dplyr::filter(SKMEL28_EP_dexseq_ANNOTATED_norm_TX_aggregated_data_avg_counts, condition == "Control") %>%
  dplyr::filter(ensembl_transcript_id_version %in% SKMEL28_EP_dexseq_ANNOTATED_low_MFE_TX_total$ensembl_transcript_id_version) %>%
  dplyr::summarise(total_counts_low_MFE_TX = sum(avg_counts))

SKMEL28_EP_dexseq_ANNOTATED_control_low_MFE_TX_counts_fraction_EP <- SKMEL28_EP_dexseq_ANNOTATED_control_low_MFE_TX_counts$total_counts_low_MFE_TX / SKMEL28_EP_dexseq_ANNOTATED_control_total_counts$total_counts
print(SKMEL28_EP_dexseq_ANNOTATED_control_low_MFE_TX_counts_fraction_EP) # This value goes in the summary table





SKMEL28_EP_dexseq_ANNOTATED_vemurafenib_total_counts <- dplyr::filter(SKMEL28_EP_dexseq_ANNOTATED_norm_TX_aggregated_data_avg_counts, condition == "vemurafenib") %>%
  dplyr::summarise(total_counts = sum(avg_counts))

SKMEL28_EP_dexseq_ANNOTATED_vemurafenib_low_MFE_TX_counts <- dplyr::filter(SKMEL28_EP_dexseq_ANNOTATED_norm_TX_aggregated_data_avg_counts, condition == "vemurafenib") %>%
  dplyr::filter(ensembl_transcript_id_version %in% SKMEL28_EP_dexseq_ANNOTATED_low_MFE_TX_total$ensembl_transcript_id_version) %>%
  dplyr::summarise(total_counts_low_MFE_TX = sum(avg_counts))

SKMEL28_EP_dexseq_ANNOTATED_vemurafenib_low_MFE_TX_counts_fraction_EP <- SKMEL28_EP_dexseq_ANNOTATED_vemurafenib_low_MFE_TX_counts$total_counts_low_MFE_TX / SKMEL28_EP_dexseq_ANNOTATED_vemurafenib_total_counts$total_counts
print(SKMEL28_EP_dexseq_ANNOTATED_vemurafenib_low_MFE_TX_counts_fraction_EP) # This value goes in the summary table





# Fisher's exact test to compare the proportions of low MFE transcripts counts between the two groups
# Create the groups
SKMEL28_EP_dexseq_ANNOTATED_low_MFE_TX_counts_control_CT <- SKMEL28_EP_dexseq_ANNOTATED_control_low_MFE_TX_counts_fraction_EP * SKMEL28_EP_dexseq_ANNOTATED_control_total_counts$total_counts
SKMEL28_EP_dexseq_ANNOTATED_low_MFE_TX_counts_control_CT <- SKMEL28_EP_dexseq_ANNOTATED_control_total_counts$total_counts - SKMEL28_EP_dexseq_ANNOTATED_control_low_MFE_TX_counts$total_counts_low_MFE_TX

SKMEL28_EP_dexseq_ANNOTATED_low_MFE_TX_counts_vemurafenib_CT <- SKMEL28_EP_dexseq_ANNOTATED_vemurafenib_low_MFE_TX_counts_fraction_EP * SKMEL28_EP_dexseq_ANNOTATED_vemurafenib_total_counts$total_counts
SKMEL28_EP_dexseq_ANNOTATED_low_MFE_TX_counts_vemurafenib_CT <- SKMEL28_EP_dexseq_ANNOTATED_vemurafenib_total_counts$total_counts - SKMEL28_EP_dexseq_ANNOTATED_vemurafenib_low_MFE_TX_counts$total_counts_low_MFE_TX

# Construct the contingency table with raw counts
SKMEL28_EP_dexseq_ANNOTATED_counts_contingency_table_low <- matrix(
  c(
    SKMEL28_EP_dexseq_ANNOTATED_control_low_MFE_TX_counts$total_counts_low_MFE_TX, 
    SKMEL28_EP_dexseq_ANNOTATED_vemurafenib_low_MFE_TX_counts$total_counts_low_MFE_TX,
    SKMEL28_EP_dexseq_ANNOTATED_control_total_counts$total_counts - SKMEL28_EP_dexseq_ANNOTATED_control_low_MFE_TX_counts$total_counts_low_MFE_TX, 
    SKMEL28_EP_dexseq_ANNOTATED_vemurafenib_total_counts$total_counts - SKMEL28_EP_dexseq_ANNOTATED_vemurafenib_low_MFE_TX_counts$total_counts_low_MFE_TX
  ),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(Category = c("Low_MFE_TX_counts", "Other_TX_counts"),
                  Group = c("Control", "vemurafenib"))
)

SKMEL28_EP_dexseq_ANNOTATED_counts_contingency_table_df_low <- as.data.frame(SKMEL28_EP_dexseq_ANNOTATED_counts_contingency_table_low)

SKMEL28_EP_dexseq_ANNOTATED_fisher_result_counts_low <- fisher.test(SKMEL28_EP_dexseq_ANNOTATED_counts_contingency_table_low)
print(SKMEL28_EP_dexseq_ANNOTATED_fisher_result_counts_low)


cat("SKMEL28_EP_dexseq_ANNOTATED Fisher P-value Fisher proportions of low MFE transcripts counts in drug-treated and control :", SKMEL28_EP_dexseq_ANNOTATED_fisher_result_counts$p.value, "\n")  # This value goes in the summary table


# Output the results
cat("SKMEL28_EP_dexseq_ANNOTATED control low MFE TX as fraction of total counts:", SKMEL28_EP_dexseq_ANNOTATED_control_low_MFE_TX_counts_fraction_EP, "\n")
cat("SKMEL28_EP_dexseq_ANNOTATED palbociclib low MFE TX as fraction of total counts:", SKMEL28_EP_dexseq_ANNOTATED_vemurafenib_low_MFE_TX_counts_fraction_EP, "\n")
cat("SKMEL28_EP_dexseq_ANNOTATED net change counts low MFE TX drug vs control", SKMEL28_EP_dexseq_ANNOTATED_vemurafenib_low_MFE_TX_counts_fraction_EP-SKMEL28_EP_dexseq_ANNOTATED_control_low_MFE_TX_counts_fraction_EP, "\n")




# Perform a two-proportion z-test instead of Fisher's exact test
# Extract the raw counts for low MFE transcripts in both groups
low_MFE_control <- SKMEL28_EP_dexseq_ANNOTATED_control_low_MFE_TX_counts$total_counts_low_MFE_TX
low_MFE_vemurafenib <- SKMEL28_EP_dexseq_ANNOTATED_vemurafenib_low_MFE_TX_counts$total_counts_low_MFE_TX

# Total counts in both groups
total_control <- SKMEL28_EP_dexseq_ANNOTATED_control_total_counts$total_counts
total_vemurafenib <- SKMEL28_EP_dexseq_ANNOTATED_vemurafenib_total_counts$total_counts

# Perform the two-proportion z-test
z_test_result_low <- prop.test(
  c(low_MFE_control, low_MFE_vemurafenib), 
  c(total_control, total_vemurafenib)
)

# Output the result
SKMEL28_EP_dexseq_ANNOTATED_z_test_counts_25_p_value <- z_test_result_low$p.value






## SKMEL28_counts_EP heat capacity changes ####

# This is calculated on the total n. of transcripts with a significant log2 fold change
# Generate a scatter plot Log2_FC vs heat capacity
ggplot(SKMEL28_counts_EP, aes(y=Specific_heat, x=log2fold_treated_untreated)) +
  geom_point(size=1.5, alpha=0.1)+
  theme_classic() +
  theme(legend.position="none") +
  geom_vline(xintercept = 0, size=0.3) +
  xlim(-12, 12) +
  ylim(0,400) +
  theme(axis.title = element_text(size = 20)) +
  theme(plot.title = element_text(size = 15)) +
  labs(y="Heat capacity (J/C)",x="Log2_FC_expression",title="SKMEL28_EP_dexseq_ANNOTATED_Log2_FC_adj.P<0.05")


# Show the data with a box plot
ggplot(SKMEL28_counts_EP, aes(x = factor(log2fold_treated_untreated > 0), y = Specific_heat)) +
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
  labs(y="Heat capacity (J/C)",x="",title="SKMEL28_EP_dexseq_ANNOTATED_Log2_FC_adj.P<0.05 heat_capacity all transcripts") +
  scale_x_discrete(labels = c("log2FoldChange < 0", "log2FoldChange > 0"))

Mean_heat_capacity_SKMEL28_EP_dexseq_ANNOTATED_UP <- mean(SKMEL28_EP_dexseq_ANNOTATED_UP_EP$Specific_heat)
Mean_heat_capacity_SKMEL28_EP_dexseq_ANNOTATED_DOWN <- mean(SKMEL28_EP_dexseq_ANNOTATED_DOWN_EP$Specific_heat)
TTEST_heat_capacity_SKMEL28_EP_dexseq_ANNOTATED_all_UP_vs_DOWN <- t.test(SKMEL28_EP_dexseq_ANNOTATED_DOWN_EP$Specific_heat, SKMEL28_EP_dexseq_ANNOTATED_UP_EP$Specific_heat, paired=FALSE)
print(TTEST_heat_capacity_SKMEL28_EP_dexseq_ANNOTATED_all_UP_vs_DOWN)
# Add K-S test to look at the differences between the up-regukated and down-regulated genes
ks_test_SKMEL28_EP_dexseq_ANNOTATED_all_UP_vs_DOWN_heat_capacity <- ks.test(SKMEL28_EP_dexseq_ANNOTATED_DOWN_EP$Specific_heat, SKMEL28_EP_dexseq_ANNOTATED_UP_EP$Specific_heat)
print(ks_test_SKMEL28_EP_dexseq_ANNOTATED_all_UP_vs_DOWN_heat_capacity)




# Calculate the changes of heat capacity only in the high MFE transcripts

# Contribution of the UP or DOWN-regulated high MFE transcripts to the heat capacity
# Look at the heat capacity of the UP-regulated high MFE transcripts
SKMEL28_EP_dexseq_ANNOTATED_high_MFE_TX_in_UP_heat_capacity <- mean(SKMEL28_EP_dexseq_ANNOTATED_high_MFE_TX_in_UP$Specific_heat)

# Look at the heat capacity of the DOWN-regulated high MFE transcripts
SKMEL28_EP_dexseq_ANNOTATED_high_MFE_TX_in_DOWN_heat_capacity <- mean(SKMEL28_EP_dexseq_ANNOTATED_high_MFE_TX_in_DOWN$Specific_heat)

# Heat capacity gain in UP-regulated high MFE transcripts vs DOWN-regulated high MFE transcripts
SKMEL28_EP_dexseq_ANNOTATED_heat_capacity_gain_UP_vs_DOWN <- SKMEL28_EP_dexseq_ANNOTATED_high_MFE_TX_in_UP_heat_capacity - SKMEL28_EP_dexseq_ANNOTATED_high_MFE_TX_in_DOWN_heat_capacity
cat("SKMEL28_EP_dexseq_ANNOTATED net heat capacity UP vs DOWN high MFE TX:", SKMEL28_EP_dexseq_ANNOTATED_heat_capacity_gain_UP_vs_DOWN, "\n")


# Extract specific heat data
SKMEL28_EP_dexseq_ANNOTATED_Specific_heat_UP <- SKMEL28_EP_dexseq_ANNOTATED_high_MFE_TX_in_UP$Specific_heat
SKMEL28_EP_dexseq_ANNOTATED_Specific_heat_DOWN <- SKMEL28_EP_dexseq_ANNOTATED_high_MFE_TX_in_DOWN$Specific_heat

# Combine data into a single data frame for easier plotting
SKMEL28_EP_dexseq_ANNOTATED_Specific_heat_data <- data.frame(
  SKMEL28_EP_dexseq_ANNOTATED_Specific_heat = c(SKMEL28_EP_dexseq_ANNOTATED_Specific_heat_UP, SKMEL28_EP_dexseq_ANNOTATED_Specific_heat_DOWN),
  regulation = factor(c(rep("UP", length(SKMEL28_EP_dexseq_ANNOTATED_Specific_heat_UP)), rep("DOWN", length(SKMEL28_EP_dexseq_ANNOTATED_Specific_heat_DOWN))))
)

# Generate box plot
ggplot(SKMEL28_EP_dexseq_ANNOTATED_Specific_heat_data, aes(x = regulation, y = SKMEL28_EP_dexseq_ANNOTATED_Specific_heat, fill = regulation)) +
  geom_boxplot(fill = "lightgrey", alpha = 0.7) +
  labs(title = "SKMEL28_EP_dexseq_ANNOTATED specific Heat of High MFE Transcripts",
       x = "",
       y = "Specific Heat (J/C)") +
  scale_x_discrete(labels = c("log2_FC < 0", "log2_FC > 0")) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 16),
        axis.title = element_text(size = 25))

Mean_heat_capacity_SKMEL28_EP_dexseq_ANNOTATED_high_MFE_TX_UP <- mean(SKMEL28_EP_dexseq_ANNOTATED_high_MFE_TX_in_UP$Specific_heat)
Mean_heat_capacity_SKMEL28_EP_dexseq_ANNOTATED_high_MFE_TX_DOWN <- mean(SKMEL28_EP_dexseq_ANNOTATED_high_MFE_TX_in_DOWN$Specific_heat)
TTEST_heat_capacity_SKMEL28_EP_dexseq_ANNOTATED_all_high_MFE_TX_UP_vs_DOWN <- t.test(SKMEL28_EP_dexseq_ANNOTATED_high_MFE_TX_in_DOWN$Specific_heat, SKMEL28_EP_dexseq_ANNOTATED_high_MFE_TX_in_UP$Specific_heat, paired=FALSE)
print(TTEST_heat_capacity_SKMEL28_EP_dexseq_ANNOTATED_all_high_MFE_TX_UP_vs_DOWN)

# Add K-S test to look at the differences between the up-regukated and down-regulated genes
ks_test_heat_capacity_SKMEL28_EP_dexseq_ANNOTATED_all_high_MFE_TX_UP_vs_DOWN <- ks.test(SKMEL28_EP_dexseq_ANNOTATED_high_MFE_TX_in_DOWN$Specific_heat, SKMEL28_EP_dexseq_ANNOTATED_high_MFE_TX_in_UP$Specific_heat)
print(ks_test_heat_capacity_SKMEL28_EP_dexseq_ANNOTATED_all_high_MFE_TX_UP_vs_DOWN)











# Analysis cell proliferation Nanopore data SKMEL28 VEM ####

SKMEL28_VEM_DESeq2 <- read.delim(file="path_to/early.phase_SKMEL28_vemurafenib.tx_counts.both_conditions.treated_untreated_deseq2.tab", header = TRUE, sep = "\t", quote = "")



# Filter the proliferation gene sets by using a function
# Use a function also to determine the mean and standard deviation of each set and the plot all of them together

# Generate the function that filters each gene set in the target data frame
# This function will be used in the following commands
filtering <- function(x) {
  Filtering <- dplyr::filter(SKMEL28_VEM_DESeq2, gene_name %in% x)
  return(Filtering)
}

# Generate the list containing all the gene sets to be used for the filtering.
CURATED_PATHWAYS.input <- list(MTOR_positive_regulators, MTOR_negative_regulators, Proliferation_Science_paper, Meta_PCNA, Proliferation_GO_genes_Magill_Sack_2018, Proliferation_STOP_genes_Magill_Sack_2018, Reactome_G1_phase, Proliferation_Carey_2021, Essential_cell_cycle_Viner_Breuer_Cell_Division_2019, Essential_S_Phase_Viner_Breuer_Cell_Division_2019, Homologous_recombination, Fanconi_anemia, High_fidelity_polymerases)

# This creates the names that will be used for the final data frames of each file resulting from the filtering of the target data frame with each gene set
CURATED_PATHWAYS.output <- c("SKMEL28_VEM_MTOR_POS_EP", "SKMEL28_VEM_MTOR_NEG_EP", "SKMEL28_VEM_Proliferation_Science_EP", "SKMEL28_VEM_Meta_PCNA_EP", "SKMEL28_VEM_Proliferation_GO_EP", "SKMEL28_VEM_Proliferation_STOP_EP", "SKMEL28_VEM_G1_phase_EP", "SKMEL28_VEM_Proliferation_Carey_EP", "SKMEL28_VEM_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019_EP", "SKMEL28_VEM_Essential_S_Phase_Viner_Breuer_Cell_Division_2019_EP", "SKMEL28_VEM_HR_EP", "SKMEL28_VEM_Fanconi_anemia_EP", "SKMEL28_VEM_HF_POL_EP")

# Assign the names of the gene sets to each entry in the list
names(CURATED_PATHWAYS.input) <- CURATED_PATHWAYS.output

# Apply the filtering function to create a list of data frames
Curated_test_df <- lapply(CURATED_PATHWAYS.input, function(x) filtering(x))

# creating a vector with desired names of filtered data frames
Curated_dfs.output <- c("df_SKMEL28_VEM_MTOR_POS_EP", "df_SKMEL28_VEM_MTOR_NEG_EP", "df_SKMEL28_VEM_Proliferation_Science_EP", "df_SKMEL28_VEM_Meta_PCNA_EP", "df_SKMEL28_VEM_Proliferation_GO_EP", "df_SKMEL28_VEM_Proliferation_STOP_EP", "df_SKMEL28_VEM_G1_phase_EP", "df_SKMEL28_VEM_Proliferation_Carey_EP", "df_SKMEL28_VEM_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019_EP", "df_SKMEL28_VEM_Essential_S_Phase_Viner_Breuer_Cell_Division_2019_EP", "df_SKMEL28_VEM_HR_EP", "df_SKMEL28_VEM_Fanconi_anemia_EP", "df_SKMEL28_VEM_HF_POL_EP")

# loop to assign data frames to the desired variable name
for (i in seq_along(Curated_dfs.output)) {
  assign(Curated_dfs.output[i], Curated_test_df[[i]])
}





# Generate the function to calculate the MEAN of the log2FoldChange of each gene set resulting from the filtering above.
Mean_curated <- function(x) {
  mean_curated <- mean(x$log2FoldChange, na.rm = TRUE)
  return(mean_curated)
}

# Create an input of the data frames generated above to calculate the mean of log2FoldChange
Curated_input_dfs <- list(df_SKMEL28_VEM_MTOR_POS_EP, df_SKMEL28_VEM_MTOR_NEG_EP, df_SKMEL28_VEM_Proliferation_Science_EP, df_SKMEL28_VEM_Meta_PCNA_EP, df_SKMEL28_VEM_Proliferation_GO_EP, df_SKMEL28_VEM_Proliferation_STOP_EP, df_SKMEL28_VEM_G1_phase_EP, df_SKMEL28_VEM_Proliferation_Carey_EP, df_SKMEL28_VEM_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019_EP, df_SKMEL28_VEM_Essential_S_Phase_Viner_Breuer_Cell_Division_2019_EP, df_SKMEL28_VEM_HR_EP, df_SKMEL28_VEM_Fanconi_anemia_EP, df_SKMEL28_VEM_HF_POL_EP)

# Apply the mean function to create a list of means
Curated_mean_test_df <- lapply(Curated_input_dfs, function(x) Mean_curated(x))

# creating a vector with desired names of the files containing the mean of each gene set log2FoldChange
Curated_Mean_SKMEL28_VEM_dfs.output <- c("Mean_SKMEL28_VEM_MTOR_POS_EP", "Mean_SKMEL28_VEM_MTOR_NEG_EP", "Mean_SKMEL28_VEM_Proliferation_Science_EP", "Mean_SKMEL28_VEM_Meta_PCNA_EP", "Mean_SKMEL28_VEM_Proliferation_GO_EP", "Mean_SKMEL28_VEM_Proliferation_STOP_EP", "Mean_SKMEL28_VEM_G1_phase_EP", "Mean_SKMEL28_VEM_Proliferation_Carey_EP", "Mean_SKMEL28_VEM_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019_EP", "Mean_SKMEL28_VEM_Essential_S_Phase_Viner_Breuer_Cell_Division_2019_EP", "Mean_SKMEL28_VEM_HR_EP", "Mean_SKMEL28_VEM_Fanconi_anemia_EP", "Mean_SKMEL28_VEM_HF_POL_EP")

# loop to assign means to the desired variable name
for (i in seq_along(Curated_Mean_SKMEL28_VEM_dfs.output)) {
  assign(Curated_Mean_SKMEL28_VEM_dfs.output[i], Curated_mean_test_df[[i]])
}




# Generate the function to calculate the STANDARD DEVIATION of the log2FoldChange of each gene set resulting from the filtering above.
STDEV_curated <- function(x) {
  stdev_curated <- sd(x$log2FoldChange, na.rm = TRUE)
  return(stdev_curated)
}

# Apply the standard deviation function to create a list of standard deviations
Curated_stdev_test_df <- lapply(Curated_input_dfs, function(x) STDEV_curated(x))

# creating a vector with desired names of the files containing the standard deviation of each gene set log2FoldChange
Curated_STDEV_SKMEL28_VEM_dfs.output <- c("STDEV_SKMEL28_VEM_MTOR_POS_EP", "STDEV_SKMEL28_VEM_MTOR_NEG_EP", "STDEV_SKMEL28_VEM_Proliferation_Science_EP", "STDEV_SKMEL28_VEM_Meta_PCNA_EP", "STDEV_SKMEL28_VEM_Proliferation_GO_EP", "STDEV_SKMEL28_VEM_Proliferation_STOP_EP", "STDEV_SKMEL28_VEM_G1_phase_EP", "STDEV_SKMEL28_VEM_Proliferation_Carey_EP", "STDEV_SKMEL28_VEM_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019_EP", "STDEV_SKMEL28_VEM_Essential_S_Phase_Viner_Breuer_Cell_Division_2019_EP", "STDEV_SKMEL28_VEM_HR_EP", "STDEV_SKMEL28_VEM_Fanconi_anemia_EP", "STDEV_SKMEL28_VEM_HF_POL_EP")
# loop to assign data.frames to the desired variable name
for (i in (seq_along(Curated_STDEV_SKMEL28_VEM_dfs.output))) {
  do.call("<-",list(Curated_STDEV_SKMEL28_VEM_dfs.output[i], Curated_stdev_test_df[[i]]))
}




# Generate the data frame for the plot

# Generate the data frame for the plot
df_778_PALBO_proliferation <- data.frame(mean=c(Mean_SKMEL28_VEM_MTOR_POS_EP, Mean_SKMEL28_VEM_MTOR_NEG_EP, Mean_SKMEL28_VEM_Proliferation_Science_EP, Mean_SKMEL28_VEM_Meta_PCNA_EP, Mean_SKMEL28_VEM_Proliferation_GO_EP, Mean_SKMEL28_VEM_Proliferation_STOP_EP, Mean_SKMEL28_VEM_G1_phase_EP, Mean_SKMEL28_VEM_Proliferation_Carey_EP, Mean_SKMEL28_VEM_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019_EP, Mean_SKMEL28_VEM_Essential_S_Phase_Viner_Breuer_Cell_Division_2019_EP, Mean_SKMEL28_VEM_HR_EP, Mean_SKMEL28_VEM_Fanconi_anemia_EP, Mean_SKMEL28_VEM_HF_POL_EP),
                                         sd=c(STDEV_SKMEL28_VEM_MTOR_POS_EP, STDEV_SKMEL28_VEM_MTOR_NEG_EP, STDEV_SKMEL28_VEM_Proliferation_Science_EP, STDEV_SKMEL28_VEM_Meta_PCNA_EP, STDEV_SKMEL28_VEM_Proliferation_GO_EP, STDEV_SKMEL28_VEM_Proliferation_STOP_EP, STDEV_SKMEL28_VEM_G1_phase_EP, STDEV_SKMEL28_VEM_Proliferation_Carey_EP, STDEV_SKMEL28_VEM_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019_EP, STDEV_SKMEL28_VEM_Essential_S_Phase_Viner_Breuer_Cell_Division_2019_EP, STDEV_SKMEL28_VEM_HR_EP, STDEV_SKMEL28_VEM_Fanconi_anemia_EP, STDEV_SKMEL28_VEM_HF_POL_EP),
                                         Category=as.factor(c("SKMEL28_VEM_MTOR_POS_EP", "SKMEL28_VEM_MTOR_NEG_EP", "SKMEL28_VEM_Proliferation_Science_EP", "SKMEL28_VEM_Meta_PCNA_EP", "SKMEL28_VEM_Proliferation_GO_EP", "SKMEL28_VEM_Proliferation_STOP_EP", "SKMEL28_VEM_G1_phase_EP", "SKMEL28_VEM_Proliferation_Carey_EP", "SKMEL28_VEM_Essential_cell_cycle_EP", "SKMEL28_VEM_Essential_S_Phase_EP", "SKMEL28_VEM_HR_EP", "SKMEL28_VEM_FA_EP", "SKMEL28_VEM_HF_POL_EP")))
df_778_PALBO_proliferation




library(ggplot2)
ggplot(df_778_PALBO_proliferation, aes(x=reorder(Category, -mean), y=mean)) +
  theme_bw() +
  geom_point(position=position_dodge(width = 1), stat="identity",colour='black') +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.4) +
  theme(axis.text.x = element_text(size=12, angle=45,hjust=0.95,vjust=0.98)) +
  theme(axis.text.y = element_text(size=14)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 10))) +
  geom_point(colour = "black", size = 3) +
  xlab("") + ylab("Relative log2 fold change") +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  ggtitle("Analysis cell proliferation 778 PALBO NANOPORE 72h vs T0") +
  scale_x_discrete(limits = c("SKMEL28_VEM_Proliferation_Science_EP", "SKMEL28_VEM_Meta_PCNA_EP", "SKMEL28_VEM_Proliferation_GO_EP", "SKMEL28_VEM_Proliferation_STOP_EP", "SKMEL28_VEM_Proliferation_Carey_EP", "SKMEL28_VEM_G1_phase_EP", "SKMEL28_VEM_MTOR_POS_EP", "SKMEL28_VEM_MTOR_NEG_EP", "SKMEL28_VEM_Essential_cell_cycle_EP", "SKMEL28_VEM_Essential_S_Phase_EP", "SKMEL28_VEM_HR_EP", "SKMEL28_VEM_FA_EP", "SKMEL28_VEM_HF_POL_EP")) +
  theme(plot.title = element_text(color="black", size=14, face="bold"))


# Merge all the data frames to generate box plots
# Add a column with the specific name of the gene set
df_SKMEL28_VEM_MTOR_POS_EP$Gene_set <- "SKMEL28_VEM_MTOR_POS_EP"
df_SKMEL28_VEM_MTOR_NEG_EP$Gene_set <- "SKMEL28_VEM_MTOR_NEG_EP"
df_SKMEL28_VEM_Proliferation_Science_EP$Gene_set <- "SKMEL28_VEM_Proliferation_Science_EP"
df_SKMEL28_VEM_Meta_PCNA_EP$Gene_set <- "SKMEL28_VEM_Meta_PCNA_EP"
df_SKMEL28_VEM_Proliferation_GO_EP$Gene_set <- "SKMEL28_VEM_Proliferation_GO_EP"
df_SKMEL28_VEM_Proliferation_STOP_EP$Gene_set <- "SKMEL28_VEM_Proliferation_STOP_EP"
df_SKMEL28_VEM_G1_phase_EP$Gene_set <- "SKMEL28_VEM_G1_phase_EP"
df_SKMEL28_VEM_Proliferation_Carey_EP$Gene_set <- "SKMEL28_VEM_Proliferation_Carey_EP"
df_SKMEL28_VEM_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019_EP$Gene_set <- "SKMEL28_VEM_Essential_cell_cycle_EP"
df_SKMEL28_VEM_Essential_S_Phase_Viner_Breuer_Cell_Division_2019_EP$Gene_set <- "SKMEL28_VEM_Essential_S_Phase_EP"
df_SKMEL28_VEM_HR_EP$Gene_set <- "SKMEL28_VEM_HR_EP"
df_SKMEL28_VEM_Fanconi_anemia_EP$Gene_set <- "SKMEL28_VEM_Fanconi_anemia_EP"
df_SKMEL28_VEM_HF_POL_EP$Gene_set <- "SKMEL28_VEM_HF_POL_EP"



# Combine the data frames
combined_df <- rbind(
  df_SKMEL28_VEM_MTOR_POS_EP,
  df_SKMEL28_VEM_MTOR_NEG_EP,
  df_SKMEL28_VEM_Proliferation_Science_EP,
  df_SKMEL28_VEM_Meta_PCNA_EP,
  df_SKMEL28_VEM_Proliferation_GO_EP,
  df_SKMEL28_VEM_Proliferation_STOP_EP,
  df_SKMEL28_VEM_G1_phase_EP,
  df_SKMEL28_VEM_Proliferation_Carey_EP,
  df_SKMEL28_VEM_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019_EP,
  df_SKMEL28_VEM_Essential_S_Phase_Viner_Breuer_Cell_Division_2019_EP,
  df_SKMEL28_VEM_HR_EP,
  df_SKMEL28_VEM_Fanconi_anemia_EP,
  df_SKMEL28_VEM_HF_POL_EP
)

# Ensure the Gene_set column is a factor with the correct order
combined_df$Gene_set <- factor(combined_df$Gene_set, levels = c("SKMEL28_VEM_Proliferation_Science_EP", "SKMEL28_VEM_Meta_PCNA_EP", "SKMEL28_VEM_Proliferation_GO_EP", "SKMEL28_VEM_Proliferation_STOP_EP", "SKMEL28_VEM_Proliferation_Carey_EP", "SKMEL28_VEM_G1_phase_EP", "SKMEL28_VEM_MTOR_POS_EP", "SKMEL28_VEM_MTOR_NEG_EP", "SKMEL28_VEM_Essential_cell_cycle_EP", "SKMEL28_VEM_Essential_S_Phase_EP", "SKMEL28_VEM_HR_EP", "SKMEL28_VEM_Fanconi_anemia_EP", "SKMEL28_VEM_HF_POL_EP"))


# Generate the box plot
library(ggplot2)
ggplot(combined_df, aes(x=Gene_set, y=log2FoldChange)) +
  geom_boxplot(fill="lightgrey") + 
  stat_summary(fun.y=mean, geom="point", shape=20, size=4, color="blue") +
  # geom_jitter(width=0.2, alpha=0.5) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle=45, hjust=0.95, vjust=0.98)) +
  theme(axis.text.y = element_text(size=14)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 10))) +
  xlab("") + ylab("Relative log2 fold change") +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  ggtitle("Analysis cell proliferation 778 palbociclib early phase") +
  scale_x_discrete(limits = levels(combined_df$Gene_set)) +
  theme(plot.title = element_text(color="black", size=14, face="bold")) 






Mean_df_MTOR_POS_SKMEL28_VEM_EP <- mean(df_SKMEL28_VEM_MTOR_POS_EP$log2FoldChange)
Mean_df_MTOR_NEG_SKMEL28_VEM_EP <- mean(df_SKMEL28_VEM_MTOR_NEG_EP$log2FoldChange)
Mean_df_Proliferation_Science_SKMEL28_VEM_EP <- mean(df_SKMEL28_VEM_Proliferation_Science_EP$log2FoldChange)
Mean_df_Meta_PCNA_SKMEL28_VEM_EP <- mean(df_SKMEL28_VEM_Meta_PCNA_EP$log2FoldChange)
Mean_df_Proliferation_GO_SKMEL28_VEM_EP <- mean(df_SKMEL28_VEM_Proliferation_GO_EP$log2FoldChange)
Mean_df_Proliferation_STOP_SKMEL28_VEM_EP <- mean(df_SKMEL28_VEM_Proliferation_STOP_EP$log2FoldChange)
Mean_df_G1_phase_SKMEL28_VEM_EP <- mean(df_SKMEL28_VEM_G1_phase_EP$log2FoldChange)
Mean_df_Proliferation_Carey_SKMEL28_VEM_EP <- mean(df_SKMEL28_VEM_Proliferation_Carey_EP$log2FoldChange)
Mean_df_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019_SKMEL28_VEM_EP <- mean(df_SKMEL28_VEM_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019_EP$log2FoldChange)
Mean_df_Essential_S_Phase_Viner_Breuer_Cell_Division_2019_SKMEL28_VEM_EP <- mean(df_SKMEL28_VEM_Essential_S_Phase_Viner_Breuer_Cell_Division_2019_EP$log2FoldChange)
Mean_df_HR_SKMEL28_VEM_EP <- mean(df_SKMEL28_VEM_HR_EP$log2FoldChange)
Mean_df_FA_SKMEL28_VEM_EP <- mean(df_SKMEL28_VEM_Fanconi_anemia_EP$log2FoldChange)
Mean_df_HF_POL_SKMEL28_VEM_EP <- mean(df_SKMEL28_VEM_HF_POL_EP$log2FoldChange)

# Perform the Wilcoxon signed-rank test
wilcox_test_MTOR_POS_SKMEL28_VEM_EP <- wilcox.test(df_SKMEL28_VEM_MTOR_POS_EP$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_MTOR_NEG_SKMEL28_VEM_EP <- wilcox.test(df_SKMEL28_VEM_MTOR_NEG_EP$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_Proliferation_Science_SKMEL28_VEM_EP <- wilcox.test(df_SKMEL28_VEM_Proliferation_Science_EP$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_Meta_PCNA_SKMEL28_VEM_EP <- wilcox.test(df_SKMEL28_VEM_Meta_PCNA_EP$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_Proliferation_GO_SKMEL28_VEM_EP <- wilcox.test(df_SKMEL28_VEM_Proliferation_GO_EP$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_Proliferation_STOP_SKMEL28_VEM_EP <- wilcox.test(df_SKMEL28_VEM_Proliferation_STOP_EP$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_G1_phase_SKMEL28_VEM_EP <- wilcox.test(df_SKMEL28_VEM_G1_phase_EP$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_Proliferation_Carey_SKMEL28_VEM_EP <- wilcox.test(df_SKMEL28_VEM_Proliferation_Carey_EP$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019_SKMEL28_VEM_EP <- wilcox.test(df_SKMEL28_VEM_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019_EP$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_Essential_S_Phase_Viner_Breuer_Cell_Division_2019_SKMEL28_VEM_EP <- wilcox.test(df_SKMEL28_VEM_Essential_S_Phase_Viner_Breuer_Cell_Division_2019_EP$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_HR_SKMEL28_VEM_EP <- wilcox.test(df_SKMEL28_VEM_HR_EP$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_FA_SKMEL28_VEM_EP <- wilcox.test(df_SKMEL28_VEM_Fanconi_anemia_EP$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_HF_POL_SKMEL28_VEM_EP <- wilcox.test(df_SKMEL28_VEM_HF_POL_EP$log2FoldChange, mu = 0, alternative = "less")
















