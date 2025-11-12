# This script generates the analysis of the immediate response of model systems sequenced at 72h with Nanopore

# Immediate 72h Nanopore Response Analysis
# Aim:
#   Integrate transcriptome-wide thermodynamic features (Normalized MFE and
#   heat capacity) with Nanopore DEXSeq/DESeq2 data to characterize the
#   immediate transcriptional response to drug treatment (72h).
#
# Workflow:
#   1) Load annotated thermodynamic table (MFE, Specific_heat, transcript length).
#   2) Join with transcript-level DEXSeq results (treated vs control, padj < 0.05).
#   3) Separate UP vs DOWN regulated transcripts; compare Normalized_MFE
#      distributions (t.test, KS test), visualize (scatter, violin, density).
#   4) Define high-MFE (>75th percentile) and low-MFE (<25th percentile).
#      - Test enrichment among UP vs DOWN (Fisher, z-test).
#      - Compute fractions of transcript counts (control vs drug).
#   5) Evaluate heat capacity:
#      - All DE transcripts (UP vs DOWN).
#      - High-MFE subsets only (UP vs DOWN).
#   6) Pathway-level validation (DESeq2 gene-level data):
#      - Proliferation and DNA repair sets; compute mean ± SD log2FC.
#      - Visualize via dot/error plots and box plots; Wilcoxon tests vs 0.


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


# Load Table_ANNOTATED_MFE_heat_capacity
Gencode_TX_norm_MFE_HC_table <- read.delim("path_to/Gencode_TX_norm_MFE_HC_table.csv",
                                           sep = ",", header = T, stringsAsFactors = F)
Table_ANNOTATED <- Gencode_TX_norm_MFE_HC_table








# 778 palbociclib ANNOTATED transcripts ####
X778_PALBO_dexseq_ANNOTATED<-read.delim(file="path_to/72h_778_palbociclib.tx_counts.both_conditions.treated_untreated_dexseq.tab", header = TRUE, sep = "\t", quote = "") %>%
  dplyr::rename(ensembl_transcript_id_version=featureID) %>%
  dplyr::inner_join(Table_ANNOTATED, by="ensembl_transcript_id_version") %>%
  dplyr::filter(padj < 0.05) # The file 'Table' contains the Normalized_MFE values of all ANNOTATED transcripts


# Generate a scatter plot Log2_FC vs Normalized_MFE
library(ggplot2)
ggplot(X778_PALBO_dexseq_ANNOTATED, aes(y=Normalized_MFE, x=log2fold_treated_untreated)) +
  geom_point(size=1.5, alpha=0.6)+
  theme_classic() +
  theme(legend.position="none") +
  geom_vline(xintercept = 0, size=0.3) +
  #xlim(-5.5,5.5) +
  theme(axis.title = element_text(size = 20)) +
  theme(plot.title = element_text(size = 15)) +
  labs(y="Normalized_MFE",x="Log2_FC_expression",title="X778_PALBO_shared_ANNOTATED_transcripts_72h_Log2_FC_adj.P<0.05")




# Show the data with a violin plot
ggplot(X778_PALBO_dexseq_ANNOTATED, aes(x = factor(log2fold_treated_untreated > 0), y = Normalized_MFE)) +
  geom_violin(fill = "lightgrey", alpha = 0.7) +
  geom_jitter(width = 0.1, size = 1, alpha = 0.3) +
  theme_classic() +
  theme(legend.position = "none", axis.title = element_text(size = 20), plot.title = element_text(size = 10),
        axis.text.x = element_text(size = 16, color = "black")) +
  geom_vline(xintercept = 0, size = 0.3) +
  labs(y = "Normalized_MFE", x = NULL, title = "778_shared_ANNOTATED_transcripts_72h_Log2_FC_adj.P<0.05") +
  scale_x_discrete(labels = c("Log2_FC < 0", "Log2_FC > 0"))



X778_PALBO_dexseq_ANNOTATED_UP <- dplyr::filter(X778_PALBO_dexseq_ANNOTATED, log2fold_treated_untreated > 0)
Mean_X778_PALBO_dexseq_ANNOTATED_UP <- mean(X778_PALBO_dexseq_ANNOTATED_UP$Normalized_MFE)

X778_PALBO_dexseq_ANNOTATED_DOWN <- dplyr::filter(X778_PALBO_dexseq_ANNOTATED, log2fold_treated_untreated < 0)
Mean_X778_PALBO_dexseq_ANNOTATED_DOWN <- mean(X778_PALBO_dexseq_ANNOTATED_DOWN$Normalized_MFE)


# Wilcoxon rank-sum test (Mann–Whitney U)
# Split data into UP and DOWN groups
X778_PALBO_dexseq_ANNOTATED_UP <- dplyr::filter(X778_PALBO_dexseq_ANNOTATED, log2fold_treated_untreated > 0)
X778_PALBO_dexseq_ANNOTATED_DOWN <- dplyr::filter(X778_PALBO_dexseq_ANNOTATED, log2fold_treated_untreated < 0)

# Wilcoxon rank-sum test
wilcox_test_MFE_UP_vs_DOWN <- wilcox.test(
  X778_PALBO_dexseq_ANNOTATED_UP$Normalized_MFE,
  X778_PALBO_dexseq_ANNOTATED_DOWN$Normalized_MFE,
  alternative = "greater" # tests if median(UP) > median(DOWN)
)

print(wilcox_test_MFE_UP_vs_DOWN)










# Show the data in a density plot
max_length_778_MFE <- max(c(length(X778_PALBO_dexseq_ANNOTATED_UP$Normalized_MFE), length(X778_PALBO_dexseq_ANNOTATED_DOWN$Normalized_MFE)))

data_778_MFE_density <- data.frame(MFE_TX_UP_778 = c(X778_PALBO_dexseq_ANNOTATED_UP$Normalized_MFE, rep(NA, max_length_778_MFE - length(X778_PALBO_dexseq_ANNOTATED_UP$Normalized_MFE))), # This is useful when we need to create a data frame and the number of values between columns are different
                                   MFE_TX_DOWN_778 = c(X778_PALBO_dexseq_ANNOTATED_DOWN$Normalized_MFE, rep(NA, max_length_778_MFE- length(X778_PALBO_dexseq_ANNOTATED_DOWN$Normalized_MFE))))
library(reshape)
data_778_MFE_density <- melt(data_778_MFE_density)

library(ggplot2)
ggplot(data= data_778_MFE_density, aes(x=value, group=variable, color=variable)) +
  geom_density (size=1) +
  theme_classic() +
  #xlim(0,1.5) +
  labs(title = "MFE_distribution_ANNOTATED_Transcripts_778_72h")+
  labs(y="Density")+
  labs(x="DGfold / nt (kcal/mol)")













## Analysis net MFE and heat capacity ####
# Determine whether within the differentially expressed genes there is a significant up-regulation of high MFE transcripts:
# Define high MFE transcripts by taking the 75th percentile of the combined MFE distributions of control and drug-treated samples.
# Count the number of transcripts in the significantly up- and down-regulated genes that have MFE values above this threshold.
# Perform a Fisher’s exact test to compare the proportions of high MFE transcripts between the two groups.

X778_PALBO_dexseq_ANNOTATED_norm_TX_counts_MFE_df <- X778_PALBO_dexseq_ANNOTATED %>%
  dplyr::select(ensembl_transcript_id_version, ensembl_transcript_id, ensembl_gene_id, external_gene_name, transcript_is_canonical, Transcript_length,
                countData.72h_778_untreated_R1, countData.72h_778_untreated_R2, countData.72h_778_untreated_R3,
                countData.72h_778_palbociclib_R1, countData.72h_778_palbociclib_R2, countData.72h_778_palbociclib_R3,
                Normalized_MFE, Cp_37C)


# Identify the 75th percentile of the combined MFE distributions of control and drug-treated samples
quantile_075_X778_PALBO_dexseq_ANNOTATED_norm_TX_counts_MFE_df <- quantile(X778_PALBO_dexseq_ANNOTATED_norm_TX_counts_MFE_df$Normalized_MFE, probs = c(0.75))

# Select all the transcripts with Normalized_MFE values above the 75th percentile
X778_PALBO_dexseq_ANNOTATED_high_MFE_TX <- dplyr::filter(X778_PALBO_dexseq_ANNOTATED_norm_TX_counts_MFE_df, Normalized_MFE > quantile_075_X778_PALBO_dexseq_ANNOTATED_norm_TX_counts_MFE_df)

# Count the number of transcripts in the significantly up- and down-regulated pools that have MFE values above this threshold.
X778_PALBO_dexseq_ANNOTATED_high_MFE_TX_in_UP <- dplyr::filter(X778_PALBO_dexseq_ANNOTATED_UP, ensembl_transcript_id_version %in% X778_PALBO_dexseq_ANNOTATED_high_MFE_TX$ensembl_transcript_id_version)
X778_PALBO_dexseq_ANNOTATED_high_MFE_TX_in_DOWN <- dplyr::filter(X778_PALBO_dexseq_ANNOTATED_DOWN, ensembl_transcript_id_version %in% X778_PALBO_dexseq_ANNOTATED_high_MFE_TX$ensembl_transcript_id_version)


# Perform a Fisher’s exact test to compare the proportions of high MFE transcripts between the two groups.
X778_PALBO_dexseq_ANNOTATED_total_UP <- nrow(X778_PALBO_dexseq_ANNOTATED_UP)
X778_PALBO_dexseq_ANNOTATED_total_DOWN <- nrow(X778_PALBO_dexseq_ANNOTATED_DOWN)
X778_PALBO_dexseq_ANNOTATED_count_high_MFE_UP <- nrow(X778_PALBO_dexseq_ANNOTATED_high_MFE_TX_in_UP)
X778_PALBO_dexseq_ANNOTATED_count_high_MFE_DOWN <- nrow(X778_PALBO_dexseq_ANNOTATED_high_MFE_TX_in_DOWN)

X778_PALBO_dexseq_ANNOTATED_contingency_table <- matrix(c(X778_PALBO_dexseq_ANNOTATED_count_high_MFE_UP, X778_PALBO_dexseq_ANNOTATED_total_UP - X778_PALBO_dexseq_ANNOTATED_count_high_MFE_UP,
                                        X778_PALBO_dexseq_ANNOTATED_count_high_MFE_DOWN, X778_PALBO_dexseq_ANNOTATED_total_DOWN - X778_PALBO_dexseq_ANNOTATED_count_high_MFE_DOWN),
                                      nrow = 2, byrow = TRUE,
                                      dimnames = list(UP = c("UP-regulated", "DOWN-regulated"),
                                                      DOWN = c("High MFE TX", "Not High MFE TX")))

X778_PALBO_dexseq_ANNOTATED_contingency_table_df <- as.data.frame(X778_PALBO_dexseq_ANNOTATED_contingency_table)

X778_PALBO_dexseq_ANNOTATED_fisher_result <- fisher.test(X778_PALBO_dexseq_ANNOTATED_contingency_table)

cat("X778_PALBO_dexseq_ANNOTATED P-value from Fisher test:", X778_PALBO_dexseq_ANNOTATED_fisher_result$p.value, "\n")


# Look at the fractions
X778_PALBO_dexseq_ANNOTATED_high_MFE_in_UP_fraction <- X778_PALBO_dexseq_ANNOTATED_count_high_MFE_UP / X778_PALBO_dexseq_ANNOTATED_total_UP
X778_PALBO_dexseq_ANNOTATED_high_MFE_in_DOWN_fraction <- X778_PALBO_dexseq_ANNOTATED_count_high_MFE_DOWN / X778_PALBO_dexseq_ANNOTATED_total_DOWN

# Output the results
cat("Fraction of X778_PALBO_dexseq_ANNOTATED high MFE TX in UP-regulated genes:", X778_PALBO_dexseq_ANNOTATED_high_MFE_in_UP_fraction, "\n") # This value goes in the summary table
cat("Fraction of X778_PALBO_dexseq_ANNOTATED high MFE TX in DOWN-regulated genes:", X778_PALBO_dexseq_ANNOTATED_high_MFE_in_DOWN_fraction, "\n") # This value goes in the summary table
cat("Mean of X778_PALBO_dexseq_ANNOTATED high MFE TX in UP-regulated genes:", mean(X778_PALBO_dexseq_ANNOTATED_high_MFE_TX_in_UP$Normalized_MFE), "\n")
cat("Mean of X778_PALBO_dexseq_ANNOTATED high MFE TX in DOWN-regulated genes:", mean(X778_PALBO_dexseq_ANNOTATED_high_MFE_TX_in_DOWN$Normalized_MFE), "\n")
cat("Mode of X778_PALBO_dexseq_ANNOTATED high MFE TX in UP-regulated genes:", calculate_mode(X778_PALBO_dexseq_ANNOTATED_high_MFE_TX_in_UP$Normalized_MFE), "\n")
cat("Mode of X778_PALBO_dexseq_ANNOTATED high MFE TX in DOWN-regulated genes:", calculate_mode(X778_PALBO_dexseq_ANNOTATED_high_MFE_TX_in_DOWN$Normalized_MFE), "\n")


# Perform a two-proportion z-test to compare proportions
# Calculate the proportions and total counts
X778_PALBO_dexseq_ANNOTATED_high_MFE_in_UP_proportion <- X778_PALBO_dexseq_ANNOTATED_high_MFE_in_UP_fraction
X778_PALBO_dexseq_ANNOTATED_high_MFE_in_DOWN_proportion <- X778_PALBO_dexseq_ANNOTATED_high_MFE_in_DOWN_fraction

# Combined proportion for pooled standard error
combined_proportion <- (X778_PALBO_dexseq_ANNOTATED_count_high_MFE_UP + X778_PALBO_dexseq_ANNOTATED_count_high_MFE_DOWN) /
  (X778_PALBO_dexseq_ANNOTATED_total_UP + X778_PALBO_dexseq_ANNOTATED_total_DOWN)

# Standard error of the difference in proportions
standard_error <- sqrt(
  combined_proportion * (1 - combined_proportion) * 
    (1 / X778_PALBO_dexseq_ANNOTATED_total_UP + 1 / X778_PALBO_dexseq_ANNOTATED_total_DOWN)
)

# Z-statistic
z_statistic <- (X778_PALBO_dexseq_ANNOTATED_high_MFE_in_UP_proportion - X778_PALBO_dexseq_ANNOTATED_high_MFE_in_DOWN_proportion) /
  standard_error

# Two-tailed p-value
X778_PALBO_dexseq_ANNOTATED_z_test_MFE_p_value <- 2 * (1 - pnorm(abs(z_statistic)))








# Show data as violin plot
# Add a new column to indicate if the value is above the 75th percentile
X778_PALBO_dexseq_ANNOTATED$above_quantile_075 <- X778_PALBO_dexseq_ANNOTATED$Normalized_MFE > quantile_075_X778_PALBO_dexseq_ANNOTATED_norm_TX_counts_MFE_df

# Generate violin plot with conditional jitter coloring
ggplot(X778_PALBO_dexseq_ANNOTATED, aes(x = factor(log2fold_treated_untreated > 0), y = Normalized_MFE)) +
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
    title = "X778_PALBO_dexseq_ANNOTATED_Log2_FC_adj.P<0.05"
  ) +
  scale_x_discrete(labels = c("log2FoldChange < 0", "log2FoldChange > 0"))













## Fraction of total counts analysis ####
# Determine the impact of these changes in terms of fraction of total counts:
# Calculate the fraction of total counts of high MFE transcripts in control and drug-treated groups. 
# This is done by summing the counts of the high MFE transcripts and dividing by the total counts in each group.
# Fisher’s exact test to compare the proportions of high MFE transcripts counts between the two groups.


# Load the data with the counts of the entire transcriptome
X778_PALBO_dexseq_ANNOTATED_whole_transcriptome <- read.delim(file="path_to/72h_778_palbociclib.tx_counts.both_conditions.treated_untreated_dexseq.tab", header = TRUE, sep = "\t", quote = "") %>%
  dplyr::rename(ensembl_transcript_id_version=featureID) %>%
  dplyr::inner_join(Table_ANNOTATED, by="ensembl_transcript_id_version") %>%
  dplyr::select(ensembl_transcript_id_version, ensembl_transcript_id, ensembl_gene_id, external_gene_name, transcript_is_canonical, Transcript_length,
                countData.72h_778_untreated_R1, countData.72h_778_untreated_R2, countData.72h_778_untreated_R3,
                countData.72h_778_palbociclib_R1, countData.72h_778_palbociclib_R2, countData.72h_778_palbociclib_R3,
                Normalized_MFE, Cp_37C)

# Convert from wide to long format (one row per transcript x replicate)
X778_PALBO_dexseq_ANNOTATED_long_data <- X778_PALBO_dexseq_ANNOTATED_whole_transcriptome %>%
  tidyr::pivot_longer(
    cols = starts_with("countData.72h"),
    names_to = "replicate",
    values_to = "counts"
  )

# Add a condition column based on the replicate
X778_PALBO_dexseq_ANNOTATED_long_data <- X778_PALBO_dexseq_ANNOTATED_long_data %>%
  dplyr::mutate(
    condition = dplyr::if_else(
      replicate %in% c(
        "countData.72h_778_untreated_R1",
        "countData.72h_778_untreated_R2",
        "countData.72h_778_untreated_R3"
      ),
      "Control", "palbociclib"
    ),
    # (optional) set factors for tidy plotting
    condition = factor(condition, levels = c("Control", "palbociclib"))
  )

X778_PALBO_dexseq_ANNOTATED_norm_TX_aggregated_data_avg_counts <- X778_PALBO_dexseq_ANNOTATED_long_data %>%
  dplyr::group_by(ensembl_transcript_id_version, condition, Normalized_MFE) %>%
  dplyr::summarise(
    avg_counts = mean(counts, na.rm = TRUE),
    .groups = "drop"
  )





# Use Fisher
# Test whether cells exposed to palbociclib have a higher fraction of transcripts in the 75th percentile of MFE values compared to control cells.

# 75th percentile threshold
mfe_threshold <- quantile(
  X778_PALBO_dexseq_ANNOTATED_norm_TX_aggregated_data_avg_counts$Normalized_MFE,
  0.75, na.rm = TRUE
)

# Label high MFE
df <- X778_PALBO_dexseq_ANNOTATED_norm_TX_aggregated_data_avg_counts %>%
  dplyr::mutate(MFE_class = ifelse(Normalized_MFE >= mfe_threshold, "High_MFE", "Other"))

# Aggregate counts per condition x class
counts_table <- df %>%
  dplyr::group_by(condition, MFE_class) %>%
  dplyr::summarise(total_counts = sum(avg_counts, na.rm = TRUE), .groups = "drop")

# 2x2 contingency (row1=palbociclib, row2=Control; col1=High_MFE, col2=Other)
contingency <- counts_table %>%
  tidyr::pivot_wider(names_from = MFE_class, values_from = total_counts, values_fill = 0) %>%
  dplyr::arrange(factor(condition, levels = c("palbociclib", "Control")))

fisher_matrix <- as.matrix(contingency[, c("High_MFE", "Other")])
rownames(fisher_matrix) <- contingency$condition

# One-sided Fisher: palbociclib > Control
Fisher_test_778_palbo_75th <- fisher.test(fisher_matrix, alternative = "greater")
print(Fisher_test_778_palbo_75th)

# Report High_MFE fractions per condition
Counts_table_778_palbo_75th <- counts_table %>%
  dplyr::group_by(condition) %>%
  dplyr::mutate(fraction = total_counts / sum(total_counts)) %>%
  dplyr::filter(MFE_class == "High_MFE")
print(Counts_table_778_palbo_75th)


# High-MFE fraction among total counts in palbociclib
X778_PALBO_dexseq_ANNOTATED_palbociclib_high_MFE_TX_counts_fraction <-
  Counts_table_778_palbo_75th %>%
  dplyr::filter(condition == "palbociclib") %>%
  dplyr::pull(fraction) %>%
  as.numeric()


X778_PALBO_dexseq_ANNOTATED_control_high_MFE_TX_counts_fraction <-
  Counts_table_778_palbo_75th %>%
  dplyr::filter(condition == "Control") %>%
  dplyr::pull(fraction) %>%
  as.numeric()






# Test whether the Control cells have a higher fraction of transcripts in the 25th percentile of MFE values compared to palbociclib-treated cells.

# 25th percentile threshold
mfe_threshold_25 <- quantile(
  X778_PALBO_dexseq_ANNOTATED_norm_TX_aggregated_data_avg_counts$Normalized_MFE,
  0.25, na.rm = TRUE
)

# Label low MFE (<= 25th percentile)
df25 <- X778_PALBO_dexseq_ANNOTATED_norm_TX_aggregated_data_avg_counts %>%
  dplyr::mutate(MFE_class = ifelse(Normalized_MFE <= mfe_threshold_25, "Low_MFE", "Other"))

# Aggregate counts per condition x class
counts_table_25 <- df25 %>%
  dplyr::group_by(condition, MFE_class) %>%
  dplyr::summarise(total_counts = sum(avg_counts, na.rm = TRUE), .groups = "drop")

# 2x2 contingency (row1=Control, row2=palbociclib; col1=Low_MFE, col2=Other)
contingency_25 <- counts_table_25 %>%
  tidyr::pivot_wider(names_from = MFE_class, values_from = total_counts, values_fill = 0) %>%
  dplyr::arrange(factor(condition, levels = c("Control", "palbociclib")))

fisher_matrix_25 <- as.matrix(contingency_25[, c("Low_MFE", "Other")])
rownames(fisher_matrix_25) <- contingency_25$condition

# One-sided Fisher: tests whether Control has a higher Low_MFE fraction than palbociclib
Fisher_test_778_palbo_25th <- fisher.test(fisher_matrix_25, alternative = "greater")
print(Fisher_test_778_palbo_25th)

# Report Low_MFE fractions per condition
Counts_table_778_palbo_25th <- counts_table_25 %>%
  dplyr::group_by(condition) %>%
  dplyr::mutate(fraction = total_counts / sum(total_counts)) %>%
  dplyr::filter(MFE_class == "Low_MFE")
print(Counts_table_778_palbo_25th)




















## X778_PALBO_dexseq_ANNOTATED heat capacity changes ####

# This is calculated on the total n. of transcripts with a significant log2 fold change
# Generate a scatter plot Log2_FC vs heat capacity
ggplot(X778_PALBO_dexseq_ANNOTATED, aes(y=Cp_37C, x=log2fold_treated_untreated)) +
  geom_point(size=1.5, alpha=0.1)+
  theme_classic() +
  theme(legend.position="none") +
  geom_vline(xintercept = 0, size=0.3) +
  xlim(-12, 12) +
  ylim(0,400) +
  theme(axis.title = element_text(size = 20)) +
  theme(plot.title = element_text(size = 15)) +
  labs(y="Heat capacity (J/C)",x="Log2_FC_expression",title="X778_PALBO_dexseq_ANNOTATED_Log2_FC_adj.P<0.05")


# Show the data with a box plot
ggplot(X778_PALBO_dexseq_ANNOTATED, aes(x = factor(log2fold_treated_untreated > 0), y = Cp_37C)) +
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
  labs(y="Heat capacity (J/C)",x="",title="X778_PALBO_dexseq_ANNOTATED_Log2_FC_adj.P<0.05 heat_capacity all transcripts") +
  scale_x_discrete(labels = c("log2FoldChange < 0", "log2FoldChange > 0"))

Mean_heat_capacity_X778_PALBO_dexseq_ANNOTATED_UP <- mean(X778_PALBO_dexseq_ANNOTATED_UP$Cp_37C, na.rm = TRUE)
Mean_heat_capacity_X778_PALBO_dexseq_ANNOTATED_DOWN <- mean(X778_PALBO_dexseq_ANNOTATED_DOWN$Cp_37C, na.rm = TRUE)
TTEST_heat_capacity_X778_PALBO_dexseq_ANNOTATED_all_UP_vs_DOWN <- t.test(X778_PALBO_dexseq_ANNOTATED_DOWN$Cp_37C, X778_PALBO_dexseq_ANNOTATED_UP$Cp_37C, paired=FALSE)
print(TTEST_heat_capacity_X778_PALBO_dexseq_ANNOTATED_all_UP_vs_DOWN)
# Add K-S test to look at the differences between the up-regukated and down-regulated genes
ks_test_X778_PALBO_dexseq_ANNOTATED_all_UP_vs_DOWN_heat_capacity <- ks.test(X778_PALBO_dexseq_ANNOTATED_DOWN$Cp_37C, X778_PALBO_dexseq_ANNOTATED_UP$Cp_37C)
print(ks_test_X778_PALBO_dexseq_ANNOTATED_all_UP_vs_DOWN_heat_capacity)




# Calculate the changes of heat capacity only in the high MFE transcripts

# Contribution of the UP or DOWN-regulated high MFE transcripts to the heat capacity
# Look at the heat capacity of the UP-regulated high MFE transcripts
X778_PALBO_dexseq_ANNOTATED_high_MFE_TX_in_UP_heat_capacity <- mean(X778_PALBO_dexseq_ANNOTATED_high_MFE_TX_in_UP$Cp_37C)

# Look at the heat capacity of the DOWN-regulated high MFE transcripts
X778_PALBO_dexseq_ANNOTATED_high_MFE_TX_in_DOWN_heat_capacity <- mean(X778_PALBO_dexseq_ANNOTATED_high_MFE_TX_in_DOWN$Cp_37C)

# Heat capacity gain in UP-regulated high MFE transcripts vs DOWN-regulated high MFE transcripts
X778_PALBO_dexseq_ANNOTATED_heat_capacity_gain_UP_vs_DOWN <- X778_PALBO_dexseq_ANNOTATED_high_MFE_TX_in_UP_heat_capacity - X778_PALBO_dexseq_ANNOTATED_high_MFE_TX_in_DOWN_heat_capacity
cat("X778_PALBO_dexseq_ANNOTATED net heat capacity UP vs DOWN high MFE TX:", X778_PALBO_dexseq_ANNOTATED_heat_capacity_gain_UP_vs_DOWN, "\n")


# Extract specific heat data
X778_PALBO_dexseq_ANNOTATED_Cp_37C_UP <- X778_PALBO_dexseq_ANNOTATED_high_MFE_TX_in_UP$Cp_37C
X778_PALBO_dexseq_ANNOTATED_Cp_37C_DOWN <- X778_PALBO_dexseq_ANNOTATED_high_MFE_TX_in_DOWN$Cp_37C

# Combine data into a single data frame for easier plotting
X778_PALBO_dexseq_ANNOTATED_Cp_37C_data <- data.frame(
  X778_PALBO_dexseq_ANNOTATED_Cp_37C = c(X778_PALBO_dexseq_ANNOTATED_Cp_37C_UP, X778_PALBO_dexseq_ANNOTATED_Cp_37C_DOWN),
  regulation = factor(c(rep("UP", length(X778_PALBO_dexseq_ANNOTATED_Cp_37C_UP)), rep("DOWN", length(X778_PALBO_dexseq_ANNOTATED_Cp_37C_DOWN))))
)

# Generate box plot
ggplot(X778_PALBO_dexseq_ANNOTATED_Cp_37C_data, aes(x = regulation, y = X778_PALBO_dexseq_ANNOTATED_Cp_37C, fill = regulation)) +
  geom_boxplot(fill = "lightgrey", alpha = 0.7) +
  labs(title = "X778_PALBO_dexseq_ANNOTATED specific Heat of High MFE Transcripts",
       x = "",
       y = "Specific Heat (J/C)") +
  scale_x_discrete(labels = c("log2_FC < 0", "log2_FC > 0")) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 16),
        axis.title = element_text(size = 25))

Mean_heat_capacity_X778_PALBO_dexseq_ANNOTATED_high_MFE_TX_UP <- mean(X778_PALBO_dexseq_ANNOTATED_high_MFE_TX_in_UP$Cp_37C)
Mean_heat_capacity_X778_PALBO_dexseq_ANNOTATED_high_MFE_TX_DOWN <- mean(X778_PALBO_dexseq_ANNOTATED_high_MFE_TX_in_DOWN$Cp_37C)
TTEST_heat_capacity_X778_PALBO_dexseq_ANNOTATED_all_high_MFE_TX_UP_vs_DOWN <- t.test(X778_PALBO_dexseq_ANNOTATED_high_MFE_TX_in_DOWN$Cp_37C, X778_PALBO_dexseq_ANNOTATED_high_MFE_TX_in_UP$Cp_37C, paired=FALSE)
print(TTEST_heat_capacity_X778_PALBO_dexseq_ANNOTATED_all_high_MFE_TX_UP_vs_DOWN)

# Add K-S test to look at the differences between the up-regukated and down-regulated genes
ks_test_heat_capacity_X778_PALBO_dexseq_ANNOTATED_all_high_MFE_TX_UP_vs_DOWN <- ks.test(X778_PALBO_dexseq_ANNOTATED_high_MFE_TX_in_DOWN$Cp_37C, X778_PALBO_dexseq_ANNOTATED_high_MFE_TX_in_UP$Cp_37C)
print(ks_test_heat_capacity_X778_PALBO_dexseq_ANNOTATED_all_high_MFE_TX_UP_vs_DOWN)










# Analysis cell proliferation Nanopore data 778 PALBO ####
source("path_to/All_gene_sets_AC.R")

# Open the DESeq2 file with the log2 fold changes at the gene level
X778_PALBO_DESeq2 <- read.delim(file="path_to/72h_778_palbociclib.tx_counts.both_conditions.treated_untreated_deseq2.tab", header = TRUE, sep = "\t", quote = "")

filtering <- function(x) {
  Filtering <- dplyr::filter(X778_PALBO_DESeq2, gene_name %in% x)
  return(Filtering)
}

# Generate the list containing all the gene sets to be used for the filtering.
CURATED_PATHWAYS.input <- list(MTOR_positive_regulators, MTOR_negative_regulators, Proliferation_Science_paper, Meta_PCNA, Proliferation_GO_genes_Magill_Sack_2018, Proliferation_STOP_genes_Magill_Sack_2018, Reactome_G1_phase, Proliferation_Carey_2021, Essential_cell_cycle_Viner_Breuer_Cell_Division_2019, Essential_S_Phase_Viner_Breuer_Cell_Division_2019, Homologous_recombination, Fanconi_anemia, High_fidelity_polymerases)

# This creates the names that will be used for the final data frames of each file resulting from the filtering of the target data frame with each gene set
CURATED_PATHWAYS.output <- c("X778_PALBO_MTOR_POS", "X778_PALBO_MTOR_NEG", "X778_PALBO_Proliferation_Science", "X778_PALBO_Meta_PCNA", "X778_PALBO_Proliferation_GO", "X778_PALBO_Proliferation_STOP", "X778_PALBO_G1_phase", "X778_PALBO_Proliferation_Carey", "X778_PALBO_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019", "X778_PALBO_Essential_S_Phase_Viner_Breuer_Cell_Division_2019", "X778_PALBO_HR", "X778_PALBO_Fanconi_anemia", "X778_PALBO_HF_POL")

# Assign the names of the gene sets to each entry in the list
names(CURATED_PATHWAYS.input) <- CURATED_PATHWAYS.output

# Apply the filtering function to create a list of data frames
Curated_test_df <- lapply(CURATED_PATHWAYS.input, function(x) filtering(x))

# creating a vector with desired names of filtered data frames
Curated_dfs.output <- c("df_X778_PALBO_MTOR_POS", "df_X778_PALBO_MTOR_NEG", "df_X778_PALBO_Proliferation_Science", "df_X778_PALBO_Meta_PCNA", "df_X778_PALBO_Proliferation_GO", "df_X778_PALBO_Proliferation_STOP", "df_X778_PALBO_G1_phase", "df_X778_PALBO_Proliferation_Carey", "df_X778_PALBO_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019", "df_X778_PALBO_Essential_S_Phase_Viner_Breuer_Cell_Division_2019", "df_X778_PALBO_HR", "df_X778_PALBO_Fanconi_anemia", "df_X778_PALBO_HF_POL")

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
Curated_input_dfs <- list(df_X778_PALBO_MTOR_POS, df_X778_PALBO_MTOR_NEG, df_X778_PALBO_Proliferation_Science, df_X778_PALBO_Meta_PCNA, df_X778_PALBO_Proliferation_GO, df_X778_PALBO_Proliferation_STOP, df_X778_PALBO_G1_phase, df_X778_PALBO_Proliferation_Carey, df_X778_PALBO_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019, df_X778_PALBO_Essential_S_Phase_Viner_Breuer_Cell_Division_2019, df_X778_PALBO_HR, df_X778_PALBO_Fanconi_anemia, df_X778_PALBO_HF_POL)

# Apply the mean function to create a list of means
Curated_mean_test_df <- lapply(Curated_input_dfs, function(x) Mean_curated(x))

# creating a vector with desired names of the files containing the mean of each gene set log2FoldChange
Curated_Mean_X778_PALBO_dfs.output <- c("Mean_X778_PALBO_MTOR_POS", "Mean_X778_PALBO_MTOR_NEG", "Mean_X778_PALBO_Proliferation_Science", "Mean_X778_PALBO_Meta_PCNA", "Mean_X778_PALBO_Proliferation_GO", "Mean_X778_PALBO_Proliferation_STOP", "Mean_X778_PALBO_G1_phase", "Mean_X778_PALBO_Proliferation_Carey", "Mean_X778_PALBO_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019", "Mean_X778_PALBO_Essential_S_Phase_Viner_Breuer_Cell_Division_2019", "Mean_X778_PALBO_HR", "Mean_X778_PALBO_Fanconi_anemia", "Mean_X778_PALBO_HF_POL")

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
Curated_STDEV_X778_PALBO_dfs.output <- c("STDEV_X778_PALBO_MTOR_POS", "STDEV_X778_PALBO_MTOR_NEG", "STDEV_X778_PALBO_Proliferation_Science", "STDEV_X778_PALBO_Meta_PCNA", "STDEV_X778_PALBO_Proliferation_GO", "STDEV_X778_PALBO_Proliferation_STOP", "STDEV_X778_PALBO_G1_phase", "STDEV_X778_PALBO_Proliferation_Carey", "STDEV_X778_PALBO_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019", "STDEV_X778_PALBO_Essential_S_Phase_Viner_Breuer_Cell_Division_2019", "STDEV_X778_PALBO_HR", "STDEV_X778_PALBO_Fanconi_anemia", "STDEV_X778_PALBO_HF_POL")
# loop to assign data.frames to the desired variable name
for (i in (seq_along(Curated_STDEV_X778_PALBO_dfs.output))) {
  do.call("<-",list(Curated_STDEV_X778_PALBO_dfs.output[i], Curated_stdev_test_df[[i]]))
}




# Generate the data frame for the plot

# Generate the data frame for the plot
df_778_PALBO_proliferation <- data.frame(mean=c(Mean_X778_PALBO_MTOR_POS, Mean_X778_PALBO_MTOR_NEG, Mean_X778_PALBO_Proliferation_Science, Mean_X778_PALBO_Meta_PCNA, Mean_X778_PALBO_Proliferation_GO, Mean_X778_PALBO_Proliferation_STOP, Mean_X778_PALBO_G1_phase, Mean_X778_PALBO_Proliferation_Carey, Mean_X778_PALBO_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019, Mean_X778_PALBO_Essential_S_Phase_Viner_Breuer_Cell_Division_2019, Mean_X778_PALBO_HR, Mean_X778_PALBO_Fanconi_anemia, Mean_X778_PALBO_HF_POL),
                                   sd=c(STDEV_X778_PALBO_MTOR_POS, STDEV_X778_PALBO_MTOR_NEG, STDEV_X778_PALBO_Proliferation_Science, STDEV_X778_PALBO_Meta_PCNA, STDEV_X778_PALBO_Proliferation_GO, STDEV_X778_PALBO_Proliferation_STOP, STDEV_X778_PALBO_G1_phase, STDEV_X778_PALBO_Proliferation_Carey, STDEV_X778_PALBO_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019, STDEV_X778_PALBO_Essential_S_Phase_Viner_Breuer_Cell_Division_2019, STDEV_X778_PALBO_HR, STDEV_X778_PALBO_Fanconi_anemia, STDEV_X778_PALBO_HF_POL),
                                   Category=as.factor(c("X778_PALBO_MTOR_POS", "X778_PALBO_MTOR_NEG", "X778_PALBO_Proliferation_Science", "X778_PALBO_Meta_PCNA", "X778_PALBO_Proliferation_GO", "X778_PALBO_Proliferation_STOP", "X778_PALBO_G1_phase", "X778_PALBO_Proliferation_Carey", "X778_PALBO_Essential_cell_cycle", "X778_PALBO_Essential_S_Phase", "X778_PALBO_HR", "X778_PALBO_FA", "X778_PALBO_HF_POL")))
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
  scale_x_discrete(limits = c("X778_PALBO_Proliferation_Science", "X778_PALBO_Meta_PCNA", "X778_PALBO_Proliferation_GO", "X778_PALBO_Proliferation_STOP", "X778_PALBO_Proliferation_Carey", "X778_PALBO_G1_phase", "X778_PALBO_MTOR_POS", "X778_PALBO_MTOR_NEG", "X778_PALBO_Essential_cell_cycle", "X778_PALBO_Essential_S_Phase", "X778_PALBO_HR", "X778_PALBO_FA", "X778_PALBO_HF_POL")) +
  theme(plot.title = element_text(color="black", size=14, face="bold"))


# Merge all the data frames to generate box plots
# Add a column with the specific name of the gene set
df_X778_PALBO_MTOR_POS$Gene_set <- "X778_PALBO_MTOR_POS"
df_X778_PALBO_MTOR_NEG$Gene_set <- "X778_PALBO_MTOR_NEG"
df_X778_PALBO_Proliferation_Science$Gene_set <- "X778_PALBO_Proliferation_Science"
df_X778_PALBO_Meta_PCNA$Gene_set <- "X778_PALBO_Meta_PCNA"
df_X778_PALBO_Proliferation_GO$Gene_set <- "X778_PALBO_Proliferation_GO"
df_X778_PALBO_Proliferation_STOP$Gene_set <- "X778_PALBO_Proliferation_STOP"
df_X778_PALBO_G1_phase$Gene_set <- "X778_PALBO_G1_phase"
df_X778_PALBO_Proliferation_Carey$Gene_set <- "X778_PALBO_Proliferation_Carey"
df_X778_PALBO_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019$Gene_set <- "X778_PALBO_Essential_cell_cycle"
df_X778_PALBO_Essential_S_Phase_Viner_Breuer_Cell_Division_2019$Gene_set <- "X778_PALBO_Essential_S_Phase"
df_X778_PALBO_HR$Gene_set <- "X778_PALBO_HR"
df_X778_PALBO_Fanconi_anemia$Gene_set <- "X778_PALBO_Fanconi_anemia"
df_X778_PALBO_HF_POL$Gene_set <- "X778_PALBO_HF_POL"



# Combine the data frames
combined_df <- rbind(
  df_X778_PALBO_MTOR_POS,
  df_X778_PALBO_MTOR_NEG,
  df_X778_PALBO_Proliferation_Science,
  df_X778_PALBO_Meta_PCNA,
  df_X778_PALBO_Proliferation_GO,
  df_X778_PALBO_Proliferation_STOP,
  df_X778_PALBO_G1_phase,
  df_X778_PALBO_Proliferation_Carey,
  df_X778_PALBO_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019,
  df_X778_PALBO_Essential_S_Phase_Viner_Breuer_Cell_Division_2019,
  df_X778_PALBO_HR,
  df_X778_PALBO_Fanconi_anemia,
  df_X778_PALBO_HF_POL
)

# Ensure the Gene_set column is a factor with the correct order
combined_df$Gene_set <- factor(combined_df$Gene_set, levels = c("X778_PALBO_Proliferation_Science", "X778_PALBO_Meta_PCNA", "X778_PALBO_Proliferation_GO", "X778_PALBO_Proliferation_STOP", "X778_PALBO_Proliferation_Carey", "X778_PALBO_G1_phase", "X778_PALBO_MTOR_POS", "X778_PALBO_MTOR_NEG", "X778_PALBO_Essential_cell_cycle", "X778_PALBO_Essential_S_Phase", "X778_PALBO_HR", "X778_PALBO_Fanconi_anemia", "X778_PALBO_HF_POL"))


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
  ggtitle("Analysis cell proliferation 778 palbociclib 72h") +
  scale_x_discrete(limits = levels(combined_df$Gene_set)) +
  theme(plot.title = element_text(color="black", size=14, face="bold")) 





Mean_df_MTOR_POS_X778_PALBO <- mean(df_X778_PALBO_MTOR_POS$log2FoldChange)
Mean_df_MTOR_NEG_X778_PALBO <- mean(df_X778_PALBO_MTOR_NEG$log2FoldChange)
Mean_df_Proliferation_Science_X778_PALBO <- mean(df_X778_PALBO_Proliferation_Science$log2FoldChange)
Mean_df_Meta_PCNA_X778_PALBO <- mean(df_X778_PALBO_Meta_PCNA$log2FoldChange)
Mean_df_Proliferation_GO_X778_PALBO <- mean(df_X778_PALBO_Proliferation_GO$log2FoldChange)
Mean_df_Proliferation_STOP_X778_PALBO <- mean(df_X778_PALBO_Proliferation_STOP$log2FoldChange)
Mean_df_G1_phase_X778_PALBO <- mean(df_X778_PALBO_G1_phase$log2FoldChange)
Mean_df_Proliferation_Carey_X778_PALBO <- mean(df_X778_PALBO_Proliferation_Carey$log2FoldChange)
Mean_df_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019_X778_PALBO <- mean(df_X778_PALBO_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019$log2FoldChange)
Mean_df_Essential_S_Phase_Viner_Breuer_Cell_Division_2019_X778_PALBO <- mean(df_X778_PALBO_Essential_S_Phase_Viner_Breuer_Cell_Division_2019$log2FoldChange)
Mean_df_HR_X778_PALBO <- mean(df_X778_PALBO_HR$log2FoldChange)
Mean_df_FA_X778_PALBO <- mean(df_X778_PALBO_Fanconi_anemia$log2FoldChange)
Mean_df_HF_POL_X778_PALBO <- mean(df_X778_PALBO_HF_POL$log2FoldChange)

# Perform the Wilcoxon signed-rank test
wilcox_test_MTOR_POS_X778_PALBO <- wilcox.test(df_X778_PALBO_MTOR_POS$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_MTOR_NEG_X778_PALBO <- wilcox.test(df_X778_PALBO_MTOR_NEG$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_Proliferation_Science_X778_PALBO <- wilcox.test(df_X778_PALBO_Proliferation_Science$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_Meta_PCNA_X778_PALBO <- wilcox.test(df_X778_PALBO_Meta_PCNA$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_Proliferation_GO_X778_PALBO <- wilcox.test(df_X778_PALBO_Proliferation_GO$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_Proliferation_STOP_X778_PALBO <- wilcox.test(df_X778_PALBO_Proliferation_STOP$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_G1_phase_X778_PALBO <- wilcox.test(df_X778_PALBO_G1_phase$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_Proliferation_Carey_X778_PALBO <- wilcox.test(df_X778_PALBO_Proliferation_Carey$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019_X778_PALBO <- wilcox.test(df_X778_PALBO_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_Essential_S_Phase_Viner_Breuer_Cell_Division_2019_X778_PALBO <- wilcox.test(df_X778_PALBO_Essential_S_Phase_Viner_Breuer_Cell_Division_2019$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_HR_X778_PALBO <- wilcox.test(df_X778_PALBO_HR$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_FA_X778_PALBO <- wilcox.test(df_X778_PALBO_Fanconi_anemia$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_HF_POL_X778_PALBO <- wilcox.test(df_X778_PALBO_HF_POL$log2FoldChange, mu = 0, alternative = "less")















# SKBR3 cell line ANNOTATED transcripts ####
SKBR3_dexseq_ANNOTATED<-read.delim(file="path_to/72h_SKBR3_lapatinib.tx_counts.both_conditions.treated_untreated_dexseq.tab", header = TRUE, sep = "\t", quote = "") %>%
  dplyr::rename(ensembl_transcript_id_version=featureID) %>%
  dplyr::inner_join(Table_ANNOTATED, by="ensembl_transcript_id_version") %>%
  dplyr::filter(padj < 0.05) # The file 'Table' contains the Normalized_MFE values of all ANNOTATED transcripts


# Generate a scatter plot Log2_FC vs Normalized_MFE
library(ggplot2)
ggplot(SKBR3_dexseq_ANNOTATED, aes(y=Normalized_MFE, x=log2fold_treated_untreated)) +
  geom_point(size=1.5, alpha=0.3)+
  theme_classic() +
  theme(legend.position="none") +
  geom_vline(xintercept = 0, size=0.3) +
  #xlim(-5.5,5.5) +
  theme(axis.title = element_text(size = 20)) +
  theme(plot.title = element_text(size = 15)) +
  labs(y="Normalized_MFE",x="Log2_FC_expression",title="SKBR3_shared_ANNOTATED_transcripts_72h_Log2_FC_adj.P<0.05")


# It is well known that ANNOTATED transcripts, in general have lower MFEs than the UNANNOTATED transcripts

SKBR3_dexseq_ANNOTATED_UP <- dplyr::filter(SKBR3_dexseq_ANNOTATED, log2fold_treated_untreated > 0)
Mean_SKBR3_dexseq_ANNOTATED_UP <- mean(SKBR3_dexseq_ANNOTATED_UP$Normalized_MFE)

SKBR3_dexseq_ANNOTATED_DOWN <- dplyr::filter(SKBR3_dexseq_ANNOTATED, log2fold_treated_untreated < 0)
Mean_SKBR3_dexseq_ANNOTATED_DOWN <- mean(SKBR3_dexseq_ANNOTATED_DOWN$Normalized_MFE)

t.test(SKBR3_dexseq_ANNOTATED_DOWN$Normalized_MFE, SKBR3_dexseq_ANNOTATED_UP$Normalized_MFE, paired=FALSE)

# There is a significant increase of the mean MFE of the UP-regulated vs the DOWN-regulated ANNOTATED transcripts
# p-value = 7.576e-08



# Show the data with a violin plot
ggplot(SKBR3_dexseq_ANNOTATED, aes(x = factor(log2fold_treated_untreated > 0), y = Normalized_MFE)) +
  geom_violin(fill = "lightgrey", alpha = 0.7) +
  geom_jitter(width = 0.1, size = 1, alpha = 0.3) +
  theme_classic() +
  theme(legend.position = "none", axis.title = element_text(size = 20), plot.title = element_text(size = 10),
        axis.text.x = element_text(size = 16, color = "black")) +  # Increase font size and make x-axis labels black
  geom_vline(xintercept = 0, size = 0.3) +
  labs(y = "Normalized_MFE", x = NULL, title = "SKBR3_shared_ANNOTATED_transcripts_72h_Log2_FC_adj.P<0.05") +
  scale_x_discrete(labels = c("Log2_FC < 0", "Log2_FC > 0"))




# Generate a data frame with all the CUB values and plot the data in a density plot
max_length_SKBR3_MFE <- max(c(length(SKBR3_dexseq_ANNOTATED_UP$Normalized_MFE), length(SKBR3_dexseq_ANNOTATED_DOWN$Normalized_MFE)))

data_SKBR3_MFE_density <- data.frame(MFE_TX_UP_SKBR3 = c(SKBR3_dexseq_ANNOTATED_UP$Normalized_MFE, rep(NA, max_length_SKBR3_MFE - length(SKBR3_dexseq_ANNOTATED_UP$Normalized_MFE))), # This is useful when we need to create a data frame and the number of values between columns are different
                                     MFE_TX_DOWN_SKBR3 = c(SKBR3_dexseq_ANNOTATED_DOWN$Normalized_MFE, rep(NA, max_length_SKBR3_MFE- length(SKBR3_dexseq_ANNOTATED_DOWN$Normalized_MFE))))
library(reshape)
data_SKBR3_MFE_density <- melt(data_SKBR3_MFE_density)

library(ggplot2)
ggplot(data= data_SKBR3_MFE_density, aes(x=value, group=variable, color=variable)) +
  geom_density (size=1) +
  theme_classic() +
  #xlim(0,1.5) +
  labs(title = "MFE_distribution_ANNOTATED_Transcripts_SKBR3_72h")+
  labs(y="Density")+
  labs(x="DGfold / nt (kcal/mol)")


# Wilcoxon rank-sum test for SKBR3
wilcox_test_MFE_SKBR3_UP_vs_DOWN <- wilcox.test(
  SKBR3_dexseq_ANNOTATED_UP$Normalized_MFE,
  SKBR3_dexseq_ANNOTATED_DOWN$Normalized_MFE,
  alternative = "greater" # one-sided: median(UP) > median(DOWN)
)

print(wilcox_test_MFE_SKBR3_UP_vs_DOWN)












## Analysis net MFE and heat capacity ####
SKBR3_dexseq_ANNOTATED_norm_TX_counts_MFE_df <- SKBR3_dexseq_ANNOTATED %>%
  dplyr::select(ensembl_transcript_id_version, ensembl_transcript_id, ensembl_gene_id, external_gene_name, transcript_is_canonical, Transcript_length,
                countData.72h_SKBR3_untreated_R1, countData.72h_SKBR3_untreated_R2, countData.72h_SKBR3_untreated_R3,
                countData.72h_SKBR3_lapatinib_R1, countData.72h_SKBR3_lapatinib_R2, countData.72h_SKBR3_lapatinib_R3,
                Normalized_MFE, Cp_37C) 



# Identify the 75th percentile of the combined MFE distributions of control and drug-treated samples
quantile_075_SKBR3_dexseq_ANNOTATED_norm_TX_counts_MFE_df <- quantile(SKBR3_dexseq_ANNOTATED_norm_TX_counts_MFE_df$Normalized_MFE, probs = c(0.75))

# Select all the transcripts with Normalized_MFE values above the 75th percentile
SKBR3_dexseq_ANNOTATED_high_MFE_TX <- dplyr::filter(SKBR3_dexseq_ANNOTATED_norm_TX_counts_MFE_df, Normalized_MFE > quantile_075_SKBR3_dexseq_ANNOTATED_norm_TX_counts_MFE_df)

# Count the number of transcripts in the significantly up- and down-regulated pools that have MFE values above this threshold.
SKBR3_dexseq_ANNOTATED_high_MFE_TX_in_UP <- dplyr::filter(SKBR3_dexseq_ANNOTATED_UP, ensembl_transcript_id_version %in% SKBR3_dexseq_ANNOTATED_high_MFE_TX$ensembl_transcript_id_version)
SKBR3_dexseq_ANNOTATED_high_MFE_TX_in_DOWN <- dplyr::filter(SKBR3_dexseq_ANNOTATED_DOWN, ensembl_transcript_id_version %in% SKBR3_dexseq_ANNOTATED_high_MFE_TX$ensembl_transcript_id_version)


# Perform a Fisher’s exact test to compare the proportions of high MFE transcripts between the two groups.
SKBR3_dexseq_ANNOTATED_total_UP <- nrow(SKBR3_dexseq_ANNOTATED_UP)
SKBR3_dexseq_ANNOTATED_total_DOWN <- nrow(SKBR3_dexseq_ANNOTATED_DOWN)
SKBR3_dexseq_ANNOTATED_count_high_MFE_UP <- nrow(SKBR3_dexseq_ANNOTATED_high_MFE_TX_in_UP)
SKBR3_dexseq_ANNOTATED_count_high_MFE_DOWN <- nrow(SKBR3_dexseq_ANNOTATED_high_MFE_TX_in_DOWN)

SKBR3_dexseq_ANNOTATED_contingency_table <- matrix(c(SKBR3_dexseq_ANNOTATED_count_high_MFE_UP, SKBR3_dexseq_ANNOTATED_total_UP - SKBR3_dexseq_ANNOTATED_count_high_MFE_UP,
                                                           SKBR3_dexseq_ANNOTATED_count_high_MFE_DOWN, SKBR3_dexseq_ANNOTATED_total_DOWN - SKBR3_dexseq_ANNOTATED_count_high_MFE_DOWN),
                                                         nrow = 2, byrow = TRUE,
                                                         dimnames = list(UP = c("UP-regulated", "DOWN-regulated"),
                                                                         DOWN = c("High MFE TX", "Not High MFE TX")))

SKBR3_dexseq_ANNOTATED_contingency_table_df <- as.data.frame(SKBR3_dexseq_ANNOTATED_contingency_table)

SKBR3_dexseq_ANNOTATED_fisher_result <- fisher.test(SKBR3_dexseq_ANNOTATED_contingency_table)

cat("SKBR3_dexseq_ANNOTATED P-value from Fisher test:", SKBR3_dexseq_ANNOTATED_fisher_result$p.value, "\n")


# Look at the fractions
SKBR3_dexseq_ANNOTATED_high_MFE_in_UP_fraction <- SKBR3_dexseq_ANNOTATED_count_high_MFE_UP / SKBR3_dexseq_ANNOTATED_total_UP
SKBR3_dexseq_ANNOTATED_high_MFE_in_DOWN_fraction <- SKBR3_dexseq_ANNOTATED_count_high_MFE_DOWN / SKBR3_dexseq_ANNOTATED_total_DOWN

# Output the results
cat("Fraction of SKBR3_dexseq_ANNOTATED high MFE TX in UP-regulated genes:", SKBR3_dexseq_ANNOTATED_high_MFE_in_UP_fraction, "\n") # This value goes in the summary table
cat("Fraction of SKBR3_dexseq_ANNOTATED high MFE TX in DOWN-regulated genes:", SKBR3_dexseq_ANNOTATED_high_MFE_in_DOWN_fraction, "\n") # This value goes in the summary table
cat("Mean of SKBR3_dexseq_ANNOTATED high MFE TX in UP-regulated genes:", mean(SKBR3_dexseq_ANNOTATED_high_MFE_TX_in_UP$Normalized_MFE), "\n")
cat("Mean of SKBR3_dexseq_ANNOTATED high MFE TX in DOWN-regulated genes:", mean(SKBR3_dexseq_ANNOTATED_high_MFE_TX_in_DOWN$Normalized_MFE), "\n")
cat("Mode of SKBR3_dexseq_ANNOTATED high MFE TX in UP-regulated genes:", calculate_mode(SKBR3_dexseq_ANNOTATED_high_MFE_TX_in_UP$Normalized_MFE), "\n")
cat("Mode of SKBR3_dexseq_ANNOTATED high MFE TX in DOWN-regulated genes:", calculate_mode(SKBR3_dexseq_ANNOTATED_high_MFE_TX_in_DOWN$Normalized_MFE), "\n")



# Perform a two-proportion z-test to compare proportions
# Calculate the proportions and total counts
SKBR3_dexseq_ANNOTATED_high_MFE_in_UP_proportion <- SKBR3_dexseq_ANNOTATED_high_MFE_in_UP_fraction
SKBR3_dexseq_ANNOTATED_high_MFE_in_DOWN_proportion <- SKBR3_dexseq_ANNOTATED_high_MFE_in_DOWN_fraction

# Combined proportion for pooled standard error
combined_proportion <- (SKBR3_dexseq_ANNOTATED_count_high_MFE_UP + SKBR3_dexseq_ANNOTATED_count_high_MFE_DOWN) /
  (SKBR3_dexseq_ANNOTATED_total_UP + SKBR3_dexseq_ANNOTATED_total_DOWN)

# Standard error of the difference in proportions
standard_error <- sqrt(
  combined_proportion * (1 - combined_proportion) * 
    (1 / SKBR3_dexseq_ANNOTATED_total_UP + 1 / SKBR3_dexseq_ANNOTATED_total_DOWN)
)

# Z-statistic
z_statistic <- (SKBR3_dexseq_ANNOTATED_high_MFE_in_UP_proportion - SKBR3_dexseq_ANNOTATED_high_MFE_in_DOWN_proportion) /
  standard_error

# Two-tailed p-value
SKBR3_dexseq_ANNOTATED_z_test_MFE_p_value <- 2 * (1 - pnorm(abs(z_statistic)))





# Show data as violin plot
# Add a new column to indicate if the value is above the 75th percentile
SKBR3_dexseq_ANNOTATED$above_quantile_075 <- SKBR3_dexseq_ANNOTATED$Normalized_MFE > quantile_075_SKBR3_dexseq_ANNOTATED_norm_TX_counts_MFE_df

# Generate violin plot with conditional jitter coloring
ggplot(SKBR3_dexseq_ANNOTATED, aes(x = factor(log2fold_treated_untreated > 0), y = Normalized_MFE)) +
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
    title = "SKBR3_dexseq_ANNOTATED_Log2_FC_adj.P<0.05"
  ) +
  scale_x_discrete(labels = c("log2FoldChange < 0", "log2FoldChange > 0"))












## Fraction of total counts analysis ####
# Determine the impact of these changes in terms of fraction of total counts:
# Calculate the fraction of total counts of high MFE transcripts in control and drug-treated groups. 
# This is done by summing the counts of the high MFE transcripts and dividing by the total counts in each group.
# Fisher’s exact test to compare the proportions of high MFE transcripts counts between the two groups.



# SKBR3 cell line ANNOTATED transcripts ####
SKBR3_dexseq_ANNOTATED_whole_transcriptome <- read.delim(file="path_to/72h_SKBR3_lapatinib.tx_counts.both_conditions.treated_untreated_dexseq.tab", header = TRUE, sep = "\t", quote = "") %>%
  dplyr::rename(ensembl_transcript_id_version=featureID) %>%
  dplyr::inner_join(Table_ANNOTATED, by="ensembl_transcript_id_version") %>%
  dplyr::select(ensembl_transcript_id_version, ensembl_transcript_id, ensembl_gene_id, external_gene_name, transcript_is_canonical, Transcript_length,
                countData.72h_SKBR3_untreated_R1, countData.72h_SKBR3_untreated_R2, countData.72h_SKBR3_untreated_R3,
                countData.72h_SKBR3_lapatinib_R1, countData.72h_SKBR3_lapatinib_R2, countData.72h_SKBR3_lapatinib_R3,
                Normalized_MFE, Cp_37C) 
  
  
# Convert from wide to long format (one row per transcript x replicate)
SKBR3_dexseq_ANNOTATED_long_data <- SKBR3_dexseq_ANNOTATED_whole_transcriptome %>%
  tidyr::pivot_longer(
    cols = starts_with("countData.72h"),
    names_to = "replicate",
    values_to = "counts"
  )

# Add a condition column based on the replicate
SKBR3_dexseq_ANNOTATED_long_data <- SKBR3_dexseq_ANNOTATED_long_data %>%
  dplyr::mutate(
    condition = dplyr::if_else(
      replicate %in% c(
        "countData.72h_SKBR3_untreated_R1",
        "countData.72h_SKBR3_untreated_R2",
        "countData.72h_SKBR3_untreated_R3"
      ),
      "Control", "lapatinib"
    ),
    condition = factor(condition, levels = c("Control", "lapatinib"))
  )

# Aggregate to average counts per transcript x condition
SKBR3_dexseq_ANNOTATED_norm_TX_aggregated_data_avg_counts <- SKBR3_dexseq_ANNOTATED_long_data %>%
  dplyr::group_by(ensembl_transcript_id_version, condition, Normalized_MFE) %>%
  dplyr::summarise(
    avg_counts = mean(counts, na.rm = TRUE),
    .groups = "drop"
  )







# High MFE (>= 75th percentile)
# Test whether cells exposed to lapatinib have a higher fraction of transcripts in the 75th percentile of MFE values compared to control cells.

# 75th percentile threshold
mfe_threshold_SKBR3 <- quantile(
  SKBR3_dexseq_ANNOTATED_norm_TX_aggregated_data_avg_counts$Normalized_MFE,
  0.75, na.rm = TRUE
)

# Label high MFE
df_SKBR3_high <- SKBR3_dexseq_ANNOTATED_norm_TX_aggregated_data_avg_counts %>%
  dplyr::mutate(MFE_class = ifelse(Normalized_MFE >= mfe_threshold_SKBR3, "High_MFE", "Other"))

# Aggregate counts per condition x class
counts_table_SKBR3_high <- df_SKBR3_high %>%
  dplyr::group_by(condition, MFE_class) %>%
  dplyr::summarise(total_counts = sum(avg_counts, na.rm = TRUE), .groups = "drop")

# 2x2 contingency (row1=lapatinib, row2=Control; col1=High_MFE, col2=Other)
contingency_SKBR3_high <- counts_table_SKBR3_high %>%
  tidyr::pivot_wider(names_from = MFE_class, values_from = total_counts, values_fill = 0) %>%
  dplyr::arrange(factor(condition, levels = c("lapatinib", "Control")))

fisher_matrix_SKBR3_high <- as.matrix(contingency_SKBR3_high[, c("High_MFE", "Other")])
rownames(fisher_matrix_SKBR3_high) <- contingency_SKBR3_high$condition

# One-sided Fisher: lapatinib > Control
fisher_high_SKBR3_result <- fisher.test(fisher_matrix_SKBR3_high, alternative = "greater")
fisher_high_SKBR3_result

# Report High_MFE fractions per condition
high_mfe_fractions_SKBR3 <- counts_table_SKBR3_high %>%
  dplyr::group_by(condition) %>%
  dplyr::mutate(fraction = total_counts / sum(total_counts)) %>%
  dplyr::filter(MFE_class == "High_MFE")
high_mfe_fractions_SKBR3






SKBR3_dexseq_ANNOTATED_lapatinib_high_MFE_TX_counts_fraction <-
  high_mfe_fractions_SKBR3 %>%
  dplyr::filter(condition == "lapatinib") %>%
  dplyr::pull(fraction) %>%
  as.numeric()

SKBR3_dexseq_ANNOTATED_control_high_MFE_TX_counts_fraction <-
  high_mfe_fractions_SKBR3 %>%
  dplyr::filter(condition == "Control") %>%
  dplyr::pull(fraction) %>%
  as.numeric()










# Low MFE (<= 25th percentile)
# Test whether control cells have a higher fraction of transcripts in the 25th percentile of MFE values compared to lapatinib-treated cells.
# 25th percentile threshold

mfe_threshold_25_SKBR3 <- quantile(
  SKBR3_dexseq_ANNOTATED_norm_TX_aggregated_data_avg_counts$Normalized_MFE,
  0.25, na.rm = TRUE
)

# Label low MFE
df_SKBR3_low <- SKBR3_dexseq_ANNOTATED_norm_TX_aggregated_data_avg_counts %>%
  dplyr::mutate(MFE_class = ifelse(Normalized_MFE <= mfe_threshold_25_SKBR3, "Low_MFE", "Other"))

# Aggregate counts per condition x class
counts_table_SKBR3_low <- df_SKBR3_low %>%
  dplyr::group_by(condition, MFE_class) %>%
  dplyr::summarise(total_counts = sum(avg_counts, na.rm = TRUE), .groups = "drop")

# 2x2 contingency (row1=Control, row2=lapatinib; col1=Low_MFE, col2=Other)
contingency_SKBR3_low <- counts_table_SKBR3_low %>%
  tidyr::pivot_wider(names_from = MFE_class, values_from = total_counts, values_fill = 0) %>%
  dplyr::arrange(factor(condition, levels = c("Control", "lapatinib")))

fisher_matrix_SKBR3_low <- as.matrix(contingency_SKBR3_low[, c("Low_MFE", "Other")])
rownames(fisher_matrix_SKBR3_low) <- contingency_SKBR3_low$condition

# One-sided Fisher: Control > lapatinib (higher Low_MFE fraction)
fisher_low_SKBR3_result <- fisher.test(fisher_matrix_SKBR3_low, alternative = "greater")
fisher_low_SKBR3_result

# Report Low_MFE fractions per condition
low_mfe_fractions_SKBR3 <- counts_table_SKBR3_low %>%
  dplyr::group_by(condition) %>%
  dplyr::mutate(fraction = total_counts / sum(total_counts)) %>%
  dplyr::filter(MFE_class == "Low_MFE")
low_mfe_fractions_SKBR3







## SKBR3_dexseq_ANNOTATED heat capacity changes ####
# This is calculated on the total n. of transcripts with a significant log2 fold change
# Generate a scatter plot Log2_FC vs heat capacity

ggplot(SKBR3_dexseq_ANNOTATED, aes(y=Cp_37C, x=log2fold_treated_untreated)) +
  geom_point(size=1.5, alpha=0.1)+
  theme_classic() +
  theme(legend.position="none") +
  geom_vline(xintercept = 0, size=0.3) +
  xlim(-12, 12) +
  ylim(0,400) +
  theme(axis.title = element_text(size = 20)) +
  theme(plot.title = element_text(size = 15)) +
  labs(y="Heat capacity (J/C)",x="Log2_FC_expression",title="SKBR3_dexseq_ANNOTATED_Log2_FC_adj.P<0.05")


# Show the data with a box plot
ggplot(SKBR3_dexseq_ANNOTATED, aes(x = factor(log2fold_treated_untreated > 0), y = Cp_37C)) +
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
  labs(y="Heat capacity (J/C)",x="",title="SKBR3_dexseq_ANNOTATED_Log2_FC_adj.P<0.05 heat_capacity all transcripts") +
  scale_x_discrete(labels = c("log2FoldChange < 0", "log2FoldChange > 0"))

Mean_heat_capacity_SKBR3_dexseq_ANNOTATED_UP <- mean(SKBR3_dexseq_ANNOTATED_UP$Cp_37C)
Mean_heat_capacity_SKBR3_dexseq_ANNOTATED_DOWN <- mean(SKBR3_dexseq_ANNOTATED_DOWN$Cp_37C)
TTEST_heat_capacity_SKBR3_dexseq_ANNOTATED_all_UP_vs_DOWN <- t.test(SKBR3_dexseq_ANNOTATED_DOWN$Cp_37C, SKBR3_dexseq_ANNOTATED_UP$Cp_37C, paired=FALSE)
print(TTEST_heat_capacity_SKBR3_dexseq_ANNOTATED_all_UP_vs_DOWN)
# Add K-S test to look at the differences between the up-regukated and down-regulated genes
ks_test_SKBR3_dexseq_ANNOTATED_all_UP_vs_DOWN_heat_capacity <- ks.test(SKBR3_dexseq_ANNOTATED_DOWN$Cp_37C, SKBR3_dexseq_ANNOTATED_UP$Cp_37C)
print(ks_test_SKBR3_dexseq_ANNOTATED_all_UP_vs_DOWN_heat_capacity)




# Calculate the changes of heat capacity only in the high MFE transcripts

# Contribution of the UP or DOWN-regulated high MFE transcripts to the heat capacity
# Look at the heat capacity of the UP-regulated high MFE transcripts
SKBR3_dexseq_ANNOTATED_high_MFE_TX_in_UP_heat_capacity <- mean(SKBR3_dexseq_ANNOTATED_high_MFE_TX_in_UP$Cp_37C)

# Look at the heat capacity of the DOWN-regulated high MFE transcripts
SKBR3_dexseq_ANNOTATED_high_MFE_TX_in_DOWN_heat_capacity <- mean(SKBR3_dexseq_ANNOTATED_high_MFE_TX_in_DOWN$Cp_37C)

# Heat capacity gain in UP-regulated high MFE transcripts vs DOWN-regulated high MFE transcripts
SKBR3_dexseq_ANNOTATED_heat_capacity_gain_UP_vs_DOWN <- SKBR3_dexseq_ANNOTATED_high_MFE_TX_in_UP_heat_capacity - SKBR3_dexseq_ANNOTATED_high_MFE_TX_in_DOWN_heat_capacity
cat("SKBR3_dexseq_ANNOTATED net heat capacity UP vs DOWN high MFE TX:", SKBR3_dexseq_ANNOTATED_heat_capacity_gain_UP_vs_DOWN, "\n")


# Extract specific heat data
SKBR3_dexseq_ANNOTATED_Cp_37C_UP <- SKBR3_dexseq_ANNOTATED_high_MFE_TX_in_UP$Cp_37C
SKBR3_dexseq_ANNOTATED_Cp_37C_DOWN <- SKBR3_dexseq_ANNOTATED_high_MFE_TX_in_DOWN$Cp_37C

# Combine data into a single data frame for easier plotting
SKBR3_dexseq_ANNOTATED_Cp_37C_data <- data.frame(
  SKBR3_dexseq_ANNOTATED_Cp_37C = c(SKBR3_dexseq_ANNOTATED_Cp_37C_UP, SKBR3_dexseq_ANNOTATED_Cp_37C_DOWN),
  regulation = factor(c(rep("UP", length(SKBR3_dexseq_ANNOTATED_Cp_37C_UP)), rep("DOWN", length(SKBR3_dexseq_ANNOTATED_Cp_37C_DOWN))))
)

# Generate box plot
ggplot(SKBR3_dexseq_ANNOTATED_Cp_37C_data, aes(x = regulation, y = SKBR3_dexseq_ANNOTATED_Cp_37C, fill = regulation)) +
  geom_boxplot(fill = "lightgrey", alpha = 0.7) +
  labs(title = "SKBR3_dexseq_ANNOTATED specific Heat of High MFE Transcripts",
       x = "",
       y = "Specific Heat (J/C)") +
  scale_x_discrete(labels = c("log2_FC < 0", "log2_FC > 0")) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 16),
        axis.title = element_text(size = 25))

Mean_heat_capacity_SKBR3_dexseq_ANNOTATED_high_MFE_TX_UP <- mean(SKBR3_dexseq_ANNOTATED_high_MFE_TX_in_UP$Cp_37C)
Mean_heat_capacity_SKBR3_dexseq_ANNOTATED_high_MFE_TX_DOWN <- mean(SKBR3_dexseq_ANNOTATED_high_MFE_TX_in_DOWN$Cp_37C)
TTEST_heat_capacity_SKBR3_dexseq_ANNOTATED_all_high_MFE_TX_UP_vs_DOWN <- t.test(SKBR3_dexseq_ANNOTATED_high_MFE_TX_in_DOWN$Cp_37C, SKBR3_dexseq_ANNOTATED_high_MFE_TX_in_UP$Cp_37C, paired=FALSE)
print(TTEST_heat_capacity_SKBR3_dexseq_ANNOTATED_all_high_MFE_TX_UP_vs_DOWN)

# Add K-S test to look at the differences between the up-regukated and down-regulated genes
ks_test_heat_capacity_SKBR3_dexseq_ANNOTATED_all_high_MFE_TX_UP_vs_DOWN <- ks.test(SKBR3_dexseq_ANNOTATED_high_MFE_TX_in_DOWN$Cp_37C, SKBR3_dexseq_ANNOTATED_high_MFE_TX_in_UP$Cp_37C)
print(ks_test_heat_capacity_SKBR3_dexseq_ANNOTATED_all_high_MFE_TX_UP_vs_DOWN)



















# Analysis cell proliferation Nanopore data SKBR3 LAP ####

SKBR3_LAP_DESeq2 <- read.delim(file="path_to/72h_SKBR3_lapatinib.tx_counts.both_conditions.treated_untreated_deseq2.tab", header = TRUE, sep = "\t", quote = "")

# Filter the proliferation gene sets by using a function
# Use a function also to determine the mean and standard deviation of each set and the plot all of them together


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
CURATED_PATHWAYS.output <- c("SKBR3_LAP_MTOR_POS", "SKBR3_LAP_MTOR_NEG", "SKBR3_LAP_Proliferation_Science", "SKBR3_LAP_Meta_PCNA", "SKBR3_LAP_Proliferation_GO", "SKBR3_LAP_Proliferation_STOP", "SKBR3_LAP_G1_phase", "SKBR3_LAP_Proliferation_Carey", "SKBR3_LAP_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019", "SKBR3_LAP_Essential_S_Phase_Viner_Breuer_Cell_Division_2019", "SKBR3_LAP_HR", "SKBR3_LAP_Fanconi_anemia", "SKBR3_LAP_HF_POL")

# Assign the names of the gene sets to each entry in the list
names(CURATED_PATHWAYS.input) <- CURATED_PATHWAYS.output

# Apply the filtering function to create a list of data frames
Curated_test_df <- lapply(CURATED_PATHWAYS.input, function(x) filtering(x))

# creating a vector with desired names of filtered data frames
Curated_dfs.output <- c("df_SKBR3_LAP_MTOR_POS", "df_SKBR3_LAP_MTOR_NEG", "df_SKBR3_LAP_Proliferation_Science", "df_SKBR3_LAP_Meta_PCNA", "df_SKBR3_LAP_Proliferation_GO", "df_SKBR3_LAP_Proliferation_STOP", "df_SKBR3_LAP_G1_phase", "df_SKBR3_LAP_Proliferation_Carey", "df_SKBR3_LAP_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019", "df_SKBR3_LAP_Essential_S_Phase_Viner_Breuer_Cell_Division_2019", "df_SKBR3_LAP_HR", "df_SKBR3_LAP_Fanconi_anemia", "df_SKBR3_LAP_HF_POL")

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
Curated_input_dfs <- list(df_SKBR3_LAP_MTOR_POS, df_SKBR3_LAP_MTOR_NEG, df_SKBR3_LAP_Proliferation_Science, df_SKBR3_LAP_Meta_PCNA, df_SKBR3_LAP_Proliferation_GO, df_SKBR3_LAP_Proliferation_STOP, df_SKBR3_LAP_G1_phase, df_SKBR3_LAP_Proliferation_Carey, df_SKBR3_LAP_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019, df_SKBR3_LAP_Essential_S_Phase_Viner_Breuer_Cell_Division_2019, df_SKBR3_LAP_HR, df_SKBR3_LAP_Fanconi_anemia, df_SKBR3_LAP_HF_POL)

# Apply the mean function to create a list of means
Curated_mean_test_df <- lapply(Curated_input_dfs, function(x) Mean_curated(x))

# creating a vector with desired names of the files containing the mean of each gene set log2FoldChange
Curated_Mean_SKBR3_LAP_dfs.output <- c("Mean_SKBR3_LAP_MTOR_POS", "Mean_SKBR3_LAP_MTOR_NEG", "Mean_SKBR3_LAP_Proliferation_Science", "Mean_SKBR3_LAP_Meta_PCNA", "Mean_SKBR3_LAP_Proliferation_GO", "Mean_SKBR3_LAP_Proliferation_STOP", "Mean_SKBR3_LAP_G1_phase", "Mean_SKBR3_LAP_Proliferation_Carey", "Mean_SKBR3_LAP_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019", "Mean_SKBR3_LAP_Essential_S_Phase_Viner_Breuer_Cell_Division_2019", "Mean_SKBR3_LAP_HR", "Mean_SKBR3_LAP_Fanconi_anemia", "Mean_SKBR3_LAP_HF_POL")

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
Curated_STDEV_SKBR3_LAP_dfs.output <- c("STDEV_SKBR3_LAP_MTOR_POS", "STDEV_SKBR3_LAP_MTOR_NEG", "STDEV_SKBR3_LAP_Proliferation_Science", "STDEV_SKBR3_LAP_Meta_PCNA", "STDEV_SKBR3_LAP_Proliferation_GO", "STDEV_SKBR3_LAP_Proliferation_STOP", "STDEV_SKBR3_LAP_G1_phase", "STDEV_SKBR3_LAP_Proliferation_Carey", "STDEV_SKBR3_LAP_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019", "STDEV_SKBR3_LAP_Essential_S_Phase_Viner_Breuer_Cell_Division_2019", "STDEV_SKBR3_LAP_HR", "STDEV_SKBR3_LAP_Fanconi_anemia", "STDEV_SKBR3_LAP_HF_POL")
# loop to assign data.frames to the desired variable name
for (i in (seq_along(Curated_STDEV_SKBR3_LAP_dfs.output))) {
  do.call("<-",list(Curated_STDEV_SKBR3_LAP_dfs.output[i], Curated_stdev_test_df[[i]]))
}




# Generate the data frame for the plot

# Generate the data frame for the plot
df_SKBR3_LAP_proliferation <- data.frame(mean=c(Mean_SKBR3_LAP_MTOR_POS, Mean_SKBR3_LAP_MTOR_NEG, Mean_SKBR3_LAP_Proliferation_Science, Mean_SKBR3_LAP_Meta_PCNA, Mean_SKBR3_LAP_Proliferation_GO, Mean_SKBR3_LAP_Proliferation_STOP, Mean_SKBR3_LAP_G1_phase, Mean_SKBR3_LAP_Proliferation_Carey, Mean_SKBR3_LAP_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019, Mean_SKBR3_LAP_Essential_S_Phase_Viner_Breuer_Cell_Division_2019, Mean_SKBR3_LAP_HR, Mean_SKBR3_LAP_Fanconi_anemia, Mean_SKBR3_LAP_HF_POL),
                                          sd=c(STDEV_SKBR3_LAP_MTOR_POS, STDEV_SKBR3_LAP_MTOR_NEG, STDEV_SKBR3_LAP_Proliferation_Science, STDEV_SKBR3_LAP_Meta_PCNA, STDEV_SKBR3_LAP_Proliferation_GO, STDEV_SKBR3_LAP_Proliferation_STOP, STDEV_SKBR3_LAP_G1_phase, STDEV_SKBR3_LAP_Proliferation_Carey, STDEV_SKBR3_LAP_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019, STDEV_SKBR3_LAP_Essential_S_Phase_Viner_Breuer_Cell_Division_2019, STDEV_SKBR3_LAP_HR, STDEV_SKBR3_LAP_Fanconi_anemia, STDEV_SKBR3_LAP_HF_POL),
                                          Category=as.factor(c("SKBR3_LAP_MTOR_POS", "SKBR3_LAP_MTOR_NEG", "SKBR3_LAP_Proliferation_Science", "SKBR3_LAP_Meta_PCNA", "SKBR3_LAP_Proliferation_GO", "SKBR3_LAP_Proliferation_STOP", "SKBR3_LAP_G1_phase", "SKBR3_LAP_Proliferation_Carey", "SKBR3_LAP_Essential_cell_cycle", "SKBR3_LAP_Essential_S_Phase", "SKBR3_LAP_HR", "SKBR3_LAP_FA", "SKBR3_LAP_HF_POL")))
df_SKBR3_LAP_proliferation




library(ggplot2)
ggplot(df_SKBR3_LAP_proliferation, aes(x=reorder(Category, -mean), y=mean)) +
  theme_bw() +
  geom_point(position=position_dodge(width = 1), stat="identity",colour='black') +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.4) +
  theme(axis.text.x = element_text(size=12, angle=45,hjust=0.95,vjust=0.98)) +
  theme(axis.text.y = element_text(size=14)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 10))) +
  geom_point(colour = "black", size = 3) +
  xlab("") + ylab("Relative log2 fold change") +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  ggtitle("Analysis cell proliferation 778 NUTLIN NANOPORE 72h vs T0") +
  scale_x_discrete(limits = c("SKBR3_LAP_Proliferation_Science", "SKBR3_LAP_Meta_PCNA", "SKBR3_LAP_Proliferation_GO", "SKBR3_LAP_Proliferation_STOP", "SKBR3_LAP_Proliferation_Carey", "SKBR3_LAP_G1_phase", "SKBR3_LAP_MTOR_POS", "SKBR3_LAP_MTOR_NEG", "SKBR3_LAP_Essential_cell_cycle", "SKBR3_LAP_Essential_S_Phase", "SKBR3_LAP_HR", "SKBR3_LAP_FA", "SKBR3_LAP_HF_POL")) +
  theme(plot.title = element_text(color="black", size=14, face="bold"))


# Merge all the data frames to generate box plots
# Add a column with the specific name of the gene set
df_SKBR3_LAP_MTOR_POS$Gene_set <- "SKBR3_LAP_MTOR_POS"
df_SKBR3_LAP_MTOR_NEG$Gene_set <- "SKBR3_LAP_MTOR_NEG"
df_SKBR3_LAP_Proliferation_Science$Gene_set <- "SKBR3_LAP_Proliferation_Science"
df_SKBR3_LAP_Meta_PCNA$Gene_set <- "SKBR3_LAP_Meta_PCNA"
df_SKBR3_LAP_Proliferation_GO$Gene_set <- "SKBR3_LAP_Proliferation_GO"
df_SKBR3_LAP_Proliferation_STOP$Gene_set <- "SKBR3_LAP_Proliferation_STOP"
df_SKBR3_LAP_G1_phase$Gene_set <- "SKBR3_LAP_G1_phase"
df_SKBR3_LAP_Proliferation_Carey$Gene_set <- "SKBR3_LAP_Proliferation_Carey"
df_SKBR3_LAP_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019$Gene_set <- "SKBR3_LAP_Essential_cell_cycle"
df_SKBR3_LAP_Essential_S_Phase_Viner_Breuer_Cell_Division_2019$Gene_set <- "SKBR3_LAP_Essential_S_Phase"
df_SKBR3_LAP_HR$Gene_set <- "SKBR3_LAP_HR"
df_SKBR3_LAP_Fanconi_anemia$Gene_set <- "SKBR3_LAP_Fanconi_anemia"
df_SKBR3_LAP_HF_POL$Gene_set <- "SKBR3_LAP_HF_POL"



# Combine the data frames
combined_df <- rbind(
  df_SKBR3_LAP_MTOR_POS,
  df_SKBR3_LAP_MTOR_NEG,
  df_SKBR3_LAP_Proliferation_Science,
  df_SKBR3_LAP_Meta_PCNA,
  df_SKBR3_LAP_Proliferation_GO,
  df_SKBR3_LAP_Proliferation_STOP,
  df_SKBR3_LAP_G1_phase,
  df_SKBR3_LAP_Proliferation_Carey,
  df_SKBR3_LAP_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019,
  df_SKBR3_LAP_Essential_S_Phase_Viner_Breuer_Cell_Division_2019,
  df_SKBR3_LAP_HR,
  df_SKBR3_LAP_Fanconi_anemia,
  df_SKBR3_LAP_HF_POL
)

# Ensure the Gene_set column is a factor with the correct order
combined_df$Gene_set <- factor(combined_df$Gene_set, levels = c("SKBR3_LAP_Proliferation_Science", "SKBR3_LAP_Meta_PCNA", "SKBR3_LAP_Proliferation_GO", "SKBR3_LAP_Proliferation_STOP", "SKBR3_LAP_Proliferation_Carey", "SKBR3_LAP_G1_phase", "SKBR3_LAP_MTOR_POS", "SKBR3_LAP_MTOR_NEG", "SKBR3_LAP_Essential_cell_cycle", "SKBR3_LAP_Essential_S_Phase", "SKBR3_LAP_HR", "SKBR3_LAP_Fanconi_anemia", "SKBR3_LAP_HF_POL"))


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
  ggtitle("Analysis cell proliferation SKBR3 lapatinib 72h") +
  scale_x_discrete(limits = levels(combined_df$Gene_set)) +
  theme(plot.title = element_text(color="black", size=14, face="bold")) 








Mean_df_MTOR_POS_SKBR3_LAP <- mean(df_SKBR3_LAP_MTOR_POS$log2FoldChange)
Mean_df_MTOR_NEG_SKBR3_LAP <- mean(df_SKBR3_LAP_MTOR_NEG$log2FoldChange)
Mean_df_Proliferation_Science_SKBR3_LAP <- mean(df_SKBR3_LAP_Proliferation_Science$log2FoldChange)
Mean_df_Meta_PCNA_SKBR3_LAP <- mean(df_SKBR3_LAP_Meta_PCNA$log2FoldChange)
Mean_df_Proliferation_GO_SKBR3_LAP <- mean(df_SKBR3_LAP_Proliferation_GO$log2FoldChange)
Mean_df_Proliferation_STOP_SKBR3_LAP <- mean(df_SKBR3_LAP_Proliferation_STOP$log2FoldChange)
Mean_df_G1_phase_SKBR3_LAP <- mean(df_SKBR3_LAP_G1_phase$log2FoldChange)
Mean_df_Proliferation_Carey_SKBR3_LAP <- mean(df_SKBR3_LAP_Proliferation_Carey$log2FoldChange)
Mean_df_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019_SKBR3_LAP <- mean(df_SKBR3_LAP_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019$log2FoldChange)
Mean_df_Essential_S_Phase_Viner_Breuer_Cell_Division_2019_SKBR3_LAP <- mean(df_SKBR3_LAP_Essential_S_Phase_Viner_Breuer_Cell_Division_2019$log2FoldChange)
Mean_df_HR_SKBR3_LAP <- mean(df_SKBR3_LAP_HR$log2FoldChange)
Mean_df_FA_SKBR3_LAP <- mean(df_SKBR3_LAP_Fanconi_anemia$log2FoldChange)
Mean_df_HF_POL_SKBR3_LAP <- mean(df_SKBR3_LAP_HF_POL$log2FoldChange)

# Perform the Wilcoxon signed-rank test
wilcox_test_MTOR_POS_SKBR3_LAP <- wilcox.test(df_SKBR3_LAP_MTOR_POS$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_MTOR_NEG_SKBR3_LAP <- wilcox.test(df_SKBR3_LAP_MTOR_NEG$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_Proliferation_Science_SKBR3_LAP <- wilcox.test(df_SKBR3_LAP_Proliferation_Science$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_Meta_PCNA_SKBR3_LAP <- wilcox.test(df_SKBR3_LAP_Meta_PCNA$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_Proliferation_GO_SKBR3_LAP <- wilcox.test(df_SKBR3_LAP_Proliferation_GO$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_Proliferation_STOP_SKBR3_LAP <- wilcox.test(df_SKBR3_LAP_Proliferation_STOP$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_G1_phase_SKBR3_LAP <- wilcox.test(df_SKBR3_LAP_G1_phase$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_Proliferation_Carey_SKBR3_LAP <- wilcox.test(df_SKBR3_LAP_Proliferation_Carey$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019_SKBR3_LAP <- wilcox.test(df_SKBR3_LAP_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_Essential_S_Phase_Viner_Breuer_Cell_Division_2019_SKBR3_LAP <- wilcox.test(df_SKBR3_LAP_Essential_S_Phase_Viner_Breuer_Cell_Division_2019$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_HR_SKBR3_LAP <- wilcox.test(df_SKBR3_LAP_HR$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_FA_SKBR3_LAP <- wilcox.test(df_SKBR3_LAP_Fanconi_anemia$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_HF_POL_SKBR3_LAP <- wilcox.test(df_SKBR3_LAP_HF_POL$log2FoldChange, mu = 0, alternative = "less")
































# SKMEL28 cell line ANNOTATED transcripts ####
SKMEL28_dexseq_ANNOTATED<-read.delim(file="path_to/72h_SKMEL28_vemurafenib.tx_counts.both_conditions.treated_untreated_dexseq.tab", header = TRUE, sep = "\t", quote = "") %>%
  dplyr::rename(ensembl_transcript_id_version=featureID) %>%
  dplyr::inner_join(Table_ANNOTATED, by="ensembl_transcript_id_version") %>%
  dplyr::filter(padj < 0.05) # The file 'Table' contains the Normalized_MFE values of all ANNOTATED transcripts


# Generate a scatter plot Log2_FC vs Normalized_MFE
library(ggplot2)
ggplot(SKMEL28_dexseq_ANNOTATED, aes(y=Normalized_MFE, x=log2fold_treated_untreated)) +
  geom_point(size=1.5, alpha=0.3)+
  theme_classic() +
  theme(legend.position="none") +
  geom_vline(xintercept = 0, size=0.3) +
  #xlim(-5.5,5.5) +
  theme(axis.title = element_text(size = 20)) +
  theme(plot.title = element_text(size = 15)) +
  labs(y="Normalized_MFE",x="Log2_FC_expression",title="SKMEL28_shared_ANNOTATED_transcripts_72h_Log2_FC_adj.P<0.05")



# Show the data with a violin plot
ggplot(SKMEL28_dexseq_ANNOTATED, aes(x = factor(log2fold_treated_untreated > 0), y = Normalized_MFE)) +
  geom_violin(fill = "lightgrey", alpha = 0.7) +
  geom_jitter(width = 0.1, size = 1, alpha = 0.3) +
  theme_classic() +
  theme(legend.position = "none", axis.title = element_text(size = 20), plot.title = element_text(size = 10),
        axis.text.x = element_text(size = 16, color = "black")) +  # Increase font size and make x-axis labels black
  geom_vline(xintercept = 0, size = 0.3) +
  labs(y = "Normalized_MFE", x = NULL, title = "SKMEL28_shared_ANNOTATED_transcripts_72h_Log2_FC_adj.P<0.05") +
  scale_x_discrete(labels = c("Log2_FC < 0", "Log2_FC > 0"))
# In this modified code, the x aesthetic is set as a factor variable factor(log2fold_treated_untreated > 0). 
# This creates two separate groups based on the condition Log2_FC > 0 or Log2_FC < 0. Each group is represented by a different violin shape.





SKMEL28_dexseq_ANNOTATED_UP <- dplyr::filter(SKMEL28_dexseq_ANNOTATED, log2fold_treated_untreated > 0)
Mean_SKMEL28_dexseq_ANNOTATED_UP <- mean(SKMEL28_dexseq_ANNOTATED_UP$Normalized_MFE)

SKMEL28_dexseq_ANNOTATED_DOWN <- dplyr::filter(SKMEL28_dexseq_ANNOTATED, log2fold_treated_untreated < 0)
Mean_SKMEL28_dexseq_ANNOTATED_DOWN <- mean(SKMEL28_dexseq_ANNOTATED_DOWN$Normalized_MFE)


# Split into UP and DOWN
SKMEL28_dexseq_ANNOTATED_UP <- dplyr::filter(SKMEL28_dexseq_ANNOTATED, log2fold_treated_untreated > 0)
SKMEL28_dexseq_ANNOTATED_DOWN <- dplyr::filter(SKMEL28_dexseq_ANNOTATED, log2fold_treated_untreated < 0)

# Wilcoxon rank-sum test
wilcox_test_MFE_SKMEL28_UP_vs_DOWN <- wilcox.test(
  SKMEL28_dexseq_ANNOTATED_UP$Normalized_MFE,
  SKMEL28_dexseq_ANNOTATED_DOWN$Normalized_MFE,
  alternative = "greater" # one-sided: median(UP) > median(DOWN)
)

print(wilcox_test_MFE_SKMEL28_UP_vs_DOWN)







# Generate a data frame with all the CUB values and plot the data in a density plot
max_length_SKMEL28_MFE <- max(c(length(SKMEL28_dexseq_ANNOTATED_UP$Normalized_MFE), length(SKMEL28_dexseq_ANNOTATED_DOWN$Normalized_MFE)))

data_SKMEL28_MFE_density <- data.frame(MFE_TX_UP_SKMEL28 = c(SKMEL28_dexseq_ANNOTATED_UP$Normalized_MFE, rep(NA, max_length_SKMEL28_MFE - length(SKMEL28_dexseq_ANNOTATED_UP$Normalized_MFE))), # This is useful when we need to create a data frame and the number of values between columns are different
                                       MFE_TX_DOWN_SKMEL28 = c(SKMEL28_dexseq_ANNOTATED_DOWN$Normalized_MFE, rep(NA, max_length_SKMEL28_MFE- length(SKMEL28_dexseq_ANNOTATED_DOWN$Normalized_MFE))))
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



SKMEL28_dexseq_ANNOTATED_norm_TX_counts_MFE_df <- SKMEL28_dexseq_ANNOTATED %>%
  dplyr::select(ensembl_transcript_id_version, ensembl_transcript_id, ensembl_gene_id, external_gene_name, transcript_is_canonical, Transcript_length,
                countData.72h_SKMEL28_untreated_R1, countData.72h_SKMEL28_untreated_R2, countData.72h_SKMEL28_untreated_R3,
                countData.72h_SKMEL28_vemurafenib_R1, countData.72h_SKMEL28_vemurafenib_R2, countData.72h_SKMEL28_vemurafenib_R3,
                Normalized_MFE, Cp_37C) 


# Identify the 75th percentile of the combined MFE distributions of control and drug-treated samples
quantile_075_SKMEL28_dexseq_ANNOTATED_norm_TX_counts_MFE_df <- quantile(SKMEL28_dexseq_ANNOTATED_norm_TX_counts_MFE_df$Normalized_MFE, probs = c(0.75))

# Select all the transcripts with Normalized_MFE values above the 75th percentile
SKMEL28_dexseq_ANNOTATED_high_MFE_TX <- dplyr::filter(SKMEL28_dexseq_ANNOTATED_norm_TX_counts_MFE_df, Normalized_MFE > quantile_075_SKMEL28_dexseq_ANNOTATED_norm_TX_counts_MFE_df)

# Count the number of transcripts in the significantly up- and down-regulated pools that have MFE values above this threshold.
SKMEL28_dexseq_ANNOTATED_high_MFE_TX_in_UP <- dplyr::filter(SKMEL28_dexseq_ANNOTATED_UP, ensembl_transcript_id_version %in% SKMEL28_dexseq_ANNOTATED_high_MFE_TX$ensembl_transcript_id_version)
SKMEL28_dexseq_ANNOTATED_high_MFE_TX_in_DOWN <- dplyr::filter(SKMEL28_dexseq_ANNOTATED_DOWN, ensembl_transcript_id_version %in% SKMEL28_dexseq_ANNOTATED_high_MFE_TX$ensembl_transcript_id_version)


# Perform a Fisher’s exact test to compare the proportions of high MFE transcripts between the two groups.
SKMEL28_dexseq_ANNOTATED_total_UP <- nrow(SKMEL28_dexseq_ANNOTATED_UP)
SKMEL28_dexseq_ANNOTATED_total_DOWN <- nrow(SKMEL28_dexseq_ANNOTATED_DOWN)
SKMEL28_dexseq_ANNOTATED_count_high_MFE_UP <- nrow(SKMEL28_dexseq_ANNOTATED_high_MFE_TX_in_UP)
SKMEL28_dexseq_ANNOTATED_count_high_MFE_DOWN <- nrow(SKMEL28_dexseq_ANNOTATED_high_MFE_TX_in_DOWN)

SKMEL28_dexseq_ANNOTATED_contingency_table <- matrix(c(SKMEL28_dexseq_ANNOTATED_count_high_MFE_UP, SKMEL28_dexseq_ANNOTATED_total_UP - SKMEL28_dexseq_ANNOTATED_count_high_MFE_UP,
                                                     SKMEL28_dexseq_ANNOTATED_count_high_MFE_DOWN, SKMEL28_dexseq_ANNOTATED_total_DOWN - SKMEL28_dexseq_ANNOTATED_count_high_MFE_DOWN),
                                                   nrow = 2, byrow = TRUE,
                                                   dimnames = list(UP = c("UP-regulated", "DOWN-regulated"),
                                                                   DOWN = c("High MFE TX", "Not High MFE TX")))

SKMEL28_dexseq_ANNOTATED_contingency_table_df <- as.data.frame(SKMEL28_dexseq_ANNOTATED_contingency_table)

SKMEL28_dexseq_ANNOTATED_fisher_result <- fisher.test(SKMEL28_dexseq_ANNOTATED_contingency_table)

cat("SKMEL28_dexseq_ANNOTATED P-value from Fisher test:", SKMEL28_dexseq_ANNOTATED_fisher_result$p.value, "\n")


# Look at the fractions
SKMEL28_dexseq_ANNOTATED_high_MFE_in_UP_fraction <- SKMEL28_dexseq_ANNOTATED_count_high_MFE_UP / SKMEL28_dexseq_ANNOTATED_total_UP
SKMEL28_dexseq_ANNOTATED_high_MFE_in_DOWN_fraction <- SKMEL28_dexseq_ANNOTATED_count_high_MFE_DOWN / SKMEL28_dexseq_ANNOTATED_total_DOWN

# Output the results
cat("Fraction of SKMEL28_dexseq_ANNOTATED high MFE TX in UP-regulated genes:", SKMEL28_dexseq_ANNOTATED_high_MFE_in_UP_fraction, "\n") # This value goes in the summary table
cat("Fraction of SKMEL28_dexseq_ANNOTATED high MFE TX in DOWN-regulated genes:", SKMEL28_dexseq_ANNOTATED_high_MFE_in_DOWN_fraction, "\n") # This value goes in the summary table
cat("Mean of SKMEL28_dexseq_ANNOTATED high MFE TX in UP-regulated genes:", mean(SKMEL28_dexseq_ANNOTATED_high_MFE_TX_in_UP$Normalized_MFE), "\n")
cat("Mean of SKMEL28_dexseq_ANNOTATED high MFE TX in DOWN-regulated genes:", mean(SKMEL28_dexseq_ANNOTATED_high_MFE_TX_in_DOWN$Normalized_MFE), "\n")
cat("Mode of SKMEL28_dexseq_ANNOTATED high MFE TX in UP-regulated genes:", calculate_mode(SKMEL28_dexseq_ANNOTATED_high_MFE_TX_in_UP$Normalized_MFE), "\n")
cat("Mode of SKMEL28_dexseq_ANNOTATED high MFE TX in DOWN-regulated genes:", calculate_mode(SKMEL28_dexseq_ANNOTATED_high_MFE_TX_in_DOWN$Normalized_MFE), "\n")





# Perform a two-proportion z-test to compare proportions
# Calculate the proportions and total counts
SKMEL28_dexseq_ANNOTATED_high_MFE_in_UP_proportion <- SKMEL28_dexseq_ANNOTATED_high_MFE_in_UP_fraction
SKMEL28_dexseq_ANNOTATED_high_MFE_in_DOWN_proportion <- SKMEL28_dexseq_ANNOTATED_high_MFE_in_DOWN_fraction

# Combined proportion for pooled standard error
combined_proportion <- (SKMEL28_dexseq_ANNOTATED_count_high_MFE_UP + SKMEL28_dexseq_ANNOTATED_count_high_MFE_DOWN) /
  (SKMEL28_dexseq_ANNOTATED_total_UP + SKMEL28_dexseq_ANNOTATED_total_DOWN)

# Standard error of the difference in proportions
standard_error <- sqrt(
  combined_proportion * (1 - combined_proportion) * 
    (1 / SKMEL28_dexseq_ANNOTATED_total_UP + 1 / SKMEL28_dexseq_ANNOTATED_total_DOWN)
)

# Z-statistic
z_statistic <- (SKMEL28_dexseq_ANNOTATED_high_MFE_in_UP_proportion - SKMEL28_dexseq_ANNOTATED_high_MFE_in_DOWN_proportion) /
  standard_error

# Two-tailed p-value
SKMEL28_dexseq_ANNOTATED_z_test_MFE_p_value <- 2 * (1 - pnorm(abs(z_statistic)))




# Show data as violin plot
# Add a new column to indicate if the value is above the 75th percentile
SKMEL28_dexseq_ANNOTATED$above_quantile_075 <- SKMEL28_dexseq_ANNOTATED$Normalized_MFE > quantile_075_SKMEL28_dexseq_ANNOTATED_norm_TX_counts_MFE_df

# Generate violin plot with conditional jitter coloring
ggplot(SKMEL28_dexseq_ANNOTATED, aes(x = factor(log2fold_treated_untreated > 0), y = Normalized_MFE)) +
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
    title = "SKMEL28_dexseq_ANNOTATED_Log2_FC_adj.P<0.05"
  ) +
  scale_x_discrete(labels = c("log2FoldChange < 0", "log2FoldChange > 0"))












## Fraction of total counts analysis ####
# Determine the impact of these changes in terms of fraction of total counts:
# Calculate the fraction of total counts of high MFE transcripts in control and drug-treated groups. 
# This is done by summing the counts of the high MFE transcripts and dividing by the total counts in each group.
# Fisher’s exact test to compare the proportions of high MFE transcripts counts between the two groups.

SKMEL28_dexseq_ANNOTATED_whole_transcriptome <- read.delim(file="path_to/72h_SKMEL28_vemurafenib.tx_counts.both_conditions.treated_untreated_dexseq.tab", header = TRUE, sep = "\t", quote = "") %>%
  dplyr::rename(ensembl_transcript_id_version=featureID) %>%
  dplyr::inner_join(Table_ANNOTATED, by="ensembl_transcript_id_version") %>%
  dplyr::select(ensembl_transcript_id_version, ensembl_transcript_id, ensembl_gene_id, external_gene_name, transcript_is_canonical, Transcript_length,
                countData.72h_SKMEL28_vemurafenib_R1, countData.72h_SKMEL28_vemurafenib_R2, countData.72h_SKMEL28_vemurafenib_R3,
                countData.72h_SKMEL28_untreated_R1, countData.72h_SKMEL28_untreated_R2, countData.72h_SKMEL28_untreated_R3,
                Normalized_MFE, Cp_37C) 

# Convert from wide to long format
SKMEL28_dexseq_ANNOTATED_long_data <- SKMEL28_dexseq_ANNOTATED_whole_transcriptome %>%
  pivot_longer(
    cols = starts_with("countData.72h"),      # Select all columns starting with "countData.72h"
    names_to = "replicate",         # Name of the new column for replicate identifiers
    values_to = "counts"            # Name of the new column for the counts
  )





# Add a condition column based on the replicate
SKMEL28_dexseq_ANNOTATED_long_data <- SKMEL28_dexseq_ANNOTATED_long_data %>%
  dplyr::mutate(
    condition = dplyr::if_else(
      replicate %in% c(
        "countData.72h_SKMEL28_untreated_R1",
        "countData.72h_SKMEL28_untreated_R2",
        "countData.72h_SKMEL28_untreated_R3"
      ),
      "Control", "vemurafenib"
    ),
    condition = factor(condition, levels = c("Control", "vemurafenib"))
  )

# Aggregate to average counts per transcript x condition
SKMEL28_dexseq_ANNOTATED_norm_TX_aggregated_data_avg_counts <- SKMEL28_dexseq_ANNOTATED_long_data %>%
  dplyr::group_by(ensembl_transcript_id_version, condition, Normalized_MFE) %>%
  dplyr::summarise(
    avg_counts = mean(counts, na.rm = TRUE),
    .groups = "drop"
  )






# High MFE (>= 75th percentile)
# Test whether cells exposed to vemurafenib have a higher fraction of transcripts in the 75th percentile of MFE values compared to control cells.
# 75th percentile threshold
mfe_threshold_SKMEL28 <- quantile(
  SKMEL28_dexseq_ANNOTATED_norm_TX_aggregated_data_avg_counts$Normalized_MFE,
  0.75, na.rm = TRUE
)

# Label high MFE
df_SKMEL28_high <- SKMEL28_dexseq_ANNOTATED_norm_TX_aggregated_data_avg_counts %>%
  dplyr::mutate(MFE_class = ifelse(Normalized_MFE >= mfe_threshold_SKMEL28, "High_MFE", "Other"))

# Aggregate counts per condition x class
counts_table_SKMEL28_high <- df_SKMEL28_high %>%
  dplyr::group_by(condition, MFE_class) %>%
  dplyr::summarise(total_counts = sum(avg_counts, na.rm = TRUE), .groups = "drop")

# 2x2 contingency (row1=vemurafenib, row2=Control; col1=High_MFE, col2=Other)
contingency_SKMEL28_high <- counts_table_SKMEL28_high %>%
  tidyr::pivot_wider(names_from = MFE_class, values_from = total_counts, values_fill = 0) %>%
  dplyr::arrange(factor(condition, levels = c("vemurafenib", "Control")))

fisher_matrix_SKMEL28_high <- as.matrix(contingency_SKMEL28_high[, c("High_MFE", "Other")])
rownames(fisher_matrix_SKMEL28_high) <- contingency_SKMEL28_high$condition

# One-sided Fisher: vemurafenib > Control (higher High_MFE fraction)
fisher_high_SKMEL28_result <- fisher.test(fisher_matrix_SKMEL28_high, alternative = "greater")
fisher_high_SKMEL28_result

# Report High_MFE fractions per condition
high_mfe_fractions_SKMEL28 <- counts_table_SKMEL28_high %>%
  dplyr::group_by(condition) %>%
  dplyr::mutate(fraction = total_counts / sum(total_counts)) %>%
  dplyr::filter(MFE_class == "High_MFE")
high_mfe_fractions_SKMEL28




SKMEL28_dexseq_ANNOTATED_vemurafenib_high_MFE_TX_counts_fraction <-
  high_mfe_fractions_SKMEL28 %>%
  dplyr::filter(condition == "vemurafenib") %>%
  dplyr::pull(fraction) %>%
  as.numeric()

SKMEL28_dexseq_ANNOTATED_control_high_MFE_TX_counts_fraction <-
  high_mfe_fractions_SKMEL28 %>%
  dplyr::filter(condition == "Control") %>%
  dplyr::pull(fraction) %>%
  as.numeric()








# Low MFE (<= 25th percentile)
# Tests whether cells exposed to vemurafenib have a lower fraction of transcripts in the 25th percentile of MFE values compared to control cells.
# 25th percentile threshold
mfe_threshold_25_SKMEL28 <- quantile(
  SKMEL28_dexseq_ANNOTATED_norm_TX_aggregated_data_avg_counts$Normalized_MFE,
  0.25, na.rm = TRUE
)

# Label low MFE
df_SKMEL28_low <- SKMEL28_dexseq_ANNOTATED_norm_TX_aggregated_data_avg_counts %>%
  dplyr::mutate(MFE_class = ifelse(Normalized_MFE <= mfe_threshold_25_SKMEL28, "Low_MFE", "Other"))

# Aggregate counts per condition x class
counts_table_SKMEL28_low <- df_SKMEL28_low %>%
  dplyr::group_by(condition, MFE_class) %>%
  dplyr::summarise(total_counts = sum(avg_counts, na.rm = TRUE), .groups = "drop")

# 2x2 contingency (row1=Control, row2=vemurafenib; col1=Low_MFE, col2=Other)
contingency_SKMEL28_low <- counts_table_SKMEL28_low %>%
  tidyr::pivot_wider(names_from = MFE_class, values_from = total_counts, values_fill = 0) %>%
  dplyr::arrange(factor(condition, levels = c("Control", "vemurafenib")))

fisher_matrix_SKMEL28_low <- as.matrix(contingency_SKMEL28_low[, c("Low_MFE", "Other")])
rownames(fisher_matrix_SKMEL28_low) <- contingency_SKMEL28_low$condition

# One-sided Fisher: Control > vemurafenib (higher Low_MFE fraction)
fisher_low_SKMEL28_result <- fisher.test(fisher_matrix_SKMEL28_low, alternative = "greater")
fisher_low_SKMEL28_result

# Report Low_MFE fractions per condition
low_mfe_fractions_SKMEL28 <- counts_table_SKMEL28_low %>%
  dplyr::group_by(condition) %>%
  dplyr::mutate(fraction = total_counts / sum(total_counts)) %>%
  dplyr::filter(MFE_class == "Low_MFE")
low_mfe_fractions_SKMEL28













## SKMEL28_dexseq_ANNOTATED heat capacity changes ####

# This is calculated on the total n. of transcripts with a significant log2 fold change
# Generate a scatter plot Log2_FC vs heat capacity
ggplot(SKMEL28_dexseq_ANNOTATED, aes(y=Cp_37C, x=log2fold_treated_untreated)) +
  geom_point(size=1.5, alpha=0.1)+
  theme_classic() +
  theme(legend.position="none") +
  geom_vline(xintercept = 0, size=0.3) +
  xlim(-12, 12) +
  ylim(0,400) +
  theme(axis.title = element_text(size = 20)) +
  theme(plot.title = element_text(size = 15)) +
  labs(y="Heat capacity (J/C)",x="Log2_FC_expression",title="SKMEL28_dexseq_ANNOTATED_Log2_FC_adj.P<0.05")


# Show the data with a box plot
ggplot(SKMEL28_dexseq_ANNOTATED, aes(x = factor(log2fold_treated_untreated > 0), y = Cp_37C)) +
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
  labs(y="Heat capacity (J/C)",x="",title="SKMEL28_dexseq_ANNOTATED_Log2_FC_adj.P<0.05 heat_capacity all transcripts") +
  scale_x_discrete(labels = c("log2FoldChange < 0", "log2FoldChange > 0"))

Mean_heat_capacity_SKMEL28_dexseq_ANNOTATED_UP <- mean(SKMEL28_dexseq_ANNOTATED_UP$Cp_37C)
Mean_heat_capacity_SKMEL28_dexseq_ANNOTATED_DOWN <- mean(SKMEL28_dexseq_ANNOTATED_DOWN$Cp_37C)
TTEST_heat_capacity_SKMEL28_dexseq_ANNOTATED_all_UP_vs_DOWN <- t.test(SKMEL28_dexseq_ANNOTATED_DOWN$Cp_37C, SKMEL28_dexseq_ANNOTATED_UP$Cp_37C, paired=FALSE)
print(TTEST_heat_capacity_SKMEL28_dexseq_ANNOTATED_all_UP_vs_DOWN)
# Add K-S test to look at the differences between the up-regukated and down-regulated genes
ks_test_SKMEL28_dexseq_ANNOTATED_all_UP_vs_DOWN_heat_capacity <- ks.test(SKMEL28_dexseq_ANNOTATED_DOWN$Cp_37C, SKMEL28_dexseq_ANNOTATED_UP$Cp_37C)
print(ks_test_SKMEL28_dexseq_ANNOTATED_all_UP_vs_DOWN_heat_capacity)




# Calculate the changes of heat capacity only in the high MFE transcripts

# Contribution of the UP or DOWN-regulated high MFE transcripts to the heat capacity
# Look at the heat capacity of the UP-regulated high MFE transcripts
SKMEL28_dexseq_ANNOTATED_high_MFE_TX_in_UP_heat_capacity <- mean(SKMEL28_dexseq_ANNOTATED_high_MFE_TX_in_UP$Cp_37C)

# Look at the heat capacity of the DOWN-regulated high MFE transcripts
SKMEL28_dexseq_ANNOTATED_high_MFE_TX_in_DOWN_heat_capacity <- mean(SKMEL28_dexseq_ANNOTATED_high_MFE_TX_in_DOWN$Cp_37C)

# Heat capacity gain in UP-regulated high MFE transcripts vs DOWN-regulated high MFE transcripts
SKMEL28_dexseq_ANNOTATED_heat_capacity_gain_UP_vs_DOWN <- SKMEL28_dexseq_ANNOTATED_high_MFE_TX_in_UP_heat_capacity - SKMEL28_dexseq_ANNOTATED_high_MFE_TX_in_DOWN_heat_capacity
cat("SKMEL28_dexseq_ANNOTATED net heat capacity UP vs DOWN high MFE TX:", SKMEL28_dexseq_ANNOTATED_heat_capacity_gain_UP_vs_DOWN, "\n")


# Extract specific heat data
SKMEL28_dexseq_ANNOTATED_Cp_37C_UP <- SKMEL28_dexseq_ANNOTATED_high_MFE_TX_in_UP$Cp_37C
SKMEL28_dexseq_ANNOTATED_Cp_37C_DOWN <- SKMEL28_dexseq_ANNOTATED_high_MFE_TX_in_DOWN$Cp_37C

# Combine data into a single data frame for easier plotting
SKMEL28_dexseq_ANNOTATED_Cp_37C_data <- data.frame(
  SKMEL28_dexseq_ANNOTATED_Cp_37C = c(SKMEL28_dexseq_ANNOTATED_Cp_37C_UP, SKMEL28_dexseq_ANNOTATED_Cp_37C_DOWN),
  regulation = factor(c(rep("UP", length(SKMEL28_dexseq_ANNOTATED_Cp_37C_UP)), rep("DOWN", length(SKMEL28_dexseq_ANNOTATED_Cp_37C_DOWN))))
)

# Generate box plot
ggplot(SKMEL28_dexseq_ANNOTATED_Cp_37C_data, aes(x = regulation, y = SKMEL28_dexseq_ANNOTATED_Cp_37C, fill = regulation)) +
  geom_boxplot(fill = "lightgrey", alpha = 0.7) +
  labs(title = "SKMEL28_dexseq_ANNOTATED specific Heat of High MFE Transcripts",
       x = "",
       y = "Specific Heat (J/C)") +
  scale_x_discrete(labels = c("log2_FC < 0", "log2_FC > 0")) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 16),
        axis.title = element_text(size = 25))

Mean_heat_capacity_SKMEL28_dexseq_ANNOTATED_high_MFE_TX_UP <- mean(SKMEL28_dexseq_ANNOTATED_high_MFE_TX_in_UP$Cp_37C)
Mean_heat_capacity_SKMEL28_dexseq_ANNOTATED_high_MFE_TX_DOWN <- mean(SKMEL28_dexseq_ANNOTATED_high_MFE_TX_in_DOWN$Cp_37C)
TTEST_heat_capacity_SKMEL28_dexseq_ANNOTATED_all_high_MFE_TX_UP_vs_DOWN <- t.test(SKMEL28_dexseq_ANNOTATED_high_MFE_TX_in_DOWN$Cp_37C, SKMEL28_dexseq_ANNOTATED_high_MFE_TX_in_UP$Cp_37C, paired=FALSE)
print(TTEST_heat_capacity_SKMEL28_dexseq_ANNOTATED_all_high_MFE_TX_UP_vs_DOWN)

# Add K-S test to look at the differences between the up-regukated and down-regulated genes
ks_test_heat_capacity_SKMEL28_dexseq_ANNOTATED_all_high_MFE_TX_UP_vs_DOWN <- ks.test(SKMEL28_dexseq_ANNOTATED_high_MFE_TX_in_DOWN$Cp_37C, SKMEL28_dexseq_ANNOTATED_high_MFE_TX_in_UP$Cp_37C)
print(ks_test_heat_capacity_SKMEL28_dexseq_ANNOTATED_all_high_MFE_TX_UP_vs_DOWN)











# Analysis cell proliferation Nanopore data SKMEL28 VEM ####

SKMEL28_VEM_DESeq2 <- read.delim(file="path_to/72h_SKMEL28_vemurafenib.tx_counts.both_conditions.treated_untreated_deseq2.tab", header = TRUE, sep = "\t", quote = "")



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
CURATED_PATHWAYS.output <- c("SKMEL28_VEM_MTOR_POS", "SKMEL28_VEM_MTOR_NEG", "SKMEL28_VEM_Proliferation_Science", "SKMEL28_VEM_Meta_PCNA", "SKMEL28_VEM_Proliferation_GO", "SKMEL28_VEM_Proliferation_STOP", "SKMEL28_VEM_G1_phase", "SKMEL28_VEM_Proliferation_Carey", "SKMEL28_VEM_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019", "SKMEL28_VEM_Essential_S_Phase_Viner_Breuer_Cell_Division_2019", "SKMEL28_VEM_HR", "SKMEL28_VEM_Fanconi_anemia", "SKMEL28_VEM_HF_POL")

# Assign the names of the gene sets to each entry in the list
names(CURATED_PATHWAYS.input) <- CURATED_PATHWAYS.output

# Apply the filtering function to create a list of data frames
Curated_test_df <- lapply(CURATED_PATHWAYS.input, function(x) filtering(x))

# creating a vector with desired names of filtered data frames
Curated_dfs.output <- c("df_SKMEL28_VEM_MTOR_POS", "df_SKMEL28_VEM_MTOR_NEG", "df_SKMEL28_VEM_Proliferation_Science", "df_SKMEL28_VEM_Meta_PCNA", "df_SKMEL28_VEM_Proliferation_GO", "df_SKMEL28_VEM_Proliferation_STOP", "df_SKMEL28_VEM_G1_phase", "df_SKMEL28_VEM_Proliferation_Carey", "df_SKMEL28_VEM_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019", "df_SKMEL28_VEM_Essential_S_Phase_Viner_Breuer_Cell_Division_2019", "df_SKMEL28_VEM_HR", "df_SKMEL28_VEM_Fanconi_anemia", "df_SKMEL28_VEM_HF_POL")

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
Curated_input_dfs <- list(df_SKMEL28_VEM_MTOR_POS, df_SKMEL28_VEM_MTOR_NEG, df_SKMEL28_VEM_Proliferation_Science, df_SKMEL28_VEM_Meta_PCNA, df_SKMEL28_VEM_Proliferation_GO, df_SKMEL28_VEM_Proliferation_STOP, df_SKMEL28_VEM_G1_phase, df_SKMEL28_VEM_Proliferation_Carey, df_SKMEL28_VEM_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019, df_SKMEL28_VEM_Essential_S_Phase_Viner_Breuer_Cell_Division_2019, df_SKMEL28_VEM_HR, df_SKMEL28_VEM_Fanconi_anemia, df_SKMEL28_VEM_HF_POL)

# Apply the mean function to create a list of means
Curated_mean_test_df <- lapply(Curated_input_dfs, function(x) Mean_curated(x))

# creating a vector with desired names of the files containing the mean of each gene set log2FoldChange
Curated_Mean_SKMEL28_VEM_dfs.output <- c("Mean_SKMEL28_VEM_MTOR_POS", "Mean_SKMEL28_VEM_MTOR_NEG", "Mean_SKMEL28_VEM_Proliferation_Science", "Mean_SKMEL28_VEM_Meta_PCNA", "Mean_SKMEL28_VEM_Proliferation_GO", "Mean_SKMEL28_VEM_Proliferation_STOP", "Mean_SKMEL28_VEM_G1_phase", "Mean_SKMEL28_VEM_Proliferation_Carey", "Mean_SKMEL28_VEM_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019", "Mean_SKMEL28_VEM_Essential_S_Phase_Viner_Breuer_Cell_Division_2019", "Mean_SKMEL28_VEM_HR", "Mean_SKMEL28_VEM_Fanconi_anemia", "Mean_SKMEL28_VEM_HF_POL")

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
Curated_STDEV_SKMEL28_VEM_dfs.output <- c("STDEV_SKMEL28_VEM_MTOR_POS", "STDEV_SKMEL28_VEM_MTOR_NEG", "STDEV_SKMEL28_VEM_Proliferation_Science", "STDEV_SKMEL28_VEM_Meta_PCNA", "STDEV_SKMEL28_VEM_Proliferation_GO", "STDEV_SKMEL28_VEM_Proliferation_STOP", "STDEV_SKMEL28_VEM_G1_phase", "STDEV_SKMEL28_VEM_Proliferation_Carey", "STDEV_SKMEL28_VEM_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019", "STDEV_SKMEL28_VEM_Essential_S_Phase_Viner_Breuer_Cell_Division_2019", "STDEV_SKMEL28_VEM_HR", "STDEV_SKMEL28_VEM_Fanconi_anemia", "STDEV_SKMEL28_VEM_HF_POL")
# loop to assign data.frames to the desired variable name
for (i in (seq_along(Curated_STDEV_SKMEL28_VEM_dfs.output))) {
  do.call("<-",list(Curated_STDEV_SKMEL28_VEM_dfs.output[i], Curated_stdev_test_df[[i]]))
}




# Generate the data frame for the plot

# Generate the data frame for the plot
df_SKMEL28_VEM_proliferation <- data.frame(mean=c(Mean_SKMEL28_VEM_MTOR_POS, Mean_SKMEL28_VEM_MTOR_NEG, Mean_SKMEL28_VEM_Proliferation_Science, Mean_SKMEL28_VEM_Meta_PCNA, Mean_SKMEL28_VEM_Proliferation_GO, Mean_SKMEL28_VEM_Proliferation_STOP, Mean_SKMEL28_VEM_G1_phase, Mean_SKMEL28_VEM_Proliferation_Carey, Mean_SKMEL28_VEM_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019, Mean_SKMEL28_VEM_Essential_S_Phase_Viner_Breuer_Cell_Division_2019, Mean_SKMEL28_VEM_HR, Mean_SKMEL28_VEM_Fanconi_anemia, Mean_SKMEL28_VEM_HF_POL),
                                          sd=c(STDEV_SKMEL28_VEM_MTOR_POS, STDEV_SKMEL28_VEM_MTOR_NEG, STDEV_SKMEL28_VEM_Proliferation_Science, STDEV_SKMEL28_VEM_Meta_PCNA, STDEV_SKMEL28_VEM_Proliferation_GO, STDEV_SKMEL28_VEM_Proliferation_STOP, STDEV_SKMEL28_VEM_G1_phase, STDEV_SKMEL28_VEM_Proliferation_Carey, STDEV_SKMEL28_VEM_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019, STDEV_SKMEL28_VEM_Essential_S_Phase_Viner_Breuer_Cell_Division_2019, STDEV_SKMEL28_VEM_HR, STDEV_SKMEL28_VEM_Fanconi_anemia, STDEV_SKMEL28_VEM_HF_POL),
                                          Category=as.factor(c("SKMEL28_VEM_MTOR_POS", "SKMEL28_VEM_MTOR_NEG", "SKMEL28_VEM_Proliferation_Science", "SKMEL28_VEM_Meta_PCNA", "SKMEL28_VEM_Proliferation_GO", "SKMEL28_VEM_Proliferation_STOP", "SKMEL28_VEM_G1_phase", "SKMEL28_VEM_Proliferation_Carey", "SKMEL28_VEM_Essential_cell_cycle", "SKMEL28_VEM_Essential_S_Phase", "SKMEL28_VEM_HR", "SKMEL28_VEM_FA", "SKMEL28_VEM_HF_POL")))
df_SKMEL28_VEM_proliferation




library(ggplot2)
ggplot(df_SKMEL28_VEM_proliferation, aes(x=reorder(Category, -mean), y=mean)) +
  theme_bw() +
  geom_point(position=position_dodge(width = 1), stat="identity",colour='black') +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.4) +
  theme(axis.text.x = element_text(size=12, angle=45,hjust=0.95,vjust=0.98)) +
  theme(axis.text.y = element_text(size=14)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 10))) +
  geom_point(colour = "black", size = 3) +
  xlab("") + ylab("Relative log2 fold change") +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  ggtitle("Analysis cell proliferation 778 NUTLIN NANOPORE 72h vs T0") +
  scale_x_discrete(limits = c("SKMEL28_VEM_Proliferation_Science", "SKMEL28_VEM_Meta_PCNA", "SKMEL28_VEM_Proliferation_GO", "SKMEL28_VEM_Proliferation_STOP", "SKMEL28_VEM_Proliferation_Carey", "SKMEL28_VEM_G1_phase", "SKMEL28_VEM_MTOR_POS", "SKMEL28_VEM_MTOR_NEG", "SKMEL28_VEM_Essential_cell_cycle", "SKMEL28_VEM_Essential_S_Phase", "SKMEL28_VEM_HR", "SKMEL28_VEM_FA", "SKMEL28_VEM_HF_POL")) +
  theme(plot.title = element_text(color="black", size=14, face="bold"))


# Merge all the data frames to generate box plots
# Add a column with the specific name of the gene set
df_SKMEL28_VEM_MTOR_POS$Gene_set <- "SKMEL28_VEM_MTOR_POS"
df_SKMEL28_VEM_MTOR_NEG$Gene_set <- "SKMEL28_VEM_MTOR_NEG"
df_SKMEL28_VEM_Proliferation_Science$Gene_set <- "SKMEL28_VEM_Proliferation_Science"
df_SKMEL28_VEM_Meta_PCNA$Gene_set <- "SKMEL28_VEM_Meta_PCNA"
df_SKMEL28_VEM_Proliferation_GO$Gene_set <- "SKMEL28_VEM_Proliferation_GO"
df_SKMEL28_VEM_Proliferation_STOP$Gene_set <- "SKMEL28_VEM_Proliferation_STOP"
df_SKMEL28_VEM_G1_phase$Gene_set <- "SKMEL28_VEM_G1_phase"
df_SKMEL28_VEM_Proliferation_Carey$Gene_set <- "SKMEL28_VEM_Proliferation_Carey"
df_SKMEL28_VEM_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019$Gene_set <- "SKMEL28_VEM_Essential_cell_cycle"
df_SKMEL28_VEM_Essential_S_Phase_Viner_Breuer_Cell_Division_2019$Gene_set <- "SKMEL28_VEM_Essential_S_Phase"
df_SKMEL28_VEM_HR$Gene_set <- "SKMEL28_VEM_HR"
df_SKMEL28_VEM_Fanconi_anemia$Gene_set <- "SKMEL28_VEM_Fanconi_anemia"
df_SKMEL28_VEM_HF_POL$Gene_set <- "SKMEL28_VEM_HF_POL"



# Combine the data frames
combined_df <- rbind(
  df_SKMEL28_VEM_MTOR_POS,
  df_SKMEL28_VEM_MTOR_NEG,
  df_SKMEL28_VEM_Proliferation_Science,
  df_SKMEL28_VEM_Meta_PCNA,
  df_SKMEL28_VEM_Proliferation_GO,
  df_SKMEL28_VEM_Proliferation_STOP,
  df_SKMEL28_VEM_G1_phase,
  df_SKMEL28_VEM_Proliferation_Carey,
  df_SKMEL28_VEM_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019,
  df_SKMEL28_VEM_Essential_S_Phase_Viner_Breuer_Cell_Division_2019,
  df_SKMEL28_VEM_HR,
  df_SKMEL28_VEM_Fanconi_anemia,
  df_SKMEL28_VEM_HF_POL
)

# Ensure the Gene_set column is a factor with the correct order
combined_df$Gene_set <- factor(combined_df$Gene_set, levels = c("SKMEL28_VEM_Proliferation_Science", "SKMEL28_VEM_Meta_PCNA", "SKMEL28_VEM_Proliferation_GO", "SKMEL28_VEM_Proliferation_STOP", "SKMEL28_VEM_Proliferation_Carey", "SKMEL28_VEM_G1_phase", "SKMEL28_VEM_MTOR_POS", "SKMEL28_VEM_MTOR_NEG", "SKMEL28_VEM_Essential_cell_cycle", "SKMEL28_VEM_Essential_S_Phase", "SKMEL28_VEM_HR", "SKMEL28_VEM_Fanconi_anemia", "SKMEL28_VEM_HF_POL"))


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
  ggtitle("Analysis cell proliferation SKMEL28 vemurafenib 72h") +
  scale_x_discrete(limits = levels(combined_df$Gene_set)) +
  theme(plot.title = element_text(color="black", size=14, face="bold")) 







Mean_df_MTOR_POS_SKMEL28_VEM <- mean(df_SKMEL28_VEM_MTOR_POS$log2FoldChange)
Mean_df_MTOR_NEG_SKMEL28_VEM <- mean(df_SKMEL28_VEM_MTOR_NEG$log2FoldChange)
Mean_df_Proliferation_Science_SKMEL28_VEM <- mean(df_SKMEL28_VEM_Proliferation_Science$log2FoldChange)
Mean_df_Meta_PCNA_SKMEL28_VEM <- mean(df_SKMEL28_VEM_Meta_PCNA$log2FoldChange)
Mean_df_Proliferation_GO_SKMEL28_VEM <- mean(df_SKMEL28_VEM_Proliferation_GO$log2FoldChange)
Mean_df_Proliferation_STOP_SKMEL28_VEM <- mean(df_SKMEL28_VEM_Proliferation_STOP$log2FoldChange)
Mean_df_G1_phase_SKMEL28_VEM <- mean(df_SKMEL28_VEM_G1_phase$log2FoldChange)
Mean_df_Proliferation_Carey_SKMEL28_VEM <- mean(df_SKMEL28_VEM_Proliferation_Carey$log2FoldChange)
Mean_df_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019_SKMEL28_VEM <- mean(df_SKMEL28_VEM_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019$log2FoldChange)
Mean_df_Essential_S_Phase_Viner_Breuer_Cell_Division_2019_SKMEL28_VEM <- mean(df_SKMEL28_VEM_Essential_S_Phase_Viner_Breuer_Cell_Division_2019$log2FoldChange)
Mean_df_HR_SKMEL28_VEM <- mean(df_SKMEL28_VEM_HR$log2FoldChange)
Mean_df_FA_SKMEL28_VEM <- mean(df_SKMEL28_VEM_Fanconi_anemia$log2FoldChange)
Mean_df_HF_POL_SKMEL28_VEM <- mean(df_SKMEL28_VEM_HF_POL$log2FoldChange)

# Perform the Wilcoxon signed-rank test
wilcox_test_MTOR_POS_SKMEL28_VEM <- wilcox.test(df_SKMEL28_VEM_MTOR_POS$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_MTOR_NEG_SKMEL28_VEM <- wilcox.test(df_SKMEL28_VEM_MTOR_NEG$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_Proliferation_Science_SKMEL28_VEM <- wilcox.test(df_SKMEL28_VEM_Proliferation_Science$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_Meta_PCNA_SKMEL28_VEM <- wilcox.test(df_SKMEL28_VEM_Meta_PCNA$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_Proliferation_GO_SKMEL28_VEM <- wilcox.test(df_SKMEL28_VEM_Proliferation_GO$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_Proliferation_STOP_SKMEL28_VEM <- wilcox.test(df_SKMEL28_VEM_Proliferation_STOP$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_G1_phase_SKMEL28_VEM <- wilcox.test(df_SKMEL28_VEM_G1_phase$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_Proliferation_Carey_SKMEL28_VEM <- wilcox.test(df_SKMEL28_VEM_Proliferation_Carey$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019_SKMEL28_VEM <- wilcox.test(df_SKMEL28_VEM_Essential_cell_cycle_Viner_Breuer_Cell_Division_2019$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_Essential_S_Phase_Viner_Breuer_Cell_Division_2019_SKMEL28_VEM <- wilcox.test(df_SKMEL28_VEM_Essential_S_Phase_Viner_Breuer_Cell_Division_2019$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_HR_SKMEL28_VEM <- wilcox.test(df_SKMEL28_VEM_HR$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_FA_SKMEL28_VEM <- wilcox.test(df_SKMEL28_VEM_Fanconi_anemia$log2FoldChange, mu = 0, alternative = "less")
wilcox_test_HF_POL_SKMEL28_VEM <- wilcox.test(df_SKMEL28_VEM_HF_POL$log2FoldChange, mu = 0, alternative = "less")















