suppressPackageStartupMessages({
  library(data.table)
  library(readr)
  library(tibble)
  library(ggplot2)
})

joined_file <- "path_to/HEK293T_icSHAPE_ViennaRNA_joined_per_position_table.tsv.gz"
output_dir <- "path_to/HEK293T_icSHAPE_ViennaRNA_followup_analysis"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

min_paired_positions <- 10L
min_unpaired_positions <- 10L

message2 <- function(...) cat(sprintf(...), "\n")

stop_if_missing <- function(x, label) {
  if (!file.exists(x)) stop(label, " not found: ", x, call. = FALSE)
}

safe_wilcox_two_group <- function(score_vec, group_vec) {
  if (length(unique(group_vec[!is.na(group_vec)])) < 2) return(NA_real_)
  out <- tryCatch(
    wilcox.test(score_vec ~ group_vec),
    error = function(e) NULL
  )
  if (is.null(out)) return(NA_real_)
  out$p.value
}

safe_median <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  median(x, na.rm = TRUE)
}

safe_mean <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

stop_if_missing(joined_file, "Joined icSHAPE-ViennaRNA per-position table")

joined_dt <- fread(joined_file)

message2("Rows in joined table: %s", format(nrow(joined_dt), big.mark = ","))
message2("Unique transcripts in joined table: %s", format(uniqueN(joined_dt$transcript_id_norm), big.mark = ","))
message2("Unique genes in joined table: %s", format(uniqueN(joined_dt$gene_id), big.mark = ","))

if (!all(c("score", "vienna_state", "transcript_id_norm", "tx_position") %in% names(joined_dt))) {
  stop("Joined table is missing required columns.", call. = FALSE)
}

global_state_summary <- joined_dt[, .(
  n = .N,
  mean_score = mean(score, na.rm = TRUE),
  median_score = median(score, na.rm = TRUE),
  sd_score = sd(score, na.rm = TRUE),
  q25 = quantile(score, 0.25, na.rm = TRUE),
  q75 = quantile(score, 0.75, na.rm = TRUE)
), by = vienna_state][order(vienna_state)]

global_wilcox_p <- safe_wilcox_two_group(joined_dt$score, joined_dt$vienna_state)

message2("Global score summary by ViennaRNA state:")
print(global_state_summary)
cat("\n")

message2("Global Wilcoxon p-value: %s", format(global_wilcox_p, scientific = TRUE))
cat("\n")

transcript_summary <- joined_dt[, .(
  gene_name = unique(gene_name)[1],
  transcript_name = unique(transcript_name)[1],
  Transcript_length = unique(Transcript_length)[1],
  MFE = unique(MFE)[1],
  n_positions = .N,
  n_paired = sum(vienna_state == "paired", na.rm = TRUE),
  n_unpaired = sum(vienna_state == "unpaired", na.rm = TRUE),
  mean_score_paired = safe_mean(score[vienna_state == "paired"]),
  mean_score_unpaired = safe_mean(score[vienna_state == "unpaired"]),
  median_score_paired = safe_median(score[vienna_state == "paired"]),
  median_score_unpaired = safe_median(score[vienna_state == "unpaired"]),
  delta_mean_unpaired_minus_paired =
    safe_mean(score[vienna_state == "unpaired"]) -
    safe_mean(score[vienna_state == "paired"]),
  delta_median_unpaired_minus_paired =
    safe_median(score[vienna_state == "unpaired"]) -
    safe_median(score[vienna_state == "paired"]),
  wilcox_p = safe_wilcox_two_group(score, vienna_state)
), by = .(transcript_id_norm, transcript_id)]

transcript_summary[, passes_coverage_filter :=
                     n_paired >= min_paired_positions & n_unpaired >= min_unpaired_positions
]

transcript_summary_filtered <- transcript_summary[passes_coverage_filter == TRUE]

message2("Rows in transcript summary: %s", format(nrow(transcript_summary), big.mark = ","))
message2("Rows passing coverage filter: %s", format(nrow(transcript_summary_filtered), big.mark = ","))
cat("\n")

direction_summary_all <- tibble(
  n_transcripts_total = nrow(transcript_summary),
  n_delta_mean_gt_0 = sum(transcript_summary$delta_mean_unpaired_minus_paired > 0, na.rm = TRUE),
  n_delta_mean_lt_0 = sum(transcript_summary$delta_mean_unpaired_minus_paired < 0, na.rm = TRUE),
  n_delta_mean_eq_0 = sum(transcript_summary$delta_mean_unpaired_minus_paired == 0, na.rm = TRUE),
  frac_delta_mean_gt_0 = mean(transcript_summary$delta_mean_unpaired_minus_paired > 0, na.rm = TRUE),
  n_delta_median_gt_0 = sum(transcript_summary$delta_median_unpaired_minus_paired > 0, na.rm = TRUE),
  n_delta_median_lt_0 = sum(transcript_summary$delta_median_unpaired_minus_paired < 0, na.rm = TRUE),
  n_delta_median_eq_0 = sum(transcript_summary$delta_median_unpaired_minus_paired == 0, na.rm = TRUE),
  frac_delta_median_gt_0 = mean(transcript_summary$delta_median_unpaired_minus_paired > 0, na.rm = TRUE)
)

direction_summary_filtered <- tibble(
  n_transcripts_total = nrow(transcript_summary_filtered),
  n_delta_mean_gt_0 = sum(transcript_summary_filtered$delta_mean_unpaired_minus_paired > 0, na.rm = TRUE),
  n_delta_mean_lt_0 = sum(transcript_summary_filtered$delta_mean_unpaired_minus_paired < 0, na.rm = TRUE),
  n_delta_mean_eq_0 = sum(transcript_summary_filtered$delta_mean_unpaired_minus_paired == 0, na.rm = TRUE),
  frac_delta_mean_gt_0 = mean(transcript_summary_filtered$delta_mean_unpaired_minus_paired > 0, na.rm = TRUE),
  n_delta_median_gt_0 = sum(transcript_summary_filtered$delta_median_unpaired_minus_paired > 0, na.rm = TRUE),
  n_delta_median_lt_0 = sum(transcript_summary_filtered$delta_median_unpaired_minus_paired < 0, na.rm = TRUE),
  n_delta_median_eq_0 = sum(transcript_summary_filtered$delta_median_unpaired_minus_paired == 0, na.rm = TRUE),
  frac_delta_median_gt_0 = mean(transcript_summary_filtered$delta_median_unpaired_minus_paired > 0, na.rm = TRUE)
)

message2("Direction summary (all transcripts):")
print(direction_summary_all)
cat("\n")

message2("Direction summary (coverage-filtered transcripts):")
print(direction_summary_filtered)
cat("\n")

transcript_summary[, wilcox_fdr := p.adjust(wilcox_p, method = "BH")]
transcript_summary_filtered[, wilcox_fdr := p.adjust(wilcox_p, method = "BH")]

signif_summary <- tibble(
  n_transcripts_all = nrow(transcript_summary),
  n_fdr_lt_0_05_all = sum(transcript_summary$wilcox_fdr < 0.05, na.rm = TRUE),
  n_fdr_lt_0_05_and_delta_gt_0_all =
    sum(transcript_summary$wilcox_fdr < 0.05 &
          transcript_summary$delta_mean_unpaired_minus_paired > 0, na.rm = TRUE),
  n_transcripts_filtered = nrow(transcript_summary_filtered),
  n_fdr_lt_0_05_filtered = sum(transcript_summary_filtered$wilcox_fdr < 0.05, na.rm = TRUE),
  n_fdr_lt_0_05_and_delta_gt_0_filtered =
    sum(transcript_summary_filtered$wilcox_fdr < 0.05 &
          transcript_summary_filtered$delta_mean_unpaired_minus_paired > 0, na.rm = TRUE)
)

message2("Transcript significance summary:")
print(signif_summary)
cat("\n")

top_positive_transcripts <- transcript_summary_filtered[
  order(-delta_mean_unpaired_minus_paired, wilcox_fdr)
][1:min(.N, 100)]

top_negative_transcripts <- transcript_summary_filtered[
  order(delta_mean_unpaired_minus_paired, wilcox_fdr)
][1:min(.N, 100)]

gene_summary <- joined_dt[, .(
  n_positions = .N,
  n_transcripts = uniqueN(transcript_id_norm),
  n_paired = sum(vienna_state == "paired", na.rm = TRUE),
  n_unpaired = sum(vienna_state == "unpaired", na.rm = TRUE),
  mean_score_paired = safe_mean(score[vienna_state == "paired"]),
  mean_score_unpaired = safe_mean(score[vienna_state == "unpaired"]),
  median_score_paired = safe_median(score[vienna_state == "paired"]),
  median_score_unpaired = safe_median(score[vienna_state == "unpaired"]),
  delta_mean_unpaired_minus_paired =
    safe_mean(score[vienna_state == "unpaired"]) -
    safe_mean(score[vienna_state == "paired"]),
  delta_median_unpaired_minus_paired =
    safe_median(score[vienna_state == "unpaired"]) -
    safe_median(score[vienna_state == "paired"])
), by = .(gene_id, gene_name)]

setorder(gene_summary, -delta_mean_unpaired_minus_paired)

write_csv(global_state_summary, file.path(output_dir, "global_state_score_summary.csv"))
write_csv(direction_summary_all, file.path(output_dir, "direction_summary_all_transcripts.csv"))
write_csv(direction_summary_filtered, file.path(output_dir, "direction_summary_filtered_transcripts.csv"))
write_csv(signif_summary, file.path(output_dir, "significance_summary.csv"))

fwrite(transcript_summary,
       file.path(output_dir, "transcript_level_followup_summary_all.tsv.gz"),
       sep = "\t")

fwrite(transcript_summary_filtered,
       file.path(output_dir, "transcript_level_followup_summary_filtered.tsv.gz"),
       sep = "\t")

fwrite(top_positive_transcripts,
       file.path(output_dir, "top_positive_transcripts_unpaired_gt_paired.tsv"),
       sep = "\t")

fwrite(top_negative_transcripts,
       file.path(output_dir, "top_negative_transcripts_unpaired_lt_paired.tsv"),
       sep = "\t")

fwrite(gene_summary,
       file.path(output_dir, "gene_level_followup_summary.tsv.gz"),
       sep = "\t")

p1 <- ggplot(as.data.frame(joined_dt), aes(x = vienna_state, y = score)) +
  geom_violin(trim = TRUE) +
  geom_boxplot(width = 0.15, outlier.alpha = 0.15) +
  theme_classic(base_size = 12) +
  labs(
    title = "icSHAPE score by ViennaRNA-predicted state",
    subtitle = paste0("Global Wilcoxon p = ", format(global_wilcox_p, scientific = TRUE)),
    x = "ViennaRNA state",
    y = "icSHAPE score"
  )
print(p1)

p2 <- ggplot(as.data.frame(transcript_summary_filtered),
             aes(x = delta_mean_unpaired_minus_paired)) +
  geom_histogram(bins = 50) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme_classic(base_size = 12) +
  labs(
    title = "Transcript-level icSHAPE difference",
    subtitle = "Coverage-filtered transcripts only",
    x = "Mean icSHAPE(unpaired) - Mean icSHAPE(paired)",
    y = "Number of transcripts"
  )
print(p2)

p3 <- ggplot(as.data.frame(transcript_summary_filtered),
             aes(x = n_positions, y = delta_mean_unpaired_minus_paired)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = 2) +
  theme_classic(base_size = 12) +
  labs(
    title = "Transcript coverage vs paired/unpaired icSHAPE difference",
    x = "Number of joined positions in transcript",
    y = "Mean icSHAPE(unpaired) - Mean icSHAPE(paired)"
  )
print(p3)

p4 <- ggplot(as.data.frame(transcript_summary_filtered),
             aes(x = MFE, y = delta_mean_unpaired_minus_paired)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = 2) +
  theme_classic(base_size = 12) +
  labs(
    title = "Transcript MFE vs paired/unpaired icSHAPE difference",
    x = "ViennaRNA MFE",
    y = "Mean icSHAPE(unpaired) - Mean icSHAPE(paired)"
  )
print(p4)

summary_lines <- c(
  "icSHAPE-MaP vs ViennaRNA follow-up analysis",
  "",
  paste0("Joined input file: ", joined_file),
  paste0("Output directory: ", output_dir),
  "",
  paste0("Global Wilcoxon p-value: ", format(global_wilcox_p, scientific = TRUE)),
  "",
  paste0("Mean score (paired): ",
         signif(global_state_summary[global_state_summary$vienna_state == "paired", "mean_score"], 5)),
  paste0("Mean score (unpaired): ",
         signif(global_state_summary[global_state_summary$vienna_state == "unpaired", "mean_score"], 5)),
  paste0("Median score (paired): ",
         signif(global_state_summary[global_state_summary$vienna_state == "paired", "median_score"], 5)),
  paste0("Median score (unpaired): ",
         signif(global_state_summary[global_state_summary$vienna_state == "unpaired", "median_score"], 5)),
  "",
  paste0("Coverage filter: >= ", min_paired_positions, " paired and >= ", min_unpaired_positions, " unpaired positions"),
  paste0("Filtered transcripts: ", nrow(transcript_summary_filtered)),
  "",
  paste0("Fraction of filtered transcripts with delta_mean > 0: ",
         signif(direction_summary_filtered$frac_delta_mean_gt_0, 4)),
  paste0("Fraction of filtered transcripts with delta_median > 0: ",
         signif(direction_summary_filtered$frac_delta_median_gt_0, 4)),
  "",
  paste0("Filtered transcripts with FDR < 0.05: ",
         signif_summary$n_fdr_lt_0_05_filtered),
  paste0("Filtered transcripts with FDR < 0.05 and delta_mean > 0: ",
         signif_summary$n_fdr_lt_0_05_and_delta_gt_0_filtered)
)

writeLines(summary_lines, file.path(output_dir, "run_summary.txt"))

message2("Follow-up analysis complete. Outputs written to: %s", output_dir)

nonredundant_gene_dt <- copy(transcript_summary_filtered)[
  order(gene_name, -abs(delta_mean_unpaired_minus_paired), wilcox_fdr, -n_positions)
][
  , .SD[1], by = gene_name
]

message2("Nonredundant genes represented by one transcript each: %s",
         format(nrow(nonredundant_gene_dt), big.mark = ","))
cat("\n")

nonredundant_direction_summary <- tibble(
  n_genes_total = nrow(nonredundant_gene_dt),
  n_delta_mean_gt_0 = sum(nonredundant_gene_dt$delta_mean_unpaired_minus_paired > 0, na.rm = TRUE),
  n_delta_mean_lt_0 = sum(nonredundant_gene_dt$delta_mean_unpaired_minus_paired < 0, na.rm = TRUE),
  frac_delta_mean_gt_0 = mean(nonredundant_gene_dt$delta_mean_unpaired_minus_paired > 0, na.rm = TRUE),
  n_fdr_lt_0_05 = sum(nonredundant_gene_dt$wilcox_fdr < 0.05, na.rm = TRUE),
  n_fdr_lt_0_05_and_delta_gt_0 =
    sum(nonredundant_gene_dt$wilcox_fdr < 0.05 &
          nonredundant_gene_dt$delta_mean_unpaired_minus_paired > 0, na.rm = TRUE)
)

message2("Nonredundant gene-level direction summary:")
print(nonredundant_direction_summary)
cat("\n")

top_positive_nonredundant <- nonredundant_gene_dt[
  order(-delta_mean_unpaired_minus_paired, wilcox_fdr)
][1:min(.N, 100)]

top_negative_nonredundant <- nonredundant_gene_dt[
  order(delta_mean_unpaired_minus_paired, wilcox_fdr)
][1:min(.N, 100)]

fwrite(
  nonredundant_gene_dt,
  file.path(output_dir, "nonredundant_one_transcript_per_gene.tsv"),
  sep = "\t"
)

fwrite(
  top_positive_nonredundant,
  file.path(output_dir, "top_positive_nonredundant.tsv"),
  sep = "\t"
)

fwrite(
  top_negative_nonredundant,
  file.path(output_dir, "top_negative_nonredundant.tsv"),
  sep = "\t"
)

p5 <- ggplot(as.data.frame(nonredundant_gene_dt),
             aes(x = delta_mean_unpaired_minus_paired)) +
  geom_histogram(bins = 50) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme_classic(base_size = 12) +
  labs(
    title = "Nonredundant gene-level icSHAPE difference",
    subtitle = "One representative transcript per gene",
    x = "Mean icSHAPE(unpaired) - Mean icSHAPE(paired)",
    y = "Number of genes"
  )
print(p5)

summary_lines_nonredundant <- c(
  "",
  "NONREDUNDANT GENE-LEVEL SUMMARY",
  paste0("One representative transcript per gene: ", nrow(nonredundant_gene_dt)),
  paste0("Genes with delta_mean > 0: ",
         nonredundant_direction_summary$n_delta_mean_gt_0),
  paste0("Genes with delta_mean < 0: ",
         nonredundant_direction_summary$n_delta_mean_lt_0),
  paste0("Fraction with delta_mean > 0: ",
         signif(nonredundant_direction_summary$frac_delta_mean_gt_0, 4)),
  paste0("Genes with FDR < 0.05: ",
         nonredundant_direction_summary$n_fdr_lt_0_05),
  paste0("Genes with FDR < 0.05 and delta_mean > 0: ",
         nonredundant_direction_summary$n_fdr_lt_0_05_and_delta_gt_0)
)

write(
  summary_lines_nonredundant,
  file = file.path(output_dir, "run_summary.txt"),
  append = TRUE
)

message2("Nonredundant gene-level analysis appended successfully.")

figure_dir <- file.path(output_dir, "paper_figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

state_counts <- joined_dt[, .N, by = vienna_state]
state_counts <- state_counts[match(c("paired", "unpaired"), vienna_state)]
x_labels <- c(
  paste0("Paired\n(n=", format(state_counts$N[state_counts$vienna_state == "paired"], big.mark = ","), ")"),
  paste0("Unpaired\n(n=", format(state_counts$N[state_counts$vienna_state == "unpaired"], big.mark = ","), ")")
)

y_max_fig1 <- quantile(joined_dt$score, 0.995, na.rm = TRUE)
y_min_fig1 <- quantile(joined_dt$score, 0.005, na.rm = TRUE)

p_paper_1 <- ggplot(as.data.frame(joined_dt), aes(x = vienna_state, y = score, fill = vienna_state)) +
  geom_violin(trim = TRUE, scale = "width", alpha = 0.85, color = NA) +
  geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white", color = "black", linewidth = 0.35) +
  stat_summary(fun = median, geom = "point", shape = 21, size = 2.4, fill = "white", color = "black") +
  scale_x_discrete(labels = x_labels, limits = c("paired", "unpaired")) +
  scale_fill_manual(values = c("paired" = "#4C78A8", "unpaired" = "#F58518")) +
  coord_cartesian(ylim = c(y_min_fig1, y_max_fig1)) +
  annotate(
    "text",
    x = 1.5,
    y = y_max_fig1,
    label = paste0(
      "Wilcoxon p = ", format(global_wilcox_p, scientific = TRUE), "\n",
      "Mean paired = ", signif(global_state_summary[global_state_summary$vienna_state == "paired", "mean_score"], 4), "\n",
      "Mean unpaired = ", signif(global_state_summary[global_state_summary$vienna_state == "unpaired", "mean_score"], 4)
    ),
    vjust = 1,
    size = 3.8
  ) +
  labs(
    title = "Experimental icSHAPE accessibility follows ViennaRNA structural state",
    subtitle = "Transcript-mapped nucleotide positions across all joined transcripts",
    x = NULL,
    y = "icSHAPE score"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold")
  )

print(p_paper_1)

ggsave(
  filename = file.path(figure_dir, "Figure_icSHAPE_ViennaRNA_global_paired_vs_unpaired.pdf"),
  plot = p_paper_1,
  width = 6.4,
  height = 5.2,
  units = "in"
)

ggsave(
  filename = file.path(figure_dir, "Figure_icSHAPE_ViennaRNA_global_paired_vs_unpaired.png"),
  plot = p_paper_1,
  width = 6.4,
  height = 5.2,
  units = "in",
  dpi = 600
)

n_genes_total_fig2 <- nrow(nonredundant_gene_dt)
n_genes_pos_fig2 <- sum(nonredundant_gene_dt$delta_mean_unpaired_minus_paired > 0, na.rm = TRUE)
frac_genes_pos_fig2 <- mean(nonredundant_gene_dt$delta_mean_unpaired_minus_paired > 0, na.rm = TRUE)
n_genes_sig_pos_fig2 <- sum(
  nonredundant_gene_dt$wilcox_fdr < 0.05 &
    nonredundant_gene_dt$delta_mean_unpaired_minus_paired > 0,
  na.rm = TRUE
)

x_max_fig2 <- max(abs(nonredundant_gene_dt$delta_mean_unpaired_minus_paired), na.rm = TRUE)

p_paper_2 <- ggplot(
  as.data.frame(nonredundant_gene_dt),
  aes(x = delta_mean_unpaired_minus_paired)
) +
  geom_histogram(bins = 36, fill = "#7A7A7A", color = "white", linewidth = 0.25) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.7) +
  annotate(
    "text",
    x = x_max_fig2,
    y = Inf,
    hjust = 1,
    vjust = 1.4,
    size = 3.8,
    label = paste0(
      "Nonredundant genes = ", n_genes_total_fig2, "\n",
      "Positive delta = ", n_genes_pos_fig2, " / ", n_genes_total_fig2,
      " (", sprintf("%.1f", 100 * frac_genes_pos_fig2), "%)\n",
      "FDR < 0.05 and positive = ", n_genes_sig_pos_fig2
    )
  ) +
  labs(
    title = "Agreement remains after collapsing transcript redundancy",
    subtitle = "One representative transcript retained per gene",
    x = "Mean icSHAPE(unpaired) - Mean icSHAPE(paired)",
    y = "Number of genes"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold")
  )

print(p_paper_2)

ggsave(
  filename = file.path(figure_dir, "Figure_icSHAPE_ViennaRNA_nonredundant_gene_distribution.pdf"),
  plot = p_paper_2,
  width = 6.4,
  height = 5.0,
  units = "in"
)

ggsave(
  filename = file.path(figure_dir, "Figure_icSHAPE_ViennaRNA_nonredundant_gene_distribution.png"),
  plot = p_paper_2,
  width = 6.4,
  height = 5.0,
  units = "in",
  dpi = 600
)

plot_dt <- copy(transcript_summary_filtered)
plot_dt[, tx_rank := rank(delta_mean_unpaired_minus_paired, ties.method = "first")]
plot_dt[, sig_group := ifelse(wilcox_fdr < 0.05, "TRUE", "FALSE")]

ggplot(plot_dt, aes(x = tx_rank, y = delta_mean_unpaired_minus_paired)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.6) +
  geom_point(
    data = subset(plot_dt, sig_group == "FALSE"),
    aes(color = sig_group),
    alpha = 0.8,
    size = 1.8
  ) +
  geom_point(
    data = subset(plot_dt, sig_group == "TRUE"),
    aes(color = sig_group),
    alpha = 0.8,
    size = 1.8
  ) +
  scale_color_manual(values = c("FALSE" = "grey50", "TRUE" = "#C62828")) +
  theme_classic(base_size = 13) +
  labs(
    title = "Transcript-level concordance between ViennaRNA and icSHAPE",
    subtitle = "Positive values indicate higher accessibility at ViennaRNA-predicted unpaired positions",
    x = "Transcript rank",
    y = "Mean icSHAPE(unpaired) - Mean icSHAPE(paired)",
    color = "FDR < 0.05"
  )

cat(
  "Transcripts above zero:", sum(plot_dt$delta_mean_unpaired_minus_paired > 0, na.rm = TRUE), "\n",
  "Total transcripts:", sum(!is.na(plot_dt$delta_mean_unpaired_minus_paired)), "\n",
  "Percent above zero:",
  round(
    100 * sum(plot_dt$delta_mean_unpaired_minus_paired > 0, na.rm = TRUE) /
      sum(!is.na(plot_dt$delta_mean_unpaired_minus_paired)),
    1
  ), "%\n",

  "Transcripts with FDR < 0.05:", sum(plot_dt$wilcox_fdr < 0.05, na.rm = TRUE), "\n",
  "Percent with FDR < 0.05:",
  round(
    100 * sum(plot_dt$wilcox_fdr < 0.05, na.rm = TRUE) /
      sum(!is.na(plot_dt$delta_mean_unpaired_minus_paired)),
    1
  ), "%\n",

  "Transcripts above zero and FDR < 0.05:",
  sum(plot_dt$delta_mean_unpaired_minus_paired > 0 & plot_dt$wilcox_fdr < 0.05, na.rm = TRUE), "\n",
  "Percent above zero that are FDR < 0.05:",
  round(
    100 * sum(plot_dt$delta_mean_unpaired_minus_paired > 0 & plot_dt$wilcox_fdr < 0.05, na.rm = TRUE) /
      sum(plot_dt$delta_mean_unpaired_minus_paired > 0, na.rm = TRUE),
    1
  ), "%\n"
)
