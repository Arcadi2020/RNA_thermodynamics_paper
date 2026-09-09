suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(ggplot2)
})

output_dir <- "path_to/Cross_cell_line_ViennaRNA_icSHAPE_summary"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

min_paired_positions <- 10L
min_unpaired_positions <- 10L

joined_files <- c(
  HEK293T = "path_to/HEK293T_icSHAPE_ViennaRNA_joined_per_position_table.tsv.gz",
  K562 = "path_to/K562_CellResearch_ViennaRNA_joined_per_position_table.tsv.gz",
  HeLa = "path_to/HeLa_CellResearch_ViennaRNA_joined_per_position_table.tsv.gz",
  HepG2 = "path_to/HepG2_CellResearch_ViennaRNA_joined_per_position_table.tsv.gz",
  H9 = "path_to/H9_CellResearch_ViennaRNA_joined_per_position_table.tsv.gz"
)

for (nm in names(joined_files)) {
  if (!file.exists(joined_files[[nm]])) {
    stop("Joined file not found for ", nm, ": ", joined_files[[nm]])
  }
}

safe_mean <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

safe_median <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  median(x, na.rm = TRUE)
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

cat("Reading joined tables...\n\n")

joined_list <- lapply(names(joined_files), function(cell_line) {
  dt <- fread(joined_files[[cell_line]])
  dt$cell_line <- cell_line
  dt
})
names(joined_list) <- names(joined_files)

for (nm in names(joined_list)) {
  cat(
    nm, ":\n",
    "  rows = ", format(nrow(joined_list[[nm]]), big.mark = ","), "\n",
    "  transcripts = ", format(data.table::uniqueN(joined_list[[nm]]$transcript_id_norm), big.mark = ","), "\n",
    "  genes = ", format(data.table::uniqueN(joined_list[[nm]]$gene_id), big.mark = ","), "\n\n",
    sep = ""
  )
}

global_state_summary_all <- bind_rows(lapply(names(joined_list), function(cell_line) {
  dt <- joined_list[[cell_line]]

  dt %>%
    group_by(vienna_state) %>%
    summarise(
      n = n(),
      mean_score = mean(score, na.rm = TRUE),
      median_score = median(score, na.rm = TRUE),
      sd_score = sd(score, na.rm = TRUE),
      q25 = quantile(score, 0.25, na.rm = TRUE),
      q75 = quantile(score, 0.75, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(cell_line = cell_line)
}))

global_wilcox_summary <- bind_rows(lapply(names(joined_list), function(cell_line) {
  dt <- joined_list[[cell_line]]

  tibble(
    cell_line = cell_line,
    global_wilcox_p = safe_wilcox_two_group(dt$score, dt$vienna_state),
    n_joined_rows = nrow(dt),
    n_unique_transcripts = dplyr::n_distinct(dt$transcript_id_norm),
    n_unique_genes = dplyr::n_distinct(dt$gene_id)
  )
}))

global_summary_wide <- global_state_summary_all %>%
  select(cell_line, vienna_state, n, mean_score, median_score, q25, q75) %>%
  pivot_wider(
    names_from = vienna_state,
    values_from = c(n, mean_score, median_score, q25, q75)
  ) %>%
  left_join(global_wilcox_summary, by = "cell_line") %>%
  mutate(
    delta_mean_unpaired_minus_paired = mean_score_unpaired - mean_score_paired,
    delta_median_unpaired_minus_paired = median_score_unpaired - median_score_paired
  ) %>%
  select(
    cell_line,
    n_joined_rows,
    n_unique_transcripts,
    n_unique_genes,
    n_paired,
    n_unpaired,
    mean_score_paired,
    mean_score_unpaired,
    delta_mean_unpaired_minus_paired,
    median_score_paired,
    median_score_unpaired,
    delta_median_unpaired_minus_paired,
    global_wilcox_p
  )

cat("Global paired vs unpaired summary across cell lines:\n")
print(global_summary_wide)
cat("\n")

transcript_summary_list <- list()
transcript_summary_filtered_list <- list()
direction_summary_filtered_all <- list()
significance_summary_all <- list()
nonredundant_gene_summary_all <- list()
rank_plot_list <- list()

for (current_cell_line in names(joined_list)) {

  cat("Processing transcript-level summaries for ", current_cell_line, "...\n", sep = "")

  dt <- joined_list[[current_cell_line]]

  transcript_summary <- dt %>%
    group_by(transcript_id_norm, transcript_id) %>%
    summarise(
      cell_line = current_cell_line,
      gene_name = dplyr::first(gene_name),
      transcript_name = dplyr::first(transcript_name),
      Transcript_length = dplyr::first(Transcript_length),
      MFE = dplyr::first(MFE),
      n_positions = n(),
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
      wilcox_p = safe_wilcox_two_group(score, vienna_state),
      .groups = "drop"
    ) %>%
    mutate(
      passes_coverage_filter =
        n_paired >= min_paired_positions &
        n_unpaired >= min_unpaired_positions,
      wilcox_fdr = p.adjust(wilcox_p, method = "BH")
    )

  transcript_summary_filtered <- transcript_summary %>%
    filter(passes_coverage_filter)

  direction_summary_filtered <- transcript_summary_filtered %>%
    summarise(
      cell_line = dplyr::first(cell_line),
      n_transcripts_total = n(),
      n_delta_mean_gt_0 = sum(delta_mean_unpaired_minus_paired > 0, na.rm = TRUE),
      n_delta_mean_lt_0 = sum(delta_mean_unpaired_minus_paired < 0, na.rm = TRUE),
      frac_delta_mean_gt_0 = mean(delta_mean_unpaired_minus_paired > 0, na.rm = TRUE),
      n_delta_median_gt_0 = sum(delta_median_unpaired_minus_paired > 0, na.rm = TRUE),
      n_delta_median_lt_0 = sum(delta_median_unpaired_minus_paired < 0, na.rm = TRUE),
      frac_delta_median_gt_0 = mean(delta_median_unpaired_minus_paired > 0, na.rm = TRUE),
      median_delta_mean = median(delta_mean_unpaired_minus_paired, na.rm = TRUE),
      median_delta_median = median(delta_median_unpaired_minus_paired, na.rm = TRUE)
    )

  significance_summary <- transcript_summary_filtered %>%
    summarise(
      cell_line = dplyr::first(cell_line),
      n_transcripts_filtered = n(),
      n_fdr_lt_0_05 = sum(wilcox_fdr < 0.05, na.rm = TRUE),
      n_fdr_lt_0_05_and_delta_gt_0 =
        sum(wilcox_fdr < 0.05 & delta_mean_unpaired_minus_paired > 0, na.rm = TRUE),
      frac_fdr_lt_0_05 = mean(wilcox_fdr < 0.05, na.rm = TRUE),
      frac_fdr_lt_0_05_and_delta_gt_0 =
        mean(wilcox_fdr < 0.05 & delta_mean_unpaired_minus_paired > 0, na.rm = TRUE)
    )

  nonredundant_gene_summary <- transcript_summary_filtered %>%
    mutate(gene_name_fixed = ifelse(is.na(gene_name) | gene_name == "", transcript_id_norm, gene_name)) %>%
    arrange(gene_name_fixed, desc(abs(delta_mean_unpaired_minus_paired)), wilcox_fdr, desc(n_positions)) %>%
    group_by(gene_name_fixed) %>%
    slice(1) %>%
    ungroup() %>%
    summarise(
      cell_line = dplyr::first(cell_line),
      n_genes_total = n(),
      n_delta_mean_gt_0 = sum(delta_mean_unpaired_minus_paired > 0, na.rm = TRUE),
      n_delta_mean_lt_0 = sum(delta_mean_unpaired_minus_paired < 0, na.rm = TRUE),
      frac_delta_mean_gt_0 = mean(delta_mean_unpaired_minus_paired > 0, na.rm = TRUE),
      n_fdr_lt_0_05 = sum(wilcox_fdr < 0.05, na.rm = TRUE),
      n_fdr_lt_0_05_and_delta_gt_0 =
        sum(wilcox_fdr < 0.05 & delta_mean_unpaired_minus_paired > 0, na.rm = TRUE)
    )

  rank_plot_dt <- transcript_summary_filtered %>%
    arrange(delta_mean_unpaired_minus_paired) %>%
    mutate(
      tx_rank = row_number(),
      sig_group = ifelse(wilcox_fdr < 0.05, "FDR < 0.05", "NS")
    )

  transcript_summary_list[[current_cell_line]] <- transcript_summary
  transcript_summary_filtered_list[[current_cell_line]] <- transcript_summary_filtered
  direction_summary_filtered_all[[current_cell_line]] <- direction_summary_filtered
  significance_summary_all[[current_cell_line]] <- significance_summary
  nonredundant_gene_summary_all[[current_cell_line]] <- nonredundant_gene_summary
  rank_plot_list[[current_cell_line]] <- rank_plot_dt
}

transcript_summary_all <- bind_rows(transcript_summary_list)
transcript_summary_filtered_all <- bind_rows(transcript_summary_filtered_list)
direction_summary_filtered_df <- bind_rows(direction_summary_filtered_all)
significance_summary_df <- bind_rows(significance_summary_all)
nonredundant_gene_summary_df <- bind_rows(nonredundant_gene_summary_all)
rank_plot_df <- bind_rows(rank_plot_list)

cat("Coverage-filtered transcript-level summary across cell lines:\n")
print(direction_summary_filtered_df)
cat("\n")

cat("Significance summary across cell lines:\n")
print(significance_summary_df)
cat("\n")

cat("Nonredundant gene-level summary across cell lines:\n")
print(nonredundant_gene_summary_df)
cat("\n")

safe_signed_rank_greater_than_zero <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_real_)
  if (all(x == 0)) return(NA_real_)

  out <- tryCatch(
    wilcox.test(x, mu = 0, alternative = "greater"),
    error = function(e) NULL
  )

  if (is.null(out)) return(NA_real_)
  out$p.value
}

plot2_summary_df <- transcript_summary_filtered_all %>%
  group_by(cell_line) %>%
  summarise(
    n_transcripts_filtered = n(),
    n_delta_mean_gt_0 = sum(delta_mean_unpaired_minus_paired > 0, na.rm = TRUE),
    frac_delta_mean_gt_0 = mean(delta_mean_unpaired_minus_paired > 0, na.rm = TRUE),
    n_delta_median_gt_0 = sum(delta_median_unpaired_minus_paired > 0, na.rm = TRUE),
    frac_delta_median_gt_0 = mean(delta_median_unpaired_minus_paired > 0, na.rm = TRUE),
    median_delta_mean = median(delta_mean_unpaired_minus_paired, na.rm = TRUE),
    median_delta_median = median(delta_median_unpaired_minus_paired, na.rm = TRUE),
    signed_rank_p_delta_mean_gt_0 =
      safe_signed_rank_greater_than_zero(delta_mean_unpaired_minus_paired),
    signed_rank_p_delta_median_gt_0 =
      safe_signed_rank_greater_than_zero(delta_median_unpaired_minus_paired),
    n_fdr_lt_0_05 = sum(wilcox_fdr < 0.05, na.rm = TRUE),
    frac_fdr_lt_0_05 = mean(wilcox_fdr < 0.05, na.rm = TRUE),
    n_fdr_lt_0_05_and_delta_gt_0 =
      sum(wilcox_fdr < 0.05 & delta_mean_unpaired_minus_paired > 0, na.rm = TRUE),
    frac_fdr_lt_0_05_and_delta_gt_0 =
      mean(wilcox_fdr < 0.05 & delta_mean_unpaired_minus_paired > 0, na.rm = TRUE),
    .groups = "drop"
  )

cat("Transcript-level summary table for histogram interpretation:\n")
print(plot2_summary_df)
cat("\n")

write_csv(
  plot2_summary_df,
  file.path(output_dir, "plot2_transcript_level_shift_summary.csv")
)

write_csv(global_state_summary_all, file.path(output_dir, "global_state_summary_all_cell_lines.csv"))
write_csv(global_summary_wide, file.path(output_dir, "global_summary_wide.csv"))
write_csv(direction_summary_filtered_df, file.path(output_dir, "direction_summary_filtered_all_cell_lines.csv"))
write_csv(significance_summary_df, file.path(output_dir, "significance_summary_all_cell_lines.csv"))
write_csv(nonredundant_gene_summary_df, file.path(output_dir, "nonredundant_gene_summary_all_cell_lines.csv"))

fwrite(
  as.data.table(transcript_summary_all),
  file.path(output_dir, "transcript_summary_all_cell_lines.tsv.gz"),
  sep = "\t"
)

fwrite(
  as.data.table(transcript_summary_filtered_all),
  file.path(output_dir, "transcript_summary_filtered_all_cell_lines.tsv.gz"),
  sep = "\t"
)

plot1_dt <- bind_rows(joined_list)

set.seed(1)

plot1_dt_plot <- plot1_dt %>%
  group_by(cell_line, vienna_state) %>%
  group_modify(~ {
    n_take <- min(nrow(.x), 10000)
    .x %>% slice_sample(n = n_take)
  }) %>%
  ungroup()

plot1_list <- list()

for (current_cell_line in unique(plot1_dt_plot$cell_line)) {

  plot_dt_cell <- plot1_dt_plot %>%
    filter(cell_line == current_cell_line)

  p_cell <- ggplot(plot_dt_cell, aes(x = vienna_state, y = score, fill = vienna_state)) +
    geom_violin(trim = TRUE, scale = "width", alpha = 0.85, color = NA) +
    geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white", color = "black", linewidth = 0.3) +
    stat_summary(fun = median, geom = "point", shape = 21, size = 1.8, fill = "white", color = "black") +
    theme_classic(base_size = 12) +
    labs(
      title = paste0("Experimental icSHAPE accessibility by ViennaRNA-predicted state: ", current_cell_line),
      subtitle = "Global nucleotide-level comparison",
      x = "ViennaRNA state",
      y = "icSHAPE score"
    ) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold")
    )

  plot1_list[[current_cell_line]] <- p_cell

  print(p_cell)
}

icshape_method_df <- tibble::tibble(
  cell_line = c("HEK293T", "K562", "HeLa", "HepG2", "H9"),
  structure_method = c(
    "icSHAPE-MaP",
    "icSHAPE",
    "icSHAPE",
    "icSHAPE",
    "icSHAPE"
  )
)

get_structure_method <- function(current_cell_line) {
  method <- icshape_method_df$structure_method[
    match(current_cell_line, icshape_method_df$cell_line)
  ]

  ifelse(is.na(method), "RNA structure probing", method)
}

plot2_ymax_non_HEK293T <- transcript_summary_filtered_all %>%
  dplyr::filter(cell_line != "HEK293T") %>%
  dplyr::group_split(cell_line) %>%
  lapply(function(df_cell) {

    p_tmp <- ggplot(
      df_cell,
      aes(x = delta_mean_unpaired_minus_paired)
    ) +
      geom_histogram(bins = 50)

    max(ggplot_build(p_tmp)$data[[1]]$count, na.rm = TRUE)
  }) %>%
  unlist() %>%
  max(na.rm = TRUE)

cat(
  "\nPlot 2 y-axis maximum calculated excluding HEK293T:",
  plot2_ymax_non_HEK293T,
  "\n"
)

plot2_list <- list()

for (current_cell_line in unique(transcript_summary_filtered_all$cell_line)) {

  current_method <- get_structure_method(current_cell_line)

  plot2_dt_cell <- transcript_summary_filtered_all %>%
    dplyr::filter(cell_line == current_cell_line)

  p2_cell <- ggplot(
    plot2_dt_cell,
    aes(x = delta_mean_unpaired_minus_paired)
  ) +
    geom_histogram(
      bins = 50,
      fill = "grey50",
      color = "white",
      linewidth = 0.2
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "solid",
      linewidth = 1
    ) +
    theme_classic(base_size = 12) +
    labs(
      title = paste0(
        current_cell_line,
        " — transcript-level concordance between ViennaRNA and ",
        current_method
      ),
      subtitle = "Coverage-filtered transcripts only",
      x = paste0(
        "Mean ",
        current_method,
        "(unpaired) - mean ",
        current_method,
        "(paired)"
      ),
      y = "Number of transcripts"
    ) +
    theme(
      plot.title = element_text(face = "bold")
    )

  if (current_cell_line != "HEK293T") {
    p2_cell <- p2_cell +
      coord_cartesian(
        ylim = c(0, plot2_ymax_non_HEK293T)
      )
  }

  plot2_list[[current_cell_line]] <- p2_cell

  print(p2_cell)

  ggsave(
    filename = file.path(
      output_dir,
      paste0("Plot2_transcript_level_delta_histogram_", current_cell_line, ".tiff")
    ),
    plot = p2_cell,
    device = "tiff",
    width = 4.5,
    height = 4.5,
    units = "in",
    dpi = 300,
    compression = "lzw"
  )
}

cat("\nTranscript-level summary table for Plot 2:\n")
print(plot2_summary_df, n = Inf, width = Inf)
cat("\n")

plot3_ylim <- range(
  rank_plot_df$delta_mean_unpaired_minus_paired,
  na.rm = TRUE
)

plot3_ymax_abs <- max(abs(plot3_ylim), na.rm = TRUE)
plot3_ylim_symmetric <- c(-plot3_ymax_abs, plot3_ymax_abs)

plot3_list <- list()

for (current_cell_line in unique(rank_plot_df$cell_line)) {

  current_method <- get_structure_method(current_cell_line)

  plot3_dt_cell <- rank_plot_df %>%
    dplyr::filter(cell_line == current_cell_line)

  p3_cell <- ggplot(
    plot3_dt_cell,
    aes(x = tx_rank, y = delta_mean_unpaired_minus_paired)
  ) +
    geom_hline(
      yintercept = 0,
      linetype = 2,
      linewidth = 0.5
    ) +
    geom_point(
      aes(color = sig_group),
      alpha = 0.8,
      size = 1.1
    ) +
    scale_color_manual(
      values = c(
        "NS" = "grey50",
        "FDR < 0.05" = "#C62828"
      )
    ) +
    coord_cartesian(
      ylim = plot3_ylim_symmetric
    ) +
    theme_classic(base_size = 12) +
    labs(
      title = paste0(
        current_cell_line,
        " — ranked transcript-level concordance between ViennaRNA and ",
        current_method
      ),
      subtitle = "Positive values indicate higher accessibility at ViennaRNA-predicted unpaired positions",
      x = "Transcript rank",
      y = paste0(
        "Mean ",
        current_method,
        "(unpaired) - mean ",
        current_method,
        "(paired)"
      ),
      color = NULL
    ) +
    theme(
      plot.title = element_text(face = "bold")
    )

  plot3_list[[current_cell_line]] <- p3_cell

  print(p3_cell)
}

plot3_summary_df <- transcript_summary_filtered_all %>%
  group_by(cell_line) %>%
  summarise(
    n_transcripts_filtered = n(),
    n_positive_delta = sum(delta_mean_unpaired_minus_paired > 0, na.rm = TRUE),
    frac_positive_delta = mean(delta_mean_unpaired_minus_paired > 0, na.rm = TRUE),
    n_negative_delta = sum(delta_mean_unpaired_minus_paired < 0, na.rm = TRUE),
    frac_negative_delta = mean(delta_mean_unpaired_minus_paired < 0, na.rm = TRUE),
    median_delta_mean = median(delta_mean_unpaired_minus_paired, na.rm = TRUE),
    n_fdr_lt_0_05 = sum(wilcox_fdr < 0.05, na.rm = TRUE),
    frac_fdr_lt_0_05 = mean(wilcox_fdr < 0.05, na.rm = TRUE),
    n_fdr_lt_0_05_and_positive =
      sum(wilcox_fdr < 0.05 & delta_mean_unpaired_minus_paired > 0, na.rm = TRUE),
    frac_fdr_lt_0_05_and_positive =
      mean(wilcox_fdr < 0.05 & delta_mean_unpaired_minus_paired > 0, na.rm = TRUE),
    frac_significant_that_are_positive = ifelse(
      sum(wilcox_fdr < 0.05, na.rm = TRUE) > 0,
      sum(wilcox_fdr < 0.05 & delta_mean_unpaired_minus_paired > 0, na.rm = TRUE) /
        sum(wilcox_fdr < 0.05, na.rm = TRUE),
      NA_real_
    ),
    signed_rank_p_delta_mean_gt_0 =
      safe_signed_rank_greater_than_zero(delta_mean_unpaired_minus_paired),
    .groups = "drop"
  ) %>%
  mutate(
    frac_positive_delta = round(100 * frac_positive_delta, 1),
    frac_negative_delta = round(100 * frac_negative_delta, 1),
    median_delta_mean = round(median_delta_mean, 3),
    frac_fdr_lt_0_05 = round(100 * frac_fdr_lt_0_05, 1),
    frac_fdr_lt_0_05_and_positive = round(100 * frac_fdr_lt_0_05_and_positive, 1),
    frac_significant_that_are_positive = round(100 * frac_significant_that_are_positive, 1),
    signed_rank_p_delta_mean_gt_0 = format(signed_rank_p_delta_mean_gt_0, scientific = TRUE)
  )

cat("\nSummary table for Plot 3:\n")
print(plot3_summary_df, n = Inf, width = Inf)
cat("\n")

plot4A_dt <- transcript_summary_all %>%
  dplyr::mutate(
    coverage_status = ifelse(
      passes_coverage_filter,
      "Passed coverage filter",
      "Failed coverage filter"
    )
  )

plot4A_list <- list()

for (current_cell_line in unique(plot4A_dt$cell_line)) {

  plot_dt_cell <- plot4A_dt %>%
    dplyr::filter(cell_line == current_cell_line)

  p4A_cell <- ggplot(
    plot_dt_cell,
    aes(x = n_positions)
  ) +
    geom_histogram(
      bins = 80,
      fill = "grey50",
      color = "white",
      linewidth = 0.2
    ) +
    scale_x_log10() +
    theme_classic(base_size = 12) +
    labs(
      title = paste0(
        "Joined nucleotide coverage per transcript: ",
        current_cell_line
      ),
      subtitle = "All transcripts before coverage filtering",
      x = "Joined nucleotide positions per transcript, log10 scale",
      y = "Number of transcripts"
    ) +
    theme(
      plot.title = element_text(face = "bold")
    )

  plot4A_list[[current_cell_line]] <- p4A_cell

  print(p4A_cell)
}

plot4A_summary_df <- plot4A_dt %>%
  dplyr::group_by(cell_line) %>%
  dplyr::summarise(
    n_joined_transcripts = dplyr::n(),
    total_joined_positions = sum(n_positions, na.rm = TRUE),
    min_joined_positions_per_tx = min(n_positions, na.rm = TRUE),
    q25_joined_positions_per_tx = quantile(n_positions, 0.25, na.rm = TRUE),
    median_joined_positions_per_tx = median(n_positions, na.rm = TRUE),
    mean_joined_positions_per_tx = mean(n_positions, na.rm = TRUE),
    q75_joined_positions_per_tx = quantile(n_positions, 0.75, na.rm = TRUE),
    max_joined_positions_per_tx = max(n_positions, na.rm = TRUE),
    n_passed_filter = sum(passes_coverage_filter, na.rm = TRUE),
    frac_passed_filter = mean(passes_coverage_filter, na.rm = TRUE),
    n_failed_filter = sum(!passes_coverage_filter, na.rm = TRUE),
    frac_failed_filter = mean(!passes_coverage_filter, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    frac_passed_filter_percent = round(100 * frac_passed_filter, 1),
    frac_failed_filter_percent = round(100 * frac_failed_filter, 1),
    dplyr::across(
      c(
        q25_joined_positions_per_tx,
        median_joined_positions_per_tx,
        mean_joined_positions_per_tx,
        q75_joined_positions_per_tx
      ),
      ~ round(.x, 1)
    )
  )

cat("\nPlot 4A coverage summary by cell line:\n")
print(plot4A_summary_df, n = Inf, width = Inf)
cat("\n")

plot4B_dt <- transcript_summary_all %>%
  dplyr::mutate(
    coverage_status = ifelse(
      passes_coverage_filter,
      "Passed paired/unpaired filter",
      "Failed paired/unpaired filter"
    ),
    coverage_status = factor(
      coverage_status,
      levels = c(
        "Failed paired/unpaired filter",
        "Passed paired/unpaired filter"
      )
    )
  )

plot4B_list <- list()

for (current_cell_line in unique(plot4B_dt$cell_line)) {

  plot_dt_cell <- plot4B_dt %>%
    dplyr::filter(cell_line == current_cell_line)

  n_total_cell <- nrow(plot_dt_cell)
  n_pass_cell <- sum(plot_dt_cell$passes_coverage_filter, na.rm = TRUE)
  frac_pass_cell <- n_pass_cell / n_total_cell

  subtitle_cell <- paste0(
    n_pass_cell,
    "/",
    n_total_cell,
    " transcripts passed paired/unpaired filter (",
    round(100 * frac_pass_cell, 1),
    "%)"
  )

  p4B_cell <- ggplot(
    plot_dt_cell,
    aes(x = n_positions, fill = coverage_status)
  ) +
    geom_histogram(
      bins = 80,
      position = "identity",
      alpha = 0.75,
      color = "white",
      linewidth = 0.2
    ) +
    scale_x_log10() +
    scale_fill_manual(
      values = c(
        "Failed paired/unpaired filter" = "grey65",
        "Passed paired/unpaired filter" = "grey25"
      )
    ) +
    theme_classic(base_size = 12) +
    labs(
      title = paste0(
        "Coverage-filter retention across joined transcripts: ",
        current_cell_line
      ),
      subtitle = subtitle_cell,
      x = "Total joined nucleotide positions per transcript, log10 scale",
      y = "Number of transcripts",
      fill = NULL
    ) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "top"
    )

  plot4B_list[[current_cell_line]] <- p4B_cell

  print(p4B_cell)
}

plot4B_summary_df <- plot4B_dt %>%
  dplyr::group_by(cell_line, coverage_status) %>%
  dplyr::summarise(
    n_transcripts = dplyr::n(),
    total_joined_positions = sum(n_positions, na.rm = TRUE),
    min_joined_positions_per_tx = min(n_positions, na.rm = TRUE),
    q25_joined_positions_per_tx = quantile(n_positions, 0.25, na.rm = TRUE),
    median_joined_positions_per_tx = median(n_positions, na.rm = TRUE),
    mean_joined_positions_per_tx = mean(n_positions, na.rm = TRUE),
    q75_joined_positions_per_tx = quantile(n_positions, 0.75, na.rm = TRUE),
    max_joined_positions_per_tx = max(n_positions, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::group_by(cell_line) %>%
  dplyr::mutate(
    total_transcripts_in_cell_line = sum(n_transcripts),
    frac_transcripts = n_transcripts / total_transcripts_in_cell_line,
    frac_transcripts_percent = round(100 * frac_transcripts, 1)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    dplyr::across(
      c(
        q25_joined_positions_per_tx,
        median_joined_positions_per_tx,
        mean_joined_positions_per_tx,
        q75_joined_positions_per_tx
      ),
      ~ round(.x, 1)
    )
  )

cat("\nPlot 4B coverage-filter retention summary by cell line:\n")
print(plot4B_summary_df, n = Inf, width = Inf)
cat("\n")

plot4B_dt <- transcript_summary_filtered_all %>%
  mutate(
    min_joined_class_n = pmin(n_paired, n_unpaired)
  ) %>%
  group_by(cell_line, min_joined_class_n) %>%
  summarise(
    n_transcripts = n(),
    .groups = "drop"
  ) %>%
  group_by(cell_line) %>%
  mutate(
    frac_transcripts = n_transcripts / sum(n_transcripts)
  ) %>%
  ungroup()

cat("Summary table for Plot 4B:\n")
print(plot4B_dt, n = Inf, width = Inf)
cat("\n")

plot4B_count_list <- list()

for (current_cell_line in unique(plot4B_dt$cell_line)) {

  plot_dt_cell <- plot4B_dt %>%
    filter(cell_line == current_cell_line)

  x_max <- max(plot_dt_cell$min_joined_class_n, na.rm = TRUE)

  p4B_count <- ggplot(
    plot_dt_cell,
    aes(x = min_joined_class_n, y = n_transcripts)
  ) +
    geom_col(
      width = 0.9,
      fill = "grey70",
      color = "black",
      linewidth = 0.2
    ) +
    geom_vline(
      xintercept = min_paired_positions,
      linetype = 2,
      linewidth = 0.6
    ) +
    scale_x_continuous(
      limits = c(min_paired_positions, x_max),
      breaks = unique(c(min_paired_positions, pretty(c(min_paired_positions, x_max)))),
      expand = c(0, 0)
    ) +
    theme_classic(base_size = 12) +
    labs(
      title = paste0(
        "Balanced paired/unpaired coverage among retained transcripts: ",
        current_cell_line
      ),
      subtitle = paste0(
        "Only transcripts with ≥",
        min_paired_positions,
        " joined nucleotides in both ViennaRNA-predicted paired and unpaired classes are shown"
      ),
      x = "Minimum joined nucleotide count in either structural class",
      y = "Number of retained transcripts"
    ) +
    theme(
      plot.title = element_text(face = "bold")
    )

  plot4B_count_list[[current_cell_line]] <- p4B_count

  print(p4B_count)
}

plot4C_stats_df <- transcript_summary_filtered_all %>%
  dplyr::mutate(
    min_joined_class_n = pmin(n_paired, n_unpaired)
  ) %>%
  dplyr::group_by(cell_line) %>%
  dplyr::summarise(
    n_retained_transcripts = dplyr::n(),

    min_limiting_class_n = min(min_joined_class_n, na.rm = TRUE),
    q25_limiting_class_n = quantile(min_joined_class_n, 0.25, na.rm = TRUE),
    median_limiting_class_n = median(min_joined_class_n, na.rm = TRUE),
    mean_limiting_class_n = mean(min_joined_class_n, na.rm = TRUE),
    q75_limiting_class_n = quantile(min_joined_class_n, 0.75, na.rm = TRUE),
    max_limiting_class_n = max(min_joined_class_n, na.rm = TRUE),

    n_at_minimum_threshold = sum(min_joined_class_n == min_paired_positions, na.rm = TRUE),
    frac_at_minimum_threshold = mean(min_joined_class_n == min_paired_positions, na.rm = TRUE),

    n_limiting_class_ge_20 = sum(min_joined_class_n >= 20, na.rm = TRUE),
    frac_limiting_class_ge_20 = mean(min_joined_class_n >= 20, na.rm = TRUE),

    n_limiting_class_ge_50 = sum(min_joined_class_n >= 50, na.rm = TRUE),
    frac_limiting_class_ge_50 = mean(min_joined_class_n >= 50, na.rm = TRUE),

    n_limiting_class_ge_100 = sum(min_joined_class_n >= 100, na.rm = TRUE),
    frac_limiting_class_ge_100 = mean(min_joined_class_n >= 100, na.rm = TRUE),

    .groups = "drop"
  ) %>%
  dplyr::mutate(
    dplyr::across(
      c(
        q25_limiting_class_n,
        median_limiting_class_n,
        mean_limiting_class_n,
        q75_limiting_class_n
      ),
      ~ round(.x, 1)
    ),
    frac_at_minimum_threshold_percent = round(100 * frac_at_minimum_threshold, 1),
    frac_limiting_class_ge_20_percent = round(100 * frac_limiting_class_ge_20, 1),
    frac_limiting_class_ge_50_percent = round(100 * frac_limiting_class_ge_50, 1),
    frac_limiting_class_ge_100_percent = round(100 * frac_limiting_class_ge_100, 1)
  )

cat("\nPlot 4C limiting paired/unpaired coverage summary by cell line:\n")
print(plot4C_stats_df, n = Inf, width = Inf)
cat("\n")

plot5_dt <- transcript_summary_filtered_all %>%
  mutate(
    min_joined_class_n = pmin(n_paired, n_unpaired)
  ) %>%
  group_by(cell_line) %>%
  mutate(
    q05_cov = quantile(min_joined_class_n, 0.05, na.rm = TRUE),
    q95_cov = quantile(min_joined_class_n, 0.95, na.rm = TRUE),
    coverage_extreme_group = case_when(
      min_joined_class_n <= q05_cov ~ "Lower 5% coverage",
      min_joined_class_n >= q95_cov ~ "Upper 5% coverage",
      TRUE ~ "Middle 90%"
    )
  ) %>%
  ungroup()

plot5_extremes_dt <- plot5_dt %>%
  filter(coverage_extreme_group != "Middle 90%")

cat("\nCoverage-extreme transcript groups:\n")
print(
  plot5_extremes_dt %>%
    count(cell_line, coverage_extreme_group),
  n = Inf
)
cat("\n")

plot5_summary_df <- plot5_extremes_dt %>%
  group_by(cell_line, coverage_extreme_group) %>%
  summarise(
    n_transcripts = n(),
    median_min_joined_class_n = median(min_joined_class_n, na.rm = TRUE),
    median_delta_mean = median(delta_mean_unpaired_minus_paired, na.rm = TRUE),
    frac_positive_delta = mean(delta_mean_unpaired_minus_paired > 0, na.rm = TRUE),
    n_fdr_lt_0_05 = sum(wilcox_fdr < 0.05, na.rm = TRUE),
    frac_fdr_lt_0_05 = mean(wilcox_fdr < 0.05, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    frac_positive_delta = round(100 * frac_positive_delta, 1),
    frac_fdr_lt_0_05 = round(100 * frac_fdr_lt_0_05, 1),
    median_delta_mean = round(median_delta_mean, 3),
    median_min_joined_class_n = round(median_min_joined_class_n, 1)
  )

cat("Summary table for coverage extremes:\n")
print(plot5_summary_df, n = Inf, width = Inf)
cat("\n")

plot5_hist_list <- list()

for (current_cell_line in unique(plot5_extremes_dt$cell_line)) {

  plot_dt_cell <- plot5_extremes_dt %>%
    filter(cell_line == current_cell_line)

  p_hist_cell <- ggplot(
    plot_dt_cell,
    aes(x = delta_mean_unpaired_minus_paired, fill = coverage_extreme_group)
  ) +
    geom_histogram(
      bins = 40,
      alpha = 0.65,
      position = "identity",
      color = "white",
      linewidth = 0.2
    ) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.6) +
    scale_fill_manual(
      values = c(
        "Lower 5% coverage" = "#2166AC",
        "Upper 5% coverage" = "#B2182B"
      )
    ) +
    theme_classic(base_size = 12) +
    labs(
      title = paste0("Coverage extremes mapped onto transcript-level delta distribution: ", current_cell_line),
      subtitle = "Lower 5% and upper 5% of transcripts ranked by minimum joined nucleotide count",
      x = "Mean icSHAPE(unpaired) - Mean icSHAPE(paired)",
      y = "Number of transcripts",
      fill = NULL
    ) +
    theme(
      plot.title = element_text(face = "bold")
    )

  plot5_hist_list[[current_cell_line]] <- p_hist_cell

  print(p_hist_cell)
}

summary_lines <- c(
  "Cross-cell-line ViennaRNA / icSHAPE comparison",
  "",
  "Cell lines included:",
  paste0("  - ", names(joined_files)),
  "",
  paste0("Coverage filter: >= ", min_paired_positions, " paired positions and >= ", min_unpaired_positions, " unpaired positions"),
  "",
  "GLOBAL SUMMARY",
  capture.output(print(global_summary_wide)),
  "",
  "TRANSCRIPT-LEVEL DIRECTION SUMMARY",
  capture.output(print(direction_summary_filtered_df)),
  "",
  "TRANSCRIPT-LEVEL SIGNIFICANCE SUMMARY",
  capture.output(print(significance_summary_df)),
  "",
  "NONREDUNDANT GENE-LEVEL SUMMARY",
  capture.output(print(nonredundant_gene_summary_df))
)

writeLines(summary_lines, file.path(output_dir, "run_summary.txt"))

cat("Cross-cell-line summary complete.\n")
cat("Outputs written to:\n", output_dir, "\n")
