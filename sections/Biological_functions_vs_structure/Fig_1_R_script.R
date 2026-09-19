suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
})

gencode_mfe_file <- "../../data/Gencode_TX_norm_MFE_HC_table.csv"
rna_comparison_rdata_file <- "../../data/Compare_ViennaRNA_vs_RNAstructure_MFE.RData"
paris2_results_rds_file <- "../../data/PARIS2_enriched_mRNA_ViennaRNA_pair_overlap_analysis_final_environment.rds"
icshape_transcript_summary_file <- "../../data/transcript_summary_filtered_all_cell_lines.tsv.gz"
icshape_plot_summary_file <- "../../data/plot2_transcript_level_shift_summary.csv"
shuffle_master_rdata_file <- "../../data/02_analyse_full_shuffle_vs_native_master_tables.RData"
shuffle_selected_objects_rds_file <- "../../data/02_analyse_full_shuffle_vs_native_selected_objects.rds"
output_dir <- "Figure_1_outputs"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

stop_if_missing <- function(file_path, label) {
  if (!file.exists(file_path)) {
    stop(label, " not found: ", file_path, call. = FALSE)
  }
}

load_rdata_objects <- function(file_path, object_names) {
  object_environment <- new.env(parent = baseenv())
  loaded_names <- load(file_path, envir = object_environment)
  missing_names <- setdiff(object_names, loaded_names)
  if (length(missing_names) > 0) {
    stop(
      "Required objects missing from ", basename(file_path), ": ",
      paste(missing_names, collapse = ", "),
      call. = FALSE
    )
  }
  mget(object_names, envir = object_environment, inherits = FALSE)
}

read_named_rds_objects <- function(file_path, object_names) {
  result <- readRDS(file_path)
  if (!is.list(result) || is.null(names(result))) {
    stop("Expected a named list in: ", file_path, call. = FALSE)
  }
  missing_names <- setdiff(object_names, names(result))
  if (length(missing_names) > 0) {
    stop(
      "Required objects missing from ", basename(file_path), ": ",
      paste(missing_names, collapse = ", "),
      call. = FALSE
    )
  }
  result[object_names]
}

input_files <- c(
  gencode_mfe_file,
  rna_comparison_rdata_file,
  paris2_results_rds_file,
  icshape_transcript_summary_file,
  icshape_plot_summary_file,
  shuffle_master_rdata_file,
  shuffle_selected_objects_rds_file
)

for (input_file in input_files) {
  stop_if_missing(input_file, "Input file")
}

Gencode_TX_norm_MFE_HC_table <- read.delim(
  gencode_mfe_file,
  sep = ",",
  header = TRUE,
  stringsAsFactors = FALSE
)
Table_ANNOTATED_gene_names <- Gencode_TX_norm_MFE_HC_table

rna_comparison_objects <- load_rdata_objects(
  rna_comparison_rdata_file,
  c("comparison_table", "pearson_test", "spearman_test")
)

comparison_table <- rna_comparison_objects[["comparison_table"]]
pearson_test <- rna_comparison_objects[["pearson_test"]]
spearman_test <- rna_comparison_objects[["spearman_test"]]

scatter_plot <- ggplot(
  comparison_table,
  aes(x = ViennaRNA_MFE, y = RNAstructure_MFE)
) +
  geom_point(alpha = 0.8, size = 1.2) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    color = "red",
    linewidth = 1
  ) +
  theme_classic() +
  labs(
    title = "ViennaRNA vs RNAstructure MFE",
    subtitle = paste0(
      "Matched transcripts = ", nrow(comparison_table),
      " | Pearson r = ", round(unname(pearson_test$estimate), 3),
      " | Spearman rho = ", round(unname(spearman_test$estimate), 3)
    ),
    x = "ViennaRNA MFE (kcal/mol)",
    y = "RNAstructure MFE (kcal/mol)"
  )

print(scatter_plot)

ggsave(
  filename = file.path(output_dir, "ViennaRNA_vs_RNAstructure_scatter.tiff"),
  plot = scatter_plot,
  width = 7,
  height = 6,
  units = "in",
  dpi = 300,
  compression = "lzw"
)

bland_altman_plot <- ggplot(
  comparison_table,
  aes(x = MFE_mean, y = MFE_difference)
) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = mean(comparison_table$MFE_difference, na.rm = TRUE), color = "red") +
  geom_hline(
    yintercept = mean(comparison_table$MFE_difference, na.rm = TRUE) + 1.96 * sd(comparison_table$MFE_difference, na.rm = TRUE),
    linetype = "dashed",
    color = "blue"
  ) +
  geom_hline(
    yintercept = mean(comparison_table$MFE_difference, na.rm = TRUE) - 1.96 * sd(comparison_table$MFE_difference, na.rm = TRUE),
    linetype = "dashed",
    color = "blue"
  ) +
  theme_classic() +
  labs(
    title = "Bland–Altman Plot: RNAstructure vs ViennaRNA",
    x = "Mean MFE",
    y = "Difference (RNAstructure - ViennaRNA)"
  )

print(bland_altman_plot)

ggsave(
  filename = file.path(output_dir, "ViennaRNA_vs_RNAstructure_BlandAltman.tiff"),
  plot = bland_altman_plot,
  width = 7,
  height = 6,
  units = "in",
  dpi = 300,
  compression = "lzw"
)

agreement_stats <- comparison_table %>%
  mutate(
    abs_MFE_difference = abs(MFE_difference),
    ViennaRNA_Normalized_MFE = ViennaRNA_MFE / length_nt,
    RNAstructure_Normalized_MFE = RNAstructure_MFE / length_nt,
    Normalized_MFE_difference = RNAstructure_Normalized_MFE - ViennaRNA_Normalized_MFE
  )

mean_diff <- mean(agreement_stats$MFE_difference, na.rm = TRUE)
sd_diff   <- sd(agreement_stats$MFE_difference, na.rm = TRUE)

loa_lower <- mean_diff - 1.96 * sd_diff
loa_upper <- mean_diff + 1.96 * sd_diff

mean_norm_diff <- mean(agreement_stats$Normalized_MFE_difference, na.rm = TRUE)

essential_agreement_stats <- c(
  paste0("Number of matched transcripts: ", nrow(agreement_stats)),
  paste0("Mean difference / bias (RNAstructure - ViennaRNA): ", round(mean_diff, 4), " kcal/mol"),
  paste0("95% limits of agreement: ",
         round(loa_lower, 4), " to ",
         round(loa_upper, 4), " kcal/mol"),
  paste0("Median absolute difference: ",
         round(median(agreement_stats$abs_MFE_difference, na.rm = TRUE), 4),
         " kcal/mol"),
  paste0("95th percentile absolute difference: ",
         round(quantile(agreement_stats$abs_MFE_difference, 0.95, na.rm = TRUE), 4),
         " kcal/mol"),
  paste0("Mean normalized difference: ",
         round(mean_norm_diff, 8),
         " kcal/mol/nt")
)

cat("\nEssential agreement statistics:\n")
cat(paste0(essential_agreement_stats, collapse = "\n"), "\n")

paris2_object_names <- c(
  "global_obs_exact",
  "global_obs_t5",
  "global_obs_t10",
  "global_obs_t15",
  "global_perm_exact",
  "global_perm_t5",
  "global_perm_t10",
  "global_perm_t15",
  "duplex_perm_summary"
)

paris2_objects <- read_named_rds_objects(
  paris2_results_rds_file,
  paris2_object_names
)

global_obs_exact <- paris2_objects[["global_obs_exact"]]
global_obs_t5 <- paris2_objects[["global_obs_t5"]]
global_obs_t10 <- paris2_objects[["global_obs_t10"]]
global_obs_t15 <- paris2_objects[["global_obs_t15"]]
global_perm_exact <- paris2_objects[["global_perm_exact"]]
global_perm_t5 <- paris2_objects[["global_perm_t5"]]
global_perm_t10 <- paris2_objects[["global_perm_t10"]]
global_perm_t15 <- paris2_objects[["global_perm_t15"]]
duplex_perm_summary <- paris2_objects[["duplex_perm_summary"]]

empirical_p_ge <- function(observed, null_values) {
  (sum(null_values >= observed, na.rm = TRUE) + 1) / (length(null_values) + 1)
}

panelA_dt <- data.table(
  tolerance = factor(
    c("Exact", "+/- 5 nt", "+/- 10 nt", "+/- 15 nt"),
    levels = c("Exact", "+/- 5 nt", "+/- 10 nt", "+/- 15 nt")
  ),
  observed_real_PARIS2 = c(
    global_obs_exact,
    global_obs_t5,
    global_obs_t10,
    global_obs_t15
  ),
  expected_random_mean = c(
    mean(global_perm_exact, na.rm = TRUE),
    mean(global_perm_t5, na.rm = TRUE),
    mean(global_perm_t10, na.rm = TRUE),
    mean(global_perm_t15, na.rm = TRUE)
  ),
  expected_random_lower = c(
    quantile(global_perm_exact, 0.025, na.rm = TRUE),
    quantile(global_perm_t5, 0.025, na.rm = TRUE),
    quantile(global_perm_t10, 0.025, na.rm = TRUE),
    quantile(global_perm_t15, 0.025, na.rm = TRUE)
  ),
  expected_random_upper = c(
    quantile(global_perm_exact, 0.975, na.rm = TRUE),
    quantile(global_perm_t5, 0.975, na.rm = TRUE),
    quantile(global_perm_t10, 0.975, na.rm = TRUE),
    quantile(global_perm_t15, 0.975, na.rm = TRUE)
  ),
  empirical_p = c(
    empirical_p_ge(global_obs_exact, global_perm_exact),
    empirical_p_ge(global_obs_t5, global_perm_t5),
    empirical_p_ge(global_obs_t10, global_perm_t10),
    empirical_p_ge(global_obs_t15, global_perm_t15)
  )
)

panelA_dt[, enrichment := observed_real_PARIS2 / expected_random_mean]

panelA_dt[, label_text := sprintf(
  "%.2fx\nP = %.3g",
  enrichment,
  empirical_p
)]

panelA_obs <- panelA_dt[, .(
  tolerance,
  group_label = "Observed PARIS2 duplexes",
  value = observed_real_PARIS2
)]

panelA_exp <- panelA_dt[, .(
  tolerance,
  group_label = "Matched random expectation",
  value = expected_random_mean,
  ymin = expected_random_lower,
  ymax = expected_random_upper
)]

panelA_obs[, x := seq_len(.N) - 0.18]
panelA_exp[, x := seq_len(.N) + 0.18]
panelA_dt[, x_mid := seq_len(.N)]

p_panelA <- ggplot() +
  geom_col(
    data = panelA_obs,
    aes(x = x, y = value, fill = group_label),
    width = 0.32,
    color = "black",
    linewidth = 0.25
  ) +
  geom_col(
    data = panelA_exp,
    aes(x = x, y = value, fill = group_label),
    width = 0.32,
    color = "black",
    linewidth = 0.25
  ) +
  geom_errorbar(
    data = panelA_exp,
    aes(x = x, ymin = ymin, ymax = ymax),
    width = 0.10,
    linewidth = 0.55,
    color = "black"
  ) +
  geom_text(
    data = panelA_dt,
    aes(
      x = x_mid,
      y = pmax(observed_real_PARIS2, expected_random_upper) + 0.055,
      label = label_text
    ),
    size = 3.2,
    lineheight = 0.9
  ) +
  scale_fill_manual(
    values = c(
      "Observed PARIS2 duplexes" = "grey50",
      "Matched random expectation" = "white"
    )
  ) +
  scale_x_continuous(
    breaks = seq_len(nrow(panelA_dt)),
    labels = levels(panelA_dt$tolerance)
  ) +
  scale_y_continuous(
    labels = function(x) paste0(round(100 * x), "%")
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    legend.title = element_blank(),
    legend.position = "top",
    legend.text = element_text(size = 11)
  ) +
  labs(
    title = "A. PARIS2 duplexes are supported by ViennaRNA above matched-random expectation",
    x = "Support criterion",
    y = "Duplexes supported by ViennaRNA",
    fill = ""
  ) +
  expand_limits(
    y = max(panelA_dt$expected_random_upper, panelA_dt$observed_real_PARIS2, na.rm = TRUE) + 0.14
  )

print(p_panelA)

panelA_tiff_file <- file.path(
  output_dir,
  "Panel_A_global_observed_vs_matched_random_support_PARIS2.tiff"
)

ggsave(
  filename = panelA_tiff_file,
  plot = p_panelA,
  device = "tiff",
  width = 4,
  height = 6,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

stopifnot(exists("duplex_perm_summary"))

panelB_dt <- duplex_perm_summary[, .(
  n_duplexes = .N,

  observed_frac_t15 = mean(obs_supported_t15, na.rm = TRUE),
  expected_frac_t15 = mean(exp_supported_t15, na.rm = TRUE),

  observed_frac_exact = mean(obs_supported_exact, na.rm = TRUE),
  expected_frac_exact = mean(exp_supported_exact, na.rm = TRUE)
), by = .(gene1, transcript_id, transcript_id_norm)]

panelB_dt[, delta_t15 := observed_frac_t15 - expected_frac_t15]
panelB_dt[, delta_exact := observed_frac_exact - expected_frac_exact]

panelB_summary <- panelB_dt[, .(
  n_transcripts = .N,

  n_above_random_t15 = sum(delta_t15 > 0, na.rm = TRUE),
  n_below_random_t15 = sum(delta_t15 < 0, na.rm = TRUE),
  n_equal_random_t15 = sum(delta_t15 == 0, na.rm = TRUE),
  percent_above_random_t15 = 100 * mean(delta_t15 > 0, na.rm = TRUE),
  median_delta_t15 = median(delta_t15, na.rm = TRUE),
  mean_delta_t15 = mean(delta_t15, na.rm = TRUE),

  n_above_random_exact = sum(delta_exact > 0, na.rm = TRUE),
  n_below_random_exact = sum(delta_exact < 0, na.rm = TRUE),
  n_equal_random_exact = sum(delta_exact == 0, na.rm = TRUE),
  percent_above_random_exact = 100 * mean(delta_exact > 0, na.rm = TRUE),
  median_delta_exact = median(delta_exact, na.rm = TRUE),
  mean_delta_exact = mean(delta_exact, na.rm = TRUE)
)]

panelB_label <- sprintf(
  "%s/%s transcripts above random expectation\n%.1f%% above random; median delta = %.2f",
  format(panelB_summary$n_above_random_t15, big.mark = ","),
  format(panelB_summary$n_transcripts, big.mark = ","),
  panelB_summary$percent_above_random_t15,
  panelB_summary$median_delta_t15
)

p_panelB <- ggplot(
  as.data.frame(panelB_dt),
  aes(x = delta_t15)
) +
  geom_histogram(
    binwidth = 0.025,
    boundary = 0,
    fill = "grey60",
    color = "black",
    linewidth = 0.20,
    alpha = 0.85
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 1.3,
    color = "black"
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = panelB_label,
    hjust = 1.05,
    vjust = 1.4,
    size = 2
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  ) +
  labs(

    x = "Observed - expected fraction supported at +/- 15 nt",
    y = "Number of transcripts"
  )

print(p_panelB)

panelB_tiff_file <- file.path(
  output_dir,
  "Panel_B_transcript_level_observed_minus_expected_t15_PARIS2.tiff"
)

ggsave(
  filename = panelB_tiff_file,
  plot = p_panelB,
  device = "tiff",
  width = 6,
  height = 6,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

transcript_summary_filtered_all <- tibble::as_tibble(
  fread(icshape_transcript_summary_file)
)

plot2_summary_df <- tibble::as_tibble(
  fread(icshape_plot_summary_file)
)

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
      linewidth = 1.2
    ) +
    theme_classic(base_size = 12) +
    labs(
      title = paste0(
        current_cell_line,
        " - ",
        current_method
      ),

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
      plot.title = element_text(face = "bold", size = 13),

      axis.text.x = element_text(size = 12, colour = "black"),
      axis.text.y = element_text(size = 12, colour = "black"),

      axis.title.x = element_text(size = 9),
      axis.title.y = element_text(size = 9),

      axis.line = element_line(linewidth = 0.5, colour = "black"),
      axis.ticks = element_line(linewidth = 0.5, colour = "black")
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

plot2_cell_line_order <- c("HeLa", "HepG2", "H9", "K562", "HEK293T")

plot2_list_ordered <- plot2_list[plot2_cell_line_order]
plot2_list_ordered <- plot2_list_ordered[!vapply(plot2_list_ordered, is.null, logical(1))]

plot2_combined_single_row <- patchwork::wrap_plots(
  plot2_list_ordered,
  nrow = 1
)

print(plot2_combined_single_row)

ggsave(
  filename = file.path(
    output_dir,
    "Plot2_transcript_level_delta_histograms_five_cell_lines.tiff"
  ),
  plot = plot2_combined_single_row,
  device = "tiff",
  width = 22.5,
  height = 4.5,
  units = "in",
  dpi = 300,
  compression = "lzw"
)

shuffle_master_objects <- load_rdata_objects(
  shuffle_master_rdata_file,
  "master_transcript_table"
)

master_transcript_table <- shuffle_master_objects[["master_transcript_table"]]

stopifnot(exists("master_transcript_table"))
stopifnot("transcript_id" %in% names(master_transcript_table))
stopifnot("transcript_biotype" %in% names(master_transcript_table))
stopifnot("native_Normalized_MFE" %in% names(master_transcript_table))
stopifnot("transcript_is_canonical" %in% names(master_transcript_table))

fig1e_biotype_levels <- c(
  "protein_coding",
  "tRNA",
  "rRNA"
)

fig1e_native_df <- master_transcript_table %>%
  dplyr::filter(
    !is.na(transcript_id),
    transcript_id != "",
    !is.na(transcript_biotype),
    transcript_biotype %in% fig1e_biotype_levels,
    !is.na(native_Normalized_MFE)
  ) %>%
  dplyr::mutate(
    Group = dplyr::case_when(
      transcript_biotype == "protein_coding" ~ "protein_coding",
      transcript_biotype == "tRNA" ~ "tRNA",
      transcript_biotype == "rRNA" ~ "rRNA",
      TRUE ~ NA_character_
    ),
    transcript_is_canonical_chr = as.character(transcript_is_canonical),
    canonical_flag = dplyr::case_when(
      transcript_is_canonical_chr %in% c("1", "TRUE", "true", "True") ~ 1L,
      TRUE ~ 0L
    ),

    gene_selection_id = dplyr::case_when(
      !is.na(ensembl_gene_id) & ensembl_gene_id != "" ~ ensembl_gene_id,
      TRUE ~ transcript_id
    )
  ) %>%
  dplyr::filter(!is.na(Group)) %>%
  dplyr::group_by(gene_selection_id) %>%
  dplyr::mutate(
    n_master_transcripts_for_gene = dplyr::n(),
    selected_transcript_rule = dplyr::case_when(
      any(canonical_flag == 1L, na.rm = TRUE) ~
        "canonical_transcript_present_in_master_table",
      TRUE ~
        "no_canonical_in_master_table_assigned_available_transcript"
    )
  ) %>%
  dplyr::arrange(
    dplyr::desc(canonical_flag),
    transcript_id,
    .by_group = TRUE
  ) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(
    transcript_id,
    ensembl_gene_id,
    external_gene_name,
    transcript_biotype,
    Group = factor(Group, levels = fig1e_biotype_levels),
    selected_transcript_rule,
    canonical_flag,
    n_master_transcripts_for_gene,
    Normalized_MFE = as.numeric(native_Normalized_MFE)
  ) %>%
  dplyr::filter(!is.na(Group))

cat("\nFig.1E native MFE table restricted to canonical/representative transcripts:\n")
cat("Rows:", nrow(fig1e_native_df), "\n")
cat("Unique transcripts:", dplyr::n_distinct(fig1e_native_df$transcript_id), "\n")

fig1e_native_df %>%
  dplyr::count(Group, name = "n_selected_transcripts") %>%
  tibble::as_tibble() %>%
  print(n = Inf)

stopifnot(
  all(fig1e_biotype_levels %in% as.character(unique(fig1e_native_df$Group)))
)

fig1e_native_summary <- fig1e_native_df %>%
  dplyr::group_by(Group) %>%
  dplyr::summarise(
    n_tx = dplyr::n(),
    mean_Normalized_MFE = mean(Normalized_MFE, na.rm = TRUE),
    median_Normalized_MFE = median(Normalized_MFE, na.rm = TRUE),
    sd_Normalized_MFE = sd(Normalized_MFE, na.rm = TRUE),
    iqr_Normalized_MFE = IQR(Normalized_MFE, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    mean_Normalized_MFE = round(mean_Normalized_MFE, 4),
    median_Normalized_MFE = round(median_Normalized_MFE, 4),
    sd_Normalized_MFE = round(sd_Normalized_MFE, 4),
    iqr_Normalized_MFE = round(iqr_Normalized_MFE, 4)
  )

cat("\nFig.1E native normalized MFE summary:\n")
fig1e_native_summary %>%
  tibble::as_tibble() %>%
  print(n = Inf, width = Inf)

protein_vals <- fig1e_native_df %>%
  dplyr::filter(Group == "protein_coding") %>%
  dplyr::pull(Normalized_MFE) %>%
  stats::na.omit()

fig1e_native_vs_protein_stats <- purrr::map_dfr(
  c("tRNA", "rRNA"),
  function(test_group) {

    test_vals <- fig1e_native_df %>%
      dplyr::filter(Group == test_group) %>%
      dplyr::pull(Normalized_MFE) %>%
      stats::na.omit()

    tibble::tibble(
      comparison = paste0(test_group, "_vs_protein_coding"),
      reference_group = "protein_coding",
      test_group = test_group,

      n_protein_coding = length(protein_vals),
      n_test_group = length(test_vals),

      mean_protein_coding = mean(protein_vals),
      mean_test_group = mean(test_vals),
      mean_difference_test_minus_protein =
        mean(test_vals) - mean(protein_vals),

      median_protein_coding = median(protein_vals),
      median_test_group = median(test_vals),
      median_difference_test_minus_protein =
        median(test_vals) - median(protein_vals),

      wilcox_p_two_sided = wilcox.test(
        protein_vals,
        test_vals,
        alternative = "two.sided",
        exact = FALSE
      )$p.value,

      wilcox_p_test_more_negative_than_protein = wilcox.test(
        test_vals,
        protein_vals,
        alternative = "less",
        exact = FALSE
      )$p.value,

      ks_p_two_sided = ks.test(
        protein_vals,
        test_vals
      )$p.value
    )
  }
) %>%
  dplyr::mutate(
    wilcox_p_two_sided_BH =
      p.adjust(wilcox_p_two_sided, method = "BH"),
    wilcox_p_test_more_negative_than_protein_BH =
      p.adjust(wilcox_p_test_more_negative_than_protein, method = "BH"),
    ks_p_two_sided_BH =
      p.adjust(ks_p_two_sided, method = "BH"),

    direction_by_median = dplyr::case_when(
      median_difference_test_minus_protein < 0 ~
        "test_group_more_negative_than_protein_coding",
      median_difference_test_minus_protein > 0 ~
        "test_group_less_negative_than_protein_coding",
      TRUE ~
        "no_median_difference"
    )
  ) %>%
  dplyr::mutate(
    dplyr::across(
      .cols = c(
        mean_protein_coding,
        mean_test_group,
        mean_difference_test_minus_protein,
        median_protein_coding,
        median_test_group,
        median_difference_test_minus_protein
      ),
      .fns = ~ round(.x, 4)
    )
  )

cat("\nFig.1E protein-coding versus tRNA/rRNA native normalized MFE statistics:\n")
fig1e_native_vs_protein_stats %>%
  tibble::as_tibble() %>%
  print(n = Inf, width = Inf)

fig1e_pairwise_wilcox <- pairwise.wilcox.test(
  x = fig1e_native_df$Normalized_MFE,
  g = fig1e_native_df$Group,
  p.adjust.method = "BH",
  exact = FALSE
)$p.value %>%
  as.data.frame() %>%
  tibble::rownames_to_column("group_1") %>%
  tidyr::pivot_longer(
    cols = -group_1,
    names_to = "group_2",
    values_to = "wilcox_p_adj"
  ) %>%
  dplyr::filter(!is.na(wilcox_p_adj)) %>%
  dplyr::arrange(wilcox_p_adj)

cat("\nFig.1E all-pairwise Wilcoxon comparisons:\n")
fig1e_pairwise_wilcox %>%
  tibble::as_tibble() %>%
  print(n = Inf, width = Inf)

fig1e_counts <- fig1e_native_df %>%
  dplyr::count(Group, name = "n_tx")

fig1e_group_label_map <- c(
  protein_coding = "Protein-coding",
  tRNA = "tRNA",
  rRNA = "rRNA"
)

fig1e_axis_labels <- paste0(
  fig1e_group_label_map[fig1e_biotype_levels],
  "\n(n=",
  scales::comma(
    fig1e_counts$n_tx[
      match(fig1e_biotype_levels, as.character(fig1e_counts$Group))
    ]
  ),
  ")"
)

names(fig1e_axis_labels) <- fig1e_biotype_levels

fig1e_biotype_levels_inverted <- rev(fig1e_biotype_levels)

MFE_boxplot_flipped <- ggplot(
  fig1e_native_df,
  aes(
    x = Group,
    y = Normalized_MFE
  )
) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    fill = "white",
    colour = "black",
    linewidth = 0.4
  ) +
  scale_x_discrete(
    limits = fig1e_biotype_levels_inverted,
    labels = fig1e_axis_labels
  ) +
  coord_flip() +
  theme_classic(base_size = 13) +
  labs(
    x = NULL,
    y = "Native normalized ΔGfold (kcal·mol⁻¹·nt⁻¹)"
  ) +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14)
  )

print(MFE_boxplot_flipped)

ggsave(
  filename = file.path(output_dir, "Fig1E_native_normalized_DGfold_boxplot.tiff"),
  plot = MFE_boxplot_flipped,
  device = "tiff",
  width = 7,
  height = 3,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

fig1e_selected_transcript_ids <- fig1e_native_df %>%
  dplyr::pull(transcript_id) %>%
  unique()

rna_class_plot_df <- master_transcript_table %>%
  dplyr::filter(
    transcript_id %in% fig1e_selected_transcript_ids,
    !is.na(transcript_id),
    transcript_id != "",
    !is.na(transcript_biotype),
    transcript_biotype %in% c("protein_coding", "tRNA", "rRNA"),
    !is.na(native_Normalized_MFE),
    !is.na(mean_shuffled_Normalized_MFE),
    !is.na(delta_Normalized_MFE)
  ) %>%
  dplyr::mutate(
    group_name = dplyr::case_when(
      transcript_biotype == "protein_coding" ~ "Protein_coding_RNAs",
      transcript_biotype == "tRNA" ~ "tRNAs",
      transcript_biotype == "rRNA" ~ "rRNAs",
      TRUE ~ NA_character_
    ),
    group_name = factor(
      group_name,
      levels = c("Protein_coding_RNAs", "tRNAs", "rRNAs")
    )
  ) %>%
  dplyr::filter(!is.na(group_name)) %>%
  dplyr::distinct(transcript_id, .keep_all = TRUE)

cat("\nMatched shuffled-control RNA-class table restricted to Fig.1E transcript universe:\n")
cat("Rows:", nrow(rna_class_plot_df), "\n")
cat("Unique transcripts:", dplyr::n_distinct(rna_class_plot_df$transcript_id), "\n")
cat("Expected selected Fig.1E transcripts:", length(fig1e_selected_transcript_ids), "\n\n")

rna_class_plot_df %>%
  dplyr::count(group_name, name = "n_transcripts") %>%
  tibble::as_tibble() %>%
  print(n = Inf)

stopifnot(
  setequal(
    rna_class_plot_df$transcript_id,
    fig1e_selected_transcript_ids
  )
)

rna_class_native_summary <- rna_class_plot_df %>%
  dplyr::group_by(group_name) %>%
  dplyr::summarise(
    n_tx = dplyr::n_distinct(transcript_id),
    n_genes = dplyr::n_distinct(ensembl_gene_id[!is.na(ensembl_gene_id)]),
    mean_native_Normalized_MFE = mean(native_Normalized_MFE, na.rm = TRUE),
    median_native_Normalized_MFE = median(native_Normalized_MFE, na.rm = TRUE),
    sd_native_Normalized_MFE = sd(native_Normalized_MFE, na.rm = TRUE),
    iqr_native_Normalized_MFE = IQR(native_Normalized_MFE, na.rm = TRUE),
    .groups = "drop"
  )

rna_class_native_pairwise_stats <- pairwise.wilcox.test(
  x = rna_class_plot_df$native_Normalized_MFE,
  g = rna_class_plot_df$group_name,
  p.adjust.method = "BH",
  exact = FALSE
)$p.value %>%
  as.data.frame() %>%
  tibble::rownames_to_column("group_1") %>%
  tidyr::pivot_longer(
    cols = -group_1,
    names_to = "group_2",
    values_to = "wilcox_p_adj"
  ) %>%
  dplyr::filter(!is.na(wilcox_p_adj)) %>%
  dplyr::mutate(metric = "native_Normalized_MFE", .before = 1) %>%
  dplyr::arrange(wilcox_p_adj)

cat("\nNative normalized MFE summary in matched shuffled-control universe:\n")
rna_class_native_summary %>%
  tibble::as_tibble() %>%
  print(n = Inf, width = Inf)

cat("\nPairwise native normalized MFE comparisons:\n")
rna_class_native_pairwise_stats %>%
  tibble::as_tibble() %>%
  print(n = Inf, width = Inf)

protein_native_vals <- rna_class_plot_df %>%
  dplyr::filter(group_name == "Protein_coding_RNAs") %>%
  dplyr::pull(native_Normalized_MFE) %>%
  stats::na.omit()

rna_class_native_vs_protein_stats <- purrr::map_dfr(
  c("tRNAs", "rRNAs"),
  function(test_group) {

    test_vals <- rna_class_plot_df %>%
      dplyr::filter(group_name == test_group) %>%
      dplyr::pull(native_Normalized_MFE) %>%
      stats::na.omit()

    tibble::tibble(
      comparison = paste0(test_group, "_vs_Protein_coding_RNAs"),
      reference_group = "Protein_coding_RNAs",
      test_group = test_group,
      n_protein_coding = length(protein_native_vals),
      n_test_group = length(test_vals),
      mean_protein_coding = mean(protein_native_vals),
      mean_test_group = mean(test_vals),
      mean_difference_test_minus_protein =
        mean(test_vals) - mean(protein_native_vals),
      median_protein_coding = median(protein_native_vals),
      median_test_group = median(test_vals),
      median_difference_test_minus_protein =
        median(test_vals) - median(protein_native_vals),
      wilcox_p_two_sided = wilcox.test(
        protein_native_vals,
        test_vals,
        alternative = "two.sided",
        exact = FALSE
      )$p.value,
      wilcox_p_test_more_negative_than_protein = wilcox.test(
        test_vals,
        protein_native_vals,
        alternative = "less",
        exact = FALSE
      )$p.value
    )
  }
) %>%
  dplyr::mutate(
    wilcox_p_two_sided_BH =
      p.adjust(wilcox_p_two_sided, method = "BH"),
    wilcox_p_test_more_negative_than_protein_BH =
      p.adjust(wilcox_p_test_more_negative_than_protein, method = "BH"),
    dplyr::across(
      .cols = c(
        mean_protein_coding,
        mean_test_group,
        mean_difference_test_minus_protein,
        median_protein_coding,
        median_test_group,
        median_difference_test_minus_protein
      ),
      .fns = ~ round(.x, 4)
    )
  )

cat("\nFocused native comparisons: protein-coding versus tRNA/rRNA\n")
rna_class_native_vs_protein_stats %>%
  tibble::as_tibble() %>%
  print(n = Inf, width = Inf)

rna_class_shuffle_summary <- rna_class_plot_df %>%
  dplyr::group_by(group_name) %>%
  dplyr::summarise(
    n_tx = dplyr::n_distinct(transcript_id),
    n_genes = dplyr::n_distinct(ensembl_gene_id[!is.na(ensembl_gene_id)]),

    mean_native_Normalized_MFE =
      mean(native_Normalized_MFE, na.rm = TRUE),
    mean_mean_shuffled_Normalized_MFE =
      mean(mean_shuffled_Normalized_MFE, na.rm = TRUE),

    median_native_Normalized_MFE =
      median(native_Normalized_MFE, na.rm = TRUE),
    median_mean_shuffled_Normalized_MFE =
      median(mean_shuffled_Normalized_MFE, na.rm = TRUE),

    mean_delta_Normalized_MFE =
      mean(delta_Normalized_MFE, na.rm = TRUE),
    median_delta_Normalized_MFE =
      median(delta_Normalized_MFE, na.rm = TRUE),
    sd_delta_Normalized_MFE =
      sd(delta_Normalized_MFE, na.rm = TRUE),
    iqr_delta_Normalized_MFE =
      IQR(delta_Normalized_MFE, na.rm = TRUE),

    pct_native_more_stable_than_shuffle =
      100 * mean(delta_Normalized_MFE < 0, na.rm = TRUE),

    n_native_more_stable_than_shuffle =
      sum(delta_Normalized_MFE < 0, na.rm = TRUE),
    n_tested =
      sum(!is.na(delta_Normalized_MFE)),

    wilcox_p_delta_less_than_zero = tryCatch(
      wilcox.test(
        delta_Normalized_MFE,
        mu = 0,
        alternative = "less",
        exact = FALSE
      )$p.value,
      error = function(e) NA_real_
    ),

    binom_p_more_than_half_less_than_zero = tryCatch(
      binom.test(
        x = sum(delta_Normalized_MFE < 0, na.rm = TRUE),
        n = sum(!is.na(delta_Normalized_MFE)),
        p = 0.5,
        alternative = "greater"
      )$p.value,
      error = function(e) NA_real_
    ),

    .groups = "drop"
  ) %>%
  dplyr::mutate(
    wilcox_p_adj_BH =
      p.adjust(wilcox_p_delta_less_than_zero, method = "BH"),
    binom_p_adj_BH =
      p.adjust(binom_p_more_than_half_less_than_zero, method = "BH"),

    shuffled_control_conclusion = dplyr::case_when(
      wilcox_p_adj_BH < 0.05 & median_delta_Normalized_MFE < 0 ~
        "Native sequences significantly more stable than matched shuffled controls",
      TRUE ~
        "No significant native excess stability detected after matched shuffling"
    )
  )

cat("\nNative-minus-shuffled summary by RNA class:\n")
rna_class_shuffle_summary %>%
  tibble::as_tibble() %>%
  print(n = Inf, width = Inf)

rna_class_delta_pairwise_stats <- pairwise.wilcox.test(
  x = rna_class_plot_df$delta_Normalized_MFE,
  g = rna_class_plot_df$group_name,
  p.adjust.method = "BH",
  exact = FALSE
)$p.value %>%
  as.data.frame() %>%
  tibble::rownames_to_column("group_1") %>%
  tidyr::pivot_longer(
    cols = -group_1,
    names_to = "group_2",
    values_to = "wilcox_p_adj"
  ) %>%
  dplyr::filter(!is.na(wilcox_p_adj)) %>%
  dplyr::mutate(metric = "delta_Normalized_MFE", .before = 1) %>%
  dplyr::arrange(wilcox_p_adj)

cat("\nPairwise comparisons of delta_Normalized_MFE:\n")
rna_class_delta_pairwise_stats %>%
  tibble::as_tibble() %>%
  print(n = Inf, width = Inf)

protein_delta_vals <- rna_class_plot_df %>%
  dplyr::filter(group_name == "Protein_coding_RNAs") %>%
  dplyr::pull(delta_Normalized_MFE) %>%
  stats::na.omit()

rna_class_delta_vs_protein_stats <- purrr::map_dfr(
  c("tRNAs", "rRNAs"),
  function(test_group) {

    test_vals <- rna_class_plot_df %>%
      dplyr::filter(group_name == test_group) %>%
      dplyr::pull(delta_Normalized_MFE) %>%
      stats::na.omit()

    tibble::tibble(
      comparison = paste0(test_group, "_vs_Protein_coding_RNAs"),
      reference_group = "Protein_coding_RNAs",
      test_group = test_group,
      n_protein_coding = length(protein_delta_vals),
      n_test_group = length(test_vals),
      mean_delta_protein_coding = mean(protein_delta_vals),
      mean_delta_test_group = mean(test_vals),
      mean_delta_difference_test_minus_protein =
        mean(test_vals) - mean(protein_delta_vals),
      median_delta_protein_coding = median(protein_delta_vals),
      median_delta_test_group = median(test_vals),
      median_delta_difference_test_minus_protein =
        median(test_vals) - median(protein_delta_vals),
      pct_protein_coding_delta_below_zero =
        100 * mean(protein_delta_vals < 0),
      pct_test_group_delta_below_zero =
        100 * mean(test_vals < 0),
      wilcox_p_two_sided = wilcox.test(
        protein_delta_vals,
        test_vals,
        alternative = "two.sided",
        exact = FALSE
      )$p.value,
      wilcox_p_test_more_negative_delta_than_protein = wilcox.test(
        test_vals,
        protein_delta_vals,
        alternative = "less",
        exact = FALSE
      )$p.value
    )
  }
) %>%
  dplyr::mutate(
    wilcox_p_two_sided_BH =
      p.adjust(wilcox_p_two_sided, method = "BH"),
    wilcox_p_test_more_negative_delta_than_protein_BH =
      p.adjust(wilcox_p_test_more_negative_delta_than_protein, method = "BH"),
    dplyr::across(
      .cols = c(
        mean_delta_protein_coding,
        mean_delta_test_group,
        mean_delta_difference_test_minus_protein,
        median_delta_protein_coding,
        median_delta_test_group,
        median_delta_difference_test_minus_protein,
        pct_protein_coding_delta_below_zero,
        pct_test_group_delta_below_zero
      ),
      .fns = ~ round(.x, 5)
    )
  )

cat("\nFocused delta comparisons: protein-coding versus tRNA/rRNA\n")
rna_class_delta_vs_protein_stats %>%
  tibble::as_tibble() %>%
  print(n = Inf, width = Inf)

rna_class_group_levels <- c(
  "Protein_coding_RNAs",
  "tRNAs",
  "rRNAs"
)

rna_class_group_levels_inverted <- rev(rna_class_group_levels)

rna_class_label_map <- c(
  Protein_coding_RNAs = "Protein-coding RNAs",
  tRNAs = "tRNAs",
  rRNAs = "rRNAs"
)

rna_class_shuffle_summary_for_plot <- rna_class_shuffle_summary %>%
  dplyr::mutate(
    group_name = factor(
      group_name,
      levels = rna_class_group_levels
    ),
    label_text = paste0(
      "tx=", scales::comma(n_tx),
      "\n", round(pct_native_more_stable_than_shuffle, 1), "% < 0"
    )
  )

rna_class_plot_df_flipped <- rna_class_plot_df %>%
  dplyr::mutate(
    group_name = factor(
      group_name,
      levels = rna_class_group_levels
    )
  )

y_min <- min(rna_class_plot_df_flipped$delta_Normalized_MFE, na.rm = TRUE)
y_max <- max(rna_class_plot_df_flipped$delta_Normalized_MFE, na.rm = TRUE)
y_range <- y_max - y_min

p_rna_class_delta_norm_flipped <- ggplot(
  rna_class_plot_df_flipped,
  aes(
    x = group_name,
    y = delta_Normalized_MFE
  )
) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    fill = "#337AB7",
    color = "black",
    linewidth = 0.4
  ) +
  geom_hline(
    yintercept = 0,
    linetype = 2,
    color = "black",
    linewidth = 0.6
  ) +
  geom_text(
    data = rna_class_shuffle_summary_for_plot,
    aes(
      x = group_name,
      y = y_max + 0.10 * y_range,
      label = label_text
    ),
    inherit.aes = FALSE,
    size = 3.4
  ) +
  scale_x_discrete(
    limits = rna_class_group_levels_inverted,
    labels = rna_class_label_map
  ) +
  coord_flip(
    ylim = c(y_min, y_max + 0.18 * y_range),
    clip = "off"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.margin = margin(t = 5.5, r = 45, b = 5.5, l = 5.5),
    legend.position = "none",
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14)
  ) +
  labs(
    title = "RNA class stability relative to matched shuffled controls",
    subtitle = "Same canonical/representative transcript universe as Fig.1E",
    x = NULL,
    y = "Native − shuffled normalized ΔGfold (kcal·mol⁻¹·nt⁻¹)"
  )

print(p_rna_class_delta_norm_flipped)

ggsave(
  filename = file.path(output_dir, "Fig1E_p_rna_class_delta_norm.tiff"),
  plot = p_rna_class_delta_norm_flipped,
  device = "tiff",
  width = 7,
  height = 3,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

ggplot(Table_ANNOTATED_gene_names, aes(x = Normalized_MFE)) +
  geom_histogram(binwidth = 0.001, color = "black") +
  labs(title = "Distribution of the Minimum Free Energy (Normalized_MFE) values of the Human transcriptome",
       x = "Normalized_MFE (kcal/mol)",
       y = "Frequency") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_reverse(limits = c(0, -0.8)) +

  theme(axis.text.y = element_text(size = 10)) +
  theme(axis.title.x = element_text(size = 12)) +
  theme(axis.title.y = element_text(size = 12)) +
  theme(legend.position = "none")

ggplot(Table_ANNOTATED_gene_names, aes(x = reorder(ensembl_transcript_id_version, Normalized_MFE), y = Normalized_MFE)) +
  theme_classic() +
  geom_point(stat = "identity", colour = "black", size = 2) +

  theme(axis.text.y = element_blank()) +
  theme(axis.title=element_text(size=20)) +
  xlab("Human annotated transcriptome") +
  ylab("Normalized_MFE (kcal/mol)") +
  ggtitle("Distribution of the Human transcriptome Normalized_MFE") +
  coord_flip()

Table_ANNOTATED_gene_filtered <- dplyr::filter(Table_ANNOTATED_gene_names, transcript_biotype == "protein_coding")

ggplot(Table_ANNOTATED_gene_filtered, aes(x = Normalized_MFE)) +
  geom_histogram(binwidth = 0.001, color = "black") +
  labs(title = "Distribution of Normalized Minimum Free Energy Values for Human Protein-Coding mRNAs",
       x = "Normalized_MFE (kcal/mol)",
       y = "Frequency") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_reverse(limits = c(0, -0.8)) +

  theme(axis.text.y = element_text(size = 10)) +
  theme(axis.title.x = element_text(size = 12)) +
  theme(axis.title.y = element_text(size = 12)) +
  theme(legend.position = "none")

Human_protein_coding_TX <- ggplot(Table_ANNOTATED_gene_filtered, aes(x = reorder(ensembl_transcript_id, Normalized_MFE), y = Normalized_MFE)) +
  theme_classic() +
  geom_point(stat = "identity", colour = "black", size = 1) +
  theme(axis.text.y = element_blank()) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.title = element_text(size = 20),
    axis.text.x = element_text(size = 18, color = "black"),
    plot.margin = margin(t = 50, r = 10, b = 10, l = 10),
    axis.line = element_line(linewidth = 0.2)
  ) +
  xlab("") +
  ylab("") +
  ggtitle("") +
  coord_flip()

print(Human_protein_coding_TX)

ggsave(
  filename = file.path(output_dir, "Human_protein_coding_TX.pdf"),
  plot = Human_protein_coding_TX,
  width = 5,
  height = 10
)

ora_object_names <- c(
  "gene_rank_mapped",
  "gene_representative_delta_table",
  "ora_canonical_signaling_lower10_df",
  "ora_canonical_signaling_upper90_df",
  "most_stable_10percent_entrez",
  "least_stable_90percent_entrez",
  "canonical_signaling_term2gene"
)

ora_objects <- read_named_rds_objects(
  shuffle_selected_objects_rds_file,
  ora_object_names
)

gene_rank_mapped <- ora_objects[["gene_rank_mapped"]]
gene_representative_delta_table <- ora_objects[["gene_representative_delta_table"]]
ora_canonical_signaling_lower10_df <- ora_objects[["ora_canonical_signaling_lower10_df"]]
ora_canonical_signaling_upper90_df <- ora_objects[["ora_canonical_signaling_upper90_df"]]
most_stable_10percent_entrez <- ora_objects[["most_stable_10percent_entrez"]]
least_stable_90percent_entrez <- ora_objects[["least_stable_90percent_entrez"]]
canonical_signaling_term2gene <- ora_objects[["canonical_signaling_term2gene"]]

entrez_native_mfe_annotation <- gene_rank_mapped %>%
  dplyr::select(
    ENTREZID,
    SYMBOL,
    ensembl_gene_id
  ) %>%
  dplyr::distinct() %>%
  dplyr::left_join(
    gene_representative_delta_table %>%
      dplyr::select(
        ensembl_gene_id,
        gene_label,
        external_gene_name,
        transcript_id,
        native_MFE,
        native_Normalized_MFE,
        mean_shuffled_MFE,
        mean_shuffled_Normalized_MFE,
        delta_MFE,
        delta_Normalized_MFE,
        stability_selection_score_normalized
      ),
    by = "ensembl_gene_id"
  ) %>%
  dplyr::mutate(
    ENTREZID = as.character(ENTREZID),
    gene_display = dplyr::case_when(
      !is.na(SYMBOL) & SYMBOL != "" ~ SYMBOL,
      !is.na(gene_label) & gene_label != "" ~ gene_label,
      !is.na(external_gene_name) & external_gene_name != "" ~ external_gene_name,
      TRUE ~ ensembl_gene_id
    )
  ) %>%
  dplyr::filter(
    !is.na(ENTREZID),
    ENTREZID != "",
    !is.na(native_Normalized_MFE),
    !is.na(delta_Normalized_MFE)
  ) %>%
  dplyr::distinct(ENTREZID, .keep_all = TRUE)

cat("\nEntrez-to-MFE annotation table:\n")
cat("Rows:", nrow(entrez_native_mfe_annotation), "\n")
cat(
  "Unique Entrez genes:",
  dplyr::n_distinct(entrez_native_mfe_annotation$ENTREZID),
  "\n"
)

build_ora_native_mfe_plot_df <- function(ora_df, contrast_label) {

  if (is.null(ora_df) || nrow(ora_df) == 0) {
    message("No ORA terms available for ", contrast_label)
    return(tibble::tibble())
  }

  required_cols <- c("ID", "Description", "geneID", "p.adjust")

  missing_cols <- setdiff(required_cols, colnames(ora_df))

  if (length(missing_cols) > 0) {
    stop(
      "The ORA result table for ",
      contrast_label,
      " is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  if (!"pathway_source" %in% colnames(ora_df)) {
    ora_df <- ora_df %>%
      dplyr::mutate(pathway_source = NA_character_)
  }

  ora_df %>%
    dplyr::filter(
      !is.na(geneID),
      geneID != ""
    ) %>%
    dplyr::mutate(
      ORA_contrast = contrast_label,
      ORA_contrast_label = dplyr::case_when(
        contrast_label == "Lower_10_percent_most_stable_vs_shuffle" ~
          "Lower 10%: most stable vs shuffled",
        contrast_label == "Upper_10_percent_least_stable_vs_shuffle" ~
          "Upper 10%: least stable vs shuffled",
        TRUE ~ contrast_label
      ),
      pathway_label = Description,
      pathway_label = stringr::str_replace_all(pathway_label, "_", " "),
      pathway_label = stringr::str_squish(pathway_label)
    ) %>%
    dplyr::select(
      ORA_contrast,
      ORA_contrast_label,
      pathway_id = ID,
      pathway_label,
      pathway_source,
      pvalue,
      p.adjust,
      qvalue,
      GeneRatio,
      BgRatio,
      Count,
      geneID
    ) %>%
    tidyr::separate_rows(geneID, sep = "/") %>%
    dplyr::rename(ENTREZID = geneID) %>%
    dplyr::mutate(
      ENTREZID = as.character(ENTREZID)
    ) %>%
    dplyr::left_join(
      entrez_native_mfe_annotation,
      by = "ENTREZID"
    ) %>%
    dplyr::filter(
      !is.na(native_Normalized_MFE),
      !is.na(delta_Normalized_MFE)
    ) %>%
    dplyr::distinct(
      ORA_contrast,
      pathway_id,
      ENTREZID,
      .keep_all = TRUE
    )
}

ora_lower10_native_mfe_plot_df <- build_ora_native_mfe_plot_df(
  ora_df = ora_canonical_signaling_lower10_df,
  contrast_label = "Lower_10_percent_most_stable_vs_shuffle"
)

ora_upper90_native_mfe_plot_df <- build_ora_native_mfe_plot_df(
  ora_df = ora_canonical_signaling_upper90_df,
  contrast_label = "Upper_10_percent_least_stable_vs_shuffle"
)

ora_combined_native_mfe_plot_df <- dplyr::bind_rows(
  ora_lower10_native_mfe_plot_df,
  ora_upper90_native_mfe_plot_df
)

cat("\nCombined ORA MFE plot table QC:\n")
cat(
  "Combined ORA pathway-gene rows:",
  nrow(ora_combined_native_mfe_plot_df),
  "\n"
)
cat(
  "Combined unique ORA pathways:",
  dplyr::n_distinct(
    paste(
      ora_combined_native_mfe_plot_df$ORA_contrast,
      ora_combined_native_mfe_plot_df$pathway_id,
      sep = "___"
    )
  ),
  "\n"
)
cat(
  "Combined unique Entrez genes:",
  dplyr::n_distinct(ora_combined_native_mfe_plot_df$ENTREZID),
  "\n"
)

cat("\nORA pathway-gene rows by contrast:\n")
ora_combined_native_mfe_plot_df %>%
  dplyr::count(ORA_contrast_label, name = "n_pathway_gene_rows") %>%
  tibble::as_tibble() %>%
  print(n = Inf)

prepare_ora_combined_plot_df <- function(plot_df, metric_col) {

  if (is.null(plot_df) || nrow(plot_df) == 0) {
    return(tibble::tibble())
  }

  metric_col <- rlang::ensym(metric_col)

  contrast_levels <- c(
    "Upper 10%: least stable vs shuffled",
    "Lower 10%: most stable vs shuffled"
  )

  plot_df <- plot_df %>%
    dplyr::mutate(
      ORA_contrast_label = factor(
        ORA_contrast_label,
        levels = contrast_levels
      )
    )

  pathway_summary <- plot_df %>%
    dplyr::group_by(
      ORA_contrast,
      ORA_contrast_label,
      pathway_id,
      pathway_label,
      pathway_source
    ) %>%
    dplyr::summarise(
      median_metric = median(!!metric_col, na.rm = TRUE),
      n_genes = dplyr::n_distinct(ENTREZID),
      min_ORA_FDR = min(p.adjust, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(
      ORA_contrast_label,
      median_metric
    ) %>%
    dplyr::mutate(
      pathway_axis_label = paste0(
        pathway_label,
        "\n(n=",
        n_genes,
        "; FDR=",
        signif(min_ORA_FDR, 3),
        ")"
      ),
      pathway_plot_id = paste(
        ORA_contrast,
        pathway_id,
        sep = "___"
      )
    )

  pathway_levels <- pathway_summary$pathway_plot_id

  axis_labels <- pathway_summary$pathway_axis_label
  names(axis_labels) <- pathway_summary$pathway_plot_id

  plot_df_prepared <- plot_df %>%
    dplyr::mutate(
      pathway_plot_id = paste(
        ORA_contrast,
        pathway_id,
        sep = "___"
      ),
      pathway_plot_id = factor(
        pathway_plot_id,
        levels = pathway_levels
      )
    )

  list(
    plot_df = plot_df_prepared,
    axis_labels = axis_labels,
    pathway_summary = pathway_summary
  )
}

native_plot_data <- prepare_ora_combined_plot_df(
  plot_df = ora_combined_native_mfe_plot_df,
  metric_col = native_Normalized_MFE
)

p_ora_combined_native_mfe <- ggplot(
  native_plot_data$plot_df,
  aes(
    x = native_Normalized_MFE,
    y = pathway_plot_id
  )
) +
  geom_vline(
    xintercept = -0.3,
    linetype = "dashed",
    linewidth = 0.35,
    colour = "grey60"
  ) +
  geom_boxplot(
    linewidth = 0.35,
    width = 0.6,
    fill = "white",
    outlier.shape = 16,
    outlier.size = 1.1,
    outlier.alpha = 0.6
  ) +
  facet_grid(
    ORA_contrast_label ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
  scale_y_discrete(
    labels = native_plot_data$axis_labels
  ) +
  scale_x_continuous(
    limits = c(-0.5, -0.2),
    breaks = seq(-0.5, -0.1, by = 0.1),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  theme_classic(base_size = 12) +
  labs(
    title = "Native MFE of genes in ORA-enriched signaling pathways",
    x = "Native normalized MFE (kcal·mol⁻¹·nt⁻¹)",
    y = ""
  ) +
  theme(
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 9),
    strip.text.y = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

print(p_ora_combined_native_mfe)

ggsave(
  filename = file.path(output_dir, "p_ora_combined_native_mfe.tiff"),
  plot = p_ora_combined_native_mfe,
  device = "tiff",
  width = 9,
  height = 5,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

delta_plot_data <- prepare_ora_combined_plot_df(
  plot_df = ora_combined_native_mfe_plot_df,
  metric_col = delta_Normalized_MFE
)

delta_x_min <- -0.1

delta_x_max <- max(
  delta_plot_data$plot_df$delta_Normalized_MFE,
  0,
  na.rm = TRUE
)

p_ora_combined_delta_mfe <- ggplot(
  delta_plot_data$plot_df,
  aes(
    x = delta_Normalized_MFE,
    y = pathway_plot_id
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.35,
    colour = "grey60"
  ) +
  geom_boxplot(
    linewidth = 0.35,
    width = 0.6,
    fill = "#337AB7",
    outlier.shape = 16,
    outlier.size = 1.1,
    outlier.alpha = 0.6
  ) +
  facet_grid(
    ORA_contrast_label ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
  scale_y_discrete(
    labels = delta_plot_data$axis_labels
  ) +
  scale_x_continuous(
    limits = c(delta_x_min, delta_x_max * 1.05),
    breaks = seq(-0.1, delta_x_max * 1.05, by = 0.02),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  theme_classic(base_size = 12) +
  labs(
    title = "Native-minus-shuffled MFE of genes in ORA-enriched signaling pathways",
    x = "Native - matched shuffled normalized MFE",
    y = ""
  ) +
  theme(
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 9),
    strip.text.y = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

print(p_ora_combined_delta_mfe)

ggsave(
  filename = file.path(output_dir, "p_ora_combined_delta_mfe.tiff"),
  plot = p_ora_combined_delta_mfe,
  device = "tiff",
  width = 9,
  height = 5,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

ora_combined_native_mfe_summary <- ora_combined_native_mfe_plot_df %>%
  dplyr::group_by(
    ORA_contrast,
    ORA_contrast_label,
    pathway_id,
    pathway_label,
    pathway_source
  ) %>%
  dplyr::summarise(
    n_genes_with_mfe = dplyr::n_distinct(ENTREZID),

    median_native_MFE = median(native_MFE, na.rm = TRUE),
    mean_native_MFE = mean(native_MFE, na.rm = TRUE),

    median_native_Normalized_MFE =
      median(native_Normalized_MFE, na.rm = TRUE),
    mean_native_Normalized_MFE =
      mean(native_Normalized_MFE, na.rm = TRUE),

    median_delta_MFE = median(delta_MFE, na.rm = TRUE),
    mean_delta_MFE = mean(delta_MFE, na.rm = TRUE),

    median_delta_Normalized_MFE =
      median(delta_Normalized_MFE, na.rm = TRUE),
    mean_delta_Normalized_MFE =
      mean(delta_Normalized_MFE, na.rm = TRUE),

    min_ORA_FDR = min(p.adjust, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    ORA_contrast_label,
    median_native_Normalized_MFE
  )

cat("\nNative and native-minus-shuffled MFE summary for all ORA-enriched pathways:\n")

ora_combined_native_mfe_summary %>%
  tibble::as_tibble() %>%
  print(n = Inf, width = Inf)

cat("\nCombined ORA-enriched signaling pathway MFE plots completed.\n")
cat("Generated plot objects:\n")
cat("  p_ora_combined_native_mfe\n")
cat("  p_ora_combined_delta_mfe\n")

required_objects <- c(
  "ora_combined_native_mfe_plot_df",
  "most_stable_10percent_entrez",
  "least_stable_90percent_entrez",
  "entrez_native_mfe_annotation"
)

missing_objects <- required_objects[!vapply(required_objects, exists, logical(1))]

if (length(missing_objects) > 0) {
  stop(
    "The following required objects are missing: ",
    paste(missing_objects, collapse = ", "),
    "\nRun the ORA pathway-gene MFE plotting block before this table block."
  )
}

required_columns <- c(
  "ORA_contrast",
  "ORA_contrast_label",
  "pathway_id",
  "pathway_label",
  "pathway_source",
  "ENTREZID",
  "native_MFE",
  "native_Normalized_MFE",
  "delta_MFE",
  "delta_Normalized_MFE",
  "GeneRatio",
  "BgRatio",
  "Count",
  "pvalue",
  "p.adjust",
  "qvalue"
)

missing_columns <- setdiff(
  required_columns,
  colnames(ora_combined_native_mfe_plot_df)
)

if (length(missing_columns) > 0) {
  stop(
    "ora_combined_native_mfe_plot_df is missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

safe_wilcox_one_sample <- function(x, mu = 0, alternative = "two.sided") {

  x <- x[is.finite(x)]

  if (length(x) < 2) {
    return(NA_real_)
  }

  if (all(x == mu, na.rm = TRUE)) {
    return(NA_real_)
  }

  tryCatch(
    stats::wilcox.test(
      x = x,
      mu = mu,
      alternative = alternative,
      exact = FALSE
    )$p.value,
    error = function(e) NA_real_
  )
}

safe_wilcox_two_sample <- function(x, y, alternative = "two.sided") {

  x <- x[is.finite(x)]
  y <- y[is.finite(y)]

  if (length(x) < 2 || length(y) < 2) {
    return(NA_real_)
  }

  tryCatch(
    stats::wilcox.test(
      x = x,
      y = y,
      alternative = alternative,
      exact = FALSE
    )$p.value,
    error = function(e) NA_real_
  )
}

foreground_reference_df <- entrez_native_mfe_annotation %>%
  dplyr::mutate(
    ENTREZID = as.character(ENTREZID),
    foreground_group = dplyr::case_when(
      ENTREZID %in% as.character(most_stable_10percent_entrez) ~
        "Lower_10_percent_most_stable_vs_shuffle",
      ENTREZID %in% as.character(least_stable_90percent_entrez) ~
        "Upper_10_percent_least_stable_vs_shuffle",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(foreground_group))

lower10_reference_native <- foreground_reference_df %>%
  dplyr::filter(foreground_group == "Lower_10_percent_most_stable_vs_shuffle") %>%
  dplyr::pull(native_Normalized_MFE)

upper90_reference_native <- foreground_reference_df %>%
  dplyr::filter(foreground_group == "Upper_10_percent_least_stable_vs_shuffle") %>%
  dplyr::pull(native_Normalized_MFE)

lower10_reference_delta <- foreground_reference_df %>%
  dplyr::filter(foreground_group == "Lower_10_percent_most_stable_vs_shuffle") %>%
  dplyr::pull(delta_Normalized_MFE)

upper90_reference_delta <- foreground_reference_df %>%
  dplyr::filter(foreground_group == "Upper_10_percent_least_stable_vs_shuffle") %>%
  dplyr::pull(delta_Normalized_MFE)

cat("\nForeground reference distributions:\n")
cat("Lower 10% reference genes:", length(lower10_reference_delta), "\n")
cat("Upper 10% reference genes:", length(upper90_reference_delta), "\n")

ora_pathway_manuscript_table <- ora_combined_native_mfe_plot_df %>%
  dplyr::group_by(
    ORA_contrast,
    ORA_contrast_label,
    pathway_id,
    pathway_label,
    pathway_source
  ) %>%
  dplyr::summarise(

    ORA_GeneRatio = dplyr::first(GeneRatio),
    ORA_BgRatio = dplyr::first(BgRatio),
    ORA_Count = dplyr::first(Count),
    ORA_pvalue = dplyr::first(pvalue),
    ORA_FDR = dplyr::first(p.adjust),
    ORA_qvalue = dplyr::first(qvalue),

    n_ORA_genes_with_MFE = dplyr::n_distinct(ENTREZID),

    median_native_MFE = median(native_MFE, na.rm = TRUE),
    mean_native_MFE = mean(native_MFE, na.rm = TRUE),
    sd_native_MFE = sd(native_MFE, na.rm = TRUE),

    median_native_Normalized_MFE =
      median(native_Normalized_MFE, na.rm = TRUE),
    mean_native_Normalized_MFE =
      mean(native_Normalized_MFE, na.rm = TRUE),
    sd_native_Normalized_MFE =
      sd(native_Normalized_MFE, na.rm = TRUE),

    median_delta_MFE = median(delta_MFE, na.rm = TRUE),
    mean_delta_MFE = mean(delta_MFE, na.rm = TRUE),
    sd_delta_MFE = sd(delta_MFE, na.rm = TRUE),

    median_delta_Normalized_MFE =
      median(delta_Normalized_MFE, na.rm = TRUE),
    mean_delta_Normalized_MFE =
      mean(delta_Normalized_MFE, na.rm = TRUE),
    sd_delta_Normalized_MFE =
      sd(delta_Normalized_MFE, na.rm = TRUE),

    pct_genes_native_more_stable_than_shuffle =
      100 * mean(delta_Normalized_MFE < 0, na.rm = TRUE),

    n_genes_native_more_stable_than_shuffle =
      sum(delta_Normalized_MFE < 0, na.rm = TRUE),

    n_genes_native_less_stable_than_shuffle =
      sum(delta_Normalized_MFE > 0, na.rm = TRUE),

    wilcox_delta_vs_zero_two_sided =
      safe_wilcox_one_sample(
        delta_Normalized_MFE,
        mu = 0,
        alternative = "two.sided"
      ),

    wilcox_delta_less_than_zero_native_more_stable =
      safe_wilcox_one_sample(
        delta_Normalized_MFE,
        mu = 0,
        alternative = "less"
      ),

    wilcox_delta_greater_than_zero_native_less_stable =
      safe_wilcox_one_sample(
        delta_Normalized_MFE,
        mu = 0,
        alternative = "greater"
      ),

    wilcox_native_vs_lower10_foreground_two_sided =
      safe_wilcox_two_sample(
        native_Normalized_MFE,
        lower10_reference_native,
        alternative = "two.sided"
      ),

    wilcox_native_less_than_lower10_foreground =
      safe_wilcox_two_sample(
        native_Normalized_MFE,
        lower10_reference_native,
        alternative = "less"
      ),

    wilcox_native_greater_than_lower10_foreground =
      safe_wilcox_two_sample(
        native_Normalized_MFE,
        lower10_reference_native,
        alternative = "greater"
      ),

    wilcox_delta_vs_lower10_foreground_two_sided =
      safe_wilcox_two_sample(
        delta_Normalized_MFE,
        lower10_reference_delta,
        alternative = "two.sided"
      ),

    wilcox_delta_less_than_lower10_foreground =
      safe_wilcox_two_sample(
        delta_Normalized_MFE,
        lower10_reference_delta,
        alternative = "less"
      ),

    wilcox_delta_greater_than_lower10_foreground =
      safe_wilcox_two_sample(
        delta_Normalized_MFE,
        lower10_reference_delta,
        alternative = "greater"
      ),

    wilcox_native_vs_upper90_foreground_two_sided =
      safe_wilcox_two_sample(
        native_Normalized_MFE,
        upper90_reference_native,
        alternative = "two.sided"
      ),

    wilcox_native_less_than_upper90_foreground =
      safe_wilcox_two_sample(
        native_Normalized_MFE,
        upper90_reference_native,
        alternative = "less"
      ),

    wilcox_native_greater_than_upper90_foreground =
      safe_wilcox_two_sample(
        native_Normalized_MFE,
        upper90_reference_native,
        alternative = "greater"
      ),

    wilcox_delta_vs_upper90_foreground_two_sided =
      safe_wilcox_two_sample(
        delta_Normalized_MFE,
        upper90_reference_delta,
        alternative = "two.sided"
      ),

    wilcox_delta_less_than_upper90_foreground =
      safe_wilcox_two_sample(
        delta_Normalized_MFE,
        upper90_reference_delta,
        alternative = "less"
      ),

    wilcox_delta_greater_than_upper90_foreground =
      safe_wilcox_two_sample(
        delta_Normalized_MFE,
        upper90_reference_delta,
        alternative = "greater"
      ),

    gene_symbols = paste(
      sort(unique(gene_display[!is.na(gene_display) & gene_display != ""])),
      collapse = "; "
    ),

    entrez_ids = paste(
      sort(unique(ENTREZID)),
      collapse = "; "
    ),

    .groups = "drop"
  ) %>%
  dplyr::mutate(
    ORA_direction = dplyr::case_when(
      ORA_contrast == "Lower_10_percent_most_stable_vs_shuffle" ~
        "Enriched_in_lower_10_percent_most_stable_native_vs_shuffle",
      ORA_contrast == "Upper_10_percent_least_stable_vs_shuffle" ~
        "Enriched_in_upper_10_percent_least_stable_native_vs_shuffle",
      TRUE ~ ORA_contrast
    ),

    pathway_native_vs_shuffle_direction = dplyr::case_when(
      median_delta_Normalized_MFE < 0 ~
        "Native_more_stable_than_matched_shuffles",
      median_delta_Normalized_MFE > 0 ~
        "Native_less_stable_than_matched_shuffles",
      TRUE ~
        "No_median_difference"
    )
  )

ora_pathway_manuscript_table <- ora_pathway_manuscript_table %>%
  dplyr::mutate(
    across(
      .cols = c(
        wilcox_delta_vs_zero_two_sided,
        wilcox_delta_less_than_zero_native_more_stable,
        wilcox_delta_greater_than_zero_native_less_stable,
        wilcox_native_vs_lower10_foreground_two_sided,
        wilcox_native_less_than_lower10_foreground,
        wilcox_native_greater_than_lower10_foreground,
        wilcox_delta_vs_lower10_foreground_two_sided,
        wilcox_delta_less_than_lower10_foreground,
        wilcox_delta_greater_than_lower10_foreground,
        wilcox_native_vs_upper90_foreground_two_sided,
        wilcox_native_less_than_upper90_foreground,
        wilcox_native_greater_than_upper90_foreground,
        wilcox_delta_vs_upper90_foreground_two_sided,
        wilcox_delta_less_than_upper90_foreground,
        wilcox_delta_greater_than_upper90_foreground
      ),
      .fns = ~ p.adjust(.x, method = "BH"),
      .names = "{.col}_BH"
    )
  ) %>%
  dplyr::arrange(
    ORA_contrast_label,
    ORA_FDR,
    median_delta_Normalized_MFE
  )

ora_pathway_manuscript_table_compact <- ora_pathway_manuscript_table %>%
  dplyr::select(
    ORA_direction,
    ORA_contrast_label,
    pathway_id,
    pathway_label,
    pathway_source,

    n_ORA_genes_with_MFE,
    ORA_GeneRatio,
    ORA_BgRatio,
    ORA_Count,
    ORA_pvalue,
    ORA_FDR,
    ORA_qvalue,

    median_native_Normalized_MFE,
    mean_native_Normalized_MFE,
    sd_native_Normalized_MFE,

    median_delta_Normalized_MFE,
    mean_delta_Normalized_MFE,
    sd_delta_Normalized_MFE,
    pct_genes_native_more_stable_than_shuffle,
    pathway_native_vs_shuffle_direction,

    wilcox_delta_vs_zero_two_sided,
    wilcox_delta_vs_zero_two_sided_BH,
    wilcox_delta_less_than_zero_native_more_stable,
    wilcox_delta_less_than_zero_native_more_stable_BH,
    wilcox_delta_greater_than_zero_native_less_stable,
    wilcox_delta_greater_than_zero_native_less_stable_BH,

    wilcox_native_vs_lower10_foreground_two_sided,
    wilcox_native_vs_lower10_foreground_two_sided_BH,
    wilcox_delta_vs_lower10_foreground_two_sided,
    wilcox_delta_vs_lower10_foreground_two_sided_BH,

    wilcox_native_vs_upper90_foreground_two_sided,
    wilcox_native_vs_upper90_foreground_two_sided_BH,
    wilcox_delta_vs_upper90_foreground_two_sided,
    wilcox_delta_vs_upper90_foreground_two_sided_BH,

    gene_symbols,
    entrez_ids
  )

cat("\nComplete manuscript-ready ORA pathway MFE statistics table:\n")

ora_pathway_manuscript_table %>%
  tibble::as_tibble() %>%
  print(n = Inf, width = Inf)

cat("\nCompact manuscript-ready ORA pathway MFE statistics table:\n")

ora_pathway_manuscript_table_compact %>%
  tibble::as_tibble() %>%
  print(n = Inf, width = Inf)

ora_enriched_pathways <- dplyr::bind_rows(
  ora_canonical_signaling_lower10_df %>%
    dplyr::mutate(ORA_contrast_label = "Lower 10%: most stable vs shuffled"),
  ora_canonical_signaling_upper90_df %>%
    dplyr::mutate(ORA_contrast_label = "Upper 10%: least stable vs shuffled")
) %>%
  dplyr::select(
    pathway_id = ID,
    pathway_label = Description,
    ORA_contrast_label,
    p.adjust
  ) %>%
  dplyr::distinct() %>%
  dplyr::mutate(
    pathway_label = stringr::str_replace_all(pathway_label, "_", " "),
    pathway_label = stringr::str_squish(pathway_label)
  )

ora_full_pathway_mfe_df <- ora_enriched_pathways %>%
  dplyr::left_join(
    canonical_signaling_term2gene %>%
      dplyr::rename(
        pathway_id = term,
        ENTREZID = gene
      ) %>%
      dplyr::mutate(ENTREZID = as.character(ENTREZID)),
    by = "pathway_id"
  ) %>%
  dplyr::left_join(
    entrez_native_mfe_annotation,
    by = "ENTREZID"
  ) %>%
  dplyr::filter(
    !is.na(native_Normalized_MFE),
    !is.na(delta_Normalized_MFE)
  ) %>%
  dplyr::distinct(
    ORA_contrast_label,
    pathway_id,
    ENTREZID,
    .keep_all = TRUE
  )

ora_full_pathway_mfe_df <- ora_full_pathway_mfe_df %>%
  dplyr::mutate(
    ORA_contrast_label = factor(
      ORA_contrast_label,
      levels = c(
        "Upper 10%: least stable vs shuffled",
        "Lower 10%: most stable vs shuffled"
      )
    )
  )

pathway_order_df <- ora_full_pathway_mfe_df %>%
  dplyr::group_by(ORA_contrast_label, pathway_id, pathway_label) %>%
  dplyr::summarise(
    n_genes = dplyr::n_distinct(ENTREZID),
    min_FDR = min(p.adjust, na.rm = TRUE),
    median_native = median(native_Normalized_MFE, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    ORA_contrast_label,
    median_native
  ) %>%
  dplyr::mutate(
    pathway_plot_id = paste(ORA_contrast_label, pathway_id, sep = "___"),
    axis_label = paste0(
      pathway_label,
      "\n(full pathway genes n=", n_genes,
      "; ORA FDR=", signif(min_FDR, 3),
      ")"
    )
  )

axis_labels <- pathway_order_df$axis_label
names(axis_labels) <- pathway_order_df$pathway_plot_id

ora_full_pathway_mfe_df <- ora_full_pathway_mfe_df %>%
  dplyr::mutate(
    pathway_plot_id = paste(ORA_contrast_label, pathway_id, sep = "___"),
    pathway_plot_id = factor(
      pathway_plot_id,
      levels = pathway_order_df$pathway_plot_id
    )
  )

p_full_pathway_native_mfe <- ggplot(
  ora_full_pathway_mfe_df,
  aes(
    x = native_Normalized_MFE,
    y = pathway_plot_id
  )
) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35) +
  geom_boxplot(
    outlier.shape = NA,
    linewidth = 0.35,
    width = 0.6,
    fill = "white"
  ) +
  geom_jitter(
    height = 0.16,
    width = 0,
    alpha = 0.45,
    size = 1.0
  ) +
  facet_grid(
    ORA_contrast_label ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
  scale_y_discrete(labels = axis_labels) +
  theme_classic(base_size = 12) +
  labs(
    title = "Native MFE of full ORA-enriched signaling pathways",
    subtitle = "Uses the full pathway gene set, not only ORA-overlap genes",
    x = "Native normalized MFE",
    y = ""
  ) +
  theme(
    axis.text.y = element_text(size = 8),
    strip.text.y = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

print(p_full_pathway_native_mfe)

pathway_native_rank <- ora_full_pathway_mfe_df %>%
  dplyr::group_by(ORA_contrast_label, pathway_id, pathway_label) %>%
  dplyr::summarise(
    n_genes = dplyr::n_distinct(ENTREZID),
    median_native_Normalized_MFE = median(native_Normalized_MFE, na.rm = TRUE),
    mean_native_Normalized_MFE = mean(native_Normalized_MFE, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(median_native_Normalized_MFE)) %>%
  dplyr::mutate(
    least_stable_native_rank = dplyr::row_number()
  )

olfactory_native <- ora_full_pathway_mfe_df %>%
  dplyr::filter(pathway_id == "REACTOME_OLFACTORY_SIGNALING_PATHWAY") %>%
  dplyr::pull(native_Normalized_MFE) %>%
  stats::na.omit()

other_native <- ora_full_pathway_mfe_df %>%
  dplyr::filter(pathway_id != "REACTOME_OLFACTORY_SIGNALING_PATHWAY") %>%
  dplyr::pull(native_Normalized_MFE) %>%
  stats::na.omit()

olfactory_native_less_stable_wilcox <- wilcox.test(
  olfactory_native,
  other_native,
  alternative = "greater",
  exact = FALSE
)

olfactory_native_rank_with_pvalue <- pathway_native_rank %>%
  dplyr::filter(pathway_id == "REACTOME_OLFACTORY_SIGNALING_PATHWAY") %>%
  dplyr::mutate(
    wilcox_p_native_less_stable_than_other_ORA_pathways =
      olfactory_native_less_stable_wilcox$p.value
  )

cat("\nREACTOME_OLFACTORY_SIGNALING_PATHWAY native MFE rank and p value:\n")

olfactory_native_rank_with_pvalue %>%
  tibble::as_tibble() %>%
  print(width = Inf)

cat("\nLeast stable pathway by median native_Normalized_MFE:\n")

pathway_native_rank %>%
  dplyr::slice(1) %>%
  tibble::as_tibble() %>%
  print(width = Inf)

cat(
  "\nWilcoxon p value testing whether REACTOME_OLFACTORY_SIGNALING_PATHWAY has less stable native MFE than the other ORA-enriched pathways:",
  olfactory_native_less_stable_wilcox$p.value,
  "\n"
)

p_full_pathway_delta_mfe <- ggplot(
  ora_full_pathway_mfe_df,
  aes(
    x = delta_Normalized_MFE,
    y = pathway_plot_id
  )
) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35) +
  geom_boxplot(
    outlier.shape = NA,
    linewidth = 0.35,
    width = 0.6,
    fill = "grey70"
  ) +
  geom_jitter(
    height = 0.16,
    width = 0,
    alpha = 0.45,
    size = 1.0
  ) +
  facet_grid(
    ORA_contrast_label ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
  scale_y_discrete(labels = axis_labels) +
  theme_classic(base_size = 12) +
  labs(
    title = "Native-minus-shuffled MFE of full ORA-enriched signaling pathways",
    subtitle = "Uses the full pathway gene set, not only ORA-overlap genes",
    x = "delta_Normalized_MFE = native - mean matched shuffled",
    y = ""
  ) +
  theme(
    axis.text.y = element_text(size = 8),
    strip.text.y = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

print(p_full_pathway_delta_mfe)

pathway_delta_pvalues <- ora_full_pathway_mfe_df %>%
  dplyr::group_by(ORA_contrast_label, pathway_id, pathway_label) %>%
  dplyr::summarise(
    n_genes = dplyr::n_distinct(ENTREZID),
    median_delta_Normalized_MFE = median(delta_Normalized_MFE, na.rm = TRUE),
    mean_delta_Normalized_MFE = mean(delta_Normalized_MFE, na.rm = TRUE),
    pct_delta_below_zero = 100 * mean(delta_Normalized_MFE < 0, na.rm = TRUE),
    wilcox_p_more_stable_than_shuffle = tryCatch(
      wilcox.test(
        delta_Normalized_MFE,
        mu = 0,
        alternative = "less",
        exact = FALSE
      )$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    wilcox_p_more_stable_than_shuffle_BH =
      p.adjust(wilcox_p_more_stable_than_shuffle, method = "BH")
  ) %>%
  dplyr::arrange(
    wilcox_p_more_stable_than_shuffle_BH,
    median_delta_Normalized_MFE
  )

pathway_delta_pvalues %>%
  tibble::as_tibble() %>%
  print(n = Inf, width = Inf)

pathway_delta_rank <- ora_full_pathway_mfe_df %>%
  dplyr::group_by(ORA_contrast_label, pathway_id, pathway_label) %>%
  dplyr::summarise(
    n_genes = dplyr::n_distinct(ENTREZID),
    median_delta_Normalized_MFE = median(delta_Normalized_MFE, na.rm = TRUE),
    mean_delta_Normalized_MFE = mean(delta_Normalized_MFE, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(median_delta_Normalized_MFE)) %>%
  dplyr::mutate(
    least_stable_rank = dplyr::row_number()
  )

olfactory_delta <- ora_full_pathway_mfe_df %>%
  dplyr::filter(pathway_id == "REACTOME_OLFACTORY_SIGNALING_PATHWAY") %>%
  dplyr::pull(delta_Normalized_MFE) %>%
  stats::na.omit()

other_delta <- ora_full_pathway_mfe_df %>%
  dplyr::filter(pathway_id != "REACTOME_OLFACTORY_SIGNALING_PATHWAY") %>%
  dplyr::pull(delta_Normalized_MFE) %>%
  stats::na.omit()

olfactory_less_stable_wilcox <- wilcox.test(
  olfactory_delta,
  other_delta,
  alternative = "greater",
  exact = FALSE
)

olfactory_rank_with_pvalue <- pathway_delta_rank %>%
  dplyr::filter(pathway_id == "REACTOME_OLFACTORY_SIGNALING_PATHWAY") %>%
  dplyr::mutate(
    wilcox_p_less_stable_than_other_ORA_pathways =
      olfactory_less_stable_wilcox$p.value
  )

cat("\nREACTOME_OLFACTORY_SIGNALING_PATHWAY rank and p value:\n")

olfactory_rank_with_pvalue %>%
  tibble::as_tibble() %>%
  print(width = Inf)

cat("\nLeast stable pathway by median delta_Normalized_MFE:\n")

pathway_delta_rank %>%
  dplyr::slice(1) %>%
  tibble::as_tibble() %>%
  print(width = Inf)

cat(
  "\nWilcoxon p value testing whether REACTOME_OLFACTORY_SIGNALING_PATHWAY is less stable than the other ORA-enriched pathways:",
  olfactory_less_stable_wilcox$p.value,
  "\n"
)

cat("\nFigure 1 analysis complete.\n")
cat("Outputs written to:\n", output_dir, "\n")
