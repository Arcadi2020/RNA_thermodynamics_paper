suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(tibble)
})

arguments <- commandArgs(trailingOnly = TRUE)

if (length(arguments) > 0L && length(arguments) < 4L) {
  stop(
    paste(
      "Provide all four arguments:",
      "GSE224463 RDS, GSE226643 RDS, GSE131992 RDS, and output directory."
    ),
    call. = FALSE
  )
}

Fig5_GSE224463_analysis_rds_file <- if (length(arguments) >= 4L) {
  arguments[[1]]
} else {
  "../../data/GSE224463_analysis_results_for_Fig_5.rds"
}

Fig5_GSE226643_analysis_rds_file <- if (length(arguments) >= 4L) {
  arguments[[2]]
} else {
  "../../data/GSE226643_analysis_results_for_Fig_5.rds"
}

Fig5_GSE131992_analysis_rds_file <- if (length(arguments) >= 4L) {
  arguments[[3]]
} else {
  "../../data/GSE131992_analysis_results_for_Fig_5.rds"
}

Fig5_output_dir <- if (length(arguments) >= 4L) {
  arguments[[4]]
} else {
  "Figure_5_outputs"
}

dir.create(Fig5_output_dir, recursive = TRUE, showWarnings = FALSE)

Fig5_input_files <- c(
  Fig5_GSE224463_analysis_rds_file,
  Fig5_GSE226643_analysis_rds_file,
  Fig5_GSE131992_analysis_rds_file
)

Fig5_missing_input_files <- Fig5_input_files[!file.exists(Fig5_input_files)]

if (length(Fig5_missing_input_files) > 0L) {
  stop(
    "Required Figure 5 input file(s) were not found: ",
    paste(Fig5_missing_input_files, collapse = ", "),
    call. = FALSE
  )
}

Fig5_GSE224463_analysis_objects <- readRDS(
  Fig5_GSE224463_analysis_rds_file
)

Fig5_GSE224463_objects <- Fig5_GSE224463_analysis_objects

Fig5_GSE226643_objects <- readRDS(
  Fig5_GSE226643_analysis_rds_file
)

Fig5_GSE131992_EXT_objects <- readRDS(
  Fig5_GSE131992_analysis_rds_file
)

if (
  !is.list(Fig5_GSE224463_objects) ||
  !is.list(Fig5_GSE226643_objects) ||
  !is.list(Fig5_GSE131992_EXT_objects)
) {
  stop("All Figure 5 input RDS files must contain named lists.", call. = FALSE)
}

Fig5_required_GSE224463_objects <- c(
  "contrast_manifest",
  "GSEA_bubble_data",
  "significant_DE_thermodynamic_transcripts",
  "significant_DE_distribution_statistics",
  "significant_DE_abundance_fractions",
  "weighted_sample_indices",
  "full_transcriptome_weighted_delta_summary",
  "full_transcriptome_weighted_gene_contributions",
  "WT_manifest",
  "DESeq2_normalised_counts_by_analysis"
)

Fig5_required_GSE226643_objects <- c(
  "GSE226643_GSEA_bubble_data",
  "GSE226643_N_sig_thermodynamic",
  "GSE226643_S_sig_thermodynamic",
  "GSE226643_significant_distribution_statistics",
  "GSE226643_significant_DE_abundance_fraction",
  "GSE226643_MFE_sample_indices_by_comparison",
  "GSE226643_MFE_delta_statistics",
  "GSE226643_MFE_comparison_gene_data"
)

Fig5_missing_GSE224463_objects <- setdiff(
  Fig5_required_GSE224463_objects,
  names(Fig5_GSE224463_objects)
)

Fig5_missing_GSE226643_objects <- setdiff(
  Fig5_required_GSE226643_objects,
  names(Fig5_GSE226643_objects)
)

if (length(Fig5_missing_GSE224463_objects) > 0L) {
  stop(
    "The GSE224463 RDS lacks required object(s): ",
    paste(Fig5_missing_GSE224463_objects, collapse = ", "),
    call. = FALSE
  )
}

if (length(Fig5_missing_GSE226643_objects) > 0L) {
  stop(
    "The GSE226643 RDS lacks required object(s): ",
    paste(Fig5_missing_GSE226643_objects, collapse = ", "),
    call. = FALSE
  )
}

Fig5_curated_GSEA_terms <- c(
  "DNA replication",
  "Cell division",
  "Cell cycle",
  "Ribosome biogenesis",
  "Translation",
  "Purine ribonucleotide biosynthesis",
  "Pyrimidine ribonucleotide biosynthesis",
  "Amino acid biosynthesis",
  "Peptidoglycan biosynthesis",
  "Cell wall organization",
  "Glycolysis",
  "TCA cycle",
  "Aerobic respiration",
  "Oxidative phosphorylation",
  "ATP metabolic process",
  "ATP biosynthetic process",
  "Stringent response",
  "Response to stress",
  "SOS response",
  "DNA damage response",
  "DNA repair",
  "Recombinational repair",
  "Error-prone translesion synthesis"
)

Fig5_GSE224463_long_timepoints <- c(
  90L,
  120L,
  180L,
  240L,
  300L
)

Fig5_GSE224463_long_time_labels <- paste0(
  Fig5_GSE224463_long_timepoints,
  " min"
)

Fig5_GSE224463_long_manifest <-
  Fig5_GSE224463_analysis_objects$contrast_manifest %>%
  dplyr::transmute(
    analysis_id = as.character(.data$analysis_id),
    series = as.character(.data$series),
    target_time_min = as.integer(.data$target_time_min)
  ) %>%
  dplyr::filter(.data$series == "long") %>%
  dplyr::arrange(.data$target_time_min) %>%
  dplyr::mutate(
    figure_order = dplyr::row_number(),
    target_time_label = factor(
      paste0(.data$target_time_min, " min"),
      levels = Fig5_GSE224463_long_time_labels
    )
  )

if (nrow(Fig5_GSE224463_long_manifest) != 5L ||
    !identical(
      Fig5_GSE224463_long_manifest$target_time_min,
      Fig5_GSE224463_long_timepoints
    )) {
  stop(
    "The GSE224463 final panel must contain 90, 120, 180, 240, and 300 min.",
    call. = FALSE
  )
}

Fig5_GSE224463_GSEA_raw <-
  Fig5_GSE224463_analysis_objects$GSEA_bubble_data %>%
  dplyr::filter(as.character(.data$series) == "long")

Fig5_missing_GSE224463_GSEA_terms <- setdiff(
  Fig5_curated_GSEA_terms,
  unique(as.character(Fig5_GSE224463_GSEA_raw$term))
)

if (length(Fig5_missing_GSE224463_GSEA_terms) > 0L) {
  stop(
    "The GSE224463 long-series curated GSEA panel is missing: ",
    paste(Fig5_missing_GSE224463_GSEA_terms, collapse = ", "),
    call. = FALSE
  )
}

Fig5_GSE224463_long_GSEA_plot_data <-
  Fig5_GSE224463_GSEA_raw %>%
  dplyr::transmute(
    analysis_id = as.character(.data$analysis_id),
    term = as.character(.data$term),
    GSEA_tested = as.logical(.data$GSEA_tested),
    NES = as.numeric(.data$NES),
    p_adjust_long_timecourse = as.numeric(
      .data$p_adjust_across_timepoints_within_series
    )
  ) %>%
  dplyr::filter(.data$term %in% Fig5_curated_GSEA_terms) %>%
  dplyr::inner_join(
    Fig5_GSE224463_long_manifest %>%
      dplyr::select(
        .data$analysis_id,
        .data$target_time_label,
        .data$figure_order
      ),
    by = "analysis_id"
  ) %>%
  dplyr::mutate(
    term = factor(
      .data$term,
      levels = rev(Fig5_curated_GSEA_terms)
    ),
    minus_log10_BH = -log10(
      pmax(
        .data$p_adjust_long_timecourse,
        .Machine$double.xmin
      )
    )
  ) %>%
  dplyr::arrange(
    .data$figure_order,
    .data$term
  )

if (nrow(Fig5_GSE224463_long_GSEA_plot_data) !=
    5L * length(Fig5_curated_GSEA_terms)) {
  stop(
    "The GSE224463 long-series GSEA table is incomplete.",
    call. = FALSE
  )
}

Fig5_GSE224463_long_GSEA <- ggplot(
  Fig5_GSE224463_long_GSEA_plot_data %>%
    dplyr::filter(.data$GSEA_tested),
  aes(
    x = .data$target_time_label,
    y = .data$term,
    size = .data$minus_log10_BH,
    colour = .data$NES
  )
) +
  geom_point(alpha = 0.95) +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  scale_colour_gradient2(
    low = "#2166ac",
    mid = "white",
    high = "#b2182b",
    midpoint = 0,
    name = "NES"
  ) +
  scale_size_continuous(
    range = c(0.5, 6),
    name = expression(-log[10]~"BH-adjusted p")
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      colour = "black",
      size = 8
    ),
    axis.text.y = element_text(
      colour = "black",
      size = 7
    )
  )

print(Fig5_GSE224463_long_GSEA)

ggsave(
  filename = file.path(Fig5_output_dir, "Fig5_GSE224463_long_GSEA.tiff"),
  plot = Fig5_GSE224463_long_GSEA,
  width = 5,
  height = 5,
  units = "in",
  dpi = 600,
  device = "tiff",
  compression = "lzw",
  bg = "white"
)

Fig5_GSE226643_GSEA_raw <-
  Fig5_GSE226643_objects$GSE226643_GSEA_bubble_data

Fig5_GSE226643_limitation_order <- c(
  "Nitrogen limitation",
  "Sulfur limitation"
)

Fig5_missing_GSE226643_GSEA_terms <- setdiff(
  Fig5_curated_GSEA_terms,
  unique(as.character(Fig5_GSE226643_GSEA_raw$term))
)

if (length(Fig5_missing_GSE226643_GSEA_terms) > 0L) {
  stop(
    "The GSE226643 curated GSEA panel is missing: ",
    paste(Fig5_missing_GSE226643_GSEA_terms, collapse = ", "),
    call. = FALSE
  )
}

Fig5_GSE226643_GSEA_plot_data <- Fig5_GSE226643_GSEA_raw %>%
  dplyr::transmute(
    limitation = as.character(.data$limitation),
    term = as.character(.data$term),
    NES = as.numeric(.data$NES),
    p_adjust = as.numeric(.data$p_adjust)
  ) %>%
  dplyr::filter(.data$term %in% Fig5_curated_GSEA_terms) %>%
  tidyr::complete(
    limitation = Fig5_GSE226643_limitation_order,
    term = Fig5_curated_GSEA_terms,
    fill = list(
      NES = 0,
      p_adjust = 1
    )
  ) %>%
  dplyr::mutate(
    limitation = factor(
      .data$limitation,
      levels = Fig5_GSE226643_limitation_order
    ),
    term = factor(
      .data$term,
      levels = rev(Fig5_curated_GSEA_terms)
    ),
    minus_log10_BH = -log10(
      pmax(.data$p_adjust, .Machine$double.xmin)
    )
  )

Fig5_GSEA_GSE226643_Nitrogen_and_Sulfur_12h_vs_3h <- ggplot(
  Fig5_GSE226643_GSEA_plot_data,
  aes(
    x = .data$limitation,
    y = .data$term,
    size = .data$minus_log10_BH,
    colour = .data$NES
  )
) +
  geom_point(alpha = 0.95) +
  scale_x_discrete(
    labels = c(
      "Nitrogen limitation" = "Nitrogen\nlimitation",
      "Sulfur limitation" = "Sulfur\nlimitation"
    )
  ) +
  scale_y_discrete(drop = FALSE) +
  scale_colour_gradient2(
    low = "#2166ac",
    mid = "white",
    high = "#b2182b",
    midpoint = 0,
    name = "NES"
  ) +
  scale_size_continuous(
    range = c(0.2, 5.5),
    name = expression(-log[10]~"BH-adjusted p")
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(
      size = 9,
      colour = "black"
    ),
    axis.text.y = element_text(
      size = 7,
      colour = "black"
    )
  )

print(Fig5_GSEA_GSE226643_Nitrogen_and_Sulfur_12h_vs_3h)

ggsave(
  filename = file.path(
    Fig5_output_dir,
    "Fig5_GSEA_GSE226643_Nitrogen_and_Sulfur_12h_vs_3h.tiff"
  ),
  plot = Fig5_GSEA_GSE226643_Nitrogen_and_Sulfur_12h_vs_3h,
  width = 4.5,
  height = 5,
  units = "in",
  dpi = 600,
  device = "tiff",
  compression = "lzw",
  bg = "white"
)

Fig5_GSE224463_long_timepoints <- c(
  90L,
  120L,
  180L,
  240L,
  300L
)

Fig5_GSE224463_long_manifest <-
  Fig5_GSE224463_objects$contrast_manifest %>%
  dplyr::transmute(
    analysis_id = as.character(.data$analysis_id),
    series = as.character(.data$series),
    target_time_min = as.integer(.data$target_time_min)
  ) %>%
  dplyr::filter(
    .data$series == "long",
    .data$target_time_min %in% Fig5_GSE224463_long_timepoints
  ) %>%
  dplyr::arrange(.data$target_time_min) %>%
  dplyr::mutate(
    figure_order = dplyr::row_number()
  )

if (nrow(Fig5_GSE224463_long_manifest) != 5L) {
  stop(
    "Expected the five long-series time points: 90, 120, 180, 240, and 300 min.",
    call. = FALSE
  )
}

Fig5_GSE224463_MFE_statistics <-
  Fig5_GSE224463_objects$significant_DE_distribution_statistics %>%
  dplyr::filter(
    as.character(.data$series) == "long",
    .data$descriptor == "Normalized_MFE"
  ) %>%
  dplyr::transmute(
    analysis_id = as.character(.data$analysis_id),
    N_target_lower = as.integer(.data$N_target_lower),
    N_target_higher = as.integer(.data$N_target_higher),
    delta_target_higher_minus_lower =
      as.numeric(.data$delta_target_higher_minus_lower),
    Wilcoxon_p = as.numeric(.data$Wilcoxon_p)
  ) %>%
  dplyr::inner_join(
    Fig5_GSE224463_long_manifest,
    by = "analysis_id"
  ) %>%
  dplyr::arrange(.data$figure_order) %>%
  dplyr::mutate(
    Wilcoxon_p_BH_long_series = stats::p.adjust(
      .data$Wilcoxon_p,
      method = "BH"
    ),
    panel_label = paste0(
      .data$target_time_min,
      " min\nΔNormalized MFE = ",
      format(
        round(.data$delta_target_higher_minus_lower, 4),
        nsmall = 4
      ),
      " | BH P = ",
      formatC(
        .data$Wilcoxon_p_BH_long_series,
        format = "e",
        digits = 2
      )
    ),
    panel_label = factor(
      .data$panel_label,
      levels = unique(.data$panel_label)
    )
  )

if (nrow(Fig5_GSE224463_MFE_statistics) != 5L) {
  stop(
    "Expected one ΔGfold statistic for each GSE224463 long-series time point.",
    call. = FALSE
  )
}

Fig5_GSE224463_MFE_violin_data <-
  Fig5_GSE224463_objects$significant_DE_thermodynamic_transcripts %>%
  dplyr::filter(
    as.character(.data$series) == "long"
  ) %>%
  dplyr::transmute(
    analysis_id = as.character(.data$analysis_id),
    FC_sign = factor(
      as.character(.data$FC_sign),
      levels = c("log2_FC < 0", "log2_FC > 0")
    ),
    Normalized_MFE = as.numeric(.data$Normalized_MFE)
  ) %>%
  dplyr::filter(
    !is.na(.data$FC_sign),
    is.finite(.data$Normalized_MFE)
  ) %>%
  dplyr::inner_join(
    Fig5_GSE224463_MFE_statistics %>%
      dplyr::select(
        .data$analysis_id,
        .data$panel_label,
        .data$figure_order
      ),
    by = "analysis_id"
  )

Fig5_GSE224463_MFE_y_limits <- c(
  min(
    Fig5_GSE224463_MFE_violin_data$Normalized_MFE,
    na.rm = TRUE
  ),
  0
)

Fig5_Violin_GSE224463_long_series_plots <- list()

for (i in seq_len(nrow(Fig5_GSE224463_MFE_statistics))) {

  current_analysis_id <-
    Fig5_GSE224463_MFE_statistics$analysis_id[i]

  current_time_min <-
    Fig5_GSE224463_MFE_statistics$target_time_min[i]

  current_panel_label <-
    as.character(Fig5_GSE224463_MFE_statistics$panel_label[i])

  current_violin_data <-
    Fig5_GSE224463_MFE_violin_data %>%
    dplyr::filter(
      .data$analysis_id == current_analysis_id
    )

  current_plot <- ggplot(
    current_violin_data,
    aes(
      x = .data$FC_sign,
      y = .data$Normalized_MFE
    )
  ) +
    geom_violin(
      trim = FALSE,
      fill = "lightgrey",
      alpha = 1,
      colour = "black",
      linewidth = 0.25
    ) +
    geom_boxplot(
      width = 0.14,
      outlier.shape = NA,
      fill = "white",
      linewidth = 0.3
    ) +
    coord_cartesian(
      ylim = Fig5_GSE224463_MFE_y_limits,
      expand = FALSE
    ) +
    scale_x_discrete(
      labels = c(
        "log2_FC < 0" = "Lower\nat target time",
        "log2_FC > 0" = "Higher\nat target time"
      )
    ) +
    theme_classic(base_size = 11) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_text(
        colour = "black",
        size = 8
      ),
      axis.text.y = element_text(
        colour = "black"
      ),
      plot.title = element_text(
        hjust = 0.5,
        size = 10
      )
    ) +
    labs(
      title = current_panel_label,
      y = "Normalized ΔGfold (ΔGfold / nt)"
    )

  Fig5_Violin_GSE224463_long_series_plots[[as.character(current_time_min)]] <-
    current_plot

  print(current_plot)

  ggsave(
    filename = file.path(
      Fig5_output_dir,
      paste0(
        "Fig5_GSE224463_long_series_significant_DE_DGfold_violin_",
        current_time_min,
        "min.tiff"
      )
    ),
    plot = current_plot,
    width = 3.2,
    height = 5,
    units = "in",
    dpi = 600,
    device = "tiff",
    compression = "lzw",
    bg = "white"
  )
}

Fig5B_GSE224463_complete_statistics <- list()

for (current_time_min in Fig5_GSE224463_long_timepoints) {

  current_analysis_id <-
    Fig5_GSE224463_MFE_statistics %>%
    dplyr::filter(
      target_time_min == current_time_min
    ) %>%
    dplyr::pull(analysis_id)

  current_data <-
    Fig5_GSE224463_MFE_violin_data %>%
    dplyr::filter(
      analysis_id == current_analysis_id
    )

  MFE_UP <-
    current_data %>%
    dplyr::filter(
      FC_sign == "log2_FC > 0"
    )

  MFE_DOWN <-
    current_data %>%
    dplyr::filter(
      FC_sign == "log2_FC < 0"
    )

  Mean_MFE_UP <-
    mean(
      MFE_UP$Normalized_MFE,
      na.rm = TRUE
    )

  Mean_MFE_DOWN <-
    mean(
      MFE_DOWN$Normalized_MFE,
      na.rm = TRUE
    )

  TTEST_all_UP_vs_DOWN <- t.test(
    MFE_DOWN$Normalized_MFE,
    MFE_UP$Normalized_MFE,
    paired = FALSE
  )

  KS_all_UP_vs_DOWN <- ks.test(
    MFE_DOWN$Normalized_MFE,
    MFE_UP$Normalized_MFE
  )

  quantile_075_MFE_UP <- quantile(
    MFE_UP$Normalized_MFE,
    probs = 0.75,
    na.rm = TRUE
  )

  quantile_075_MFE_DOWN <- quantile(
    MFE_DOWN$Normalized_MFE,
    probs = 0.75,
    na.rm = TRUE
  )

  KS_075_above <- ks.test(
    MFE_UP$Normalized_MFE[
      MFE_UP$Normalized_MFE >= quantile_075_MFE_UP
    ],
    MFE_DOWN$Normalized_MFE[
      MFE_DOWN$Normalized_MFE >= quantile_075_MFE_DOWN
    ]
  )

  quantile_025_MFE_UP <- quantile(
    MFE_UP$Normalized_MFE,
    probs = 0.25,
    na.rm = TRUE
  )

  quantile_025_MFE_DOWN <- quantile(
    MFE_DOWN$Normalized_MFE,
    probs = 0.25,
    na.rm = TRUE
  )

  KS_025_below <- ks.test(
    MFE_UP$Normalized_MFE[
      MFE_UP$Normalized_MFE <= quantile_025_MFE_UP
    ],
    MFE_DOWN$Normalized_MFE[
      MFE_DOWN$Normalized_MFE <= quantile_025_MFE_DOWN
    ]
  )

  Wilcoxon_all_UP_vs_DOWN <- wilcox.test(
    MFE_DOWN$Normalized_MFE,
    MFE_UP$Normalized_MFE,
    paired = FALSE,
    alternative = "two.sided"
  )

  Fig5B_GSE224463_complete_statistics[[as.character(current_time_min)]] <- list(

    time_min = current_time_min,

    N_UP = nrow(MFE_UP),
    N_DOWN = nrow(MFE_DOWN),

    Mean_MFE_UP = Mean_MFE_UP,
    Mean_MFE_DOWN = Mean_MFE_DOWN,

    TTEST = TTEST_all_UP_vs_DOWN,

    KS_all = KS_all_UP_vs_DOWN,

    quantile_075_MFE_UP = quantile_075_MFE_UP,
    quantile_075_MFE_DOWN = quantile_075_MFE_DOWN,
    KS_075_above = KS_075_above,

    quantile_025_MFE_UP = quantile_025_MFE_UP,
    quantile_025_MFE_DOWN = quantile_025_MFE_DOWN,
    KS_025_below = KS_025_below,

    Wilcoxon = Wilcoxon_all_UP_vs_DOWN
  )
}

Fig5B_GSE224463_stats <- dplyr::bind_rows(
  lapply(
    Fig5_GSE224463_long_timepoints,
    function(x) {

      s <- Fig5B_GSE224463_complete_statistics[[as.character(x)]]

      data.frame(
        Time_min = x,
        N_DOWN = s$N_DOWN,
        N_UP = s$N_UP,
        Mean_MFE_DOWN = s$Mean_MFE_DOWN,
        Mean_MFE_UP = s$Mean_MFE_UP,
        T_test_P = s$TTEST$p.value,
        KS_all_P = s$KS_all$p.value,
        KS_075_P = s$KS_075_above$p.value,
        KS_025_P = s$KS_025_below$p.value,
        Wilcoxon_P = s$Wilcoxon$p.value
      )
    }
  )
)

print(Fig5B_GSE224463_stats)

Fig5_GSE224463_long_MFE_timecourse_data <-
  Fig5_GSE224463_MFE_statistics %>%
  dplyr::transmute(
    target_time_min = .data$target_time_min,
    target_time_factor = factor(
      .data$target_time_min,
      levels = Fig5_GSE224463_long_timepoints
    ),
    delta_target_higher_minus_lower =
      .data$delta_target_higher_minus_lower,
    FDR_status = dplyr::if_else(
      .data$Wilcoxon_p_BH_long_series < 0.05,
      "BH FDR < 0.05",
      "Not BH-significant"
    )
  ) %>%
  dplyr::arrange(.data$target_time_min)

Fig5_GSE224463_long_MFE_timecourse_y_limit <- max(
  abs(
    Fig5_GSE224463_long_MFE_timecourse_data$
      delta_target_higher_minus_lower
  ),
  na.rm = TRUE
) * 1.15

Fig5_TimeCourse_GSE224463_LongSeries_significant_DE_DGfold <- ggplot(
  Fig5_GSE224463_long_MFE_timecourse_data,
  aes(
    x = .data$target_time_factor,
    y = .data$delta_target_higher_minus_lower,
    group = 1
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.4,
    colour = "grey40"
  ) +
  geom_line(
    colour = "grey10",
    linewidth = 0.7
  ) +
  geom_point(
    aes(shape = .data$FDR_status),
    colour = "grey10",
    size = 2
  ) +
  scale_shape_manual(
    values = c(
      "BH FDR < 0.05" = 16,
      "Not BH-significant" = 1
    ),
    name = "Directional\nWilcoxon test"
  ) +
  scale_x_discrete(
    labels = paste0(Fig5_GSE224463_long_timepoints, " min")
  ) +
  coord_cartesian(
    ylim = c(
      -Fig5_GSE224463_long_MFE_timecourse_y_limit,
      Fig5_GSE224463_long_MFE_timecourse_y_limit
    )
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black")
  ) +
  labs(
    x = "Time after exponential-phase baseline",
    y = "ΔNormalized ΔGfold\n(higher minus lower, ΔGfold / nt)"
  )

print(Fig5_TimeCourse_GSE224463_LongSeries_significant_DE_DGfold)

ggsave(
  filename = file.path(
    Fig5_output_dir,
    "Fig5_GSE224463_long_series_significant_DE_DGfold_time_course.tiff"
  ),
  plot = Fig5_TimeCourse_GSE224463_LongSeries_significant_DE_DGfold,
  width = 3,
  height = 5,
  units = "in",
  dpi = 600,
  device = "tiff",
  compression = "lzw",
  bg = "white"
)

Fig5C_Spearman_trend <- cor.test(
  Fig5_GSE224463_long_MFE_timecourse_data$target_time_min,
  Fig5_GSE224463_long_MFE_timecourse_data$delta_target_higher_minus_lower,
  method = "spearman",
  alternative = "greater",
  exact = TRUE
)

print(Fig5C_Spearman_trend)

Fig5C_linear_trend <- lm(
  delta_target_higher_minus_lower ~ target_time_min,
  data = Fig5_GSE224463_long_MFE_timecourse_data
)

summary(Fig5C_linear_trend)

Fig5C_trend_statistics <- data.frame(
  Spearman_rho = unname(Fig5C_Spearman_trend$estimate),
  Spearman_P = Fig5C_Spearman_trend$p.value,
  Linear_slope_per_min = coef(Fig5C_linear_trend)[["target_time_min"]],
  Linear_slope_per_hour = coef(Fig5C_linear_trend)[["target_time_min"]] * 60,
  Linear_P_two_sided = summary(Fig5C_linear_trend)$coefficients[
    "target_time_min", "Pr(>|t|)"
  ]
)

print(Fig5C_trend_statistics)

Fig5C_P <- Fig5C_Spearman_trend$p.value

Fig5C_report <- data.frame(
  n_timepoints = nrow(Fig5_GSE224463_long_MFE_timecourse_data),
  
  Spearman_rho = round(
    unname(Fig5C_Spearman_trend$estimate),
    2
  ),
  
  P_value = Fig5C_P,
  
  P_report = ifelse(
    Fig5C_P < 0.0001,
    "P < 0.0001",
    paste0(
      "P = ",
      trimws(base::format.pval(Fig5C_P, digits = 4))
    )
  ),
  
  Test = "One-sided exact Spearman rank correlation"
)

print(Fig5C_report, row.names = FALSE)


Fig5_GSE226643_limitation_order <- c(
  "Nitrogen limitation",
  "Sulfur limitation"
)

Fig5_GSE226643_MFE_statistics <-
  Fig5_GSE226643_objects$GSE226643_significant_distribution_statistics %>%
  dplyr::filter(.data$descriptor == "Normalized_MFE") %>%
  dplyr::mutate(
    limitation = factor(
      as.character(.data$limitation),
      levels = Fig5_GSE226643_limitation_order
    ),
    panel_label = paste0(
      as.character(.data$limitation),
      "\nΔNormalized MFE = ",
      format(
        round(.data$delta_12h_higher_minus_lower, 4),
        nsmall = 4
      ),
      " | BH P = ",
      formatC(
        .data$Wilcoxon_p_BH,
        format = "e",
        digits = 2
      )
    )
  ) %>%
  dplyr::arrange(.data$limitation)

if (nrow(Fig5_GSE226643_MFE_statistics) != 2L) {
  stop(
    "Expected exactly two GSE226643 ΔGfold violin panels.",
    call. = FALSE
  )
}

Fig5_GSE226643_MFE_violin_data <- dplyr::bind_rows(
  Fig5_GSE226643_objects$GSE226643_N_sig_thermodynamic,
  Fig5_GSE226643_objects$GSE226643_S_sig_thermodynamic
) %>%
  dplyr::transmute(
    limitation = as.character(.data$limitation),
    FC_sign = factor(
      as.character(.data$FC_sign),
      levels = c("log2_FC < 0", "log2_FC > 0")
    ),
    Normalized_MFE = as.numeric(.data$Normalized_MFE)
  ) %>%
  dplyr::filter(
    !is.na(.data$FC_sign),
    is.finite(.data$Normalized_MFE)
  ) %>%
  dplyr::left_join(
    Fig5_GSE226643_MFE_statistics %>%
      dplyr::transmute(
        limitation = as.character(.data$limitation),
        panel_label = .data$panel_label
      ),
    by = "limitation"
  ) %>%
  dplyr::mutate(
    panel_label = factor(
      .data$panel_label,
      levels = Fig5_GSE226643_MFE_statistics$panel_label
    )
  )

Fig5_GSE226643_MFE_y_limits <- c(
  min(
    Fig5_GSE226643_MFE_violin_data$Normalized_MFE,
    na.rm = TRUE
  ),
  0
)

Fig5_GSE226643_Nitrogen_violin_data <-
  Fig5_GSE226643_MFE_violin_data %>%
  dplyr::filter(
    .data$limitation == "Nitrogen limitation"
  )

Fig5_Violin_GSE226643_Nitrogen <- ggplot(
  Fig5_GSE226643_Nitrogen_violin_data,
  aes(
    x = .data$FC_sign,
    y = .data$Normalized_MFE
  )
) +
  geom_violin(
    trim = FALSE,
    fill = "lightgrey",
    alpha = 1,
    colour = "black",
    linewidth = 0.25
  ) +
  geom_boxplot(
    width = 0.14,
    outlier.shape = NA,
    fill = "white",
    linewidth = 0.3
  ) +
  coord_cartesian(
    ylim = Fig5_GSE226643_MFE_y_limits,
    expand = FALSE
  ) +
  scale_x_discrete(
    labels = c(
      "log2_FC < 0" = "Lower\nat 12 h",
      "log2_FC > 0" = "Higher\nat 12 h"
    )
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(
      colour = "black",
      size = 9
    ),
    axis.text.y = element_text(
      colour = "black"
    )
  ) +
  labs(
    y = "Normalized ΔGfold (ΔGfold / nt)"
  )

print(Fig5_Violin_GSE226643_Nitrogen)

ggsave(
  filename = file.path(Fig5_output_dir, "Fig5_Violin_GSE226643_Nitrogen.tiff"),
  plot = Fig5_Violin_GSE226643_Nitrogen,
  width = 4,
  height = 5,
  units = "in",
  dpi = 600,
  device = "tiff",
  compression = "lzw",
  bg = "white"
)

Fig5_GSE226643_Sulfur_violin_data <-
  Fig5_GSE226643_MFE_violin_data[
    Fig5_GSE226643_MFE_violin_data$limitation ==
      "Sulfur limitation",
    ,
    drop = FALSE
  ] %>%
  tibble::as_tibble()

stopifnot(nrow(Fig5_GSE226643_Sulfur_violin_data) == 297L)

Fig5_Violin_GSE226643_Sulfur <- ggplot(
  Fig5_GSE226643_Sulfur_violin_data,
  aes(
    x = .data$FC_sign,
    y = .data$Normalized_MFE
  )
) +
  geom_violin(
    trim = FALSE,
    fill = "lightgrey",
    alpha = 1,
    colour = "black",
    linewidth = 0.25
  ) +
  geom_boxplot(
    width = 0.14,
    outlier.shape = NA,
    fill = "white",
    linewidth = 0.3
  ) +
  coord_cartesian(
    ylim = Fig5_GSE226643_MFE_y_limits,
    expand = FALSE
  ) +
  scale_x_discrete(
    labels = c(
      "log2_FC < 0" = "Lower\nat 12 h",
      "log2_FC > 0" = "Higher\nat 12 h"
    )
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(
      colour = "black",
      size = 9
    ),
    axis.text.y = element_text(
      colour = "black"
    )
  ) +
  labs(
    y = "Normalized ΔGfold (ΔGfold / nt)"
  )

print(Fig5_Violin_GSE226643_Sulfur)

ggsave(
  filename = file.path(Fig5_output_dir, "Fig5_Violin_GSE226643_Sulfur.tiff"),
  plot = Fig5_Violin_GSE226643_Sulfur,
  width = 4,
  height = 5,
  units = "in",
  dpi = 600,
  device = "tiff",
  compression = "lzw",
  bg = "white"
)

MFE_GSE226643_Nitrogen_UP <-
  dplyr::filter(
    Fig5_GSE226643_Nitrogen_violin_data,
    FC_sign == "log2_FC > 0"
  )

Mean_MFE_GSE226643_Nitrogen_UP <-
  mean(MFE_GSE226643_Nitrogen_UP$Normalized_MFE, na.rm = TRUE)

MFE_GSE226643_Nitrogen_DOWN <-
  dplyr::filter(
    Fig5_GSE226643_Nitrogen_violin_data,
    FC_sign == "log2_FC < 0"
  )

Mean_MFE_GSE226643_Nitrogen_DOWN <-
  mean(MFE_GSE226643_Nitrogen_DOWN$Normalized_MFE, na.rm = TRUE)

TTEST_MFE_GSE226643_Nitrogen_all_UP_vs_DOWN <- t.test(
  MFE_GSE226643_Nitrogen_DOWN$Normalized_MFE,
  MFE_GSE226643_Nitrogen_UP$Normalized_MFE,
  paired = FALSE
)

print(TTEST_MFE_GSE226643_Nitrogen_all_UP_vs_DOWN)

ks_test_GSE226643_Nitrogen_all_UP_vs_DOWN <- ks.test(
  MFE_GSE226643_Nitrogen_DOWN$Normalized_MFE,
  MFE_GSE226643_Nitrogen_UP$Normalized_MFE
)

print(ks_test_GSE226643_Nitrogen_all_UP_vs_DOWN)

quantile_075_MFE_GSE226643_Nitrogen_UP <- quantile(
  MFE_GSE226643_Nitrogen_UP$Normalized_MFE,
  probs = 0.75,
  na.rm = TRUE
)

quantile_075_MFE_GSE226643_Nitrogen_DOWN <- quantile(
  MFE_GSE226643_Nitrogen_DOWN$Normalized_MFE,
  probs = 0.75,
  na.rm = TRUE
)

ks_test_MFE_GSE226643_Nitrogen_075_above <- ks.test(
  MFE_GSE226643_Nitrogen_UP$Normalized_MFE[
    MFE_GSE226643_Nitrogen_UP$Normalized_MFE >=
      quantile_075_MFE_GSE226643_Nitrogen_UP
  ],
  MFE_GSE226643_Nitrogen_DOWN$Normalized_MFE[
    MFE_GSE226643_Nitrogen_DOWN$Normalized_MFE >=
      quantile_075_MFE_GSE226643_Nitrogen_DOWN
  ]
)

print(ks_test_MFE_GSE226643_Nitrogen_075_above)

quantile_025_MFE_GSE226643_Nitrogen_UP <- quantile(
  MFE_GSE226643_Nitrogen_UP$Normalized_MFE,
  probs = 0.25,
  na.rm = TRUE
)

quantile_025_MFE_GSE226643_Nitrogen_DOWN <- quantile(
  MFE_GSE226643_Nitrogen_DOWN$Normalized_MFE,
  probs = 0.25,
  na.rm = TRUE
)

ks_test_MFE_GSE226643_Nitrogen_025_below <- ks.test(
  MFE_GSE226643_Nitrogen_UP$Normalized_MFE[
    MFE_GSE226643_Nitrogen_UP$Normalized_MFE <=
      quantile_025_MFE_GSE226643_Nitrogen_UP
  ],
  MFE_GSE226643_Nitrogen_DOWN$Normalized_MFE[
    MFE_GSE226643_Nitrogen_DOWN$Normalized_MFE <=
      quantile_025_MFE_GSE226643_Nitrogen_DOWN
  ]
)

print(ks_test_MFE_GSE226643_Nitrogen_025_below)

wilcox_test_GSE226643_Nitrogen_all_UP_vs_DOWN <- wilcox.test(
  MFE_GSE226643_Nitrogen_DOWN$Normalized_MFE,
  MFE_GSE226643_Nitrogen_UP$Normalized_MFE,
  paired = FALSE,
  alternative = "two.sided"
)

print(wilcox_test_GSE226643_Nitrogen_all_UP_vs_DOWN)

MFE_GSE226643_Sulfur_UP <-
  dplyr::filter(
    Fig5_GSE226643_Sulfur_violin_data,
    FC_sign == "log2_FC > 0"
  )

Mean_MFE_GSE226643_Sulfur_UP <-
  mean(MFE_GSE226643_Sulfur_UP$Normalized_MFE, na.rm = TRUE)

MFE_GSE226643_Sulfur_DOWN <-
  dplyr::filter(
    Fig5_GSE226643_Sulfur_violin_data,
    FC_sign == "log2_FC < 0"
  )

Mean_MFE_GSE226643_Sulfur_DOWN <-
  mean(MFE_GSE226643_Sulfur_DOWN$Normalized_MFE, na.rm = TRUE)

TTEST_MFE_GSE226643_Sulfur_all_UP_vs_DOWN <- t.test(
  MFE_GSE226643_Sulfur_DOWN$Normalized_MFE,
  MFE_GSE226643_Sulfur_UP$Normalized_MFE,
  paired = FALSE
)

print(TTEST_MFE_GSE226643_Sulfur_all_UP_vs_DOWN)

ks_test_GSE226643_Sulfur_all_UP_vs_DOWN <- ks.test(
  MFE_GSE226643_Sulfur_DOWN$Normalized_MFE,
  MFE_GSE226643_Sulfur_UP$Normalized_MFE
)

print(ks_test_GSE226643_Sulfur_all_UP_vs_DOWN)

quantile_075_MFE_GSE226643_Sulfur_UP <- quantile(
  MFE_GSE226643_Sulfur_UP$Normalized_MFE,
  probs = 0.75,
  na.rm = TRUE
)

quantile_075_MFE_GSE226643_Sulfur_DOWN <- quantile(
  MFE_GSE226643_Sulfur_DOWN$Normalized_MFE,
  probs = 0.75,
  na.rm = TRUE
)

ks_test_MFE_GSE226643_Sulfur_075_above <- ks.test(
  MFE_GSE226643_Sulfur_UP$Normalized_MFE[
    MFE_GSE226643_Sulfur_UP$Normalized_MFE >=
      quantile_075_MFE_GSE226643_Sulfur_UP
  ],
  MFE_GSE226643_Sulfur_DOWN$Normalized_MFE[
    MFE_GSE226643_Sulfur_DOWN$Normalized_MFE >=
      quantile_075_MFE_GSE226643_Sulfur_DOWN
  ]
)

print(ks_test_MFE_GSE226643_Sulfur_075_above)

quantile_025_MFE_GSE226643_Sulfur_UP <- quantile(
  MFE_GSE226643_Sulfur_UP$Normalized_MFE,
  probs = 0.25,
  na.rm = TRUE
)

quantile_025_MFE_GSE226643_Sulfur_DOWN <- quantile(
  MFE_GSE226643_Sulfur_DOWN$Normalized_MFE,
  probs = 0.25,
  na.rm = TRUE
)

ks_test_MFE_GSE226643_Sulfur_025_below <- ks.test(
  MFE_GSE226643_Sulfur_UP$Normalized_MFE[
    MFE_GSE226643_Sulfur_UP$Normalized_MFE <=
      quantile_025_MFE_GSE226643_Sulfur_UP
  ],
  MFE_GSE226643_Sulfur_DOWN$Normalized_MFE[
    MFE_GSE226643_Sulfur_DOWN$Normalized_MFE <=
      quantile_025_MFE_GSE226643_Sulfur_DOWN
  ]
)

print(ks_test_MFE_GSE226643_Sulfur_025_below)

wilcox_test_GSE226643_Sulfur_all_UP_vs_DOWN <- wilcox.test(
  MFE_GSE226643_Sulfur_DOWN$Normalized_MFE,
  MFE_GSE226643_Sulfur_UP$Normalized_MFE,
  paired = FALSE,
  alternative = "two.sided"
)

print(wilcox_test_GSE226643_Sulfur_all_UP_vs_DOWN)

Fig5D_GSE226643_stats <- data.frame(
  Limitation = c("Nitrogen", "Sulfur"),

  Mean_DOWN = c(
    Mean_MFE_GSE226643_Nitrogen_DOWN,
    Mean_MFE_GSE226643_Sulfur_DOWN
  ),

  Mean_UP = c(
    Mean_MFE_GSE226643_Nitrogen_UP,
    Mean_MFE_GSE226643_Sulfur_UP
  ),

  T_test_P = c(
    TTEST_MFE_GSE226643_Nitrogen_all_UP_vs_DOWN$p.value,
    TTEST_MFE_GSE226643_Sulfur_all_UP_vs_DOWN$p.value
  ),

  KS_all_P = c(
    ks_test_GSE226643_Nitrogen_all_UP_vs_DOWN$p.value,
    ks_test_GSE226643_Sulfur_all_UP_vs_DOWN$p.value
  ),

  KS_075_P = c(
    ks_test_MFE_GSE226643_Nitrogen_075_above$p.value,
    ks_test_MFE_GSE226643_Sulfur_075_above$p.value
  ),

  KS_025_P = c(
    ks_test_MFE_GSE226643_Nitrogen_025_below$p.value,
    ks_test_MFE_GSE226643_Sulfur_025_below$p.value
  ),

  Wilcoxon_P = c(
    wilcox_test_GSE226643_Nitrogen_all_UP_vs_DOWN$p.value,
    wilcox_test_GSE226643_Sulfur_all_UP_vs_DOWN$p.value
  )
)

print(Fig5D_GSE226643_stats)

if (!exists("Fig5_GSE224463_objects", inherits = FALSE)) {
  stop(
    "Run the GSE224463 violin block first to load the saved analysis RDS.",
    call. = FALSE
  )
}
Fig5_GSE224463_fraction_timepoints <- c(
  90L,
  120L,
  180L,
  240L,
  300L
)

Fig5_GSE224463_fraction_manifest <-
  Fig5_GSE224463_objects$contrast_manifest %>%
  dplyr::transmute(
    analysis_id = as.character(.data$analysis_id),
    series = as.character(.data$series),
    target_time_min = as.integer(.data$target_time_min)
  ) %>%
  dplyr::filter(
    .data$series == "long",
    .data$target_time_min %in% Fig5_GSE224463_fraction_timepoints
  ) %>%
  dplyr::arrange(.data$target_time_min) %>%
  dplyr::mutate(
    target_time_label = factor(
      paste0(.data$target_time_min, " min"),
      levels = paste0(
        Fig5_GSE224463_fraction_timepoints,
        " min"
      )
    )
  )

if (nrow(Fig5_GSE224463_fraction_manifest) != 5L) {
  stop(
    "Expected five long-series GSE224463 contrasts.",
    call. = FALSE
  )
}

Fig5_GSE224463_fraction_source <-
  Fig5_GSE224463_objects$significant_DE_abundance_fractions %>%
  dplyr::filter(
    as.character(.data$series) == "long",
    .data$abundance_fraction_eligible %in% TRUE
  ) %>%
  dplyr::transmute(
    analysis_id = as.character(.data$analysis_id),
    Control_low_fraction = as.numeric(.data$Control_low_fraction),
    Treated_low_fraction = as.numeric(.data$Treated_low_fraction),
    Control_high_fraction = as.numeric(.data$Control_high_fraction),
    Treated_high_fraction = as.numeric(.data$Treated_high_fraction)
  ) %>%
  dplyr::inner_join(
    Fig5_GSE224463_fraction_manifest,
    by = "analysis_id"
  ) %>%
  dplyr::arrange(.data$target_time_min)

if (nrow(Fig5_GSE224463_fraction_source) != 5L) {
  stop(
    "One or more long-series contrasts lack eligible significant-DE abundance fractions.",
    call. = FALSE
  )
}

Fig5_GSE224463_fraction_plot_data <- dplyr::bind_rows(
  Fig5_GSE224463_fraction_source %>%
    dplyr::transmute(
      target_time_label = .data$target_time_label,
      DGfold_quartile = "High ΔGfold quartile",
      baseline_0min = .data$Control_high_fraction,
      target_time = .data$Treated_high_fraction
    ),
  Fig5_GSE224463_fraction_source %>%
    dplyr::transmute(
      target_time_label = .data$target_time_label,
      DGfold_quartile = "Low ΔGfold quartile",
      baseline_0min = .data$Control_low_fraction,
      target_time = .data$Treated_low_fraction
    )
) %>%
  tidyr::pivot_longer(
    cols = c(baseline_0min, target_time),
    names_to = "stage",
    values_to = "fraction"
  ) %>%
  dplyr::filter(is.finite(.data$fraction)) %>%
  dplyr::mutate(
    DGfold_quartile = factor(
      .data$DGfold_quartile,
      levels = c(
        "High ΔGfold quartile",
        "Low ΔGfold quartile"
      )
    ),
    stage = factor(
      .data$stage,
      levels = c("baseline_0min", "target_time"),
      labels = c("0 min baseline", "Target time")
    )
  )

Fig5_GSE224463_fraction_y_max <- max(
  Fig5_GSE224463_fraction_plot_data$fraction,
  na.rm = TRUE
) * 1.20

Fig5_Fraction_GSE224463_LongSeries_block <- ggplot(
  Fig5_GSE224463_fraction_plot_data,
  aes(
    x = .data$target_time_label,
    y = .data$fraction,
    fill = .data$stage
  )
) +
  geom_col(
    position = position_dodge(width = 0.78),
    width = 0.68,
    colour = "black",
    linewidth = 0.25
  ) +
  geom_text(
    aes(label = sprintf("%.3f", .data$fraction)),
    position = position_dodge(width = 0.78),
    vjust = -0.35,
    size = 1.5
  ) +
  facet_wrap(
    ~ DGfold_quartile,
    nrow = 1
  ) +
  coord_cartesian(
    ylim = c(0, Fig5_GSE224463_fraction_y_max)
  ) +
  scale_fill_manual(
    values = c(
      "0 min baseline" = "#8cc5e3",
      "Target time" = "#2066a8"
    )
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "top",
    axis.title.x = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      colour = "black",
      size = 8
    ),
    axis.text.y = element_text(colour = "black"),
    strip.text = element_text(
      size = 5,
      face = "bold"
    )
  ) +
  labs(
    y = "Fraction of total significant-DE counts",
    fill = NULL
  )

print(Fig5_Fraction_GSE224463_LongSeries_block)

ggsave(
  filename = file.path(Fig5_output_dir, "Fig5_GSE224463_long_series_significant_DE_DGfold_abundance_fractions.tiff"),
  plot = Fig5_Fraction_GSE224463_LongSeries_block,
  width = 6,
  height = 5,
  units = "in",
  dpi = 600,
  device = "tiff",
  compression = "lzw",
  bg = "white"
)


Fig5_GSE226643_fraction_limitation_order <- c(
  "Nitrogen limitation",
  "Sulfur limitation"
)

Fig5_GSE226643_fraction_source <-
  Fig5_GSE226643_objects$GSE226643_significant_DE_abundance_fraction %>%
  dplyr::transmute(
    limitation = as.character(.data$limitation),
    Control_low_fraction = as.numeric(.data$Control_low_fraction),
    Treated_low_fraction = as.numeric(.data$Treated_low_fraction),
    Control_high_fraction = as.numeric(.data$Control_high_fraction),
    Treated_high_fraction = as.numeric(.data$Treated_high_fraction)
  ) %>%
  dplyr::filter(
    .data$limitation %in% Fig5_GSE226643_fraction_limitation_order
  ) %>%
  dplyr::arrange(
    match(
      .data$limitation,
      Fig5_GSE226643_fraction_limitation_order
    )
  )

if (nrow(Fig5_GSE226643_fraction_source) != 2L) {
  stop(
    "Expected one eligible significant-DE fraction summary for nitrogen and sulfur limitation.",
    call. = FALSE
  )
}

Fig5_GSE226643_fraction_plot_data <- dplyr::bind_rows(
  Fig5_GSE226643_fraction_source %>%
    dplyr::transmute(
      limitation = .data$limitation,
      DGfold_quartile = "High ΔGfold quartile",
      early_3h = .data$Control_high_fraction,
      stationary_12h = .data$Treated_high_fraction
    ),
  Fig5_GSE226643_fraction_source %>%
    dplyr::transmute(
      limitation = .data$limitation,
      DGfold_quartile = "Low ΔGfold quartile",
      early_3h = .data$Control_low_fraction,
      stationary_12h = .data$Treated_low_fraction
    )
) %>%
  tidyr::pivot_longer(
    cols = c(early_3h, stationary_12h),
    names_to = "stage",
    values_to = "fraction"
  ) %>%
  dplyr::filter(is.finite(.data$fraction)) %>%
  dplyr::mutate(
    limitation = factor(
      .data$limitation,
      levels = Fig5_GSE226643_fraction_limitation_order
    ),
    DGfold_quartile = factor(
      .data$DGfold_quartile,
      levels = c(
        "High ΔGfold quartile",
        "Low ΔGfold quartile"
      )
    ),
    stage = factor(
      .data$stage,
      levels = c("early_3h", "stationary_12h"),
      labels = c("3 h early", "12 h stationary")
    )
  )

Fig5_GSE226643_fraction_y_max <- 0.6

Fig5_Fraction_GSE226643_Nitrogen_Sulfur_block <- ggplot(
  Fig5_GSE226643_fraction_plot_data,
  aes(
    x = .data$limitation,
    y = .data$fraction,
    fill = .data$stage
  )
) +
  geom_col(
    position = position_dodge(width = 0.78),
    width = 0.68,
    colour = "black",
    linewidth = 0.25
  ) +
  geom_text(
    aes(label = sprintf("%.3f", .data$fraction)),
    position = position_dodge(width = 0.78),
    vjust = -0.35,
    size = 3
  ) +
  facet_wrap(
    ~ DGfold_quartile,
    nrow = 1
  ) +
  scale_y_continuous(
    breaks = seq(0, 0.6, by = 0.2)
  ) +
  coord_cartesian(
    ylim = c(0, 0.6)
  ) +
  scale_x_discrete(
    labels = c(
      "Nitrogen limitation" = "Nitrogen\nlimitation",
      "Sulfur limitation" = "Sulfur\nlimitation"
    )
  ) +
  scale_fill_manual(
    values = c(
      "3 h early" = "#8cc5e3",
      "12 h stationary" = "#2066a8"
    )
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "top",
    axis.title.x = element_blank(),
    axis.text.x = element_text(
      colour = "black",
      size = 9
    ),
    axis.text.y = element_text(colour = "black"),
    strip.text = element_text(
      size = 9,
      face = "bold"
    )
  ) +
  labs(
    y = "Fraction of total significant-DE counts",
    fill = NULL
  )

print(Fig5_Fraction_GSE226643_Nitrogen_Sulfur_block)

ggsave(
  filename = file.path(Fig5_output_dir, "Fig5_GSE226643_N_S_significant_DE_DGfold_abundance_fractions.tiff"),
  plot = Fig5_Fraction_GSE226643_Nitrogen_Sulfur_block,
  width = 3.5,
  height = 5,
  units = "in",
  dpi = 600,
  device = "tiff",
  compression = "lzw",
  bg = "white"
)


Fig5E_perm_B <- 10000L
Fig5E_perm_seed <- 20260918L

Fig5E_prepare_permutation <- function(
    dat, dataset, comparison, reference_col, target_col, old
) {
  
  fail <- function(msg) stop(dataset, " / ", comparison, ": ", msg)
  
  required <- c("Normalized_MFE", reference_col, target_col)
  
  if (!all(required %in% names(dat)) || nrow(dat) < 4L) {
    fail("Missing transcript-level inputs.")
  }
  
  mfe <- dat$Normalized_MFE
  counts <- cbind(dat[[reference_col]], dat[[target_col]])
  
  if (
    any(!is.finite(mfe)) ||
    any(!is.finite(counts)) ||
    any(counts < 0) ||
    any(!is.finite(colSums(counts))) ||
    any(colSums(counts) <= 0)
  ) {
    fail("Invalid folding values or normalized counts.")
  }
  
  cutoffs <- quantile(mfe, c(0.25, 0.75), names = FALSE)
  
  category <- ifelse(
    mfe <= cutoffs[1], 2L,
    ifelse(mfe >= cutoffs[2], 1L, 0L)
  )
  
  if (!any(category == 1L) || !any(category == 2L)) {
    fail("One of the quartile groups is empty.")
  }
  
  weights <- sweep(counts, 2, colSums(counts), "/")
  
  fractions <- vapply(1:2, function(k) {
    colSums(weights[category == k, , drop = FALSE])
  }, numeric(2))
  
  if (nrow(old) != 1L) {
    fail("Expected one existing fraction summary.")
  }
  
  expected <- c(
    old$Control_high_fraction, old$Treated_high_fraction,
    old$Control_low_fraction, old$Treated_low_fraction
  )
  
  if (
    length(expected) != 4L ||
    any(!is.finite(expected)) ||
    max(abs(as.vector(fractions) - expected)) > 1e-10
  ) {
    fail("Recalculated fractions do not match the existing bars.")
  }
  
  list(
    Dataset = dataset,
    Comparison = comparison,
    category = category,
    fractions = fractions,
    delta = weights[, 2] - weights[, 1]
  )
}

stopifnot(
  nrow(Fig5_GSE224463_fraction_manifest) == 5L,
  !anyDuplicated(Fig5_GSE224463_fraction_manifest$analysis_id)
)

Fig5E_perm_inputs <- lapply(
  seq_len(nrow(Fig5_GSE224463_fraction_manifest)),
  function(i) {
    
    info <- Fig5_GSE224463_fraction_manifest[i, ]
    
    source <-
      Fig5_GSE224463_objects$significant_DE_thermodynamic_transcripts
    
    dat <- source[which(
      as.character(source$analysis_id) == as.character(info$analysis_id) &
        as.character(source$series) == "long"
    ), ]
    
    old <- Fig5_GSE224463_fraction_source[which(
      as.character(Fig5_GSE224463_fraction_source$analysis_id) ==
        as.character(info$analysis_id)
    ), ]
    
    Fig5E_prepare_permutation(
      dat,
      "GSE224463",
      paste0(info$target_time_min, " min vs 0 min"),
      "reference_0min_mean_normalised_counts",
      "target_time_mean_normalised_counts",
      old
    )
  }
)

Fig5E_perm_limitation_objects <- c(
  "Nitrogen limitation" = "GSE226643_N_sig_thermodynamic",
  "Sulfur limitation" = "GSE226643_S_sig_thermodynamic"
)

Fig5E_perm_inputs <- c(
  Fig5E_perm_inputs,
  lapply(names(Fig5E_perm_limitation_objects), function(limitation) {
    
    dat <- Fig5_GSE226643_objects[[
      Fig5E_perm_limitation_objects[[limitation]]
    ]]
    
    old <- Fig5_GSE226643_fraction_source[which(
      as.character(Fig5_GSE226643_fraction_source$limitation) == limitation
    ), ]
    
    Fig5E_prepare_permutation(
      dat,
      "GSE226643",
      paste0(limitation, ": 12 h vs 3 h"),
      "early_3h_mean_normalised_counts",
      "stationary_12h_mean_normalised_counts",
      old
    )
  })
)

stopifnot(length(Fig5E_perm_inputs) == 7L)

message("All 28 Fig. 5E bar fractions match the existing source data.")

set.seed(Fig5E_perm_seed)

Fig5E_perm_runs <- lapply(Fig5E_perm_inputs, function(x) {
  
  message("Permutations: ", x$Dataset, " / ", x$Comparison)
  
  observed <- x$fractions[2, ] - x$fractions[1, ]
  
  null <- replicate(Fig5E_perm_B, {
    shuffled <- x$category[sample.int(length(x$category))]
    c(
      sum(x$delta[shuffled == 1L]),
      sum(x$delta[shuffled == 2L])
    )
  })
  
  b <- rowSums(abs(null) >= abs(observed) - 1e-12)
  
  list(
    results = data.frame(
      Dataset = x$Dataset,
      Comparison = x$Comparison,
      Quartile = c("Upper / less negative", "Lower / more negative"),
      Fraction_reference = x$fractions[1, ],
      Fraction_target = x$fractions[2, ],
      Delta_fraction = observed,
      N_transcripts = length(x$category),
      N_quartile = c(sum(x$category == 1L), sum(x$category == 2L)),
      P_permutation = (b + 1) / (Fig5E_perm_B + 1),
      stringsAsFactors = FALSE
    ),
    null_delta = null
  )
})

Fig5E_quartile_permutation_results <- do.call(
  rbind,
  lapply(Fig5E_perm_runs, function(x) x$results)
)

stopifnot(nrow(Fig5E_quartile_permutation_results) == 14L)

Fig5E_quartile_permutation_results$P_adj_BH <- p.adjust(
  Fig5E_quartile_permutation_results$P_permutation,
  method = "BH"
)

Fig5E_quartile_permutation_results$Figure_label <- ifelse(
  Fig5E_quartile_permutation_results$P_adj_BH >= 0.05,
  "NS",
  sprintf("P = %.4f", Fig5E_quartile_permutation_results$P_adj_BH)
)

print(
  Fig5E_quartile_permutation_results[, c(
    "Dataset", "Comparison", "Quartile", "Fraction_reference",
    "Fraction_target", "Delta_fraction", "P_permutation",
    "P_adj_BH", "Figure_label"
  )],
  row.names = FALSE,
  digits = 7
)


Fig5E_rep_times <- c(90, 120, 180, 240, 300)
Fig5E_rep_ids <- c("LSP_A", "LSP_C", "LSP_E")

Fig5E_rep_meta <- as.data.frame(dplyr::filter(
  Fig5_GSE224463_objects$WT_manifest,
  series == "long",
  time_min %in% c(0, Fig5E_rep_times)
))

stopifnot(
  nrow(Fig5E_rep_meta) == 18L,
  !anyDuplicated(Fig5E_rep_meta$file_name),
  setequal(
    as.character(Fig5E_rep_meta$author_replicate),
    Fig5E_rep_ids
  ),
  all(table(
    Fig5E_rep_meta$author_replicate,
    Fig5E_rep_meta$time_min
  ) == 1L)
)

Fig5E_rep_count_tables <-
  Fig5_GSE224463_objects$DESeq2_normalised_counts_by_analysis

if (!is.list(Fig5E_rep_count_tables) || is.null(names(Fig5E_rep_count_tables))) {
  stop(
    "The GSE224463 RDS lacks the named normalized-count tables required for Figure 5E.",
    call. = FALSE
  )
}

Fig5E_rep_fraction_data <- dplyr::bind_rows(lapply(
  Fig5E_rep_times, function(tm) {
    
    old <- dplyr::filter(
      Fig5_GSE224463_fraction_source,
      target_time_min == tm
    )
    stopifnot(nrow(old) == 1L)
    id <- as.character(old$analysis_id[[1]])
    
    genes <- as.data.frame(dplyr::filter(
      Fig5_GSE224463_objects$significant_DE_thermodynamic_transcripts,
      analysis_id == id,
      series == "long",
      target_time_min == tm
    ))
    
    tab <- Fig5E_rep_count_tables[[id]]

    if (is.null(tab)) {
      stop("Missing normalized-count table for analysis: ", id, call. = FALSE)
    }

    tab <- as.data.frame(tab)
    
    tab$Gene_id <- tolower(as.character(tab$Gene_id))
    genes$Gene_id <- tolower(as.character(genes$Gene_id))
    
    meta <- Fig5E_rep_meta[
      Fig5E_rep_meta$time_min %in% c(0, tm), ,
      drop = FALSE
    ]
    
    gene_rows <- match(genes$Gene_id, tab$Gene_id)
    
    stopifnot(
      nrow(genes) >= 4L,
      !anyNA(genes$Gene_id),
      !anyDuplicated(genes$Gene_id),
      !anyDuplicated(tab$Gene_id),
      !anyNA(gene_rows),
      all(meta$file_name %in% names(tab)),
      all(is.finite(genes$Normalized_MFE))
    )
    
    counts <- as.matrix(
      tab[gene_rows, meta$file_name, drop = FALSE]
    )
    
    stopifnot(all(is.finite(counts)), all(counts >= 0))
    
    total <- colSums(counts)
    stopifnot(all(is.finite(total)), all(total > 0))
    
    q <- quantile(
      genes$Normalized_MFE,
      c(0.25, 0.75),
      names = FALSE
    )
    
    lower <- genes$Normalized_MFE <= q[1]
    upper <- !lower & genes$Normalized_MFE >= q[2]
    stopifnot(any(lower), any(upper))
    
    upper_counts <- colSums(counts[upper, , drop = FALSE])
    lower_counts <- colSums(counts[lower, , drop = FALSE])
    ref <- meta$time_min == 0
    
    reconstructed <- c(
      sum(upper_counts[ref]) / sum(total[ref]),
      sum(lower_counts[ref]) / sum(total[ref]),
      sum(upper_counts[!ref]) / sum(total[!ref]),
      sum(lower_counts[!ref]) / sum(total[!ref])
    )
    
    expected <- c(
      old$Control_high_fraction,
      old$Control_low_fraction,
      old$Treated_high_fraction,
      old$Treated_low_fraction
    )
    
    if (any(!is.finite(expected)) ||
        max(abs(reconstructed - expected)) > 1e-10) {
      stop(
        "The saved counts do not reproduce the bars at ",
        tm, " min."
      )
    }
    
    data.frame(
      Time_min = tm,
      Replicate = rep(as.character(meta$author_replicate), 2),
      State = rep(ifelse(ref, "Baseline", "Target"), 2),
      Quartile = rep(
        c("Upper / less negative", "Lower / more negative"),
        each = nrow(meta)
      ),
      Fraction = c(upper_counts / total, lower_counts / total)
    )
  }
))

message("All 20 original GSE224463 bar fractions reproduced correctly.")

Fig5E_rep_fraction_data$Quartile <- factor(
  Fig5E_rep_fraction_data$Quartile,
  levels = c("Upper / less negative", "Lower / more negative")
)

Fig5E_rep_target <- dplyr::filter(
  Fig5E_rep_fraction_data,
  State == "Target"
)

Fig5E_rep_fraction_table <- Fig5E_rep_target |>
  dplyr::select(Time_min, Quartile, Replicate, Fraction) |>
  tidyr::pivot_wider(
    names_from = Replicate,
    values_from = Fraction
  ) |>
  dplyr::arrange(Quartile, Time_min)

Fig5E_rep_slopes <- Fig5E_rep_target |>
  dplyr::group_by(Quartile, Replicate) |>
  dplyr::summarise(
    Slope_fraction_per_hour = unname(
      coef(lm(Fraction ~ I(Time_min / 60)))[2]
    ),
    .groups = "drop"
  )

print(
  as.data.frame(Fig5E_rep_fraction_table),
  row.names = FALSE,
  digits = 6
)

print(
  as.data.frame(Fig5E_rep_slopes),
  row.names = FALSE,
  digits = 6
)

Fig5E_rep_plot <- ggplot2::ggplot(
  Fig5E_rep_target,
  ggplot2::aes(
    Time_min, Fraction,
    colour = Replicate,
    group = Replicate
  )
) +
  ggplot2::geom_line() +
  ggplot2::geom_point(size = 2) +
  ggplot2::facet_wrap(~ Quartile) +
  ggplot2::scale_x_continuous(breaks = Fig5E_rep_times) +
  ggplot2::theme_classic() +
  ggplot2::labs(
    x = "Time (min)",
    y = "Fraction of counts in the comparison-specific DE pool"
  )

print(Fig5E_rep_plot)


Fig5E_rep_trend_results <- Fig5E_rep_slopes |>
  dplyr::group_by(Quartile) |>
  dplyr::group_modify(~ {
    
    slopes <- .x$Slope_fraction_per_hour
    
    stopifnot(
      length(slopes) == 3L,
      !anyDuplicated(.x$Replicate),
      all(is.finite(slopes))
    )
    
    test <- stats::t.test(
      slopes,
      mu = 0,
      alternative = "two.sided"
    )
    
    data.frame(
      N_experiments = length(slopes),
      Mean_slope_per_hour = mean(slopes),
      CI95_lower = test$conf.int[1],
      CI95_upper = test$conf.int[2],
      P_t_test = test$p.value
    )
  }) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    P_adj_BH = p.adjust(P_t_test, method = "BH")
  )

print(
  as.data.frame(Fig5E_rep_trend_results),
  row.names = FALSE,
  digits = 7
)


if (!exists("Fig5_GSE224463_objects", inherits = FALSE)) {
  stop(
    "Run the GSE224463 violin block first to load the saved analysis RDS.",
    call. = FALSE
  )
}
Fig5_GSE224463_weighted_bar_colour <- "#D55E00"

Fig5_GSE224463_weighted_timepoints <- c(
  90L,
  120L,
  180L,
  240L,
  300L
)

Fig5_GSE224463_weighted_manifest <-
  Fig5_GSE224463_objects$contrast_manifest %>%
  dplyr::transmute(
    analysis_id = as.character(.data$analysis_id),
    series = as.character(.data$series),
    target_time_min = as.integer(.data$target_time_min)
  ) %>%
  dplyr::filter(
    .data$series == "long",
    .data$target_time_min %in% Fig5_GSE224463_weighted_timepoints
  ) %>%
  dplyr::arrange(.data$target_time_min) %>%
  dplyr::mutate(
    figure_order = dplyr::row_number(),
    target_time_label = factor(
      paste0(.data$target_time_min, " min"),
      levels = paste0(
        Fig5_GSE224463_weighted_timepoints,
        " min"
      )
    )
  )

Fig5_GSE224463_weighted_sample_indices <-
  Fig5_GSE224463_objects$weighted_sample_indices %>%
  dplyr::filter(as.character(.data$series) == "long") %>%
  dplyr::transmute(
    series = as.character(.data$series),
    time_min = as.integer(.data$time_min),
    author_replicate = as.character(.data$author_replicate),
    abundance_weighted_Normalized_MFE =
      as.numeric(.data$abundance_weighted_Normalized_MFE)
  )

Fig5_GSE224463_weighted_reference_indices <-
  Fig5_GSE224463_weighted_sample_indices %>%
  dplyr::filter(.data$time_min == 0L) %>%
  dplyr::transmute(
    series = .data$series,
    author_replicate = .data$author_replicate,
    reference_0min_abundance_weighted_Normalized_MFE =
      .data$abundance_weighted_Normalized_MFE
  )

Fig5_GSE224463_weighted_target_indices <-
  Fig5_GSE224463_weighted_sample_indices %>%
  dplyr::filter(
    .data$time_min %in% Fig5_GSE224463_weighted_timepoints
  ) %>%
  dplyr::transmute(
    series = .data$series,
    author_replicate = .data$author_replicate,
    target_time_min = .data$time_min,
    target_time_abundance_weighted_Normalized_MFE =
      .data$abundance_weighted_Normalized_MFE
  )

Fig5_GSE224463_weighted_paired_deltas <-
  Fig5_GSE224463_weighted_target_indices %>%
  dplyr::inner_join(
    Fig5_GSE224463_weighted_reference_indices,
    by = c("series", "author_replicate")
  ) %>%
  dplyr::inner_join(
    Fig5_GSE224463_weighted_manifest,
    by = c("series", "target_time_min")
  ) %>%
  dplyr::mutate(
    paired_delta_target_time_minus_0min =
      .data$target_time_abundance_weighted_Normalized_MFE -
      .data$reference_0min_abundance_weighted_Normalized_MFE
  )

Fig5_GSE224463_weighted_MFE_plot_data <-
  Fig5_GSE224463_weighted_paired_deltas %>%
  dplyr::group_by(
    .data$analysis_id,
    .data$target_time_min,
    .data$target_time_label,
    .data$figure_order
  ) %>%
  dplyr::summarise(
    N_paired_biological_replicates = dplyr::n(),
    mean_paired_delta_target_time_minus_0min = mean(
      .data$paired_delta_target_time_minus_0min
    ),
    sd_paired_delta_target_time_minus_0min = stats::sd(
      .data$paired_delta_target_time_minus_0min
    ),
    sem_paired_delta_target_time_minus_0min =
      .data$sd_paired_delta_target_time_minus_0min /
      sqrt(.data$N_paired_biological_replicates),
    .groups = "drop"
  ) %>%
  dplyr::left_join(
    Fig5_GSE224463_objects$full_transcriptome_weighted_delta_summary %>%
      dplyr::filter(
        as.character(.data$series) == "long",
        .data$descriptor == "Normalized_MFE"
      ) %>%
      dplyr::transmute(
        analysis_id = as.character(.data$analysis_id),
        master_delta_target_time_minus_0min = as.numeric(
          .data$delta_target_time_minus_0min
        )
      ),
    by = "analysis_id"
  ) %>%
  dplyr::arrange(.data$figure_order)

if (nrow(Fig5_GSE224463_weighted_MFE_plot_data) != 5L ||
    any(
      Fig5_GSE224463_weighted_MFE_plot_data$
      N_paired_biological_replicates != 3L
    )) {
  stop(
    "Expected three paired biological replicates at each of five long-series time points.",
    call. = FALSE
  )
}

if (max(
  abs(
    Fig5_GSE224463_weighted_MFE_plot_data$
    mean_paired_delta_target_time_minus_0min -
    Fig5_GSE224463_weighted_MFE_plot_data$
    master_delta_target_time_minus_0min
  ),
  na.rm = TRUE
) > 1e-12) {
  stop(
    "The GSE224463 weighted ΔGfold bars do not match the saved analysis values.",
    call. = FALSE
  )
}

Fig5_GSE224463_weighted_MFE_y_limit <- 0.05

Fig5_GSE224463_weighted_MFE_label_offset <-
  Fig5_GSE224463_weighted_MFE_y_limit * 0.06

Fig5_GSE224463_weighted_MFE_plot_data <-
  Fig5_GSE224463_weighted_MFE_plot_data %>%
  dplyr::mutate(
    errorbar_ymin = dplyr::if_else(
      .data$mean_paired_delta_target_time_minus_0min >= 0,
      .data$mean_paired_delta_target_time_minus_0min,
      .data$mean_paired_delta_target_time_minus_0min -
        .data$sem_paired_delta_target_time_minus_0min
    ),
    errorbar_ymax = dplyr::if_else(
      .data$mean_paired_delta_target_time_minus_0min >= 0,
      .data$mean_paired_delta_target_time_minus_0min +
        .data$sem_paired_delta_target_time_minus_0min,
      .data$mean_paired_delta_target_time_minus_0min
    ),
    label_y = dplyr::if_else(
      .data$mean_paired_delta_target_time_minus_0min >= 0,
      .data$mean_paired_delta_target_time_minus_0min +
        .data$sem_paired_delta_target_time_minus_0min +
        Fig5_GSE224463_weighted_MFE_label_offset,
      .data$mean_paired_delta_target_time_minus_0min -
        .data$sem_paired_delta_target_time_minus_0min -
        Fig5_GSE224463_weighted_MFE_label_offset
    ),
    label_vjust = dplyr::if_else(
      .data$mean_paired_delta_target_time_minus_0min >= 0,
      0,
      1
    )
  )

Fig5_Weighted_GSE224463_LongSeries_DGfold <- ggplot(
  Fig5_GSE224463_weighted_MFE_plot_data,
  aes(
    x = .data$target_time_label,
    y = .data$mean_paired_delta_target_time_minus_0min
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.4,
    colour = "black"
  ) +
  geom_col(
    width = 0.75,
    fill = Fig5_GSE224463_weighted_bar_colour,
    colour = "black",
    linewidth = 0.3
  ) +
  geom_errorbar(
    aes(
      ymin = .data$errorbar_ymin,
      ymax = .data$errorbar_ymax
    ),
    width = 0.16,
    linewidth = 0.4
  ) +
  geom_text(
    aes(
      y = .data$label_y,
      label = sprintf(
        "%+.4f",
        .data$mean_paired_delta_target_time_minus_0min
      ),
      vjust = .data$label_vjust
    ),
    size = 1
  ) +
  scale_y_continuous(
    breaks = seq(-0.05, 0.05, by = 0.025)
  ) +
  coord_cartesian(
    ylim = c(-0.05, 0.05)
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      colour = "black",
      size = 8
    )
  ) +
  labs(
    y = expression(
      "Paired abundance-weighted " * Delta * G[fold] * " / nt"
    )
  )

print(Fig5_Weighted_GSE224463_LongSeries_DGfold)

ggsave(
  filename = file.path(
    Fig5_output_dir,
    "Fig5_GSE224463_long_series_full_transcriptome_weighted_DGfold.tiff"
  ),
  plot = Fig5_Weighted_GSE224463_LongSeries_DGfold,
  width = 2.5,
  height = 5,
  units = "in",
  dpi = 600,
  device = "tiff",
  compression = "lzw",
  bg = "white"
)

Fig5F_GSE224463_gene_contributions <-
  Fig5_GSE224463_objects$full_transcriptome_weighted_gene_contributions

Fig5F_GSE224463_weighted_stats <- function(df, seed) {

  df <- df %>%
    dplyr::filter(
      is.finite(Normalized_MFE),
      is.finite(reference_0min_mean_relative_RNA_abundance),
      is.finite(target_time_mean_relative_RNA_abundance),
      is.finite(contrast_weight_target_time_minus_0min),
      is.finite(abundance_weighted_MFE_contribution)
    )

  x <- df$Normalized_MFE

  w0 <- df$reference_0min_mean_relative_RNA_abundance
  wt <- df$target_time_mean_relative_RNA_abundance

  contrast_weight <- df$contrast_weight_target_time_minus_0min
  contribution <- df$abundance_weighted_MFE_contribution

  baseline_weighted_MFE <- sum(w0 * x)
  target_weighted_MFE <- sum(wt * x)
  delta <- sum(contribution)

  stopifnot(
    abs(
      delta -
        (target_weighted_MFE - baseline_weighted_MFE)
    ) < 1e-12
  )

  set.seed(seed)

  boot <- replicate(
    5000L,
    sum(
      sample(
        contribution,
        size = length(contribution),
        replace = TRUE
      )
    )
  )

  perm <- replicate(
    5000L,
    sum(
      contrast_weight *
        sample(
          x,
          size = length(x),
          replace = FALSE
        )
    )
  )

  tibble::tibble(
    baseline_weighted_MFE = baseline_weighted_MFE,
    target_weighted_MFE = target_weighted_MFE,
    delta_target_minus_0min = delta,
    percent_change =
      100 * delta / abs(baseline_weighted_MFE),
    bootstrap_CI_low =
      unname(stats::quantile(boot, 0.025)),
    bootstrap_CI_high =
      unname(stats::quantile(boot, 0.975)),
    permutation_p_value =
      (
        sum(abs(perm) >= abs(delta)) + 1L
      ) / 5001L,
    bootstrap_p_sign =
      2 * min(
        mean(boot <= 0),
        mean(boot >= 0)
      ),
    N_transcripts = length(x)
  )
}

Fig5F_GSE224463_contrast_order <-
  Fig5_GSE224463_objects$contrast_manifest %>%
  dplyr::mutate(
    contrast_index = dplyr::row_number()
  ) %>%
  dplyr::filter(
    as.character(series) == "long",
    target_time_min %in% c(90, 120, 180, 240, 300)
  ) %>%
  dplyr::select(
    analysis_id,
    target_time_min,
    contrast_index
  )

Fig5F_GSE224463_stats <- dplyr::bind_rows(

  lapply(
    seq_len(nrow(Fig5F_GSE224463_contrast_order)),
    function(i) {

      current_analysis_id <-
        Fig5F_GSE224463_contrast_order$analysis_id[i]

      current_time <-
        Fig5F_GSE224463_contrast_order$target_time_min[i]

      current_contrast_index <-
        Fig5F_GSE224463_contrast_order$contrast_index[i]

      current_data <-
        Fig5F_GSE224463_gene_contributions %>%
        dplyr::filter(
          analysis_id == current_analysis_id,
          series == "long",
          target_time_min == current_time
        )

      Fig5F_GSE224463_weighted_stats(
        current_data,
        seed = 20260807L + current_contrast_index
      ) %>%
        dplyr::mutate(
          analysis_id = current_analysis_id,
          target_time_min = current_time,
          .before = 1
        )
    }
  )

) %>%
  dplyr::mutate(
    permutation_p_adj_BH =
      p.adjust(
        permutation_p_value,
        method = "BH"
      ),
    bootstrap_p_sign_adj_BH =
      p.adjust(
        bootstrap_p_sign,
        method = "BH"
      )
  ) %>%
  dplyr::arrange(target_time_min)

Fig5F_GSE224463_master_check <-
  Fig5_GSE224463_objects$full_transcriptome_weighted_delta_summary %>%
  dplyr::filter(
    as.character(series) == "long",
    descriptor == "Normalized_MFE",
    target_time_min %in% c(90, 120, 180, 240, 300)
  ) %>%
  dplyr::select(
    analysis_id,
    master_delta = delta_target_time_minus_0min
  )

Fig5F_GSE224463_stats <-
  Fig5F_GSE224463_stats %>%
  dplyr::left_join(
    Fig5F_GSE224463_master_check,
    by = "analysis_id"
  )

stopifnot(
  max(
    abs(
      Fig5F_GSE224463_stats$delta_target_minus_0min -
        Fig5F_GSE224463_stats$master_delta
    )
  ) < 1e-12
)

print(Fig5F_GSE224463_stats)

Fig5_GSE226643_weighted_bar_colour <- "#D55E00"

Fig5_GSE226643_weighted_limitation_order <- c(
  "Nitrogen limitation",
  "Sulfur limitation"
)

Fig5_GSE226643_MFE_sample_indices <-
  Fig5_GSE226643_objects$GSE226643_MFE_sample_indices_by_comparison

Fig5_GSE226643_MFE_master_delta_statistics <-
  Fig5_GSE226643_objects$GSE226643_MFE_delta_statistics

Fig5_GSE226643_weighted_MFE_plot_data <-
  Fig5_GSE226643_MFE_sample_indices %>%
  dplyr::group_by(.data$analysis_id, .data$limitation) %>%
  dplyr::summarise(
    N_early_replicates = sum(
      as.character(.data$arm) == "3 h early exponential"
    ),
    N_stationary_replicates = sum(
      as.character(.data$arm) == "12 h stationary"
    ),
    early_mean = mean(
      .data$abundance_weighted_Normalized_MFE[
        as.character(.data$arm) == "3 h early exponential"
      ]
    ),
    stationary_mean = mean(
      .data$abundance_weighted_Normalized_MFE[
        as.character(.data$arm) == "12 h stationary"
      ]
    ),
    early_SD = stats::sd(
      .data$abundance_weighted_Normalized_MFE[
        as.character(.data$arm) == "3 h early exponential"
      ]
    ),
    stationary_SD = stats::sd(
      .data$abundance_weighted_Normalized_MFE[
        as.character(.data$arm) == "12 h stationary"
      ]
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    delta_stationary_minus_early =
      .data$stationary_mean - .data$early_mean,
    early_SEM =
      .data$early_SD / sqrt(.data$N_early_replicates),
    stationary_SEM =
      .data$stationary_SD / sqrt(.data$N_stationary_replicates),
    delta_SEM = sqrt(
      .data$early_SEM^2 +
        .data$stationary_SEM^2
    )
  ) %>%
  dplyr::left_join(
    Fig5_GSE226643_MFE_master_delta_statistics %>%
      dplyr::transmute(
        analysis_id = as.character(.data$analysis_id),
        master_delta =
          as.numeric(.data$delta_abundance_weighted_Normalized_MFE)
      ),
    by = "analysis_id"
  ) %>%
  dplyr::mutate(
    limitation = factor(
      .data$limitation,
      levels = Fig5_GSE226643_weighted_limitation_order
    )
  ) %>%
  dplyr::arrange(.data$limitation)

if (nrow(Fig5_GSE226643_weighted_MFE_plot_data) != 2L ||
    any(
      Fig5_GSE226643_weighted_MFE_plot_data$N_early_replicates != 2L |
      Fig5_GSE226643_weighted_MFE_plot_data$N_stationary_replicates != 2L
    )) {
  stop(
    "Expected two biological replicates per arm for both GSE226643 contrasts.",
    call. = FALSE
  )
}

if (max(
  abs(
    Fig5_GSE226643_weighted_MFE_plot_data$
    delta_stationary_minus_early -
    Fig5_GSE226643_weighted_MFE_plot_data$master_delta
  ),
  na.rm = TRUE
) > 1e-12) {
  stop(
    "The GSE226643 weighted ΔGfold bars do not match the saved analysis values.",
    call. = FALSE
  )
}

Fig5_GSE226643_weighted_MFE_y_limit <- 0.05

Fig5_GSE226643_weighted_MFE_label_offset <-
  Fig5_GSE226643_weighted_MFE_y_limit * 0.06

Fig5_GSE226643_weighted_MFE_plot_data <-
  Fig5_GSE226643_weighted_MFE_plot_data %>%
  dplyr::mutate(
    errorbar_ymin = dplyr::if_else(
      .data$delta_stationary_minus_early >= 0,
      .data$delta_stationary_minus_early,
      .data$delta_stationary_minus_early - .data$delta_SEM
    ),
    errorbar_ymax = dplyr::if_else(
      .data$delta_stationary_minus_early >= 0,
      .data$delta_stationary_minus_early + .data$delta_SEM,
      .data$delta_stationary_minus_early
    ),
    label_y = dplyr::if_else(
      .data$delta_stationary_minus_early >= 0,
      .data$delta_stationary_minus_early +
        .data$delta_SEM +
        Fig5_GSE226643_weighted_MFE_label_offset,
      .data$delta_stationary_minus_early -
        .data$delta_SEM -
        Fig5_GSE226643_weighted_MFE_label_offset
    ),
    label_vjust = dplyr::if_else(
      .data$delta_stationary_minus_early >= 0,
      0,
      1
    )
  )

Fig5_Weighted_GSE226643_Nitrogen_Sulfur_DGfold <- ggplot(
  Fig5_GSE226643_weighted_MFE_plot_data,
  aes(
    x = .data$limitation,
    y = .data$delta_stationary_minus_early
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.4,
    colour = "black"
  ) +
  geom_col(
    width = 0.65,
    fill = Fig5_GSE226643_weighted_bar_colour,
    colour = "black",
    linewidth = 0.3
  ) +
  geom_errorbar(
    aes(
      ymin = .data$errorbar_ymin,
      ymax = .data$errorbar_ymax
    ),
    width = 0.16,
    linewidth = 0.4
  ) +
  geom_text(
    aes(
      y = .data$label_y,
      label = sprintf(
        "%+.4f",
        .data$delta_stationary_minus_early
      ),
      vjust = .data$label_vjust
    ),
    size = 3
  ) +
  scale_x_discrete(
    labels = c(
      "Nitrogen limitation" = "Nitrogen\nlimitation",
      "Sulfur limitation" = "Sulfur\nlimitation"
    )
  ) +
  scale_y_continuous(
    breaks = seq(-0.05, 0.05, by = 0.025)
  ) +
  coord_cartesian(
    ylim = c(
      -Fig5_GSE226643_weighted_MFE_y_limit,
      Fig5_GSE226643_weighted_MFE_y_limit
    )
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(colour = "black")
  ) +
  labs(
    y = expression(
      Delta * " normalized ΔGfold (ΔGfold / nt)"
    )
  )

print(Fig5_Weighted_GSE226643_Nitrogen_Sulfur_DGfold)

ggsave(
  filename = file.path(
    Fig5_output_dir,
    "Fig5_GSE226643_full_transcriptome_weighted_DGfold.tiff"
  ),
  plot = Fig5_Weighted_GSE226643_Nitrogen_Sulfur_DGfold,
  width = 1.8,
  height = 5,
  units = "in",
  dpi = 600,
  device = "tiff",
  compression = "lzw",
  bg = "white"
)

Fig5F_GSE226643_gene_data <-
  Fig5_GSE226643_objects$GSE226643_MFE_comparison_gene_data

set.seed(20260807L)

Fig5F_GSE226643_bootstrap_sign <-
  Fig5F_GSE226643_gene_data %>%
  dplyr::group_by(analysis_id, limitation) %>%
  dplyr::group_modify(~ {

    boot <- replicate(
      5000L,
      sum(
        sample(
          .x$abundance_weighted_MFE_contribution,
          replace = TRUE
        )
      )
    )

    tibble::tibble(
      bootstrap_sign_P =
        2 * min(
          mean(boot <= 0),
          mean(boot >= 0)
        )
    )
  }) %>%
  dplyr::ungroup()

Fig5F_GSE226643_stats <-
  Fig5_GSE226643_objects$GSE226643_MFE_delta_statistics %>%
  dplyr::left_join(
    Fig5F_GSE226643_bootstrap_sign,
    by = c("analysis_id", "limitation")
  ) %>%
  dplyr::mutate(
    permutation_P_BH =
      p.adjust(descriptor_label_permutation_p, method = "BH"),
    bootstrap_sign_P_BH =
      p.adjust(bootstrap_sign_P, method = "BH")
  ) %>%
  dplyr::transmute(
    limitation,
    early_weighted_MFE =
      early_3h_abundance_weighted_Normalized_MFE,
    stationary_weighted_MFE =
      stationary_12h_abundance_weighted_Normalized_MFE,
    delta_stationary_minus_early =
      delta_abundance_weighted_Normalized_MFE,
    bootstrap_CI_low = gene_bootstrap_CI_low,
    bootstrap_CI_high = gene_bootstrap_CI_high,
    permutation_P = descriptor_label_permutation_p,
    permutation_P_BH,
    bootstrap_sign_P,
    bootstrap_sign_P_BH,
    N_transcripts = N_full_transcriptome_genes
  )

print(Fig5F_GSE226643_stats)
Fig5_GSE131992_EXT_required_objects <- c(
  "sample_manifest",
  "extended_sample_manifest",
  "contrast_manifest",
  "strong_change_features",
  "strong_change_MFE_statistics",
  "fixed_MFE_quartiles",
  "fixed_MFE_quartile_sample_fractions",
  "relative_TPM_weight_matrices",
  "weighted_sample_indices",
  "weighted_contrast_summary",
  "GSEA_bubble_data",
  "exact_MFE",
  "TPM_mass_coverage"
)

Fig5_GSE131992_EXT_missing_objects <- setdiff(
  Fig5_GSE131992_EXT_required_objects,
  names(Fig5_GSE131992_EXT_objects)
)

if (length(Fig5_GSE131992_EXT_missing_objects) > 0L) {
  stop(
    "The GSE131992 analysis RDS lacks required object(s): ",
    paste(Fig5_GSE131992_EXT_missing_objects, collapse = ", "),
    call. = FALSE
  )
}

Fig5_GSE131992_EXT_required_manifest_columns <- c(
  "sample_id",
  "geo_timepoint",
  "biological_replicate"
)

Fig5_GSE131992_EXT_required_weighted_columns <- c(
  "reference_n_biological_replicates",
  "target_n_biological_replicates",
  "reference_SEM",
  "target_SEM",
  "propagated_SEM"
)

if (
  length(
    setdiff(
      Fig5_GSE131992_EXT_required_manifest_columns,
      colnames(Fig5_GSE131992_EXT_objects$extended_sample_manifest)
    )
  ) > 0L ||
  length(
    setdiff(
      Fig5_GSE131992_EXT_required_weighted_columns,
      colnames(Fig5_GSE131992_EXT_objects$weighted_contrast_summary)
    )
  ) > 0L
) {
  stop(
    paste0(
      "The GSE131992 RDS predates the biological-replicate analysis. ",
      "Rerun scripts 00, 01 and 03 before generating Fig. 5."
    ),
    call. = FALSE
  )
}

Fig5_GSE131992_EXT_expected_contrasts <- c(
  "Established_starvation_t19_minus_pre_starvation",
  "Early_recovery_t23_minus_established_starvation_t19",
  "Intermediate_recovery_t26_minus_established_starvation_t19",
  "Late_recovery_t29_minus_established_starvation_t19"
)

Fig5_GSE131992_EXT_selected_timepoints <- c(
  "t01",
  "t19",
  "t23",
  "t26",
  "t29"
)

Fig5_GSE131992_EXT_contrast_labels <- c(
  "Established_starvation_t19_minus_pre_starvation" =
    "Established starvation\nvs pre-starvation",
  "Early_recovery_t23_minus_established_starvation_t19" =
    "Recovery 25 min\nvs starvation",
  "Intermediate_recovery_t26_minus_established_starvation_t19" =
    "Recovery 65 min\nvs starvation",
  "Late_recovery_t29_minus_established_starvation_t19" =
    "Recovery 155 min\nvs starvation"
)

Fig5_GSE131992_EXT_selected_manifest <-
  Fig5_GSE131992_EXT_objects$extended_sample_manifest %>%
  dplyr::arrange(
    .data$timepoint_number,
    .data$biological_replicate
  )

Fig5_GSE131992_EXT_selected_sample_counts <- table(
  factor(
    Fig5_GSE131992_EXT_selected_manifest$geo_timepoint,
    levels = Fig5_GSE131992_EXT_selected_timepoints
  )
)

if (
  nrow(Fig5_GSE131992_EXT_selected_manifest) != 12L ||
  !identical(
    as.integer(Fig5_GSE131992_EXT_selected_sample_counts),
    c(4L, 2L, 2L, 2L, 2L)
  ) ||
  !setequal(
    Fig5_GSE131992_EXT_objects$contrast_manifest$contrast_id,
    Fig5_GSE131992_EXT_expected_contrasts
  )
) {
  stop(
    paste0(
      "The GSE131992 RDS does not contain the required ",
      "t01/t19/t23/t26/t29 analysis with 4/2/2/2/2 biological replicates and four ",
      "anchor-based contrasts."
    ),
    call. = FALSE
  )
}

Fig5_GSE131992_EXT_contrast_colours <- stats::setNames(
  c("#D55E00", "#0072B2", "#0072B2", "#0072B2"),
  unname(Fig5_GSE131992_EXT_contrast_labels)
)


if (!exists("Fig5_curated_GSEA_terms", inherits = FALSE)) {
  stop(
    paste0(
      "Fig5_curated_GSEA_terms is absent. Run the existing Fig. 5A block ",
      "before the GSE131992 addition."
    ),
    call. = FALSE
  )
}

Fig5G_GSE131992_EXT_missing_GSEA_terms <- setdiff(
  Fig5_curated_GSEA_terms,
  unique(
    as.character(
      Fig5_GSE131992_EXT_objects$GSEA_bubble_data$term
    )
  )
)

if (length(Fig5G_GSE131992_EXT_missing_GSEA_terms) > 0L) {
  stop(
    "The GSE131992 GSEA object lacks curated Fig. 5 term(s): ",
    paste(Fig5G_GSE131992_EXT_missing_GSEA_terms, collapse = ", "),
    call. = FALSE
  )
}

Fig5G_GSE131992_EXT_GSEA_plot_data <-
  Fig5_GSE131992_EXT_objects$GSEA_bubble_data %>%
  dplyr::transmute(
    contrast_id = as.character(.data$contrast_id),
    term = as.character(.data$term),
    NES = as.numeric(.data$NES),
    p_adjust = as.numeric(.data$p_adjust)
  ) %>%
  dplyr::filter(
    .data$contrast_id %in% Fig5_GSE131992_EXT_expected_contrasts,
    .data$term %in% Fig5_curated_GSEA_terms
  ) %>%
  dplyr::mutate(
    comparison = factor(
      .data$contrast_id,
      levels = Fig5_GSE131992_EXT_expected_contrasts,
      labels = unname(Fig5_GSE131992_EXT_contrast_labels)
    ),
    term = factor(
      .data$term,
      levels = rev(Fig5_curated_GSEA_terms)
    ),
    significant =
      is.finite(.data$NES) &
      is.finite(.data$p_adjust) &
      .data$p_adjust < 0.05,
    minus_log10_BH = dplyr::if_else(
      .data$significant,
      -log10(
        pmax(
          .data$p_adjust,
          .Machine$double.xmin
        )
      ),
      NA_real_
    )
  )

if (
  nrow(Fig5G_GSE131992_EXT_GSEA_plot_data) !=
  4L * length(Fig5_curated_GSEA_terms) ||
  anyNA(Fig5G_GSE131992_EXT_GSEA_plot_data$comparison) ||
  anyNA(Fig5G_GSE131992_EXT_GSEA_plot_data$term)
) {
  stop(
    "The GSE131992 matched GSEA table is incomplete.",
    call. = FALSE
  )
}

Fig5G_GSE131992_EXT_GSEA <- ggplot(
  Fig5G_GSE131992_EXT_GSEA_plot_data,
  aes(
    x = .data$comparison,
    y = .data$term
  )
) +
  geom_blank() +
  geom_point(
    data = Fig5G_GSE131992_EXT_GSEA_plot_data %>%
      dplyr::filter(.data$significant),
    aes(
      size = .data$minus_log10_BH,
      colour = .data$NES
    ),
    alpha = 0.95
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  scale_colour_gradient2(
    low = "#2166ac",
    mid = "white",
    high = "#b2182b",
    midpoint = 0,
    name = "NES"
  ) +
  scale_size_continuous(
    range = c(0.5, 6),
    name = expression(-log[10]~"BH-adjusted p")
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(
      colour = "black",
      size = 3
    ),
    axis.text.y = element_text(
      colour = "black",
      size = 3
    ),
    plot.caption = element_text(
      hjust = 0,
      size = 3
    )
  ) +
  labs(
    caption = paste0(
      "GSEA of DESeq2 Wald-statistic rankings calculated from independent ",
      "biological replicates"
    )
  )

print(Fig5G_GSE131992_EXT_GSEA)

ggsave(
  filename = file.path(
    Fig5_output_dir,
    "Fig5G_GSE131992_EXT_extended_recovery_GSEA.tiff"
  ),
  plot = Fig5G_GSE131992_EXT_GSEA,
  width = 4,
  height = 5,
  units = "in",
  dpi = 600,
  device = "tiff",
  compression = "lzw",
  bg = "white"
)


Fig5H_GSE131992_EXT_violin_statistics <-
  Fig5_GSE131992_EXT_objects$strong_change_MFE_statistics %>%
  dplyr::filter(
    .data$universe_id == "exclude_rRNA",
    .data$contrast_id %in% Fig5_GSE131992_EXT_expected_contrasts
  ) %>%
  dplyr::mutate(
    contrast_id = factor(
      as.character(.data$contrast_id),
      levels = Fig5_GSE131992_EXT_expected_contrasts
    )
  ) %>%
  dplyr::arrange(.data$contrast_id)

if (
  nrow(Fig5H_GSE131992_EXT_violin_statistics) != 4L ||
  any(
    !is.finite(
      Fig5H_GSE131992_EXT_violin_statistics$
      higher_minus_lower_mean_MFE
    )
  ) ||
  any(
    !is.finite(
      Fig5H_GSE131992_EXT_violin_statistics$
      exploratory_Wilcoxon_BH
    )
  )
) {
  stop(
    "The GSE131992 strong-change MFE statistics are incomplete or inconsistent.",
    call. = FALSE
  )
}

Fig5H_GSE131992_EXT_violin_data <-
  Fig5_GSE131992_EXT_objects$strong_change_features %>%
  dplyr::filter(
    !.data$is_rRNA,
    .data$contrast_id %in% Fig5_GSE131992_EXT_expected_contrasts,
    is.finite(.data$Normalized_MFE),
    !is.na(.data$expression_direction)
  ) %>%
  dplyr::mutate(
    contrast_id = factor(
      as.character(.data$contrast_id),
      levels = Fig5_GSE131992_EXT_expected_contrasts
    ),
    expression_direction = factor(
      as.character(.data$expression_direction),
      levels = c(
        "Lower in target state",
        "Higher in target state"
      ),
      labels = c(
        "Lower\nin target state",
        "Higher\nin target state"
      )
    )
  )

Fig5H_GSE131992_EXT_violin_counts <-
  Fig5H_GSE131992_EXT_violin_data %>%
  dplyr::count(
    .data$contrast_id,
    .data$expression_direction,
    name = "N_genes"
  )

if (
  nrow(Fig5H_GSE131992_EXT_violin_counts) != 8L ||
  any(Fig5H_GSE131992_EXT_violin_counts$N_genes < 2L)
) {
  stop(
    paste0(
      "The GSE131992 violin data do not contain both expression ",
      "directions for all four extended contrasts."
    ),
    call. = FALSE
  )
}

Fig5H_GSE131992_EXT_MFE_y_limits <- c(-0.5, 0)

Fig5H_GSE131992_EXT_facet_table <-
  Fig5H_GSE131992_EXT_violin_statistics %>%
  dplyr::mutate(
    contrast_id = as.character(.data$contrast_id),
    comparison_label = unname(
      Fig5_GSE131992_EXT_contrast_labels[.data$contrast_id]
    ),
    facet_label = paste0(
      .data$comparison_label,
      "\nΔNormalized MFE = ",
      sprintf("%+.4f", .data$higher_minus_lower_mean_MFE),
      " | MFE Wilcoxon BH P = ",
      formatC(
        .data$exploratory_Wilcoxon_BH,
        format = "e",
        digits = 2
      )
    )
  )

Fig5H_GSE131992_EXT_facet_labels <- stats::setNames(
  Fig5H_GSE131992_EXT_facet_table$facet_label,
  Fig5H_GSE131992_EXT_facet_table$contrast_id
)

Fig5H_GSE131992_EXT_violin_output_stubs <- stats::setNames(
  c(
    "established_starvation",
    "recovery_25min",
    "recovery_65min",
    "recovery_155min"
  ),
  Fig5_GSE131992_EXT_expected_contrasts
)

Fig5H_GSE131992_EXT_violin_plots <- list()

for (
  Fig5H_GSE131992_EXT_contrast_id in
  Fig5_GSE131992_EXT_expected_contrasts
) {
  Fig5H_GSE131992_EXT_violin_fill <- if (
    Fig5H_GSE131992_EXT_contrast_id ==
    "Established_starvation_t19_minus_pre_starvation"
  ) {
    "lightgrey"
  } else {
    "grey45"
  }
  
  Fig5H_GSE131992_EXT_current_data <-
    Fig5H_GSE131992_EXT_violin_data %>%
    dplyr::filter(
      as.character(.data$contrast_id) ==
        Fig5H_GSE131992_EXT_contrast_id
    )
  
  Fig5H_GSE131992_EXT_current_plot <- ggplot(
    Fig5H_GSE131992_EXT_current_data,
    aes(
      x = .data$expression_direction,
      y = .data$Normalized_MFE
    )
  ) +
    geom_violin(
      trim = FALSE,
      fill = Fig5H_GSE131992_EXT_violin_fill,
      alpha = 1,
      colour = "black",
      linewidth = 0.25
    ) +
    geom_boxplot(
      width = 0.14,
      outlier.shape = NA,
      fill = "white",
      linewidth = 0.3
    ) +
    coord_cartesian(
      ylim = Fig5H_GSE131992_EXT_MFE_y_limits,
      expand = FALSE
    ) +
    theme_classic(base_size = 11) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_text(
        colour = "black",
        size = 8
      ),
      axis.text.y = element_text(colour = "black"),
      plot.title = element_text(
        hjust = 0.5,
        size = 10
      )
    ) +
    labs(
      title = unname(
        Fig5H_GSE131992_EXT_facet_labels[
          Fig5H_GSE131992_EXT_contrast_id
        ]
      ),
      y = "Normalized ΔGfold (ΔGfold / nt)"
    )
  
  Fig5H_GSE131992_EXT_violin_plots[[
    Fig5H_GSE131992_EXT_contrast_id
  ]] <- Fig5H_GSE131992_EXT_current_plot
  
  print(Fig5H_GSE131992_EXT_current_plot)
  
  ggsave(
    filename = file.path(
      Fig5_output_dir,
      paste0(
        "Fig5H_GSE131992_EXT_DGfold_violin_",
        Fig5H_GSE131992_EXT_violin_output_stubs[
          Fig5H_GSE131992_EXT_contrast_id
        ],
        ".tiff"
      )
    ),
    plot = Fig5H_GSE131992_EXT_current_plot,
    width = 3.2,
    height = 5,
    units = "in",
    dpi = 600,
    device = "tiff",
    compression = "lzw",
    bg = "white"
  )
}

print(
  as.data.frame(
    Fig5H_GSE131992_EXT_violin_statistics
  ),
  row.names = FALSE,
  digits = 10
)


Fig5I_GSE131992_EXT_quartile_contrast_states <-
  Fig5_GSE131992_EXT_objects$contrast_manifest %>%
  dplyr::filter(
    .data$contrast_id %in% Fig5_GSE131992_EXT_expected_contrasts
  ) %>%
  dplyr::transmute(
    contrast_id = as.character(.data$contrast_id),
    reference_timepoint = as.character(.data$reference_timepoint),
    target_timepoint = as.character(.data$target_timepoint)
  ) %>%
  dplyr::arrange(
    match(.data$contrast_id, Fig5_GSE131992_EXT_expected_contrasts)
  ) %>%
  tidyr::pivot_longer(
    cols = c(reference_timepoint, target_timepoint),
    names_to = "contrast_role",
    values_to = "geo_timepoint"
  ) %>%
  dplyr::mutate(
    comparison = factor(
      .data$contrast_id,
      levels = Fig5_GSE131992_EXT_expected_contrasts,
      labels = unname(Fig5_GSE131992_EXT_contrast_labels)
    ),
    contrast_role = factor(
      .data$contrast_role,
      levels = c("reference_timepoint", "target_timepoint"),
      labels = c("Reference state", "Target state")
    )
  )

if (
  nrow(Fig5I_GSE131992_EXT_quartile_contrast_states) != 8L ||
  anyNA(Fig5I_GSE131992_EXT_quartile_contrast_states$comparison) ||
  anyNA(Fig5I_GSE131992_EXT_quartile_contrast_states$contrast_role)
) {
  stop(
    "The GSE131992 fixed-quartile contrast-state map is incomplete.",
    call. = FALSE
  )
}

Fig5I_GSE131992_EXT_quartile_plot_data <-
  Fig5_GSE131992_EXT_objects$fixed_MFE_quartile_sample_fractions %>%
  dplyr::filter(
    .data$geo_timepoint %in% Fig5_GSE131992_EXT_selected_timepoints,
    .data$fixed_MFE_group != "Middle 50%"
  ) %>%
  dplyr::transmute(
    sample_id = as.character(.data$sample_id),
    geo_timepoint = as.character(.data$geo_timepoint),
    DGfold_quartile = factor(
      as.character(.data$fixed_MFE_group),
      levels = c(
        "Less-negative MFE (>75th percentile)",
        "More-negative MFE (<25th percentile)"
      ),
      labels = c(
        "High ΔGfold quartile",
        "Low ΔGfold quartile"
      )
    ),
    TPM_fraction = as.numeric(.data$TPM_fraction)
  ) %>%
  dplyr::inner_join(
    Fig5I_GSE131992_EXT_quartile_contrast_states,
    by = "geo_timepoint"
  ) %>%
  dplyr::arrange(
    .data$DGfold_quartile,
    .data$comparison,
    .data$contrast_role,
    .data$sample_id
  )

Fig5I_GSE131992_EXT_quartile_group_counts <-
  Fig5I_GSE131992_EXT_quartile_plot_data %>%
  dplyr::count(
    .data$comparison,
    .data$contrast_role,
    .data$DGfold_quartile,
    name = "N_biological_replicates"
  )

if (
  nrow(Fig5I_GSE131992_EXT_quartile_group_counts) != 16L ||
  !setequal(
    Fig5I_GSE131992_EXT_quartile_group_counts$N_biological_replicates,
    c(2L, 4L)
  ) ||
  any(!is.finite(Fig5I_GSE131992_EXT_quartile_plot_data$TPM_fraction))
) {
  stop(
    "The GSE131992 extended-recovery fixed-quartile data are incomplete.",
    call. = FALSE
  )
}

Fig5I_GSE131992_EXT_quartile_summary <-
  Fig5I_GSE131992_EXT_quartile_plot_data %>%
  dplyr::group_by(
    .data$comparison,
    .data$contrast_role,
    .data$DGfold_quartile
  ) %>%
  dplyr::summarise(
    N_biological_replicates = dplyr::n(),
    mean_fraction = mean(.data$TPM_fraction),
    SD_fraction = stats::sd(.data$TPM_fraction),
    SEM_fraction = .data$SD_fraction /
      sqrt(.data$N_biological_replicates),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    errorbar_low = pmax(
      0,
      .data$mean_fraction - .data$SEM_fraction
    ),
    errorbar_high = pmin(
      1,
      .data$mean_fraction + .data$SEM_fraction
    )
  )

if (
  nrow(Fig5I_GSE131992_EXT_quartile_summary) != 16L ||
  any(Fig5I_GSE131992_EXT_quartile_summary$N_biological_replicates < 2L) ||
  any(!is.finite(Fig5I_GSE131992_EXT_quartile_summary$SEM_fraction))
) {
  stop(
    "Panel I could not calculate all biological-replicate SEM values.",
    call. = FALSE
  )
}

Fig5I_GSE131992_EXT_quartile_y_max <- 1.20 * max(
  Fig5I_GSE131992_EXT_quartile_summary$errorbar_high,
  na.rm = TRUE
)

Fig5I_GSE131992_EXT_quartile_fractions <- ggplot(
  Fig5I_GSE131992_EXT_quartile_summary,
  aes(
    x = .data$comparison,
    y = .data$mean_fraction,
    fill = .data$contrast_role
  )
) +
  geom_col(
    position = position_dodge(width = 0.78),
    width = 0.68,
    colour = "black",
    linewidth = 0.25
  ) +
  geom_errorbar(
    aes(
      ymin = .data$mean_fraction,
      ymax = .data$errorbar_high
    ),
    position = position_dodge(width = 0.78),
    width = 0.12,
    linewidth = 0.35,
    colour = "black"
  ) +
  geom_text(
    aes(
      y = .data$errorbar_high,
      label = sprintf("%.3f", .data$mean_fraction)
    ),
    position = position_dodge(width = 0.78),
    vjust = -0.65,
    size = 2
  ) +
  facet_wrap(
    ~ DGfold_quartile,
    nrow = 1
  ) +
  coord_cartesian(
    ylim = c(0, Fig5I_GSE131992_EXT_quartile_y_max)
  ) +
  scale_fill_manual(
    values = c(
      "Reference state" = "#8cc5e3",
      "Target state" = "#2066a8"
    )
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "top",
    axis.title.x = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      colour = "black",
      size = 8
    ),
    axis.text.y = element_text(colour = "black"),
    strip.text = element_text(
      size = 5,
      face = "bold"
    )
  ) +
  labs(
    y = "Fraction of total retained TPM",
    fill = NULL
  )

print(Fig5I_GSE131992_EXT_quartile_fractions)

ggsave(
  filename = file.path(
    Fig5_output_dir,
    "Fig5I_GSE131992_EXT_fixed_MFE_quartile_TPM_fractions.tiff"
  ),
  plot = Fig5I_GSE131992_EXT_quartile_fractions,
  width = 6,
  height = 5,
  units = "in",
  dpi = 600,
  device = "tiff",
  compression = "lzw",
  bg = "white"
)


Fig5I_perm_B <- 10000L
Fig5I_perm_seed <- 20260918L

stopifnot(all(c(
  "fixed_MFE_quartiles", "relative_TPM_weight_matrices"
) %in% names(Fig5_GSE131992_EXT_objects)))

Fig5I_perm_weights <-
  Fig5_GSE131992_EXT_objects$relative_TPM_weight_matrices[[
    "all_finite_MFE"
  ]]

Fig5I_perm_groups <- Fig5_GSE131992_EXT_objects$fixed_MFE_quartiles

Fig5I_perm_plot <- as.data.frame(
  Fig5I_GSE131992_EXT_quartile_plot_data
)

Fig5I_perm_tail_labels <- c(
  "High ΔGfold quartile", "Low ΔGfold quartile"
)

stopifnot(
  is.matrix(Fig5I_perm_weights),
  is.numeric(Fig5I_perm_weights),
  !is.null(rownames(Fig5I_perm_weights)),
  !is.null(colnames(Fig5I_perm_weights)),
  !anyDuplicated(rownames(Fig5I_perm_weights)),
  !anyDuplicated(colnames(Fig5I_perm_weights)),
  !anyDuplicated(Fig5I_perm_groups$feature_id),
  setequal(rownames(Fig5I_perm_weights), Fig5I_perm_groups$feature_id),
  all(Fig5I_perm_plot$sample_id %in% colnames(Fig5I_perm_weights))
)

Fig5I_perm_weights <- Fig5I_perm_weights[
  , unique(Fig5I_perm_plot$sample_id), drop = FALSE
]

Fig5I_perm_category <- match(
  as.character(Fig5I_perm_groups$fixed_MFE_group[
    match(rownames(Fig5I_perm_weights), Fig5I_perm_groups$feature_id)
  ]),
  c(
    "Less-negative MFE (>75th percentile)",
    "More-negative MFE (<25th percentile)",
    "Middle 50%"
  )
)

stopifnot(
  !anyNA(Fig5I_perm_category),
  setequal(Fig5I_perm_category, 1:3),
  all(is.finite(Fig5I_perm_weights)),
  all(Fig5I_perm_weights >= 0),
  max(abs(colSums(Fig5I_perm_weights) - 1)) < 1e-10
)

Fig5I_perm_sample_fractions <- rbind(
  colSums(Fig5I_perm_weights[Fig5I_perm_category == 1L, , drop = FALSE]),
  colSums(Fig5I_perm_weights[Fig5I_perm_category == 2L, , drop = FALSE])
)

Fig5I_perm_sample_check <- Fig5I_perm_sample_fractions[cbind(
  match(
    as.character(Fig5I_perm_plot$DGfold_quartile),
    Fig5I_perm_tail_labels
  ),
  match(Fig5I_perm_plot$sample_id, colnames(Fig5I_perm_weights))
)]

stopifnot(
  all(is.finite(Fig5I_perm_sample_check)),
  all(is.finite(Fig5I_perm_plot$TPM_fraction)),
  max(abs(Fig5I_perm_sample_check - Fig5I_perm_plot$TPM_fraction)) < 1e-10,
  length(Fig5_GSE131992_EXT_expected_contrasts) == 4L,
  !anyDuplicated(Fig5_GSE131992_EXT_expected_contrasts)
)

Fig5I_perm_inputs <- lapply(
  Fig5_GSE131992_EXT_expected_contrasts, function(id) {
    
    dat <- Fig5I_perm_plot[Fig5I_perm_plot$contrast_id == id, ]
    
    reference <- unique(
      dat$sample_id[dat$contrast_role == "Reference state"]
    )
    target <- unique(
      dat$sample_id[dat$contrast_role == "Target state"]
    )
    comparison <- unique(as.character(dat$comparison))
    
    stopifnot(
      length(reference) >= 2L,
      length(target) >= 2L,
      length(intersect(reference, target)) == 0L,
      length(comparison) == 1L
    )
    
    state_weights <- cbind(
      rowMeans(Fig5I_perm_weights[, reference, drop = FALSE]),
      rowMeans(Fig5I_perm_weights[, target, drop = FALSE])
    )
    
    fractions <- vapply(1:2, function(k) {
      colSums(state_weights[Fig5I_perm_category == k, , drop = FALSE])
    }, numeric(2))
    
    old <- as.data.frame(Fig5I_GSE131992_EXT_quartile_summary)
    old <- old[as.character(old$comparison) == comparison, ]
    stopifnot(nrow(old) == 4L)
    
    expected <- vapply(Fig5I_perm_tail_labels, function(q) {
      rows <- old[as.character(old$DGfold_quartile) == q, ]
      stopifnot(nrow(rows) == 2L)
      
      rows$mean_fraction[match(
        c("Reference state", "Target state"),
        as.character(rows$contrast_role)
      )]
    }, numeric(2))
    
    if (any(!is.finite(expected)) ||
        max(abs(fractions - expected)) > 1e-10) {
      stop("Recalculated bar fractions do not match: ", id)
    }
    
    list(
      contrast_id = id,
      Comparison = gsub("\n", " ", comparison, fixed = TRUE),
      Reference_timepoint = unique(as.character(
        dat$geo_timepoint[dat$contrast_role == "Reference state"]
      )),
      Target_timepoint = unique(as.character(
        dat$geo_timepoint[dat$contrast_role == "Target state"]
      )),
      fractions = fractions,
      delta = state_weights[, 2] - state_weights[, 1]
    )
  }
)

message("All sample fractions and all 16 Fig. 5I bars match the source data.")

Fig5I_perm_recovery <- Filter(function(x) {
  identical(x$Reference_timepoint, "t19") &&
    x$Target_timepoint %in% c("t23", "t26", "t29")
}, Fig5I_perm_inputs)

stopifnot(
  length(Fig5I_perm_recovery) == 3L,
  setequal(vapply(Fig5I_perm_recovery, function(x) {
    x$Target_timepoint
  }, character(1)), c("t23", "t26", "t29"))
)

Fig5I_perm_inputs[[5L]] <- list(
  contrast_id = "Aggregate_recovery_minus_starvation",
  Comparison = "Mean recovery (25, 65, 155 min) vs starvation",
  fractions = Reduce("+", lapply(Fig5I_perm_recovery, function(x) {
    x$fractions
  })) / 3,
  delta = Reduce("+", lapply(Fig5I_perm_recovery, function(x) {
    x$delta
  })) / 3
)

Fig5I_perm_delta_weights <- do.call(
  cbind, lapply(Fig5I_perm_inputs, function(x) x$delta)
)

Fig5I_perm_observed <- vapply(Fig5I_perm_inputs, function(x) {
  x$fractions[2, ] - x$fractions[1, ]
}, numeric(2))

Fig5I_perm_extreme <- matrix(0L, nrow = 2L, ncol = 5L)

set.seed(Fig5I_perm_seed)

for (Fig5I_perm_b in seq_len(Fig5I_perm_B)) {
  
  Fig5I_perm_shuffled <- sample(Fig5I_perm_category)
  
  Fig5I_perm_null <- rbind(
    colSums(Fig5I_perm_delta_weights[
      Fig5I_perm_shuffled == 1L, , drop = FALSE
    ]),
    colSums(Fig5I_perm_delta_weights[
      Fig5I_perm_shuffled == 2L, , drop = FALSE
    ])
  )
  
  Fig5I_perm_extreme <- Fig5I_perm_extreme +
    (abs(Fig5I_perm_null) >= abs(Fig5I_perm_observed) - 1e-12)
}

Fig5I_perm_all_results <- do.call(rbind, lapply(
  seq_along(Fig5I_perm_inputs), function(i) {
    x <- Fig5I_perm_inputs[[i]]
    
    data.frame(
      Contrast = x$contrast_id,
      Comparison = x$Comparison,
      Quartile = c("Upper / less negative", "Lower / more negative"),
      Fraction_reference = x$fractions[1, ],
      Fraction_target = x$fractions[2, ],
      Delta_fraction = Fig5I_perm_observed[, i],
      P_permutation = (Fig5I_perm_extreme[, i] + 1) / (Fig5I_perm_B + 1)
    )
  }
))

Fig5I_format_results <- function(tab, prefix = "P") {
  tab$P_adj_BH <- p.adjust(tab$P_permutation, method = "BH")
  
  tab$Figure_label <- ifelse(
    tab$P_adj_BH >= 0.05,
    if (prefix == "P") "NS" else "Aggregate: NS",
    sprintf("%s = %.4f", prefix, tab$P_adj_BH)
  )
  
  rownames(tab) <- NULL
  tab
}

Fig5I_quartile_permutation_results <- Fig5I_format_results(
  Fig5I_perm_all_results[1:8, ]
)

Fig5I_aggregate_permutation_results <- Fig5I_format_results(
  Fig5I_perm_all_results[9:10, ], prefix = "Aggregate P"
)

Fig5I_print_columns <- c(
  "Comparison", "Quartile", "Fraction_reference", "Fraction_target",
  "Delta_fraction", "P_permutation", "P_adj_BH", "Figure_label"
)

print(
  Fig5I_quartile_permutation_results[, Fig5I_print_columns],
  row.names = FALSE,
  digits = 7
)

print(
  Fig5I_aggregate_permutation_results[, Fig5I_print_columns],
  row.names = FALSE,
  digits = 7
)


Fig5I_sampling_times <- c("t19", "t23", "t26", "t29")

Fig5I_sampling_groups <- c(
  "Upper / less negative" = "Less-negative MFE (>75th percentile)",
  "Lower / more negative" = "More-negative MFE (<25th percentile)"
)

Fig5I_sampling_source <- as.data.frame(
  Fig5_GSE131992_EXT_objects$fixed_MFE_quartile_sample_fractions
)

Fig5I_sampling_aggregate_results <- do.call(rbind, lapply(
  seq_along(Fig5I_sampling_groups), function(i) {
    
    dat <- Fig5I_sampling_source[
      Fig5I_sampling_source$fixed_MFE_group ==
        Fig5I_sampling_groups[[i]] &
        Fig5I_sampling_source$geo_timepoint %in%
        Fig5I_sampling_times,
    ]
    
    stopifnot(
      nrow(dat) == 8L,
      !anyNA(dat$sample_id),
      !anyDuplicated(dat$sample_id)
    )
    
    values <- lapply(Fig5I_sampling_times, function(tm) {
      x <- dat$TPM_fraction[dat$geo_timepoint == tm]
      
      stopifnot(
        length(x) == 2L,
        all(is.finite(x)),
        all(x >= 0 & x <= 1)
      )
      
      x
    })
    
    n <- lengths(values)
    means <- vapply(values, mean, numeric(1))
    variances <- vapply(values, stats::var, numeric(1))
    
    coefficients <- c(-1, 1/3, 1/3, 1/3)
    
    delta <- sum(coefficients * means)
    variance_components <- coefficients^2 * variances / n
    total_variance <- sum(variance_components)
    
    if (!is.finite(total_variance) || total_variance <= 0) {
      stop("Cannot estimate sampling uncertainty for ",
           names(Fig5I_sampling_groups)[i])
    }
    
    standard_error <- sqrt(total_variance)
    
    degrees_freedom <- total_variance^2 /
      sum(variance_components^2 / (n - 1))
    
    t_statistic <- delta / standard_error
    
    p_value <- 2 * stats::pt(
      -abs(t_statistic),
      df = degrees_freedom
    )
    
    margin <- stats::qt(0.975, df = degrees_freedom) *
      standard_error
    
    data.frame(
      Quartile = names(Fig5I_sampling_groups)[i],
      N_samples_per_timepoint = 2L,
      Mean_starvation = means[1],
      Mean_recovery = mean(means[2:4]),
      Delta_fraction = delta,
      SE_difference = standard_error,
      Degrees_freedom = degrees_freedom,
      CI95_lower = delta - margin,
      CI95_upper = delta + margin,
      P_sampling = p_value
    )
  }
))

Fig5I_sampling_aggregate_results$P_adj_BH <- p.adjust(
  Fig5I_sampling_aggregate_results$P_sampling,
  method = "BH"
)

Fig5I_sampling_match <- match(
  Fig5I_sampling_aggregate_results$Quartile,
  Fig5I_aggregate_permutation_results$Quartile
)

stopifnot(
  !anyNA(Fig5I_sampling_match),
  max(abs(
    Fig5I_sampling_aggregate_results$Delta_fraction -
      Fig5I_aggregate_permutation_results$Delta_fraction[
        Fig5I_sampling_match
      ]
  )) < 1e-10
)

print(
  Fig5I_sampling_aggregate_results,
  row.names = FALSE,
  digits = 7
)


Fig5I_recovery_times <- c(
  t23 = 25,
  t26 = 65,
  t29 = 155
)

Fig5I_recovery_hours <- unname(Fig5I_recovery_times) / 60

Fig5I_recovery_coefficients <-
  (Fig5I_recovery_hours - mean(Fig5I_recovery_hours)) /
  sum(
    (Fig5I_recovery_hours - mean(Fig5I_recovery_hours))^2
  )

Fig5I_recovery_trend_results <- do.call(
  rbind,
  lapply(seq_along(Fig5I_sampling_groups), function(i) {
    
    dat <- Fig5I_sampling_source[
      which(
        Fig5I_sampling_source$fixed_MFE_group ==
          Fig5I_sampling_groups[[i]] &
          Fig5I_sampling_source$geo_timepoint %in%
          names(Fig5I_recovery_times)
      ),
      ,
      drop = FALSE
    ]
    
    stopifnot(
      nrow(dat) == 6L,
      !anyNA(dat$sample_id),
      !anyDuplicated(dat$sample_id)
    )
    
    values <- lapply(names(Fig5I_recovery_times), function(tm) {
      x <- dat$TPM_fraction[dat$geo_timepoint == tm]
      
      stopifnot(
        length(x) == 2L,
        all(is.finite(x)),
        all(x >= 0 & x <= 1)
      )
      
      x
    })
    
    n <- lengths(values)
    means <- vapply(values, mean, numeric(1))
    variances <- vapply(values, stats::var, numeric(1))
    
    slope <- sum(Fig5I_recovery_coefficients * means)
    
    variance_components <-
      Fig5I_recovery_coefficients^2 * variances / n
    
    total_variance <- sum(variance_components)
    
    if (!is.finite(total_variance) || total_variance <= 0) {
      stop(
        "Cannot estimate recovery-trend uncertainty for ",
        names(Fig5I_sampling_groups)[i]
      )
    }
    
    standard_error <- sqrt(total_variance)
    
    degrees_freedom <- total_variance^2 /
      sum(variance_components^2 / (n - 1))
    
    p_value <- 2 * stats::pt(
      -abs(slope / standard_error),
      df = degrees_freedom
    )
    
    margin <- stats::qt(0.975, df = degrees_freedom) *
      standard_error
    
    data.frame(
      Quartile = names(Fig5I_sampling_groups)[i],
      N_samples_per_timepoint = 2L,
      Mean_25min = means[1],
      Mean_65min = means[2],
      Mean_155min = means[3],
      Slope_fraction_per_hour = slope,
      SE_slope = standard_error,
      Degrees_freedom = degrees_freedom,
      CI95_lower = slope - margin,
      CI95_upper = slope + margin,
      P_trend = p_value
    )
  })
)

Fig5I_recovery_trend_results$P_adj_BH <- stats::p.adjust(
  Fig5I_recovery_trend_results$P_trend,
  method = "BH"
)

Fig5I_recovery_trend_results$Figure_label <- with(
  Fig5I_recovery_trend_results,
  ifelse(
    P_adj_BH < 0.0001,
    "Trend P < 0.0001",
    ifelse(
      P_adj_BH < 0.05,
      sprintf("Trend P = %.4f", P_adj_BH),
      "Trend: NS"
    )
  )
)

print(
  Fig5I_recovery_trend_results,
  row.names = FALSE,
  digits = 7
)


Fig5J_GSE131992_EXT_weighted_sample_data <-
  Fig5_GSE131992_EXT_objects$weighted_sample_indices %>%
  dplyr::filter(
    .data$universe_id == "all_finite_MFE",
    .data$geo_timepoint %in% Fig5_GSE131992_EXT_selected_timepoints
  ) %>%
  dplyr::mutate(
    geo_timepoint = factor(
      as.character(.data$geo_timepoint),
      levels = Fig5_GSE131992_EXT_selected_timepoints
    )
  )

Fig5J_GSE131992_EXT_weighted_state_summary <-
  Fig5J_GSE131992_EXT_weighted_sample_data %>%
  dplyr::group_by(.data$geo_timepoint) %>%
  dplyr::summarise(
    N_biological_replicates = dplyr::n(),
    mean_index = mean(
      .data$TPM_weighted_Normalized_MFE
    ),
    SD_index = stats::sd(
      .data$TPM_weighted_Normalized_MFE
    ),
    SEM_index = .data$SD_index /
      sqrt(.data$N_biological_replicates),
    .groups = "drop"
  )

if (
  nrow(Fig5J_GSE131992_EXT_weighted_state_summary) != 5L ||
  !identical(
    Fig5J_GSE131992_EXT_weighted_state_summary$N_biological_replicates,
    c(4L, 2L, 2L, 2L, 2L)
  ) ||
  any(
    !is.finite(
      Fig5J_GSE131992_EXT_weighted_state_summary$SEM_index
    )
  )
) {
  stop(
    paste0(
      "Panel J requires 4/2/2/2/2 t01/t19/t23/t26/t29 biological replicates ",
      "with finite biological-replicate SEM values."
    ),
    call. = FALSE
  )
}

Fig5J_GSE131992_EXT_weighted_state_statistics <-
  Fig5J_GSE131992_EXT_weighted_state_summary %>%
  dplyr::mutate(
    geo_timepoint = as.character(.data$geo_timepoint)
  )

Fig5J_GSE131992_EXT_weighted_plot_data <-
  Fig5_GSE131992_EXT_objects$weighted_contrast_summary %>%
  dplyr::filter(
    .data$universe_id == "all_finite_MFE",
    .data$contrast_id %in% Fig5_GSE131992_EXT_expected_contrasts
  ) %>%
  dplyr::mutate(
    comparison = factor(
      as.character(.data$contrast_id),
      levels = Fig5_GSE131992_EXT_expected_contrasts,
      labels = unname(Fig5_GSE131992_EXT_contrast_labels)
    ),
    reference_SEM =
      Fig5J_GSE131992_EXT_weighted_state_statistics$SEM_index[
        match(
          as.character(.data$reference_timepoint),
          Fig5J_GSE131992_EXT_weighted_state_statistics$geo_timepoint
        )
      ],
    target_SEM =
      Fig5J_GSE131992_EXT_weighted_state_statistics$SEM_index[
        match(
          as.character(.data$target_timepoint),
          Fig5J_GSE131992_EXT_weighted_state_statistics$geo_timepoint
        )
      ],
    propagated_SEM = sqrt(
      .data$reference_SEM^2 + .data$target_SEM^2
    )
  ) %>%
  dplyr::arrange(.data$comparison)

if (
  nrow(Fig5J_GSE131992_EXT_weighted_plot_data) != 4L ||
  anyNA(Fig5J_GSE131992_EXT_weighted_plot_data$comparison) ||
  any(
    !is.finite(
      Fig5J_GSE131992_EXT_weighted_plot_data$
      delta_TPM_weighted_Normalized_MFE
    )
  ) ||
  any(
    !is.finite(
      Fig5J_GSE131992_EXT_weighted_plot_data$propagated_SEM
    )
  )
) {
  stop(
    "The GSE131992 complete-pool weighted contrasts are incomplete.",
    call. = FALSE
  )
}

Fig5J_GSE131992_EXT_master_delta_check <- max(
  abs(
    Fig5J_GSE131992_EXT_weighted_plot_data$
      delta_TPM_weighted_Normalized_MFE -
      (
        Fig5J_GSE131992_EXT_weighted_plot_data$target_mean_index -
          Fig5J_GSE131992_EXT_weighted_plot_data$reference_mean_index
      )
  ),
  na.rm = TRUE
)

if (Fig5J_GSE131992_EXT_master_delta_check > 1e-12) {
  stop(
    "Panel J deltas disagree with the saved condition means.",
    call. = FALSE
  )
}

Fig5J_GSE131992_EXT_weighted_plot_data <-
  Fig5J_GSE131992_EXT_weighted_plot_data %>%
  dplyr::mutate(
    errorbar_low = dplyr::if_else(
      .data$delta_TPM_weighted_Normalized_MFE >= 0,
      .data$delta_TPM_weighted_Normalized_MFE,
      .data$delta_TPM_weighted_Normalized_MFE -
        .data$propagated_SEM
    ),
    errorbar_high = dplyr::if_else(
      .data$delta_TPM_weighted_Normalized_MFE >= 0,
      .data$delta_TPM_weighted_Normalized_MFE +
        .data$propagated_SEM,
      .data$delta_TPM_weighted_Normalized_MFE
    )
  )

Fig5J_GSE131992_EXT_weighted_y_limit <- 1.40 * max(
  abs(
    c(
      Fig5J_GSE131992_EXT_weighted_plot_data$errorbar_low,
      Fig5J_GSE131992_EXT_weighted_plot_data$errorbar_high
    )
  ),
  na.rm = TRUE
)

Fig5J_GSE131992_EXT_weighted_label_offset <-
  0.05 * Fig5J_GSE131992_EXT_weighted_y_limit

Fig5J_GSE131992_EXT_weighted_plot_data <-
  Fig5J_GSE131992_EXT_weighted_plot_data %>%
  dplyr::mutate(
    label_y = dplyr::if_else(
      .data$delta_TPM_weighted_Normalized_MFE >= 0,
      .data$errorbar_high +
        Fig5J_GSE131992_EXT_weighted_label_offset,
      .data$errorbar_low -
        Fig5J_GSE131992_EXT_weighted_label_offset
    ),
    label_vjust = dplyr::if_else(
      .data$delta_TPM_weighted_Normalized_MFE >= 0,
      0,
      1
    )
  )

Fig5J_GSE131992_EXT_weighted_DeltaDeltaGfold <- ggplot(
  Fig5J_GSE131992_EXT_weighted_plot_data,
  aes(
    x = .data$comparison,
    y = .data$delta_TPM_weighted_Normalized_MFE,
    fill = .data$comparison
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.4,
    colour = "black"
  ) +
  geom_col(
    width = 0.65,
    colour = "black",
    linewidth = 0.3
  ) +
  geom_errorbar(
    aes(
      ymin = .data$errorbar_low,
      ymax = .data$errorbar_high
    ),
    width = 0.16,
    linewidth = 0.4,
    colour = "black"
  ) +
  geom_text(
    aes(
      y = .data$label_y,
      label = sprintf(
        "%+.5f",
        .data$delta_TPM_weighted_Normalized_MFE
      ),
      vjust = .data$label_vjust
    ),
    size = 2.8
  ) +
  scale_fill_manual(
    values = Fig5_GSE131992_EXT_contrast_colours
  ) +
  scale_y_continuous(
    labels = scales::number_format(accuracy = 0.001)
  ) +
  coord_cartesian(
    ylim = c(
      -Fig5J_GSE131992_EXT_weighted_y_limit,
      Fig5J_GSE131992_EXT_weighted_y_limit
    ),
    clip = "off"
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      colour = "black",
      size = 8
    ),
    axis.text.y = element_text(colour = "black"),
    plot.caption = element_text(
      hjust = 0,
      size = 7
    )
  ) +
  labs(
    y = expression(
      "Abundance-weighted " * Delta * Delta * G[fold] * " / nt"
    ),
    caption = paste0(
      "Recovery bars use the common t19 anchor. Error bars are propagated ",
      "SEM across independent biological replicates"
    )
  )

print(Fig5J_GSE131992_EXT_weighted_DeltaDeltaGfold)

ggsave(
  filename = file.path(
    Fig5_output_dir,
    "Fig5J_GSE131992_EXT_complete_pool_weighted_DeltaDeltaGfold.tiff"
  ),
  plot = Fig5J_GSE131992_EXT_weighted_DeltaDeltaGfold,
  width = 4,
  height = 5,
  units = "in",
  dpi = 600,
  device = "tiff",
  compression = "lzw",
  bg = "white"
)
Fig5J_GSE131992_EXT_required_aggregate_objects <- c(
  "sample_manifest",
  "exact_MFE",
  "relative_TPM_weight_matrices",
  "TPM_mass_coverage",
  "weighted_sample_indices"
)

Fig5J_GSE131992_EXT_missing_aggregate_objects <- setdiff(
  Fig5J_GSE131992_EXT_required_aggregate_objects,
  names(Fig5_GSE131992_EXT_objects)
)

if (length(Fig5J_GSE131992_EXT_missing_aggregate_objects) > 0L) {
  stop(
    paste0(
      "The GSE131992 RDS lacks object(s) required for the aggregate ",
      "descriptor permutation: ",
      paste(
        Fig5J_GSE131992_EXT_missing_aggregate_objects,
        collapse = ", "
      )
    ),
    call. = FALSE
  )
}


Fig5J_GSE131992_EXT_recovery_timepoints <- c(
  "t23",
  "t26",
  "t29"
)

Fig5J_GSE131992_EXT_aggregate_timepoints <- c(
  "t19",
  Fig5J_GSE131992_EXT_recovery_timepoints
)

Fig5J_GSE131992_EXT_aggregate_manifest <-
  Fig5_GSE131992_EXT_objects$sample_manifest %>%
  dplyr::filter(
    .data$geo_timepoint %in%
      Fig5J_GSE131992_EXT_aggregate_timepoints
  ) %>%
  dplyr::transmute(
    sample_id = as.character(.data$sample_id),
    geo_timepoint = as.character(.data$geo_timepoint),
    biological_replicate = as.character(
      .data$biological_replicate
    )
  )

Fig5J_GSE131992_EXT_timepoint_counts <- table(
  factor(
    Fig5J_GSE131992_EXT_aggregate_manifest$geo_timepoint,
    levels = Fig5J_GSE131992_EXT_aggregate_timepoints
  )
)

if (
  !identical(
    as.integer(Fig5J_GSE131992_EXT_timepoint_counts),
    c(2L, 2L, 2L, 2L)
  ) ||
  dplyr::n_distinct(
    Fig5J_GSE131992_EXT_aggregate_manifest$biological_replicate
  ) != 2L ||
  anyDuplicated(
    paste(
      Fig5J_GSE131992_EXT_aggregate_manifest$biological_replicate,
      Fig5J_GSE131992_EXT_aggregate_manifest$geo_timepoint,
      sep = "__"
    )
  ) > 0L
) {
  stop(
    paste0(
      "The aggregate test requires two matched biological replicates ",
      "at each of t19, t23, t26 and t29."
    ),
    call. = FALSE
  )
}

Fig5J_GSE131992_EXT_replicates <- sort(
  unique(
    Fig5J_GSE131992_EXT_aggregate_manifest$biological_replicate
  )
)


Fig5J_GSE131992_EXT_weight_matrix <-
  Fig5_GSE131992_EXT_objects$relative_TPM_weight_matrices[[
    "all_finite_MFE"
  ]]

if (
  is.null(Fig5J_GSE131992_EXT_weight_matrix) ||
  is.null(rownames(Fig5J_GSE131992_EXT_weight_matrix)) ||
  is.null(colnames(Fig5J_GSE131992_EXT_weight_matrix)) ||
  any(
    !Fig5J_GSE131992_EXT_aggregate_manifest$sample_id %in%
    colnames(Fig5J_GSE131992_EXT_weight_matrix)
  )
) {
  stop(
    "The all-finite-MFE relative-TPM weight matrix is invalid.",
    call. = FALSE
  )
}

Fig5J_GSE131992_EXT_delta_weights_by_replicate <- vapply(
  Fig5J_GSE131992_EXT_replicates,
  function(Fig5J_GSE131992_EXT_replicate_id) {
    
    Fig5J_GSE131992_EXT_t19_sample <-
      Fig5J_GSE131992_EXT_aggregate_manifest$sample_id[
        Fig5J_GSE131992_EXT_aggregate_manifest$
          biological_replicate ==
          Fig5J_GSE131992_EXT_replicate_id &
          Fig5J_GSE131992_EXT_aggregate_manifest$
          geo_timepoint == "t19"
      ]
    
    Fig5J_GSE131992_EXT_recovery_samples <-
      Fig5J_GSE131992_EXT_aggregate_manifest$sample_id[
        Fig5J_GSE131992_EXT_aggregate_manifest$
          biological_replicate ==
          Fig5J_GSE131992_EXT_replicate_id &
          Fig5J_GSE131992_EXT_aggregate_manifest$
          geo_timepoint %in%
          Fig5J_GSE131992_EXT_recovery_timepoints
      ]
    
    if (
      length(Fig5J_GSE131992_EXT_t19_sample) != 1L ||
      length(Fig5J_GSE131992_EXT_recovery_samples) != 3L
    ) {
      stop(
        "Incomplete trajectory for biological replicate ",
        Fig5J_GSE131992_EXT_replicate_id,
        ".",
        call. = FALSE
      )
    }
    
    rowMeans(
      Fig5J_GSE131992_EXT_weight_matrix[
        ,
        Fig5J_GSE131992_EXT_recovery_samples,
        drop = FALSE
      ]
    ) -
      Fig5J_GSE131992_EXT_weight_matrix[
        ,
        Fig5J_GSE131992_EXT_t19_sample
      ]
  },
  numeric(nrow(Fig5J_GSE131992_EXT_weight_matrix))
)

Fig5J_GSE131992_EXT_aggregate_delta_weight <- rowMeans(
  Fig5J_GSE131992_EXT_delta_weights_by_replicate
)


Fig5J_GSE131992_EXT_feature_match <- match(
  rownames(Fig5J_GSE131992_EXT_weight_matrix),
  Fig5_GSE131992_EXT_objects$exact_MFE$feature_id
)

if (anyNA(Fig5J_GSE131992_EXT_feature_match)) {
  stop(
    "The weight matrix cannot be aligned to exact_MFE.",
    call. = FALSE
  )
}

Fig5J_GSE131992_EXT_normalized_MFE <- as.numeric(
  Fig5_GSE131992_EXT_objects$exact_MFE$Normalized_MFE[
    Fig5J_GSE131992_EXT_feature_match
  ]
)

Fig5J_GSE131992_EXT_feature_length <- as.numeric(
  Fig5_GSE131992_EXT_objects$exact_MFE$feature_length[
    Fig5J_GSE131992_EXT_feature_match
  ]
)

if (
  any(!is.finite(Fig5J_GSE131992_EXT_normalized_MFE)) ||
  any(!is.finite(Fig5J_GSE131992_EXT_feature_length)) ||
  any(Fig5J_GSE131992_EXT_feature_length <= 0) ||
  abs(sum(Fig5J_GSE131992_EXT_aggregate_delta_weight)) > 1e-10
) {
  stop(
    "Invalid feature annotation or aggregate contrast weights.",
    call. = FALSE
  )
}

Fig5J_GSE131992_EXT_observed_aggregate_delta <- sum(
  Fig5J_GSE131992_EXT_normalized_MFE *
    Fig5J_GSE131992_EXT_aggregate_delta_weight
)


Fig5J_GSE131992_EXT_TPM_coverage <-
  Fig5_GSE131992_EXT_objects$TPM_mass_coverage %>%
  dplyr::filter(
    .data$universe_id == "all_finite_MFE"
  ) %>%
  dplyr::transmute(
    sample_id = as.character(.data$sample_id),
    retained_TPM = as.numeric(.data$retained_TPM)
  )

Fig5J_GSE131992_EXT_retained_TPM <- stats::setNames(
  Fig5J_GSE131992_EXT_TPM_coverage$retained_TPM,
  Fig5J_GSE131992_EXT_TPM_coverage$sample_id
)

Fig5J_GSE131992_EXT_retained_TPM <- unname(
  Fig5J_GSE131992_EXT_retained_TPM[
    colnames(Fig5J_GSE131992_EXT_weight_matrix)
  ]
)

if (
  length(Fig5J_GSE131992_EXT_retained_TPM) !=
  ncol(Fig5J_GSE131992_EXT_weight_matrix) ||
  any(!is.finite(Fig5J_GSE131992_EXT_retained_TPM)) ||
  any(Fig5J_GSE131992_EXT_retained_TPM <= 0)
) {
  stop(
    "The retained-TPM totals cannot be aligned to the weight matrix.",
    call. = FALSE
  )
}

Fig5J_GSE131992_EXT_reconstructed_TPM <- sweep(
  Fig5J_GSE131992_EXT_weight_matrix,
  2L,
  Fig5J_GSE131992_EXT_retained_TPM,
  "*"
)

Fig5J_GSE131992_EXT_mean_TPM <- rowMeans(
  Fig5J_GSE131992_EXT_reconstructed_TPM
)

Fig5J_GSE131992_EXT_permutation_strata <- paste(
  dplyr::ntile(
    log10(Fig5J_GSE131992_EXT_feature_length),
    5L
  ),
  dplyr::ntile(
    log10(Fig5J_GSE131992_EXT_mean_TPM + 1e-8),
    5L
  ),
  sep = "_"
)

Fig5J_GSE131992_EXT_stratum_indices <- split(
  seq_along(Fig5J_GSE131992_EXT_normalized_MFE),
  Fig5J_GSE131992_EXT_permutation_strata
)


Fig5J_GSE131992_EXT_N_permutations <- 10000L
set.seed(140995L)

Fig5J_GSE131992_EXT_permutation_null <- numeric(
  Fig5J_GSE131992_EXT_N_permutations
)

Fig5J_GSE131992_EXT_permuted_MFE <-
  Fig5J_GSE131992_EXT_normalized_MFE

for (
  Fig5J_GSE131992_EXT_permutation_index in
  seq_len(Fig5J_GSE131992_EXT_N_permutations)
) {
  for (
    Fig5J_GSE131992_EXT_one_stratum in
    Fig5J_GSE131992_EXT_stratum_indices
  ) {
    Fig5J_GSE131992_EXT_permuted_MFE[
      Fig5J_GSE131992_EXT_one_stratum
    ] <- Fig5J_GSE131992_EXT_normalized_MFE[
      Fig5J_GSE131992_EXT_one_stratum
    ][
      sample.int(
        length(Fig5J_GSE131992_EXT_one_stratum)
      )
    ]
  }
  
  Fig5J_GSE131992_EXT_permutation_null[
    Fig5J_GSE131992_EXT_permutation_index
  ] <- sum(
    Fig5J_GSE131992_EXT_permuted_MFE *
      Fig5J_GSE131992_EXT_aggregate_delta_weight
  )
}

Fig5J_GSE131992_EXT_permutation_null_mean <- mean(
  Fig5J_GSE131992_EXT_permutation_null
)

Fig5J_GSE131992_EXT_descriptor_permutation_two_sided_P <-
  (
    1L + sum(
      abs(
        Fig5J_GSE131992_EXT_permutation_null -
          Fig5J_GSE131992_EXT_permutation_null_mean
      ) >=
        abs(
          Fig5J_GSE131992_EXT_observed_aggregate_delta -
            Fig5J_GSE131992_EXT_permutation_null_mean
        )
    )
  ) /
  (Fig5J_GSE131992_EXT_N_permutations + 1L)


Fig5J_GSE131992_EXT_paired_data <-
  Fig5_GSE131992_EXT_objects$weighted_sample_indices %>%
  dplyr::filter(
    .data$universe_id == "all_finite_MFE",
    .data$geo_timepoint %in%
      Fig5J_GSE131992_EXT_aggregate_timepoints
  ) %>%
  dplyr::transmute(
    biological_replicate = as.character(
      .data$biological_replicate
    ),
    geo_timepoint = as.character(.data$geo_timepoint),
    weighted_index = as.numeric(
      .data$TPM_weighted_Normalized_MFE
    )
  ) %>%
  tidyr::pivot_wider(
    names_from = geo_timepoint,
    values_from = weighted_index
  ) %>%
  dplyr::arrange(.data$biological_replicate)

if (
  nrow(Fig5J_GSE131992_EXT_paired_data) != 2L ||
  anyNA(
    Fig5J_GSE131992_EXT_paired_data[
      Fig5J_GSE131992_EXT_aggregate_timepoints
    ]
  )
) {
  stop(
    "The two matched t19/t23/t26/t29 trajectories are incomplete.",
    call. = FALSE
  )
}

Fig5J_GSE131992_EXT_paired_data$mean_recovery_index <- rowMeans(
  as.matrix(
    Fig5J_GSE131992_EXT_paired_data[
      Fig5J_GSE131992_EXT_recovery_timepoints
    ]
  )
)

Fig5J_GSE131992_EXT_paired_data$
  paired_delta_recovery_minus_t19 <-
  Fig5J_GSE131992_EXT_paired_data$mean_recovery_index -
  Fig5J_GSE131992_EXT_paired_data$t19

Fig5J_GSE131992_EXT_paired_test <- stats::t.test(
  Fig5J_GSE131992_EXT_paired_data$mean_recovery_index,
  Fig5J_GSE131992_EXT_paired_data$t19,
  paired = TRUE,
  alternative = "two.sided"
)


Fig5J_GSE131992_EXT_paired_test_one_sided <- stats::t.test(
  Fig5J_GSE131992_EXT_paired_data$mean_recovery_index,
  Fig5J_GSE131992_EXT_paired_data$t19,
  paired = TRUE,
  alternative = "less"
)


Fig5J_GSE131992_EXT_delta_check <- abs(
  Fig5J_GSE131992_EXT_observed_aggregate_delta -
    mean(
      Fig5J_GSE131992_EXT_paired_data$
        paired_delta_recovery_minus_t19
    )
)

if (Fig5J_GSE131992_EXT_delta_check > 1e-10) {
  stop(
    "The gene-weight and paired-sample aggregate deltas disagree.",
    call. = FALSE
  )
}


Fig5J_GSE131992_EXT_paired_delta_SD <- stats::sd(
  Fig5J_GSE131992_EXT_paired_data$
    paired_delta_recovery_minus_t19
)

Fig5J_GSE131992_EXT_paired_delta_SEM <-
  Fig5J_GSE131992_EXT_paired_delta_SD /
  sqrt(nrow(Fig5J_GSE131992_EXT_paired_data))

Fig5J_GSE131992_EXT_aggregate_statistics <- tibble::tibble(
  contrast_id =
    "Aggregate_recovery_t23_t26_t29_minus_t19",
  reference_timepoint = "t19",
  recovery_timepoints = "t23,t26,t29",
  N_independent_biological_trajectories =
    nrow(Fig5J_GSE131992_EXT_paired_data),
  N_recovery_measurements_per_trajectory =
    length(Fig5J_GSE131992_EXT_recovery_timepoints),
  reference_mean_index = mean(
    Fig5J_GSE131992_EXT_paired_data$t19
  ),
  aggregate_recovery_mean_index = mean(
    Fig5J_GSE131992_EXT_paired_data$mean_recovery_index
  ),
  aggregate_delta_recovery_minus_t19 =
    Fig5J_GSE131992_EXT_observed_aggregate_delta,
  paired_delta_SD =
    Fig5J_GSE131992_EXT_paired_delta_SD,
  paired_delta_SEM =
    Fig5J_GSE131992_EXT_paired_delta_SEM,
  paired_t = unname(
    Fig5J_GSE131992_EXT_paired_test$statistic
  ),
  paired_df = unname(
    Fig5J_GSE131992_EXT_paired_test$parameter
  ),
  paired_two_sided_P =
    Fig5J_GSE131992_EXT_paired_test$p.value,
  paired_one_sided_P =
    Fig5J_GSE131992_EXT_paired_test_one_sided$p.value,
  descriptor_permutations =
    Fig5J_GSE131992_EXT_N_permutations,
  descriptor_null_mean =
    Fig5J_GSE131992_EXT_permutation_null_mean,
  descriptor_null_SD = stats::sd(
    Fig5J_GSE131992_EXT_permutation_null
  ),
  descriptor_permutation_two_sided_P =
    Fig5J_GSE131992_EXT_descriptor_permutation_two_sided_P,
  descriptor_permutation_significant_at_0.05 =
    Fig5J_GSE131992_EXT_descriptor_permutation_two_sided_P < 0.05
)

cat(
  paste0(
    "\nGSE131992 Fig. 5J aggregate recovery versus t19:\n",
    "Recovery is the within-trajectory mean of t23, t26 and t29.\n",
    "The paired P value evaluates two biological trajectories.\n",
    "The descriptor-permutation P value evaluates transcriptome-wide ",
    "alignment with normalized MFE.\n"
  )
)

print(
  as.data.frame(Fig5J_GSE131992_EXT_paired_data),
  row.names = FALSE,
  digits = 12
)

print(
  as.data.frame(Fig5J_GSE131992_EXT_aggregate_statistics),
  row.names = FALSE,
  digits = 12
)


message(
  "\nFigure 5 analysis complete.\n",
  "GSE224463 input: ", Fig5_GSE224463_analysis_rds_file, "\n",
  "GSE226643 input: ", Fig5_GSE226643_analysis_rds_file, "\n",
  "GSE131992 input: ", Fig5_GSE131992_analysis_rds_file, "\n",
  "Output directory: ", normalizePath(Fig5_output_dir, mustWork = FALSE)
)
