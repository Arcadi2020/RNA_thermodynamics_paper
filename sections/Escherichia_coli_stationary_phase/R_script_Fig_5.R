suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(tibble)
})

arguments <- commandArgs(trailingOnly = TRUE)

if (length(arguments) > 0L && length(arguments) < 3L) {
  stop(
    paste(
      "Provide all three arguments:",
      "GSE224463 RDS, GSE226643 RDS, and output directory."
    ),
    call. = FALSE
  )
}

Fig5_GSE224463_analysis_rds_file <- if (length(arguments) >= 3L) {
  arguments[[1]]
} else {
  "../../data/GSE224463_analysis_results_for_Fig_5.rds"
}

Fig5_GSE226643_analysis_rds_file <- if (length(arguments) >= 3L) {
  arguments[[2]]
} else {
  "../../data/GSE226643_analysis_results_for_Fig_5.rds"
}

Fig5_output_dir <- if (length(arguments) >= 3L) {
  arguments[[3]]
} else {
  "Figure_5_outputs"
}

dir.create(Fig5_output_dir, recursive = TRUE, showWarnings = FALSE)

Fig5_missing_input_files <- c(
  Fig5_GSE224463_analysis_rds_file,
  Fig5_GSE226643_analysis_rds_file
)[!file.exists(c(
  Fig5_GSE224463_analysis_rds_file,
  Fig5_GSE226643_analysis_rds_file
))]

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

if (!is.list(Fig5_GSE224463_objects) || !is.list(Fig5_GSE226643_objects)) {
  stop("Both Figure 5 input RDS files must contain named lists.", call. = FALSE)
}

Fig5_required_GSE224463_objects <- c(
  "contrast_manifest",
  "GSEA_bubble_data",
  "significant_DE_thermodynamic_transcripts",
  "significant_DE_distribution_statistics",
  "significant_DE_abundance_fractions",
  "weighted_sample_indices",
  "full_transcriptome_weighted_delta_summary",
  "full_transcriptome_weighted_gene_contributions"
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
  Fig5_GSE226643_MFE_violin_data %>%
  dplyr::filter(
    .data$limitation == "Sulfur limitation"
  )

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
  filename = file.path(
    Fig5_output_dir,
    "Fig5_GSE224463_long_series_significant_DE_DGfold_abundance_fractions.tiff"
  ),
  plot = Fig5_Fraction_GSE224463_LongSeries_block,
  width = 6,
  height = 5,
  units = "in",
  dpi = 600,
  device = "tiff",
  compression = "lzw",
  bg = "white"
)

Fig5_GSE224463_fraction_z_tests <-
  Fig5_GSE224463_objects$significant_DE_thermodynamic_transcripts %>%
  dplyr::filter(
    as.character(.data$series) == "long",
    .data$target_time_min %in% Fig5_GSE224463_fraction_timepoints
  ) %>%
  dplyr::group_by(.data$analysis_id) %>%
  dplyr::group_modify(~ {

    q25 <- quantile(.x$Normalized_MFE, 0.25, na.rm = TRUE)
    q75 <- quantile(.x$Normalized_MFE, 0.75, na.rm = TRUE)

    x <- .x %>%
      dplyr::mutate(
        MFE_group = dplyr::case_when(
          Normalized_MFE <= q25 ~ "Low ΔGfold quartile",
          Normalized_MFE >= q75 ~ "High ΔGfold quartile",
          TRUE ~ "Intermediate"
        )
      )

    n1 <- sum(x$reference_0min_mean_normalised_counts, na.rm = TRUE)
    n2 <- sum(x$target_time_mean_normalised_counts, na.rm = TRUE)

    dplyr::bind_rows(
      lapply(
        c("High ΔGfold quartile", "Low ΔGfold quartile"),
        function(q) {

          x1 <- sum(
            x$reference_0min_mean_normalised_counts[x$MFE_group == q],
            na.rm = TRUE
          )

          x2 <- sum(
            x$target_time_mean_normalised_counts[x$MFE_group == q],
            na.rm = TRUE
          )

          p1 <- x1 / n1
          p2 <- x2 / n2
          p_pool <- (x1 + x2) / (n1 + n2)

          se <- sqrt(
            p_pool * (1 - p_pool) *
              ((1 / n1) + (1 / n2))
          )

          z <- (p2 - p1) / se

          data.frame(
            DGfold_quartile = q,
            Fraction_0min = p1,
            Fraction_target = p2,
            Delta_target_minus_0min = p2 - p1,
            Z = z,
            P = 2 * (1 - pnorm(abs(z)))
          )
        }
      )
    )
  }) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(
    Fig5_GSE224463_fraction_manifest %>%
      dplyr::select(
        .data$analysis_id,
        .data$target_time_min
      ),
    by = "analysis_id"
  ) %>%
  dplyr::mutate(
    P_BH = p.adjust(P, method = "BH")
  ) %>%
  dplyr::select(
    target_time_min,
    DGfold_quartile,
    Fraction_0min,
    Fraction_target,
    Delta_target_minus_0min,
    Z,
    P,
    P_BH
  ) %>%
  dplyr::arrange(
    target_time_min,
    factor(
      DGfold_quartile,
      levels = c(
        "High ΔGfold quartile",
        "Low ΔGfold quartile"
      )
    )
  )

print(Fig5_GSE224463_fraction_z_tests)

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
  filename = file.path(
    Fig5_output_dir,
    "Fig5_GSE226643_N_S_significant_DE_DGfold_abundance_fractions.tiff"
  ),
  plot = Fig5_Fraction_GSE226643_Nitrogen_Sulfur_block,
  width = 3.5,
  height = 5,
  units = "in",
  dpi = 600,
  device = "tiff",
  compression = "lzw",
  bg = "white"
)

Fig5E_GSE226643_fraction_z_tests <- dplyr::bind_rows(

  Fig5_GSE226643_objects$GSE226643_N_sig_thermodynamic %>%
    dplyr::mutate(limitation = "Nitrogen limitation"),

  Fig5_GSE226643_objects$GSE226643_S_sig_thermodynamic %>%
    dplyr::mutate(limitation = "Sulfur limitation")

) %>%
  dplyr::group_by(limitation) %>%
  dplyr::mutate(
    q25 = quantile(Normalized_MFE, 0.25, na.rm = TRUE),
    q75 = quantile(Normalized_MFE, 0.75, na.rm = TRUE),
    MFE_group = dplyr::case_when(
      Normalized_MFE <= q25 ~ "Low ΔGfold quartile",
      Normalized_MFE >= q75 ~ "High ΔGfold quartile",
      TRUE ~ "Intermediate"
    )
  ) %>%
  dplyr::group_modify(~ {

    total_early <- sum(.x$early_3h_mean_normalised_counts, na.rm = TRUE)
    total_stationary <- sum(.x$stationary_12h_mean_normalised_counts, na.rm = TRUE)

    dplyr::bind_rows(
      lapply(
        c("High ΔGfold quartile", "Low ΔGfold quartile"),
        function(q) {

          x1 <- sum(
            .x$early_3h_mean_normalised_counts[.x$MFE_group == q],
            na.rm = TRUE
          )

          x2 <- sum(
            .x$stationary_12h_mean_normalised_counts[.x$MFE_group == q],
            na.rm = TRUE
          )

          p1 <- x1 / total_early
          p2 <- x2 / total_stationary

          p_pool <- (x1 + x2) / (total_early + total_stationary)

          se <- sqrt(
            p_pool * (1 - p_pool) *
              ((1 / total_early) + (1 / total_stationary))
          )

          z <- (p2 - p1) / se

          data.frame(
            DGfold_quartile = q,
            Fraction_3h = p1,
            Fraction_12h = p2,
            Delta_12h_minus_3h = p2 - p1,
            Z = z,
            P = 2 * (1 - pnorm(abs(z)))
          )
        }
      )
    )
  }) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    P_BH = p.adjust(P, method = "BH")
  )

print(Fig5E_GSE226643_fraction_z_tests)

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

message(
  "\nFigure 5 analysis complete.\n",
  "GSE224463 input: ", Fig5_GSE224463_analysis_rds_file, "\n",
  "GSE226643 input: ", Fig5_GSE226643_analysis_rds_file, "\n",
  "Output directory: ", normalizePath(Fig5_output_dir, mustWork = FALSE)
)
