suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(pheatmap)
})

arguments <- commandArgs(trailingOnly = TRUE)

figure4_data_file <- if (length(arguments) >= 1L) {
  arguments[[1]]
} else {
  "path_to/iLINCS_analysis_results_for_Fig_4.rds"
}

output_dir <- if (length(arguments) >= 2L) {
  arguments[[2]]
} else {
  "path_to/Figure_4_outputs"
}

if (!file.exists(figure4_data_file)) {
  stop("Figure 4 data file not found: ", figure4_data_file, call. = FALSE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

figure4_data <- readRDS(figure4_data_file)

required_objects <- c(
  "analysis_3.33uM",
  "gsea_3.33uM",
  "violin_3.33uM",
  "analysis_10uM",
  "gsea_10uM",
  "violin_10uM",
  "zweighted"
)

if (!is.list(figure4_data) || is.null(names(figure4_data))) {
  stop("The Figure 4 RDS must contain a named list.", call. = FALSE)
}

missing_objects <- setdiff(required_objects, names(figure4_data))

if (length(missing_objects) > 0L) {
  stop(
    "The Figure 4 RDS is missing: ",
    paste(missing_objects, collapse = ", "),
    call. = FALSE
  )
}

validate_table <- function(data, label, required_columns) {
  if (!is.data.frame(data) || nrow(data) == 0L) {
    stop(label, " must be a non-empty data frame.", call. = FALSE)
  }

  missing_columns <- setdiff(required_columns, colnames(data))

  if (length(missing_columns) > 0L) {
    stop(
      label,
      " is missing: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

analysis_3.33uM <- figure4_data[["analysis_3.33uM"]]
gsea_3.33uM <- figure4_data[["gsea_3.33uM"]]
violin_3.33uM <- figure4_data[["violin_3.33uM"]]
analysis_10uM <- figure4_data[["analysis_10uM"]]
gsea_10uM <- figure4_data[["gsea_10uM"]]
violin_10uM <- figure4_data[["violin_10uM"]]
zweighted <- figure4_data[["zweighted"]]

analysis_columns <- c(
  "SampleName",
  "cell_line",
  "Drug",
  "MeanMFE_UP",
  "MeanMFE_DOWN",
  "wilcox_test_all_UP_vs_DOWN_pValue",
  "Log10_MFE_change"
)

gsea_columns <- c("ID", "Description", "Sample", "NES")
violin_columns <- c("SampleName", "z_score", "Normalized_MFE")
zweighted_columns <- c(
  "SampleName",
  "dose",
  "cell_line",
  "Drug",
  "delta_z_weighted_Normalized_MFE_UP_minus_DOWN",
  "permutation_p_delta_UP_minus_DOWN"
)

validate_table(analysis_3.33uM, "3.33 µM analysis", analysis_columns)
validate_table(gsea_3.33uM, "3.33 µM GSEA", gsea_columns)
validate_table(violin_3.33uM, "3.33 µM violin data", violin_columns)
validate_table(analysis_10uM, "10 µM analysis", analysis_columns)
validate_table(gsea_10uM, "10 µM GSEA", gsea_columns)
validate_table(violin_10uM, "10 µM violin data", violin_columns)
validate_table(zweighted, "z-score-weighted data", zweighted_columns)

make_heatmap_matrix <- function(data, row_order) {
  data <- data %>%
    dplyr::filter(
      is.finite(NES),
      !is.na(Description),
      Description != "",
      !is.na(Sample),
      Sample != ""
    ) %>%
    dplyr::distinct(Description, Sample, .keep_all = TRUE)

  sample_order <- unique(data$Sample)

  if (length(sample_order) == 0L) {
    stop("No GSEA samples were available for the heatmap.", call. = FALSE)
  }

  matrix_data <- matrix(
    0,
    nrow = length(row_order),
    ncol = length(sample_order),
    dimnames = list(row_order, sample_order)
  )

  row_position <- match(data$Description, row_order)
  column_position <- match(data$Sample, sample_order)
  keep <- !is.na(row_position) & !is.na(column_position)

  matrix_data[cbind(row_position[keep], column_position[keep])] <-
    data$NES[keep]

  matrix_data
}

ordered_pathways_3.33uM <- c(
  "AMPK signaling",
  "Apelin signaling",
  "Apoptosis",
  "ATP-chromatin remodel.",
  "Base excision repair",
  "C-type lectin receptor signaling",
  "Cellular senescence",
  "cGMP-PKG signaling",
  "Chemokine signaling",
  "DNA replication",
  "ECM-receptor interaction",
  "ErbB signaling",
  "Fanconi anemia",
  "FoxO signaling",
  "HIF-1 signaling",
  "Hippo signaling",
  "Homologous recombination",
  "IL-17 signaling",
  "JAK-STAT signaling",
  "MAPK signaling",
  "Mismatch repair",
  "mTOR signaling",
  "NF-kappa B signaling",
  "Notch signaling",
  "Nucleotide excision repair",
  "p53 signaling",
  "PI3K-Akt signaling",
  "PPAR signaling",
  "Rap1 signaling",
  "Relaxin signaling",
  "RIG-I-like receptor signaling",
  "TGF-beta signaling",
  "TNF signaling",
  "Wnt signaling",
  "Cell cycle",
  "Signt_1",
  "Signt_2",
  "Signt_3"
)

heatmap_matrix_3.33uM <- make_heatmap_matrix(
  gsea_3.33uM,
  ordered_pathways_3.33uM
)

Heatmap_cell_cycle_3.33uM <- pheatmap::pheatmap(
  heatmap_matrix_3.33uM,
  color = colorRampPalette(c("blue", "white"))(50),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_format = "%.2g",
  fontsize_row = 11,
  fontsize_col = 2,
  fontsize_number = 0,
  border_color = "black",
  filename = file.path(
    output_dir,
    "Fig_4A_GSEA_heatmap_3.33uM.pdf"
  ),
  width = 10,
  height = 8,
  silent = TRUE
)

ordered_pathways_10uM <- c(
  setdiff(
    unique(gsea_10uM$Description),
    c("Signt_1", "Signt_2", "Signt_3")
  ),
  "Signt_1",
  "Signt_2",
  "Signt_3"
)

heatmap_matrix_10uM <- make_heatmap_matrix(
  gsea_10uM,
  ordered_pathways_10uM
)

Heatmap_cell_cycle_10uM <- pheatmap::pheatmap(
  heatmap_matrix_10uM,
  color = colorRampPalette(c("blue", "white"))(50),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_format = "%.2g",
  fontsize_row = 11,
  fontsize_col = 2,
  fontsize_number = 0,
  border_color = "black",
  labels_col = rep("", ncol(heatmap_matrix_10uM)),
  filename = file.path(
    output_dir,
    "Fig_4D_GSEA_heatmap_10uM.pdf"
  ),
  width = 12,
  height = 14,
  silent = TRUE
)

add_mfe_class <- function(data) {
  data %>%
    dplyr::mutate(
      MFE_class_plot = dplyr::if_else(
        Log10_MFE_change > 0 &
          wilcox_test_all_UP_vs_DOWN_pValue < 0.05,
        "Significant positive MFE shift",
        "Other sample"
      )
    )
}

summarise_mfe_shift <- function(data, dose_label) {
  n_total <- nrow(data)
  n_positive <- sum(data$Log10_MFE_change > 0, na.rm = TRUE)
  n_positive_significant <- sum(
    data$Log10_MFE_change > 0 &
      data$wilcox_test_all_UP_vs_DOWN_pValue < 0.05,
    na.rm = TRUE
  )
  n_negative <- sum(data$Log10_MFE_change < 0, na.rm = TRUE)
  n_negative_significant <- sum(
    data$Log10_MFE_change < 0 &
      data$wilcox_test_all_UP_vs_DOWN_pValue < 0.05,
    na.rm = TRUE
  )

  data.frame(
    dose = dose_label,
    n_total = n_total,
    n_positive = n_positive,
    percent_positive = 100 * n_positive / n_total,
    n_positive_significant = n_positive_significant,
    percent_positive_significant_total =
      100 * n_positive_significant / n_total,
    n_negative = n_negative,
    percent_negative = 100 * n_negative / n_total,
    n_negative_significant = n_negative_significant,
    percent_negative_significant_total =
      100 * n_negative_significant / n_total,
    stringsAsFactors = FALSE
  )
}

run_one_sample_tests <- function(data, dose_label) {
  values <- data$Log10_MFE_change
  values <- values[is.finite(values)]

  if (length(values) < 2L || stats::sd(values) == 0) {
    stop("Insufficient variation at ", dose_label, call. = FALSE)
  }

  data.frame(
    dose = dose_label,
    n = length(values),
    mean_Log10_MFE_change = mean(values),
    median_Log10_MFE_change = median(values),
    one_sample_t_greater_p_value =
      t.test(values, mu = 0, alternative = "greater")$p.value,
    one_sample_wilcox_greater_p_value =
      wilcox.test(
        values,
        mu = 0,
        alternative = "greater",
        exact = FALSE
      )$p.value,
    one_sample_t_two_sided_p_value =
      t.test(values, mu = 0, alternative = "two.sided")$p.value,
    one_sample_wilcox_two_sided_p_value =
      wilcox.test(
        values,
        mu = 0,
        alternative = "two.sided",
        exact = FALSE
      )$p.value,
    stringsAsFactors = FALSE
  )
}

make_volcano_panel <- function(
  data,
  dose_label,
  highlighted_conditions,
  output_filename,
  label_separator
) {
  data <- add_mfe_class(data)
  data$condition_key <- paste(data$cell_line, data$Drug, sep = "_")
  data$plot_label <- ifelse(
    data$Log10_MFE_change > 0 &
      data$wilcox_test_all_UP_vs_DOWN_pValue < 0.05 &
      data$condition_key %in% highlighted_conditions,
    paste(data$cell_line, data$Drug, sep = label_separator),
    ""
  )

  volcano_plot <- ggplot(
    data,
    aes(
      x = Log10_MFE_change,
      y = -log10(wilcox_test_all_UP_vs_DOWN_pValue)
    )
  ) +
    theme_classic() +
    geom_point(
      aes(
        fill = MFE_class_plot,
        alpha = MFE_class_plot
      ),
      size = 4,
      shape = 21,
      color = "black",
      stroke = 0.6
    ) +
    scale_fill_manual(
      values = c(
        "Significant positive MFE shift" = "blue",
        "Other sample" = "lightgrey"
      )
    ) +
    scale_alpha_manual(
      values = c(
        "Significant positive MFE shift" = 1,
        "Other sample" = 0.55
      )
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "black"
    ) +
    geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed",
      color = "black"
    ) +
    ggrepel::geom_text_repel(
      aes(label = plot_label),
      size = ifelse(dose_label == "10 µM", 1.5, 3),
      box.padding = 1,
      point.padding = 0.8,
      force = 6,
      max.overlaps = Inf,
      min.segment.length = 0,
      segment.color = "grey50",
      segment.size = 0.5,
      seed = 42
    ) +
    coord_cartesian(
      xlim = c(-0.06, 0.06),
      ylim = c(0, 8)
    ) +
    labs(
      x = expression(Log[10]~normalized~Delta*G[fold]~"change"),
      y = expression(-Log[10]~p.value),
      title = ""
    ) +
    theme(
      legend.position = "none",
      plot.title = element_text(
        color = "black",
        size = 14,
        face = "bold"
      )
    )

  boxplot_horizontal <- ggplot(
    data,
    aes(y = "", x = Log10_MFE_change)
  ) +
    theme_classic() +
    geom_boxplot(
      fill = "lightgrey",
      color = "black",
      width = 0.8,
      outlier.shape = NA
    ) +
    geom_point(
      aes(fill = MFE_class_plot),
      position = position_jitter(
        width = 0,
        height = 0.18,
        seed = 42
      ),
      shape = 21,
      size = 2,
      stroke = 0.4,
      color = "black",
      alpha = 0.8
    ) +
    scale_fill_manual(
      values = c(
        "Significant positive MFE shift" = "blue",
        "Other sample" = "lightgrey"
      )
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.8,
      color = "black"
    ) +
    coord_cartesian(xlim = c(-0.06, 0.06)) +
    labs(y = "", x = "", title = "") +
    theme(
      legend.position = "none",
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank()
    )

  combined_plot <- boxplot_horizontal / volcano_plot +
    patchwork::plot_layout(heights = c(0.07, 0.93))

  ggsave(
    filename = output_filename,
    plot = combined_plot,
    path = output_dir,
    width = 5,
    height = 8.5,
    dpi = 300,
    device = "tiff",
    compression = "lzw"
  )

  list(
    plot = combined_plot,
    summary = summarise_mfe_shift(data, dose_label),
    tests = run_one_sample_tests(data, dose_label)
  )
}

Fig4B <- make_volcano_panel(
  analysis_3.33uM,
  "3.33 µM",
  c(
    "HCC515_neratinib",
    "HT29_trametinib",
    "PC3_palbociclib"
  ),
  "Fig_4B_iLINCS_analysis_3.33uM.tiff",
  "_"
)

Fig4E <- make_volcano_panel(
  analysis_10uM,
  "10 µM",
  c(
    "ES2_lapatinib",
    "PC3_torin-1",
    "A375_tozasertib"
  ),
  "Fig_4E_iLINCS_analysis_10uM.tiff",
  " "
)

print(Fig4B$summary)
print(Fig4B$tests)
print(Fig4E$summary)
print(Fig4E$tests)

make_violin_panel <- function(
  data,
  target_conditions,
  dose_label,
  output_filename
) {
  data$condition_key <- sub(
    "^([^_]+_[^_]+).*",
    "\\1",
    data$SampleName
  )

  plots <- vector("list", length(target_conditions))
  result_tables <- vector("list", length(target_conditions))

  for (i in seq_along(target_conditions)) {
    condition <- target_conditions[[i]]
    sample_names <- unique(data$SampleName[data$condition_key == condition])

    if (length(sample_names) != 1L) {
      stop(
        dose_label,
        " violin data must contain exactly one signature for ",
        condition,
        call. = FALSE
      )
    }

    sample_name <- sample_names[[1]]
    sample_data <- data[data$SampleName == sample_name, , drop = FALSE]
    mfe_up <- sample_data$Normalized_MFE[sample_data$z_score > 0]
    mfe_down <- sample_data$Normalized_MFE[sample_data$z_score < 0]

    if (length(mfe_up) < 2L || length(mfe_down) < 2L) {
      stop("Insufficient violin data for ", sample_name, call. = FALSE)
    }

    t_test <- t.test(mfe_down, mfe_up, paired = FALSE)
    wilcox_test <- wilcox.test(
      mfe_down,
      mfe_up,
      alternative = "two.sided",
      paired = FALSE,
      exact = FALSE
    )

    net_mfe_change <- sum(
      (sample_data$Normalized_MFE * sample_data$z_score) * (-1),
      na.rm = TRUE
    )

    title_text <- if (dose_label == "10 µM") condition else ""

    plots[[i]] <- ggplot(
      sample_data,
      aes(
        x = factor(z_score > 0, levels = c(FALSE, TRUE)),
        y = Normalized_MFE
      )
    ) +
      geom_violin(fill = "lightgrey", alpha = 0.7) +
      geom_boxplot(
        width = 0.15,
        fill = "white",
        color = "black",
        outlier.shape = NA
      ) +
      theme_classic() +
      scale_x_discrete(
        labels = c("z_score < 0", "z_score > 0")
      ) +
      labs(
        x = NULL,
        y = ifelse(dose_label == "3.33 µM", "Normalized MFE", ""),
        title = title_text
      ) +
      theme(
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 20),
        axis.title.y = if (dose_label == "3.33 µM") {
          element_blank()
        } else {
          element_text(size = 20)
        },
        plot.title = if (dose_label == "3.33 µM") {
          element_blank()
        } else {
          element_text(size = 10)
        }
      )

    result_tables[[i]] <- data.frame(
      SampleName = sample_name,
      Net_MFE_change = net_mfe_change,
      MeanMFE_UP = mean(mfe_up, na.rm = TRUE),
      MeanMFE_DOWN = mean(mfe_down, na.rm = TRUE),
      TTEST_pValue = t_test$p.value,
      wilcox_test_all_UP_vs_DOWN_pValue = wilcox_test$p.value,
      stringsAsFactors = FALSE
    )
  }

  combined_plot <-
    plots[[1]] /
    patchwork::plot_spacer() /
    plots[[2]] /
    patchwork::plot_spacer() /
    plots[[3]] +
    patchwork::plot_layout(
      ncol = 1,
      heights = c(1, 0.1, 1, 0.1, 1)
    ) &
    theme(plot.margin = margin(20, 10, 20, 10))

  ggsave(
    filename = output_filename,
    plot = combined_plot,
    path = output_dir,
    width = 5,
    height = ifelse(dose_label == "3.33 µM", 16, 13)
  )

  list(
    plot = combined_plot,
    results = dplyr::bind_rows(result_tables)
  )
}

Fig4C <- make_violin_panel(
  violin_3.33uM,
  c(
    "HCC515_neratinib",
    "HT29_trametinib",
    "PC3_palbociclib"
  ),
  "3.33 µM",
  "Fig_4C_combined_violin_plots_3.33uM.pdf"
)

Fig4F <- make_violin_panel(
  violin_10uM,
  c(
    "ES2_lapatinib",
    "PC3_torin-1",
    "A375_tozasertib"
  ),
  "10 µM",
  "Fig_4F_combined_violin_plots_10uM.pdf"
)

print(Fig4C$results)
print(Fig4F$results)

zweighted <- zweighted %>%
  dplyr::mutate(
    dose = factor(dose, levels = c("3.33 µM", "10 µM")),
    condition_label = paste(cell_line, Drug, sep = "_"),
    permutation_significance = dplyr::case_when(
      is.finite(permutation_p_delta_UP_minus_DOWN) &
        permutation_p_delta_UP_minus_DOWN < 0.05 &
        delta_z_weighted_Normalized_MFE_UP_minus_DOWN > 0 ~
        "Significant positive",
      is.finite(permutation_p_delta_UP_minus_DOWN) &
        permutation_p_delta_UP_minus_DOWN < 0.05 &
        delta_z_weighted_Normalized_MFE_UP_minus_DOWN < 0 ~
        "Significant negative",
      TRUE ~ "Not significant"
    ),
    permutation_significance = factor(
      permutation_significance,
      levels = c(
        "Significant positive",
        "Significant negative",
        "Not significant"
      )
    )
  ) %>%
  dplyr::group_by(dose, condition_label) %>%
  dplyr::mutate(
    duplicated_condition_n = dplyr::n(),
    duplicated_condition_index = dplyr::row_number(),
    condition_label_unique = dplyr::if_else(
      duplicated_condition_n > 1L,
      paste0(condition_label, "_", duplicated_condition_index),
      condition_label
    )
  ) %>%
  dplyr::ungroup()

axis_limit_zweighted <- max(
  abs(zweighted$delta_z_weighted_Normalized_MFE_UP_minus_DOWN),
  na.rm = TRUE
) * 1.10

if (!is.finite(axis_limit_zweighted) || axis_limit_zweighted == 0) {
  axis_limit_zweighted <- 0.01
}

make_zweighted_dotplot <- function(data, dose_label, output_filename) {
  plot_data <- data %>%
    dplyr::filter(dose == dose_label) %>%
    dplyr::arrange(delta_z_weighted_Normalized_MFE_UP_minus_DOWN) %>%
    dplyr::mutate(
      condition_label_unique = factor(
        condition_label_unique,
        levels = condition_label_unique
      )
    )

  if (nrow(plot_data) == 0L) {
    stop("No z-score-weighted rows for ", dose_label, call. = FALSE)
  }

  plot_object <- ggplot(
    plot_data,
    aes(
      x = condition_label_unique,
      y = delta_z_weighted_Normalized_MFE_UP_minus_DOWN
    )
  ) +
    theme_classic() +
    geom_point(
      aes(fill = permutation_significance),
      shape = 21,
      color = "black",
      size = ifelse(dose_label == "3.33 µM", 3.2, 2.7),
      stroke = 0.4,
      alpha = 0.9
    ) +
    scale_fill_manual(
      values = c(
        "Significant positive" = "blue",
        "Significant negative" = "grey40",
        "Not significant" = "white"
      ),
      drop = FALSE
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = 0.6,
      color = "black"
    ) +
    coord_cartesian(
      ylim = c(-axis_limit_zweighted, axis_limit_zweighted)
    ) +
    labs(
      x = "",
      y = expression(
        Delta~"z-score-weighted normalized "*Delta*G[fold]
      ),
      fill = "",
      title = dose_label
    ) +
    theme(
      plot.title = element_text(
        size = 13,
        face = "bold",
        hjust = 0.5
      ),
      axis.text.x = element_text(
        size = ifelse(dose_label == "3.33 µM", 8, 5.5),
        color = "black",
        angle = 45,
        hjust = 1
      ),
      axis.text.y = element_text(size = 10, color = "black"),
      axis.title.y = element_text(size = 12, color = "black"),
      axis.ticks.x = element_blank(),
      legend.position = "right"
    )

  ggsave(
    filename = output_filename,
    plot = plot_object,
    path = output_dir,
    width = ifelse(dose_label == "3.33 µM", 8, 12),
    height = 5
  )

  plot_object
}

Fig4G <- make_zweighted_dotplot(
  zweighted,
  "3.33 µM",
  "Fig_4G_iLINCS_zweighted_DGfold_3.33uM.pdf"
)

Fig4H <- make_zweighted_dotplot(
  zweighted,
  "10 µM",
  "Fig_4H_iLINCS_zweighted_DGfold_10uM.pdf"
)

safe_t_greater <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L || sd(x) == 0) return(NA_real_)
  t.test(x, mu = 0, alternative = "greater")$p.value
}

safe_wilcox_greater <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(NA_real_)
  wilcox.test(
    x,
    mu = 0,
    alternative = "greater",
    exact = FALSE
  )$p.value
}

safe_binomial_greater <- function(x) {
  x <- x[is.finite(x) & x != 0]
  if (length(x) == 0L) return(NA_real_)
  binom.test(
    sum(x > 0),
    length(x),
    p = 0.5,
    alternative = "greater"
  )$p.value
}

zweighted_summary <- zweighted %>%
  dplyr::group_by(dose) %>%
  dplyr::summarise(
    n_conditions = dplyr::n(),
    n_nonzero_conditions = sum(
      delta_z_weighted_Normalized_MFE_UP_minus_DOWN != 0,
      na.rm = TRUE
    ),
    mean_delta_z_weighted_DGfold = mean(
      delta_z_weighted_Normalized_MFE_UP_minus_DOWN,
      na.rm = TRUE
    ),
    median_delta_z_weighted_DGfold = median(
      delta_z_weighted_Normalized_MFE_UP_minus_DOWN,
      na.rm = TRUE
    ),
    n_positive = sum(
      delta_z_weighted_Normalized_MFE_UP_minus_DOWN > 0,
      na.rm = TRUE
    ),
    n_negative = sum(
      delta_z_weighted_Normalized_MFE_UP_minus_DOWN < 0,
      na.rm = TRUE
    ),
    n_significant_positive = sum(
      permutation_significance == "Significant positive",
      na.rm = TRUE
    ),
    n_significant_negative = sum(
      permutation_significance == "Significant negative",
      na.rm = TRUE
    ),
    one_sample_t_test_greater_p_value = safe_t_greater(
      delta_z_weighted_Normalized_MFE_UP_minus_DOWN
    ),
    one_sample_wilcox_greater_p_value = safe_wilcox_greater(
      delta_z_weighted_Normalized_MFE_UP_minus_DOWN
    ),
    binomial_positive_greater_than_half_p_value = safe_binomial_greater(
      delta_z_weighted_Normalized_MFE_UP_minus_DOWN
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    percent_positive = 100 * n_positive / n_nonzero_conditions
  )

y_limit_zweighted <- max(
  abs(zweighted$delta_z_weighted_Normalized_MFE_UP_minus_DOWN),
  na.rm = TRUE
) * 1.35

if (!is.finite(y_limit_zweighted) || y_limit_zweighted == 0) {
  y_limit_zweighted <- 0.01
}

format_exact_p_value <- function(p, digits = 6) {
  ifelse(
    is.na(p),
    "NA",
    format(
      p,
      scientific = TRUE,
      digits = digits,
      trim = TRUE
    )
  )
}

zweighted_summary <- zweighted_summary %>%
  dplyr::mutate(
    dose = factor(dose, levels = c("3.33 µM", "10 µM")),
    plot_label = paste0(
      "Wilcoxon p = ",
      format_exact_p_value(
        one_sample_wilcox_greater_p_value,
        digits = 6
      ),
      "\npositive = ",
      n_positive,
      "/",
      n_nonzero_conditions
    ),
    label_y = y_limit_zweighted * 0.92
  )

GREIN_orange <- "#D55E00"

Fig4I <- ggplot(
  zweighted,
  aes(
    x = dose,
    y = delta_z_weighted_Normalized_MFE_UP_minus_DOWN
  )
) +
  theme_classic() +
  geom_boxplot(
    fill = GREIN_orange,
    color = "black",
    width = 0.55,
    outlier.shape = NA,
    alpha = 0.75
  ) +
  geom_point(
    aes(fill = permutation_significance),
    position = position_jitter(
      width = 0.14,
      height = 0,
      seed = 42
    ),
    size = 2.6,
    shape = 21,
    color = "black",
    stroke = 0.35,
    alpha = 0.85
  ) +
  scale_fill_manual(
    values = c(
      "Significant positive" = "blue",
      "Significant negative" = "grey90",
      "Not significant" = "grey90"
    ),
    drop = FALSE
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.6,
    color = "black"
  ) +
  geom_text(
    data = zweighted_summary,
    aes(
      x = dose,
      y = label_y,
      label = plot_label
    ),
    inherit.aes = FALSE,
    size = 2,
    color = "black"
  ) +
  coord_cartesian(
    ylim = c(-y_limit_zweighted, y_limit_zweighted),
    clip = "off"
  ) +
  labs(
    x = "",
    y = expression(
      Delta~"z-score-weighted normalized "*Delta*G[fold]
    ),
    fill = "",
    title = ""
  ) +
  theme(
    axis.text.x = element_text(size = 11, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 12, color = "black"),
    axis.ticks.x = element_blank(),
    legend.position = "right",
    plot.margin = margin(10, 10, 10, 10)
  )

ggsave(
  filename = "Fig_4I_iLINCS_zweighted_DGfold_summary_boxplot_3.33uM_10uM.tiff",
  plot = Fig4I,
  path = output_dir,
  width = 5,
  height = 6,
  units = "in",
  dpi = 300,
  compression = "lzw"
)

print(zweighted_summary)

cat("\nFigure 4 analysis complete.\n")
cat("Input RDS: ", figure4_data_file, "\n", sep = "")
cat("Output directory: ", output_dir, "\n", sep = "")
