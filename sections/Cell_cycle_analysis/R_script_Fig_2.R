suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(scales)
})

figure2_results_file <- "../../data/cell_cycle_analysis_results_for_Fig_2.rds"
output_dir <- "Figure_2_outputs"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

stop_if_missing <- function(file_path, label) {
  if (!file.exists(file_path)) {
    stop(label, " not found: ", file_path, call. = FALSE)
  }
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

figure2_required_objects <- c(
  "GSE94479_combined_kk2_df_filtered_plot",
  "GSE104736_combined_kk2_df_filtered_plot",
  "GSE94479_G0G1vsG2M",
  "GSE94479_SvsG0G1",
  "MFE_GSE104736_G1_vs_G2",
  "MFE_GSE104736_S_vs_G1",
  "Cell_cycle_MCF7_HeLa_DGfold_quartile_DE_fraction",
  "GSE94479_MCF7_TPM_long",
  "GSE94479_MCF7_TPM_weight_validation",
  "GSE104736_HeLa_TPM_long",
  "GSE104736_HeLa_TPM_weight_validation"
)

stop_if_missing(figure2_results_file, "Figure 2 analysis input")
figure2_objects <- read_named_rds_objects(
  figure2_results_file,
  figure2_required_objects
)

GSE94479_combined_kk2_df_filtered_plot <-
  figure2_objects[["GSE94479_combined_kk2_df_filtered_plot"]]
GSE104736_combined_kk2_df_filtered_plot <-
  figure2_objects[["GSE104736_combined_kk2_df_filtered_plot"]]
GSE94479_G0G1vsG2M <- figure2_objects[["GSE94479_G0G1vsG2M"]]
GSE94479_SvsG0G1 <- figure2_objects[["GSE94479_SvsG0G1"]]
MFE_GSE104736_G1_vs_G2 <-
  figure2_objects[["MFE_GSE104736_G1_vs_G2"]]
MFE_GSE104736_S_vs_G1 <-
  figure2_objects[["MFE_GSE104736_S_vs_G1"]]
Cell_cycle_MCF7_HeLa_DGfold_quartile_DE_fraction <-
  figure2_objects[["Cell_cycle_MCF7_HeLa_DGfold_quartile_DE_fraction"]]
GSE94479_MCF7_TPM_long <-
  figure2_objects[["GSE94479_MCF7_TPM_long"]]
GSE94479_MCF7_TPM_weight_validation <-
  figure2_objects[["GSE94479_MCF7_TPM_weight_validation"]]
GSE104736_HeLa_TPM_long <-
  figure2_objects[["GSE104736_HeLa_TPM_long"]]
GSE104736_HeLa_TPM_weight_validation <-
  figure2_objects[["GSE104736_HeLa_TPM_weight_validation"]]

GSE94479_combined_kk2_df_filtered_plot_edit <- GSE94479_combined_kk2_df_filtered_plot %>%
  dplyr::mutate(cell_line = "MCF7") %>%
  dplyr::filter(!Sample %in% c("GSE94479_G2MvsS"))

GSE94479_combined_kk2_df_filtered_plot_edit <- GSE94479_combined_kk2_df_filtered_plot_edit %>%
  bind_rows(
    tibble(
      ID = c("hsa03410", "hsa01230", "hsa04370"),
      Description = c(
        "Base excision repair",
        "Biosynthesis of amino acids",
        "VEGF signaling pathway"
      ),
      cell_line = "MCF7",
      Sample    = "GSE94479_SvsG0G1",
      !!!setNames(
        rep(list(NA),
            ncol(GSE94479_combined_kk2_df_filtered_plot_edit) - 4),
        setdiff(
          names(GSE94479_combined_kk2_df_filtered_plot_edit),
          c("ID", "Description", "cell_line", "Sample")
        )
      )
    )
  ) %>%
  dplyr::mutate(Sample = gsub("GSE94479_G0G1vsG2M", "GSE94479_G1vsG2M", Sample)) %>%
  dplyr::mutate(Sample = gsub("GSE94479_SvsG0G1", "GSE94479_SvsG1", Sample))

GSE104736_combined_kk2_df_filtered_plot_edit <- GSE104736_combined_kk2_df_filtered_plot %>%
  mutate(cell_line = "HeLa") %>%
  dplyr::filter(!Sample %in% c("GSE104736_G2vsS"))

new_descriptions <- c(
  "Calcium signaling pathway", "Ras signaling pathway", "PI3K-Akt signaling pathway",
  "VEGF signaling pathway", "Rap1 signaling pathway", "MAPK signaling pathway",
  "p53 signaling pathway", "Central carbon metabolism in cancer"
)

new_ids <- c(
  "hsa04020", "hsa04014", "hsa04151", "hsa04370", "hsa04015",
  "hsa04010", "hsa04115", "hsa05230"
)

new_rows <- tibble(
  ID = new_ids,
  Description = new_descriptions,
  cell_line = "HeLa"
) %>%
  mutate(
    !!!setNames(rep(list(NA), ncol(GSE104736_combined_kk2_df_filtered_plot_edit) - 3),
                setdiff(names(GSE104736_combined_kk2_df_filtered_plot_edit), c("ID", "Description", "cell_line")))
  )

GSE104736_combined_kk2_df_filtered_plot_edit <- GSE104736_combined_kk2_df_filtered_plot_edit %>%
  bind_rows(new_rows) %>%
  dplyr::mutate(Sample = if_else(is.na(Sample) & Description == "Ras signaling pathway", "GSE104736_G0G1vsG2M", Sample)) %>%
  dplyr::mutate(Sample = if_else(is.na(Sample) & Description == "Calcium signaling pathway", "GSE104736_G0G1vsG2M", Sample)) %>%
  dplyr::mutate(Sample = if_else(is.na(Sample) & Description == "PI3K-Akt signaling pathway", "GSE104736_G0G1vsG2M", Sample)) %>%
  dplyr::mutate(Sample = if_else(is.na(Sample) & Description == "VEGF signaling pathway", "GSE104736_G0G1vsG2M", Sample)) %>%
  dplyr::mutate(Sample = if_else(is.na(Sample) & Description == "Rap1 signaling pathway", "GSE104736_G0G1vsG2M", Sample)) %>%
  dplyr::mutate(Sample = if_else(is.na(Sample) & Description == "MAPK signaling pathway", "GSE104736_G0G1vsG2M", Sample)) %>%
  dplyr::mutate(Sample = if_else(is.na(Sample) & Description == "p53 signaling pathway", "GSE104736_G0G1vsG2M", Sample)) %>%
  dplyr::mutate(Sample = if_else(is.na(Sample) & Description == "Central carbon metabolism in cancer", "GSE104736_SvsG0G1", Sample)) %>%
  dplyr::mutate(Sample = gsub("GSE104736_G0G1vsG2M", "GSE104736_G1vsG2M", Sample))

Data_cell_cycle_MCF7_HeLa <- rbind(GSE94479_combined_kk2_df_filtered_plot_edit,
                                   GSE104736_combined_kk2_df_filtered_plot_edit) %>%
  mutate(
    Sample = stringr::str_replace(Sample, "GSE94479_", "MCF7_"),
    Sample = stringr::str_replace(Sample, "GSE104736_", "HeLa_"),
    Sample = stringr::str_replace(Sample, "HeLa_SvsG0G1", "HeLa_SvsG1")
  )

desired_order <- c(     "Central carbon metabolism in cancer",
                        "Biosynthesis of amino acids",
                        "Wnt signaling pathway",
                        "Ras signaling pathway",
                        "PI3K-Akt signaling pathway",
                        "MAPK signaling pathway",
                        "Rap1 signaling pathway",
                        "Calcium signaling pathway",
                        "p53 signaling pathway",
                        "DNA replication",
                        "Base excision repair",
                        "Homologous recombination",
                        "Fanconi anemia pathway"
)

Data_cell_cycle_MCF7_HeLa_filtered_plot <- Data_cell_cycle_MCF7_HeLa %>%

  dplyr::filter(Description %in% desired_order) %>%
  dplyr::mutate(
    Description = if_else(
      is.na(Description) & ID == "hsa03410",
      "Base excision repair",
      Description
    ),
    Description = factor(Description, levels = desired_order),
    neg_log10_p_adjust = -log10(p.adjust),
    Sample = factor(
      Sample,
      levels = c("MCF7_G1vsG2M", "MCF7_SvsG1",
                 "HeLa_G1vsG2M", "HeLa_SvsG1")
    )
  )

Plot_cell_proliferation <- ggplot(Data_cell_cycle_MCF7_HeLa_filtered_plot, aes(x = Sample, y = Description, fill = NES, size = neg_log10_p_adjust)) +
  geom_point(shape = 21, color = "black", stroke = 0.5) +
  scale_fill_gradient2(low = "dodgerblue4", mid = "white", high = "firebrick3", midpoint = 0, limits = c(-4, 4)) +
  scale_size_continuous(range = c(2, 10)) +
  theme_bw() +
  labs(title = "",
       x = "",
       y = "",
       fill = "NES",
       size = "-log10(p.adjust)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 10))

print(Plot_cell_proliferation)

ggsave("Plot_cell_proliferation.pdf", plot = Plot_cell_proliferation, width = 6.5, height = 8, path = output_dir)

MCF7_G0G1_vs_G2M_plot <- ggplot(GSE94479_G0G1vsG2M, aes(x = factor(log2FoldChange > 0), y = Normalized_MFE)) +
  geom_violin(fill = "lightgrey", color = "black", size = 0.5) +
  geom_boxplot(width = 0.1, color = "black", fill = "white", outlier.shape = NA) +
  theme_classic() +
  ylim(-0.6, 0) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 20),
    plot.title = element_text(size = 10),
    axis.text.x = element_text(size = 16, color = "black")
  ) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  labs(y = expression(Delta*G[fold]~"/ nt (kcal/mol)"), x = "", title = "") +
  scale_x_discrete(labels = c("log2_FC < 0", "log2_FC > 0"))

print(MCF7_G0G1_vs_G2M_plot)

ggsave("MCF7_G0G1_vs_G2M_plot.pdf", plot = MCF7_G0G1_vs_G2M_plot, width = 5, height = 6, path = output_dir)

MFE_MCF7_G0G1vsG2M_UP <- dplyr::filter(GSE94479_G0G1vsG2M, log2FoldChange > 0)
Mean_MFE_MCF7_G0G1vsG2M_UP <- mean(MFE_MCF7_G0G1vsG2M_UP$Normalized_MFE, na.rm = TRUE)

MFE_MCF7_G0G1vsG2M_DOWN <- dplyr::filter(GSE94479_G0G1vsG2M, log2FoldChange < 0)
Mean_MFE_MCF7_G0G1vsG2M_DOWN <- mean(MFE_MCF7_G0G1vsG2M_DOWN$Normalized_MFE, na.rm = TRUE)

TTEST_MFE_MCF7_G0G1vsG2M_all_UP_vs_DOWN <- t.test(MFE_MCF7_G0G1vsG2M_DOWN$Normalized_MFE, MFE_MCF7_G0G1vsG2M_UP$Normalized_MFE, paired=FALSE)
print(TTEST_MFE_MCF7_G0G1vsG2M_all_UP_vs_DOWN)

ks_test_MCF7_G0G1vsG2M_all_UP_vs_DOWN <- ks.test(MFE_MCF7_G0G1vsG2M_DOWN$Normalized_MFE, MFE_MCF7_G0G1vsG2M_UP$Normalized_MFE)
print(ks_test_MCF7_G0G1vsG2M_all_UP_vs_DOWN)

quantile_075_MFE_MCF7_G0G1vsG2M_UP <- quantile(MFE_MCF7_G0G1vsG2M_UP$Normalized_MFE, probs = c(0.75))
quantile_075_MFE_MCF7_G0G1vsG2M_DOWN <- quantile(MFE_MCF7_G0G1vsG2M_DOWN$Normalized_MFE, probs = c(0.75))

ks_test_MFE_MCF7_G0G1vsG2M_075_above <- ks.test(MFE_MCF7_G0G1vsG2M_UP$Normalized_MFE[MFE_MCF7_G0G1vsG2M_UP$Normalized_MFE >= quantile_075_MFE_MCF7_G0G1vsG2M_UP],
                                                MFE_MCF7_G0G1vsG2M_DOWN$Normalized_MFE[MFE_MCF7_G0G1vsG2M_DOWN$Normalized_MFE >= quantile_075_MFE_MCF7_G0G1vsG2M_DOWN])
print(ks_test_MFE_MCF7_G0G1vsG2M_075_above)

quantile_025_MFE_MCF7_G0G1vsG2M_UP <- quantile(MFE_MCF7_G0G1vsG2M_UP$Normalized_MFE, probs = c(0.25))
quantile_025_MFE_MCF7_G0G1vsG2M_DOWN <- quantile(MFE_MCF7_G0G1vsG2M_DOWN$Normalized_MFE, probs = c(0.25))

ks_test_MF_MCF7_G0G1vsG2M_025_below <- ks.test(MFE_MCF7_G0G1vsG2M_UP$Normalized_MFE[MFE_MCF7_G0G1vsG2M_UP$Normalized_MFE <= quantile_025_MFE_MCF7_G0G1vsG2M_UP],
                                               MFE_MCF7_G0G1vsG2M_DOWN$Normalized_MFE[MFE_MCF7_G0G1vsG2M_DOWN$Normalized_MFE <= quantile_025_MFE_MCF7_G0G1vsG2M_DOWN])
print(ks_test_MF_MCF7_G0G1vsG2M_025_below)

wilcox_test_MCF7_G0G1vsG2M_all_UP_vs_DOWN <- wilcox.test(
  MFE_MCF7_G0G1vsG2M_DOWN$Normalized_MFE,
  MFE_MCF7_G0G1vsG2M_UP$Normalized_MFE,
  paired = FALSE,
  alternative = "two.sided"
)
print(wilcox_test_MCF7_G0G1vsG2M_all_UP_vs_DOWN)

MCF7_SvsG0G1_plot <- ggplot(GSE94479_SvsG0G1, aes(x = factor(log2FoldChange > 0), y = Normalized_MFE)) +
  geom_violin(fill = "lightgrey", color = "black", size = 0.5) +
  geom_boxplot(width = 0.1, color = "black", fill = "white", outlier.shape = NA) +
  theme_classic() +
  ylim(-0.6, 0) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 20),
    plot.title = element_text(size = 10),
    axis.text.x = element_text(size = 16, color = "black")
  ) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  labs(y = expression(Delta*G[fold]~"/ nt (kcal/mol)"), x = "", title = "") +
  scale_x_discrete(labels = c("log2_FC < 0", "log2_FC > 0"))

print(MCF7_SvsG0G1_plot)

ggsave("MCF7_SvsG0G1_plot.pdf", plot = MCF7_SvsG0G1_plot, width = 5, height = 6, path = output_dir)

MFE_MCF7_SvsG0G1_UP <- dplyr::filter(GSE94479_SvsG0G1, log2FoldChange > 0)
Mean_MFE_MCF7_SvsG0G1_UP <- mean(MFE_MCF7_SvsG0G1_UP$Normalized_MFE)

MFE_MCF7_SvsG0G1_DOWN <- dplyr::filter(GSE94479_SvsG0G1, log2FoldChange < 0)
Mean_MFE_MCF7_SvsG0G1_DOWN <- mean(MFE_MCF7_SvsG0G1_DOWN$Normalized_MFE)

TTEST_MFE_MCF7_SvsG0G1_all_UP_vs_DOWN <- t.test(MFE_MCF7_SvsG0G1_DOWN$Normalized_MFE, MFE_MCF7_SvsG0G1_UP$Normalized_MFE, paired=FALSE)
print(TTEST_MFE_MCF7_SvsG0G1_all_UP_vs_DOWN)

ks_test_MCF7_SvsG0G1_all_UP_vs_DOWN <- ks.test(MFE_MCF7_SvsG0G1_DOWN$Normalized_MFE, MFE_MCF7_SvsG0G1_UP$Normalized_MFE)
print(ks_test_MCF7_SvsG0G1_all_UP_vs_DOWN)

quantile_075_MFE_MCF7_SvsG0G1_UP <- quantile(MFE_MCF7_SvsG0G1_UP$Normalized_MFE, probs = c(0.75))
quantile_075_MFE_MCF7_SvsG0G1_DOWN <- quantile(MFE_MCF7_SvsG0G1_DOWN$Normalized_MFE, probs = c(0.75))

ks_test_MFE_MCF7_SvsG0G1_075_above <- ks.test(MFE_MCF7_SvsG0G1_UP$Normalized_MFE[MFE_MCF7_SvsG0G1_UP$Normalized_MFE >= quantile_075_MFE_MCF7_SvsG0G1_UP],
                                              MFE_MCF7_SvsG0G1_DOWN$Normalized_MFE[MFE_MCF7_SvsG0G1_DOWN$Normalized_MFE >= quantile_075_MFE_MCF7_SvsG0G1_DOWN])
print(ks_test_MFE_MCF7_SvsG0G1_075_above)

quantile_025_MFE_MCF7_SvsG0G1_UP <- quantile(MFE_MCF7_SvsG0G1_UP$Normalized_MFE, probs = c(0.25))
quantile_025_MFE_MCF7_SvsG0G1_DOWN <- quantile(MFE_MCF7_SvsG0G1_DOWN$Normalized_MFE, probs = c(0.25))

ks_test_MF_MCF7_SvsG0G1_025_below <- ks.test(MFE_MCF7_SvsG0G1_UP$Normalized_MFE[MFE_MCF7_SvsG0G1_UP$Normalized_MFE <= quantile_025_MFE_MCF7_SvsG0G1_UP],
                                             MFE_MCF7_SvsG0G1_DOWN$Normalized_MFE[MFE_MCF7_SvsG0G1_DOWN$Normalized_MFE <= quantile_025_MFE_MCF7_SvsG0G1_DOWN])
print(ks_test_MF_MCF7_SvsG0G1_025_below)

wilcox_test_MCF7_SvsG0G1_all_UP_vs_DOWN <- wilcox.test(
  MFE_MCF7_SvsG0G1_DOWN$Normalized_MFE,
  MFE_MCF7_SvsG0G1_UP$Normalized_MFE,
  paired = FALSE,
  alternative = "two.sided"
)
print(wilcox_test_MCF7_SvsG0G1_all_UP_vs_DOWN)

GSE104736_HeLa_G1_vs_G2_plot <- ggplot(MFE_GSE104736_G1_vs_G2, aes(x = factor(log2FoldChange > 0), y = Normalized_MFE)) +
  geom_violin(fill = "lightgrey") +
  geom_boxplot(width = 0.1, position = position_dodge(0.9), color = "black", fill = "white", outlier.shape = NA) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 20),
    plot.title = element_text(size = 10),
    axis.text.x = element_text(size = 16, color = "black")
  ) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  ylim(-0.6, 0) +
  labs(y = expression(Delta*G[fold]~"/ nt (kcal/mol)"), x = "", title = "") +
  scale_x_discrete(labels = c("log2FoldChange < 0", "log2FoldChange > 0"))

print(GSE104736_HeLa_G1_vs_G2_plot)

ggsave("GSE104736_HeLa_G1_vs_G2_plot.pdf", plot = GSE104736_HeLa_G1_vs_G2_plot, width = 5, height = 6, path = output_dir)

MFE_GSE104736_G1_vs_G2_UP <- dplyr::filter(MFE_GSE104736_G1_vs_G2, log2FoldChange > 0)
Mean_MFE_GSE104736_G1_vs_G2_UP <- mean(MFE_GSE104736_G1_vs_G2_UP$Normalized_MFE)

MFE_GSE104736_G1_vs_G2_DOWN <- dplyr::filter(MFE_GSE104736_G1_vs_G2, log2FoldChange < 0)
Mean_MFE_GSE104736_G1_vs_G2_DOWN <- mean(MFE_GSE104736_G1_vs_G2_DOWN$Normalized_MFE)

TTEST_MFE_GSE104736_G1_vs_G2_all_UP_vs_DOWN <- t.test(MFE_GSE104736_G1_vs_G2_DOWN$Normalized_MFE, MFE_GSE104736_G1_vs_G2_UP$Normalized_MFE, paired=FALSE)
print(TTEST_MFE_GSE104736_G1_vs_G2_all_UP_vs_DOWN)

ks_test_GSE104736_all_UP_vs_DOWN <- ks.test(MFE_GSE104736_G1_vs_G2_DOWN$Normalized_MFE, MFE_GSE104736_G1_vs_G2_UP$Normalized_MFE)
print(ks_test_GSE104736_all_UP_vs_DOWN)

quantile_075_MFE_GSE104736_G1_vs_G2_UP <- quantile(MFE_GSE104736_G1_vs_G2_UP$Normalized_MFE, probs = c(0.75))
quantile_075_MFE_GSE104736_G1_vs_G2_DOWN <- quantile(MFE_GSE104736_G1_vs_G2_DOWN$Normalized_MFE, probs = c(0.75))

ks_test_MFE_GSE104736_G1_vs_G2_075_above <- ks.test(MFE_GSE104736_G1_vs_G2_UP$Normalized_MFE[MFE_GSE104736_G1_vs_G2_UP$Normalized_MFE >= quantile_075_MFE_GSE104736_G1_vs_G2_UP],
                                                    MFE_GSE104736_G1_vs_G2_DOWN$Normalized_MFE[MFE_GSE104736_G1_vs_G2_DOWN$Normalized_MFE >= quantile_075_MFE_GSE104736_G1_vs_G2_DOWN])
print(ks_test_MFE_GSE104736_G1_vs_G2_075_above)

quantile_025_MFE_GSE104736_G1_vs_G2_UP <- quantile(MFE_GSE104736_G1_vs_G2_UP$Normalized_MFE, probs = c(0.25))
quantile_025_MFE_GSE104736_G1_vs_G2_DOWN <- quantile(MFE_GSE104736_G1_vs_G2_DOWN$Normalized_MFE, probs = c(0.25))

ks_test_MF_GSE104736_025_below <- ks.test(MFE_GSE104736_G1_vs_G2_UP$Normalized_MFE[MFE_GSE104736_G1_vs_G2_UP$Normalized_MFE <= quantile_025_MFE_GSE104736_G1_vs_G2_UP],
                                          MFE_GSE104736_G1_vs_G2_DOWN$Normalized_MFE[MFE_GSE104736_G1_vs_G2_DOWN$Normalized_MFE <= quantile_025_MFE_GSE104736_G1_vs_G2_DOWN])
print(ks_test_MF_GSE104736_025_below)

wilcox_test_GSE104736_G1_vs_G2_all_UP_vs_DOWN <- wilcox.test(
  MFE_GSE104736_G1_vs_G2_DOWN$Normalized_MFE,
  MFE_GSE104736_G1_vs_G2_UP$Normalized_MFE,
  paired = FALSE,
  alternative = "two.sided"
)
print(wilcox_test_GSE104736_G1_vs_G2_all_UP_vs_DOWN)

GSE104736_HeLa_S_vs_G1_plot <- ggplot(MFE_GSE104736_S_vs_G1, aes(x = factor(log2FoldChange > 0), y = Normalized_MFE)) +
  geom_violin(fill = "lightgrey") +
  geom_boxplot(width = 0.1, position = position_dodge(0.9), color = "black", fill = "white", outlier.shape = NA) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 20),
    plot.title = element_text(size = 10),
    axis.text.x = element_text(size = 16, color = "black")
  ) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  ylim(-0.6, 0) +
  labs(y = expression(Delta*G[fold]~"/ nt (kcal/mol)"), x = "", title = "") +
  scale_x_discrete(labels = c("log2FoldChange < 0", "log2FoldChange > 0"))

print(GSE104736_HeLa_S_vs_G1_plot)

ggsave("GSE104736_HeLa_S_vs_G1_plot.pdf", plot = GSE104736_HeLa_S_vs_G1_plot, width = 5, height = 6, path = output_dir)

MFE_GSE104736_S_vs_G1_UP <- dplyr::filter(MFE_GSE104736_S_vs_G1, log2FoldChange > 0)
Mean_MFE_GSE104736_S_vs_G1_UP <- mean(MFE_GSE104736_S_vs_G1_UP$Normalized_MFE)

MFE_GSE104736_S_vs_G1_DOWN <- dplyr::filter(MFE_GSE104736_S_vs_G1, log2FoldChange < 0)
Mean_MFE_GSE104736_S_vs_G1_DOWN <- mean(MFE_GSE104736_S_vs_G1_DOWN$Normalized_MFE)

TTEST_MFE_GSE104736_S_vs_G1_all_UP_vs_DOWN <- t.test(MFE_GSE104736_S_vs_G1_DOWN$Normalized_MFE, MFE_GSE104736_S_vs_G1_UP$Normalized_MFE, paired=FALSE)
print(TTEST_MFE_GSE104736_S_vs_G1_all_UP_vs_DOWN)

ks_test_GSE104736_all_UP_vs_DOWN <- ks.test(MFE_GSE104736_S_vs_G1_DOWN$Normalized_MFE, MFE_GSE104736_S_vs_G1_UP$Normalized_MFE)
print(ks_test_GSE104736_all_UP_vs_DOWN)

quantile_075_MFE_GSE104736_S_vs_G1_UP <- quantile(MFE_GSE104736_S_vs_G1_UP$Normalized_MFE, probs = c(0.75))
quantile_075_MFE_GSE104736_S_vs_G1_DOWN <- quantile(MFE_GSE104736_S_vs_G1_DOWN$Normalized_MFE, probs = c(0.75))

ks_test_MFE_GSE104736_S_vs_G1_075_above <- ks.test(MFE_GSE104736_S_vs_G1_UP$Normalized_MFE[MFE_GSE104736_S_vs_G1_UP$Normalized_MFE >= quantile_075_MFE_GSE104736_S_vs_G1_UP],
                                                   MFE_GSE104736_S_vs_G1_DOWN$Normalized_MFE[MFE_GSE104736_S_vs_G1_DOWN$Normalized_MFE >= quantile_075_MFE_GSE104736_S_vs_G1_DOWN])
print(ks_test_MFE_GSE104736_S_vs_G1_075_above)

quantile_025_MFE_GSE104736_S_vs_G1_UP <- quantile(MFE_GSE104736_S_vs_G1_UP$Normalized_MFE, probs = c(0.25))
quantile_025_MFE_GSE104736_S_vs_G1_DOWN <- quantile(MFE_GSE104736_S_vs_G1_DOWN$Normalized_MFE, probs = c(0.25))

ks_test_MF_GSE104736_025_below <- ks.test(MFE_GSE104736_S_vs_G1_UP$Normalized_MFE[MFE_GSE104736_S_vs_G1_UP$Normalized_MFE <= quantile_025_MFE_GSE104736_S_vs_G1_UP],
                                          MFE_GSE104736_S_vs_G1_DOWN$Normalized_MFE[MFE_GSE104736_S_vs_G1_DOWN$Normalized_MFE <= quantile_025_MFE_GSE104736_S_vs_G1_DOWN])
print(ks_test_MF_GSE104736_025_below)

wilcox_test_GSE104736_S_vs_G1_all_UP_vs_DOWN <- wilcox.test(
  MFE_GSE104736_S_vs_G1_DOWN$Normalized_MFE,
  MFE_GSE104736_S_vs_G1_UP$Normalized_MFE,
  paired = FALSE,
  alternative = "two.sided"
)
print(wilcox_test_GSE104736_S_vs_G1_all_UP_vs_DOWN)

stopifnot(
  exists("Cell_cycle_MCF7_HeLa_DGfold_quartile_DE_fraction")
)

stopifnot(
  all(
    c(
      "Dataset",
      "Transition",
      "Phase",
      "DGfold_quartile",
      "normalized_count_fraction"
    ) %in%
      colnames(Cell_cycle_MCF7_HeLa_DGfold_quartile_DE_fraction)
  )
)

Fig2D_DGfold_quartile_DE_fraction_plot <-
  Cell_cycle_MCF7_HeLa_DGfold_quartile_DE_fraction %>%
  dplyr::mutate(
    Dataset = factor(
      Dataset,
      levels = c("MCF7", "HeLa")
    ),
    Transition = factor(
      Transition,
      levels = c("G1 vs G2/M", "S vs G1")
    ),
    Phase = factor(
      Phase,
      levels = c("G2/M", "G1", "S")
    ),
    DGfold_quartile_plot = dplyr::case_when(
      DGfold_quartile == "Least stable\nupper DGfold quartile" ~
        "Upper DGfold\nquartile",
      DGfold_quartile == "Most stable\nlower DGfold quartile" ~
        "Lower DGfold\nquartile",
      TRUE ~ as.character(DGfold_quartile)
    ),
    DGfold_quartile_plot = factor(
      DGfold_quartile_plot,
      levels = c(
        "Upper DGfold\nquartile",
        "Lower DGfold\nquartile"
      )
    )
  )

Data_Fig2D_MCF7_more_negative_G1_vs_G2M <-
  Fig2D_DGfold_quartile_DE_fraction_plot %>%
  dplyr::filter(
    Dataset == "MCF7",
    DGfold_quartile_plot == "Lower DGfold\nquartile",
    Transition == "G1 vs G2/M"
  ) %>%
  droplevels()

Plot_Fig2D_MCF7_more_negative_G1_vs_G2M <- ggplot(
  Data_Fig2D_MCF7_more_negative_G1_vs_G2M,
  aes(
    x = Phase,
    y = normalized_count_fraction,
    fill = Phase
  )
) +
  geom_col(
    width = 0.85,
    color = "black",
    linewidth = 0.3
  ) +
  geom_text(
    aes(label = round(normalized_count_fraction, 3)),
    vjust = -0.4,
    size = 2,
    color = "black"
  ) +
  labs(

    subtitle = "G1 vs G2/M",
    x = "",
    y = "Fraction of total counts",
    fill = ""
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 0,
      hjust = 1,
      size = 8,
      color = "black"
    ),
    axis.text.y = element_text(
      size = 8,
      color = "black"
    ),
    axis.title.y = element_text(size = 8),
    plot.title = element_text(size = 8, hjust = 0.5, color = "black"),
    plot.subtitle = element_text(size = 8, hjust = 0.5, color = "black"),
    legend.position = "none",
    plot.margin = margin(5.5, 20, 5.5, 5.5)
  ) +
  coord_cartesian(
    ylim = c(0, 0.45),
    clip = "off"
  ) +
  scale_fill_manual(
    values = c(
      "G2/M" = "grey85",
      "G1" = "#8cc5e3",
      "S" = "#2066a8"
    ),
    drop = TRUE
  )

print(Plot_Fig2D_MCF7_more_negative_G1_vs_G2M)

ggsave(
  filename = "Plot_Fig2D_MCF7_more_negative_G1_vs_G2M.tiff",
  plot = Plot_Fig2D_MCF7_more_negative_G1_vs_G2M,
  path = output_dir,
  device = "tiff",
  width = 1.5,
  height = 5,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

Data_Fig2D_MCF7_more_negative_S_vs_G1 <-
  Fig2D_DGfold_quartile_DE_fraction_plot %>%
  dplyr::filter(
    Dataset == "MCF7",
    DGfold_quartile_plot == "Lower DGfold\nquartile",
    Transition == "S vs G1"
  ) %>%
  droplevels()

Plot_Fig2D_MCF7_more_negative_S_vs_G1 <- ggplot(
  Data_Fig2D_MCF7_more_negative_S_vs_G1,
  aes(
    x = Phase,
    y = normalized_count_fraction,
    fill = Phase
  )
) +
  geom_col(
    width = 0.85,
    color = "black",
    linewidth = 0.3
  ) +
  geom_text(
    aes(label = round(normalized_count_fraction, 3)),
    vjust = -0.4,
    size = 2,
    color = "black"
  ) +
  labs(

    subtitle = "S vs G1",
    x = "",
    y = "Fraction of total counts",
    fill = ""
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 0,
      hjust = 1,
      size = 8,
      color = "black"
    ),
    axis.text.y = element_text(
      size = 8,
      color = "black"
    ),
    axis.title.y = element_text(size = 8),
    plot.title = element_text(size = 8, hjust = 0.5, color = "black"),
    plot.subtitle = element_text(size = 8, hjust = 0.5, color = "black"),
    legend.position = "none",
    plot.margin = margin(5.5, 20, 5.5, 5.5)
  ) +
  coord_cartesian(
    ylim = c(0, 0.45),
    clip = "off"
  ) +
  scale_fill_manual(
    values = c(
      "G2/M" = "grey85",
      "G1" = "#8cc5e3",
      "S" = "#2066a8"
    ),
    drop = TRUE
  )

print(Plot_Fig2D_MCF7_more_negative_S_vs_G1)

ggsave(
  filename = "Plot_Fig2D_MCF7_more_negative_S_vs_G1.tiff",
  plot = Plot_Fig2D_MCF7_more_negative_S_vs_G1,
  path = output_dir,
  device = "tiff",
  width = 1.5,
  height = 5,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

Data_Fig2D_MCF7_less_negative_G1_vs_G2M <-
  Fig2D_DGfold_quartile_DE_fraction_plot %>%
  dplyr::filter(
    Dataset == "MCF7",
    DGfold_quartile_plot == "Upper DGfold\nquartile",
    Transition == "G1 vs G2/M"
  ) %>%
  droplevels()

Plot_Fig2D_MCF7_less_negative_G1_vs_G2M <- ggplot(
  Data_Fig2D_MCF7_less_negative_G1_vs_G2M,
  aes(
    x = Phase,
    y = normalized_count_fraction,
    fill = Phase
  )
) +
  geom_col(width = 0.85, color = "black", linewidth = 0.3) +
  geom_text(
    aes(label = round(normalized_count_fraction, 3)),
    vjust = -0.4,
    size = 2,
    color = "black"
  ) +
  labs(

    subtitle = "G1 vs G2/M",
    x = "",
    y = "Fraction of total counts",
    fill = ""
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1, size = 8, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.title.y = element_text(size = 8),
    plot.title = element_text(size = 8, hjust = 0.5, color = "black"),
    plot.subtitle = element_text(size = 8, hjust = 0.5, color = "black"),
    legend.position = "none",
    plot.margin = margin(5.5, 20, 5.5, 5.5)
  ) +
  coord_cartesian(
    ylim = c(0, 0.45),
    clip = "off"
  ) +
  scale_fill_manual(
    values = c("G2/M" = "grey85", "G1" = "#8cc5e3", "S" = "#2066a8"),
    drop = TRUE
  )

print(Plot_Fig2D_MCF7_less_negative_G1_vs_G2M)

ggsave(
  filename = "Plot_Fig2D_MCF7_less_negative_G1_vs_G2M.tiff",
  plot = Plot_Fig2D_MCF7_less_negative_G1_vs_G2M,
  path = output_dir,
  device = "tiff",
  width = 1.5,
  height = 5,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

Data_Fig2D_MCF7_less_negative_S_vs_G1 <-
  Fig2D_DGfold_quartile_DE_fraction_plot %>%
  dplyr::filter(
    Dataset == "MCF7",
    DGfold_quartile_plot == "Upper DGfold\nquartile",
    Transition == "S vs G1"
  ) %>%
  droplevels()

Plot_Fig2D_MCF7_less_negative_S_vs_G1 <- ggplot(
  Data_Fig2D_MCF7_less_negative_S_vs_G1,
  aes(
    x = Phase,
    y = normalized_count_fraction,
    fill = Phase
  )
) +
  geom_col(width = 0.85, color = "black", linewidth = 0.3) +
  geom_text(
    aes(label = round(normalized_count_fraction, 3)),
    vjust = -0.4,
    size = 2,
    color = "black"
  ) +
  labs(

    subtitle = "S vs G1",
    x = "",
    y = "Fraction of total counts",
    fill = ""
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1, size = 8, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.title.y = element_text(size = 8),
    plot.title = element_text(size = 8, hjust = 0.5, color = "black"),
    plot.subtitle = element_text(size = 8, hjust = 0.5, color = "black"),
    legend.position = "none",
    plot.margin = margin(5.5, 20, 5.5, 5.5)
  ) +
  coord_cartesian(
    ylim = c(0, 0.45),
    clip = "off"
  ) +
  scale_fill_manual(
    values = c("G2/M" = "grey85", "G1" = "#8cc5e3", "S" = "#2066a8"),
    drop = TRUE
  )

print(Plot_Fig2D_MCF7_less_negative_S_vs_G1)

ggsave(
  filename = "Plot_Fig2D_MCF7_less_negative_S_vs_G1.tiff",
  plot = Plot_Fig2D_MCF7_less_negative_S_vs_G1,
  path = output_dir,
  device = "tiff",
  width = 1.5,
  height = 5,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

Data_Fig2D_HeLa_more_negative_G1_vs_G2M <-
  Fig2D_DGfold_quartile_DE_fraction_plot %>%
  dplyr::filter(
    Dataset == "HeLa",
    DGfold_quartile_plot == "Lower DGfold\nquartile",
    Transition == "G1 vs G2/M"
  ) %>%
  droplevels()

Plot_Fig2D_HeLa_more_negative_G1_vs_G2M <- ggplot(
  Data_Fig2D_HeLa_more_negative_G1_vs_G2M,
  aes(
    x = Phase,
    y = normalized_count_fraction,
    fill = Phase
  )
) +
  geom_col(width = 0.85, color = "black", linewidth = 0.3) +
  geom_text(
    aes(label = round(normalized_count_fraction, 3)),
    vjust = -0.4,
    size = 2,
    color = "black"
  ) +
  labs(

    subtitle = "G1 vs G2/M",
    x = "",
    y = "Fraction of total counts",
    fill = ""
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1, size = 8, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.title.y = element_text(size = 8),
    plot.title = element_text(size = 8, hjust = 0.5, color = "black"),
    plot.subtitle = element_text(size = 8, hjust = 0.5, color = "black"),
    legend.position = "none",
    plot.margin = margin(5.5, 20, 5.5, 5.5)
  ) +
  coord_cartesian(
    ylim = c(0, 0.45),
    clip = "off"
  ) +
  scale_fill_manual(
    values = c("G2/M" = "grey85", "G1" = "#8cc5e3", "S" = "#2066a8"),
    drop = TRUE
  )

print(Plot_Fig2D_HeLa_more_negative_G1_vs_G2M)

ggsave(
  filename = "Plot_Fig2D_HeLa_more_negative_G1_vs_G2M.tiff",
  plot = Plot_Fig2D_HeLa_more_negative_G1_vs_G2M,
  path = output_dir,
  device = "tiff",
  width = 1.5,
  height = 5,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

Data_Fig2D_HeLa_more_negative_S_vs_G1 <-
  Fig2D_DGfold_quartile_DE_fraction_plot %>%
  dplyr::filter(
    Dataset == "HeLa",
    DGfold_quartile_plot == "Lower DGfold\nquartile",
    Transition == "S vs G1"
  ) %>%
  droplevels()

Plot_Fig2D_HeLa_more_negative_S_vs_G1 <- ggplot(
  Data_Fig2D_HeLa_more_negative_S_vs_G1,
  aes(
    x = Phase,
    y = normalized_count_fraction,
    fill = Phase
  )
) +
  geom_col(width = 0.85, color = "black", linewidth = 0.3) +
  geom_text(
    aes(label = round(normalized_count_fraction, 3)),
    vjust = -0.4,
    size = 2,
    color = "black"
  ) +
  labs(

    subtitle = "S vs G1",
    x = "",
    y = "Fraction of total counts",
    fill = ""
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1, size = 8, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.title.y = element_text(size = 8),
    plot.title = element_text(size = 8, hjust = 0.5, color = "black"),
    plot.subtitle = element_text(size = 8, hjust = 0.5, color = "black"),
    legend.position = "none",
    plot.margin = margin(5.5, 20, 5.5, 5.5)
  ) +
  coord_cartesian(
    ylim = c(0, 0.45),
    clip = "off"
  ) +
  scale_fill_manual(
    values = c("G2/M" = "grey85", "G1" = "#8cc5e3", "S" = "#2066a8"),
    drop = TRUE
  )

print(Plot_Fig2D_HeLa_more_negative_S_vs_G1)

ggsave(
  filename = "Plot_Fig2D_HeLa_more_negative_S_vs_G1.tiff",
  plot = Plot_Fig2D_HeLa_more_negative_S_vs_G1,
  path = output_dir,
  device = "tiff",
  width = 1.5,
  height = 5,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

Data_Fig2D_HeLa_less_negative_G1_vs_G2M <-
  Fig2D_DGfold_quartile_DE_fraction_plot %>%
  dplyr::filter(
    Dataset == "HeLa",
    DGfold_quartile_plot == "Upper DGfold\nquartile",
    Transition == "G1 vs G2/M"
  ) %>%
  droplevels()

Plot_Fig2D_HeLa_less_negative_G1_vs_G2M <- ggplot(
  Data_Fig2D_HeLa_less_negative_G1_vs_G2M,
  aes(
    x = Phase,
    y = normalized_count_fraction,
    fill = Phase
  )
) +
  geom_col(width = 0.85, color = "black", linewidth = 0.3) +
  geom_text(
    aes(label = round(normalized_count_fraction, 3)),
    vjust = -0.4,
    size = 2,
    color = "black"
  ) +
  labs(

    subtitle = "G1 vs G2/M",
    x = "",
    y = "Fraction of total counts",
    fill = ""
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1, size = 8, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.title.y = element_text(size = 8),
    plot.title = element_text(size = 8, hjust = 0.5, color = "black"),
    plot.subtitle = element_text(size = 8, hjust = 0.5, color = "black"),
    legend.position = "none",
    plot.margin = margin(5.5, 20, 5.5, 5.5)
  ) +
  coord_cartesian(
    ylim = c(0, 0.45),
    clip = "off"
  ) +
  scale_fill_manual(
    values = c("G2/M" = "grey85", "G1" = "#8cc5e3", "S" = "#2066a8"),
    drop = TRUE
  )

print(Plot_Fig2D_HeLa_less_negative_G1_vs_G2M)

ggsave(
  filename = "Plot_Fig2D_HeLa_less_negative_G1_vs_G2M.tiff",
  plot = Plot_Fig2D_HeLa_less_negative_G1_vs_G2M,
  path = output_dir,
  device = "tiff",
  width = 1.5,
  height = 5,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

Data_Fig2D_HeLa_less_negative_S_vs_G1 <-
  Fig2D_DGfold_quartile_DE_fraction_plot %>%
  dplyr::filter(
    Dataset == "HeLa",
    DGfold_quartile_plot == "Upper DGfold\nquartile",
    Transition == "S vs G1"
  ) %>%
  droplevels()

Plot_Fig2D_HeLa_less_negative_S_vs_G1 <- ggplot(
  Data_Fig2D_HeLa_less_negative_S_vs_G1,
  aes(
    x = Phase,
    y = normalized_count_fraction,
    fill = Phase
  )
) +
  geom_col(width = 0.85, color = "black", linewidth = 0.3) +
  geom_text(
    aes(label = round(normalized_count_fraction, 3)),
    vjust = -0.4,
    size = 2,
    color = "black"
  ) +
  labs(

    subtitle = "S vs G1",
    x = "",
    y = "Fraction of total counts",
    fill = ""
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1, size = 8, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.title.y = element_text(size = 8),
    plot.title = element_text(size = 8, hjust = 0.5, color = "black"),
    plot.subtitle = element_text(size = 8, hjust = 0.5, color = "black"),
    legend.position = "none",
    plot.margin = margin(5.5, 20, 5.5, 5.5)
  ) +
  coord_cartesian(
    ylim = c(0, 0.45),
    clip = "off"
  )+
  scale_fill_manual(
    values = c("G2/M" = "grey85", "G1" = "#8cc5e3", "S" = "#2066a8"),
    drop = TRUE
  )

print(Plot_Fig2D_HeLa_less_negative_S_vs_G1)

ggsave(
  filename = "Plot_Fig2D_HeLa_less_negative_S_vs_G1.tiff",
  plot = Plot_Fig2D_HeLa_less_negative_S_vs_G1,
  path = output_dir,
  device = "tiff",
  width = 1.5,
  height = 5,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

stopifnot(
  all(
    c(
      "Dataset",
      "Transition",
      "Phase",
      "DGfold_quartile",
      "normalized_count_fraction",
      "normalized_counts_in_quartile",
      "total_normalized_counts_DE_universe",
      "phase_1",
      "phase_2"
    ) %in%
      colnames(Cell_cycle_MCF7_HeLa_DGfold_quartile_DE_fraction)
  )
)

Fig2D_DGfold_quartile_DE_z_tests <-
  Cell_cycle_MCF7_HeLa_DGfold_quartile_DE_fraction %>%
  dplyr::mutate(
    Dataset = as.character(Dataset),
    Transition = as.character(Transition),
    Phase = as.character(Phase),
    DGfold_quartile = as.character(DGfold_quartile)
  ) %>%
  dplyr::group_by(
    Dataset,
    Transition,
    DGfold_quartile
  ) %>%
  dplyr::group_modify(~ {

    phase_1 <- unique(.x$phase_1)
    phase_2 <- unique(.x$phase_2)

    stopifnot(length(phase_1) == 1)
    stopifnot(length(phase_2) == 1)

    row_1 <- .x %>%
      dplyr::filter(Phase == phase_1)

    row_2 <- .x %>%
      dplyr::filter(Phase == phase_2)

    stopifnot(nrow(row_1) == 1)
    stopifnot(nrow(row_2) == 1)

    x1 <- row_1$normalized_counts_in_quartile
    n1 <- row_1$total_normalized_counts_DE_universe

    x2 <- row_2$normalized_counts_in_quartile
    n2 <- row_2$total_normalized_counts_DE_universe

    p1 <- x1 / n1
    p2 <- x2 / n2

    p_pool <- (x1 + x2) / (n1 + n2)

    se <- sqrt(
      p_pool *
        (1 - p_pool) *
        ((1 / n1) + (1 / n2))
    )

    z_value <- (p2 - p1) / se

    p_value <- 2 * (1 - pnorm(abs(z_value)))

    tibble::tibble(
      phase_1 = phase_1,
      phase_2 = phase_2,
      comparison = paste(phase_2, "vs", phase_1),
      fraction_phase_1 = p1,
      fraction_phase_2 = p2,
      delta_phase2_minus_phase1 = p2 - p1,
      counts_phase_1 = x1,
      counts_phase_2 = x2,
      total_counts_phase_1 = n1,
      total_counts_phase_2 = n2,
      z_value = z_value,
      p_value = p_value
    )
  }) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    p_adj_BH = p.adjust(p_value, method = "BH"),
    DGfold_quartile_plot = dplyr::case_when(
      DGfold_quartile == "Least stable\nupper DGfold quartile" ~
        "Upper DGfold quartile / less negative DGfold",
      DGfold_quartile == "Most stable\nlower DGfold quartile" ~
        "Lower DGfold quartile / more negative DGfold",
      TRUE ~ DGfold_quartile
    ),
    p_label = dplyr::case_when(
      is.na(p_value) ~ NA_character_,
      p_value < 0.001 ~ "p < 0.001",
      TRUE ~ paste0("p = ", signif(p_value, 2))
    ),
    p_adj_label = dplyr::case_when(
      is.na(p_adj_BH) ~ NA_character_,
      p_adj_BH < 0.001 ~ "BH-adjusted p < 0.001",
      TRUE ~ paste0("BH-adjusted p = ", signif(p_adj_BH, 2))
    )
  ) %>%
  dplyr::select(
    Dataset,
    Transition,
    DGfold_quartile_plot,
    phase_1,
    phase_2,
    comparison,
    fraction_phase_1,
    fraction_phase_2,
    delta_phase2_minus_phase1,
    counts_phase_1,
    counts_phase_2,
    total_counts_phase_1,
    total_counts_phase_2,
    z_value,
    p_value,
    p_adj_BH,
    p_label,
    p_adj_label
  ) %>%
  dplyr::arrange(
    Dataset,
    Transition,
    DGfold_quartile_plot
  )

print(Fig2D_DGfold_quartile_DE_z_tests)

Fig2E_B_perm <- 5000L
Fig2E_requested_cores <- 8L
Fig2E_output_dir <- output_dir

Fig2E_available_cores <- parallel::detectCores(
  logical = TRUE
)

if (is.na(Fig2E_available_cores)) {
  Fig2E_available_cores <- 1L
}

Fig2E_n_cores <- min(
  Fig2E_requested_cores,
  max(1L, as.integer(Fig2E_available_cores))
)

if (Fig2E_n_cores < Fig2E_requested_cores) {
  message(
    "Using ", Fig2E_n_cores,
    " permutation worker(s); ", Fig2E_requested_cores,
    " were requested."
  )
}

Fig2E_required_objects <- c(
  "GSE94479_MCF7_TPM_long",
  "GSE94479_MCF7_TPM_weight_validation",
  "GSE104736_HeLa_TPM_long",
  "GSE104736_HeLa_TPM_weight_validation"
)

Fig2E_missing_objects <- setdiff(
  Fig2E_required_objects,
  names(figure2_objects)
)

if (length(Fig2E_missing_objects) > 0) {
  stop(
    paste0(
      "Run the completed cell-cycle master script first. ",
      "The following TPM objects are missing: ",
      paste(
        Fig2E_missing_objects,
        collapse = ", "
      )
    )
  )
}

Fig2E_standardize_phase <- function(phase_value) {

  phase_value <- as.character(phase_value)

  phase_value[
    phase_value == "G2"
  ] <- "G2/M"

  phase_value
}

Fig2E_analysis_config <- tibble::tribble(
  ~Analysis_ID,
  ~Dataset,
  ~Comparison,
  ~Phase_1,
  ~Phase_2,
  ~TPM_long_object,
  ~TPM_sample_object,
  ~Seed,

  "MCF7_G1_minus_G2M",
  "MCF7",
  "G1 - G2/M",
  "G1",
  "G2/M",
  "GSE94479_MCF7_TPM_long",
  "GSE94479_MCF7_TPM_weight_validation",
  201L,

  "MCF7_S_minus_G1",
  "MCF7",
  "S - G1",
  "S",
  "G1",
  "GSE94479_MCF7_TPM_long",
  "GSE94479_MCF7_TPM_weight_validation",
  202L,

  "HeLa_G1_minus_G2M",
  "HeLa",
  "G1 - G2/M",
  "G1",
  "G2/M",
  "GSE104736_HeLa_TPM_long",
  "GSE104736_HeLa_TPM_weight_validation",
  203L,

  "HeLa_S_minus_G1",
  "HeLa",
  "S - G1",
  "S",
  "G1",
  "GSE104736_HeLa_TPM_long",
  "GSE104736_HeLa_TPM_weight_validation",
  204L
)

Fig2E_run_one_TPM_contrast <- function(
    TPM_long,
    TPM_sample,
    analysis_id,
    dataset_label,
    comparison_label,
    phase_1,
    phase_2,
    B_perm = 5000L,
    n_cores = 8L,
    seed = 123L
) {

  required_long_columns <- c(
    "ensembl_transcript_id",
    "Sample",
    "Phase",
    "TPM",
    "Normalized_MFE",
    "relative_transcript_abundance_MFE"
  )

  required_sample_columns <- c(
    "Sample",
    "Phase",
    "annotated_transcripts",
    "expressed_annotated_transcripts",
    "retained_TPM_MFE",
    "sum_MFE_weights",
    "TPM_weighted_Normalized_MFE"
  )

  missing_long_columns <- setdiff(
    required_long_columns,
    colnames(TPM_long)
  )

  missing_sample_columns <- setdiff(
    required_sample_columns,
    colnames(TPM_sample)
  )

  if (length(missing_long_columns) > 0) {
    stop(
      analysis_id,
      ": missing TPM-long columns: ",
      paste(
        missing_long_columns,
        collapse = ", "
      )
    )
  }

  if (length(missing_sample_columns) > 0) {
    stop(
      analysis_id,
      ": missing TPM-summary columns: ",
      paste(
        missing_sample_columns,
        collapse = ", "
      )
    )
  }

  sample_values <- TPM_sample %>%
    dplyr::transmute(
      Analysis_ID = analysis_id,
      Dataset = dataset_label,
      Sample = as.character(Sample),
      Phase = Fig2E_standardize_phase(Phase),

      annotated_transcripts =
        as.integer(annotated_transcripts),

      expressed_annotated_transcripts =
        as.integer(
          expressed_annotated_transcripts
        ),

      retained_TPM_MFE =
        as.numeric(retained_TPM_MFE),

      sum_MFE_weights =
        as.numeric(sum_MFE_weights),

      TPM_weighted_Normalized_MFE =
        as.numeric(
          TPM_weighted_Normalized_MFE
        )
    ) %>%
    dplyr::filter(
      Phase %in% c(
        phase_1,
        phase_2
      )
    )

  if (
    !all(
      c(
        phase_1,
        phase_2
      ) %in% sample_values$Phase
    )
  ) {
    stop(
      analysis_id,
      ": one or both required phases are absent."
    )
  }

  if (
    any(
      abs(
        sample_values$sum_MFE_weights - 1
      ) > 1e-8
    )
  ) {
    stop(
      analysis_id,
      ": one or more TPM-weight sums differ from one."
    )
  }

  phase_stats <- sample_values %>%
    dplyr::group_by(
      Analysis_ID,
      Dataset,
      Phase
    ) %>%
    dplyr::summarise(
      n_samples = dplyr::n(),

      mean_TPM_weighted_Normalized_MFE =
        mean(
          TPM_weighted_Normalized_MFE,
          na.rm = TRUE
        ),

      sd_TPM_weighted_Normalized_MFE =
        if (dplyr::n() > 1L) {
          stats::sd(
            TPM_weighted_Normalized_MFE,
            na.rm = TRUE
          )
        } else {
          NA_real_
        },

      sem_TPM_weighted_Normalized_MFE =
        if (dplyr::n() > 1L) {
          stats::sd(
            TPM_weighted_Normalized_MFE,
            na.rm = TRUE
          ) /
            sqrt(dplyr::n())
        } else {
          NA_real_
        },

      .groups = "drop"
    )

  phase_1_row <- phase_stats %>%
    dplyr::filter(
      Phase == phase_1
    )

  phase_2_row <- phase_stats %>%
    dplyr::filter(
      Phase == phase_2
    )

  if (
    nrow(phase_1_row) != 1L ||
    nrow(phase_2_row) != 1L
  ) {
    stop(
      analysis_id,
      ": phase summaries could not be identified uniquely."
    )
  }

  phase_1_value <- getElement(
    phase_1_row,
    "mean_TPM_weighted_Normalized_MFE"
  )[[1]]

  phase_2_value <- getElement(
    phase_2_row,
    "mean_TPM_weighted_Normalized_MFE"
  )[[1]]

  phase_1_SE <- getElement(
    phase_1_row,
    "sem_TPM_weighted_Normalized_MFE"
  )[[1]]

  phase_2_SE <- getElement(
    phase_2_row,
    "sem_TPM_weighted_Normalized_MFE"
  )[[1]]

  observed_delta <-
    phase_1_value - phase_2_value

  delta_SE <- if (
    is.finite(phase_1_SE) &&
    is.finite(phase_2_SE)
  ) {
    sqrt(
      phase_1_SE^2 +
        phase_2_SE^2
    )
  } else {
    NA_real_
  }

  phase_profile <- TPM_long %>%
    dplyr::transmute(
      ensembl_transcript_id =
        as.character(
          ensembl_transcript_id
        ),

      Sample =
        as.character(Sample),

      Phase =
        Fig2E_standardize_phase(Phase),

      Normalized_MFE =
        as.numeric(Normalized_MFE),

      relative_TPM_weight =
        as.numeric(
          relative_transcript_abundance_MFE
        )
    ) %>%
    dplyr::filter(
      Phase %in% c(
        phase_1,
        phase_2
      ),

      is.finite(Normalized_MFE),
      is.finite(relative_TPM_weight),
      relative_TPM_weight >= 0
    ) %>%
    dplyr::group_by(
      ensembl_transcript_id,
      Normalized_MFE,
      Phase
    ) %>%
    dplyr::summarise(
      mean_relative_TPM_weight =
        mean(
          relative_TPM_weight,
          na.rm = TRUE
        ),

      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(
      names_from = Phase,
      values_from =
        mean_relative_TPM_weight,
      values_fill = 0
    )

  if (
    !(phase_1 %in% colnames(phase_profile)) ||
    !(phase_2 %in% colnames(phase_profile))
  ) {
    stop(
      analysis_id,
      ": failed to construct both phase-level TPM profiles."
    )
  }

  w_phase_1 <- as.numeric(
    getElement(
      phase_profile,
      phase_1
    )
  )

  w_phase_2 <- as.numeric(
    getElement(
      phase_profile,
      phase_2
    )
  )

  x_mfe <- as.numeric(
    phase_profile$Normalized_MFE
  )

  if (
    abs(sum(w_phase_1) - 1) > 1e-8 ||
    abs(sum(w_phase_2) - 1) > 1e-8
  ) {
    stop(
      analysis_id,
      ": phase-level mean TPM weights do not sum to one."
    )
  }

  profile_phase_1_value <- sum(
    w_phase_1 * x_mfe
  )

  profile_phase_2_value <- sum(
    w_phase_2 * x_mfe
  )

  profile_delta <-
    profile_phase_1_value -
    profile_phase_2_value

  if (
    abs(
      profile_phase_1_value -
      phase_1_value
    ) > 1e-10 ||

    abs(
      profile_phase_2_value -
      phase_2_value
    ) > 1e-10 ||

    abs(
      profile_delta -
      observed_delta
    ) > 1e-10
  ) {
    stop(
      analysis_id,
      paste(
        ": sample-summary and transcript-profile",
        "calculations do not agree."
      )
    )
  }

  contrast_weight <-
    w_phase_1 - w_phase_2

  if (
    abs(
      sum(contrast_weight)
    ) > 1e-8
  ) {
    stop(
      analysis_id,
      ": phase-contrast TPM weights do not sum to zero."
    )
  }

  n_transcripts <- length(x_mfe)

  active_indices <- which(
    contrast_weight != 0
  )

  active_contrast_weight <-
    contrast_weight[
      active_indices
    ]

  n_active_weights <- length(
    active_contrast_weight
  )

  if (n_active_weights == 0L) {
    stop(
      analysis_id,
      ": all phase-contrast TPM weights are zero."
    )
  }

  run_one_permutation <- function(iteration) {

    set.seed(
      seed + 10000L + iteration
    )

    sampled_MFE <- x_mfe[
      sample.int(
        n = n_transcripts,
        size = n_active_weights,
        replace = FALSE
      )
    ]

    sum(
      active_contrast_weight *
        sampled_MFE
    )
  }

  if (
    .Platform$OS.type != "windows" &&
    n_cores > 1L
  ) {

    permutation_delta <- unlist(
      parallel::mclapply(
        X = seq_len(B_perm),
        FUN = run_one_permutation,
        mc.cores = n_cores,
        mc.preschedule = TRUE,
        mc.set.seed = FALSE
      ),
      use.names = FALSE
    )

  } else {

    permutation_delta <- vapply(
      X = seq_len(B_perm),
      FUN = run_one_permutation,
      FUN.VALUE = numeric(1)
    )
  }

  if (
    length(permutation_delta) != B_perm ||
    any(!is.finite(permutation_delta))
  ) {
    stop(
      analysis_id,
      ": one or more permutation statistics are invalid."
    )
  }

  permutation_p_value <- (
    sum(
      abs(permutation_delta) >=
        abs(observed_delta)
    ) + 1
  ) / (B_perm + 1)

  contrast_stats <- tibble::tibble(
    Analysis_ID = analysis_id,
    Dataset = dataset_label,
    Comparison = comparison_label,
    phase_1 = phase_1,
    phase_2 = phase_2,

    n_phase_1 =
      getElement(
        phase_1_row,
        "n_samples"
      )[[1]],

    n_phase_2 =
      getElement(
        phase_2_row,
        "n_samples"
      )[[1]],

    phase_1_value =
      phase_1_value,

    phase_2_value =
      phase_2_value,

    delta_TPM_weighted_Normalized_MFE =
      observed_delta,

    phase_1_SE =
      phase_1_SE,

    phase_2_SE =
      phase_2_SE,

    delta_SE =
      delta_SE,

    errorbar_low =
      ifelse(
        is.finite(delta_SE),
        observed_delta - delta_SE,
        NA_real_
      ),

    errorbar_high =
      ifelse(
        is.finite(delta_SE),
        observed_delta + delta_SE,
        NA_real_
      ),

    percent_change_phase1_minus_phase2 =
      100 *
      observed_delta /
      abs(phase_2_value),

    abundance_measure =
      "Transcript-level TPM",

    transcript_universe =
      paste(
        "Transcripts matched to Table_ANNOTATED",
        "with finite Normalized_MFE"
      ),

    n_transcripts_used =
      n_transcripts,

    n_transcripts_expressed_in_either_phase =
      sum(
        w_phase_1 > 0 |
          w_phase_2 > 0
      ),

    n_nonzero_contrast_weights =
      n_active_weights,

    permutation_p_value =
      permutation_p_value,

    permutation_test =
      paste(
        "Normalized_MFE permuted across transcripts;",
        "phase-specific TPM contrast weights held fixed"
      ),

    B_perm =
      B_perm,

    permutation_cores =
      n_cores,

    direction =
      dplyr::case_when(
        observed_delta < 0 ~
          paste0(
            phase_1,
            " is more negative than ",
            phase_2
          ),

        observed_delta > 0 ~
          paste0(
            phase_1,
            " is less negative than ",
            phase_2
          ),

        TRUE ~
          "No change"
      )
  )

  list(
    sample_values =
      sample_values,

    phase_stats =
      phase_stats,

    contrast_stats =
      contrast_stats
  )
}

Fig2E_TPM_results <- vector(
  mode = "list",
  length = nrow(
    Fig2E_analysis_config
  )
)

for (
  analysis_index in seq_len(
    nrow(
      Fig2E_analysis_config
    )
  )
) {

  current_config <-
    Fig2E_analysis_config[
      analysis_index,
      ,
      drop = FALSE
    ]

  message(
    "Running Fig. 2E TPM analysis: ",
    current_config$Analysis_ID[[1]]
  )

  Fig2E_TPM_results[[analysis_index]] <-
    Fig2E_run_one_TPM_contrast(
      TPM_long = figure2_objects[[
        current_config$TPM_long_object[[1]]
      ]],

      TPM_sample = figure2_objects[[
        current_config$TPM_sample_object[[1]]
      ]],

      analysis_id =
        current_config$
        Analysis_ID[[1]],

      dataset_label =
        current_config$
        Dataset[[1]],

      comparison_label =
        current_config$
        Comparison[[1]],

      phase_1 =
        current_config$
        Phase_1[[1]],

      phase_2 =
        current_config$
        Phase_2[[1]],

      B_perm =
        Fig2E_B_perm,

      n_cores =
        Fig2E_n_cores,

      seed =
        current_config$
        Seed[[1]]
    )

  invisible(gc())
}

names(Fig2E_TPM_results) <-
  Fig2E_analysis_config$Analysis_ID

Fig2E_TPM_sample_values <- dplyr::bind_rows(
  lapply(
    Fig2E_TPM_results,
    function(current_result) {
      current_result$sample_values
    }
  )
) %>%
  dplyr::distinct()

Fig2E_TPM_phase_stats <- dplyr::bind_rows(
  lapply(
    Fig2E_TPM_results,
    function(current_result) {
      current_result$phase_stats
    }
  )
) %>%
  dplyr::distinct()

Fig2E_TPM_delta_stats <- dplyr::bind_rows(
  lapply(
    Fig2E_TPM_results,
    function(current_result) {
      current_result$contrast_stats
    }
  )
) %>%
  dplyr::mutate(
    permutation_p_adj_BH =
      stats::p.adjust(
        permutation_p_value,
        method = "BH"
      ),

    permutation_p_label =
      dplyr::case_when(
        permutation_p_value < 0.001 ~
          "P < 0.001",

        TRUE ~
          paste0(
            "P = ",
            format(
              signif(
                permutation_p_value,
                2
              ),
              scientific = FALSE,
              trim = TRUE
            )
          )
      ),

    permutation_p_adj_BH_label =
      dplyr::case_when(
        permutation_p_adj_BH < 0.001 ~
          "BH-adjusted P < 0.001",

        TRUE ~
          paste0(
            "BH-adjusted P = ",
            format(
              signif(
                permutation_p_adj_BH,
                2
              ),
              scientific = FALSE,
              trim = TRUE
            )
          )
      ),

    Dataset = factor(
      Dataset,
      levels = c(
        "MCF7",
        "HeLa"
      )
    ),

    Comparison = factor(
      Comparison,
      levels = c(
        "G1 - G2/M",
        "S - G1"
      )
    ),

    delta_label = sprintf(
      "%+.4f",
      delta_TPM_weighted_Normalized_MFE
    )
  ) %>%
  dplyr::arrange(
    Dataset,
    Comparison
  )

stopifnot(
  nrow(
    Fig2E_TPM_delta_stats
  ) == 4,

  all(
    is.finite(
      Fig2E_TPM_delta_stats$
        delta_TPM_weighted_Normalized_MFE
    )
  ),

  all(
    is.finite(
      Fig2E_TPM_delta_stats$
        permutation_p_value
    )
  ),

  all(
    Fig2E_TPM_delta_stats$
      permutation_p_value > 0 &
      Fig2E_TPM_delta_stats$
      permutation_p_value <= 1
  ),

  all(
    Fig2E_TPM_delta_stats$
      abundance_measure ==
      "Transcript-level TPM"
  )
)

Fig2E_errorbar_values <- c(
  Fig2E_TPM_delta_stats$
    delta_TPM_weighted_Normalized_MFE,

  Fig2E_TPM_delta_stats$
    errorbar_low,

  Fig2E_TPM_delta_stats$
    errorbar_high
)

Fig2E_delta_y_limit <- max(
  abs(
    Fig2E_errorbar_values[
      is.finite(
        Fig2E_errorbar_values
      )
    ]
  ),
  na.rm = TRUE
) * 1.35

if (
  !is.finite(Fig2E_delta_y_limit) ||
  Fig2E_delta_y_limit <= 0
) {
  Fig2E_delta_y_limit <- 0.01
}

make_Fig2E_TPM_plot <- function(
    plot_data,
    plot_title
) {

  ggplot(
    plot_data,
    aes(
      x = Comparison,
      y =
        delta_TPM_weighted_Normalized_MFE
    )
  ) +
    geom_hline(
      yintercept = 0,
      linewidth = 0.4,
      linetype = "dashed",
      color = "black"
    ) +
    geom_col(
      width = 0.75,
      fill = "#D55E00",
      color = "black",
      linewidth = 0.3
    ) +
    geom_errorbar(
      data = plot_data %>%
        dplyr::filter(
          is.finite(delta_SE)
        ),

      aes(
        x = Comparison,
        ymin = errorbar_low,
        ymax = delta_TPM_weighted_Normalized_MFE
      ),

      inherit.aes = FALSE,
      width = 0.18,
      linewidth = 0.4,
      color = "black"
    ) +
    geom_text(
      aes(
        y = ifelse(
          is.finite(errorbar_low),
          errorbar_low - 0.0006,
          delta_TPM_weighted_Normalized_MFE
        ),

        label = delta_label,

        vjust = ifelse(
          is.finite(errorbar_low),
          1,
          ifelse(
            delta_TPM_weighted_Normalized_MFE >= 0,
            -0.45,
            1.25
          )
        )
      ),

      size = 3.5,
      color = "black"
    ) +
    scale_y_continuous(
      labels = scales::number_format(
        accuracy = 0.001
      ),

      limits = c(
        -Fig2E_delta_y_limit,
        Fig2E_delta_y_limit
      ),

      expand = expansion(
        mult = c(
          0.05,
          0.05
        )
      )
    ) +
    theme_classic() +
    theme(
      legend.position = "none",

      axis.text.x = element_text(
        size = 11,
        color = "black",
        angle = 45,
        hjust = 1
      ),

      axis.text.y = element_text(
        size = 11,
        color = "black"
      ),

      axis.title.y = element_text(
        size = 14
      ),

      plot.title = element_text(
        size = 14,
        hjust = 0.5
      )
    ) +
    labs(
      title = plot_title,
      subtitle = "",
      x = "",

      y = expression(
        Delta * Delta * G[fold] ~
          (kcal %.% mol^{-1} %.% nt^{-1})
      )
    )
}

Fig2E_MCF7_TPM_plot_data <-
  Fig2E_TPM_delta_stats %>%
  dplyr::filter(
    Dataset == "MCF7"
  ) %>%
  droplevels()

Plot_Cell_cycle_MCF7_TPM_weighted_DGfold_delta <-
  make_Fig2E_TPM_plot(
    plot_data =
      Fig2E_MCF7_TPM_plot_data,

    plot_title =
      "MCF7"
  )

print(
  Plot_Cell_cycle_MCF7_TPM_weighted_DGfold_delta
)

ggsave(
  filename =
    "Plot_Cell_cycle_MCF7_TPM_weighted_DGfold_delta.tiff",

  plot =
    Plot_Cell_cycle_MCF7_TPM_weighted_DGfold_delta,

  path =
    Fig2E_output_dir,

  device = "tiff",
  width = 2.5,
  height = 4.5,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

Fig2E_HeLa_TPM_plot_data <-
  Fig2E_TPM_delta_stats %>%
  dplyr::filter(
    Dataset == "HeLa"
  ) %>%
  droplevels()

Plot_Cell_cycle_HeLa_TPM_weighted_DGfold_delta <-
  make_Fig2E_TPM_plot(
    plot_data =
      Fig2E_HeLa_TPM_plot_data,

    plot_title =
      "HeLa"
  )

print(
  Plot_Cell_cycle_HeLa_TPM_weighted_DGfold_delta
)

ggsave(
  filename =
    "Plot_Cell_cycle_HeLa_TPM_weighted_DGfold_delta.tiff",

  plot =
    Plot_Cell_cycle_HeLa_TPM_weighted_DGfold_delta,

  path =
    Fig2E_output_dir,

  device = "tiff",
  width = 2.5,
  height = 4.5,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

Fig2E_TPM_delta_relative_stats <-
  Fig2E_TPM_delta_stats %>%
  dplyr::select(
    Analysis_ID,
    Dataset,
    Comparison,
    phase_1,
    phase_2,
    n_phase_1,
    n_phase_2,
    phase_1_value,
    phase_2_value,
    delta_TPM_weighted_Normalized_MFE,
    phase_1_SE,
    phase_2_SE,
    delta_SE,
    errorbar_low,
    errorbar_high,
    percent_change_phase1_minus_phase2,
    abundance_measure,
    transcript_universe,
    n_transcripts_used,
    n_transcripts_expressed_in_either_phase,
    n_nonzero_contrast_weights,
    permutation_p_value,
    permutation_p_adj_BH,
    permutation_p_label,
    permutation_p_adj_BH_label,
    permutation_test,
    direction,
    B_perm,
    permutation_cores
  ) %>%
  dplyr::arrange(
    Dataset,
    Comparison
  )

cat(
  "\nFig. 2E TPM-weighted phase-transition statistics:\n"
)

print(
  as.data.frame(
    Fig2E_TPM_delta_relative_stats
  ),
  row.names = FALSE,
  digits = 12
)

cat("\nFigure 2 analysis complete.\n")
cat("Outputs written to:\n", output_dir, "\n")
