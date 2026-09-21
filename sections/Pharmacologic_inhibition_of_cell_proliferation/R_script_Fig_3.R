input_dir <- "../../data"
scripts_dir <- "."
figure_output_dir <- "Figure_3_outputs"

configured_paths <- c(input_dir, scripts_dir, figure_output_dir)
if (any(startsWith(configured_paths, "path_to_"))) {
  stop("Set input_dir, scripts_dir, and figure_output_dir before running this script.")
}

dir.create(figure_output_dir, showWarnings = FALSE, recursive = TRUE)

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
library(stringr)
library(clusterProfiler)
library(pathview)
library(enrichplot)
library(grid)
library(scales)

master_script_path <- file.path(scripts_dir, "Cross_validation_MFE.R")
analysis_results_file <- file.path(
  input_dir,
  "pharmacologic_inhibition_analysis_results_for_Fig_3.rds"
)

if (!file.exists(master_script_path)) {
  stop("Cross_validation_MFE.R was not found in scripts_dir.")
}

if (!file.exists(analysis_results_file)) {
  stop("The precomputed Figure 3 results file was not found in input_dir.")
}

figure3_objects <- readRDS(analysis_results_file)

if (!is.list(figure3_objects) || is.null(names(figure3_objects))) {
  stop("The precomputed Figure 3 results file must contain a named list.")
}

invisible(
  list2env(
    figure3_objects,
    envir = .GlobalEnv
  )
)

gsea_output_dir <- file.path(figure_output_dir, "GSEA_KEGG_clusterProfiler")
dir.create(gsea_output_dir, showWarnings = FALSE, recursive = TRUE)

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

Expression_changes_GSE241735_ICC13_7_Futibatinib_4h <- res_GSE241735_Futibatinib_4h_gene_df %>%
  dplyr::select(Gene_symbol, log2FoldChange) %>%
  dplyr::rename(
    GSE241735_ICC13_7_Futibatinib_4h_Log_FC = log2FoldChange
  )

Expression_changes_GSE241735_ICC13_7_Futibatinib_12h <- res_GSE241735_Futibatinib_12h_gene_df %>%
  dplyr::select(Gene_symbol, log2FoldChange) %>%
  dplyr::rename(
    GSE241735_ICC13_7_Futibatinib_12h_Log_FC = log2FoldChange
  )

Expression_changes_GSE241735_ICC13_7_Futibatinib_24h <- res_GSE241735_Futibatinib_24h_gene_df %>%
  dplyr::select(Gene_symbol, log2FoldChange) %>%
  dplyr::rename(
    GSE241735_ICC13_7_Futibatinib_24h_Log_FC = log2FoldChange
  )

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
  dplyr::full_join(Expression_changes_GSE241735_ICC13_7_Futibatinib_4h, by = "Gene_symbol") %>%
  dplyr::full_join(Expression_changes_GSE241735_ICC13_7_Futibatinib_12h, by = "Gene_symbol") %>%
  dplyr::full_join(Expression_changes_GSE241735_ICC13_7_Futibatinib_24h, by = "Gene_symbol")

master_script_lines <- readLines(
  master_script_path,
  warn = FALSE
)

dataset_heading_lines <- grep(
  "^# GSE[0-9]+\\b.*####",
  master_script_lines
)

GSE163249_start <- grep(
  "^# GSE163249\\b.*####",
  master_script_lines
)[1]

next_dataset_start <- dataset_heading_lines[
  dataset_heading_lines > GSE163249_start
][1]

if (is.na(next_dataset_start)) {
  cancer_cell_line_script <- master_script_lines
} else {
  cancer_cell_line_script <- master_script_lines[
    seq_len(next_dataset_start - 1)
  ]
}

Fig3A_screened_datasets <- cancer_cell_line_script %>%
  stringr::str_subset("^# GSE[0-9]+\\b.*####") %>%
  stringr::str_extract("GSE[0-9]+") %>%
  unique()

n_Fig3A_screened_datasets <- length(
  Fig3A_screened_datasets
)

Fig3A_screened_profiles <- cancer_cell_line_script %>%
  stringr::str_subset(
    "^MFE(_res)?_GSE[0-9]+.*<- dplyr::inner_join"
  ) %>%
  stringr::str_extract(
    "^MFE(_res)?_GSE[[:alnum:]_]+"
  ) %>%
  unique()

n_Fig3A_screened_profiles <- length(
  Fig3A_screened_profiles
)

Fig3A_design_selected_datasets <- c(
  "GSE105083",
  "GSE129221",
  "GSE135514",
  "GSE137535",
  "GSE143462",
  "GSE148096",
  "GSE177054",
  "GSE178755",
  "GSE103021",
  "GSE118548",
  "GSE93124",
  "GSE121105",
  "GSE241735"
)

Fig3A_design_selected_datasets <- intersect(
  Fig3A_design_selected_datasets,
  Fig3A_screened_datasets
)

n_Fig3A_design_selected_datasets <- length(
  Fig3A_design_selected_datasets
)

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host="https://useast.ensembl.org")
gene_ids <- (Expression_changes_targeted_agents$Gene_symbol)

res_gene_entrez_id <- getBM(attributes = c('entrezgene_id',
                                           'external_gene_name'),
                            filters = 'external_gene_name',
                            values = gene_ids,
                            mart = mart)
res_gene_entrez_id <- dplyr::filter(res_gene_entrez_id, !is.na(entrezgene_id)) %>%
  dplyr::rename(Gene_symbol = external_gene_name)

Expression_changes_targeted_agents_genes <- dplyr::left_join(Expression_changes_targeted_agents, res_gene_entrez_id, by = "Gene_symbol")

write.csv(Expression_changes_targeted_agents_genes, file = file.path(gsea_output_dir, "Expression_changes_targeted_agents_genes.csv"))

process_and_save_plots <- function(data_frame, log_fc_cols, plot_titles, output_dir) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  plot_list <- list()
  kk2_list <- list()

  for (i in seq_along(log_fc_cols)) {
    log_fc_col <- log_fc_cols[i]
    plot_title <- plot_titles[i]
    file_name <- file.path(output_dir, paste0(plot_title, "_GSEA_KEGG.pdf"))

    all_values_entrez <- data_frame %>%
      dplyr::select(entrezgene_id, !!sym(log_fc_col)) %>%
      dplyr::filter(!is.na(!!sym(log_fc_col))) %>%
      dplyr::filter(!is.na(entrezgene_id))

    all_values_Log2_FC_entrez <- all_values_entrez[[log_fc_col]]
    names(all_values_Log2_FC_entrez) <- all_values_entrez$entrezgene_id

    all_values_Log2_FC_entrez <- all_values_Log2_FC_entrez[!duplicated(names(all_values_Log2_FC_entrez))]

    all_values_Log2_FC_entrez <- na.omit(all_values_Log2_FC_entrez)
    all_values_Log2_FC_entrez <- sort(all_values_Log2_FC_entrez, decreasing = TRUE)

    cat("Sample:", plot_titles[i], "\n")
    cat("Number of genes:", length(all_values_Log2_FC_entrez), "\n")

    if (length(all_values_Log2_FC_entrez) == 0) {
      cat("No genes available for", plot_title, "\n")
      next
    }

    kk2 <- gseKEGG(all_values_Log2_FC_entrez,
                   organism = "hsa",
                   keyType = "kegg",
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "BH",
                   minGSSize = 10,
                   maxGSSize = 500)

    if (is.null(kk2) || nrow(kk2) == 0) {
      cat("No significant pathways found for", plot_title, "\n")
      next
    }

    kk2_list[[plot_title]] <- kk2

    p <- dotplot(kk2,
                 x = "GeneRatio",
                 color = "p.adjust",
                 showCategory = 50,
                 size = NULL,
                 split = ".sign",
                 font.size = 8,
                 title = plot_title,
                 orderBy = "x",
                 label_format = 100) +
      facet_grid(~.sign) +
      theme(panel.spacing = unit(0.5, "cm"))

    ggsave(file_name, plot = p, width = 8, height = 12)
  }

  combined_kk2 <- do.call(rbind, lapply(names(kk2_list), function(name) {
    df <- as.data.frame(kk2_list[[name]])
    df$Sample <- name
    return(df)
  }))

  return(combined_kk2)
}
data_frame <- Expression_changes_targeted_agents_genes

log_fc_cols <- c(
  "GSE93124_SW480_crizotinib_Log_FC",
  "GSE178755_HCC827_erlotinib_Log_FC",
  "GSE137535_SH454_Log_FC",
  "GSE105083_PANC1_CM03_24h_Log_FC",
  "GSE118548_T84_Trametinib_Log_FC",
  "GSE148096_TNBC_erlotinib_Log_FC",
  "GSE143462_KYSE150_omipalisib_Log_FC",
  "GSE121105_BT474_trastuzumab_Log_FC",
  "GSE103021_LU65_ARS_1620_Log_FC",

  "GSE135514_MDAMB468_Das_24h_Log_FC",
  "GSE135514_MDAMB468_Das_72h_Log_FC",

  "GSE241735_ICC13_7_Futibatinib_4h_Log_FC",
  "GSE241735_ICC13_7_Futibatinib_12h_Log_FC",
  "GSE241735_ICC13_7_Futibatinib_24h_Log_FC",

  "GSE105083_MIA_PACA2_CM03_6h_Log_FC",
  "GSE105083_MIA_PACA2_CM03_24h_Log_FC",
  "GSE129221_PC9GR_apatinib_Log_FC"
)

plot_titles <- c(
  "GSE93124",
  "GSE178755",
  "GSE137535_SH454",
  "GSE105083_PANC1_CM03_24h",
  "GSE118548_tra",
  "GSE148096",
  "GSE143462",
  "GSE121105",
  "GSE103021",

  "GSE135514_24h",
  "GSE135514_72h",

  "GSE241735_Futibatinib_4h",
  "GSE241735_Futibatinib_12h",
  "GSE241735_Futibatinib_24h",

  "GSE105083_MIA_PACA2_6h",
  "GSE105083_MIA_PACA2_24h",
  "GSE129221_PC9GR_apatinib"
)

output_dir <- gsea_output_dir

combined_kk2_df <- process_and_save_plots(data_frame, log_fc_cols, plot_titles, output_dir)

GSE241735_timecourse_samples <- c(
  "GSE241735_Futibatinib_4h",
  "GSE241735_Futibatinib_12h",
  "GSE241735_Futibatinib_24h"
)

write.csv(combined_kk2_df, file = file.path(figure_output_dir, "combined_kk2_df.csv"))

Selected_KEGG_pathways <- c("hsa04150", "hsa04668", "hsa04657", "hsa03030", "hsa03410", "hsa03430", "hsa03420",
                            "hsa03320", "hsa03460", "hsa03440", "hsa04014", "hsa04066", "hsa04070",
                            "hsa04310", "hsa03082", "hsa04015", "hsa04625", "hsa04371", "hsa04390", "hsa04010", "hsa04622", "hsa04350",
                            "hsa04152", "hsa04926", "hsa04630", "hsa04115", "hsa04151", "hsa04910", "hsa04330", "hsa04012",
                            "hsa04110")

combined_kk2_df_filtered <- combined_kk2_df %>%
  dplyr::filter(
    ID %in% Selected_KEGG_pathways
  ) %>%
  dplyr::mutate(
    NES_raw = NES
  )

missing_GSE241735_samples <- setdiff(
  GSE241735_timecourse_samples,
  unique(combined_kk2_df_filtered$Sample)
)

if (length(missing_GSE241735_samples) > 0) {

  KEGG_pathway_descriptions <- combined_kk2_df %>%
    dplyr::select(ID, Description) %>%
    dplyr::distinct(ID, .keep_all = TRUE)

  GSE241735_placeholder_rows <- tidyr::expand_grid(
    Sample = missing_GSE241735_samples,
    ID = Selected_KEGG_pathways
  ) %>%
    dplyr::left_join(
      KEGG_pathway_descriptions,
      by = "ID"
    ) %>%
    dplyr::mutate(
      Description = ifelse(
        is.na(Description),
        ID,
        Description
      ),
      NES_raw = 0,
      NES = 0,
      p.adjust = 1
    )

  combined_kk2_df_filtered <- dplyr::bind_rows(
    combined_kk2_df_filtered,
    GSE241735_placeholder_rows
  )
}

combined_kk2_df_filtered <- combined_kk2_df_filtered %>%
  dplyr::mutate(
    Description = dplyr::case_when(
      ID == "hsa04070" ~ "Phosphatidylinositol signaling system",
      ID == "hsa04622" ~ "RIG-I-like receptor signaling pathway",
      ID == "hsa04350" ~ "TGF-beta signaling pathway",
      ID == "hsa04926" ~ "Relaxin signaling pathway",
      ID == "hsa04630" ~ "JAK-STAT signaling pathway",
      ID == "hsa04151" ~ "PI3K-Akt signaling pathway",
      ID == "hsa04910" ~ "Insulin signaling pathway",
      ID == "hsa04330" ~ "Notch signaling pathway",
      ID == "hsa04012" ~ "ErbB signaling pathway",
      TRUE ~ Description
    )
  )

generate_heatmap <- function(kk2_df) {

  kk2_df <- kk2_df %>%
    dplyr::select(
      Description,
      Sample,
      NES
    ) %>%
    tidyr::spread(
      Sample,
      NES
    )

  kk2_df[is.na(kk2_df)] <- 0

  rownames(kk2_df) <- kk2_df$Description
  kk2_df <- kk2_df[, -1]

  display_matrix <- as.matrix(kk2_df)

  maximum_absolute_NES <- max(
    abs(display_matrix),
    na.rm = TRUE
  )

  if (
    !is.finite(maximum_absolute_NES) ||
    maximum_absolute_NES == 0
  ) {
    maximum_absolute_NES <- 1
  }

  heatmap_breaks <- seq(
    -maximum_absolute_NES,
    maximum_absolute_NES,
    length.out = 51
  )

  display_numbers <- ifelse(
    display_matrix == 0,
    "",
    sprintf("%.2g", display_matrix)
  )

  pheatmap(
    display_matrix,
    color = colorRampPalette(
      c("blue", "white", "red")
    )(50),
    breaks = heatmap_breaks,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    display_numbers = display_numbers,
    number_format = "%.2g",
    main = "Heatmap of KEGG signaling pathways"
  )
}

Heatmap_cell_cycle <- generate_heatmap(combined_kk2_df_filtered)
print(Heatmap_cell_cycle)

write.csv(combined_kk2_df_filtered, file = file.path(figure_output_dir, "combined_kk2_results.csv"), row.names = FALSE)

Fig3B_original_selected_samples <- c(
  "GSE103021",
  "GSE121105",
  "GSE143462",
  "GSE135514_24h",
  "GSE135514_72h"
)

Fig3B_sample_order <- c(
  Fig3B_original_selected_samples,
  "GSE241735_Futibatinib_4h",
  "GSE241735_Futibatinib_12h",
  "GSE241735_Futibatinib_24h"
)

Cell_cycle_arrest_samples <- combined_kk2_df_filtered %>%
  dplyr::filter(
    ID == "hsa04110",
    NES_raw < 0
  ) %>%
  dplyr::select(Sample) %>%
  dplyr::distinct()

Cell_cycle_arrest_samples <- dplyr::bind_rows(
  Cell_cycle_arrest_samples,
  data.frame(
    Sample = GSE241735_timecourse_samples
  )
) %>%
  dplyr::distinct()

Seleted_non_proliferating_samples <- combined_kk2_df_filtered %>%
  dplyr::filter(
    Sample %in% Cell_cycle_arrest_samples$Sample
  )

generate_heatmap <- function(kk2_df) {

  kk2_df <- kk2_df %>%
    dplyr::select(
      Description,
      Sample,
      NES
    ) %>%
    tidyr::spread(
      Sample,
      NES
    )

  kk2_df[is.na(kk2_df)] <- 0

  rownames(kk2_df) <- kk2_df$Description
  kk2_df <- kk2_df[, -1]

  display_matrix <- as.matrix(kk2_df)

  maximum_absolute_NES <- max(
    abs(display_matrix),
    na.rm = TRUE
  )

  if (
    !is.finite(maximum_absolute_NES) ||
    maximum_absolute_NES == 0
  ) {
    maximum_absolute_NES <- 1
  }

  heatmap_colours <- colorRampPalette(
    c("blue", "white", "red")
  )(51)

  heatmap_breaks <- seq(
    -maximum_absolute_NES,
    maximum_absolute_NES,
    length.out = 52
  )

  display_numbers <- ifelse(
    display_matrix == 0,
    "",
    sprintf("%.2g", display_matrix)
  )

  pheatmap(
    display_matrix,
    color = heatmap_colours,
    breaks = heatmap_breaks,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    display_numbers = display_numbers,
    number_format = "%.2g",
    main = "Heatmap of KEGG signaling pathways"
  )
}

Heatmap_cell_cycle <- generate_heatmap(Seleted_non_proliferating_samples)
print(Heatmap_cell_cycle)

ids_to_check <- c(
  "hsa04150", "hsa04668", "hsa04657", "hsa03030", "hsa03410", "hsa03430", "hsa03420",
  "hsa03320", "hsa03460", "hsa03440", "hsa04014", "hsa04066", "hsa04070",
  "hsa04310", "hsa03082", "hsa04015", "hsa04625", "hsa04371", "hsa04390", "hsa04010", "hsa04622", "hsa04350",
  "hsa04152", "hsa04926", "hsa04630", "hsa04115", "hsa04151", "hsa04910", "hsa04330", "hsa04012"
)

samples_to_exclude <- Seleted_non_proliferating_samples %>%
  dplyr::filter(
    ID %in% ids_to_check,
    NES_raw > 0,
    !Sample %in% GSE241735_timecourse_samples
  ) %>%
  dplyr::select(Sample) %>%
  dplyr::distinct()

Cell_cycle_arrest_mitogenic_signaling <-
  Seleted_non_proliferating_samples %>%
  dplyr::filter(
    !Sample %in% samples_to_exclude$Sample
  )

Fig3B_selected_samples <- dplyr::bind_rows(

  Cell_cycle_arrest_mitogenic_signaling %>%
    dplyr::filter(
      Sample %in% Fig3B_original_selected_samples
    ),

  combined_kk2_df_filtered %>%
    dplyr::filter(
      Sample %in% GSE241735_timecourse_samples
    )

) %>%
  dplyr::distinct()

setdiff(
  Fig3B_sample_order,
  unique(Fig3B_selected_samples$Sample)
)

Fig3A_final_selected_datasets <- Fig3B_selected_samples %>%
  dplyr::transmute(
    Dataset = stringr::str_extract(
      as.character(Sample),
      "GSE[0-9]+"
    )
  ) %>%
  dplyr::filter(
    !is.na(Dataset),
    Dataset %in% Fig3A_design_selected_datasets
  ) %>%
  dplyr::distinct(Dataset) %>%
  dplyr::pull(Dataset)

n_Fig3A_final_selected_datasets <- length(
  Fig3A_final_selected_datasets
)

Fig3A_selection_numbers <- data.frame(
  Selection_step = c(
    "GREIN datasets screened",
    "Datasets meeting design criteria",
    "Cell lines with repression of proliferation pathways"
  ),
  Number = c(
    n_Fig3A_screened_datasets,
    n_Fig3A_design_selected_datasets,
    n_Fig3A_final_selected_datasets
  )
)

print(Fig3A_selection_numbers)

cat(
  "\nIndividual treatment profiles screened:",
  n_Fig3A_screened_profiles,
  "\n\n"
)

cat("Datasets screened:\n")
print(Fig3A_screened_datasets)

cat("\nDatasets meeting design criteria:\n")
print(Fig3A_design_selected_datasets)

cat("\nFinal selected datasets:\n")
print(Fig3A_final_selected_datasets)

generate_heatmap <- function(kk2_df) {

  kk2_df <- kk2_df %>%
    dplyr::select(
      Description,
      Sample,
      NES
    ) %>%
    tidyr::spread(
      Sample,
      NES
    )

  kk2_df[is.na(kk2_df)] <- 0

  kk2_df <- kk2_df %>%
    dplyr::select(
      Description,
      dplyr::all_of(Fig3B_sample_order)
    )

  rownames(kk2_df) <- kk2_df$Description

  display_matrix <- as.matrix(
    kk2_df[, -1, drop = FALSE]
  )

  maximum_absolute_NES <- max(
    abs(display_matrix),
    na.rm = TRUE
  )

  if (
    !is.finite(maximum_absolute_NES) ||
    maximum_absolute_NES == 0
  ) {
    maximum_absolute_NES <- 1
  }

  heatmap_colours <- colorRampPalette(
    c("blue", "white", "red")
  )(51)

  heatmap_breaks <- seq(
    -maximum_absolute_NES,
    maximum_absolute_NES,
    length.out = 52
  )

  pheatmap(
    display_matrix,
    color = heatmap_colours,
    breaks = heatmap_breaks,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    display_numbers = FALSE,
    main = "Heatmap of KEGG pathways",
    fontsize_row = 12,
    fontsize_col = 12
  )
}

Heatmap_cell_cycle <- generate_heatmap(
  Fig3B_selected_samples
)
print(Heatmap_cell_cycle)

Plot_GREIN_selected_samples <- Fig3B_selected_samples %>%
  dplyr::select(
    ID,
    Description,
    NES,
    p.adjust,
    Sample
  )

all_ids <- c("hsa04110", "hsa03030", "hsa04010", "hsa03082", "hsa04625",
             "hsa03460", "hsa03440", "hsa03420", "hsa03430", "hsa03410",
             "hsa04115", "hsa04014", "hsa03320", "hsa04150", "hsa04390",
             "hsa04310", "hsa04657", "hsa04015", "hsa04668", "hsa04371",
             "hsa04151", "hsa04910", "hsa04330", "hsa04350",
             "hsa04630", "hsa04012", "hsa04152")

df_complete <- Plot_GREIN_selected_samples %>%
  dplyr::select(Sample, ID) %>%
  dplyr::distinct() %>%
  complete(Sample, ID = all_ids) %>%
  dplyr::left_join(Plot_GREIN_selected_samples, by = c("Sample", "ID")) %>%
  tidyr::replace_na(
    list(
      Description = "",
      NES = 0,
      p.adjust = 1
    )
  )

Description_names <- combined_kk2_df_filtered %>%
  dplyr::select(
    ID,
    Description
  ) %>%
  dplyr::distinct(
    ID,
    .keep_all = TRUE
  )

Plot_GREIN_selected_samples_complete <- dplyr::left_join(df_complete, Description_names, by = "ID") %>%
  dplyr::select(-Description.x) %>%
  dplyr::rename(Description = Description.y) %>%
  dplyr::mutate(Description = gsub("signaling ", "", Description)) %>%
  filter(!Sample %in% c("X778_PALBO_72h", "X778_PALBO_EP", "X778_nutlin_72h", "X778_TUNIC_72h"))

ggplot(Plot_GREIN_selected_samples_complete, aes(x = Sample, y = Description, fill = NES)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low = "blue", high = "white") +
  theme_minimal() +
  labs(title = "",
       x = "",
       y = "",
       fill = "NES") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 14))

proliferation_science_genes <- c(
  "TPX2", "PRC1", "FOXM1", "CDK1", "C20orf24", "MCM2", "H2AFZ", "TOP2A", "PCNA", "UBE2C",
  "MELK", "TRIP13", "NCAPD2", "MCM7", "RNASEH2A", "RAD51AP1", "KIF20A", "CDC45", "MAD2L1", "ESPL1",
  "CCNB2", "FEN1", "TTK", "CCT5", "RFC4", "ATAD2", "CKAP5", "NUP205", "CDC20", "CKS2",
  "RRM2", "ELAVL1", "CCNB1", "RRM1", "AURKB", "MSH6", "EZH2", "CTPS", "DKC1", "OIP5",
  "CDCA8", "PTTG1", "CEP55", "H2AFX", "CMAS", "NCAPH", "MCM10", "LSM4", "NCAPG2", "ASF1B",
  "ZWINT", "PBK", "ZWILCH", "CDCA3", "ECT2", "CDC6", "UNG", "MTCH2", "RAD21", "ACTL6A",
  "GPI", "SRSF2", "HDGF", "NXT1", "NEK2", "DHCR7", "AURKA", "NDUFAB1", "TMEM194A", "KIF4A"
)

meta_pcna_genes <- c(
  "PCNA", "CKS2", "NUSAP1", "RRM2", "ZWINT", "PRC1", "TFDP1", "CCNA2", "CCNB1", "MELK",
  "TPX2", "BIRC5", "NCAPG2", "RFWD3", "TACC3", "CDK1", "PCLAF", "MCM2", "MCM7", "MCM5",
  "ASF1B", "CCNB2", "GTPBP2", "KIF20A", "PTTG1", "AURKA", "CDC20", "FOXM1", "GINS2", "MAD2L1",
  "UBE2C", "GINS1", "MCM6", "NCAPD2", "NUP37", "CKS1B", "MCM7", "ARID3A", "AURKB", "CDCA8",
  "SNRPD1", "SNRPD2", "TROAP", "MIS18A", "DDX39A", "AHSP", "ESPL1", "HMBS", "LSM6", "MCM4",
  "NCAPD3", "PPIH", "CDC45", "DTL", "EPB42", "HMGN2", "MCM3", "PF4", "PPBP", "RFC4",
  "TOP2A", "BUB1B", "CDT1", "FEN1", "GATA1", "GYPB", "HMGB2", "KIF22", "KLF1", "CENPU",
  "RFC3", "RHCE", "TAL1", "TCF3", "TCF7L1", "ALAS2", "TSPO2", "CDCA3", "CENPA", "CKLF",
  "GTSE1", "NFE2", "OIP5", "SHCBP1", "SNF8", "SNRPB", "SPTA1", "KEL", "KIF2C", "LYL1",
  "PSMD9", "RACGAP1", "RPIA", "TIMELESS", "TRMT5", "TYMS", "VRK1", "FBXO7", "H3-3A", "NUDT1",
  "PLEK", "POLE2", "RHAG", "NSD2", "APOBEC3B", "BPGM", "CHAF1A", "DNAJC9", "FBXO5", "KIF4A",
  "KLF15", "MKI67", "PGD", "RPA3", "SRSF2", "TRIM58", "ADAMTS13", "CDCA4", "CDKN3", "FECH",
  "LBR", "LIG1", "LMNB1", "MICB", "NUP210", "ORC6", "RAD51AP1", "RHD", "SMC4", "GYPA",
  "RPP30", "TRIM10"
)

essential_cell_cycle_genes <- c(
  "ANAPC2", "ATP6V1B2", "ATR", "AURKA", "BIRC5", "BUB3", "CCDC101", "CDC16", "CDC23", "CDC45",
  "CDC5L", "CDC6", "CDK1", "CDT1", "CHAF1A", "CHEK1", "DAD1", "DCTN2", "DDB1", "DDX54",
  "DGCR14", "DHX30", "DONSON", "DTL", "EEF2", "EIF3C", "EP300", "ESPL1", "FAU", "GINS2",
  "HINFP", "INTS3", "KIF11", "KIF23", "NDC80", "NELFB", "NOL9", "NOP9", "NSF", "ORAOV1",
  "ORC3", "PCID2", "PLK1", "POLA1", "POLE2", "PPP1CB", "PPP2CA", "PRMT1", "PRPF19", "PSMA3",
  "PSMB7", "RACGAP1", "RFC2", "RPA1", "RPS13", "RPS2", "RPS6", "RRM1", "RRM2", "SEC61A1",
  "SLC7A6OS", "SMARCA4", "SRSF3", "TARS", "TICRR", "TIMELESS", "TONSL", "TOP2A", "TPX2", "U2AF2",
  "UPF1", "WDR74", "WEE1", "YME1L1", "ZNF830"
)

TERM2GENE <- data.frame(
  term = c(rep("Proliferation_Science_paper", length(proliferation_science_genes)),
           rep("Meta_PCNA", length(meta_pcna_genes)),
           rep("Essential_cell_cycle_Viner_Breuer_Cell_Division_2019", length(essential_cell_cycle_genes))),
  Gene_symbol = c(proliferation_science_genes, meta_pcna_genes, essential_cell_cycle_genes)
)

TERM2GENE <- TERM2GENE %>%
  mutate(term = recode(term,
                       "Proliferation_Science_paper" = "Signt_1",
                       "Meta_PCNA" = "Signt_2",
                       "Essential_cell_cycle_Viner_Breuer_Cell_Division_2019" = "Signt_3"))
print(TERM2GENE)

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host="https://useast.ensembl.org")
gene_ids <- (TERM2GENE$Gene_symbol)

res_gene_entrez_id <- getBM(attributes = c('entrezgene_id',
                                           'external_gene_name'),
                            filters = 'external_gene_name',
                            values = gene_ids,
                            mart = mart)
res_gene_entrez_id <- dplyr::filter(res_gene_entrez_id, !is.na(entrezgene_id)) %>%
  dplyr::rename(Gene_symbol = external_gene_name) %>%
  distinct(Gene_symbol, .keep_all = TRUE)

TERM2GENE_entrez <- dplyr::left_join(TERM2GENE, res_gene_entrez_id, by = "Gene_symbol") %>%
  dplyr::filter(!is.na(entrezgene_id)) %>%
  dplyr::select(-Gene_symbol)
  TERM2GENE <- TERM2GENE_entrez
  print(TERM2GENE)

  process_and_get_gsea_results <- function(data_frame, log_fc_cols, TERM2GENE) {
    gsea_results_list <- list()

    for (i in seq_along(log_fc_cols)) {
      log_fc_col <- log_fc_cols[i]
      plot_title <- plot_titles[i]

      all_values_entrez <- data_frame %>%
        dplyr::select(entrezgene_id, !!sym(log_fc_col)) %>%
        dplyr::filter(!is.na(!!sym(log_fc_col))) %>%
        dplyr::filter(!is.na(entrezgene_id))

      all_values_Log2_FC_entrez <- all_values_entrez[[log_fc_col]]
      names(all_values_Log2_FC_entrez) <- all_values_entrez$entrezgene_id

      all_values_Log2_FC_entrez <- all_values_Log2_FC_entrez[!duplicated(names(all_values_Log2_FC_entrez))]

      all_values_Log2_FC_entrez <- na.omit(all_values_Log2_FC_entrez)
      all_values_Log2_FC_entrez <- sort(all_values_Log2_FC_entrez, decreasing = TRUE)

      cat("Sample:", plot_titles[i], "\n")
      cat("Number of genes:", length(all_values_Log2_FC_entrez), "\n")

      if (length(all_values_Log2_FC_entrez) == 0) {
        cat("No genes available for", plot_title, "\n")
        next
      }

      gsea_res <- GSEA(geneList = all_values_Log2_FC_entrez,
                       TERM2GENE = TERM2GENE,
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "BH",
                       minGSSize = 10,
                       maxGSSize = 1000)

      if (is.null(gsea_res) || nrow(gsea_res) == 0) {
        cat("No significant pathways found for", plot_title, "\n")
        next
      }

      gsea_results_list[[plot_title]] <- gsea_res
    }

    combined_gsea <- do.call(rbind, lapply(names(gsea_results_list), function(name) {
      df <- as.data.frame(gsea_results_list[[name]])
      df$Sample <- name
      return(df)
    }))

    return(combined_gsea)
  }

  data_frame <- Expression_changes_targeted_agents_genes
  log_fc_cols <- c(
    "GSE93124_SW480_crizotinib_Log_FC",
    "GSE178755_HCC827_erlotinib_Log_FC",
    "GSE137535_SH454_Log_FC",
    "GSE105083_PANC1_CM03_24h_Log_FC",
    "GSE118548_T84_Trametinib_Log_FC",
    "GSE148096_TNBC_erlotinib_Log_FC",
    "GSE143462_KYSE150_omipalisib_Log_FC",
    "GSE121105_BT474_trastuzumab_Log_FC",
    "GSE103021_LU65_ARS_1620_Log_FC",
    "GSE135514_MDAMB468_Das_24h_Log_FC",
    "GSE135514_MDAMB468_Das_72h_Log_FC",

    "GSE241735_ICC13_7_Futibatinib_4h_Log_FC",
    "GSE241735_ICC13_7_Futibatinib_12h_Log_FC",
    "GSE241735_ICC13_7_Futibatinib_24h_Log_FC",

    "GSE105083_MIA_PACA2_CM03_6h_Log_FC",
    "GSE105083_MIA_PACA2_CM03_24h_Log_FC"
  )

  plot_titles <- c(
    "GSE93124",
    "GSE178755",
    "GSE137535_SH454",
    "GSE105083_PANC1_CM03_24h",
    "GSE118548_tra",
    "GSE148096",
    "GSE143462",
    "GSE121105",
    "GSE103021",
    "GSE135514_24h",
    "GSE135514_72h",

    "GSE241735_Futibatinib_4h",
    "GSE241735_Futibatinib_12h",
    "GSE241735_Futibatinib_24h",

    "GSE105083_MIA_PACA2_6h",
    "GSE105083_MIA_PACA2_24h"
  )

  stopifnot(
    length(log_fc_cols) == length(plot_titles)
  )

  combined_gsea <- process_and_get_gsea_results(data_frame, log_fc_cols, TERM2GENE)
  write.csv(combined_gsea, file.path(figure_output_dir, "combined_gsea_results.csv"), row.names = FALSE)

combined_gsea_filtered <- dplyr::select(combined_gsea, Sample, ID, NES, p.adjust, Description)
rownames(combined_gsea_filtered) <- NULL
combined_gsea_seleted <- tidyr::expand_grid(
  Sample = Fig3B_sample_order,
  ID = c(
    "Signt_1",
    "Signt_2",
    "Signt_3"
  )
) %>%
  dplyr::left_join(
    combined_gsea_filtered,
    by = c(
      "Sample",
      "ID"
    )
  ) %>%
  dplyr::mutate(
    Description = dplyr::coalesce(
      Description,
      ID
    ),
    NES = tidyr::replace_na(
      NES,
      0
    ),
    p.adjust = tidyr::replace_na(
      p.adjust,
      1
    )
  )

Plot_GREIN_selected_samples_complete_plot <- dplyr::bind_rows(combined_gsea_seleted, Plot_GREIN_selected_samples_complete)

desired_order <- c("Signt_1", "Signt_2", "Signt_3", "Cell cycle", "DNA replication",
                   "Homologous recombination", "Mismatch repair", "Fanconi anemia pathway", "Base excision repair", "Nucleotide excision repair", "ATP-dependent chromatin remodeling", "p53 pathway",
                   "PPAR pathway", "Apelin pathway", "C-type lectin receptor pathway", "ErbB pathway", "Hippo pathway",
                   "IL-17 pathway", "Insulin pathway", "JAK-STAT pathway", "MAPK pathway",  "AMPK pathway", "mTOR pathway", "PI3K-Akt pathway", "Ras pathway",
                   "Rap1 pathway", "TGF-beta pathway", "TNF pathway", "Wnt pathway", "Notch pathway")

Plot_GREIN_selected_samples_complete_plot <- Plot_GREIN_selected_samples_complete_plot %>%
  mutate(Description = factor(Description, levels = desired_order),
         neg_log10_p_adjust = -log10(p.adjust))

Plot_GREIN_selected_samples_complete_plot <-
  Plot_GREIN_selected_samples_complete_plot %>%
  dplyr::mutate(
    Sample = factor(
      Sample,
      levels = Fig3B_sample_order
    )
  ) %>%
  dplyr::filter(
    !is.na(Sample),
    !is.na(Description)
  )

Plot_cell_proliferation <- ggplot(Plot_GREIN_selected_samples_complete_plot, aes(x = Sample, y = Description, fill = NES, size = neg_log10_p_adjust)) +
  geom_point(
    data = dplyr::filter(Plot_GREIN_selected_samples_complete_plot, !is.na(NES), NES != 0),
    shape = 21,
    color = "black",
    stroke = 0.5
  ) +
  scale_fill_gradient2(low = "dodgerblue4", mid = "white", high = "firebrick3", midpoint = 0, limits = c(-3, 3)) +
  scale_size_continuous(range = c(2, 10)) +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
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

Fig3B_size_max <- max(
  Plot_GREIN_selected_samples_complete_plot$neg_log10_p_adjust[
    is.finite(
      Plot_GREIN_selected_samples_complete_plot$
        neg_log10_p_adjust
    ) &
      Plot_GREIN_selected_samples_complete_plot$NES != 0
  ],
  na.rm = TRUE
)

Plot_cell_proliferation_1_data <-
  Plot_GREIN_selected_samples_complete_plot %>%
  dplyr::filter(
    Sample %in% c(
      "GSE143462",
      "GSE121105",
      "GSE103021"
    )
  ) %>%
  dplyr::mutate(
    Sample = factor(
      as.character(Sample),
      levels = c(
        "GSE143462",
        "GSE121105",
        "GSE103021"
      )
    )
  )

Plot_cell_proliferation_1 <- ggplot(
  Plot_cell_proliferation_1_data,
  aes(
    x = Sample,
    y = Description,
    fill = NES,
    size = neg_log10_p_adjust
  )
) +
  geom_point(
    data = dplyr::filter(
      Plot_cell_proliferation_1_data,
      !is.na(NES),
      NES != 0
    ),
    shape = 21,
    color = "black",
    stroke = 0.5
  ) +
  scale_fill_gradient2(
    low = "dodgerblue4",
    mid = "white",
    high = "firebrick3",
    midpoint = 0,
    limits = c(-3, 3)
  ) +
  scale_size_continuous(
    range = c(2, 10),
    limits = c(0, Fig3B_size_max)
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  theme_bw() +
  labs(
    title = "",
    x = "",
    y = "",
    fill = "NES",
    size = "-log10(p.adjust)"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 16
    ),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 10)
  )

print(Plot_cell_proliferation_1)

Plot_cell_proliferation_2_data <-
  Plot_GREIN_selected_samples_complete_plot %>%
  dplyr::filter(
    Sample %in% c(
      "GSE135514_24h",
      "GSE135514_72h"
    )
  ) %>%
  dplyr::mutate(
    Sample = factor(
      as.character(Sample),
      levels = c(
        "GSE135514_24h",
        "GSE135514_72h"
      )
    )
  )

Plot_cell_proliferation_2 <- ggplot(
  Plot_cell_proliferation_2_data,
  aes(
    x = Sample,
    y = Description,
    fill = NES,
    size = neg_log10_p_adjust
  )
) +
  geom_point(
    data = dplyr::filter(
      Plot_cell_proliferation_2_data,
      !is.na(NES),
      NES != 0
    ),
    shape = 21,
    color = "black",
    stroke = 0.5
  ) +
  scale_fill_gradient2(
    low = "dodgerblue4",
    mid = "white",
    high = "firebrick3",
    midpoint = 0,
    limits = c(-3, 3)
  ) +
  scale_size_continuous(
    range = c(2, 10),
    limits = c(0, Fig3B_size_max)
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  theme_bw() +
  labs(
    title = "",
    x = "",
    y = "",
    fill = "NES",
    size = "-log10(p.adjust)"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 16
    ),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 10)
  )

print(Plot_cell_proliferation_2)

Plot_cell_proliferation_3_data <-
  Plot_GREIN_selected_samples_complete_plot %>%
  dplyr::filter(
    Sample %in% c(
      "GSE241735_Futibatinib_4h",
      "GSE241735_Futibatinib_12h",
      "GSE241735_Futibatinib_24h"
    )
  ) %>%
  dplyr::mutate(
    Sample = factor(
      as.character(Sample),
      levels = c(
        "GSE241735_Futibatinib_4h",
        "GSE241735_Futibatinib_12h",
        "GSE241735_Futibatinib_24h"
      )
    )
  )

Plot_cell_proliferation_3 <- ggplot(
  Plot_cell_proliferation_3_data,
  aes(
    x = Sample,
    y = Description,
    fill = NES,
    size = neg_log10_p_adjust
  )
) +
  geom_point(
    data = dplyr::filter(
      Plot_cell_proliferation_3_data,
      !is.na(NES),
      NES != 0
    ),
    shape = 21,
    color = "black",
    stroke = 0.5
  ) +
  scale_fill_gradient2(
    low = "dodgerblue4",
    mid = "white",
    high = "firebrick3",
    midpoint = 0,
    limits = c(-3, 3)
  ) +
  scale_size_continuous(
    range = c(2, 10),
    limits = c(0, Fig3B_size_max)
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  theme_bw() +
  labs(
    title = "",
    x = "",
    y = "",
    fill = "NES",
    size = "-log10(p.adjust)"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 16
    ),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 10)
  )

print(Plot_cell_proliferation_3)

ggsave(
  filename = "Plot_cell_proliferation_single_timepoints.tiff",
  plot = Plot_cell_proliferation_1,
  width = 7,
  height = 12,
  units = "in",
  dpi = 600,
  device = "tiff",
  compression = "lzw",
  bg = "white",
  path = figure_output_dir
)

ggsave(
  filename = "Plot_cell_proliferation_GSE135514_timecourse.tiff",
  plot = Plot_cell_proliferation_2,
  width = 6.5,
  height = 12,
  units = "in",
  dpi = 600,
  device = "tiff",
  compression = "lzw",
  bg = "white",
  path = figure_output_dir
)

ggsave(
  filename = "Plot_cell_proliferation_GSE241735_timecourse.tiff",
  plot = Plot_cell_proliferation_3,
  width = 6.5,
  height = 12,
  units = "in",
  dpi = 600,
  device = "tiff",
  compression = "lzw",
  bg = "white",
  path = figure_output_dir
)

Plot_MFE_GSE103021_ARS_1620_LU65 <- ggplot(
  MFE_res_GSE103021_ARS_1620_LU65,
  aes(
    x = factor(log2FoldChange > 0),
    y = Normalized_MFE
  )
) +
  geom_violin(
    fill = "white",
    size = 1.5
  ) +
  theme_classic() +
  ylim(-0.6, 0) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 20),
    plot.title = element_text(size = 10),
    axis.text.x = element_text(
      size = 16,
      color = "black"
    )
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.3
  ) +
  labs(
    y = "Normalized_MFE",
    x = NULL,
    title = "GSE103021_ARS_1620_LU65_Log2_FC_adj.P<0.05"
  ) +
  scale_x_discrete(
    labels = c(
      "log2_FC < 0",
      "log2_FC > 0"
    )
  ) +
  geom_violin(
    fill = "lightgrey",
    color = "black",
    size = 0.1
  ) +
  geom_boxplot(
    width = 0.1,
    color = "black",
    outlier.size = 0
  ) +
  theme(
    plot.margin = unit(
      c(4, 0.2, 0.2, 0.2),
      "cm"
    ),
    axis.title = element_blank(),
    axis.text.x = element_text(size = 14),
    plot.title = element_blank()
  )

Plot_MFE_GSE121105 <- ggplot(
  MFE_res_GSE121105_trastuzumab,
  aes(
    x = factor(log2FoldChange > 0),
    y = Normalized_MFE
  )
) +
  geom_violin(
    fill = "white",
    size = 1.5
  ) +
  theme_classic() +
  ylim(-0.6, 0) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 20),
    plot.title = element_text(size = 10),
    axis.text.x = element_text(
      size = 16,
      color = "black"
    )
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.3
  ) +
  labs(
    y = "Normalized_MFE",
    x = NULL,
    title = "GSE121105_trastuzumab_Log2_FC_adj.P<0.05"
  ) +
  scale_x_discrete(
    labels = c(
      "log2_FC < 0",
      "log2_FC > 0"
    )
  ) +
  geom_violin(
    fill = "lightgrey",
    color = "black",
    size = 0.1
  ) +
  geom_boxplot(
    width = 0.1,
    color = "black",
    outlier.size = 0
  ) +
  theme(
    plot.margin = unit(
      c(4, 0.2, 0.2, 0.2),
      "cm"
    ),
    axis.title = element_blank(),
    axis.text.x = element_text(size = 14),
    plot.title = element_blank()
  )

Plot_MFE_GSE143462 <- ggplot(
  MFE_res_GSE143462,
  aes(
    x = factor(log2FoldChange > 0),
    y = Normalized_MFE
  )
) +
  geom_violin(
    fill = "white",
    size = 1.5
  ) +
  theme_classic() +
  ylim(-0.6, 0) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 20),
    plot.title = element_text(size = 10),
    axis.text.x = element_text(
      size = 16,
      color = "black"
    )
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.3
  ) +
  labs(
    y = "Normalized_MFE",
    x = NULL,
    title = "GSE143462_Log2_FC_adj.P<0.05"
  ) +
  scale_x_discrete(
    labels = c(
      "log2_FC < 0",
      "log2_FC > 0"
    )
  ) +
  geom_violin(
    fill = "lightgrey",
    color = "black",
    size = 0.1
  ) +
  geom_boxplot(
    width = 0.1,
    color = "black",
    outlier.size = 0
  ) +
  theme(
    plot.margin = unit(
      c(4, 0.2, 0.2, 0.2),
      "cm"
    ),
    axis.title = element_blank(),
    axis.text.x = element_text(size = 14),
    plot.title = element_blank()
  )

Fig3C_panel_1 <- gridExtra::arrangeGrob(
  Plot_MFE_GSE103021_ARS_1620_LU65,
  Plot_MFE_GSE121105,
  Plot_MFE_GSE143462,
  ncol = 3,
  nrow = 1
)

grid::grid.newpage()
grid::grid.draw(Fig3C_panel_1)

Plot_MFE_GSE135514_24h <- ggplot(
  MFE_res_GSE135514_Das_24h,
  aes(
    x = factor(log2FoldChange > 0),
    y = Normalized_MFE
  )
) +
  geom_violin(
    fill = "white",
    size = 1.5
  ) +
  theme_classic() +
  ylim(-0.6, 0) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 20),
    plot.title = element_text(size = 10),
    axis.text.x = element_text(
      size = 16,
      color = "black"
    )
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.3
  ) +
  labs(
    y = "Normalized_MFE",
    x = NULL,
    title = "GSE135514_Das_24h_Log2_FC_adj.P<0.05"
  ) +
  scale_x_discrete(
    labels = c(
      "log2_FC < 0",
      "log2_FC > 0"
    )
  ) +
  geom_violin(
    fill = "lightgrey",
    color = "black",
    size = 0.1
  ) +
  geom_boxplot(
    width = 0.1,
    color = "black",
    outlier.size = 0
  ) +
  theme(
    plot.margin = unit(
      c(4, 0.2, 0.2, 0.2),
      "cm"
    ),
    axis.title = element_blank(),
    axis.text.x = element_text(size = 14),
    plot.title = element_blank()
  )

Plot_MFE_GSE135514_72h <- ggplot(
  MFE_res_GSE135514_Das_72h,
  aes(
    x = factor(log2FoldChange > 0),
    y = Normalized_MFE
  )
) +
  geom_violin(
    fill = "white",
    size = 1.5
  ) +
  theme_classic() +
  ylim(-0.6, 0) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 20),
    plot.title = element_text(size = 10),
    axis.text.x = element_text(
      size = 16,
      color = "black"
    )
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.3
  ) +
  labs(
    y = "Normalized_MFE",
    x = NULL,
    title = "GSE135514_Das_72h_Log2_FC_adj.P<0.05"
  ) +
  scale_x_discrete(
    labels = c(
      "log2_FC < 0",
      "log2_FC > 0"
    )
  ) +
  geom_violin(
    fill = "lightgrey",
    color = "black",
    size = 0.1
  ) +
  geom_boxplot(
    width = 0.1,
    color = "black",
    outlier.size = 0
  ) +
  theme(
    plot.margin = unit(
      c(4, 0.2, 0.2, 0.2),
      "cm"
    ),
    axis.title = element_blank(),
    axis.text.x = element_text(size = 14),
    plot.title = element_blank()
  )

Fig3C_panel_2 <- gridExtra::arrangeGrob(
  Plot_MFE_GSE135514_24h,
  Plot_MFE_GSE135514_72h,
  ncol = 2,
  nrow = 1
)

grid::grid.newpage()
grid::grid.draw(Fig3C_panel_2)

Plot_MFE_GSE241735_Futibatinib_4h <- ggplot(
  MFE_res_GSE241735_Futibatinib_4h,
  aes(
    x = factor(log2FoldChange > 0),
    y = Normalized_MFE
  )
) +
  geom_violin(
    fill = "white",
    size = 1.5
  ) +
  theme_classic() +
  ylim(-0.6, 0) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 20),
    plot.title = element_text(size = 10),
    axis.text.x = element_text(
      size = 16,
      color = "black"
    )
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.3
  ) +
  labs(
    y = "Normalized_MFE",
    x = NULL,
    title = "GSE241735_Futibatinib_4h_Log2_FC_adj.P<0.05"
  ) +
  scale_x_discrete(
    labels = c(
      "log2_FC < 0",
      "log2_FC > 0"
    )
  ) +
  geom_violin(
    fill = "lightgrey",
    color = "black",
    size = 0.1
  ) +
  geom_boxplot(
    width = 0.1,
    color = "black",
    outlier.size = 0
  ) +
  theme(
    plot.margin = unit(
      c(4, 0.2, 0.2, 0.2),
      "cm"
    ),
    axis.title = element_blank(),
    axis.text.x = element_text(size = 14),
    plot.title = element_blank()
  )

Plot_MFE_GSE241735_Futibatinib_12h <- ggplot(
  MFE_res_GSE241735_Futibatinib_12h,
  aes(
    x = factor(log2FoldChange > 0),
    y = Normalized_MFE
  )
) +
  geom_violin(
    fill = "white",
    size = 1.5
  ) +
  theme_classic() +
  ylim(-0.6, 0) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 20),
    plot.title = element_text(size = 10),
    axis.text.x = element_text(
      size = 16,
      color = "black"
    )
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.3
  ) +
  labs(
    y = "Normalized_MFE",
    x = NULL,
    title = "GSE241735_Futibatinib_12h_Log2_FC_adj.P<0.05"
  ) +
  scale_x_discrete(
    labels = c(
      "log2_FC < 0",
      "log2_FC > 0"
    )
  ) +
  geom_violin(
    fill = "lightgrey",
    color = "black",
    size = 0.1
  ) +
  geom_boxplot(
    width = 0.1,
    color = "black",
    outlier.size = 0
  ) +
  theme(
    plot.margin = unit(
      c(4, 0.2, 0.2, 0.2),
      "cm"
    ),
    axis.title = element_blank(),
    axis.text.x = element_text(size = 14),
    plot.title = element_blank()
  )

Plot_MFE_GSE241735_Futibatinib_24h <- ggplot(
  MFE_res_GSE241735_Futibatinib_24h,
  aes(
    x = factor(log2FoldChange > 0),
    y = Normalized_MFE
  )
) +
  geom_violin(
    fill = "white",
    size = 1.5
  ) +
  theme_classic() +
  ylim(-0.6, 0) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 20),
    plot.title = element_text(size = 10),
    axis.text.x = element_text(
      size = 16,
      color = "black"
    )
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.3
  ) +
  labs(
    y = "Normalized_MFE",
    x = NULL,
    title = "GSE241735_Futibatinib_24h_Log2_FC_adj.P<0.05"
  ) +
  scale_x_discrete(
    labels = c(
      "log2_FC < 0",
      "log2_FC > 0"
    )
  ) +
  geom_violin(
    fill = "lightgrey",
    color = "black",
    size = 0.1
  ) +
  geom_boxplot(
    width = 0.1,
    color = "black",
    outlier.size = 0
  ) +
  theme(
    plot.margin = unit(
      c(4, 0.2, 0.2, 0.2),
      "cm"
    ),
    axis.title = element_blank(),
    axis.text.x = element_text(size = 14),
    plot.title = element_blank()
  )

Fig3C_panel_3 <- gridExtra::arrangeGrob(
  Plot_MFE_GSE241735_Futibatinib_4h,
  Plot_MFE_GSE241735_Futibatinib_12h,
  Plot_MFE_GSE241735_Futibatinib_24h,
  ncol = 3,
  nrow = 1
)

grid::grid.newpage()
grid::grid.draw(Fig3C_panel_3)

ggsave(
  filename = "Fig3C_panel_1_single_timepoint_datasets.tiff",
  plot = Fig3C_panel_1,
  width = 12,
  height = 7,
  units = "in",
  dpi = 600,
  device = "tiff",
  compression = "lzw",
  bg = "white",
  path = figure_output_dir
)

ggsave(
  filename = "Fig3C_panel_2_GSE135514_timecourse.tiff",
  plot = Fig3C_panel_2,
  width = 8,
  height = 7,
  units = "in",
  dpi = 600,
  device = "tiff",
  compression = "lzw",
  bg = "white",
  path = figure_output_dir
)

ggsave(
  filename = "Fig3C_panel_3_GSE241735_timecourse.tiff",
  plot = Fig3C_panel_3,
  width = 12,
  height = 7,
  units = "in",
  dpi = 600,
  device = "tiff",
  compression = "lzw",
  bg = "white",
  path = figure_output_dir
)

format_P_value <- function(p) {

  label <- rep(NA_character_, length(p))

  below_threshold <- !is.na(p) & p < 0.0001
  exact_value <- !is.na(p) & p >= 0.0001

  label[below_threshold] <- "P < 0.0001"

  label[exact_value] <- paste0(
    "P = ",
    trimws(
      base::format.pval(
        p[exact_value],
        digits = 4
      )
    )
  )

  label
}

results <- data.frame(

  Dataset = c(
    "GSE103021_ARS_1620_LU65",
    "GSE121105_trastuzumab",
    "GSE143462",
    "GSE135514_Das_24h",
    "GSE135514_Das_72h",
    "GSE241735_Futibatinib_4h",
    "GSE241735_Futibatinib_12h",
    "GSE241735_Futibatinib_24h"
  ),

  Mean_MFE_UP = c(
    mean(
      MFE_res_GSE103021_ARS_1620_LU65_UP$Normalized_MFE,
      na.rm = TRUE
    ),
    mean(
      MFE_res_GSE121105_trastuzumab_UP$Normalized_MFE,
      na.rm = TRUE
    ),
    mean(
      MFE_res_GSE143462_UP$Normalized_MFE,
      na.rm = TRUE
    ),
    mean(
      MFE_res_GSE135514_Das_24h_UP$Normalized_MFE,
      na.rm = TRUE
    ),
    mean(
      MFE_res_GSE135514_Das_72h_UP$Normalized_MFE,
      na.rm = TRUE
    ),
    mean(
      MFE_res_GSE241735_Futibatinib_4h_UP$Normalized_MFE,
      na.rm = TRUE
    ),
    mean(
      MFE_res_GSE241735_Futibatinib_12h_UP$Normalized_MFE,
      na.rm = TRUE
    ),
    mean(
      MFE_res_GSE241735_Futibatinib_24h_UP$Normalized_MFE,
      na.rm = TRUE
    )
  ),

  Mean_MFE_DOWN = c(
    mean(
      MFE_res_GSE103021_ARS_1620_LU65_DOWN$Normalized_MFE,
      na.rm = TRUE
    ),
    mean(
      MFE_res_GSE121105_trastuzumab_DOWN$Normalized_MFE,
      na.rm = TRUE
    ),
    mean(
      MFE_res_GSE143462_DOWN$Normalized_MFE,
      na.rm = TRUE
    ),
    mean(
      MFE_res_GSE135514_Das_24h_DOWN$Normalized_MFE,
      na.rm = TRUE
    ),
    mean(
      MFE_res_GSE135514_Das_72h_DOWN$Normalized_MFE,
      na.rm = TRUE
    ),
    mean(
      MFE_res_GSE241735_Futibatinib_4h_DOWN$Normalized_MFE,
      na.rm = TRUE
    ),
    mean(
      MFE_res_GSE241735_Futibatinib_12h_DOWN$Normalized_MFE,
      na.rm = TRUE
    ),
    mean(
      MFE_res_GSE241735_Futibatinib_24h_DOWN$Normalized_MFE,
      na.rm = TRUE
    )
  ),

  KS_Test_p_value = c(
    ks_test_res_GSE103021_ARS_1620_LU65_all_UP_vs_DOWN$p.value,
    ks_test_res_GSE121105_trastuzumab_all_UP_vs_DOWN$p.value,
    ks_test_res_GSE143462_all_UP_vs_DOWN$p.value,
    ks_test_res_GSE135514_Das_24h_all_UP_vs_DOWN$p.value,
    ks_test_res_GSE135514_Das_72h_all_UP_vs_DOWN$p.value,
    ks_test_res_GSE241735_Futibatinib_4h_all_UP_vs_DOWN$p.value,
    ks_test_res_GSE241735_Futibatinib_12h_all_UP_vs_DOWN$p.value,
    ks_test_res_GSE241735_Futibatinib_24h_all_UP_vs_DOWN$p.value
  ),

  Wilcoxon_p_value = c(
    stats::wilcox.test(
      MFE_res_GSE103021_ARS_1620_LU65_DOWN$Normalized_MFE,
      MFE_res_GSE103021_ARS_1620_LU65_UP$Normalized_MFE,
      paired = FALSE,
      exact = FALSE
    )$p.value,

    stats::wilcox.test(
      MFE_res_GSE121105_trastuzumab_DOWN$Normalized_MFE,
      MFE_res_GSE121105_trastuzumab_UP$Normalized_MFE,
      paired = FALSE,
      exact = FALSE
    )$p.value,

    stats::wilcox.test(
      MFE_res_GSE143462_DOWN$Normalized_MFE,
      MFE_res_GSE143462_UP$Normalized_MFE,
      paired = FALSE,
      exact = FALSE
    )$p.value,

    stats::wilcox.test(
      MFE_res_GSE135514_Das_24h_DOWN$Normalized_MFE,
      MFE_res_GSE135514_Das_24h_UP$Normalized_MFE,
      paired = FALSE,
      exact = FALSE
    )$p.value,

    stats::wilcox.test(
      MFE_res_GSE135514_Das_72h_DOWN$Normalized_MFE,
      MFE_res_GSE135514_Das_72h_UP$Normalized_MFE,
      paired = FALSE,
      exact = FALSE
    )$p.value,

    stats::wilcox.test(
      MFE_res_GSE241735_Futibatinib_4h_DOWN$Normalized_MFE,
      MFE_res_GSE241735_Futibatinib_4h_UP$Normalized_MFE,
      paired = FALSE,
      exact = FALSE
    )$p.value,

    stats::wilcox.test(
      MFE_res_GSE241735_Futibatinib_12h_DOWN$Normalized_MFE,
      MFE_res_GSE241735_Futibatinib_12h_UP$Normalized_MFE,
      paired = FALSE,
      exact = FALSE
    )$p.value,

    stats::wilcox.test(
      MFE_res_GSE241735_Futibatinib_24h_DOWN$Normalized_MFE,
      MFE_res_GSE241735_Futibatinib_24h_UP$Normalized_MFE,
      paired = FALSE,
      exact = FALSE
    )$p.value
  )
)

results_report <- results %>%
  dplyr::mutate(
    KS_Test_P = format_P_value(KS_Test_p_value),
    Wilcoxon_P = format_P_value(Wilcoxon_p_value)
  ) %>%
  dplyr::select(
    Dataset,
    Mean_MFE_UP,
    Mean_MFE_DOWN,
    KS_Test_P,
    Wilcoxon_P
  )

print(results_report, row.names = FALSE)

library(ggplot2)
library(tidyr)

sample_order <- c(
  "GSE103021",
  "GSE121105",
  "GSE143462",
  "Spacer_1",
  "GSE135514_24h",
  "GSE135514_72h",
  "Spacer_2",
  "GSE241735_4h",
  "GSE241735_12h",
  "GSE241735_24h"
)

sample_labels <- c(
  "GSE103021" = "GSE103021",
  "GSE121105" = "GSE121105",
  "GSE143462" = "GSE143462",
  "Spacer_1" = "",
  "GSE135514_24h" = "GSE135514_24h",
  "GSE135514_72h" = "GSE135514_72h",
  "Spacer_2" = "",
  "GSE241735_4h" = "GSE241735_4h",
  "GSE241735_12h" = "GSE241735_12h",
  "GSE241735_24h" = "GSE241735_24h"
)

Data_high_MFE_fraction_total_counts <- data.frame(
  Samples = c(
    "GSE103021",
    "GSE121105",
    "GSE143462",
    "GSE135514_24h",
    "GSE135514_72h",
    "GSE241735_4h",
    "GSE241735_12h",
    "GSE241735_24h"
  ),

  High_MFE_of_total_control = c(
    GSE103021_ARS_1620_LU65_control_high_MFE_TX_counts_fraction,
    GSE121105_trastuzumab_control_high_MFE_TX_counts_fraction,
    GSE143462_control_high_MFE_TX_counts_fraction,
    GSE135514_Das_24h_control_high_MFE_TX_counts_fraction,
    GSE135514_Das_72h_control_high_MFE_TX_counts_fraction,
    GSE241735_Futibatinib_4h_control_high_MFE_TX_counts_fraction,
    GSE241735_Futibatinib_12h_control_high_MFE_TX_counts_fraction,
    GSE241735_Futibatinib_24h_control_high_MFE_TX_counts_fraction
  ),

  High_MFE_of_total_drug_treated = c(
    GSE103021_ARS_1620_LU65_ARS_1620_LU65_high_MFE_TX_counts_fraction,
    GSE121105_trastuzumab_trastuzumab_high_MFE_TX_counts_fraction,
    GSE143462_Omipalisib_high_MFE_TX_counts_fraction,
    GSE135514_Das_24h_Das_24h_high_MFE_TX_counts_fraction,
    GSE135514_Das_72h_Das_72h_high_MFE_TX_counts_fraction,
    GSE241735_Futibatinib_4h_Futibatinib_4h_high_MFE_TX_counts_fraction,
    GSE241735_Futibatinib_12h_Futibatinib_12h_high_MFE_TX_counts_fraction,
    GSE241735_Futibatinib_24h_Futibatinib_24h_high_MFE_TX_counts_fraction
  )
)

print(Data_high_MFE_fraction_total_counts)

Data_high_MFE_fraction_total_counts$Samples <- factor(
  Data_high_MFE_fraction_total_counts$Samples,
  levels = sample_order
)

df_long_high <- pivot_longer(
  Data_high_MFE_fraction_total_counts,
  cols = c(
    High_MFE_of_total_control,
    High_MFE_of_total_drug_treated
  ),
  names_to = "Condition",
  values_to = "Value"
)

Plot_high_fraction_tot_counts <- ggplot(
  df_long_high,
  aes(
    x = Samples,
    y = Value,
    fill = Condition
  )
) +
  geom_bar(
    stat = "identity",
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  geom_text(
    aes(label = round(Value, 3)),
    position = position_dodge(width = 0.8),
    hjust = -0.3,
    vjust = 0.4,
    angle = 90,
    size = 4
  ) +
  labs(
    title = "",
    x = "",
    y = "Fraction of total counts",
    fill = ""
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 12
    ),
    axis.title.y = element_text(size = 14),
    legend.position = "none"
  ) +
  ylim(0, 0.4) +
  scale_x_discrete(
    drop = FALSE,
    labels = sample_labels
  ) +
  scale_fill_manual(
    values = c(
      "High_MFE_of_total_control" = "#8cc5e3",
      "High_MFE_of_total_drug_treated" = "#2066a8"
    ),
    labels = c(
      "High_MFE_of_total_control" = "Untreated",
      "High_MFE_of_total_drug_treated" = "Drug_treated"
    )
  )

print(Plot_high_fraction_tot_counts)

ggsave(
  filename = "Plot_high_fraction_tot_counts.tiff",
  plot = Plot_high_fraction_tot_counts,
  width = 5,
  height = 4.5,
  units = "in",
  dpi = 600,
  device = "tiff",
  compression = "lzw",
  bg = "white",
  path = figure_output_dir
)

Data_low_MFE_fraction_total_counts <- data.frame(
  Samples = c(
    "GSE103021",
    "GSE121105",
    "GSE143462",
    "GSE135514_24h",
    "GSE135514_72h",
    "GSE241735_4h",
    "GSE241735_12h",
    "GSE241735_24h"
  ),

  low_MFE_of_total_control = c(
    GSE103021_ARS_1620_LU65_control_low_MFE_TX_counts_fraction,
    GSE121105_trastuzumab_control_low_MFE_TX_counts_fraction,
    GSE143462_control_low_MFE_TX_counts_fraction,
    GSE135514_Das_24h_control_low_MFE_TX_counts_fraction,
    GSE135514_Das_72h_control_low_MFE_TX_counts_fraction,
    GSE241735_Futibatinib_4h_control_low_MFE_TX_counts_fraction,
    GSE241735_Futibatinib_12h_control_low_MFE_TX_counts_fraction,
    GSE241735_Futibatinib_24h_control_low_MFE_TX_counts_fraction
  ),

  low_MFE_of_total_drug_treated = c(
    GSE103021_ARS_1620_LU65_low_MFE_TX_counts_fraction,
    GSE121105_trastuzumab_trastuzumab_low_MFE_TX_counts_fraction,
    GSE143462_Omipalisib_low_MFE_TX_counts_fraction,
    GSE135514_Das_24h_low_MFE_TX_counts_fraction,
    GSE135514_Das_72h_low_MFE_TX_counts_fraction,
    GSE241735_Futibatinib_4h_Futibatinib_4h_low_MFE_TX_counts_fraction,
    GSE241735_Futibatinib_12h_Futibatinib_12h_low_MFE_TX_counts_fraction,
    GSE241735_Futibatinib_24h_Futibatinib_24h_low_MFE_TX_counts_fraction
  )
)

print(Data_low_MFE_fraction_total_counts)

Data_low_MFE_fraction_total_counts$Samples <- factor(
  Data_low_MFE_fraction_total_counts$Samples,
  levels = sample_order
)

df_long_low <- pivot_longer(
  Data_low_MFE_fraction_total_counts,
  cols = c(
    low_MFE_of_total_control,
    low_MFE_of_total_drug_treated
  ),
  names_to = "Condition",
  values_to = "Value"
)

Plot_low_fraction_tot_counts <- ggplot(
  df_long_low,
  aes(
    x = Samples,
    y = Value,
    fill = Condition
  )
) +
  geom_bar(
    stat = "identity",
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  geom_text(
    aes(label = round(Value, 3)),
    position = position_dodge(width = 0.8),
    hjust = -0.3,
    vjust = 0.4,
    angle = 90,
    size = 4
  ) +
  labs(
    title = "",
    x = "",
    y = "Fraction of total counts",
    fill = ""
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 12
    ),
    axis.title.y = element_text(size = 14),
    legend.position = "none"
  ) +
  ylim(0, 0.4) +
  scale_x_discrete(
    drop = FALSE,
    labels = sample_labels
  ) +
  scale_fill_manual(
    values = c(
      "low_MFE_of_total_control" = "#8cc5e3",
      "low_MFE_of_total_drug_treated" = "#2066a8"
    ),
    labels = c(
      "low_MFE_of_total_control" = "Untreated",
      "low_MFE_of_total_drug_treated" = "Drug_treated"
    )
  )

print(Plot_low_fraction_tot_counts)

ggsave(
  filename = "Plot_low_fraction_tot_counts.tiff",
  plot = Plot_low_fraction_tot_counts,
  width = 5,
  height = 4.5,
  units = "in",
  dpi = 600,
  device = "tiff",
  compression = "lzw",
  bg = "white",
  path = figure_output_dir
)

Fig3D_complete_stats <- Data_high_MFE_fraction_total_counts %>%
  dplyr::select(
    Samples,
    High_MFE_control = High_MFE_of_total_control,
    High_MFE_drug_treated = High_MFE_of_total_drug_treated
  ) %>%
  dplyr::left_join(
    Data_low_MFE_fraction_total_counts %>%
      dplyr::select(
        Samples,
        Low_MFE_control = low_MFE_of_total_control,
        Low_MFE_drug_treated = low_MFE_of_total_drug_treated
      ),
    by = "Samples"
  ) %>%
  dplyr::mutate(
    High_MFE_drug_minus_control =
      High_MFE_drug_treated - High_MFE_control,
    Low_MFE_drug_minus_control =
      Low_MFE_drug_treated - Low_MFE_control
  ) %>%
  dplyr::arrange(Samples)

print(Fig3D_complete_stats)

Fig3D_perm_B <- 10000L
Fig3D_perm_seed <- 20260918L

Fig3D_perm_spec <- data.frame(
  Samples = c(
    "GSE103021", "GSE121105", "GSE143462",
    "GSE135514_24h", "GSE135514_72h",
    "GSE241735_4h", "GSE241735_12h", "GSE241735_24h"
  ),
  Prefix = c(
    "GSE103021_ARS_1620_LU65", "GSE121105_trastuzumab", "GSE143462",
    "GSE135514_Das_24h", "GSE135514_Das_72h",
    "GSE241735_Futibatinib_4h", "GSE241735_Futibatinib_12h",
    "GSE241735_Futibatinib_24h"
  ),
  Treatment = c(
    "ARS_1620_LU65", "trastuzumab", "Omipalisib",
    "Das_24h", "Das_72h",
    "Futibatinib_4h", "Futibatinib_12h", "Futibatinib_24h"
  ),
  Upper_suffix = c(
    "_high_MFE_TX_total", "_high_MFE_TX_total", "_high_MFE_TX_total",
    "_high_MFE_TX", "_high_MFE_TX",
    "_high_MFE_TX_total", "_high_MFE_TX_total", "_high_MFE_TX_total"
  ),
  stringsAsFactors = FALSE
)

Fig3D_perm_required <- c(
  "Fig3D_complete_stats",
  paste0(Fig3D_perm_spec$Prefix, "_norm_TX_aggregated_data_avg_counts"),
  paste0(Fig3D_perm_spec$Prefix, Fig3D_perm_spec$Upper_suffix),
  paste0(Fig3D_perm_spec$Prefix, "_low_MFE_TX_total")
)

Fig3D_perm_missing <- Fig3D_perm_required[
  !vapply(
    Fig3D_perm_required, exists, logical(1),
    envir = .GlobalEnv, inherits = FALSE
  )
]

if (length(Fig3D_perm_missing)) {
  stop(
    "Missing Fig. 3D source objects: ",
    paste(Fig3D_perm_missing, collapse = ", ")
  )
}

Fig3D_prepare_permutation <- function(i) {

  spec <- Fig3D_perm_spec[i, ]
  fail <- function(msg) stop(spec$Samples, ": ", msg)

  fetch <- function(suffix) {
    get(
      paste0(spec$Prefix, suffix),
      envir = .GlobalEnv, inherits = FALSE
    )
  }

  dat <- fetch("_norm_TX_aggregated_data_avg_counts")
  id_col <- "ensembl_transcript_id_version"
  required <- c(id_col, "condition", "Normalized_MFE", "avg_counts")

  if (!all(required %in% names(dat))) {
    fail("Missing source columns.")
  }

  if (
    anyNA(dat[, required]) ||
    any(!is.finite(dat$avg_counts)) ||
    any(dat$avg_counts < 0) ||
    any(!is.finite(dat$Normalized_MFE)) ||
    !setequal(as.character(dat$condition), c("Control", spec$Treatment))
  ) {
    fail("Check source counts, folding values or condition names.")
  }

  ref <- dat[dat$condition == "Control", ]
  cmp <- dat[dat$condition == spec$Treatment, ]

  if (
    anyDuplicated(ref[[id_col]]) ||
    anyDuplicated(cmp[[id_col]]) ||
    !setequal(ref[[id_col]], cmp[[id_col]])
  ) {
    fail("Transcript IDs are duplicated or differ between conditions.")
  }

  cmp <- cmp[match(ref[[id_col]], cmp[[id_col]]), ]

  if (any(abs(ref$Normalized_MFE - cmp$Normalized_MFE) > 1e-12)) {
    fail("Folding values differ between conditions.")
  }

  upper <- fetch(spec$Upper_suffix)
  lower <- fetch("_low_MFE_TX_total")

  if (!(id_col %in% names(upper)) || !(id_col %in% names(lower))) {
    fail("Missing quartile transcript IDs.")
  }

  is_upper <- ref[[id_col]] %in% upper[[id_col]]
  is_lower <- ref[[id_col]] %in% lower[[id_col]]

  if (any(is_upper & is_lower) || !any(is_upper) || !any(is_lower)) {
    fail("Invalid quartile memberships.")
  }

  category <- ifelse(is_upper, 1L, ifelse(is_lower, 2L, 0L))

  counts <- cbind(ref$avg_counts, cmp$avg_counts)

  if (any(!is.finite(colSums(counts))) || any(colSums(counts) <= 0)) {
    fail("Invalid total counts.")
  }

  weights <- sweep(counts, 2, colSums(counts), "/")

  fractions <- vapply(1:2, function(k) {
    colSums(weights[category == k, , drop = FALSE])
  }, numeric(2))

  old <- Fig3D_complete_stats[
    as.character(Fig3D_complete_stats$Samples) == spec$Samples, ]

  if (nrow(old) != 1L) {
    fail("Expected one existing summary row.")
  }

  expected <- c(
    old$High_MFE_control, old$High_MFE_drug_treated,
    old$Low_MFE_control, old$Low_MFE_drug_treated
  )

  if (
    length(expected) != 4L ||
    any(!is.finite(expected)) ||
    max(abs(as.vector(fractions) - expected)) > 1e-10
  ) {
    fail("Recalculated fractions do not match the existing bars.")
  }

  list(
    Samples = spec$Samples,
    category = category,
    fractions = fractions,
    delta = weights[, 2] - weights[, 1]
  )
}

Fig3D_perm_inputs <- lapply(
  seq_len(nrow(Fig3D_perm_spec)),
  Fig3D_prepare_permutation
)

message("All Fig. 3D fractions match the existing bars.")

set.seed(Fig3D_perm_seed)

Fig3D_perm_runs <- lapply(Fig3D_perm_inputs, function(x) {

  message("Permutations: ", x$Samples)

  observed <- x$fractions[2, ] - x$fractions[1, ]

  null <- replicate(Fig3D_perm_B, {

    shuffled <- x$category[sample.int(length(x$category))]
    c(
      sum(x$delta[shuffled == 1L]),
      sum(x$delta[shuffled == 2L])
    )
  })

  b <- rowSums(abs(null) >= abs(observed) - 1e-12)

  list(
    results = data.frame(
      Samples = x$Samples,
      Quartile = c("Upper / less negative", "Lower / more negative"),
      Fraction_control = x$fractions[1, ],
      Fraction_treated = x$fractions[2, ],
      Delta_fraction = observed,
      N_transcripts = length(x$category),
      N_quartile = c(sum(x$category == 1L), sum(x$category == 2L)),
      P_permutation = (b + 1) / (Fig3D_perm_B + 1),
      stringsAsFactors = FALSE
    ),
    null_delta = null
  )
})

Fig3D_quartile_permutation_results <- do.call(
  rbind,
  lapply(Fig3D_perm_runs, function(x) x$results)
)

stopifnot(nrow(Fig3D_quartile_permutation_results) == 16L)

Fig3D_quartile_permutation_results$P_adj_BH <- p.adjust(
  Fig3D_quartile_permutation_results$P_permutation,
  method = "BH"
)

Fig3D_quartile_permutation_results$Figure_label <- ifelse(
  Fig3D_quartile_permutation_results$P_adj_BH >= 0.05,
  "NS",
  sprintf("P = %.4f", Fig3D_quartile_permutation_results$P_adj_BH)
)

print(
  Fig3D_quartile_permutation_results[, c(
    "Samples", "Quartile", "Fraction_control", "Fraction_treated",
    "Delta_fraction", "P_permutation", "P_adj_BH", "Figure_label"
  )],
  row.names = FALSE,
  digits = 7
)

Fig3D_time_names <- c(
  "GSE135514_24h", "GSE135514_72h",
  "GSE241735_4h", "GSE241735_12h", "GSE241735_24h"
)

Fig3D_time_spec <- Fig3D_perm_spec[
  match(Fig3D_time_names, Fig3D_perm_spec$Samples),
  ,
  drop = FALSE
]

stopifnot(!anyNA(Fig3D_time_spec))

Fig3D_time_spec$Dataset <- c(
  "GSE135514", "GSE135514",
  "GSE241735", "GSE241735", "GSE241735"
)

Fig3D_time_spec$Time_h <- c(24, 72, 4, 12, 24)

Fig3D_time_quartiles <- c(
  "Upper / less negative",
  "Lower / more negative"
)

Fig3D_time_fraction_data <- dplyr::bind_rows(lapply(
  seq_len(nrow(Fig3D_time_spec)),
  function(i) {

    spec <- Fig3D_time_spec[i, ]

    fetch <- function(suffix) {
      object_name <- paste0(spec$Prefix, suffix)
      if (!exists(object_name, envir = .GlobalEnv, inherits = FALSE)) {
        stop("Missing source object: ", object_name)
      }
      get(object_name, envir = .GlobalEnv, inherits = FALSE)
    }

    dat <- as.data.frame(fetch("_long_data"))
    upper <- fetch(spec$Upper_suffix)
    lower <- fetch("_low_MFE_TX_total")
    id_col <- "ensembl_transcript_id_version"

    stopifnot(
      all(c(id_col, "replicate", "condition", "counts") %in% names(dat)),
      !anyNA(dat[, c(id_col, "replicate", "condition", "counts")]),
      !anyDuplicated(dat[, c(id_col, "replicate")]),
      all(is.finite(dat$counts)),
      all(dat$counts >= 0),
      setequal(unique(dat$condition), c("Control", spec$Treatment))
    )

    dat$is_upper <- dat[[id_col]] %in% upper[[id_col]]
    dat$is_lower <- dat[[id_col]] %in% lower[[id_col]]

    stopifnot(
      any(dat$is_upper),
      any(dat$is_lower),
      !any(dat$is_upper & dat$is_lower)
    )

    sample_counts <- dat %>%
      dplyr::group_by(replicate, condition) %>%
      dplyr::summarise(
        N_transcripts = dplyr::n(),
        Total_counts = sum(counts),
        Upper_counts = sum(counts[is_upper]),
        Lower_counts = sum(counts[is_lower]),
        .groups = "drop"
      )

    stopifnot(
      nrow(sample_counts) == 6L,
      !anyDuplicated(sample_counts$replicate),
      all(table(sample_counts$condition) == 3L),
      all(sample_counts$Total_counts > 0),
      all(
        sample_counts$N_transcripts ==
          length(unique(dat[[id_col]]))
      )
    )

    ref <- sample_counts$condition == "Control"

    reconstructed <- c(
      sum(sample_counts$Upper_counts[ref]) /
        sum(sample_counts$Total_counts[ref]),
      sum(sample_counts$Lower_counts[ref]) /
        sum(sample_counts$Total_counts[ref]),
      sum(sample_counts$Upper_counts[!ref]) /
        sum(sample_counts$Total_counts[!ref]),
      sum(sample_counts$Lower_counts[!ref]) /
        sum(sample_counts$Total_counts[!ref])
    )

    old <- Fig3D_quartile_permutation_results[
      Fig3D_quartile_permutation_results$Samples == spec$Samples,
      ,
      drop = FALSE
    ]

    stopifnot(nrow(old) == 2L)

    old <- old[
      match(Fig3D_time_quartiles, old$Quartile),
      ,
      drop = FALSE
    ]

    expected <- c(old$Fraction_control, old$Fraction_treated)

    if (
      any(!is.finite(expected)) ||
      max(abs(reconstructed - expected)) > 1e-10
    ) {
      stop(spec$Samples, ": sample counts do not reproduce the plotted bars.")
    }

    data.frame(
      Dataset = spec$Dataset,
      Comparison = spec$Samples,
      Time_h = spec$Time_h,
      Sample = rep(as.character(sample_counts$replicate), 2),
      Condition = rep(ifelse(ref, "Control", "Treated"), 2),
      Quartile = rep(
        Fig3D_time_quartiles,
        each = nrow(sample_counts)
      ),
      Fraction = c(
        sample_counts$Upper_counts / sample_counts$Total_counts,
        sample_counts$Lower_counts / sample_counts$Total_counts
      )
    )
  }
))

message("All original time-course bar fractions reproduced correctly.")

Fig3D_time_fraction_summary <- Fig3D_time_fraction_data %>%
  dplyr::group_by(Dataset, Time_h, Quartile, Condition) %>%
  dplyr::summarise(
    Mean_fraction = mean(Fraction),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    names_from = Condition,
    values_from = Mean_fraction
  ) %>%
  dplyr::mutate(Delta_fraction = Treated - Control) %>%
  dplyr::arrange(Dataset, Quartile, Time_h)

print(
  as.data.frame(Fig3D_time_fraction_summary),
  row.names = FALSE,
  digits = 7
)

Fig3D_time_test <- function(dataset, quartile) {

  dat <- Fig3D_time_fraction_data[
    Fig3D_time_fraction_data$Dataset == dataset &
      Fig3D_time_fraction_data$Quartile == quartile,
    ,
    drop = FALSE
  ]

  times <- sort(unique(dat$Time_h))
  coefficients <- (times - mean(times)) /
    sum((times - mean(times))^2)

  dat$Time_coefficient <- coefficients[match(dat$Time_h, times)]

  dat$Sampling_group <- paste(
    dat$Condition, dat$Time_h, sep = "_"
  )

  if (dataset == "GSE241735") {

    control_sets <- lapply(times, function(tm) {
      dat$Sample[
        dat$Condition == "Control" & dat$Time_h == tm
      ]
    })

    stopifnot(all(vapply(
      control_sets,
      function(x) setequal(x, control_sets[[1]]),
      logical(1)
    )))

    dat$Sampling_group[dat$Condition == "Control"] <- "Shared_control"
  }

  sample_groups <- unique(dat[, c("Sample", "Sampling_group")])
  stopifnot(!anyDuplicated(sample_groups$Sample))

  dat$Weighted_fraction <- dat$Time_coefficient *
    ifelse(dat$Condition == "Treated", 1, -1) *
    dat$Fraction

  sample_scores <- stats::aggregate(
    Weighted_fraction ~ Sampling_group + Sample,
    data = dat,
    FUN = sum
  )

  components <- do.call(rbind, lapply(
    split(
      sample_scores$Weighted_fraction,
      sample_scores$Sampling_group
    ),
    function(x) {
      stopifnot(length(x) >= 2L, all(is.finite(x)))
      c(
        Mean = mean(x),
        Variance = stats::var(x) / length(x),
        DF = length(x) - 1
      )
    }
  ))

  slope <- sum(components[, "Mean"])
  variance <- sum(components[, "Variance"])

  if (!is.finite(variance) || variance <= 0) {
    stop(dataset, ": cannot estimate temporal-contrast uncertainty.")
  }

  se <- sqrt(variance)

  df <- variance^2 / sum(
    components[, "Variance"]^2 / components[, "DF"]
  )

  p_value <- 2 * stats::pt(-abs(slope / se), df = df)
  ci <- slope + c(-1, 1) * stats::qt(0.975, df = df) * se

  data.frame(
    Dataset = dataset,
    Quartile = quartile,
    Temporal_test = if (length(times) == 2L) {
      "24-72 h change"
    } else {
      "4-24 h linear trend"
    },
    Slope_delta_fraction_per_hour = slope,
    SE_slope = se,
    Degrees_freedom = df,
    CI95_lower = ci[1],
    CI95_upper = ci[2],
    P_temporal = p_value
  )
}

Fig3D_time_results <- dplyr::bind_rows(lapply(
  c("GSE135514", "GSE241735"),
  function(dataset) {
    dplyr::bind_rows(lapply(
      Fig3D_time_quartiles,
      function(quartile) Fig3D_time_test(dataset, quartile)
    ))
  }
))

Fig3D_time_results$P_adj_BH <- stats::p.adjust(
  Fig3D_time_results$P_temporal,
  method = "BH"
)

Fig3D_time_results$Figure_label <- with(
  Fig3D_time_results,
  ifelse(
    P_adj_BH >= 0.05,
    "NS",
    ifelse(
      P_adj_BH < 0.0001,
      "P < 0.0001",
      sprintf("P = %.4f", P_adj_BH)
    )
  )
)

print(
  as.data.frame(Fig3D_time_results),
  row.names = FALSE,
  digits = 7
)

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(scales)

Fig3E_B_perm <- 5000L
Fig3E_n_cores <- 8L

Fig3E_available_cores <- parallel::detectCores(
  logical = TRUE
)

if (
  is.na(Fig3E_available_cores) ||
  Fig3E_available_cores < Fig3E_n_cores
) {
  stop(
    paste0(
      "Eight CPU workers were requested, but only ",
      Fig3E_available_cores,
      " logical cores were detected."
    )
  )
}

Fig3E_output_dir <- figure_output_dir

Fig3E_required_TPM_objects <- c(
  "GSE103021_ARS_1620_LU65_TPM_long",
  "GSE103021_ARS_1620_LU65_TPM_weight_validation",
  "GSE121105_trastuzumab_72h_TPM_long",
  "GSE121105_trastuzumab_72h_TPM_weight_validation",
  "GSE143462_omipalisib_12h_TPM_long",
  "GSE143462_omipalisib_12h_TPM_weight_validation",
  "GSE135514_Das_24h_TPM_long",
  "GSE135514_Das_24h_TPM_weight_validation",
  "GSE135514_Das_72h_TPM_long",
  "GSE135514_Das_72h_TPM_weight_validation",
  "GSE241735_Futibatinib_4h_TPM_long",
  "GSE241735_Futibatinib_4h_TPM_weight_validation",
  "GSE241735_Futibatinib_12h_TPM_long",
  "GSE241735_Futibatinib_12h_TPM_weight_validation",
  "GSE241735_Futibatinib_24h_TPM_long",
  "GSE241735_Futibatinib_24h_TPM_weight_validation"
)

Fig3E_missing_TPM_objects <- Fig3E_required_TPM_objects[
  !vapply(
    Fig3E_required_TPM_objects,
    exists,
    logical(1),
    envir = .GlobalEnv,
    inherits = FALSE
  )
]

if (length(Fig3E_missing_TPM_objects) > 0) {
  stop(
    paste0(
      "Run the completed master script before running Fig. 3E.\n",
      "The following TPM objects are missing:\n",
      paste(Fig3E_missing_TPM_objects, collapse = "\n")
    )
  )
}

Fig3E_analysis_config <- tibble::tribble(
  ~Analysis_ID,
  ~Sample_label,
  ~Dataset,
  ~Treatment,
  ~Analysis_label,
  ~Treatment_condition,
  ~Comparison,
  ~TPM_long_object,
  ~TPM_sample_object,
  ~Seed,

  "GSE103021_ARS_1620_LU65",
  "GSE103021",
  "GSE103021",
  "ARS-1620 LU65",
  "GSE103021\nARS-1620 LU65",
  "ARS-1620",
  "ARS-1620 minus Control",
  "GSE103021_ARS_1620_LU65_TPM_long",
  "GSE103021_ARS_1620_LU65_TPM_weight_validation",
  124L,

  "GSE121105_trastuzumab_72h",
  "GSE121105",
  "GSE121105",
  "Trastuzumab 72h",
  "GSE121105\nTrastuzumab 72h",
  "Trastuzumab",
  "Trastuzumab 72h minus Control",
  "GSE121105_trastuzumab_72h_TPM_long",
  "GSE121105_trastuzumab_72h_TPM_weight_validation",
  125L,

  "GSE143462_omipalisib_12h",
  "GSE143462",
  "GSE143462",
  "Omipalisib 12h",
  "GSE143462\nOmipalisib 12h",
  "Omipalisib",
  "Omipalisib 12h minus Control",
  "GSE143462_omipalisib_12h_TPM_long",
  "GSE143462_omipalisib_12h_TPM_weight_validation",
  126L,

  "GSE135514_Das_24h",
  "GSE135514_24h",
  "GSE135514",
  "Dasatinib 24h",
  "GSE135514\nDasatinib 24h",
  "Dasatinib 24h",
  "Dasatinib 24h minus Control",
  "GSE135514_Das_24h_TPM_long",
  "GSE135514_Das_24h_TPM_weight_validation",
  127L,

  "GSE135514_Das_72h",
  "GSE135514_72h",
  "GSE135514",
  "Dasatinib 72h",
  "GSE135514\nDasatinib 72h",
  "Dasatinib 72h",
  "Dasatinib 72h minus Control",
  "GSE135514_Das_72h_TPM_long",
  "GSE135514_Das_72h_TPM_weight_validation",
  128L,

  "GSE241735_Futibatinib_4h",
  "GSE241735_4h",
  "GSE241735",
  "Futibatinib 4h",
  "GSE241735\nFutibatinib 4h",
  "Futibatinib 4h",
  "Futibatinib 4h minus Control",
  "GSE241735_Futibatinib_4h_TPM_long",
  "GSE241735_Futibatinib_4h_TPM_weight_validation",
  129L,

  "GSE241735_Futibatinib_12h",
  "GSE241735_12h",
  "GSE241735",
  "Futibatinib 12h",
  "GSE241735\nFutibatinib 12h",
  "Futibatinib 12h",
  "Futibatinib 12h minus Control",
  "GSE241735_Futibatinib_12h_TPM_long",
  "GSE241735_Futibatinib_12h_TPM_weight_validation",
  130L,

  "GSE241735_Futibatinib_24h",
  "GSE241735_24h",
  "GSE241735",
  "Futibatinib 24h",
  "GSE241735\nFutibatinib 24h",
  "Futibatinib 24h",
  "Futibatinib 24h minus Control",
  "GSE241735_Futibatinib_24h_TPM_long",
  "GSE241735_Futibatinib_24h_TPM_weight_validation",
  131L
)

Fig3E_run_one_TPM_profile <- function(
    TPM_long,
    TPM_sample,
    analysis_id,
    sample_label,
    dataset_label,
    treatment_label,
    analysis_label,
    treatment_condition,
    comparison_label,
    B_perm = 5000L,
    n_cores = 8L,
    seed = 123L
) {

  required_long_columns <- c(
    "ensembl_transcript_id",
    "Sample",
    "Condition",
    "TPM",
    "Normalized_MFE",
    "relative_transcript_abundance_MFE"
  )

  required_sample_columns <- c(
    "Sample",
    "Condition",
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
      ": missing columns in the TPM-long table: ",
      paste(missing_long_columns, collapse = ", ")
    )
  }

  if (length(missing_sample_columns) > 0) {
    stop(
      analysis_id,
      ": missing columns in the TPM sample table: ",
      paste(missing_sample_columns, collapse = ", ")
    )
  }

  sample_values <- TPM_sample %>%
    dplyr::transmute(
      Analysis_ID = analysis_id,
      Sample_label = sample_label,
      Dataset = dataset_label,
      Treatment = treatment_label,
      Analysis_label = analysis_label,
      Sample = as.character(Sample),
      Condition = as.character(Condition),
      annotated_transcripts = as.integer(annotated_transcripts),
      expressed_annotated_transcripts =
        as.integer(expressed_annotated_transcripts),
      retained_TPM_MFE = as.numeric(retained_TPM_MFE),
      sum_MFE_weights = as.numeric(sum_MFE_weights),
      TPM_weighted_Normalized_MFE =
        as.numeric(TPM_weighted_Normalized_MFE)
    )

  if (
    !all(
      sample_values$Condition %in%
      c("Control", treatment_condition)
    )
  ) {
    stop(
      analysis_id,
      ": unexpected condition labels were found."
    )
  }

  if (
    sum(sample_values$Condition == "Control") < 2 ||
    sum(sample_values$Condition == treatment_condition) < 2
  ) {
    stop(
      analysis_id,
      ": at least two biological replicates are required per condition."
    )
  }

  if (
    any(
      abs(sample_values$sum_MFE_weights - 1) > 1e-8
    )
  ) {
    stop(
      analysis_id,
      ": one or more sample-level TPM weights do not sum to one."
    )
  }

  condition_stats <- sample_values %>%
    dplyr::group_by(
      Analysis_ID,
      Sample_label,
      Dataset,
      Treatment,
      Analysis_label,
      Condition
    ) %>%
    dplyr::summarise(
      n_samples = dplyr::n(),
      mean_TPM_weighted_Normalized_MFE = mean(
        TPM_weighted_Normalized_MFE,
        na.rm = TRUE
      ),
      sd_TPM_weighted_Normalized_MFE = stats::sd(
        TPM_weighted_Normalized_MFE,
        na.rm = TRUE
      ),
      sem_TPM_weighted_Normalized_MFE =
        sd_TPM_weighted_Normalized_MFE /
        sqrt(n_samples),
      .groups = "drop"
    )

  control_row <- condition_stats %>%
    dplyr::filter(Condition == "Control")

  treatment_row <- condition_stats %>%
    dplyr::filter(Condition == treatment_condition)

  if (nrow(control_row) != 1 || nrow(treatment_row) != 1) {
    stop(
      analysis_id,
      ": control or treatment condition could not be identified uniquely."
    )
  }

  control_value <-
    control_row$mean_TPM_weighted_Normalized_MFE[[1]]

  treatment_value <-
    treatment_row$mean_TPM_weighted_Normalized_MFE[[1]]

  control_SE <-
    control_row$sem_TPM_weighted_Normalized_MFE[[1]]

  treatment_SE <-
    treatment_row$sem_TPM_weighted_Normalized_MFE[[1]]

  observed_delta <- treatment_value - control_value

  delta_SE <- sqrt(
    control_SE^2 + treatment_SE^2
  )

  condition_profile <- TPM_long %>%
    dplyr::transmute(
      ensembl_transcript_id =
        as.character(ensembl_transcript_id),
      Sample = as.character(Sample),
      Condition = as.character(Condition),
      Normalized_MFE = as.numeric(Normalized_MFE),
      relative_TPM_weight =
        as.numeric(relative_transcript_abundance_MFE)
    ) %>%
    dplyr::filter(
      Condition %in% c("Control", treatment_condition),
      is.finite(Normalized_MFE),
      is.finite(relative_TPM_weight),
      relative_TPM_weight >= 0
    ) %>%
    dplyr::group_by(
      ensembl_transcript_id,
      Normalized_MFE,
      Condition
    ) %>%
    dplyr::summarise(
      mean_relative_TPM_weight = mean(
        relative_TPM_weight,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(
      names_from = Condition,
      values_from = mean_relative_TPM_weight,
      values_fill = 0
    )

  if (
    !("Control" %in% colnames(condition_profile)) ||
    !(treatment_condition %in% colnames(condition_profile))
  ) {
    stop(
      analysis_id,
      ": failed to construct both condition-level TPM profiles."
    )
  }

  w_control <- as.numeric(
    condition_profile[["Control"]]
  )

  w_treatment <- as.numeric(
    condition_profile[[treatment_condition]]
  )

  x_mfe <- as.numeric(
    condition_profile$Normalized_MFE
  )

  if (
    abs(sum(w_control) - 1) > 1e-8 ||
    abs(sum(w_treatment) - 1) > 1e-8
  ) {
    stop(
      analysis_id,
      ": condition-level mean TPM weights do not sum to one."
    )
  }

  profile_control_value <- sum(
    w_control * x_mfe
  )

  profile_treatment_value <- sum(
    w_treatment * x_mfe
  )

  profile_delta <-
    profile_treatment_value - profile_control_value

  if (
    abs(profile_control_value - control_value) > 1e-10 ||
    abs(profile_treatment_value - treatment_value) > 1e-10 ||
    abs(profile_delta - observed_delta) > 1e-10
  ) {
    stop(
      analysis_id,
      ": transcript-profile and sample-summary calculations do not agree."
    )
  }

  contrast_weight <- w_treatment - w_control

  if (abs(sum(contrast_weight)) > 1e-8) {
    stop(
      analysis_id,
      ": treatment-minus-control TPM contrast weights do not sum to zero."
    )
  }

  n_transcripts <- length(x_mfe)

  active_indices <- which(
    contrast_weight != 0
  )

  active_contrast_weight <- contrast_weight[
    active_indices
  ]

  n_active_weights <- length(
    active_contrast_weight
  )

  if (n_active_weights == 0L) {
    stop(
      analysis_id,
      ": all treatment-minus-control TPM contrast weights are zero."
    )
  }

  B_perm <- as.integer(B_perm)
  n_cores <- as.integer(n_cores)

  if (length(B_perm) != 1L || is.na(B_perm) || B_perm < 1L) {
    stop(
      analysis_id,
      ": B_perm must be one positive integer."
    )
  }

  if (length(n_cores) != 1L || is.na(n_cores) || n_cores < 1L) {
    stop(
      analysis_id,
      ": n_cores must be one positive integer."
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
      active_contrast_weight * sampled_MFE
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
      ": one or more permutation statistics were not generated correctly."
    )
  }

  permutation_p_value <- (
    sum(
      abs(permutation_delta) >= abs(observed_delta)
    ) + 1
  ) / (B_perm + 1)

  profile_stats <- tibble::tibble(
    Analysis_ID = analysis_id,
    Sample = sample_label,
    Dataset = dataset_label,
    Treatment = treatment_label,
    Analysis_label = analysis_label,
    comparison = comparison_label,
    abundance_measure = "GREIN transcript-level TPM",
    transcript_universe = paste(
      "Transcripts matched to Table_ANNOTATED",
      "with finite Normalized_MFE"
    ),
    n_control = control_row$n_samples[[1]],
    n_treatment = treatment_row$n_samples[[1]],
    n_transcripts_used = n_transcripts,
    n_transcripts_expressed_in_either_condition = sum(
      w_control > 0 | w_treatment > 0
    ),
    n_nonzero_contrast_weights = n_active_weights,
    control_value = control_value,
    treatment_value = treatment_value,
    delta_abundance_weighted_Normalized_MFE =
      observed_delta,
    control_SE = control_SE,
    treatment_SE = treatment_SE,
    delta_SE = delta_SE,
    errorbar_low = observed_delta - delta_SE,
    errorbar_high = observed_delta + delta_SE,
    percent_change_treatment_minus_control =
      100 * observed_delta / abs(control_value),
    permutation_p_value = permutation_p_value,
    permutation_test = paste(
      "Normalized_MFE permuted across transcripts;",
      "treatment-minus-control TPM weights held fixed"
    ),
    B_perm = B_perm,
    permutation_cores = n_cores,
    direction = dplyr::case_when(
      observed_delta > 0 ~
        "Shift toward less-negative predicted DeltaGfold",
      observed_delta < 0 ~
        "Shift toward more-negative predicted DeltaGfold",
      TRUE ~
        "No change"
    )
  )

  list(
    sample_values = sample_values,
    condition_stats = condition_stats,
    profile_stats = profile_stats
  )
}

Fig3E_TPM_results <- vector(
  mode = "list",
  length = nrow(Fig3E_analysis_config)
)

for (analysis_index in seq_len(nrow(Fig3E_analysis_config))) {

  current_config <- Fig3E_analysis_config[
    analysis_index,
    ,
    drop = FALSE
  ]

  message(
    "Running Fig. 3E TPM analysis: ",
    current_config$Analysis_ID[[1]]
  )

  Fig3E_TPM_results[[analysis_index]] <-
    Fig3E_run_one_TPM_profile(
      TPM_long = get(
        current_config$TPM_long_object[[1]],
        envir = .GlobalEnv,
        inherits = FALSE
      ),
      TPM_sample = get(
        current_config$TPM_sample_object[[1]],
        envir = .GlobalEnv,
        inherits = FALSE
      ),
      analysis_id = current_config$Analysis_ID[[1]],
      sample_label = current_config$Sample_label[[1]],
      dataset_label = current_config$Dataset[[1]],
      treatment_label = current_config$Treatment[[1]],
      analysis_label = current_config$Analysis_label[[1]],
      treatment_condition =
        current_config$Treatment_condition[[1]],
      comparison_label = current_config$Comparison[[1]],
      B_perm = Fig3E_B_perm,
      n_cores = Fig3E_n_cores,
      seed = current_config$Seed[[1]]
    )

  invisible(gc())
}

names(Fig3E_TPM_results) <-
  Fig3E_analysis_config$Analysis_ID

Fig3E_TPM_sample_values <- dplyr::bind_rows(
  lapply(
    Fig3E_TPM_results,
    function(result) result$sample_values
  )
)

Fig3E_TPM_condition_stats <- dplyr::bind_rows(
  lapply(
    Fig3E_TPM_results,
    function(result) result$condition_stats
  )
)

Fig3_GREIN_abundance_weighted_DGfold_delta_stats <-
  dplyr::bind_rows(
    lapply(
      Fig3E_TPM_results,
      function(result) result$profile_stats
    )
  ) %>%
  dplyr::mutate(
    permutation_p_adj_BH = stats::p.adjust(
      permutation_p_value,
      method = "BH"
    ),
    permutation_p_adj_BH_label = dplyr::case_when(
      permutation_p_adj_BH < 0.001 ~
        "P < 0.001",
      TRUE ~ paste0(
        "P = ",
        format(
          signif(permutation_p_adj_BH, 2),
          scientific = FALSE,
          trim = TRUE
        )
      )
    ),
    delta_label = dplyr::case_when(
      abs(delta_abundance_weighted_Normalized_MFE) < 0.001 ~
        sprintf(
          "%.5f",
          delta_abundance_weighted_Normalized_MFE
        ),
      TRUE ~
        sprintf(
          "%.4f",
          delta_abundance_weighted_Normalized_MFE
        )
    )
  )

Fig3E_sample_order <- c(
  "GSE103021",
  "GSE121105",
  "GSE143462",
  "Spacer_1",
  "GSE135514_24h",
  "GSE135514_72h",
  "Spacer_2",
  "GSE241735_4h",
  "GSE241735_12h",
  "GSE241735_24h"
)

Fig3E_sample_labels <- c(
  "GSE103021" = "GSE103021",
  "GSE121105" = "GSE121105",
  "GSE143462" = "GSE143462",
  "Spacer_1" = "",
  "GSE135514_24h" = "GSE135514_24h",
  "GSE135514_72h" = "GSE135514_72h",
  "Spacer_2" = "",
  "GSE241735_4h" = "GSE241735_4h",
  "GSE241735_12h" = "GSE241735_12h",
  "GSE241735_24h" = "GSE241735_24h"
)

Fig3_GREIN_abundance_weighted_DGfold_delta_stats <-
  Fig3_GREIN_abundance_weighted_DGfold_delta_stats %>%
  dplyr::mutate(
    Sample = factor(
      Sample,
      levels = Fig3E_sample_order
    )
  ) %>%
  dplyr::arrange(Sample)

cat("\nFig. 3E TPM-based sample-level values:\n")

print(
  as.data.frame(Fig3E_TPM_sample_values),
  row.names = FALSE,
  digits = 12
)

cat("\nFig. 3E TPM-based condition-level statistics:\n")

print(
  as.data.frame(Fig3E_TPM_condition_stats),
  row.names = FALSE,
  digits = 12
)

cat("\nFig. 3E TPM-based treatment-minus-control statistics:\n")

print(
  as.data.frame(
    Fig3_GREIN_abundance_weighted_DGfold_delta_stats
  ),
  row.names = FALSE,
  digits = 12
)

Fig3E_max_errorbar_extent <- max(
  abs(
    c(
      Fig3_GREIN_abundance_weighted_DGfold_delta_stats$errorbar_low,
      Fig3_GREIN_abundance_weighted_DGfold_delta_stats$errorbar_high
    )
  ),
  na.rm = TRUE
)

if (
  !is.finite(Fig3E_max_errorbar_extent) ||
  Fig3E_max_errorbar_extent == 0
) {
  Fig3E_max_errorbar_extent <- 0.001
}

Fig3E_delta_label_offset <-
  0.035 * Fig3E_max_errorbar_extent

Fig3E_p_label_offset <-
  0.23 * Fig3E_max_errorbar_extent

Fig3E_y_limit <-
  1.55 * Fig3E_max_errorbar_extent

Fig3_GREIN_abundance_weighted_DGfold_delta_stats <-
  Fig3_GREIN_abundance_weighted_DGfold_delta_stats %>%
  dplyr::mutate(
    delta_label_y = dplyr::if_else(
      delta_abundance_weighted_Normalized_MFE >= 0,
      errorbar_high + Fig3E_delta_label_offset,
      errorbar_low - Fig3E_delta_label_offset
    ),
    p_label_y = dplyr::if_else(
      delta_abundance_weighted_Normalized_MFE >= 0,
      errorbar_high + Fig3E_p_label_offset,
      errorbar_low - Fig3E_p_label_offset
    ),
    text_hjust = dplyr::if_else(
      delta_abundance_weighted_Normalized_MFE >= 0,
      0,
      1
    )
  )

Plot_Fig3_GREIN_abundance_weighted_DGfold_delta <- ggplot(
  Fig3_GREIN_abundance_weighted_DGfold_delta_stats,
  aes(
    x = Sample,
    y = delta_abundance_weighted_Normalized_MFE
  )
) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.4,
    linetype = "dashed",
    color = "black"
  ) +
  geom_col(
    width = 0.65,
    fill = "#D55E00",
    color = "black",
    linewidth = 0.3
  ) +
  geom_errorbar(
    aes(
      ymin = delta_abundance_weighted_Normalized_MFE,
      ymax = errorbar_high
    ),
    width = 0.18,
    linewidth = 0.4,
    color = "black"
  ) +
  geom_text(
    aes(
      y = delta_label_y,
      label = delta_label,
      hjust = text_hjust
    ),
    angle = 90,
    size = 1,
    color = "black"
  ) +
  geom_text(
    aes(
      y = p_label_y,
      label = permutation_p_adj_BH_label,
      hjust = text_hjust
    ),
    angle = 90,
    size = 1,
    color = "black"
  ) +
  scale_x_discrete(
    drop = FALSE,
    labels = Fig3E_sample_labels,
    expand = expansion(
      add = c(1.0, 0.6)
    )
  ) +
  scale_y_continuous(
    labels = scales::number_format(
      accuracy = 0.001
    ),
    expand = expansion(
      mult = c(0.02, 0.02)
    )
  ) +
  coord_cartesian(
    ylim = c(-0.02, 0.02),
    clip = "off"
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
    plot.margin = margin(
      t = 8,
      r = 8,
      b = 8,
      l = 8
    )
  ) +
  labs(
    title = "",
    subtitle = "",
    x = "",
    y = expression(
      Delta * Delta * G[fold] ~
        (kcal %.% mol^{-1} %.% nt^{-1})
    )
  )

print(
  Plot_Fig3_GREIN_abundance_weighted_DGfold_delta
)

ggsave(
  filename =
    "Plot_Fig3_GREIN_TPM_weighted_DGfold_delta.tiff",
  plot =
    Plot_Fig3_GREIN_abundance_weighted_DGfold_delta,
  path = Fig3E_output_dir,
  device = "tiff",
  width = 5.5,
  height = 4.5,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

Fig3_GREIN_abundance_weighted_DGfold_delta_stats <-
  dplyr::bind_rows(
    lapply(
      Fig3E_TPM_results,
      function(result) result$profile_stats
    )
  ) %>%
  dplyr::mutate(

    permutation_p_adj_BH = stats::p.adjust(
      permutation_p_value,
      method = "BH"
    ),

    permutation_p_adj_BH_label = dplyr::case_when(
      is.na(permutation_p_adj_BH) ~ NA_character_,

      permutation_p_adj_BH < 0.0001 ~
        "P < 0.0001",

      TRUE ~ paste0(
        "P = ",
        trimws(
          base::format.pval(
            permutation_p_adj_BH,
            digits = 4
          )
        )
      )
    ),

    delta_label = dplyr::case_when(
      abs(delta_abundance_weighted_Normalized_MFE) < 0.001 ~
        sprintf(
          "%.5f",
          delta_abundance_weighted_Normalized_MFE
        ),

      TRUE ~
        sprintf(
          "%.4f",
          delta_abundance_weighted_Normalized_MFE
        )
    )
  )

print(Fig3_GREIN_abundance_weighted_DGfold_delta_stats)

library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(ggpubr)
library(scales)
library(grid)

Fig3F_I_output_dir <- figure_output_dir

GSE97681_Fig3_environment <- .GlobalEnv

GSE97681_Fig3_loaded_objects <- ls(
  envir = .GlobalEnv,
  all.names = TRUE
)

GSE97681_Fig3_required_objects <- c(
  "GSE97681_sample_map",
  "GSE97681_excluded_samples",
  "GSE97681_acute_48h_growth_arrest_eligibility",
  "GSE97681_Fig3_dot_plot_single_panel_data",
  "MFE_res_GSE97681_WM989_acute_vs_untreated",
  "MFE_res_GSE97681_WM989_resistant_vs_acute",
  "MFE_res_GSE97681_WM983B_acute_vs_untreated",
  "MFE_res_GSE97681_WM983B_resistant_vs_acute",
  "GSE97681_WM989_acute_vs_untreated_permutation_p_value",
  "GSE97681_WM989_resistant_vs_acute_permutation_p_value",
  "GSE97681_WM983B_acute_vs_untreated_permutation_p_value",
  "GSE97681_WM983B_resistant_vs_acute_permutation_p_value",
  "GSE97681_Fig3D_complete_stats",
  "GSE97681_primary_DGfold_contrast_stats",
  "GSE97681_Fig3D_transcript_count_data",
  "GSE97681_founder_balanced_trajectory_wide"
)

GSE97681_Fig3_missing_objects <- setdiff(
  GSE97681_Fig3_required_objects,
  GSE97681_Fig3_loaded_objects
)

if (length(GSE97681_Fig3_missing_objects) > 0) {
  stop(
    paste0(
      "The GSE97681 environment is missing objects required for Fig. 3F-I:\n",
      paste(
        GSE97681_Fig3_missing_objects,
        collapse = "\n"
      )
    )
  )
}

Fig3F_GSE97681_GSEA_data <- get(
  "GSE97681_Fig3_dot_plot_single_panel_data",
  envir = GSE97681_Fig3_environment,
  inherits = FALSE
) %>%
  tibble::as_tibble() %>%
  dplyr::select(
    cell_line,
    comparison,
    ID,
    Description,
    NES,
    pvalue,
    p.adjust,
    significant_GSEA_result,
    neg_log10_p_adjust
  )

Fig3G_GSE97681_violin_data <- dplyr::bind_rows(
  get(
    "MFE_res_GSE97681_WM989_acute_vs_untreated",
    envir = GSE97681_Fig3_environment,
    inherits = FALSE
  ) %>%
    tibble::as_tibble() %>%
    dplyr::transmute(
      cell_line = "WM989-A6",
      comparison = "Acute_48h versus Untreated",
      log2FoldChange,
      Normalized_MFE
    ),
  get(
    "MFE_res_GSE97681_WM989_resistant_vs_acute",
    envir = GSE97681_Fig3_environment,
    inherits = FALSE
  ) %>%
    tibble::as_tibble() %>%
    dplyr::transmute(
      cell_line = "WM989-A6",
      comparison = "Resistant versus Acute_48h",
      log2FoldChange,
      Normalized_MFE
    ),
  get(
    "MFE_res_GSE97681_WM983B_acute_vs_untreated",
    envir = GSE97681_Fig3_environment,
    inherits = FALSE
  ) %>%
    tibble::as_tibble() %>%
    dplyr::transmute(
      cell_line = "WM983B-E9",
      comparison = "Acute_48h versus Untreated",
      log2FoldChange,
      Normalized_MFE
    ),
  get(
    "MFE_res_GSE97681_WM983B_resistant_vs_acute",
    envir = GSE97681_Fig3_environment,
    inherits = FALSE
  ) %>%
    tibble::as_tibble() %>%
    dplyr::transmute(
      cell_line = "WM983B-E9",
      comparison = "Resistant versus Acute_48h",
      log2FoldChange,
      Normalized_MFE
    )
)

Fig3G_GSE97681_permutation_statistics <- tibble::tibble(
  cell_line = c(
    "WM989-A6",
    "WM989-A6",
    "WM983B-E9",
    "WM983B-E9"
  ),
  comparison = c(
    "Acute_48h versus Untreated",
    "Resistant versus Acute_48h",
    "Acute_48h versus Untreated",
    "Resistant versus Acute_48h"
  ),
  permutation_p_value = c(
    as.numeric(
      get(
        "GSE97681_WM989_acute_vs_untreated_permutation_p_value",
        envir = GSE97681_Fig3_environment,
        inherits = FALSE
      )
    ),
    as.numeric(
      get(
        "GSE97681_WM989_resistant_vs_acute_permutation_p_value",
        envir = GSE97681_Fig3_environment,
        inherits = FALSE
      )
    ),
    as.numeric(
      get(
        "GSE97681_WM983B_acute_vs_untreated_permutation_p_value",
        envir = GSE97681_Fig3_environment,
        inherits = FALSE
      )
    ),
    as.numeric(
      get(
        "GSE97681_WM983B_resistant_vs_acute_permutation_p_value",
        envir = GSE97681_Fig3_environment,
        inherits = FALSE
      )
    )
  )
)

Fig3H_GSE97681_quartile_statistics <- get(
  "GSE97681_Fig3D_complete_stats",
  envir = GSE97681_Fig3_environment,
  inherits = FALSE
) %>%
  tibble::as_tibble() %>%
  dplyr::select(
    Samples,
    reference_high_MFE_fraction,
    comparison_high_MFE_fraction,
    reference_low_MFE_fraction,
    comparison_low_MFE_fraction
  )

Fig3I_GSE97681_DeltaDeltaGfold_statistics <- get(
  "GSE97681_primary_DGfold_contrast_stats",
  envir = GSE97681_Fig3_environment,
  inherits = FALSE
) %>%
  tibble::as_tibble()

Fig3F_I_GSE97681_sample_map <- get(
  "GSE97681_sample_map",
  envir = GSE97681_Fig3_environment,
  inherits = FALSE
) %>%
  tibble::as_tibble() %>%
  dplyr::select(
    Sample,
    cell_line,
    condition,
    primary_state
  )

Fig3F_I_GSE97681_excluded_samples <- get(
  "GSE97681_excluded_samples",
  envir = GSE97681_Fig3_environment,
  inherits = FALSE
)

Fig3F_I_GSE97681_growth_arrest_eligibility <- get(
  "GSE97681_acute_48h_growth_arrest_eligibility",
  envir = GSE97681_Fig3_environment,
  inherits = FALSE
) %>%
  tibble::as_tibble()

rm(
  GSE97681_Fig3_environment,
  GSE97681_Fig3_loaded_objects,
  GSE97681_Fig3_missing_objects
)

invisible(gc())

Fig3F_I_excluded_samples <- c(
  "GSM2575347",
  "GSM2575397",
  "GSM2575399"
)

if (
  !setequal(
    as.character(
      Fig3F_I_GSE97681_excluded_samples
    ),
    Fig3F_I_excluded_samples
  )
) {
  stop(
    "The saved GSE97681 environment does not contain the expected three-sample exclusion set."
  )
}

if (
  length(
    intersect(
      Fig3F_I_excluded_samples,
      Fig3F_I_GSE97681_sample_map$Sample
    )
  ) > 0
) {
  stop(
    "At least one of the three excluded WM983B-E9 samples remains in the GSE97681 environment."
  )
}

if (
  any(
    Fig3F_I_GSE97681_sample_map$primary_state &
    grepl(
      "^Withdrawal_",
      as.character(
        Fig3F_I_GSE97681_sample_map$condition
      )
    )
  )
) {
  stop(
    "A withdrawal sample is marked as a primary GSE97681 sample."
  )
}

Fig3F_I_observed_sample_counts <-
  Fig3F_I_GSE97681_sample_map %>%
  dplyr::filter(primary_state) %>%
  dplyr::mutate(
    condition = as.character(condition)
  ) %>%
  dplyr::count(
    cell_line,
    condition,
    name = "observed_n"
  )

Fig3F_I_expected_sample_counts <- tibble::tribble(
  ~cell_line, ~condition, ~expected_n,
  "WM989-A6", "Untreated", 8L,
  "WM989-A6", "Acute_48h", 8L,
  "WM989-A6", "Resistant", 10L,
  "WM983B-E9", "Untreated", 20L,
  "WM983B-E9", "Acute_48h", 20L,
  "WM983B-E9", "Resistant", 35L
)

Fig3F_I_sample_count_check <-
  Fig3F_I_expected_sample_counts %>%
  dplyr::full_join(
    Fig3F_I_observed_sample_counts,
    by = c(
      "cell_line",
      "condition"
    )
  ) %>%
  dplyr::filter(
    is.na(expected_n) |
      is.na(observed_n) |
      expected_n != observed_n
  )

if (nrow(Fig3F_I_sample_count_check) > 0) {
  stop(
    paste0(
      "The GSE97681 primary sample groups do not contain the expected ",
      "8/8/10 and 20/20/35 samples."
    )
  )
}

Fig3F_I_growth_arrest_eligibility_check <-
  Fig3F_I_GSE97681_growth_arrest_eligibility %>%
  dplyr::filter(
    cell_line %in% c(
      "WM989-A6",
      "WM983B-E9"
    )
  )

if (
  nrow(Fig3F_I_growth_arrest_eligibility_check) != 2 ||
  any(
    is.na(
      Fig3F_I_growth_arrest_eligibility_check[[
        "mandatory_growth_arrest_eligible"
      ]]
    )
  ) ||
  any(
    !Fig3F_I_growth_arrest_eligibility_check[[
      "mandatory_growth_arrest_eligible"
    ]]
  )
) {
  stop(
    paste0(
      "Both GSE97681 Acute_48h groups must pass the prespecified ",
      "growth-arrest GSEA screen before panels F-I are generated."
    )
  )
}

Fig3F_description_order <- c(
  "Signt_1",
  "Signt_2",
  "Signt_3",
  "Cell cycle",
  "DNA replication",
  "Homologous recombination",
  "Mismatch repair",
  "Fanconi anemia pathway",
  "Base excision repair",
  "Nucleotide excision repair",
  "ATP-dependent chromatin remodeling",
  "p53 pathway",
  "PPAR pathway",
  "Apelin pathway",
  "C-type lectin receptor pathway",
  "ErbB pathway",
  "Hippo pathway",
  "IL-17 pathway",
  "Insulin pathway",
  "JAK-STAT pathway",
  "MAPK pathway",
  "AMPK pathway",
  "mTOR pathway",
  "PI3K-Akt pathway",
  "Ras pathway",
  "Rap1 pathway",
  "TGF-beta pathway",
  "TNF pathway",
  "Wnt pathway",
  "Notch pathway"
)

Fig3F_required_columns <- c(
  "cell_line",
  "comparison",
  "ID",
  "Description",
  "NES",
  "pvalue",
  "p.adjust",
  "significant_GSEA_result"
)

Fig3F_missing_columns <- setdiff(
  Fig3F_required_columns,
  colnames(Fig3F_GSE97681_GSEA_data)
)

if (length(Fig3F_missing_columns) > 0) {
  stop(
    paste0(
      "Panel F is missing the following columns:\n",
      paste(Fig3F_missing_columns, collapse = "\n")
    )
  )
}

if (
  nrow(Fig3F_GSE97681_GSEA_data) != 120 ||
  anyDuplicated(
    paste(
      Fig3F_GSE97681_GSEA_data$cell_line,
      Fig3F_GSE97681_GSEA_data$comparison,
      Fig3F_GSE97681_GSEA_data$ID,
      sep = "::"
    )
  ) > 0
) {
  stop(
    paste0(
      "Panel F must contain 120 unique rows: 2 cell lines x ",
      "2 comparisons x 30 signatures/pathways."
    )
  )
}

if (
  anyNA(
    Fig3F_GSE97681_GSEA_data$significant_GSEA_result
  ) ||
  any(
    Fig3F_GSE97681_GSEA_data$significant_GSEA_result &
    (
      !is.finite(
        Fig3F_GSE97681_GSEA_data$NES
      ) |
      !is.finite(
        Fig3F_GSE97681_GSEA_data$p.adjust
      )
    )
  )
) {
  stop(
    "Panel F contains an invalid significant GSEA result."
  )
}

Fig3F_GSE97681_GSEA_data <-
  Fig3F_GSE97681_GSEA_data %>%
  dplyr::mutate(
    neg_log10_p_adjust = dplyr::if_else(
      significant_GSEA_result,
      -log10(
        pmax(
          p.adjust,
          .Machine$double.xmin
        )
      ),
      NA_real_
    ),
    Description = factor(
      as.character(Description),
      levels = Fig3F_description_order
    ),
    cell_line = factor(
      as.character(cell_line),
      levels = c(
        "WM989-A6",
        "WM983B-E9"
      )
    ),
    comparison = factor(
      as.character(comparison),
      levels = c(
        "Acute_48h versus Untreated",
        "Resistant versus Acute_48h"
      ),
      labels = c(
        "Acute 48 h\nversus untreated",
        "Stable resistant (drug on)\nversus acute 48 h"
      )
    )
  )

if (
  anyNA(Fig3F_GSE97681_GSEA_data$Description) ||
  anyNA(Fig3F_GSE97681_GSEA_data$cell_line) ||
  anyNA(Fig3F_GSE97681_GSEA_data$comparison)
) {
  stop(
    "Panel F contains an unrecognized pathway, cell line or comparison."
  )
}

Plot_Fig3F_GSE97681_GSEA <- ggplot(
  Fig3F_GSE97681_GSEA_data,
  aes(
    x = comparison,
    y = Description,
    fill = NES,
    size = neg_log10_p_adjust
  )
) +
  geom_blank(
    data = Fig3F_GSE97681_GSEA_data,
    aes(
      x = comparison,
      y = Description
    ),
    inherit.aes = FALSE
  ) +
  geom_point(
    data = dplyr::filter(
      Fig3F_GSE97681_GSEA_data,
      significant_GSEA_result
    ),
    shape = 21,
    color = "black",
    stroke = 0.5
  ) +
  facet_grid(
    cols = vars(cell_line),
    scales = "free_x",
    space = "free_x"
  ) +
  scale_fill_gradient2(
    low = "dodgerblue4",
    mid = "white",
    high = "firebrick3",
    midpoint = 0,
    limits = c(-3, 3),
    oob = scales::squish
  ) +
  scale_size_continuous(
    range = c(2, 10)
  ) +
  scale_x_discrete(
    drop = FALSE
  ) +
  scale_y_discrete(
    drop = FALSE
  ) +
  theme_bw() +
  labs(
    title = "",
    x = "",
    y = "",
    fill = "NES",
    size = "-log10(p.adjust)"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 12,
      color = "black"
    ),
    axis.text.y = element_text(
      size = 14,
      color = "black"
    ),
    strip.background = element_blank(),
    strip.text = element_text(
      size = 14,
      color = "black"
    ),
    panel.spacing.x = grid::unit(
      0.8,
      "lines"
    )
  )

print(Plot_Fig3F_GSE97681_GSEA)

ggsave(
  filename = "Fig3F_GSE97681_GSEA.tiff",
  plot = Plot_Fig3F_GSE97681_GSEA,
  path = Fig3F_I_output_dir,
  device = "tiff",
  width = 8,
  height = 12,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

Fig3G_cell_line_levels <- c(
  "WM989-A6",
  "WM983B-E9"
)

Fig3G_comparison_keys <- c(
  "Acute_48h versus Untreated",
  "Resistant versus Acute_48h"
)

Fig3G_comparison_labels <- c(
  "Acute 48 h versus untreated",
  "Stable resistant (drug on) versus acute 48 h"
)

Fig3G_standardize_comparison <- function(x) {

  x <- as.character(x)

  dplyr::case_when(
    x %in% c(
      "Acute_48h versus Untreated",
      "Acute 48 h versus untreated"
    ) ~
      "Acute_48h versus Untreated",

    x %in% c(
      "Resistant versus Acute_48h",
      "Stable resistant (drug on) versus acute 48 h"
    ) ~
      "Resistant versus Acute_48h",

    TRUE ~ NA_character_
  )
}

Fig3G_GSE97681_violin_plot_data <-
  Fig3G_GSE97681_violin_data %>%
  dplyr::filter(
    is.finite(log2FoldChange),
    is.finite(Normalized_MFE),
    log2FoldChange != 0
  ) %>%
  dplyr::mutate(

    expression_group = factor(
      dplyr::if_else(
        log2FoldChange > 0,
        "LFC > 0",
        "LFC < 0"
      ),
      levels = c(
        "LFC < 0",
        "LFC > 0"
      )
    ),

    cell_line = factor(
      as.character(cell_line),
      levels = Fig3G_cell_line_levels
    ),

    comparison_key =
      Fig3G_standardize_comparison(comparison),

    comparison_label = factor(
      comparison_key,
      levels = Fig3G_comparison_keys,
      labels = Fig3G_comparison_labels
    )
  )

if (
  anyNA(Fig3G_GSE97681_violin_plot_data$cell_line) ||
  anyNA(Fig3G_GSE97681_violin_plot_data$comparison_key) ||
  anyNA(Fig3G_GSE97681_violin_plot_data$comparison_label) ||
  anyNA(Fig3G_GSE97681_violin_plot_data$expression_group)
) {
  stop(
    paste0(
      "Panel G contains missing or unrecognized grouping labels. ",
      "Rerun the preceding GSE97681 environment-loading and ",
      "data-extraction section before running Panel G."
    )
  )
}

if (
  any(
    abs(
      Fig3G_GSE97681_violin_plot_data$log2FoldChange
    ) >= 10
  )
) {
  stop(
    paste0(
      "Panel G contains a transcript outside the required ",
      "-10 to +10 log2-fold-change display range."
    )
  )
}

Fig3G_GSE97681_group_counts <-
  Fig3G_GSE97681_violin_plot_data %>%
  dplyr::count(
    cell_line,
    comparison_key,
    expression_group,
    name = "n_transcripts"
  )

if (
  nrow(Fig3G_GSE97681_group_counts) != 8 ||
  any(Fig3G_GSE97681_group_counts$n_transcripts < 2)
) {

  print(
    as.data.frame(Fig3G_GSE97681_group_counts),
    row.names = FALSE
  )

  stop(
    paste0(
      "Panel G does not contain both expression directions ",
      "for all four contrasts."
    )
  )
}

Fig3G_GSE97681_permutation_plot_statistics <-
  Fig3G_GSE97681_permutation_statistics %>%
  dplyr::transmute(

    cell_line = factor(
      as.character(cell_line),
      levels = Fig3G_cell_line_levels
    ),

    comparison_key =
      Fig3G_standardize_comparison(comparison),

    permutation_p_value =
      as.numeric(permutation_p_value)
  )

if (
  nrow(Fig3G_GSE97681_permutation_plot_statistics) != 4 ||
  nrow(
    dplyr::distinct(
      Fig3G_GSE97681_permutation_plot_statistics,
      cell_line,
      comparison_key
    )
  ) != 4 ||
  anyNA(
    Fig3G_GSE97681_permutation_plot_statistics$cell_line
  ) ||
  anyNA(
    Fig3G_GSE97681_permutation_plot_statistics$comparison_key
  ) ||
  any(
    !is.finite(
      Fig3G_GSE97681_permutation_plot_statistics$
      permutation_p_value
    )
  ) ||
  any(
    Fig3G_GSE97681_permutation_plot_statistics$
    permutation_p_value < 0 |
    Fig3G_GSE97681_permutation_plot_statistics$
    permutation_p_value > 1
  )
) {
  stop(
    "Panel G does not contain four unique and valid permutation P values."
  )
}

Fig3G_GSE97681_violin_statistics <-
  Fig3G_GSE97681_violin_plot_data %>%
  dplyr::group_by(
    cell_line,
    comparison_key
  ) %>%
  dplyr::summarise(

    n_LFC_below_zero = sum(
      expression_group == "LFC < 0"
    ),

    n_LFC_above_zero = sum(
      expression_group == "LFC > 0"
    ),

    mean_MFE_LFC_below_zero = mean(
      Normalized_MFE[
        expression_group == "LFC < 0"
      ]
    ),

    mean_MFE_LFC_above_zero = mean(
      Normalized_MFE[
        expression_group == "LFC > 0"
      ]
    ),

    .groups = "drop"
  ) %>%
  dplyr::left_join(
    Fig3G_GSE97681_permutation_plot_statistics,
    by = c(
      "cell_line",
      "comparison_key"
    )
  ) %>%
  dplyr::mutate(

    comparison_label = factor(
      comparison_key,
      levels = Fig3G_comparison_keys,
      labels = Fig3G_comparison_labels
    ),

    p_label = dplyr::case_when(
      is.na(permutation_p_value) ~
        NA_character_,

      permutation_p_value < 0.0001 ~
        "P < 0.0001",

      TRUE ~
        paste0(
          "P = ",
          trimws(
            base::format.pval(
              permutation_p_value,
              digits = 4
            )
          )
        )
    )
  )

if (
  nrow(Fig3G_GSE97681_violin_statistics) != 4 ||
  anyNA(
    Fig3G_GSE97681_violin_statistics$permutation_p_value
  ) ||
  anyNA(
    Fig3G_GSE97681_violin_statistics$comparison_label
  )
) {
  stop(
    "Panel G statistics could not be matched to all four violin plots."
  )
}

Fig3G_GSE97681_mean_labels <-
  dplyr::bind_rows(

    Fig3G_GSE97681_violin_statistics %>%
      dplyr::transmute(
        cell_line,
        comparison_label,
        x_position = 1,
        mean_label = sprintf(
          "%.3f",
          mean_MFE_LFC_below_zero
        )
      ),

    Fig3G_GSE97681_violin_statistics %>%
      dplyr::transmute(
        cell_line,
        comparison_label,
        x_position = 2,
        mean_label = sprintf(
          "%.3f",
          mean_MFE_LFC_above_zero
        )
      )
  )

Plot_Fig3G_GSE97681_violin <- ggplot(
  Fig3G_GSE97681_violin_plot_data,
  aes(
    x = as.numeric(expression_group),
    y = Normalized_MFE,
    group = expression_group
  )
) +
  geom_violin(
    aes(fill = comparison_label),
    color = "black",
    linewidth = 0.3,
    scale = "width",
    trim = TRUE
  ) +
  geom_boxplot(
    width = 0.10,
    fill = "white",
    color = "black",
    linewidth = 0.3,
    outlier.shape = NA
  ) +
  geom_text(
    data = Fig3G_GSE97681_violin_statistics,
    aes(
      x = 1.5,
      y = 0.060,
      label = p_label
    ),
    inherit.aes = FALSE,
    size = 3.2,
    color = "black"
  ) +
  geom_text(
    data = Fig3G_GSE97681_mean_labels,
    aes(
      x = x_position,
      y = 0.018,
      label = mean_label
    ),
    inherit.aes = FALSE,
    size = 3.0,
    color = "black"
  ) +
  scale_fill_manual(
    values = c(
      "Acute 48 h versus untreated" = "lightgrey",
      "Stable resistant (drug on) versus acute 48 h" =
        "grey50"
    )
  ) +
  facet_grid(
    rows = vars(cell_line),
    cols = vars(comparison_label),
    drop = FALSE
  ) +
  scale_x_continuous(
    breaks = c(1, 2),
    labels = c(
      "LFC < 0",
      "LFC > 0"
    ),
    limits = c(0.5, 2.5)
  ) +
  scale_y_continuous(
    breaks = seq(
      -0.6,
      0,
      by = 0.2
    ),
    expand = expansion(
      mult = c(0, 0)
    )
  ) +
  coord_cartesian(
    ylim = c(-0.6, 0.08),
    clip = "on"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",

    axis.text.x = element_text(
      size = 11,
      color = "black"
    ),

    axis.text.y = element_text(
      size = 11,
      color = "black"
    ),

    axis.title.y = element_text(
      size = 14
    ),

    axis.title.x = element_blank(),

    strip.background = element_blank(),

    strip.text = element_text(
      size = 11,
      color = "black"
    ),

    panel.spacing = grid::unit(
      0.8,
      "lines"
    ),

    plot.margin = margin(
      t = 8,
      r = 8,
      b = 8,
      l = 8
    )
  ) +
  labs(
    title = "",
    x = "",
    y = expression(
      Delta * G[fold] ~
        (kcal %.% mol^{-1} %.% nt^{-1})
    )
  )

print(Plot_Fig3G_GSE97681_violin)

print(
  as.data.frame(
    Fig3G_GSE97681_violin_statistics
  ),
  row.names = FALSE,
  digits = 10
)

ggsave(
  filename = "Fig3G_GSE97681_violin_plots.tiff",
  plot = Plot_Fig3G_GSE97681_violin,
  path = Fig3F_I_output_dir,
  device = "tiff",
  width = 6,
  height = 8,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

Fig3H_expected_profiles <- c(
  "WM989_Acute_vs_Untreated",
  "WM989_Resistant_vs_Acute",
  "WM983B_Acute_vs_Untreated",
  "WM983B_Resistant_vs_Acute"
)

Fig3H_required_columns <- c(
  "Samples",
  "reference_high_MFE_fraction",
  "comparison_high_MFE_fraction",
  "reference_low_MFE_fraction",
  "comparison_low_MFE_fraction"
)

Fig3H_missing_columns <- setdiff(
  Fig3H_required_columns,
  colnames(Fig3H_GSE97681_quartile_statistics)
)

if (length(Fig3H_missing_columns) > 0) {
  stop(
    paste0(
      "Panel H is missing the following columns:\n",
      paste(Fig3H_missing_columns, collapse = "\n")
    )
  )
}

if (
  nrow(Fig3H_GSE97681_quartile_statistics) != 4 ||
  !setequal(
    Fig3H_GSE97681_quartile_statistics$Samples,
    Fig3H_expected_profiles
  )
) {
  stop(
    "Panel H does not contain exactly the four required GSE97681 contrasts."
  )
}

Fig3H_sample_labels <- c(
  "WM989_Acute_vs_Untreated" =
    "WM989-A6\nAcute 48 h versus untreated",
  "WM989_Resistant_vs_Acute" =
    "WM989-A6\nStable resistant (drug on) versus acute 48 h",
  "WM983B_Acute_vs_Untreated" =
    "WM983B-E9\nAcute 48 h versus untreated",
  "WM983B_Resistant_vs_Acute" =
    "WM983B-E9\nStable resistant (drug on) versus acute 48 h"
)

Fig3H_GSE97681_plot_data <- dplyr::bind_rows(
  Fig3H_GSE97681_quartile_statistics %>%
    dplyr::transmute(
      Samples,
      MFE_group = "Transcripts with less-negative DeltaGfold values",
      Reference = reference_high_MFE_fraction,
      Comparison = comparison_high_MFE_fraction
    ) %>%
    tidyr::pivot_longer(
      cols = c(
        Reference,
        Comparison
      ),
      names_to = "Condition",
      values_to = "Value"
    ),
  Fig3H_GSE97681_quartile_statistics %>%
    dplyr::transmute(
      Samples,
      MFE_group = "Transcripts with more-negative DeltaGfold values",
      Reference = reference_low_MFE_fraction,
      Comparison = comparison_low_MFE_fraction
    ) %>%
    tidyr::pivot_longer(
      cols = c(
        Reference,
        Comparison
      ),
      names_to = "Condition",
      values_to = "Value"
    )
) %>%
  dplyr::mutate(
    Samples = factor(
      Samples,
      levels = Fig3H_expected_profiles
    ),
    MFE_group = factor(
      MFE_group,
      levels = c(
        "Transcripts with less-negative DeltaGfold values",
        "Transcripts with more-negative DeltaGfold values"
      ),
      labels = c(
        paste0(
          "Transcripts with less-negative ΔGfold values\n",
          "(>75th percentile)"
        ),
        paste0(
          "Transcripts with more-negative ΔGfold values\n",
          "(<25th percentile)"
        )
      )
    ),
    Condition = factor(
      Condition,
      levels = c(
        "Reference",
        "Comparison"
      ),
      labels = c(
        "Reference state",
        "Comparison state"
      )
    )
  )

if (
  any(!is.finite(Fig3H_GSE97681_plot_data$Value)) ||
  any(Fig3H_GSE97681_plot_data$Value < 0) ||
  any(Fig3H_GSE97681_plot_data$Value > 1)
) {
  stop("Panel H contains an invalid fraction.")
}

Fig3H_y_max <- max(
  0.4,
  ceiling(
    20 *
      1.35 *
      max(
        Fig3H_GSE97681_plot_data$Value,
        na.rm = TRUE
      )
  ) / 20
)

Plot_Fig3H_GSE97681_quartiles <- ggplot(
  Fig3H_GSE97681_plot_data,
  aes(
    x = Samples,
    y = Value,
    fill = Condition
  )
) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  geom_text(
    aes(
      label = sprintf(
        "%.3f",
        Value
      )
    ),
    position = position_dodge(width = 0.8),
    hjust = -0.30,
    vjust = 0.40,
    angle = 90,
    size = 3.3,
    color = "black"
  ) +
  facet_grid(
    cols = vars(MFE_group)
  ) +
  scale_x_discrete(
    drop = FALSE,
    labels = Fig3H_sample_labels
  ) +
  scale_y_continuous(
    limits = c(
      0,
      Fig3H_y_max
    ),
    expand = expansion(
      mult = c(0, 0)
    )
  ) +
  scale_fill_manual(
    values = c(
      "Reference state" = "#8cc5e3",
      "Comparison state" = "#2066a8"
    )
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 10,
      color = "black"
    ),
    axis.text.y = element_text(
      size = 11,
      color = "black"
    ),
    axis.title.y = element_text(
      size = 14
    ),
    strip.background = element_blank(),
    strip.text = element_text(
      size = 12,
      color = "black"
    ),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(
      size = 10
    ),
    panel.spacing.x = grid::unit(
      1,
      "lines"
    )
  ) +
  labs(
    title = "",
    x = "",
    y = paste(
      "Fraction of normalized counts",
      "among significant transcripts",
      sep = "\n"
    ),
    fill = ""
  )

print(Plot_Fig3H_GSE97681_quartiles)

print(
  as.data.frame(
    Fig3H_GSE97681_quartile_statistics
  ),
  row.names = FALSE,
  digits = 10
)

ggsave(
  filename = "Fig3H_GSE97681_quartile_fractions.tiff",
  plot = Plot_Fig3H_GSE97681_quartiles,
  path = Fig3F_I_output_dir,
  device = "tiff",
  width = 5,
  height = 6,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

Fig3H_perm_B <- 10000L
Fig3H_perm_seed <- 20260918L

Fig3H_perm_source <- get(
  "GSE97681_Fig3D_transcript_count_data",
  envir = .GlobalEnv,
  inherits = FALSE
)

Fig3H_perm_required <- c(
  "Samples", "GREIN_transcript_id", "Normalized_MFE",
  "mean_normalized_count_reference",
  "mean_normalized_count_comparison",
  "high_MFE_transcript", "low_MFE_transcript"
)

stopifnot(
  all(Fig3H_perm_required %in% names(Fig3H_perm_source)),
  length(Fig3H_expected_profiles) == 4L,
  !anyDuplicated(Fig3H_expected_profiles),
  !anyNA(Fig3H_perm_source$Samples),
  setequal(Fig3H_perm_source$Samples, Fig3H_expected_profiles)
)

Fig3H_prepare_permutation <- function(profile) {

  fail <- function(msg) stop(profile, ": ", msg)

  dat <- Fig3H_perm_source[
    as.character(Fig3H_perm_source$Samples) == profile, ]

  if (
    anyNA(dat[, Fig3H_perm_required]) ||
    anyDuplicated(dat$GREIN_transcript_id) ||
    any(!is.finite(dat$Normalized_MFE))
  ) {
    fail("Missing values, duplicated IDs or invalid folding predictions.")
  }

  upper <- dat$high_MFE_transcript
  lower <- dat$low_MFE_transcript

  if (
    !is.logical(upper) || !is.logical(lower) ||
    any(upper & lower) || !any(upper) || !any(lower)
  ) {
    fail("Invalid quartile memberships.")
  }

  category <- ifelse(upper, 1L, ifelse(lower, 2L, 0L))

  counts <- cbind(
    dat$mean_normalized_count_reference,
    dat$mean_normalized_count_comparison
  )

  if (
    any(!is.finite(counts)) || any(counts < 0) ||
    any(!is.finite(colSums(counts))) || any(colSums(counts) <= 0)
  ) {
    fail("Invalid normalized counts.")
  }

  weights <- sweep(counts, 2, colSums(counts), "/")

  fractions <- vapply(1:2, function(k) {
    colSums(weights[category == k, , drop = FALSE])
  }, numeric(2))

  old <- Fig3H_GSE97681_quartile_statistics[
    as.character(Fig3H_GSE97681_quartile_statistics$Samples) == profile, ]

  if (nrow(old) != 1L) {
    fail("Expected one existing summary row.")
  }

  expected <- c(
    old$reference_high_MFE_fraction,
    old$comparison_high_MFE_fraction,
    old$reference_low_MFE_fraction,
    old$comparison_low_MFE_fraction
  )

  if (
    length(expected) != 4L ||
    any(!is.finite(expected)) ||
    max(abs(as.vector(fractions) - expected)) > 1e-10
  ) {
    fail("Recalculated fractions do not match the existing bars.")
  }

  list(
    Samples = profile,
    category = category,
    fractions = fractions,
    delta = weights[, 2] - weights[, 1]
  )
}

Fig3H_perm_inputs <- lapply(
  Fig3H_expected_profiles,
  Fig3H_prepare_permutation
)

message("All Fig. 3H fractions match the existing bars.")

set.seed(Fig3H_perm_seed)

Fig3H_perm_runs <- lapply(Fig3H_perm_inputs, function(x) {

  message("Permutations: ", x$Samples)

  observed <- x$fractions[2, ] - x$fractions[1, ]

  null <- replicate(Fig3H_perm_B, {

    shuffled <- x$category[sample.int(length(x$category))]
    c(
      sum(x$delta[shuffled == 1L]),
      sum(x$delta[shuffled == 2L])
    )
  })

  b <- rowSums(abs(null) >= abs(observed) - 1e-12)

  list(
    results = data.frame(
      Samples = x$Samples,
      Quartile = c("Upper / less negative", "Lower / more negative"),
      Fraction_reference = x$fractions[1, ],
      Fraction_comparison = x$fractions[2, ],
      Delta_fraction = observed,
      N_transcripts = length(x$category),
      N_quartile = c(sum(x$category == 1L), sum(x$category == 2L)),
      P_permutation = (b + 1) / (Fig3H_perm_B + 1),
      stringsAsFactors = FALSE
    ),
    null_delta = null
  )
})

Fig3H_quartile_permutation_results <- do.call(
  rbind,
  lapply(Fig3H_perm_runs, function(x) x$results)
)

stopifnot(nrow(Fig3H_quartile_permutation_results) == 8L)

Fig3H_quartile_permutation_results$P_adj_BH <- p.adjust(
  Fig3H_quartile_permutation_results$P_permutation,
  method = "BH"
)

Fig3H_quartile_permutation_results$Figure_label <- ifelse(
  Fig3H_quartile_permutation_results$P_adj_BH >= 0.05,
  "NS",
  sprintf("P = %.4f", Fig3H_quartile_permutation_results$P_adj_BH)
)

print(
  Fig3H_quartile_permutation_results[, c(
    "Samples", "Quartile", "Fraction_reference", "Fraction_comparison",
    "Delta_fraction", "P_permutation", "P_adj_BH", "Figure_label"
  )],
  row.names = FALSE,
  digits = 7
)

Fig3I_required_contrasts <- c(
  "WM989_Acute_minus_Untreated",
  "WM989_Resistant_minus_Acute",
  "WM983B_Acute_minus_Untreated",
  "WM983B_Resistant_minus_Acute"
)

Fig3I_required_columns <- c(
  "contrast_id",
  "cell_line",
  "delta_TPM_weighted_Normalized_MFE",
  "propagated_sem",
  "permutation_p_value",
  "permutation_p_adj_BH"
)

Fig3I_missing_columns <- setdiff(
  Fig3I_required_columns,
  colnames(Fig3I_GSE97681_DeltaDeltaGfold_statistics)
)

if (length(Fig3I_missing_columns) > 0) {
  stop(
    paste0(
      "Panel I is missing the following columns:\n",
      paste(Fig3I_missing_columns, collapse = "\n")
    )
  )
}

Fig3I_GSE97681_DeltaDeltaGfold_statistics <-
  Fig3I_GSE97681_DeltaDeltaGfold_statistics %>%
  dplyr::filter(
    contrast_id %in% Fig3I_required_contrasts
  )

if (
  nrow(Fig3I_GSE97681_DeltaDeltaGfold_statistics) != 4 ||
  !setequal(
    Fig3I_GSE97681_DeltaDeltaGfold_statistics$contrast_id,
    Fig3I_required_contrasts
  ) ||
  any(
    grepl(
      "Resistant_minus_Untreated",
      Fig3I_GSE97681_DeltaDeltaGfold_statistics$contrast_id
    )
  )
) {
  stop(
    "Panel I does not contain exactly the four required sequential contrasts."
  )
}

Fig3I_GSE97681_DeltaDeltaGfold_statistics <-
  Fig3I_GSE97681_DeltaDeltaGfold_statistics %>%
  dplyr::mutate(
    cell_line = factor(
      cell_line,
      levels = c(
        "WM989-A6",
        "WM983B-E9"
      )
    ),
    comparison = dplyr::case_when(
      contrast_id %in% c(
        "WM989_Acute_minus_Untreated",
        "WM983B_Acute_minus_Untreated"
      ) ~
        "Acute 48 h versus untreated",
      contrast_id %in% c(
        "WM989_Resistant_minus_Acute",
        "WM983B_Resistant_minus_Acute"
      ) ~
        "Stable resistant (drug on) versus acute 48 h",
      TRUE ~
        NA_character_
    ),
    comparison = factor(
      comparison,
      levels = c(
        "Acute 48 h versus untreated",
        "Stable resistant (drug on) versus acute 48 h"
      )
    ),
    errorbar_low =
      delta_TPM_weighted_Normalized_MFE -
      propagated_sem,
    errorbar_high =
      delta_TPM_weighted_Normalized_MFE +
      propagated_sem,
    plot_errorbar_low = dplyr::if_else(
      delta_TPM_weighted_Normalized_MFE >= 0,
      delta_TPM_weighted_Normalized_MFE,
      errorbar_low
    ),
    plot_errorbar_high = dplyr::if_else(
      delta_TPM_weighted_Normalized_MFE >= 0,
      errorbar_high,
      delta_TPM_weighted_Normalized_MFE
    ),
    delta_label = dplyr::case_when(
      abs(
        delta_TPM_weighted_Normalized_MFE
      ) < 0.001 ~
        sprintf(
          "%.5f",
          delta_TPM_weighted_Normalized_MFE
        ),
      TRUE ~
        sprintf(
          "%.4f",
          delta_TPM_weighted_Normalized_MFE
        )
    ),
    p_label = dplyr::case_when(
      permutation_p_adj_BH < 0.001 ~
        "BH-adjusted P < 0.001",
      permutation_p_adj_BH >= 0.05 ~
        "NS",
      TRUE ~ paste0(
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
    )
  ) %>%
  dplyr::arrange(
    cell_line,
    comparison
  )

if (
  anyNA(Fig3I_GSE97681_DeltaDeltaGfold_statistics$cell_line) ||
  anyNA(Fig3I_GSE97681_DeltaDeltaGfold_statistics$comparison) ||
  any(
    !is.finite(
      as.matrix(
        Fig3I_GSE97681_DeltaDeltaGfold_statistics[
          c(
            "delta_TPM_weighted_Normalized_MFE",
            "propagated_sem",
            "permutation_p_value",
            "permutation_p_adj_BH"
          )
        ]
      )
    )
  ) ||
  any(
    Fig3I_GSE97681_DeltaDeltaGfold_statistics[[
      "propagated_sem"
    ]] < 0
  ) ||
  any(
    Fig3I_GSE97681_DeltaDeltaGfold_statistics[[
      "permutation_p_value"
    ]] < 0 |
    Fig3I_GSE97681_DeltaDeltaGfold_statistics[[
      "permutation_p_value"
    ]] > 1 |
    Fig3I_GSE97681_DeltaDeltaGfold_statistics[[
      "permutation_p_adj_BH"
    ]] < 0 |
    Fig3I_GSE97681_DeltaDeltaGfold_statistics[[
      "permutation_p_adj_BH"
    ]] > 1
  )
) {
  stop(
    "Panel I contains an invalid cell line, comparison or statistic."
  )
}

Fig3I_y_limit <- max(
  0.02,
  1.45 *
    max(
      abs(
        c(
          Fig3I_GSE97681_DeltaDeltaGfold_statistics$errorbar_low,
          Fig3I_GSE97681_DeltaDeltaGfold_statistics$errorbar_high
        )
      ),
      na.rm = TRUE
    )
)

Fig3I_GSE97681_DeltaDeltaGfold_statistics <-
  Fig3I_GSE97681_DeltaDeltaGfold_statistics %>%
  dplyr::mutate(
    delta_label_y = dplyr::if_else(
      delta_TPM_weighted_Normalized_MFE >= 0,
      errorbar_high + 0.025 * Fig3I_y_limit,
      errorbar_low - 0.025 * Fig3I_y_limit
    ),
    p_label_y = dplyr::if_else(
      delta_TPM_weighted_Normalized_MFE >= 0,
      errorbar_high + 0.18 * Fig3I_y_limit,
      errorbar_low - 0.18 * Fig3I_y_limit
    ),
    text_hjust = dplyr::if_else(
      delta_TPM_weighted_Normalized_MFE >= 0,
      0,
      1
    )
  )

Plot_Fig3I_GSE97681_DeltaDeltaGfold <- ggplot(
  Fig3I_GSE97681_DeltaDeltaGfold_statistics,
  aes(
    x = comparison,
    y = delta_TPM_weighted_Normalized_MFE,
    fill = comparison
  )
) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.4,
    linetype = "dashed",
    color = "black"
  ) +
  geom_col(
    width = 0.85,
    color = "black",
    linewidth = 0.3
  ) +
  geom_errorbar(
    aes(
      ymin = plot_errorbar_low,
      ymax = plot_errorbar_high
    ),
    width = 0.15,
    linewidth = 0.5,
    color = "black"
  ) +
  geom_text(
    aes(
      y = delta_label_y,
      label = delta_label,
      hjust = text_hjust
    ),
    angle = 90,
    size = 2.8,
    color = "black"
  ) +
  facet_grid(
    cols = vars(cell_line)
  ) +
  scale_fill_manual(
    values = c(
      "Acute 48 h versus untreated" = "#D55E00",
      "Stable resistant (drug on) versus acute 48 h" = "#0072B2"
    )
  ) +
  scale_x_discrete(
    labels = c(
      "Acute 48 h versus untreated" =
        "Acute 48 h\nversus untreated",
      "Stable resistant (drug on) versus acute 48 h" =
        "Stable resistant (drug on)\nversus acute 48 h"
    )
  ) +
  scale_y_continuous(
    labels = scales::number_format(
      accuracy = 0.001
    ),
    expand = expansion(
      mult = c(0.02, 0.02)
    )
  ) +
  coord_cartesian(
    ylim = c(
      -Fig3I_y_limit,
      Fig3I_y_limit
    ),
    clip = "off"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      size = 10,
      color = "black"
    ),
    axis.text.y = element_text(
      size = 11,
      color = "black"
    ),
    axis.title.y = element_text(
      size = 14
    ),
    strip.background = element_blank(),
    strip.text = element_text(
      size = 12,
      color = "black"
    ),
    plot.margin = margin(
      t = 8,
      r = 8,
      b = 8,
      l = 8
    )
  ) +
  labs(
    title = "",
    x = "",
    y = expression(
      Delta * Delta * G[fold] ~
        (kcal %.% mol^{-1} %.% nt^{-1})
    )
  )

print(Plot_Fig3I_GSE97681_DeltaDeltaGfold)

print(
  as.data.frame(
    Fig3I_GSE97681_DeltaDeltaGfold_statistics
  ),
  row.names = FALSE,
  digits = 12
)

ggsave(
  filename = "Fig3I_GSE97681_TPM_weighted_DeltaDeltaGfold.tiff",
  plot = Plot_Fig3I_GSE97681_DeltaDeltaGfold,
  path = Fig3F_I_output_dir,
  device = "tiff",
  width = 4,
  height = 6,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

Fig3I_GSE97681_paired_statistics <-
  get(
    "GSE97681_founder_balanced_trajectory_wide",
    envir = .GlobalEnv,
    inherits = FALSE
  ) %>%
  tibble::as_tibble() %>%
  dplyr::filter(
    is.finite(.data$Resistant_minus_Acute)
  ) %>%
  dplyr::group_by(.data$cell_line) %>%
  dplyr::summarise(
    N_founders = dplyr::n(),
    mean_Resistant_minus_Acute =
      mean(.data$Resistant_minus_Acute),
    paired_one_sided_P = stats::t.test(
      .data$Resistant_minus_Acute,
      mu = 0,
      alternative = "less"
    )$p.value,
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    paired_BH_P = stats::p.adjust(
      .data$paired_one_sided_P,
      method = "BH"
    )
  )

print(
  as.data.frame(Fig3I_GSE97681_paired_statistics),
  row.names = FALSE,
  digits = 12
)
