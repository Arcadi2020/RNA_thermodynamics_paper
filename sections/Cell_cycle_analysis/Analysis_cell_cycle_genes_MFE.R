suppressPackageStartupMessages({
  library(dplyr)
  library(biomaRt)
  library(clusterProfiler)
  library(DESeq2)
})

gencode_mfe_file <- "../../data/Gencode_TX_norm_MFE_HC_table.csv"
mcf7_differential_expression_rdata_file <- "path_to/Analysis_RNAseq_PNAS_2017_environment.RData"
mcf7_transcript_tpm_file <- "path_to/GSE94479_TranscriptLevel_TPM_data.csv"
mcf7_transcript_raw_counts_file <- "path_to/GSE94479_TranscriptLevel_Raw_data.csv"
hela_transcript_raw_counts_file <- "path_to/GSE104736_TranscriptLevel_Raw_data.csv"
hela_transcript_tpm_file <- "path_to/GSE104736_TranscriptLevel_TPM_data.csv"
hela_gene_raw_counts_file <- "path_to/GSE104736_GeneLevel_Raw_data.csv"
analysis_output_dir <- "path_to/cell_cycle_analysis_outputs"
figure2_results_file <- file.path(
  analysis_output_dir,
  "cell_cycle_analysis_results_for_Fig_2.rds"
)

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

input_files <- c(
  gencode_mfe_file,
  mcf7_differential_expression_rdata_file,
  mcf7_transcript_tpm_file,
  mcf7_transcript_raw_counts_file,
  hela_transcript_raw_counts_file,
  hela_transcript_tpm_file,
  hela_gene_raw_counts_file
)

for (input_file in input_files) {
  stop_if_missing(input_file, "Input file")
}

dir.create(analysis_output_dir, recursive = TRUE, showWarnings = FALSE)

Gencode_TX_norm_MFE_HC_table <- read.delim(
  gencode_mfe_file,
  sep = ",",
  header = TRUE,
  stringsAsFactors = FALSE
)
Table_ANNOTATED <- Gencode_TX_norm_MFE_HC_table

mcf7_differential_expression_objects <- load_rdata_objects(
  mcf7_differential_expression_rdata_file,
  c(
    "GSE94479_Deseq_G0G1vsG2M_IDs",
    "GSE94479_Deseq_SvsG0G1_IDs",
    "GSE94479_Deseq_G2MvsS_IDs"
  )
)

GSE94479_Deseq_G0G1vsG2M_IDs <-
  mcf7_differential_expression_objects[["GSE94479_Deseq_G0G1vsG2M_IDs"]]
GSE94479_Deseq_SvsG0G1_IDs <-
  mcf7_differential_expression_objects[["GSE94479_Deseq_SvsG0G1_IDs"]]
GSE94479_Deseq_G2MvsS_IDs <-
  mcf7_differential_expression_objects[["GSE94479_Deseq_G2MvsS_IDs"]]

GSE94479_G0G1vsG2M <- GSE94479_Deseq_G0G1vsG2M_IDs %>%
  dplyr::filter(padj < 0.05) %>%
  dplyr::left_join(Table_ANNOTATED, by = "ensembl_transcript_id") %>%
  dplyr::filter(!is.na(Normalized_MFE))

GSE94479_SvsG0G1 <- GSE94479_Deseq_SvsG0G1_IDs %>%
  dplyr::filter(padj < 0.05) %>%
  dplyr::left_join(Table_ANNOTATED, by = "ensembl_transcript_id") %>%
  dplyr::filter(!is.na(Normalized_MFE))

aggregate_transcript_results_to_gene <- function(differential_expression_table) {
  differential_expression_table %>%
    dplyr::filter(
      !is.na(log2FoldChange),
      !is.na(baseMean),
      baseMean > 0,
      !is.na(ensembl_gene_id)
    ) %>%
    dplyr::group_by(ensembl_gene_id) %>%
    dplyr::summarise(
      weighted_log2FC = sum(log2FoldChange * baseMean) / sum(baseMean),
      min_pvalue = min(pvalue, na.rm = TRUE),
      mean_baseMean = mean(baseMean, na.rm = TRUE),
      n_transcripts = dplyr::n(),
      .groups = "drop"
    )
}

map_ensembl_gene_to_entrez <- function(ensembl_gene_ids, ensembl_mart) {
  biomaRt::getBM(
    attributes = c("entrezgene_id", "ensembl_gene_id"),
    filters = "ensembl_gene_id",
    values = ensembl_gene_ids,
    mart = ensembl_mart
  ) %>%
    dplyr::distinct(ensembl_gene_id, .keep_all = TRUE) %>%
    dplyr::filter(!is.na(entrezgene_id))
}

select_one_gene_per_entrez <- function(expression_table) {
  expression_table %>%
    dplyr::group_by(entrezgene_id) %>%
    dplyr::slice_max(
      order_by = abs(weighted_log2FC),
      n = 1,
      with_ties = FALSE
    ) %>%
    dplyr::ungroup()
}

run_kegg_gsea <- function(data_frame, log_fc_columns, analysis_labels) {
  if (length(log_fc_columns) != length(analysis_labels)) {
    stop("GSEA column and label counts do not match.", call. = FALSE)
  }

  gsea_results <- list()

  for (analysis_index in seq_along(log_fc_columns)) {
    log_fc_column <- log_fc_columns[[analysis_index]]
    analysis_label <- analysis_labels[[analysis_index]]

    ranked_table <- data_frame %>%
      dplyr::select(entrezgene_id, !!rlang::sym(log_fc_column)) %>%
      dplyr::filter(
        !is.na(entrezgene_id),
        !is.na(!!rlang::sym(log_fc_column))
      )

    ranked_values <- ranked_table[[log_fc_column]]
    names(ranked_values) <- ranked_table$entrezgene_id
    ranked_values <- ranked_values[!duplicated(names(ranked_values))]
    ranked_values <- sort(ranked_values, decreasing = TRUE)

    if (length(ranked_values) == 0) {
      stop("No ranked genes available for ", analysis_label, ".", call. = FALSE)
    }

    gsea_result <- clusterProfiler::gseKEGG(
      ranked_values,
      organism = "hsa",
      keyType = "kegg",
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      minGSSize = 10,
      maxGSSize = 500
    )

    if (!is.null(gsea_result) && nrow(gsea_result) > 0) {
      gsea_results[[analysis_label]] <- gsea_result
    }
  }

  if (length(gsea_results) == 0) {
    stop("No significant KEGG GSEA results were returned.", call. = FALSE)
  }

  do.call(
    rbind,
    lapply(names(gsea_results), function(analysis_label) {
      result_table <- as.data.frame(gsea_results[[analysis_label]])
      result_table$Sample <- analysis_label
      result_table
    })
  )
}

GSE94479_Deseq_G0G1vsG2M_IDs_gene_level_df <-
  aggregate_transcript_results_to_gene(GSE94479_Deseq_G0G1vsG2M_IDs)
GSE94479_Deseq_SvsG0G1_IDs_gene_level_df <-
  aggregate_transcript_results_to_gene(GSE94479_Deseq_SvsG0G1_IDs)
GSE94479_Deseq_G2MvsS_IDs_gene_level_df <-
  aggregate_transcript_results_to_gene(GSE94479_Deseq_G2MvsS_IDs)

mart <- biomaRt::useMart(
  biomart = "ensembl",
  dataset = "hsapiens_gene_ensembl",
  host = "https://useast.ensembl.org"
)

GSE94479_Deseq_G0G1vsG2M_entrez_IDs <-
  GSE94479_Deseq_G0G1vsG2M_IDs_gene_level_df %>%
  dplyr::left_join(
    map_ensembl_gene_to_entrez(
      GSE94479_Deseq_G0G1vsG2M_IDs_gene_level_df$ensembl_gene_id,
      mart
    ),
    by = "ensembl_gene_id"
  ) %>%
  select_one_gene_per_entrez()

GSE94479_Deseq_SvsG0G1_entrez_IDs <-
  GSE94479_Deseq_SvsG0G1_IDs_gene_level_df %>%
  dplyr::left_join(
    map_ensembl_gene_to_entrez(
      GSE94479_Deseq_SvsG0G1_IDs_gene_level_df$ensembl_gene_id,
      mart
    ),
    by = "ensembl_gene_id"
  ) %>%
  select_one_gene_per_entrez()

GSE94479_Deseq_G2MvsS_entrez_IDs <-
  GSE94479_Deseq_G2MvsS_IDs_gene_level_df %>%
  dplyr::left_join(
    map_ensembl_gene_to_entrez(
      GSE94479_Deseq_G2MvsS_IDs_gene_level_df$ensembl_gene_id,
      mart
    ),
    by = "ensembl_gene_id"
  ) %>%
  select_one_gene_per_entrez()

GSE94479_gene_level_Log2_FC_df <-
  GSE94479_Deseq_G0G1vsG2M_entrez_IDs %>%
  dplyr::select(ensembl_gene_id, entrezgene_id, weighted_log2FC) %>%
  dplyr::rename(G0G1vsG2M_Log_FC = weighted_log2FC) %>%
  dplyr::left_join(
    GSE94479_Deseq_SvsG0G1_entrez_IDs %>%
      dplyr::select(ensembl_gene_id, entrezgene_id, weighted_log2FC) %>%
      dplyr::rename(SvsG0G1_Log_FC = weighted_log2FC),
    by = c("ensembl_gene_id", "entrezgene_id")
  ) %>%
  dplyr::left_join(
    GSE94479_Deseq_G2MvsS_entrez_IDs %>%
      dplyr::select(ensembl_gene_id, entrezgene_id, weighted_log2FC) %>%
      dplyr::rename(G2MvsS_Log_FC = weighted_log2FC),
    by = c("ensembl_gene_id", "entrezgene_id")
  )

combined_kk2_df <- run_kegg_gsea(
  GSE94479_gene_level_Log2_FC_df,
  c("G0G1vsG2M_Log_FC", "SvsG0G1_Log_FC", "G2MvsS_Log_FC"),
  c("GSE94479_G0G1vsG2M", "GSE94479_SvsG0G1", "GSE94479_G2MvsS")
)

Selected_KEGG_pathways <- c(
  "hsa04150", "hsa04668", "hsa04657", "hsa03030", "hsa03410",
  "hsa03430", "hsa03420", "hsa04370", "hsa4072", "hsa03320",
  "hsa03460", "hsa03440", "hsa04014", "hsa04070", "hsa04020",
  "hsa01230", "hsa04310", "hsa03082", "hsa04015", "hsa04625",
  "hsa04371", "hsa04390", "hsa04010", "hsa04622", "hsa04350",
  "hsa04152", "hsa04926", "hsa04630", "hsa04115", "hsa04151",
  "hsa04910", "hsa04330", "hsa04012", "hsa03440", "hsa03460",
  "hsa05230"
)

combined_kk2_df_filtered <- combined_kk2_df %>%
  dplyr::filter(ID %in% Selected_KEGG_pathways) %>%
  dplyr::mutate(NES = scales::rescale(NES, to = c(-4, 4)))

desired_order <- c(
  "Wnt signaling pathway",
  "Ras signaling pathway",
  "PI3K-Akt signaling pathway",
  "MAPK signaling pathway",
  "Rap1 signaling pathway",
  "VEGF signaling pathway",
  "Calcium signaling pathway",
  "p53 signaling pathway",
  "DNA replication",
  "Homologous recombination",
  "Fanconi anemia pathway",
  "Central carbon metabolism in cancer"
)

GSE94479_combined_kk2_df_filtered_plot <- combined_kk2_df_filtered %>%
  dplyr::mutate(
    Description = factor(Description, levels = desired_order),
    neg_log10_p_adjust = -log10(p.adjust),
    Sample = factor(
      Sample,
      levels = c(
        "GSE94479_G0G1vsG2M",
        "GSE94479_SvsG0G1",
        "GSE94479_G2MvsS"
      )
    )
  )

GSE94479_TPM_file <- mcf7_transcript_tpm_file

GSE94479_MCF7_TPM_required_samples <- c(
  "GSM2476255",
  "GSM2476256",
  "GSM2476257"
)

GSE94479_MCF7_TPM_sample_metadata <- tibble::tibble(
  Sample = GSE94479_MCF7_TPM_required_samples,
  Phase = factor(
    c("G1", "S", "G2/M"),
    levels = c("G1", "S", "G2/M")
  )
)

if (!file.exists(GSE94479_TPM_file)) {
  stop("The GSE94479 transcript-level TPM file was not found.")
}

GSE94479_MCF7_TranscriptLevel_TPM <- read.csv(
  GSE94479_TPM_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

colnames(GSE94479_MCF7_TranscriptLevel_TPM)[1] <-
  "GREIN_transcript_id"

GSE94479_MCF7_TPM_missing_samples <- setdiff(
  GSE94479_MCF7_TPM_required_samples,
  colnames(GSE94479_MCF7_TranscriptLevel_TPM)
)

if (length(GSE94479_MCF7_TPM_missing_samples) > 0) {
  stop(
    paste0(
      "The following GSE94479 TPM samples are missing: ",
      paste(
        GSE94479_MCF7_TPM_missing_samples,
        collapse = ", "
      )
    )
  )
}

GSE94479_MCF7_TranscriptLevel_TPM <-
  GSE94479_MCF7_TranscriptLevel_TPM %>%
  dplyr::mutate(
    GREIN_transcript_id =
      as.character(GREIN_transcript_id),

    ensembl_transcript_id = sub(
      "\\..*$",
      "",
      GREIN_transcript_id
    ),

    dplyr::across(
      dplyr::all_of(
        GSE94479_MCF7_TPM_required_samples
      ),
      as.numeric
    )
  ) %>%
  dplyr::select(
    GREIN_transcript_id,
    ensembl_transcript_id,
    dplyr::all_of(
      GSE94479_MCF7_TPM_required_samples
    )
  )

if (
  anyDuplicated(
    GSE94479_MCF7_TranscriptLevel_TPM$
    ensembl_transcript_id
  ) > 0
) {
  stop(
    paste(
      "Duplicated versionless transcript IDs were found",
      "in the GSE94479 TPM file."
    )
  )
}

GSE94479_MCF7_TPM_values_are_valid <- vapply(
  GSE94479_MCF7_TPM_required_samples,

  function(current_sample) {

    current_TPM <- getElement(
      GSE94479_MCF7_TranscriptLevel_TPM,
      current_sample
    )

    all(
      is.finite(current_TPM) &
        current_TPM >= 0
    )
  },

  logical(1)
)

if (!all(GSE94479_MCF7_TPM_values_are_valid)) {
  stop(
    paste(
      "At least one GSE94479 TPM sample contains",
      "missing, non-finite, or negative values."
    )
  )
}

if (
  !all(
    c(
      "ensembl_transcript_id",
      "Normalized_MFE"
    ) %in% colnames(Table_ANNOTATED)
  )
) {
  stop(
    paste(
      "Table_ANNOTATED lacks ensembl_transcript_id",
      "or Normalized_MFE."
    )
  )
}

Table_ANNOTATED_for_GSE94479_TPM <-
  Table_ANNOTATED %>%
  dplyr::transmute(
    ensembl_transcript_id = sub(
      "\\..*$",
      "",
      as.character(ensembl_transcript_id)
    ),

    Normalized_MFE =
      as.numeric(Normalized_MFE)
  ) %>%
  dplyr::filter(
    !is.na(ensembl_transcript_id),
    ensembl_transcript_id != "",
    is.finite(Normalized_MFE)
  ) %>%
  dplyr::distinct(
    ensembl_transcript_id,
    .keep_all = TRUE
  )

GSE94479_MCF7_TPM_thermodynamic <-
  GSE94479_MCF7_TranscriptLevel_TPM %>%
  dplyr::inner_join(
    Table_ANNOTATED_for_GSE94479_TPM,
    by = "ensembl_transcript_id"
  )

if (
  nrow(GSE94479_MCF7_TPM_thermodynamic) == 0 ||
  anyDuplicated(
    GSE94479_MCF7_TPM_thermodynamic$
    ensembl_transcript_id
  ) > 0
) {
  stop(
    paste(
      "The GSE94479 TPM-to-annotation join is empty",
      "or generated duplicated transcripts."
    )
  )
}

GSE94479_MCF7_TPM_original_totals <-
  tibble::tibble(
    Sample =
      GSE94479_MCF7_TPM_required_samples,

    Original_total_TPM = vapply(
      GSE94479_MCF7_TPM_required_samples,

      function(current_sample) {
        sum(
          getElement(
            GSE94479_MCF7_TranscriptLevel_TPM,
            current_sample
          )
        )
      },

      numeric(1)
    )
  )

GSE94479_MCF7_TPM_long <-
  GSE94479_MCF7_TPM_thermodynamic %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(
      GSE94479_MCF7_TPM_required_samples
    ),
    names_to = "Sample",
    values_to = "TPM"
  ) %>%
  dplyr::left_join(
    GSE94479_MCF7_TPM_sample_metadata,
    by = "Sample"
  ) %>%
  dplyr::mutate(
    Dataset = "GSE94479",
    Cell_line = "MCF7",
    TPM = as.numeric(TPM)
  ) %>%
  dplyr::group_by(
    Sample,
    Phase
  ) %>%
  dplyr::mutate(
    retained_TPM_MFE =
      sum(TPM),

    relative_transcript_abundance_MFE =
      TPM / retained_TPM_MFE,

    TPM_weighted_DGfold_contribution =
      relative_transcript_abundance_MFE *
      Normalized_MFE
  ) %>%
  dplyr::ungroup()

GSE94479_MCF7_TPM_weight_validation <-
  GSE94479_MCF7_TPM_long %>%
  dplyr::group_by(
    Dataset,
    Cell_line,
    Sample,
    Phase
  ) %>%
  dplyr::summarise(
    annotated_transcripts =
      dplyr::n_distinct(
        ensembl_transcript_id
      ),

    expressed_annotated_transcripts =
      dplyr::n_distinct(
        ensembl_transcript_id[TPM > 0]
      ),

    retained_TPM_MFE =
      dplyr::first(retained_TPM_MFE),

    sum_MFE_weights = sum(
      relative_transcript_abundance_MFE,
      na.rm = TRUE
    ),

    TPM_weighted_Normalized_MFE = sum(
      TPM_weighted_DGfold_contribution,
      na.rm = TRUE
    ),

    .groups = "drop"
  ) %>%
  dplyr::left_join(
    GSE94479_MCF7_TPM_original_totals,
    by = "Sample"
  ) %>%
  dplyr::mutate(
    MFE_coverage_percent =
      100 *
      retained_TPM_MFE /
      Original_total_TPM,

    abundance_measure =
      "Transcript-level TPM"
  ) %>%
  dplyr::select(
    Dataset,
    Cell_line,
    Sample,
    Phase,
    annotated_transcripts,
    expressed_annotated_transcripts,
    Original_total_TPM,
    retained_TPM_MFE,
    MFE_coverage_percent,
    sum_MFE_weights,
    TPM_weighted_Normalized_MFE,
    abundance_measure
  )

if (
  any(
    abs(
      GSE94479_MCF7_TPM_weight_validation$
      sum_MFE_weights - 1
    ) > 1e-8
  )
) {
  stop(
    paste(
      "The GSE94479 TPM-derived ΔGfold weights",
      "do not sum to one."
    )
  )
}

cat(
  "\nGSE94479 transcriptome-wide TPM-weighted ΔGfold:\n"
)

print(
  as.data.frame(
    GSE94479_MCF7_TPM_weight_validation
  ),
  row.names = FALSE,
  digits = 12
)

Counts_GSE104736_GREIN <- read.delim(
  hela_transcript_raw_counts_file,
  sep = ",",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
) %>%
  tibble::column_to_rownames(var = "X") %>%
  as.matrix()

Expdata_GSE104736 <- data.frame(
  condition = rep(c("HeLa_G1", "HeLa_S", "HeLa_G2"), each = 3),
  Cell_cycle_phase = rep(c("G1", "S", "G2"), each = 3),
  Cell_line = rep("HeLa", 9),
  row.names = c(
    "GSM2806897", "GSM2806898", "GSM2806899",
    "GSM2806900", "GSM2806901", "GSM2806902",
    "GSM2806903", "GSM2806904", "GSM2806905"
  ),
  stringsAsFactors = FALSE
)

missing_hela_transcript_samples <- setdiff(
  rownames(Expdata_GSE104736),
  colnames(Counts_GSE104736_GREIN)
)

if (length(missing_hela_transcript_samples) > 0) {
  stop(
    "Missing HeLa transcript-count samples: ",
    paste(missing_hela_transcript_samples, collapse = ", "),
    call. = FALSE
  )
}

Counts_GSE104736_GREIN <- floor(
  Counts_GSE104736_GREIN[
    ,
    rownames(Expdata_GSE104736),
    drop = FALSE
  ]
)

stopifnot(
  identical(
    colnames(Counts_GSE104736_GREIN),
    rownames(Expdata_GSE104736)
  )
)

dds_GSE104736 <- DESeq2::DESeqDataSetFromMatrix(
  countData = Counts_GSE104736_GREIN,
  colData = Expdata_GSE104736,
  design = ~condition
)

dds_GSE104736 <- dds_GSE104736[
  rowSums(DESeq2::counts(dds_GSE104736)) >= 3,
]
dds_GSE104736 <- DESeq2::DESeq(dds_GSE104736)

res_GSE104736_G1_vs_G2 <- DESeq2::results(
  dds_GSE104736,
  contrast = c("condition", "HeLa_G1", "HeLa_G2")
) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("ensembl_transcript_id")

MFE_GSE104736_G1_vs_G2 <- Table_ANNOTATED %>%
  dplyr::inner_join(
    res_GSE104736_G1_vs_G2,
    by = "ensembl_transcript_id"
  ) %>%
  dplyr::filter(padj < 0.05)

res_GSE104736_S_vs_G1 <- DESeq2::results(
  dds_GSE104736,
  contrast = c("condition", "HeLa_S", "HeLa_G1")
) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("ensembl_transcript_id")

MFE_GSE104736_S_vs_G1 <- Table_ANNOTATED %>%
  dplyr::inner_join(
    res_GSE104736_S_vs_G1,
    by = "ensembl_transcript_id"
  ) %>%
  dplyr::filter(
    padj < 0.05,
    log2FoldChange < 15,
    log2FoldChange > -15
  )

GSE104736_HeLa_TX_norm_counts <- read.delim(
  hela_transcript_tpm_file,
  sep = ",",
  header = TRUE,
  stringsAsFactors = FALSE
)

GSE104736_HeLa_TX_norm_counts <- dplyr::select(GSE104736_HeLa_TX_norm_counts,
                                               X,
                                               GSM2806897,
                                               GSM2806898,
                                               GSM2806899,
                                               GSM2806900,
                                               GSM2806901,
                                               GSM2806902,
                                               GSM2806903,
                                               GSM2806904,
                                               GSM2806905) %>%
  dplyr::rename(ensembl_transcript_id = X)

GSE104736_HeLa_TranscriptLevel_TPM <- GSE104736_HeLa_TX_norm_counts

GSE104736_HeLa_TPM_required_samples <- c(
  "GSM2806897",
  "GSM2806898",
  "GSM2806899",
  "GSM2806900",
  "GSM2806901",
  "GSM2806902",
  "GSM2806903",
  "GSM2806904",
  "GSM2806905"
)

GSE104736_HeLa_TPM_sample_metadata <-
  tibble::tibble(
    Sample =
      GSE104736_HeLa_TPM_required_samples,

    Phase = factor(
      c(
        "G1", "G1", "G1",
        "S", "S", "S",
        "G2", "G2", "G2"
      ),
      levels = c(
        "G1",
        "S",
        "G2"
      )
    )
  )

GSE104736_HeLa_TPM_missing_samples <- setdiff(
  GSE104736_HeLa_TPM_required_samples,
  colnames(
    GSE104736_HeLa_TranscriptLevel_TPM
  )
)

if (
  length(
    GSE104736_HeLa_TPM_missing_samples
  ) > 0
) {
  stop(
    paste0(
      "The following HeLa TPM samples are missing: ",
      paste(
        GSE104736_HeLa_TPM_missing_samples,
        collapse = ", "
      )
    )
  )
}

GSE104736_HeLa_TranscriptLevel_TPM <-
  GSE104736_HeLa_TranscriptLevel_TPM %>%
  dplyr::mutate(
    GREIN_transcript_id =
      as.character(
        ensembl_transcript_id
      ),

    ensembl_transcript_id = sub(
      "\\..*$",
      "",
      GREIN_transcript_id
    ),

    dplyr::across(
      dplyr::all_of(
        GSE104736_HeLa_TPM_required_samples
      ),
      as.numeric
    )
  ) %>%
  dplyr::select(
    GREIN_transcript_id,
    ensembl_transcript_id,
    dplyr::all_of(
      GSE104736_HeLa_TPM_required_samples
    )
  )

if (
  anyDuplicated(
    GSE104736_HeLa_TranscriptLevel_TPM$
    ensembl_transcript_id
  ) > 0
) {
  stop(
    paste(
      "Duplicated versionless transcript identifiers",
      "were found in the GSE104736 TPM table."
    )
  )
}

GSE104736_HeLa_TPM_values_are_valid <- vapply(
  GSE104736_HeLa_TPM_required_samples,

  function(current_sample) {

    current_TPM <- getElement(
      GSE104736_HeLa_TranscriptLevel_TPM,
      current_sample
    )

    all(
      is.finite(current_TPM) &
        current_TPM >= 0
    )
  },

  logical(1)
)

if (
  !all(
    GSE104736_HeLa_TPM_values_are_valid
  )
) {
  stop(
    paste(
      "At least one GSE104736 sample contains",
      "missing, non-finite, or negative TPM values."
    )
  )
}

if (
  !all(
    c(
      "ensembl_transcript_id",
      "Normalized_MFE"
    ) %in% colnames(Table_ANNOTATED)
  )
) {
  stop(
    paste(
      "Table_ANNOTATED lacks ensembl_transcript_id",
      "or Normalized_MFE."
    )
  )
}

Table_ANNOTATED_for_GSE104736_TPM <-
  Table_ANNOTATED %>%
  dplyr::transmute(
    ensembl_transcript_id = sub(
      "\\..*$",
      "",
      as.character(
        ensembl_transcript_id
      )
    ),

    Normalized_MFE =
      as.numeric(Normalized_MFE)
  ) %>%
  dplyr::filter(
    !is.na(
      ensembl_transcript_id
    ),

    ensembl_transcript_id != "",

    is.finite(
      Normalized_MFE
    )
  ) %>%
  dplyr::distinct(
    ensembl_transcript_id,
    .keep_all = TRUE
  )

GSE104736_HeLa_TPM_thermodynamic <-
  GSE104736_HeLa_TranscriptLevel_TPM %>%
  dplyr::inner_join(
    Table_ANNOTATED_for_GSE104736_TPM,
    by = "ensembl_transcript_id"
  )

if (
  nrow(
    GSE104736_HeLa_TPM_thermodynamic
  ) == 0
) {
  stop(
    "The GSE104736 TPM-to-annotation join produced no transcripts."
  )
}

if (
  anyDuplicated(
    GSE104736_HeLa_TPM_thermodynamic$
    ensembl_transcript_id
  ) > 0
) {
  stop(
    paste(
      "The GSE104736 annotation join generated",
      "duplicated transcript identifiers."
    )
  )
}

GSE104736_HeLa_TPM_original_totals <-
  tibble::tibble(
    Sample =
      GSE104736_HeLa_TPM_required_samples,

    Original_total_TPM = vapply(
      GSE104736_HeLa_TPM_required_samples,

      function(current_sample) {
        sum(
          getElement(
            GSE104736_HeLa_TranscriptLevel_TPM,
            current_sample
          )
        )
      },

      numeric(1)
    )
  )

GSE104736_HeLa_TPM_long <-
  GSE104736_HeLa_TPM_thermodynamic %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(
      GSE104736_HeLa_TPM_required_samples
    ),
    names_to = "Sample",
    values_to = "TPM"
  ) %>%
  dplyr::left_join(
    GSE104736_HeLa_TPM_sample_metadata,
    by = "Sample"
  ) %>%
  dplyr::mutate(
    Dataset = "GSE104736",

    Cell_line = "HeLa",

    TPM =
      as.numeric(TPM)
  ) %>%
  dplyr::group_by(
    Sample,
    Phase
  ) %>%
  dplyr::mutate(
    retained_TPM_MFE =
      sum(TPM),

    relative_transcript_abundance_MFE =
      TPM / retained_TPM_MFE,

    TPM_weighted_DGfold_contribution =
      relative_transcript_abundance_MFE *
      Normalized_MFE
  ) %>%
  dplyr::ungroup()

GSE104736_HeLa_TPM_weight_validation <-
  GSE104736_HeLa_TPM_long %>%
  dplyr::group_by(
    Dataset,
    Cell_line,
    Sample,
    Phase
  ) %>%
  dplyr::summarise(
    annotated_transcripts =
      dplyr::n_distinct(
        ensembl_transcript_id
      ),

    expressed_annotated_transcripts =
      dplyr::n_distinct(
        ensembl_transcript_id[TPM > 0]
      ),

    retained_TPM_MFE =
      dplyr::first(
        retained_TPM_MFE
      ),

    sum_MFE_weights = sum(
      relative_transcript_abundance_MFE,
      na.rm = TRUE
    ),

    TPM_weighted_Normalized_MFE = sum(
      TPM_weighted_DGfold_contribution,
      na.rm = TRUE
    ),

    .groups = "drop"
  ) %>%
  dplyr::left_join(
    GSE104736_HeLa_TPM_original_totals,
    by = "Sample"
  ) %>%
  dplyr::mutate(
    MFE_coverage_percent =
      100 *
      retained_TPM_MFE /
      Original_total_TPM,

    abundance_measure =
      "Transcript-level TPM"
  ) %>%
  dplyr::select(
    Dataset,
    Cell_line,
    Sample,
    Phase,
    annotated_transcripts,
    expressed_annotated_transcripts,
    Original_total_TPM,
    retained_TPM_MFE,
    MFE_coverage_percent,
    sum_MFE_weights,
    TPM_weighted_Normalized_MFE,
    abundance_measure
  )

if (
  any(
    abs(
      GSE104736_HeLa_TPM_weight_validation$
      sum_MFE_weights - 1
    ) > 1e-8
  )
) {
  stop(
    paste(
      "The GSE104736 TPM-derived ΔGfold weights",
      "do not sum to one."
    )
  )
}

cat(
  "\nGSE104736 sample-level TPM-weighted ΔGfold:\n"
)

print(
  as.data.frame(
    GSE104736_HeLa_TPM_weight_validation
  ),
  row.names = FALSE,
  digits = 12
)

Counts_GSE104736_GREIN_gene_level <- read.delim(
  hela_gene_raw_counts_file,
  sep = ",",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
) %>%
  tibble::column_to_rownames(var = "X") %>%
  dplyr::select(
    -c(
      gene_symbol,
      GSM2806906,
      GSM2806907,
      GSM2806908,
      GSM2806909,
      GSM2806910,
      GSM2806911
    )
  ) %>%
  as.matrix()

missing_hela_gene_samples <- setdiff(
  rownames(Expdata_GSE104736),
  colnames(Counts_GSE104736_GREIN_gene_level)
)

if (length(missing_hela_gene_samples) > 0) {
  stop(
    "Missing HeLa gene-count samples: ",
    paste(missing_hela_gene_samples, collapse = ", "),
    call. = FALSE
  )
}

Counts_GSE104736_GREIN_gene_level <- floor(
  Counts_GSE104736_GREIN_gene_level[
    ,
    rownames(Expdata_GSE104736),
    drop = FALSE
  ]
)

stopifnot(
  identical(
    colnames(Counts_GSE104736_GREIN_gene_level),
    rownames(Expdata_GSE104736)
  )
)

dds_GSE104736_gene <- DESeq2::DESeqDataSetFromMatrix(
  countData = Counts_GSE104736_GREIN_gene_level,
  colData = Expdata_GSE104736,
  design = ~condition
)

dds_GSE104736_gene <- dds_GSE104736_gene[
  rowSums(DESeq2::counts(dds_GSE104736_gene)) >= 3,
]
dds_GSE104736_gene <- DESeq2::DESeq(dds_GSE104736_gene)

res_GSE104736_G1_vs_G2_gene_df <- DESeq2::results(
  dds_GSE104736_gene,
  contrast = c("condition", "HeLa_G1", "HeLa_G2")
) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("ensembl_gene_id") %>%
  dplyr::select(ensembl_gene_id, log2FoldChange)

res_GSE104736_S_vs_G1_gene_df <- DESeq2::results(
  dds_GSE104736_gene,
  contrast = c("condition", "HeLa_S", "HeLa_G1")
) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("ensembl_gene_id") %>%
  dplyr::select(ensembl_gene_id, log2FoldChange)

res_GSE104736_G2_vs_S_gene_df <- DESeq2::results(
  dds_GSE104736_gene,
  contrast = c("condition", "HeLa_G2", "HeLa_S")
) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("ensembl_gene_id") %>%
  dplyr::select(ensembl_gene_id, log2FoldChange)

GSE104736_gene_entrez_id_G1vsG2_entrez_IDs <-
  res_GSE104736_G1_vs_G2_gene_df %>%
  dplyr::left_join(
    map_ensembl_gene_to_entrez(
      res_GSE104736_G1_vs_G2_gene_df$ensembl_gene_id,
      mart
    ),
    by = "ensembl_gene_id"
  ) %>%
  dplyr::rename(G1vsG2_log2FoldChange = log2FoldChange)

GSE104736_gene_entrez_id_SvsG1_entrez_IDs <-
  res_GSE104736_S_vs_G1_gene_df %>%
  dplyr::left_join(
    map_ensembl_gene_to_entrez(
      res_GSE104736_S_vs_G1_gene_df$ensembl_gene_id,
      mart
    ),
    by = "ensembl_gene_id"
  ) %>%
  dplyr::rename(SvsG1_log2FoldChange = log2FoldChange)

GSE104736_gene_entrez_id_G2vsS_entrez_IDs <-
  res_GSE104736_G2_vs_S_gene_df %>%
  dplyr::left_join(
    map_ensembl_gene_to_entrez(
      res_GSE104736_G2_vs_S_gene_df$ensembl_gene_id,
      mart
    ),
    by = "ensembl_gene_id"
  ) %>%
  dplyr::rename(G2vsS_log2FoldChange = log2FoldChange)

GSE104736_gene_level_Log2_FC_df <-
  GSE104736_gene_entrez_id_G1vsG2_entrez_IDs %>%
  dplyr::select(
    ensembl_gene_id,
    entrezgene_id,
    G1vsG2_log2FoldChange
  ) %>%
  dplyr::left_join(
    GSE104736_gene_entrez_id_SvsG1_entrez_IDs %>%
      dplyr::select(
        ensembl_gene_id,
        entrezgene_id,
        SvsG1_log2FoldChange
      ),
    by = c("ensembl_gene_id", "entrezgene_id")
  ) %>%
  dplyr::left_join(
    GSE104736_gene_entrez_id_G2vsS_entrez_IDs %>%
      dplyr::select(
        ensembl_gene_id,
        entrezgene_id,
        G2vsS_log2FoldChange
      ),
    by = c("ensembl_gene_id", "entrezgene_id")
  )

combined_kk2_df <- run_kegg_gsea(
  GSE104736_gene_level_Log2_FC_df,
  c(
    "G1vsG2_log2FoldChange",
    "SvsG1_log2FoldChange",
    "G2vsS_log2FoldChange"
  ),
  c("GSE104736_G1vsG2", "GSE104736_SvsG1", "GSE104736_G2vsS")
)

Selected_KEGG_pathways <- c(
  "hsa04150", "hsa04668", "hsa04657", "hsa03030", "hsa03410",
  "hsa03430", "hsa03420", "hsa04370", "hsa4072", "hsa03320",
  "hsa03460", "hsa03440", "hsa04014", "hsa04070", "hsa04020",
  "hsa01230", "hsa04310", "hsa03082", "hsa04015", "hsa04625",
  "hsa04371", "hsa04390", "hsa04010", "hsa04622", "hsa04350",
  "hsa04152", "hsa04926", "hsa04630", "hsa04115", "hsa04151",
  "hsa04910", "hsa04330", "hsa04012", "hsa03440", "hsa03460"
)

combined_kk2_df_filtered <- combined_kk2_df %>%
  dplyr::filter(ID %in% Selected_KEGG_pathways) %>%
  dplyr::mutate(NES = scales::rescale(NES, to = c(-4, 4)))

desired_order <- c(
  "Biosynthesis of amino acids",
  "DNA replication",
  "Homologous recombination",
  "Fanconi anemia pathway",
  "Base excision repair"
)

GSE104736_combined_kk2_df_filtered_plot <- combined_kk2_df_filtered %>%
  dplyr::mutate(
    Description = factor(Description, levels = desired_order),
    neg_log10_p_adjust = -log10(p.adjust),
    Sample = factor(
      Sample,
      levels = c("GSE104736_SvsG1", "GSE104736_G2vsS")
    )
  )

prepare_Fig2D_normalized_phase_counts <- function(
    raw_count_df,
    sample_metadata,
    dataset_label
) {

  required_metadata_columns <- c(
    "Sample",
    "Phase",
    "Replicate"
  )

  stopifnot(
    all(
      required_metadata_columns %in%
        colnames(sample_metadata)
    )
  )

  sample_metadata <- sample_metadata %>%
    dplyr::mutate(
      Sample = as.character(Sample),
      Phase = as.character(Phase),
      Replicate = as.character(Replicate)
    )

  missing_samples <- setdiff(
    sample_metadata$Sample,
    colnames(raw_count_df)
  )

  if (length(missing_samples) > 0) {
    stop(
      paste0(
        dataset_label,
        " is missing the following raw-count samples: ",
        paste(
          missing_samples,
          collapse = ", "
        )
      )
    )
  }

  raw_count_df <- raw_count_df %>%
    dplyr::mutate(
      ensembl_transcript_id = sub(
        "\\..*$",
        "",
        as.character(ensembl_transcript_id)
      )
    ) %>%
    dplyr::select(
      ensembl_transcript_id,
      dplyr::all_of(
        sample_metadata$Sample
      )
    ) %>%
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(
          sample_metadata$Sample
        ),
        as.numeric
      )
    ) %>%
    dplyr::group_by(
      ensembl_transcript_id
    ) %>%
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(
          sample_metadata$Sample
        ),
        ~ sum(.x, na.rm = TRUE)
      ),
      .groups = "drop"
    )

  count_matrix <- raw_count_df %>%
    tibble::column_to_rownames(
      "ensembl_transcript_id"
    ) %>%
    as.matrix()

  if (
    any(!is.finite(count_matrix)) ||
    any(count_matrix < 0)
  ) {
    stop(
      paste0(
        dataset_label,
        " contains invalid raw-count values."
      )
    )
  }

  count_matrix <- floor(count_matrix)
  storage.mode(count_matrix) <- "integer"

  sample_metadata <- sample_metadata[
    match(
      colnames(count_matrix),
      sample_metadata$Sample
    ),
    ,
    drop = FALSE
  ]

  sample_information <- data.frame(
    Phase = factor(
      sample_metadata$Phase
    ),
    row.names = sample_metadata$Sample
  )

  stopifnot(
    identical(
      colnames(count_matrix),
      rownames(sample_information)
    )
  )

  dds_Fig2D <- DESeq2::DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = sample_information,
    design = ~ 1
  )

  dds_Fig2D <- DESeq2::estimateSizeFactors(
    dds_Fig2D
  )

  normalized_count_df <- DESeq2::counts(
    dds_Fig2D,
    normalized = TRUE
  ) %>%
    as.data.frame(
      check.names = FALSE
    ) %>%
    tibble::rownames_to_column(
      "ensembl_transcript_id"
    )

  normalized_count_long <- normalized_count_df %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(
        sample_metadata$Sample
      ),
      names_to = "Sample",
      values_to = "normalized_counts"
    ) %>%
    dplyr::left_join(
      sample_metadata,
      by = "Sample"
    ) %>%
    dplyr::mutate(
      Dataset = dataset_label,
      normalized_counts =
        as.numeric(normalized_counts)
    )

  normalized_phase_counts <- normalized_count_long %>%
    dplyr::group_by(
      Dataset,
      ensembl_transcript_id,
      Phase
    ) %>%
    dplyr::summarise(
      normalized_counts =
        mean(
          normalized_counts,
          na.rm = TRUE
        ),

      n_biological_samples =
        dplyr::n_distinct(Sample),

      .groups = "drop"
    ) %>%
    dplyr::mutate(
      Sample =
        paste0(
          Dataset,
          "_",
          Phase,
          "_phase_mean"
        ),

      Replicate = dplyr::if_else(
        n_biological_samples == 1,
        "single_sample",
        paste0(
          "mean_of_",
          n_biological_samples
        )
      )
    ) %>%
    dplyr::select(
      Dataset,
      ensembl_transcript_id,
      Sample,
      Phase,
      Replicate,
      n_biological_samples,
      normalized_counts
    )

  return(
    list(
      sample_level = normalized_count_long,
      phase_level = normalized_phase_counts
    )
  )
}

GSE94479_Fig2D_raw_count_file <- mcf7_transcript_raw_counts_file

GSE94479_Fig2D_raw_counts <- read.delim(
  GSE94479_Fig2D_raw_count_file,
  sep = ",",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

colnames(
  GSE94479_Fig2D_raw_counts
)[1] <- "ensembl_transcript_id"

GSE94479_Fig2D_sample_metadata <-
  tibble::tribble(
    ~Sample,                    ~Phase,  ~Replicate,
    "SRR5227986_GSM2476255",    "G1",    "1",
    "SRR5227987_GSM2476256",    "S",     "1",
    "SRR5227988_GSM2476257",    "G2/M",  "1"
  )

GSE94479_Fig2D_normalized_counts <-
  prepare_Fig2D_normalized_phase_counts(
    raw_count_df =
      GSE94479_Fig2D_raw_counts,

    sample_metadata =
      GSE94479_Fig2D_sample_metadata,

    dataset_label =
      "MCF7"
  )

if (!exists("Counts_GSE104736_GREIN")) {
  stop(
    paste(
      "Counts_GSE104736_GREIN is missing.",
      "Run the GSE104736 differential-expression section first."
    )
  )
}

GSE104736_Fig2D_raw_counts <-
  as.data.frame(
    Counts_GSE104736_GREIN,
    check.names = FALSE
  ) %>%
  tibble::rownames_to_column(
    "ensembl_transcript_id"
  )

GSE104736_Fig2D_sample_metadata <-
  tibble::tribble(
    ~Sample,       ~Phase,  ~Replicate,
    "GSM2806897",  "G1",    "1",
    "GSM2806898",  "G1",    "2",
    "GSM2806899",  "G1",    "3",
    "GSM2806900",  "S",     "1",
    "GSM2806901",  "S",     "2",
    "GSM2806902",  "S",     "3",
    "GSM2806903",  "G2/M",  "1",
    "GSM2806904",  "G2/M",  "2",
    "GSM2806905",  "G2/M",  "3"
  )

GSE104736_Fig2D_normalized_counts <-
  prepare_Fig2D_normalized_phase_counts(
    raw_count_df =
      GSE104736_Fig2D_raw_counts,

    sample_metadata =
      GSE104736_Fig2D_sample_metadata,

    dataset_label =
      "HeLa"
  )

calculate_Fig2D_quartile_fraction <- function(
    normalized_phase_counts,
    de_df,
    dataset_label,
    transition_label,
    phase_1,
    phase_2
) {

  required_DE_columns <- c(
    "ensembl_transcript_id",
    "padj",
    "Normalized_MFE"
  )

  stopifnot(
    all(
      required_DE_columns %in%
        colnames(de_df)
    )
  )

  de_universe <- de_df %>%
    dplyr::transmute(
      ensembl_transcript_id = sub(
        "\\..*$",
        "",
        as.character(ensembl_transcript_id)
      ),

      padj = as.numeric(padj),

      Normalized_MFE =
        as.numeric(Normalized_MFE)
    ) %>%
    dplyr::filter(
      !is.na(padj),
      padj < 0.05,
      is.finite(Normalized_MFE)
    ) %>%
    dplyr::distinct(
      ensembl_transcript_id,
      .keep_all = TRUE
    )

  if (nrow(de_universe) == 0) {
    stop(
      paste(
        "No significant-DE transcripts were retained for",
        dataset_label,
        transition_label
      )
    )
  }

  q25 <- as.numeric(
    stats::quantile(
      de_universe$Normalized_MFE,
      probs = 0.25,
      na.rm = TRUE
    )
  )

  q75 <- as.numeric(
    stats::quantile(
      de_universe$Normalized_MFE,
      probs = 0.75,
      na.rm = TRUE
    )
  )

  de_universe <- de_universe %>%
    dplyr::mutate(
      DGfold_quartile =
        dplyr::case_when(
          Normalized_MFE >= q75 ~
            "Least stable\nupper DGfold quartile",

          Normalized_MFE <= q25 ~
            "Most stable\nlower DGfold quartile",

          TRUE ~
            "Middle 50%"
        )
    )

  de_count_data <- normalized_phase_counts %>%
    dplyr::filter(
      Phase %in% c(
        phase_1,
        phase_2
      )
    ) %>%
    dplyr::inner_join(
      de_universe,
      by = "ensembl_transcript_id"
    ) %>%
    dplyr::mutate(
      Dataset = dataset_label,
      Transition = transition_label,
      normalized_counts =
        as.numeric(normalized_counts)
    ) %>%
    dplyr::filter(
      is.finite(normalized_counts),
      normalized_counts >= 0
    )

  phase_totals <- de_count_data %>%
    dplyr::group_by(
      Dataset,
      Transition,
      Phase
    ) %>%
    dplyr::summarise(
      total_normalized_counts_DE_universe =
        sum(
          normalized_counts,
          na.rm = TRUE
        ),

      n_significant_DE_transcripts =
        dplyr::n_distinct(
          ensembl_transcript_id
        ),

      n_significant_DE_transcripts_with_counts =
        dplyr::n_distinct(
          ensembl_transcript_id[
            normalized_counts > 0
          ]
        ),

      .groups = "drop"
    )

  quartile_fraction <- de_count_data %>%
    dplyr::filter(
      DGfold_quartile != "Middle 50%"
    ) %>%
    dplyr::group_by(
      Dataset,
      Transition,
      Phase,
      DGfold_quartile
    ) %>%
    dplyr::summarise(
      normalized_counts_in_quartile =
        sum(
          normalized_counts,
          na.rm = TRUE
        ),

      n_transcript_species =
        dplyr::n_distinct(
          ensembl_transcript_id
        ),

      n_transcript_species_with_counts =
        dplyr::n_distinct(
          ensembl_transcript_id[
            normalized_counts > 0
          ]
        ),

      .groups = "drop"
    ) %>%
    dplyr::left_join(
      phase_totals,
      by = c(
        "Dataset",
        "Transition",
        "Phase"
      )
    ) %>%
    dplyr::mutate(
      normalized_count_fraction =
        normalized_counts_in_quartile /
        total_normalized_counts_DE_universe,

      q25 = q25,
      q75 = q75,
      phase_1 = phase_1,
      phase_2 = phase_2,

      abundance_measure =
        "DESeq2-normalized transcript-level read counts"
    )

  return(quartile_fraction)
}

MCF7_G1_vs_G2M_DGfold_quartile_fraction <-
  calculate_Fig2D_quartile_fraction(
    normalized_phase_counts =
      GSE94479_Fig2D_normalized_counts$
      phase_level,

    de_df =
      GSE94479_G0G1vsG2M,

    dataset_label =
      "MCF7",

    transition_label =
      "G1 vs G2/M",

    phase_1 =
      "G2/M",

    phase_2 =
      "G1"
  )

MCF7_S_vs_G1_DGfold_quartile_fraction <-
  calculate_Fig2D_quartile_fraction(
    normalized_phase_counts =
      GSE94479_Fig2D_normalized_counts$
      phase_level,

    de_df =
      GSE94479_SvsG0G1,

    dataset_label =
      "MCF7",

    transition_label =
      "S vs G1",

    phase_1 =
      "G1",

    phase_2 =
      "S"
  )

HeLa_G1_vs_G2M_DGfold_quartile_fraction <-
  calculate_Fig2D_quartile_fraction(
    normalized_phase_counts =
      GSE104736_Fig2D_normalized_counts$
      phase_level,

    de_df =
      MFE_GSE104736_G1_vs_G2,

    dataset_label =
      "HeLa",

    transition_label =
      "G1 vs G2/M",

    phase_1 =
      "G2/M",

    phase_2 =
      "G1"
  )

HeLa_S_vs_G1_DGfold_quartile_fraction <-
  calculate_Fig2D_quartile_fraction(
    normalized_phase_counts =
      GSE104736_Fig2D_normalized_counts$
      phase_level,

    de_df =
      MFE_GSE104736_S_vs_G1,

    dataset_label =
      "HeLa",

    transition_label =
      "S vs G1",

    phase_1 =
      "G1",

    phase_2 =
      "S"
  )

Cell_cycle_MCF7_HeLa_DGfold_quartile_DE_fraction <-
  dplyr::bind_rows(
    MCF7_G1_vs_G2M_DGfold_quartile_fraction,
    MCF7_S_vs_G1_DGfold_quartile_fraction,
    HeLa_G1_vs_G2M_DGfold_quartile_fraction,
    HeLa_S_vs_G1_DGfold_quartile_fraction
  ) %>%
  dplyr::mutate(
    Dataset = factor(
      Dataset,
      levels = c(
        "MCF7",
        "HeLa"
      )
    ),

    Transition = factor(
      Transition,
      levels = c(
        "G1 vs G2/M",
        "S vs G1"
      )
    ),

    Phase = factor(
      Phase,
      levels = c(
        "G2/M",
        "G1",
        "S"
      )
    ),

    DGfold_quartile = factor(
      DGfold_quartile,
      levels = c(
        "Least stable\nupper DGfold quartile",
        "Most stable\nlower DGfold quartile"
      )
    )
  ) %>%
  dplyr::arrange(
    Dataset,
    Transition,
    DGfold_quartile,
    Phase
  )

Fig2D_required_columns <- c(
  "Dataset",
  "Transition",
  "Phase",
  "DGfold_quartile",
  "normalized_count_fraction",
  "normalized_counts_in_quartile",
  "total_normalized_counts_DE_universe",
  "phase_1",
  "phase_2",
  "abundance_measure"
)

stopifnot(
  all(
    Fig2D_required_columns %in%
      colnames(
        Cell_cycle_MCF7_HeLa_DGfold_quartile_DE_fraction
      )
  ),

  nrow(
    Cell_cycle_MCF7_HeLa_DGfold_quartile_DE_fraction
  ) == 16,

  all(
    is.finite(
      Cell_cycle_MCF7_HeLa_DGfold_quartile_DE_fraction$
        normalized_count_fraction
    )
  ),

  all(
    dplyr::between(
      Cell_cycle_MCF7_HeLa_DGfold_quartile_DE_fraction$
        normalized_count_fraction,
      0,
      1
    )
  )
)

cat(
  "\nFig. 2D DESeq2-normalized read-count fractions:\n"
)

print(
  as.data.frame(
    Cell_cycle_MCF7_HeLa_DGfold_quartile_DE_fraction
  ),
  row.names = FALSE,
  digits = 10
)

figure2_results <- list(
  GSE94479_combined_kk2_df_filtered_plot =
    GSE94479_combined_kk2_df_filtered_plot,
  GSE104736_combined_kk2_df_filtered_plot =
    GSE104736_combined_kk2_df_filtered_plot,
  GSE94479_G0G1vsG2M = GSE94479_G0G1vsG2M,
  GSE94479_SvsG0G1 = GSE94479_SvsG0G1,
  MFE_GSE104736_G1_vs_G2 = MFE_GSE104736_G1_vs_G2,
  MFE_GSE104736_S_vs_G1 = MFE_GSE104736_S_vs_G1,
  GSE94479_Fig2D_normalized_counts =
    GSE94479_Fig2D_normalized_counts,
  GSE104736_Fig2D_normalized_counts =
    GSE104736_Fig2D_normalized_counts,
  Cell_cycle_MCF7_HeLa_DGfold_quartile_DE_fraction =
    Cell_cycle_MCF7_HeLa_DGfold_quartile_DE_fraction,
  GSE94479_MCF7_TPM_long = GSE94479_MCF7_TPM_long,
  GSE94479_MCF7_TPM_weight_validation =
    GSE94479_MCF7_TPM_weight_validation,
  GSE104736_HeLa_TPM_long = GSE104736_HeLa_TPM_long,
  GSE104736_HeLa_TPM_weight_validation =
    GSE104736_HeLa_TPM_weight_validation
)

Fig2_GSEA_objects <- c(
  "GSE94479_combined_kk2_df_filtered_plot",
  "GSE104736_combined_kk2_df_filtered_plot"
)

Fig2_GSEA_required_columns <- c(
  "ID",
  "Description",
  "Sample",
  "NES",
  "p.adjust"
)

for (current_object in Fig2_GSEA_objects) {
  current_data <- figure2_results[[current_object]]
  current_missing_columns <- setdiff(
    Fig2_GSEA_required_columns,
    colnames(current_data)
  )

  if (length(current_missing_columns) > 0) {
    stop(
      current_object,
      " lacks the following columns: ",
      paste(current_missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
}

Fig2_DE_objects <- c(
  "GSE94479_G0G1vsG2M",
  "GSE94479_SvsG0G1",
  "MFE_GSE104736_G1_vs_G2",
  "MFE_GSE104736_S_vs_G1"
)

Fig2_DE_required_columns <- c(
  "ensembl_transcript_id",
  "log2FoldChange",
  "padj",
  "Normalized_MFE"
)

for (current_object in Fig2_DE_objects) {
  current_data <- figure2_results[[current_object]]

  current_missing_columns <- setdiff(
    Fig2_DE_required_columns,
    colnames(current_data)
  )

  if (length(current_missing_columns) > 0) {
    stop(
      paste0(
        current_object,
        " lacks the following columns: ",
        paste(
          current_missing_columns,
          collapse = ", "
        )
      )
    )
  }
}

Fig2_TPM_long_objects <- c(
  "GSE94479_MCF7_TPM_long",
  "GSE104736_HeLa_TPM_long"
)

Fig2_TPM_long_required_columns <- c(
  "ensembl_transcript_id",
  "Sample",
  "Phase",
  "TPM",
  "Normalized_MFE",
  "relative_transcript_abundance_MFE",
  "TPM_weighted_DGfold_contribution"
)

for (current_object in Fig2_TPM_long_objects) {
  current_data <- figure2_results[[current_object]]

  current_missing_columns <- setdiff(
    Fig2_TPM_long_required_columns,
    colnames(current_data)
  )

  if (length(current_missing_columns) > 0) {
    stop(
      paste0(
        current_object,
        " lacks the following columns: ",
        paste(
          current_missing_columns,
          collapse = ", "
        )
      )
    )
  }
}

Fig2_TPM_summary_required_columns <- c(
  "Sample",
  "Phase",
  "annotated_transcripts",
  "expressed_annotated_transcripts",
  "retained_TPM_MFE",
  "sum_MFE_weights",
  "TPM_weighted_Normalized_MFE"
)

for (
  current_object in c(
    "GSE94479_MCF7_TPM_weight_validation",
    "GSE104736_HeLa_TPM_weight_validation"
  )
) {
  current_data <- figure2_results[[current_object]]

  current_missing_columns <- setdiff(
    Fig2_TPM_summary_required_columns,
    colnames(current_data)
  )

  if (length(current_missing_columns) > 0) {
    stop(
      paste0(
        current_object,
        " lacks the following columns: ",
        paste(
          current_missing_columns,
          collapse = ", "
        )
      )
    )
  }
}

stopifnot(
  nrow(
    GSE94479_MCF7_TPM_weight_validation
  ) == 3,

  nrow(
    GSE104736_HeLa_TPM_weight_validation
  ) == 9,

  setequal(
    as.character(
      GSE94479_MCF7_TPM_weight_validation$Phase
    ),
    c("G1", "S", "G2/M")
  ),

  setequal(
    as.character(
      GSE104736_HeLa_TPM_weight_validation$Phase
    ),
    c("G1", "S", "G2")
  ),

  all(
    abs(
      GSE94479_MCF7_TPM_weight_validation$
        sum_MFE_weights - 1
    ) < 1e-8
  ),

  all(
    abs(
      GSE104736_HeLa_TPM_weight_validation$
        sum_MFE_weights - 1
    ) < 1e-8
  )
)

saveRDS(figure2_results, figure2_results_file)

cat("\nCell-cycle analysis complete.\n")
cat("Figure 2 input objects written to:\n", figure2_results_file, "\n")
cat("Other outputs written to:\n", analysis_output_dir, "\n")
