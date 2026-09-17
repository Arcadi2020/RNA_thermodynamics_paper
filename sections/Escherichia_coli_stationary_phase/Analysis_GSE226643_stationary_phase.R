options(timeout = max(1200, getOption("timeout")))

suppressPackageStartupMessages({
  library(GEOquery)
  library(DESeq2)
  library(SummarizedExperiment)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
  library(clusterProfiler)
  library(httr)
  library(jsonlite)
})

GSE226643_project_dir <- "path_to/GSE226643"
GSE226643_prediction_file <- "../../data/Prediction_Ecoli_MFE_HC_df.csv"
GSE226643_dir <- GSE226643_project_dir

GSE226643_metadata_dir <- file.path(GSE226643_dir, "metadata")
GSE226643_geo_dir <- file.path(GSE226643_dir, "GEO_download")
GSE226643_github_dir <- file.path(GSE226643_dir, "author_GitHub_repository")
GSE226643_sra_dir <- file.path(GSE226643_dir, "SRA_selected_samples")
GSE226643_fastq_dir <- file.path(GSE226643_sra_dir, "FASTQ")

for (one_dir in c(
  GSE226643_dir,
  GSE226643_metadata_dir,
  GSE226643_geo_dir,
  GSE226643_github_dir,
  GSE226643_sra_dir,
  GSE226643_fastq_dir
)) {
  dir.create(one_dir, recursive = TRUE, showWarnings = FALSE)
}

download_selected_fastq <- FALSE

as_one_string <- function(x) {
  if (is.null(x) || length(x) == 0L) {
    return(NA_character_)
  }
  paste(as.character(unlist(x, use.names = FALSE)), collapse = " | ")
}

write_tsv <- function(x, file) {
  utils::write.table(
    x,
    file = file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = ""
  )
}

download_if_needed <- function(url, destination, expected_size = NA_real_) {
  have_valid_file <- file.exists(destination) && file.info(destination)$size > 0L
  if (have_valid_file && is.finite(expected_size)) {
    have_valid_file <- file.info(destination)$size == expected_size
  }

  if (!have_valid_file) {
    message("Downloading: ", basename(destination))
    utils::download.file(url, destination, mode = "wb", quiet = FALSE)
  }

  if (!file.exists(destination) || file.info(destination)$size == 0L) {
    stop("Download failed or produced an empty file: ", destination, call. = FALSE)
  }
}

extract_SRX <- function(x) {
  hit <- regmatches(x, regexpr("SRX[0-9]+", x, perl = TRUE))
  ifelse(length(hit) == 0L || identical(hit, ""), NA_character_, hit)
}

classify_data_file <- function(file_name) {
  if (grepl("raw.*count|count.*matrix|featurecount", file_name, ignore.case = TRUE)) {
    return("candidate_raw_count_matrix")
  }
  if (grepl("tpm|fpkm|rpkm|normaliz|log", file_name, ignore.case = TRUE)) {
    return("normalised_expression_or_log_expression")
  }
  if (grepl("metadata|sample|design|annotation|gene_info", file_name, ignore.case = TRUE)) {
    return("metadata_or_annotation")
  }
  if (grepl("\\.R$|script|code", file_name, ignore.case = TRUE)) {
    return("code_or_script")
  }
  "other"
}

GSE226643_selected_manifest <- data.frame(
  GSM = c(
    "GSM7080903", "GSM7080904", "GSM7080905", "GSM7080906",
    "GSM7080907", "GSM7080908", "GSM7080909", "GSM7080910",
    "GSM7080911", "GSM7080912", "GSM7080913", "GSM7080914",
    "GSM7080915", "GSM7080916", "GSM7080917", "GSM7080918",
    "GSM7080919", "GSM7080920", "GSM7080921", "GSM7080922",
    "GSM7080923", "GSM7080924", "GSM7080925", "GSM7080926"
  ),
  nutrient_limitation = rep(c("Nitrogen", "Sulfur"), each = 12L),
  time_h = rep(rep(c(3L, 4L, 6L, 8L, 12L, 24L), each = 2L), 2L),
  biological_replicate = rep(c(1L, 2L), 12L),
  stringsAsFactors = FALSE
)

GSE226643_selected_manifest$condition <- paste0(
  GSE226643_selected_manifest$nutrient_limitation,
  "_limitation_",
  GSE226643_selected_manifest$time_h,
  "h"
)
GSE226643_selected_manifest$analysis_role <- "time_course_only"
GSE226643_selected_manifest$analysis_role[
  GSE226643_selected_manifest$time_h %in% c(3L, 12L)
] <- "primary_12h_vs_3h"
GSE226643_selected_manifest$contrast_id <- NA_character_
GSE226643_selected_manifest$contrast_id[
  GSE226643_selected_manifest$nutrient_limitation == "Nitrogen" &
    GSE226643_selected_manifest$time_h %in% c(3L, 12L)
] <- "N_limitation_12h_vs_3h"
GSE226643_selected_manifest$contrast_id[
  GSE226643_selected_manifest$nutrient_limitation == "Sulfur" &
    GSE226643_selected_manifest$time_h %in% c(3L, 12L)
] <- "S_limitation_12h_vs_3h"

stopifnot(
  nrow(GSE226643_selected_manifest) == 24L,
  sum(GSE226643_selected_manifest$contrast_id == "N_limitation_12h_vs_3h", na.rm = TRUE) == 4L,
  sum(GSE226643_selected_manifest$contrast_id == "S_limitation_12h_vs_3h", na.rm = TRUE) == 4L
)

write_tsv(
  GSE226643_selected_manifest,
  file.path(GSE226643_metadata_dir, "GSE226643_N_S_limitation_selected_manifest.tsv")
)

GSE226643_accession <- "GSE226643"

GSE226643_GEO <- GEOquery::getGEO(
  GSE226643_accession,
  GSEMatrix = FALSE,
  getGPL = FALSE
)

GSE226643_GSM_list <- GEOquery::GSMList(GSE226643_GEO)

GSE226643_GEO_metadata <- do.call(
  rbind,
  lapply(names(GSE226643_GSM_list), function(one_GSM) {
    one_meta <- GEOquery::Meta(GSE226643_GSM_list[[one_GSM]])
    data.frame(
      GSM = one_GSM,
      title = as_one_string(one_meta$title),
      source_name_ch1 = as_one_string(one_meta$source_name_ch1),
      characteristics_ch1 = as_one_string(one_meta$characteristics_ch1),
      relation = as_one_string(one_meta$relation),
      stringsAsFactors = FALSE
    )
  })
)

write_tsv(
  GSE226643_GEO_metadata,
  file.path(GSE226643_metadata_dir, "GSE226643_GEO_sample_metadata_all.tsv")
)

GSE226643_selected_metadata <- merge(
  GSE226643_selected_manifest,
  GSE226643_GEO_metadata,
  by = "GSM",
  all.x = TRUE,
  sort = FALSE
)
GSE226643_selected_metadata <- GSE226643_selected_metadata[
  match(GSE226643_selected_manifest$GSM, GSE226643_selected_metadata$GSM),
  , drop = FALSE
]

if (any(is.na(GSE226643_selected_metadata$title))) {
  stop("At least one selected GSM accession was not found in GEO metadata.", call. = FALSE)
}

write_tsv(
  GSE226643_selected_metadata,
  file.path(GSE226643_metadata_dir, "GSE226643_N_S_limitation_selected_metadata.tsv")
)

GSE226643_supplementary_status <- tryCatch(
  {
    GEOquery::getGEOSuppFiles(
      GSE226643_accession,
      baseDir = GSE226643_geo_dir,
      makeDirectory = FALSE
    )
    "GEO supplementary-file download completed."
  },
  error = function(e) paste("GEO supplementary-file download failed:", conditionMessage(e))
)
writeLines(
  GSE226643_supplementary_status,
  con = file.path(GSE226643_metadata_dir, "GSE226643_GEO_supplementary_download_status.txt")
)

GSE226643_GitHub_archive <- file.path(
  GSE226643_github_dir,
  "Ecoli_transition_stationary_GitHub_default_branch.zip"
)
GSE226643_GitHub_url <- paste0(
  "https://api.github.com/repos/hyungyulim/",
  "Ecoli_transition_stationary/zipball"
)

GSE226643_GitHub_status <- tryCatch(
  {
    download_if_needed(
      url = GSE226643_GitHub_url,
      destination = GSE226643_GitHub_archive
    )

    GSE226643_GitHub_unpacked_dir <- file.path(GSE226643_github_dir, "unpacked")
    if (length(list.files(GSE226643_GitHub_unpacked_dir, all.files = TRUE, no.. = TRUE)) == 0L) {
      dir.create(GSE226643_GitHub_unpacked_dir, recursive = TRUE, showWarnings = FALSE)
      utils::unzip(GSE226643_GitHub_archive, exdir = GSE226643_GitHub_unpacked_dir)
    }
    "GitHub repository archive download completed."
  },
  error = function(e) paste("GitHub repository download failed:", conditionMessage(e))
)
writeLines(
  GSE226643_GitHub_status,
  con = file.path(GSE226643_metadata_dir, "GSE226643_GitHub_download_status.txt")
)

GSE226643_downloaded_files <- list.files(
  c(GSE226643_geo_dir, GSE226643_github_dir),
  recursive = TRUE,
  full.names = TRUE
)
GSE226643_file_inventory <- data.frame(
  local_file = GSE226643_downloaded_files,
  file_name = basename(GSE226643_downloaded_files),
  bytes = file.info(GSE226643_downloaded_files)$size,
  stringsAsFactors = FALSE
)
GSE226643_file_inventory$likely_content <- vapply(
  GSE226643_file_inventory$file_name,
  classify_data_file,
  character(1)
)

write_tsv(
  GSE226643_file_inventory,
  file.path(GSE226643_metadata_dir, "GSE226643_downloaded_file_inventory.tsv")
)

GSE226643_selected_metadata$SRX <- vapply(
  GSE226643_selected_metadata$relation,
  extract_SRX,
  character(1)
)

fetch_ena_runs <- function(SRX, GSM) {
  if (is.na(SRX) || SRX == "") {
    return(data.frame(
      GSM = GSM,
      run_accession = NA_character_,
      experiment_accession = NA_character_,
      fastq_ftp = NA_character_,
      fastq_md5 = NA_character_,
      library_layout = NA_character_,
      stringsAsFactors = FALSE
    ))
  }

  ena_url <- paste0(
    "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=", SRX,
    "&result=read_run&fields=run_accession,experiment_accession,fastq_ftp,",
    "fastq_md5,library_layout&format=tsv&download=false"
  )
  ena_file <- tempfile(fileext = ".tsv")
  on.exit(unlink(ena_file), add = TRUE)

  tryCatch(
    utils::download.file(ena_url, ena_file, mode = "wb", quiet = TRUE),
    error = function(e) NULL
  )

  if (!file.exists(ena_file) || file.info(ena_file)$size == 0L) {
    return(data.frame(
      GSM = GSM,
      run_accession = NA_character_,
      experiment_accession = SRX,
      fastq_ftp = NA_character_,
      fastq_md5 = NA_character_,
      library_layout = NA_character_,
      stringsAsFactors = FALSE
    ))
  }

  ena_tbl <- tryCatch(
    utils::read.delim(ena_file, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )

  if (is.null(ena_tbl) || nrow(ena_tbl) == 0L) {
    return(data.frame(
      GSM = GSM,
      run_accession = NA_character_,
      experiment_accession = SRX,
      fastq_ftp = NA_character_,
      fastq_md5 = NA_character_,
      library_layout = NA_character_,
      stringsAsFactors = FALSE
    ))
  }

  ena_tbl$GSM <- GSM
  ena_tbl[, c(
    "GSM", "run_accession", "experiment_accession", "fastq_ftp",
    "fastq_md5", "library_layout"
  ), drop = FALSE]
}

GSE226643_SRA_manifest <- do.call(
  rbind,
  Map(
    fetch_ena_runs,
    SRX = GSE226643_selected_metadata$SRX,
    GSM = GSE226643_selected_metadata$GSM
  )
)

write_tsv(
  GSE226643_SRA_manifest,
  file.path(GSE226643_sra_dir, "GSE226643_selected_ENA_SRA_run_manifest.tsv")
)

if (all(is.na(GSE226643_SRA_manifest$run_accession))) {
  warning(
    "No ENA runs were resolved automatically. Use the GEO sample metadata ",
    "and SRA Run Selector if raw FASTQ reprocessing is required.",
    call. = FALSE
  )
}

if (download_selected_fastq) {
  all_fastq_paths <- unlist(
    strsplit(
      GSE226643_SRA_manifest$fastq_ftp[!is.na(GSE226643_SRA_manifest$fastq_ftp)],
      split = ";",
      fixed = TRUE
    ),
    use.names = FALSE
  )
  all_fastq_paths <- unique(all_fastq_paths[nzchar(all_fastq_paths)])

  if (length(all_fastq_paths) == 0L) {
    stop("No FASTQ URLs are available in the ENA run manifest.", call. = FALSE)
  }

  for (one_fastq_path in all_fastq_paths) {
    download_if_needed(
      url = paste0("https://", one_fastq_path),
      destination = file.path(GSE226643_fastq_dir, basename(one_fastq_path))
    )
  }
}

writeLines(
  capture.output(sessionInfo()),
  con = file.path(GSE226643_metadata_dir, "GSE226643_download_sessionInfo.txt")
)

message(
  "GSE226643 download and manifests completed.\n",
  "Next, inspect metadata/GSE226643_downloaded_file_inventory.tsv and verify ",
  "the raw gene-by-sample count matrix before starting DE."
)

GSE226643_unpacked_dir <- file.path(
  GSE226643_project_dir,
  "author_GitHub_repository",
  "unpacked"
)

GSE226643_audit_dir <- file.path(
  GSE226643_project_dir,
  "analysis_manifests",
  "00_source_audit"
)

dir.create(
  GSE226643_audit_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!dir.exists(GSE226643_project_dir)) {
  stop(
    "GSE226643_project_dir does not exist: ",
    GSE226643_project_dir,
    call. = FALSE
  )
}

if (!dir.exists(GSE226643_unpacked_dir)) {
  stop(
    "The unpacked author-repository directory does not exist: ",
    GSE226643_unpacked_dir,
    call. = FALSE
  )
}

GSE226643_author_root_candidates <- list.files(
  GSE226643_unpacked_dir,
  full.names = TRUE,
  recursive = FALSE,
  no.. = TRUE
)

GSE226643_author_root_candidates <- GSE226643_author_root_candidates[
  file.info(GSE226643_author_root_candidates)$isdir
]

GSE226643_author_root_candidates <- GSE226643_author_root_candidates[
  file.exists(
    file.path(
      GSE226643_author_root_candidates,
      "raw_data",
      "iM",
      "counts.csv"
    )
  )
]

if (length(GSE226643_author_root_candidates) != 1L) {
  stop(
    "Exactly one unpacked repository containing raw_data/iM/counts.csv must ",
    "be present. Found ",
    length(GSE226643_author_root_candidates),
    ".",
    call. = FALSE
  )
}

GSE226643_author_root <- GSE226643_author_root_candidates[[1]]

GSE226643_counts_file <- file.path(
  GSE226643_author_root,
  "raw_data",
  "iM",
  "counts.csv"
)

GSE226643_sample_table_file <- file.path(
  GSE226643_author_root,
  "raw_data",
  "iM",
  "stationary_phase_sample_table_2.csv"
)

if (!file.exists(GSE226643_counts_file)) {
  stop("Raw-count file not found: ", GSE226643_counts_file, call. = FALSE)
}

if (!file.exists(GSE226643_sample_table_file)) {
  stop(
    "Author sample table not found: ",
    GSE226643_sample_table_file,
    call. = FALSE
  )
}

GSE226643_sample_table_raw <- readr::read_csv(
  GSE226643_sample_table_file,
  show_col_types = FALSE,
  progress = FALSE,
  name_repair = "unique"
)

if (ncol(GSE226643_sample_table_raw) < 2L) {
  stop(
    "The author sample table must contain the leading p1k-ID column plus ",
    "metadata columns.",
    call. = FALSE
  )
}

GSE226643_leading_source_column_name <- colnames(
  GSE226643_sample_table_raw
)[1]

colnames(GSE226643_sample_table_raw)[1] <- "raw_count_column"

GSE226643_required_author_columns <- c(
  "raw_count_column",
  "sample_id",
  "study",
  "project",
  "condition",
  "rep_id",
  "Strain",
  "Temperature (C)"
)

GSE226643_missing_author_columns <- setdiff(
  GSE226643_required_author_columns,
  colnames(GSE226643_sample_table_raw)
)

if (length(GSE226643_missing_author_columns) > 0L) {
  stop(
    "The author sample table lacks required column(s): ",
    paste(GSE226643_missing_author_columns, collapse = ", "),
    call. = FALSE
  )
}

GSE226643_sample_table <- GSE226643_sample_table_raw %>%
  dplyr::rename(
    author_sample_id = .data$sample_id,
    author_study = .data$study,
    author_project = .data$project,
    author_condition = .data$condition,
    author_replicate_id = .data$rep_id
  ) %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::where(is.character),
      trimws
    ),
    raw_count_column = as.character(.data$raw_count_column),
    author_sample_id = as.character(.data$author_sample_id),
    author_study = as.character(.data$author_study),
    author_project = as.character(.data$author_project),
    author_condition = as.character(.data$author_condition),
    author_replicate_id = suppressWarnings(
      as.integer(.data$author_replicate_id)
    ),
    temperature_C = suppressWarnings(
      as.numeric(.data[["Temperature (C)"]])
    )
  )

if (anyNA(GSE226643_sample_table$raw_count_column) ||
    any(GSE226643_sample_table$raw_count_column == "")) {
  stop(
    "The leading author-table column contains missing raw-count IDs.",
    call. = FALSE
  )
}

if (anyDuplicated(GSE226643_sample_table$raw_count_column) > 0L) {
  stop(
    "The leading author-table column contains duplicated raw-count IDs.",
    call. = FALSE
  )
}

if (!all(grepl("^p1k_[0-9]{5}$", GSE226643_sample_table$raw_count_column))) {
  stop(
    "The leading author-table column is not composed of expected p1k IDs. ",
    "Do not continue until the source-file structure has been checked.",
    call. = FALSE
  )
}

GSE226643_counts_header <- readr::read_csv(
  GSE226643_counts_file,
  n_max = 0L,
  show_col_types = FALSE,
  progress = FALSE,
  name_repair = "minimal"
)

if (!identical(colnames(GSE226643_counts_header)[1], "Geneid")) {
  stop(
    "The first column in counts.csv must be Geneid.",
    call. = FALSE
  )
}

GSE226643_raw_count_columns <- colnames(GSE226643_counts_header)[-1]

if (anyDuplicated(GSE226643_raw_count_columns) > 0L) {
  stop("counts.csv contains duplicated p1k column names.", call. = FALSE)
}

GSE226643_metadata_ids_missing_from_counts <- setdiff(
  GSE226643_sample_table$raw_count_column,
  GSE226643_raw_count_columns
)

GSE226643_count_columns_missing_from_metadata <- setdiff(
  GSE226643_raw_count_columns,
  GSE226643_sample_table$raw_count_column
)

if (length(GSE226643_metadata_ids_missing_from_counts) > 0L ||
    length(GSE226643_count_columns_missing_from_metadata) > 0L) {
  stop(
    "The p1k IDs in author metadata and counts.csv do not match exactly.\n",
    "Metadata IDs absent from counts.csv: ",
    paste(GSE226643_metadata_ids_missing_from_counts, collapse = ", "),
    "\ncounts.csv columns absent from metadata: ",
    paste(GSE226643_count_columns_missing_from_metadata, collapse = ", "),
    call. = FALSE
  )
}

GSE226643_complete_sample_mapping <- GSE226643_sample_table %>%
  dplyr::select(
    .data$raw_count_column,
    .data$author_sample_id,
    .data$author_study,
    .data$author_project,
    .data$author_condition,
    .data$author_replicate_id,
    .data$Strain,
    .data$temperature_C,
    dplyr::any_of(c(
      "Base Media",
      "pH",
      "Carbon Source (g/L)",
      "Nitrogen Source (g/L)",
      "Trace Element Mixture",
      "Supplement",
      "Additional Details",
      "GEO",
      "SRX",
      "Run",
      "full_name"
    ))
  ) %>%
  dplyr::arrange(
    .data$author_condition,
    .data$author_replicate_id,
    .data$raw_count_column
  )

GSE226643_condition_audit <- GSE226643_complete_sample_mapping %>%
  dplyr::group_by(
    .data$author_study,
    .data$author_project,
    .data$author_condition,
    .data$Strain,
    .data$temperature_C
  ) %>%
  dplyr::summarise(
    N_samples = dplyr::n(),
    replicate_ids = paste(
      sort(.data$author_replicate_id),
      collapse = ";"
    ),
    p1k_raw_count_columns = paste(
      .data$raw_count_column,
      collapse = ";"
    ),
    author_sample_ids = paste(
      .data$author_sample_id,
      collapse = ";"
    ),
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    .data$author_condition
  )

GSE226643_primary_candidate_audit <- GSE226643_complete_sample_mapping %>%
  dplyr::filter(
    grepl("(N|S).*GlcM9|GlcM9.*(N|S)", .data$author_sample_id),
    grepl("(3h|12h)", .data$author_sample_id)
  ) %>%
  dplyr::arrange(
    .data$author_sample_id,
    .data$author_replicate_id
  )

readr::write_csv(
  GSE226643_complete_sample_mapping,
  file.path(
    GSE226643_audit_dir,
    "GSE226643_complete_author_metadata_to_p1k_mapping.csv"
  )
)

readr::write_csv(
  GSE226643_condition_audit,
  file.path(
    GSE226643_audit_dir,
    "GSE226643_available_conditions_and_replicates.csv"
  )
)

readr::write_csv(
  GSE226643_primary_candidate_audit,
  file.path(
    GSE226643_audit_dir,
    "GSE226643_N_and_S_3h_12h_candidate_samples.csv"
  )
)

GSE226643_audit_text <- capture.output({
  cat("GSE226643 source-metadata / raw-count-column audit\n")
  cat("Generated: ", format(Sys.time()), "\n\n", sep = "")
  cat("The original leading blank header was repaired as: ",
      GSE226643_leading_source_column_name, "\n", sep = "")
  cat("It is retained as raw_count_column and must match counts.csv p1k IDs.\n\n")
  cat("Author repository root:\n", GSE226643_author_root, "\n\n", sep = "")
  cat("Available author conditions:\n")
  print(GSE226643_condition_audit, n = Inf)
  cat("\nCandidate N/S 3 h and 12 h sample rows:\n")
  print(GSE226643_primary_candidate_audit, n = Inf)
})

writeLines(
  GSE226643_audit_text,
  con = file.path(
    GSE226643_audit_dir,
    "GSE226643_source_metadata_audit.txt"
  )
)

print(GSE226643_condition_audit, n = Inf)

print(GSE226643_primary_candidate_audit, n = Inf)

message(
  "\nAudit completed successfully. Review these two files before creating the final manifests:\n",
  file.path(
    GSE226643_audit_dir,
    "GSE226643_available_conditions_and_replicates.csv"
  ),
  "\n",
  file.path(
    GSE226643_audit_dir,
    "GSE226643_N_and_S_3h_12h_candidate_samples.csv"
  )
)

GSE226643_unpacked_dir <- file.path(
  GSE226643_project_dir,
  "author_GitHub_repository",
  "unpacked"
)

GSE226643_manifest_dir <- file.path(
  GSE226643_project_dir,
  "analysis_manifests",
  "01_verified"
)

dir.create(
  GSE226643_manifest_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!dir.exists(GSE226643_project_dir)) {
  stop(
    "GSE226643_project_dir does not exist: ",
    GSE226643_project_dir,
    call. = FALSE
  )
}

if (!dir.exists(GSE226643_unpacked_dir)) {
  stop(
    "The unpacked author-repository directory does not exist: ",
    GSE226643_unpacked_dir,
    call. = FALSE
  )
}

GSE226643_author_root_candidates <- list.files(
  GSE226643_unpacked_dir,
  full.names = TRUE,
  recursive = FALSE,
  no.. = TRUE
)

GSE226643_author_root_candidates <- GSE226643_author_root_candidates[
  file.info(GSE226643_author_root_candidates)$isdir
]

GSE226643_author_root_candidates <- GSE226643_author_root_candidates[
  file.exists(
    file.path(
      GSE226643_author_root_candidates,
      "raw_data",
      "iM",
      "counts.csv"
    )
  )
]

if (length(GSE226643_author_root_candidates) != 1L) {
  stop(
    "Exactly one unpacked author repository containing raw_data/iM/counts.csv ",
    "must be present. Found ",
    length(GSE226643_author_root_candidates),
    ".",
    call. = FALSE
  )
}

GSE226643_author_root <- GSE226643_author_root_candidates[[1]]

GSE226643_counts_file <- file.path(
  GSE226643_author_root,
  "raw_data",
  "iM",
  "counts.csv"
)

GSE226643_sample_table_file <- file.path(
  GSE226643_author_root,
  "raw_data",
  "iM",
  "stationary_phase_sample_table_2.csv"
)

GSE226643_gene_info_file <- file.path(
  GSE226643_author_root,
  "raw_data",
  "MG1655_gene_info.csv"
)

GSE226643_required_source_files <- c(
  raw_count_matrix = GSE226643_counts_file,
  author_sample_table = GSE226643_sample_table_file,
  MG1655_gene_information = GSE226643_gene_info_file
)

GSE226643_missing_source_files <- GSE226643_required_source_files[
  !file.exists(GSE226643_required_source_files)
]

if (length(GSE226643_missing_source_files) > 0L) {
  stop(
    "Required source file(s) are missing:\n",
    paste(
      paste0(
        names(GSE226643_missing_source_files),
        ": ",
        GSE226643_missing_source_files
      ),
      collapse = "\n"
    ),
    call. = FALSE
  )
}

GSE226643_primary_manifest_expected <- tibble::tribble(
  ~raw_count_column, ~author_sample_id,       ~author_condition, ~author_replicate_id, ~analysis_id,             ~limitation,            ~time_h, ~DESeq2_stage,     ~stage_label,                            ~sample_role,
  "p1k_01068",      "N-GlcM9_1_3h__1",  "N-GlcM9_3h",    1L,                   "GSE226643_N_12h_vs_3h", "Nitrogen limitation", 3L,      "Early_3h",       "Early exponential reference (3 h)",    "reference",
  "p1k_01069",      "N-GlcM9_1_3h__2",  "N-GlcM9_3h",    2L,                   "GSE226643_N_12h_vs_3h", "Nitrogen limitation", 3L,      "Early_3h",       "Early exponential reference (3 h)",    "reference",
  "p1k_01076",      "N-GlcM9_5_12h__1", "N-GlcM9_12h",   1L,                   "GSE226643_N_12h_vs_3h", "Nitrogen limitation", 12L,     "Stationary_12h", "Established stationary phase (12 h)",  "stationary_phase",
  "p1k_01077",      "N-GlcM9_5_12h__2", "N-GlcM9_12h",   2L,                   "GSE226643_N_12h_vs_3h", "Nitrogen limitation", 12L,     "Stationary_12h", "Established stationary phase (12 h)",  "stationary_phase",
  "p1k_01080",      "S-GlcM9_1_3h__1",  "S-GlcM9_3h",    1L,                   "GSE226643_S_12h_vs_3h", "Sulfur limitation",   3L,      "Early_3h",       "Early exponential reference (3 h)",    "reference",
  "p1k_01081",      "S-GlcM9_1_3h__2",  "S-GlcM9_3h",    2L,                   "GSE226643_S_12h_vs_3h", "Sulfur limitation",   3L,      "Early_3h",       "Early exponential reference (3 h)",    "reference",
  "p1k_01088",      "S-GlcM9_5_12h__1", "S-GlcM9_12h",   1L,                   "GSE226643_S_12h_vs_3h", "Sulfur limitation",   12L,     "Stationary_12h", "Established stationary phase (12 h)",  "stationary_phase",
  "p1k_01089",      "S-GlcM9_5_12h__2", "S-GlcM9_12h",   2L,                   "GSE226643_S_12h_vs_3h", "Sulfur limitation",   12L,     "Stationary_12h", "Established stationary phase (12 h)",  "stationary_phase"
) %>%
  dplyr::mutate(
    contrast_numerator = "Stationary_12h",
    contrast_denominator = "Early_3h",
    DESeq2_design = "~ DESeq2_stage",
    primary_analysis = TRUE,
    biological_interpretation = paste0(
      "Established stationary phase (12 h) minus matched early-exponential ",
      "reference (3 h) under the same nutrient-limitation condition"
    ),
    stage_order = dplyr::case_when(
      .data$DESeq2_stage == "Early_3h" ~ 1L,
      .data$DESeq2_stage == "Stationary_12h" ~ 2L,
      TRUE ~ NA_integer_
    )
  ) %>%
  dplyr::arrange(
    .data$analysis_id,
    .data$stage_order,
    .data$author_replicate_id
  )

if (nrow(GSE226643_primary_manifest_expected) != 8L) {
  stop("The hard-locked primary design must contain eight samples.", call. = FALSE)
}

if (anyDuplicated(GSE226643_primary_manifest_expected$raw_count_column) > 0L) {
  stop("The hard-locked primary design contains duplicated p1k IDs.", call. = FALSE)
}

if (anyDuplicated(
  paste(
    GSE226643_primary_manifest_expected$analysis_id,
    GSE226643_primary_manifest_expected$DESeq2_stage,
    GSE226643_primary_manifest_expected$author_replicate_id,
    sep = "__"
  )
) > 0L) {
  stop(
    "The hard-locked primary design contains duplicated analysis/stage/replicate ",
    "combinations.",
    call. = FALSE
  )
}

readr::write_csv(
  GSE226643_primary_manifest_expected,
  file.path(
    GSE226643_manifest_dir,
    "GSE226643_primary_design_hard_locked.csv"
  )
)

GSE226643_sample_table_with_header <- readr::read_csv(
  GSE226643_sample_table_file,
  col_names = FALSE,
  show_col_types = FALSE,
  progress = FALSE
)

if (nrow(GSE226643_sample_table_with_header) < 2L ||
    ncol(GSE226643_sample_table_with_header) < 2L) {
  stop(
    "The author sample table must contain a p1k-ID column plus metadata columns.",
    call. = FALSE
  )
}

GSE226643_source_header <- as.character(
  unlist(
    GSE226643_sample_table_with_header[1, , drop = FALSE],
    use.names = FALSE
  )
)

GSE226643_source_header <- trimws(GSE226643_source_header)

GSE226643_original_leading_header <- GSE226643_source_header[1]

GSE226643_source_header[1] <- "raw_count_column"

if (anyNA(GSE226643_source_header) ||
    any(GSE226643_source_header == "") ||
    anyDuplicated(GSE226643_source_header) > 0L) {
  stop(
    "The author-table header could not be converted to unique, non-empty ",
    "column names.",
    call. = FALSE
  )
}

GSE226643_sample_table_raw <- GSE226643_sample_table_with_header[-1, , drop = FALSE]

colnames(GSE226643_sample_table_raw) <- GSE226643_source_header

GSE226643_required_metadata_columns <- c(
  "raw_count_column",
  "sample_id",
  "study",
  "project",
  "condition",
  "rep_id",
  "Strain",
  "Temperature (C)"
)

GSE226643_missing_metadata_columns <- setdiff(
  GSE226643_required_metadata_columns,
  colnames(GSE226643_sample_table_raw)
)

if (length(GSE226643_missing_metadata_columns) > 0L) {
  stop(
    "The author sample table lacks required column(s): ",
    paste(GSE226643_missing_metadata_columns, collapse = ", "),
    call. = FALSE
  )
}

GSE226643_sample_table <- GSE226643_sample_table_raw %>%
  dplyr::rename(
    author_sample_id = sample_id,
    author_study = study,
    author_project = project,
    author_condition = condition,
    author_replicate_id = rep_id
  ) %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::where(is.character),
      trimws
    ),
    raw_count_column = as.character(.data$raw_count_column),
    author_sample_id = as.character(.data$author_sample_id),
    author_study = as.character(.data$author_study),
    author_project = as.character(.data$author_project),
    author_condition = as.character(.data$author_condition),
    author_replicate_id = suppressWarnings(
      as.integer(.data$author_replicate_id)
    ),
    temperature_C = suppressWarnings(
      as.numeric(.data[["Temperature (C)"]])
    )
  )

if (anyNA(GSE226643_sample_table$raw_count_column) ||
    any(GSE226643_sample_table$raw_count_column == "")) {
  stop("The author metadata contains missing p1k IDs.", call. = FALSE)
}

if (anyDuplicated(GSE226643_sample_table$raw_count_column) > 0L) {
  stop("The author metadata contains duplicated p1k IDs.", call. = FALSE)
}

if (!all(grepl("^p1k_[0-9]{5}$", GSE226643_sample_table$raw_count_column))) {
  stop(
    "The leading author-metadata column does not contain the expected p1k IDs.",
    call. = FALSE
  )
}

GSE226643_source_primary_rows <- GSE226643_sample_table %>%
  dplyr::filter(
    .data$raw_count_column %in%
      GSE226643_primary_manifest_expected$raw_count_column
  ) %>%
  dplyr::select(
    .data$raw_count_column,
    .data$author_sample_id,
    .data$author_study,
    .data$author_project,
    .data$author_condition,
    .data$author_replicate_id,
    .data$Strain,
    .data$temperature_C,
    dplyr::any_of(c(
      "Base Media",
      "pH",
      "Carbon Source (g/L)",
      "Nitrogen Source (g/L)",
      "Trace Element Mixture",
      "Supplement",
      "Additional Details",
      "GEO",
      "SRX",
      "Run",
      "full_name"
    ))
  )

if (nrow(GSE226643_source_primary_rows) != 8L) {
  stop(
    "The author metadata did not return all eight hard-locked primary samples.",
    call. = FALSE
  )
}

GSE226643_metadata_match_audit <- GSE226643_primary_manifest_expected %>%
  dplyr::select(
    .data$raw_count_column,
    expected_author_sample_id = .data$author_sample_id,
    expected_author_condition = .data$author_condition,
    expected_author_replicate_id = .data$author_replicate_id
  ) %>%
  dplyr::left_join(
    GSE226643_source_primary_rows %>%
      dplyr::select(
        .data$raw_count_column,
        actual_author_sample_id = .data$author_sample_id,
        actual_author_condition = .data$author_condition,
        actual_author_replicate_id = .data$author_replicate_id,
        actual_author_study = .data$author_study,
        actual_author_project = .data$author_project,
        actual_strain = .data$Strain,
        actual_temperature_C = .data$temperature_C
      ),
    by = "raw_count_column"
  ) %>%
  dplyr::arrange(.data$raw_count_column)

GSE226643_metadata_mismatch <- GSE226643_metadata_match_audit %>%
  dplyr::filter(
    is.na(.data$actual_author_sample_id) |
      .data$expected_author_sample_id != .data$actual_author_sample_id |
      .data$expected_author_condition != .data$actual_author_condition |
      .data$expected_author_replicate_id != .data$actual_author_replicate_id |
      .data$actual_author_study != "Stationary_phase" |
      .data$actual_author_project != "stationary" |
      .data$actual_strain != "MG1655" |
      .data$actual_temperature_C != 37
  )

if (nrow(GSE226643_metadata_mismatch) > 0L) {
  print(GSE226643_metadata_mismatch)
  stop(
    "At least one hard-locked primary sample does not exactly match the audited ",
    "author metadata. Do not continue until this discrepancy is resolved.",
    call. = FALSE
  )
}

GSE226643_counts_header <- readr::read_csv(
  GSE226643_counts_file,
  n_max = 0L,
  show_col_types = FALSE,
  progress = FALSE,
  name_repair = "minimal"
)

if (!identical(colnames(GSE226643_counts_header)[1], "Geneid")) {
  stop(
    "The first counts.csv column must be named Geneid.",
    call. = FALSE
  )
}

GSE226643_count_columns <- colnames(GSE226643_counts_header)[-1]

GSE226643_primary_p1k_missing_from_counts <- setdiff(
  GSE226643_primary_manifest_expected$raw_count_column,
  GSE226643_count_columns
)

if (length(GSE226643_primary_p1k_missing_from_counts) > 0L) {
  stop(
    "The following hard-locked p1k ID(s) are absent from counts.csv: ",
    paste(GSE226643_primary_p1k_missing_from_counts, collapse = ", "),
    call. = FALSE
  )
}

GSE226643_primary_manifest_verified <- GSE226643_source_primary_rows %>%
  dplyr::inner_join(
    GSE226643_primary_manifest_expected %>%
      dplyr::select(
        .data$raw_count_column,
        .data$analysis_id,
        .data$limitation,
        .data$time_h,
        .data$DESeq2_stage,
        .data$stage_label,
        .data$sample_role,
        .data$contrast_numerator,
        .data$contrast_denominator,
        .data$DESeq2_design,
        .data$primary_analysis,
        .data$biological_interpretation,
        .data$stage_order
      ),
    by = "raw_count_column"
  ) %>%
  dplyr::mutate(
    DESeq2_stage = factor(
      .data$DESeq2_stage,
      levels = c("Early_3h", "Stationary_12h")
    ),
    manifest_status = "verified_from_audited_author_metadata",
    manifest_frozen_date = as.character(Sys.Date()),
    raw_count_file = GSE226643_counts_file,
    source_sample_table = GSE226643_sample_table_file
  ) %>%
  dplyr::arrange(
    .data$analysis_id,
    .data$stage_order,
    .data$author_replicate_id
  )

if (nrow(GSE226643_primary_manifest_verified) != 8L) {
  stop("The verified primary manifest must contain eight samples.", call. = FALSE)
}

GSE226643_replicate_audit <- GSE226643_primary_manifest_verified %>%
  dplyr::group_by(
    .data$analysis_id,
    .data$limitation,
    .data$author_condition,
    .data$time_h,
    .data$DESeq2_stage
  ) %>%
  dplyr::summarise(
    N_samples = dplyr::n(),
    N_unique_replicates = dplyr::n_distinct(.data$author_replicate_id),
    replicate_ids = paste(
      sort(.data$author_replicate_id),
      collapse = ";"
    ),
    raw_count_columns = paste(
      .data$raw_count_column,
      collapse = ";"
    ),
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    .data$analysis_id,
    .data$time_h
  )

if (any(GSE226643_replicate_audit$N_samples != 2L) ||
    any(GSE226643_replicate_audit$N_unique_replicates != 2L) ||
    any(GSE226643_replicate_audit$replicate_ids != "1;2")) {
  stop(
    "Each primary condition must contain exactly biological replicates 1 and 2.",
    call. = FALSE
  )
}

GSE226643_N_manifest_verified <- GSE226643_primary_manifest_verified %>%
  dplyr::filter(
    .data$analysis_id == "GSE226643_N_12h_vs_3h"
  )

GSE226643_S_manifest_verified <- GSE226643_primary_manifest_verified %>%
  dplyr::filter(
    .data$analysis_id == "GSE226643_S_12h_vs_3h"
  )

if (nrow(GSE226643_N_manifest_verified) != 4L ||
    nrow(GSE226643_S_manifest_verified) != 4L) {
  stop(
    "Each verified analysis-specific manifest must contain four samples.",
    call. = FALSE
  )
}

GSE226643_counts <- readr::read_csv(
  GSE226643_counts_file,
  show_col_types = FALSE,
  progress = FALSE,
  name_repair = "minimal"
)

if (!identical(colnames(GSE226643_counts)[1], "Geneid")) {
  stop("The first counts.csv column must be Geneid.", call. = FALSE)
}

GSE226643_gene_ids <- as.character(GSE226643_counts$Geneid)

if (anyNA(GSE226643_gene_ids) || any(GSE226643_gene_ids == "")) {
  stop("counts.csv contains missing Geneid values.", call. = FALSE)
}

if (anyDuplicated(GSE226643_gene_ids) > 0L) {
  stop("counts.csv contains duplicated Geneid values.", call. = FALSE)
}

GSE226643_selected_count_matrix <- as.matrix(
  GSE226643_counts[
    ,
    GSE226643_primary_manifest_verified$raw_count_column,
    drop = FALSE
  ]
)

storage.mode(GSE226643_selected_count_matrix) <- "numeric"

rownames(GSE226643_selected_count_matrix) <- GSE226643_gene_ids

if (any(!is.finite(GSE226643_selected_count_matrix))) {
  stop("The selected count matrix contains non-finite values.", call. = FALSE)
}

if (any(GSE226643_selected_count_matrix < 0)) {
  stop("The selected count matrix contains negative values.", call. = FALSE)
}

if (any(abs(
  GSE226643_selected_count_matrix - round(GSE226643_selected_count_matrix)
) > 1e-8)) {
  stop(
    "counts.csv does not contain integer raw counts for every selected sample.",
    call. = FALSE
  )
}

GSE226643_library_size_audit <- tibble::tibble(
  raw_count_column = colnames(GSE226643_selected_count_matrix),
  raw_library_size = colSums(GSE226643_selected_count_matrix),
  detected_genes_raw_count_gt_0 = colSums(
    GSE226643_selected_count_matrix > 0
  ),
  proportion_zero_counts = colMeans(
    GSE226643_selected_count_matrix == 0
  )
) %>%
  dplyr::left_join(
    GSE226643_primary_manifest_verified %>%
      dplyr::select(
        .data$raw_count_column,
        .data$analysis_id,
        .data$limitation,
        .data$author_condition,
        .data$time_h,
        .data$DESeq2_stage,
        .data$author_replicate_id
      ),
    by = "raw_count_column"
  ) %>%
  dplyr::arrange(
    .data$analysis_id,
    .data$DESeq2_stage,
    .data$author_replicate_id
  )

if (any(GSE226643_library_size_audit$raw_library_size <= 0)) {
  stop(
    "At least one selected sample has a zero raw library size.",
    call. = FALSE
  )
}

GSE226643_log_count_correlation <- stats::cor(
  log2(GSE226643_selected_count_matrix + 1),
  method = "pearson"
)

GSE226643_log_count_correlation_long <- as.data.frame(
  as.table(GSE226643_log_count_correlation),
  stringsAsFactors = FALSE
) %>%
  tibble::as_tibble() %>%
  dplyr::rename(
    raw_count_column_row = Var1,
    raw_count_column_column = Var2,
    Pearson_correlation_log2_count_plus_1 = Freq
  )

GSE226643_source_file_audit <- tibble::tibble(
  source_name = names(GSE226643_required_source_files),
  file = unname(GSE226643_required_source_files),
  size_bytes = file.info(unname(GSE226643_required_source_files))$size,
  md5 = unname(tools::md5sum(unname(GSE226643_required_source_files)))
)

readr::write_csv(
  GSE226643_primary_manifest_verified,
  file.path(
    GSE226643_manifest_dir,
    "GSE226643_primary_manifest_VERIFIED.csv"
  )
)

readr::write_csv(
  GSE226643_N_manifest_verified,
  file.path(
    GSE226643_manifest_dir,
    "GSE226643_N_12h_vs_3h_manifest_VERIFIED.csv"
  )
)

readr::write_csv(
  GSE226643_S_manifest_verified,
  file.path(
    GSE226643_manifest_dir,
    "GSE226643_S_12h_vs_3h_manifest_VERIFIED.csv"
  )
)

readr::write_csv(
  GSE226643_metadata_match_audit,
  file.path(
    GSE226643_manifest_dir,
    "GSE226643_hard_locked_metadata_match_audit.csv"
  )
)

readr::write_csv(
  GSE226643_replicate_audit,
  file.path(
    GSE226643_manifest_dir,
    "GSE226643_verified_replicate_audit.csv"
  )
)

readr::write_csv(
  GSE226643_library_size_audit,
  file.path(
    GSE226643_manifest_dir,
    "GSE226643_verified_selected_library_size_audit.csv"
  )
)

readr::write_csv(
  GSE226643_log_count_correlation_long,
  file.path(
    GSE226643_manifest_dir,
    "GSE226643_verified_selected_log_count_correlation.csv"
  )
)

readr::write_csv(
  GSE226643_source_file_audit,
  file.path(
    GSE226643_manifest_dir,
    "GSE226643_verified_source_file_audit.csv"
  )
)

GSE226643_verified_manifest_audit_text <- capture.output({
  cat("GSE226643 verified primary-manifest audit\n")
  cat("Generated: ", format(Sys.time()), "\n\n", sep = "")
  cat("Original repaired leading author-table header: ",
      GSE226643_original_leading_header, "\n", sep = "")
  cat("Author repository root:\n", GSE226643_author_root, "\n\n", sep = "")
  cat("Verified primary manifest:\n")
  print(
    GSE226643_primary_manifest_verified %>%
      dplyr::select(
        .data$analysis_id,
        .data$limitation,
        .data$raw_count_column,
        .data$author_sample_id,
        .data$author_condition,
        .data$time_h,
        .data$DESeq2_stage,
        .data$author_replicate_id,
        .data$Strain,
        .data$temperature_C
      ),
    n = Inf
  )
  cat("\nReplicate audit:\n")
  print(GSE226643_replicate_audit, n = Inf)
  cat("\nRaw-count library-size audit:\n")
  print(GSE226643_library_size_audit, n = Inf)
  cat("\nSource-file checksums:\n")
  print(GSE226643_source_file_audit, n = Inf)
})

writeLines(
  GSE226643_verified_manifest_audit_text,
  con = file.path(
    GSE226643_manifest_dir,
    "GSE226643_verified_manifest_audit.txt"
  )
)

GSE226643_library_size_plot <- ggplot(
  GSE226643_library_size_audit,
  aes(
    x = .data$raw_count_column,
    y = .data$raw_library_size,
    fill = .data$DESeq2_stage
  )
) +
  geom_col(
    colour = "black",
    linewidth = 0.2
  ) +
  facet_wrap(~ analysis_id, scales = "free_x") +
  scale_fill_manual(
    values = c(
      "Early_3h" = "#8cc5e3",
      "Stationary_12h" = "#2066a8"
    )
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      colour = "black",
      size = 9
    ),
    legend.position = "top"
  ) +
  labs(
    x = "Raw-count p1k column",
    y = "Raw library size",
    fill = "Stage"
  )

ggsave(
  filename = file.path(
    GSE226643_manifest_dir,
    "GSE226643_verified_selected_library_size_audit.png"
  ),
  plot = GSE226643_library_size_plot,
  width = 10,
  height = 5,
  units = "in",
  dpi = 300
)

GSE226643_log_count_correlation_plot <- ggplot(
  GSE226643_log_count_correlation_long,
  aes(
    x = .data$raw_count_column_column,
    y = .data$raw_count_column_row,
    fill = .data$Pearson_correlation_log2_count_plus_1
  )
) +
  geom_tile(
    colour = "white",
    linewidth = 0.2
  ) +
  geom_text(
    aes(
      label = format(
        round(.data$Pearson_correlation_log2_count_plus_1, 2),
        nsmall = 2
      )
    ),
    size = 3
  ) +
  scale_fill_gradient(
    low = "white",
    high = "#2066a8",
    limits = c(0, 1),
    name = "Pearson r"
  ) +
  coord_equal() +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      colour = "black",
      size = 9
    ),
    axis.text.y = element_text(
      colour = "black",
      size = 9
    )
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "Verified primary samples: correlation of log2(raw count + 1)"
  )

ggsave(
  filename = file.path(
    GSE226643_manifest_dir,
    "GSE226643_verified_selected_log_count_correlation.png"
  ),
  plot = GSE226643_log_count_correlation_plot,
  width = 7,
  height = 6,
  units = "in",
  dpi = 300
)

saveRDS(
  list(
    primary_manifest = GSE226643_primary_manifest_verified,
    nitrogen_manifest = GSE226643_N_manifest_verified,
    sulfur_manifest = GSE226643_S_manifest_verified,
    hard_locked_design = GSE226643_primary_manifest_expected,
    metadata_match_audit = GSE226643_metadata_match_audit,
    replicate_audit = GSE226643_replicate_audit,
    library_size_audit = GSE226643_library_size_audit,
    source_file_audit = GSE226643_source_file_audit,
    raw_count_file = GSE226643_counts_file,
    author_sample_table_file = GSE226643_sample_table_file,
    MG1655_gene_info_file = GSE226643_gene_info_file
  ),
  file = file.path(
    GSE226643_manifest_dir,
    "GSE226643_verified_analysis_manifests.rds"
  )
)

message(
  "Verified GSE226643 manifests were written to:\n",
  GSE226643_manifest_dir,
  "\n\nThe next script can now perform the N- and S-limitation DESeq2 and ",
  "thermodynamic analyses using only GSE226643_verified_analysis_manifests.rds."
)

GSE226643_B_BOOT <- 5000L
GSE226643_B_PERM <- 5000L

GSE226643_refresh_GO_cache <- FALSE

GSE226643_manifest_file <- file.path(
  GSE226643_project_dir,
  "analysis_manifests",
  "01_verified",
  "GSE226643_verified_analysis_manifests.rds"
)

GSE226643_output_dir <- file.path(
  GSE226643_project_dir,
  "analysis_results",
  "02_N_S_stationary_phase_thermodynamic_analysis"
)

GSE226643_table_dir <- file.path(GSE226643_output_dir, "tables")
GSE226643_plot_dir <- file.path(GSE226643_output_dir, "plots")
GSE226643_qc_dir <- file.path(GSE226643_output_dir, "quality_control")
GSE226643_gsea_dir <- file.path(GSE226643_output_dir, "GSEA")

dir.create(GSE226643_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(GSE226643_plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(GSE226643_qc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(GSE226643_gsea_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(GSE226643_manifest_file)) {
  stop(
    "The verified manifest RDS is missing. Run the v2 manifest script first: ",
    GSE226643_manifest_file,
    call. = FALSE
  )
}

if (!file.exists(GSE226643_prediction_file)) {
  stop(
    "The MG1655 ΔGfold prediction file is missing: ",
    GSE226643_prediction_file,
    call. = FALSE
  )
}

GSE226643_manifest_bundle <- readRDS(GSE226643_manifest_file)

GSE226643_required_manifest_objects <- c(
  "primary_manifest",
  "nitrogen_manifest",
  "sulfur_manifest",
  "raw_count_file",
  "MG1655_gene_info_file"
)

GSE226643_missing_manifest_objects <- setdiff(
  GSE226643_required_manifest_objects,
  names(GSE226643_manifest_bundle)
)

if (length(GSE226643_missing_manifest_objects) > 0L) {
  stop(
    "The verified manifest RDS lacks required object(s): ",
    paste(GSE226643_missing_manifest_objects, collapse = ", "),
    call. = FALSE
  )
}

GSE226643_primary_manifest <- GSE226643_manifest_bundle$primary_manifest %>%
  tibble::as_tibble() %>%
  dplyr::mutate(
    raw_count_column = as.character(.data$raw_count_column),
    analysis_id = as.character(.data$analysis_id),
    limitation = as.character(.data$limitation),
    author_condition = as.character(.data$author_condition),
    author_replicate_id = as.integer(.data$author_replicate_id),
    DESeq2_stage = as.character(.data$DESeq2_stage),
    time_h = as.integer(.data$time_h)
  )

GSE226643_N_manifest <- GSE226643_manifest_bundle$nitrogen_manifest %>%
  tibble::as_tibble() %>%
  dplyr::mutate(
    raw_count_column = as.character(.data$raw_count_column),
    DESeq2_stage = as.character(.data$DESeq2_stage),
    author_replicate_id = as.integer(.data$author_replicate_id)
  ) %>%
  dplyr::arrange(.data$stage_order, .data$author_replicate_id)

GSE226643_S_manifest <- GSE226643_manifest_bundle$sulfur_manifest %>%
  tibble::as_tibble() %>%
  dplyr::mutate(
    raw_count_column = as.character(.data$raw_count_column),
    DESeq2_stage = as.character(.data$DESeq2_stage),
    author_replicate_id = as.integer(.data$author_replicate_id)
  ) %>%
  dplyr::arrange(.data$stage_order, .data$author_replicate_id)

GSE226643_required_manifest_columns <- c(
  "raw_count_column",
  "analysis_id",
  "limitation",
  "author_condition",
  "author_replicate_id",
  "time_h",
  "DESeq2_stage",
  "stage_order"
)

if (!all(GSE226643_required_manifest_columns %in% names(GSE226643_primary_manifest))) {
  stop(
    "The primary manifest is missing required analysis columns.",
    call. = FALSE
  )
}

if (nrow(GSE226643_primary_manifest) != 8L ||
    nrow(GSE226643_N_manifest) != 4L ||
    nrow(GSE226643_S_manifest) != 4L) {
  stop(
    "Expected an eight-sample combined manifest and two four-sample manifests.",
    call. = FALSE
  )
}

GSE226643_expected_primary_p1k <- c(
  "p1k_01068", "p1k_01069", "p1k_01076", "p1k_01077",
  "p1k_01080", "p1k_01081", "p1k_01088", "p1k_01089"
)

if (!setequal(
  GSE226643_primary_manifest$raw_count_column,
  GSE226643_expected_primary_p1k
)) {
  stop(
    "The verified manifest does not contain the eight audited GSE226643 p1k ",
    "raw-count columns.",
    call. = FALSE
  )
}

GSE226643_manifest_stage_audit <- GSE226643_primary_manifest %>%
  dplyr::count(
    .data$analysis_id,
    .data$limitation,
    .data$author_condition,
    .data$time_h,
    .data$DESeq2_stage,
    name = "N_biological_replicates"
  ) %>%
  dplyr::arrange(.data$analysis_id, .data$time_h)

if (any(GSE226643_manifest_stage_audit$N_biological_replicates != 2L)) {
  stop(
    "Each analysis/stage combination must contain exactly two biological replicates.",
    call. = FALSE
  )
}

print(GSE226643_manifest_stage_audit)

readr::write_tsv(
  GSE226643_primary_manifest,
  file.path(GSE226643_table_dir, "GSE226643_primary_manifest_used.tsv")
)

readr::write_tsv(
  GSE226643_manifest_stage_audit,
  file.path(GSE226643_table_dir, "GSE226643_manifest_stage_audit.tsv")
)

GSE226643_counts_file <- GSE226643_manifest_bundle$raw_count_file
GSE226643_gene_info_file <- GSE226643_manifest_bundle$MG1655_gene_info_file

if (!file.exists(GSE226643_counts_file)) {
  stop("The verified raw-count file is missing: ", GSE226643_counts_file, call. = FALSE)
}

if (!file.exists(GSE226643_gene_info_file)) {
  stop("The MG1655 gene-information file is missing: ", GSE226643_gene_info_file, call. = FALSE)
}

GSE226643_counts <- readr::read_csv(
  GSE226643_counts_file,
  show_col_types = FALSE,
  progress = FALSE,
  name_repair = "minimal"
)

if (!identical(colnames(GSE226643_counts)[1], "Geneid")) {
  stop("The first counts.csv column must be Geneid.", call. = FALSE)
}

GSE226643_count_sample_columns <- colnames(GSE226643_counts)[-1]
GSE226643_gene_ids <- tolower(trimws(as.character(GSE226643_counts$Geneid)))

if (anyNA(GSE226643_gene_ids) || any(GSE226643_gene_ids == "")) {
  stop("counts.csv contains missing Geneid values.", call. = FALSE)
}

if (anyDuplicated(GSE226643_gene_ids) > 0L) {
  stop("counts.csv contains duplicated Geneid values.", call. = FALSE)
}

if (!all(GSE226643_primary_manifest$raw_count_column %in% GSE226643_count_sample_columns)) {
  stop(
    "At least one verified p1k sample is absent from counts.csv.",
    call. = FALSE
  )
}

GSE226643_raw_count_matrix <- as.matrix(
  GSE226643_counts[, GSE226643_primary_manifest$raw_count_column, drop = FALSE]
)

storage.mode(GSE226643_raw_count_matrix) <- "numeric"
rownames(GSE226643_raw_count_matrix) <- GSE226643_gene_ids

if (any(!is.finite(GSE226643_raw_count_matrix)) ||
    any(GSE226643_raw_count_matrix < 0) ||
    any(abs(GSE226643_raw_count_matrix - round(GSE226643_raw_count_matrix)) > 1e-8)) {
  stop(
    "The selected GSE226643 expression values are not non-negative integer raw counts.",
    call. = FALSE
  )
}

storage.mode(GSE226643_raw_count_matrix) <- "integer"

GSE226643_gene_info <- readr::read_csv(
  GSE226643_gene_info_file,
  show_col_types = FALSE,
  progress = FALSE
)

GSE226643_required_gene_info_columns <- c(
  "locus_tag",
  "gene_name",
  "start",
  "end"
)

GSE226643_missing_gene_info_columns <- setdiff(
  GSE226643_required_gene_info_columns,
  colnames(GSE226643_gene_info)
)

if (length(GSE226643_missing_gene_info_columns) > 0L) {
  stop(
    "MG1655_gene_info.csv lacks required column(s): ",
    paste(GSE226643_missing_gene_info_columns, collapse = ", "),
    call. = FALSE
  )
}

GSE226643_gene_annotation <- GSE226643_gene_info %>%
  dplyr::transmute(
    Gene_id = tolower(trimws(as.character(.data$locus_tag))),
    gene_symbol = toupper(trimws(as.character(.data$gene_name))),
    gene_length_nt = abs(
      as.numeric(.data$end) - as.numeric(.data$start)
    ) + 1
  ) %>%
  dplyr::filter(!is.na(.data$Gene_id), .data$Gene_id != "") %>%
  dplyr::distinct(.data$Gene_id, .keep_all = TRUE)

if (anyDuplicated(GSE226643_gene_annotation$Gene_id) > 0L) {
  stop("The MG1655 gene annotation contains duplicated locus tags.", call. = FALSE)
}

GSE226643_prediction_raw <- readr::read_csv(
  GSE226643_prediction_file,
  show_col_types = FALSE,
  progress = FALSE
)

GSE226643_required_prediction_columns <- c("Gene_id", "Normalized_MFE")
GSE226643_missing_prediction_columns <- setdiff(
  GSE226643_required_prediction_columns,
  colnames(GSE226643_prediction_raw)
)

if (length(GSE226643_missing_prediction_columns) > 0L) {
  stop(
    "The prediction table lacks required column(s): ",
    paste(GSE226643_missing_prediction_columns, collapse = ", "),
    call. = FALSE
  )
}

GSE226643_prediction <- GSE226643_prediction_raw %>%
  dplyr::transmute(
    Gene_id = tolower(trimws(as.character(.data$Gene_id))),
    Normalized_MFE = as.numeric(.data$Normalized_MFE)
  ) %>%
  dplyr::filter(!is.na(.data$Gene_id), .data$Gene_id != "") %>%
  dplyr::distinct(.data$Gene_id, .keep_all = TRUE)

GSE226643_thermodynamic_annotation_audit <- tibble::tibble(
  Gene_id = rownames(GSE226643_raw_count_matrix)
) %>%
  dplyr::left_join(GSE226643_gene_annotation, by = "Gene_id") %>%
  dplyr::left_join(GSE226643_prediction, by = "Gene_id") %>%
  dplyr::mutate(
    finite_Normalized_MFE = is.finite(.data$Normalized_MFE),
    positive_gene_length = is.finite(.data$gene_length_nt) &
      .data$gene_length_nt > 0,
    included_in_DE_descriptor_analysis = .data$finite_Normalized_MFE,
    included_in_full_transcriptome_weighted_analysis =
      .data$finite_Normalized_MFE &
      .data$positive_gene_length,
    exclusion_reason_weighted_analysis = dplyr::case_when(
      .data$included_in_full_transcriptome_weighted_analysis ~ "",
      !.data$positive_gene_length &
        !.data$finite_Normalized_MFE ~
        "Invalid gene length; missing/non-finite Normalized_MFE",
      !.data$positive_gene_length ~ "Invalid or missing gene length",
      !.data$finite_Normalized_MFE ~ "Missing/non-finite Normalized_MFE"
    )
  )

GSE226643_DE_descriptor_annotation <- GSE226643_thermodynamic_annotation_audit %>%
  dplyr::filter(.data$included_in_DE_descriptor_analysis)

GSE226643_weighted_annotation <- GSE226643_thermodynamic_annotation_audit %>%
  dplyr::filter(.data$included_in_full_transcriptome_weighted_analysis) %>%
  dplyr::mutate(
    raw_matrix_row = match(.data$Gene_id, rownames(GSE226643_raw_count_matrix))
  ) %>%
  dplyr::arrange(.data$raw_matrix_row)

if (nrow(GSE226643_DE_descriptor_annotation) < 3500L) {
  stop(
    "Too few genes have finite ΔGfold predictions for the descriptor analyses.",
    call. = FALSE
  )
}

if (nrow(GSE226643_weighted_annotation) < 3500L) {
  stop(
    "Too few genes have valid length and ΔGfold annotation for the ",
    "full-transcriptome weighted analyses.",
    call. = FALSE
  )
}

readr::write_tsv(
  GSE226643_thermodynamic_annotation_audit,
  file.path(GSE226643_table_dir, "GSE226643_thermodynamic_annotation_audit.tsv")
)

GSE226643_N_count_matrix <- GSE226643_raw_count_matrix[
  ,
  GSE226643_N_manifest$raw_count_column,
  drop = FALSE
]

GSE226643_N_nonzero_genes <- rowSums(GSE226643_N_count_matrix) > 0L

if (sum(GSE226643_N_nonzero_genes) < 1000L) {
  stop("Too few non-zero genes are available for the nitrogen DESeq2 analysis.", call. = FALSE)
}

GSE226643_N_count_matrix <- GSE226643_N_count_matrix[
  GSE226643_N_nonzero_genes,
  ,
  drop = FALSE
]

GSE226643_N_coldata <- GSE226643_N_manifest %>%
  dplyr::mutate(
    DESeq2_stage = factor(
      .data$DESeq2_stage,
      levels = c("Early_3h", "Stationary_12h")
    )
  ) %>%
  dplyr::arrange(match(.data$raw_count_column, colnames(GSE226643_N_count_matrix))) %>%
  dplyr::select(
    .data$raw_count_column,
    .data$DESeq2_stage,
    .data$author_replicate_id,
    .data$author_condition
  ) %>%
  as.data.frame()

rownames(GSE226643_N_coldata) <- GSE226643_N_coldata$raw_count_column
GSE226643_N_coldata$raw_count_column <- NULL

if (!identical(colnames(GSE226643_N_count_matrix), rownames(GSE226643_N_coldata))) {
  stop("Nitrogen count-matrix columns and DESeq2 colData rows are not aligned.", call. = FALSE)
}

GSE226643_N_dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = GSE226643_N_count_matrix,
  colData = GSE226643_N_coldata,
  design = ~ DESeq2_stage
)

GSE226643_N_dds <- DESeq2::DESeq(
  GSE226643_N_dds,
  quiet = TRUE
)

GSE226643_N_DESeq2_results <- DESeq2::results(
  GSE226643_N_dds,
  contrast = c("DESeq2_stage", "Stationary_12h", "Early_3h"),
  alpha = 0.05,
  independentFiltering = TRUE,
  cooksCutoff = TRUE
)

GSE226643_N_normalised_counts <- DESeq2::counts(
  GSE226643_N_dds,
  normalized = TRUE
)

GSE226643_N_DE_all <- tibble::as_tibble(
  as.data.frame(GSE226643_N_DESeq2_results),
  rownames = "Gene_id"
) %>%
  dplyr::transmute(
    Gene_id = tolower(as.character(.data$Gene_id)),
    base_mean = as.numeric(.data$baseMean),
    log2FC_stationary_12h_minus_early_3h = as.numeric(.data$log2FoldChange),
    lfc_standard_error = as.numeric(.data$lfcSE),
    wald_statistic = as.numeric(.data$stat),
    p_value = as.numeric(.data$pvalue),
    padj = as.numeric(.data$padj),
    early_3h_mean_normalised_counts = rowMeans(
      GSE226643_N_normalised_counts[
        ,
        GSE226643_N_coldata$DESeq2_stage == "Early_3h",
        drop = FALSE
      ]
    ),
    stationary_12h_mean_normalised_counts = rowMeans(
      GSE226643_N_normalised_counts[
        ,
        GSE226643_N_coldata$DESeq2_stage == "Stationary_12h",
        drop = FALSE
      ]
    ),
    analysis_id = "GSE226643_N_12h_vs_3h",
    limitation = "Nitrogen limitation",
    comparison_label = "N limitation: 12 h stationary minus 3 h early exponential"
  )

GSE226643_N_DE_thermodynamic <- GSE226643_N_DE_all %>%
  dplyr::left_join(
    GSE226643_thermodynamic_annotation_audit %>%
      dplyr::select(
        .data$Gene_id,
        .data$gene_symbol,
        .data$gene_length_nt,
        .data$Normalized_MFE,
        .data$included_in_DE_descriptor_analysis,
        .data$included_in_full_transcriptome_weighted_analysis
      ),
    by = "Gene_id"
  )

GSE226643_S_count_matrix <- GSE226643_raw_count_matrix[
  ,
  GSE226643_S_manifest$raw_count_column,
  drop = FALSE
]

GSE226643_S_nonzero_genes <- rowSums(GSE226643_S_count_matrix) > 0L

if (sum(GSE226643_S_nonzero_genes) < 1000L) {
  stop("Too few non-zero genes are available for the sulfur DESeq2 analysis.", call. = FALSE)
}

GSE226643_S_count_matrix <- GSE226643_S_count_matrix[
  GSE226643_S_nonzero_genes,
  ,
  drop = FALSE
]

GSE226643_S_coldata <- GSE226643_S_manifest %>%
  dplyr::mutate(
    DESeq2_stage = factor(
      .data$DESeq2_stage,
      levels = c("Early_3h", "Stationary_12h")
    )
  ) %>%
  dplyr::arrange(match(.data$raw_count_column, colnames(GSE226643_S_count_matrix))) %>%
  dplyr::select(
    .data$raw_count_column,
    .data$DESeq2_stage,
    .data$author_replicate_id,
    .data$author_condition
  ) %>%
  as.data.frame()

rownames(GSE226643_S_coldata) <- GSE226643_S_coldata$raw_count_column
GSE226643_S_coldata$raw_count_column <- NULL

if (!identical(colnames(GSE226643_S_count_matrix), rownames(GSE226643_S_coldata))) {
  stop("Sulfur count-matrix columns and DESeq2 colData rows are not aligned.", call. = FALSE)
}

GSE226643_S_dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = GSE226643_S_count_matrix,
  colData = GSE226643_S_coldata,
  design = ~ DESeq2_stage
)

GSE226643_S_dds <- DESeq2::DESeq(
  GSE226643_S_dds,
  quiet = TRUE
)

GSE226643_S_DESeq2_results <- DESeq2::results(
  GSE226643_S_dds,
  contrast = c("DESeq2_stage", "Stationary_12h", "Early_3h"),
  alpha = 0.05,
  independentFiltering = TRUE,
  cooksCutoff = TRUE
)

GSE226643_S_normalised_counts <- DESeq2::counts(
  GSE226643_S_dds,
  normalized = TRUE
)

GSE226643_S_DE_all <- tibble::as_tibble(
  as.data.frame(GSE226643_S_DESeq2_results),
  rownames = "Gene_id"
) %>%
  dplyr::transmute(
    Gene_id = tolower(as.character(.data$Gene_id)),
    base_mean = as.numeric(.data$baseMean),
    log2FC_stationary_12h_minus_early_3h = as.numeric(.data$log2FoldChange),
    lfc_standard_error = as.numeric(.data$lfcSE),
    wald_statistic = as.numeric(.data$stat),
    p_value = as.numeric(.data$pvalue),
    padj = as.numeric(.data$padj),
    early_3h_mean_normalised_counts = rowMeans(
      GSE226643_S_normalised_counts[
        ,
        GSE226643_S_coldata$DESeq2_stage == "Early_3h",
        drop = FALSE
      ]
    ),
    stationary_12h_mean_normalised_counts = rowMeans(
      GSE226643_S_normalised_counts[
        ,
        GSE226643_S_coldata$DESeq2_stage == "Stationary_12h",
        drop = FALSE
      ]
    ),
    analysis_id = "GSE226643_S_12h_vs_3h",
    limitation = "Sulfur limitation",
    comparison_label = "S limitation: 12 h stationary minus 3 h early exponential"
  )

GSE226643_S_DE_thermodynamic <- GSE226643_S_DE_all %>%
  dplyr::left_join(
    GSE226643_thermodynamic_annotation_audit %>%
      dplyr::select(
        .data$Gene_id,
        .data$gene_symbol,
        .data$gene_length_nt,
        .data$Normalized_MFE,
        .data$included_in_DE_descriptor_analysis,
        .data$included_in_full_transcriptome_weighted_analysis
      ),
    by = "Gene_id"
  )

GSE226643_N_sig_thermodynamic <- GSE226643_N_DE_thermodynamic %>%
  dplyr::filter(
    !is.na(.data$padj),
    .data$padj < 0.05,
    is.finite(.data$log2FC_stationary_12h_minus_early_3h),
    .data$log2FC_stationary_12h_minus_early_3h != 0,
    .data$included_in_DE_descriptor_analysis %in% TRUE
  ) %>%
  dplyr::mutate(
    FC_sign = dplyr::case_when(
      .data$log2FC_stationary_12h_minus_early_3h < 0 ~ "log2_FC < 0",
      .data$log2FC_stationary_12h_minus_early_3h > 0 ~ "log2_FC > 0"
    ),
    FC_sign = factor(
      .data$FC_sign,
      levels = c("log2_FC < 0", "log2_FC > 0")
    )
  ) %>%
  dplyr::filter(!is.na(.data$FC_sign))

GSE226643_S_sig_thermodynamic <- GSE226643_S_DE_thermodynamic %>%
  dplyr::filter(
    !is.na(.data$padj),
    .data$padj < 0.05,
    is.finite(.data$log2FC_stationary_12h_minus_early_3h),
    .data$log2FC_stationary_12h_minus_early_3h != 0,
    .data$included_in_DE_descriptor_analysis %in% TRUE
  ) %>%
  dplyr::mutate(
    FC_sign = dplyr::case_when(
      .data$log2FC_stationary_12h_minus_early_3h < 0 ~ "log2_FC < 0",
      .data$log2FC_stationary_12h_minus_early_3h > 0 ~ "log2_FC > 0"
    ),
    FC_sign = factor(
      .data$FC_sign,
      levels = c("log2_FC < 0", "log2_FC > 0")
    )
  ) %>%
  dplyr::filter(!is.na(.data$FC_sign))

if (nrow(GSE226643_N_sig_thermodynamic) == 0L ||
    nrow(GSE226643_S_sig_thermodynamic) == 0L) {
  stop(
    "At least one primary contrast has no FDR-significant transcripts with ",
    "finite ΔGfold predictions.",
    call. = FALSE
  )
}

GSE226643_DE_summary <- dplyr::bind_rows(
  GSE226643_N_DE_all %>%
    dplyr::summarise(
      analysis_id = dplyr::first(.data$analysis_id),
      limitation = dplyr::first(.data$limitation),
      N_genes_tested = dplyr::n(),
      N_significant_padj_lt_0_05 = sum(!is.na(.data$padj) & .data$padj < 0.05),
      N_stationary_up = sum(
        !is.na(.data$padj) & .data$padj < 0.05 &
          .data$log2FC_stationary_12h_minus_early_3h > 0,
        na.rm = TRUE
      ),
      N_stationary_down = sum(
        !is.na(.data$padj) & .data$padj < 0.05 &
          .data$log2FC_stationary_12h_minus_early_3h < 0,
        na.rm = TRUE
      )
    ),
  GSE226643_S_DE_all %>%
    dplyr::summarise(
      analysis_id = dplyr::first(.data$analysis_id),
      limitation = dplyr::first(.data$limitation),
      N_genes_tested = dplyr::n(),
      N_significant_padj_lt_0_05 = sum(!is.na(.data$padj) & .data$padj < 0.05),
      N_stationary_up = sum(
        !is.na(.data$padj) & .data$padj < 0.05 &
          .data$log2FC_stationary_12h_minus_early_3h > 0,
        na.rm = TRUE
      ),
      N_stationary_down = sum(
        !is.na(.data$padj) & .data$padj < 0.05 &
          .data$log2FC_stationary_12h_minus_early_3h < 0,
        na.rm = TRUE
      )
    )
)

print(GSE226643_DE_summary)

readr::write_tsv(
  GSE226643_N_DE_all,
  file.path(GSE226643_table_dir, "GSE226643_N_DESeq2_all_gene_results.tsv")
)

readr::write_tsv(
  GSE226643_S_DE_all,
  file.path(GSE226643_table_dir, "GSE226643_S_DESeq2_all_gene_results.tsv")
)

readr::write_tsv(
  GSE226643_N_DE_thermodynamic,
  file.path(GSE226643_table_dir, "GSE226643_N_DESeq2_MFE_results.tsv")
)

readr::write_tsv(
  GSE226643_S_DE_thermodynamic,
  file.path(GSE226643_table_dir, "GSE226643_S_DESeq2_MFE_results.tsv")
)

readr::write_tsv(
  GSE226643_N_sig_thermodynamic,
  file.path(GSE226643_table_dir, "GSE226643_N_significant_DE_MFE_transcripts.tsv")
)

readr::write_tsv(
  GSE226643_S_sig_thermodynamic,
  file.path(GSE226643_table_dir, "GSE226643_S_significant_DE_MFE_transcripts.tsv")
)

readr::write_tsv(
  GSE226643_DE_summary,
  file.path(GSE226643_table_dir, "GSE226643_DESeq2_primary_contrast_summary.tsv")
)

GSE226643_N_vst <- DESeq2::vst(GSE226643_N_dds, blind = TRUE)
GSE226643_S_vst <- DESeq2::vst(GSE226643_S_dds, blind = TRUE)

GSE226643_N_pca_data_raw <- DESeq2::plotPCA(
  GSE226643_N_vst,
  intgroup = "DESeq2_stage",
  returnData = TRUE
)

GSE226643_N_pca_percent_variance <- round(
  100 * attr(GSE226643_N_pca_data_raw, "percentVar")
)

GSE226643_N_pca_data <- GSE226643_N_pca_data_raw %>%
  tibble::as_tibble() %>%
  dplyr::mutate(
    limitation = "Nitrogen limitation",
    DESeq2_stage = factor(
      .data$DESeq2_stage,
      levels = c("Early_3h", "Stationary_12h")
    )
  )

GSE226643_S_pca_data_raw <- DESeq2::plotPCA(
  GSE226643_S_vst,
  intgroup = "DESeq2_stage",
  returnData = TRUE
)

GSE226643_S_pca_percent_variance <- round(
  100 * attr(GSE226643_S_pca_data_raw, "percentVar")
)

GSE226643_S_pca_data <- GSE226643_S_pca_data_raw %>%
  tibble::as_tibble() %>%
  dplyr::mutate(
    limitation = "Sulfur limitation",
    DESeq2_stage = factor(
      .data$DESeq2_stage,
      levels = c("Early_3h", "Stationary_12h")
    )
  )

GSE226643_N_PCA_plot <- ggplot(
  GSE226643_N_pca_data,
  aes(
    x = .data$PC1,
    y = .data$PC2,
    colour = .data$DESeq2_stage
  )
) +
  geom_point(size = 3) +
  geom_text(
    aes(label = .data$name),
    vjust = -0.8,
    size = 3,
    show.legend = FALSE
  ) +
  scale_colour_manual(
    values = c(
      "Early_3h" = "#8cc5e3",
      "Stationary_12h" = "#2066a8"
    ),
    labels = c(
      "Early_3h" = "3 h",
      "Stationary_12h" = "12 h"
    )
  ) +
  theme_classic() +
  labs(
    title = "Nitrogen limitation",
    x = paste0("PC1: ", GSE226643_N_pca_percent_variance[1], "% variance"),
    y = paste0("PC2: ", GSE226643_N_pca_percent_variance[2], "% variance"),
    colour = "Time point"
  )

GSE226643_S_PCA_plot <- ggplot(
  GSE226643_S_pca_data,
  aes(
    x = .data$PC1,
    y = .data$PC2,
    colour = .data$DESeq2_stage
  )
) +
  geom_point(size = 3) +
  geom_text(
    aes(label = .data$name),
    vjust = -0.8,
    size = 3,
    show.legend = FALSE
  ) +
  scale_colour_manual(
    values = c(
      "Early_3h" = "#8cc5e3",
      "Stationary_12h" = "#2066a8"
    ),
    labels = c(
      "Early_3h" = "3 h",
      "Stationary_12h" = "12 h"
    )
  ) +
  theme_classic() +
  labs(
    title = "Sulfur limitation",
    x = paste0("PC1: ", GSE226643_S_pca_percent_variance[1], "% variance"),
    y = paste0("PC2: ", GSE226643_S_pca_percent_variance[2], "% variance"),
    colour = "Time point"
  )

GSE226643_PCA_panel <- GSE226643_N_PCA_plot +
  GSE226643_S_PCA_plot +
  patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

print(GSE226643_PCA_panel)

GSE226643_N_sample_correlation <- cor(
  SummarizedExperiment::assay(GSE226643_N_vst)
)

GSE226643_N_sample_correlation_long <- as.data.frame(
  GSE226643_N_sample_correlation,
  check.names = FALSE
) %>%
  tibble::rownames_to_column("sample_x") %>%
  tidyr::pivot_longer(
    cols = -sample_x,
    names_to = "sample_y",
    values_to = "Pearson_r"
  )

GSE226643_S_sample_correlation <- cor(
  SummarizedExperiment::assay(GSE226643_S_vst)
)

GSE226643_S_sample_correlation_long <- as.data.frame(
  GSE226643_S_sample_correlation,
  check.names = FALSE
) %>%
  tibble::rownames_to_column("sample_x") %>%
  tidyr::pivot_longer(
    cols = -sample_x,
    names_to = "sample_y",
    values_to = "Pearson_r"
  )

GSE226643_N_correlation_plot <- ggplot(
  GSE226643_N_sample_correlation_long,
  aes(
    x = .data$sample_x,
    y = .data$sample_y,
    fill = .data$Pearson_r
  )
) +
  geom_tile(colour = "white", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%.2f", .data$Pearson_r)), size = 3) +
  scale_fill_gradient(limits = c(0, 1), low = "white", high = "#2066a8") +
  coord_equal() +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.ticks = element_blank()
  ) +
  labs(title = "Nitrogen limitation: VST sample correlation", fill = "Pearson r")

GSE226643_S_correlation_plot <- ggplot(
  GSE226643_S_sample_correlation_long,
  aes(
    x = .data$sample_x,
    y = .data$sample_y,
    fill = .data$Pearson_r
  )
) +
  geom_tile(colour = "white", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%.2f", .data$Pearson_r)), size = 3) +
  scale_fill_gradient(limits = c(0, 1), low = "white", high = "#2066a8") +
  coord_equal() +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.ticks = element_blank()
  ) +
  labs(title = "Sulfur limitation: VST sample correlation", fill = "Pearson r")

GSE226643_correlation_panel <- GSE226643_N_correlation_plot +
  GSE226643_S_correlation_plot +
  patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

print(GSE226643_correlation_panel)

GSE226643_N_MA_data <- GSE226643_N_DE_all %>%
  dplyr::filter(
    is.finite(.data$base_mean),
    .data$base_mean > 0,
    is.finite(.data$log2FC_stationary_12h_minus_early_3h)
  ) %>%
  dplyr::mutate(
    significance = dplyr::if_else(
      !is.na(.data$padj) & .data$padj < 0.05,
      "FDR < 0.05",
      "Not significant"
    )
  )

GSE226643_S_MA_data <- GSE226643_S_DE_all %>%
  dplyr::filter(
    is.finite(.data$base_mean),
    .data$base_mean > 0,
    is.finite(.data$log2FC_stationary_12h_minus_early_3h)
  ) %>%
  dplyr::mutate(
    significance = dplyr::if_else(
      !is.na(.data$padj) & .data$padj < 0.05,
      "FDR < 0.05",
      "Not significant"
    )
  )

GSE226643_N_MA_plot <- ggplot(
  GSE226643_N_MA_data,
  aes(
    x = .data$base_mean,
    y = .data$log2FC_stationary_12h_minus_early_3h,
    colour = .data$significance
  )
) +
  geom_hline(yintercept = 0, linewidth = 0.25) +
  geom_point(alpha = 0.45, size = 0.75) +
  scale_x_log10() +
  scale_colour_manual(values = c("Not significant" = "grey70", "FDR < 0.05" = "#b2182b")) +
  theme_classic() +
  labs(
    title = "Nitrogen limitation",
    x = "DESeq2 base mean (log10 scale)",
    y = expression(log[2]~"FC (12 h / 3 h)"),
    colour = ""
  )

GSE226643_S_MA_plot <- ggplot(
  GSE226643_S_MA_data,
  aes(
    x = .data$base_mean,
    y = .data$log2FC_stationary_12h_minus_early_3h,
    colour = .data$significance
  )
) +
  geom_hline(yintercept = 0, linewidth = 0.25) +
  geom_point(alpha = 0.45, size = 0.75) +
  scale_x_log10() +
  scale_colour_manual(values = c("Not significant" = "grey70", "FDR < 0.05" = "#b2182b")) +
  theme_classic() +
  labs(
    title = "Sulfur limitation",
    x = "DESeq2 base mean (log10 scale)",
    y = expression(log[2]~"FC (12 h / 3 h)"),
    colour = ""
  )

GSE226643_MA_panel <- GSE226643_N_MA_plot +
  GSE226643_S_MA_plot +
  patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

print(GSE226643_MA_panel)

readr::write_tsv(
  GSE226643_N_pca_data,
  file.path(GSE226643_qc_dir, "GSE226643_N_vst_PCA_coordinates.tsv")
)

readr::write_tsv(
  GSE226643_S_pca_data,
  file.path(GSE226643_qc_dir, "GSE226643_S_vst_PCA_coordinates.tsv")
)

readr::write_tsv(
  GSE226643_N_sample_correlation_long,
  file.path(GSE226643_qc_dir, "GSE226643_N_vst_sample_correlations.tsv")
)

readr::write_tsv(
  GSE226643_S_sample_correlation_long,
  file.path(GSE226643_qc_dir, "GSE226643_S_vst_sample_correlations.tsv")
)

GSE226643_N_MFE_down <- GSE226643_N_sig_thermodynamic %>%
  dplyr::filter(.data$FC_sign == "log2_FC < 0") %>%
  dplyr::pull(.data$Normalized_MFE)

GSE226643_N_MFE_up <- GSE226643_N_sig_thermodynamic %>%
  dplyr::filter(.data$FC_sign == "log2_FC > 0") %>%
  dplyr::pull(.data$Normalized_MFE)

GSE226643_S_MFE_down <- GSE226643_S_sig_thermodynamic %>%
  dplyr::filter(.data$FC_sign == "log2_FC < 0") %>%
  dplyr::pull(.data$Normalized_MFE)

GSE226643_S_MFE_up <- GSE226643_S_sig_thermodynamic %>%
  dplyr::filter(.data$FC_sign == "log2_FC > 0") %>%
  dplyr::pull(.data$Normalized_MFE)

if (min(
  length(GSE226643_N_MFE_down),
  length(GSE226643_N_MFE_up),
  length(GSE226643_S_MFE_down),
  length(GSE226643_S_MFE_up)
) < 2L) {
  stop(
    "Each primary comparison must retain at least two significant transcripts ",
    "in both fold-change directions for the ΔGfold distribution plots.",
    call. = FALSE
  )
}

GSE226643_N_MFE_Wilcoxon_p <- stats::wilcox.test(
  GSE226643_N_MFE_up,
  GSE226643_N_MFE_down,
  exact = FALSE
)$p.value

GSE226643_S_MFE_Wilcoxon_p <- stats::wilcox.test(
  GSE226643_S_MFE_up,
  GSE226643_S_MFE_down,
  exact = FALSE
)$p.value

GSE226643_N_MFE_delta <- mean(GSE226643_N_MFE_up) - mean(GSE226643_N_MFE_down)
GSE226643_S_MFE_delta <- mean(GSE226643_S_MFE_up) - mean(GSE226643_S_MFE_down)

GSE226643_significant_distribution_statistics <- dplyr::bind_rows(
  tibble::tibble(
    analysis_id = "GSE226643_N_12h_vs_3h",
    limitation = "Nitrogen limitation",
    descriptor = "Normalized_MFE",
    N_12h_lower = length(GSE226643_N_MFE_down),
    N_12h_higher = length(GSE226643_N_MFE_up),
    mean_12h_lower = mean(GSE226643_N_MFE_down),
    mean_12h_higher = mean(GSE226643_N_MFE_up),
    delta_12h_higher_minus_lower = mean(GSE226643_N_MFE_up) - mean(GSE226643_N_MFE_down),
    Wilcoxon_p = GSE226643_N_MFE_Wilcoxon_p
  ),
  tibble::tibble(
    analysis_id = "GSE226643_S_12h_vs_3h",
    limitation = "Sulfur limitation",
    descriptor = "Normalized_MFE",
    N_12h_lower = length(GSE226643_S_MFE_down),
    N_12h_higher = length(GSE226643_S_MFE_up),
    mean_12h_lower = mean(GSE226643_S_MFE_down),
    mean_12h_higher = mean(GSE226643_S_MFE_up),
    delta_12h_higher_minus_lower = mean(GSE226643_S_MFE_up) - mean(GSE226643_S_MFE_down),
    Wilcoxon_p = GSE226643_S_MFE_Wilcoxon_p
  )
) %>%
  dplyr::mutate(
    Wilcoxon_p_BH = stats::p.adjust(
      .data$Wilcoxon_p,
      method = "BH",
      n = 4L
    )
  )

print(GSE226643_significant_distribution_statistics)

readr::write_tsv(
  GSE226643_significant_distribution_statistics,
  file.path(
    GSE226643_table_dir,
    "GSE226643_significant_DE_MFE_directional_statistics.tsv"
  )
)

GSE226643_MFE_y_limits <- range(
  c(
    GSE226643_N_sig_thermodynamic$Normalized_MFE,
    GSE226643_S_sig_thermodynamic$Normalized_MFE
  ),
  finite = TRUE
)

GSE226643_N_MFE_plot <- ggplot(
  GSE226643_N_sig_thermodynamic,
  aes(x = .data$FC_sign, y = .data$Normalized_MFE, fill = .data$FC_sign)
) +
  geom_violin(trim = FALSE, alpha = 0.75, colour = "black", linewidth = 0.2) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", linewidth = 0.25) +
  coord_cartesian(ylim = GSE226643_MFE_y_limits) +
  scale_fill_manual(values = c("log2_FC < 0" = "#8cc5e3", "log2_FC > 0" = "#2066a8")) +
  scale_x_discrete(labels = c("log2_FC < 0" = "Lower at 12 h", "log2_FC > 0" = "Higher at 12 h")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1, colour = "black")
  ) +
  labs(
    title = paste0(
      "Nitrogen limitation\nΔMFE = ",
      format(round(GSE226643_N_MFE_delta, 3), nsmall = 3),
      " (Wilcoxon p = ",
      formatC(GSE226643_N_MFE_Wilcoxon_p, format = "e", digits = 2),
      ")"
    ),
    y = "Normalized ΔGfold (ΔGfold / nt)"
  )

GSE226643_S_MFE_plot <- ggplot(
  GSE226643_S_sig_thermodynamic,
  aes(x = .data$FC_sign, y = .data$Normalized_MFE, fill = .data$FC_sign)
) +
  geom_violin(trim = FALSE, alpha = 0.75, colour = "black", linewidth = 0.2) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", linewidth = 0.25) +
  coord_cartesian(ylim = GSE226643_MFE_y_limits) +
  scale_fill_manual(values = c("log2_FC < 0" = "#8cc5e3", "log2_FC > 0" = "#2066a8")) +
  scale_x_discrete(labels = c("log2_FC < 0" = "Lower at 12 h", "log2_FC > 0" = "Higher at 12 h")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1, colour = "black")
  ) +
  labs(
    title = paste0(
      "Sulfur limitation\nΔMFE = ",
      format(round(GSE226643_S_MFE_delta, 3), nsmall = 3),
      " (Wilcoxon p = ",
      formatC(GSE226643_S_MFE_Wilcoxon_p, format = "e", digits = 2),
      ")"
    ),
    y = "Normalized ΔGfold (ΔGfold / nt)"
  )

GSE226643_MFE_distribution_panel <- GSE226643_N_MFE_plot + GSE226643_S_MFE_plot

print(GSE226643_MFE_distribution_panel)

GSE226643_N_MFE_quartiles <- stats::quantile(
  GSE226643_N_sig_thermodynamic$Normalized_MFE,
  probs = c(0.25, 0.75),
  na.rm = TRUE,
  names = FALSE
)

GSE226643_S_MFE_quartiles <- stats::quantile(
  GSE226643_S_sig_thermodynamic$Normalized_MFE,
  probs = c(0.25, 0.75),
  na.rm = TRUE,
  names = FALSE
)

GSE226643_N_sig_MFE_abundance_groups <- GSE226643_N_sig_thermodynamic %>%
  dplyr::mutate(
    MFE_group = dplyr::case_when(
      .data$Normalized_MFE <= GSE226643_N_MFE_quartiles[1] ~ "Low_MFE",
      .data$Normalized_MFE >= GSE226643_N_MFE_quartiles[2] ~ "High_MFE",
      TRUE ~ "Intermediate_MFE"
    )
  )

GSE226643_S_sig_MFE_abundance_groups <- GSE226643_S_sig_thermodynamic %>%
  dplyr::mutate(
    MFE_group = dplyr::case_when(
      .data$Normalized_MFE <= GSE226643_S_MFE_quartiles[1] ~ "Low_MFE",
      .data$Normalized_MFE >= GSE226643_S_MFE_quartiles[2] ~ "High_MFE",
      TRUE ~ "Intermediate_MFE"
    )
  )

GSE226643_N_significant_DE_abundance_fraction <- GSE226643_N_sig_MFE_abundance_groups %>%
  dplyr::summarise(
    analysis_id = dplyr::first(.data$analysis_id),
    limitation = dplyr::first(.data$limitation),
    comparison_label = dplyr::first(.data$comparison_label),
    Control = "3 h early exponential",
    Treated = "12 h stationary",
    N_significant_transcripts = dplyr::n(),
    N_low_MFE = sum(.data$MFE_group == "Low_MFE"),
    N_high_MFE = sum(.data$MFE_group == "High_MFE"),
    low_MFE_quartile_cutoff = GSE226643_N_MFE_quartiles[1],
    high_MFE_quartile_cutoff = GSE226643_N_MFE_quartiles[2],
    Control_low_fraction = sum(
      .data$early_3h_mean_normalised_counts[.data$MFE_group == "Low_MFE"]
    ) / sum(.data$early_3h_mean_normalised_counts),
    Treated_low_fraction = sum(
      .data$stationary_12h_mean_normalised_counts[.data$MFE_group == "Low_MFE"]
    ) / sum(.data$stationary_12h_mean_normalised_counts),
    Control_high_fraction = sum(
      .data$early_3h_mean_normalised_counts[.data$MFE_group == "High_MFE"]
    ) / sum(.data$early_3h_mean_normalised_counts),
    Treated_high_fraction = sum(
      .data$stationary_12h_mean_normalised_counts[.data$MFE_group == "High_MFE"]
    ) / sum(.data$stationary_12h_mean_normalised_counts)
  )

GSE226643_S_significant_DE_abundance_fraction <- GSE226643_S_sig_MFE_abundance_groups %>%
  dplyr::summarise(
    analysis_id = dplyr::first(.data$analysis_id),
    limitation = dplyr::first(.data$limitation),
    comparison_label = dplyr::first(.data$comparison_label),
    Control = "3 h early exponential",
    Treated = "12 h stationary",
    N_significant_transcripts = dplyr::n(),
    N_low_MFE = sum(.data$MFE_group == "Low_MFE"),
    N_high_MFE = sum(.data$MFE_group == "High_MFE"),
    low_MFE_quartile_cutoff = GSE226643_S_MFE_quartiles[1],
    high_MFE_quartile_cutoff = GSE226643_S_MFE_quartiles[2],
    Control_low_fraction = sum(
      .data$early_3h_mean_normalised_counts[.data$MFE_group == "Low_MFE"]
    ) / sum(.data$early_3h_mean_normalised_counts),
    Treated_low_fraction = sum(
      .data$stationary_12h_mean_normalised_counts[.data$MFE_group == "Low_MFE"]
    ) / sum(.data$stationary_12h_mean_normalised_counts),
    Control_high_fraction = sum(
      .data$early_3h_mean_normalised_counts[.data$MFE_group == "High_MFE"]
    ) / sum(.data$early_3h_mean_normalised_counts),
    Treated_high_fraction = sum(
      .data$stationary_12h_mean_normalised_counts[.data$MFE_group == "High_MFE"]
    ) / sum(.data$stationary_12h_mean_normalised_counts)
  )

GSE226643_significant_DE_abundance_fraction <- dplyr::bind_rows(
  GSE226643_N_significant_DE_abundance_fraction,
  GSE226643_S_significant_DE_abundance_fraction
) %>%
  dplyr::mutate(
    comparison_label = factor(
      .data$comparison_label,
      levels = c(
        "N limitation: 12 h stationary minus 3 h early exponential",
        "S limitation: 12 h stationary minus 3 h early exponential"
      )
    )
  ) %>%
  dplyr::arrange(.data$comparison_label)

print(GSE226643_significant_DE_abundance_fraction)

readr::write_tsv(
  GSE226643_significant_DE_abundance_fraction,
  file.path(
    GSE226643_table_dir,
    "GSE226643_significant_DE_MFE_abundance_fractions.tsv"
  )
)

GSE226643_high_MFE_fraction_plot_data <- GSE226643_significant_DE_abundance_fraction %>%
  dplyr::select(
    .data$comparison_label,
    .data$Control,
    .data$Treated,
    .data$Control_high_fraction,
    .data$Treated_high_fraction
  ) %>%
  tidyr::pivot_longer(
    cols = c(.data$Control_high_fraction, .data$Treated_high_fraction),
    names_to = "Condition_type",
    values_to = "Fraction"
  ) %>%
  dplyr::mutate(
    Sample_group = dplyr::if_else(
      .data$Condition_type == "Control_high_fraction",
      "Control",
      "Treated"
    ),
    Sample_group = factor(.data$Sample_group, levels = c("Control", "Treated"))
  )

GSE226643_low_MFE_fraction_plot_data <- GSE226643_significant_DE_abundance_fraction %>%
  dplyr::select(
    .data$comparison_label,
    .data$Control,
    .data$Treated,
    .data$Control_low_fraction,
    .data$Treated_low_fraction
  ) %>%
  tidyr::pivot_longer(
    cols = c(.data$Control_low_fraction, .data$Treated_low_fraction),
    names_to = "Condition_type",
    values_to = "Fraction"
  ) %>%
  dplyr::mutate(
    Sample_group = dplyr::if_else(
      .data$Condition_type == "Control_low_fraction",
      "Control",
      "Treated"
    ),
    Sample_group = factor(.data$Sample_group, levels = c("Control", "Treated"))
  )

GSE226643_high_MFE_fraction_plot <- ggplot(
  GSE226643_high_MFE_fraction_plot_data,
  aes(x = .data$comparison_label, y = .data$Fraction, fill = .data$Sample_group)
) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7,
    colour = "black",
    linewidth = 0.2
  ) +
  geom_text(
    aes(label = sprintf("%.3f", .data$Fraction)),
    position = position_dodge(width = 0.8),
    vjust = -0.35,
    size = 3
  ) +
  coord_cartesian(
    ylim = c(
      0,
      max(GSE226643_high_MFE_fraction_plot_data$Fraction, na.rm = TRUE) * 1.18
    )
  ) +
  scale_fill_manual(
    values = c("Control" = "#8cc5e3", "Treated" = "#2066a8"),
    labels = c("Control" = "3 h early exponential", "Treated" = "12 h stationary")
  ) +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 20, hjust = 1, colour = "black"),
    legend.position = "top"
  ) +
  labs(
    title = "High-ΔGfold transcript abundance fraction",
    y = "Fraction of significant-DE transcript abundance",
    fill = ""
  )

GSE226643_low_MFE_fraction_plot <- ggplot(
  GSE226643_low_MFE_fraction_plot_data,
  aes(x = .data$comparison_label, y = .data$Fraction, fill = .data$Sample_group)
) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7,
    colour = "black",
    linewidth = 0.2
  ) +
  geom_text(
    aes(label = sprintf("%.3f", .data$Fraction)),
    position = position_dodge(width = 0.8),
    vjust = -0.35,
    size = 3
  ) +
  coord_cartesian(
    ylim = c(
      0,
      max(GSE226643_low_MFE_fraction_plot_data$Fraction, na.rm = TRUE) * 1.18
    )
  ) +
  scale_fill_manual(
    values = c("Control" = "#8cc5e3", "Treated" = "#2066a8"),
    labels = c("Control" = "3 h early exponential", "Treated" = "12 h stationary")
  ) +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 20, hjust = 1, colour = "black"),
    legend.position = "top"
  ) +
  labs(
    title = "Low-ΔGfold transcript abundance fraction",
    y = "Fraction of significant-DE transcript abundance",
    fill = ""
  )

print(GSE226643_high_MFE_fraction_plot)
print(GSE226643_low_MFE_fraction_plot)

GSE226643_GO_terms <- c(
  "GO:0006260", "GO:0051301", "GO:0007049", "GO:0042254", "GO:0006412",
  "GO:0009150", "GO:0009219", "GO:0008652", "GO:0003700", "GO:0017000",
  "GO:0009432", "GO:0006950", "GO:0006270", "GO:0071555", "GO:0015986",
  "GO:0009060", "GO:0006096", "GO:0006413", "GO:0046677", "GO:0015968",
  "GO:0043038", "GO:0009252", "GO:0006275", "GO:0006265", "GO:0044038",
  "GO:0006414", "GO:0006450", "GO:0007059", "GO:0000910", "GO:0006094",
  "GO:0006099", "GO:0006119", "GO:0006113", "GO:0046034", "GO:0006633",
  "GO:0006635", "GO:0004812", "GO:1990904", "GO:0009103", "GO:0009245",
  "GO:0009234", "GO:0008614", "GO:0071897", "GO:0006754", "GO:0006281",
  "GO:0006974", "GO:0000719", "GO:0006284", "GO:0006285", "GO:0006287",
  "GO:0006289", "GO:0006296", "GO:0006297", "GO:0006298", "GO:0000725",
  "GO:0006302", "GO:0000724", "GO:0006303", "GO:0006307", "GO:0036297",
  "GO:0019985", "GO:0042276", "GO:0070987", "GO:0000731"
)

GSE226643_GO_names <- c(
  "DNA replication", "Cell division", "Cell cycle", "Ribosome biogenesis", "Translation",
  "Purine ribonucleotide biosynthesis", "Pyrimidine ribonucleotide biosynthesis", "Amino acid biosynthesis", "DNA-binding TF activity", "Antibiotic metabolic process",
  "SOS response", "Response to stress", "DNA replication initiation", "Cell wall organization", "ATP synthesis coupled proton transport",
  "Aerobic respiration", "Glycolysis", "Translation initiation", "Response to antibiotic", "Stringent response",
  "tRNA aminoacylation", "Peptidoglycan biosynthesis", "Regulation of DNA replication", "DNA topology", "Cell wall macromolecule biosynthetic process",
  "Translation elongation", "Regulation of translation", "Chromosome segregation", "Cytokinesis", "Gluconeogenesis",
  "TCA cycle", "Oxidative phosphorylation", "Fermentation", "ATP metabolic process", "Fatty acid biosynthesis",
  "Fatty acid beta-oxidation", "Aminoacyl-tRNA ligase activity", "Ribonucleoprotein complex", "Lipopolysaccharide biosynthesis", "Lipid A biosynthesis",
  "Menaquinone biosynthesis", "Pyridoxine biosynthesis", "DNA biosynthetic process", "ATP biosynthetic process", "DNA repair",
  "DNA damage response", "Photoreactive repair", "Base-excision repair", "BER: AP site formation", "BER: gap-filling",
  "Nucleotide-excision repair", "NER: DNA incision", "NER: DNA gap filling", "Mismatch repair", "Recombinational repair",
  "Double-strand break repair", "DSB repair via homologous recombination", "DSB repair via nonhomologous end joining", "DNA alkylation repair", "Interstrand cross-link repair",
  "Translesion synthesis", "Error-prone translesion synthesis", "Error-free translesion synthesis", "DNA synthesis involved in DNA repair"
)

if (length(GSE226643_GO_terms) != length(GSE226643_GO_names)) {
  stop("The curated GO-ID and GO-name vectors have different lengths.", call. = FALSE)
}

GSE226643_GO_panel <- tibble::tibble(
  GO = GSE226643_GO_terms,
  term = GSE226643_GO_names
)

GSE226643_TERM2GENE_cache <- file.path(
  GSE226643_gsea_dir,
  "GSE226643_GO_panel_TERM2GENE_E_coli_K12.tsv"
)

if (file.exists(GSE226643_TERM2GENE_cache) && !GSE226643_refresh_GO_cache) {
  GSE226643_GO_genes_all <- readr::read_tsv(
    GSE226643_TERM2GENE_cache,
    show_col_types = FALSE,
    progress = FALSE
  )
} else {

  GSE226643_GO_gene_tables <- vector("list", length(GSE226643_GO_terms))

  for (GSE226643_GO_index in seq_along(GSE226643_GO_terms)) {
    GSE226643_one_GO <- GSE226643_GO_terms[GSE226643_GO_index]
    GSE226643_GO_url <- paste0(
      "https://api.geneontology.org/api/bioentity/function/",
      GSE226643_one_GO,
      "/genes?taxon=NCBITaxon:83333"
    )

    GSE226643_GO_response <- tryCatch(
      httr::GET(GSE226643_GO_url),
      error = function(error_message) NULL
    )

    if (is.null(GSE226643_GO_response) ||
        httr::status_code(GSE226643_GO_response) != 200L) {
      warning("No GO response retrieved for ", GSE226643_one_GO, call. = FALSE)
      GSE226643_GO_gene_tables[[GSE226643_GO_index]] <- tibble::tibble(
        SYMBOL = character(0),
        GO = GSE226643_one_GO
      )
      next
    }

    GSE226643_GO_json <- tryCatch(
      jsonlite::fromJSON(
        httr::content(
          GSE226643_GO_response,
          as = "text",
          encoding = "UTF-8"
        )
      ),
      error = function(error_message) NULL
    )

    if (is.null(GSE226643_GO_json) ||
        is.null(GSE226643_GO_json$associations) ||
        length(GSE226643_GO_json$associations) == 0L) {
      GSE226643_GO_gene_tables[[GSE226643_GO_index]] <- tibble::tibble(
        SYMBOL = character(0),
        GO = GSE226643_one_GO
      )
      next
    }

    GSE226643_GO_symbols <- unique(
      as.character(GSE226643_GO_json$associations$subject$label)
    )

    GSE226643_GO_gene_tables[[GSE226643_GO_index]] <- tibble::tibble(
      SYMBOL = GSE226643_GO_symbols,
      GO = GSE226643_one_GO
    )
  }

  GSE226643_GO_genes_all <- dplyr::bind_rows(GSE226643_GO_gene_tables) %>%
    dplyr::left_join(GSE226643_GO_panel, by = "GO") %>%
    dplyr::mutate(
      SYMBOL = stringr::str_trim(as.character(.data$SYMBOL)),
      SYMBOL = stringr::str_replace(.data$SYMBOL, "_ECOLI$", ""),
      SYMBOL = toupper(.data$SYMBOL)
    ) %>%
    dplyr::filter(!is.na(.data$SYMBOL), .data$SYMBOL != "") %>%
    dplyr::distinct(.data$term, .data$SYMBOL, .keep_all = TRUE) %>%
    dplyr::arrange(.data$term, .data$SYMBOL)

  readr::write_tsv(
    GSE226643_GO_genes_all,
    GSE226643_TERM2GENE_cache
  )
}

GSE226643_TERM2GENE <- GSE226643_GO_genes_all %>%
  dplyr::transmute(
    term = as.character(.data$term),
    gene = toupper(trimws(as.character(.data$SYMBOL)))
  ) %>%
  dplyr::filter(!is.na(.data$term), .data$term != "", !is.na(.data$gene), .data$gene != "") %>%
  dplyr::distinct()

if (nrow(GSE226643_TERM2GENE) == 0L) {
  stop(
    "The cached/retrieved GO panel has no E. coli K-12 gene annotations.",
    call. = FALSE
  )
}

print(
  GSE226643_TERM2GENE %>%
    dplyr::count(.data$term, name = "N_genes") %>%
    dplyr::arrange(.data$term)
)

GSE226643_N_GSEA_rank_table <- GSE226643_N_DE_all %>%
  dplyr::left_join(
    GSE226643_gene_annotation %>%
      dplyr::select(.data$Gene_id, .data$gene_symbol),
    by = "Gene_id"
  ) %>%
  dplyr::mutate(gene_symbol = toupper(trimws(as.character(.data$gene_symbol)))) %>%
  dplyr::filter(
    is.finite(.data$log2FC_stationary_12h_minus_early_3h),
    !is.na(.data$gene_symbol),
    .data$gene_symbol != ""
  ) %>%
  dplyr::group_by(.data$gene_symbol) %>%
  dplyr::slice_max(
    order_by = abs(.data$log2FC_stationary_12h_minus_early_3h),
    n = 1,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(dplyr::desc(.data$log2FC_stationary_12h_minus_early_3h))

GSE226643_S_GSEA_rank_table <- GSE226643_S_DE_all %>%
  dplyr::left_join(
    GSE226643_gene_annotation %>%
      dplyr::select(.data$Gene_id, .data$gene_symbol),
    by = "Gene_id"
  ) %>%
  dplyr::mutate(gene_symbol = toupper(trimws(as.character(.data$gene_symbol)))) %>%
  dplyr::filter(
    is.finite(.data$log2FC_stationary_12h_minus_early_3h),
    !is.na(.data$gene_symbol),
    .data$gene_symbol != ""
  ) %>%
  dplyr::group_by(.data$gene_symbol) %>%
  dplyr::slice_max(
    order_by = abs(.data$log2FC_stationary_12h_minus_early_3h),
    n = 1,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(dplyr::desc(.data$log2FC_stationary_12h_minus_early_3h))

GSE226643_N_GSEA_ranks <- GSE226643_N_GSEA_rank_table$log2FC_stationary_12h_minus_early_3h
names(GSE226643_N_GSEA_ranks) <- GSE226643_N_GSEA_rank_table$gene_symbol
GSE226643_N_GSEA_ranks <- sort(GSE226643_N_GSEA_ranks, decreasing = TRUE)

GSE226643_S_GSEA_ranks <- GSE226643_S_GSEA_rank_table$log2FC_stationary_12h_minus_early_3h
names(GSE226643_S_GSEA_ranks) <- GSE226643_S_GSEA_rank_table$gene_symbol
GSE226643_S_GSEA_ranks <- sort(GSE226643_S_GSEA_ranks, decreasing = TRUE)

if (length(GSE226643_N_GSEA_ranks) < 30L || length(GSE226643_S_GSEA_ranks) < 30L) {
  stop("Too few ranked gene symbols are available for GSEA.", call. = FALSE)
}

GSE226643_N_GSEA_object <- tryCatch(
  clusterProfiler::GSEA(
    geneList = GSE226643_N_GSEA_ranks,
    TERM2GENE = GSE226643_TERM2GENE,
    pvalueCutoff = 0.10,
    pAdjustMethod = "BH",
    minGSSize = 5,
    maxGSSize = 1000,
    verbose = FALSE
  ),
  error = function(error_message) NULL
)

GSE226643_S_GSEA_object <- tryCatch(
  clusterProfiler::GSEA(
    geneList = GSE226643_S_GSEA_ranks,
    TERM2GENE = GSE226643_TERM2GENE,
    pvalueCutoff = 0.10,
    pAdjustMethod = "BH",
    minGSSize = 5,
    maxGSSize = 1000,
    verbose = FALSE
  ),
  error = function(error_message) NULL
)

GSE226643_N_GSEA_table <- if (is.null(GSE226643_N_GSEA_object)) {
  data.frame()
} else {
  as.data.frame(GSE226643_N_GSEA_object)
}

GSE226643_S_GSEA_table <- if (is.null(GSE226643_S_GSEA_object)) {
  data.frame()
} else {
  as.data.frame(GSE226643_S_GSEA_object)
}

GSE226643_N_GSEA_results <- if (nrow(GSE226643_N_GSEA_table) == 0L) {
  tibble::tibble(
    analysis_id = character(0),
    limitation = character(0),
    term = character(0),
    NES = numeric(0),
    p.adjust = numeric(0)
  )
} else {
  tibble::as_tibble(GSE226643_N_GSEA_table) %>%
    dplyr::mutate(
      analysis_id = "GSE226643_N_12h_vs_3h",
      limitation = "Nitrogen limitation",
      term = as.character(.data$ID)
    )
}

GSE226643_S_GSEA_results <- if (nrow(GSE226643_S_GSEA_table) == 0L) {
  tibble::tibble(
    analysis_id = character(0),
    limitation = character(0),
    term = character(0),
    NES = numeric(0),
    p.adjust = numeric(0)
  )
} else {
  tibble::as_tibble(GSE226643_S_GSEA_table) %>%
    dplyr::mutate(
      analysis_id = "GSE226643_S_12h_vs_3h",
      limitation = "Sulfur limitation",
      term = as.character(.data$ID)
    )
}

GSE226643_GSEA_results <- dplyr::bind_rows(
  GSE226643_N_GSEA_results,
  GSE226643_S_GSEA_results
)

readr::write_tsv(
  GSE226643_N_GSEA_rank_table,
  file.path(GSE226643_gsea_dir, "GSE226643_N_GSEA_ranked_genes.tsv")
)

readr::write_tsv(
  GSE226643_S_GSEA_rank_table,
  file.path(GSE226643_gsea_dir, "GSE226643_S_GSEA_ranked_genes.tsv")
)

readr::write_tsv(
  GSE226643_GSEA_results,
  file.path(GSE226643_gsea_dir, "GSE226643_N_S_curated_GO_GSEA_results.tsv")
)

GSE226643_GSEA_bubble_data <- GSE226643_GSEA_results %>%
  dplyr::transmute(
    limitation = .data$limitation,
    term = .data$term,
    NES = as.numeric(.data$NES),
    p_adjust = as.numeric(.data$p.adjust)
  ) %>%
  tidyr::complete(
    limitation = c("Nitrogen limitation", "Sulfur limitation"),
    term = GSE226643_GO_names,
    fill = list(NES = 0, p_adjust = 1)
  ) %>%
  dplyr::mutate(
    negative_log10_BH = -log10(pmax(.data$p_adjust, .Machine$double.xmin)),
    term = factor(.data$term, levels = rev(GSE226643_GO_names)),
    limitation = factor(
      .data$limitation,
      levels = c("Nitrogen limitation", "Sulfur limitation")
    )
  )

GSE226643_GSEA_bubble_plot <- ggplot(
  GSE226643_GSEA_bubble_data,
  aes(
    x = .data$limitation,
    y = .data$term,
    size = .data$negative_log10_BH,
    colour = .data$NES
  )
) +
  geom_point(alpha = 0.9) +
  scale_colour_gradient2(
    low = "#2166ac",
    mid = "white",
    high = "#b2182b",
    midpoint = 0,
    name = "NES"
  ) +
  scale_size_continuous(name = expression(-log[10]~"BH-adjusted p")) +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 20, hjust = 1, colour = "black"),
    axis.text.y = element_text(size = 7, colour = "black")
  ) +
  labs(title = "Curated GO-panel GSEA: 12 h stationary versus 3 h")

print(GSE226643_GSEA_bubble_plot)

GSE226643_counts_for_weighted_analysis <- GSE226643_raw_count_matrix[
  GSE226643_weighted_annotation$Gene_id,
  GSE226643_primary_manifest$raw_count_column,
  drop = FALSE
]

if (!identical(
  rownames(GSE226643_counts_for_weighted_analysis),
  GSE226643_weighted_annotation$Gene_id
)) {
  stop(
    "The full-transcriptome weighted count matrix is not aligned to the ",
    "validated annotation table.",
    call. = FALSE
  )
}

GSE226643_length_adjusted_counts <- sweep(
  GSE226643_counts_for_weighted_analysis,
  MARGIN = 1L,
  STATS = GSE226643_weighted_annotation$gene_length_nt,
  FUN = "/"
)

GSE226643_length_adjusted_library_totals <- colSums(
  GSE226643_length_adjusted_counts
)

if (any(!is.finite(GSE226643_length_adjusted_library_totals)) ||
    any(GSE226643_length_adjusted_library_totals <= 0)) {
  stop(
    "At least one sample has an invalid length-adjusted library total.",
    call. = FALSE
  )
}

GSE226643_relative_RNA_abundance <- sweep(
  GSE226643_length_adjusted_counts,
  MARGIN = 2L,
  STATS = GSE226643_length_adjusted_library_totals,
  FUN = "/"
)

if (any(!is.finite(GSE226643_relative_RNA_abundance)) ||
    any(abs(colSums(GSE226643_relative_RNA_abundance) - 1) > 1e-10)) {
  stop(
    "The full-transcriptome relative-RNA-abundance matrix is invalid.",
    call. = FALSE
  )
}

GSE226643_weighted_sample_metadata <- GSE226643_primary_manifest %>%
  dplyr::mutate(
    total_raw_counts = as.numeric(
      colSums(GSE226643_raw_count_matrix)[.data$raw_count_column]
    ),
    length_adjusted_library_total = as.numeric(
      GSE226643_length_adjusted_library_totals[.data$raw_count_column]
    )
  ) %>%
  dplyr::arrange(.data$analysis_id, .data$time_h, .data$author_replicate_id)

GSE226643_abundance_long <- tibble::as_tibble(
  as.data.frame(GSE226643_relative_RNA_abundance),
  rownames = "Gene_id"
) %>%
  tidyr::pivot_longer(
    cols = -Gene_id,
    names_to = "raw_count_column",
    values_to = "relative_RNA_molecule_abundance"
  ) %>%
  dplyr::left_join(
    GSE226643_weighted_annotation %>%
      dplyr::select(
        .data$Gene_id,
        .data$gene_symbol,
        .data$gene_length_nt,
        .data$Normalized_MFE
      ),
    by = "Gene_id"
  ) %>%
  dplyr::left_join(
    GSE226643_weighted_sample_metadata,
    by = "raw_count_column"
  )

GSE226643_weighted_sample_gene_audit <- GSE226643_abundance_long %>%
  dplyr::group_by(
    .data$raw_count_column,
    .data$analysis_id,
    .data$author_condition,
    .data$DESeq2_stage
  ) %>%
  dplyr::summarise(
    N_full_transcriptome_genes = dplyr::n(),
    relative_abundance_total = sum(.data$relative_RNA_molecule_abundance),
    .groups = "drop"
  )

if (any(GSE226643_weighted_sample_gene_audit$N_full_transcriptome_genes !=
        nrow(GSE226643_weighted_annotation)) ||
    any(abs(GSE226643_weighted_sample_gene_audit$relative_abundance_total - 1) > 1e-10)) {
  stop(
    "The full-transcriptome weighted-abundance audit failed.",
    call. = FALSE
  )
}

print(GSE226643_weighted_sample_gene_audit)

readr::write_tsv(
  GSE226643_weighted_annotation,
  file.path(
    GSE226643_table_dir,
    "GSE226643_full_transcriptome_weighted_analysis_annotation.tsv"
  )
)

readr::write_tsv(
  GSE226643_weighted_sample_metadata,
  file.path(
    GSE226643_table_dir,
    "GSE226643_full_transcriptome_weighted_analysis_sample_metadata.tsv"
  )
)

readr::write_tsv(
  GSE226643_weighted_sample_gene_audit,
  file.path(
    GSE226643_table_dir,
    "GSE226643_full_transcriptome_weighted_analysis_sample_gene_audit.tsv"
  )
)

readr::write_tsv(
  GSE226643_abundance_long %>%
    dplyr::select(
      .data$Gene_id,
      .data$gene_symbol,
      .data$gene_length_nt,
      .data$Normalized_MFE,
      .data$raw_count_column,
      .data$analysis_id,
      .data$limitation,
      .data$author_condition,
      .data$DESeq2_stage,
      .data$time_h,
      .data$author_replicate_id,
      .data$relative_RNA_molecule_abundance
    ),
  file.path(
    GSE226643_table_dir,
    "GSE226643_full_transcriptome_gene_length_adjusted_relative_RNA_abundance.tsv.gz"
  )
)

GSE226643_MFE_sample_indices <- GSE226643_abundance_long %>%
  dplyr::group_by(
    .data$raw_count_column,
    .data$analysis_id,
    .data$limitation,
    .data$author_condition,
    .data$DESeq2_stage,
    .data$time_h,
    .data$author_replicate_id
  ) %>%
  dplyr::summarise(
    N_full_transcriptome_genes = dplyr::n(),
    abundance_weighted_Normalized_MFE = sum(
      .data$relative_RNA_molecule_abundance * .data$Normalized_MFE
    ),
    .groups = "drop"
  )

GSE226643_MFE_sample_indices_by_comparison <- GSE226643_MFE_sample_indices %>%
  dplyr::mutate(
    arm = dplyr::if_else(
      .data$DESeq2_stage == "Early_3h",
      "3 h early exponential",
      "12 h stationary"
    ),
    arm = factor(
      .data$arm,
      levels = c("3 h early exponential", "12 h stationary")
    )
  )

GSE226643_MFE_replicate_audit <- GSE226643_MFE_sample_indices_by_comparison %>%
  dplyr::count(.data$analysis_id, .data$limitation, .data$arm, name = "N_replicates")

if (any(GSE226643_MFE_replicate_audit$N_replicates != 2L)) {
  stop(
    "The abundance-weighted ΔGfold analysis does not have two replicates per arm.",
    call. = FALSE
  )
}

GSE226643_MFE_condition_indices <- GSE226643_MFE_sample_indices_by_comparison %>%
  dplyr::group_by(.data$analysis_id, .data$limitation, .data$arm) %>%
  dplyr::summarise(
    N_biological_replicates = dplyr::n(),
    mean_abundance_weighted_Normalized_MFE = mean(
      .data$abundance_weighted_Normalized_MFE
    ),
    sd_abundance_weighted_Normalized_MFE = stats::sd(
      .data$abundance_weighted_Normalized_MFE
    ),
    sem_abundance_weighted_Normalized_MFE =
      .data$sd_abundance_weighted_Normalized_MFE /
      sqrt(.data$N_biological_replicates),
    .groups = "drop"
  )

GSE226643_MFE_condition_profiles <- GSE226643_abundance_long %>%
  dplyr::group_by(
    .data$analysis_id,
    .data$limitation,
    .data$author_condition,
    .data$DESeq2_stage,
    .data$Gene_id,
    .data$Normalized_MFE
  ) %>%
  dplyr::summarise(
    mean_relative_RNA_molecule_abundance = mean(
      .data$relative_RNA_molecule_abundance
    ),
    .groups = "drop"
  )

GSE226643_MFE_early_profiles <- GSE226643_MFE_condition_profiles %>%
  dplyr::filter(.data$DESeq2_stage == "Early_3h") %>%
  dplyr::transmute(
    analysis_id = .data$analysis_id,
    limitation = .data$limitation,
    early_condition = .data$author_condition,
    Gene_id = .data$Gene_id,
    Normalized_MFE = .data$Normalized_MFE,
    early_3h_mean_relative_RNA_abundance =
      .data$mean_relative_RNA_molecule_abundance
  )

GSE226643_MFE_stationary_profiles <- GSE226643_MFE_condition_profiles %>%
  dplyr::filter(.data$DESeq2_stage == "Stationary_12h") %>%
  dplyr::transmute(
    analysis_id = .data$analysis_id,
    limitation = .data$limitation,
    stationary_condition = .data$author_condition,
    Gene_id = .data$Gene_id,
    Normalized_MFE = .data$Normalized_MFE,
    stationary_12h_mean_relative_RNA_abundance =
      .data$mean_relative_RNA_molecule_abundance
  )

GSE226643_MFE_comparison_gene_data <- GSE226643_MFE_early_profiles %>%
  dplyr::inner_join(
    GSE226643_MFE_stationary_profiles,
    by = c("analysis_id", "limitation", "Gene_id", "Normalized_MFE")
  ) %>%
  dplyr::mutate(
    contrast_weight_stationary_12h_minus_early_3h =
      .data$stationary_12h_mean_relative_RNA_abundance -
      .data$early_3h_mean_relative_RNA_abundance,
    abundance_weighted_MFE_contribution =
      .data$contrast_weight_stationary_12h_minus_early_3h *
      .data$Normalized_MFE
  )

GSE226643_MFE_delta_statistics <- GSE226643_MFE_comparison_gene_data %>%
  dplyr::group_by(
    .data$analysis_id,
    .data$limitation,
    .data$early_condition,
    .data$stationary_condition
  ) %>%
  dplyr::summarise(
    N_full_transcriptome_genes = dplyr::n(),
    early_3h_abundance_weighted_Normalized_MFE = sum(
      .data$early_3h_mean_relative_RNA_abundance * .data$Normalized_MFE
    ),
    stationary_12h_abundance_weighted_Normalized_MFE = sum(
      .data$stationary_12h_mean_relative_RNA_abundance * .data$Normalized_MFE
    ),
    delta_abundance_weighted_Normalized_MFE = sum(
      .data$abundance_weighted_MFE_contribution
    ),
    .groups = "drop"
  )

GSE226643_MFE_replicate_statistics <- GSE226643_MFE_sample_indices_by_comparison %>%
  dplyr::group_by(.data$analysis_id, .data$limitation) %>%
  dplyr::summarise(
    biological_replicate_Welch_t_p = tryCatch(
      stats::t.test(
        .data$abundance_weighted_Normalized_MFE[
          .data$arm == "12 h stationary"
        ],
        .data$abundance_weighted_Normalized_MFE[
          .data$arm == "3 h early exponential"
        ]
      )$p.value,
      error = function(error_message) NA_real_
    ),
    .groups = "drop"
  )

set.seed(20260807L)

GSE226643_MFE_resampling_statistics <- vector(
  "list",
  length(unique(GSE226643_MFE_comparison_gene_data$analysis_id))
)

GSE226643_MFE_analysis_ids <- unique(
  GSE226643_MFE_comparison_gene_data$analysis_id
)

for (GSE226643_MFE_index in seq_along(GSE226643_MFE_analysis_ids)) {
  GSE226643_one_analysis_id <- GSE226643_MFE_analysis_ids[GSE226643_MFE_index]

  GSE226643_one_MFE_gene_data <- GSE226643_MFE_comparison_gene_data %>%
    dplyr::filter(.data$analysis_id == GSE226643_one_analysis_id)

  GSE226643_observed_MFE_delta <- sum(
    GSE226643_one_MFE_gene_data$abundance_weighted_MFE_contribution
  )

  GSE226643_MFE_bootstrap_delta <- vapply(
    seq_len(GSE226643_B_BOOT),
    function(iteration) {
      sum(
        sample(
          GSE226643_one_MFE_gene_data$abundance_weighted_MFE_contribution,
          size = nrow(GSE226643_one_MFE_gene_data),
          replace = TRUE
        )
      )
    },
    numeric(1)
  )

  GSE226643_MFE_permutation_delta <- vapply(
    seq_len(GSE226643_B_PERM),
    function(iteration) {
      sum(
        GSE226643_one_MFE_gene_data$contrast_weight_stationary_12h_minus_early_3h *
          sample(
            GSE226643_one_MFE_gene_data$Normalized_MFE,
            size = nrow(GSE226643_one_MFE_gene_data),
            replace = FALSE
          )
      )
    },
    numeric(1)
  )

  GSE226643_MFE_resampling_statistics[[GSE226643_MFE_index]] <- tibble::tibble(
    analysis_id = GSE226643_one_analysis_id,
    gene_bootstrap_CI_low = unname(
      stats::quantile(GSE226643_MFE_bootstrap_delta, 0.025)
    ),
    gene_bootstrap_CI_high = unname(
      stats::quantile(GSE226643_MFE_bootstrap_delta, 0.975)
    ),
    descriptor_label_permutation_p = (
      sum(
        abs(GSE226643_MFE_permutation_delta) >=
          abs(GSE226643_observed_MFE_delta)
      ) + 1
    ) / (GSE226643_B_PERM + 1L)
  )
}

GSE226643_MFE_resampling_statistics <- dplyr::bind_rows(
  GSE226643_MFE_resampling_statistics
)

GSE226643_MFE_delta_statistics <- GSE226643_MFE_delta_statistics %>%
  dplyr::left_join(
    GSE226643_MFE_replicate_statistics,
    by = c("analysis_id", "limitation")
  ) %>%
  dplyr::left_join(
    GSE226643_MFE_resampling_statistics,
    by = "analysis_id"
  ) %>%
  dplyr::mutate(
    direction = dplyr::case_when(
      .data$delta_abundance_weighted_Normalized_MFE < 0 ~
        "12 h RNA pool is more negative / more stably folded",
      .data$delta_abundance_weighted_Normalized_MFE > 0 ~
        "12 h RNA pool is less negative / less stably folded",
      TRUE ~ "No change"
    ),
    B_boot = GSE226643_B_BOOT,
    B_perm = GSE226643_B_PERM
  )

print(GSE226643_MFE_delta_statistics)

readr::write_tsv(
  GSE226643_MFE_sample_indices_by_comparison,
  file.path(
    GSE226643_table_dir,
    "GSE226643_full_transcriptome_abundance_weighted_DGfold_sample_indices.tsv"
  )
)

readr::write_tsv(
  GSE226643_MFE_condition_indices,
  file.path(
    GSE226643_table_dir,
    "GSE226643_full_transcriptome_abundance_weighted_DGfold_condition_indices.tsv"
  )
)

readr::write_tsv(
  GSE226643_MFE_comparison_gene_data,
  file.path(
    GSE226643_table_dir,
    "GSE226643_full_transcriptome_abundance_weighted_DGfold_gene_contributions.tsv.gz"
  )
)

readr::write_tsv(
  GSE226643_MFE_delta_statistics,
  file.path(
    GSE226643_table_dir,
    "GSE226643_full_transcriptome_abundance_weighted_DGfold_delta_statistics.tsv"
  )
)

GSE226643_stationary_marker_panel <- tibble::tribble(
  ~Marker_module, ~Gene_symbol, ~Interpretive_role,

  "RpoS regulatory axis", "RPOS",
  "Master stationary-phase sigma factor; transcript abundance is not a direct activity readout",

  "RpoS regulatory axis", "IRAP",
  "Stabilises RpoS under starvation/stress",

  "RpoS regulatory axis", "RSSB",
  "Adaptor controlling RpoS degradation",

  "ppGpp / stringent-response axis", "RELA",
  "ppGpp synthetase activated by nutrient limitation",

  "ppGpp / stringent-response axis", "SPOT",
  "Bifunctional ppGpp synthetase/hydrolase",

  "ppGpp / stringent-response axis", "DKSA",
  "RNA-polymerase cofactor acting with ppGpp",

  "Stationary-phase / RpoS-associated outputs", "DPS",
  "DNA protection during starvation and stationary phase",

  "Stationary-phase / RpoS-associated outputs", "OSMY",
  "Canonical RpoS-responsive osmotic/stationary-phase marker",

  "Stationary-phase / RpoS-associated outputs", "KATE",
  "Catalase associated with the RpoS general-stress response",

  "Stationary-phase / RpoS-associated outputs", "CSID",
  "Carbon-starvation/stationary-phase associated gene",

  "Stationary-phase / RpoS-associated outputs", "HDEA",
  "Acid-stress survival protein induced in stationary phase",

  "Stationary-phase / RpoS-associated outputs", "HDEB",
  "Acid-stress survival protein induced in stationary phase",

  "Stationary-phase / RpoS-associated outputs", "GADA",
  "Glutamate-dependent acid-resistance system component",

  "Stationary-phase / RpoS-associated outputs", "GADB",
  "Glutamate-dependent acid-resistance system component",

  "Stationary-phase / RpoS-associated outputs", "GADC",
  "Glutamate antiporter associated with acid resistance",

  "Stationary-phase / RpoS-associated outputs", "BOLA",
  "Morphology and stationary-phase-associated regulator",

  "Stationary-phase / RpoS-associated outputs", "GLGS",
  "Stationary-phase glycogen-associated gene",

  "Stationary-phase / RpoS-associated outputs", "OTSA",
  "Trehalose-biosynthesis pathway component",

  "Translation hibernation", "RMF",
  "Ribosome modulation factor; stationary-phase ribosome hibernation",

  "Translation hibernation", "HPF",
  "Hibernation-promoting factor",

  "Translation hibernation", "RAIA",
  "Ribosome-associated inhibitor during stress/starvation",

  "SOS response", "RECA",
  "Central SOS recombinase/coprotase",

  "SOS response", "LEXA",
  "SOS transcriptional repressor",

  "SOS response", "SULA",
  "SOS-induced cell-division inhibitor",

  "SOS response", "DINB",
  "SOS-associated DNA polymerase IV",

  "SOS response", "RECN",
  "DNA-damage repair/SOS-associated protein",

  "SOS response", "UVRA",
  "Nucleotide-excision repair protein",

  "SOS response", "UVRB",
  "Nucleotide-excision repair protein",

  "SOS response", "UMUC",
  "SOS-associated translesion DNA-polymerase component",

  "SOS response", "UMUD",
  "SOS-associated translesion DNA-polymerase component",

  "SOS response", "TISB",
  "SOS-induced toxin involved in growth inhibition",

  "Growth / DNA-replication anchors", "DNAA",
  "Replication-initiation protein",

  "Growth / DNA-replication anchors", "DNAE",
  "DNA polymerase III alpha subunit",

  "Growth / DNA-replication anchors", "DNAN",
  "DNA polymerase III beta clamp",

  "Growth / DNA-replication anchors", "GYRA",
  "DNA gyrase subunit A",

  "Growth / DNA-replication anchors", "GYRB",
  "DNA gyrase subunit B",

  "Growth / cell-division anchors", "FTSA",
  "Cell-division protein",

  "Growth / cell-division anchors", "FTSI",
  "Cell-division transpeptidase",

  "Growth / cell-division anchors", "FTSZ",
  "Tubulin-like cell-division protein",

  "Translation / ribosome anchors", "RPLA",
  "Large ribosomal-subunit protein",

  "Translation / ribosome anchors", "RPLB",
  "Large ribosomal-subunit protein",

  "Translation / ribosome anchors", "RPLC",
  "Large ribosomal-subunit protein",

  "Translation / ribosome anchors", "RPLD",
  "Large ribosomal-subunit protein",

  "Translation / ribosome anchors", "RPLE",
  "Large ribosomal-subunit protein",

  "Translation / ribosome anchors", "RPSA",
  "Small ribosomal-subunit protein",

  "Translation / ribosome anchors", "RPSB",
  "Small ribosomal-subunit protein",

  "Translation / ribosome anchors", "RPSC",
  "Small ribosomal-subunit protein",

  "Translation / ribosome anchors", "RPSD",
  "Small ribosomal-subunit protein",

  "Translation / ribosome anchors", "RPSE",
  "Small ribosomal-subunit protein"
) %>%
  dplyr::distinct(
    .data$Marker_module,
    .data$Gene_symbol,
    .keep_all = TRUE
  ) %>%
  dplyr::mutate(
    Gene_symbol = toupper(trimws(.data$Gene_symbol))
  )

GSE226643_all_gene_DE_with_symbols <- dplyr::bind_rows(
  GSE226643_N_DE_all,
  GSE226643_S_DE_all
) %>%
  dplyr::left_join(
    GSE226643_gene_annotation %>%
      dplyr::select(
        .data$Gene_id,
        .data$gene_symbol
      ),
    by = "Gene_id"
  ) %>%
  dplyr::mutate(
    gene_symbol = toupper(trimws(as.character(.data$gene_symbol)))
  ) %>%
  dplyr::filter(
    !is.na(.data$gene_symbol),
    .data$gene_symbol != ""
  ) %>%
  dplyr::group_by(
    .data$analysis_id,
    .data$limitation,
    .data$gene_symbol
  ) %>%
  dplyr::slice_max(
    order_by = abs(.data$log2FC_stationary_12h_minus_early_3h),
    n = 1,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup()

GSE226643_marker_comparisons <- GSE226643_all_gene_DE_with_symbols %>%
  dplyr::distinct(
    .data$analysis_id,
    .data$limitation,
    .data$comparison_label
  )

GSE226643_targeted_marker_table <- tidyr::crossing(
  GSE226643_marker_comparisons,
  GSE226643_stationary_marker_panel
) %>%
  dplyr::left_join(
    GSE226643_all_gene_DE_with_symbols %>%
      dplyr::select(
        .data$analysis_id,
        .data$limitation,
        gene_symbol = .data$gene_symbol,
        .data$Gene_id,
        .data$base_mean,
        .data$log2FC_stationary_12h_minus_early_3h,
        .data$lfc_standard_error,
        .data$wald_statistic,
        .data$p_value,
        .data$padj,
        .data$early_3h_mean_normalised_counts,
        .data$stationary_12h_mean_normalised_counts
      ),
    by = c(
      "analysis_id",
      "limitation",
      "Gene_symbol" = "gene_symbol"
    )
  ) %>%
  dplyr::mutate(
    marker_present_in_DESeq2_results = !is.na(.data$Gene_id),

    DE_status = dplyr::case_when(
      !.data$marker_present_in_DESeq2_results ~
        "Marker absent from DESeq2 result table",

      is.na(.data$padj) ~
        "No adjusted P value",

      .data$padj < 0.05 &
        .data$log2FC_stationary_12h_minus_early_3h > 0 ~
        "Higher at 12 h (FDR < 0.05)",

      .data$padj < 0.05 &
        .data$log2FC_stationary_12h_minus_early_3h < 0 ~
        "Lower at 12 h (FDR < 0.05)",

      TRUE ~
        "Not FDR-significant"
    ),

    FC_direction = dplyr::case_when(
      !.data$marker_present_in_DESeq2_results ~ NA_character_,

      .data$log2FC_stationary_12h_minus_early_3h > 0 ~
        "Higher at 12 h",

      .data$log2FC_stationary_12h_minus_early_3h < 0 ~
        "Lower at 12 h",

      TRUE ~ "No estimated change"
    )
  ) %>%
  dplyr::arrange(
    .data$analysis_id,
    .data$Marker_module,
    .data$Gene_symbol
  )

print(GSE226643_targeted_marker_table)

GSE226643_targeted_marker_module_summary <- GSE226643_targeted_marker_table %>%
  dplyr::group_by(
    .data$analysis_id,
    .data$limitation,
    .data$comparison_label,
    .data$Marker_module
  ) %>%
  dplyr::summarise(
    N_pre_specified_markers = dplyr::n(),

    N_markers_present_in_DESeq2 = sum(
      .data$marker_present_in_DESeq2_results
    ),

    N_higher_at_12h_FDR_lt_0_05 = sum(
      .data$DE_status == "Higher at 12 h (FDR < 0.05)"
    ),

    N_lower_at_12h_FDR_lt_0_05 = sum(
      .data$DE_status == "Lower at 12 h (FDR < 0.05)"
    ),

    median_log2FC_among_present_markers = dplyr::if_else(
      all(is.na(.data$log2FC_stationary_12h_minus_early_3h)),
      NA_real_,
      stats::median(
        .data$log2FC_stationary_12h_minus_early_3h,
        na.rm = TRUE
      )
    ),

    higher_at_12h_FDR_marker_symbols = paste(
      .data$Gene_symbol[
        .data$DE_status == "Higher at 12 h (FDR < 0.05)"
      ],
      collapse = "; "
    ),

    lower_at_12h_FDR_marker_symbols = paste(
      .data$Gene_symbol[
        .data$DE_status == "Lower at 12 h (FDR < 0.05)"
      ],
      collapse = "; "
    ),

    .groups = "drop"
  ) %>%
  dplyr::arrange(
    .data$analysis_id,
    .data$Marker_module
  )

print(GSE226643_targeted_marker_module_summary)

readr::write_tsv(
  GSE226643_stationary_marker_panel,
  file.path(
    GSE226643_table_dir,
    "GSE226643_pre_specified_RpoS_ppGpp_stationary_marker_panel.tsv"
  )
)

readr::write_tsv(
  GSE226643_targeted_marker_table,
  file.path(
    GSE226643_table_dir,
    "GSE226643_targeted_RpoS_ppGpp_stationary_phase_marker_table.tsv"
  )
)

readr::write_tsv(
  GSE226643_targeted_marker_module_summary,
  file.path(
    GSE226643_table_dir,
    "GSE226643_targeted_RpoS_ppGpp_stationary_phase_marker_module_summary.tsv"
  )
)

GSE226643_analysis_objects <- list(
  primary_manifest = GSE226643_primary_manifest,
  annotation_audit = GSE226643_thermodynamic_annotation_audit,
  DE_summary = GSE226643_DE_summary,
  N_DE_all = GSE226643_N_DE_all,
  S_DE_all = GSE226643_S_DE_all,
  GSEA_results = GSE226643_GSEA_results,
  GSE226643_GSEA_bubble_data = GSE226643_GSEA_bubble_data,
  GSE226643_N_sig_thermodynamic = GSE226643_N_sig_thermodynamic,
  GSE226643_S_sig_thermodynamic = GSE226643_S_sig_thermodynamic,
  GSE226643_significant_distribution_statistics =
    GSE226643_significant_distribution_statistics,
  GSE226643_significant_DE_abundance_fraction =
    GSE226643_significant_DE_abundance_fraction,
  GSE226643_MFE_sample_indices_by_comparison =
    GSE226643_MFE_sample_indices_by_comparison,
  GSE226643_MFE_delta_statistics = GSE226643_MFE_delta_statistics,
  GSE226643_MFE_comparison_gene_data = GSE226643_MFE_comparison_gene_data,
  pre_specified_stationary_marker_panel = GSE226643_stationary_marker_panel,
  targeted_RpoS_ppGpp_stationary_marker_table = GSE226643_targeted_marker_table,
  targeted_RpoS_ppGpp_stationary_marker_module_summary =
    GSE226643_targeted_marker_module_summary
)

saveRDS(
  GSE226643_analysis_objects,
  file.path(
    GSE226643_output_dir,
    "GSE226643_analysis_results_for_Fig_5.rds"
  )
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(
    GSE226643_output_dir,
    "GSE226643_nutrient_limitation_analysis_sessionInfo.txt"
  )
)

message(
  "\nGSE226643 nutrient-limitation analysis is complete.\n",
  "Results directory: ", GSE226643_output_dir, "\n",
  "The full-transcriptome ΔGfold results use no DE filter."
)
