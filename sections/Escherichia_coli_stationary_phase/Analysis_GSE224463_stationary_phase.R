options(timeout = max(1200, getOption("timeout")))

suppressPackageStartupMessages({
  library(GEOquery)
  library(jsonlite)
  library(DESeq2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
})

GSE224463_project_dir <- "path_to/GSE224463"
GSE224463_prediction_file <- "path_to/Prediction_Ecoli_MFE.csv"
GSE224463_gene_info_file <- "path_to/MG1655_gene_info.csv"
GSE224463_dir <- GSE224463_project_dir

GSE224463_metadata_dir <- file.path(GSE224463_dir, "metadata")
GSE224463_geo_dir <- file.path(GSE224463_dir, "GEO_download")
GSE224463_zenodo_dir <- file.path(GSE224463_dir, "Zenodo_8083605")
GSE224463_sra_dir <- file.path(GSE224463_dir, "SRA_selected_samples")
GSE224463_fastq_dir <- file.path(GSE224463_sra_dir, "FASTQ")

for (one_dir in c(
  GSE224463_dir,
  GSE224463_metadata_dir,
  GSE224463_geo_dir,
  GSE224463_zenodo_dir,
  GSE224463_sra_dir,
  GSE224463_fastq_dir
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

GSE224463_selected_manifest <- data.frame(
  GSM = c(
    "GSM7024604", "GSM7024605", "GSM7024606", "GSM7024607", "GSM7024608", "GSM7024609",
    "GSM7024610", "GSM7024611", "GSM7024612", "GSM7024613", "GSM7024614", "GSM7024615",
    "GSM7024616", "GSM7024617", "GSM7024618", "GSM7024619", "GSM7024620", "GSM7024621",
    "GSM7024640", "GSM7024641", "GSM7024642", "GSM7024643", "GSM7024644", "GSM7024645",
    "GSM7024646", "GSM7024647", "GSM7024648", "GSM7024649", "GSM7024650", "GSM7024651",
    "GSM7024652", "GSM7024653", "GSM7024654", "GSM7024655", "GSM7024656", "GSM7024657"
  ),
  series = c(
    rep("SP_short", 18L),
    rep("SP_long", 18L)
  ),
  biological_replicate = c(
    rep("SP_A", 6L), rep("SP_B", 6L), rep("SP_C", 6L),
    rep("LSP_A", 6L), rep("LSP_C", 6L), rep("LSP_E", 6L)
  ),
  time_min = c(
    rep(c(0L, 30L, 60L, 90L, 120L, 150L), 3L),
    rep(c(0L, 90L, 120L, 180L, 240L, 300L), 3L)
  ),
  stringsAsFactors = FALSE
)

GSE224463_selected_manifest$analysis_role <- "time_course_only"
GSE224463_selected_manifest$analysis_role[
  GSE224463_selected_manifest$series == "SP_short" &
    GSE224463_selected_manifest$time_min %in% c(0L, 150L)
] <- "primary_SP_150min_vs_0min"
GSE224463_selected_manifest$analysis_role[
  GSE224463_selected_manifest$series == "SP_long" &
    GSE224463_selected_manifest$time_min %in% c(0L, 300L)
] <- "confirmatory_LSP_300min_vs_0min"

GSE224463_selected_manifest$condition <- NA_character_
GSE224463_selected_manifest$condition[
  GSE224463_selected_manifest$time_min == 0L
] <- "Exponential_reference"
GSE224463_selected_manifest$condition[
  GSE224463_selected_manifest$time_min > 0L
] <- "Stationary_phase_entry"
GSE224463_selected_manifest$genotype <- "wild_type"

stopifnot(
  nrow(GSE224463_selected_manifest) == 36L,
  sum(GSE224463_selected_manifest$analysis_role == "primary_SP_150min_vs_0min") == 6L,
  sum(GSE224463_selected_manifest$analysis_role == "confirmatory_LSP_300min_vs_0min") == 6L
)

write_tsv(
  GSE224463_selected_manifest,
  file.path(GSE224463_metadata_dir, "GSE224463_WT_stationary_selected_manifest.tsv")
)

GSE224463_accession <- "GSE224463"

GSE224463_GEO <- GEOquery::getGEO(
  GSE224463_accession,
  GSEMatrix = FALSE,
  getGPL = FALSE
)

GSE224463_GSM_list <- GEOquery::GSMList(GSE224463_GEO)

GSE224463_GEO_metadata <- do.call(
  rbind,
  lapply(names(GSE224463_GSM_list), function(one_GSM) {
    one_meta <- GEOquery::Meta(GSE224463_GSM_list[[one_GSM]])
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
  GSE224463_GEO_metadata,
  file.path(GSE224463_metadata_dir, "GSE224463_GEO_sample_metadata_all.tsv")
)

GSE224463_selected_metadata <- merge(
  GSE224463_selected_manifest,
  GSE224463_GEO_metadata,
  by = "GSM",
  all.x = TRUE,
  sort = FALSE
)
GSE224463_selected_metadata <- GSE224463_selected_metadata[
  match(GSE224463_selected_manifest$GSM, GSE224463_selected_metadata$GSM),
  , drop = FALSE
]

if (any(is.na(GSE224463_selected_metadata$title))) {
  stop("At least one selected GSM accession was not found in GEO metadata.", call. = FALSE)
}

write_tsv(
  GSE224463_selected_metadata,
  file.path(GSE224463_metadata_dir, "GSE224463_WT_stationary_selected_metadata.tsv")
)

GSE224463_supplementary_status <- tryCatch(
  {
    GEOquery::getGEOSuppFiles(
      GSE224463_accession,
      baseDir = GSE224463_geo_dir,
      makeDirectory = FALSE
    )
    "GEO supplementary-file download completed."
  },
  error = function(e) paste("GEO supplementary-file download failed:", conditionMessage(e))
)
writeLines(
  GSE224463_supplementary_status,
  con = file.path(GSE224463_metadata_dir, "GSE224463_GEO_supplementary_download_status.txt")
)

download_zenodo_record <- function(record_id, destination_dir) {
  record_api <- paste0("https://zenodo.org/api/records/", record_id)
  record_json_file <- file.path(destination_dir, paste0("Zenodo_", record_id, "_record.json"))
  download_if_needed(record_api, record_json_file)

  record <- jsonlite::fromJSON(record_json_file, simplifyVector = FALSE)
  if (is.null(record$files) || length(record$files) == 0L) {
    stop("The Zenodo record contains no downloadable files: ", record_id, call. = FALSE)
  }

  file_manifest <- do.call(rbind, lapply(record$files, function(one_file) {
    all_links <- unlist(one_file$links, use.names = FALSE)
    download_url <- all_links[grepl("^https?://", all_links)][1L]
    data.frame(
      filename = as.character(one_file$key),
      expected_bytes = as.numeric(one_file$size),
      checksum = as_one_string(one_file$checksum),
      download_url = as.character(download_url),
      stringsAsFactors = FALSE
    )
  }))

  if (any(is.na(file_manifest$download_url) | file_manifest$download_url == "")) {
    stop("A Zenodo-file download URL could not be identified.", call. = FALSE)
  }

  file_manifest$local_file <- file.path(destination_dir, basename(file_manifest$filename))

  for (i in seq_len(nrow(file_manifest))) {
    download_if_needed(
      url = file_manifest$download_url[i],
      destination = file_manifest$local_file[i],
      expected_size = file_manifest$expected_bytes[i]
    )
  }

  write_tsv(
    file_manifest,
    file.path(destination_dir, paste0("Zenodo_", record_id, "_download_manifest.tsv"))
  )

  invisible(file_manifest)
}

GSE224463_Zenodo_manifest <- download_zenodo_record(
  record_id = "8083605",
  destination_dir = GSE224463_zenodo_dir
)

GSE224463_Zenodo_inventory <- data.frame(
  local_file = list.files(GSE224463_zenodo_dir, recursive = TRUE, full.names = TRUE),
  stringsAsFactors = FALSE
)
GSE224463_Zenodo_inventory$file_name <- basename(GSE224463_Zenodo_inventory$local_file)
GSE224463_Zenodo_inventory$bytes <- file.info(GSE224463_Zenodo_inventory$local_file)$size
GSE224463_Zenodo_inventory$likely_content <- ifelse(
  grepl("raw.*count|count.*matrix|featurecount", GSE224463_Zenodo_inventory$file_name, ignore.case = TRUE),
  "candidate_raw_count_matrix",
  ifelse(
    grepl("metadata|sample|design", GSE224463_Zenodo_inventory$file_name, ignore.case = TRUE),
    "candidate_metadata",
    ifelse(
      grepl("script|code|\\.R$", GSE224463_Zenodo_inventory$file_name, ignore.case = TRUE),
      "code_or_script",
      "other"
    )
  )
)
write_tsv(
  GSE224463_Zenodo_inventory,
  file.path(GSE224463_metadata_dir, "GSE224463_Zenodo_file_inventory.tsv")
)

GSE224463_selected_metadata$SRX <- vapply(
  GSE224463_selected_metadata$relation,
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

GSE224463_SRA_manifest <- do.call(
  rbind,
  Map(
    fetch_ena_runs,
    SRX = GSE224463_selected_metadata$SRX,
    GSM = GSE224463_selected_metadata$GSM
  )
)

write_tsv(
  GSE224463_SRA_manifest,
  file.path(GSE224463_sra_dir, "GSE224463_selected_ENA_SRA_run_manifest.tsv")
)

if (all(is.na(GSE224463_SRA_manifest$run_accession))) {
  warning(
    "No ENA runs were resolved automatically. Use the GEO sample metadata ",
    "and SRA Run Selector if raw FASTQ reprocessing is required.",
    call. = FALSE
  )
}

if (download_selected_fastq) {
  all_fastq_paths <- unlist(
    strsplit(
      GSE224463_SRA_manifest$fastq_ftp[!is.na(GSE224463_SRA_manifest$fastq_ftp)],
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
      destination = file.path(GSE224463_fastq_dir, basename(one_fastq_path))
    )
  }
}

writeLines(
  capture.output(sessionInfo()),
  con = file.path(GSE224463_metadata_dir, "GSE224463_download_sessionInfo.txt")
)

message(
  "GSE224463 download and manifests completed.\n",
  "Next, inspect metadata/GSE224463_Zenodo_file_inventory.tsv and identify ",
  "the author-supplied raw gene-by-sample count matrix before starting DE."
)

GSE224463_metadata_dir <- file.path(
  GSE224463_project_dir,
  "metadata"
)

GSE224463_audit_dir <- file.path(
  GSE224463_project_dir,
  "analysis_manifests",
  "00_stationary_phase_count_audit"
)

dir.create(
  GSE224463_audit_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!dir.exists(GSE224463_project_dir)) {
  stop(
    "GSE224463_project_dir does not exist:\n",
    GSE224463_project_dir,
    "\nUpdate the path in Section 1 before continuing.",
    call. = FALSE
  )
}

GSE224463_geo_metadata_file <- file.path(
  GSE224463_metadata_dir,
  "GSE224463_GEO_sample_metadata_all.tsv"
)

if (!file.exists(GSE224463_geo_metadata_file)) {
  stop(
    "Required GEO metadata file was not found:\n",
    GSE224463_geo_metadata_file,
    "\nRun the GSE224463 download script first, or place the file in metadata/.",
    call. = FALSE
  )
}

GSE224463_geo_metadata <- readr::read_tsv(
  GSE224463_geo_metadata_file,
  show_col_types = FALSE,
  progress = FALSE,
  name_repair = "unique"
)

names(GSE224463_geo_metadata) <- make.names(
  names(GSE224463_geo_metadata),
  unique = TRUE
)

required_GEO_columns <- c("GSM", "title")

missing_GEO_columns <- setdiff(
  required_GEO_columns,
  names(GSE224463_geo_metadata)
)

if (length(missing_GEO_columns) > 0L) {
  stop(
    "The GEO metadata table lacks required column(s): ",
    paste(missing_GEO_columns, collapse = ", "),
    call. = FALSE
  )
}

GSE224463_geo_metadata <- GSE224463_geo_metadata %>%
  dplyr::mutate(
    GSM = trimws(as.character(.data$GSM)),
    title = trimws(as.character(.data$title))
  )

GSE224463_WT_stationary_GEO <- GSE224463_geo_metadata %>%
  dplyr::filter(
    stringr::str_detect(
      .data$title,
      stringr::regex(
        "^wild-type,\\s*stationary phase,",
        ignore_case = TRUE
      )
    )
  ) %>%
  dplyr::mutate(
    time_min = as.integer(
      stringr::str_match(
        .data$title,
        stringr::regex(
          ",\\s*([0-9]+)min,",
          ignore_case = TRUE
        )
      )[, 2]
    ),
    author_replicate = stringr::str_match(
      .data$title,
      stringr::regex(
        "replicate\\s+(SP_[A-Z]|LSP_[A-Z])\\s*$",
        ignore_case = TRUE
      )
    )[, 2],
    author_replicate = toupper(.data$author_replicate),
    series = dplyr::case_when(
      stringr::str_detect(.data$author_replicate, "^SP_") ~ "short",
      stringr::str_detect(.data$author_replicate, "^LSP_") ~ "long",
      TRUE ~ NA_character_
    ),
    condition = dplyr::if_else(
      .data$time_min == 0L,
      "Exponential reference (0 min)",
      "Stationary-phase entry"
    ),
    genotype = "wild_type"
  ) %>%
  dplyr::select(
    GSM,
    title,
    genotype,
    series,
    author_replicate,
    time_min,
    condition,
    dplyr::everything()
  ) %>%
  dplyr::arrange(
    .data$series,
    .data$author_replicate,
    .data$time_min
  )

if (nrow(GSE224463_WT_stationary_GEO) != 36L) {
  stop(
    "Expected exactly 36 wild-type stationary-phase GEO samples, but found ",
    nrow(GSE224463_WT_stationary_GEO),
    ". Inspect the title parsing before continuing.",
    call. = FALSE
  )
}

if (anyNA(GSE224463_WT_stationary_GEO$time_min) ||
    anyNA(GSE224463_WT_stationary_GEO$author_replicate) ||
    anyNA(GSE224463_WT_stationary_GEO$series)) {
  stop(
    "At least one wild-type stationary-phase GEO title could not be parsed ",
    "into a time point, replicate, and series.",
    call. = FALSE
  )
}

expected_short_times <- c(0L, 30L, 60L, 90L, 120L, 150L)
expected_long_times <- c(0L, 90L, 120L, 180L, 240L, 300L)

expected_WT_stationary_geometry <- dplyr::bind_rows(
  tidyr::expand_grid(
    series = "short",
    author_replicate = c("SP_A", "SP_B", "SP_C"),
    time_min = expected_short_times
  ),
  tidyr::expand_grid(
    series = "long",
    author_replicate = c("LSP_A", "LSP_C", "LSP_E"),
    time_min = expected_long_times
  )
)

observed_WT_stationary_geometry <- GSE224463_WT_stationary_GEO %>%
  dplyr::select(
    series,
    author_replicate,
    time_min
  )

missing_WT_stationary_GEO_rows <- dplyr::anti_join(
  expected_WT_stationary_geometry,
  observed_WT_stationary_geometry,
  by = c("series", "author_replicate", "time_min")
)

unexpected_WT_stationary_GEO_rows <- dplyr::anti_join(
  observed_WT_stationary_geometry,
  expected_WT_stationary_geometry,
  by = c("series", "author_replicate", "time_min")
)

if (nrow(missing_WT_stationary_GEO_rows) > 0L ||
    nrow(unexpected_WT_stationary_GEO_rows) > 0L) {
  stop(
    "The parsed wild-type stationary-phase GEO design does not match the ",
    "expected author design. Inspect the audit tables written before this error ",
    "or the GEO metadata source.",
    call. = FALSE
  )
}

readr::write_csv(
  GSE224463_WT_stationary_GEO,
  file.path(
    GSE224463_audit_dir,
    "GSE224463_WT_stationary_phase_GEO_manifest_author_labels.csv"
  ),
  na = ""
)

print(GSE224463_WT_stationary_GEO)

GSE224463_all_count_files <- list.files(
  GSE224463_project_dir,
  pattern = "_counts\\.txt$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(GSE224463_all_count_files) == 0L) {
  stop(
    "No *_counts.txt files were found under:\n",
    GSE224463_project_dir,
    call. = FALSE
  )
}

GSE224463_all_count_files <- normalizePath(
  GSE224463_all_count_files,
  winslash = "/",
  mustWork = TRUE
)

GSE224463_project_dir_normalized <- normalizePath(
  GSE224463_project_dir,
  winslash = "/",
  mustWork = TRUE
)

GSE224463_stationary_count_files <- GSE224463_all_count_files[
  stringr::str_detect(
    GSE224463_all_count_files,
    "/moc1533/"
  ) &
    stringr::str_detect(
      basename(GSE224463_all_count_files),
      "^JH(01|02)_(?:LSP_)?[A-Z]_?[0-9]{2}_counts\\.txt$"
    )
]

if (length(GSE224463_stationary_count_files) != 72L) {
  stop(
    "Expected 72 stationary-phase count files in moc1533, but found ",
    length(GSE224463_stationary_count_files),
    ". Do not continue until the Zenodo archive structure has been checked.",
    call. = FALSE
  )
}

GSE224463_stationary_file_inventory <- tibble::tibble(
  count_file = GSE224463_stationary_count_files,
  file_name = basename(GSE224463_stationary_count_files),
  bytes = as.numeric(file.info(GSE224463_stationary_count_files)$size)
) %>%
  dplyr::mutate(
    relative_count_file = sub(
      paste0(GSE224463_project_dir_normalized, "/"),
      "",
      .data$count_file,
      fixed = TRUE
    ),
    JH_prefix = stringr::str_match(
      .data$file_name,
      "^(JH[0-9]{2})_"
    )[, 2],
    is_long_series = stringr::str_detect(
      .data$file_name,
      "^JH[0-9]{2}_LSP_"
    ),
    series = dplyr::if_else(
      .data$is_long_series,
      "long",
      "short"
    ),
    replicate_letter = dplyr::case_when(
      .data$is_long_series ~ stringr::str_match(
        .data$file_name,
        "^JH[0-9]{2}_LSP_([A-Z])_[0-9]{2}_counts\\.txt$"
      )[, 2],
      TRUE ~ stringr::str_match(
        .data$file_name,
        "^JH[0-9]{2}_([A-Z])[0-9]{2}_counts\\.txt$"
      )[, 2]
    ),
    time_index = as.integer(
      stringr::str_match(
        .data$file_name,
        "([0-9]{2})_counts\\.txt$"
      )[, 2]
    ),
    author_replicate = dplyr::case_when(
      .data$series == "short" ~ paste0("SP_", .data$replicate_letter),
      .data$series == "long" ~ paste0("LSP_", .data$replicate_letter),
      TRUE ~ NA_character_
    ),
    time_min = dplyr::case_when(
      .data$series == "short" ~ expected_short_times[.data$time_index],
      .data$series == "long" ~ expected_long_times[.data$time_index],
      TRUE ~ NA_integer_
    )
  ) %>%
  dplyr::arrange(
    .data$series,
    .data$JH_prefix,
    .data$author_replicate,
    .data$time_min
  )

if (anyNA(GSE224463_stationary_file_inventory$JH_prefix) ||
    anyNA(GSE224463_stationary_file_inventory$replicate_letter) ||
    anyNA(GSE224463_stationary_file_inventory$time_index) ||
    anyNA(GSE224463_stationary_file_inventory$time_min)) {
  stop(
    "At least one stationary count filename could not be parsed completely.",
    call. = FALSE
  )
}

if (any(!GSE224463_stationary_file_inventory$time_index %in% 1:6)) {
  stop(
    "A stationary count filename has a time index outside 01--06.",
    call. = FALSE
  )
}

expected_stationary_file_geometry <- dplyr::bind_rows(
  tidyr::expand_grid(
    JH_prefix = c("JH01", "JH02"),
    series = "short",
    author_replicate = c("SP_A", "SP_B", "SP_C"),
    time_index = 1:6
  ) %>%
    dplyr::mutate(
      time_min = expected_short_times[.data$time_index]
    ),
  tidyr::expand_grid(
    JH_prefix = c("JH01", "JH02"),
    series = "long",
    author_replicate = c("LSP_A", "LSP_C", "LSP_E"),
    time_index = 1:6
  ) %>%
    dplyr::mutate(
      time_min = expected_long_times[.data$time_index]
    )
)

observed_stationary_file_geometry <- GSE224463_stationary_file_inventory %>%
  dplyr::select(
    JH_prefix,
    series,
    author_replicate,
    time_index,
    time_min
  )

missing_stationary_count_files <- dplyr::anti_join(
  expected_stationary_file_geometry,
  observed_stationary_file_geometry,
  by = c(
    "JH_prefix",
    "series",
    "author_replicate",
    "time_index",
    "time_min"
  )
)

unexpected_stationary_count_files <- dplyr::anti_join(
  observed_stationary_file_geometry,
  expected_stationary_file_geometry,
  by = c(
    "JH_prefix",
    "series",
    "author_replicate",
    "time_index",
    "time_min"
  )
)

if (nrow(missing_stationary_count_files) > 0L ||
    nrow(unexpected_stationary_count_files) > 0L) {
  stop(
    "The stationary count filenames do not form the expected 2 x 2 x 3 x 6 ",
    "design. Inspect the parsed file inventory before continuing.",
    call. = FALSE
  )
}

GSE224463_stationary_count_tables <- lapply(
  seq_len(nrow(GSE224463_stationary_file_inventory)),
  function(i) {
    one_file <- GSE224463_stationary_file_inventory$count_file[i]
    one_file_name <- GSE224463_stationary_file_inventory$file_name[i]

    one_count_table <- readr::read_tsv(
      one_file,
      col_names = FALSE,
      show_col_types = FALSE,
      progress = FALSE,
      name_repair = "minimal"
    )

    if (ncol(one_count_table) != 2L) {
      stop(
        "Expected exactly two columns in raw count file: ",
        one_file_name,
        call. = FALSE
      )
    }

    colnames(one_count_table) <- c("Gene_id", "raw_count")

    one_count_table <- one_count_table %>%
      dplyr::transmute(
        Gene_id = trimws(as.character(.data$Gene_id)),
        raw_count = suppressWarnings(as.numeric(.data$raw_count))
      )

    if (anyNA(one_count_table$Gene_id) ||
        any(one_count_table$Gene_id == "")) {
      stop(
        "A blank or missing gene identifier was found in: ",
        one_file_name,
        call. = FALSE
      )
    }

    if (anyDuplicated(one_count_table$Gene_id) > 0L) {
      stop(
        "Duplicate gene identifiers were found in: ",
        one_file_name,
        call. = FALSE
      )
    }

    if (any(!is.finite(one_count_table$raw_count)) ||
        any(one_count_table$raw_count < 0) ||
        any(abs(one_count_table$raw_count - round(one_count_table$raw_count)) > 1e-8)) {
      stop(
        "Counts must be finite, non-negative integers. Validation failed for: ",
        one_file_name,
        call. = FALSE
      )
    }

    one_count_table
  }
)

reference_gene_ids <- GSE224463_stationary_count_tables[[1L]]$Gene_id

GSE224463_stationary_gene_set_audit <- dplyr::bind_rows(
  lapply(
    seq_along(GSE224463_stationary_count_tables),
    function(i) {
      one_count_table <- GSE224463_stationary_count_tables[[i]]

      tibble::tibble(
        file_name = GSE224463_stationary_file_inventory$file_name[i],
        n_genes = nrow(one_count_table),
        gene_order_matches_first_file = identical(
          one_count_table$Gene_id,
          reference_gene_ids
        ),
        gene_set_matches_first_file = setequal(
          one_count_table$Gene_id,
          reference_gene_ids
        )
      )
    }
  )
)

if (any(!GSE224463_stationary_gene_set_audit$gene_set_matches_first_file)) {
  stop(
    "At least one stationary count table has a different gene set. Do not ",
    "construct a count matrix until this author-data inconsistency is resolved.",
    call. = FALSE
  )
}

GSE224463_stationary_count_matrix <- vapply(
  GSE224463_stationary_count_tables,
  function(one_count_table) {
    one_count_table$raw_count[
      match(
        reference_gene_ids,
        one_count_table$Gene_id
      )
    ]
  },
  numeric(length(reference_gene_ids))
)

rownames(GSE224463_stationary_count_matrix) <- reference_gene_ids
colnames(GSE224463_stationary_count_matrix) <-
  GSE224463_stationary_file_inventory$file_name

if (anyNA(GSE224463_stationary_count_matrix)) {
  stop(
    "The stationary count matrix contains missing values after gene alignment.",
    call. = FALSE
  )
}

GSE224463_stationary_library_QC <- GSE224463_stationary_file_inventory %>%
  dplyr::mutate(
    raw_library_size = as.numeric(
      colSums(GSE224463_stationary_count_matrix)
    ),
    genes_detected = as.integer(
      colSums(GSE224463_stationary_count_matrix > 0)
    ),
    percent_zero_counts = round(
      100 * colMeans(GSE224463_stationary_count_matrix == 0),
      3
    )
  )

if (any(GSE224463_stationary_library_QC$raw_library_size <= 0)) {
  stop(
    "At least one stationary-phase raw count file has a zero total library ",
    "size. Do not interpret the rpoS CPM audit until this is resolved.",
    call. = FALSE
  )
}

rpoS_locus_tag <- "b2741"

if (!rpoS_locus_tag %in% rownames(GSE224463_stationary_count_matrix)) {
  stop(
    "The expected rpoS locus tag (", rpoS_locus_tag,
    ") was not found in the stationary count matrix.",
    call. = FALSE
  )
}

GSE224463_stationary_library_QC <- GSE224463_stationary_library_QC %>%
  dplyr::mutate(
    rpoS_raw_count = as.numeric(
      GSE224463_stationary_count_matrix[
        rpoS_locus_tag,
        .data$file_name
      ]
    ),
    rpoS_counts_per_million = 1e6 * .data$rpoS_raw_count / .data$raw_library_size
  )

GSE224463_rpoS_prefix_audit <- GSE224463_stationary_library_QC %>%
  dplyr::group_by(
    .data$JH_prefix,
    .data$series
  ) %>%
  dplyr::summarise(
    N_libraries = dplyr::n(),
    median_rpoS_raw_count = median(.data$rpoS_raw_count),
    minimum_rpoS_raw_count = min(.data$rpoS_raw_count),
    maximum_rpoS_raw_count = max(.data$rpoS_raw_count),
    median_rpoS_CPM = median(.data$rpoS_counts_per_million),
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    .data$series,
    .data$JH_prefix
  )

readr::write_csv(
  GSE224463_stationary_library_QC,
  file.path(
    GSE224463_audit_dir,
    "GSE224463_stationary_phase_raw_count_library_QC.csv"
  ),
  na = ""
)

readr::write_csv(
  GSE224463_stationary_gene_set_audit,
  file.path(
    GSE224463_audit_dir,
    "GSE224463_stationary_phase_gene_set_audit.csv"
  ),
  na = ""
)

readr::write_csv(
  GSE224463_rpoS_prefix_audit,
  file.path(
    GSE224463_audit_dir,
    "GSE224463_stationary_phase_rpoS_JH_prefix_audit.csv"
  ),
  na = ""
)

print(GSE224463_rpoS_prefix_audit)

GSE224463_stationary_log_count_correlation <- stats::cor(
  log2(GSE224463_stationary_count_matrix + 1),
  method = "pearson"
)

GSE224463_stationary_log_count_correlation_long <- as.data.frame(
  GSE224463_stationary_log_count_correlation,
  check.names = FALSE
) %>%
  tibble::rownames_to_column("file_name_1") %>%
  tidyr::pivot_longer(
    cols = -dplyr::all_of("file_name_1"),
    names_to = "file_name_2",
    values_to = "Pearson_r"
  ) %>%
  dplyr::left_join(
    GSE224463_stationary_library_QC %>%
      dplyr::select(
        file_name,
        JH_prefix,
        series,
        author_replicate,
        time_min
      ) %>%
      dplyr::rename_with(
        ~ paste0(.x, "_1"),
        -file_name
      ) %>%
      dplyr::rename(file_name_1 = file_name),
    by = "file_name_1"
  ) %>%
  dplyr::left_join(
    GSE224463_stationary_library_QC %>%
      dplyr::select(
        file_name,
        JH_prefix,
        series,
        author_replicate,
        time_min
      ) %>%
      dplyr::rename_with(
        ~ paste0(.x, "_2"),
        -file_name
      ) %>%
      dplyr::rename(file_name_2 = file_name),
    by = "file_name_2"
  ) %>%
  dplyr::arrange(
    .data$file_name_1,
    .data$file_name_2
  )

readr::write_tsv(
  GSE224463_stationary_log_count_correlation_long,
  file.path(
    GSE224463_audit_dir,
    "GSE224463_stationary_phase_log_count_correlations.tsv"
  ),
  na = ""
)

GSE224463_stationary_raw_PCA <- stats::prcomp(
  t(log2(GSE224463_stationary_count_matrix + 1)),
  center = TRUE,
  scale. = FALSE
)

GSE224463_stationary_raw_PCA_variance <- 100 *
  GSE224463_stationary_raw_PCA$sdev^2 /
  sum(GSE224463_stationary_raw_PCA$sdev^2)

GSE224463_stationary_raw_PCA_coordinates <- as.data.frame(
  GSE224463_stationary_raw_PCA$x[, 1:2, drop = FALSE]
) %>%
  tibble::rownames_to_column("file_name") %>%
  dplyr::as_tibble() %>%
  dplyr::left_join(
    GSE224463_stationary_library_QC %>%
      dplyr::select(
        file_name,
        JH_prefix,
        series,
        author_replicate,
        time_min,
        raw_library_size,
        rpoS_counts_per_million
      ),
    by = "file_name"
  ) %>%
  dplyr::arrange(
    .data$series,
    .data$JH_prefix,
    .data$author_replicate,
    .data$time_min
  )

readr::write_tsv(
  GSE224463_stationary_raw_PCA_coordinates,
  file.path(
    GSE224463_audit_dir,
    "GSE224463_stationary_phase_raw_log_count_PCA_coordinates.tsv"
  ),
  na = ""
)

GSE224463_stationary_library_plot <- ggplot(
  GSE224463_stationary_library_QC,
  aes(
    x = .data$author_replicate,
    y = .data$raw_library_size,
    colour = .data$JH_prefix,
    shape = .data$JH_prefix
  )
) +
  geom_point(
    size = 2.6,
    position = position_jitter(width = 0.08, height = 0)
  ) +
  facet_grid(
    rows = vars(series),
    cols = vars(time_min),
    scales = "free_x",
    space = "free_x"
  ) +
  scale_y_continuous(
    trans = "log10"
  ) +
  theme_classic() +
  theme(
    legend.position = "top",
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      colour = "black"
    ),
    strip.text = element_text(size = 8)
  ) +
  labs(
    title = "GSE224463 stationary-phase raw-library audit",
    x = "Author biological replicate",
    y = "Raw library size (log10 scale)",
    colour = "Filename prefix",
    shape = "Filename prefix"
  )

GSE224463_stationary_rpoS_plot <- ggplot(
  GSE224463_stationary_library_QC,
  aes(
    x = .data$time_min,
    y = .data$rpoS_counts_per_million,
    colour = .data$JH_prefix,
    group = interaction(
      .data$JH_prefix,
      .data$author_replicate
    )
  )
) +
  geom_line(
    linewidth = 0.45,
    alpha = 0.75
  ) +
  geom_point(
    size = 1.9
  ) +
  facet_wrap(
    vars(series),
    nrow = 1,
    scales = "free_x"
  ) +
  theme_classic() +
  theme(
    legend.position = "top"
  ) +
  labs(
    title = expression(italic(rpoS)~"count audit: JH01/JH02 mapping evidence"),
    x = "Author time point (min)",
    y = expression(italic(rpoS)~"raw counts per million"),
    colour = "Filename prefix"
  )

GSE224463_stationary_raw_PCA_plot <- ggplot(
  GSE224463_stationary_raw_PCA_coordinates,
  aes(
    x = .data$PC1,
    y = .data$PC2,
    colour = .data$JH_prefix,
    shape = .data$series
  )
) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.25,
    colour = "grey70"
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.25,
    colour = "grey70"
  ) +
  geom_point(
    size = 2.4
  ) +
  geom_text(
    aes(
      label = paste0(
        .data$author_replicate,
        "_",
        .data$time_min
      )
    ),
    size = 2.5,
    vjust = -0.7,
    show.legend = FALSE
  ) +
  theme_classic() +
  theme(
    legend.position = "top"
  ) +
  labs(
    title = "GSE224463 stationary-phase raw-count PCA",
    subtitle = "Descriptive QC only; not batch-corrected and not used for DESeq2",
    x = paste0(
      "PC1: ",
      format(round(GSE224463_stationary_raw_PCA_variance[1], 1), nsmall = 1),
      "% variance"
    ),
    y = paste0(
      "PC2: ",
      format(round(GSE224463_stationary_raw_PCA_variance[2], 1), nsmall = 1),
      "% variance"
    ),
    colour = "Filename prefix",
    shape = "Series"
  )

print(GSE224463_stationary_library_plot)
print(GSE224463_stationary_rpoS_plot)
print(GSE224463_stationary_raw_PCA_plot)

ggsave(
  filename = file.path(
    GSE224463_audit_dir,
    "GSE224463_stationary_phase_raw_library_QC.tiff"
  ),
  plot = GSE224463_stationary_library_plot,
  width = 12,
  height = 6,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  filename = file.path(
    GSE224463_audit_dir,
    "GSE224463_stationary_phase_rpoS_JH_prefix_audit.tiff"
  ),
  plot = GSE224463_stationary_rpoS_plot,
  width = 10,
  height = 5,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  filename = file.path(
    GSE224463_audit_dir,
    "GSE224463_stationary_phase_raw_log_count_PCA.tiff"
  ),
  plot = GSE224463_stationary_raw_PCA_plot,
  width = 10,
  height = 7,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

GSE224463_WT_stationary_count_file_candidates <- GSE224463_WT_stationary_GEO %>%
  dplyr::select(
    GSM,
    title,
    series,
    author_replicate,
    time_min,
    condition
  ) %>%
  dplyr::left_join(
    GSE224463_stationary_library_QC %>%
      dplyr::select(
        count_file,
        relative_count_file,
        file_name,
        JH_prefix,
        series,
        author_replicate,
        time_min,
        raw_library_size,
        genes_detected,
        percent_zero_counts,
        rpoS_raw_count,
        rpoS_counts_per_million
      ),
    by = c("series", "author_replicate", "time_min")
  ) %>%
  dplyr::mutate(
    candidate_mapping_status =
      "UNRESOLVED: verify whether this JH prefix is wild type before DESeq2",
    intended_analysis = dplyr::case_when(
      .data$series == "short" & .data$time_min %in% c(0L, 150L) ~
        "primary_short_150min_vs_0min",
      .data$series == "long" & .data$time_min %in% c(0L, 300L) ~
        "confirmatory_long_300min_vs_0min",
      TRUE ~ "time_course_context_only"
    )
  ) %>%
  dplyr::arrange(
    .data$series,
    .data$author_replicate,
    .data$time_min,
    .data$JH_prefix
  )

if (nrow(GSE224463_WT_stationary_count_file_candidates) != 72L) {
  stop(
    "Each of the 36 GEO wild-type samples should have two candidate JH files. ",
    "Expected 72 candidate rows but found ",
    nrow(GSE224463_WT_stationary_count_file_candidates),
    ".",
    call. = FALSE
  )
}

readr::write_csv(
  GSE224463_stationary_file_inventory,
  file.path(
    GSE224463_audit_dir,
    "GSE224463_stationary_phase_count_file_inventory.csv"
  ),
  na = ""
)

readr::write_csv(
  GSE224463_WT_stationary_count_file_candidates,
  file.path(
    GSE224463_audit_dir,
    "GSE224463_WT_stationary_phase_count_file_candidate_mapping.csv"
  ),
  na = ""
)

GSE224463_author_text_files <- list.files(
  GSE224463_project_dir,
  pattern = "\\.(R|Rmd|md|txt|csv|tsv|yaml|yml)$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

GSE224463_author_text_files <- GSE224463_author_text_files[
  !stringr::str_detect(
    basename(GSE224463_author_text_files),
    "_counts\\.txt$"
  )
]

GSE224463_author_text_file_info <- file.info(GSE224463_author_text_files)

GSE224463_author_text_files <- GSE224463_author_text_files[
  !is.na(GSE224463_author_text_file_info$size) &
    GSE224463_author_text_file_info$size <= 5 * 1024^2
]

GSE224463_author_text_file_inventory <- tibble::tibble(
  text_file = normalizePath(
    GSE224463_author_text_files,
    winslash = "/",
    mustWork = TRUE
  ),
  bytes = as.numeric(file.info(GSE224463_author_text_files)$size)
) %>%
  dplyr::mutate(
    relative_text_file = sub(
      paste0(GSE224463_project_dir_normalized, "/"),
      "",
      .data$text_file,
      fixed = TRUE
    )
  ) %>%
  dplyr::arrange(.data$relative_text_file)

author_mapping_keywords <- c(
  "JH01",
  "JH02",
  "LSP",
  "stationary phase",
  "delta-rpoS",
  "wild-type",
  "SP_A",
  "SP_B",
  "SP_C"
)

GSE224463_author_mapping_text_hits <- dplyr::bind_rows(
  lapply(
    seq_len(nrow(GSE224463_author_text_file_inventory)),
    function(i) {
      one_file <- GSE224463_author_text_file_inventory$text_file[i]
      one_relative_file <-
        GSE224463_author_text_file_inventory$relative_text_file[i]

      one_text <- tryCatch(
        readLines(
          one_file,
          warn = FALSE,
          encoding = "UTF-8",
          n = 50000L
        ),
        error = function(e) character(0)
      )

      if (length(one_text) == 0L) {
        return(tibble::tibble())
      }

      dplyr::bind_rows(
        lapply(
          author_mapping_keywords,
          function(one_keyword) {
            hit_lines <- which(
              stringr::str_detect(
                one_text,
                stringr::fixed(one_keyword, ignore_case = TRUE)
              )
            )

            if (length(hit_lines) == 0L) {
              return(tibble::tibble())
            }

            tibble::tibble(
              relative_text_file = one_relative_file,
              keyword = one_keyword,
              line_number = hit_lines,
              line_text = trimws(one_text[hit_lines])
            )
          }
        )
      )
    }
  )
) %>%
  dplyr::arrange(
    .data$relative_text_file,
    .data$line_number,
    .data$keyword
  )

readr::write_csv(
  GSE224463_author_text_file_inventory,
  file.path(
    GSE224463_audit_dir,
    "GSE224463_author_text_file_inventory.csv"
  ),
  na = ""
)

readr::write_tsv(
  GSE224463_author_mapping_text_hits,
  file.path(
    GSE224463_audit_dir,
    "GSE224463_author_source_JH_mapping_text_hits.tsv"
  ),
  na = ""
)

GSE224463_stationary_audit_objects <- list(
  WT_stationary_GEO_manifest = GSE224463_WT_stationary_GEO,
  stationary_count_file_inventory = GSE224463_stationary_file_inventory,
  stationary_gene_set_audit = GSE224463_stationary_gene_set_audit,
  stationary_library_QC = GSE224463_stationary_library_QC,
  rpoS_prefix_audit = GSE224463_rpoS_prefix_audit,
  stationary_log_count_correlation = GSE224463_stationary_log_count_correlation,
  stationary_log_count_correlation_long =
    GSE224463_stationary_log_count_correlation_long,
  stationary_raw_PCA_coordinates = GSE224463_stationary_raw_PCA_coordinates,
  WT_stationary_count_file_candidates =
    GSE224463_WT_stationary_count_file_candidates,
  author_mapping_text_hits = GSE224463_author_mapping_text_hits,
  source_count_matrix = GSE224463_stationary_count_matrix,
  source_gene_ids = reference_gene_ids,
  rpoS_locus_tag = rpoS_locus_tag,
  audit_directory = GSE224463_audit_dir
)

saveRDS(
  GSE224463_stationary_audit_objects,
  file = file.path(
    GSE224463_audit_dir,
    "GSE224463_stationary_phase_count_audit_objects.rds"
  )
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(
    GSE224463_audit_dir,
    "GSE224463_stationary_phase_count_audit_sessionInfo.txt"
  )
)

message(
  "GSE224463 stationary-phase raw-count audit completed.\n\n",
  "Important: no JH prefix was assigned to wild type by this script.\n",
  "Review these first:\n",
  "  1. GSE224463_stationary_phase_rpoS_JH_prefix_audit.csv\n",
  "  2. GSE224463_author_source_JH_mapping_text_hits.tsv\n",
  "  3. GSE224463_WT_stationary_phase_count_file_candidate_mapping.csv\n\n",
  "After the JH01/JH02 genotype mapping is verified, the next script will write ",
  "hard-locked short- and long-series manifests and then run the separate ",
  "DESeq2/thermodynamic analyses.\n\n",
  "Audit output directory:\n",
  GSE224463_audit_dir
)

GSE224463_B_BOOT <- 5000L
GSE224463_B_PERM <- 5000L
GSE224463_resampling_seed <- 20260807L

GSE224463_extreme_depth_threshold <- 100000L
GSE224463_low_depth_threshold <- 1000000L
GSE224463_extreme_detected_gene_threshold <- 1000L
GSE224463_low_detected_gene_threshold <- 3000L

GSE224463_audit_dir <- file.path(
  GSE224463_project_dir,
  "analysis_manifests",
  "00_stationary_phase_count_audit"
)

GSE224463_audit_rds_file <- file.path(
  GSE224463_audit_dir,
  "GSE224463_stationary_phase_count_audit_objects.rds"
)

GSE224463_manifest_dir <- file.path(
  GSE224463_project_dir,
  "analysis_manifests",
  "01_WT_stationary_phase_time_course_verified"
)

GSE224463_output_dir <- file.path(
  GSE224463_project_dir,
  "analysis_results",
  "01_WT_stationary_phase_time_course_thermodynamic_analysis"
)

GSE224463_table_dir <- file.path(GSE224463_output_dir, "tables")
GSE224463_DE_table_dir <- file.path(GSE224463_table_dir, "DESeq2_by_contrast")
GSE224463_plot_dir <- file.path(GSE224463_output_dir, "plots")
GSE224463_qc_dir <- file.path(GSE224463_output_dir, "quality_control")

for (one_dir in c(
  GSE224463_manifest_dir,
  GSE224463_table_dir,
  GSE224463_DE_table_dir,
  GSE224463_plot_dir,
  GSE224463_qc_dir
)) {
  dir.create(one_dir, recursive = TRUE, showWarnings = FALSE)
}

if (!dir.exists(GSE224463_project_dir)) {
  stop(
    "GSE224463_project_dir does not exist:\n",
    GSE224463_project_dir,
    "\nUpdate Section 0 before continuing.",
    call. = FALSE
  )
}

if (!file.exists(GSE224463_audit_rds_file)) {
  stop(
    "The completed raw-count audit RDS is missing:\n",
    GSE224463_audit_rds_file,
    "\nRun 00_GSE224463_audit_stationary_phase_counts_and_candidate_manifest.R ",
    "first.",
    call. = FALSE
  )
}

if (!file.exists(GSE224463_prediction_file)) {
  stop(
    "The E. coli ΔGfold prediction file is missing:\n",
    GSE224463_prediction_file,
    call. = FALSE
  )
}

GSE224463_audit_objects <- readRDS(GSE224463_audit_rds_file)

GSE224463_required_audit_objects <- c(
  "WT_stationary_GEO_manifest",
  "stationary_library_QC",
  "rpoS_prefix_audit",
  "source_count_matrix"
)

GSE224463_missing_audit_objects <- setdiff(
  GSE224463_required_audit_objects,
  names(GSE224463_audit_objects)
)

if (length(GSE224463_missing_audit_objects) > 0L) {
  stop(
    "The audit RDS lacks required object(s): ",
    paste(GSE224463_missing_audit_objects, collapse = ", "),
    call. = FALSE
  )
}

GSE224463_confirmed_WT_prefix <- "JH01"

GSE224463_rpoS_prefix_audit <- GSE224463_audit_objects$rpoS_prefix_audit %>%
  tibble::as_tibble() %>%
  dplyr::mutate(JH_prefix = as.character(.data$JH_prefix))

if (!GSE224463_confirmed_WT_prefix %in%
    GSE224463_rpoS_prefix_audit$JH_prefix) {
  stop(
    "The audit RDS does not contain the confirmed wild-type prefix ",
    GSE224463_confirmed_WT_prefix,
    ".",
    call. = FALSE
  )
}

GSE224463_raw_count_matrix <- as.matrix(
  GSE224463_audit_objects$source_count_matrix
)

if (is.null(rownames(GSE224463_raw_count_matrix)) ||
    is.null(colnames(GSE224463_raw_count_matrix))) {
  stop(
    "The audited source count matrix must have gene IDs as row names and ",
    "author count-file names as column names.",
    call. = FALSE
  )
}

storage.mode(GSE224463_raw_count_matrix) <- "numeric"

if (any(!is.finite(GSE224463_raw_count_matrix)) ||
    any(GSE224463_raw_count_matrix < 0) ||
    any(abs(
      GSE224463_raw_count_matrix - round(GSE224463_raw_count_matrix)
    ) > 1e-8)) {
  stop(
    "The audited matrix contains values that are not non-negative integer raw ",
    "counts.",
    call. = FALSE
  )
}

rownames(GSE224463_raw_count_matrix) <- tolower(
  trimws(as.character(rownames(GSE224463_raw_count_matrix)))
)

if (anyNA(rownames(GSE224463_raw_count_matrix)) ||
    any(rownames(GSE224463_raw_count_matrix) == "") ||
    anyDuplicated(rownames(GSE224463_raw_count_matrix)) > 0L) {
  stop(
    "The audited count matrix has missing or duplicated gene IDs after ",
    "case-standardisation.",
    call. = FALSE
  )
}

storage.mode(GSE224463_raw_count_matrix) <- "integer"

GSE224463_stationary_library_QC <- GSE224463_audit_objects$stationary_library_QC %>%
  tibble::as_tibble() %>%
  dplyr::mutate(
    file_name = as.character(.data$file_name),
    JH_prefix = as.character(.data$JH_prefix),
    series = as.character(.data$series),
    author_replicate = as.character(.data$author_replicate),
    time_min = as.integer(.data$time_min),
    raw_library_size = as.numeric(.data$raw_library_size),
    genes_detected = as.integer(.data$genes_detected),
    percent_zero_counts = as.numeric(.data$percent_zero_counts)
  )

if (!setequal(
  colnames(GSE224463_raw_count_matrix),
  GSE224463_stationary_library_QC$file_name
)) {
  stop(
    "The audited count-matrix columns and library-QC file names do not match.",
    call. = FALSE
  )
}

GSE224463_stationary_library_QC <- GSE224463_stationary_library_QC %>%
  dplyr::arrange(
    match(.data$file_name, colnames(GSE224463_raw_count_matrix))
  )

GSE224463_WT_GEO_manifest <- GSE224463_audit_objects$WT_stationary_GEO_manifest %>%
  tibble::as_tibble() %>%
  dplyr::transmute(
    GSM = as.character(.data$GSM),
    GEO_title = as.character(.data$title),
    series = as.character(.data$series),
    author_replicate = as.character(.data$author_replicate),
    time_min = as.integer(.data$time_min)
  )

GSE224463_WT_manifest <- GSE224463_stationary_library_QC %>%
  dplyr::filter(.data$JH_prefix == GSE224463_confirmed_WT_prefix) %>%
  dplyr::left_join(
    GSE224463_WT_GEO_manifest,
    by = c("series", "author_replicate", "time_min")
  ) %>%
  dplyr::mutate(
    genotype = "wild type",
    filename_prefix_mapping =
      "JH01 = wild type; verified from author source code and audit",
    reference_or_target = dplyr::if_else(
      .data$time_min == 0L,
      "Exponential reference",
      "Nutrient-depletion time-course sample"
    ),
    condition_label = dplyr::if_else(
      .data$time_min == 0L,
      "Exponential phase, 0 min",
      paste0("Nutrient depletion, ", .data$time_min, " min")
    ),
    library_QC_flag = dplyr::case_when(
      .data$raw_library_size < GSE224463_extreme_depth_threshold |
        .data$genes_detected < GSE224463_extreme_detected_gene_threshold ~
        "Extreme shallow library: retained; interpret results cautiously",
      .data$raw_library_size < GSE224463_low_depth_threshold |
        .data$genes_detected < GSE224463_low_detected_gene_threshold ~
        "Lower-depth library: retained; inspect contrast-level QC",
      TRUE ~
        "No automatic-exclusion QC flag"
    ),
    sample_order_within_series = dplyr::case_when(
      .data$series == "short" ~ match(
        .data$time_min,
        c(0L, 30L, 60L, 90L, 120L, 150L)
      ),
      .data$series == "long" ~ match(
        .data$time_min,
        c(0L, 90L, 120L, 180L, 240L, 300L)
      )
    )
  ) %>%
  dplyr::arrange(
    factor(.data$series, levels = c("short", "long")),
    .data$author_replicate,
    .data$sample_order_within_series
  )

if (nrow(GSE224463_WT_manifest) != 36L ||
    anyNA(GSE224463_WT_manifest$GSM) ||
    anyNA(GSE224463_WT_manifest$GEO_title)) {
  stop(
    "The hard-locked WT manifest must contain 36 mapped GEO/count-file rows.",
    call. = FALSE
  )
}

GSE224463_expected_WT_geometry <- dplyr::bind_rows(
  tidyr::expand_grid(
    series = "short",
    author_replicate = c("SP_A", "SP_B", "SP_C"),
    time_min = c(0L, 30L, 60L, 90L, 120L, 150L)
  ),
  tidyr::expand_grid(
    series = "long",
    author_replicate = c("LSP_A", "LSP_C", "LSP_E"),
    time_min = c(0L, 90L, 120L, 180L, 240L, 300L)
  )
)

GSE224463_observed_WT_geometry <- GSE224463_WT_manifest %>%
  dplyr::select(.data$series, .data$author_replicate, .data$time_min)

GSE224463_missing_WT_geometry <- dplyr::anti_join(
  GSE224463_expected_WT_geometry,
  GSE224463_observed_WT_geometry,
  by = c("series", "author_replicate", "time_min")
)

GSE224463_unexpected_WT_geometry <- dplyr::anti_join(
  GSE224463_observed_WT_geometry,
  GSE224463_expected_WT_geometry,
  by = c("series", "author_replicate", "time_min")
)

if (nrow(GSE224463_missing_WT_geometry) > 0L ||
    nrow(GSE224463_unexpected_WT_geometry) > 0L) {
  stop(
    "The hard-locked WT sample geometry does not match the author design.",
    call. = FALSE
  )
}

GSE224463_contrast_manifest <- dplyr::bind_rows(
  tibble::tibble(
    series = "short",
    target_time_min = c(30L, 60L, 90L, 120L, 150L)
  ),
  tibble::tibble(
    series = "long",
    target_time_min = c(90L, 120L, 180L, 240L, 300L)
  )
) %>%
  dplyr::mutate(
    reference_time_min = 0L,
    analysis_id = paste0(
      "GSE224463_",
      .data$series,
      "_",
      .data$target_time_min,
      "min_vs_0min"
    ),
    comparison_label = paste0(
      dplyr::if_else(.data$series == "short", "Short", "Long"),
      ": ",
      .data$target_time_min,
      " min minus 0 min"
    ),
    target_label = paste0(.data$target_time_min, " min"),
    reference_label = "0 min exponential phase",
    biological_context = dplyr::case_when(
      .data$series == "short" ~
        "Nutrient-depletion / entry-into-stationary-phase time course",
      .data$series == "long" & .data$target_time_min %in% c(90L, 120L) ~
        "Extended nutrient-depletion time course before the late stationary-phase window",
      .data$series == "long" & .data$target_time_min == 180L ~
        "Late nutrient-depletion transition after the median RpoS-regulated gene-onset time",
      .data$series == "long" & .data$target_time_min == 240L ~
        "Late stationary-phase time point near the reported RpoS peak",
      .data$series == "long" & .data$target_time_min == 300L ~
        "Latest sampled time point; retained with its raw-depth QC flag"
    ),
    comparison_order = dplyr::row_number()
  )

GSE224463_reference_QC_by_series <- GSE224463_WT_manifest %>%
  dplyr::filter(.data$time_min == 0L) %>%
  dplyr::group_by(.data$series) %>%
  dplyr::summarise(
    N_reference_replicates = dplyr::n(),
    reference_min_raw_library_size = min(.data$raw_library_size),
    reference_min_genes_detected = min(.data$genes_detected),
    .groups = "drop"
  )

GSE224463_target_QC_by_time <- GSE224463_WT_manifest %>%
  dplyr::filter(.data$time_min > 0L) %>%
  dplyr::group_by(.data$series, target_time_min = .data$time_min) %>%
  dplyr::summarise(
    N_target_replicates = dplyr::n(),
    target_min_raw_library_size = min(.data$raw_library_size),
    target_median_raw_library_size = median(.data$raw_library_size),
    target_min_genes_detected = min(.data$genes_detected),
    target_max_percent_zero_counts = max(.data$percent_zero_counts),
    target_library_QC_flags = paste(
      unique(.data$library_QC_flag),
      collapse = " | "
    ),
    .groups = "drop"
  )

GSE224463_contrast_manifest <- GSE224463_contrast_manifest %>%
  dplyr::left_join(
    GSE224463_reference_QC_by_series,
    by = "series"
  ) %>%
  dplyr::left_join(
    GSE224463_target_QC_by_time,
    by = c("series", "target_time_min")
  ) %>%
  dplyr::mutate(
    contrast_QC_status = dplyr::case_when(
      .data$target_min_raw_library_size < GSE224463_extreme_depth_threshold |
        .data$target_min_genes_detected <
          GSE224463_extreme_detected_gene_threshold ~
        "Extreme shallow target library present: results retained but require caution",
      .data$target_min_raw_library_size < GSE224463_low_depth_threshold |
        .data$target_min_genes_detected <
          GSE224463_low_detected_gene_threshold ~
        "Lower-depth target library present: results retained with QC note",
      TRUE ~
        "No automatic-exclusion QC flag"
    )
  ) %>%
  dplyr::arrange(.data$comparison_order)

if (nrow(GSE224463_contrast_manifest) != 10L ||
    any(GSE224463_contrast_manifest$N_reference_replicates != 3L) ||
    any(GSE224463_contrast_manifest$N_target_replicates != 3L)) {
  stop(
    "Each of the ten time-point contrasts must have three reference and three ",
    "target biological replicates.",
    call. = FALSE
  )
}

print(GSE224463_WT_manifest)
print(GSE224463_contrast_manifest)

readr::write_tsv(
  GSE224463_WT_manifest,
  file.path(
    GSE224463_manifest_dir,
    "GSE224463_WT_stationary_phase_time_course_manifest_VERIFIED.tsv"
  )
)

readr::write_tsv(
  GSE224463_contrast_manifest,
  file.path(
    GSE224463_manifest_dir,
    "GSE224463_WT_stationary_phase_time_course_contrasts_VERIFIED.tsv"
  )
)

saveRDS(
  list(
    WT_JH_prefix = GSE224463_confirmed_WT_prefix,
    WT_manifest = GSE224463_WT_manifest,
    contrast_manifest = GSE224463_contrast_manifest,
    rpoS_prefix_audit = GSE224463_rpoS_prefix_audit,
    source_audit_rds = GSE224463_audit_rds_file
  ),
  file.path(
    GSE224463_manifest_dir,
    "GSE224463_WT_stationary_phase_time_course_manifest_objects.rds"
  )
)

GSE224463_prediction_raw <- readr::read_csv(
  GSE224463_prediction_file,
  show_col_types = FALSE,
  progress = FALSE,
  name_repair = "unique"
)

GSE224463_required_prediction_columns <- c(
  "Gene_id",
  "Normalized_MFE",
  "Length"
)

GSE224463_missing_prediction_columns <- setdiff(
  GSE224463_required_prediction_columns,
  colnames(GSE224463_prediction_raw)
)

if (length(GSE224463_missing_prediction_columns) > 0L) {
  stop(
    "The prediction table lacks required column(s): ",
    paste(GSE224463_missing_prediction_columns, collapse = ", "),
    ". The Length column is required for the full-transcriptome count/length ",
    "abundance weighting.",
    call. = FALSE
  )
}

GSE224463_prediction <- GSE224463_prediction_raw %>%
  dplyr::transmute(
    Gene_id = tolower(trimws(as.character(.data$Gene_id))),
    transcript_length_nt = as.numeric(.data$Length),
    Normalized_MFE = as.numeric(.data$Normalized_MFE)
  ) %>%
  dplyr::filter(!is.na(.data$Gene_id), .data$Gene_id != "") %>%
  dplyr::distinct(.data$Gene_id, .keep_all = TRUE)

if (anyDuplicated(GSE224463_prediction$Gene_id) > 0L) {
  stop(
    "The prediction table contains duplicated Gene_id values after ",
    "case-standardisation.",
    call. = FALSE
  )
}

GSE224463_thermodynamic_annotation_audit <- tibble::tibble(
  Gene_id = rownames(GSE224463_raw_count_matrix)
) %>%
  dplyr::left_join(
    GSE224463_prediction,
    by = "Gene_id"
  ) %>%
  dplyr::mutate(
    finite_Normalized_MFE = is.finite(.data$Normalized_MFE),
    positive_transcript_length = is.finite(.data$transcript_length_nt) &
      .data$transcript_length_nt > 0,
    included_in_DE_descriptor_analysis = .data$finite_Normalized_MFE,
    included_in_full_transcriptome_weighted_analysis =
      .data$finite_Normalized_MFE &
      .data$positive_transcript_length,
    exclusion_reason_weighted_analysis = dplyr::case_when(
      .data$included_in_full_transcriptome_weighted_analysis ~ "",
      !.data$positive_transcript_length &
        !.data$finite_Normalized_MFE ~
        "Invalid transcript length; missing/non-finite Normalized_MFE",
      !.data$positive_transcript_length ~ "Invalid or missing transcript length",
      !.data$finite_Normalized_MFE ~ "Missing/non-finite Normalized_MFE"
    )
  )

GSE224463_DE_descriptor_annotation <-
  GSE224463_thermodynamic_annotation_audit %>%
  dplyr::filter(.data$included_in_DE_descriptor_analysis)

GSE224463_weighted_annotation <-
  GSE224463_thermodynamic_annotation_audit %>%
  dplyr::filter(.data$included_in_full_transcriptome_weighted_analysis) %>%
  dplyr::mutate(
    raw_matrix_row = match(
      .data$Gene_id,
      rownames(GSE224463_raw_count_matrix)
    )
  ) %>%
  dplyr::arrange(.data$raw_matrix_row)

if (nrow(GSE224463_DE_descriptor_annotation) < 3500L) {
  stop(
    "Too few genes have finite ΔGfold predictions for the DE descriptor ",
    "analyses.",
    call. = FALSE
  )
}

if (nrow(GSE224463_weighted_annotation) < 3500L) {
  stop(
    "Too few genes have a positive transcript length and finite ΔGfold ",
    "values for the full-transcriptome weighted analyses.",
    call. = FALSE
  )
}

readr::write_tsv(
  GSE224463_thermodynamic_annotation_audit,
  file.path(
    GSE224463_table_dir,
    "GSE224463_thermodynamic_annotation_audit.tsv"
  )
)

GSE224463_WT_library_QC_plot <- ggplot(
  GSE224463_WT_manifest,
  aes(
    x = .data$time_min,
    y = .data$raw_library_size,
    group = .data$author_replicate,
    colour = .data$library_QC_flag,
    shape = .data$author_replicate
  )
) +
  geom_line(colour = "grey65", linewidth = 0.4) +
  geom_point(size = 2.8) +
  facet_wrap(~ series, nrow = 1, scales = "free_x") +
  scale_y_log10() +
  scale_colour_manual(
    values = c(
      "No automatic-exclusion QC flag" = "#2066a8",
      "Lower-depth library: retained; inspect contrast-level QC" = "#e69f00",
      "Extreme shallow library: retained; interpret results cautiously" = "#d55e00"
    )
  ) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank()
  ) +
  labs(
    title = "GSE224463 WT stationary-phase time course: raw-library QC",
    subtitle = "All author-provided JH01 samples are retained",
    x = "Time after exponential-phase baseline (min)",
    y = "Raw library size (log10 scale)",
    shape = "Author replicate"
  )

print(GSE224463_WT_library_QC_plot)

ggsave(
  filename = file.path(
    GSE224463_qc_dir,
    "GSE224463_WT_stationary_phase_time_course_library_QC.tiff"
  ),
  plot = GSE224463_WT_library_QC_plot,
  width = 11,
  height = 5,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

GSE224463_series_PCA_data_list <- vector("list", 2L)
GSE224463_series_correlation_list <- vector("list", 2L)
GSE224463_series_PCA_plot_list <- vector("list", 2L)
GSE224463_series_correlation_plot_list <- vector("list", 2L)

GSE224463_series_order <- c("short", "long")

for (GSE224463_series_index in seq_along(GSE224463_series_order)) {
  GSE224463_one_series <- GSE224463_series_order[GSE224463_series_index]

  GSE224463_one_series_manifest <- GSE224463_WT_manifest %>%
    dplyr::filter(.data$series == GSE224463_one_series) %>%
    dplyr::arrange(.data$author_replicate, .data$sample_order_within_series)

  GSE224463_one_series_counts <- GSE224463_raw_count_matrix[
    ,
    GSE224463_one_series_manifest$file_name,
    drop = FALSE
  ]

  GSE224463_one_series_nonzero_genes <-
    rowSums(GSE224463_one_series_counts) > 0L

  if (sum(GSE224463_one_series_nonzero_genes) < 1000L) {
    stop(
      "Too few non-zero genes are available for the ",
      GSE224463_one_series,
      "-series VST QC.",
      call. = FALSE
    )
  }

  GSE224463_one_series_counts <- GSE224463_one_series_counts[
    GSE224463_one_series_nonzero_genes,
    ,
    drop = FALSE
  ]

  GSE224463_one_series_coldata <- GSE224463_one_series_manifest %>%
    dplyr::mutate(
      author_replicate = factor(.data$author_replicate),
      timepoint = factor(
        paste0("T", .data$time_min, "min"),
        levels = paste0(
          "T",
          sort(unique(.data$time_min)),
          "min"
        )
      )
    ) %>%
    dplyr::select(
      .data$file_name,
      .data$author_replicate,
      .data$timepoint
    ) %>%
    as.data.frame()

  rownames(GSE224463_one_series_coldata) <-
    GSE224463_one_series_coldata$file_name
  GSE224463_one_series_coldata$file_name <- NULL

  if (!identical(
    colnames(GSE224463_one_series_counts),
    rownames(GSE224463_one_series_coldata)
  )) {
    stop(
      "VST-QC count-matrix columns and colData rows are not aligned for the ",
      GSE224463_one_series,
      " series.",
      call. = FALSE
    )
  }

  GSE224463_one_series_dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = GSE224463_one_series_counts,
    colData = GSE224463_one_series_coldata,
    design = ~ author_replicate + timepoint
  )

  GSE224463_one_series_dds <- DESeq2::estimateSizeFactors(
    GSE224463_one_series_dds
  )

  GSE224463_one_series_vst <- DESeq2::vst(
    GSE224463_one_series_dds,
    blind = TRUE
  )

  GSE224463_one_series_vst_matrix <- SummarizedExperiment::assay(
    GSE224463_one_series_vst
  )

  GSE224463_one_series_PCA <- stats::prcomp(
    t(GSE224463_one_series_vst_matrix),
    scale. = FALSE
  )

  GSE224463_one_series_PCA_percent_variance <-
    100 * GSE224463_one_series_PCA$sdev^2 /
    sum(GSE224463_one_series_PCA$sdev^2)

  GSE224463_one_series_PCA_data <- tibble::tibble(
    file_name = rownames(GSE224463_one_series_PCA$x),
    PC1 = GSE224463_one_series_PCA$x[, 1L],
    PC2 = GSE224463_one_series_PCA$x[, 2L],
    series = GSE224463_one_series,
    PC1_percent_variance = GSE224463_one_series_PCA_percent_variance[1L],
    PC2_percent_variance = GSE224463_one_series_PCA_percent_variance[2L]
  ) %>%
    dplyr::left_join(
      GSE224463_one_series_manifest %>%
        dplyr::select(
          .data$file_name,
          .data$author_replicate,
          .data$time_min,
          .data$raw_library_size,
          .data$library_QC_flag
        ),
      by = "file_name"
    )

  GSE224463_one_series_correlation_matrix <- stats::cor(
    GSE224463_one_series_vst_matrix,
    method = "pearson"
  )

  GSE224463_one_series_correlation_long <- as.data.frame(
    as.table(GSE224463_one_series_correlation_matrix),
    stringsAsFactors = FALSE
  ) %>%
    tibble::as_tibble() %>%
    dplyr::transmute(
      file_name_1 = as.character(.data$Var1),
      file_name_2 = as.character(.data$Var2),
      Pearson_r = as.numeric(.data$Freq),
      series = GSE224463_one_series
    ) %>%
    dplyr::left_join(
      GSE224463_one_series_manifest %>%
        dplyr::transmute(
          file_name_1 = .data$file_name,
          time_min_1 = .data$time_min,
          author_replicate_1 = .data$author_replicate
        ),
      by = "file_name_1"
    ) %>%
    dplyr::left_join(
      GSE224463_one_series_manifest %>%
        dplyr::transmute(
          file_name_2 = .data$file_name,
          time_min_2 = .data$time_min,
          author_replicate_2 = .data$author_replicate
        ),
      by = "file_name_2"
    ) %>%
    dplyr::mutate(
      file_name_1 = factor(
        .data$file_name_1,
        levels = GSE224463_one_series_manifest$file_name
      ),
      file_name_2 = factor(
        .data$file_name_2,
        levels = rev(GSE224463_one_series_manifest$file_name)
      )
    )

  GSE224463_one_series_PCA_plot <- ggplot(
    GSE224463_one_series_PCA_data,
    aes(
      x = .data$PC1,
      y = .data$PC2,
      colour = factor(.data$time_min),
      shape = .data$author_replicate
    )
  ) +
    geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.3) +
    geom_point(size = 3) +
    geom_text(
      aes(label = paste0(.data$author_replicate, "_", .data$time_min)),
      vjust = -0.8,
      size = 2.7,
      show.legend = FALSE
    ) +
    scale_colour_brewer(palette = "Dark2") +
    theme_classic() +
    theme(legend.position = "bottom") +
    labs(
      title = paste0(
        ifelse(
          GSE224463_one_series == "short",
          "Short series",
          "Long series"
        ),
        ": VST sample PCA"
      ),
      x = paste0(
        "PC1: ",
        format(
          round(GSE224463_one_series_PCA_percent_variance[1L], 1),
          nsmall = 1
        ),
        "% variance"
      ),
      y = paste0(
        "PC2: ",
        format(
          round(GSE224463_one_series_PCA_percent_variance[2L], 1),
          nsmall = 1
        ),
        "% variance"
      ),
      colour = "Time (min)",
      shape = "Author replicate"
    )

  GSE224463_one_series_correlation_plot <- ggplot(
    GSE224463_one_series_correlation_long,
    aes(
      x = .data$file_name_1,
      y = .data$file_name_2,
      fill = .data$Pearson_r
    )
  ) +
    geom_tile(colour = "white", linewidth = 0.2) +
    geom_text(
      aes(label = sprintf("%.2f", .data$Pearson_r)),
      size = 2.2
    ) +
    scale_fill_gradient(
      low = "white",
      high = "#2066a8",
      limits = c(0, 1)
    ) +
    theme_classic() +
    theme(
      axis.title = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "bottom"
    ) +
    labs(
      title = paste0(
        ifelse(
          GSE224463_one_series == "short",
          "Short series",
          "Long series"
        ),
        ": VST sample correlation"
      ),
      fill = "Pearson r"
    )

  GSE224463_series_PCA_data_list[[GSE224463_series_index]] <-
    GSE224463_one_series_PCA_data
  GSE224463_series_correlation_list[[GSE224463_series_index]] <-
    GSE224463_one_series_correlation_long
  GSE224463_series_PCA_plot_list[[GSE224463_series_index]] <-
    GSE224463_one_series_PCA_plot
  GSE224463_series_correlation_plot_list[[GSE224463_series_index]] <-
    GSE224463_one_series_correlation_plot
}

GSE224463_series_PCA_data <- dplyr::bind_rows(
  GSE224463_series_PCA_data_list
)

GSE224463_series_correlation_long <- dplyr::bind_rows(
  GSE224463_series_correlation_list
)

GSE224463_VST_PCA_panel <- GSE224463_series_PCA_plot_list[[1L]] +
  GSE224463_series_PCA_plot_list[[2L]] +
  patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

GSE224463_VST_correlation_panel <-
  GSE224463_series_correlation_plot_list[[1L]] +
  GSE224463_series_correlation_plot_list[[2L]] +
  patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

print(GSE224463_VST_PCA_panel)
print(GSE224463_VST_correlation_panel)

readr::write_tsv(
  GSE224463_series_PCA_data,
  file.path(
    GSE224463_qc_dir,
    "GSE224463_WT_stationary_phase_time_course_VST_PCA_coordinates.tsv"
  )
)

readr::write_tsv(
  GSE224463_series_correlation_long,
  file.path(
    GSE224463_qc_dir,
    "GSE224463_WT_stationary_phase_time_course_VST_sample_correlations.tsv"
  )
)

ggsave(
  filename = file.path(
    GSE224463_qc_dir,
    "GSE224463_WT_stationary_phase_time_course_VST_PCA.tiff"
  ),
  plot = GSE224463_VST_PCA_panel,
  width = 12,
  height = 6,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  filename = file.path(
    GSE224463_qc_dir,
    "GSE224463_WT_stationary_phase_time_course_VST_sample_correlations.tiff"
  ),
  plot = GSE224463_VST_correlation_panel,
  width = 11,
  height = 5,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

GSE224463_DE_all_list <- vector(
  "list",
  nrow(GSE224463_contrast_manifest)
)

GSE224463_DE_thermodynamic_list <- vector(
  "list",
  nrow(GSE224463_contrast_manifest)
)

GSE224463_significant_DE_thermodynamic_list <- vector(
  "list",
  nrow(GSE224463_contrast_manifest)
)

GSE224463_DE_summary_list <- vector(
  "list",
  nrow(GSE224463_contrast_manifest)
)

for (GSE224463_contrast_index in
     seq_len(nrow(GSE224463_contrast_manifest))) {

  GSE224463_one_contrast <-
    GSE224463_contrast_manifest[GSE224463_contrast_index, ]

  GSE224463_one_series <- GSE224463_one_contrast$series[[1L]]
  GSE224463_one_target_time <-
    GSE224463_one_contrast$target_time_min[[1L]]
  GSE224463_one_analysis_id <-
    GSE224463_one_contrast$analysis_id[[1L]]

  GSE224463_one_manifest <- GSE224463_WT_manifest %>%
    dplyr::filter(
      .data$series == GSE224463_one_series,
      .data$time_min %in% c(0L, GSE224463_one_target_time)
    ) %>%
    dplyr::arrange(.data$author_replicate, .data$time_min)

  GSE224463_one_replicate_count <- GSE224463_one_manifest %>%
    dplyr::count(.data$time_min, name = "N_biological_replicates")

  if (nrow(GSE224463_one_manifest) != 6L ||
      nrow(GSE224463_one_replicate_count) != 2L ||
      any(GSE224463_one_replicate_count$N_biological_replicates != 3L)) {
    stop(
      "The contrast ",
      GSE224463_one_analysis_id,
      " does not contain exactly three biological replicates at 0 min and ",
      "the target time point.",
      call. = FALSE
    )
  }

  GSE224463_one_count_matrix <- GSE224463_raw_count_matrix[
    ,
    GSE224463_one_manifest$file_name,
    drop = FALSE
  ]

  GSE224463_one_nonzero_gene <- rowSums(
    GSE224463_one_count_matrix
  ) > 0L

  if (sum(GSE224463_one_nonzero_gene) < 1000L) {
    stop(
      "Too few non-zero genes are available for DESeq2 contrast ",
      GSE224463_one_analysis_id,
      ".",
      call. = FALSE
    )
  }

  GSE224463_one_count_matrix <- GSE224463_one_count_matrix[
    GSE224463_one_nonzero_gene,
    ,
    drop = FALSE
  ]

  GSE224463_one_coldata <- GSE224463_one_manifest %>%
    dplyr::mutate(
      author_replicate = factor(.data$author_replicate),
      DESeq2_timepoint = factor(
        dplyr::if_else(
          .data$time_min == 0L,
          "Baseline_0min",
          "Target_timepoint"
        ),
        levels = c("Baseline_0min", "Target_timepoint")
      )
    ) %>%
    dplyr::select(
      .data$file_name,
      .data$author_replicate,
      .data$DESeq2_timepoint
    ) %>%
    as.data.frame()

  rownames(GSE224463_one_coldata) <-
    GSE224463_one_coldata$file_name
  GSE224463_one_coldata$file_name <- NULL

  if (!identical(
    colnames(GSE224463_one_count_matrix),
    rownames(GSE224463_one_coldata)
  )) {
    stop(
      "Count-matrix columns and DESeq2 colData rows are not aligned for ",
      GSE224463_one_analysis_id,
      ".",
      call. = FALSE
    )
  }

  GSE224463_one_dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = GSE224463_one_count_matrix,
    colData = GSE224463_one_coldata,
    design = ~ author_replicate + DESeq2_timepoint
  )

  GSE224463_one_dds <- DESeq2::DESeq(
    GSE224463_one_dds,
    quiet = TRUE
  )

  GSE224463_one_DESeq2_result <- DESeq2::results(
    GSE224463_one_dds,
    contrast = c(
      "DESeq2_timepoint",
      "Target_timepoint",
      "Baseline_0min"
    ),
    alpha = 0.05,
    independentFiltering = TRUE,
    cooksCutoff = TRUE
  )

  GSE224463_one_normalised_counts <- DESeq2::counts(
    GSE224463_one_dds,
    normalized = TRUE
  )

  GSE224463_one_reference_columns <- rownames(
    GSE224463_one_coldata
  )[
    GSE224463_one_coldata$DESeq2_timepoint == "Baseline_0min"
  ]

  GSE224463_one_target_columns <- rownames(
    GSE224463_one_coldata
  )[
    GSE224463_one_coldata$DESeq2_timepoint == "Target_timepoint"
  ]

  GSE224463_one_DE_all <- tibble::as_tibble(
    as.data.frame(GSE224463_one_DESeq2_result),
    rownames = "Gene_id"
  ) %>%
    dplyr::transmute(
      Gene_id = tolower(as.character(.data$Gene_id)),
      base_mean = as.numeric(.data$baseMean),
      log2FC_target_time_minus_0min = as.numeric(.data$log2FoldChange),
      lfc_standard_error = as.numeric(.data$lfcSE),
      wald_statistic = as.numeric(.data$stat),
      p_value = as.numeric(.data$pvalue),
      padj = as.numeric(.data$padj),
      reference_0min_mean_normalised_counts = rowMeans(
        GSE224463_one_normalised_counts[
          ,
          GSE224463_one_reference_columns,
          drop = FALSE
        ]
      ),
      target_time_mean_normalised_counts = rowMeans(
        GSE224463_one_normalised_counts[
          ,
          GSE224463_one_target_columns,
          drop = FALSE
        ]
      ),
      analysis_id = GSE224463_one_analysis_id,
      series = GSE224463_one_series,
      target_time_min = GSE224463_one_target_time,
      comparison_label = GSE224463_one_contrast$comparison_label[[1L]],
      biological_context = GSE224463_one_contrast$biological_context[[1L]],
      contrast_QC_status = GSE224463_one_contrast$contrast_QC_status[[1L]]
    )

  GSE224463_one_DE_thermodynamic <- GSE224463_one_DE_all %>%
    dplyr::left_join(
      GSE224463_thermodynamic_annotation_audit %>%
        dplyr::select(
          .data$Gene_id,
          .data$transcript_length_nt,
          .data$Normalized_MFE,
          .data$included_in_DE_descriptor_analysis,
          .data$included_in_full_transcriptome_weighted_analysis
        ),
      by = "Gene_id"
    )

  GSE224463_one_significant_DE_thermodynamic <-
    GSE224463_one_DE_thermodynamic %>%
    dplyr::filter(
      !is.na(.data$padj),
      .data$padj < 0.05,
      is.finite(.data$log2FC_target_time_minus_0min),
      .data$log2FC_target_time_minus_0min != 0,
      .data$included_in_DE_descriptor_analysis %in% TRUE
    ) %>%
    dplyr::mutate(
      FC_sign = dplyr::case_when(
        .data$log2FC_target_time_minus_0min < 0 ~ "log2_FC < 0",
        .data$log2FC_target_time_minus_0min > 0 ~ "log2_FC > 0"
      ),
      FC_sign = factor(
        .data$FC_sign,
        levels = c("log2_FC < 0", "log2_FC > 0")
      )
    ) %>%
    dplyr::filter(!is.na(.data$FC_sign))

  GSE224463_one_DE_summary <- GSE224463_one_DE_all %>%
    dplyr::summarise(
      analysis_id = GSE224463_one_analysis_id,
      series = GSE224463_one_series,
      target_time_min = GSE224463_one_target_time,
      comparison_label = GSE224463_one_contrast$comparison_label[[1L]],
      biological_context = GSE224463_one_contrast$biological_context[[1L]],
      contrast_QC_status = GSE224463_one_contrast$contrast_QC_status[[1L]],
      N_genes_tested = dplyr::n(),
      N_significant_padj_lt_0_05 = sum(
        !is.na(.data$padj) & .data$padj < 0.05
      ),
      N_target_higher = sum(
        !is.na(.data$padj) &
          .data$padj < 0.05 &
          .data$log2FC_target_time_minus_0min > 0,
        na.rm = TRUE
      ),
      N_target_lower = sum(
        !is.na(.data$padj) &
          .data$padj < 0.05 &
          .data$log2FC_target_time_minus_0min < 0,
        na.rm = TRUE
      ),
      N_significant_with_finite_MFE = nrow(
        GSE224463_one_significant_DE_thermodynamic
      )
    )

  GSE224463_one_normalised_count_table <- tibble::as_tibble(
    as.data.frame(GSE224463_one_normalised_counts),
    rownames = "Gene_id"
  ) %>%
    dplyr::mutate(
      Gene_id = tolower(as.character(.data$Gene_id))
    )

  readr::write_tsv(
    GSE224463_one_DE_all,
    file.path(
      GSE224463_DE_table_dir,
      paste0(GSE224463_one_analysis_id, "_DESeq2_all_gene_results.tsv")
    )
  )

  readr::write_tsv(
    GSE224463_one_DE_thermodynamic,
    file.path(
      GSE224463_DE_table_dir,
      paste0(GSE224463_one_analysis_id, "_DESeq2_MFE_results.tsv")
    )
  )

  readr::write_tsv(
    GSE224463_one_significant_DE_thermodynamic,
    file.path(
      GSE224463_DE_table_dir,
      paste0(
        GSE224463_one_analysis_id,
        "_significant_DE_MFE_transcripts.tsv"
      )
    )
  )

  readr::write_tsv(
    GSE224463_one_normalised_count_table,
    file.path(
      GSE224463_DE_table_dir,
      paste0(
        GSE224463_one_analysis_id,
        "_DESeq2_normalised_counts.tsv.gz"
      )
    )
  )

  GSE224463_DE_all_list[[GSE224463_contrast_index]] <-
    GSE224463_one_DE_all
  GSE224463_DE_thermodynamic_list[[GSE224463_contrast_index]] <-
    GSE224463_one_DE_thermodynamic
  GSE224463_significant_DE_thermodynamic_list[[GSE224463_contrast_index]] <-
    GSE224463_one_significant_DE_thermodynamic
  GSE224463_DE_summary_list[[GSE224463_contrast_index]] <-
    GSE224463_one_DE_summary
}

GSE224463_DE_all <- dplyr::bind_rows(GSE224463_DE_all_list)
GSE224463_DE_thermodynamic <- dplyr::bind_rows(
  GSE224463_DE_thermodynamic_list
)
GSE224463_significant_DE_thermodynamic <- dplyr::bind_rows(
  GSE224463_significant_DE_thermodynamic_list
)
GSE224463_DE_summary <- dplyr::bind_rows(GSE224463_DE_summary_list) %>%
  dplyr::arrange(
    factor(.data$analysis_id, levels = GSE224463_contrast_manifest$analysis_id)
  )

print(GSE224463_DE_summary)

readr::write_tsv(
  GSE224463_DE_summary,
  file.path(
    GSE224463_table_dir,
    "GSE224463_WT_stationary_phase_time_course_DESeq2_summary.tsv"
  )
)

readr::write_tsv(
  GSE224463_DE_all,
  file.path(
    GSE224463_table_dir,
    "GSE224463_WT_stationary_phase_time_course_DESeq2_all_results.tsv.gz"
  )
)

readr::write_tsv(
  GSE224463_DE_thermodynamic,
  file.path(
    GSE224463_table_dir,
    "GSE224463_WT_stationary_phase_time_course_DESeq2_MFE_results.tsv.gz"
  )
)

readr::write_tsv(
  GSE224463_significant_DE_thermodynamic,
  file.path(
    GSE224463_table_dir,
    "GSE224463_WT_stationary_phase_time_course_significant_DE_MFE_transcripts.tsv.gz"
  )
)

GSE224463_distribution_statistics_list <- vector(
  "list",
  nrow(GSE224463_contrast_manifest)
)

for (GSE224463_contrast_index in
     seq_len(nrow(GSE224463_contrast_manifest))) {

  GSE224463_one_contrast <-
    GSE224463_contrast_manifest[GSE224463_contrast_index, ]

  GSE224463_one_analysis_id <-
    GSE224463_one_contrast$analysis_id[[1L]]

  GSE224463_one_sig <- GSE224463_significant_DE_thermodynamic %>%
    dplyr::filter(.data$analysis_id == GSE224463_one_analysis_id)

  GSE224463_one_MFE_lower <- GSE224463_one_sig$Normalized_MFE[
    GSE224463_one_sig$FC_sign == "log2_FC < 0"
  ]

  GSE224463_one_MFE_higher <- GSE224463_one_sig$Normalized_MFE[
    GSE224463_one_sig$FC_sign == "log2_FC > 0"
  ]

  GSE224463_one_MFE_Wilcoxon_p <- if (
    length(GSE224463_one_MFE_lower) > 0L &&
      length(GSE224463_one_MFE_higher) > 0L
  ) {
    tryCatch(
      stats::wilcox.test(
        GSE224463_one_MFE_higher,
        GSE224463_one_MFE_lower,
        alternative = "two.sided",
        exact = FALSE
      )$p.value,
      error = function(error_message) NA_real_
    )
  } else {
    NA_real_
  }

  GSE224463_distribution_statistics_list[[GSE224463_contrast_index]] <-
    tibble::tibble(
      analysis_id = GSE224463_one_analysis_id,
      series = GSE224463_one_contrast$series[[1L]],
      target_time_min = GSE224463_one_contrast$target_time_min[[1L]],
      comparison_label = GSE224463_one_contrast$comparison_label[[1L]],
      biological_context = GSE224463_one_contrast$biological_context[[1L]],
      contrast_QC_status = GSE224463_one_contrast$contrast_QC_status[[1L]],
      descriptor = "Normalized_MFE",
      N_target_lower = length(GSE224463_one_MFE_lower),
      N_target_higher = length(GSE224463_one_MFE_higher),
      mean_target_lower = if (
        length(GSE224463_one_MFE_lower) > 0L
      ) mean(GSE224463_one_MFE_lower) else NA_real_,
      mean_target_higher = if (
        length(GSE224463_one_MFE_higher) > 0L
      ) mean(GSE224463_one_MFE_higher) else NA_real_,
      delta_target_higher_minus_lower = if (
        length(GSE224463_one_MFE_lower) > 0L &&
          length(GSE224463_one_MFE_higher) > 0L
      ) {
        mean(GSE224463_one_MFE_higher) -
          mean(GSE224463_one_MFE_lower)
      } else {
        NA_real_
      },
      Wilcoxon_p = GSE224463_one_MFE_Wilcoxon_p
    )
}

GSE224463_significant_DE_distribution_statistics <-
  dplyr::bind_rows(GSE224463_distribution_statistics_list) %>%
  dplyr::group_by(.data$descriptor) %>%
  dplyr::mutate(
    Wilcoxon_p_BH = stats::p.adjust(.data$Wilcoxon_p, method = "BH")
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(
    .data$descriptor,
    factor(
      .data$analysis_id,
      levels = GSE224463_contrast_manifest$analysis_id
    )
  )

print(GSE224463_significant_DE_distribution_statistics)

readr::write_tsv(
  GSE224463_significant_DE_distribution_statistics,
  file.path(
    GSE224463_table_dir,
    "GSE224463_WT_stationary_phase_time_course_significant_DE_MFE_directional_statistics.tsv"
  )
)

GSE224463_significant_DE_thermodynamic <-
  GSE224463_significant_DE_thermodynamic %>%
  dplyr::mutate(
    comparison_label = factor(
      .data$comparison_label,
      levels = GSE224463_contrast_manifest$comparison_label
    )
  )

if (nrow(GSE224463_significant_DE_thermodynamic) > 0L) {

  GSE224463_significant_MFE_y_limits <- range(
    GSE224463_significant_DE_thermodynamic$Normalized_MFE,
    finite = TRUE
  )

  GSE224463_significant_MFE_violin_plot <- ggplot(
    GSE224463_significant_DE_thermodynamic,
    aes(
      x = .data$FC_sign,
      y = .data$Normalized_MFE,
      fill = .data$FC_sign
    )
  ) +
    geom_violin(
      trim = FALSE,
      alpha = 0.75,
      colour = "black",
      linewidth = 0.2
    ) +
    geom_boxplot(
      width = 0.15,
      outlier.shape = NA,
      fill = "white",
      linewidth = 0.25
    ) +
    facet_wrap(
      ~ comparison_label,
      ncol = 5,
      drop = FALSE
    ) +
    coord_cartesian(ylim = GSE224463_significant_MFE_y_limits) +
    scale_fill_manual(
      values = c(
        "log2_FC < 0" = "#8cc5e3",
        "log2_FC > 0" = "#2066a8"
      )
    ) +
    scale_x_discrete(
      labels = c(
        "log2_FC < 0" = "Lower",
        "log2_FC > 0" = "Higher"
      )
    ) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.text.x = element_text(
        angle = 30,
        hjust = 1,
        colour = "black",
        size = 8
      ),
      axis.text.y = element_text(colour = "black", size = 8),
      strip.text = element_text(size = 7)
    ) +
    labs(
      title = "Significant-DE transcript ΔGfold distributions",
      subtitle = "Target time point versus matched 0-min exponential reference",
      y = "Normalized ΔGfold (ΔGfold / nt)"
    )

  print(GSE224463_significant_MFE_violin_plot)

  ggsave(
    filename = file.path(
      GSE224463_plot_dir,
      "GSE224463_WT_stationary_phase_time_course_significant_DE_MFE_violin_plots.tiff"
    ),
    plot = GSE224463_significant_MFE_violin_plot,
    width = 15,
    height = 7,
    units = "in",
    dpi = 600,
    compression = "lzw"
  )

} else {
  warning(
    "No FDR-significant DE transcripts with finite ΔGfold values were ",
    "retained across the ten contrasts; significant-DE distribution plots were ",
    "not written.",
    call. = FALSE
  )
}

GSE224463_significant_DE_abundance_fraction_list <- vector(
  "list",
  nrow(GSE224463_contrast_manifest)
)

for (GSE224463_contrast_index in
     seq_len(nrow(GSE224463_contrast_manifest))) {

  GSE224463_one_contrast <-
    GSE224463_contrast_manifest[GSE224463_contrast_index, ]

  GSE224463_one_analysis_id <-
    GSE224463_one_contrast$analysis_id[[1L]]

  GSE224463_one_sig <- GSE224463_significant_DE_thermodynamic %>%
    dplyr::filter(.data$analysis_id == GSE224463_one_analysis_id)

  if (nrow(GSE224463_one_sig) >= 4L &&
      sum(
        GSE224463_one_sig$reference_0min_mean_normalised_counts,
        na.rm = TRUE
      ) > 0 &&
      sum(
        GSE224463_one_sig$target_time_mean_normalised_counts,
        na.rm = TRUE
      ) > 0) {

    GSE224463_one_MFE_quartiles <- stats::quantile(
      GSE224463_one_sig$Normalized_MFE,
      probs = c(0.25, 0.75),
      na.rm = TRUE,
      names = FALSE
    )

    GSE224463_one_sig <- GSE224463_one_sig %>%
      dplyr::mutate(
        MFE_group = dplyr::case_when(
          .data$Normalized_MFE <= GSE224463_one_MFE_quartiles[1L] ~
            "Low_MFE",
          .data$Normalized_MFE >= GSE224463_one_MFE_quartiles[2L] ~
            "High_MFE",
          TRUE ~ "Intermediate_MFE"
        )
      )

    GSE224463_significant_DE_abundance_fraction_list[[GSE224463_contrast_index]] <-
      GSE224463_one_sig %>%
      dplyr::summarise(
        analysis_id = GSE224463_one_analysis_id,
        series = GSE224463_one_contrast$series[[1L]],
        target_time_min = GSE224463_one_contrast$target_time_min[[1L]],
        comparison_label = GSE224463_one_contrast$comparison_label[[1L]],
        contrast_QC_status = GSE224463_one_contrast$contrast_QC_status[[1L]],
        Control = "0 min exponential phase",
        Treated = paste0(
          GSE224463_one_contrast$target_time_min[[1L]],
          " min"
        ),
        N_significant_transcripts = dplyr::n(),
        N_low_MFE = sum(.data$MFE_group == "Low_MFE"),
        N_high_MFE = sum(.data$MFE_group == "High_MFE"),
        low_MFE_quartile_cutoff = GSE224463_one_MFE_quartiles[1L],
        high_MFE_quartile_cutoff = GSE224463_one_MFE_quartiles[2L],
        Control_low_fraction = sum(
          .data$reference_0min_mean_normalised_counts[
            .data$MFE_group == "Low_MFE"
          ]
        ) / sum(.data$reference_0min_mean_normalised_counts),
        Treated_low_fraction = sum(
          .data$target_time_mean_normalised_counts[
            .data$MFE_group == "Low_MFE"
          ]
        ) / sum(.data$target_time_mean_normalised_counts),
        Control_high_fraction = sum(
          .data$reference_0min_mean_normalised_counts[
            .data$MFE_group == "High_MFE"
          ]
        ) / sum(.data$reference_0min_mean_normalised_counts),
        Treated_high_fraction = sum(
          .data$target_time_mean_normalised_counts[
            .data$MFE_group == "High_MFE"
          ]
        ) / sum(.data$target_time_mean_normalised_counts),
        abundance_fraction_eligible = TRUE
      )

  } else {

    GSE224463_significant_DE_abundance_fraction_list[[GSE224463_contrast_index]] <-
      tibble::tibble(
      analysis_id = GSE224463_one_analysis_id,
      series = GSE224463_one_contrast$series[[1L]],
      target_time_min = GSE224463_one_contrast$target_time_min[[1L]],
      comparison_label = GSE224463_one_contrast$comparison_label[[1L]],
      contrast_QC_status = GSE224463_one_contrast$contrast_QC_status[[1L]],
      Control = "0 min exponential phase",
      Treated = paste0(GSE224463_one_contrast$target_time_min[[1L]], " min"),
      N_significant_transcripts = nrow(GSE224463_one_sig),
      N_low_MFE = NA_integer_,
      N_high_MFE = NA_integer_,
      low_MFE_quartile_cutoff = NA_real_,
      high_MFE_quartile_cutoff = NA_real_,
      Control_low_fraction = NA_real_,
      Treated_low_fraction = NA_real_,
      Control_high_fraction = NA_real_,
      Treated_high_fraction = NA_real_,
      abundance_fraction_eligible = FALSE
    )
  }
}

GSE224463_significant_DE_abundance_fraction <- dplyr::bind_rows(
  GSE224463_significant_DE_abundance_fraction_list
) %>%
  dplyr::mutate(
    comparison_label = factor(
      .data$comparison_label,
      levels = GSE224463_contrast_manifest$comparison_label
    )
  ) %>%
  dplyr::arrange(.data$comparison_label)

print(GSE224463_significant_DE_abundance_fraction)

readr::write_tsv(
  GSE224463_significant_DE_abundance_fraction,
  file.path(
    GSE224463_table_dir,
    "GSE224463_WT_stationary_phase_time_course_significant_DE_MFE_abundance_fractions.tsv"
  )
)

GSE224463_high_MFE_fraction_plot_data <-
  GSE224463_significant_DE_abundance_fraction %>%
  dplyr::filter(.data$abundance_fraction_eligible) %>%
  dplyr::select(
    .data$comparison_label,
    .data$Control,
    .data$Treated,
    .data$Control_high_fraction,
    .data$Treated_high_fraction
  ) %>%
  tidyr::pivot_longer(
    cols = c(
      .data$Control_high_fraction,
      .data$Treated_high_fraction
    ),
    names_to = "Condition_type",
    values_to = "Fraction"
  ) %>%
  dplyr::mutate(
    Sample_group = dplyr::if_else(
      .data$Condition_type == "Control_high_fraction",
      "Control",
      "Target"
    ),
    Sample_group = factor(.data$Sample_group, levels = c("Control", "Target"))
  )

GSE224463_low_MFE_fraction_plot_data <-
  GSE224463_significant_DE_abundance_fraction %>%
  dplyr::filter(.data$abundance_fraction_eligible) %>%
  dplyr::select(
    .data$comparison_label,
    .data$Control,
    .data$Treated,
    .data$Control_low_fraction,
    .data$Treated_low_fraction
  ) %>%
  tidyr::pivot_longer(
    cols = c(
      .data$Control_low_fraction,
      .data$Treated_low_fraction
    ),
    names_to = "Condition_type",
    values_to = "Fraction"
  ) %>%
  dplyr::mutate(
    Sample_group = dplyr::if_else(
      .data$Condition_type == "Control_low_fraction",
      "Control",
      "Target"
    ),
    Sample_group = factor(.data$Sample_group, levels = c("Control", "Target"))
  )

if (nrow(GSE224463_high_MFE_fraction_plot_data) > 0L &&
    nrow(GSE224463_low_MFE_fraction_plot_data) > 0L) {

  GSE224463_high_MFE_fraction_plot <- ggplot(
    GSE224463_high_MFE_fraction_plot_data,
    aes(
      x = .data$comparison_label,
      y = .data$Fraction,
      fill = .data$Sample_group
    )
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
      size = 2.7
    ) +
    coord_cartesian(
      ylim = c(
        0,
        max(GSE224463_high_MFE_fraction_plot_data$Fraction, na.rm = TRUE) *
          1.18
      )
    ) +
    scale_fill_manual(
      values = c("Control" = "#8cc5e3", "Target" = "#2066a8"),
      labels = c("Control" = "0 min exponential", "Target" = "Target time")
    ) +
    theme_classic() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        colour = "black",
        size = 8
      ),
      legend.position = "top"
    ) +
    labs(
      title = "High-ΔGfold significant-DE transcript abundance fraction",
      y = "Fraction of significant-DE transcript abundance",
      fill = ""
    )

  GSE224463_low_MFE_fraction_plot <- ggplot(
    GSE224463_low_MFE_fraction_plot_data,
    aes(
      x = .data$comparison_label,
      y = .data$Fraction,
      fill = .data$Sample_group
    )
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
      size = 2.7
    ) +
    coord_cartesian(
      ylim = c(
        0,
        max(GSE224463_low_MFE_fraction_plot_data$Fraction, na.rm = TRUE) *
          1.18
      )
    ) +
    scale_fill_manual(
      values = c("Control" = "#8cc5e3", "Target" = "#2066a8"),
      labels = c("Control" = "0 min exponential", "Target" = "Target time")
    ) +
    theme_classic() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        colour = "black",
        size = 8
      ),
      legend.position = "top"
    ) +
    labs(
      title = "Low-ΔGfold significant-DE transcript abundance fraction",
      y = "Fraction of significant-DE transcript abundance",
      fill = ""
    )

  print(GSE224463_high_MFE_fraction_plot)
  print(GSE224463_low_MFE_fraction_plot)

  ggsave(
    filename = file.path(
      GSE224463_plot_dir,
      "GSE224463_WT_stationary_phase_time_course_significant_DE_high_MFE_abundance_fraction.tiff"
    ),
    plot = GSE224463_high_MFE_fraction_plot,
    width = 14,
    height = 6,
    units = "in",
    dpi = 600,
    compression = "lzw"
  )

  ggsave(
    filename = file.path(
      GSE224463_plot_dir,
      "GSE224463_WT_stationary_phase_time_course_significant_DE_low_MFE_abundance_fraction.tiff"
    ),
    plot = GSE224463_low_MFE_fraction_plot,
    width = 14,
    height = 6,
    units = "in",
    dpi = 600,
    compression = "lzw"
  )

} else {
  warning(
    "No contrast retained enough significant thermodynamic transcripts to ",
    "calculate high-/low-ΔGfold abundance fractions.",
    call. = FALSE
  )
}

GSE224463_counts_for_weighted_analysis <- GSE224463_raw_count_matrix[
  GSE224463_weighted_annotation$Gene_id,
  GSE224463_WT_manifest$file_name,
  drop = FALSE
]

if (!identical(
  rownames(GSE224463_counts_for_weighted_analysis),
  GSE224463_weighted_annotation$Gene_id
)) {
  stop(
    "The full-transcriptome weighted count matrix is not aligned to the ",
    "validated annotation table.",
    call. = FALSE
  )
}

GSE224463_length_adjusted_counts <- sweep(
  GSE224463_counts_for_weighted_analysis,
  MARGIN = 1L,
  STATS = GSE224463_weighted_annotation$transcript_length_nt,
  FUN = "/"
)

GSE224463_length_adjusted_library_totals <- colSums(
  GSE224463_length_adjusted_counts
)

if (any(!is.finite(GSE224463_length_adjusted_library_totals)) ||
    any(GSE224463_length_adjusted_library_totals <= 0)) {
  stop(
    "At least one sample has an invalid length-adjusted library total.",
    call. = FALSE
  )
}

GSE224463_relative_RNA_abundance <- sweep(
  GSE224463_length_adjusted_counts,
  MARGIN = 2L,
  STATS = GSE224463_length_adjusted_library_totals,
  FUN = "/"
)

if (any(!is.finite(GSE224463_relative_RNA_abundance)) ||
    any(abs(colSums(GSE224463_relative_RNA_abundance) - 1) > 1e-10)) {
  stop(
    "The full-transcriptome relative-RNA-abundance matrix is invalid.",
    call. = FALSE
  )
}

GSE224463_weighted_sample_metadata <- GSE224463_WT_manifest %>%
  dplyr::mutate(
    total_raw_counts_all_genes = as.numeric(
      colSums(GSE224463_raw_count_matrix)[.data$file_name]
    ),
    total_raw_counts_full_valid_transcriptome = as.numeric(
      colSums(GSE224463_counts_for_weighted_analysis)[.data$file_name]
    ),
    length_adjusted_library_total = as.numeric(
      GSE224463_length_adjusted_library_totals[.data$file_name]
    )
  ) %>%
  dplyr::arrange(
    factor(.data$series, levels = c("short", "long")),
    .data$author_replicate,
    .data$sample_order_within_series
  )

GSE224463_abundance_long <- tibble::as_tibble(
  as.data.frame(GSE224463_relative_RNA_abundance),
  rownames = "Gene_id"
) %>%
  tidyr::pivot_longer(
    cols = -Gene_id,
    names_to = "file_name",
    values_to = "relative_RNA_molecule_abundance"
  ) %>%
  dplyr::left_join(
    GSE224463_weighted_annotation %>%
      dplyr::select(
        .data$Gene_id,
        .data$transcript_length_nt,
        .data$Normalized_MFE
      ),
    by = "Gene_id"
  ) %>%
  dplyr::left_join(
    GSE224463_weighted_sample_metadata,
    by = "file_name"
  )

GSE224463_weighted_sample_gene_audit <- GSE224463_abundance_long %>%
  dplyr::group_by(
    .data$file_name,
    .data$series,
    .data$author_replicate,
    .data$time_min
  ) %>%
  dplyr::summarise(
    N_full_transcriptome_genes = dplyr::n(),
    relative_abundance_total = sum(.data$relative_RNA_molecule_abundance),
    .groups = "drop"
  )

if (any(
  GSE224463_weighted_sample_gene_audit$N_full_transcriptome_genes !=
    nrow(GSE224463_weighted_annotation)
) ||
    any(
      abs(
        GSE224463_weighted_sample_gene_audit$relative_abundance_total - 1
      ) > 1e-10
    )) {
  stop(
    "The full-transcriptome weighted-abundance audit failed.",
    call. = FALSE
  )
}

print(GSE224463_weighted_sample_gene_audit)

readr::write_tsv(
  GSE224463_weighted_annotation,
  file.path(
    GSE224463_table_dir,
    "GSE224463_full_transcriptome_weighted_analysis_annotation.tsv"
  )
)

readr::write_tsv(
  GSE224463_weighted_sample_metadata,
  file.path(
    GSE224463_table_dir,
    "GSE224463_full_transcriptome_weighted_analysis_sample_metadata.tsv"
  )
)

readr::write_tsv(
  GSE224463_weighted_sample_gene_audit,
  file.path(
    GSE224463_table_dir,
    "GSE224463_full_transcriptome_weighted_analysis_sample_gene_audit.tsv"
  )
)

readr::write_tsv(
  GSE224463_abundance_long %>%
    dplyr::select(
      .data$Gene_id,
      .data$transcript_length_nt,
      .data$Normalized_MFE,
      .data$file_name,
      .data$series,
      .data$author_replicate,
      .data$time_min,
      .data$GSM,
      .data$condition_label,
      .data$relative_RNA_molecule_abundance
    ),
  file.path(
    GSE224463_table_dir,
    "GSE224463_full_transcriptome_gene_length_adjusted_relative_RNA_abundance.tsv.gz"
  )
)

GSE224463_weighted_sample_indices <- GSE224463_abundance_long %>%
  dplyr::group_by(
    .data$file_name,
    .data$series,
    .data$author_replicate,
    .data$time_min,
    .data$GSM,
    .data$condition_label,
    .data$library_QC_flag
  ) %>%
  dplyr::summarise(
    N_full_transcriptome_genes = dplyr::n(),
    abundance_weighted_Normalized_MFE = sum(
      .data$relative_RNA_molecule_abundance * .data$Normalized_MFE
    ),
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    factor(.data$series, levels = c("short", "long")),
    .data$author_replicate,
    .data$time_min
  )

if (nrow(GSE224463_weighted_sample_indices) != 36L) {
  stop(
    "The full-transcriptome sample-index table must contain one row per ",
    "retained WT sample.",
    call. = FALSE
  )
}

GSE224463_weighted_timecourse_summary <- GSE224463_weighted_sample_indices %>%
  dplyr::group_by(.data$series, .data$time_min) %>%
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
  ) %>%
  dplyr::arrange(
    factor(.data$series, levels = c("short", "long")),
    .data$time_min
  )

print(GSE224463_weighted_timecourse_summary)

readr::write_tsv(
  GSE224463_weighted_sample_indices,
  file.path(
    GSE224463_table_dir,
    "GSE224463_full_transcriptome_abundance_weighted_MFE_sample_indices.tsv"
  )
)

readr::write_tsv(
  GSE224463_weighted_timecourse_summary,
  file.path(
    GSE224463_table_dir,
    "GSE224463_full_transcriptome_abundance_weighted_MFE_time_course_summary.tsv"
  )
)

GSE224463_weighted_MFE_timecourse_plot <- ggplot(
  GSE224463_weighted_timecourse_summary,
  aes(
    x = .data$time_min,
    y = .data$mean_abundance_weighted_Normalized_MFE,
    group = .data$series,
    colour = .data$series
  )
) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_errorbar(
    aes(
      ymin =
        .data$mean_abundance_weighted_Normalized_MFE -
        .data$sem_abundance_weighted_Normalized_MFE,
      ymax =
        .data$mean_abundance_weighted_Normalized_MFE +
        .data$sem_abundance_weighted_Normalized_MFE
    ),
    width = 7,
    linewidth = 0.4
  ) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 2.4) +
  facet_wrap(~ series, nrow = 1, scales = "free_x") +
  scale_colour_manual(values = c("short" = "#8cc5e3", "long" = "#2066a8")) +
  theme_classic() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 10)
  ) +
  labs(
    title = "Full-transcriptome abundance-weighted RNA folding stability",
    subtitle = "Mean ± SEM across author biological replicates",
    x = "Time after exponential-phase baseline (min)",
    y = "Abundance-weighted Normalized ΔGfold"
  )

print(GSE224463_weighted_MFE_timecourse_plot)

ggsave(
  filename = file.path(
    GSE224463_plot_dir,
    "GSE224463_WT_stationary_phase_time_course_full_transcriptome_weighted_MFE_time_course.tiff"
  ),
  plot = GSE224463_weighted_MFE_timecourse_plot,
  width = 11,
  height = 5,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

GSE224463_weighted_condition_profiles <- GSE224463_abundance_long %>%
  dplyr::group_by(
    .data$series,
    .data$time_min,
    .data$Gene_id,
    .data$transcript_length_nt,
    .data$Normalized_MFE
  ) %>%
  dplyr::summarise(
    mean_relative_RNA_molecule_abundance = mean(
      .data$relative_RNA_molecule_abundance
    ),
    .groups = "drop"
  )

GSE224463_weighted_metric_statistics_list <- vector(
  "list",
  nrow(GSE224463_contrast_manifest)
)

GSE224463_weighted_gene_contribution_list <- vector(
  "list",
  nrow(GSE224463_contrast_manifest)
)

for (GSE224463_contrast_index in
     seq_len(nrow(GSE224463_contrast_manifest))) {

  GSE224463_one_contrast <-
    GSE224463_contrast_manifest[GSE224463_contrast_index, ]

  GSE224463_one_series <- GSE224463_one_contrast$series[[1L]]
  GSE224463_one_target_time <-
    GSE224463_one_contrast$target_time_min[[1L]]
  GSE224463_one_analysis_id <-
    GSE224463_one_contrast$analysis_id[[1L]]

  GSE224463_one_reference_profile <-
    GSE224463_weighted_condition_profiles %>%
    dplyr::filter(
      .data$series == GSE224463_one_series,
      .data$time_min == 0L
    ) %>%
    dplyr::transmute(
      Gene_id = .data$Gene_id,
      transcript_length_nt = .data$transcript_length_nt,
      Normalized_MFE = .data$Normalized_MFE,
      reference_0min_mean_relative_RNA_abundance =
        .data$mean_relative_RNA_molecule_abundance
    )

  GSE224463_one_target_profile <-
    GSE224463_weighted_condition_profiles %>%
    dplyr::filter(
      .data$series == GSE224463_one_series,
      .data$time_min == GSE224463_one_target_time
    ) %>%
    dplyr::transmute(
      Gene_id = .data$Gene_id,
      transcript_length_nt = .data$transcript_length_nt,
      Normalized_MFE = .data$Normalized_MFE,
      target_time_mean_relative_RNA_abundance =
        .data$mean_relative_RNA_molecule_abundance
    )

  GSE224463_one_weighted_gene_data <-
    GSE224463_one_reference_profile %>%
    dplyr::inner_join(
      GSE224463_one_target_profile,
      by = c(
        "Gene_id",
        "transcript_length_nt",
        "Normalized_MFE"
      )
    ) %>%
    dplyr::mutate(
      contrast_weight_target_time_minus_0min =
        .data$target_time_mean_relative_RNA_abundance -
        .data$reference_0min_mean_relative_RNA_abundance,
      abundance_weighted_MFE_contribution =
        .data$contrast_weight_target_time_minus_0min *
        .data$Normalized_MFE
    )

  if (nrow(GSE224463_one_weighted_gene_data) !=
      nrow(GSE224463_weighted_annotation)) {
    stop(
      "The full-transcriptome profiles could not be aligned for ",
      GSE224463_one_analysis_id,
      ".",
      call. = FALSE
    )
  }

  GSE224463_one_reference_sample_indices <-
    GSE224463_weighted_sample_indices %>%
    dplyr::filter(
      .data$series == GSE224463_one_series,
      .data$time_min == 0L
    ) %>%
    dplyr::transmute(
      author_replicate = .data$author_replicate,
      reference_abundance_weighted_Normalized_MFE =
        .data$abundance_weighted_Normalized_MFE
    )

  GSE224463_one_target_sample_indices <-
    GSE224463_weighted_sample_indices %>%
    dplyr::filter(
      .data$series == GSE224463_one_series,
      .data$time_min == GSE224463_one_target_time
    ) %>%
    dplyr::transmute(
      author_replicate = .data$author_replicate,
      target_abundance_weighted_Normalized_MFE =
        .data$abundance_weighted_Normalized_MFE
    )

  GSE224463_one_paired_sample_indices <-
    GSE224463_one_reference_sample_indices %>%
    dplyr::inner_join(
      GSE224463_one_target_sample_indices,
      by = "author_replicate"
    ) %>%
    dplyr::arrange(.data$author_replicate)

  if (nrow(GSE224463_one_paired_sample_indices) != 3L) {
    stop(
      "The full-transcriptome paired replicate table does not contain three ",
      "matched samples for ",
      GSE224463_one_analysis_id,
      ".",
      call. = FALSE
    )
  }

  GSE224463_one_MFE_paired_t_p <- tryCatch(
    stats::t.test(
      GSE224463_one_paired_sample_indices$target_abundance_weighted_Normalized_MFE,
      GSE224463_one_paired_sample_indices$reference_abundance_weighted_Normalized_MFE,
      paired = TRUE
    )$p.value,
    error = function(error_message) NA_real_
  )

  GSE224463_one_MFE_paired_Wilcoxon_p <- tryCatch(
    stats::wilcox.test(
      GSE224463_one_paired_sample_indices$target_abundance_weighted_Normalized_MFE,
      GSE224463_one_paired_sample_indices$reference_abundance_weighted_Normalized_MFE,
      paired = TRUE,
      exact = FALSE
    )$p.value,
    error = function(error_message) NA_real_
  )

  GSE224463_one_observed_MFE_delta <- sum(
    GSE224463_one_weighted_gene_data$abundance_weighted_MFE_contribution
  )

  set.seed(GSE224463_resampling_seed + GSE224463_contrast_index)

  GSE224463_one_MFE_bootstrap_delta <- vapply(
    seq_len(GSE224463_B_BOOT),
    function(iteration) {
      sum(
        sample(
          GSE224463_one_weighted_gene_data$abundance_weighted_MFE_contribution,
          size = nrow(GSE224463_one_weighted_gene_data),
          replace = TRUE
        )
      )
    },
    numeric(1)
  )

  GSE224463_one_MFE_permutation_delta <- vapply(
    seq_len(GSE224463_B_PERM),
    function(iteration) {
      sum(
        GSE224463_one_weighted_gene_data$contrast_weight_target_time_minus_0min *
          sample(
            GSE224463_one_weighted_gene_data$Normalized_MFE,
            size = nrow(GSE224463_one_weighted_gene_data),
            replace = FALSE
          )
      )
    },
    numeric(1)
  )

  GSE224463_weighted_metric_statistics_list[[GSE224463_contrast_index]] <-
    tibble::tibble(
    analysis_id = GSE224463_one_analysis_id,
    series = GSE224463_one_series,
    target_time_min = GSE224463_one_target_time,
    comparison_label = GSE224463_one_contrast$comparison_label[[1L]],
    biological_context = GSE224463_one_contrast$biological_context[[1L]],
    contrast_QC_status = GSE224463_one_contrast$contrast_QC_status[[1L]],
    descriptor = "Normalized_MFE",
    N_full_transcriptome_genes = nrow(GSE224463_one_weighted_gene_data),
    reference_0min_abundance_weighted_value = sum(
      GSE224463_one_weighted_gene_data$reference_0min_mean_relative_RNA_abundance *
        GSE224463_one_weighted_gene_data$Normalized_MFE
    ),
    target_time_abundance_weighted_value = sum(
      GSE224463_one_weighted_gene_data$target_time_mean_relative_RNA_abundance *
        GSE224463_one_weighted_gene_data$Normalized_MFE
    ),
    delta_target_time_minus_0min = GSE224463_one_observed_MFE_delta,
    gene_bootstrap_CI_low = unname(
      stats::quantile(GSE224463_one_MFE_bootstrap_delta, 0.025)
    ),
    gene_bootstrap_CI_high = unname(
      stats::quantile(GSE224463_one_MFE_bootstrap_delta, 0.975)
    ),
    descriptor_label_permutation_p = (
      sum(
        abs(GSE224463_one_MFE_permutation_delta) >=
          abs(GSE224463_one_observed_MFE_delta)
      ) + 1L
    ) / (GSE224463_B_PERM + 1L),
    paired_biological_replicate_t_p = GSE224463_one_MFE_paired_t_p,
    paired_biological_replicate_Wilcoxon_p =
      GSE224463_one_MFE_paired_Wilcoxon_p,
    direction = dplyr::case_when(
      GSE224463_one_observed_MFE_delta > 0 ~
        "Target RNA pool is less negative / less stably folded",
      GSE224463_one_observed_MFE_delta < 0 ~
        "Target RNA pool is more negative / more stably folded",
      TRUE ~ "No change"
    ),
    B_boot = GSE224463_B_BOOT,
    B_perm = GSE224463_B_PERM
  )

  GSE224463_weighted_gene_contribution_list[[GSE224463_contrast_index]] <-
    GSE224463_one_weighted_gene_data %>%
    dplyr::transmute(
      analysis_id = GSE224463_one_analysis_id,
      series = GSE224463_one_series,
      target_time_min = GSE224463_one_target_time,
      comparison_label = GSE224463_one_contrast$comparison_label[[1L]],
      Gene_id = .data$Gene_id,
      transcript_length_nt = .data$transcript_length_nt,
      Normalized_MFE = .data$Normalized_MFE,
      reference_0min_mean_relative_RNA_abundance =
        .data$reference_0min_mean_relative_RNA_abundance,
      target_time_mean_relative_RNA_abundance =
        .data$target_time_mean_relative_RNA_abundance,
      contrast_weight_target_time_minus_0min =
        .data$contrast_weight_target_time_minus_0min,
      abundance_weighted_MFE_contribution =
        .data$abundance_weighted_MFE_contribution
    )
}

GSE224463_full_transcriptome_weighted_delta_summary <- dplyr::bind_rows(
  GSE224463_weighted_metric_statistics_list
) %>%
  dplyr::group_by(.data$descriptor) %>%
  dplyr::mutate(
    descriptor_label_permutation_p_BH = stats::p.adjust(
      .data$descriptor_label_permutation_p,
      method = "BH"
    ),
    paired_biological_replicate_t_p_BH = stats::p.adjust(
      .data$paired_biological_replicate_t_p,
      method = "BH"
    ),
    paired_biological_replicate_Wilcoxon_p_BH = stats::p.adjust(
      .data$paired_biological_replicate_Wilcoxon_p,
      method = "BH"
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(
    .data$descriptor,
    factor(
      .data$analysis_id,
      levels = GSE224463_contrast_manifest$analysis_id
    )
  )

GSE224463_full_transcriptome_weighted_gene_contributions <-
  dplyr::bind_rows(GSE224463_weighted_gene_contribution_list)

print(GSE224463_full_transcriptome_weighted_delta_summary)

readr::write_tsv(
  GSE224463_full_transcriptome_weighted_delta_summary,
  file.path(
    GSE224463_table_dir,
    "GSE224463_full_transcriptome_abundance_weighted_MFE_delta_summary.tsv"
  )
)

readr::write_tsv(
  GSE224463_full_transcriptome_weighted_gene_contributions,
  file.path(
    GSE224463_table_dir,
    "GSE224463_full_transcriptome_abundance_weighted_MFE_gene_contributions.tsv.gz"
  )
)

GSE224463_required_GSEA_packages <- c(
  "clusterProfiler",
  "httr",
  "jsonlite",
  "stringr"
)

GSE224463_missing_GSEA_packages <- GSE224463_required_GSEA_packages[
  !vapply(
    GSE224463_required_GSEA_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(GSE224463_missing_GSEA_packages) > 0L) {
  stop(
    "Install the required GSEA package(s) before continuing: ",
    paste(GSE224463_missing_GSEA_packages, collapse = ", "),
    call. = FALSE
  )
}

GSE224463_gsea_dir <- file.path(
  GSE224463_output_dir,
  "GSEA"
)

GSE224463_gsea_rank_dir <- file.path(
  GSE224463_gsea_dir,
  "ranked_gene_lists"
)

GSE224463_gsea_result_dir <- file.path(
  GSE224463_gsea_dir,
  "results_by_contrast"
)

for (one_dir in c(
  GSE224463_gsea_dir,
  GSE224463_gsea_rank_dir,
  GSE224463_gsea_result_dir
)) {
  dir.create(one_dir, recursive = TRUE, showWarnings = FALSE)
}

GSE224463_refresh_GO_cache <- FALSE
GSE224463_GSEA_seed <- 20260808L

if (!file.exists(GSE224463_gene_info_file)) {
  stop(
    "The MG1655 gene-information file is missing:\n",
    GSE224463_gene_info_file,
    "\nUpdate GSE224463_gene_info_file before continuing.",
    call. = FALSE
  )
}

GSE224463_gene_info <- readr::read_csv(
  GSE224463_gene_info_file,
  show_col_types = FALSE,
  progress = FALSE,
  name_repair = "unique"
)

GSE224463_b_number_column_candidates <- c(
  "...1",
  "locus_tag",
  "b_number",
  "bnum"
)

GSE224463_gene_symbol_column_candidates <- c(
  "gene_name",
  "Gene_name",
  "Gene_Name",
  "symbol",
  "Symbol",
  "gene"
)

GSE224463_available_b_number_columns <- intersect(
  GSE224463_b_number_column_candidates,
  colnames(GSE224463_gene_info)
)

GSE224463_available_gene_symbol_columns <- intersect(
  GSE224463_gene_symbol_column_candidates,
  colnames(GSE224463_gene_info)
)

if (length(GSE224463_available_b_number_columns) == 0L ||
    length(GSE224463_available_gene_symbol_columns) == 0L) {
  stop(
    "Could not identify b-number and gene-symbol columns in:\n",
    GSE224463_gene_info_file,
    "\nAvailable columns:\n",
    paste(colnames(GSE224463_gene_info), collapse = ", "),
    call. = FALSE
  )
}

GSE224463_b_number_column <-
  GSE224463_available_b_number_columns[[1L]]

GSE224463_gene_symbol_column <-
  GSE224463_available_gene_symbol_columns[[1L]]

GSE224463_b_number_to_symbol <- GSE224463_gene_info %>%
  dplyr::transmute(
    b_number = stringr::str_extract(
      tolower(trimws(as.character(.data[[GSE224463_b_number_column]]))),
      "b\\d{4}"
    ),
    gene_symbol = toupper(
      stringr::str_replace(
        trimws(as.character(.data[[GSE224463_gene_symbol_column]])),
        "_ECOLI$",
        ""
      )
    )
  ) %>%
  dplyr::filter(
    !is.na(.data$b_number),
    !is.na(.data$gene_symbol),
    .data$gene_symbol != ""
  ) %>%
  dplyr::distinct(
    .data$b_number,
    .keep_all = TRUE
  ) %>%
  dplyr::arrange(.data$b_number)

if (nrow(GSE224463_b_number_to_symbol) < 3500L) {
  stop(
    "Too few MG1655 b-number-to-symbol mappings were recovered. ",
    "Do not run GSEA until the annotation source is checked.",
    call. = FALSE
  )
}

readr::write_tsv(
  GSE224463_b_number_to_symbol,
  file.path(
    GSE224463_gsea_dir,
    "GSE224463_MG1655_b_number_to_gene_symbol_mapping.tsv"
  )
)

GSE224463_GO_terms <- c(
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

GSE224463_GO_names <- c(
  "DNA replication", "Cell division", "Cell cycle", "Ribosome biogenesis",
  "Translation", "Purine ribonucleotide biosynthesis",
  "Pyrimidine ribonucleotide biosynthesis", "Amino acid biosynthesis",
  "DNA-binding TF activity", "Antibiotic metabolic process", "SOS response",
  "Response to stress", "DNA replication initiation",
  "Cell wall organization", "ATP synthesis coupled proton transport",
  "Aerobic respiration", "Glycolysis", "Translation initiation",
  "Response to antibiotic", "Stringent response", "tRNA aminoacylation",
  "Peptidoglycan biosynthesis", "Regulation of DNA replication",
  "DNA topology", "Cell wall macromolecule biosynthetic process",
  "Translation elongation", "Regulation of translation",
  "Chromosome segregation", "Cytokinesis", "Gluconeogenesis", "TCA cycle",
  "Oxidative phosphorylation", "Fermentation", "ATP metabolic process",
  "Fatty acid biosynthesis", "Fatty acid beta-oxidation",
  "Aminoacyl-tRNA ligase activity", "Ribonucleoprotein complex",
  "Lipopolysaccharide biosynthesis", "Lipid A biosynthesis",
  "Menaquinone biosynthesis", "Pyridoxine biosynthesis",
  "DNA biosynthetic process", "ATP biosynthetic process", "DNA repair",
  "DNA damage response", "Photoreactive repair", "Base-excision repair",
  "BER: AP site formation", "BER: gap-filling",
  "Nucleotide-excision repair", "NER: DNA incision", "NER: DNA gap filling",
  "Mismatch repair", "Recombinational repair",
  "Double-strand break repair",
  "DSB repair via homologous recombination",
  "DSB repair via nonhomologous end joining",
  "DNA alkylation repair", "Interstrand cross-link repair",
  "Translesion synthesis", "Error-prone translesion synthesis",
  "Error-free translesion synthesis",
  "DNA synthesis involved in DNA repair"
)

if (length(GSE224463_GO_terms) != length(GSE224463_GO_names)) {
  stop(
    "The curated GO-ID and GO-name vectors have different lengths.",
    call. = FALSE
  )
}

GSE224463_GO_panel <- tibble::tibble(
  GO = GSE224463_GO_terms,
  term = GSE224463_GO_names
) %>%
  dplyr::mutate(
    process_class = dplyr::case_when(
      .data$term %in% c(
        "DNA replication",
        "DNA replication initiation",
        "DNA biosynthetic process",
        "Regulation of DNA replication",
        "DNA topology",
        "Cell division",
        "Cell cycle",
        "Chromosome segregation",
        "Cytokinesis"
      ) ~ "Growth and chromosome biology",

      .data$term %in% c(
        "Ribosome biogenesis",
        "Translation",
        "Translation initiation",
        "Translation elongation",
        "Regulation of translation",
        "tRNA aminoacylation",
        "Aminoacyl-tRNA ligase activity",
        "Ribonucleoprotein complex",
        "DNA-binding TF activity"
      ) ~ "Translation and transcriptional regulation",

      .data$term %in% c(
        "Purine ribonucleotide biosynthesis",
        "Pyrimidine ribonucleotide biosynthesis",
        "Amino acid biosynthesis",
        "Fatty acid biosynthesis",
        "Fatty acid beta-oxidation",
        "Peptidoglycan biosynthesis",
        "Cell wall organization",
        "Cell wall macromolecule biosynthetic process",
        "Lipopolysaccharide biosynthesis",
        "Lipid A biosynthesis",
        "Menaquinone biosynthesis",
        "Pyridoxine biosynthesis"
      ) ~ "Biosynthesis and envelope metabolism",

      .data$term %in% c(
        "Aerobic respiration",
        "Glycolysis",
        "Gluconeogenesis",
        "TCA cycle",
        "Oxidative phosphorylation",
        "Fermentation",
        "ATP metabolic process",
        "ATP biosynthetic process",
        "ATP synthesis coupled proton transport"
      ) ~ "Central and energy metabolism",

      .data$term %in% c(
        "SOS response",
        "Response to stress",
        "Response to antibiotic",
        "Stringent response",
        "Antibiotic metabolic process"
      ) ~ "Stress-response programmes",

      TRUE ~ "DNA repair and genome maintenance"
    )
  )

readr::write_tsv(
  GSE224463_GO_panel,
  file.path(
    GSE224463_gsea_dir,
    "GSE224463_curated_GO_panel.tsv"
  )
)

GSE224463_TERM2GENE_cache <- file.path(
  GSE224463_gsea_dir,
  "GSE224463_GO_panel_TERM2GENE_E_coli_K12.tsv"
)

if (file.exists(GSE224463_TERM2GENE_cache) &&
    !GSE224463_refresh_GO_cache) {

  GSE224463_GO_mapping_raw <- readr::read_tsv(
    GSE224463_TERM2GENE_cache,
    show_col_types = FALSE,
    progress = FALSE
  ) %>%
    dplyr::transmute(
      GO = as.character(.data$GO),
      SYMBOL = as.character(.data$SYMBOL)
    )

  GSE224463_GO_retrieval_status <- GSE224463_GO_panel %>%
    dplyr::left_join(
      GSE224463_GO_mapping_raw %>%
        dplyr::count(.data$GO, name = "N_retrieved_symbols"),
      by = "GO"
    ) %>%
    dplyr::mutate(
      N_retrieved_symbols = dplyr::coalesce(
        .data$N_retrieved_symbols,
        0L
      ),
      retrieval_status = "Loaded from local GO cache"
    )

} else {

  GSE224463_GO_gene_tables <- vector(
    "list",
    length(GSE224463_GO_terms)
  )

  GSE224463_GO_status_tables <- vector(
    "list",
    length(GSE224463_GO_terms)
  )

  for (GSE224463_GO_index in seq_along(GSE224463_GO_terms)) {

    GSE224463_one_GO <- GSE224463_GO_terms[GSE224463_GO_index]

    GSE224463_one_GO_url <- paste0(
      "https://api.geneontology.org/api/bioentity/function/",
      GSE224463_one_GO,
      "/genes?taxon=NCBITaxon:83333"
    )

    GSE224463_one_GO_response <- tryCatch(
      httr::GET(
        GSE224463_one_GO_url,
        httr::timeout(60)
      ),
      error = function(error_message) NULL
    )

    if (is.null(GSE224463_one_GO_response)) {

      GSE224463_GO_gene_tables[[GSE224463_GO_index]] <- tibble::tibble(
        GO = character(0),
        SYMBOL = character(0)
      )

      GSE224463_GO_status_tables[[GSE224463_GO_index]] <- tibble::tibble(
        GO = GSE224463_one_GO,
        N_retrieved_symbols = 0L,
        retrieval_status = "HTTP request failed"
      )

      next
    }

    if (httr::status_code(GSE224463_one_GO_response) != 200L) {

      GSE224463_GO_gene_tables[[GSE224463_GO_index]] <- tibble::tibble(
        GO = character(0),
        SYMBOL = character(0)
      )

      GSE224463_GO_status_tables[[GSE224463_GO_index]] <- tibble::tibble(
        GO = GSE224463_one_GO,
        N_retrieved_symbols = 0L,
        retrieval_status = paste0(
          "HTTP status ",
          httr::status_code(GSE224463_one_GO_response)
        )
      )

      next
    }

    GSE224463_one_GO_json <- tryCatch(
      jsonlite::fromJSON(
        httr::content(
          GSE224463_one_GO_response,
          as = "text",
          encoding = "UTF-8"
        )
      ),
      error = function(error_message) NULL
    )

    if (is.null(GSE224463_one_GO_json) ||
        is.null(GSE224463_one_GO_json$associations) ||
        length(GSE224463_one_GO_json$associations) == 0L) {

      GSE224463_GO_gene_tables[[GSE224463_GO_index]] <- tibble::tibble(
        GO = character(0),
        SYMBOL = character(0)
      )

      GSE224463_GO_status_tables[[GSE224463_GO_index]] <- tibble::tibble(
        GO = GSE224463_one_GO,
        N_retrieved_symbols = 0L,
        retrieval_status = "No E. coli K-12 associations returned"
      )

      next
    }

    GSE224463_one_GO_symbols <- unique(
      as.character(
        GSE224463_one_GO_json$associations$subject$label
      )
    )

    GSE224463_GO_gene_tables[[GSE224463_GO_index]] <- tibble::tibble(
      GO = GSE224463_one_GO,
      SYMBOL = GSE224463_one_GO_symbols
    )

    GSE224463_GO_status_tables[[GSE224463_GO_index]] <- tibble::tibble(
      GO = GSE224463_one_GO,
      N_retrieved_symbols = length(GSE224463_one_GO_symbols),
      retrieval_status = "Retrieved from GO API"
    )
  }

  GSE224463_GO_mapping_raw <- dplyr::bind_rows(
    GSE224463_GO_gene_tables
  ) %>%
    dplyr::mutate(
      SYMBOL = stringr::str_trim(as.character(.data$SYMBOL)),
      SYMBOL = stringr::str_replace(.data$SYMBOL, "_ECOLI$", ""),
      SYMBOL = toupper(.data$SYMBOL)
    ) %>%
    dplyr::filter(
      !is.na(.data$GO),
      .data$GO != "",
      !is.na(.data$SYMBOL),
      .data$SYMBOL != ""
    ) %>%
    dplyr::distinct(
      .data$GO,
      .data$SYMBOL
    ) %>%
    dplyr::arrange(
      .data$GO,
      .data$SYMBOL
    )

  GSE224463_GO_retrieval_status <- dplyr::bind_rows(
    GSE224463_GO_status_tables
  )

  readr::write_tsv(
    GSE224463_GO_mapping_raw,
    GSE224463_TERM2GENE_cache
  )
}

GSE224463_GO_genes_all <- GSE224463_GO_mapping_raw %>%
  dplyr::left_join(
    GSE224463_GO_panel,
    by = "GO"
  ) %>%
  dplyr::filter(
    !is.na(.data$GO),
    !is.na(.data$term),
    !is.na(.data$SYMBOL),
    .data$SYMBOL != ""
  ) %>%
  dplyr::distinct(
    .data$GO,
    .data$term,
    .data$SYMBOL
  ) %>%
  dplyr::arrange(
    .data$term,
    .data$SYMBOL
  )

GSE224463_GO_term_coverage <- GSE224463_GO_panel %>%
  dplyr::left_join(
    GSE224463_GO_genes_all %>%
      dplyr::count(.data$GO, name = "N_genes_in_TERM2GENE"),
    by = "GO"
  ) %>%
  dplyr::left_join(
    GSE224463_GO_retrieval_status %>%
      dplyr::select(
        .data$GO,
        .data$N_retrieved_symbols,
        .data$retrieval_status
      ),
    by = "GO"
  ) %>%
  dplyr::mutate(
    N_genes_in_TERM2GENE = dplyr::coalesce(
      .data$N_genes_in_TERM2GENE,
      0L
    ),
    N_retrieved_symbols = dplyr::coalesce(
      .data$N_retrieved_symbols,
      0L
    ),
    retrieval_status = dplyr::coalesce(
      .data$retrieval_status,
      "No retrieval-status record"
    )
  ) %>%
  dplyr::arrange(.data$term)

GSE224463_TERM2GENE <- GSE224463_GO_genes_all %>%
  dplyr::transmute(
    term = as.character(.data$term),
    gene = toupper(trimws(as.character(.data$SYMBOL)))
  ) %>%
  dplyr::filter(
    !is.na(.data$term),
    .data$term != "",
    !is.na(.data$gene),
    .data$gene != ""
  ) %>%
  dplyr::distinct()

if (nrow(GSE224463_TERM2GENE) == 0L ||
    dplyr::n_distinct(GSE224463_TERM2GENE$term) < 10L) {
  stop(
    "The cached/retrieved GO panel has insufficient E. coli K-12 annotation.",
    call. = FALSE
  )
}

print(GSE224463_GO_term_coverage)

readr::write_tsv(
  GSE224463_GO_term_coverage,
  file.path(
    GSE224463_gsea_dir,
    "GSE224463_GO_panel_term_annotation_coverage.tsv"
  )
)

GSE224463_GSEA_rank_table_list <- vector(
  "list",
  nrow(GSE224463_contrast_manifest)
)

GSE224463_GSEA_result_list <- vector(
  "list",
  nrow(GSE224463_contrast_manifest)
)

GSE224463_GSEA_mapping_audit_list <- vector(
  "list",
  nrow(GSE224463_contrast_manifest)
)

GSE224463_GSEA_execution_audit_list <- vector(
  "list",
  nrow(GSE224463_contrast_manifest)
)

GSE224463_empty_GSEA_result <- tibble::tibble(
  analysis_id = character(0),
  series = character(0),
  target_time_min = integer(0),
  comparison_label = character(0),
  comparison_order = integer(0),
  contrast_QC_status = character(0),
  term = character(0),
  description = character(0),
  set_size = integer(0),
  enrichment_score = numeric(0),
  NES = numeric(0),
  p_value = numeric(0),
  p_adjust_within_contrast = numeric(0),
  q_value_within_contrast = numeric(0),
  rank_at_maximum_enrichment = integer(0),
  leading_edge = character(0),
  core_enrichment = character(0)
)

for (GSE224463_contrast_index in
     seq_len(nrow(GSE224463_contrast_manifest))) {

  GSE224463_one_contrast <-
    GSE224463_contrast_manifest[GSE224463_contrast_index, ]

  GSE224463_one_analysis_id <-
    GSE224463_one_contrast$analysis_id[[1L]]

  GSE224463_one_DE_table <- GSE224463_DE_all %>%
    dplyr::filter(.data$analysis_id == GSE224463_one_analysis_id)

  GSE224463_one_rank_table_before_mapping <-
    GSE224463_one_DE_table %>%
    dplyr::transmute(
      Gene_id = tolower(trimws(as.character(.data$Gene_id))),
      log2FC_target_time_minus_0min =
        as.numeric(.data$log2FC_target_time_minus_0min),
      b_number = stringr::str_extract(
        tolower(trimws(as.character(.data$Gene_id))),
        "b\\d{4}"
      )
    ) %>%
    dplyr::filter(
      is.finite(.data$log2FC_target_time_minus_0min),
      !is.na(.data$b_number)
    )

  GSE224463_one_rank_table <-
    GSE224463_one_rank_table_before_mapping %>%
    dplyr::left_join(
      GSE224463_b_number_to_symbol,
      by = "b_number"
    ) %>%
    dplyr::filter(
      !is.na(.data$gene_symbol),
      .data$gene_symbol != ""
    ) %>%
    dplyr::group_by(.data$gene_symbol) %>%
    dplyr::slice_max(
      order_by = abs(.data$log2FC_target_time_minus_0min),
      n = 1,
      with_ties = FALSE
    ) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(
      dplyr::desc(.data$log2FC_target_time_minus_0min)
    )

  GSE224463_GSEA_mapping_audit_list[[GSE224463_contrast_index]] <-
    tibble::tibble(
      analysis_id = GSE224463_one_analysis_id,
      series = GSE224463_one_contrast$series[[1L]],
      target_time_min = GSE224463_one_contrast$target_time_min[[1L]],
      N_DESeq2_rows = nrow(GSE224463_one_DE_table),
      N_finite_log2FC_b_numbers = nrow(
        GSE224463_one_rank_table_before_mapping
      ),
      N_b_numbers_mapped_to_gene_symbols = sum(
        GSE224463_one_rank_table_before_mapping$b_number %in%
          GSE224463_b_number_to_symbol$b_number
      ),
      N_unique_gene_symbols_ranked = nrow(GSE224463_one_rank_table),
      fraction_b_numbers_mapped_to_symbols =
        sum(
          GSE224463_one_rank_table_before_mapping$b_number %in%
            GSE224463_b_number_to_symbol$b_number
        ) /
        nrow(GSE224463_one_rank_table_before_mapping)
    )

  readr::write_tsv(
    GSE224463_one_rank_table,
    file.path(
      GSE224463_gsea_rank_dir,
      paste0(
        GSE224463_one_analysis_id,
        "_GSEA_ranked_genes.tsv"
      )
    )
  )

  GSE224463_one_gene_ranks <-
    GSE224463_one_rank_table$log2FC_target_time_minus_0min

  names(GSE224463_one_gene_ranks) <-
    GSE224463_one_rank_table$gene_symbol

  GSE224463_one_gene_ranks <- sort(
    GSE224463_one_gene_ranks,
    decreasing = TRUE
  )

  if (length(GSE224463_one_gene_ranks) < 30L) {

    GSE224463_GSEA_result_list[[GSE224463_contrast_index]] <-
      GSE224463_empty_GSEA_result

    GSE224463_GSEA_execution_audit_list[[GSE224463_contrast_index]] <-
      tibble::tibble(
        analysis_id = GSE224463_one_analysis_id,
        series = GSE224463_one_contrast$series[[1L]],
        target_time_min = GSE224463_one_contrast$target_time_min[[1L]],
        N_ranked_gene_symbols = length(GSE224463_one_gene_ranks),
        GSEA_status = "Not run: fewer than 30 ranked gene symbols",
        GSEA_error_message = NA_character_
      )

    next
  }

  set.seed(GSE224463_GSEA_seed + GSE224463_contrast_index)

  GSE224463_one_GSEA_run <- tryCatch(
    list(
      GSEA_object = clusterProfiler::GSEA(
        geneList = GSE224463_one_gene_ranks,
        TERM2GENE = GSE224463_TERM2GENE,
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        minGSSize = 5,
        maxGSSize = 1000,
        verbose = FALSE
      ),
      GSEA_error_message = NA_character_
    ),
    error = function(error_message) {
      list(
        GSEA_object = NULL,
        GSEA_error_message = conditionMessage(error_message)
      )
    }
  )

  GSE224463_one_GSEA_table <- if (
    is.null(GSE224463_one_GSEA_run$GSEA_object)
  ) {
    data.frame()
  } else {
    as.data.frame(GSE224463_one_GSEA_run$GSEA_object)
  }

  if (nrow(GSE224463_one_GSEA_table) == 0L) {

    GSE224463_one_GSEA_result <-
      GSE224463_empty_GSEA_result

    GSE224463_one_GSEA_status <- if (
      is.na(GSE224463_one_GSEA_run$GSEA_error_message)
    ) {
      "Completed: no GO terms returned"
    } else {
      "Failed"
    }

  } else {

    GSE224463_required_GSEA_output_columns <- c(
      "ID",
      "setSize",
      "enrichmentScore",
      "NES",
      "pvalue",
      "p.adjust",
      "rank"
    )

    GSE224463_missing_GSEA_output_columns <- setdiff(
      GSE224463_required_GSEA_output_columns,
      colnames(GSE224463_one_GSEA_table)
    )

    if (length(GSE224463_missing_GSEA_output_columns) > 0L) {
      stop(
        "clusterProfiler returned a GSEA table lacking required column(s): ",
        paste(
          GSE224463_missing_GSEA_output_columns,
          collapse = ", "
        ),
        call. = FALSE
      )
    }

    GSE224463_one_GSEA_description <- if (
      "Description" %in% colnames(GSE224463_one_GSEA_table)
    ) {
      as.character(GSE224463_one_GSEA_table$Description)
    } else {
      as.character(GSE224463_one_GSEA_table$ID)
    }

    GSE224463_one_GSEA_q_values <- if (
      "qvalues" %in% colnames(GSE224463_one_GSEA_table)
    ) {
      as.numeric(GSE224463_one_GSEA_table$qvalues)
    } else {
      rep(NA_real_, nrow(GSE224463_one_GSEA_table))
    }

    GSE224463_one_GSEA_leading_edge <- if (
      "leading_edge" %in% colnames(GSE224463_one_GSEA_table)
    ) {
      as.character(GSE224463_one_GSEA_table$leading_edge)
    } else {
      rep(NA_character_, nrow(GSE224463_one_GSEA_table))
    }

    GSE224463_one_GSEA_core_enrichment <- if (
      "core_enrichment" %in% colnames(GSE224463_one_GSEA_table)
    ) {
      as.character(GSE224463_one_GSEA_table$core_enrichment)
    } else {
      rep(NA_character_, nrow(GSE224463_one_GSEA_table))
    }

    GSE224463_one_GSEA_result <-
      tibble::tibble(
        analysis_id = GSE224463_one_analysis_id,
        series = GSE224463_one_contrast$series[[1L]],
        target_time_min =
          GSE224463_one_contrast$target_time_min[[1L]],
        comparison_label =
          GSE224463_one_contrast$comparison_label[[1L]],
        comparison_order =
          GSE224463_one_contrast$comparison_order[[1L]],
        contrast_QC_status =
          GSE224463_one_contrast$contrast_QC_status[[1L]],
        term = as.character(GSE224463_one_GSEA_table$ID),
        description = GSE224463_one_GSEA_description,
        set_size = as.integer(GSE224463_one_GSEA_table$setSize),
        enrichment_score = as.numeric(
          GSE224463_one_GSEA_table$enrichmentScore
        ),
        NES = as.numeric(GSE224463_one_GSEA_table$NES),
        p_value = as.numeric(GSE224463_one_GSEA_table$pvalue),
        p_adjust_within_contrast = as.numeric(
          GSE224463_one_GSEA_table$p.adjust
        ),
        q_value_within_contrast = GSE224463_one_GSEA_q_values,
        rank_at_maximum_enrichment = as.integer(
          GSE224463_one_GSEA_table$rank
        ),
        leading_edge = GSE224463_one_GSEA_leading_edge,
        core_enrichment = GSE224463_one_GSEA_core_enrichment
      ) %>%
      dplyr::arrange(
        dplyr::desc(abs(.data$NES)),
        .data$p_adjust_within_contrast
      )

    GSE224463_one_GSEA_status <- "Completed"
  }

  GSE224463_GSEA_result_list[[GSE224463_contrast_index]] <-
    GSE224463_one_GSEA_result

  GSE224463_GSEA_execution_audit_list[[GSE224463_contrast_index]] <-
    tibble::tibble(
      analysis_id = GSE224463_one_analysis_id,
      series = GSE224463_one_contrast$series[[1L]],
      target_time_min = GSE224463_one_contrast$target_time_min[[1L]],
      N_ranked_gene_symbols = length(GSE224463_one_gene_ranks),
      GSEA_status = GSE224463_one_GSEA_status,
      GSEA_error_message =
        GSE224463_one_GSEA_run$GSEA_error_message
    )

  readr::write_tsv(
    GSE224463_one_GSEA_result,
    file.path(
      GSE224463_gsea_result_dir,
      paste0(
        GSE224463_one_analysis_id,
        "_curated_GO_GSEA_results.tsv"
      )
    )
  )
}

GSE224463_GSEA_results <- dplyr::bind_rows(
  GSE224463_GSEA_result_list
)

GSE224463_GSEA_mapping_audit <- dplyr::bind_rows(
  GSE224463_GSEA_mapping_audit_list
)

GSE224463_GSEA_execution_audit <- dplyr::bind_rows(
  GSE224463_GSEA_execution_audit_list
)

if (nrow(GSE224463_GSEA_results) > 0L) {

  GSE224463_GSEA_results <- GSE224463_GSEA_results %>%
    dplyr::group_by(.data$series) %>%
    dplyr::mutate(
      p_adjust_across_timepoints_within_series = stats::p.adjust(
        .data$p_value,
        method = "BH"
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      p_adjust_across_all_timepoints = stats::p.adjust(
        .data$p_value,
        method = "BH"
      )
    ) %>%
    dplyr::arrange(
      .data$comparison_order,
      .data$p_adjust_within_contrast,
      dplyr::desc(abs(.data$NES))
    )
} else {
  GSE224463_GSEA_results <-
    GSE224463_empty_GSEA_result %>%
    dplyr::mutate(
      p_adjust_across_timepoints_within_series = numeric(0),
      p_adjust_across_all_timepoints = numeric(0)
    )
}

print(GSE224463_GSEA_execution_audit)
print(GSE224463_GSEA_mapping_audit)

readr::write_tsv(
  GSE224463_GSEA_mapping_audit,
  file.path(
    GSE224463_gsea_dir,
    "GSE224463_GSEA_gene_symbol_mapping_audit.tsv"
  )
)

readr::write_tsv(
  GSE224463_GSEA_execution_audit,
  file.path(
    GSE224463_gsea_dir,
    "GSE224463_GSEA_execution_audit.tsv"
  )
)

readr::write_tsv(
  GSE224463_GSEA_results,
  file.path(
    GSE224463_gsea_dir,
    "GSE224463_WT_stationary_phase_time_course_curated_GO_GSEA_results.tsv"
  )
)

GSE224463_GSEA_bubble_grid <- tidyr::expand_grid(
  analysis_id = GSE224463_contrast_manifest$analysis_id,
  term = GSE224463_GO_panel$term
) %>%
  dplyr::left_join(
    GSE224463_contrast_manifest %>%
      dplyr::select(
        .data$analysis_id,
        .data$series,
        .data$target_time_min,
        .data$comparison_label,
        .data$comparison_order,
        .data$contrast_QC_status
      ),
    by = "analysis_id"
  ) %>%
  dplyr::left_join(
    GSE224463_GO_panel %>%
      dplyr::select(
        .data$term,
        .data$process_class
      ),
    by = "term"
  )

GSE224463_GSEA_bubble_data <- GSE224463_GSEA_bubble_grid %>%
  dplyr::left_join(
    GSE224463_GSEA_results %>%
      dplyr::select(
        .data$analysis_id,
        .data$term,
        .data$set_size,
        .data$NES,
        .data$p_value,
        .data$p_adjust_within_contrast,
        .data$p_adjust_across_timepoints_within_series,
        .data$p_adjust_across_all_timepoints
      ),
    by = c("analysis_id", "term")
  ) %>%
  dplyr::mutate(
    GSEA_tested = !is.na(.data$NES),
    negative_log10_series_timecourse_FDR = dplyr::if_else(
      .data$GSEA_tested,
      -log10(
        pmax(
          .data$p_adjust_across_timepoints_within_series,
          .Machine$double.xmin
        )
      ),
      NA_real_
    ),
    series_timecourse_FDR_status = dplyr::case_when(
      !.data$GSEA_tested ~ "Not testable / not returned",
      .data$p_adjust_across_timepoints_within_series < 0.05 ~
        "FDR < 0.05",
      TRUE ~ "Not FDR significant"
    ),
    target_time_label = factor(
      paste0(.data$target_time_min, " min"),
      levels = c(
        "30 min",
        "60 min",
        "90 min",
        "120 min",
        "150 min",
        "180 min",
        "240 min",
        "300 min"
      )
    ),
    term = factor(
      .data$term,
      levels = rev(GSE224463_GO_panel$term)
    )
  ) %>%
  dplyr::arrange(
    .data$series,
    .data$comparison_order,
    .data$term
  )

readr::write_tsv(
  GSE224463_GSEA_bubble_data,
  file.path(
    GSE224463_gsea_dir,
    "GSE224463_WT_stationary_phase_time_course_curated_GO_GSEA_bubble_data.tsv"
  )
)

GSE224463_GSEA_bubble_plot <- ggplot(
  GSE224463_GSEA_bubble_data %>%
    dplyr::filter(.data$GSEA_tested),
  aes(
    x = .data$target_time_label,
    y = .data$term,
    size = .data$negative_log10_series_timecourse_FDR,
    colour = .data$NES
  )
) +
  geom_point(alpha = 0.90) +
  facet_grid(
    . ~ series,
    scales = "free_x",
    space = "free_x"
  ) +
  scale_colour_gradient2(
    low = "#2166ac",
    mid = "white",
    high = "#b2182b",
    midpoint = 0,
    name = "NES"
  ) +
  scale_size_continuous(
    range = c(0.7, 7),
    name = expression(
      -log[10]~"BH FDR across time points within series"
    )
  ) +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      colour = "black",
      size = 9
    ),
    axis.text.y = element_text(
      colour = "black",
      size = 7
    ),
    strip.text = element_text(
      size = 11,
      face = "bold"
    ),
    panel.spacing.x = grid::unit(1.0, "lines")
  ) +
  labs(
    title = "Curated GO-panel GSEA across the nutrient-depletion time course",
    subtitle = paste(
      "Positive NES: enrichment among transcripts higher than matched 0 min;",
      "negative NES: enrichment among transcripts lower than matched 0 min."
    )
  )

print(GSE224463_GSEA_bubble_plot)

ggsave(
  filename = file.path(
    GSE224463_gsea_dir,
    "GSE224463_WT_stationary_phase_time_course_curated_GO_GSEA_bubble_plot.tiff"
  ),
  plot = GSE224463_GSEA_bubble_plot,
  width = 15,
  height = 18,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

GSE224463_growth_transcription_terms <- c(
  "DNA replication",
  "Cell division",
  "Ribosome biogenesis",
  "Translation",
  "Translation initiation",
  "Translation elongation",
  "Regulation of translation",
  "tRNA aminoacylation",
  "Amino acid biosynthesis",
  "DNA-binding TF activity",
  "Stringent response",
  "Response to stress",
  "SOS response"
)

GSE224463_growth_transcription_trajectory_data <-
  GSE224463_GSEA_bubble_data %>%
  dplyr::filter(
    .data$GSEA_tested,
    as.character(.data$term) %in%
      GSE224463_growth_transcription_terms
  ) %>%
  dplyr::mutate(
    term = factor(
      as.character(.data$term),
      levels = rev(GSE224463_growth_transcription_terms)
    )
  )

GSE224463_growth_transcription_GSEA_trajectory_plot <- ggplot(
  GSE224463_growth_transcription_trajectory_data,
  aes(
    x = .data$target_time_min,
    y = .data$NES,
    group = .data$series,
    colour = .data$series
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "grey55",
    linewidth = 0.35
  ) +
  geom_line(
    linewidth = 0.55,
    alpha = 0.75
  ) +
  geom_point(
    aes(shape = .data$series_timecourse_FDR_status),
    size = 2.1
  ) +
  facet_grid(
    term ~ series,
    scales = "free_x",
    space = "free_x",
    drop = FALSE
  ) +
  scale_colour_manual(
    values = c(
      "short" = "#8cc5e3",
      "long" = "#2066a8"
    ),
    guide = "none"
  ) +
  scale_shape_manual(
    values = c(
      "FDR < 0.05" = 16,
      "Not FDR significant" = 1,
      "Not testable / not returned" = 4
    ),
    name = "Series-time-course FDR"
  ) +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(
      colour = "black",
      size = 8
    ),
    axis.text.y = element_text(
      colour = "black",
      size = 8
    ),
    strip.text.y = element_text(
      size = 7,
      angle = 0
    ),
    strip.text.x = element_text(
      size = 10,
      face = "bold"
    ),
    legend.position = "bottom"
  ) +
  labs(
    title = "Progressive growth, translation, and stress-response transitions",
    subtitle = paste(
      "Each series is analysed independently relative to its own 0-min",
      "exponential reference."
    ),
    y = "Normalised enrichment score (NES)"
  )

print(GSE224463_growth_transcription_GSEA_trajectory_plot)

ggsave(
  filename = file.path(
    GSE224463_gsea_dir,
    paste0(
      "GSE224463_WT_stationary_phase_time_course_",
      "growth_translation_stress_GSEA_trajectories.tiff"
    )
  ),
  plot = GSE224463_growth_transcription_GSEA_trajectory_plot,
  width = 15,
  height = 16,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

GSE224463_metabolic_terms <- c(
  "Purine ribonucleotide biosynthesis",
  "Pyrimidine ribonucleotide biosynthesis",
  "Glycolysis",
  "Gluconeogenesis",
  "TCA cycle",
  "Aerobic respiration",
  "Oxidative phosphorylation",
  "Fermentation",
  "ATP synthesis coupled proton transport",
  "ATP metabolic process",
  "Fatty acid biosynthesis",
  "Fatty acid beta-oxidation"
)

GSE224463_metabolic_trajectory_data <-
  GSE224463_GSEA_bubble_data %>%
  dplyr::filter(
    .data$GSEA_tested,
    as.character(.data$term) %in% GSE224463_metabolic_terms
  ) %>%
  dplyr::mutate(
    term = factor(
      as.character(.data$term),
      levels = rev(GSE224463_metabolic_terms)
    )
  )

GSE224463_metabolic_GSEA_trajectory_plot <- ggplot(
  GSE224463_metabolic_trajectory_data,
  aes(
    x = .data$target_time_min,
    y = .data$NES,
    group = .data$series,
    colour = .data$series
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "grey55",
    linewidth = 0.35
  ) +
  geom_line(
    linewidth = 0.55,
    alpha = 0.75
  ) +
  geom_point(
    aes(shape = .data$series_timecourse_FDR_status),
    size = 2.1
  ) +
  facet_grid(
    term ~ series,
    scales = "free_x",
    space = "free_x",
    drop = FALSE
  ) +
  scale_colour_manual(
    values = c(
      "short" = "#8cc5e3",
      "long" = "#2066a8"
    ),
    guide = "none"
  ) +
  scale_shape_manual(
    values = c(
      "FDR < 0.05" = 16,
      "Not FDR significant" = 1,
      "Not testable / not returned" = 4
    ),
    name = "Series-time-course FDR"
  ) +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(
      colour = "black",
      size = 8
    ),
    axis.text.y = element_text(
      colour = "black",
      size = 8
    ),
    strip.text.y = element_text(
      size = 7,
      angle = 0
    ),
    strip.text.x = element_text(
      size = 10,
      face = "bold"
    ),
    legend.position = "bottom"
  ) +
  labs(
    title = "Progressive central-metabolic and energy transitions",
    subtitle = paste(
      "Each series is analysed independently relative to its own 0-min",
      "exponential reference."
    ),
    y = "Normalised enrichment score (NES)"
  )

print(GSE224463_metabolic_GSEA_trajectory_plot)

ggsave(
  filename = file.path(
    GSE224463_gsea_dir,
    paste0(
      "GSE224463_WT_stationary_phase_time_course_",
      "metabolic_energy_GSEA_trajectories.tiff"
    )
  ),
  plot = GSE224463_metabolic_GSEA_trajectory_plot,
  width = 15,
  height = 15,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

GSE224463_analysis_objects <- list(
  source_audit_rds = GSE224463_audit_rds_file,
  WT_JH_prefix = GSE224463_confirmed_WT_prefix,
  WT_manifest = GSE224463_WT_manifest,
  contrast_manifest = GSE224463_contrast_manifest,
  thermodynamic_annotation_audit = GSE224463_thermodynamic_annotation_audit,
  DESeq2_time_course_summary = GSE224463_DE_summary,
  DESeq2_all_results = GSE224463_DE_all,
  DESeq2_thermodynamic_results = GSE224463_DE_thermodynamic,
  significant_DE_thermodynamic_transcripts =
    GSE224463_significant_DE_thermodynamic,
  significant_DE_distribution_statistics =
    GSE224463_significant_DE_distribution_statistics,
  significant_DE_abundance_fractions =
    GSE224463_significant_DE_abundance_fraction,
  weighted_sample_indices = GSE224463_weighted_sample_indices,
  weighted_timecourse_summary = GSE224463_weighted_timecourse_summary,
  full_transcriptome_weighted_delta_summary =
    GSE224463_full_transcriptome_weighted_delta_summary,
  full_transcriptome_weighted_gene_contributions =
    GSE224463_full_transcriptome_weighted_gene_contributions,
  GSEA_GO_panel = GSE224463_GO_panel,
  GSEA_TERM2GENE = GSE224463_TERM2GENE,
  GSEA_results = GSE224463_GSEA_results,
  GSEA_bubble_data = GSE224463_GSEA_bubble_data,
  GSEA_mapping_audit = GSE224463_GSEA_mapping_audit,
  GSEA_execution_audit = GSE224463_GSEA_execution_audit
)

saveRDS(
  GSE224463_analysis_objects,
  file.path(
    GSE224463_output_dir,
    "GSE224463_analysis_results_for_Fig_5.rds"
  )
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(
    GSE224463_output_dir,
    "GSE224463_WT_stationary_phase_time_course_analysis_sessionInfo.txt"
  )
)

message(
  "\nGSE224463 WT stationary-phase time-course thermodynamic analysis is complete.\n",
  "Results directory: ", GSE224463_output_dir, "\n",
  "All ten target-timepoint-minus-0-min contrasts were analysed within their ",
  "own series.\n",
  "The full-transcriptome ΔGfold results use no DE filter.\n",
  "Depth flags are descriptive only; no author-provided sample was removed."
)
