runtime_arguments <- commandArgs(trailingOnly = TRUE)

if (length(runtime_arguments) != 2L) {
  stop(
    paste(
      "Usage:",
      "Rscript create_GREIN_counts_RDS.R <GREIN_count_directory> <output_directory>",
      sep = "\n"
    ),
    call. = FALSE
  )
}

source_dir <- normalizePath(runtime_arguments[[1L]], mustWork = TRUE)
output_dir <- normalizePath(runtime_arguments[[2L]], mustWork = TRUE)
output_file <- file.path(output_dir, "GREIN_counts_for_Fig_3.rds")

if (!dir.exists(output_dir)) {
  stop("Output directory does not exist: ", output_dir, call. = FALSE)
}

if (file.exists(output_file)) {
  stop("Output already exists and will not be overwritten: ", output_file, call. = FALSE)
}

all_files <- list.files(
  source_dir,
  recursive = TRUE,
  full.names = TRUE,
  include.dirs = FALSE,
  no.. = TRUE
)

count_pattern <- paste0(
  "^GSE[0-9]+_(GeneLevel|TranscriptLevel)_",
  "(Raw|TPM|Normalized).*\\.csv(\\.gz)?$"
)

count_files <- all_files[
  grepl(count_pattern, basename(all_files), ignore.case = TRUE, perl = TRUE)
]

if (!length(count_files)) {
  stop("No GREIN gene- or transcript-level count files were found.", call. = FALSE)
}

count_files <- sort(normalizePath(count_files, mustWork = TRUE))
file_names <- basename(count_files)

grein_datasets <- c(
  "GSE100327", "GSE103021", "GSE105083", "GSE115550", "GSE118548",
  "GSE121105", "GSE126779", "GSE129221", "GSE129944", "GSE130161",
  "GSE135514", "GSE137535", "GSE143462", "GSE143798", "GSE148096",
  "GSE163249", "GSE165325", "GSE177054", "GSE178755", "GSE182018",
  "GSE184125", "GSE185141", "GSE241600", "GSE241735", "GSE62056",
  "GSE63582", "GSE73318", "GSE93124", "GSE99060"
)

accession <- sub("^(GSE[0-9]+)_.*$", "\\1", file_names, perl = TRUE)
level <- ifelse(
  grepl("_TranscriptLevel_", file_names, ignore.case = TRUE),
  "transcript",
  "gene"
)
measure <- ifelse(
  grepl("_Raw", file_names, ignore.case = TRUE),
  "raw_counts",
  ifelse(
    grepl("_TPM", file_names, ignore.case = TRUE),
    "tpm",
    "normalized"
  )
)
keys <- paste(accession, level, measure, sep = "__")
checksums <- unname(tools::md5sum(count_files))

candidate_manifest <- data.frame(
  key = keys,
  accession = accession,
  level = level,
  measure = measure,
  file_name = file_names,
  source_file = count_files,
  bytes = unname(file.info(count_files)$size),
  md5 = checksums,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

candidate_manifest <- candidate_manifest[
  candidate_manifest$accession %in% grein_datasets,
  ,
  drop = FALSE
]

if (!nrow(candidate_manifest)) {
  stop("No count files belonging to the GREIN master-script datasets were found.", call. = FALSE)
}

selected_rows <- integer()
duplicate_keys <- unique(candidate_manifest$key[duplicated(candidate_manifest$key)])

for (current_key in unique(candidate_manifest$key)) {
  group_rows <- which(candidate_manifest$key == current_key)

  if (length(group_rows) == 1L) {
    selected_rows <- c(selected_rows, group_rows)
    next
  }

  group_md5 <- unique(candidate_manifest$md5[group_rows])

  if (length(group_md5) != 1L) {
    details <- paste(candidate_manifest$source_file[group_rows], collapse = "\n")
    stop(
      "Different files represent the same GREIN table (", current_key,
      "). Select the correct version before continuing:\n", details,
      call. = FALSE
    )
  }

  relative_paths <- substring(
    candidate_manifest$source_file[group_rows],
    nchar(source_dir) + 2L
  )
  ordering <- order(nchar(relative_paths), relative_paths)
  selected_rows <- c(selected_rows, group_rows[ordering[[1L]]])
}

selected <- candidate_manifest[sort(selected_rows), , drop = FALSE]
selected <- selected[order(selected$accession, selected$level, selected$measure), , drop = FALSE]
rownames(selected) <- NULL

required_products <- c(
  "gene__raw_counts",
  "gene__normalized",
  "transcript__raw_counts",
  "transcript__tpm"
)

required_keys <- as.vector(
  outer(grein_datasets, required_products, paste, sep = "__")
)

missing_required <- setdiff(required_keys, selected$key)

if (length(missing_required)) {
  stop(
    paste(
      "Four-product GREIN tables are missing for master-script datasets:",
      paste(missing_required, collapse = "\n"),
      sep = "\n"
    ),
    call. = FALSE
  )
}

selected$included_in_complete_bundle <- selected$key %in% required_keys

if (nrow(selected) != length(required_keys) ||
    !setequal(selected$key, required_keys)) {
  stop("The selected GREIN data do not form a complete four-table set.", call. = FALSE)
}

cat("GREIN COUNT BUNDLE CREATOR\n\n")
cat("Source directory: ", source_dir, "\n", sep = "")
cat("New output:      ", output_file, "\n\n", sep = "")
cat("The source files are read only. They are not edited, moved, renamed, or deleted.\n")
cat("An existing output file will never be overwritten.\n\n")
cat("GREIN datasets:          ", length(grein_datasets), "\n", sep = "")
cat("Tables per dataset:      4\n")
cat("Count tables selected:   ", nrow(selected), "\n", sep = "")

if (length(duplicate_keys)) {
  cat("Identical duplicate categories resolved: ", length(duplicate_keys), "\n", sep = "")
}

read_count_table <- function(path) {
  connection <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    gzfile(path, open = "rt")
  } else {
    file(path, open = "rt")
  }
  on.exit(close(connection), add = TRUE)

  utils::read.csv(
    connection,
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    comment.char = "",
    quote = "\""
  )
}

tables <- vector("list", nrow(selected))
names(tables) <- selected$key
selected$rows <- NA_integer_
selected$columns <- NA_integer_

cat("\nReading and validating tables. This may take considerable time and memory.\n\n")

for (i in seq_len(nrow(selected))) {
  cat(sprintf("[%d/%d] %s\n", i, nrow(selected), selected$file_name[[i]]))
  current_table <- read_count_table(selected$source_file[[i]])

  if (!is.data.frame(current_table) || nrow(current_table) < 1L || ncol(current_table) < 2L) {
    stop("Invalid or empty count table: ", selected$source_file[[i]], call. = FALSE)
  }

  selected$rows[[i]] <- nrow(current_table)
  selected$columns[[i]] <- ncol(current_table)
  tables[[i]] <- current_table
}

public_manifest <- selected[
  , c(
    "key", "accession", "level", "measure", "file_name", "bytes", "md5",
    "rows", "columns", "included_in_complete_bundle"
  ),
  drop = FALSE
]

bundle <- list(
  format_version = 1L,
  source_database = "GREIN",
  grein_master_script_datasets = grein_datasets,
  products_per_dataset = c(
    gene_raw_counts = "GREIN gene-level raw counts",
    gene_normalized = "GREIN gene-level normalized CPM/TMM values",
    transcript_raw_counts = "GREIN transcript-level raw counts",
    transcript_tpm = "GREIN transcript-level TPM values"
  ),
  manifest = public_manifest,
  tables = tables
)

temporary_output <- file.path(
  output_dir,
  paste0(".", basename(output_file), ".partial_", Sys.getpid())
)

if (file.exists(temporary_output)) {
  stop("Temporary output already exists: ", temporary_output, call. = FALSE)
}

on.exit({
  if (file.exists(temporary_output)) {
    unlink(temporary_output)
  }
}, add = TRUE)

cat("\nWriting compressed RDS. This may take considerable time.\n")
saveRDS(bundle, temporary_output, compress = "gzip", version = 3L)

validation <- readRDS(temporary_output)

if (!is.list(validation) ||
    !identical(names(validation$tables), names(bundle$tables)) ||
    !identical(validation$manifest$key, bundle$manifest$key) ||
    !identical(validation$manifest$rows, bundle$manifest$rows) ||
    !identical(validation$manifest$columns, bundle$manifest$columns)) {
  stop("The new RDS failed validation.", call. = FALSE)
}

rm(validation)
gc(verbose = FALSE)

if (!file.rename(temporary_output, output_file)) {
  stop("The validated temporary RDS could not be renamed to the final output.", call. = FALSE)
}

output_size <- unname(file.info(output_file)$size)

cat("\nValidation successful.\n")
cat("Tables in new RDS: ", length(bundle$tables), "\n", sep = "")
cat("Output size:       ", format(output_size, big.mark = ",", scientific = FALSE),
    " bytes\n", sep = "")
cat("Created: ", output_file, "\n", sep = "")
cat("\nNo source file was changed or deleted.\n")
