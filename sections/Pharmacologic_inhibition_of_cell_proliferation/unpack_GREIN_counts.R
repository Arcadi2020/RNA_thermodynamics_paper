# Expand GREIN_counts_for_Fig_3.rds back into the CSV files the analysis reads.
#
# create_GREIN_counts_RDS.R bundles the 116 GREIN count tables into one object
# for deposit. Nothing read it back: Cross_validation_MFE.R opens 57 of those
# tables as individual CSVs from input_dir, so a reader who downloads the
# bundle has the data but no way to feed it to the script. This is that step.
#
#   Rscript unpack_GREIN_counts.R <GREIN_counts_for_Fig_3.rds> <output_directory>
#
# Each table is written under the file name recorded in the bundle's manifest.
# The first column header is empty in the originals, which is why the analysis
# reads it back as "X"; write.csv reproduces that, so the round trip is exact.

runtime_arguments <- commandArgs(trailingOnly = TRUE)

if (length(runtime_arguments) != 2L) {
  stop("Usage: Rscript unpack_GREIN_counts.R <bundle.rds> <output_directory>",
       call. = FALSE)
}

bundle_file <- normalizePath(runtime_arguments[[1L]], mustWork = TRUE)
output_dir  <- runtime_arguments[[2L]]
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

bundle <- readRDS(bundle_file)
if (!is.list(bundle) || is.null(bundle$tables) || is.null(bundle$manifest)) {
  stop("That does not look like a GREIN count bundle.", call. = FALSE)
}

manifest <- bundle$manifest
cat("Tables in bundle: ", nrow(manifest), "\n", sep = "")
cat("Writing to:       ", normalizePath(output_dir, mustWork = TRUE), "\n\n", sep = "")

written <- 0L
for (i in seq_len(nrow(manifest))) {
  key <- manifest$key[[i]]
  target <- file.path(output_dir, manifest$file_name[[i]])
  if (file.exists(target)) {
    cat("  exists, skipping  ", manifest$file_name[[i]], "\n", sep = "")
    next
  }
  table_data <- bundle$tables[[key]]
  if (is.null(table_data)) {
    stop("Manifest names a table the bundle does not contain: ", key, call. = FALSE)
  }
  write.csv(table_data, target, row.names = FALSE, quote = TRUE, na = "")
  written <- written + 1L
  if (written %% 20L == 0L) cat("  written ", written, " of ", nrow(manifest), "\n", sep = "")
}

cat("\nWrote ", written, " file(s).\n", sep = "")
cat("The bundle was read only; no source file was changed.\n")
