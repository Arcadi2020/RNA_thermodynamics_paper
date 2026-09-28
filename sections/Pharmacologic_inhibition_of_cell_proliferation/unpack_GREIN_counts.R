# Expand the deposited GREIN archives into the flat directory the analysis reads.
#
# Figure 3's GREIN inputs are deposited on Zenodo as one zip per GEO accession
#
#   GSE<accession>.zip -> GSE<accession>/GSE<accession>_{GeneLevel,TranscriptLevel}_*.csv
#
# plus GREIN_supporting_inputs.zip, which holds the per-series metadata objects
# (Data_GSE*.rds and Expdata_GSE182018_raw_gene.rds) and the one accession whose
# counts are a raw NCBI table (GSE192446_raw_counts_GRCh38.p13_NCBI.tsv).
#
# Cross_validation_MFE.R reads every one of these as a flat file from input_dir
# -- file.path(input_dir, "GSE182018_TranscriptLevel_Raw_data.csv"),
# file.path(input_dir, "Data_GSE182018.rds"), and so on. So this expands the
# archives and flattens them: each member is written to the output directory
# under its base name, dropping the per-accession folder and ignoring the
# __MACOSX / .DS_Store entries some of the archives carry.
#
#   Rscript unpack_GREIN_counts.R <zip_directory> <output_directory>
#
# <zip_directory> is where download_zenodo_data.sh placed the archives (its
# data/ directory). Files that already exist in the output directory are left
# untouched, so re-running is cheap and never overwrites. The archives are read
# only; nothing in <zip_directory> is changed.

runtime_arguments <- commandArgs(trailingOnly = TRUE)

if (length(runtime_arguments) != 2L) {
  stop("Usage: Rscript unpack_GREIN_counts.R <zip_directory> <output_directory>",
       call. = FALSE)
}

zip_dir    <- normalizePath(runtime_arguments[[1L]], mustWork = TRUE)
output_dir <- runtime_arguments[[2L]]
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
output_dir <- normalizePath(output_dir, mustWork = TRUE)

archives <- list.files(zip_dir, pattern = "\\.zip$", full.names = TRUE)
grein_pattern <- "^(GSE[0-9]+|GREIN_supporting_inputs)\\.zip$"
archives <- archives[grepl(grein_pattern, basename(archives), perl = TRUE)]
archives <- sort(archives)

if (!length(archives)) {
  stop("No GREIN archives (GSE*.zip or GREIN_supporting_inputs.zip) in ", zip_dir,
       call. = FALSE)
}

cat("GREIN archives found: ", length(archives), "\n", sep = "")
cat("Writing to:           ", output_dir, "\n\n", sep = "")

is_cruft <- function(paths) {
  grepl("(^|/)__MACOSX(/|$)", paths) | grepl("(^|/)\\.DS_Store$", paths)
}

written <- 0L
skipped <- 0L
for (archive in archives) {
  members <- utils::unzip(archive, list = TRUE)$Name
  members <- members[!grepl("/$", members)]     # drop directory entries
  members <- members[!is_cruft(members)]        # drop macOS cruft

  for (member in members) {
    target <- file.path(output_dir, basename(member))
    if (file.exists(target)) {
      skipped <- skipped + 1L
      next
    }
    tmp <- tempfile()
    dir.create(tmp)
    utils::unzip(archive, files = member, exdir = tmp, junkpaths = TRUE)
    extracted <- file.path(tmp, basename(member))
    if (!file.exists(extracted)) {
      unlink(tmp, recursive = TRUE)
      stop("Failed to extract ", member, " from ", basename(archive), call. = FALSE)
    }
    if (!file.copy(extracted, target)) {
      unlink(tmp, recursive = TRUE)
      stop("Failed to write ", target, call. = FALSE)
    }
    unlink(tmp, recursive = TRUE)
    written <- written + 1L
    if (written %% 20L == 0L) cat("  written ", written, " file(s)\n", sep = "")
  }
}

cat("\nWrote ", written, " file(s); ", skipped, " already present and left as-is.\n",
    sep = "")
cat("The archives were read only; nothing in the zip directory was changed.\n")
