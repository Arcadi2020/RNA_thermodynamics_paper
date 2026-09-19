suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

outdir <- "path_to/ViennaRNA_position_table_for_icSHAPE"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

rnafold_job1 <- "../../../data/Prediction_MFE_Gencode_37C.txt.gz"
rnafold_job2 <- "../../../data/GENCODE_filtered_transcripts.txt.gz"
rnafold_job3 <- "../../../data/GENCODE_filtered_transcripts_the_last_one.txt.gz"

message2 <- function(...) cat(sprintf(...), "\n")

stop_if_missing <- function(x, label) {
  if (!file.exists(x)) stop(label, " not found: ", x, call. = FALSE)
}

# The three RNAfold outputs are distributed gzipped (they compress ~87%).
# read.delim() decompresses .gz transparently, so no explicit gunzip is
# needed and the parsed result is identical to the uncompressed input.
parse_rnafold_triplets <- function(file_path) {
  x <- read.delim(
    file_path,
    header = FALSE,
    col.names = "raw",
    stringsAsFactors = FALSE
  )

  if (nrow(x) %% 3 != 0) {
    stop("RNAfold output is not a multiple of 3 rows: ", file_path, call. = FALSE)
  }

  mat <- matrix(x$raw, ncol = 3, byrow = TRUE)
  df <- as.data.frame(mat, stringsAsFactors = FALSE)
  colnames(df) <- c("header", "sequence", "structure_line")

  df$ensembl_transcript_id_version <- sub("^>", "", df$header)
  df$ensembl_transcript_id <- sub("\\.\\d+$", "", df$ensembl_transcript_id_version)

  df$structure_dotbracket <- sub("\\s+\\(.*$", "", df$structure_line)

  df$MFE <- sub("^.*\\(([-0-9\\.]+)\\)\\s*$", "\\1", df$structure_line)
  df$MFE <- suppressWarnings(as.numeric(df$MFE))

  df$Transcript_length <- nchar(df$sequence)

  bad_len <- which(nchar(df$sequence) != nchar(df$structure_dotbracket))
  if (length(bad_len) > 0) {
    message("Warning: ", length(bad_len),
            " rows have sequence/structure length mismatch in ",
            basename(file_path))
  }

  df %>%
    dplyr::select(
      ensembl_transcript_id_version,
      ensembl_transcript_id,
      sequence,
      Transcript_length,
      structure_dotbracket,
      MFE
    )
}

expand_structure_to_positions <- function(tx_id_version, tx_id_norm, sequence, structure_dotbracket, tx_length, mfe) {
  seq_vec <- strsplit(sequence, "", fixed = TRUE)[[1]]
  struct_vec <- strsplit(structure_dotbracket, "", fixed = TRUE)[[1]]

  if (length(seq_vec) != length(struct_vec)) {
    stop("Sequence and structure lengths do not match for transcript: ", tx_id_version, call. = FALSE)
  }

  data.table(
    ensembl_transcript_id_version = tx_id_version,
    ensembl_transcript_id = tx_id_norm,
    tx_position = seq_along(seq_vec),
    base = seq_vec,
    dot_bracket_symbol = struct_vec,
    vienna_state = ifelse(struct_vec == ".", "unpaired", "paired"),
    Transcript_length = tx_length,
    MFE = mfe
  )
}

stop_if_missing(rnafold_job1, "RNAfold job 1")
stop_if_missing(rnafold_job2, "RNAfold job 2")
stop_if_missing(rnafold_job3, "RNAfold job 3")

message2("Output directory: %s", outdir)

message2("Parsing RNAfold outputs to recover ViennaRNA transcript structures...")

job1_struct <- parse_rnafold_triplets(rnafold_job1)
job2_struct <- parse_rnafold_triplets(rnafold_job2)
job3_struct <- parse_rnafold_triplets(rnafold_job3)

vienna_df <- bind_rows(job1_struct, job2_struct, job3_struct) %>%
  distinct(ensembl_transcript_id_version, .keep_all = TRUE)

vienna_dt <- as.data.table(vienna_df)

message2("Recovered ViennaRNA structures for %s transcripts.",
         format(nrow(vienna_dt), big.mark = ","))

vienna_dt[, sequence_length := nchar(sequence)]
vienna_dt[, structure_length := nchar(structure_dotbracket)]

length_qc <- vienna_dt[, .(
  n_transcripts = .N,
  n_equal_length = sum(sequence_length == structure_length),
  n_mismatched_length = sum(sequence_length != structure_length)
)]

message2("Length QC:")
print(length_qc)

if (vienna_dt[sequence_length != structure_length, .N] > 0) {
  stop("Some ViennaRNA rows have mismatched sequence and structure lengths.", call. = FALSE)
}

fwrite(
  vienna_dt[, .(
    ensembl_transcript_id_version,
    ensembl_transcript_id,
    sequence,
    Transcript_length,
    structure_dotbracket,
    MFE
  )],
  file.path(outdir, "ViennaRNA_structures_recovered_for_icSHAPE.tsv"),
  sep = "\t"
)

message2("Expanding ViennaRNA structures to transcript-position state table...")

vienna_pos_list <- lapply(seq_len(nrow(vienna_dt)), function(i) {
  expand_structure_to_positions(
    tx_id_version = vienna_dt$ensembl_transcript_id_version[i],
    tx_id_norm = vienna_dt$ensembl_transcript_id[i],
    sequence = vienna_dt$sequence[i],
    structure_dotbracket = vienna_dt$structure_dotbracket[i],
    tx_length = vienna_dt$Transcript_length[i],
    mfe = vienna_dt$MFE[i]
  )
})

vienna_pos_dt <- rbindlist(vienna_pos_list, use.names = TRUE)

setorder(vienna_pos_dt, ensembl_transcript_id, tx_position)

message2("Rows in ViennaRNA position table: %s", format(nrow(vienna_pos_dt), big.mark = ","))
message2("Unique transcripts in ViennaRNA position table: %s", format(uniqueN(vienna_pos_dt$ensembl_transcript_id), big.mark = ","))

state_summary <- vienna_pos_dt[, .N, by = vienna_state][order(vienna_state)]
base_summary <- vienna_pos_dt[, .N, by = base][order(base)]

transcript_summary <- vienna_pos_dt[, .(
  Transcript_length = unique(Transcript_length)[1],
  MFE = unique(MFE)[1],
  n_paired = sum(vienna_state == "paired"),
  n_unpaired = sum(vienna_state == "unpaired"),
  frac_paired = mean(vienna_state == "paired"),
  frac_unpaired = mean(vienna_state == "unpaired")
), by = .(ensembl_transcript_id, ensembl_transcript_id_version)]

message2("ViennaRNA state summary:")
print(state_summary)
cat("\n")

message2("ViennaRNA base summary:")
print(base_summary)
cat("\n")

fwrite(
  vienna_pos_dt,
  file.path(outdir, "ViennaRNA_per_position_state_table_for_icSHAPE.tsv.gz"),
  sep = "\t"
)

fwrite(
  state_summary,
  file.path(outdir, "ViennaRNA_state_summary_for_icSHAPE.tsv"),
  sep = "\t"
)

fwrite(
  base_summary,
  file.path(outdir, "ViennaRNA_base_summary_for_icSHAPE.tsv"),
  sep = "\t"
)

fwrite(
  transcript_summary,
  file.path(outdir, "ViennaRNA_transcript_summary_for_icSHAPE.tsv"),
  sep = "\t"
)

message2("ViennaRNA position-state table written.")
