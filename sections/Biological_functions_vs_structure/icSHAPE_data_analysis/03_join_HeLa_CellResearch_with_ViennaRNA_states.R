suppressPackageStartupMessages({
  library(data.table)
  library(readr)
  library(tibble)
})

hela_mapped_file <- "path_to/HeLa_CellResearch_icSHAPE_positions_mapped_to_transcripts.tsv.gz"
vienna_state_file <- "path_to/ViennaRNA_per_position_state_table_for_icSHAPE.tsv.gz"
outdir <- "path_to/HeLa_CellResearch_ViennaRNA_joined_analysis"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

message2 <- function(...) cat(sprintf(...), "\n")

stop_if_missing <- function(x, label) {
  if (!file.exists(x)) stop(label, " not found: ", x, call. = FALSE)
}

normalize_tx_id <- function(x) sub("\\..*$", "", x)

safe_wilcox <- function(x, g) {
  out <- tryCatch(
    wilcox.test(x ~ g),
    error = function(e) NULL
  )
  if (is.null(out)) return(NA_real_)
  out$p.value
}

stop_if_missing(hela_mapped_file, "Mapped HeLa icSHAPE transcript table")
stop_if_missing(vienna_state_file, "ViennaRNA per-position state table")

message2("Reading mapped HeLa icSHAPE transcript table...")
icshape_dt <- fread(hela_mapped_file)

message2("Reading ViennaRNA per-position state table...")
vienna_dt <- fread(vienna_state_file)

message2("Rows in mapped HeLa icSHAPE table: %s", format(nrow(icshape_dt), big.mark = ","))
message2("Rows in ViennaRNA position-state table: %s", format(nrow(vienna_dt), big.mark = ","))

icshape_dt[, transcript_id_norm := normalize_tx_id(transcript_id)]

if ("ensembl_transcript_id" %in% names(vienna_dt)) {
  vienna_dt[, transcript_id_norm := normalize_tx_id(ensembl_transcript_id)]
} else if ("ensembl_transcript_id_version" %in% names(vienna_dt)) {
  vienna_dt[, transcript_id_norm := normalize_tx_id(ensembl_transcript_id_version)]
} else {
  stop("ViennaRNA state table lacks recognizable transcript ID column.", call. = FALSE)
}

common_tx <- intersect(
  unique(icshape_dt$transcript_id_norm),
  unique(vienna_dt$transcript_id_norm)
)

message2("Common transcripts between HeLa icSHAPE and ViennaRNA: %s",
         format(length(common_tx), big.mark = ","))

icshape_dt <- icshape_dt[transcript_id_norm %in% common_tx]
vienna_dt  <- vienna_dt[transcript_id_norm %in% common_tx]

message2("Rows in mapped HeLa icSHAPE table after transcript restriction: %s",
         format(nrow(icshape_dt), big.mark = ","))
message2("Rows in ViennaRNA state table after transcript restriction: %s",
         format(nrow(vienna_dt), big.mark = ","))

vienna_join <- unique(
  vienna_dt[, .(
    transcript_id_norm,
    tx_position,
    ensembl_transcript_id_version,
    ensembl_transcript_id,
    base,
    dot_bracket_symbol,
    vienna_state,
    Transcript_length,
    MFE
  )],
  by = c("transcript_id_norm", "tx_position")
)

icshape_join <- icshape_dt[, .(
  transcript_id,
  transcript_name,
  transcript_id_norm,
  gene_id,
  gene_name,
  chrom,
  genomic_pos_1based,
  strand,
  score,
  tx_position
)]

joined_dt <- merge(
  icshape_join,
  vienna_join,
  by = c("transcript_id_norm", "tx_position"),
  all = FALSE
)

setorder(joined_dt, transcript_id_norm, tx_position)

message2("Rows in joined HeLa icSHAPE-ViennaRNA table: %s",
         format(nrow(joined_dt), big.mark = ","))

if (nrow(joined_dt) == 0) {
  stop("Join produced zero rows. Check transcript IDs and tx_position definitions.", call. = FALSE)
}

joined_summary <- tibble(
  n_joined_rows = nrow(joined_dt),
  n_unique_transcripts = uniqueN(joined_dt$transcript_id_norm),
  n_unique_genes = uniqueN(joined_dt$gene_id),
  n_unique_transcript_positions = uniqueN(joined_dt[, .(transcript_id_norm, tx_position)]),
  n_paired = sum(joined_dt$vienna_state == "paired", na.rm = TRUE),
  n_unpaired = sum(joined_dt$vienna_state == "unpaired", na.rm = TRUE)
)

message2("Joined table summary:")
print(joined_summary)
cat("\n")

state_score_summary <- joined_dt[, .(
  n = .N,
  mean_score = mean(score, na.rm = TRUE),
  median_score = median(score, na.rm = TRUE),
  sd_score = sd(score, na.rm = TRUE),
  q25 = quantile(score, 0.25, na.rm = TRUE),
  q75 = quantile(score, 0.75, na.rm = TRUE)
), by = vienna_state][order(vienna_state)]

message2("HeLa icSHAPE score summary by ViennaRNA state:")
print(state_score_summary)
cat("\n")

global_wilcox_p <- safe_wilcox(joined_dt$score, joined_dt$vienna_state)

message2("Global Wilcoxon p-value (HeLa score ~ paired/unpaired): %s",
         format(global_wilcox_p, scientific = TRUE))
cat("\n")

fwrite(
  joined_dt,
  file.path(outdir, "HeLa_CellResearch_ViennaRNA_joined_per_position_table.tsv.gz"),
  sep = "\t"
)

write_csv(
  joined_summary,
  file.path(outdir, "joined_table_summary.csv")
)

fwrite(
  state_score_summary,
  file.path(outdir, "HeLa_icSHAPE_score_summary_by_Vienna_state.tsv"),
  sep = "\t"
)

writeLines(
  c(
    paste0("Common transcripts: ", format(length(common_tx), big.mark = ",")),
    paste0("Joined rows: ", format(nrow(joined_dt), big.mark = ",")),
    paste0("Global Wilcoxon p-value: ", format(global_wilcox_p, scientific = TRUE))
  ),
  file.path(outdir, "run_summary.txt")
)

message2("HeLa joined analysis complete. Outputs written to: %s", outdir)
