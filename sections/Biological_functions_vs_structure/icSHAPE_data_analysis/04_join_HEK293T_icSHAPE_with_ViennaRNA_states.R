suppressPackageStartupMessages({
  library(data.table)
  library(readr)
  library(tibble)
  library(ggplot2)
})

icshape_mapped_file <- "path_to/HEK293T_icSHAPE_positions_mapped_to_transcripts.tsv.gz"
vienna_state_file <- "path_to/ViennaRNA_per_position_state_table_for_icSHAPE.tsv.gz"
outdir <- "path_to/HEK293T_icSHAPE_ViennaRNA_joined_analysis"
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

stop_if_missing(icshape_mapped_file, "Mapped icSHAPE transcript table")
stop_if_missing(vienna_state_file, "ViennaRNA per-position state table")

message2("Reading mapped icSHAPE transcript table...")
icshape_dt <- fread(icshape_mapped_file)

message2("Reading ViennaRNA per-position state table...")
vienna_dt <- fread(vienna_state_file)

message2("Rows in mapped icSHAPE table: %s", format(nrow(icshape_dt), big.mark = ","))
message2("Rows in ViennaRNA position-state table: %s", format(nrow(vienna_dt), big.mark = ","))

if (!"transcript_id" %in% names(icshape_dt)) {
  stop("Mapped icSHAPE table does not contain transcript_id.", call. = FALSE)
}
icshape_dt[, transcript_id_norm := normalize_tx_id(transcript_id)]

if ("ensembl_transcript_id" %in% names(vienna_dt)) {
  vienna_dt[, transcript_id_norm := normalize_tx_id(ensembl_transcript_id)]
} else if ("ensembl_transcript_id_version" %in% names(vienna_dt)) {
  vienna_dt[, transcript_id_norm := normalize_tx_id(ensembl_transcript_id_version)]
} else {
  stop("ViennaRNA state table lacks recognizable transcript ID column.", call. = FALSE)
}

if (!"tx_position" %in% names(icshape_dt)) {
  stop("Mapped icSHAPE table does not contain tx_position.", call. = FALSE)
}

if (!"tx_position" %in% names(vienna_dt)) {
  stop("ViennaRNA state table does not contain tx_position.", call. = FALSE)
}

common_tx <- intersect(
  unique(icshape_dt$transcript_id_norm),
  unique(vienna_dt$transcript_id_norm)
)

message2("Common transcripts between icSHAPE and ViennaRNA: %s",
         format(length(common_tx), big.mark = ","))

icshape_dt <- icshape_dt[transcript_id_norm %in% common_tx]
vienna_dt  <- vienna_dt[transcript_id_norm %in% common_tx]

message2("Rows in mapped icSHAPE table after transcript restriction: %s",
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

message2("Rows in joined icSHAPE-ViennaRNA table: %s",
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

message2("icSHAPE score summary by ViennaRNA state:")
print(state_score_summary)
cat("\n")

global_wilcox_p <- safe_wilcox(joined_dt$score, joined_dt$vienna_state)

message2("Global Wilcoxon p-value (score ~ paired/unpaired): %s",
         format(global_wilcox_p, scientific = TRUE))
cat("\n")

transcript_state_summary <- joined_dt[, .(
  gene_name = unique(gene_name)[1],
  transcript_name = unique(transcript_name)[1],
  Transcript_length = unique(Transcript_length)[1],
  MFE = unique(MFE)[1],
  n_positions = .N,
  n_paired = sum(vienna_state == "paired", na.rm = TRUE),
  n_unpaired = sum(vienna_state == "unpaired", na.rm = TRUE),
  mean_score_paired = mean(score[vienna_state == "paired"], na.rm = TRUE),
  mean_score_unpaired = mean(score[vienna_state == "unpaired"], na.rm = TRUE),
  median_score_paired = median(score[vienna_state == "paired"], na.rm = TRUE),
  median_score_unpaired = median(score[vienna_state == "unpaired"], na.rm = TRUE),
  delta_mean_unpaired_minus_paired =
    mean(score[vienna_state == "unpaired"], na.rm = TRUE) -
    mean(score[vienna_state == "paired"], na.rm = TRUE),
  delta_median_unpaired_minus_paired =
    median(score[vienna_state == "unpaired"], na.rm = TRUE) -
    median(score[vienna_state == "paired"], na.rm = TRUE),
  wilcox_p = safe_wilcox(score, vienna_state)
), by = .(transcript_id_norm, transcript_id)]

setorder(transcript_state_summary, -delta_mean_unpaired_minus_paired)

message2("Rows in transcript-level summary: %s",
         format(nrow(transcript_state_summary), big.mark = ","))

gene_state_summary <- joined_dt[, .(
  n_positions = .N,
  n_transcripts = uniqueN(transcript_id_norm),
  n_paired = sum(vienna_state == "paired", na.rm = TRUE),
  n_unpaired = sum(vienna_state == "unpaired", na.rm = TRUE),
  mean_score_paired = mean(score[vienna_state == "paired"], na.rm = TRUE),
  mean_score_unpaired = mean(score[vienna_state == "unpaired"], na.rm = TRUE),
  median_score_paired = median(score[vienna_state == "paired"], na.rm = TRUE),
  median_score_unpaired = median(score[vienna_state == "unpaired"], na.rm = TRUE),
  delta_mean_unpaired_minus_paired =
    mean(score[vienna_state == "unpaired"], na.rm = TRUE) -
    mean(score[vienna_state == "paired"], na.rm = TRUE)
), by = .(gene_id, gene_name)]

setorder(gene_state_summary, -delta_mean_unpaired_minus_paired)

fwrite(
  joined_dt,
  file.path(outdir, "HEK293T_icSHAPE_ViennaRNA_joined_per_position_table.tsv.gz"),
  sep = "\t"
)

write_csv(
  joined_summary,
  file.path(outdir, "joined_table_summary.csv")
)

fwrite(
  state_score_summary,
  file.path(outdir, "icSHAPE_score_summary_by_Vienna_state.tsv"),
  sep = "\t"
)

fwrite(
  transcript_state_summary,
  file.path(outdir, "transcript_level_icSHAPE_Vienna_state_summary.tsv.gz"),
  sep = "\t"
)

fwrite(
  gene_state_summary,
  file.path(outdir, "gene_level_icSHAPE_Vienna_state_summary.tsv.gz"),
  sep = "\t"
)

p1 <- ggplot(as.data.frame(joined_dt), aes(x = vienna_state, y = score)) +
  geom_boxplot(outlier.alpha = 0.2) +
  theme_classic(base_size = 12) +
  labs(
    title = "icSHAPE score by ViennaRNA-predicted state",
    x = "ViennaRNA state",
    y = "icSHAPE score"
  )
print(p1)

p2 <- ggplot(as.data.frame(joined_dt), aes(x = score, fill = vienna_state)) +
  geom_histogram(bins = 50, alpha = 0.6, position = "identity") +
  theme_classic(base_size = 12) +
  labs(
    title = "Distribution of icSHAPE scores by ViennaRNA-predicted state",
    x = "icSHAPE score",
    y = "Number of positions",
    fill = "ViennaRNA state"
  )
print(p2)

p3 <- ggplot(as.data.frame(transcript_state_summary),
             aes(x = delta_mean_unpaired_minus_paired)) +
  geom_histogram(bins = 50) +
  theme_classic(base_size = 12) +
  labs(
    title = "Transcript-level difference in icSHAPE score",
    subtitle = "Unpaired mean score minus paired mean score",
    x = "Delta mean (unpaired - paired)",
    y = "Number of transcripts"
  )
print(p3)

summary_lines <- c(
  "icSHAPE-MaP vs ViennaRNA joined transcript-position analysis",
  "",
  paste0("Mapped icSHAPE file: ", icshape_mapped_file),
  paste0("ViennaRNA state file: ", vienna_state_file),
  paste0("Output directory: ", outdir),
  "",
  paste0("Common transcripts: ", format(length(common_tx), big.mark = ",")),
  paste0("Joined rows: ", format(nrow(joined_dt), big.mark = ",")),
  paste0("Unique joined transcripts: ", format(uniqueN(joined_dt$transcript_id_norm), big.mark = ",")),
  paste0("Unique joined genes: ", format(uniqueN(joined_dt$gene_id), big.mark = ",")),
  "",
  paste0("Paired positions: ", sum(joined_dt$vienna_state == "paired", na.rm = TRUE)),
  paste0("Unpaired positions: ", sum(joined_dt$vienna_state == "unpaired", na.rm = TRUE)),
  "",
  paste0("Mean score (paired): ",
         signif(state_score_summary[vienna_state == "paired", mean_score], 5)),
  paste0("Mean score (unpaired): ",
         signif(state_score_summary[vienna_state == "unpaired", mean_score], 5)),
  paste0("Median score (paired): ",
         signif(state_score_summary[vienna_state == "paired", median_score], 5)),
  paste0("Median score (unpaired): ",
         signif(state_score_summary[vienna_state == "unpaired", median_score], 5)),
  "",
  paste0("Global Wilcoxon p-value: ", format(global_wilcox_p, scientific = TRUE))
)

writeLines(
  summary_lines,
  con = file.path(outdir, "run_summary.txt")
)

message2("Joined analysis complete. Outputs written to: %s", outdir)
