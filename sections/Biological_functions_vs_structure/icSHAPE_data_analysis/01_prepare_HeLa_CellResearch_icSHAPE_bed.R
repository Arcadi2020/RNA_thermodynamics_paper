suppressPackageStartupMessages({
  library(data.table)
  library(readr)
  library(tibble)
})

tar_file <- "path_to/RASP_files_HeLa_icSHAPE.tar"
output_dir <- "path_to/prepared_HeLa_CellResearch_icSHAPE_bed"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(tar_file)) {
  stop("Tar file not found: ", tar_file)
}

message("Tar file found: ", tar_file)

tar_contents <- utils::untar(tar_file, list = TRUE)

if (length(tar_contents) != 1) {
  stop("Expected exactly one file in tar archive, found: ", length(tar_contents))
}

bed_member <- tar_contents[[1]]

extract_dir <- file.path(output_dir, "extracted")
dir.create(extract_dir, recursive = TRUE, showWarnings = FALSE)

bed_file <- file.path(extract_dir, bed_member)

if (!file.exists(bed_file)) {
  utils::untar(tar_file, exdir = extract_dir)
}

if (!file.exists(bed_file)) {
  stop("Extracted BED file not found: ", bed_file)
}

cat("Extracted BED file:\n", bed_file, "\n\n")

bed_dt <- fread(
  cmd = paste("gzip -dc", shQuote(bed_file)),
  header = FALSE,
  sep = "\t",
  col.names = c("chrom", "start", "end", "name", "score", "strand"),
  na.strings = c("NA", ".", "")
)

cat("Dimensions of BED table:\n")
print(dim(bed_dt))
cat("\n")

width_summary <- bed_dt[, .N, by = .(width = end - start)][order(width)]
cat("Interval width summary:\n")
print(width_summary)
cat("\n")

dup_summary <- tibble(
  n_duplicated_chr_start_end_strand =
    bed_dt[, .N, by = .(chrom, start, end, strand)][N > 1, sum(N)],
  n_unique_positions =
    nrow(unique(bed_dt[, .(chrom, start, end, strand)]))
)

cat("Duplicate-position summary:\n")
print(dup_summary)
cat("\n")

basic_summary <- tibble(
  n_rows = nrow(bed_dt),
  n_cols = ncol(bed_dt),
  n_unique_chrom = dplyr::n_distinct(bed_dt$chrom),
  n_unique_strand = dplyr::n_distinct(bed_dt$strand),
  n_missing_name = sum(is.na(bed_dt$name)),
  n_missing_score = sum(is.na(bed_dt$score)),
  n_negative_score = sum(bed_dt$score < 0, na.rm = TRUE),
  n_zero_score = sum(bed_dt$score == 0, na.rm = TRUE),
  n_positive_score = sum(bed_dt$score > 0, na.rm = TRUE)
)

cat("Basic summary:\n")
print(basic_summary)
cat("\n")

score_summary <- tibble(
  min_score = min(bed_dt$score, na.rm = TRUE),
  q01 = quantile(bed_dt$score, 0.01, na.rm = TRUE),
  q05 = quantile(bed_dt$score, 0.05, na.rm = TRUE),
  q25 = quantile(bed_dt$score, 0.25, na.rm = TRUE),
  median = median(bed_dt$score, na.rm = TRUE),
  mean = mean(bed_dt$score, na.rm = TRUE),
  q75 = quantile(bed_dt$score, 0.75, na.rm = TRUE),
  q95 = quantile(bed_dt$score, 0.95, na.rm = TRUE),
  q99 = quantile(bed_dt$score, 0.99, na.rm = TRUE),
  max_score = max(bed_dt$score, na.rm = TRUE)
)

cat("Score summary:\n")
print(score_summary)
cat("\n")

strand_counts <- as_tibble(bed_dt[, .N, by = strand][order(-N)])
cat("Strand counts:\n")
print(strand_counts, n = Inf)
cat("\n")

chrom_counts <- as_tibble(bed_dt[, .N, by = chrom][order(-N)])
cat("Chromosome counts:\n")
print(chrom_counts, n = Inf)
cat("\n")

bed_dt[, genomic_pos_1based := start + 1L]

bed_clean <- bed_dt[, .(
  chrom,
  start,
  end,
  genomic_pos_1based,
  strand,
  score
)]

setorder(bed_clean, chrom, start, end, strand)

fwrite(
  bed_clean,
  file.path(output_dir, "HeLa_CellResearch_icSHAPE_bed_clean_sorted.tsv.gz"),
  sep = "\t"
)

write_csv(basic_summary, file.path(output_dir, "basic_summary.csv"))
write_csv(score_summary, file.path(output_dir, "score_summary.csv"))
write_csv(strand_counts, file.path(output_dir, "strand_counts.csv"))
write_csv(chrom_counts, file.path(output_dir, "chrom_counts.csv"))
write_csv(as_tibble(width_summary), file.path(output_dir, "width_summary.csv"))
write_csv(dup_summary, file.path(output_dir, "dup_summary.csv"))

cat("Saved cleaned BED table to:\n",
    file.path(output_dir, "HeLa_CellResearch_icSHAPE_bed_clean_sorted.tsv.gz"), "\n\n")

cat("HeLa Cell Research BED preparation complete.\n")
