suppressPackageStartupMessages({
  library(data.table)
  library(readr)
  library(tibble)
  library(GenomicRanges)
  library(IRanges)
  library(rtracklayer)
})

icshape_clean_file <- "path_to/K562_CellResearch_icSHAPE_bed_clean_sorted.tsv.gz"
gencode_gtf <- "../../../data/upstream/gencode/gencode.v47.annotation.gtf"
output_dir <- "path_to/K562_CellResearch_icSHAPE_mapped_to_transcripts"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(icshape_clean_file)) {
  stop("Cleaned icSHAPE file not found: ", icshape_clean_file)
}

if (!file.exists(gencode_gtf)) {
  stop("GENCODE GTF not found: ", gencode_gtf)
}

icshape_dt <- fread(icshape_clean_file)

cat("Rows in cleaned icSHAPE table:", nrow(icshape_dt), "\n")

icshape_gr <- GRanges(
  seqnames = icshape_dt$chrom,
  ranges   = IRanges(
    start = icshape_dt$genomic_pos_1based,
    end   = icshape_dt$genomic_pos_1based
  ),
  strand   = icshape_dt$strand
)

mcols(icshape_gr)$score <- icshape_dt$score
mcols(icshape_gr)$bed_start <- icshape_dt$start
mcols(icshape_gr)$bed_end <- icshape_dt$end

gtf_gr <- import(gencode_gtf)

cat("Rows in full GTF:", length(gtf_gr), "\n")

gtf_seqlevels <- unique(as.character(seqnames(gtf_gr)))
icshape_seqlevels <- unique(icshape_dt$chrom)

gtf_has_chr <- any(grepl("^chr", gtf_seqlevels))
icshape_has_chr <- any(grepl("^chr", icshape_seqlevels))

if (gtf_has_chr && !icshape_has_chr) {
  seqlevels(gtf_gr) <- paste0("chr", seqlevels(gtf_gr))
}

if (!gtf_has_chr && icshape_has_chr) {
  seqlevels(gtf_gr) <- sub("^chr", "", seqlevels(gtf_gr))
}

exon_gr <- gtf_gr[gtf_gr$type == "exon"]

cat("Rows in exon subset:", length(exon_gr), "\n")

exon_gr <- exon_gr[!is.na(mcols(exon_gr)$transcript_id)]

cat("Rows in exon subset with transcript_id:", length(exon_gr), "\n")

exon_dt <- as.data.table(as.data.frame(exon_gr))

wanted_cols <- c(
  "seqnames", "start", "end", "strand",
  "gene_id", "gene_name",
  "transcript_id", "transcript_name",
  "exon_number"
)

present_cols <- intersect(wanted_cols, names(exon_dt))
exon_dt <- exon_dt[, ..present_cols]

setnames(exon_dt, old = "seqnames", new = "chrom")

if ("exon_number" %in% names(exon_dt)) {
  suppressWarnings(exon_dt[, exon_number := as.integer(exon_number)])
}

plus_dt <- exon_dt[strand == "+"]
minus_dt <- exon_dt[strand == "-"]

if (nrow(plus_dt) > 0) {
  setorder(plus_dt, transcript_id, start, end)
}

if (nrow(minus_dt) > 0) {
  setorder(minus_dt, transcript_id, -start, -end)
}

exon_dt <- rbindlist(list(plus_dt, minus_dt), use.names = TRUE)

exon_dt[, exon_width := end - start + 1L]
exon_dt[, exon_tx_start := cumsum(c(0L, head(exon_width, -1))) + 1L, by = transcript_id]
exon_dt[, exon_tx_end := exon_tx_start + exon_width - 1L]

cat("Rows in transcript exon table:", nrow(exon_dt), "\n")
cat("Unique transcripts in exon table:", uniqueN(exon_dt$transcript_id), "\n")

fwrite(
  exon_dt,
  file.path(output_dir, "transcript_exon_coordinate_table.tsv.gz"),
  sep = "\t"
)

tx_exon_gr <- GRanges(
  seqnames = exon_dt$chrom,
  ranges   = IRanges(start = exon_dt$start, end = exon_dt$end),
  strand   = exon_dt$strand
)

mcols(tx_exon_gr)$gene_id         <- if ("gene_id" %in% names(exon_dt)) exon_dt$gene_id else NA_character_
mcols(tx_exon_gr)$gene_name       <- if ("gene_name" %in% names(exon_dt)) exon_dt$gene_name else NA_character_
mcols(tx_exon_gr)$transcript_id   <- exon_dt$transcript_id
mcols(tx_exon_gr)$transcript_name <- if ("transcript_name" %in% names(exon_dt)) exon_dt$transcript_name else NA_character_
mcols(tx_exon_gr)$exon_number     <- if ("exon_number" %in% names(exon_dt)) exon_dt$exon_number else NA_integer_
mcols(tx_exon_gr)$exon_start      <- exon_dt$start
mcols(tx_exon_gr)$exon_end        <- exon_dt$end
mcols(tx_exon_gr)$exon_tx_start   <- exon_dt$exon_tx_start
mcols(tx_exon_gr)$exon_tx_end     <- exon_dt$exon_tx_end

hits <- findOverlaps(icshape_gr, tx_exon_gr, ignore.strand = FALSE)

cat("Number of overlaps between icSHAPE positions and exons:", length(hits), "\n")

if (length(hits) == 0) {
  stop("No overlaps found. Check chromosome naming and annotation build.")
}

q_idx <- queryHits(hits)
s_idx <- subjectHits(hits)

ic_hits_dt <- as.data.table(as.data.frame(icshape_gr[q_idx]))
ex_hits_dt <- as.data.table(as.data.frame(mcols(tx_exon_gr))[s_idx, , drop = FALSE])

mapped_dt <- data.table(
  chrom = as.character(ic_hits_dt$seqnames),
  genomic_pos_1based = ic_hits_dt$start,
  strand = as.character(ic_hits_dt$strand),
  score = ic_hits_dt$score,
  bed_start = ic_hits_dt$bed_start,
  bed_end = ic_hits_dt$bed_end,
  gene_id = ex_hits_dt$gene_id,
  gene_name = ex_hits_dt$gene_name,
  transcript_id = ex_hits_dt$transcript_id,
  transcript_name = ex_hits_dt$transcript_name,
  exon_number = ex_hits_dt$exon_number,
  exon_start = ex_hits_dt$exon_start,
  exon_end = ex_hits_dt$exon_end,
  exon_tx_start = ex_hits_dt$exon_tx_start,
  exon_tx_end = ex_hits_dt$exon_tx_end
)

mapped_dt[strand == "+",
          tx_position := exon_tx_start + (genomic_pos_1based - exon_start)
]

mapped_dt[strand == "-",
          tx_position := exon_tx_start + (exon_end - genomic_pos_1based)
]

mapped_dt <- mapped_dt[
  tx_position >= exon_tx_start & tx_position <= exon_tx_end
]

setorder(mapped_dt, transcript_id, tx_position, chrom, genomic_pos_1based)

mapping_summary <- tibble(
  n_icshape_positions = nrow(icshape_dt),
  n_overlap_rows = nrow(mapped_dt),
  n_unique_transcripts = uniqueN(mapped_dt$transcript_id),
  n_unique_genes = uniqueN(mapped_dt$gene_id),
  n_unique_genomic_positions = uniqueN(mapped_dt[, .(chrom, genomic_pos_1based, strand)]),
  n_unique_transcript_positions = uniqueN(mapped_dt[, .(transcript_id, tx_position)])
)

cat("Mapping summary:\n")
print(mapping_summary)
cat("\n")

multi_map_summary <- mapped_dt[
  , .(n_transcripts = uniqueN(transcript_id)),
  by = .(chrom, genomic_pos_1based, strand)
][
  , .N, by = n_transcripts
][order(n_transcripts)]

cat("Multi-transcript mapping summary:\n")
print(as_tibble(multi_map_summary), n = 30)
cat("\n")

fwrite(
  mapped_dt,
  file.path(output_dir, "K562_CellResearch_icSHAPE_positions_mapped_to_transcripts.tsv.gz"),
  sep = "\t"
)

write_csv(mapping_summary, file.path(output_dir, "mapping_summary.csv"))
fwrite(
  multi_map_summary,
  file.path(output_dir, "multi_transcript_mapping_summary.tsv"),
  sep = "\t"
)

cat("Saved mapped transcript-level icSHAPE table to:\n",
    file.path(output_dir, "K562_CellResearch_icSHAPE_positions_mapped_to_transcripts.tsv.gz"), "\n\n")

cat("K562 transcript mapping complete.\n")
