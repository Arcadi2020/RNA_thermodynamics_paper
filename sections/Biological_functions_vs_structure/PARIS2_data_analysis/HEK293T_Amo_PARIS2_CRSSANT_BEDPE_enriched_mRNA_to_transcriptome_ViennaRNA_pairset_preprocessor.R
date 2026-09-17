suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(rtracklayer)
  library(ggplot2)
})

bedpe_file <- "../../../data/HEK293T_Amo_mRNA_enriched_CRSSANT_combined_128shards_dg.bedpe"
gtf_file <- "../../../data/upstream/gencode/gencode.v47.annotation.gtf"
vienna_master_pairs_file <- "path_to/ViennaRNA_transcript_level_pairs.tsv.gz"
outdir <- "path_to/PARIS2_enriched_mRNA_ViennaRNA_precomparison_output_HEK293T"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

min_support <- 2L
min_support_strict <- 5L
min_overlap_bp <- 5L

cat("\nSelected BEDPE file:\n")
print(bedpe_file)

message2 <- function(...) cat(sprintf(...), "\n")

stop_if_missing <- function(x, label) {
  if (!file.exists(x)) stop(label, " not found: ", x, call. = FALSE)
}

normalize_tx_id <- function(x) sub("\\..*$", "", x)

read_crssant_bedpe <- function(path) {
  stop_if_missing(path, "CRSSANT BEDPE")

  dt <- fread(
    path,
    header = FALSE,
    sep = "\t"
  )

  if (ncol(dt) != 10) {
    stop("BEDPE does not have 10 columns. Found: ", ncol(dt), call. = FALSE)
  }

  setnames(dt, c(
    "chrom1", "start1", "end1",
    "chrom2", "start2", "end2",
    "dg_info", "support",
    "strand1", "strand2"
  ))

  dg_split <- tstrsplit(dt$dg_info, ",", fixed = TRUE)

  if (length(dg_split) != 4) {
    stop("dg_info does not split into 4 fields consistently.", call. = FALSE)
  }

  dt[, gene1 := dg_split[[1]]]
  dt[, gene2 := dg_split[[2]]]
  dt[, dg_id := suppressWarnings(as.integer(dg_split[[3]]))]
  dt[, dg_score := suppressWarnings(as.numeric(dg_split[[4]]))]

  dt[, start1 := as.integer(start1)]
  dt[, end1 := as.integer(end1)]
  dt[, start2 := as.integer(start2)]
  dt[, end2 := as.integer(end2)]
  dt[, support := as.integer(support)]

  dt[, arm1_len := end1 - start1 + 1L]
  dt[, arm2_len := end2 - start2 + 1L]
  dt[, genomic_span := pmax(end1, end2) - pmin(start1, start2) + 1L]

  dt[, same_chr := chrom1 == chrom2]
  dt[, same_gene := gene1 == gene2]
  dt[, same_strand := strand1 == strand2]

  dt
}

build_exon_models <- function(gtf_path) {
  stop_if_missing(gtf_path, "GTF")

  gtf <- import(gtf_path)
  exons <- gtf[mcols(gtf)$type == "exon"]

  if (is.null(mcols(exons)$transcript_id)) {
    stop("The GTF does not contain transcript_id in exon rows.", call. = FALSE)
  }

  if (is.null(mcols(exons)$gene_name)) {
    stop("The GTF does not contain gene_name in exon rows.", call. = FALSE)
  }

  exons_dt <- as.data.table(exons)
  exons_dt[, transcript_id := normalize_tx_id(transcript_id)]

  exons_dt <- exons_dt[, .(
    seqnames = as.character(seqnames),
    start,
    end,
    strand = as.character(strand),
    transcript_id,
    gene_name = as.character(gene_name)
  )]

  setorder(exons_dt, transcript_id, start, end)

  tx_models <- exons_dt[, {
    strand_val <- unique(strand)
    if (length(strand_val) != 1L) {
      strand_val <- strand_val[1]
    }

    gene_name_val <- unique(gene_name)
    if (length(gene_name_val) != 1L) {
      gene_name_val <- gene_name_val[1]
    }

    if (strand_val == "+") {
      ord <- order(start, end)
      s <- start[ord]
      e <- end[ord]
    } else {
      ord <- order(start, end, decreasing = TRUE)
      s <- start[ord]
      e <- end[ord]
    }

    exon_len <- e - s + 1L
    tx_start <- cumsum(c(1L, head(exon_len, -1L)))
    tx_end <- tx_start + exon_len - 1L

    .(
      seqnames = seqnames[1],
      strand = strand_val,
      gene_name = gene_name_val,
      exon_g_start = s,
      exon_g_end = e,
      exon_tx_start = tx_start,
      exon_tx_end = tx_end
    )
  }, by = transcript_id]

  setkey(tx_models, transcript_id)
  tx_models
}

map_genomic_interval_to_tx <- function(tx_exons, g_start, g_end) {
  ov <- tx_exons[!(exon_g_end < g_start | exon_g_start > g_end)]
  if (!nrow(ov)) return(c(NA_integer_, NA_integer_))

  parts <- ov[, {
    os <- pmax(exon_g_start, g_start)
    oe <- pmin(exon_g_end, g_end)

    if (strand[1] == "+") {
      ts <- exon_tx_start + (os - exon_g_start)
      te <- exon_tx_start + (oe - exon_g_start)
    } else {
      ts <- exon_tx_start + (exon_g_end - oe)
      te <- exon_tx_start + (exon_g_end - os)
    }

    .(ts = ts, te = te)
  }]

  coords <- c(parts$ts, parts$te)
  c(min(coords), max(coords))
}

assign_duplexes_to_transcripts <- function(dg, tx_models, min_overlap = 5L) {
  dg <- copy(dg)
  dg[, dg_row_id := .I]

  tx_models_by_tx <- split(tx_models, by = "transcript_id", keep.by = TRUE)

  tx_index <- unique(tx_models[, .(gene_name, seqnames, strand, transcript_id)])
  tx_index[, group_key := paste(gene_name, seqnames, strand, sep = "___")]
  tx_candidates_by_group <- split(tx_index$transcript_id, f = tx_index$group_key)

  results_list <- vector("list", nrow(dg))
  n_total <- nrow(dg)

  for (k in seq_len(n_total)) {
    if (k %% 1000L == 0L) {
      message(sprintf("Processed %d / %d duplexes", k, n_total))
    }

    row <- dg[k]

    group_key <- paste(row$gene1, row$chrom1, row$strand1, sep = "___")
    candidate_txs <- tx_candidates_by_group[[group_key]]

    if (is.null(candidate_txs) || length(candidate_txs) == 0L) {
      next
    }

    mapped_candidates <- vector("list", length(candidate_txs))
    keep_n <- 0L

    for (tx in candidate_txs) {
      tx_ex <- tx_models_by_tx[[tx]]

      if (is.null(tx_ex) || nrow(tx_ex) == 0L) next

      ov_left <- tx_ex[
        !(exon_g_end < row$start1 | exon_g_start > row$end1)
      ]
      ov_right <- tx_ex[
        !(exon_g_end < row$start2 | exon_g_start > row$end2)
      ]

      if (nrow(ov_left) == 0 || nrow(ov_right) == 0) {
        next
      }

      left_bp <- sum(
        pmax(
          0L,
          pmin(ov_left$exon_g_end, row$end1) -
            pmax(ov_left$exon_g_start, row$start1) + 1L
        )
      )

      right_bp <- sum(
        pmax(
          0L,
          pmin(ov_right$exon_g_end, row$end2) -
            pmax(ov_right$exon_g_start, row$start2) + 1L
        )
      )

      if (left_bp < min_overlap || right_bp < min_overlap) {
        next
      }

      left_span <- map_genomic_interval_to_tx(tx_ex, row$start1, row$end1)
      right_span <- map_genomic_interval_to_tx(tx_ex, row$start2, row$end2)

      if (any(is.na(c(left_span, right_span)))) {
        next
      }

      keep_n <- keep_n + 1L
      mapped_candidates[[keep_n]] <- data.table(
        dg_row_id = row$dg_row_id,
        transcript_id = tx,
        gene1 = row$gene1,
        gene2 = row$gene2,
        chrom1 = row$chrom1,
        start1 = row$start1,
        end1 = row$end1,
        chrom2 = row$chrom2,
        start2 = row$start2,
        end2 = row$end2,
        strand1 = row$strand1,
        strand2 = row$strand2,
        support = row$support,
        dg_id = row$dg_id,
        dg_score = row$dg_score,
        arm1_len = row$arm1_len,
        arm2_len = row$arm2_len,
        genomic_span = row$genomic_span,
        tx_left_start = left_span[1],
        tx_left_end = left_span[2],
        tx_right_start = right_span[1],
        tx_right_end = right_span[2],
        left_overlap_bp = left_bp,
        right_overlap_bp = right_bp
      )
    }

    if (keep_n == 0L) {
      next
    }

    cand_dt <- rbindlist(
      mapped_candidates[seq_len(keep_n)],
      use.names = TRUE,
      fill = TRUE
    )

    cand_dt[, mapped_nt := (tx_left_end - tx_left_start + 1L) +
              (tx_right_end - tx_right_start + 1L)]

    setorder(cand_dt, -mapped_nt, -left_overlap_bp, -right_overlap_bp, -support)

    results_list[[k]] <- cand_dt[1]
  }

  out <- rbindlist(results_list, use.names = TRUE, fill = TRUE)
  out
}

stop_if_missing(bedpe_file, "CRSSANT BEDPE")
stop_if_missing(gtf_file, "GTF")
stop_if_missing(vienna_master_pairs_file, "Complete master ViennaRNA pair table")

message2("Output directory: %s", outdir)

message2("Loading complete master ViennaRNA transcript-level base-pair table...")

vienna_pairs_df <- fread(vienna_master_pairs_file)

message2("Columns in master ViennaRNA pair table:")
print(names(vienna_pairs_df))

if (!"ensembl_transcript_id" %in% names(vienna_pairs_df)) {
  if ("ensembl_transcript_id_version" %in% names(vienna_pairs_df)) {
    vienna_pairs_df[, ensembl_transcript_id := sub("\\.\\d+$", "", ensembl_transcript_id_version)]
  } else if ("transcript_id" %in% names(vienna_pairs_df)) {
    vienna_pairs_df[, ensembl_transcript_id := sub("\\.\\d+$", "", transcript_id)]
  } else {
    stop("Cannot find transcript ID column in master ViennaRNA pair table.", call. = FALSE)
  }
}

if (!"ensembl_transcript_id_version" %in% names(vienna_pairs_df)) {
  vienna_pairs_df[, ensembl_transcript_id_version := ensembl_transcript_id]
}

message2(
  "Loaded %s ViennaRNA base-pair rows from complete master table.",
  format(nrow(vienna_pairs_df), big.mark = ",")
)

message2(
  "Master ViennaRNA pair table contains %s unique transcripts.",
  format(length(unique(vienna_pairs_df$ensembl_transcript_id)), big.mark = ",")
)

fwrite(
  vienna_pairs_df,
  file.path(outdir, "ViennaRNA_transcript_level_pairs.tsv.gz"),
  sep = "\t"
)

fwrite(
  vienna_pairs_df,
  file.path(outdir, "ViennaRNA_transcript_level_pairs_from_MASTER.tsv.gz"),
  sep = "\t"
)

message2("Reading CRSSANT BEDPE...")
bedpe_dt <- read_crssant_bedpe(bedpe_file)

message2("Total BEDPE rows: %s", format(nrow(bedpe_dt), big.mark = ","))
message2("Same chromosome rows: %s", format(sum(bedpe_dt$same_chr), big.mark = ","))
message2("Same gene rows: %s", format(sum(bedpe_dt$same_gene), big.mark = ","))
message2("Same strand rows: %s", format(sum(bedpe_dt$same_strand), big.mark = ","))

fwrite(
  bedpe_dt,
  file.path(outdir, "CRSSANT_BEDPE_clean_all.tsv.gz"),
  sep = "\t"
)

bedpe_intra <- bedpe_dt[
  same_chr == TRUE &
    same_gene == TRUE &
    same_strand == TRUE &
    support >= min_support
]

bedpe_intra_strict <- bedpe_dt[
  same_chr == TRUE &
    same_gene == TRUE &
    same_strand == TRUE &
    support >= min_support_strict
]

message2("Rows after intra-gene / same-strand / support >= %d filtering: %s",
         min_support, format(nrow(bedpe_intra), big.mark = ","))

message2("Rows after stricter support >= %d filtering: %s",
         min_support_strict, format(nrow(bedpe_intra_strict), big.mark = ","))

fwrite(
  bedpe_intra,
  file.path(outdir, "CRSSANT_BEDPE_intragene_support2.tsv.gz"),
  sep = "\t"
)

fwrite(
  bedpe_intra_strict,
  file.path(outdir, "CRSSANT_BEDPE_intragene_support5.tsv.gz"),
  sep = "\t"
)

message2("Building transcript exon models from GENCODE v47...")
tx_models <- build_exon_models(gtf_file)

message2("Transcript exon models built for %s transcripts.",
         format(length(unique(tx_models$transcript_id)), big.mark = ","))

fwrite(
  tx_models,
  file.path(outdir, "Gencode_v47_transcript_exon_models.tsv.gz"),
  sep = "\t"
)

message2("Assigning filtered CRSSANT duplexes to transcripts...")

mapped_duplexes <- assign_duplexes_to_transcripts(
  dg = bedpe_intra,
  tx_models = tx_models,
  min_overlap = min_overlap_bp
)

message2("Mapped intra-transcript duplexes retained: %s",
         format(nrow(mapped_duplexes), big.mark = ","))

if (nrow(mapped_duplexes) == 0) {
  stop("No duplexes were assigned to the same transcript. Check GTF compatibility and coordinates.",
       call. = FALSE)
}

fwrite(
  mapped_duplexes,
  file.path(outdir, "mapped_intra_transcript_duplexes.tsv.gz"),
  sep = "\t"
)

tx_summary <- mapped_duplexes[, .(
  n_duplexes = .N,
  total_support = sum(support, na.rm = TRUE),
  median_support = median(support, na.rm = TRUE),
  mean_support = mean(support, na.rm = TRUE)
), by = .(transcript_id)][order(-n_duplexes, -total_support)]

fwrite(
  tx_summary,
  file.path(outdir, "mapped_duplexes_transcript_summary.tsv"),
  sep = "\t"
)

gene_summary <- mapped_duplexes[, .(
  n_duplexes = .N,
  total_support = sum(support, na.rm = TRUE),
  median_support = median(support, na.rm = TRUE),
  mean_support = mean(support, na.rm = TRUE)
), by = .(gene1)][order(-n_duplexes, -total_support)]

fwrite(
  gene_summary,
  file.path(outdir, "mapped_duplexes_gene_summary.tsv"),
  sep = "\t"
)

mapped_compact <- mapped_duplexes[, .(
  gene1,
  transcript_id,
  chrom1, start1, end1,
  chrom2, start2, end2,
  tx_left_start, tx_left_end,
  tx_right_start, tx_right_end,
  support,
  dg_id,
  dg_score,
  strand1,
  strand2
)]

fwrite(
  mapped_compact,
  file.path(outdir, "mapped_intra_transcript_duplexes_compact.tsv"),
  sep = "\t"
)

p_support <- ggplot(as.data.frame(bedpe_intra), aes(x = support)) +
  geom_histogram(bins = 40) +
  theme_classic(base_size = 12) +
  labs(
    title = "CRSSANT intra-gene duplex support distribution",
    x = "Support",
    y = "Number of duplexes"
  )

ggsave(
  file.path(outdir, "CRSSANT_intragene_support_histogram.pdf"),
  p_support,
  width = 7,
  height = 5
)

p_tx_top <- ggplot(as.data.frame(head(tx_summary, 25)), aes(x = reorder(transcript_id, n_duplexes), y = n_duplexes)) +
  geom_col() +
  coord_flip() +
  theme_classic(base_size = 11) +
  labs(
    title = "Top transcripts by mapped duplex count",
    x = "Transcript",
    y = "Number of mapped duplexes"
  )

ggsave(
  file.path(outdir, "Top_transcripts_by_mapped_duplex_count.pdf"),
  p_tx_top,
  width = 8,
  height = 7
)

p_gene_top <- ggplot(as.data.frame(head(gene_summary, 25)), aes(x = reorder(gene1, n_duplexes), y = n_duplexes)) +
  geom_col() +
  coord_flip() +
  theme_classic(base_size = 11) +
  labs(
    title = "Top genes by mapped duplex count",
    x = "Gene",
    y = "Number of mapped duplexes"
  )

ggsave(
  file.path(outdir, "Top_genes_by_mapped_duplex_count.pdf"),
  p_gene_top,
  width = 8,
  height = 7
)

summary_lines <- c(
  "PARIS2 enriched mRNA / CRSSANT vs ViennaRNA pre-comparison pipeline",
  "",
  paste0("BEDPE input: ", bedpe_file),
  paste0("GTF: ", gtf_file),
  paste0("Output directory: ", outdir),
  "",
  paste0("ViennaRNA master pair table: ", vienna_master_pairs_file),
  paste0("ViennaRNA base-pair rows loaded from master table: ", format(nrow(vienna_pairs_df), big.mark = ",")),
  paste0("Unique ViennaRNA transcripts in master pair table: ", format(length(unique(vienna_pairs_df$ensembl_transcript_id)), big.mark = ",")),
  "",
  paste0("Raw BEDPE rows: ", format(nrow(bedpe_dt), big.mark = ",")),
  paste0("Same chromosome rows: ", format(sum(bedpe_dt$same_chr), big.mark = ",")),
  paste0("Same gene rows: ", format(sum(bedpe_dt$same_gene), big.mark = ",")),
  paste0("Same strand rows: ", format(sum(bedpe_dt$same_strand), big.mark = ",")),
  "",
  paste0("Filtered intragene rows (support >= ", min_support, "): ", format(nrow(bedpe_intra), big.mark = ",")),
  paste0("Filtered intragene rows (support >= ", min_support_strict, "): ", format(nrow(bedpe_intra_strict), big.mark = ",")),
  paste0("Mapped intra-transcript duplexes: ", format(nrow(mapped_duplexes), big.mark = ",")),
  paste0("Unique mapped transcripts: ", format(length(unique(mapped_duplexes$transcript_id)), big.mark = ",")),
  paste0("Unique mapped genes: ", format(length(unique(mapped_duplexes$gene1)), big.mark = ","))
)

writeLines(
  summary_lines,
  con = file.path(outdir, "run_summary.txt")
)

message2("Done. Outputs written to: %s", outdir)
