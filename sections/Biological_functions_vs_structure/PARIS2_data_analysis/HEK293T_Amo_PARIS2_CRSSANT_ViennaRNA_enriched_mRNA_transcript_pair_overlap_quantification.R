suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

mapped_duplex_file <- "path_to/mapped_intra_transcript_duplexes.tsv.gz"
vienna_pairs_file <- "path_to/ViennaRNA_transcript_level_pairs.tsv.gz"
vienna_struct_file <- "path_to/ViennaRNA_structures_recovered.tsv"
outdir <- "path_to/PARIS2_enriched_mRNA_ViennaRNA_pair_overlap_analysis_output_HEK293T"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

RUN_PERMUTATIONS <- TRUE
n_perm <- 1000L
set.seed(123)

message2 <- function(...) cat(sprintf(...), "\n")
stop_if_missing <- function(x, label) if (!file.exists(x)) stop(label, " not found: ", x, call. = FALSE)
normalize_tx_id <- function(x) sub("\\..*$", "", x)
safe_fraction <- function(num, den) ifelse(den > 0, num / den, NA_real_)

count_pairs_linking_regions <- function(pair_dt, left_start, left_end, right_start, right_end, tol = 0L) {
  if (nrow(pair_dt) == 0) return(0L)
  ls <- left_start - tol; le <- left_end + tol
  rs <- right_start - tol; re <- right_end + tol
  hits <- pair_dt[((i >= ls & i <= le) & (j >= rs & j <= re)) |
                    ((j >= ls & j <= le) & (i >= rs & i <= re))]
  nrow(hits)
}

has_pairs_linking_regions <- function(pair_dt, left_start, left_end, right_start, right_end, tol = 0L) {
  as.integer(count_pairs_linking_regions(pair_dt, left_start, left_end, right_start, right_end, tol) > 0L)
}

empirical_p_ge <- function(obs, null_vec) {
  if (length(null_vec) == 0L) return(NA_real_)
  (sum(null_vec >= obs, na.rm = TRUE) + 1) / (length(null_vec) + 1)
}

stop_if_missing(mapped_duplex_file, "Mapped duplex table")
stop_if_missing(vienna_pairs_file, "ViennaRNA pair table")
stop_if_missing(vienna_struct_file, "ViennaRNA structure table")

message2("Reading mapped duplexes from: %s", mapped_duplex_file)
mapped_duplexes <- fread(file = mapped_duplex_file)
vienna_pairs    <- fread(file = vienna_pairs_file)
vienna_struct   <- fread(file = vienna_struct_file)

message2("ViennaRNA structure rows loaded: %s", format(nrow(vienna_struct), big.mark = ","))
message2("Unique ViennaRNA structure transcripts: %s", format(uniqueN(vienna_struct$ensembl_transcript_id), big.mark = ","))

if (!"transcript_id" %in% names(mapped_duplexes)) stop("mapped_duplexes does not contain transcript_id", call. = FALSE)
mapped_duplexes[, transcript_id_norm := normalize_tx_id(transcript_id)]

if ("ensembl_transcript_id" %in% names(vienna_pairs)) {
  vienna_pairs[, transcript_id_norm := normalize_tx_id(ensembl_transcript_id)]
} else if ("ensembl_transcript_id_version" %in% names(vienna_pairs)) {
  vienna_pairs[, transcript_id_norm := normalize_tx_id(ensembl_transcript_id_version)]
} else stop("ViennaRNA pair table lacks transcript ID column.", call. = FALSE)

if ("ensembl_transcript_id" %in% names(vienna_struct)) {
  vienna_struct[, transcript_id_norm := normalize_tx_id(ensembl_transcript_id)]
} else if ("ensembl_transcript_id_version" %in% names(vienna_struct)) {
  vienna_struct[, transcript_id_norm := normalize_tx_id(ensembl_transcript_id_version)]
} else stop("ViennaRNA structure table lacks transcript ID column.", call. = FALSE)

vienna_meta <- unique(vienna_struct[, .(transcript_id_norm, ensembl_transcript_id_version, Transcript_length, MFE)], by = "transcript_id_norm")

message2("Mapped duplex rows loaded: %s", format(nrow(mapped_duplexes), big.mark = ","))
message2("ViennaRNA pair rows loaded: %s", format(nrow(vienna_pairs), big.mark = ","))

common_tx <- intersect(unique(mapped_duplexes$transcript_id_norm), unique(vienna_pairs$transcript_id_norm))
mapped_duplexes <- mapped_duplexes[transcript_id_norm %in% common_tx]
vienna_pairs <- vienna_pairs[transcript_id_norm %in% common_tx]
message2("Transcripts present in both mapped duplexes and ViennaRNA pairs: %s", format(length(common_tx), big.mark = ","))
message2("Mapped duplex rows retained: %s", format(nrow(mapped_duplexes), big.mark = ","))
if (nrow(mapped_duplexes) == 0) stop("No mapped duplexes remain after transcript intersection.", call. = FALSE)

if (!all(c("pos1", "pos2") %in% names(vienna_pairs))) stop("ViennaRNA pair table lacks pos1/pos2.", call. = FALSE)
setnames(vienna_pairs, c("pos1", "pos2"), c("i", "j"))
vienna_pairs_by_tx <- split(
  vienna_pairs[, .(transcript_id_norm, i, j)],
  by = "transcript_id_norm",
  keep.by = FALSE
)

message2("Computing duplex-level ViennaRNA support...")
duplex_results_list <- vector("list", nrow(mapped_duplexes))
for (k in seq_len(nrow(mapped_duplexes))) {
  if (k %% 1000L == 0L) message2("Processed %s / %s duplexes", format(k, big.mark = ","), format(nrow(mapped_duplexes), big.mark = ","))
  row <- mapped_duplexes[k]
  pred_pairs <- vienna_pairs_by_tx[[row$transcript_id_norm]]
  if (is.null(pred_pairs)) pred_pairs <- data.table(i = integer(), j = integer())
  n_exact <- count_pairs_linking_regions(pred_pairs, row$tx_left_start, row$tx_left_end, row$tx_right_start, row$tx_right_end, 0L)
  n_t5 <- count_pairs_linking_regions(pred_pairs, row$tx_left_start, row$tx_left_end, row$tx_right_start, row$tx_right_end, 5L)
  n_t10 <- count_pairs_linking_regions(pred_pairs, row$tx_left_start, row$tx_left_end, row$tx_right_start, row$tx_right_end, 10L)
  n_t15 <- count_pairs_linking_regions(pred_pairs, row$tx_left_start, row$tx_left_end, row$tx_right_start, row$tx_right_end, 15L)
  duplex_results_list[[k]] <- data.table(
    gene1 = row$gene1, gene2 = row$gene2, transcript_id = row$transcript_id,
    transcript_id_norm = row$transcript_id_norm,
    chrom1 = row$chrom1, start1 = row$start1, end1 = row$end1,
    chrom2 = row$chrom2, start2 = row$start2, end2 = row$end2,
    tx_left_start = row$tx_left_start, tx_left_end = row$tx_left_end,
    tx_right_start = row$tx_right_start, tx_right_end = row$tx_right_end,
    support = row$support, dg_id = row$dg_id, dg_score = row$dg_score,
    mapped_nt = row$mapped_nt, left_overlap_bp = row$left_overlap_bp, right_overlap_bp = row$right_overlap_bp,
    n_pred_pairs_transcript = nrow(pred_pairs),
    n_pairs_exact = n_exact, n_pairs_t5 = n_t5, n_pairs_t10 = n_t10, n_pairs_t15 = n_t15,
    supported_exact = as.integer(n_exact > 0), supported_t5 = as.integer(n_t5 > 0),
    supported_t10 = as.integer(n_t10 > 0), supported_t15 = as.integer(n_t15 > 0)
  )
}

duplex_support <- rbindlist(duplex_results_list, use.names = TRUE, fill = TRUE)
duplex_support <- merge(duplex_support, vienna_meta, by = "transcript_id_norm", all.x = TRUE)
fwrite(duplex_support, file.path(outdir, "duplex_level_PARIS2_Amo_ViennaRNA_support.tsv.gz"), sep = "\t")

transcript_summary <- duplex_support[, .(
  gene_label = paste(unique(gene1), collapse = ";"),
  n_duplexes = as.integer(.N),
  total_support_reads = as.integer(sum(support, na.rm = TRUE)),
  median_duplex_support = as.numeric(median(support, na.rm = TRUE)),
  mean_duplex_support = as.numeric(mean(support, na.rm = TRUE)),
  Transcript_length = as.integer(unique(Transcript_length)[1]),
  MFE = as.numeric(unique(MFE)[1]),
  n_pred_pairs_transcript = as.integer(unique(n_pred_pairs_transcript)[1]),
  n_supported_exact = as.integer(sum(supported_exact, na.rm = TRUE)),
  n_supported_t5 = as.integer(sum(supported_t5, na.rm = TRUE)),
  n_supported_t10 = as.integer(sum(supported_t10, na.rm = TRUE)),
  n_supported_t15 = as.integer(sum(supported_t15, na.rm = TRUE)),
  frac_supported_exact = as.numeric(safe_fraction(sum(supported_exact, na.rm = TRUE), .N)),
  frac_supported_t5 = as.numeric(safe_fraction(sum(supported_t5, na.rm = TRUE), .N)),
  frac_supported_t10 = as.numeric(safe_fraction(sum(supported_t10, na.rm = TRUE), .N)),
  frac_supported_t15 = as.numeric(safe_fraction(sum(supported_t15, na.rm = TRUE), .N))
), by = .(transcript_id_norm, transcript_id, ensembl_transcript_id_version)]

transcript_summary[, Normalized_MFE := as.numeric(MFE) / as.numeric(Transcript_length)]

fwrite(
  transcript_summary,
  file.path(outdir, "transcript_level_PARIS2_Amo_ViennaRNA_support_summary.tsv"),
  sep = "\t"
)

p1 <- ggplot(as.data.frame(duplex_support), aes(x = factor(supported_exact))) + geom_bar() + theme_classic(base_size = 12) + labs(title = "Exact ViennaRNA support of mapped PARIS2 Amo duplexes", x = "Exact support", y = "Number of duplexes")
ggsave(file.path(outdir, "duplex_exact_support_barplot.pdf"), p1, width = 6, height = 5)
p2 <- ggplot(as.data.frame(duplex_support), aes(x = factor(supported_t15))) + geom_bar() + theme_classic(base_size = 12) + labs(title = "ViennaRNA support of mapped PARIS2 Amo duplexes at +/- 15 nt", x = "Support at +/- 15 nt", y = "Number of duplexes")
ggsave(file.path(outdir, "duplex_t15_support_barplot.pdf"), p2, width = 6, height = 5)
p3 <- ggplot(as.data.frame(duplex_support), aes(x = support, fill = factor(supported_t15))) + geom_histogram(bins = 30, position = "identity", alpha = 0.6) + theme_classic(base_size = 12) + labs(title = "Experimental duplex support vs ViennaRNA support at +/- 15 nt", x = "CRSSANT support", y = "Number of duplexes", fill = "Supported")
ggsave(file.path(outdir, "duplex_support_vs_t15_support_histogram.pdf"), p3, width = 7, height = 5)
if (all(!is.na(transcript_summary$Normalized_MFE))) {
  p4 <- ggplot(as.data.frame(transcript_summary), aes(x = Normalized_MFE, y = frac_supported_t15)) + geom_point(alpha = 0.6) + geom_smooth(method = "lm", se = TRUE) + theme_classic(base_size = 12) + labs(title = "Normalized transcript MFE vs fraction of PARIS2 Amo duplexes supported by ViennaRNA", x = "ViennaRNA normalized MFE (kcal/mol/nt)", y = "Fraction supported at +/- 15 nt")
  ggsave(file.path(outdir, "normalized_MFE_vs_fraction_supported_t15.pdf"), p4, width = 7, height = 5)
}

perm_outdir <- file.path(outdir, "matched_within_transcript_randomization")

if (RUN_PERMUTATIONS) {

  dir.create(perm_outdir, recursive = TRUE, showWarnings = FALSE)

  message2("Starting matched within-transcript randomization: %d permutations per duplex", n_perm)

  duplex_perm_input <- copy(duplex_support)
  duplex_perm_input[, duplex_row := .I]

  duplex_perm_input[, armA_start := fifelse(
    tx_left_start <= tx_right_start,
    tx_left_start,
    tx_right_start
  )]

  duplex_perm_input[, armA_end := fifelse(
    tx_left_start <= tx_right_start,
    tx_left_end,
    tx_right_end
  )]

  duplex_perm_input[, armB_start := fifelse(
    tx_left_start <= tx_right_start,
    tx_right_start,
    tx_left_start
  )]

  duplex_perm_input[, armB_end := fifelse(
    tx_left_start <= tx_right_start,
    tx_right_end,
    tx_left_end
  )]

  duplex_perm_input[, armA_len := armA_end - armA_start + 1L]
  duplex_perm_input[, armB_len := armB_end - armB_start + 1L]
  duplex_perm_input[, gap_len := armB_start - armA_end - 1L]

  duplex_perm_input <- duplex_perm_input[
    !is.na(Transcript_length) &
      !is.na(armA_len) &
      !is.na(armB_len) &
      !is.na(gap_len) &
      Transcript_length > 0L &
      armA_len > 0L &
      armB_len > 0L &
      gap_len >= 0L &
      (armA_len + gap_len + armB_len) <= Transcript_length
  ]

  message2(
    "Duplexes entering matched randomization: %s",
    format(nrow(duplex_perm_input), big.mark = ",")
  )

  if (nrow(duplex_perm_input) == 0) {
    stop("No duplexes available for matched randomization.", call. = FALSE)
  }

  obs_exact <- duplex_perm_input$supported_exact
  obs_t5    <- duplex_perm_input$supported_t5
  obs_t10   <- duplex_perm_input$supported_t10
  obs_t15   <- duplex_perm_input$supported_t15

  generate_matched_random_duplex_ordered <- function(tx_length, armA_len, armB_len, gap_len) {

    total_span <- armA_len + gap_len + armB_len

    if (
      is.na(tx_length) ||
      is.na(armA_len) ||
      is.na(armB_len) ||
      is.na(gap_len)
    ) {
      return(NULL)
    }

    if (
      tx_length <= 0L ||
      armA_len <= 0L ||
      armB_len <= 0L ||
      gap_len < 0L ||
      total_span > tx_length
    ) {
      return(NULL)
    }

    max_armA_start <- tx_length - total_span + 1L

    if (max_armA_start < 1L) {
      return(NULL)
    }

    rand_armA_start <- sample.int(max_armA_start, size = 1L)
    rand_armA_end   <- rand_armA_start + armA_len - 1L

    rand_armB_start <- rand_armA_end + gap_len + 1L
    rand_armB_end   <- rand_armB_start + armB_len - 1L

    list(
      armA_start = rand_armA_start,
      armA_end   = rand_armA_end,
      armB_start = rand_armB_start,
      armB_end   = rand_armB_end
    )
  }

  perm_exact_mat <- matrix(NA_integer_, nrow = nrow(duplex_perm_input), ncol = n_perm)
  perm_t5_mat    <- matrix(NA_integer_, nrow = nrow(duplex_perm_input), ncol = n_perm)
  perm_t10_mat   <- matrix(NA_integer_, nrow = nrow(duplex_perm_input), ncol = n_perm)
  perm_t15_mat   <- matrix(NA_integer_, nrow = nrow(duplex_perm_input), ncol = n_perm)

  message2("Generating matched random duplexes and scoring ViennaRNA support...")

  for (k in seq_len(nrow(duplex_perm_input))) {

    if (k %% 25L == 0L || k == 1L || k == nrow(duplex_perm_input)) {
      message2(
        "Permutation scoring duplex %s / %s",
        format(k, big.mark = ","),
        format(nrow(duplex_perm_input), big.mark = ",")
      )
    }

    row <- duplex_perm_input[k]
    tx  <- row$transcript_id_norm

    pred_pairs <- vienna_pairs_by_tx[[tx]]

    if (is.null(pred_pairs)) {
      pred_pairs <- data.table(i = integer(), j = integer())
    }

    for (p in seq_len(n_perm)) {

      rand <- generate_matched_random_duplex_ordered(
        tx_length = as.integer(row$Transcript_length),
        armA_len  = as.integer(row$armA_len),
        armB_len  = as.integer(row$armB_len),
        gap_len   = as.integer(row$gap_len)
      )

      if (is.null(rand)) {
        next
      }

      perm_exact_mat[k, p] <- has_pairs_linking_regions(
        pred_pairs,
        rand$armA_start,
        rand$armA_end,
        rand$armB_start,
        rand$armB_end,
        tol = 0L
      )

      perm_t5_mat[k, p] <- has_pairs_linking_regions(
        pred_pairs,
        rand$armA_start,
        rand$armA_end,
        rand$armB_start,
        rand$armB_end,
        tol = 5L
      )

      perm_t10_mat[k, p] <- has_pairs_linking_regions(
        pred_pairs,
        rand$armA_start,
        rand$armA_end,
        rand$armB_start,
        rand$armB_end,
        tol = 10L
      )

      perm_t15_mat[k, p] <- has_pairs_linking_regions(
        pred_pairs,
        rand$armA_start,
        rand$armA_end,
        rand$armB_start,
        rand$armB_end,
        tol = 15L
      )
    }
  }

  duplex_perm_summary <- duplex_perm_input[, .(
    duplex_row,
    gene1,
    gene2,
    transcript_id,
    transcript_id_norm,
    support,
    dg_id,
    dg_score,
    Transcript_length,
    MFE,
    armA_start,
    armA_end,
    armB_start,
    armB_end,
    armA_len,
    armB_len,
    gap_len,
    obs_supported_exact = supported_exact,
    obs_supported_t5    = supported_t5,
    obs_supported_t10   = supported_t10,
    obs_supported_t15   = supported_t15
  )]

  duplex_perm_summary[, exp_supported_exact := rowMeans(perm_exact_mat, na.rm = TRUE)]
  duplex_perm_summary[, exp_supported_t5    := rowMeans(perm_t5_mat, na.rm = TRUE)]
  duplex_perm_summary[, exp_supported_t10   := rowMeans(perm_t10_mat, na.rm = TRUE)]
  duplex_perm_summary[, exp_supported_t15   := rowMeans(perm_t15_mat, na.rm = TRUE)]

  duplex_perm_summary[, enrichment_exact := ifelse(exp_supported_exact > 0, obs_supported_exact / exp_supported_exact, NA_real_)]
  duplex_perm_summary[, enrichment_t5    := ifelse(exp_supported_t5 > 0, obs_supported_t5 / exp_supported_t5, NA_real_)]
  duplex_perm_summary[, enrichment_t10   := ifelse(exp_supported_t10 > 0, obs_supported_t10 / exp_supported_t10, NA_real_)]
  duplex_perm_summary[, enrichment_t15   := ifelse(exp_supported_t15 > 0, obs_supported_t15 / exp_supported_t15, NA_real_)]

  duplex_perm_summary[, empirical_p_exact := vapply(
    seq_len(.N),
    function(i) empirical_p_ge(obs_supported_exact[i], perm_exact_mat[i, ]),
    numeric(1)
  )]

  duplex_perm_summary[, empirical_p_t5 := vapply(
    seq_len(.N),
    function(i) empirical_p_ge(obs_supported_t5[i], perm_t5_mat[i, ]),
    numeric(1)
  )]

  duplex_perm_summary[, empirical_p_t10 := vapply(
    seq_len(.N),
    function(i) empirical_p_ge(obs_supported_t10[i], perm_t10_mat[i, ]),
    numeric(1)
  )]

  duplex_perm_summary[, empirical_p_t15 := vapply(
    seq_len(.N),
    function(i) empirical_p_ge(obs_supported_t15[i], perm_t15_mat[i, ]),
    numeric(1)
  )]

  fwrite(
    duplex_perm_summary,
    file.path(perm_outdir, "duplex_level_observed_vs_random_support.tsv.gz"),
    sep = "\t"
  )

  global_obs_exact <- mean(obs_exact, na.rm = TRUE)
  global_obs_t5    <- mean(obs_t5, na.rm = TRUE)
  global_obs_t10   <- mean(obs_t10, na.rm = TRUE)
  global_obs_t15   <- mean(obs_t15, na.rm = TRUE)

  global_perm_exact <- colMeans(perm_exact_mat, na.rm = TRUE)
  global_perm_t5    <- colMeans(perm_t5_mat, na.rm = TRUE)
  global_perm_t10   <- colMeans(perm_t10_mat, na.rm = TRUE)
  global_perm_t15   <- colMeans(perm_t15_mat, na.rm = TRUE)

  null_distribution <- data.table(
    permutation = seq_len(n_perm),
    fraction_supported_exact = global_perm_exact,
    fraction_supported_t5 = global_perm_t5,
    fraction_supported_t10 = global_perm_t10,
    fraction_supported_t15 = global_perm_t15
  )

  fwrite(
    null_distribution,
    file.path(perm_outdir, "global_random_support_null_distribution.tsv.gz"),
    sep = "\t"
  )

  global_randomization_summary <- data.table(
    tolerance = c("exact", "plus_minus_5_nt", "plus_minus_10_nt", "plus_minus_15_nt"),
    n_duplexes_tested = nrow(duplex_perm_input),
    n_permutations = n_perm,
    observed_fraction_supported = c(
      global_obs_exact,
      global_obs_t5,
      global_obs_t10,
      global_obs_t15
    ),
    expected_fraction_supported = c(
      mean(global_perm_exact, na.rm = TRUE),
      mean(global_perm_t5, na.rm = TRUE),
      mean(global_perm_t10, na.rm = TRUE),
      mean(global_perm_t15, na.rm = TRUE)
    ),
    sd_random_fraction_supported = c(
      sd(global_perm_exact, na.rm = TRUE),
      sd(global_perm_t5, na.rm = TRUE),
      sd(global_perm_t10, na.rm = TRUE),
      sd(global_perm_t15, na.rm = TRUE)
    ),
    empirical_p_observed_ge_random = c(
      empirical_p_ge(global_obs_exact, global_perm_exact),
      empirical_p_ge(global_obs_t5, global_perm_t5),
      empirical_p_ge(global_obs_t10, global_perm_t10),
      empirical_p_ge(global_obs_t15, global_perm_t15)
    )
  )

  global_randomization_summary[, enrichment_observed_over_expected :=
                                 ifelse(
                                   expected_fraction_supported > 0,
                                   observed_fraction_supported / expected_fraction_supported,
                                   NA_real_
                                 )
  ]

  global_randomization_summary[, observed_percent_supported :=
                                 100 * observed_fraction_supported]

  global_randomization_summary[, expected_percent_supported :=
                                 100 * expected_fraction_supported]

  fwrite(
    global_randomization_summary,
    file.path(perm_outdir, "global_observed_vs_random_support_summary.tsv"),
    sep = "\t"
  )

  save(
    duplex_perm_input,
    duplex_perm_summary,
    global_randomization_summary,
    null_distribution,
    perm_exact_mat,
    perm_t5_mat,
    perm_t10_mat,
    perm_t15_mat,
    file = file.path(perm_outdir, "matched_randomization_results_all_tolerances.RData")
  )



  null_long <- melt(
    null_distribution,
    id.vars = "permutation",
    variable.name = "tolerance",
    value.name = "random_fraction_supported"
  )

  observed_long <- data.table(
    tolerance = c(
      "fraction_supported_exact",
      "fraction_supported_t5",
      "fraction_supported_t10",
      "fraction_supported_t15"
    ),
    observed_fraction_supported = c(
      global_obs_exact,
      global_obs_t5,
      global_obs_t10,
      global_obs_t15
    )
  )

  p_null <- ggplot(as.data.frame(null_long), aes(x = random_fraction_supported)) +
    geom_histogram(bins = 40) +
    geom_vline(
      data = as.data.frame(observed_long),
      aes(xintercept = observed_fraction_supported),
      linewidth = 0.8
    ) +
    facet_wrap(~ tolerance, scales = "free_x") +
    theme_classic(base_size = 12) +
    labs(
      title = "Matched within-transcript randomization: observed vs random ViennaRNA support",
      subtitle = "Vertical lines show observed PARIS2 support fractions",
      x = "Random fraction supported by ViennaRNA",
      y = "Number of permutations"
    )

  ggsave(
    file.path(perm_outdir, "matched_randomization_observed_vs_random_support.pdf"),
    p_null,
    width = 10,
    height = 7
  )

  message2("Matched randomization complete. Outputs written to: %s", perm_outdir)
  print(global_randomization_summary)
}

if (!RUN_PERMUTATIONS) {
  randomization_rdata_file <- file.path(
    perm_outdir,
    "matched_randomization_results_all_tolerances.RData"
  )
  randomization_summary_file <- file.path(
    perm_outdir,
    "global_observed_vs_random_support_summary.tsv"
  )
  stop_if_missing(randomization_rdata_file, "Saved matched-randomization results")
  stop_if_missing(randomization_summary_file, "Saved matched-randomization summary")
  load(randomization_rdata_file)
  global_randomization_summary <- fread(randomization_summary_file)
}

global_obs_exact <- global_randomization_summary[tolerance == "exact", observed_fraction_supported]
global_obs_t5    <- global_randomization_summary[tolerance == "plus_minus_5_nt", observed_fraction_supported]
global_obs_t10   <- global_randomization_summary[tolerance == "plus_minus_10_nt", observed_fraction_supported]
global_obs_t15   <- global_randomization_summary[tolerance == "plus_minus_15_nt", observed_fraction_supported]

global_perm_exact <- null_distribution$fraction_supported_exact
global_perm_t5    <- null_distribution$fraction_supported_t5
global_perm_t10   <- null_distribution$fraction_supported_t10
global_perm_t15   <- null_distribution$fraction_supported_t15

figure_outdir <- file.path(
  outdir,
  "manuscript_style_support_figure"
)

dir.create(figure_outdir, recursive = TRUE, showWarnings = FALSE)

message2("Generating manuscript-style PARIS2 enriched mRNA support panels...")

panelA_dt <- data.table(
  tolerance = factor(
    c("Exact", "+/- 5 nt", "+/- 10 nt", "+/- 15 nt"),
    levels = c("Exact", "+/- 5 nt", "+/- 10 nt", "+/- 15 nt")
  ),
  observed_real_PARIS2 = c(
    global_obs_exact,
    global_obs_t5,
    global_obs_t10,
    global_obs_t15
  ),
  expected_random_mean = c(
    mean(global_perm_exact, na.rm = TRUE),
    mean(global_perm_t5, na.rm = TRUE),
    mean(global_perm_t10, na.rm = TRUE),
    mean(global_perm_t15, na.rm = TRUE)
  ),
  expected_random_lower = c(
    quantile(global_perm_exact, 0.025, na.rm = TRUE),
    quantile(global_perm_t5, 0.025, na.rm = TRUE),
    quantile(global_perm_t10, 0.025, na.rm = TRUE),
    quantile(global_perm_t15, 0.025, na.rm = TRUE)
  ),
  expected_random_upper = c(
    quantile(global_perm_exact, 0.975, na.rm = TRUE),
    quantile(global_perm_t5, 0.975, na.rm = TRUE),
    quantile(global_perm_t10, 0.975, na.rm = TRUE),
    quantile(global_perm_t15, 0.975, na.rm = TRUE)
  ),
  empirical_p = c(
    empirical_p_ge(global_obs_exact, global_perm_exact),
    empirical_p_ge(global_obs_t5, global_perm_t5),
    empirical_p_ge(global_obs_t10, global_perm_t10),
    empirical_p_ge(global_obs_t15, global_perm_t15)
  )
)

panelA_dt[, enrichment := observed_real_PARIS2 / expected_random_mean]

panelA_dt[, label_text := sprintf(
  "%.2fx\nP = %.3g",
  enrichment,
  empirical_p
)]

panelA_obs <- panelA_dt[, .(
  tolerance,
  group_label = "Observed PARIS2 duplexes",
  value = observed_real_PARIS2
)]

panelA_exp <- panelA_dt[, .(
  tolerance,
  group_label = "Matched random expectation",
  value = expected_random_mean,
  ymin = expected_random_lower,
  ymax = expected_random_upper
)]

panelA_obs[, x := seq_len(.N) - 0.18]
panelA_exp[, x := seq_len(.N) + 0.18]
panelA_dt[, x_mid := seq_len(.N)]

p_panelA <- ggplot() +
  geom_col(
    data = panelA_obs,
    aes(x = x, y = value, fill = group_label),
    width = 0.32,
    color = "black",
    linewidth = 0.25
  ) +
  geom_col(
    data = panelA_exp,
    aes(x = x, y = value, fill = group_label),
    width = 0.32,
    color = "black",
    linewidth = 0.25
  ) +
  geom_errorbar(
    data = panelA_exp,
    aes(x = x, ymin = ymin, ymax = ymax),
    width = 0.10,
    linewidth = 0.55,
    color = "black"
  ) +
  geom_text(
    data = panelA_dt,
    aes(
      x = x_mid,
      y = pmax(observed_real_PARIS2, expected_random_upper) + 0.055,
      label = label_text
    ),
    size = 3.2,
    lineheight = 0.9
  ) +
  scale_fill_manual(
    values = c(
      "Observed PARIS2 duplexes" = "#0072B2",
      "Matched random expectation" = "#BDBDBD"
    )
  ) +
  scale_x_continuous(
    breaks = seq_len(nrow(panelA_dt)),
    labels = levels(panelA_dt$tolerance)
  ) +
  scale_y_continuous(
    labels = function(x) paste0(round(100 * x), "%")
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    legend.title = element_blank(),
    legend.position = "top",
    legend.text = element_text(size = 11)
  ) +
  labs(
    title = "A. PARIS2 duplexes are supported by ViennaRNA above matched-random expectation",
    x = "Support criterion",
    y = "Duplexes supported by ViennaRNA",
    fill = ""
  ) +
  expand_limits(
    y = max(panelA_dt$expected_random_upper, panelA_dt$observed_real_PARIS2, na.rm = TRUE) + 0.14
  )

print(p_panelA)

ggsave(
  filename = file.path(figure_outdir, "Panel_A_global_observed_vs_matched_random_support_PARIS2.tiff"),
  plot = p_panelA,
  device = "tiff",
  width = 4,
  height = 6,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

stopifnot(exists("duplex_perm_summary"))

panelB_dt <- duplex_perm_summary[, .(
  n_duplexes = .N,

  observed_frac_t15 = mean(obs_supported_t15, na.rm = TRUE),
  expected_frac_t15 = mean(exp_supported_t15, na.rm = TRUE),

  observed_frac_exact = mean(obs_supported_exact, na.rm = TRUE),
  expected_frac_exact = mean(exp_supported_exact, na.rm = TRUE)
), by = .(gene1, transcript_id, transcript_id_norm)]

panelB_dt[, delta_t15 := observed_frac_t15 - expected_frac_t15]
panelB_dt[, delta_exact := observed_frac_exact - expected_frac_exact]

panelB_summary <- panelB_dt[, .(
  n_transcripts = .N,

  n_above_random_t15 = sum(delta_t15 > 0, na.rm = TRUE),
  n_below_random_t15 = sum(delta_t15 < 0, na.rm = TRUE),
  n_equal_random_t15 = sum(delta_t15 == 0, na.rm = TRUE),
  percent_above_random_t15 = 100 * mean(delta_t15 > 0, na.rm = TRUE),
  median_delta_t15 = median(delta_t15, na.rm = TRUE),
  mean_delta_t15 = mean(delta_t15, na.rm = TRUE),

  n_above_random_exact = sum(delta_exact > 0, na.rm = TRUE),
  n_below_random_exact = sum(delta_exact < 0, na.rm = TRUE),
  n_equal_random_exact = sum(delta_exact == 0, na.rm = TRUE),
  percent_above_random_exact = 100 * mean(delta_exact > 0, na.rm = TRUE),
  median_delta_exact = median(delta_exact, na.rm = TRUE),
  mean_delta_exact = mean(delta_exact, na.rm = TRUE)
)]

panelB_label <- sprintf(
  "%s/%s transcripts above random expectation\n%.1f%% above random; median delta = %.2f",
  format(panelB_summary$n_above_random_t15, big.mark = ","),
  format(panelB_summary$n_transcripts, big.mark = ","),
  panelB_summary$percent_above_random_t15,
  panelB_summary$median_delta_t15
)

p_panelB <- ggplot(
  as.data.frame(panelB_dt),
  aes(x = delta_t15)
) +
  geom_histogram(
    binwidth = 0.025,
    boundary = 0,
    fill = "grey30",
    color = "black",
    linewidth = 0.20,
    alpha = 0.85
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 1.3,
    color = "black"
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = panelB_label,
    hjust = 1.05,
    vjust = 1.4,
    size = 3.4
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  ) +
  labs(
    title = "B. Most transcripts show greater observed than random PARIS2 support",
    x = "Observed - expected fraction supported at +/- 15 nt",
    y = "Number of transcripts"
  )

print(p_panelB)

ggsave(
  filename = file.path(figure_outdir, "Panel_B_transcript_level_observed_minus_expected_t15_PARIS2.tiff"),
  plot = p_panelB,
  device = "tiff",
  width = 6,
  height = 6,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

panelC_dt <- copy(panelA_dt)

p_panelC <- ggplot(
  as.data.frame(panelC_dt),
  aes(x = tolerance, y = enrichment)
) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = 1, linetype = 2, linewidth = 0.6) +
  geom_text(
    aes(label = sprintf("%.2fx", enrichment)),
    vjust = -0.5,
    size = 4
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  ) +
  labs(
    title = "Observed PARIS2 support is enriched over matched-random expectation",
    x = "Support criterion",
    y = "Observed / expected support"
  ) +
  expand_limits(y = max(panelC_dt$enrichment, na.rm = TRUE) + 0.4)

print(p_panelC)

ggsave(
  filename = file.path(figure_outdir, "Supplementary_global_enrichment_over_matched_random.pdf"),
  plot = p_panelC,
  width = 7,
  height = 5
)

fwrite(
  panelA_dt,
  file.path(figure_outdir, "Panel_A_global_observed_vs_matched_random_support_data.tsv"),
  sep = "\t"
)

fwrite(
  panelB_dt,
  file.path(figure_outdir, "Panel_B_transcript_level_observed_minus_expected_t15_data.tsv"),
  sep = "\t"
)

fwrite(
  panelB_summary,
  file.path(figure_outdir, "Panel_B_transcript_level_observed_minus_expected_t15_summary.tsv"),
  sep = "\t"
)

fwrite(
  panelC_dt,
  file.path(figure_outdir, "Supplementary_global_enrichment_over_matched_random_data.tsv"),
  sep = "\t"
)

message2("Manuscript-style figure panels written to: %s", figure_outdir)

summary_lines <- c(
  "PARIS2 Amo / CRSSANT transcriptome duplex vs ViennaRNA pairset overlap analysis",
  "", paste0("Input mapped duplex table: ", mapped_duplex_file), paste0("Input ViennaRNA pair table: ", vienna_pairs_file), paste0("Input ViennaRNA structure table: ", vienna_struct_file), paste0("Output directory: ", outdir),
  "", paste0("Mapped duplex rows analyzed: ", format(nrow(duplex_support), big.mark = ",")), paste0("Unique transcripts analyzed: ", format(uniqueN(duplex_support$transcript_id_norm), big.mark = ",")), paste0("Unique genes analyzed: ", format(uniqueN(duplex_support$gene1), big.mark = ",")),
  "", paste0("Duplexes with exact support: ", sum(duplex_support$supported_exact, na.rm = TRUE)), paste0("Duplexes with support at +/- 5 nt: ", sum(duplex_support$supported_t5, na.rm = TRUE)), paste0("Duplexes with support at +/- 10 nt: ", sum(duplex_support$supported_t10, na.rm = TRUE)), paste0("Duplexes with support at +/- 15 nt: ", sum(duplex_support$supported_t15, na.rm = TRUE)),
  "", paste0("Fraction supported exactly: ", signif(mean(duplex_support$supported_exact, na.rm = TRUE), 4)), paste0("Fraction supported at +/- 5 nt: ", signif(mean(duplex_support$supported_t5, na.rm = TRUE), 4)), paste0("Fraction supported at +/- 10 nt: ", signif(mean(duplex_support$supported_t10, na.rm = TRUE), 4)), paste0("Fraction supported at +/- 15 nt: ", signif(mean(duplex_support$supported_t15, na.rm = TRUE), 4)),
  "", paste0("Matched randomization enabled: ", RUN_PERMUTATIONS), paste0("n_perm setting: ", n_perm)
)
writeLines(summary_lines, con = file.path(outdir, "run_summary.txt"))
message2("Done. Outputs written to: %s", outdir)

perm_outdir <- file.path(outdir, "matched_within_transcript_randomization")

duplex_file <- file.path(
  outdir,
  "duplex_level_PARIS2_Amo_ViennaRNA_support.tsv.gz"
)

random_summary_file <- file.path(
  perm_outdir,
  "global_observed_vs_random_support_summary.tsv"
)

rdata_file <- file.path(
  perm_outdir,
  "matched_randomization_results_all_tolerances.RData"
)

stopifnot(file.exists(duplex_file))
stopifnot(file.exists(random_summary_file))
stopifnot(file.exists(rdata_file))

duplex_support <- fread(file = duplex_file)
global_randomization_summary <- fread(file = random_summary_file)

load(rdata_file)
stopifnot(exists("duplex_perm_summary"))

direct_PARIS2_ViennaRNA_summary <- data.table(
  metric = c(
    "mapped_duplexes_evaluated",
    "unique_transcripts_evaluated",
    "unique_genes_evaluated",
    "duplexes_supported_exact",
    "duplexes_supported_t5",
    "duplexes_supported_t10",
    "duplexes_supported_t15",
    "fraction_supported_exact",
    "fraction_supported_t5",
    "fraction_supported_t10",
    "fraction_supported_t15",
    "percent_supported_exact",
    "percent_supported_t5",
    "percent_supported_t10",
    "percent_supported_t15"
  ),
  value = c(
    nrow(duplex_support),
    uniqueN(duplex_support$transcript_id_norm),
    uniqueN(duplex_support$gene1),
    sum(duplex_support$supported_exact, na.rm = TRUE),
    sum(duplex_support$supported_t5, na.rm = TRUE),
    sum(duplex_support$supported_t10, na.rm = TRUE),
    sum(duplex_support$supported_t15, na.rm = TRUE),
    mean(duplex_support$supported_exact, na.rm = TRUE),
    mean(duplex_support$supported_t5, na.rm = TRUE),
    mean(duplex_support$supported_t10, na.rm = TRUE),
    mean(duplex_support$supported_t15, na.rm = TRUE),
    100 * mean(duplex_support$supported_exact, na.rm = TRUE),
    100 * mean(duplex_support$supported_t5, na.rm = TRUE),
    100 * mean(duplex_support$supported_t10, na.rm = TRUE),
    100 * mean(duplex_support$supported_t15, na.rm = TRUE)
  )
)

randomization_PARIS2_ViennaRNA_summary <- global_randomization_summary[, .(
  tolerance,
  n_duplexes_tested,
  n_permutations,
  observed_fraction_supported,
  expected_fraction_supported,
  observed_percent_supported,
  expected_percent_supported,
  enrichment_observed_over_expected,
  empirical_p = empirical_p_observed_ge_random
)]

transcript_random_summary <- duplex_perm_summary[, .(
  n_duplexes = .N,
  observed_frac_exact = mean(obs_supported_exact, na.rm = TRUE),
  expected_frac_exact = mean(exp_supported_exact, na.rm = TRUE),
  observed_frac_t5 = mean(obs_supported_t5, na.rm = TRUE),
  expected_frac_t5 = mean(exp_supported_t5, na.rm = TRUE),
  observed_frac_t10 = mean(obs_supported_t10, na.rm = TRUE),
  expected_frac_t10 = mean(exp_supported_t10, na.rm = TRUE),
  observed_frac_t15 = mean(obs_supported_t15, na.rm = TRUE),
  expected_frac_t15 = mean(exp_supported_t15, na.rm = TRUE)
), by = .(gene1, transcript_id, transcript_id_norm)]

transcript_random_summary[, delta_exact := observed_frac_exact - expected_frac_exact]
transcript_random_summary[, delta_t15 := observed_frac_t15 - expected_frac_t15]

transcript_level_PARIS2_ViennaRNA_summary <- transcript_random_summary[, .(
  n_transcripts_evaluated = .N,
  n_above_random_exact = sum(delta_exact > 0, na.rm = TRUE),
  n_below_random_exact = sum(delta_exact < 0, na.rm = TRUE),
  n_equal_random_exact = sum(delta_exact == 0, na.rm = TRUE),
  percent_above_random_exact = 100 * mean(delta_exact > 0, na.rm = TRUE),
  n_above_random_t15 = sum(delta_t15 > 0, na.rm = TRUE),
  n_below_random_t15 = sum(delta_t15 < 0, na.rm = TRUE),
  n_equal_random_t15 = sum(delta_t15 == 0, na.rm = TRUE),
  percent_above_random_t15 = 100 * mean(delta_t15 > 0, na.rm = TRUE),
  median_delta_t15 = median(delta_t15, na.rm = TRUE),
  mean_delta_t15 = mean(delta_t15, na.rm = TRUE)
)]

summary_outdir <- file.path(outdir, "final_PARIS2_ViennaRNA_agreement_summary")
dir.create(summary_outdir, recursive = TRUE, showWarnings = FALSE)

fwrite(
  direct_PARIS2_ViennaRNA_summary,
  file.path(summary_outdir, "direct_PARIS2_ViennaRNA_support_summary.tsv"),
  sep = "\t"
)

fwrite(
  randomization_PARIS2_ViennaRNA_summary,
  file.path(summary_outdir, "matched_randomization_PARIS2_ViennaRNA_summary.tsv"),
  sep = "\t"
)

fwrite(
  transcript_level_PARIS2_ViennaRNA_summary,
  file.path(summary_outdir, "transcript_level_PARIS2_ViennaRNA_observed_vs_random_summary.tsv"),
  sep = "\t"
)

fwrite(
  transcript_random_summary,
  file.path(summary_outdir, "transcript_level_PARIS2_ViennaRNA_observed_vs_random_full_table.tsv"),
  sep = "\t"
)

cat("\n===== FINAL PARIS2 / VIENNARNA AGREEMENT SUMMARY =====\n")

cat(
  "\nDirect support: ",
  format(direct_PARIS2_ViennaRNA_summary[metric == "mapped_duplexes_evaluated", value], big.mark = ","),
  " mapped PARIS2 duplexes across ",
  format(direct_PARIS2_ViennaRNA_summary[metric == "unique_transcripts_evaluated", value], big.mark = ","),
  " transcripts and ",
  format(direct_PARIS2_ViennaRNA_summary[metric == "unique_genes_evaluated", value], big.mark = ","),
  " genes.\n",
  sep = ""
)

cat(
  "Exact direct support: ",
  format(direct_PARIS2_ViennaRNA_summary[metric == "duplexes_supported_exact", value], big.mark = ","),
  " / ",
  format(direct_PARIS2_ViennaRNA_summary[metric == "mapped_duplexes_evaluated", value], big.mark = ","),
  " = ",
  round(direct_PARIS2_ViennaRNA_summary[metric == "percent_supported_exact", value], 1),
  "%.\n",
  sep = ""
)

cat(
  "±15 nt direct support: ",
  format(direct_PARIS2_ViennaRNA_summary[metric == "duplexes_supported_t15", value], big.mark = ","),
  " / ",
  format(direct_PARIS2_ViennaRNA_summary[metric == "mapped_duplexes_evaluated", value], big.mark = ","),
  " = ",
  round(direct_PARIS2_ViennaRNA_summary[metric == "percent_supported_t15", value], 1),
  "%.\n",
  sep = ""
)

cat("\nMatched within-transcript randomization:\n")

for (tol in randomization_PARIS2_ViennaRNA_summary$tolerance) {
  row <- randomization_PARIS2_ViennaRNA_summary[tolerance == tol]

  cat(
    tol,
    ": observed = ",
    round(row$observed_percent_supported, 2),
    "%, expected = ",
    round(row$expected_percent_supported, 2),
    "%, enrichment = ",
    round(row$enrichment_observed_over_expected, 2),
    "x, empirical P = ",
    signif(row$empirical_p, 3),
    ".\n",
    sep = ""
  )
}

cat(
  "\nTranscript-level ±15 nt: ",
  format(transcript_level_PARIS2_ViennaRNA_summary$n_above_random_t15, big.mark = ","),
  " / ",
  format(transcript_level_PARIS2_ViennaRNA_summary$n_transcripts_evaluated, big.mark = ","),
  " transcripts above random expectation = ",
  round(transcript_level_PARIS2_ViennaRNA_summary$percent_above_random_t15, 1),
  "%.\n",
  sep = ""
)

cat(
  "\nFinal summary tables saved to: ",
  summary_outdir,
  "\n",
  sep = ""
)
