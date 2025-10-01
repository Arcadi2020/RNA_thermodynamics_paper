# Script: Generate normalized MFE inputs for the human transcriptome (GENCODE v47)
# Author: Arcadi Cipponi
#
# What this script does (current functionality):
# 1) Downloads the GENCODE v47 GTF (basic annotation) and imports it as a GFF table,
#    filtering to `type == "transcript"`. It removes version suffixes from gene_id and
#    transcript_id and creates a combined gene–transcript key (geneIDtxID).
# 2) Downloads the GENCODE v47 transcriptome FASTA (.fa.gz), decompresses it, reads all
#    transcript sequences, and standardizes FASTA headers to keep only the first field
#    (the Ensembl transcript stable ID with version).
# 3) Loads an Ensembl BioMart export that includes transcript lengths. It then filters
#    that table to the set of transcripts present in the GENCODE FASTA, so lengths and
#    sequences refer to the same transcript universe.
#
# Notes / assumptions:
# - Genome build: GRCh38; annotations: GENCODE v47 (basic).
# - Transcript IDs in the FASTA are expected to match the BioMart "Transcript stable ID
#   version" field (e.g., ENST00000335137.4).
# - The script prepares inputs for downstream MFE/heat-capacity calculations with
#   ViennaRNA (RNAfold / RNAheat). RNAfold’s length limit ~40,000 nt is noted below.
#
# What is left to do (planned):
# - Run ViennaRNA on the FASTA to compute raw MFE per transcript (see provided Bash
#   snippets/files: RNA_MFE_prediction.rtf; RNA_heat_capacity_prediction.rtf).
# - Join RNAfold outputs to the BioMart length table, then compute a normalized MFE
#   metric (e.g., MFE per nucleotide: ΔG_total / length) and, if needed, z-scores
#   across transcripts or per-biotype centering.
# - Exclude or handle transcripts >40,000 nt before RNAfold (split or skip).
# - Save tidy outputs (e.g., CSV/Parquet) with columns: transcript_id_version,
#   gene_id, gene_name, length_nt, MFE_kcal_mol, norm_MFE_kcal_mol_per_nt, etc.
#
# Output (intended):
# - A tidy table aligning each transcript with its sequence, length, and (after the
#   ViennaRNA step) MFE and normalized MFE, ready for downstream analyses and plots.


library(httr)
library(glue)
library(rtracklayer)  
library(dplyr)        
library(stringr)     
library(glue)     

# Set download directory
download_dir <- "Downloads/GTF_annotations"
dir.create(download_dir, recursive = TRUE, showWarnings = FALSE)

# Define correct HTTP-based URL and file name
file_name <- "gencode.v47.basic.annotation.gtf.gz"
gencode_url <- glue("http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/{file_name}")
gencode_path <- file.path(download_dir, file_name)

# Always download the file (overwrite if already exists)
message("Downloading GENCODE GTF file...")
response <- GET(gencode_url, write_disk(gencode_path, overwrite = TRUE))

# Check for successful download
if (response$status_code == 200) {
  message("Download complete.")
} else {
  stop("Download failed with status code: ", response$status_code)
}




# Define path to downloaded GTF file
download_dir <- "Downloads/GTF_annotations"
file_name <- "gencode.v47.basic.annotation.gtf.gz"
gencode_path <- file.path(download_dir, file_name)

# Import GTF file and filter for transcript entries
gencode <- readGFF(gencode_path) %>%
  filter(type == "transcript")

# Remove version numbers from gene_id and transcript_id
gencode <- gencode %>%
  mutate(
    gene_id = str_remove(gene_id, "\\..*$"),
    transcript_id = str_remove(transcript_id, "\\..*$")
  )

# Create combined gene-transcript ID column
gencode <- gencode %>%
  mutate(geneIDtxID = glue("{gene_id}-{transcript_id}"))

# Sort by combined ID
gencode <- arrange(gencode, geneIDtxID)





library(R.utils)       
library(seqinr)        

# Create download directory
download_dir <- "Downloads/GTF_annotations"
dir.create(download_dir, recursive = TRUE, showWarnings = FALSE)

# Define file name and URL
file_name <- "gencode.v47.transcripts.fa.gz"
gencode_url <- "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.transcripts.fa.gz"
gz_path <- file.path(download_dir, file_name)

# Download the transcript FASTA file (overwrite any existing file)
message("Downloading GENCODE transcriptome FASTA...")
GET(gencode_url, write_disk(gz_path, overwrite = TRUE))

# Decompress the file (create a copy and keep the original .gz)
decompressed_path <- gunzip(gz_path, remove = FALSE, overwrite = TRUE)

# Read the decompressed FASTA
message("Reading FASTA...")
Gencode_transcripts_FASTA <- read.fasta(file = decompressed_path, seqtype = "DNA")

# Clean sequence names: keep only the part before the first '|'
cleaned_names <- sapply(strsplit(names(Gencode_transcripts_FASTA), "\\|"), `[`, 1)
names(Gencode_transcripts_FASTA) <- cleaned_names

# Report number of transcripts
message("Number of transcripts loaded: ", length(Gencode_transcripts_FASTA))











# This file has been downloaded manually from Ensembl website
# This file contains the entire list of transcripts present in the Genecode data and unlike to the Gencode transcripts annotation, it contains the transcripts length

# I have followed these steps:
# Via Web: https://www.ensembl.org --> Human --> BioMart at top --> Ensembl Genes --> Human Genes --> Select "Attributes" in left column --> 
# Expand "Gene" (+ sign) in right column --> Select (Check mark) Gene Stable ID/Transcript Stable ID/Gene name/Transctipt length --> Click on count button in top left corner --> 86402 as of today. --> 
# Click on "Results" --> Export All results to file --> Select CSV format (the file will have th extension .txt) --> Click "Go"
Human_transcriptome_GRCh38_attributes <- read.csv("/Volumes/AC_HD_8TB/Exp.\ 277\ Analysis\ of\ Gibbs\ free\ energy\ in\ response\ to\ drug\ selection/R_scripts/Prediction_transcript_MFE/Ensembl_Human_transcriptome_GRCh38/Human_transcriptome_mart_export.txt",
                                                  stringsAsFactors = FALSE) %>%
  dplyr::rename(ensembl_transcript_id_version = Transcript.stable.ID.version) %>%
  dplyr::rename(Transcript_length = Transcript.length..including.UTRs.and.CDS.)
# The file contains 412,243 transcripts

# Check whether all Gencode transcripts are present in the ensembl annotation file
# # Filter the file for the transcripts present in the FASTA file
# # Extract the names from the FASTA headers
fasta_transcript_ids <- names(Gencode_transcripts_FASTA)

# # Filter the data frame based on the transcript IDs from the Gencode FASTA file
Human_transcriptome_GRCh38_attributes_Genecode <- Human_transcriptome_GRCh38_attributes[
  Human_transcriptome_GRCh38_attributes$ensembl_transcript_id_version %in% fasta_transcript_ids, ]







# CALCULATE MFE FOR THE FULL GENCODE TRANSCRIPTOME
#
# Overview:
# We compute minimum free energy (MFE) for all GENCODE v47 transcripts using ViennaRNA.
# Results will later be merged with transcript lengths to derive a normalized metric
# (e.g., ΔG per nucleotide).
#
# Important constraint:
# RNAfold reliably handles transcripts up to ~40,000 nt. Transcripts exceeding this
# length should be excluded or split prior to computation.
#
# How to run the external steps:
# - Use the provided ViennaRNA command wrappers to generate:
#   1) Per-transcript MFE (kcal/mol)
#   2) Heat-capacity profiles (Cp) for downstream thermodynamic analyses
#
# Refer to these scripts:
# - RNA_MFE_prediction.rtf
# - RNA_heat_capacity_prediction.rtf
#
# Note:
# If these are documentation/notes saved as .rtf, convert the commands into executable
# shell scripts (e.g., RNA_MFE_prediction.sh, RNA_heat_capacity_prediction.sh) and
# ensure ViennaRNA is on PATH. After running, import the outputs here and compute
# normalized MFE (ΔG_total / length_nt).











# The predicted MFE and heat capacity values for the human transcriptome 
# have already been computed externally with ViennaRNA and stored in:
#
#   Gencode_TX_norm_MFE_HC_table.csv
#
# Contents:
# - Each row corresponds to one GENCODE transcript.
# - Columns include transcript ID, gene ID, transcript length, raw MFE,
#   normalized MFE (ΔG per nucleotide), and heat capacity (Cp).
# - Normalization at this stage is ONLY by transcript length.
#
# Next steps:
# - This table can be used directly for downstream thermodynamic analyses
#   or further normalized (e.g., z-scores, per-biotype scaling).

Gencode_TX_norm_MFE_HC_table <- read.csv(
  "path_to/Gencode_TX_norm_MFE_HC_table.csv",
  stringsAsFactors = FALSE
)
















# Analysis normalized MFE values ####
# Analysis RNA secondary structures normalized MFE values

Table_ANNOTATED_gene_names <- Gencode_TX_norm_MFE_HC_table
library(ggplot2)

# Scatter plot: Normalized_MFE vs Transcript_length
ggplot(Table_ANNOTATED_gene_names, aes(x = Transcript_length, y = Normalized_MFE)) +
  geom_point(alpha = 0.5, size = 1.5) +  # semi-transparent points
  theme_classic() +
  labs(
    x = "Transcript length (nt)",
    y = expression("Normalized MFE (kcal/mol/nt)"),
    title = "Transcript Length vs Normalized MFE"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )

# The scatter plot of normalized MFE versus transcript length reveals that short transcripts display a wide spread of per-nucleotide MFE values, while longer transcripts converge toward a narrower range. 
# This indicates that although length normalization reduces the direct scaling of folding free energy with transcript size, residual length-dependent bias remains. Such bias likely arises from nonlinear scaling of thermodynamic stability and compositional differences across transcript length ranges.
# Further correction, for example via GC- and length-adjusted residuals or shuffle-based z-scores, would be needed for fully length-independent comparisons.


# Fit and plot the LOESS line so the residual length bias is visually obvious
library(dplyr)
library(ggplot2)

df <- Table_ANNOTATED_gene_names %>%
  filter(is.finite(Transcript_length), is.finite(Normalized_MFE)) %>%
  mutate(L = pmax(Transcript_length, 1))  # guard zeros

# Compute a smooth trend via running medians on log-length bins
trend <- df %>%
  mutate(bin = cut(log1p(L), breaks = 80)) %>%              # ~80 bins; adjust if needed
  group_by(bin) %>%
  summarise(
    x = exp(mean(log1p(L), na.rm = TRUE)) - 1,               # bin center back on original scale
    y = median(Normalized_MFE, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(is.finite(x), is.finite(y))

ggplot(df, aes(x = L, y = Normalized_MFE)) +
  stat_bin2d(bins = 120) +                                   # density instead of points
  geom_line(data = trend, aes(x = x, y = y),                 # smooth running-median trend
            linewidth = 1.2, color = "white") +
  scale_fill_continuous(trans = "log10") +                   # highlight dense areas
  scale_x_log10() +                                          # spreads short/long transcripts
  theme_classic() +
  labs(
    x = "Transcript length (nt, log scale)",
    y = expression("Normalized MFE (kcal/mol/nt)"),
    fill = "Count",
    title = "Transcript Length vs Normalized MFE (binned density + running-median trend)"
  )


# Short RNAs (<200 nt): The median normalized MFE is highly variable and shifts abruptly. This volatility likely reflects both the small number of base pairs possible and compositional differences among very short RNAs (e.g., small ncRNAs, partial transcripts).
# Intermediate lengths (~200 nt–1 kb): The median rapidly stabilizes toward ~–0.28 kcal/mol/nt.
# Long RNAs (>1 kb): The trend flattens further, but there is still a slight upward drift (less negative per-nt MFE) as length increases, suggesting residual length-dependent bias even after per-nt normalization.
# The flattening in the middle range indicates that normalization removes much of the first-order scaling, but the small slope for long transcripts suggests systematic effects remain—possibly due to:
# Composition shifts with transcript size (e.g., GC enrichment or depletion in long transcripts).
# Structural modularity in large RNAs that alters per-nt folding energetics.
# ViennaRNA’s algorithmic limitations for very large molecules (ensemble approximations, energy model scaling).


# This plot confirms that simple length normalization reduces but does not fully eliminate length-dependent structure in normalized MFE values. 
# Short RNAs display inherently greater variance in folding stability per nt, while long RNAs converge toward a narrower range with a mild upward bias. 
# This issue has been addressed in the manuscript of E. Trotta:
# Trotta, E., On the normalization of the minimum free energy of RNAs by sequence length, PLOS ONE 2014


# Here below, I normalize the transcripts that are , 600 nt using the table and the formula provided in the manuscript of E. Trotta
# This code generates the normalized MFE values for the entire Human transcriptome using the Trotta approach: 
# Trotta, E., On the normalization of the minimum free energy of RNAs by sequence length, PLOS ONE 2014

# Background:
# The raw folding energy of a transcript of length L (nt) is the minimum free energy (MFE), measured in kcal/mol.
#
# A common normalization is the "average MFE per 100 nt" (AMFE):
#
#     AMFE = (100 * MFE) / L
#
# However, this naive normalization introduces a systematic bias
# for short transcripts (< ~600 nt), which appear artificially
# more stable simply because of scaling with length.
#
# Trotta (PLOS ONE, 2014) proposed a correction to remove this
# length-dependent bias by referencing an empirical dataset of
# expected folding energies for different transcript lengths.
#
# Formula:
# For transcripts with 40 <= L <= 600 nt:
#
#     MFE_den = 100 * ( MFE - MFE_ref(L) ) / (L - 8)
#
# where:
#   - MFE           : raw folding energy of the transcript
#   - L             : transcript length (nt)
#   - MFE_ref(L)    : expected folding energy for length L
#                     (from Trotta's empirical reference curve)
#   - (L - 8)       : adjustment because the first ~8 nt
#                     contribute little to folding stability
#
# Key points:
# - For L > 600 nt, transcripts are left unchanged (correction was validated only for shorter RNAs).
# - Subtracting MFE_ref(L) removes the baseline length effect.
# - Dividing by (L - 8) produces a corrected, per-length metric that is comparable across transcripts.
#
# Interpretation:
#   MFE_den quantifies how much more or less stable a transcript is relative to the expected baseline at that length. 
# It is a length-bias–corrected measure of folding energy.



library(dplyr)
library(ggplot2)
library(scales)
suppressPackageStartupMessages(require(mgcv))      # for GAM in ggpubr/geom_smooth
use_hex <- requireNamespace("hexbin", quietly = TRUE)





## Load the Gencode transcript MFE table

# No need to escape spaces when the path is inside quotes
# Gencode_TX_norm_MFE_HC_table <- read.csv(
#   "/Volumes/AC_HD_8TB/Exp. 277 Analysis of Gibbs free energy in response to drug selection/R_scripts/Prediction_transcript_MFE/Gencode_TX_norm_MFE_HC_table.csv",
#   stringsAsFactors = FALSE
# )





## Helpers

# Safe fetch of a numeric column; returns NULL if absent
get_col <- function(df, nm) if (nm %in% names(df)) df[[nm]] else NULL

# Derive MFE if missing (from Normalized_MFE and length)
derive_MFE_if_needed <- function(df, len_col = "Transcript_length",
                                 nmfe_col = "Normalized_MFE", mfe_col = "MFE") {
  L   <- get_col(df, len_col)
  NM  <- get_col(df, nmfe_col)
  MFE <- get_col(df, mfe_col)
  
  if (!is.null(MFE)) return(df)  # already present
  
  if (!is.null(L) && !is.null(NM)) {
    # NM is kcal/mol per 100 nt (AMFE). So MFE = NM * L / 100
    df[[mfe_col]] <- as.numeric(NM) * as.numeric(L) / 100
    return(df)
  }
  
  stop("Neither 'MFE' nor a combination of 'Normalized_MFE' and 'Transcript_length' is available to compute MFE.")
}

# Median trend on bins with a cap on number of bins
median_trend <- function(df, x = "L", y = "Y",
                         from = NULL, to = NULL, target_bin_width = 10, max_bins = 2000) {
  xv <- df[[x]]
  yv <- df[[y]]
  ok <- is.finite(xv) & is.finite(yv)
  xv <- xv[ok]; yv <- yv[ok]
  
  if (is.null(from)) from <- floor(min(xv, na.rm = TRUE))
  if (is.null(to))   to   <- ceiling(max(xv, na.rm = TRUE))
  
  n_bins <- ceiling((to - from) / target_bin_width)
  n_bins <- max(70, min(n_bins, max_bins))  # keep at least 70, cap at max_bins
  
  brks <- seq(from, to, length.out = n_bins + 1)
  cutx <- cut(xv, breaks = brks, include.lowest = TRUE)
  out <- dplyr::tibble(bin = cutx, x = xv, y = yv) |>
    dplyr::group_by(bin) |>
    dplyr::summarise(x = mean(x, na.rm = TRUE),
                     y = median(y, na.rm = TRUE),
                     .groups = "drop") |>
    dplyr::filter(is.finite(x), is.finite(y))
  out
}

## Prep main DF

Gencode_TX_norm_MFE_HC_table <- derive_MFE_if_needed(Gencode_TX_norm_MFE_HC_table)

df_full <- Gencode_TX_norm_MFE_HC_table %>%
  filter(is.finite(Transcript_length),
         is.finite(Normalized_MFE)) %>%
  mutate(L = pmax(Transcript_length, 0))

maxL <- max(df_full$L, na.rm = TRUE)









## Full-range plot (binned density + running-median trend)

trend_full <- median_trend(
  df_full %>% transmute(L = L, Y = Normalized_MFE),
  x = "L", y = "Y",
  from = 0, to = maxL,
  target_bin_width = 10,  # ~10-nt bins
  max_bins = 2000
)

p_full <- ggplot(df_full, aes(x = L, y = Normalized_MFE)) +
  { if (use_hex) ggplot2::geom_hex(bins = 150) else ggplot2::stat_bin2d(bins = 150) } +
  geom_line(data = trend_full, aes(x = x, y = y),
            linewidth = 0.6, color = "white") +
  scale_fill_continuous(trans = "log10", labels = label_number(accuracy = 1), name = "Count") +
  coord_cartesian(xlim = c(0, maxL), ylim = c(-1.3, 0)) +
  theme_classic() +
  labs(
    x = "Transcript length (nt)",
    y = expression("Normalized MFE (kcal/mol/100 nt)"),
    title = "Transcript length vs Normalized MFE (full transcriptome)"
  )
print(p_full)











## Short-range zoom (0–700 nt)

df_short <- df_full %>% filter(L <= 700)

trend_short <- median_trend(
  df_short %>% transmute(L = L, Y = Normalized_MFE),
  x = "L", y = "Y",
  from = 0, to = 700,
  target_bin_width = 10
)

p_short <- ggplot(df_short, aes(x = L, y = Normalized_MFE)) +
  { if (use_hex) ggplot2::geom_hex(bins = 120) else ggplot2::stat_bin2d(bins = 120) } +
  geom_line(data = trend_short, aes(x = x, y = y),
            linewidth = 0.6, color = "white") +
  geom_vline(xintercept = 600, linetype = "dashed") +
  scale_fill_continuous(trans = "log10", name = "Count") +
  coord_cartesian(xlim = c(0, 700), ylim = c(-1.3, 0)) +
  theme_classic() +
  labs(
    x = "Transcript length (nt)",
    y = expression("Normalized MFE (kcal/mol/100 nt)"),
    title = "Transcript length vs Normalized MFE (0–700 nt)"
  )
print(p_short)














## Trotta correction (40–600 nt only; unchanged >600)

# Load Trotta reference and interpolator
trotta_ref <- read.csv(
  "path_to/Dataset_S1_Edoardo_Trotta.csv",
  stringsAsFactors = FALSE
) %>%
  dplyr::rename(
    Transcript_length = L,
    MFEref = MFE
  ) %>%
  arrange(Transcript_length)


mfe_ref_fun <- stats::approxfun(
  x = trotta_ref$Transcript_length,
  y = trotta_ref$MFEref,
  rule = 1  # NA outside [min,max]; keeps >600 uncorrected
)

Complete_table_MFEden <- Gencode_TX_norm_MFE_HC_table %>%
  mutate(
    L = Transcript_length,
    AMFE = 100 * MFE / L,                        # baseline "normalized" (kcal/mol/100 nt)
    MFEref_interp = mfe_ref_fun(L),
    MFEden = dplyr::case_when(
      !is.finite(L) | !is.finite(MFE) ~ NA_real_,
      L <= 8                          ~ NA_real_,
      L < 40                          ~ NA_real_,
      is.na(MFEref_interp)            ~ NA_real_,  # e.g., L > 600
      TRUE ~ 100 * (MFE - MFEref_interp) / (L - 8) # Trotta formula
    )
  ) %>%
  dplyr::select(ensembl_transcript_id_version, Transcript_length, MFE, Normalized_MFE, AMFE, MFEden)

# Valid plotting range for MFEden
df_plot <- Complete_table_MFEden %>%
  filter(is.finite(Transcript_length), is.finite(MFEden), Transcript_length <= 700)

trend_plot <- median_trend(
  df_plot %>% transmute(L = Transcript_length, Y = MFEden),
  x = "L", y = "Y", from = 40, to = 600, target_bin_width = 10
)

p_den <- ggplot(df_plot, aes(x = Transcript_length, y = MFEden)) +
  { if (use_hex) ggplot2::geom_hex(bins = 120) else ggplot2::stat_bin2d(bins = 120) } +
  geom_line(data = trend_plot, aes(x = x, y = y), linewidth = 0.6, color = "white") +
  geom_vline(xintercept = 600, linetype = "dashed") +
  scale_fill_continuous(trans = "log10", name = "Count") +
  coord_cartesian(xlim = c(0, 700)) +
  theme_classic() +
  labs(
    x = "Transcript length (nt)",
    y = expression(MFE[den]~"(kcal/mol/100 nt)"),
    title = "Transcript length vs MFEden (valid range; binned density + running median)"
  )
print(p_den)







# Quick correlation within valid range
cor_test_mfeden <- with(
  subset(Complete_table_MFEden, !is.na(MFEden)),
  cor.test(Transcript_length, MFEden, method = "pearson")
)
print(cor_test_mfeden)

## Anchored short-length adjustment (≤600 adjusted; >600 unchanged)

set.seed(123)
df_all <- Gencode_TX_norm_MFE_HC_table %>%
  filter(is.finite(Transcript_length), Transcript_length > 8,
         is.finite(Normalized_MFE)) %>%
  transmute(
    ensembl_transcript_id_version,
    L   = Transcript_length,
    NM0 = Normalized_MFE       # source metric to preserve >600
  )

# Fit GAM on ≤600 only
df_fit <- df_all %>% filter(L <= 600)
N_FIT   <- 10000
K_BASIS <- 40
if (nrow(df_fit) > N_FIT) df_fit <- dplyr::slice_sample(df_fit, n = N_FIT)

gam_short <- mgcv::gam(NM0 ~ s(L, k = K_BASIS), data = df_fit, method = "REML")

anchor_val <- as.numeric(predict(gam_short, newdata = data.frame(L = 600)))
trend_short_full <- as.numeric(predict(gam_short, newdata = data.frame(L = df_all$L)))

df_all <- df_all %>%
  mutate(
    trend_short = trend_short_full,
    delta       = dplyr::if_else(L <= 600, trend_short - anchor_val, 0.0),
    NM_adj      = NM0 - delta
  )

# Hard invariant check: >600 unchanged
idx_gt600 <- df_all$L > 600
stopifnot(identical(df_all$NM_adj[idx_gt600], df_all$NM0[idx_gt600]))

# Lightweight backgrounds
BINS_COARSE <- 40
MAX_POINTS  <- 20000

bg_before <- if (use_hex) list(geom_hex(bins = BINS_COARSE)) else {
  set.seed(1)
  df_sub <- dplyr::slice_sample(df_all, n = min(MAX_POINTS, nrow(df_all)))
  list(geom_point(data = df_sub, aes(L, NM0), shape = 16, size = 0.4, alpha = 0.15))
}

bg_after <- if (use_hex) list(geom_hex(bins = BINS_COARSE)) else {
  set.seed(2)
  df_sub <- dplyr::slice_sample(df_all, n = min(MAX_POINTS, nrow(df_all)))
  list(geom_point(data = df_sub, aes(L, NM_adj), shape = 16, size = 0.4, alpha = 0.15))
}

# Trend lines for before/after
Lmin <- max(40, floor(min(df_all$L, na.rm = TRUE) / 10) * 10)
Lmax <- ceiling(max(df_all$L, na.rm = TRUE) / 10) * 10

L_grid_short <- data.frame(L = seq(Lmin, 600, by = 5))
trend_before_short <- transform(L_grid_short,
                                y = as.numeric(predict(gam_short, newdata = L_grid_short)))
trend_before_long <- median_trend(df_all %>% transmute(L, Y = NM0),
                                  x = "L", y = "Y", from = 600, to = Lmax, target_bin_width = 25)
trend_before_full <- dplyr::bind_rows(dplyr::rename(trend_before_short, x = L), trend_before_long)

trend_after_short <- data.frame(x = c(Lmin, 600), y = rep(anchor_val, 2))
trend_after_long  <- median_trend(df_all %>% transmute(L, Y = NM_adj),
                                  x = "L", y = "Y", from = 600, to = Lmax, target_bin_width = 25)
trend_after_full  <- dplyr::bind_rows(trend_after_short, trend_after_long)





# Plot the normalized values before/after the Trotta normalization on transcripts < 600 nts with full-width trend lines
p_before <- ggplot(df_all, aes(L, NM0)) +
  bg_before +
  geom_line(data = trend_before_full, aes(x, y), color = "white", linewidth = 0.6) +
  geom_vline(xintercept = 600, linetype = "dashed") +
  { if (use_hex) scale_fill_continuous(trans = "log10", name = "Count") else scale_y_continuous() } +
  coord_cartesian(xlim = c(Lmin, Lmax)) +
  labs(x = "Sequence length (nt)", y = "Normalized MFE (original)",
       title = "Before: original normalized values + full-width trend") +
  theme_classic()
print(p_before)




p_after <- ggplot(df_all, aes(L, NM_adj)) +
  bg_after +
  geom_line(data = trend_after_full, aes(x, y), color = "white", linewidth = 0.6) +
  geom_vline(xintercept = 600, linetype = "dashed") +
  { if (use_hex) scale_fill_continuous(trans = "log10", name = "Count") else scale_y_continuous() } +
  coord_cartesian(xlim = c(Lmin, Lmax)) +
  labs(x = "Sequence length (nt)", y = "Normalized MFE (≤600 adjusted)",
       title = "After: anchored correction (≤600 adjusted; >600 unchanged) + full-width trend") +
  theme_classic()
print(p_after)






# Zoomed plots for transcripts < 600 nt (before vs after)

# Lightweight backgrounds restricted to L <= 600
bg_before_short <- if (use_hex) {
  list(geom_hex(bins = 60))
} else {
  set.seed(11)
  df_sub <- dplyr::slice_sample(dplyr::filter(df_all, L <= 600),
                                n = min(MAX_POINTS, sum(df_all$L <= 600)))
  list(geom_point(data = df_sub, aes(L, NM0), shape = 16, size = 0.4, alpha = 0.15))
}

bg_after_short <- if (use_hex) {
  list(geom_hex(bins = 60))
} else {
  set.seed(12)
  df_sub <- dplyr::slice_sample(dplyr::filter(df_all, L <= 600),
                                n = min(MAX_POINTS, sum(df_all$L <= 600)))
  list(geom_point(data = df_sub, aes(L, NM_adj), shape = 16, size = 0.4, alpha = 0.15))
}

# BEFORE (zoom < 600 nt): use the short-range trend you already computed
p_before_zoom <- ggplot(dplyr::filter(df_all, L <= 600), aes(L, NM0)) +
  bg_before_short +
  geom_line(data = dplyr::rename(trend_before_short, x = L), aes(x, y),
            color = "white", linewidth = 0.6) +
  geom_vline(xintercept = 600, linetype = "dashed") +
  { if (use_hex) scale_fill_continuous(trans = "log10", name = "Count") else scale_y_continuous() } +
  coord_cartesian(xlim = c(Lmin, 600)) +
  labs(x = "Sequence length (nt)", y = "Normalized MFE (original)",
       title = "Before (zoom < 600 nt): original normalized values + trend") +
  theme_classic()
print(p_before_zoom)

# AFTER (zoom < 600 nt): flat anchored trend between [Lmin, 600]
p_after_zoom <- ggplot(dplyr::filter(df_all, L <= 600), aes(L, NM_adj)) +
  bg_after_short +
  geom_line(data = trend_after_short, aes(x, y),
            color = "white", linewidth = 0.6) +
  geom_vline(xintercept = 600, linetype = "dashed") +
  { if (use_hex) scale_fill_continuous(trans = "log10", name = "Count") else scale_y_continuous() } +
  coord_cartesian(xlim = c(Lmin, 600)) +
  labs(x = "Sequence length (nt)", y = "Normalized MFE (≤600 adjusted)",
       title = "After (zoom < 600 nt): anchored correction + trend") +
  theme_classic()
print(p_after_zoom)







# Map adjusted values back to the original table
adj_map <- df_all %>%
  transmute(ensembl_transcript_id_version, Transcript_length = L, Normalized_MFE_new = NM_adj)


# 1) Diagnose key uniqueness on both sides
x_dups <- Gencode_TX_norm_MFE_HC_table %>%
  count(ensembl_transcript_id_version, Transcript_length, name = "n") %>%
  filter(n > 1)

y_dups <- adj_map %>%
  count(ensembl_transcript_id_version, Transcript_length, name = "n") %>%
  filter(n > 1)

if (nrow(y_dups) > 0) {
  message("De-duplicating adj_map: found ", nrow(y_dups), " duplicated keys")
}

# 2) Force adj_map to one row per (id, length)
adj_map_dedup <- adj_map %>%
  arrange(ensembl_transcript_id_version, Transcript_length) %>%
  distinct(ensembl_transcript_id_version, Transcript_length, .keep_all = TRUE)

# (Optional, stronger) Or aggregate if duplicates exist with identical NM_adj:
# adj_map_dedup <- adj_map %>%
#   group_by(ensembl_transcript_id_version, Transcript_length) %>%
#   summarise(Normalized_MFE_new = first(Normalized_MFE_new), .groups = "drop")




Gencode_TX_norm_MFE_HC_table_final <- Gencode_TX_norm_MFE_HC_table %>%
  left_join(
    adj_map_dedup,
    by = c("ensembl_transcript_id_version", "Transcript_length"),
    relationship = "many-to-one"    # will error if still not unique in adj_map
  ) %>%
  mutate(
    Normalized_MFE = if_else(
      Transcript_length <= 600 & !is.na(Normalized_MFE_new),
      Normalized_MFE_new,
      Normalized_MFE
    )
  ) %>%
  dplyr::select(-Normalized_MFE_new)


# NM_adj really equals NM0 for >600 in the generator df
stopifnot(max(abs(df_all$NM_adj[df_all$L > 600] - df_all$NM0[df_all$L > 600])) == 0)

# Final table unchanged for >600
idx_gt600_final <- Gencode_TX_norm_MFE_HC_table$Transcript_length > 600
stopifnot(identical(
  Gencode_TX_norm_MFE_HC_table_final$Normalized_MFE[idx_gt600_final],
  Gencode_TX_norm_MFE_HC_table$Normalized_MFE[idx_gt600_final]
))





library(ggplot2)

# Build a comparison table
compare_df <- Gencode_TX_norm_MFE_HC_table %>%
  dplyr::select(ensembl_transcript_id_version, Transcript_length, Normalized_MFE_orig = Normalized_MFE) %>%
  left_join(
    Gencode_TX_norm_MFE_HC_table_final %>%
      dplyr::select(ensembl_transcript_id_version, Transcript_length, Normalized_MFE_new = Normalized_MFE),
    by = c("ensembl_transcript_id_version", "Transcript_length")
  )

# Scatter plot: old vs new
ggplot(compare_df, aes(x = Normalized_MFE_orig, y = Normalized_MFE_new,
                       color = Transcript_length <= 600)) +
  geom_point(alpha = 0.4, size = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  scale_color_manual(values = c("TRUE" = "blue", "FALSE" = "red"),
                     name = "Transcript length ≤ 600") +
  theme_classic() +
  labs(
    x = "Original Normalized MFE",
    y = "Trotta-corrected Normalized MFE",
    title = "Comparison of Normalized MFE (original vs. Trotta correction)"
  )

# The plot shows the changes of Normalized MFE after applying the Trotta correction.
# The blues dots refer to the transcripts with length ≤ 600 nts, which have been adjusted, while the red dots refer to the transcripts with length > 600 nts, which remain unchanged.


# What fraction of the transcriptome are the transcripts < 600 nts?
library(dplyr)

# Assuming your working table is Gencode_TX_norm_MFE_HC_table_final
transcript_stats <- Gencode_TX_norm_MFE_HC_table_final %>%
  mutate(Group = if_else(Transcript_length <= 600, "<=600 nt", ">600 nt")) %>%
  count(Group) %>%
  mutate(Percent = 100 * n / sum(n))

print(transcript_stats)


# The transcripts with length ≤ 600 nts are 29% of the total transcriptome, while the transcripts with length > 600 nts are 71% of the total transcriptome.


write.csv(Gencode_TX_norm_MFE_HC_table_final, "path_to/Gencode_TX_norm_MFE_HC_table.csv", row.names = FALSE)


# This is the final table normalized only by transcript length and by Trotta approach (for transcripts < 600 nts)










































