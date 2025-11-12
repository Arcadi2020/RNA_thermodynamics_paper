# RNA folding stability & heat-capacity (Cp) integration (GENCODE)
# Input files: predicted MFE + Cp tables, integrate Cp (35–39 °C),
#              diagnose residual length bias, apply Trotta normalization
#              (40–600 nt) and an anchored ≤600-nt correction, export final tables.

# 1) Load input tables:
#    - Heat_capacity_complete_table: per-transcript Normalized_MFE, Transcript_length, and Cp columns (Cp_35C..Cp_39C).
#    - Dataset_S1_Edoardo_Trotta.csv: reference MFE_ref(L) curve (Trotta 2014).
#
# 2) Integrate Cp over 35–39 °C (trapezoidal rule):
#    - Adds `integrated_Cp_trapz` per transcript using pracma::trapz on
#      Cp_35C..Cp_39C; robust to missing temperatures.
#
# 3) Exploratory diagnostics of Normalized_MFE vs length:
#    - Binned density + running-median trend (0–all nt; zoom 0–700 nt).
#    - Confirms residual length-dependence after per-nt normalization.
#
# 4) Trotta normalization (valid for 40–600 nt):
#    - Interpolates MFE_ref(L); computes MFEden = 100 * (MFE - MFEref) / (L - 8).
#
# OUTPUTS
# - "Gencode_TX_norm_MFE_HC_table.csv": final table containing:
#     * Normalized_MFE (anchored correction applied for L ≤ 600 only)
#     * integrated_Cp_trapz (35–39 °C)
#     * Trotta normalization for transcripts 40–600 nt


library(seqinr)
library(dplyr)
library(devtools)
library(tidyverse)
library(Biostrings)
library(FastaUtils)
library(biomaRt)
library(cbioportalR)
library(tibble)
library(R.utils)
library(rtracklayer)




# Integrate Cp over 35–39 °C (trapezoid)
library(pracma)

temps <- 35:39
cp_cols <- paste0("Cp_", temps, "C")

# Safety: ensure Cp columns exist (in case a file had no values at some T)
missing_cols <- setdiff(cp_cols, names(Heat_capacity_complete_table))
if (length(missing_cols)) {
  # add missing columns as NA so the mutate below doesn't error
  Heat_capacity_complete_table[missing_cols] <- NA_real_
}

Heat_capacity_complete_table_IHC <- Heat_capacity_complete_table %>%
  rowwise() %>%
  mutate(
    integrated_Cp_trapz = {
      xi <- temps
      yi <- as.numeric(c_across(all_of(cp_cols)))
      ok <- is.finite(yi)
      if (sum(ok) >= 2) pracma::trapz(xi[ok], yi[ok]) else NA_real_
    }
  ) %>%
  ungroup()

# Optional: quick QC
cat("Rows with integrated_Cp_trapz = NA:",
    sum(is.na(Heat_capacity_complete_table_IHC$integrated_Cp_trapz)), "\n")


Gencode_TX_norm_MFE_HC_table <- Heat_capacity_complete_table


write.csv(Gencode_TX_norm_MFE_HC_table, "path_to/Gencode_TX_norm_MFE_HC_table.csv", row.names = FALSE)
#  This is the table with MFE values normalized by transcript length and with the heat capacity values across the T range 35 - 39 ˚C range













# Analysis normalized MFE values ####

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
# This indicates that although length normalization reduces the direct scaling of folding free energy with transcript size, residual length-dependent bias remains. 


# Fit and plot the LOESS line so the residual length bias is visually obvious
library(dplyr)
library(ggplot2)

df <- Table_ANNOTATED_gene_names %>%
  filter(is.finite(Transcript_length), is.finite(Normalized_MFE)) %>%
  mutate(L = pmax(Transcript_length, 1))  # guard zeros

# Compute a smooth trend via running medians on log-length bins
trend <- df %>%
  mutate(bin = cut(log1p(L), breaks = 80)) %>%              
  group_by(bin) %>%
  summarise(
    x = exp(mean(log1p(L), na.rm = TRUE)) - 1,               
    y = median(Normalized_MFE, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(is.finite(x), is.finite(y))

ggplot(df, aes(x = L, y = Normalized_MFE)) +
  stat_bin2d(bins = 120) +                                  
  geom_line(data = trend, aes(x = x, y = y),                
            linewidth = 1.2, color = "white") +
  scale_fill_continuous(trans = "log10") +                   
  scale_x_log10() +                                          
  theme_classic() +
  labs(
    x = "Transcript length (nt, log scale)",
    y = expression("Normalized MFE (kcal/mol/nt)"),
    fill = "Count",
    title = "Transcript Length vs Normalized MFE (binned density + running-median trend)"
  )



# The plot confirms that simple length normalization reduces but does not fully eliminate length-dependent structure in normalized MFE values. 
# Short RNAs display inherently greater variance in folding stability per nt, while long RNAs converge toward a narrower range with a mild upward bias. 
# This issue has been addressed in the manuscript of E. Trotta:
# Trotta, E., On the normalization of the minimum free energy of RNAs by sequence length, PLOS ONE 2014

# Here below, I normalize the transcripts that are , 600 nt using the table and the formula provided in the manuscript of E. Trotta
# This code generates the normalized MFE values for the entire Human transcriptome using the Trotta approach: 
# Trotta, E., On the normalization of the minimum free energy of RNAs by sequence length, PLOS ONE 2014



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







# Correlation within valid range
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



write.csv(Gencode_TX_norm_MFE_HC_table_final, "path_to/Gencode_TX_norm_MFE_HC_table.csv", row.names = FALSE)
# This is the final table normalized only by transcript length and by Trotta approach (for transcripts < 600 nts)







