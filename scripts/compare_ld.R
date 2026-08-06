#!/usr/bin/env Rscript
#
# compare_ld.R - compare a safeld LD matrix against the original one.
#
# Reads two plink2 --r2-unphased outputs (.vcor) and joins them on the SNP pair,
# keyed by CHROM/POS_A/POS_B rather than by ID. Position keying matters here:
# after a liftover with allele swaps, the ID column can encode the pre-swap
# allele order while REF/ALT hold the post-swap one, so IDs are not a safe join
# key. r2 is sign-invariant, so a swapped pair still compares correctly.
#
# Usage:
#   Rscript compare_ld.R SAFELD.vcor ORIGINAL.vcor [OUT_PREFIX] [EXTRA.vcor EXTRA_NAME]
#
# Prints correlation, attenuation slope and banding diagnostics, and writes
# <OUT_PREFIX>_scatter.png if ggplot2 is available.

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("usage: Rscript compare_ld.R SAFELD.vcor ORIGINAL.vcor [OUT_PREFIX] [EXTRA.vcor EXTRA_NAME]\n")
  quit(status = 1)
}
safeld_file <- args[1]
orig_file   <- args[2]
out_prefix  <- if (length(args) >= 3) args[3] else "ld_compare"
extra_file  <- if (length(args) >= 5) args[4] else NULL
extra_name  <- if (length(args) >= 5) args[5] else NULL

prep <- function(path, name) {
  dt <- fread(path)
  need <- c("#CHROM_A", "POS_A", "POS_B", "UNPHASED_R2")
  miss <- setdiff(need, names(dt))
  if (length(miss)) {
    stop(path, ": missing column(s) ", paste(miss, collapse = ", "),
         ". Expected plink2 --r2-unphased output.")
  }
  out <- data.table(
    pair = paste(dt[["#CHROM_A"]], dt$POS_A, dt$POS_B, sep = "_"),
    r2   = pmin(pmax(dt$UNPHASED_R2, 0), 1)
  )
  setnames(out, "r2", name)
  unique(out, by = "pair")
}

s <- prep(safeld_file, "SAFELD")
o <- prep(orig_file,   "ORIGINAL")
cat(sprintf("pairs: safeld=%d  original=%d\n", nrow(s), nrow(o)))

df <- merge(s, o, by = "pair")
if (!is.null(extra_file)) {
  e <- prep(extra_file, extra_name)
  df <- merge(df, e, by = "pair")
  cat(sprintf("pairs: %s=%d\n", extra_name, nrow(e)))
}
cat(sprintf("pairs in common          : %d\n", nrow(df)))
if (nrow(df) == 0) {
  stop("no shared SNP pairs. The two runs kept different variants; check the ",
       "safeld preprocess log for how many variants it dropped.")
}
only_orig <- nrow(o) - nrow(df)
if (only_orig > 0) {
  cat(sprintf("pairs only in original   : %d  (safeld dropped those variants)\n", only_orig))
}

x <- df$ORIGINAL; y <- df$SAFELD
cat("\n--- SAFELD vs ORIGINAL ---\n")
cat(sprintf("  Pearson r                : %.4f\n", cor(x, y)))
cat(sprintf("  Pearson r2               : %.4f\n", cor(x, y)^2))
cat(sprintf("  Spearman rho             : %.4f\n", cor(x, y, method = "spearman")))
cat(sprintf("  RMSE                     : %.4f\n", sqrt(mean((y - x)^2))))
cat(sprintf("  mean(SAFELD) - mean(ORIG): %+.4f\n", mean(y) - mean(x)))

# Slope through the origin. A faithful reconstruction sits at ~1.0; systematic
# mean-imputation of genotypes drags it below 1 by the squared call rate.
slope <- sum(x * y) / sum(x * x)
cat(sprintf("  attenuation slope (y~0+x): %.4f", slope))
if (slope < 0.9) cat("   <- LD is being lost")
cat("\n")

# Banding check: if the deficit y/x clusters at separated values rather than
# spreading smoothly, distinct groups of variants are attenuated by different
# constant factors, which is what a scatter shows as parallel lines.
strong <- df[x > 0.2]
if (nrow(strong) > 50) {
  ratio <- strong$SAFELD / strong$ORIGINAL
  qs <- quantile(ratio, c(0.05, 0.25, 0.5, 0.75, 0.95))
  cat("\n  r2 ratio (safeld/original) on pairs with original r2 > 0.2:\n")
  cat(sprintf("    5%%=%.3f  25%%=%.3f  median=%.3f  75%%=%.3f  95%%=%.3f\n",
              qs[1], qs[2], qs[3], qs[4], qs[5]))
  h <- hist(ratio, breaks = seq(0, 1.2, by = 0.05), plot = FALSE)
  peaks <- 0
  for (i in seq_along(h$counts)) {
    lo <- if (i > 1) h$counts[i - 1] else 0
    hi <- if (i < length(h$counts)) h$counts[i + 1] else 0
    if (h$counts[i] > lo && h$counts[i] > hi && h$counts[i] > nrow(strong) * 0.03) peaks <- peaks + 1
  }
  cat(sprintf("    distinct modes in the ratio: %d", peaks))
  if (peaks > 1) cat("   <- separated bands, not smooth noise")
  cat("\n")
}

fwrite(df, paste0(out_prefix, "_pairs.tsv.gz"), sep = "\t")
cat(sprintf("\nwrote %s_pairs.tsv.gz\n", out_prefix))

if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  p <- ggplot(df, aes(ORIGINAL, SAFELD)) +
    geom_point(size = 0.35, alpha = 0.25) +
    geom_abline(slope = 1, intercept = 0, colour = "red", linewidth = 0.4) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = "original r2", y = "safeld r2",
         subtitle = sprintf("n=%d pairs   r=%.4f   slope=%.3f", nrow(df), cor(x, y), slope)) +
    theme_classic(base_size = 12)
  ggsave(paste0(out_prefix, "_scatter.png"), p, width = 6, height = 6, dpi = 200)
  cat(sprintf("wrote %s_scatter.png\n", out_prefix))
} else {
  cat("ggplot2 not installed; skipping the scatter plot\n")
}
