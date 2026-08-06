#!/usr/bin/env Rscript
#
# check_cor.R - compare the LD of safeld output against the original genotypes.
#
# Adapted from the three-way ggpairs version. The R implementation is dropped
# from the comparison because it cannot read these VCFs at all: 100% of records
# carry no INFO/AF, so its MAF filter turns every frequency into NA and keeps
# zero variants. Pass a third file explicitly if you do have one.
#
# Pairs are keyed on CHROM/POS_A/POS_B rather than on variant ID. After a
# liftover with allele swaps the ID column can hold the pre-swap allele order
# while REF/ALT hold the post-swap one, so IDs are not a safe join key. r2 is
# sign-invariant, so a swapped pair still compares correctly.
#
# Usage:
#   Rscript check_cor.R                                   # uses the defaults below
#   Rscript check_cor.R SAFELD.vcor ORIGINAL.vcor PREFIX
#   Rscript check_cor.R SAFELD.vcor ORIGINAL.vcor PREFIX EXTRA.vcor EXTRA_NAME

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args        <- commandArgs(trailingOnly = TRUE)
safeld_file <- if (length(args) >= 1) args[1] else "SAFELD_pge.vcor"
orig_file   <- if (length(args) >= 2) args[2] else "ORIGINAL_pge.vcor"
out_prefix  <- if (length(args) >= 3) args[3] else "ld_compare"
extra_file  <- if (length(args) >= 5) args[4] else NULL
extra_name  <- if (length(args) >= 5) args[5] else NULL

prep <- function(path, name) {
  if (!file.exists(path)) stop("no such file: ", path)
  dt   <- fread(path, showProgress = FALSE)
  need <- c("#CHROM_A", "POS_A", "POS_B", "UNPHASED_R2")
  miss <- setdiff(need, names(dt))
  if (length(miss)) {
    stop(path, ": missing column(s) ", paste(miss, collapse = ", "),
         "\nExpected plink2 --r2-unphased output. Columns found: ",
         paste(names(dt), collapse = ", "))
  }
  out <- data.table(
    snp_pair = paste(dt[["#CHROM_A"]], dt$POS_A, dt$POS_B, sep = "_"),
    r2       = pmin(pmax(dt$UNPHASED_R2, 0), 1)
  )
  setnames(out, "r2", name)
  unique(out, by = "snp_pair")
}

cat("reading...\n")
s <- prep(safeld_file, "SAFELD")
o <- prep(orig_file,   "ORIGINAL")
cat(sprintf("  SAFELD   : %9d pairs   (%s)\n", nrow(s), basename(safeld_file)))
cat(sprintf("  ORIGINAL : %9d pairs   (%s)\n", nrow(o), basename(orig_file)))

df <- merge(s, o, by = "snp_pair")
if (!is.null(extra_file)) {
  e  <- prep(extra_file, extra_name)
  cat(sprintf("  %-9s: %9d pairs\n", extra_name, nrow(e)))
  df <- merge(df, e, by = "snp_pair")
}

cat(sprintf("  in common: %9d pairs\n\n", nrow(df)))
if (nrow(df) == 0) {
  stop("no shared SNP pairs. The two runs kept different variants; check how ",
       "many variants safeld preprocess reported emitting.")
}
lost <- nrow(o) - nrow(df)
if (lost > 0) {
  cat(sprintf("NOTE: %d pair(s) present in ORIGINAL but not in SAFELD (%.1f%%).\n",
              lost, 100 * lost / nrow(o)))
  cat("      If you did not pass --ld-window-r2 0 to both plink2 calls, the\n")
  cat("      pairs safeld attenuated below the default threshold have deleted\n")
  cat("      themselves from this comparison.\n\n")
}

x <- df$ORIGINAL
y <- df$SAFELD

# Slope through the origin: 1.0 means faithful, below 1 means LD is being lost.
slope <- sum(x * y) / sum(x * x)

cat("--- SAFELD vs ORIGINAL ---\n")
cat(sprintf("  Pearson r                 : %.4f\n", cor(x, y)))
cat(sprintf("  Pearson r2                : %.4f\n", cor(x, y)^2))
cat(sprintf("  Spearman rho              : %.4f\n", cor(x, y, method = "spearman")))
cat(sprintf("  RMSE                      : %.4f\n", sqrt(mean((y - x)^2))))
cat(sprintf("  mean(SAFELD) - mean(ORIG) : %+.4f\n", mean(y) - mean(x)))
cat(sprintf("  attenuation slope (y~0+x) : %.4f%s\n", slope,
            if (slope < 0.9) "   <- LD is being lost" else ""))

# Attenuation across the range: uniform shrinkage looks flat here, whereas a
# subset of damaged variants drags only some deciles down.
cat("\n  by original r2 decile:\n")
cat(sprintf("    %-14s %8s %10s %10s %8s\n", "orig r2 range", "n", "mean orig", "mean safeld", "ratio"))
brk <- seq(0, 1, by = 0.1)
grp <- cut(x, brk, include.lowest = TRUE)
for (lv in levels(grp)) {
  i <- which(grp == lv)
  if (length(i) < 20) next
  mo <- mean(x[i]); ms <- mean(y[i])
  cat(sprintf("    %-14s %8d %10.4f %10.4f %8.3f\n", lv, length(i), mo, ms, ms / mo))
}

# Banding check: separated modes in the per-pair ratio mean distinct groups of
# variants are each attenuated by their own constant factor, which is what a
# scatter renders as parallel lines.
strong <- df[x > 0.2]
if (nrow(strong) > 50) {
  ratio <- strong$SAFELD / strong$ORIGINAL
  qs <- quantile(ratio, c(0.05, 0.25, 0.5, 0.75, 0.95))
  cat("\n  ratio safeld/original on pairs with original r2 > 0.2:\n")
  cat(sprintf("    5%%=%.3f  25%%=%.3f  median=%.3f  75%%=%.3f  95%%=%.3f\n",
              qs[1], qs[2], qs[3], qs[4], qs[5]))
  h <- hist(ratio, breaks = seq(0, 1.5, by = 0.05), plot = FALSE)
  peaks <- 0
  for (i in seq_along(h$counts)) {
    lo <- if (i > 1) h$counts[i - 1] else 0
    hi <- if (i < length(h$counts)) h$counts[i + 1] else 0
    if (h$counts[i] > lo && h$counts[i] > hi && h$counts[i] > nrow(strong) * 0.03) peaks <- peaks + 1
  }
  cat(sprintf("    distinct modes: %d%s\n", peaks,
              if (peaks > 1) "   <- separated bands, not smooth noise" else "   <- single mode"))
}

fwrite(df, paste0(out_prefix, "_pairs.tsv.gz"), sep = "\t")
cat(sprintf("\nwrote %s_pairs.tsv.gz\n", out_prefix))

# ---- plots ----------------------------------------------------------------
sub <- sprintf("n=%s pairs   r=%.4f   slope=%.3f",
               format(nrow(df), big.mark = ","), cor(x, y), slope)

p <- ggplot(df, aes(ORIGINAL, SAFELD)) +
  geom_point(size = 0.3, alpha = 0.15) +
  geom_abline(slope = 1, intercept = 0, colour = "red", linewidth = 0.4) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "original r2", y = "safeld r2", subtitle = sub) +
  theme_classic(base_size = 12)
ggsave(paste0(out_prefix, "_scatter.png"), p, width = 6, height = 6, dpi = 200)
cat(sprintf("wrote %s_scatter.png\n", out_prefix))

# Density version: with ~10^6 pairs an alpha scatter saturates and hides where
# the mass actually sits.
if (requireNamespace("hexbin", quietly = TRUE)) {
  ph <- ggplot(df, aes(ORIGINAL, SAFELD)) +
    geom_hex(bins = 100) +
    scale_fill_viridis_c(trans = "log10") +
    geom_abline(slope = 1, intercept = 0, colour = "red", linewidth = 0.4) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = "original r2", y = "safeld r2", subtitle = sub) +
    theme_classic(base_size = 12)
  ggsave(paste0(out_prefix, "_hex.png"), ph, width = 6.5, height = 6, dpi = 200)
  cat(sprintf("wrote %s_hex.png\n", out_prefix))
}

if (!is.null(extra_file) && requireNamespace("GGally", quietly = TRUE)) {
  library(GGally)
  pp <- ggpairs(
    as.data.frame(df[, c("SAFELD", extra_name, "ORIGINAL"), with = FALSE]),
    lower = list(continuous = wrap("points", alpha = 0.15, size = 0.2)),
    upper = list(continuous = wrap("cor", method = "pearson", use = "complete.obs", size = 4)),
    diag  = list(continuous = wrap("densityDiag", alpha = 0.6)),
    progress = FALSE
  ) + theme_classic(base_size = 12)
  ggsave(paste0(out_prefix, "_pairs.png"), pp, width = 10, height = 10, dpi = 200)
  cat(sprintf("wrote %s_pairs.png\n", out_prefix))
}
