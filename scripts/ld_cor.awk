#!/usr/bin/awk -f
#
# ld_cor.awk - correlation between two plink2 --r2-unphased outputs, no R needed.
#
# Streams both .vcor files, joins pairs on CHROM/POS_A/POS_B, and prints Pearson
# r, the attenuation slope, and the ratio profile across the r2 range. The first
# file is held in memory (about 60 bytes per pair, so ~60 MB per million pairs);
# the second streams past it.
#
# Pairs are keyed on positions rather than variant IDs: after a liftover with
# allele swaps the ID column can carry the pre-swap allele order while REF/ALT
# carry the post-swap one. r2 is sign-invariant, so a swapped pair still
# compares correctly.
#
# Usage:
#   awk -f ld_cor.awk SAFELD.vcor ORIGINAL.vcor
#
# The FIRST file is treated as safeld, the SECOND as the reference.

BEGIN {
    FS = "\t"
    nbin = 10
}

function key(a, b, c) { return a "_" b "_" c }

# ---- file 1: safeld ----
FNR == NR {
    if (FNR == 1) {                     # locate columns from the header
        for (i = 1; i <= NF; i++) {
            if ($i == "#CHROM_A")    cA = i
            else if ($i == "POS_A")  pA = i
            else if ($i == "POS_B")  pB = i
            else if ($i == "UNPHASED_R2") rc = i
        }
        if (!cA || !pA || !pB || !rc) {
            print "ERROR: " FILENAME " is not plink2 --r2-unphased output" > "/dev/stderr"
            print "       header was: " $0 > "/dev/stderr"
            bad = 1; exit 1
        }
        next
    }
    v = $rc + 0
    if (v < 0) v = 0; if (v > 1) v = 1
    S[key($cA, $pA, $pB)] = v
    n1++
    next
}

# ---- file 2: reference ----
{
    if (FNR == 1) {
        for (i = 1; i <= NF; i++) {
            if ($i == "#CHROM_A")    cA2 = i
            else if ($i == "POS_A")  pA2 = i
            else if ($i == "POS_B")  pB2 = i
            else if ($i == "UNPHASED_R2") rc2 = i
        }
        if (!cA2 || !pA2 || !pB2 || !rc2) {
            print "ERROR: " FILENAME " is not plink2 --r2-unphased output" > "/dev/stderr"
            bad = 1; exit 1
        }
        next
    }
    n2++
    k = key($cA2, $pA2, $pB2)
    if (!(k in S)) { only2++; next }

    y = S[k]                     # safeld
    x = $rc2 + 0                 # original
    if (x < 0) x = 0; if (x > 1) x = 1

    n++
    sx += x; sy += y; sxx += x * x; syy += y * y; sxy += x * y
    sd += (y - x) * (y - x)

    b = int(x * nbin); if (b >= nbin) b = nbin - 1
    BN[b]++; BX[b] += x; BY[b] += y

    if (x > 0.2) {
        st++
        r = (x > 0 ? y / x : 0)
        rb = int(r / 0.05); if (rb > 29) rb = 29
        RB[rb]++
        RA[st] = r
    }
}

END {
    if (bad) exit 1
    printf "\n"
    printf "  safeld pairs      : %d\n", n1 + 0
    printf "  reference pairs   : %d\n", n2 + 0
    printf "  in common         : %d\n", n + 0
    if (only2 > 0)
        printf "  reference only    : %d  (%.1f%% - if you omitted --ld-window-r2 0,\n                      the pairs safeld attenuated below the default\n                      threshold have removed themselves from this join)\n", only2, 100 * only2 / n2

    if (n < 2) { print "\n  too few shared pairs to correlate"; exit 1 }

    mx = sx / n; my = sy / n
    cov = sxy / n - mx * my
    vx  = sxx / n - mx * mx
    vy  = syy / n - my * my
    r   = (vx > 0 && vy > 0) ? cov / sqrt(vx * vy) : 0
    slope = (sxx > 0) ? sxy / sxx : 0

    printf "\n  --- safeld vs reference ---\n"
    printf "  Pearson r         : %.4f\n", r
    printf "  Pearson r2        : %.4f\n", r * r
    printf "  RMSE              : %.4f\n", sqrt(sd / n)
    printf "  mean safeld - ref : %+.4f\n", my - mx
    printf "  attenuation slope : %.4f%s\n", slope, (slope < 0.9 ? "   <- LD is being lost" : "")

    printf "\n  by reference r2 decile:\n"
    printf "    %-12s %9s %10s %10s %8s\n", "range", "n", "mean ref", "mean safeld", "ratio"
    for (b = 0; b < nbin; b++) {
        if (BN[b] < 20) continue
        mo = BX[b] / BN[b]; ms = BY[b] / BN[b]
        printf "    %.1f - %.1f    %9d %10.4f %10.4f %8.3f\n", b / nbin, (b + 1) / nbin, BN[b], mo, ms, (mo > 0 ? ms / mo : 0)
    }

    if (st > 50) {
        # median and quartiles of the per-pair ratio, on pairs with ref r2 > 0.2
        nq = asort_local(RA, st)
        printf "\n  ratio safeld/reference where reference r2 > 0.2  (n=%d):\n", st
        printf "    5%%=%.3f  25%%=%.3f  median=%.3f  75%%=%.3f  95%%=%.3f\n",
               RA[int(st * 0.05) + 1], RA[int(st * 0.25) + 1], RA[int(st * 0.50) + 1],
               RA[int(st * 0.75) + 1], RA[int(st * 0.95) + 1]

        peaks = 0
        for (b = 0; b < 30; b++) {
            lo = (b > 0 ? RB[b - 1] : 0); hi = RB[b + 1]
            if (RB[b] > lo && RB[b] > hi && RB[b] > st * 0.03) peaks++
        }
        printf "    distinct modes: %d%s\n", peaks,
               (peaks > 1 ? "   <- separated bands, not smooth noise" : "   <- single mode")
    }
    printf "\n"
}

# insertion-free simple sort (gawk asort is not portable to mawk/BSD awk)
function asort_local(A, m,   i, j, t, gap) {
    gap = int(m / 2)
    while (gap > 0) {
        for (i = gap + 1; i <= m; i++) {
            t = A[i]; j = i
            while (j > gap && A[j - gap] > t) { A[j] = A[j - gap]; j -= gap }
            A[j] = t
        }
        gap = int(gap / 2)
    }
    return m
}
