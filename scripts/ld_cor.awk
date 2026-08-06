#!/usr/bin/awk -f
#
# ld_cor.awk - Pearson correlation between two plink2 --r2-unphased outputs.
#
# Two ways to run it:
#
#   awk -f ld_cor.awk ORIGINAL.vcor SAFELD.vcor
#       Reads both files and joins pairs on CHROM/POS_A/POS_B. Safe even when
#       the two runs kept different variants. The first file is held in memory,
#       about 60 bytes per pair.
#
#   paste ORIGINAL.vcor SAFELD.vcor | awk -f ld_cor.awk -v paste_mode=1
#       Streams, constant memory, but assumes both files list exactly the same
#       pairs in exactly the same order. That assumption is checked here rather
#       than trusted: paste silently pairs row 1 with row 1, so if one run
#       dropped a variant every subsequent row is compared against the wrong
#       pair and the correlation is meaningless. This mode aborts if the
#       positions on a row disagree.
#
# Pairs are keyed on positions rather than variant IDs: after a liftover with
# allele swaps the ID column can carry the pre-swap allele order while REF/ALT
# carry the post-swap one. r2 is sign-invariant, so a swapped pair still
# compares correctly.

BEGIN { FS = "\t" }

function hdr(prefix,    i) {
    for (i = 1; i <= NF; i++) {
        if ($i == "#CHROM_A")         C[prefix] = i
        else if ($i == "POS_A")       A[prefix] = i
        else if ($i == "POS_B")       B[prefix] = i
        else if ($i == "UNPHASED_R2") R[prefix] = i
    }
    if (!C[prefix] || !A[prefix] || !B[prefix] || !R[prefix]) {
        print "ERROR: not plink2 --r2-unphased output: " FILENAME > "/dev/stderr"
        bad = 1
        exit 1
    }
}

function clamp(v) { v = v + 0; return (v < 0 ? 0 : (v > 1 ? 1 : v)) }

function accumulate(x, y) {
    n++
    sx += x; sy += y; sxx += x * x; syy += y * y; sxy += x * y
    sdd += (y - x) * (y - x)
    b = int(x * 10); if (b > 9) b = 9
    BN[b]++; BX[b] += x; BY[b] += y
}

# ---------------------------------------------------------------- paste mode
paste_mode && FNR == 1 {
    # Header of the pasted stream: the left file's columns, then the right's.
    nleft = 0
    for (i = 1; i <= NF; i++) if ($i == "#CHROM_A") { if (!first_seen) first_seen = i; else nleft = i - 1 }
    if (!nleft) {
        print "ERROR: could not find two #CHROM_A columns in the pasted header." > "/dev/stderr"
        print "       Are both inputs plink2 --r2-unphased files?" > "/dev/stderr"
        bad = 1; exit 1
    }
    for (i = 1; i <= nleft; i++) {
        if ($i == "#CHROM_A") cL = i; else if ($i == "POS_A") aL = i
        else if ($i == "POS_B") bL = i; else if ($i == "UNPHASED_R2") rL = i
    }
    for (i = nleft + 1; i <= NF; i++) {
        if ($i == "#CHROM_A") cR = i; else if ($i == "POS_A") aR = i
        else if ($i == "POS_B") bR = i; else if ($i == "UNPHASED_R2") rR = i
    }
    next
}
paste_mode {
    # The check that makes this mode safe.
    if ($cL != $cR || $aL != $aR || $bL != $bR) {
        printf "ERROR: rows are misaligned at line %d.\n", FNR > "/dev/stderr"
        printf "       left  = %s:%s-%s\n", $cL, $aL, $bL > "/dev/stderr"
        printf "       right = %s:%s-%s\n", $cR, $aR, $bR > "/dev/stderr"
        print  "       The two .vcor files do not contain the same pairs in the" > "/dev/stderr"
        print  "       same order, so paste has lined up unrelated rows. Re-run" > "/dev/stderr"
        print  "       without paste_mode to join on position instead." > "/dev/stderr"
        bad = 1
        exit 1
    }
    accumulate(clamp($rL), clamp($rR))
    next
}

# ------------------------------------------------------------------ join mode
FNR == 1 { hdr(FILENAME == ARGV[1] ? "L" : "R"); next }

FNR == NR {                                   # first file: the reference
    R2[$C["L"] "_" $A["L"] "_" $B["L"]] = clamp($R["L"])
    n1++
    next
}
{                                             # second file: safeld
    n2++
    k = $C["R"] "_" $A["R"] "_" $B["R"]
    if (!(k in R2)) { only2++; next }
    accumulate(R2[k], clamp($R["R"]))
}

END {
    if (bad) exit 1
    if (n < 2) { print "no shared pairs to correlate" > "/dev/stderr"; exit 1 }

    if (!paste_mode) {
        printf "\n  reference pairs : %d\n", n1 + 0
        printf "  safeld pairs    : %d\n", n2 + 0
        printf "  in common       : %d\n", n
        if (only2 > 0)
            printf "  safeld only     : %d\n", only2
    } else {
        printf "\n  pairs           : %d  (row alignment verified)\n", n
    }

    denom = sqrt((n * sxx - sx * sx) * (n * syy - sy * sy))
    if (denom == 0) {
        print "\n  Pearson r: undefined (zero variance)"
        exit 0
    }
    r = (n * sxy - sx * sy) / denom
    slope = (sxx > 0) ? sxy / sxx : 0

    printf "\n  Pearson r       : %.6f\n", r
    printf "  Pearson r2      : %.6f\n", r * r
    printf "  RMSE            : %.6f\n", sqrt(sdd / n)
    printf "  mean difference : %+.6f  (safeld - reference)\n", (sy - sx) / n
    printf "  slope (y~0+x)   : %.4f%s\n", slope, (slope < 0.9 ? "   <- LD is being lost" : "")

    printf "\n  by reference r2 decile:\n"
    printf "    %-12s %9s %10s %10s %8s\n", "range", "n", "mean ref", "mean safeld", "ratio"
    for (b = 0; b < 10; b++) {
        if (BN[b] < 20) continue
        mo = BX[b] / BN[b]; ms = BY[b] / BN[b]
        printf "    %.1f - %.1f    %9d %10.4f %10.4f %8.3f\n",
               b / 10, (b + 1) / 10, BN[b], mo, ms, (mo > 0 ? ms / mo : 0)
    }
    printf "\n"
}
