#!/usr/bin/env bash
#
# safeld_check.sh - single-pass report on everything that makes safeld's output
# differ from the original R implementation, or from the true LD.
#
# WHAT IT PRINTS
#   Aggregate numbers only: counts, means, histograms over variants. No sample
#   identifiers, no per-sample genotypes, no per-variant lines. Roughly 70 lines,
#   safe to copy off a secured host by hand.
#
# WHY EACH SECTION EXISTS
#   [2] variant properties  - things the R code and safeld handle differently
#   [3] FIELD COMPLETENESS  - the one that matters most. A VCF can carry GT for
#                             every sample but omit the DS subfield for many of
#                             them ("0|1:0.97" next to a bare "0|0"). Checking
#                             GT alone reports 0% missing and hides it entirely.
#                             Mean-imputing absent dosages attenuates pairwise r2
#                             in proportion to the DS call rate, which appears as
#                             discrete lower bands against the original LD.
#   [4] DS vs GT agreement  - whether hard calls are a usable substitute
#   [5] predicted attenuation
#   [6] WHAT SAFELD WILL DO - variants surviving under -dosage-field auto/DS/GT
#
# REQUIREMENTS: bcftools and awk. Nothing else.
#
# USAGE
#   ./safeld_check.sh INPUT.vcf[.gz] [N_RECORDS] [REGION] [MAF] [MAX_MISSING]
#
#   N_RECORDS    variants to sample (default 20000; 0 = all). Duplicate
#                detection only sees the sampled window.
#   REGION       optional bcftools region, e.g. 1:113000000-114000000
#   MAF          MAF threshold to simulate (default 0.01, safeld's -maf)
#   MAX_MISSING  missingness threshold  (default 0.1,  safeld's -max-missing)
#
set -uo pipefail

VCF=${1:?usage: safeld_check.sh INPUT.vcf[.gz] [N_RECORDS] [REGION] [MAF] [MAX_MISSING]}
NREC=${2:-20000}
REGION=${3:-}
MAF=${4:-0.01}
MAXMISS=${5:-0.1}

command -v bcftools >/dev/null || { echo "bcftools not found in PATH" >&2; exit 1; }

TMPQ=$(mktemp) || exit 1
trap 'rm -f "$TMPQ"' EXIT

echo "=================================================================="
echo " safeld input check"
echo " file    : $(basename "$VCF")"
echo " sampled : ${NREC} records${REGION:+   region: $REGION}"
echo " filters : maf=$MAF  max-missing=$MAXMISS"
echo " date    : $(date -u '+%Y-%m-%dT%H:%M:%SZ')"
echo "=================================================================="

# ------------------------------------------------------------------ [1] -----
echo
echo "[1] HEADER"
NSAMP=$(bcftools query -l "$VCF" 2>/dev/null | wc -l | tr -d ' ')
echo "    samples in file          : $NSAMP"
SRC=$(bcftools view -h "$VCF" 2>/dev/null | sed -n 's/^##source=//p' | head -1)
echo "    ##source                 : ${SRC:-<none>}"
for tag in GT DS GP HDS; do
  line=$(bcftools view -h "$VCF" 2>/dev/null | grep -m1 "FORMAT=<ID=$tag," || true)
  if [ -n "$line" ]; then
    num=$(echo "$line" | sed -n 's/.*Number=\([^,]*\).*/\1/p')
    echo "    FORMAT $tag declared      : yes (Number=$num)"
  fi
done
INFO_IDS=$(bcftools view -h "$VCF" 2>/dev/null | sed -n 's/^##INFO=<ID=\([^,]*\).*/\1/p' | tr '\n' ' ')
echo "    INFO fields declared     : ${INFO_IDS:-<none>}"
case " $INFO_IDS " in
  *" R2 "*|*" INFO "*|*" ER2 "*) ;;
  *) echo "    note: no imputation-quality field (R2/INFO/ER2). The chr22 reference"
     echo "          test used data pre-filtered at R2>=0.8; this file was not." ;;
esac

# ------------------------------------------------------------------ query ---
# %DS prints "." when the subfield is absent for a sample, which is exactly the
# signal section [3] needs. No arrays here: `set -u` with an empty array is an
# error on bash 3.x and would silently produce an all-zero report.
# Only ask bcftools for fields the header actually declares; querying an absent
# tag is a hard error. Missing tags are filled with a literal "." so the awk sees
# a uniform "GT=DS" shape either way.
HAS_GT=$(bcftools view -h "$VCF" 2>/dev/null | grep -c 'FORMAT=<ID=GT,' || true)
HAS_DS=$(bcftools view -h "$VCF" 2>/dev/null | grep -c 'FORMAT=<ID=DS,' || true)
if [ "$HAS_GT" -gt 0 ] && [ "$HAS_DS" -gt 0 ]; then SFMT='[\t%GT=%DS]'
elif [ "$HAS_GT" -gt 0 ];                       then SFMT='[\t%GT=.]'
elif [ "$HAS_DS" -gt 0 ];                       then SFMT='[\t.=%DS]'
else
  echo >&2
  echo "  !! the header declares neither FORMAT/GT nor FORMAT/DS." >&2
  echo "     safeld has no genotype source in this file." >&2
  exit 1
fi

STREAM() {
  if [ -n "$REGION" ]; then
    bcftools query -r "$REGION" -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO$SFMT\n" "$VCF" 2>/dev/null
  else
    bcftools query -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO$SFMT\n" "$VCF" 2>/dev/null
  fi
}

if [ "$NREC" = "0" ]; then
  STREAM > "$TMPQ"
else
  STREAM 2>/dev/null | head -"$NREC" > "$TMPQ"
fi

if [ ! -s "$TMPQ" ]; then
  echo >&2
  echo "  !! bcftools produced no output. Check the file and region." >&2
  echo "     Aborting rather than reporting zeros." >&2
  exit 1
fi

awk -v MAF="$MAF" -v MAXMISS="$MAXMISS" -v NSAMP="$NSAMP" '
function fbin(x,   i) {                      # fraction -> 8 bins
  if (x >= 0.999) return 8
  if (x >= 0.99)  return 7
  if (x >= 0.95)  return 6
  if (x >= 0.90)  return 5
  if (x >= 0.75)  return 4
  if (x >= 0.50)  return 3
  if (x >= 0.25)  return 2
  return 1
}
function bar(pc,   s, k, nb) { s = ""; nb = int(pc / 2 + 0.5); for (k = 0; k < nb; k++) s = s "#"; return s }
BEGIN { FS = "\t" }
# pass 1: count locus/ID occurrences so duplicates can be dropped the way safeld
# drops them, i.e. every copy including the first.
NR == FNR {
  if (index($5, ",")) next
  CNT_LOC[$1 ":" $2 ":" $4 ":" $5]++
  if ($3 != "." && $3 != "") CNT_ID[$3]++
  next
}
{
  nrec++
  chrom = $1; pos = $2; id = $3; ref = $4; alt = $5; info = $6

  # ---- variant-level -------------------------------------------------------
  if (index(alt, ",")) { multi++; next }          # safeld skips non-biallelic
  biallelic++

  if (id == "." || id == "") noid++
  locus = chrom ":" pos ":" ref ":" alt
  has_real_id = (id != "." && id != "")
  is_dup = (CNT_LOC[locus] > 1) || (has_real_id && CNT_ID[id] > 1)
  if (has_real_id && CNT_ID[id] > 1) dupid++
  if (CNT_LOC[locus] > 1) duploc++

  has_af = (info ~ /(^|;)AF=/)
  if (!has_af) noaf++
  else {
    afs = info; sub(/^.*(^|;)AF=/, "", afs); sub(/;.*$/, "", afs)
    if (index(afs, ",")) multiaf++
  }

  # ---- per-sample field scan ----------------------------------------------
  nds = 0; ngt = 0; nboth = 0; na = 0
  sd = 0; sd2 = 0; sg = 0; sg2 = 0; sgd = 0; sa = 0; sa2 = 0
  for (i = 7; i <= NF; i++) {
    p = index($i, "=")
    if (p == 0) continue
    gt = substr($i, 1, p - 1)
    ds = substr($i, p + 1)

    # --- DS present and usable?
    ds_ok = 0; d = 0
    if (ds != "." && ds != "" && !index(ds, ",")) {
      d = ds + 0
      if (d >= 0 && d == d) { ds_ok = 1 }          # d==d rejects nan
      else neg_ds++
    }
    # --- GT present and fully called?
    gt_ok = 0; h = 0
    if (gt != "." && gt != "") {
      a = gt; gsub(/\|/, "/", a)
      m = split(a, AL, "/"); bad = 0; nalt = 0
      for (k = 1; k <= m; k++) {
        if (AL[k] == ".") { bad = 1; break }
        if (AL[k] + 0 > 0) nalt++
      }
      if (!bad && m > 0) { gt_ok = 1; h = 2.0 * nalt / m }
      else if (bad) half_or_missing_gt++
    }

    if (ds_ok) { nds++; sd += d; sd2 += d * d }
    if (gt_ok) { ngt++; sg += h; sg2 += h * h }
    if (ds_ok && gt_ok) { nboth++; sgd += d * h }
    # merged vector: DS when written, else the hard call
    if (ds_ok)      { na++; sa += d; sa2 += d * d }
    else if (gt_ok) { na++; sa += h; sa2 += h * h }
    ncell++
  }

  ns = NF - 6
  if (ns <= 0) next
  ds_rate = nds / ns
  gt_rate = ngt / ns
  DSB[fbin(ds_rate)]++
  GTB[fbin(gt_rate)]++
  sum_ds_rate += ds_rate; sum_gt_rate += gt_rate

  if (nds == 0) ds_absent_var++
  if (ds_rate >= 0.999) ds_full_var++
  if (gt_rate >= 0.999) gt_full_var++

  # ---- DS vs GT agreement (samples carrying both) -------------------------
  if (nboth >= 10) {
    md = sd / nboth; mg = sg / nboth
    vd = sd2 / nboth - md * md
    vg = sg2 / nboth - mg * mg
    cv = sgd / nboth - md * mg
    if (vd > 0 && vg > 0) {
      r2 = (cv * cv) / (vd * vg); if (r2 > 1) r2 = 1
      sum_r2 += r2; nr2++
      if (r2 >= 0.999) r2_perfect++
      R2B[fbin(r2)]++
    }
  }

  # ---- what safeld would keep, per -dosage-field mode ---------------------
  if (is_dup) { ndup_dropped++; next }            # safeld drops every copy
  keep_mode("DS", nds, ns, sd, sd2)
  keep_mode("GT", ngt, ns, sg, sg2)
  # auto: per sample, the real dosage where one was written, else the matching
  # GT hard call. Only a sample lacking BOTH is missing.
  keep_mode("AUTO", na, ns, sa, sa2)
  if (na > nds) { auto_filled_var++; auto_filled_calls += na - nds }
}
function keep_mode(tag, nobs, ns, s, s2,    mean, var, af, maf) {
  if (nobs == 0) { DROP_MISS[tag]++; return }
  if (1 - nobs / ns > MAXMISS) { DROP_MISS[tag]++; return }
  mean = s / nobs
  af = mean / 2.0
  maf = (af < 0.5 ? af : 1 - af)
  if (maf < MAF) { DROP_MAF[tag]++; return }
  var = s2 / nobs - mean * mean
  if (var <= 0) { DROP_VAR[tag]++; return }
  KEEP[tag]++
}
END {
  if (nrec == 0) { print "no records"; exit 1 }
  nb = biallelic + 0

  printf "\n[2] VARIANT PROPERTIES   (%d records sampled)\n", nrec
  printf "    multiallelic (skipped)   : %7d  %5.1f%%\n", multi + 0, 100 * multi / nrec
  printf "    biallelic (usable)       : %7d  %5.1f%%\n", nb, 100 * nb / nrec
  if (nb == 0) { print "    nothing usable"; exit 1 }
  printf "    ID is \".\"                : %7d  %5.1f%%\n", noid + 0, 100 * noid / nb
  printf "    variants w/ repeated ID  : %7d  %5.1f%%  (all copies dropped)\n", dupid + 0, 100 * dupid / nb
  printf "    variants w/ repeated locus: %6d  %5.1f%%  (all copies dropped)\n", duploc + 0, 100 * duploc / nb
  printf "    no INFO/AF               : %7d  %5.1f%%", noaf + 0, 100 * noaf / nb
  if (noaf == nb) printf "   <- R code cannot run: AF=NA drops all"
  printf "\n"
  printf "    INFO/AF multi-valued     : %7d  %5.1f%%\n", multiaf + 0, 100 * multiaf / nb

  lab[1] = "  < 25%"; lab[2] = " 25 - 50%"; lab[3] = " 50 - 75%"; lab[4] = " 75 - 90%"
  lab[5] = " 90 - 95%"; lab[6] = " 95 - 99%"; lab[7] = " 99 - 99.9%"; lab[8] = ">= 99.9%"

  printf "\n[3] FIELD COMPLETENESS   <-- the section that matters most\n"
  printf "    mean DS presence rate    : %.4f\n", sum_ds_rate / nb
  printf "    mean GT call rate        : %.4f\n", sum_gt_rate / nb
  printf "    variants with DS on 100%% of samples : %d (%.1f%%)\n", ds_full_var + 0, 100 * ds_full_var / nb
  printf "    variants with GT on 100%% of samples : %d (%.1f%%)\n", gt_full_var + 0, 100 * gt_full_var / nb
  printf "    variants with NO DS at all          : %d (%.1f%%)\n", ds_absent_var + 0, 100 * ds_absent_var / nb
  printf "    negative/NaN DS values (treated missing): %d\n", neg_ds + 0
  printf "    half-called or missing GT calls         : %d\n", half_or_missing_gt + 0
  printf "\n    per-variant DS presence rate:\n"
  for (i = 1; i <= 8; i++) { c = DSB[i] + 0; pc = 100 * c / nb; printf "      %-12s %7d  %5.1f%%  %s\n", lab[i], c, pc, bar(pc) }
  printf "\n    per-variant GT call rate:\n"
  for (i = 1; i <= 8; i++) { c = GTB[i] + 0; pc = 100 * c / nb; printf "      %-12s %7d  %5.1f%%  %s\n", lab[i], c, pc, bar(pc) }

  printf "\n[4] DS vs GT AGREEMENT   rho2 = cor(DS,GT)^2 on samples carrying both\n"
  if (nr2 == 0) printf "    not computable (too few samples carry both fields)\n"
  else {
    printf "    variants scored          : %d\n", nr2
    printf "    mean rho2                : %.4f\n", sum_r2 / nr2
    printf "    rho2 >= 0.999            : %d (%.1f%%)  <- DS is effectively the hard call\n", r2_perfect + 0, 100 * r2_perfect / nr2
    for (i = 1; i <= 8; i++) { c = R2B[i] + 0; pc = 100 * c / nr2; printf "      %-12s %7d  %5.1f%%  %s\n", lab[i], c, pc, bar(pc) }
  }

  printf "\n[5] PREDICTED LD ATTENUATION UNDER -dosage-field DS\n"
  mdr = sum_ds_rate / nb
  printf "    Mean-imputing absent dosages attenuates a pair roughly by the\n"
  printf "    product of the two call rates:\n"
  printf "      mean DS presence %.3f  ->  r2=0.80 typically reads as %.3f\n", mdr, 0.80 * mdr * mdr
  printf "      variants at 100%% DS    ->  r2=0.80 reads as 0.800\n"
  if (mdr < 0.99)
    printf "    The gap between those lines IS the banding seen against original LD.\n"
  else
    printf "    DS is essentially complete; this is not your problem.\n"

  printf "\n[6] WHAT SAFELD WILL DO   (of %d biallelic variants; %d dropped as duplicates)\n", nb, ndup_dropped + 0
  printf "    %-6s %10s %10s %10s %10s\n", "mode", "kept", "drop:miss", "drop:maf", "drop:novar"
  split("AUTO DS GT", M, " ")
  for (i = 1; i <= 3; i++) {
    t = M[i]
    printf "    %-6s %10d %10d %10d %10d\n", t, KEEP[t] + 0, DROP_MISS[t] + 0, DROP_MAF[t] + 0, DROP_VAR[t] + 0
  }
  if (auto_filled_var > 0) {
    printf "    auto filled %d absent DS call(s) across %d variant(s) from the\n", auto_filled_calls + 0, auto_filled_var + 0
    printf "    matching GT hard call (filled, not imputed).\n"
  }

  printf "\n[7] RECOMMENDATION\n"
  if (mdr < 0.99 && sum_gt_rate / nb > 0.99) {
    printf "    DS is incomplete (%.1f%% mean) while GT is complete (%.1f%%).\n", 100 * mdr, 100 * sum_gt_rate / nb
    printf "    Absent DS is NOT a missing genotype here - GT knows it. safeld\n"
    printf "    default (-dosage-field auto) now fills each gap from the matching\n"
    printf "    GT hard call instead of imputing, so no attenuation.\n"
    printf "    For a matrix from one field throughout use -dosage-field GT, and\n"
    printf "    then compare against a HARD-CALL reference LD.\n"
  } else if (mdr >= 0.99) {
    printf "    Both fields look complete. Parsing is not the problem here.\n"
    printf "    Report the -ntraits used for the run; that is the next suspect.\n"
  } else {
    printf "    Both DS and GT are incomplete. Raise -max-missing deliberately or\n"
    printf "    filter the input first; do not let variants be mean-imputed silently.\n"
  }

  printf "\n[SUMMARY - copy this line back]\n"
  printf "  nsamp=%s nrec=%d biallelic=%d multi=%d noid=%d dupid=%d duploc=%d noaf=%d", NSAMP, nrec, nb, multi+0, noid+0, dupid+0, duploc+0, noaf+0
  printf " ds_rate=%.4f gt_rate=%.4f ds_full=%.1f%% gt_full=%.1f%% ds_absent=%d", mdr, sum_gt_rate/nb, 100*ds_full_var/nb, 100*gt_full_var/nb, ds_absent_var+0
  printf " rho2=%.4f keepAUTO=%d keepDS=%d keepGT=%d\n", (nr2? sum_r2/nr2 : -1), KEEP["AUTO"]+0, KEEP["DS"]+0, KEEP["GT"]+0
}
' "$TMPQ" "$TMPQ"

echo
echo "=================================================================="
echo " Also report, since it is not in the VCF:"
echo "   * the -ntraits used for the safeld run"
echo "   * whether the reference LD was computed from GT or from DS"
echo "     (plink2 --r2-unphased reads hard calls unless you pass 'dosage'=DS)"
echo "=================================================================="
