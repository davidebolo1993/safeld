#!/usr/bin/env python3
"""Generate the safeld test fixtures.

Writes the same genotypes twice, as a VCF and as a plink1 .bed/.bim/.fam, so
the VCF reader and the pgenlib reader can be required to produce byte-identical
standardized matrices. Any divergence between them is then a reader bug rather
than a difference in the data.

Also writes the awkward inputs that caused real failures:
  dirty.vcf       ID-less records, duplicate loci and IDs, a multiallelic site,
                  an all-missing variant, a monomorphic variant, a half call
  sparse_ds.vcf   GT complete but DS written for only a fraction of samples,
                  which is what plink2 emits when exporting a hard-call pgen
                  with vcf-dosage=DS

Usage: python3 make_testdata.py OUTDIR
"""
import os
import random
import struct
import sys


def write_vcf(path, samples, variants, with_ds=None):
    """variants: list of (chrom, pos, vid, ref, alt, [alt_count or None]).

    with_ds: None for GT only, else a function(v_index, s_index) -> bool saying
    whether that call carries a DS subfield.
    """
    with open(path, "w") as f:
        f.write("##fileformat=VCFv4.2\n##contig=<ID=1>\n")
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        if with_ds is not None:
            f.write('##FORMAT=<ID=DS,Number=A,Type=Float,Description="Dosage">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t")
        f.write("\t".join(samples) + "\n")
        for vi, (chrom, pos, vid, ref, alt, calls) in enumerate(variants):
            fmt = "GT:DS" if with_ds is not None else "GT"
            cells = []
            for si, g in enumerate(calls):
                if g is None:
                    cells.append("./." if with_ds is None else "./.")
                    continue
                gt = "0/0" if g == 0 else "0/1" if g == 1 else "1/1"
                if with_ds is not None and with_ds(vi, si):
                    cells.append(f"{gt}:{float(g):.3f}")
                elif with_ds is not None:
                    cells.append(gt)          # DS subfield simply omitted
                else:
                    cells.append(gt)
            f.write(f"{chrom}\t{pos}\t{vid}\t{ref}\t{alt}\t99\tPASS\t.\t{fmt}\t")
            f.write("\t".join(cells) + "\n")


def write_plink1(stem, samples, variants):
    """.bed is variant-major 2-bit: 00 hom A1, 01 missing, 10 het, 11 hom A2.

    .bim column 5 is A1 and column 6 is A2, and plink2 reads A1 as ALT, so a
    count of 2 ALT alleles is hom A1 -> 00.
    """
    with open(stem + ".fam", "w") as f:
        for s in samples:
            f.write(f"{s}\t{s}\t0\t0\t0\t-9\n")

    with open(stem + ".bim", "w") as f:
        for chrom, pos, vid, ref, alt, _ in variants:
            f.write(f"{chrom}\t{vid}\t0\t{pos}\t{alt}\t{ref}\n")

    code = {2: 0b00, None: 0b01, 1: 0b10, 0: 0b11}
    with open(stem + ".bed", "wb") as f:
        f.write(struct.pack("BBB", 0x6C, 0x1B, 0x01))   # magic + variant-major
        for _, _, _, _, _, calls in variants:
            buf = bytearray((len(samples) + 3) // 4)
            for si, g in enumerate(calls):
                buf[si // 4] |= code[g] << (2 * (si % 4))
            f.write(bytes(buf))


def main():
    outdir = sys.argv[1] if len(sys.argv) > 1 else "testdata"
    os.makedirs(outdir, exist_ok=True)
    random.seed(20260807)

    # ---- matched pair: identical genotypes as VCF and as bed/bim/fam --------
    n_samples, n_variants = 150, 60
    samples = [f"S{i:03d}" for i in range(n_samples)]
    variants = []
    for v in range(n_variants):
        af = random.uniform(0.12, 0.5)
        calls = [sum(1 for _ in range(2) if random.random() < af) for _ in range(n_samples)]
        # A few genuinely missing calls, which both readers must treat alike.
        if v % 17 == 0:
            calls[3] = None
            calls[11] = None
        pos = (v + 1) * 1000
        variants.append(("1", pos, f"1:{pos}:A:G", "A", "G", calls))

    write_vcf(os.path.join(outdir, "matched.vcf"), samples, variants)
    write_plink1(os.path.join(outdir, "matched"), samples, variants)

    with open(os.path.join(outdir, "extract.txt"), "w") as f:
        for _, pos, vid, _, _, _ in variants[::3]:
            f.write(vid + "\n")
    with open(os.path.join(outdir, "samples.txt"), "w") as f:
        for s in samples[:60]:
            f.write(s + "\n")

    # ---- sparse DS: GT complete, DS on a fraction (the plink2 export shape) --
    def ds_present(vi, si):
        return random.random() < [1.0, 0.5, 0.05][vi % 3]

    random.seed(7)
    write_vcf(os.path.join(outdir, "sparse_ds.vcf"), samples, variants, with_ds=ds_present)

    # ---- dirty VCF: every pathology that caused a real failure --------------
    ds = [f"S{i:02d}" for i in range(10)]

    def g(*vals):
        return list(vals)

    dirty = [
        ("1", 100, ".", "A", "G", g(1, 0, 0, 1, 0, 0, 0, 0, 0, 0)),      # ID-less, keep
        ("1", 200, ".", "C", "T", g(1, 1, 0, 0, 0, 0, 0, 0, 0, 0)),      # ID-less, keep
        ("1", 300, "rs1", "A", "T", g(1, 0, 1, 0, 0, 0, 0, 0, 0, 0)),    # dup ID, drop both
        ("1", 400, "rs1", "A", "C", g(1, 0, 0, 0, 1, 0, 0, 0, 0, 0)),
        ("1", 500, ".", "A", "G", g(None,) * 10),                        # all missing, drop
        ("1", 600, ".", "A", "G", g(0, 0, 0, 0, 0, 0, 0, 0, 0, 0)),      # monomorphic, drop
        ("1", 700, ".", "A", "G", g(1, 1, 0, 0, 0, 0, 0, 0, 0, 0)),      # dup locus, drop both
        ("1", 700, ".", "A", "G", g(1, 0, 0, 0, 0, 0, 0, 0, 0, 0)),
        ("1", 800, ".", "A", "G", g(1, 1, 1, 0, 0, 0, 0, 0, 0, 0)),      # keep
    ]
    write_vcf(os.path.join(outdir, "dirty.vcf"), ds, dirty)
    # A multiallelic record cannot be produced by the writer above, so append it.
    with open(os.path.join(outdir, "dirty.vcf"), "a") as f:
        cells = "\t".join(["0/1", "0/2"] + ["0/0"] * 8)
        f.write(f"1\t900\t.\tA\tG,T\t99\tPASS\t.\tGT\t{cells}\n")

    print(f"wrote fixtures to {outdir}/")
    print(f"  matched.vcf + matched.bed/.bim/.fam  {n_variants} variants x {n_samples} samples")
    print(f"  sparse_ds.vcf                        DS on 100%/50%/5% of samples by variant")
    print(f"  dirty.vcf                            10 records, 8 pathologies")
    print(f"  extract.txt, samples.txt             subset lists")


if __name__ == "__main__":
    main()
