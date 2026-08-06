#pragma once

// Optional native reader for plink2 .pgen/.pvar/.psam triples, so safeld can
// start from a pgen without an intermediate VCF export.
//
// Why this is worth having beyond convenience: exporting a hard-call pgen with
// "--export vcf vcf-dosage=DS" writes the DS subfield only for the entries that
// happen to carry a dosage track, producing a VCF where DS covers a few percent
// of calls while GT covers all of them. Reading the pgen directly removes that
// failure mode at the source, because the format records dosage presence
// explicitly per sample rather than leaving it to be inferred from an absent
// subfield.
//
// Built only when SAFELD_PGEN=ON. pgenlib is LGPL-3.0; safeld remains MIT and
// links against it. Keep it a shared library, or ship the object files, so the
// LGPL relinking requirement is satisfied.

#include <cstdint>
#include <string>
#include <vector>

#include "vcf_processor.h"  // Variant, DosageStats

struct PgenVariantRecord {
    std::string chrom;
    std::string id;
    int pos = 0;
    std::string ref;
    std::string alt;
    bool biallelic = true;
};

// Summary of one .pgen, mirroring InputScan so preprocess can report the same
// things it reports for a VCF.
struct PgenScan {
    bool ok = false;
    int n_samples = 0;
    long long variants = 0;
    long long multiallelic = 0;
    // Fraction of calls carrying an explicit dosage, averaged over variants.
    // Unlike a VCF this is exact and free: the format stores it.
    double dosage_presence = 0.0;
    double hardcall_rate = 0.0;
};

class PgenSource {
public:
    PgenSource();
    ~PgenSource();

    PgenSource(const PgenSource&) = delete;
    PgenSource& operator=(const PgenSource&) = delete;

    // pgen_path is the .pgen; .pvar and .psam are derived from the same stem
    // unless given explicitly. Returns false with a message in `error`.
    bool open(const std::string& pgen_path, std::string& error,
              const std::string& pvar_path = "", const std::string& psam_path = "");

    int nSamples() const { return n_samples_; }
    long long nVariants() const { return static_cast<long long>(variants_.size()); }
    const std::vector<std::string>& sampleIds() const { return sample_ids_; }
    const PgenVariantRecord& variant(long long vidx) const { return variants_[vidx]; }

    // Scans up to max_variants records for the same report preprocess prints
    // for a VCF. Cheap: it reads dosage presence counts only.
    PgenScan scan(long long max_variants = 5000);

    // Fills `dosages` with one ALT dosage per selected sample on the diploid
    // 0..2 scale, using MISSING_DOSAGE for samples with neither a dosage nor a
    // hard call, and reports observed totals in `stats`. `prefer_dosage`
    // selects DS-with-hardcall-fill; when false only hard calls are used.
    bool readVariant(long long vidx, bool prefer_dosage,
                     std::vector<double>& dosages, DosageStats& stats,
                     long long* filled_from_hardcall = nullptr);

    // Restricts subsequent reads to these sample indices, in this order.
    void setSampleSubset(const std::vector<int>& indices);

private:
    struct Impl;
    Impl* impl_;

    int n_samples_ = 0;
    std::vector<std::string> sample_ids_;
    std::vector<PgenVariantRecord> variants_;
    std::vector<int> sample_indices_;
    // Per-variant scratch: sample index -> position in dosage_main, or UINT32_MAX.
    std::vector<uint32_t> dosage_rank_;

    bool loadPvar(const std::string& path, std::string& error);
    bool loadPsam(const std::string& path, std::string& error);
};
