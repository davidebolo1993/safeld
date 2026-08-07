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
#include <unordered_set>
#include <vector>

#include "genotype_source.h"
#include "vcf_processor.h"  // Variant, DosageStats, DosageField

struct PgenVariantRecord {
    // Index of this record in the .pgen, which is not its position in the
    // vector once an extract list has been applied.
    long long vidx = 0;
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

class PgenSource : public GenotypeSource {
public:
    // metadata_style picks which sidecar pair to read: Pvar for .pvar/.psam
    // (plink 2), Bim for .bim/.fam (plink 1). A .bed is read through the same
    // pgenlib entry points as a .pgen; it simply never carries dosages.
    enum class Metadata { Pvar, Bim };

    PgenSource(const std::string& genotype_path, Metadata metadata,
               double maf_filter, double max_missing_rate, DosageField dosage_field);
    ~PgenSource() override;

    PgenSource(const PgenSource&) = delete;
    PgenSource& operator=(const PgenSource&) = delete;

    // GenotypeSource
    bool initialize(const std::string& sample_list) override;
    void streamVariants(VariantCallback callback) override;
    const std::vector<std::string>& getTargetSamples() const override { return target_samples_; }
    std::vector<std::string> getContigNames() const override;
    const SourceCounts& counts() const override { return counts_; }
    void setExtractIds(std::vector<std::string> ids) override;
    std::string describe() const override {
        return metadata_ == Metadata::Bim ? "bed/bim/fam" : "pgen/pvar/psam";
    }

    // Explicit sidecar paths; otherwise derived from the genotype file stem.
    void setMetadataPaths(const std::string& variants, const std::string& samples) {
        variants_path_ = variants;
        samples_path_ = samples;
    }

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

private:
    std::string genotype_path_;
    Metadata metadata_;
    double maf_filter_;
    double max_missing_rate_;
    DosageField dosage_field_;
    DosageField effective_field_ = DosageField::Auto;
    std::string variants_path_;
    std::string samples_path_;
    std::vector<std::string> target_samples_;
    SourceCounts counts_;
    std::vector<std::string> extract_ids_;

    bool openFile(std::string& error);
    // Non-empty once setExtractIds has run; consulted while loading the table.
    std::unordered_set<std::string> extract_set_;
    void chooseField();

    struct Impl;
    Impl* impl_;

    int n_samples_ = 0;
    // Records in the file, which exceeds variants_.size() once an extract list
    // has been applied. pgenlib must be initialised with the true count.
    long long total_in_file_ = 0;
    std::vector<std::string> sample_ids_;
    std::vector<PgenVariantRecord> variants_;
    std::vector<int> sample_indices_;
    // Per-variant scratch: sample index -> position in dosage_main, or UINT32_MAX.
    std::vector<uint32_t> dosage_rank_;

    bool loadPvar(const std::string& path, std::string& error);
    bool loadPsam(const std::string& path, std::string& error);
    bool loadBim(const std::string& path, std::string& error);
    bool loadFam(const std::string& path, std::string& error);
};
