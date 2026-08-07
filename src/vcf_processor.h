#pragma once

#include <string>
#include <vector>
#include <memory>
#include <functional>
#include <htslib/vcf.h>
#include <htslib/hts.h>
#include <unordered_map>
#include <unordered_set>

#include "genotype_source.h"

struct Variant {
    std::string id;
    std::string chrom;
    int pos;
    std::string ref;
    std::string alt;
    double af;
    std::vector<double> dosages;
};

// Which FORMAT field supplies genotype dosages.
//   Auto - prefer DS, but fall back to GT for records where DS is absent or too
//          sparse. Files exported by plink2 routinely omit the DS subfield for
//          many samples while GT stays complete; mean-imputing those absences
//          silently attenuates LD, so using the complete field is better.
//   DS   - dosages only, never fall back.
//   GT   - hard calls only, ignore DS entirely.
enum class DosageField { Auto, DS, GT };

// What an up-front scan of the first records found. Reported at the start of
// preprocessing and used to choose the dosage source in Auto mode, so the tool
// explains its own input rather than needing a companion script.
struct InputScan {
    bool ok = false;
    int n_samples = 0;
    long long records = 0;
    long long biallelic = 0;
    long long multiallelic = 0;
    long long with_info_af = 0;
    long long no_id = 0;

    bool ds_declared = false;
    bool gt_declared = false;
    // Mean fraction of the selected samples carrying a usable value, over the
    // scanned biallelic records.
    double ds_presence = 0.0;
    double gt_call_rate = 0.0;
    long long ds_absent_records = 0;
};

// Summary of the observed (non-missing) entries of one variant's dosage vector.
// Collected while the missing sentinels are still in place so that allele
// frequency and missingness are derived from real calls only.
struct DosageStats {
    double observed_sum = 0.0;
    int observed_count = 0;
};

class VCFProcessor : public GenotypeSource {
private:
    std::string vcf_file_;
    double maf_filter_;
    double max_missing_rate_;
    DosageField dosage_field_;
    // What the run will actually read, after Auto has consulted the scan.
    DosageField effective_field_;
    InputScan scan_;
    std::string temp_dir_;
    htsFile* vcf_fp_;
    bcf_hdr_t* hdr_;
    bcf1_t* rec_;

    std::vector<std::string> target_samples_;
    std::vector<int> sample_indices_;
    // Whether to trust INFO/AF instead of computing the frequency from the
    // genotypes. Off by default: the field is only correct if it describes
    // exactly the samples present, and tools that subset samples routinely
    // recompute AC and AN while leaving AF stale.
    bool use_info_af_requested_;
    bool use_info_af_;

    int total_variants_;
    int filtered_variants_;
    int duplicate_variants_;
    int multiallelic_variants_;
    int missing_filtered_variants_;
    int gt_fallback_variants_;
    long long gt_filled_calls_;
    SourceCounts counts_;
    long long not_extracted_ = 0;
    long long maf_filtered_ = 0;
    std::unordered_set<std::string> extract_ids_;
    std::unordered_set<std::string> extract_seen_;

    bool openVCF();
    void closeVCF();
    bool parseHeader();
    void setupTargetSamples(const std::string& sample_list_str);

    bool extractDosages(bcf1_t* rec, std::vector<double>& dosages, DosageStats& stats);
    // Chooses effective_field_ from the scan and reports the decision.
    void reportScanAndChooseField();
    bool extractDosagesFromGT(bcf1_t* rec, std::vector<double>& dosages, DosageStats& stats);
    std::string makeSpoolPath() const;

public:
    // max_missing_rate: variants whose fraction of missing calls exceeds this are
    // dropped. temp_dir: where the deduplication spool is written (empty = $TMPDIR).
    VCFProcessor(const std::string& vcf_file, double maf_filter,
                 double max_missing_rate = 0.1, const std::string& temp_dir = "",
                 DosageField dosage_field = DosageField::Auto,
                 bool use_info_af = false);
    ~VCFProcessor();

    bool initialize(const std::string& sample_list = "") override;
    void streamVariants(VariantCallback callback) override;
    std::vector<std::string> getContigNames() const override;
    const SourceCounts& counts() const override { return counts_; }
    void setExtractIds(std::vector<std::string> ids) override;
    std::string describe() const override { return "VCF"; }

    // Reads up to max_records from a second handle on the same file and
    // summarises it. Called by initialize(); exposed for callers that want the
    // numbers without streaming.
    InputScan scanInput(int max_records = 5000);
    const InputScan& inputScan() const { return scan_; }

    const std::vector<std::string>& getTargetSamples() const override { return target_samples_; }
    int getTotalVariants() const { return total_variants_; }
    int getFilteredVariants() const { return filtered_variants_; }
    int getDuplicateVariants() const { return duplicate_variants_; }
    int getMultiallelicVariants() const { return multiallelic_variants_; }
    int getMissingFilteredVariants() const { return missing_filtered_variants_; }
    int getGtFallbackVariants() const { return gt_fallback_variants_; }
    long long getGtFilledCalls() const { return gt_filled_calls_; }
};
