#pragma once

#include <string>
#include <vector>
#include <memory>
#include <functional>
#include <htslib/vcf.h>
#include <htslib/hts.h>
#include <unordered_map>

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

// Summary of the observed (non-missing) entries of one variant's dosage vector.
// Collected while the missing sentinels are still in place so that allele
// frequency and missingness are derived from real calls only.
struct DosageStats {
    double observed_sum = 0.0;
    int observed_count = 0;
};

class VCFProcessor {
private:
    std::string vcf_file_;
    double maf_filter_;
    double max_missing_rate_;
    DosageField dosage_field_;
    std::string temp_dir_;
    htsFile* vcf_fp_;
    bcf_hdr_t* hdr_;
    bcf1_t* rec_;

    std::vector<std::string> target_samples_;
    std::vector<int> sample_indices_;
    // INFO/AF describes every sample in the file, so it is only a valid filter
    // when the whole cohort is being used.
    bool use_info_af_;

    int total_variants_;
    int filtered_variants_;
    int duplicate_variants_;
    int multiallelic_variants_;
    int missing_filtered_variants_;
    int gt_fallback_variants_;
    long long gt_filled_calls_;

    bool openVCF();
    void closeVCF();
    bool parseHeader();
    void setupTargetSamples(const std::string& sample_list_str);

    bool extractDosages(bcf1_t* rec, std::vector<double>& dosages, DosageStats& stats);
    bool extractDosagesFromGT(bcf1_t* rec, std::vector<double>& dosages, DosageStats& stats);
    std::string makeSpoolPath() const;

public:
    // max_missing_rate: variants whose fraction of missing calls exceeds this are
    // dropped. temp_dir: where the deduplication spool is written (empty = $TMPDIR).
    VCFProcessor(const std::string& vcf_file, double maf_filter,
                 double max_missing_rate = 0.1, const std::string& temp_dir = "",
                 DosageField dosage_field = DosageField::Auto);
    ~VCFProcessor();

    bool initialize(const std::string& sample_list = "");

    using VariantCallback = std::function<void(std::unique_ptr<Variant>)>;
    void streamVariants(VariantCallback callback);
    std::vector<std::string> getContigNames() const;

    const std::vector<std::string>& getTargetSamples() const { return target_samples_; }
    int getTotalVariants() const { return total_variants_; }
    int getFilteredVariants() const { return filtered_variants_; }
    int getDuplicateVariants() const { return duplicate_variants_; }
    int getMultiallelicVariants() const { return multiallelic_variants_; }
    int getMissingFilteredVariants() const { return missing_filtered_variants_; }
    int getGtFallbackVariants() const { return gt_fallback_variants_; }
    long long getGtFilledCalls() const { return gt_filled_calls_; }
};
