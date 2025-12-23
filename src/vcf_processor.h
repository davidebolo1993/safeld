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

class VCFProcessor {
private:
    std::string vcf_file_;
    double maf_filter_;
    htsFile* vcf_fp_;
    bcf_hdr_t* hdr_;
    bcf1_t* rec_;

    std::vector<std::string> target_samples_;
    std::vector<int> sample_indices_;

    int total_variants_;
    int filtered_variants_;
    int duplicate_variants_;

    bool openVCF();
    void closeVCF();
    bool parseHeader();
    void setupTargetSamples(const std::string& sample_list_str);

    double extractAlleleFrequency(bcf1_t* rec);
    bool extractDosages(bcf1_t* rec, std::vector<double>& dosages);

public:
    VCFProcessor(const std::string& vcf_file, double maf_filter);
    ~VCFProcessor();

    bool initialize(const std::string& sample_list = "");

    // OLD: Loads entire VCF into memory (kept for backward compatibility)
    std::vector<std::unique_ptr<Variant>> processVariants();

    // NEW: Streaming callback-based processing
    using VariantCallback = std::function<void(std::unique_ptr<Variant>)>;
    void streamVariants(VariantCallback callback);

    const std::vector<std::string>& getTargetSamples() const { return target_samples_; }
    int getTotalVariants() const { return total_variants_; }
    int getFilteredVariants() const { return filtered_variants_; }
    int getDuplicateVariants() const { return duplicate_variants_; }
};
