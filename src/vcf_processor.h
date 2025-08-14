#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <memory>
#include <htslib/vcf.h>
#include <htslib/hts.h>

struct VariantRecord {
    std::string id;
    std::string chrom;
    int pos;
    std::string ref;
    std::string alt;
    double af;
    std::vector<double> dosages;
    
    VariantRecord() = default;
    VariantRecord(const VariantRecord&) = delete;
    VariantRecord& operator=(const VariantRecord&) = delete;
    VariantRecord(VariantRecord&&) = default;
    VariantRecord& operator=(VariantRecord&&) = default;
};

class VCFProcessor {
private:
    std::string vcf_file_;
    std::vector<std::string> target_samples_;
    std::vector<int> sample_indices_;
    double maf_filter_;
    
    // HTSlib objects
    htsFile* vcf_fp_;
    bcf_hdr_t* hdr_;
    bcf1_t* rec_;
    
    // Processing stats
    size_t total_variants_;
    size_t filtered_variants_;
    size_t duplicate_variants_;
    
    // Internal methods
    bool openVCF();
    void closeVCF();
    bool parseHeader();
    bool extractDosages(bcf1_t* rec, std::vector<double>& dosages);
    double extractAlleleFrequency(bcf1_t* rec);
    void setupTargetSamples(const std::string& sample_list_str);
    
public:
    explicit VCFProcessor(const std::string& vcf_file, double maf_filter = 0.01);
    ~VCFProcessor();
    
    bool initialize(const std::string& sample_list_str = "");
    std::vector<std::unique_ptr<VariantRecord>> processVariants();
    
    // Getters
    const std::vector<std::string>& getTargetSamples() const { return target_samples_; }
    size_t getTotalVariants() const { return total_variants_; }
    size_t getFilteredVariants() const { return filtered_variants_; }
    size_t getDuplicateVariants() const { return duplicate_variants_; }
};

