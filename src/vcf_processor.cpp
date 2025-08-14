#include "vcf_processor.h"
#include "utils.h"
#include <iostream>
#include <sstream>
#include <algorithm>
#include <cstring>
#include <cmath>          // <- ADD THIS LINE


VCFProcessor::VCFProcessor(const std::string& vcf_file, double maf_filter)
    : vcf_file_(vcf_file), maf_filter_(maf_filter), vcf_fp_(nullptr), 
      hdr_(nullptr), rec_(nullptr), total_variants_(0), 
      filtered_variants_(0), duplicate_variants_(0) {
}

VCFProcessor::~VCFProcessor() {
    closeVCF();
}

bool VCFProcessor::openVCF() {
    vcf_fp_ = hts_open(vcf_file_.c_str(), "r");
    if (!vcf_fp_) {
        logError("Failed to open VCF file: " + vcf_file_);
        return false;
    }
    
    hdr_ = bcf_hdr_read(vcf_fp_);
    if (!hdr_) {
        logError("Failed to read VCF header");
        return false;
    }
    
    rec_ = bcf_init();
    if (!rec_) {
        logError("Failed to initialize BCF record");
        return false;
    }
    
    return true;
}

void VCFProcessor::closeVCF() {
    if (rec_) {
        bcf_destroy(rec_);
        rec_ = nullptr;
    }
    
    if (hdr_) {
        bcf_hdr_destroy(hdr_);
        hdr_ = nullptr;
    }
    
    if (vcf_fp_) {
        hts_close(vcf_fp_);
        vcf_fp_ = nullptr;
    }
}

bool VCFProcessor::parseHeader() {
    int n_samples = bcf_hdr_nsamples(hdr_);
    if (n_samples == 0) {
        logError("No samples found in VCF header");
        return false;
    }
    
    logInfo("Found " + std::to_string(n_samples) + " samples in VCF");
    return true;
}

void VCFProcessor::setupTargetSamples(const std::string& sample_list_str) {
    int n_samples = bcf_hdr_nsamples(hdr_);
    
    if (sample_list_str.empty()) {
        // Use all samples
        target_samples_.reserve(n_samples);
        sample_indices_.reserve(n_samples);
        
        for (int i = 0; i < n_samples; ++i) {
            target_samples_.emplace_back(hdr_->samples[i]);
            sample_indices_.push_back(i);
        }
    } else {
        // Parse requested samples
        auto requested_samples = split(sample_list_str, ',');
        std::unordered_map<std::string, int> sample_map;
        
        for (int i = 0; i < n_samples; ++i) {
            sample_map[hdr_->samples[i]] = i;
        }
        
        for (const auto& sample : requested_samples) {
            std::string trimmed = sample;
            // Trim whitespace
            trimmed.erase(0, trimmed.find_first_not_of(" \t"));
            trimmed.erase(trimmed.find_last_not_of(" \t") + 1);
            
            auto it = sample_map.find(trimmed);
            if (it != sample_map.end()) {
                target_samples_.push_back(trimmed);
                sample_indices_.push_back(it->second);
            }
        }
        
        logInfo("Found " + std::to_string(target_samples_.size()) + " of " + 
                std::to_string(requested_samples.size()) + " requested samples");
    }
}

bool VCFProcessor::initialize(const std::string& sample_list_str) {
    Timer timer("VCF initialization");
    
    if (!openVCF()) {
        return false;
    }
    
    if (!parseHeader()) {
        return false;
    }
    
    setupTargetSamples(sample_list_str);
    
    if (target_samples_.empty()) {
        logError("No target samples found");
        return false;
    }
    
    logInfo("VCF processor initialized with " + std::to_string(target_samples_.size()) + " samples");
    return true;
}

double VCFProcessor::extractAlleleFrequency(bcf1_t* rec) {
    // Try to get AF from INFO field
    int n_values = 0;
    float* af_values = nullptr;
    
    if (bcf_get_info_float(hdr_, rec, "AF", &af_values, &n_values) > 0 && n_values > 0) {
        double af = static_cast<double>(af_values[0]);
        free(af_values);
        // Ensure consistent precision with Go (6 decimal places)
        return std::round(af * 1000000.0) / 1000000.0;
    }
    
    // If AF not available, calculate from dosages
    std::vector<double> dosages;
    if (extractDosages(rec, dosages)) {
        double sum = 0.0;
        int count = 0;
        for (double d : dosages) {
            if (d >= 0.0) {  // Valid dosage
                sum += d;
                count++;
            }
        }
        if (count > 0) {
            double af = sum / (2.0 * count);
            return std::round(af * 1000000.0) / 1000000.0;
        }
    }
    
    return -1.0;  // Unable to determine AF
}

bool VCFProcessor::extractDosages(bcf1_t* rec, std::vector<double>& dosages) {
    int n_values = 0;
    float* ds_values = nullptr;
    
    int ret = bcf_get_format_float(hdr_, rec, "DS", &ds_values, &n_values);
    if (ret <= 0) {
        return false;
    }
    
    int n_samples = bcf_hdr_nsamples(hdr_);
    dosages.clear();
    dosages.reserve(target_samples_.size());
    
    for (int idx : sample_indices_) {
        if (idx < n_samples && idx < n_values) {
            float ds_val = ds_values[idx];
            if (bcf_float_is_missing(ds_val)) {
                dosages.push_back(0.0);
            } else {
                dosages.push_back(static_cast<double>(ds_val));
            }
        } else {
            dosages.push_back(0.0);
        }
    }
    
    free(ds_values);
    return true;
}

std::vector<std::unique_ptr<VariantRecord>> VCFProcessor::processVariants() {
    Timer timer("VCF processing");
    
    std::vector<std::unique_ptr<VariantRecord>> variants;
    std::unordered_map<std::string, int> id_counts;
    
    variants.reserve(100000);  // Pre-allocate
    
    logInfo("Starting variant processing...");
    
    while (bcf_read(vcf_fp_, hdr_, rec_) == 0) {
        total_variants_++;
        
        // Unpack record for processing
        bcf_unpack(rec_, BCF_UN_ALL);
        
        // Extract basic info
        std::string id = rec_->d.id ? rec_->d.id : ".";
        std::string chrom = bcf_hdr_id2name(hdr_, rec_->rid);
        int pos = rec_->pos + 1;  // Convert to 1-based
        
        // Get REF and ALT
        std::string ref = rec_->d.allele[0];
        std::string alt = rec_->n_allele > 1 ? rec_->d.allele[1] : ".";
        
        // Extract allele frequency
        double af = extractAlleleFrequency(rec_);
        if (af < 0) {
            continue;  // Skip if AF cannot be determined
        }
        
        // Apply MAF filter FIRST (exactly like Go)
        if (af < maf_filter_ || af > (1.0 - maf_filter_)) {
            continue;
        }
        
        // Extract dosages
        std::vector<double> dosages;
        if (!extractDosages(rec_, dosages)) {
            continue;  // Skip if dosages cannot be extracted
        }
        
        // Create variant record
        auto variant = std::make_unique<VariantRecord>();
        variant->id = std::move(id);
        variant->chrom = std::move(chrom);
        variant->pos = pos;
        variant->ref = std::move(ref);
        variant->alt = std::move(alt);
        variant->af = af;
        variant->dosages = std::move(dosages);
        
        // Count IDs for MAF-filtered variants only (exactly like Go)
        id_counts[variant->id]++;
        variants.push_back(std::move(variant));
        
        if (variants.size() % 10000 == 0) {
            logInfo("Processed " + std::to_string(variants.size()) + " variants...");
        }
    }
    
    // Remove duplicates - EXACT GO LOGIC MATCHING
    std::vector<std::unique_ptr<VariantRecord>> final_variants;
    final_variants.reserve(variants.size());
    
    size_t original_size = variants.size();
    for (auto& variant : variants) {
        if (id_counts[variant->id] == 1) {
            final_variants.push_back(std::move(variant));
        }
    }
    
    duplicate_variants_ = original_size - final_variants.size();
    filtered_variants_ = final_variants.size();
    
    logInfo("Processed " + std::to_string(total_variants_) + " total variants");
    logInfo("Filtered to " + std::to_string(filtered_variants_) + " variants after MAF filtering");
    logInfo("Removed " + std::to_string(duplicate_variants_) + " duplicate variants");
    
    return final_variants;
}

