#include "vcf_processor.h"
#include "utils.h"
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstring>
#include <climits>
#include <unistd.h>

namespace {
struct SpoolRecordHeader {
    uint8_t keep;
    int32_t pos;
    double af;
    uint32_t id_len;
    uint32_t chrom_len;
    uint32_t ref_len;
    uint32_t alt_len;
    uint32_t dosages_count;
};

struct IdState {
    std::streamoff keep_offset = -1;
    int count = 0;
    bool has_spooled_record = false;
};

struct TempFileGuard {
    std::string path;
    ~TempFileGuard() {
        if (!path.empty()) {
            std::remove(path.c_str());
        }
    }
};

void writeSpoolRecord(std::fstream& file, const Variant& variant, std::streamoff& keep_offset) {
    keep_offset = file.tellp();

    SpoolRecordHeader header{};
    header.keep = 1;
    header.pos = static_cast<int32_t>(variant.pos);
    header.af = variant.af;
    header.id_len = static_cast<uint32_t>(variant.id.size());
    header.chrom_len = static_cast<uint32_t>(variant.chrom.size());
    header.ref_len = static_cast<uint32_t>(variant.ref.size());
    header.alt_len = static_cast<uint32_t>(variant.alt.size());
    header.dosages_count = static_cast<uint32_t>(variant.dosages.size());

    file.write(reinterpret_cast<const char*>(&header), sizeof(header));
    file.write(variant.id.data(), static_cast<std::streamsize>(variant.id.size()));
    file.write(variant.chrom.data(), static_cast<std::streamsize>(variant.chrom.size()));
    file.write(variant.ref.data(), static_cast<std::streamsize>(variant.ref.size()));
    file.write(variant.alt.data(), static_cast<std::streamsize>(variant.alt.size()));
    file.write(reinterpret_cast<const char*>(variant.dosages.data()),
               static_cast<std::streamsize>(variant.dosages.size() * sizeof(double)));

    if (!file) {
        throw std::runtime_error("Failed writing temporary deduplication spool record");
    }
}

bool readSpoolRecord(std::fstream& file, uint8_t& keep, Variant& variant) {
    SpoolRecordHeader header{};
    file.read(reinterpret_cast<char*>(&header), sizeof(header));

    if (file.eof()) {
        return false;
    }
    if (!file) {
        throw std::runtime_error("Failed reading temporary deduplication spool header");
    }

    keep = header.keep;
    variant.pos = header.pos;
    variant.af = header.af;
    variant.id.resize(header.id_len);
    variant.chrom.resize(header.chrom_len);
    variant.ref.resize(header.ref_len);
    variant.alt.resize(header.alt_len);
    variant.dosages.resize(header.dosages_count);

    file.read(variant.id.data(), static_cast<std::streamsize>(variant.id.size()));
    file.read(variant.chrom.data(), static_cast<std::streamsize>(variant.chrom.size()));
    file.read(variant.ref.data(), static_cast<std::streamsize>(variant.ref.size()));
    file.read(variant.alt.data(), static_cast<std::streamsize>(variant.alt.size()));
    file.read(reinterpret_cast<char*>(variant.dosages.data()),
              static_cast<std::streamsize>(variant.dosages.size() * sizeof(double)));

    if (!file) {
        throw std::runtime_error("Failed reading temporary deduplication spool record body");
    }

    return true;
}
}  // namespace

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

std::vector<std::string> VCFProcessor::getContigNames() const {
    std::vector<std::string> contigs;
    if (!hdr_) {
        return contigs;
    }

    int nseq = 0;
    const char** seqnames = bcf_hdr_seqnames(hdr_, &nseq);
    if (!seqnames || nseq <= 0) {
        free(const_cast<char**>(seqnames));
        return contigs;
    }

    contigs.reserve(nseq);
    for (int i = 0; i < nseq; ++i) {
        contigs.emplace_back(seqnames[i]);
    }
    free(const_cast<char**>(seqnames));
    return contigs;
}

double VCFProcessor::extractAlleleFrequency(bcf1_t* rec) {
    // Try to get AF from INFO field
    int n_values = 0;
    float* af_values = nullptr;
    if (bcf_get_info_float(hdr_, rec, "AF", &af_values, &n_values) > 0 && n_values > 0) {
        double af = static_cast<double>(af_values[0]);
        free(af_values);
        // Ensure consistent precision (6 decimal places)
        return std::round(af * 1000000.0) / 1000000.0;
    }

    // If AF not available, calculate from dosages
    std::vector<double> dosages;
    if (extractDosages(rec, dosages)) {
        double sum = 0.0;
        int count = 0;
        for (double d : dosages) {
            if (d >= 0.0) {
                sum += d;
                count++;
            }
        }
        if (count > 0) {
            double af = sum / (2.0 * count);
            return std::round(af * 1000000.0) / 1000000.0;
        }
    }

    return -1.0; // Unable to determine AF
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

// Streaming with duplicate detection.
void VCFProcessor::streamVariants(VariantCallback callback) {
    Timer timer("VCF streaming");
    std::unordered_map<std::string, IdState> id_states;
    std::unordered_map<std::string, int> contig_rank;

    total_variants_ = 0;
    filtered_variants_ = 0;
    duplicate_variants_ = 0;

    {
        int rank = 0;
        for (const auto& contig : getContigNames()) {
            contig_rank[contig] = rank++;
        }
    }

    auto tryExtractAfFromInfo = [&](bcf1_t* rec, double& af) -> bool {
        int n_values = 0;
        float* af_values = nullptr;
        int ret = bcf_get_info_float(hdr_, rec, "AF", &af_values, &n_values);
        if (ret > 0 && n_values > 0) {
            af = std::round(static_cast<double>(af_values[0]) * 1000000.0) / 1000000.0;
            free(af_values);
            return true;
        }
        free(af_values);
        return false;
    };

    auto computeAfFromDosages = [](const std::vector<double>& dosages) -> double {
        double sum = 0.0;
        int count = 0;
        for (double d : dosages) {
            if (d >= 0.0) {
                sum += d;
                count++;
            }
        }
        if (count == 0) {
            return -1.0;
        }
        return std::round((sum / (2.0 * count)) * 1000000.0) / 1000000.0;
    };

    char temp_path_template[] = "/tmp/safeld_dedup_XXXXXX";
    int temp_fd = mkstemp(temp_path_template);
    if (temp_fd < 0) {
        throw std::runtime_error("Failed to create temporary deduplication spool file");
    }
    close(temp_fd);
    TempFileGuard temp_guard{temp_path_template};

    std::fstream spool(temp_guard.path, std::ios::in | std::ios::out | std::ios::binary | std::ios::trunc);
    if (!spool) {
        throw std::runtime_error("Failed to open temporary deduplication spool file");
    }

    logInfo("Starting single-pass variant scan with disk-backed duplicate handling...");

    // Single pass on VCF: spool first occurrences, invalidate if duplicates appear later.
    bool has_prev_coord = false;
    std::string prev_chrom;
    int prev_pos = 0;

    while (bcf_read(vcf_fp_, hdr_, rec_) == 0) {
        total_variants_++;
        bcf_unpack(rec_, BCF_UN_ALL);
        if (total_variants_ % 50000 == 0) {
            logInfo("Scanned " + std::to_string(total_variants_) + " variants...");
        }

        std::string chrom = bcf_hdr_id2name(hdr_, rec_->rid);
        int pos = rec_->pos + 1;
        if (has_prev_coord) {
            int prev_rank = contig_rank.contains(prev_chrom) ? contig_rank[prev_chrom] : INT_MAX;
            int curr_rank = contig_rank.contains(chrom) ? contig_rank[chrom] : INT_MAX;

            bool out_of_order = false;
            if (prev_rank != INT_MAX && curr_rank != INT_MAX) {
                out_of_order = (curr_rank < prev_rank) || (curr_rank == prev_rank && pos < prev_pos);
            } else if (chrom != prev_chrom) {
                out_of_order = chrom < prev_chrom;
            } else {
                out_of_order = pos < prev_pos;
            }

            if (out_of_order) {
                throw std::runtime_error(
                    "Input VCF is not coordinate-sorted. Offending record: " + chrom + ":" +
                    std::to_string(pos) + " after " + prev_chrom + ":" + std::to_string(prev_pos));
            }
        }
        prev_chrom = chrom;
        prev_pos = pos;
        has_prev_coord = true;

        std::string id = rec_->d.id ? rec_->d.id : ".";
        double af = -1.0;
        std::vector<double> dosages;

        if (!tryExtractAfFromInfo(rec_, af)) {
            if (!extractDosages(rec_, dosages)) {
                continue;
            }
            af = computeAfFromDosages(dosages);
        }

        if (af < 0 || af < maf_filter_ || af > (1.0 - maf_filter_)) {
            continue;
        }

        auto& state = id_states[id];
        state.count++;
        if (state.count > 1) {
            // This ID is duplicated among MAF-passing variants.
            // Invalidate the first spooled record once, then skip all later duplicates.
            if (state.count == 2 && state.has_spooled_record) {
                spool.seekp(state.keep_offset);
                uint8_t keep = 0;
                spool.write(reinterpret_cast<const char*>(&keep), sizeof(keep));
                if (!spool) {
                    throw std::runtime_error("Failed to invalidate duplicate spool record");
                }
                spool.seekp(0, std::ios::end);
            }
            continue;
        }

        if (dosages.empty()) {
            if (!extractDosages(rec_, dosages)) {
                // Keep ID state for duplicate accounting semantics, but don't spool.
                state.has_spooled_record = false;
                continue;
            }
        }

        Variant variant;
        variant.id = std::move(id);
        variant.chrom = std::move(chrom);
        variant.pos = pos;
        variant.ref = rec_->d.allele[0];
        variant.alt = rec_->n_allele > 1 ? rec_->d.allele[1] : ".";
        variant.af = af;
        variant.dosages = std::move(dosages);

        std::streamoff keep_offset = -1;
        writeSpoolRecord(spool, variant, keep_offset);
        state.keep_offset = keep_offset;
        state.has_spooled_record = true;
    }

    duplicate_variants_ = 0;
    for (const auto& [_, state] : id_states) {
        if (state.count > 1) {
            duplicate_variants_ += state.count;
        }
    }

    // Second phase on spool only (not on VCF): emit only records still marked as kept.
    spool.flush();
    spool.clear();
    spool.seekg(0, std::ios::beg);

    int emitted = 0;
    uint8_t keep = 0;
    Variant spooled_variant;
    while (readSpoolRecord(spool, keep, spooled_variant)) {
        if (keep == 0) {
            continue;
        }

        auto variant = std::make_unique<Variant>();
        variant->id = std::move(spooled_variant.id);
        variant->chrom = std::move(spooled_variant.chrom);
        variant->pos = spooled_variant.pos;
        variant->ref = std::move(spooled_variant.ref);
        variant->alt = std::move(spooled_variant.alt);
        variant->af = spooled_variant.af;
        variant->dosages = std::move(spooled_variant.dosages);
        callback(std::move(variant));

        filtered_variants_++;
        emitted++;

        if (emitted % 10000 == 0) {
            logInfo("Emitted " + std::to_string(emitted) + " deduplicated variants...");
        }
    }

    logInfo("Streaming complete!");
    logInfo("Total variants scanned: " + std::to_string(total_variants_));
    logInfo("Variants after MAF filter: " + std::to_string(filtered_variants_));
    logInfo("Duplicate variants skipped: " + std::to_string(duplicate_variants_));
}
