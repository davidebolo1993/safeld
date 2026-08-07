#include "vcf_processor.h"
#include "utils.h"
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <cerrno>
#include <climits>
#include <unordered_set>
#include <unistd.h>

namespace {
// Sentinel marking a not-yet-imputed missing genotype. Real dosages are always
// >= 0, so any negative value is unambiguously "missing" until imputeMissingWithMean runs.
constexpr double MISSING_DOSAGE = -1.0;

// Replace missing entries (negative sentinels) with the mean of the observed
// dosages. Mean-imputation keeps the variant's mean unchanged and avoids the
// downward bias that filling missing calls with 0.0 introduced.
//
// The all-missing case is only reachable when the caller has chosen to keep such
// a variant (max_missing_rate >= 1.0); the resulting constant row is dropped
// later by standardize(), which refuses zero-variance input.
void imputeMissingWithMean(std::vector<double>& dosages, const DosageStats& stats) {
    if (stats.observed_count == 0) {
        std::fill(dosages.begin(), dosages.end(), 0.0);
        return;
    }
    const double mean = stats.observed_sum / stats.observed_count;
    for (double& d : dosages) {
        if (d < 0.0) {
            d = mean;
        }
    }
}

// Fraction of samples with no usable call, computed before imputation.
double missingRate(const DosageStats& stats, size_t n_samples) {
    if (n_samples == 0) {
        return 1.0;
    }
    return 1.0 - static_cast<double>(stats.observed_count) / static_cast<double>(n_samples);
}

// ALT allele frequency over the observed calls only. Dosages are on the diploid
// 0..2 scale, so the denominator is 2 alleles per observed sample.
double alleleFrequencyFromObserved(const DosageStats& stats) {
    if (stats.observed_count == 0) {
        return -1.0;
    }
    const double af = stats.observed_sum / (2.0 * stats.observed_count);
    return std::round(af * 1000000.0) / 1000000.0;
}

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

struct DupState {
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

// Flip a spooled record's keep flag to 0 so the emit pass skips it.
void invalidateSpoolRecord(std::fstream& file, std::streamoff offset) {
    file.seekp(offset);
    const uint8_t keep = 0;
    file.write(reinterpret_cast<const char*>(&keep), sizeof(keep));
    if (!file) {
        throw std::runtime_error("Failed to invalidate duplicate spool record");
    }
    file.seekp(0, std::ios::end);
}
}  // namespace

VCFProcessor::VCFProcessor(const std::string& vcf_file, double maf_filter,
                           double max_missing_rate, const std::string& temp_dir,
                           DosageField dosage_field, bool use_info_af)
    : vcf_file_(vcf_file), maf_filter_(maf_filter), max_missing_rate_(max_missing_rate),
      dosage_field_(dosage_field), effective_field_(dosage_field), temp_dir_(temp_dir), vcf_fp_(nullptr),
      hdr_(nullptr), rec_(nullptr), use_info_af_requested_(use_info_af),
      use_info_af_(false), total_variants_(0),
      filtered_variants_(0), duplicate_variants_(0), multiallelic_variants_(0),
      missing_filtered_variants_(0), gt_fallback_variants_(0), gt_filled_calls_(0) {
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

    logDebug("Found " + std::to_string(n_samples) + " samples in VCF");
    return true;
}

void VCFProcessor::setupTargetSamples(const std::string& sample_list_str) {
    int n_samples = bcf_hdr_nsamples(hdr_);

    if (sample_list_str.empty()) {
        target_samples_.reserve(n_samples);
        sample_indices_.reserve(n_samples);
        for (int i = 0; i < n_samples; ++i) {
            target_samples_.emplace_back(hdr_->samples[i]);
            sample_indices_.push_back(i);
        }
    } else {
        auto requested_samples = split(sample_list_str, ',');
        std::unordered_map<std::string, int> sample_map;
        for (int i = 0; i < n_samples; ++i) {
            sample_map[hdr_->samples[i]] = i;
        }

        for (const auto& sample : requested_samples) {
            std::string trimmed = sample;
            trimmed.erase(0, trimmed.find_first_not_of(" \t"));
            trimmed.erase(trimmed.find_last_not_of(" \t") + 1);

            auto it = sample_map.find(trimmed);
            if (it != sample_map.end()) {
                target_samples_.push_back(trimmed);
                sample_indices_.push_back(it->second);
            }
        }

        logDebug("Found " + std::to_string(target_samples_.size()) + " of " +
                 std::to_string(requested_samples.size()) + " requested samples");
    }

    // INFO/AF is only a valid filter when it describes exactly these samples.
    // Even without -samples that is not guaranteed: subsetting a VCF upstream
    // typically recomputes AC and AN while leaving AF at its original value, so
    // trusting it silently filters on a frequency the data does not have.
    const bool whole_cohort = (sample_indices_.size() == static_cast<size_t>(n_samples));
    use_info_af_ = use_info_af_requested_ && whole_cohort;
    if (use_info_af_requested_ && !whole_cohort) {
        logWarning("-use-info-af ignored: a sample subset is in use, so INFO/AF does "
                   "not describe the selected samples.");
    }
    if (use_info_af_) {
        logWarning("-use-info-af: filtering on INFO/AF without checking it. This is "
                   "wrong if the field does not describe exactly the samples in this "
                   "file, which upstream subsetting commonly breaks.");
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

    scan_ = scanInput();
    reportScanAndChooseField();

    logDebug("VCF processor initialized with " + std::to_string(target_samples_.size()) + " samples");
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

// Fill `dosages` with one entry per target sample, using MISSING_DOSAGE for
// samples with no usable call, and report the observed totals in `stats`.
// Imputation is deliberately left to the caller so that allele frequency and
// missingness can be measured against real calls first.
bool VCFProcessor::extractDosages(bcf1_t* rec, std::vector<double>& dosages, DosageStats& stats) {
    if (effective_field_ == DosageField::GT) {
        return extractDosagesFromGT(rec, dosages, stats);
    }

    int n_values = 0;
    float* ds_values = nullptr;
    int ret = bcf_get_format_float(hdr_, rec, "DS", &ds_values, &n_values);

    if (ret <= 0) {
        // No DS field on this record.
        free(ds_values);
        // Forcing DS on a record that has none: the genotype is still known
        // from GT, so read that rather than dropping the variant.
        return extractDosagesFromGT(rec, dosages, stats);
    }

    int n_samples = bcf_hdr_nsamples(hdr_);
    // DS is one value per ALT allele; records are biallelic here, so anything
    // other than one value per sample is a shape we cannot index by sample.
    if (n_samples <= 0 || n_values / n_samples != 1) {
        free(ds_values);
        return extractDosagesFromGT(rec, dosages, stats);
    }

    dosages.clear();
    dosages.reserve(target_samples_.size());
    stats = DosageStats{};

    for (int idx : sample_indices_) {
        if (idx >= n_samples) {
            dosages.push_back(MISSING_DOSAGE);
            continue;
        }

        float ds_val = ds_values[idx];
        // The htslib sentinels are NaN bit patterns, and the MISSING_DOSAGE
        // scheme relies on stored dosages being non-negative, so reject
        // anything that is not a real value in [0, inf).
        if (bcf_float_is_missing(ds_val) || bcf_float_is_vector_end(ds_val) ||
            !(ds_val >= 0.0f)) {
            dosages.push_back(MISSING_DOSAGE);
            continue;
        }

        double d = static_cast<double>(ds_val);
        dosages.push_back(d);
        stats.observed_sum += d;
        stats.observed_count++;
    }

    free(ds_values);

    // A VCF may carry GT for every sample while omitting the DS subfield for
    // many of them: plink2 exports look exactly like "0|1:0.97" sitting next to
    // a bare "0|0". Such a sample is NOT missing - its genotype is known from
    // GT - so mean-imputing it would invent data and attenuate every pairwise r2
    // in proportion to the DS presence rate. Fill each gap from that sample's own
    // hard call instead, keeping the real dosage wherever one was written.
    if (stats.observed_count < static_cast<int>(dosages.size())) {
        std::vector<double> gt_dosages;
        DosageStats gt_stats;
        if (extractDosagesFromGT(rec, gt_dosages, gt_stats) &&
            gt_dosages.size() == dosages.size()) {
            int filled = 0;
            for (size_t i = 0; i < dosages.size(); ++i) {
                if (dosages[i] < 0.0 && gt_dosages[i] >= 0.0) {
                    dosages[i] = gt_dosages[i];
                    filled++;
                }
            }
            if (filled > 0) {
                // Recompute the observed totals over the merged vector.
                stats = DosageStats{};
                for (double d : dosages) {
                    if (d >= 0.0) {
                        stats.observed_sum += d;
                        stats.observed_count++;
                    }
                }
                gt_fallback_variants_++;
                gt_filled_calls_ += filled;
            }
        }
    }

    return true;
}

// Derive an ALT-allele dosage on the diploid 0..2 scale from GT hard calls when
// DS is absent. Each sample is normalised by its own ploidy, so a hemizygous ALT
// call (chrX/chrY in a male, mitochondria) scores 2.0 and cannot be mistaken for
// a heterozygote. A genotype with any missing allele (./1) is treated as missing
// outright rather than being silently counted as a reference call.
bool VCFProcessor::extractDosagesFromGT(bcf1_t* rec, std::vector<double>& dosages,
                                        DosageStats& stats) {
    int n_gt = 0;
    int32_t* gt_arr = nullptr;
    int ret = bcf_get_genotypes(hdr_, rec, &gt_arr, &n_gt);
    if (ret <= 0) {
        free(gt_arr);
        return false;
    }

    int n_samples = bcf_hdr_nsamples(hdr_);
    int max_ploidy = n_samples > 0 ? n_gt / n_samples : 0;
    if (max_ploidy <= 0) {
        free(gt_arr);
        return false;
    }

    dosages.clear();
    dosages.reserve(target_samples_.size());
    stats = DosageStats{};

    for (int idx : sample_indices_) {
        if (idx >= n_samples) {
            dosages.push_back(MISSING_DOSAGE);
            continue;
        }

        const int32_t* ptr = gt_arr + static_cast<size_t>(idx) * max_ploidy;
        int alt_count = 0;
        int ploidy = 0;
        bool any_missing = false;
        for (int p = 0; p < max_ploidy; ++p) {
            if (ptr[p] == bcf_int32_vector_end) {
                break;  // sample has lower ploidy than the record max
            }
            if (bcf_gt_is_missing(ptr[p])) {
                any_missing = true;
                break;
            }
            ploidy++;
            // Records are biallelic here, so any non-reference allele is ALT.
            if (bcf_gt_allele(ptr[p]) > 0) {
                alt_count++;
            }
        }

        if (any_missing || ploidy == 0) {
            dosages.push_back(MISSING_DOSAGE);
            continue;
        }

        double dosage = 2.0 * alt_count / ploidy;
        dosages.push_back(dosage);
        stats.observed_sum += dosage;
        stats.observed_count++;
    }

    free(gt_arr);
    return true;
}

// Location of the deduplication spool. It holds the full dosage vector of every
// kept variant, so it must land on a filesystem with room for a second copy of
// the genotype matrix, not on a small /tmp.
std::string VCFProcessor::makeSpoolPath() const {
    std::string dir = temp_dir_;
    if (dir.empty()) {
        const char* env_tmp = std::getenv("TMPDIR");
        dir = (env_tmp && *env_tmp) ? env_tmp : "/tmp";
    }
    while (dir.size() > 1 && dir.back() == '/') {
        dir.pop_back();
    }
    return dir + "/safeld_dedup_XXXXXX";
}

// Reads the head of the file on a separate handle and summarises what it finds.
// Cheap relative to the run and it removes the need for an external script.
InputScan VCFProcessor::scanInput(int max_records) {
    InputScan scan;
    if (sample_indices_.empty()) {
        return scan;
    }

    htsFile* fp = hts_open(vcf_file_.c_str(), "r");
    if (!fp) {
        return scan;
    }
    bcf_hdr_t* hdr = bcf_hdr_read(fp);
    if (!hdr) {
        hts_close(fp);
        return scan;
    }
    bcf1_t* rec = bcf_init();
    if (!rec) {
        bcf_hdr_destroy(hdr);
        hts_close(fp);
        return scan;
    }

    scan.n_samples = static_cast<int>(sample_indices_.size());
    scan.gt_declared = bcf_hdr_id2int(hdr, BCF_DT_ID, "GT") >= 0;
    scan.ds_declared = bcf_hdr_id2int(hdr, BCF_DT_ID, "DS") >= 0;

    const int file_samples = bcf_hdr_nsamples(hdr);
    double ds_rate_sum = 0.0;
    double gt_rate_sum = 0.0;

    float* ds_values = nullptr;
    int32_t* gt_arr = nullptr;
    int n_ds = 0, n_gt = 0;

    while (scan.records < max_records && bcf_read(fp, hdr, rec) == 0) {
        scan.records++;
        bcf_unpack(rec, BCF_UN_ALL);

        if (rec->n_allele != 2) {
            scan.multiallelic++;
            continue;
        }
        scan.biallelic++;

        if (!rec->d.id || rec->d.id[0] == '.') {
            scan.no_id++;
        }

        int n_af = 0;
        float* af_values = nullptr;
        if (bcf_get_info_float(hdr, rec, "AF", &af_values, &n_af) > 0 && n_af > 0) {
            scan.with_info_af++;
        }
        free(af_values);

        // DS presence: htslib fills the missing sentinel for samples whose DS
        // subfield was simply omitted, which is precisely what we want to count.
        int ds_present = 0;
        if (scan.ds_declared &&
            bcf_get_format_float(hdr, rec, "DS", &ds_values, &n_ds) > 0 &&
            file_samples > 0 && n_ds / file_samples == 1) {
            for (int idx : sample_indices_) {
                if (idx >= file_samples) continue;
                float v = ds_values[idx];
                if (!bcf_float_is_missing(v) && !bcf_float_is_vector_end(v) && v >= 0.0f) {
                    ds_present++;
                }
            }
        }
        if (ds_present == 0) scan.ds_absent_records++;
        ds_rate_sum += static_cast<double>(ds_present) / scan.n_samples;

        int gt_called = 0;
        if (scan.gt_declared && bcf_get_genotypes(hdr, rec, &gt_arr, &n_gt) > 0 && file_samples > 0) {
            const int max_ploidy = n_gt / file_samples;
            if (max_ploidy > 0) {
                for (int idx : sample_indices_) {
                    if (idx >= file_samples) continue;
                    const int32_t* ptr = gt_arr + static_cast<size_t>(idx) * max_ploidy;
                    bool ok = false;
                    for (int p = 0; p < max_ploidy; ++p) {
                        if (ptr[p] == bcf_int32_vector_end) break;
                        if (bcf_gt_is_missing(ptr[p])) { ok = false; break; }
                        ok = true;
                    }
                    if (ok) gt_called++;
                }
            }
        }
        gt_rate_sum += static_cast<double>(gt_called) / scan.n_samples;
    }

    free(ds_values);
    free(gt_arr);
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    hts_close(fp);

    if (scan.biallelic > 0) {
        scan.ds_presence  = ds_rate_sum / static_cast<double>(scan.biallelic);
        scan.gt_call_rate = gt_rate_sum / static_cast<double>(scan.biallelic);
    }
    scan.ok = scan.records > 0;
    return scan;
}

// Report what the scan found, then settle which FORMAT field the run will read.
void VCFProcessor::reportScanAndChooseField() {
    if (!scan_.ok) {
        logWarning("Could not scan the input; proceeding without a summary");
        effective_field_ = dosage_field_;
        return;
    }

    std::ostringstream head;
    head << "Input: " << formatCount(scan_.n_samples) << " samples, scanned "
         << formatCount(scan_.records) << " records";
    logInfo(head.str());

    if (scan_.multiallelic > 0) {
        logInfo("  " + formatCount(scan_.multiallelic) +
                " non-biallelic (will be skipped; split with bcftools norm -m -any)");
    }
    if (scan_.biallelic == 0) {
        logWarning("  no biallelic records in the scanned window");
        effective_field_ = dosage_field_;
        return;
    }
    if (scan_.no_id > 0) {
        logDebug("  " + formatCount(scan_.no_id) + " records without an ID (deduplicated by locus)");
    }
    if (scan_.with_info_af == 0) {
        logDebug("  no INFO/AF: allele frequencies will be computed from the genotypes");
    }

    auto pct = [](double v) {
        std::ostringstream os;
        os << std::fixed << std::setprecision(1) << (v * 100.0) << "%";
        return os.str();
    };

    if (scan_.gt_declared) {
        logInfo("  GT present for " + pct(scan_.gt_call_rate) + " of calls");
    }
    if (scan_.ds_declared) {
        logInfo("  DS present for " + pct(scan_.ds_presence) + " of calls");
    }

    // Resolve Auto. The decision that matters: a file can carry GT for every
    // sample while writing DS for only a fraction of them, and treating those
    // absences as missing genotypes attenuates every pairwise r2 in proportion
    // to the DS presence rate.
    effective_field_ = dosage_field_;
    if (dosage_field_ == DosageField::Auto) {
        if (!scan_.ds_declared || scan_.ds_presence <= 0.0) {
            effective_field_ = DosageField::GT;
        } else if (scan_.ds_presence >= 1.0 - max_missing_rate_) {
            effective_field_ = DosageField::DS;
        } else if (scan_.gt_call_rate > scan_.ds_presence) {
            effective_field_ = DosageField::GT;
            logWarning("DS covers only " + pct(scan_.ds_presence) + " of calls while GT covers " +
                       pct(scan_.gt_call_rate) + "; reading GT so the matrix comes from one "
                       "complete field. Override with -dosage-field DS.");
        } else {
            effective_field_ = DosageField::DS;
        }
    }

    if (dosage_field_ == DosageField::DS && !scan_.ds_declared) {
        logWarning("-dosage-field DS was requested but the file declares no DS; "
                   "reading GT hard calls instead.");
        effective_field_ = DosageField::GT;
    }
    if (dosage_field_ == DosageField::GT && !scan_.gt_declared) {
        logWarning("-dosage-field GT was requested but the file declares no GT; "
                   "reading DS instead.");
        effective_field_ = DosageField::DS;
    }

    const char* chosen = (effective_field_ == DosageField::GT) ? "GT" : "DS";
    const char* how = (dosage_field_ == DosageField::Auto) ? " (auto)" : " (forced)";
    logInfo(std::string("Dosage source: ") + chosen + how);
}

void VCFProcessor::setExtractIds(std::vector<std::string> ids) {
    extract_ids_.clear();
    extract_ids_.reserve(ids.size());
    for (auto& id : ids) {
        extract_ids_.insert(std::move(id));
    }
}

// stream variants with one-pass duplicate tracking.
void VCFProcessor::streamVariants(VariantCallback callback) {
    Timer timer("VCF streaming");
    std::unordered_map<std::string, DupState> dup_states;
    std::unordered_set<std::streamoff> retracted_offsets;
    std::unordered_map<std::string, int> contig_rank;

    total_variants_ = 0;
    filtered_variants_ = 0;
    duplicate_variants_ = 0;
    multiallelic_variants_ = 0;
    missing_filtered_variants_ = 0;
    gt_fallback_variants_ = 0;
    gt_filled_calls_ = 0;
    not_extracted_ = 0;
    maf_filtered_ = 0;
    extract_seen_.clear();

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

    auto failsMafFilter = [this](double af) {
        return af < 0.0 || af < maf_filter_ || af > (1.0 - maf_filter_);
    };

    std::string spool_template = makeSpoolPath();
    std::vector<char> spool_path_buf(spool_template.begin(), spool_template.end());
    spool_path_buf.push_back('\0');
    int temp_fd = mkstemp(spool_path_buf.data());
    if (temp_fd < 0) {
        throw std::runtime_error("Failed to create temporary deduplication spool file at " +
                                 spool_template + ": " + std::string(strerror(errno)));
    }
    close(temp_fd);
    TempFileGuard temp_guard{std::string(spool_path_buf.data())};
    logDebug("Deduplication spool: " + temp_guard.path);

    std::fstream spool(temp_guard.path, std::ios::in | std::ios::out | std::ios::binary | std::ios::trunc);
    if (!spool) {
        throw std::runtime_error("Failed to open temporary deduplication spool file");
    }

    logDebug("Starting single-pass variant scan with disk-backed duplicate handling...");

    // make one pass over the VCF, spooling first occurrences and invalidating them if duplicates appear.
    bool has_prev_coord = false;
    std::string prev_chrom;
    int prev_pos = 0;

    ProgressCounter scan("Scanning", "variants");
    while (bcf_read(vcf_fp_, hdr_, rec_) == 0) {
        total_variants_++;
        scan.increment();
        bcf_unpack(rec_, BCF_UN_ALL);

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

        // Only biallelic records are meaningful downstream: a single ALT column,
        // a single INFO/AF value and one dosage per sample. A multiallelic record
        // would collapse several ALT alleles into one column, so require it to be
        // split upstream (bcftools norm -m -any).
        if (rec_->n_allele != 2) {
            multiallelic_variants_++;
            continue;
        }

        std::string id = rec_->d.id ? rec_->d.id : ".";
        std::string ref = rec_->d.allele[0];
        std::string alt = rec_->d.allele[1];

        if (!extract_ids_.empty()) {
            if (!extract_ids_.contains(id)) {
                not_extracted_++;
                continue;
            }
            extract_seen_.insert(id);
        }

        double af = -1.0;
        std::vector<double> dosages;
        DosageStats stats;
        bool have_dosages = false;

        // INFO/AF is the cheap path: it lets whole-cohort runs reject a variant
        // before touching per-sample data.
        if (use_info_af_ && tryExtractAfFromInfo(rec_, af)) {
            if (failsMafFilter(af)) {
                maf_filtered_++;
                continue;
            }
        } else {
            if (!extractDosages(rec_, dosages, stats)) {
                continue;
            }
            have_dosages = true;
            // With nothing observed there is no frequency to test. Report that
            // as missingness rather than letting it fall through the MAF filter,
            // where it would be counted as a frequency exclusion it never had.
            if (stats.observed_count == 0) {
                missing_filtered_variants_++;
                continue;
            }
            af = alleleFrequencyFromObserved(stats);
            if (failsMafFilter(af)) {
                maf_filtered_++;
                continue;
            }
        }

        if (!have_dosages && !extractDosages(rec_, dosages, stats)) {
            continue;
        }

        // Drop variants that are mostly uncalled. Without this a variant whose
        // calls are all missing passes the INFO/AF filter, is imputed to a
        // constant vector, and surfaces as a synthetic dosage of exactly 1.0 for
        // every trait with no warning anywhere.
        if (missingRate(stats, dosages.size()) > max_missing_rate_) {
            missing_filtered_variants_++;
            continue;
        }
        imputeMissingWithMean(dosages, stats);

        // Deduplicate on the locus, plus on the ID when the record actually has
        // one. Keying on the ID alone silently discarded every ID-less record,
        // because they all share the "." placeholder and so looked like copies of
        // each other.
        std::string locus_key = "L:" + chrom + ":" + std::to_string(pos) + ":" + ref + ":" + alt;
        std::vector<std::string> dup_keys{locus_key};
        if (id != ".") {
            dup_keys.push_back("I:" + id);
        }

        bool is_duplicate = false;
        for (const auto& key : dup_keys) {
            auto& state = dup_states[key];
            state.count++;
            if (state.count == 1) {
                continue;
            }
            is_duplicate = true;
            // Retract the first copy once, then skip every later one.
            if (state.has_spooled_record) {
                if (retracted_offsets.insert(state.keep_offset).second) {
                    invalidateSpoolRecord(spool, state.keep_offset);
                    duplicate_variants_++;
                }
                state.has_spooled_record = false;
            }
        }
        if (is_duplicate) {
            duplicate_variants_++;
            continue;
        }

        Variant variant;
        variant.id = std::move(id);
        variant.chrom = std::move(chrom);
        variant.pos = pos;
        variant.ref = std::move(ref);
        variant.alt = std::move(alt);
        variant.af = af;
        variant.dosages = std::move(dosages);

        std::streamoff keep_offset = -1;
        writeSpoolRecord(spool, variant, keep_offset);
        for (const auto& key : dup_keys) {
            auto& state = dup_states[key];
            state.keep_offset = keep_offset;
            state.has_spooled_record = true;
        }
    }

    scan.finish();

    // read the spool and emit only records still marked as kept.
    spool.flush();
    spool.clear();
    spool.seekg(0, std::ios::beg);

    ProgressBar emit_bar("Emitting variants",
                         total_variants_ - duplicate_variants_ - multiallelic_variants_ -
                         missing_filtered_variants_ - not_extracted_);
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
        emit_bar.increment();
    }
    emit_bar.finish();

    // Report the totals, then only the exclusions that actually happened: a wall
    // of zeroes buries the one line that matters.
    logInfo("Scanned " + formatCount(total_variants_) + " variants, kept " +
            formatCount(filtered_variants_));
    if (not_extracted_ > 0) {
        logInfo("  excluded " + formatCount(not_extracted_) + " not in the extract list");
    }
    if (multiallelic_variants_ > 0) {
        logInfo("  excluded " + formatCount(multiallelic_variants_) + " non-biallelic");
    }
    if (maf_filtered_ > 0) {
        std::ostringstream mf;
        mf << std::fixed << std::setprecision(4) << maf_filter_;
        logInfo("  excluded " + formatCount(maf_filtered_) + " below the MAF threshold (" + mf.str() + ")");
    }
    if (missing_filtered_variants_ > 0) {
        std::ostringstream mm;
        mm << std::fixed << std::setprecision(3) << max_missing_rate_;
        logInfo("  excluded " + formatCount(missing_filtered_variants_) +
                " over the missingness limit (" + mm.str() + ")");
    }
    if (duplicate_variants_ > 0) {
        logInfo("  excluded " + formatCount(duplicate_variants_) +
                " duplicated by locus or ID (every copy)");
    }
    if (gt_fallback_variants_ > 0) {
        logWarning("DS was absent for " + std::to_string(gt_filled_calls_) +
                   " call(s) across " + std::to_string(gt_fallback_variants_) +
                   " variant(s); filled from each sample's own GT hard call "
                   "rather than imputing. Use -dosage-field GT for a matrix built "
                   "from one field throughout.");
    }

    counts_.total = total_variants_;
    counts_.emitted = filtered_variants_;
    counts_.multiallelic = multiallelic_variants_;
    counts_.missing_filtered = missing_filtered_variants_;
    counts_.duplicates = duplicate_variants_;
    counts_.not_extracted = not_extracted_;
    counts_.maf_filtered = maf_filtered_;
    counts_.extract_unmatched = static_cast<long long>(extract_ids_.size()) -
                                static_cast<long long>(extract_seen_.size());

    // Silence here would be the worst outcome: an ID list that matches nothing
    // looks exactly like a successful run that happened to keep no variants.
    if (!extract_ids_.empty() && counts_.extract_unmatched > 0) {
        logWarning(formatCount(counts_.extract_unmatched) + " of " +
                   formatCount(extract_ids_.size()) + " extract IDs matched no variant. "
                   "The list must use the same IDs as the input's ID column.");
    }
    counts_.hardcall_filled = gt_filled_calls_;
    counts_.hardcall_filled_variants = gt_fallback_variants_;

    if (filtered_variants_ == 0) {
        logWarning("No variants passed filtering; downstream stages will have nothing to process");
    }
}
