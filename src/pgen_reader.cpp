#include "pgen_reader.h"

#include "utils.h"

#include <cstring>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <iomanip>

#include "pgenlib_read.h"
#include "pgenlib_misc.h"

namespace {

constexpr double kMissingDosage = -1.0;

// pgenlib stores a dosage as a uint16 where 16384 == 1.0 ALT allele, so the
// diploid 0..2 scale is val / 16384.0.
constexpr double kDosageScale = 1.0 / 16384.0;

std::string stripSuffix(const std::string& path, const std::string& suffix) {
    if (path.size() >= suffix.size() &&
        path.compare(path.size() - suffix.size(), suffix.size(), suffix) == 0) {
        return path.substr(0, path.size() - suffix.size());
    }
    return path;
}

// 2-bit hard call codes: 0 = hom ref, 1 = het, 2 = hom alt, 3 = missing.
inline uint32_t genoAt(const uintptr_t* genovec, uint32_t idx) {
    constexpr uint32_t kPerWord = plink2::kBitsPerWordD2;
    return (genovec[idx / kPerWord] >> (2 * (idx % kPerWord))) & 3;
}

inline bool bitAt(const uintptr_t* bitarr, uint32_t idx) {
    return (bitarr[idx / plink2::kBitsPerWord] >> (idx % plink2::kBitsPerWord)) & 1;
}

}  // namespace

struct PgenSource::Impl {
    plink2::PgenFileInfo pgfi;
    plink2::PgenReader pgr;
    unsigned char* pgfi_alloc = nullptr;
    unsigned char* pgr_alloc = nullptr;
    bool pgfi_ready = false;
    bool pgr_ready = false;

    // Per-variant scratch, sized once at open().
    std::vector<uintptr_t> genovec;
    std::vector<uintptr_t> dosage_present;
    std::vector<uint16_t> dosage_main;

    Impl() {
        plink2::PreinitPgfi(&pgfi);
        plink2::PreinitPgr(&pgr);
    }

    ~Impl() {
        plink2::PglErr reterr = plink2::kPglRetSuccess;
        if (pgr_ready) plink2::CleanupPgr(&pgr, &reterr);
        if (pgfi_ready) plink2::CleanupPgfi(&pgfi, &reterr);
        plink2::aligned_free_cond(pgr_alloc);
        plink2::aligned_free_cond(pgfi_alloc);
    }
};

PgenSource::PgenSource(const std::string& genotype_path, Metadata metadata,
                       double maf_filter, double max_missing_rate, DosageField dosage_field)
    : impl_(new Impl()), genotype_path_(genotype_path), metadata_(metadata),
      maf_filter_(maf_filter), max_missing_rate_(max_missing_rate),
      dosage_field_(dosage_field) {}

PgenSource::~PgenSource() { delete impl_; }

void PgenSource::setExtractIds(std::vector<std::string> ids) {
    extract_ids_ = std::move(ids);
}

std::vector<std::string> PgenSource::getContigNames() const {
    // Contig order as first seen in the variant table.
    std::vector<std::string> contigs;
    std::string last;
    for (const auto& v : variants_) {
        if (v.chrom != last) {
            if (std::find(contigs.begin(), contigs.end(), v.chrom) == contigs.end()) {
                contigs.push_back(v.chrom);
            }
            last = v.chrom;
        }
    }
    return contigs;
}

bool PgenSource::loadPsam(const std::string& path, std::string& error) {
    std::ifstream in(path);
    if (!in) {
        error = "cannot open " + path;
        return false;
    }
    sample_ids_.clear();
    std::string line;
    int iid_col = -1;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '#') {
            // Header names the columns; IID is the sample identifier.
            std::istringstream hs(line);
            std::string tok;
            int col = 0;
            while (hs >> tok) {
                if (tok == "IID" || tok == "#IID") iid_col = col;
                col++;
            }
            continue;
        }
        std::istringstream ls(line);
        std::vector<std::string> fields;
        std::string tok;
        while (ls >> tok) fields.push_back(tok);
        if (fields.empty()) continue;
        // Without a header, plink writes FID IID ...; a single column is IID.
        int use = iid_col;
        if (use < 0) use = (fields.size() >= 2) ? 1 : 0;
        if (use >= static_cast<int>(fields.size())) use = static_cast<int>(fields.size()) - 1;
        sample_ids_.push_back(fields[use]);
    }
    if (sample_ids_.empty()) {
        error = "no samples in " + path;
        return false;
    }
    return true;
}

// .fam: FID IID PID MID SEX PHENO, no header. IID is column 2.
bool PgenSource::loadFam(const std::string& path, std::string& error) {
    std::ifstream in(path);
    if (!in) { error = "cannot open " + path; return false; }
    sample_ids_.clear();
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::istringstream ls(line);
        std::string fid, iid;
        if (!(ls >> fid >> iid)) { error = "malformed .fam record: " + line; return false; }
        sample_ids_.push_back(iid);
    }
    if (sample_ids_.empty()) { error = "no samples in " + path; return false; }
    return true;
}

// .bim: CHROM ID CM POS ALT REF, no header. Note the allele order: plink 1
// writes A1 (usually minor/ALT) before A2 (usually REF), the opposite of a VCF.
bool PgenSource::loadBim(const std::string& path, std::string& error) {
    std::ifstream in(path);
    if (!in) { error = "cannot open " + path; return false; }
    variants_.clear();
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::istringstream ls(line);
        PgenVariantRecord rec;
        std::string cm;
        if (!(ls >> rec.chrom >> rec.id >> cm >> rec.pos >> rec.alt >> rec.ref)) {
            error = "malformed .bim record: " + line;
            return false;
        }
        rec.biallelic = true;  // .bed is biallelic by construction
        variants_.push_back(std::move(rec));
    }
    if (variants_.empty()) { error = "no variants in " + path; return false; }
    return true;
}

bool PgenSource::loadPvar(const std::string& path, std::string& error) {
    std::ifstream in(path);
    if (!in) {
        error = "cannot open " + path;
        return false;
    }
    variants_.clear();
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ls(line);
        PgenVariantRecord rec;
        if (!(ls >> rec.chrom >> rec.pos >> rec.id >> rec.ref >> rec.alt)) {
            error = "malformed record in " + path + ": " + line;
            return false;
        }
        rec.biallelic = (rec.alt.find(',') == std::string::npos);
        variants_.push_back(std::move(rec));
    }
    if (variants_.empty()) {
        error = "no variants in " + path;
        return false;
    }
    return true;
}

bool PgenSource::openFile(std::string& error) {
    const std::string& pgen_path = genotype_path_;
    std::string stem = stripSuffix(pgen_path, ".pgen");
    stem = stripSuffix(stem, ".bed");

    const bool bim_style = (metadata_ == Metadata::Bim);
    const std::string vpath = !variants_path_.empty() ? variants_path_
                                                      : stem + (bim_style ? ".bim" : ".pvar");
    const std::string spath = !samples_path_.empty() ? samples_path_
                                                     : stem + (bim_style ? ".fam" : ".psam");

    if (bim_style) {
        if (!loadFam(spath, error)) return false;
        if (!loadBim(vpath, error)) return false;
    } else {
        if (!loadPsam(spath, error)) return false;
        if (!loadPvar(vpath, error)) return false;
    }

    n_samples_ = static_cast<int>(sample_ids_.size());
    const uint32_t raw_sample_ct = static_cast<uint32_t>(n_samples_);
    const uint32_t raw_variant_ct = static_cast<uint32_t>(variants_.size());

    if (sample_indices_.empty()) {
        sample_indices_.resize(n_samples_);
        for (int i = 0; i < n_samples_; ++i) sample_indices_[i] = i;
    }

    char errbuf[plink2::kPglErrstrBufBlen];
    plink2::PgenHeaderCtrl header_ctrl;
    uintptr_t pgfi_alloc_cacheline_ct = 0;

    plink2::PglErr reterr = plink2::PgfiInitPhase1(
        pgen_path.c_str(), nullptr, raw_variant_ct, raw_sample_ct,
        &header_ctrl, &impl_->pgfi, &pgfi_alloc_cacheline_ct, errbuf);
    if (reterr != plink2::kPglRetSuccess) {
        error = std::string("PgfiInitPhase1: ") + errbuf;
        return false;
    }
    impl_->pgfi_ready = true;

    if (pgfi_alloc_cacheline_ct) {
        if (plink2::cachealigned_malloc(pgfi_alloc_cacheline_ct * plink2::kCacheline,
                                        &impl_->pgfi_alloc)) {
            error = "out of memory allocating pgfi";
            return false;
        }
    }

    uint32_t max_vrec_width = 0;
    uintptr_t pgr_alloc_cacheline_ct = 0;
    reterr = plink2::PgfiInitPhase2(header_ctrl, 1, 1, 0, 0, raw_variant_ct,
                                    &max_vrec_width, &impl_->pgfi, impl_->pgfi_alloc,
                                    &pgr_alloc_cacheline_ct, errbuf);
    if (reterr != plink2::kPglRetSuccess) {
        error = std::string("PgfiInitPhase2: ") + errbuf;
        return false;
    }

    if (pgr_alloc_cacheline_ct) {
        if (plink2::cachealigned_malloc(pgr_alloc_cacheline_ct * plink2::kCacheline,
                                        &impl_->pgr_alloc)) {
            error = "out of memory allocating pgr";
            return false;
        }
    }

    reterr = plink2::PgrInit(pgen_path.c_str(), max_vrec_width, &impl_->pgfi,
                             &impl_->pgr, impl_->pgr_alloc);
    if (reterr != plink2::kPglRetSuccess) {
        error = "PgrInit failed on " + pgen_path;
        return false;
    }
    impl_->pgr_ready = true;

    const uint32_t word_ct = plink2::DivUp(raw_sample_ct, plink2::kBitsPerWordD2);
    impl_->genovec.assign(word_ct, 0);
    impl_->dosage_present.assign(plink2::DivUp(raw_sample_ct, plink2::kBitsPerWord), 0);
    impl_->dosage_main.assign(raw_sample_ct, 0);

    return true;
}

bool PgenSource::readVariant(long long vidx, bool prefer_dosage,
                             std::vector<double>& dosages, DosageStats& stats,
                             long long* filled_from_hardcall) {
    if (vidx < 0 || vidx >= static_cast<long long>(variants_.size())) return false;

    const uint32_t sample_ct = static_cast<uint32_t>(n_samples_);
    plink2::PgrSampleSubsetIndex pssi;
    plink2::PgrSetSampleSubsetIndex(nullptr, &impl_->pgr, &pssi);

    uint32_t dosage_ct = 0;
    plink2::PglErr reterr;
    if (prefer_dosage) {
        reterr = plink2::PgrGetD(nullptr, pssi, sample_ct, static_cast<uint32_t>(vidx),
                                 &impl_->pgr, impl_->genovec.data(),
                                 impl_->dosage_present.data(), impl_->dosage_main.data(),
                                 &dosage_ct);
    } else {
        reterr = plink2::PgrGet(nullptr, pssi, sample_ct, static_cast<uint32_t>(vidx),
                                &impl_->pgr, impl_->genovec.data());
    }
    if (reterr != plink2::kPglRetSuccess) return false;

    // dosage_main is packed over the samples whose dosage_present bit is set, so
    // reaching sample u needs its rank among those. Computing that per sample by
    // popcounting the prefix is quadratic; build the mapping once per variant
    // instead, which matters at UK Biobank sample counts.
    if (prefer_dosage && dosage_ct) {
        dosage_rank_.assign(sample_ct, UINT32_MAX);
        uint32_t rank = 0;
        for (uint32_t u = 0; u < sample_ct; ++u) {
            if (bitAt(impl_->dosage_present.data(), u)) {
                dosage_rank_[u] = rank++;
            }
        }
    }

    dosages.clear();
    dosages.reserve(sample_indices_.size());
    stats = DosageStats{};
    long long filled = 0;

    for (int idx : sample_indices_) {
        const uint32_t u = static_cast<uint32_t>(idx);
        if (u >= sample_ct) { dosages.push_back(kMissingDosage); continue; }

        // An explicit dosage wins where the format says one exists.
        if (prefer_dosage && dosage_ct && dosage_rank_[u] != UINT32_MAX) {
            const double d = static_cast<double>(impl_->dosage_main[dosage_rank_[u]]) * kDosageScale;
            dosages.push_back(d);
            stats.observed_sum += d;
            stats.observed_count++;
            continue;
        }

        // No stored dosage: the hard call is still known unless it is missing.
        const uint32_t g = genoAt(impl_->genovec.data(), u);
        if (g == 3) {
            dosages.push_back(kMissingDosage);
            continue;
        }
        const double d = static_cast<double>(g);
        dosages.push_back(d);
        stats.observed_sum += d;
        stats.observed_count++;
        if (prefer_dosage) filled++;
    }

    if (filled_from_hardcall) *filled_from_hardcall = filled;
    return true;
}

bool PgenSource::initialize(const std::string& sample_list) {
    std::string error;
    if (!openFile(error)) {
        logError("Failed to open " + genotype_path_ + ": " + error);
        return false;
    }

    // Resolve the sample subset against the .psam/.fam IDs.
    sample_indices_.clear();
    target_samples_.clear();
    if (sample_list.empty()) {
        sample_indices_.resize(n_samples_);
        target_samples_ = sample_ids_;
        for (int i = 0; i < n_samples_; ++i) sample_indices_[i] = i;
    } else {
        std::unordered_map<std::string, int> index;
        for (int i = 0; i < n_samples_; ++i) index[sample_ids_[i]] = i;
        for (auto& raw : split(sample_list, ',')) {
            std::string id = raw;
            id.erase(0, id.find_first_not_of(" \t"));
            id.erase(id.find_last_not_of(" \t") + 1);
            auto it = index.find(id);
            if (it != index.end()) {
                target_samples_.push_back(id);
                sample_indices_.push_back(it->second);
            }
        }
        if (target_samples_.empty()) {
            logError("None of the requested samples are present in " + genotype_path_);
            return false;
        }
        logInfo("Sample subset: " + formatCount(target_samples_.size()) + "/" +
                formatCount(n_samples_));
    }

    long long multiallelic = 0;
    for (const auto& v : variants_) if (!v.biallelic) multiallelic++;

    logInfo("Input: " + formatCount(target_samples_.size()) + " samples, " +
            formatCount(variants_.size()) + " variants (" + describe() + ")");
    if (multiallelic > 0) {
        logInfo("  " + formatCount(multiallelic) + " non-biallelic (will be skipped)");
    }

    chooseField();
    return true;
}

// A .bed never stores dosages, and a .pgen records presence explicitly, so the
// choice here needs no sampling heuristic: ask the format.
void PgenSource::chooseField() {
    effective_field_ = dosage_field_;

    PgenScan s = scan(std::min<long long>(5000, static_cast<long long>(variants_.size())));
    if (s.ok) {
        std::ostringstream os;
        os << std::fixed << std::setprecision(1);
        os << "  hard calls for " << (s.hardcall_rate * 100.0) << "% of calls";
        if (s.dosage_presence > 0.0) {
            os << ", dosages for " << (s.dosage_presence * 100.0) << "%";
        }
        logInfo(os.str());
    }

    if (dosage_field_ == DosageField::Auto) {
        // Dosages where the format has them, hard calls elsewhere. There is no
        // ambiguity to resolve: presence is recorded per sample.
        effective_field_ = (s.dosage_presence > 0.0) ? DosageField::DS : DosageField::GT;
    }
    if (effective_field_ == DosageField::DS && s.dosage_presence <= 0.0) {
        if (dosage_field_ == DosageField::DS) {
            logWarning("-dosage-field DS requested but this input stores no dosages; "
                       "reading hard calls.");
        }
        effective_field_ = DosageField::GT;
    }

    logInfo(std::string("Dosage source: ") +
            (effective_field_ == DosageField::GT ? "hard calls" : "dosages, hard calls where absent") +
            (dosage_field_ == DosageField::Auto ? " (auto)" : " (forced)"));
}

void PgenSource::streamVariants(VariantCallback callback) {
    counts_ = SourceCounts{};
    counts_.total = static_cast<long long>(variants_.size());

    std::unordered_set<std::string> extract;
    for (const auto& id : extract_ids_) extract.insert(id);
    std::unordered_set<std::string> seen;

    // Random access plus the full variant table means duplicates can be found up
    // front, so no spool file is needed here: every copy is simply skipped.
    std::unordered_map<std::string, int> locus_count;
    std::unordered_map<std::string, int> id_count;
    for (const auto& v : variants_) {
        if (!v.biallelic) continue;
        locus_count[v.chrom + ":" + std::to_string(v.pos) + ":" + v.ref + ":" + v.alt]++;
        if (v.id != "." && !v.id.empty()) id_count[v.id]++;
    }

    const bool prefer_dosage = (effective_field_ != DosageField::GT);
    std::vector<double> dosages;
    DosageStats stats;

    ProgressBar bar("Reading variants", static_cast<long long>(variants_.size()));
    for (long long v = 0; v < static_cast<long long>(variants_.size()); ++v) {
        bar.increment();
        const PgenVariantRecord& rec = variants_[v];

        if (!rec.biallelic) { counts_.multiallelic++; continue; }
        if (!extract.empty()) {
            if (!extract.contains(rec.id)) { counts_.not_extracted++; continue; }
            seen.insert(rec.id);
        }

        const std::string locus = rec.chrom + ":" + std::to_string(rec.pos) + ":" +
                                  rec.ref + ":" + rec.alt;
        const bool dup = locus_count[locus] > 1 ||
                         (rec.id != "." && !rec.id.empty() && id_count[rec.id] > 1);
        if (dup) { counts_.duplicates++; continue; }

        long long filled = 0;
        if (!readVariant(v, prefer_dosage, dosages, stats, &filled)) continue;

        const double n = static_cast<double>(dosages.size());
        if (n <= 0) continue;
        if (1.0 - static_cast<double>(stats.observed_count) / n > max_missing_rate_) {
            counts_.missing_filtered++;
            continue;
        }
        if (stats.observed_count == 0) { counts_.missing_filtered++; continue; }

        const double af = stats.observed_sum / (2.0 * stats.observed_count);
        if (af < maf_filter_ || af > 1.0 - maf_filter_) { counts_.maf_filtered++; continue; }

        // Mean-impute whatever is left, matching the VCF path.
        const double mean = stats.observed_sum / stats.observed_count;
        for (double& d : dosages) if (d < 0.0) d = mean;

        if (filled > 0) {
            counts_.hardcall_filled += filled;
            counts_.hardcall_filled_variants++;
        }

        auto out = std::make_unique<Variant>();
        out->id = rec.id;
        out->chrom = rec.chrom;
        out->pos = rec.pos;
        out->ref = rec.ref;
        out->alt = rec.alt;
        out->af = af;
        out->dosages = dosages;
        callback(std::move(out));
        counts_.emitted++;
    }
    bar.finish();

    logInfo("Read " + formatCount(counts_.total) + " variants, kept " +
            formatCount(counts_.emitted));
    if (counts_.not_extracted > 0)
        logInfo("  excluded " + formatCount(counts_.not_extracted) + " not in the extract list");
    if (counts_.multiallelic > 0)
        logInfo("  excluded " + formatCount(counts_.multiallelic) + " non-biallelic");
    if (counts_.maf_filtered > 0)
        logInfo("  excluded " + formatCount(counts_.maf_filtered) + " below the MAF threshold");
    if (counts_.missing_filtered > 0)
        logInfo("  excluded " + formatCount(counts_.missing_filtered) + " over the missingness limit");
    if (counts_.duplicates > 0)
        logInfo("  excluded " + formatCount(counts_.duplicates) + " duplicated by locus or ID (every copy)");
    counts_.extract_unmatched = static_cast<long long>(extract.size()) -
                                static_cast<long long>(seen.size());
    if (!extract.empty() && counts_.extract_unmatched > 0) {
        logWarning(formatCount(counts_.extract_unmatched) + " of " + formatCount(extract.size()) +
                   " extract IDs matched no variant. The list must use the same IDs as the "
                   "ID column of the .pvar/.bim.");
    }
    if (counts_.hardcall_filled > 0)
        logInfo("  " + formatCount(counts_.hardcall_filled) + " call(s) taken from hard calls where no dosage was stored");
}

PgenScan PgenSource::scan(long long max_variants) {
    PgenScan s;
    s.n_samples = n_samples_;
    const long long limit = std::min<long long>(max_variants,
                                                static_cast<long long>(variants_.size()));
    if (limit <= 0 || n_samples_ <= 0) return s;

    double dosage_sum = 0.0;
    double hardcall_sum = 0.0;
    long long counted = 0;

    for (long long v = 0; v < limit; ++v) {
        if (!variants_[v].biallelic) { s.multiallelic++; continue; }

        const uint32_t sample_ct = static_cast<uint32_t>(n_samples_);
        plink2::PgrSampleSubsetIndex pssi;
        plink2::PgrSetSampleSubsetIndex(nullptr, &impl_->pgr, &pssi);
        uint32_t dosage_ct = 0;
        if (plink2::PgrGetD(nullptr, pssi, sample_ct, static_cast<uint32_t>(v), &impl_->pgr,
                            impl_->genovec.data(), impl_->dosage_present.data(),
                            impl_->dosage_main.data(), &dosage_ct) != plink2::kPglRetSuccess) {
            continue;
        }

        uint32_t called = 0;
        for (uint32_t u = 0; u < sample_ct; ++u) {
            if (genoAt(impl_->genovec.data(), u) != 3) called++;
        }
        dosage_sum += static_cast<double>(dosage_ct) / sample_ct;
        hardcall_sum += static_cast<double>(called) / sample_ct;
        counted++;
    }

    s.variants = limit;
    if (counted > 0) {
        s.dosage_presence = dosage_sum / static_cast<double>(counted);
        s.hardcall_rate = hardcall_sum / static_cast<double>(counted);
    }
    s.ok = true;
    return s;
}
