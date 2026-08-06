#include "pgen_reader.h"

#include "utils.h"

#include <cstring>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <algorithm>

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

PgenSource::PgenSource() : impl_(new Impl()) {}

PgenSource::~PgenSource() { delete impl_; }

void PgenSource::setSampleSubset(const std::vector<int>& indices) {
    sample_indices_ = indices;
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

bool PgenSource::open(const std::string& pgen_path, std::string& error,
                      const std::string& pvar_path, const std::string& psam_path) {
    const std::string stem = stripSuffix(pgen_path, ".pgen");
    const std::string pvar = pvar_path.empty() ? stem + ".pvar" : pvar_path;
    const std::string psam = psam_path.empty() ? stem + ".psam" : psam_path;

    if (!loadPsam(psam, error)) return false;
    if (!loadPvar(pvar, error)) return false;

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
