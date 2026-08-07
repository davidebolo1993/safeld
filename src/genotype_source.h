#pragma once

// Common interface over the genotype inputs safeld can preprocess.
//
// Both implementations are responsible for applying the same variant-level
// policy before emitting: biallelic only, MAF and missingness thresholds,
// deduplication on locus (and on ID where one exists), and mean imputation of
// whatever missingness survives. They differ only in how they get there, since
// a VCF is a one-pass stream while a pgen offers random access and its variant
// table up front.

#include <functional>
#include <memory>
#include <string>
#include <vector>

struct Variant;

// Counts reported after streaming, so preprocess can print one summary
// regardless of which source produced the variants.
struct SourceCounts {
    long long total = 0;
    long long emitted = 0;
    long long multiallelic = 0;
    long long missing_filtered = 0;
    long long duplicates = 0;
    long long not_extracted = 0;
    long long maf_filtered = 0;
    // Extract-list IDs that matched nothing in the input.
    long long extract_unmatched = 0;
    // Calls that had no stored dosage and were taken from the hard call.
    long long hardcall_filled = 0;
    long long hardcall_filled_variants = 0;
};

class GenotypeSource {
public:
    using VariantCallback = std::function<void(std::unique_ptr<Variant>)>;

    virtual ~GenotypeSource() = default;

    // Opens the input and resolves the sample set. Returns false after logging.
    virtual bool initialize(const std::string& sample_list) = 0;

    virtual void streamVariants(VariantCallback callback) = 0;

    virtual const std::vector<std::string>& getTargetSamples() const = 0;
    virtual std::vector<std::string> getContigNames() const = 0;
    virtual const SourceCounts& counts() const = 0;

    // Restrict to these variant IDs. Empty means keep everything.
    virtual void setExtractIds(std::vector<std::string> ids) = 0;

    // Human-readable description for the log ("VCF", "pgen", "bed").
    virtual std::string describe() const = 0;
};

// Reads one ID per line, ignoring blanks and '#' comments. Throws on failure.
std::vector<std::string> readIdList(const std::string& path);
