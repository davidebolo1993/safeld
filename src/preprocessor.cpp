#include "preprocessor.h"
#include "vcf_processor.h"
#ifdef SAFELD_HAVE_PGEN
#include "pgen_reader.h"
#endif
#include "utils.h"
#include <fstream>
#include <sstream>
#include <random>
#include <sys/stat.h>
#include <sys/types.h>
#include <cerrno>
#include <cstring>
#include <algorithm>

Preprocessor::Preprocessor(const PreprocessConfig& config) 
    : config_(config), n_samples_(0) {
}

std::string Preprocessor::getTraitsDir() const {
    return config_.output_dir + "/traits";
}

std::string Preprocessor::getChunksDir() const {
    return config_.output_dir + "/chunks";
}

std::string Preprocessor::getHeaderMetaFile() const {
    return config_.output_dir + "/header_contigs.txt";
}

void Preprocessor::saveHeaderMetadata(const std::vector<std::string>& contig_names) {
    std::ofstream out(getHeaderMetaFile());
    if (!out) {
        throw std::runtime_error("Failed to open header metadata file");
    }

    for (const auto& contig : contig_names) {
        out << "##contig=<ID=" << contig << ">\n";
    }
    logDebug("Saved " + std::to_string(contig_names.size()) + " contig header records");
}

void Preprocessor::createOutputDirectories() {
    if (mkdir(config_.output_dir.c_str(), 0755) != 0 && errno != EEXIST) {
        throw std::runtime_error("Failed to create output directory: " + 
                               std::string(strerror(errno)));
    }

    std::string traits_dir = getTraitsDir();
    if (mkdir(traits_dir.c_str(), 0755) != 0 && errno != EEXIST) {
        throw std::runtime_error("Failed to create traits directory: " + 
                               std::string(strerror(errno)));
    }

    std::string chunks_dir = getChunksDir();
    if (mkdir(chunks_dir.c_str(), 0755) != 0 && errno != EEXIST) {
        throw std::runtime_error("Failed to create chunks directory: " + 
                               std::string(strerror(errno)));
    }

    logDebug("Created output directories in " + config_.output_dir);
}

void Preprocessor::generateAndSaveTraits() {
    Timer timer("Traits matrix generation and serialization");
    saveTraitsTiled();
}

int Preprocessor::calculateTraitsPerTile() const {
    if (config_.traits_per_tile > 0) {
        logDebug("Using user-specified traits per tile: " +
                 std::to_string(config_.traits_per_tile));
        return config_.traits_per_tile;
    }

    const size_t target_tile_size = 1024 * 1024 * 1024;
    size_t bytes_per_trait = n_samples_ * sizeof(double);
    int traits_per_tile = std::max(1, static_cast<int>(target_tile_size / bytes_per_trait));

    traits_per_tile = std::min(traits_per_tile, config_.n_traits);

    logDebug("Traits per tile: " + std::to_string(traits_per_tile) +
             " (auto, target ~1GB tiles)");

    return traits_per_tile;
}

void Preprocessor::saveTraitsTiled() {
    Timer timer("Saving traits to disk");

    int traits_per_tile = calculateTraitsPerTile();
    int n_tiles = (config_.n_traits + traits_per_tile - 1) / traits_per_tile;

    size_t bytes_per_trait = n_samples_ * sizeof(double);
    size_t tile_size_mb = (traits_per_tile * bytes_per_trait) / (1024 * 1024);

    logInfo("Trait matrix: " + formatCount(config_.n_traits) + " x " +
            formatCount(n_samples_) + " in " + std::to_string(n_tiles) + " tile(s) of " +
            formatBytes(static_cast<unsigned long long>(traits_per_tile) * bytes_per_trait));
    (void)tile_size_mb;

    TraitsMetadata meta;
    meta.n_traits = config_.n_traits;
    meta.n_samples = n_samples_;
    meta.n_tiles = n_tiles;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> dist(0.0, 1.0);
    std::vector<double> trait_row(n_samples_);
    int generated_traits = 0;
    ProgressBar traits_bar("Generating traits", config_.n_traits);

    for (int tile_id = 0; tile_id < n_tiles; ++tile_id) {
        int start_trait = tile_id * traits_per_tile;
        int end_trait = std::min(start_trait + traits_per_tile, config_.n_traits);
        int traits_in_tile = end_trait - start_trait;

        std::string tile_file = getTraitsDir() + "/W_tile_" + std::to_string(tile_id) + ".bin";
        std::ofstream out(tile_file, std::ios::binary);
        if (!out) {
            throw std::runtime_error("Failed to open trait tile file: " + tile_file);
        }

        for (int t = 0; t < traits_in_tile; ++t) {
            for (int j = 0; j < n_samples_; ++j) {
                trait_row[j] = dist(gen);
            }
            out.write(reinterpret_cast<const char*>(trait_row.data()),
                      static_cast<std::streamsize>(n_samples_ * sizeof(double)));

            generated_traits++;
            traits_bar.update(generated_traits);
        }

        if (!out) {
            throw std::runtime_error("Failed while writing trait tile file: " + tile_file);
        }

        meta.tile_trait_counts.push_back(traits_in_tile);

        logDebug("Wrote tile " + std::to_string(tile_id) + " (" +
                 std::to_string(traits_in_tile) + " traits)");
    }

    traits_bar.finish();
    saveTraitsMetadata(meta);
}

void Preprocessor::saveTraitsMetadata(const TraitsMetadata& meta) {
    std::string meta_file = getTraitsDir() + "/metadata.txt";
    std::ofstream out(meta_file);
    if (!out) {
        throw std::runtime_error("Failed to open traits metadata file");
    }

    out << "n_traits=" << meta.n_traits << "\n";
    out << "n_samples=" << meta.n_samples << "\n";
    out << "n_tiles=" << meta.n_tiles << "\n";
    out << "tile_trait_counts=";
    for (size_t i = 0; i < meta.tile_trait_counts.size(); ++i) {
        if (i > 0) out << ",";
        out << meta.tile_trait_counts[i];
    }
    out << "\n";

    logDebug("Saved traits metadata");
}

void Preprocessor::processAndChunkVCF(GenotypeSource& source) {
    Timer timer("VCF processing and chunking");

    logDebug("Chunk size: " + formatCount(config_.chunk_size) + " variants");

    int chunk_id = 0;
    int chunks_written = 0;
    int monomorphic_variants = 0;
    std::vector<std::vector<double>> current_chunk_genotypes;
    ChunkMetadata current_meta;

    current_chunk_genotypes.reserve(config_.chunk_size);
    current_meta.variant_ids.reserve(config_.chunk_size);
    current_meta.chroms.reserve(config_.chunk_size);
    current_meta.positions.reserve(config_.chunk_size);
    current_meta.refs.reserve(config_.chunk_size);
    current_meta.alts.reserve(config_.chunk_size);

    source.streamVariants([&](std::unique_ptr<Variant> variant) {
        std::vector<double> standardized;
        if (!standardize(variant->dosages, standardized)) {
            // Zero variance across the selected samples: the variant carries no
            // LD information, and a constant row would be emitted downstream as a
            // synthetic dosage of 1.0 for every trait.
            monomorphic_variants++;
            return;
        }

        current_chunk_genotypes.push_back(std::move(standardized));
        current_meta.variant_ids.push_back(std::move(variant->id));
        current_meta.chroms.push_back(std::move(variant->chrom));
        current_meta.positions.push_back(variant->pos);
        current_meta.refs.push_back(std::move(variant->ref));
        current_meta.alts.push_back(std::move(variant->alt));

        if (current_chunk_genotypes.size() >= static_cast<size_t>(config_.chunk_size)) {
            current_meta.chunk_id = chunk_id;
            current_meta.n_variants = current_chunk_genotypes.size();
            current_meta.n_samples = n_samples_;

            saveChunk(chunk_id, current_chunk_genotypes, current_meta);
            chunks_written++;
            logDebug("Wrote chunk " + std::to_string(chunk_id) + " (" +
                     formatCount(current_chunk_genotypes.size()) + " variants)");

            chunk_id++;
            current_chunk_genotypes.clear();
            current_meta.variant_ids.clear();
            current_meta.chroms.clear();
            current_meta.positions.clear();
            current_meta.refs.clear();
            current_meta.alts.clear();
        }
    });

    if (!current_chunk_genotypes.empty()) {
        current_meta.chunk_id = chunk_id;
        current_meta.n_variants = current_chunk_genotypes.size();
        current_meta.n_samples = n_samples_;

        saveChunk(chunk_id, current_chunk_genotypes, current_meta);
        chunks_written++;
        logDebug("Wrote final chunk " + std::to_string(chunk_id) + " (" +
                 formatCount(current_chunk_genotypes.size()) + " variants)");
    }

    if (monomorphic_variants > 0) {
        logInfo("  excluded " + formatCount(monomorphic_variants) +
                " with no variance across the selected samples");
    }

    logInfo("Wrote " + std::to_string(chunks_written) + " chunk(s)");
    if (chunks_written == 0) {
        logWarning("No chunks were created; check the MAF, missingness and "
                   "deduplication counts above");
    }
}

void Preprocessor::saveChunk(int chunk_id, const std::vector<std::vector<double>>& genotypes,
                             const ChunkMetadata& meta) {
    std::string bin_file = getChunksDir() + "/chunk_" + std::to_string(chunk_id) + ".bin";
    std::ofstream out(bin_file, std::ios::binary);
    if (!out) {
        throw std::runtime_error("Failed to open chunk file: " + bin_file);
    }

    for (const auto& row : genotypes) {
        if (row.size() != static_cast<size_t>(n_samples_)) {
            throw std::runtime_error("Genotype row length " + std::to_string(row.size()) +
                                     " does not match sample count " + std::to_string(n_samples_));
        }
        out.write(reinterpret_cast<const char*>(row.data()),
                 n_samples_ * sizeof(double));
    }

    out.flush();
    if (!out) {
        throw std::runtime_error("Failed while writing chunk file: " + bin_file);
    }

    saveChunkMetadata(meta);
}

void Preprocessor::saveChunkMetadata(const ChunkMetadata& meta) {
    std::string meta_file = getChunksDir() + "/chunk_" + 
                           std::to_string(meta.chunk_id) + ".meta";
    std::ofstream out(meta_file);
    if (!out) {
        throw std::runtime_error("Failed to open chunk metadata file");
    }

    out << "chunk_id=" << meta.chunk_id << "\n";
    out << "n_variants=" << meta.n_variants << "\n";
    out << "n_samples=" << meta.n_samples << "\n";

    for (int i = 0; i < meta.n_variants; ++i) {
        out << meta.chroms[i] << "\t" 
            << meta.positions[i] << "\t"
            << meta.variant_ids[i] << "\t"
            << meta.refs[i] << "\t"
            << meta.alts[i] << "\n";
    }
}

// Builds the reader for whichever input was given. The deduplication spool used
// by the VCF path lives next to the output, since it holds a full copy of the
// genotype matrix; the pgen path needs no spool because it can see the whole
// variant table up front.
std::unique_ptr<GenotypeSource> Preprocessor::makeSource() {
    if (!config_.genotype_file.empty()) {
#ifdef SAFELD_HAVE_PGEN
        auto src = std::make_unique<PgenSource>(
            config_.genotype_file,
            config_.plink1_metadata ? PgenSource::Metadata::Bim : PgenSource::Metadata::Pvar,
            config_.maf_filter, config_.max_missing_rate, config_.dosage_field);
        if (!config_.variants_file.empty() || !config_.samples_file.empty()) {
            src->setMetadataPaths(config_.variants_file, config_.samples_file);
        }
        return src;
#else
        throw std::runtime_error(
            "This build has no .pgen/.bed support. Rebuild with "
            "-DSAFELD_PGEN=ON -DPLINK_NG_DIR=/path/to/plink-ng, or convert to VCF first.");
#endif
    }
    return std::make_unique<VCFProcessor>(config_.vcf_file, config_.maf_filter,
                                          config_.max_missing_rate, config_.output_dir,
                                          config_.dosage_field, config_.use_info_af);
}

void Preprocessor::run() {
    LogModule module("preprocess");
    Timer timer("preprocess stage");

    createOutputDirectories();

    std::unique_ptr<GenotypeSource> source = makeSource();

    // Before initialize(), so a source with random access can skip records it
    // was never asked for instead of loading and indexing the whole file.
    if (!config_.extract_file.empty()) {
        auto ids = readIdList(config_.extract_file);
        logInfo("Extract list: " + formatCount(ids.size()) + " variant IDs");
        source->setExtractIds(std::move(ids));
    }

    if (!source->initialize(config_.sample_list)) {
        throw std::runtime_error("Failed to open the genotype input");
    }

    n_samples_ = static_cast<int>(source->getTargetSamples().size());
    saveHeaderMetadata(source->getContigNames());

    generateAndSaveTraits();

    processAndChunkVCF(*source);

    logInfo("Done in " + formatDuration(timer.elapsed()) + " -> " + config_.output_dir);
}
