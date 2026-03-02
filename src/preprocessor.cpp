#include "preprocessor.h"
#include "vcf_processor.h"
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

    logInfo("Created output directories in: " + config_.output_dir);
}

void Preprocessor::generateAndSaveTraits() {
    Timer timer("Traits matrix generation and serialization");
    logInfo("Generating " + std::to_string(config_.n_traits) + " traits for " +
           std::to_string(n_samples_) + " samples");
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

    logInfo("Traits per tile: " + std::to_string(traits_per_tile) + 
            " (auto, target ~1GB tiles)");

    return traits_per_tile;
}

void Preprocessor::saveTraitsTiled() {
    Timer timer("Saving traits to disk");

    int traits_per_tile = calculateTraitsPerTile();
    int n_tiles = (config_.n_traits + traits_per_tile - 1) / traits_per_tile;

    size_t bytes_per_trait = n_samples_ * sizeof(double);
    size_t tile_size_mb = (traits_per_tile * bytes_per_trait) / (1024 * 1024);

    logInfo("Saving traits matrix in " + std::to_string(n_tiles) + " tiles (" +
           std::to_string(traits_per_tile) + " traits per tile, ~" +
           std::to_string(tile_size_mb) + " MB per tile)");

    TraitsMetadata meta;
    meta.n_traits = config_.n_traits;
    meta.n_samples = n_samples_;
    meta.n_tiles = n_tiles;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> dist(0.0, 1.0);
    std::vector<double> trait_row(n_samples_);
    int generated_traits = 0;

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
            if (generated_traits % 1000 == 0 || generated_traits == config_.n_traits) {
                logDebug("Generated " + std::to_string(generated_traits) + "/" +
                         std::to_string(config_.n_traits) + " traits");
            }
        }

        if (!out) {
            throw std::runtime_error("Failed while writing trait tile file: " + tile_file);
        }

        meta.tile_trait_counts.push_back(traits_in_tile);

        if (n_tiles > 50 && (tile_id + 1) % 10 == 0) {
            logDebug("Saved " + std::to_string(tile_id + 1) + "/" + 
                     std::to_string(n_tiles) + " tiles");
        } else if (n_tiles <= 50) {
            logDebug("Saved tile " + std::to_string(tile_id) + " (" + 
                     std::to_string(traits_in_tile) + " traits)");
        }
    }

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

void Preprocessor::processAndChunkVCF() {
    Timer timer("VCF processing and chunking");

    VCFProcessor processor(config_.vcf_file, config_.maf_filter);
    if (!processor.initialize(config_.sample_list)) {
        throw std::runtime_error("Failed to initialize VCF processor");
    }

    logInfo("Starting streaming VCF processing with chunk size: " + 
           std::to_string(config_.chunk_size));
    logDebug("Memory-efficient mode: processing variants one at a time");

    int chunk_id = 0;
    std::vector<std::vector<double>> current_chunk_genotypes;
    ChunkMetadata current_meta;

    current_chunk_genotypes.reserve(config_.chunk_size);
    current_meta.variant_ids.reserve(config_.chunk_size);
    current_meta.chroms.reserve(config_.chunk_size);
    current_meta.positions.reserve(config_.chunk_size);
    current_meta.refs.reserve(config_.chunk_size);
    current_meta.alts.reserve(config_.chunk_size);

    processor.streamVariants([&](std::unique_ptr<Variant> variant) {
        current_chunk_genotypes.push_back(standardize(variant->dosages));
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
            logDebug("Saved chunk " + std::to_string(chunk_id) + " (" + 
                     std::to_string(current_chunk_genotypes.size()) + " variants)");

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
        logDebug("Saved final chunk " + std::to_string(chunk_id) + " (" + 
                 std::to_string(current_chunk_genotypes.size()) + " variants)");
    }

    logInfo("VCF streaming completed: " + std::to_string(chunk_id + 1) + " chunks created");
}

void Preprocessor::saveChunk(int chunk_id, const std::vector<std::vector<double>>& genotypes,
                             const ChunkMetadata& meta) {
    std::string bin_file = getChunksDir() + "/chunk_" + std::to_string(chunk_id) + ".bin";
    std::ofstream out(bin_file, std::ios::binary);
    if (!out) {
        throw std::runtime_error("Failed to open chunk file: " + bin_file);
    }

    for (const auto& row : genotypes) {
        out.write(reinterpret_cast<const char*>(row.data()), 
                 n_samples_ * sizeof(double));
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

void Preprocessor::run() {
    logInfo("Starting preprocessing...");

    createOutputDirectories();

    VCFProcessor temp_processor(config_.vcf_file, config_.maf_filter);
    if (!temp_processor.initialize(config_.sample_list)) {
        throw std::runtime_error("Failed to initialize VCF processor");
    }
    n_samples_ = temp_processor.getTargetSamples().size();
    saveHeaderMetadata(temp_processor.getContigNames());

    generateAndSaveTraits();

    processAndChunkVCF();

    logInfo("Preprocessing completed successfully!");
    logInfo("Output directory: " + config_.output_dir);
}
