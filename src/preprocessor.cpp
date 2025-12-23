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

Preprocessor::Preprocessor(const PreprocessConfig& config) 
    : config_(config), n_samples_(0) {
}

std::string Preprocessor::getTraitsDir() const {
    return config_.output_dir + "/traits";
}

std::string Preprocessor::getChunksDir() const {
    return config_.output_dir + "/chunks";
}

void Preprocessor::createOutputDirectories() {
    // Create base directory
    if (mkdir(config_.output_dir.c_str(), 0755) != 0 && errno != EEXIST) {
        throw std::runtime_error("Failed to create output directory: " + 
                               std::string(strerror(errno)));
    }

    // Create traits subdirectory
    std::string traits_dir = getTraitsDir();
    if (mkdir(traits_dir.c_str(), 0755) != 0 && errno != EEXIST) {
        throw std::runtime_error("Failed to create traits directory: " + 
                               std::string(strerror(errno)));
    }

    // Create chunks subdirectory
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

    std::vector<std::vector<double>> traits_matrix(config_.n_traits);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> dist(0.0, 1.0);

    for (int i = 0; i < config_.n_traits; ++i) {
        traits_matrix[i].resize(n_samples_);
        for (int j = 0; j < n_samples_; ++j) {
            traits_matrix[i][j] = dist(gen);
        }

        if ((i + 1) % 1000 == 0) {
            logInfo("Generated " + std::to_string(i + 1) + "/" + 
                   std::to_string(config_.n_traits) + " traits");
        }
    }

    saveTraitsTiled(traits_matrix);
}

int Preprocessor::calculateTraitsPerTile() const {
    if (config_.traits_per_tile > 0) {
        // User-specified value
        logInfo("Using user-specified traits per tile: " + 
               std::to_string(config_.traits_per_tile));
        return config_.traits_per_tile;
    }

    // Auto-calculate to aim for ~1GB per tile
    const size_t target_tile_size = 1024 * 1024 * 1024;  // 1 GB
    size_t bytes_per_trait = n_samples_ * sizeof(double);
    int traits_per_tile = std::max(1, static_cast<int>(target_tile_size / bytes_per_trait));

    // But don't exceed total traits
    traits_per_tile = std::min(traits_per_tile, config_.n_traits);

    logInfo("Auto-calculated traits per tile: " + std::to_string(traits_per_tile) + 
           " (aiming for ~1GB tiles)");

    return traits_per_tile;
}

void Preprocessor::saveTraitsTiled(const std::vector<std::vector<double>>& traits_matrix) {
    Timer timer("Saving traits to disk");

    // Calculate traits per tile (auto or user-specified)
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

    for (int tile_id = 0; tile_id < n_tiles; ++tile_id) {
        int start_trait = tile_id * traits_per_tile;
        int end_trait = std::min(start_trait + traits_per_tile, config_.n_traits);
        int traits_in_tile = end_trait - start_trait;

        std::string tile_file = getTraitsDir() + "/W_tile_" + std::to_string(tile_id) + ".bin";
        std::ofstream out(tile_file, std::ios::binary);
        if (!out) {
            throw std::runtime_error("Failed to open trait tile file: " + tile_file);
        }

        // Write traits row by row
        for (int t = start_trait; t < end_trait; ++t) {
            out.write(reinterpret_cast<const char*>(traits_matrix[t].data()),
                     n_samples_ * sizeof(double));
        }

        meta.tile_trait_counts.push_back(traits_in_tile);

        // Log progress for large number of tiles
        if (n_tiles > 50 && (tile_id + 1) % 10 == 0) {
            logInfo("Saved " + std::to_string(tile_id + 1) + "/" + 
                   std::to_string(n_tiles) + " tiles");
        } else if (n_tiles <= 50) {
            logInfo("Saved tile " + std::to_string(tile_id) + " (" + 
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

    logInfo("Saved traits metadata");
}

void Preprocessor::processAndChunkVCF() {
    Timer timer("VCF processing and chunking");

    // Initialize VCF processor
    VCFProcessor processor(config_.vcf_file, config_.maf_filter);
    if (!processor.initialize(config_.sample_list)) {
        throw std::runtime_error("Failed to initialize VCF processor");
    }

    logInfo("Starting streaming VCF processing with chunk size: " + 
           std::to_string(config_.chunk_size));
    logInfo("Memory-efficient mode: processing variants one at a time");

    // Streaming state
    int chunk_id = 0;
    std::vector<std::vector<double>> current_chunk_genotypes;
    ChunkMetadata current_meta;

    // Pre-allocate for one chunk
    current_chunk_genotypes.reserve(config_.chunk_size);
    current_meta.variant_ids.reserve(config_.chunk_size);
    current_meta.chroms.reserve(config_.chunk_size);
    current_meta.positions.reserve(config_.chunk_size);
    current_meta.refs.reserve(config_.chunk_size);
    current_meta.alts.reserve(config_.chunk_size);

    // Stream variants one at a time
    processor.streamVariants([&](std::unique_ptr<Variant> variant) {
        // Standardize and add to current chunk
        current_chunk_genotypes.push_back(standardize(variant->dosages));
        current_meta.variant_ids.push_back(std::move(variant->id));
        current_meta.chroms.push_back(std::move(variant->chrom));
        current_meta.positions.push_back(variant->pos);
        current_meta.refs.push_back(std::move(variant->ref));
        current_meta.alts.push_back(std::move(variant->alt));

        // When chunk is full, save it and clear
        if (current_chunk_genotypes.size() >= static_cast<size_t>(config_.chunk_size)) {
            current_meta.chunk_id = chunk_id;
            current_meta.n_variants = current_chunk_genotypes.size();
            current_meta.n_samples = n_samples_;

            saveChunk(chunk_id, current_chunk_genotypes, current_meta);
            logInfo("Saved chunk " + std::to_string(chunk_id) + " (" + 
                   std::to_string(current_chunk_genotypes.size()) + " variants)");

            // Reset for next chunk
            chunk_id++;
            current_chunk_genotypes.clear();
            current_meta.variant_ids.clear();
            current_meta.chroms.clear();
            current_meta.positions.clear();
            current_meta.refs.clear();
            current_meta.alts.clear();
        }
    });

    // Save final partial chunk if any
    if (!current_chunk_genotypes.empty()) {
        current_meta.chunk_id = chunk_id;
        current_meta.n_variants = current_chunk_genotypes.size();
        current_meta.n_samples = n_samples_;

        saveChunk(chunk_id, current_chunk_genotypes, current_meta);
        logInfo("Saved final chunk " + std::to_string(chunk_id) + " (" + 
               std::to_string(current_chunk_genotypes.size()) + " variants)");
    }

    logInfo("VCF streaming completed: " + std::to_string(chunk_id + 1) + " chunks created");
}

void Preprocessor::saveChunk(int chunk_id, const std::vector<std::vector<double>>& genotypes,
                             const ChunkMetadata& meta) {
    // Save binary genotype data
    std::string bin_file = getChunksDir() + "/chunk_" + std::to_string(chunk_id) + ".bin";
    std::ofstream out(bin_file, std::ios::binary);
    if (!out) {
        throw std::runtime_error("Failed to open chunk file: " + bin_file);
    }

    for (const auto& row : genotypes) {
        out.write(reinterpret_cast<const char*>(row.data()), 
                 n_samples_ * sizeof(double));
    }

    // Save metadata
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

    // Create output structure
    createOutputDirectories();

    // Process VCF to get sample count
    VCFProcessor temp_processor(config_.vcf_file, config_.maf_filter);
    if (!temp_processor.initialize(config_.sample_list)) {
        throw std::runtime_error("Failed to initialize VCF processor");
    }
    n_samples_ = temp_processor.getTargetSamples().size();

    // Generate and save traits
    generateAndSaveTraits();

    // Process and chunk VCF (now using streaming!)
    processAndChunkVCF();

    logInfo("Preprocessing completed successfully!");
    logInfo("Output directory: " + config_.output_dir);
}
