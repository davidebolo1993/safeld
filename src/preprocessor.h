#pragma once

#include <string>
#include <vector>
#include <memory>

struct PreprocessConfig {
    std::string vcf_file;
    std::string output_dir;
    double maf_filter = 0.01;
    int n_traits = 10;
    int chunk_size = 10000;  // variants per chunk
    int traits_per_tile = 0;  // 0 = auto (aim for ~1GB per tile)
    std::string sample_list;
};

struct ChunkMetadata {
    int chunk_id;
    int n_variants;
    int n_samples;
    std::vector<std::string> variant_ids;
    std::vector<std::string> chroms;
    std::vector<int> positions;
    std::vector<std::string> refs;
    std::vector<std::string> alts;
};

struct TraitsMetadata {
    int n_traits;
    int n_samples;
    int n_tiles;
    std::vector<int> tile_trait_counts;  // traits per tile
};

class Preprocessor {
private:
    PreprocessConfig config_;
    int n_samples_;
    
    std::string getTraitsDir() const;
    std::string getChunksDir() const;
    
    void createOutputDirectories();
    void generateAndSaveTraits();
    void processAndChunkVCF();
    
    void saveTraitsTiled(const std::vector<std::vector<double>>& traits_matrix);
    void saveTraitsMetadata(const TraitsMetadata& meta);
    
    void saveChunk(int chunk_id, const std::vector<std::vector<double>>& genotypes,
                   const ChunkMetadata& meta);
    void saveChunkMetadata(const ChunkMetadata& meta);
    
    int calculateTraitsPerTile() const;

public:
    Preprocessor(const PreprocessConfig& config);
    
    void run();
};

