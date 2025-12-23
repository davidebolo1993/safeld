#pragma once

#include <string>
#include <vector>
#include <memory>
#include "preprocessor.h"

struct SimulationConfig {
    std::string preprocessed_dir;
    std::string output_dir;
    int n_workers = 0;  // 0 = auto-detect
    bool compress_output = false;
    int start_chunk = -1;  // -1 = process all chunks
    int end_chunk = -1;
};

class ChunkedSimulator {
private:
    SimulationConfig config_;
    TraitsMetadata traits_meta_;
    std::vector<std::vector<double>> traits_matrix_;  // Full or tiled loading
    
    void loadTraitsMetadata();
    void loadTraitsMatrix();
    void loadTraitsTile(int tile_id, int start_trait, int n_traits);
    
    ChunkMetadata loadChunkMetadata(int chunk_id);
    std::vector<std::vector<double>> loadChunkGenotypes(int chunk_id, int n_variants);
    
    void simulateChunk(int chunk_id);
    void writeChunkVCF(int chunk_id, const ChunkMetadata& meta,
                      const std::vector<std::vector<double>>& synthetic_dosages);

public:
    ChunkedSimulator(const SimulationConfig& config);
    
    void run();
};

