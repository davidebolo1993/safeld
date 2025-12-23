#pragma once

#include <string>
#include <vector>

struct MergerConfig {
    std::string input_dir;
    std::string output_file;
    bool compress_output = true;
};

class VCFMerger {
private:
    MergerConfig config_;
    
    std::vector<std::string> findChunkFiles();
    void mergeChunks(const std::vector<std::string>& chunk_files);

public:
    VCFMerger(const MergerConfig& config);
    
    void run();
};

