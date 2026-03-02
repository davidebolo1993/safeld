#pragma once

#include <string>
#include <vector>

struct MergerConfig {
    std::string input_dir;
    std::string output_file;
    bool compress_output = true;
    bool write_index = true;
    bool enforce_sort = true;
};

class VCFMerger {
private:
    MergerConfig config_;
    
    std::vector<std::string> findChunkFiles();
    bool mergeChunks(const std::vector<std::string>& chunk_files);
    void finalizeOutput(bool already_sorted);

public:
    VCFMerger(const MergerConfig& config);
    
    void run();
};
