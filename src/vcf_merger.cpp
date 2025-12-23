#include "vcf_merger.h"
#include "utils.h"
#include <fstream>
#include <algorithm>
#include <htslib/bgzf.h>
#include <htslib/kstring.h>
#include <dirent.h>
#include <cstring>

VCFMerger::VCFMerger(const MergerConfig& config) : config_(config) {
}

std::vector<std::string> VCFMerger::findChunkFiles() {
    std::vector<std::string> files;
    
    DIR* dir = opendir(config_.input_dir.c_str());
    if (!dir) {
        throw std::runtime_error("Failed to open input directory: " + config_.input_dir);
    }
    
    struct dirent* entry;
    while ((entry = readdir(dir)) != nullptr) {
        std::string filename = entry->d_name;
        if (filename.find("chunk_") == 0 && filename.find(".vcf") != std::string::npos) {
            files.push_back(config_.input_dir + "/" + filename);
        }
    }
    closedir(dir);
    
    // Sort by chunk number
    std::sort(files.begin(), files.end());
    
    return files;
}

void VCFMerger::mergeChunks(const std::vector<std::string>& chunk_files) {
    Timer timer("Merging VCF chunks");
    
    if (chunk_files.empty()) {
        throw std::runtime_error("No chunk files found");
    }
    
    logInfo("Merging " + std::to_string(chunk_files.size()) + " chunk files using fast block I/O");
    
    // Open output file
    BGZF* out_fp = nullptr;
    std::ofstream out_file;
    
    if (config_.compress_output) {
        out_fp = bgzf_open(config_.output_file.c_str(), "w6");
        if (!out_fp) {
            throw std::runtime_error("Failed to open output file: " + config_.output_file);
        }
    } else {
        out_file.open(config_.output_file);
        if (!out_file) {
            throw std::runtime_error("Failed to open output file: " + config_.output_file);
        }
    }
    
    bool header_written = false;
    kstring_t line = KS_INITIALIZE;
    
    for (size_t chunk_idx = 0; chunk_idx < chunk_files.size(); ++chunk_idx) {
        const auto& chunk_file = chunk_files[chunk_idx];
        logInfo("Processing chunk " + std::to_string(chunk_idx + 1) + "/" + 
               std::to_string(chunk_files.size()) + ": " + chunk_file);
        
        BGZF* in_fp = nullptr;
        std::ifstream in_file;
        bool is_compressed = chunk_file.find(".gz") != std::string::npos;
        
        if (is_compressed) {
            in_fp = bgzf_open(chunk_file.c_str(), "r");
            if (!in_fp) {
                logWarning("Failed to open chunk file: " + chunk_file);
                continue;
            }
            
            // Fast path for compressed files: use kstring with bgzf_getline
            int line_count = 0;
            while (bgzf_getline(in_fp, '\n', &line) >= 0) {
                line_count++;
                
                // Skip empty lines
                if (line.l == 0) continue;
                
                // Handle header
                if (line.s[0] == '#') {
                    if (!header_written) {
                        if (out_fp) {
                            bgzf_write(out_fp, line.s, line.l);
                            bgzf_write(out_fp, "\n", 1);
                        } else {
                            out_file.write(line.s, line.l);
                            out_file << "\n";
                        }
                    }
                    continue;
                }
                
                header_written = true;
                
                // Write variant line
                if (out_fp) {
                    bgzf_write(out_fp, line.s, line.l);
                    bgzf_write(out_fp, "\n", 1);
                } else {
                    out_file.write(line.s, line.l);
                    out_file << "\n";
                }
                
                // Progress update for large files
                if (line_count % 50000 == 0) {
                    logInfo("  Processed " + std::to_string(line_count) + " lines...");
                }
            }
            
            bgzf_close(in_fp);
            
        } else {
            // Uncompressed files: use std::getline
            in_file.open(chunk_file);
            if (!in_file) {
                logWarning("Failed to open chunk file: " + chunk_file);
                continue;
            }
            
            std::string str_line;
            str_line.reserve(10000);
            int line_count = 0;
            
            while (std::getline(in_file, str_line)) {
                line_count++;
                
                if (str_line.empty()) continue;
                
                // Handle header
                if (str_line[0] == '#') {
                    if (!header_written) {
                        if (out_fp) {
                            bgzf_write(out_fp, str_line.c_str(), str_line.length());
                            bgzf_write(out_fp, "\n", 1);
                        } else {
                            out_file << str_line << "\n";
                        }
                    }
                    continue;
                }
                
                header_written = true;
                
                // Write variant line
                if (out_fp) {
                    bgzf_write(out_fp, str_line.c_str(), str_line.length());
                    bgzf_write(out_fp, "\n", 1);
                } else {
                    out_file << str_line << "\n";
                }
                
                if (line_count % 50000 == 0) {
                    logInfo("  Processed " + std::to_string(line_count) + " lines...");
                }
            }
        }
    }
    
    ks_free(&line);
    
    if (out_fp) {
        bgzf_close(out_fp);
    }
    
    logInfo("Merged VCF written to: " + config_.output_file);
}

void VCFMerger::run() {
    logInfo("Starting VCF merge...");
    
    auto chunk_files = findChunkFiles();
    
    if (chunk_files.empty()) {
        throw std::runtime_error("No chunk files found in: " + config_.input_dir);
    }
    
    logInfo("Found " + std::to_string(chunk_files.size()) + " chunk files");
    
    mergeChunks(chunk_files);
    
    logInfo("VCF merge completed!");
}

