#include "vcf_merger.h"
#include "utils.h"
#include <fstream>
#include <algorithm>
#include <htslib/bgzf.h>
#include <htslib/kstring.h>
#include <htslib/tbx.h>
#include <dirent.h>
#include <cstring>
#include <climits>
#include <unordered_map>
#include <string_view>
#include <cstdlib>
#include <cstdio>

namespace {
std::string parseContigId(std::string_view line) {
    constexpr std::string_view prefix = "##contig=<ID=";
    if (!line.starts_with(prefix)) {
        return "";
    }
    size_t start = prefix.size();
    size_t end = line.find_first_of(",>", start);
    if (end == std::string_view::npos || end <= start) {
        return "";
    }
    return std::string(line.substr(start, end - start));
}

bool parseChromPos(const char* line, size_t len, std::string& chrom, int& pos) {
    const char* tab1 = static_cast<const char*>(memchr(line, '\t', len));
    if (!tab1) {
        return false;
    }
    const char* pos_start = tab1 + 1;
    size_t remaining = len - static_cast<size_t>(pos_start - line);
    const char* tab2 = static_cast<const char*>(memchr(pos_start, '\t', remaining));
    if (!tab2) {
        return false;
    }

    chrom.assign(line, static_cast<size_t>(tab1 - line));
    char* end_ptr = nullptr;
    long parsed_pos = std::strtol(pos_start, &end_ptr, 10);
    if (end_ptr != tab2 || parsed_pos < 0 || parsed_pos > INT_MAX) {
        return false;
    }
    pos = static_cast<int>(parsed_pos);
    return true;
}

std::string shellEscape(const std::string& value) {
    std::string escaped = "'";
    for (char c : value) {
        if (c == '\'') {
            escaped += "'\\''";
        } else {
            escaped.push_back(c);
        }
    }
    escaped.push_back('\'');
    return escaped;
}

bool commandExists(const std::string& command) {
    std::string probe = "command -v " + command + " >/dev/null 2>&1";
    return std::system(probe.c_str()) == 0;
}

void buildTabixIndex(const std::string& vcf_gz_path) {
    if (tbx_index_build(vcf_gz_path.c_str(), 0, &tbx_conf_vcf) != 0) {
        throw std::runtime_error("Failed to build tabix index for: " + vcf_gz_path);
    }
}
}  // namespace

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
    
    // Sort numerically by chunk id to avoid lexical ordering issues
    // (e.g., chunk_10 before chunk_2).
    auto chunkId = [](const std::string& path) {
        size_t slash = path.find_last_of('/');
        std::string filename = (slash == std::string::npos) ? path : path.substr(slash + 1);
        if (filename.rfind("chunk_", 0) != 0) {
            return INT_MAX;
        }

        size_t start = std::string("chunk_").size();
        size_t end = filename.find('.', start);
        std::string id_str = filename.substr(start, end - start);

        try {
            return std::stoi(id_str);
        } catch (...) {
            return INT_MAX;
        }
    };

    std::sort(files.begin(), files.end(), [&](const std::string& a, const std::string& b) {
        int id_a = chunkId(a);
        int id_b = chunkId(b);
        if (id_a != id_b) {
            return id_a < id_b;
        }
        return a < b;
    });
    
    return files;
}

bool VCFMerger::mergeChunks(const std::vector<std::string>& chunk_files) {
    Timer timer("Merging VCF chunks");
    
    if (chunk_files.empty()) {
        throw std::runtime_error("No chunk files found");
    }
    
    logDebug("Merging " + std::to_string(chunk_files.size()) + " chunk files using fast block I/O");
    
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
    bool is_sorted = true;
    bool has_prev_coord = false;
    std::string prev_chrom;
    int prev_pos = 0;
    std::unordered_map<std::string, int> contig_rank;
    int next_contig_rank = 0;
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
                        std::string contig_id = parseContigId(std::string_view(line.s, line.l));
                        if (!contig_id.empty() && !contig_rank.contains(contig_id)) {
                            contig_rank[contig_id] = next_contig_rank++;
                        }
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

                if (config_.enforce_sort) {
                    std::string chrom;
                    int pos = 0;
                    if (parseChromPos(line.s, line.l, chrom, pos)) {
                        if (has_prev_coord) {
                            auto rankOf = [&](const std::string& contig) {
                                auto it = contig_rank.find(contig);
                                return it == contig_rank.end() ? INT_MAX : it->second;
                            };
                            int prev_rank = rankOf(prev_chrom);
                            int curr_rank = rankOf(chrom);

                            bool out_of_order = false;
                            if (prev_rank != INT_MAX && curr_rank != INT_MAX) {
                                out_of_order = (curr_rank < prev_rank) ||
                                               (curr_rank == prev_rank && pos < prev_pos);
                            } else if (chrom != prev_chrom) {
                                out_of_order = chrom < prev_chrom;
                            } else {
                                out_of_order = pos < prev_pos;
                            }
                            if (out_of_order) {
                                is_sorted = false;
                            }
                        }
                        prev_chrom = std::move(chrom);
                        prev_pos = pos;
                        has_prev_coord = true;
                    }
                }
                
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
                    logDebug("  Processed " + std::to_string(line_count) + " lines...");
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
                        std::string contig_id = parseContigId(str_line);
                        if (!contig_id.empty() && !contig_rank.contains(contig_id)) {
                            contig_rank[contig_id] = next_contig_rank++;
                        }
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

                if (config_.enforce_sort) {
                    std::string chrom;
                    int pos = 0;
                    if (parseChromPos(str_line.c_str(), str_line.size(), chrom, pos)) {
                        if (has_prev_coord) {
                            auto rankOf = [&](const std::string& contig) {
                                auto it = contig_rank.find(contig);
                                return it == contig_rank.end() ? INT_MAX : it->second;
                            };
                            int prev_rank = rankOf(prev_chrom);
                            int curr_rank = rankOf(chrom);

                            bool out_of_order = false;
                            if (prev_rank != INT_MAX && curr_rank != INT_MAX) {
                                out_of_order = (curr_rank < prev_rank) ||
                                               (curr_rank == prev_rank && pos < prev_pos);
                            } else if (chrom != prev_chrom) {
                                out_of_order = chrom < prev_chrom;
                            } else {
                                out_of_order = pos < prev_pos;
                            }
                            if (out_of_order) {
                                is_sorted = false;
                            }
                        }
                        prev_chrom = std::move(chrom);
                        prev_pos = pos;
                        has_prev_coord = true;
                    }
                }
                
                // Write variant line
                if (out_fp) {
                    bgzf_write(out_fp, str_line.c_str(), str_line.length());
                    bgzf_write(out_fp, "\n", 1);
                } else {
                    out_file << str_line << "\n";
                }
                
                if (line_count % 50000 == 0) {
                    logDebug("  Processed " + std::to_string(line_count) + " lines...");
                }
            }
        }
    }
    
    ks_free(&line);
    
    if (out_fp) {
        bgzf_close(out_fp);
    }
    
    logInfo("Merged VCF written to: " + config_.output_file);
    return is_sorted;
}

void VCFMerger::finalizeOutput(bool already_sorted) {
    bool want_index = config_.write_index;
    if (!config_.compress_output && want_index) {
        logWarning("Indexing is only supported for compressed output; skipping index generation");
        want_index = false;
    }

    bool need_sort = config_.enforce_sort && !already_sorted;
    bool need_index = want_index && config_.compress_output;
    if (!need_sort && !need_index) {
        return;
    }

    if (need_sort) {
        if (!commandExists("bcftools")) {
            throw std::runtime_error(
                "Merged output is unsorted and bcftools was not found in PATH for fallback sorting");
        }

        std::string tmp_sorted = config_.output_file + ".tmp.sorted" +
                                 (config_.compress_output ? ".vcf.gz" : ".vcf");
        std::string cmd;
        if (config_.compress_output) {
            cmd = "bcftools sort -Oz ";
            cmd += "-o " + shellEscape(tmp_sorted) + " " + shellEscape(config_.output_file);
        } else {
            cmd = "bcftools sort -Ov -o " + shellEscape(tmp_sorted) + " " +
                  shellEscape(config_.output_file);
        }

        logInfo("Detected unsorted merged output; sorting with bcftools");
        logDebug("Running: " + cmd);
        if (std::system(cmd.c_str()) != 0) {
            throw std::runtime_error("Failed to sort merged VCF with bcftools");
        }

        if (std::remove(config_.output_file.c_str()) != 0) {
            throw std::runtime_error("Failed to replace merged output after sorting");
        }
        if (std::rename(tmp_sorted.c_str(), config_.output_file.c_str()) != 0) {
            throw std::runtime_error("Failed to move sorted output into final location");
        }
    }

    if (need_index) {
        logInfo("Indexing merged output with HTSlib tabix");
        buildTabixIndex(config_.output_file);
    }
}

void VCFMerger::run() {
    logInfo("Starting VCF merge...");
    
    auto chunk_files = findChunkFiles();
    
    if (chunk_files.empty()) {
        throw std::runtime_error("No chunk files found in: " + config_.input_dir);
    }
    
    logInfo("Found " + std::to_string(chunk_files.size()) + " chunk files");
    
    bool is_sorted = mergeChunks(chunk_files);
    if (config_.enforce_sort && !is_sorted) {
        logWarning("Merged records are not coordinate-sorted; sorting will be applied");
    }
    finalizeOutput(is_sorted);
    
    logInfo("VCF merge completed!");
}
