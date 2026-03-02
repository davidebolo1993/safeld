#include "chunked_simulator.h"
#include "utils.h"

#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <htslib/bgzf.h>
#include <sys/stat.h>
#include <cerrno>

extern "C" {
    #include <cblas.h>
}

ChunkedSimulator::ChunkedSimulator(const SimulationConfig& config)
    : config_(config) {
}

void ChunkedSimulator::loadTraitsMetadata() {
    Timer timer("Loading traits metadata");
    traits_meta_ = TraitsMetadata{};
    std::string meta_file = config_.preprocessed_dir + "/traits/metadata.txt";
    std::ifstream in(meta_file);
    if (!in) {
        throw std::runtime_error("Failed to open traits metadata: " + meta_file);
    }

    std::string line;
    while (std::getline(in, line)) {
        auto pos = line.find('=');
        if (pos == std::string::npos) continue;
        std::string key = line.substr(0, pos);
        std::string value = line.substr(pos + 1);
        if (key == "n_traits") {
            traits_meta_.n_traits = std::stoi(value);
        } else if (key == "n_samples") {
            traits_meta_.n_samples = std::stoi(value);
        } else if (key == "n_tiles") {
            traits_meta_.n_tiles = std::stoi(value);
        } else if (key == "tile_trait_counts") {
            std::istringstream ss(value);
            std::string token;
            while (std::getline(ss, token, ',')) {
                traits_meta_.tile_trait_counts.push_back(std::stoi(token));
            }
        }
    }

    if (traits_meta_.n_tiles <= 0 || traits_meta_.n_traits <= 0 || traits_meta_.n_samples <= 0) {
        throw std::runtime_error("Invalid traits metadata values in: " + meta_file);
    }
    if (traits_meta_.tile_trait_counts.size() != static_cast<size_t>(traits_meta_.n_tiles)) {
        throw std::runtime_error("Trait metadata mismatch: tile count does not match n_tiles");
    }
    int total_tile_traits = 0;
    for (int count : traits_meta_.tile_trait_counts) {
        total_tile_traits += count;
    }
    if (total_tile_traits != traits_meta_.n_traits) {
        throw std::runtime_error("Trait metadata mismatch: sum(tile_trait_counts) != n_traits");
    }

    logInfo("Traits metadata: " + std::to_string(traits_meta_.n_traits) + " traits, " +
            std::to_string(traits_meta_.n_samples) + " samples, " +
            std::to_string(traits_meta_.n_tiles) + " tiles");
}

void ChunkedSimulator::loadHeaderMetadata() {
    contig_header_lines_.clear();
    std::string header_file = config_.preprocessed_dir + "/header_contigs.txt";
    std::ifstream in(header_file);
    if (!in) {
        logWarning("Contig header metadata not found: " + header_file);
        return;
    }

    std::string line;
    while (std::getline(in, line)) {
        if (line.rfind("##contig=<ID=", 0) == 0) {
            contig_header_lines_.push_back(line);
        }
    }

    logInfo("Loaded " + std::to_string(contig_header_lines_.size()) + " contig header records");
}

std::vector<double> ChunkedSimulator::loadTraitsTileData(int tile_id, int n_traits) {
    std::string tile_file = config_.preprocessed_dir + "/traits/W_tile_" +
                            std::to_string(tile_id) + ".bin";
    std::ifstream in(tile_file, std::ios::binary);
    if (!in) {
        throw std::runtime_error("Failed to open trait tile: " + tile_file);
    }

    std::vector<double> tile_data(static_cast<size_t>(n_traits) * traits_meta_.n_samples);
    if (!tile_data.empty()) {
        in.read(reinterpret_cast<char*>(tile_data.data()),
                static_cast<std::streamsize>(tile_data.size() * sizeof(double)));
    }

    if (!in) {
        throw std::runtime_error("Failed while reading trait tile: " + tile_file);
    }

    return tile_data;
}

ChunkMetadata ChunkedSimulator::loadChunkMetadata(int chunk_id) {
    std::string meta_file = config_.preprocessed_dir + "/chunks/chunk_" +
                           std::to_string(chunk_id) + ".meta";
    std::ifstream in(meta_file);
    if (!in) {
        throw std::runtime_error("Failed to open chunk metadata: " + meta_file);
    }

    ChunkMetadata meta;
    meta.chunk_id = chunk_id;
    std::string line;
    // Read header lines
    while (std::getline(in, line) && line.find('=') != std::string::npos) {
        auto pos = line.find('=');
        std::string key = line.substr(0, pos);
        std::string value = line.substr(pos + 1);
        if (key == "n_variants") meta.n_variants = std::stoi(value);
        else if (key == "n_samples") meta.n_samples = std::stoi(value);
    }

    // Read variant records (already read first line in loop)
    do {
        std::istringstream ss(line);
        std::string chrom, id, ref, alt;
        int pos;
        ss >> chrom >> pos >> id >> ref >> alt;
        meta.chroms.push_back(chrom);
        meta.positions.push_back(pos);
        meta.variant_ids.push_back(id);
        meta.refs.push_back(ref);
        meta.alts.push_back(alt);
    } while (std::getline(in, line));
    return meta;
}

std::vector<std::vector<double>> ChunkedSimulator::loadChunkGenotypes(int chunk_id, int n_variants) {
    std::string bin_file = config_.preprocessed_dir + "/chunks/chunk_" +
                          std::to_string(chunk_id) + ".bin";
    std::ifstream in(bin_file, std::ios::binary);
    if (!in) {
        throw std::runtime_error("Failed to open chunk genotypes: " + bin_file);
    }

    std::vector<std::vector<double>> genotypes(n_variants);
    for (int v = 0; v < n_variants; ++v) {
        genotypes[v].resize(traits_meta_.n_samples);
        in.read(reinterpret_cast<char*>(genotypes[v].data()),
                traits_meta_.n_samples * sizeof(double));
    }
    return genotypes;
}

void ChunkedSimulator::simulateChunk(int chunk_id) {
    logInfo("Starting simulation for chunk " + std::to_string(chunk_id));
    
    // Load chunk data
    auto meta = loadChunkMetadata(chunk_id);
    auto genotypes = loadChunkGenotypes(chunk_id, meta.n_variants);
    
    logInfo("Processing chunk " + std::to_string(chunk_id) + ": " +
            std::to_string(meta.n_variants) + " variants");

    // Process variants in batches to avoid exceeding vector max_size
    static constexpr int VARIANT_BATCH_SIZE = 4000;  // Process 4000 variants at a time
    std::vector<std::vector<double>> all_synthetic_dosages;
    all_synthetic_dosages.reserve(meta.n_variants);
    
    int num_variant_batches = (meta.n_variants + VARIANT_BATCH_SIZE - 1) / VARIANT_BATCH_SIZE;
    
    for (int vbatch = 0; vbatch < num_variant_batches; ++vbatch) {
        int vbatch_start = vbatch * VARIANT_BATCH_SIZE;
        int vbatch_size = std::min(VARIANT_BATCH_SIZE, meta.n_variants - vbatch_start);
        
        logInfo("Processing variant batch " + std::to_string(vbatch + 1) + "/" +
                std::to_string(num_variant_batches) + " (" + std::to_string(vbatch_size) + " variants)");
        
        // Flatten genotypes for this batch: G is V_batch x S (row-major)
        logInfo("  Flattening genotypes...");
        std::vector<double> G_flat(vbatch_size * traits_meta_.n_samples);
        for (int v = 0; v < vbatch_size; ++v) {
            std::copy(genotypes[vbatch_start + v].begin(), genotypes[vbatch_start + v].end(),
                      G_flat.begin() + v * traits_meta_.n_samples);
        }

        // Allocate result matrix: Y is V_batch x T (row-major)
        logInfo("  Allocating result matrix...");
        std::vector<double> Y_flat(vbatch_size * traits_meta_.n_traits, 0.0);

        // Process traits tile-by-tile to avoid loading the full trait matrix.
        logInfo("  Starting GEMM computation...");
        Timer gemm_timer("BLAS GEMM for variant batch " + std::to_string(vbatch + 1));

        int global_trait_start = 0;
        for (int tile_id = 0; tile_id < traits_meta_.n_tiles; ++tile_id) {
            int traits_in_tile = traits_meta_.tile_trait_counts[tile_id];
            auto tile_data = loadTraitsTileData(tile_id, traits_in_tile);
            logInfo("    Loaded trait tile " + std::to_string(tile_id + 1) + "/" +
                    std::to_string(traits_meta_.n_tiles) + " (" +
                    std::to_string(traits_in_tile) + " traits)");

            int num_trait_batches = (traits_in_tile + TRAIT_BATCH_SIZE - 1) / TRAIT_BATCH_SIZE;
            for (int tbatch = 0; tbatch < num_trait_batches; ++tbatch) {
                int local_start = tbatch * TRAIT_BATCH_SIZE;
                int tbatch_size = std::min(TRAIT_BATCH_SIZE, traits_in_tile - local_start);

                const double* w_ptr = tile_data.data() +
                    static_cast<size_t>(local_start) * traits_meta_.n_samples;
                double* y_ptr = Y_flat.data() + global_trait_start + local_start;

                cblas_dgemm(CblasRowMajor,
                            CblasNoTrans,
                            CblasTrans,
                            vbatch_size,
                            tbatch_size,
                            traits_meta_.n_samples,
                            1.0 / traits_meta_.n_samples,
                            G_flat.data(),
                            traits_meta_.n_samples,
                            w_ptr,
                            traits_meta_.n_samples,
                            0.0,
                            y_ptr,
                            traits_meta_.n_traits);
            }

            global_trait_start += traits_in_tile;
        }

        if (global_trait_start != traits_meta_.n_traits) {
            throw std::runtime_error("Trait tile metadata mismatch during simulation");
        }
        
        logInfo("  GEMM completed for " + std::to_string(vbatch_size) + " variants");

        // Convert flat result to scaled dosages for this batch
        logInfo("  Scaling results...");
        for (int v = 0; v < vbatch_size; ++v) {
            std::vector<double> row;
            row.reserve(traits_meta_.n_traits);
            for (int t = 0; t < traits_meta_.n_traits; ++t) {
                row.push_back(Y_flat[v * traits_meta_.n_traits + t]);
            }
            
            auto scaled = scaleToDosageRange(row);
            all_synthetic_dosages.push_back(std::move(scaled));
        }
    }
    
    logInfo("All variants processed and scaled");

    // Write output
    logInfo("Starting to write output for chunk " + std::to_string(chunk_id));
    writeChunkVCF(chunk_id, meta, all_synthetic_dosages);
    logInfo("Chunk " + std::to_string(chunk_id) + " completed successfully");
}

void ChunkedSimulator::writeChunkVCF(int chunk_id, const ChunkMetadata& meta,
                                     const std::vector<std::vector<double>>& synthetic_dosages) {
    Timer timer("Writing chunk " + std::to_string(chunk_id));
    std::string output_file = config_.output_dir + "/chunk_" + std::to_string(chunk_id) + ".vcf";
    if (config_.compress_output) {
        output_file += ".gz";
    }

    if (config_.compress_output) {
        BGZF* fp = bgzf_open(output_file.c_str(), "w6");
        if (!fp) {
            throw std::runtime_error("Failed to open output file: " + output_file);
        }

        // Write header in smaller pieces
        std::string header = "##fileformat=VCFv4.1\n";
        header += "##source=safeld\n";
        for (const auto& contig_line : contig_header_lines_) {
            header += contig_line + "\n";
        }
        header += "##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Dosage\">\n";
        bgzf_write(fp, header.c_str(), header.length());
        
        // Write column header line
        std::string col_header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
        bgzf_write(fp, col_header.c_str(), col_header.length());
        
        // Write trait column names one at a time
        for (int i = 1; i <= traits_meta_.n_traits; i++) {
            std::string trait_col = "\t" + std::to_string(i);
            bgzf_write(fp, trait_col.c_str(), trait_col.length());
        }
        bgzf_write(fp, "\n", 1);

        // Write variants
        for (int v = 0; v < meta.n_variants; v++) {
            std::string line;
            line.reserve(150000);
            line = meta.chroms[v] + "\t" +
                   std::to_string(meta.positions[v]) + "\t" +
                   meta.variant_ids[v] + "\t" +
                   meta.refs[v] + "\t" +
                   meta.alts[v] + "\t.\t.\t.\tDS";

            char dosage_buf[16];
            for (double dosage : synthetic_dosages[v]) {
                snprintf(dosage_buf, sizeof(dosage_buf), "\t%.4f", dosage);
                line += dosage_buf;
            }
            line += "\n";

            bgzf_write(fp, line.c_str(), line.length());

            if ((v + 1) % 1000 == 0 || v == meta.n_variants - 1) {
                logInfo("  Written " + std::to_string(v + 1) + "/" +
                        std::to_string(meta.n_variants) + " variants");
            }
        }
        bgzf_close(fp);
    } else {
        std::ofstream out(output_file);
        if (!out) {
            throw std::runtime_error("Failed to open output file: " + output_file);
        }

        // Write header
        out << "##fileformat=VCFv4.1\n";
        out << "##source=safeld\n";
        for (const auto& contig_line : contig_header_lines_) {
            out << contig_line << "\n";
        }
        out << "##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Dosage\">\n";
        out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
        for (int i = 1; i <= traits_meta_.n_traits; i++) {
            out << "\t" << i;
        }
        out << "\n";

        // Write variants
        for (int v = 0; v < meta.n_variants; v++) {
            out << meta.chroms[v] << "\t"
                << meta.positions[v] << "\t"
                << meta.variant_ids[v] << "\t"
                << meta.refs[v] << "\t"
                << meta.alts[v] << "\t.\t.\t.\tDS";
            for (double dosage : synthetic_dosages[v]) {
                out << "\t" << std::fixed << std::setprecision(4) << dosage;
            }
            out << "\n";
            if ((v + 1) % 1000 == 0 || v == meta.n_variants - 1) {
                logInfo("  Written " + std::to_string(v + 1) + "/" +
                        std::to_string(meta.n_variants) + " variants");
            }
        }
    }
    logInfo("Written chunk " + std::to_string(chunk_id) + " to: " + output_file);
}

void ChunkedSimulator::run() {
    logInfo("Starting chunked simulation...");
    
    // Load metadata once; trait tiles are streamed on demand.
    loadTraitsMetadata();
    loadHeaderMetadata();
    logInfo("Trait tiles will be streamed on demand during simulation");

    // Create output directory
    if (mkdir(config_.output_dir.c_str(), 0755) != 0 && errno != EEXIST) {
        throw std::runtime_error("Failed to create output directory");
    }

    // Determine which chunks to process
    std::string chunks_dir = config_.preprocessed_dir + "/chunks";
    std::vector<int> chunk_ids;
    if (config_.start_chunk >= 0 && config_.end_chunk >= 0) {
        // Process specified range
        for (int i = config_.start_chunk; i <= config_.end_chunk; ++i) {
            chunk_ids.push_back(i);
        }
    } else {
        // Find all available chunks
        int chunk_id = 0;
        while (true) {
            std::string meta_file = chunks_dir + "/chunk_" +
                                   std::to_string(chunk_id) + ".meta";
            std::ifstream test(meta_file);
            if (!test) break;
            chunk_ids.push_back(chunk_id++);
        }
    }
    logInfo("Processing " + std::to_string(chunk_ids.size()) + " chunks");

    // Process chunks sequentially
    for (int chunk_id : chunk_ids) {
        simulateChunk(chunk_id);
    }

    logInfo("Chunked simulation completed!");
}
