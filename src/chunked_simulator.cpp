#include "chunked_simulator.h"
#include "utils.h"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <thread>
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
    
    logInfo("Traits metadata: " + std::to_string(traits_meta_.n_traits) + " traits, " +
           std::to_string(traits_meta_.n_samples) + " samples, " +
           std::to_string(traits_meta_.n_tiles) + " tiles");
}

void ChunkedSimulator::loadTraitsMatrix() {
    Timer timer("Loading traits matrix");
    
    traits_matrix_.resize(traits_meta_.n_traits);
    
    int current_trait = 0;
    for (int tile_id = 0; tile_id < traits_meta_.n_tiles; ++tile_id) {
        int traits_in_tile = traits_meta_.tile_trait_counts[tile_id];
        loadTraitsTile(tile_id, current_trait, traits_in_tile);
        current_trait += traits_in_tile;
        
        logInfo("Loaded tile " + std::to_string(tile_id + 1) + "/" +
               std::to_string(traits_meta_.n_tiles));
    }
}

void ChunkedSimulator::loadTraitsTile(int tile_id, int start_trait, int n_traits) {
    std::string tile_file = config_.preprocessed_dir + "/traits/W_tile_" +
                           std::to_string(tile_id) + ".bin";
    std::ifstream in(tile_file, std::ios::binary);
    if (!in) {
        throw std::runtime_error("Failed to open trait tile: " + tile_file);
    }
    
    for (int t = 0; t < n_traits; ++t) {
        traits_matrix_[start_trait + t].resize(traits_meta_.n_samples);
        in.read(reinterpret_cast<char*>(traits_matrix_[start_trait + t].data()),
               traits_meta_.n_samples * sizeof(double));
    }
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
    Timer timer("Simulating chunk " + std::to_string(chunk_id));
    
    // Load chunk data
    auto meta = loadChunkMetadata(chunk_id);
    auto genotypes = loadChunkGenotypes(chunk_id, meta.n_variants);
    
    logInfo("Processing chunk " + std::to_string(chunk_id) + ": " +
           std::to_string(meta.n_variants) + " variants");
    
    // Flatten genotypes matrix for BLAS: G is V x S (row-major)
    std::vector<double> G_flat(meta.n_variants * traits_meta_.n_samples);
    for (int v = 0; v < meta.n_variants; ++v) {
        std::copy(genotypes[v].begin(), genotypes[v].end(), 
                 G_flat.begin() + v * traits_meta_.n_samples);
    }
    
    // Flatten traits matrix for BLAS: W is T x S (row-major)
    std::vector<double> W_flat(traits_meta_.n_traits * traits_meta_.n_samples);
    for (int t = 0; t < traits_meta_.n_traits; ++t) {
        std::copy(traits_matrix_[t].begin(), traits_matrix_[t].end(),
                 W_flat.begin() + t * traits_meta_.n_samples);
    }
    
    // Allocate result matrix: Y is V x T (row-major)
    std::vector<double> Y_flat(meta.n_variants * traits_meta_.n_traits);
    
    // Compute Y = G * W^T using GEMM
    // G: V x S, W^T: S x T, Y: V x T
    // cblas_dgemm(Layout, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc)
    // Y = alpha * op(G) * op(W) + beta * Y
    // We want: Y = (1/S) * G * W^T
    Timer gemm_timer("BLAS GEMM computation");
    cblas_dgemm(CblasRowMajor,           // Row-major layout
                CblasNoTrans,            // Don't transpose G
                CblasTrans,              // Transpose W (to get W^T)
                meta.n_variants,         // M = number of rows in G
                traits_meta_.n_traits,   // N = number of columns in W^T
                traits_meta_.n_samples,  // K = shared dimension
                1.0 / traits_meta_.n_samples,  // alpha = 1/S for normalization
                G_flat.data(),           // Matrix G
                traits_meta_.n_samples,  // Leading dimension of G
                W_flat.data(),           // Matrix W
                traits_meta_.n_samples,  // Leading dimension of W
                0.0,                     // beta = 0 (don't add to existing Y)
                Y_flat.data(),           // Result matrix Y
                traits_meta_.n_traits);  // Leading dimension of Y
    
    logInfo("GEMM completed for " + std::to_string(meta.n_variants) + " variants x " +
           std::to_string(traits_meta_.n_traits) + " traits");
    
    // Convert flat result back to 2D and scale to [0, 2] range
    std::vector<std::vector<double>> synthetic_dosages(meta.n_variants,
                                                       std::vector<double>(traits_meta_.n_traits));
    
    #pragma omp parallel for if(config_.n_workers > 0) num_threads(config_.n_workers)
    for (int v = 0; v < meta.n_variants; ++v) {
        // Extract row for this variant
        std::vector<double> row(traits_meta_.n_traits);
        for (int t = 0; t < traits_meta_.n_traits; ++t) {
            row[t] = Y_flat[v * traits_meta_.n_traits + t];
        }
        
        // Scale to [0, 2] dosage range
        auto scaled = scaleToDosageRange(row);
        synthetic_dosages[v] = scaled;
    }
    
    // Write output
    writeChunkVCF(chunk_id, meta, synthetic_dosages);
}

void ChunkedSimulator::writeChunkVCF(int chunk_id, const ChunkMetadata& meta,
                                     const std::vector<std::vector<double>>& synthetic_dosages) {
    Timer timer("Writing chunk " + std::to_string(chunk_id));
    
    std::string output_file = config_.output_dir + "/chunk_" + 
                             std::to_string(chunk_id) + ".vcf";
    if (config_.compress_output) {
        output_file += ".gz";
    }
    
    if (config_.compress_output) {
        BGZF* fp = bgzf_open(output_file.c_str(), "w");
        if (!fp) {
            throw std::runtime_error("Failed to open output file: " + output_file);
        }
        
        // Write header
        std::ostringstream header;
        header << "##fileformat=VCFv4.1\n";
        header << "##source=safeld-chunked\n";
        header << "##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Dosage\">\n";
        header << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
        for (int i = 1; i <= traits_meta_.n_traits; ++i) {
            header << "\tT" << i;
        }
        header << "\n";
        std::string header_str = header.str();
        
        ssize_t written = bgzf_write(fp, header_str.c_str(), header_str.length());
        if (written != static_cast<ssize_t>(header_str.length())) {
            bgzf_close(fp);
            throw std::runtime_error("Failed to write header to: " + output_file);
        }
        
        // Pre-allocate buffer for batch writing
        std::ostringstream buffer;
        const int BATCH_SIZE = 1000;
        
        for (int v = 0; v < meta.n_variants; ++v) {
            buffer << meta.chroms[v] << '\t'
                   << meta.positions[v] << '\t'
                   << meta.variant_ids[v] << '\t'
                   << meta.refs[v] << '\t'
                   << meta.alts[v] << '\t'
                   << ".\tPASS\t.\tDS";
            
            for (double dosage : synthetic_dosages[v]) {
                buffer << '\t' << std::fixed << std::setprecision(4) << dosage;
            }
            buffer << '\n';
            
            // Flush buffer periodically
            if ((v + 1) % BATCH_SIZE == 0 || v == meta.n_variants - 1) {
                std::string batch_str = buffer.str();
                written = bgzf_write(fp, batch_str.c_str(), batch_str.length());
                if (written != static_cast<ssize_t>(batch_str.length())) {
                    bgzf_close(fp);
                    throw std::runtime_error("Failed to write variants to: " + output_file);
                }
                buffer.str("");  // Clear buffer
                buffer.clear();
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
        out << "##source=safeld-chunked\n";
        out << "##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Dosage\">\n";
        out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
        for (int i = 1; i <= traits_meta_.n_traits; ++i) {
            out << "\tT" << i;
        }
        out << "\n";
        
        // Batch write for better I/O performance
        std::ostringstream buffer;
        const int BATCH_SIZE = 1000;
        
        for (int v = 0; v < meta.n_variants; ++v) {
            buffer << meta.chroms[v] << '\t'
                   << meta.positions[v] << '\t'
                   << meta.variant_ids[v] << '\t'
                   << meta.refs[v] << '\t'
                   << meta.alts[v] << '\t'
                   << ".\tPASS\t.\tDS";
            
            for (double dosage : synthetic_dosages[v]) {
                buffer << '\t' << std::fixed << std::setprecision(4) << dosage;
            }
            buffer << '\n';
            
            // Flush buffer periodically
            if ((v + 1) % BATCH_SIZE == 0 || v == meta.n_variants - 1) {
                out << buffer.str();
                buffer.str("");
                buffer.clear();
            }
        }
    }
    
    logInfo("Written chunk " + std::to_string(chunk_id) + " to: " + output_file);
}

void ChunkedSimulator::run() {
    logInfo("Starting chunked simulation...");
    
    // Load traits once into memory
    loadTraitsMetadata();
    loadTraitsMatrix();
    
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
    
    // Process chunks sequentially (each chunk uses all cores via BLAS)
    for (int chunk_id : chunk_ids) {
        simulateChunk(chunk_id);
    }
    
    logInfo("Chunked simulation completed!");
}

