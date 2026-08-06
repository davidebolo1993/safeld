#include "chunked_simulator.h"
#include "utils.h"

#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <htslib/bgzf.h>
#include <sys/stat.h>
#include <cerrno>
#include <dlfcn.h>

extern "C" {
    #include <cblas.h>
}

#ifdef _OPENMP
#include <omp.h>
#endif

// glibc only exposes RTLD_DEFAULT under _GNU_SOURCE. This is its definition
// there, and the macro is already provided on the BSDs and macOS.
#ifndef RTLD_DEFAULT
#define RTLD_DEFAULT ((void*)0)
#endif

ChunkedSimulator::ChunkedSimulator(const SimulationConfig& config)
    : config_(config) {
}

// Apply -workers to the threads that actually do the work. Nearly all of the
// runtime is inside cblas_dgemm, and the CBLAS interface has no standard way to
// set a thread count, so the vendor entry point is resolved at runtime instead
// of being linked against.
void ChunkedSimulator::configureThreads() {
    if (config_.n_workers <= 0) {
        logDebug("Thread count: auto-detected by BLAS/OpenMP");
        return;
    }

    const int n_workers = config_.n_workers;
#ifdef _OPENMP
    omp_set_num_threads(n_workers);
#endif

    using SetThreadsFn = void (*)(int);
    SetThreadsFn set_blas_threads = nullptr;
    for (const char* symbol : {"openblas_set_num_threads", "goto_set_num_threads",
                               "MKL_Set_Num_Threads"}) {
        set_blas_threads = reinterpret_cast<SetThreadsFn>(dlsym(RTLD_DEFAULT, symbol));
        if (set_blas_threads) {
            break;
        }
    }

    if (set_blas_threads) {
        set_blas_threads(n_workers);
        logInfo("Using " + std::to_string(n_workers) + " threads");
    } else {
        logWarning("This BLAS does not expose a runtime thread-count setter; set "
                   "OPENBLAS_NUM_THREADS=" + std::to_string(n_workers) +
                   " (or OMP_NUM_THREADS) in the environment instead");
    }
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

    logInfo("Traits: " + formatCount(traits_meta_.n_traits) + " x " +
            formatCount(traits_meta_.n_samples) + " samples in " +
            std::to_string(traits_meta_.n_tiles) + " tile(s)");
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

    logDebug("Loaded " + std::to_string(contig_header_lines_.size()) + " contig header records");
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
    meta.n_variants = -1;
    meta.n_samples = -1;

    std::string line;
    bool in_key_value_header = true;
    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }

        // key/value lines come first, variant records follow.
        auto sep = line.find('=');
        if (in_key_value_header && sep != std::string::npos) {
            std::string key = line.substr(0, sep);
            std::string value = line.substr(sep + 1);
            if (key == "n_variants") meta.n_variants = std::stoi(value);
            else if (key == "n_samples") meta.n_samples = std::stoi(value);
            continue;
        }
        in_key_value_header = false;

        std::istringstream ss(line);
        std::string chrom, id, ref, alt;
        int pos = 0;
        if (!(ss >> chrom >> pos >> id >> ref >> alt)) {
            throw std::runtime_error("Malformed variant record in " + meta_file + ": " + line);
        }
        meta.chroms.push_back(std::move(chrom));
        meta.positions.push_back(pos);
        meta.variant_ids.push_back(std::move(id));
        meta.refs.push_back(std::move(ref));
        meta.alts.push_back(std::move(alt));
    }

    if (meta.n_variants < 0 || meta.n_samples < 0) {
        throw std::runtime_error("Missing n_variants/n_samples in " + meta_file);
    }
    if (static_cast<size_t>(meta.n_variants) != meta.chroms.size()) {
        throw std::runtime_error("Chunk metadata mismatch in " + meta_file + ": n_variants=" +
                                 std::to_string(meta.n_variants) + " but " +
                                 std::to_string(meta.chroms.size()) + " variant records");
    }
    // The chunk and the trait matrix must describe the same samples, otherwise
    // the genotype rows below would be read at the wrong stride.
    if (meta.n_samples != traits_meta_.n_samples) {
        throw std::runtime_error("Sample count mismatch: chunk " + std::to_string(chunk_id) +
                                 " has " + std::to_string(meta.n_samples) +
                                 " samples but the trait matrix has " +
                                 std::to_string(traits_meta_.n_samples) +
                                 "; the preprocessed directory is inconsistent");
    }

    return meta;
}

std::vector<std::vector<double>> ChunkedSimulator::loadChunkGenotypes(int chunk_id, int n_variants) {
    std::string bin_file = config_.preprocessed_dir + "/chunks/chunk_" +
                          std::to_string(chunk_id) + ".bin";
    std::ifstream in(bin_file, std::ios::binary);
    if (!in) {
        throw std::runtime_error("Failed to open chunk genotypes: " + bin_file);
    }

    const size_t row_bytes = static_cast<size_t>(traits_meta_.n_samples) * sizeof(double);
    std::vector<std::vector<double>> genotypes(n_variants);
    for (int v = 0; v < n_variants; ++v) {
        genotypes[v].resize(traits_meta_.n_samples);
        in.read(reinterpret_cast<char*>(genotypes[v].data()),
                static_cast<std::streamsize>(row_bytes));
        // An unchecked short read leaves the row zero-filled, which would sail
        // through the rest of the pipeline and produce silently wrong output.
        if (!in || static_cast<size_t>(in.gcount()) != row_bytes) {
            throw std::runtime_error("Truncated chunk genotype file " + bin_file +
                                     ": expected " + std::to_string(n_variants) +
                                     " variants x " + std::to_string(traits_meta_.n_samples) +
                                     " samples, short read at variant " + std::to_string(v));
        }
    }
    return genotypes;
}

void ChunkedSimulator::simulateChunk(int chunk_id) {
    auto meta = loadChunkMetadata(chunk_id);
    if (meta.n_variants == 0) {
        logWarning("Chunk " + std::to_string(chunk_id) + " contains no variants; skipping");
        return;
    }

    auto genotypes = loadChunkGenotypes(chunk_id, meta.n_variants);

    const int variant_batch_size = std::max(1, config_.variant_batch_size);
    std::vector<std::vector<double>> all_synthetic_dosages;
    all_synthetic_dosages.reserve(meta.n_variants);
    
    int num_variant_batches = (meta.n_variants + variant_batch_size - 1) / variant_batch_size;
    
    ProgressBar chunk_bar("Chunk " + std::to_string(chunk_id), meta.n_variants);
    for (int vbatch = 0; vbatch < num_variant_batches; ++vbatch) {
        int vbatch_start = vbatch * variant_batch_size;
        int vbatch_size = std::min(variant_batch_size, meta.n_variants - vbatch_start);
        
        logDebug("  Flattening genotypes...");
        std::vector<double> G_flat(vbatch_size * traits_meta_.n_samples);
        for (int v = 0; v < vbatch_size; ++v) {
            std::copy(genotypes[vbatch_start + v].begin(), genotypes[vbatch_start + v].end(),
                      G_flat.begin() + v * traits_meta_.n_samples);
        }

        logDebug("  Allocating result matrix...");
        std::vector<double> Y_flat(vbatch_size * traits_meta_.n_traits, 0.0);

        logDebug("  Starting GEMM computation...");
        Timer gemm_timer("BLAS GEMM for variant batch " + std::to_string(vbatch + 1));

        int global_trait_start = 0;
        for (int tile_id = 0; tile_id < traits_meta_.n_tiles; ++tile_id) {
            int traits_in_tile = traits_meta_.tile_trait_counts[tile_id];
            auto tile_data = loadTraitsTileData(tile_id, traits_in_tile);
            logDebug("    Loaded trait tile " + std::to_string(tile_id + 1) + "/" +
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
        
        logDebug("  GEMM completed for " + std::to_string(vbatch_size) + " variants");

        logDebug("  Scaling results...");
        for (int v = 0; v < vbatch_size; ++v) {
            std::vector<double> row;
            row.reserve(traits_meta_.n_traits);
            for (int t = 0; t < traits_meta_.n_traits; ++t) {
                row.push_back(Y_flat[v * traits_meta_.n_traits + t]);
            }
            
            auto scaled = scaleToDosageRange(row);
            all_synthetic_dosages.push_back(std::move(scaled));
        }
        chunk_bar.update(vbatch_start + vbatch_size);
    }
    chunk_bar.finish();
    
    writeChunkVCF(chunk_id, meta, all_synthetic_dosages);
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

        std::string header = "##fileformat=VCFv4.1\n";
        header += "##source=safeld\n";
        for (const auto& contig_line : contig_header_lines_) {
            header += contig_line + "\n";
        }
        header += "##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Dosage\">\n";
        bgzf_write(fp, header.c_str(), header.length());
        
        std::string col_header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
        bgzf_write(fp, col_header.c_str(), col_header.length());
        
        for (int i = 1; i <= traits_meta_.n_traits; i++) {
            std::string trait_col = "\t" + std::to_string(i);
            bgzf_write(fp, trait_col.c_str(), trait_col.length());
        }
        bgzf_write(fp, "\n", 1);

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

        }
        bgzf_close(fp);
    } else {
        std::ofstream out(output_file);
        if (!out) {
            throw std::runtime_error("Failed to open output file: " + output_file);
        }

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
        }
    }
    logDebug("Wrote " + output_file);
}

void ChunkedSimulator::run() {
    LogModule module("simulate");
    Timer timer("simulate stage");
    if (config_.variant_batch_size <= 0) {
        throw std::runtime_error("variant_batch_size must be > 0");
    }
    logDebug("Variant batch size: " + formatCount(config_.variant_batch_size));

    configureThreads();
    loadTraitsMetadata();
    loadHeaderMetadata();
    if (traits_meta_.n_traits < 2) {
        // Output dosages are min-max scaled across traits within each variant, so
        // a single trait leaves nothing to scale against and every dosage
        // collapses to the same value.
        logWarning("Only " + std::to_string(traits_meta_.n_traits) + " trait(s): per-variant "
                   "scaling has no spread to work with and every output dosage will be identical");
    }

    if (mkdir(config_.output_dir.c_str(), 0755) != 0 && errno != EEXIST) {
        throw std::runtime_error("Failed to create output directory");
    }

    std::string chunks_dir = config_.preprocessed_dir + "/chunks";
    std::vector<int> chunk_ids;
    if (config_.start_chunk >= 0 && config_.end_chunk >= 0) {
        for (int i = config_.start_chunk; i <= config_.end_chunk; ++i) {
            chunk_ids.push_back(i);
        }
    } else {
        int chunk_id = 0;
        while (true) {
            std::string meta_file = chunks_dir + "/chunk_" +
                                   std::to_string(chunk_id) + ".meta";
            std::ifstream test(meta_file);
            if (!test) break;
            chunk_ids.push_back(chunk_id++);
        }
    }
    logInfo("Processing " + std::to_string(chunk_ids.size()) + " chunk(s)");

    for (int chunk_id : chunk_ids) {
        simulateChunk(chunk_id);
    }

    logInfo("Done in " + formatDuration(timer.elapsed()) + " -> " + config_.output_dir);
}
