#include <iostream>
#include "preprocessor.h"
#include "chunked_simulator.h"
#include "vcf_merger.h"
#include "utils.h"

void printUsage(const char* program_name) {
    std::cout << "SAFELD Chunked Mode - Three-stage workflow\n\n";
    std::cout << "Stage 1 - Preprocessing:\n";
    std::cout << "  " << program_name << " preprocess [OPTIONS]\n";
    std::cout << "    -vcf FILE            Input VCF file (required)\n";
    std::cout << "    -out DIR             Output directory for preprocessed data\n";
    std::cout << "    -samples LIST        Comma-separated sample IDs\n";
    std::cout << "    -maf FLOAT           MAF filter (default: 0.01)\n";
    std::cout << "    -ntraits INT         Number of traits (default: 10)\n";
    std::cout << "    -chunk-size INT      Variants per chunk (default: 10000)\n";
    std::cout << "    -traits-per-tile INT Traits per tile (default: auto, ~1GB tiles)\n";
    std::cout << "    -h, --help           Show this help message\n\n";
    
    std::cout << "Stage 2 - Simulation:\n";
    std::cout << "  " << program_name << " simulate [OPTIONS]\n";
    std::cout << "    -prep DIR            Preprocessed data directory (required)\n";
    std::cout << "    -out DIR             Output directory for results\n";
    std::cout << "    -workers INT         Number of threads (default: auto)\n";
    std::cout << "    -compress            Compress output\n";
    std::cout << "    -start-chunk INT     First chunk to process (default: all)\n";
    std::cout << "    -end-chunk INT       Last chunk to process (default: all)\n";
    std::cout << "    -h, --help           Show this help message\n\n";
    
    std::cout << "Stage 3 - Merge:\n";
    std::cout << "  " << program_name << " merge [OPTIONS]\n";
    std::cout << "    -in DIR              Directory with chunk VCF files (required)\n";
    std::cout << "    -out FILE            Output merged VCF file\n";
    std::cout << "    -compress            Compress output (default: true)\n";
    std::cout << "    -h, --help           Show this help message\n\n";
}

int main(int argc, char* argv[]) {
    try {
        if (argc < 2) {
            printUsage(argv[0]);
            return 1;
        }
        
        std::string mode = argv[1];
        
        // Handle help before initializing anything
        if (mode == "-h" || mode == "--help" || mode == "help") {
            printUsage(argv[0]);
            return 0;
        }
        
        // Now that we're past help checks, log initialization
        logInfo("Vector pool initialized with max size: 2000");
        
        if (mode == "preprocess") {
            PreprocessConfig config;
            
            for (int i = 2; i < argc; ++i) {
                std::string arg = argv[i];
                if (arg == "-h" || arg == "--help") {
                    std::cout << "SAFELD Preprocessing\n\n";
                    std::cout << "Usage: " << argv[0] << " preprocess [OPTIONS]\n\n";
                    std::cout << "Options:\n";
                    std::cout << "  -vcf FILE            Input VCF file (required)\n";
                    std::cout << "  -out DIR             Output directory for preprocessed data (required)\n";
                    std::cout << "  -samples LIST        Comma-separated sample IDs\n";
                    std::cout << "  -maf FLOAT           MAF filter (default: 0.01)\n";
                    std::cout << "  -ntraits INT         Number of traits (default: 10)\n";
                    std::cout << "  -chunk-size INT      Variants per chunk (default: 10000)\n";
                    std::cout << "  -traits-per-tile INT Traits per tile (default: auto, ~1GB tiles)\n";
                    std::cout << "  -h, --help           Show this help message\n\n";
                    std::cout << "Examples:\n";
                    std::cout << "  # Auto tile size (recommended)\n";
                    std::cout << "  " << argv[0] << " preprocess -vcf input.vcf.gz -out prep -ntraits 10000\n\n";
                    std::cout << "  # Custom tile size (250 traits per tile)\n";
                    std::cout << "  " << argv[0] << " preprocess -vcf input.vcf.gz -out prep -ntraits 10000 -traits-per-tile 250\n\n";
                    return 0;
                }
                
                if (arg == "-vcf" && i + 1 < argc) {
                    config.vcf_file = argv[++i];
                } else if (arg == "-out" && i + 1 < argc) {
                    config.output_dir = argv[++i];
                } else if (arg == "-samples" && i + 1 < argc) {
                    config.sample_list = argv[++i];
                } else if (arg == "-maf" && i + 1 < argc) {
                    config.maf_filter = std::stod(argv[++i]);
                } else if (arg == "-ntraits" && i + 1 < argc) {
                    config.n_traits = std::stoi(argv[++i]);
                } else if (arg == "-chunk-size" && i + 1 < argc) {
                    config.chunk_size = std::stoi(argv[++i]);
                } else if (arg == "-traits-per-tile" && i + 1 < argc) {
                    config.traits_per_tile = std::stoi(argv[++i]);
                }
            }
            
            if (config.vcf_file.empty() || config.output_dir.empty()) {
                logError("VCF file and output directory are required");
                std::cout << "\nUse: " << argv[0] << " preprocess --help for usage information\n";
                return 1;
            }
            
            Preprocessor preprocessor(config);
            preprocessor.run();
            
        } else if (mode == "simulate") {
            SimulationConfig config;
            
            for (int i = 2; i < argc; ++i) {
                std::string arg = argv[i];
                if (arg == "-h" || arg == "--help") {
                    std::cout << "SAFELD Simulation\n\n";
                    std::cout << "Usage: " << argv[0] << " simulate [OPTIONS]\n\n";
                    std::cout << "Options:\n";
                    std::cout << "  -prep DIR          Preprocessed data directory (required)\n";
                    std::cout << "  -out DIR           Output directory for results (required)\n";
                    std::cout << "  -workers INT       Number of threads (default: auto-detect)\n";
                    std::cout << "  -compress          Compress output VCF chunks\n";
                    std::cout << "  -start-chunk INT   First chunk to process (default: all)\n";
                    std::cout << "  -end-chunk INT     Last chunk to process (default: all)\n";
                    std::cout << "  -h, --help         Show this help message\n\n";
                    return 0;
                }
                
                if (arg == "-prep" && i + 1 < argc) {
                    config.preprocessed_dir = argv[++i];
                } else if (arg == "-out" && i + 1 < argc) {
                    config.output_dir = argv[++i];
                } else if (arg == "-workers" && i + 1 < argc) {
                    config.n_workers = std::stoi(argv[++i]);
                } else if (arg == "-compress") {
                    config.compress_output = true;
                } else if (arg == "-start-chunk" && i + 1 < argc) {
                    config.start_chunk = std::stoi(argv[++i]);
                } else if (arg == "-end-chunk" && i + 1 < argc) {
                    config.end_chunk = std::stoi(argv[++i]);
                }
            }
            
            if (config.preprocessed_dir.empty() || config.output_dir.empty()) {
                logError("Preprocessed directory and output directory are required");
                std::cout << "\nUse: " << argv[0] << " simulate --help for usage information\n";
                return 1;
            }
            
            ChunkedSimulator simulator(config);
            simulator.run();
            
        } else if (mode == "merge") {
            MergerConfig config;
            config.compress_output = true;  // default
            
            for (int i = 2; i < argc; ++i) {
                std::string arg = argv[i];
                if (arg == "-h" || arg == "--help") {
                    std::cout << "SAFELD Merge\n\n";
                    std::cout << "Usage: " << argv[0] << " merge [OPTIONS]\n\n";
                    std::cout << "Options:\n";
                    std::cout << "  -in DIR            Directory with chunk VCF files (required)\n";
                    std::cout << "  -out FILE          Output merged VCF file (required)\n";
                    std::cout << "  -no-compress       Don't compress output (default: compressed)\n";
                    std::cout << "  -h, --help         Show this help message\n\n";
                    return 0;
                }
                
                if (arg == "-in" && i + 1 < argc) {
                    config.input_dir = argv[++i];
                } else if (arg == "-out" && i + 1 < argc) {
                    config.output_file = argv[++i];
                } else if (arg == "-no-compress") {
                    config.compress_output = false;
                }
            }
            
            if (config.input_dir.empty() || config.output_file.empty()) {
                logError("Input directory and output file are required");
                std::cout << "\nUse: " << argv[0] << " merge --help for usage information\n";
                return 1;
            }
            
            VCFMerger merger(config);
            merger.run();
            
        } else {
            logError("Unknown mode: " + mode);
            std::cout << "\nUse: " << argv[0] << " --help for usage information\n";
            return 1;
        }
        
        return 0;
        
    } catch (const std::exception& e) {
        logError("Error: " + std::string(e.what()));
        return 1;
    }
}

