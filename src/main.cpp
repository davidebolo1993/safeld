#include <iostream>
#include <fstream>
#include <string>
#include "genotype_source.h"
#include "preprocessor.h"
#include "chunked_simulator.h"
#include "vcf_merger.h"
#include "utils.h"

void printUsage(const char* program_name) {
    std::cout << "SAFELD - Three-stage workflow\n\n";
    std::cout << "Stage 1 - Preprocessing:\n";
    std::cout << "  " << program_name << " preprocess [OPTIONS]\n";
    std::cout << "    -vcf FILE            Input VCF file\n";
    std::cout << "    -pfile PREFIX        Input plink2 .pgen/.pvar/.psam (requires SAFELD_PGEN)\n";
    std::cout << "    -bfile PREFIX        Input plink1 .bed/.bim/.fam  (requires SAFELD_PGEN)\n";
    std::cout << "    -out DIR             Output directory for preprocessed data\n";
    std::cout << "    -samples LIST        Comma-separated sample IDs, or a file with one per line\n";
    std::cout << "    -extract FILE        Keep only these variant IDs, one per line\n";
    std::cout << "    -maf FLOAT           MAF filter (default: 0.01)\n";
    std::cout << "    -max-missing FLOAT   Max fraction of missing calls per variant (default: 0.1)\n";
    std::cout << "    -dosage-field FIELD  auto|DS|GT: which FORMAT field to read (default: auto)\n";
    std::cout << "    -ntraits INT         Number of traits (default: 10)\n";
    std::cout << "    -chunk-size INT      Variants per chunk (default: 10000)\n";
    std::cout << "    -traits-per-tile INT Traits per tile (default: auto, ~1GB tiles)\n";
    std::cout << "    -h, --help           Show this help message\n\n";
    std::cout << "Global options:\n";
    std::cout << "    -verbose             Enable detailed debug logging\n\n";
    std::cout << "  Note: input VCF must be coordinate-sorted and biallelic\n";
    std::cout << "        (split with: bcftools norm -m -any).\n\n";

    std::cout << "Stage 2 - Simulation:\n";
    std::cout << "  " << program_name << " simulate [OPTIONS]\n";
    std::cout << "    -prep DIR            Preprocessed data directory (required)\n";
    std::cout << "    -out DIR             Output directory for results\n";
    std::cout << "    -workers INT         Number of threads (default: auto)\n";
    std::cout << "    -variant-batch-size INT Variants per simulation batch (default: 4000)\n";
    std::cout << "    -compress            Compress output\n";
    std::cout << "    -start-chunk INT     First chunk to process (default: all)\n";
    std::cout << "    -end-chunk INT       Last chunk to process (default: all)\n";
    std::cout << "    -h, --help           Show this help message\n\n";

    std::cout << "Stage 3 - Merge:\n";
    std::cout << "  " << program_name << " merge [OPTIONS]\n";
    std::cout << "    -in DIR              Directory with chunk VCF files (required)\n";
    std::cout << "    -out FILE            Output merged VCF file\n";
    std::cout << "    -no-compress         Don't compress output (default: compressed)\n";
    std::cout << "    -no-index            Don't create tabix index for compressed output\n";
    std::cout << "    -no-sort             Skip sortedness enforcement during merge\n";
    std::cout << "    -h, --help           Show this help message\n\n";
}

int main(int argc, char* argv[]) {
    try {
        if (argc < 2) {
            printUsage(argv[0]);
            return 1;
        }

        std::string mode = argv[1];

        if (mode == "-h" || mode == "--help" || mode == "help") {
            printUsage(argv[0]);
            return 0;
        }

        bool verbose_logs = false;
        for (int i = 1; i < argc; ++i) {
            std::string arg = argv[i];
            if (arg == "-verbose") {
                verbose_logs = true;
            }
        }
        setVerboseLogging(verbose_logs);

        if (mode == "preprocess") {
            PreprocessConfig config;

            for (int i = 2; i < argc; ++i) {
                std::string arg = argv[i];
                if (arg == "-h" || arg == "--help") {
                    std::cout << "SAFELD Preprocessing\n\n";
                    std::cout << "Usage: " << argv[0] << " preprocess [OPTIONS]\n\n";
                    std::cout << "Options:\n";
                    std::cout << "  -vcf FILE            Input VCF file\n";
                    std::cout << "  -pfile PREFIX        Input plink2 .pgen/.pvar/.psam\n";
                    std::cout << "  -bfile PREFIX        Input plink1 .bed/.bim/.fam\n";
                    std::cout << "                       (exactly one of -vcf/-pfile/-bfile)\n";
                    std::cout << "  -out DIR             Output directory for preprocessed data (required)\n";
                    std::cout << "  -samples LIST        Comma-separated sample IDs, or a file with one per line\n";
                    std::cout << "  -extract FILE        Keep only these variant IDs, one per line\n";
                    std::cout << "  -maf FLOAT           MAF filter (default: 0.01)\n";
                    std::cout << "  -max-missing FLOAT   Max fraction of missing calls per variant (default: 0.1)\n";
                    std::cout << "  -dosage-field FIELD  auto|DS|GT: which FORMAT field to read (default: auto)\n";
                    std::cout << "                       auto prefers DS but falls back to GT where DS is sparse\n";
                    std::cout << "  -ntraits INT         Number of traits (default: 10)\n";
                    std::cout << "  -chunk-size INT      Variants per chunk (default: 10000)\n";
                    std::cout << "  -traits-per-tile INT Traits per tile (default: auto, ~1GB tiles)\n";
                    std::cout << "  -h, --help           Show this help message\n\n";
                    std::cout << "Note: input VCF must be coordinate-sorted and biallelic\n";
                    std::cout << "      (split with: bcftools norm -m -any).\n\n";
                    return 0;
                }

                if (arg == "-vcf" && i + 1 < argc) {
                    config.vcf_file = argv[++i];
                } else if (arg == "-pfile" && i + 1 < argc) {
                    config.genotype_file = std::string(argv[++i]) + ".pgen";
                    config.plink1_metadata = false;
                } else if (arg == "-bfile" && i + 1 < argc) {
                    config.genotype_file = std::string(argv[++i]) + ".bed";
                    config.plink1_metadata = true;
                } else if (arg == "-pgen" && i + 1 < argc) {
                    config.genotype_file = argv[++i];
                } else if (arg == "-pvar" && i + 1 < argc) {
                    config.variants_file = argv[++i];
                } else if (arg == "-psam" && i + 1 < argc) {
                    config.samples_file = argv[++i];
                } else if (arg == "-extract" && i + 1 < argc) {
                    config.extract_file = argv[++i];
                } else if (arg == "-out" && i + 1 < argc) {
                    config.output_dir = argv[++i];
                } else if (arg == "-samples" && i + 1 < argc) {
                    // Accept either a comma-separated list or a file of IDs, so
                    // a subset does not have to fit on a command line.
                    std::string value = argv[++i];
                    std::ifstream probe(value);
                    if (probe.good() && value.find(',') == std::string::npos) {
                        auto ids = readIdList(value);
                        std::string joined;
                        for (size_t k = 0; k < ids.size(); ++k) {
                            if (k) joined += ",";
                            joined += ids[k];
                        }
                        config.sample_list = joined;
                    } else {
                        config.sample_list = value;
                    }
                } else if (arg == "-maf" && i + 1 < argc) {
                    config.maf_filter = std::stod(argv[++i]);
                } else if (arg == "-max-missing" && i + 1 < argc) {
                    config.max_missing_rate = std::stod(argv[++i]);
                } else if (arg == "-dosage-field" && i + 1 < argc) {
                    std::string field = argv[++i];
                    if (field == "auto") config.dosage_field = DosageField::Auto;
                    else if (field == "DS" || field == "ds") config.dosage_field = DosageField::DS;
                    else if (field == "GT" || field == "gt") config.dosage_field = DosageField::GT;
                    else { logError("-dosage-field must be auto, DS or GT"); return 1; }
                } else if (arg == "-ntraits" && i + 1 < argc) {
                    config.n_traits = std::stoi(argv[++i]);
                } else if (arg == "-chunk-size" && i + 1 < argc) {
                    config.chunk_size = std::stoi(argv[++i]);
                } else if (arg == "-traits-per-tile" && i + 1 < argc) {
                    config.traits_per_tile = std::stoi(argv[++i]);
                }
            }

            const int n_inputs = (!config.vcf_file.empty() ? 1 : 0) +
                                 (!config.genotype_file.empty() ? 1 : 0);
            if (n_inputs == 0 || config.output_dir.empty()) {
                logError("An input (-vcf, -pfile or -bfile) and an output directory are required");
                std::cout << "\nUse: " << argv[0] << " preprocess --help for usage information\n";
                return 1;
            }
            if (n_inputs > 1) {
                logError("Give exactly one of -vcf, -pfile or -bfile");
                return 1;
            }
            if (config.maf_filter < 0.0 || config.maf_filter >= 0.5) {
                logError("-maf must be in [0, 0.5)");
                return 1;
            }
            if (config.max_missing_rate < 0.0 || config.max_missing_rate > 1.0) {
                logError("-max-missing must be in [0, 1]");
                return 1;
            }
            if (config.n_traits < 1) {
                logError("-ntraits must be >= 1");
                return 1;
            }
            if (config.chunk_size < 1) {
                logError("-chunk-size must be >= 1");
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
                    std::cout << "  -variant-batch-size INT Variants per simulation batch (default: 4000)\n";
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
                } else if (arg == "-variant-batch-size" && i + 1 < argc) {
                    config.variant_batch_size = std::stoi(argv[++i]);
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
            if (config.variant_batch_size <= 0) {
                logError("variant-batch-size must be > 0");
                return 1;
            }
            if (config.n_workers < 0) {
                logError("-workers must be >= 0 (0 = auto-detect)");
                return 1;
            }
            if ((config.start_chunk >= 0) != (config.end_chunk >= 0)) {
                logError("-start-chunk and -end-chunk must be given together");
                return 1;
            }
            if (config.start_chunk >= 0 && config.end_chunk < config.start_chunk) {
                logError("-end-chunk must be >= -start-chunk");
                return 1;
            }

            ChunkedSimulator simulator(config);
            simulator.run();

        } else if (mode == "merge") {
            MergerConfig config;
            config.compress_output = true;

            for (int i = 2; i < argc; ++i) {
                std::string arg = argv[i];
                if (arg == "-h" || arg == "--help") {
                    std::cout << "SAFELD Merge\n\n";
                    std::cout << "Usage: " << argv[0] << " merge [OPTIONS]\n\n";
                    std::cout << "Options:\n";
                    std::cout << "  -in DIR            Directory with chunk VCF files (required)\n";
                    std::cout << "  -out FILE          Output merged VCF file (required)\n";
                    std::cout << "  -no-compress       Don't compress output (default: compressed)\n";
                    std::cout << "  -no-index          Don't create tabix index for compressed output\n";
                    std::cout << "  -no-sort           Skip sortedness enforcement during merge\n";
                    std::cout << "  -h, --help         Show this help message\n\n";
                    return 0;
                }

                if (arg == "-in" && i + 1 < argc) {
                    config.input_dir = argv[++i];
                } else if (arg == "-out" && i + 1 < argc) {
                    config.output_file = argv[++i];
                } else if (arg == "-no-compress") {
                    config.compress_output = false;
                } else if (arg == "-no-index") {
                    config.write_index = false;
                } else if (arg == "-no-sort") {
                    config.enforce_sort = false;
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
