#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <htslib/bgzf.h>
#include <cstring>
#include <cerrno>
#include "vcf_processor.h"
#include "simulation_engine.h"
#include "utils.h"

struct CommandLineArgs {
    std::string vcf_file;
    std::string output_file = "SAFE_LD.vcf";
    std::string sample_list;
    double maf_filter = 0.01;
    int n_traits = 10;
    int n_workers = 0; // 0 = auto-detect
    bool compress_output = false;
};

void printUsage(const char* program_name) {
    std::cout << "Usage: " << program_name << " [OPTIONS]\n\n"
              << "Options:\n"
              << "  -vcf FILE        Input VCF file (required, supports .vcf, .vcf.gz)\n"
              << "  -out FILE        Output VCF file (default: SAFE_LD.vcf)\n"
              << "  -samples LIST    Comma-separated sample IDs to include\n"
              << "  -maf FLOAT       Minimum allele frequency filter (default: 0.01)\n"
              << "  -ntraits INT     Number of synthetic traits (default: 10)\n"
              << "  -workers INT     Number of worker threads (default: auto-detect)\n"
              << "  -compress        Compress output (bgzip)\n"
              << "  -h, --help       Show this help message\n\n";
}

CommandLineArgs parseCommandLine(int argc, char* argv[]) {
    CommandLineArgs args;
    
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        
        if (arg == "-h" || arg == "--help") {
            printUsage(argv[0]);
            exit(0);
        } else if (arg == "-vcf" && i + 1 < argc) {
            args.vcf_file = argv[++i];
        } else if (arg == "-out" && i + 1 < argc) {
            args.output_file = argv[++i];
        } else if (arg == "-samples" && i + 1 < argc) {
            args.sample_list = argv[++i];
        } else if (arg == "-maf" && i + 1 < argc) {
            args.maf_filter = std::stod(argv[++i]);
        } else if (arg == "-ntraits" && i + 1 < argc) {
            args.n_traits = std::stoi(argv[++i]);
        } else if (arg == "-workers" && i + 1 < argc) {
            args.n_workers = std::stoi(argv[++i]);
        } else if (arg == "-compress") {
            args.compress_output = true;
        }
    }
    
    return args;
}

void writeVCFHeader(std::ostream& out, int n_traits) {
    out << "##fileformat=VCFv4.1\n";
    out << "##source=safeld-cpp\n";
    out << "##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Dosage\">\n";
    out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    for (int i = 1; i <= n_traits; ++i) {
        out << "\tT" << i;
    }
    out << "\n";
}

void writeVariant(std::ostream& out, const ProcessedVariant& variant) {
    out << variant.chrom << '\t'
        << variant.pos << '\t'
        << variant.id << '\t'
        << variant.ref << '\t'
        << variant.alt << '\t'
        << ".\tPASS\t.\tDS";
    
    for (double dosage : variant.synthetic_dosages) {
        out << '\t' << std::fixed << std::setprecision(4) << dosage;
    }
    out << '\n';
}

bool writeResults(const std::string& output_file, const std::vector<ProcessedVariant>& results,
                 int n_traits, bool compress) {
    Timer timer("Writing results");
    logInfo("Writing " + std::to_string(results.size()) + " variants to: " + output_file);
    
    if (compress || output_file.ends_with(".gz")) {
        // Use bgzf for compressed output
        BGZF* fp = bgzf_open(output_file.c_str(), "w");
        if (!fp) {
            logError("Failed to open output file for writing: " + output_file);
            logError("Check disk space and file permissions");
            return false;
        }
        
        // Write header
        std::ostringstream header;
        writeVCFHeader(header, n_traits);
        std::string header_str = header.str();
        ssize_t written = bgzf_write(fp, header_str.c_str(), header_str.length());
        if (written != static_cast<ssize_t>(header_str.length())) {
            logError("Failed to write VCF header: " + std::string(strerror(errno)));
            bgzf_close(fp);
            return false;
        }
        
        // Write variants with error checking
        for (size_t i = 0; i < results.size(); ++i) {
            std::ostringstream line;
            writeVariant(line, results[i]);
            std::string line_str = line.str();
            written = bgzf_write(fp, line_str.c_str(), line_str.length());
            if (written != static_cast<ssize_t>(line_str.length())) {
                logError("Failed to write variant " + std::to_string(i) + ": " + std::string(strerror(errno)));
                bgzf_close(fp);
                return false;
            }
            
            if (i % 10000 == 0 && i > 0) {
                logInfo("Written " + std::to_string(i) + " variants...");
            }
        }
        
        if (bgzf_close(fp) != 0) {
            logError("Failed to close output file properly: " + std::string(strerror(errno)));
            return false;
        }
        
    } else {
        // Regular file output with better error handling
        std::ofstream out(output_file);
        if (!out) {
            logError("Failed to open output file for writing: " + output_file);
            logError("Error: " + std::string(strerror(errno)));
            return false;
        }
        
        writeVCFHeader(out, n_traits);
        if (!out.good()) {
            logError("Failed to write VCF header");
            return false;
        }
        
        for (size_t i = 0; i < results.size(); ++i) {
            writeVariant(out, results[i]);
            if (!out.good()) {
                logError("Failed to write variant " + std::to_string(i));
                return false;
            }
            
            if (i % 10000 == 0 && i > 0) {
                logInfo("Written " + std::to_string(i) + " variants...");
            }
        }
    }
    
    logInfo("Results written successfully to: " + output_file);
    return true;
}

int main(int argc, char* argv[]) {
    try {
        // Parse command line arguments (handles help before any initialization)
        CommandLineArgs args = parseCommandLine(argc, argv);
        
        if (args.vcf_file.empty()) {
            logError("VCF file is required. Use -h for help.");
            return 1;
        }
        
        // Now that we're past help, initialize logging
        logInfo("Vector pool initialized with max size: 2000");
        logInfo("Starting SAFE_LD simulation...");
        logInfo("Input VCF: " + args.vcf_file);
        logInfo("Output file: " + args.output_file);
        logInfo("MAF filter: " + std::to_string(args.maf_filter));
        logInfo("Number of traits: " + std::to_string(args.n_traits));
        
        // Step 1: Process VCF
        VCFProcessor processor(args.vcf_file, args.maf_filter);
        if (!processor.initialize(args.sample_list)) {
            logError("Failed to initialize VCF processor");
            return 1;
        }
        
        auto variants = processor.processVariants();
        if (variants.empty()) {
            logError("No variants passed filtering");
            return 1;
        }
        
        int n_samples = processor.getTargetSamples().size();
        
        // Step 2: Initialize simulation engine
        SimulationEngine engine(args.n_traits, n_samples, args.n_workers);
        engine.initialize();
        
        // Step 3: Run simulation (using original method)
        auto results = engine.simulateVariants(variants);
        
        // Step 4: Write results
        if (!writeResults(args.output_file, results, args.n_traits, args.compress_output)) {
            return 1;
        }
        
        logInfo("SAFE_LD simulation completed successfully!");
        logInfo("Processed " + std::to_string(variants.size()) + " variants with " +
                std::to_string(n_samples) + " samples");
        
        return 0;
        
    } catch (const std::exception& e) {
        logError("Error: " + std::string(e.what()));
        return 1;
    }
}

