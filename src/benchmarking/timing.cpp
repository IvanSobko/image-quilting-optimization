#include "timing.h"

#include <dirent.h>
#include <algorithm>
#include <cstring>
#include <string>
#include <ctime>
#include <iostream>
#include <fstream>

#include "PngReader.h"
#include "tsc_x86.h"

//TODO: should we modify cycles required? 1e8 is value from homeworks
#define CYCLES_REQUIRED 1e8
#define RDTSC_LATENCY 26
#define REP 1

// Original timing
void timing::run_timing(int inputBlockRatio) {
    std::vector<std::string> files = read_files("./gallery", "input0_");
    if (files.empty()) {
        perror("Failed to collect files for timing");
        return;
    }

    FILE* results_txt = NULL;
    std::string filename = "timing_results_default_flags.txt";

    #ifdef _CompileFlags // get variable _CompileFlags from CmakeLists.txt
        filename = std::string("timing_results_") + _CompileFlags + ".txt";
    #endif
    
    std::string results_path = std::string("results/blockdiv_") + std::to_string(inputBlockRatio) + '/' + filename;
    results_txt = fopen(results_path.c_str(), "w");

    for (int i = 0; i < files.size(); i++) {
        printf("Timing for %s\n", files[i].c_str());
        ImgData img_data;
        file::read_png_file(files[i].c_str(), img_data);
        img_data.output_w = img_data.width * 2;
        img_data.output_h = img_data.height * 2;
        img_data.block_w = img_data.width / inputBlockRatio;
        img_data.block_h = img_data.height / inputBlockRatio;
        img_data.AllocateOutput();

        ImageQuilting quilting(&img_data);
        double cycles = rdtsc(&quilting);
        printf("Done. Quilting for %s took %lli flops and %f cycles\n\n", files[i].c_str(),
               quilting.getFlopCount(), cycles);

        int64_t data_size = img_data.width * img_data.height;
        int64_t flops = quilting.getFlopCount();
        double flopC = ((double)flops) / cycles;
        fprintf(results_txt, "size=%ix%i, n=%lli, performance=%f, flops=%lli, cycles=%f\n", img_data.width,
                img_data.height,
                data_size, flopC, flops, cycles);
        img_data.FreeOutput();
        img_data.FreeInput();
    }
    fclose(results_txt);
}

double timing::rdtsc(ImageQuilting* quilting) {
    double cycles = 0;
    int num_runs = 1;
    double multiplier = 1;
    myInt64 start, end;

    // Warm-up phase: we determine a number of executions that allows
    // the code to be executed for at least CYCLES_REQUIRED cycles.
    // This helps excluding timing overhead when measuring small runtimes.
    printf("Doing warmup phase...");
    do {
        num_runs = num_runs * multiplier;
        start = start_tsc();
        for (size_t i = 0; i < num_runs; i++) {
            quilting->Synthesis();
        }
        end = stop_tsc(start);

        cycles = (double)end;
        multiplier = (CYCLES_REQUIRED) / (cycles);

    } while (multiplier > 2);

    // Actual performance measurements repeated REP times.
    // We simply store all results and compute medians during post-processing.
    printf("actually measuring performance (%i times).\n", num_runs * REP);
    double total_cycles = 0;
    for (size_t j = 0; j < REP; j++) {
        start = start_tsc();
        for (size_t i = 0; i < num_runs; ++i) {
            quilting->Synthesis();
        }
        end = stop_tsc(start) - RDTSC_LATENCY;

        cycles = (double)end / (double)num_runs;
        total_cycles += cycles;
    }
    total_cycles /= REP;
    return total_cycles;
}

std::vector<std::string> timing::read_files(const std::string& directory,
                                            const std::string& filename_filter) {
    DIR* folder = opendir(directory.c_str());
    if (folder == NULL) {
        perror("Unable to read directory");
        return {};
    }

    std::vector<std::string> result;
    struct dirent* entry;
    while ((entry = readdir(folder))) {
        if (strstr(entry->d_name, filename_filter.c_str()) != NULL) {
            std::string filepath = directory + '/' + entry->d_name;
            result.push_back(filepath);
        }
    }
    closedir(folder);
    return result;
}


// Call the image quilting variant on the specified input and return the number of cycles and flops
timing::TimingData timing::rdtsc_functional(const ImageQuiltingFunction & imageQuiltingFunction, ImgData* imgData, int seed)
{
    double cycles = 0;
    int num_runs = 1;
    double multiplier = 1;
    myInt64 start, end;

    // Warm-up phase: we determine a number of executions that allows
    // the code to be executed for at least CYCLES_REQUIRED cycles.
    // This helps excluding timing overhead when measuring small runtimes.
    do {
        num_runs = num_runs * multiplier;
        start = start_tsc();
        for (size_t i = 0; i < num_runs; i++) {
            imageQuiltingFunction(imgData, seed);
        }
        end = stop_tsc(start);

        cycles = (double)end;
        multiplier = (CYCLES_REQUIRED) / (cycles);

    } while (multiplier > 2);

    // Actual performance measurements repeated REP times.
    // We simply store all results and compute medians during post-processing.
    double total_cycles = 0;
    double total_flops = 0;
    for (size_t j = 0; j < REP; j++) {
        start = start_tsc();
        for (size_t i = 0; i < num_runs; ++i) {
            total_flops += imageQuiltingFunction(imgData, seed);
        }
        end = stop_tsc(start) - RDTSC_LATENCY;

        cycles = (double)end / (double)num_runs;
        total_cycles += cycles;
        total_flops = (total_flops - num_runs) / num_runs;
    }
    total_cycles /= REP;
    total_flops /= REP;
    return TimingData(total_cycles, total_flops);
}

// Get the current date and time as a string
// From ChatGPT
std::string getCurrentDateTime() {
    std::time_t now = std::time(nullptr);
    std::tm timeinfo;
    localtime_r(&now, &timeinfo);

    char buffer[20];
    std::strftime(buffer, sizeof(buffer), "%m-%d-%H-%M-%S", &timeinfo);

    return buffer;
}

// Empty image quilting function for debugging purposes
double timing::EmptyImageQuiltingFunction(ImgData* imgData, int seed)
{
    return -1;
}

// Get the relative paths of files from a specified directory starting with the given filter
// From ChatGPT
std::vector<std::string> getRelativePaths(const std::string & directory, const std::string & filter)
{
    std::vector<std::string> relativePaths;
    for (const auto& entry : std::filesystem::directory_iterator(directory)) {
        if (std::filesystem::is_regular_file(entry.path())) {
            std::string filename = entry.path().filename().string();
            if (filename.compare(0, filter.length(), filter) == 0) {
                relativePaths.push_back(entry.path().relative_path().string());
            }
        }
    }
    return relativePaths;
}

// Run the timing for the image quilting variant and write the results to a file with the given label (and current date and time)
void timing::run_timing_functional(
    const std::string & label,
    const std::string & inputDirectory, const std::string & filenameFilter, const std::string & outputDirectory,
    const ImageQuiltingFunction & imageQuiltingFunction, int outputScale, int blockDivisor, int seed)
{
    // Intro
    std::cout << "Running functional timing for:" << std::endl << "\t" << label << std::endl;

    // Read the input files
    std::vector<std::string> inputFiles = getRelativePaths(inputDirectory, filenameFilter);
    std::cout << "Input files: " << std::endl;
    for (const auto & inputFile : inputFiles) std::cout << "\t" << inputFile << std::endl;

    // WARNING: _CompileFlags must be defined!
    std::string compileFlagsString = std::string(_CompileFlags);
    std::string timeString = getCurrentDateTime();
    std::string filename =
        outputDirectory + "/" + label + "_"
        + compileFlagsString + "_"
        + timeString + ".txt";
    std::cout << "Writing results to: " << std::endl << "\t" << filename << std::endl;

    // Perform the timing
    std::cout << "Performing timing and writing results..." << std::endl;

    // Open the output file
    std::ofstream outputFile(filename);
    if (!outputFile.is_open()) {
        std::cerr << "ERROR: failed to open the file: \"" << filename << "\"" << std::endl;
        return;
    }
    // Output file header
    outputFile << "label=" << label;
    outputFile << ", compileFlags=" << compileFlagsString;
    outputFile << ", time=" << timeString;
    outputFile << ", outputScale=" << outputScale;
    outputFile << ", blockDivisor=" << blockDivisor;
    outputFile << ", seed=" << seed << std::endl;

    // Run the image quilting variant on each of the input files
    for (const auto & inputFile : inputFiles) {

        // Read the input
        ImgData imgData;
        file::read_png_file(inputFile.c_str(), imgData);

        // Set the image quilting parameters and allocate the output
        imgData.output_h = outputScale * imgData.height;
        imgData.output_w = outputScale * imgData.width;
        imgData.block_h = imgData.height / blockDivisor;
        imgData.block_w = imgData.width / blockDivisor;
        imgData.AllocateOutput();

        // Write to the output file
        auto timingData = rdtsc_functional(imageQuiltingFunction, &imgData, seed);
        double performance = (double)timingData.flops / timingData.cycles;
        outputFile << "size=" << imgData.height << "x" << imgData.width;
        outputFile << ", n=" << imgData.height * imgData.width;
        outputFile << ", performance=" << performance;
        outputFile << ", flops=" << timingData.flops;
        outputFile << ", cycles=" << timingData.cycles << std::endl;

        // Clean up
        imgData.FreeInput();
        imgData.FreeOutput();
    }

    // Clean up
    outputFile.close();
    std::cout << "Done!" << std::endl << std::endl;
}