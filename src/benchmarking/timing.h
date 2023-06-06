#pragma once

#include "../ImageQuilting.h"

#include <string>
#include <vector>
#include <functional>
#include <iostream>
#include <fstream>
#include <chrono>
#include <ctime>
#include <filesystem>

namespace timing {

// Original timing
void run_timing(int inputBlockRatio);
double rdtsc(ImageQuilting* quilting);
std::vector<std::string> read_files(const std::string &directory, const std::string &filename_filter);

// Function that runs the image quilting variant on the specified input and returns the number of flops
typedef std::function<double(ImgData* imgData, int seed)> ImageQuiltingFunction;
double EmptyImageQuiltingFunction(ImgData* imgData, int seed);

// Function that sets the image quilting parameters
typedef std::function<void(ImgData* imgData)> ParameterFunction;

// Struct to keep track of the cycles and flops for a given function call
struct TimingData {
    double cycles;
    double flops;
    TimingData(double cycles, double flops) : cycles(cycles), flops(flops) {}
};
// Call the image quilting variant on the specified input and return the number of cycles and flops
TimingData rdtsc_functional(const ImageQuiltingFunction & imageQuiltingFunction, ImgData* imgData, int seed);

// Run the timing for the image quilting variant and write the results to a file with the given label (and current date and time)
void run_timing_functional(
    const std::string & label,
    const std::string & inputDirectory, const std::string & filenameFilter, const std::string & outputDirectory,
    const ImageQuiltingFunction & imageQuiltingFunction, int outputScale, int blockDivisor);

}  // namespace timing