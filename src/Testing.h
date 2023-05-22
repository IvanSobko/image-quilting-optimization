#pragma once

#include <iostream>
#include <vector>
#include <filesystem>
#include <functional>
#include <cmath>
#include <tuple>

#include "PngReader.h"
#include "ImageQuilting.h"
#include "benchmarking/tsc_x86.h"

// Defines for timing
#define CYCLES_REQUIRED 1e8
#define RDTSC_LATENCY 26
#define REP 5

class Testing {
   public:
    // Construct a Testing object with the given input files
    Testing(int seed);
    // Seed getter and setter
    int GetSeed();
    void SetSeed(int seed);

    // TestFunction applies the image quilting algorithm to imgData
    typedef std::function<void(ImgData* imgData, int seed)> TestFunction;

    // Run the image quilting algorithm on all the input files to generate the output files
    void GenerateOutputFiles();
    // Functional wrapper for the image quilting algorithm
    static void ImageQuiltingFunction(ImgData* imgData, int seed);
    // Functional wrapper for the empty image quilting algorithm
    static void EmptyImageQuiltingFunction(ImgData* imgData, int seed);
    // Register a function to test
    void RegisterTestFunction(const TestFunction & testFunction, std::string label);
    // Register a component function to test
    void RegisterComponentTestFunction(
        const TestFunction & baseFunction, const TestFunction & testFunction, std::string label);
    // Test the correctness of all the registered functions
    void TestCorrectness();
    // Test the correctness and timing of all the registered test functions
    void TestCorrectnessAndTiming(bool stabilize);
    // Test the correctness and timing of all the registered component test functions
    void TestComponentsTiming(bool stabilize);

   private:
    const std::string inputDirectory = "./testing/input/";
    const std::string outputDirectory = "./testing/output/";
    int seed = 0;
    std::vector<std::filesystem::path> inputPaths;
    std::vector<std::pair<TestFunction, std::string>> testFunctions;
    std::vector<std::tuple<TestFunction, TestFunction, std::string>> testComponentFunctions;

    // Get the output file path
    std::string GetOutputfile(const std::filesystem::path & input);
    // Set the image quilting algorithm parameters
    static void SetImageQuiltingParameters(ImgData* imgData);
    // Compute the l2 error between two images
    static double ComputeError(unsigned char ** image0, unsigned char ** image1, int height, int width);
    // Count the number of cycles of a testFunction call
    static double rdtsc(const TestFunction& testFunction, ImgData* inputData, int seed);
};