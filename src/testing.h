#pragma once

#include <cmath>
#include <filesystem>
#include <functional>
#include <iostream>
#include <map>
#include <tuple>
#include <vector>

#include "benchmarking/tsc_x86.h"
#include "image_quilting.h"
#include "png_reader.h"

// Defines for timing
#define CYCLES_REQUIRED 1e8
#define RDTSC_LATENCY 26
#define REP 1

class Testing {
public:
    // Construct a Testing object with the given input files
    Testing(int seed);
    // Seed getter and setter
    int GetSeed();
    void SetSeed(int seed);

    // TestFunction applies the image quilting algorithm to imgData
    typedef std::function<void(ImgData* imgData, int seed)> TestFunction;
    typedef std::function<void(ImgData* imgData)> ImgDataFunction;

    // Set the image quilting parameters function
    void SetParameterFunction(const ImgDataFunction& parameterFunction);
    // Run the image quilting algorithm on all the input files to generate the output files
    void GenerateOutputFiles();
    // Functional wrapper for the image quilting algorithm
    static double ImageQuiltingFunction(ImgData* imgData, int seed);
    // Functional wrapper for the empty image quilting algorithm
    static void EmptyImageQuiltingFunction(ImgData* imgData, int seed);
    // Register a function to test
    void RegisterTestFunction(const TestFunction& testFunction, std::string label);
    // Register a component function to test
    void RegisterComponentTestFunction(const TestFunction& baseFunction, const TestFunction& testFunction,
                                       std::string label);
    // Test the correctness of all the registered functions on every input
    void TestCorrectness();
    // Set the input for the TestCorrectnessAndTiming test suite
    void SetCorrectnessAndTimingInput(const std::string& label);
    // Test the correctness and timing of all the registered test functions against the specified input
    void TestCorrectnessAndTiming(bool stabilize);
    // Test the correctness and timing of all the registered component test functions against the specified input
    void TestComponentsTiming(bool stabilize);
    // Runs a given collection of functions to be timed by Intel Advisor
    void ComponentsTimingAdvisor();
    // Adds a function we wanted tested by Intel Advisor in ComponentsTimingAdvisor()
    void RegisterTestingComponentAdvisor(const TestFunction& testFunction, const std::string label);

private:
    const std::string inputDirectory = "./testing/input/";
    const std::string outputDirectory = "./testing/output/";
    int seed = 0;
    std::vector<std::filesystem::path> inputPaths;
    std::map<std::string, int> inputLabelsToIndices;
    std::vector<std::pair<TestFunction, std::string>> testFunctions;
    std::vector<std::tuple<TestFunction, TestFunction, std::string>> testComponentFunctions;
    std::vector<std::pair<TestFunction, std::string>> toTestAdvisor;
    ImgDataFunction parameterFunction;
    int testingInputIndex = 0;

    // Get the output file path
    std::string GetOutputfile(const std::filesystem::path& input);
    // Default image quilting algorithm parameters
    static void DefaultParameters(ImgData* imgData);
    // Call the parameter function to set the image quilting algorithm parameters
    void SetImageQuiltingParameters(ImgData* imgData);
    // Compute the l2 error between two images
    static double ComputeError(unsigned char** image0, unsigned char** image1, int height, int width);
    // Count the number of cycles of a testFunction call
    static double rdtsc(const TestFunction& testFunction, ImgData* inputData, int seed);
};