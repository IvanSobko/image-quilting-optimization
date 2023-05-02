#pragma once

#include <iostream>
#include <vector>
#include <filesystem>
#include <functional>
#include <cmath>

#include "PngReader.h"
#include "ImageQuilting.h"

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
    // Register a function to test
    void RegisterTestFunction(const TestFunction& testFunction, std::string label);
    // Test all the registered functions
    void TestRegisteredTestFunctions();

   private:
    const std::string inputDirectory = "./testing/input/";
    const std::string outputDirectory = "./testing/output/";
    int seed = 0;
    std::vector<std::filesystem::path> inputPaths;
    std::vector<std::pair<TestFunction, std::string>> testFunctions;

    // Get the output file path
    std::string GetOutputfile(const std::filesystem::path & input);
    // Compute the l2 error between two images
    static double ComputeError(unsigned char ** image0, unsigned char ** image1, int height, int width);
};