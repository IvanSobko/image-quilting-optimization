#include "Testing.h"

// Construct a Testing object with the given input files
Testing::Testing(int seed) {
    // Set the seed
    this->seed = seed;

    // https://stackoverflow.com/questions/612097/how-can-i-get-the-list-of-files-in-a-directory-using-c-or-c
    // Construct the vector of input files
    std::cout << "Input files:" << std::endl;
    for (const auto & inputPath : std::filesystem::directory_iterator(inputDirectory)){
        std::cout << inputPath.path() << std::endl;
        inputPaths.push_back(inputPath);
    }
}

// Seed getter
int Testing::GetSeed() { return seed; }

// Seed setter
void Testing::SetSeed(int seed) { this->seed = seed; }

// Get the output file path
std::string Testing::GetOutputfile(const std::filesystem::path & input) {
    return outputDirectory + input.stem().string() + input.extension().string();
}

// Run the image quilting algorithm on all the input files to generate the output files
void Testing::GenerateOutputFiles() {

    std::cout << "Output files:" << std::endl;
    for (const auto & input : inputPaths){
        // Construct the output path
        auto output = GetOutputfile(input);

        // Run the image quilting algorithm and write the output
        ImgData imgData;
        file::read_png_file(input.c_str(), imgData);

        // Image quilting parameters
        imgData.output_h = 2 * imgData.height;
        imgData.output_w = 2 * imgData.width;
        imgData.block_h = imgData.height / 2;
        imgData.block_w = imgData.width / 2;

        // Image quilting algorithm
        ImageQuiltingFunction(&imgData, seed);

        // Write the output
        file::write_png_file(output.c_str(), imgData);
        std::cout << output << std::endl;
    }
}

// Functional wrapper for the image quilting algorithm
void Testing::ImageQuiltingFunction(ImgData* imgData, int seed) {
    ImageQuilting imageQuilting(imgData);
    imageQuilting.Synthesis(seed);
}

// Register a function to test
void Testing::RegisterTestFunction(const TestFunction& testFunction, const std::string label) {
    testFunctions.emplace_back(testFunction, label);
    std::cout << "Registered function: " << label << std::endl;
}

// Compute the l2 error between two images
double Testing::ComputeError(unsigned char** image0, unsigned char** image1, int height, int width) {
    double l2norm = 0;
    for (int i = 0; i < height; i++){
        for (int j = 0; j < width; j++){
            for (int k = 0; k < CHANNEL_NUM; k++){
                double x0 = image0[i][CHANNEL_NUM*j+k];
                double x1 = image0[i][CHANNEL_NUM*j+k];
                double norm = x0 - x1;
                l2norm += norm*norm;
            }
        }
    }
    return std::sqrt(l2norm);
}

// Test all the registered functions
void Testing::TestRegisteredTestFunctions() {
    for (const auto & pair : testFunctions) {
        TestFunction testFunction;
        std::string label;
        std::tie(testFunction, label) = pair;

        bool correct = true;
        std::cout << "Testing function: " << label << std::endl;

        // Test against all output images
        for (const auto & input : inputPaths) {
            // Read the input and output images
            ImgData inputImgData;
            file::read_png_file(input.c_str(), inputImgData);
            auto output = GetOutputfile(input);
            ImgData outputImgData;
            file::read_png_file(output.c_str(), outputImgData);

            // Run the test function
            inputImgData.output_h = 2 * inputImgData.height;
            inputImgData.output_w = 2 * inputImgData.width;
            inputImgData.block_h = inputImgData.height / 2;
            inputImgData.block_w = inputImgData.width / 2;
            testFunction(&inputImgData, seed);

            // Compute the error
            double error = ComputeError(
                inputImgData.output_d, outputImgData.data,
                outputImgData.height, outputImgData.width);

            // Print the error
            std::cout << input << ", error: " << error << std::endl;
            if (error != 0) {
                correct = false;
                break;
            }
        }

        if (correct) std::cout << label << " is correct" << std::endl;
        else std::cout << label << " is incorrect" << std::endl;
    }
}

