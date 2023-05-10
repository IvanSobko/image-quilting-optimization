#include "Testing.h"

// Construct a Testing object with the given input files
Testing::Testing(int seed) {
    // Set the seed
    this->seed = seed;

    // https://stackoverflow.com/questions/612097/how-can-i-get-the-list-of-files-in-a-directory-using-c-or-c
    // Construct the vector of input files
    for (const auto & inputPath : std::filesystem::directory_iterator(inputDirectory)){
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

// Set the image quilting algorithm parameters
void Testing::SetImageQuiltingParameters(ImgData* imgData) {
    imgData->output_h = 2 * imgData->height;
    imgData->output_w = 2 * imgData->width;
    imgData->block_h = imgData->height / 4;
    imgData->block_w = imgData->width / 4;
}

// Run the image quilting algorithm on all the input files to generate the output files
void Testing::GenerateOutputFiles() {

    std::cout << "Input files:" << std::endl;
    for (const auto & input : inputPaths) {
        std::cout << input << std::endl;
    }

    std::cout << "Output files:" << std::endl;
    for (const auto & input : inputPaths){
        // Construct the output path
        auto output = GetOutputfile(input);

        // Run the image quilting algorithm and write the output
        ImgData imgData;
        file::read_png_file(input.c_str(), imgData);

        // Image quilting algorithm
        SetImageQuiltingParameters(&imgData);
        imgData.AllocateOutput();
        ImageQuiltingFunction(&imgData, seed);

        // Write the output
        file::write_png_file(output.c_str(), imgData);
        std::cout << output << std::endl;

        // Cleanup
        imgData.FreeInput();
        imgData.FreeOutput();
    }
}

// Functional wrapper for the image quilting algorithm
void Testing::ImageQuiltingFunction(ImgData* imgData, int seed) {
    ImageQuilting imageQuilting(imgData);
    imageQuilting.Synthesis(seed);
}

// Functional wrapper for the empty image quilting algorithm
void Testing::EmptyImageQuiltingFunction(ImgData* imgData, int seed) { }

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
                double x1 = image1[i][CHANNEL_NUM*j+k];
                double norm = x0 - x1;
                l2norm += norm*norm;
            }
        }
    }
    return std::sqrt(l2norm);
}

// Test the correctness of all the registered functions
void Testing::TestCorrectness() {
    for (const auto & pair : testFunctions) {
        TestFunction testFunction;
        std::string label;
        std::tie(testFunction, label) = pair;

        bool correct = true;
        std::cout << "Testing function: " << label << std::endl;

        // Test against all output images
        for (const auto & input : inputPaths) {
            // Read the input image
            ImgData inputImgData;
            file::read_png_file(input.c_str(), inputImgData);

            // Read the already computed output image
            auto output = GetOutputfile(input);
            ImgData outputImgData;
            file::read_png_file(output.c_str(), outputImgData);

            // Run the test function
            SetImageQuiltingParameters(&inputImgData);
            inputImgData.AllocateOutput();
            testFunction(&inputImgData, seed);

            // Compute the error
            double error = ComputeError(
                inputImgData.output_d, outputImgData.data,
                outputImgData.height, outputImgData.width);

            // Print the error
            std::cout << input << ", error: " << error << std::endl;

            // Clean up
            inputImgData.FreeOutput();
            outputImgData.FreeInput();
            outputImgData.FreeOutput();

            // Terminate early
            if (error != 0) {
                correct = false;
                break;
            }
        }

        if (correct) std::cout << label << " is correct" << std::endl;
        else std::cout << label << " is incorrect" << std::endl;
    }
}

// Count the number of cycles of a testFunction call
double Testing::rdtsc(const TestFunction& testFunction, ImgData* inputData, int seed) {
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
            testFunction(inputData, seed);
        }
        end = stop_tsc(start);

        cycles = (double)end;
        multiplier = (CYCLES_REQUIRED) / (cycles);

    } while (multiplier > 2);

    // Actual performance measurements repeated REP times.
    // We simply store all results and compute medians during post-processing.
    double total_cycles = 0;
    for (size_t j = 0; j < REP; j++) {
        start = start_tsc();
        for (size_t i = 0; i < num_runs; ++i) {
            testFunction(inputData, seed);
        }
        end = stop_tsc(start) - RDTSC_LATENCY;

        cycles = (double)end / (double)num_runs;
        total_cycles += cycles;
    }
    total_cycles /= REP;
    return total_cycles;
}

// Test the correctness and timing of all the registered test functions
void Testing::TestCorrectnessAndTiming() {
    // Use the first input
    auto input = inputPaths[0];

    // Read the input image
    ImgData inputImgData;
    file::read_png_file(input.c_str(), inputImgData);

    // Set the image quilting parameters and allocate the output image
    SetImageQuiltingParameters(&inputImgData);
    inputImgData.AllocateOutput();

    // Read the already computed output image
    auto output = GetOutputfile(input);
    ImgData outputImgData;
    file::read_png_file(output.c_str(), outputImgData);

    // Compute the base number of cycles
    double baseCycles = rdtsc(ImageQuiltingFunction, &inputImgData, seed);

    // Test the correctness and timing of all the registered test functions
    for (const auto & pair : testFunctions) {
        TestFunction testFunction;
        std::string label;
        std::tie(testFunction, label) = pair;

        std::cout << "Testing function: " << label << std::endl;

        // Run the test function
        inputImgData.AllocateOutput();
        testFunction(&inputImgData, seed);

        // Compute the error
        double error = ComputeError(
            inputImgData.output_d, outputImgData.data,
            outputImgData.height, outputImgData.width);

        // Print the correctness
        if (error == 0)  std::cout << label << " is correct" << std::endl;
        else std::cout<< label << " is incorrect" << std::endl;

        // Compute the timing
        double cycles = rdtsc(testFunction, &inputImgData, seed);

        // Print the timing
        std::cout << "Speedup: " << cycles / baseCycles << std::endl;
        std::cout << std::endl;

        // Clean up
        inputImgData.FreeOutput();
    }

    // Clean up
    inputImgData.FreeInput();
    outputImgData.FreeInput();
    outputImgData.FreeOutput();
}