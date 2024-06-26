#include "testing.h"

// Construct a Testing object with the given input files
Testing::Testing(int seed) {
    // Set the seed
    this->seed = seed;

    // Initialize the parameter function
    parameterFunction = DefaultParameters;

    // https://stackoverflow.com/questions/612097/how-can-i-get-the-list-of-files-in-a-directory-using-c-or-c
    // Construct the vector of input files
    for (const std::filesystem::path& inputPath : std::filesystem::directory_iterator(inputDirectory)) {
        inputPaths.push_back(inputPath);
        inputLabelsToIndices.emplace(inputPath.filename().string(), inputPaths.size() - 1);
    }
}

// Seed getter
int Testing::GetSeed() {
    return seed;
}

// Seed setter
void Testing::SetSeed(int seed) {
    this->seed = seed;
}

// Get the output file path
std::string Testing::GetOutputfile(const std::filesystem::path& input) {
    return outputDirectory + input.stem().string() + input.extension().string();
}

// Default image quilting algorithm parameters
void Testing::DefaultParameters(ImgData* imgData) {
    imgData->output_h = 2 * imgData->height;
    imgData->output_w = 2 * imgData->width;
    imgData->block_h = imgData->height / 4;
    imgData->block_w = imgData->width / 4;
}

// Call the parameter function to set the image quilting algorithm parameters
void Testing::SetImageQuiltingParameters(ImgData* imgData) {
    parameterFunction(imgData);
}

// Set the image quilting parameters function
void Testing::SetParameterFunction(const ImgDataFunction& parameterFunction) {
    this->parameterFunction = parameterFunction;
}

// Run the image quilting algorithm on all the input files to generate the output files
void Testing::GenerateOutputFiles() {

    std::cout << "Input files:" << std::endl;
    for (const auto& input : inputPaths) {
        std::cout << input << std::endl;
    }

    if (!std::filesystem::exists(outputDirectory)) {
        std::filesystem::create_directory(outputDirectory);
    }

    std::cout << "Output files:" << std::endl;
    for (const auto& input : inputPaths) {
        // Construct the output path
        auto output = GetOutputfile(input);

        // Run the image quilting algorithm and write the output
        ImgData imgData;
        file::read_png_file(input.string().c_str(), imgData);

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
double Testing::ImageQuiltingFunction(ImgData* imgData, int seed) {
    ImageQuilting imageQuilting(imgData);
    imageQuilting.Synthesis(seed);
    return static_cast<double>(imageQuilting.getFlopCount());
}

// Functional wrapper for the empty image quilting algorithm
void Testing::EmptyImageQuiltingFunction(ImgData* imgData, int seed) {}

// Register a function to test
void Testing::RegisterTestFunction(const TestFunction& testFunction, const std::string label) {
    testFunctions.emplace_back(testFunction, label);
    std::cout << "Registered test function: " << label << std::endl;
}

// Register a component function to test
void Testing::RegisterComponentTestFunction(const TestFunction& baseFunction,
                                            const TestFunction& testFunction, const std::string label) {
    testComponentFunctions.emplace_back(baseFunction, testFunction, label);
    std::cout << "Registered component function: " << label << std::endl;
}

// Compute the l2 error between two images
double Testing::ComputeError(unsigned char** image0, unsigned char** image1, int height, int width) {
    double l2norm = 0;
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            for (int k = 0; k < CHANNEL_NUM; k++) {
                double x0 = image0[i][CHANNEL_NUM * j + k];
                double x1 = image1[i][CHANNEL_NUM * j + k];
                double norm = x0 - x1;
                l2norm += norm * norm;
            }
        }
    }
    return std::sqrt(l2norm);
}

// Test the correctness of all the registered functions
void Testing::TestCorrectness() {
    for (const auto& pair : testFunctions) {
        TestFunction testFunction;
        std::string label;
        std::tie(testFunction, label) = pair;

        bool correct = true;
        std::cout << "Testing function: " << label << std::endl;

        // Test against all output images
        for (const auto& input : inputPaths) {
            // Read the input image
            ImgData inputImgData;
            file::read_png_file(input.string().c_str(), inputImgData);

            // Read the already computed output image
            auto output = GetOutputfile(input);
            ImgData outputImgData;
            file::read_png_file(output.c_str(), outputImgData);

            // Run the test function
            SetImageQuiltingParameters(&inputImgData);
            inputImgData.AllocateOutput();
            testFunction(&inputImgData, seed);

            // Compute the error
            double error = ComputeError(inputImgData.output_d, outputImgData.data, outputImgData.height,
                                        outputImgData.width);

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

        if (correct) {
            std::cout << label << " is correct" << std::endl;
        } else {
            std::cout << label << " is incorrect" << std::endl;
        }
    }
}

// Set the input for the TestCorrectnessAndTiming test suite
void Testing::SetCorrectnessAndTimingInput(const std::string& label) {
    if (auto search = inputLabelsToIndices.find(label); search != inputLabelsToIndices.end()) {
        testingInputIndex = search->second;
        std::cout << "Set TestCorrectnessAndTiming input to \"" << label << "\"" << std::endl;
    } else {
        std::cout << "Input with label \"" << label << "\" not found" << std::endl;
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
void Testing::TestCorrectnessAndTiming(bool stabilize) {
    // Use the specified input
    auto input = inputPaths[testingInputIndex];

    // Read the input image
    ImgData inputImgData;
    file::read_png_file(input.string().c_str(), inputImgData);

    // Set the image quilting parameters and allocate the output image
    SetImageQuiltingParameters(&inputImgData);

    // Read the already computed output image
    auto output = GetOutputfile(input);
    ImgData outputImgData;
    file::read_png_file(output.c_str(), outputImgData);

    // Compute the base number of cycles
    if (stabilize) {
        std::cout << "Stabilizing the rdtsc function calls..." << std::endl;
        int baselineIterations = 5;  // TODO magic number
        double speedup = 0.0;
        do {
            double computedCycles[baselineIterations];
            for (int j = 0; j < baselineIterations; j++) {
                inputImgData.AllocateOutput();
                computedCycles[j] = rdtsc(ImageQuiltingFunction, &inputImgData, seed);
                inputImgData.FreeOutput();
            }
            double averageCycles = 0;
            for (int j = 0; j < baselineIterations; j++)
                averageCycles += computedCycles[j];
            averageCycles /= baselineIterations;
            speedup = averageCycles / computedCycles[baselineIterations - 1];
            std::cout << "Speedup: " << speedup << std::endl;
        } while (std::abs(speedup - 1.0) > .01);  // TODO magic number
        std::cout << "Stabilized! Computing the base number of cycles..." << std::endl;
    }

    inputImgData.AllocateOutput();
    double baseCycles = rdtsc(ImageQuiltingFunction, &inputImgData, seed);
    inputImgData.FreeOutput();
    std::cout << "Base number of cycles: " << baseCycles << std::endl << std::endl;

    // Test the correctness and timing of all the registered test functions
    for (const auto& pair : testFunctions) {
        TestFunction testFunction;
        std::string label;
        std::tie(testFunction, label) = pair;

        std::cout << "Testing function: " << label << std::endl;

        // Run the test function
        inputImgData.AllocateOutput();
        testFunction(&inputImgData, seed);

        // Compute the error
        double error = ComputeError(inputImgData.output_d, outputImgData.data, outputImgData.height,
                                    outputImgData.width);

        // Print the correctness
        if (std::abs(error) < 1e-8f) {
            std::cout << label << " is correct" << std::endl;
        } else {
            std::cout << label << " is incorrect" << std::endl;
        }

        // Compute the timing
        double cycles = rdtsc(testFunction, &inputImgData, seed);

        // Print the timing
        std::cout << "Speedup: " << baseCycles / cycles << std::endl;
        std::cout << std::endl;

        // Clean up
        inputImgData.FreeOutput();
    }

    // Clean up
    inputImgData.FreeInput();
    outputImgData.FreeInput();
    outputImgData.FreeOutput();
}

// Test the correctness and timing of all the registered component test functions
void Testing::TestComponentsTiming(bool stabilize) {
    // Seed the random number generator with the system time
    srand(time(nullptr));

    // Generate a random image
    ImgData imgData;
    imgData.height = 256;
    imgData.width = 256;
    SetImageQuiltingParameters(&imgData);
    imgData.AllocateInput();
    imgData.RandomizeInput();
    imgData.AllocateOutput();

    // Test the timing of all the registered component functions
    for (const auto& tuple : testComponentFunctions) {
        TestFunction baseFunction;
        TestFunction testFunction;
        std::string label;
        std::tie(baseFunction, testFunction, label) = tuple;

        std::cout << "Testing function: " << label << std::endl;

        // Compute the speedup
        // TODO add stabilization
        double baseCycles = rdtsc(baseFunction, &imgData, seed);
        double cycles = rdtsc(testFunction, &imgData, seed);
        std::cout << "Speedup: " << baseCycles / cycles << std::endl;
        std::cout << std::endl;
    }

    imgData.FreeInput();
    imgData.FreeOutput();
}

// Runs a given collection of functions to be timed by Intel Advisor
void Testing::ComponentsTimingAdvisor() {
    // Seed the random number generator with the system time
    srand(time(nullptr));

    // Generate a random image
    ImgData imgData;
    imgData.height = 256;
    imgData.width = 256;
    SetImageQuiltingParameters(&imgData);
    imgData.AllocateInput();
    imgData.RandomizeInput();
    imgData.AllocateOutput();

    // Test the timing of all the registered component functions
    for (const auto& pair : toTestAdvisor) {
        TestFunction testFunction;
        std::string label;
        std::tie(testFunction, label) = pair;

        std::cout << "Running function: " << label << std::endl;

        // Time the given function (just to run it)
        // TODO add stabilization
        double cycles = rdtsc(testFunction, &imgData, seed);
        std::cout << "Number of cycles: " << cycles << std::endl;
        std::cout << std::endl;
    }

    imgData.FreeInput();
    imgData.FreeOutput();
}

void Testing::RegisterTestingComponentAdvisor(const TestFunction& testFunction, const std::string label) {
    toTestAdvisor.emplace_back(testFunction, label);
    std::cout << "Registered Advisor component function: " << label << std::endl;
}
