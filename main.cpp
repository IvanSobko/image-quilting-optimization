#include <cstdlib>
#include <string>
#include "ImageQuilting.h"
#include "PngReader.h"

#include "Testing.h"
#include "src/CompOverlapOptimiz.h"
#include "src/AdvanceAlgOptimiz.h"
#include "src/benchmarking/timing.h"
#include "Blocking.h"

// input parameters
std::string input_file = "./gallery/input0.png";
std::string output_file = "./gallery/output0.png";
ImgData img_data;

bool runTiming = false;
bool generate = false;
bool test = false;
bool testCorrectnessAndTiming = false;
bool advisor = true;
bool stabilize = false;

void parse_args(int argc, char* argv[]) {
    std::string delimiter = "=";
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        std::string par = arg.substr(0, arg.find(delimiter));
        std::string val = arg.substr(arg.find(delimiter) + 1, std::string::npos);
        if (par == "--input") {
            input_file = "./gallery/" + val;
        } else if (par == "--output") {
            output_file = "./gallery/" + val;
        } else if (par == "--width") {
            img_data.output_w = std::stoi(val);
        } else if (par == "--height") {
            img_data.output_h = std::stoi(val);
        } else if (par == "--blockW") {
            img_data.block_w = std::stoi(val);
        } else if (par == "--blockH") {
            img_data.block_h = std::stoi(val);
        } else if (par == "--timing") {
            runTiming = true;
        } else if (par == "--generate") {
            generate = true;
        } else if (par == "--test") {
            test = true;
        } else if (par == "--testCorrectnessAndTiming") {
            testCorrectnessAndTiming = true;
        } else if (par == "--stabilize") {
            stabilize = true;
        } else if (par == "--advisor") {
            advisor = true;
        } 
    }
}

// if invalid or were not set before
void set_default() {
    if (img_data.output_w == 0) {
        img_data.output_w = img_data.width * 2;
    }
    if (img_data.output_h == 0) {
        img_data.output_h = img_data.height * 2;
    }
    uint32_t max_block_w = img_data.width / 3;
    uint32_t max_block_h = img_data.height / 3;
    if (img_data.block_w == 0 || img_data.block_w > max_block_w) {
        img_data.block_w = img_data.width / 4;
    }
    if (img_data.block_h == 0 || img_data.block_h > max_block_h) {
        img_data.block_h = img_data.height / 4;
    }
}

int main(int argc, char* argv[]) {

    // Parse the command line arguments
    parse_args(argc, argv);

    // Benchmarking
    if (runTiming) {
        timing::run_timing(2);
        timing::run_timing(4);
        timing::run_timing(8);
    }
    // Generate output for testing
    else if (generate) {
        Testing testing = Testing(2);
        testing.SetParameterFunction(AdvanceAlgOptimiz::CustomParameters);
        testing.GenerateOutputFiles();
    }
    // Test the correctness of our base implementation
    else if (test) {
        Testing testing = Testing(2);
        testing.SetParameterFunction(AdvanceAlgOptimiz::CustomParameters);
        //testing.RegisterTestFunction(Testing::ImageQuiltingFunction, "default");
        //testing.RegisterTestFunction(AdvanceAlgOptimiz::StdC_KUnroll_BoundsRefactor, "refactor");
        //testing.RegisterTestFunction(AdvanceAlgOptimiz::StdC_KUnroll_BoundsRefactor_LoopReorder, "loop reorder");
        testing.RegisterTestFunction(AdvanceAlgOptimiz::StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking32, "blocking 32x32");
        testing.RegisterTestFunction(AdvanceAlgOptimiz::StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking48, "blocking 48x48");
        testing.RegisterTestFunction(AdvanceAlgOptimiz::StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking64, "blocking 64x64");
        testing.RegisterTestFunction(AdvanceAlgOptimiz::StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking96, "blocking 96x96");
        testing.RegisterTestFunction(AdvanceAlgOptimiz::StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking128, "blocking 128x128");
        testing.TestCorrectness();
    }
    // Test the correctness and timing of the variants of our implementation
    else if (testCorrectnessAndTiming) {
        Testing testing = Testing(2);
        testing.SetParameterFunction(AdvanceAlgOptimiz::CustomParameters);
        testing.SetCorrectnessAndTimingInput("input_256x256.png");

        testing.RegisterTestFunction(Testing::ImageQuiltingFunction, "default");
        testing.RegisterTestFunction(AdvanceAlgOptimiz::StdC_KUnroll_BoundsRefactor, "refactor");
        testing.RegisterTestFunction(AdvanceAlgOptimiz::StdC_KUnroll_BoundsRefactor_LoopReorder, "loop reorder");
        testing.RegisterTestFunction(AdvanceAlgOptimiz::StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking32, "blocking 32x32");
        testing.RegisterTestFunction(AdvanceAlgOptimiz::StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking48, "blocking 48x48");
        testing.RegisterTestFunction(AdvanceAlgOptimiz::StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking64, "blocking 64x64");
        testing.RegisterTestFunction(AdvanceAlgOptimiz::StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking96, "blocking 96x96");
        testing.RegisterTestFunction(AdvanceAlgOptimiz::StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking128, "blocking 128x128");

        std::cout << std::endl;
        testing.TestCorrectnessAndTiming(stabilize);
    }
    else if (advisor) {
        // Read the input and allocate output
        file::read_png_file("./testing/input/input_256x256.png", img_data);
        AdvanceAlgOptimiz::CustomParameters(&img_data);
        img_data.AllocateOutput();

        int seed = 2;

        // Default image quilting algorithm
        Testing::ImageQuiltingFunction(&img_data, seed);

        // Std c, k unroll, and bounds refactor optimizations
        AdvanceAlgOptimiz::StdC_KUnroll_BoundsRefactor(&img_data, seed);

        // Std c, k unroll, bounds refactor, and loop reorder optimizations
        AdvanceAlgOptimiz::StdC_KUnroll_BoundsRefactor_LoopReorder(&img_data, seed);

        // Std c, k unroll, bounds refactor, loop reorder, and blocking 32x23 optimizations
        AdvanceAlgOptimiz::StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking32(&img_data, seed);

        // Std c, k unroll, bounds refactor, loop reorder, and blocking 48x48 optimizations
        AdvanceAlgOptimiz::StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking48(&img_data, seed);

        // Std c, k unroll, bounds refactor, loop reorder, and blocking 64x64 optimizations
        AdvanceAlgOptimiz::StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking64(&img_data, seed);

        // Std c, k unroll, bounds refactor, loop reorder, and blocking 96x96 optimizations
        AdvanceAlgOptimiz::StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking96(&img_data, seed);

        // Std c, k unroll, bounds refactor, loop reorder, and blocking 128x128 optimizations
        AdvanceAlgOptimiz::StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking128(&img_data, seed);

        // Clean up
        img_data.FreeOutput();
        img_data.FreeInput();
    }
    // Main
    else {
        // Read the input data
        file::read_png_file(input_file.c_str(), img_data);
        set_default();

        // Allocate the output data and run the image quilting algorithm
        img_data.AllocateOutput();
        ImageQuilting imageQuilting(&img_data);
        imageQuilting.Synthesis();

        // Write the output file and free the members of img_data
        file::write_png_file(output_file.c_str(), img_data);
        img_data.FreeOutput();
        img_data.FreeInput();
    }
    return EXIT_SUCCESS;
}