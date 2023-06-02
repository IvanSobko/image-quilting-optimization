#include <cstdlib>
#include <string>
#include "ImageQuilting.h"
#include "PngReader.h"

#include "Testing.h"
#include "src/CompOverlapOptimiz.h"
#include "src/AdvanceAlgOptimiz.h"
#include "src/benchmarking/timing.h"
#include "Blocking.h"
#include "AdvisorHelper.h"

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
        testing.GenerateOutputFiles();
    }
    // Test the correctness of our base implementation
    else if (test) {
        Testing testing = Testing(0);
        testing.RegisterTestFunction(Testing::ImageQuiltingFunction, "default");
        testing.TestCorrectness();
    }
    // Test the correctness and timing of the variants of our implementation
    else if (testCorrectnessAndTiming) {
        Testing testing = Testing(2);

        testing.RegisterTestFunction(Testing::ImageQuiltingFunction, "default");
        // testing.RegisterTestFunction(CompOverlapOptimiz::BasicOpt, "compBasic");
        testing.RegisterTestFunction(CompOverlapOptimiz::AlgOpt, "compBasic+AlgImpr");
        // testing.RegisterTestFunction(CompOverlapOptimiz::UnrollOpt, "compBasic+AlgImpr+Unroll");
        testing.RegisterTestFunction(CompOverlapOptimiz::UnrollChnls, "compBasic+AlgImpr+UnrollChnls");
        testing.RegisterTestFunction(AdvanceAlgOptimiz::DividedFuncOpt, "Unroll+DividedFunctions");
        // testing.RegisterTestFunction(CompOverlapOptimiz::UnrollMaxOpt, "compBasic+AlgImpr+UnrollTheoreticalMax");
        // testing.RegisterTestFunction(CompOverlapOptimiz::VectorizeOpt, "compBasic+AlgImpr+Unroll+Vectorize");

        std::cout << std::endl;
        testing.TestCorrectnessAndTiming(stabilize);

        // testing.RegisterComponentTestFunction(CompOverlapOptimiz::BaseComponent,
        //                                       CompOverlapOptimiz::BaseComponent, "default");

        // testing.RegisterComponentTestFunction(CompOverlapOptimiz::BaseComponent,
        //                                       CompOverlapOptimiz::BasicOptComponent, "compBasic");

        // testing.RegisterComponentTestFunction(CompOverlapOptimiz::BaseComponent,
        //                                       CompOverlapOptimiz::AlgoOptComponent, "compBasic+AlgImpr");

        // testing.RegisterComponentTestFunction(CompOverlapOptimiz::BaseComponent,
        //                                       CompOverlapOptimiz::UnrollOptComponent,
        //                                       "compBasic+AlgImpr+Unroll");

        // testing.RegisterComponentTestFunction(CompOverlapOptimiz::BaseComponent,
        //                                       CompOverlapOptimiz::UnrollMaxOptComponent,
        //                                       "compBasic+AlgImpr+UnrollTheoreticalMax");

        // testing.RegisterComponentTestFunction(CompOverlapOptimiz::BaseComponent,
        //                                       CompOverlapOptimiz::VectorizeOptComponent,
        //                                       "compBasic+AlgImpr+Unroll+Vectorize");

        std::cout << std::endl;
        testing.TestComponentsTiming(stabilize);
    } else if (advisor) {
        file::read_png_file(input_file.c_str(), img_data);
        set_default();

        // Allocate the output data and run the image quilting algorithm
        img_data.AllocateOutput();
        Advisor::baseline(&img_data, 0);
        Advisor::basicOpt(&img_data, 0);
        Advisor::unrollChnls(&img_data, 0);
        Advisor::unrollMemory(&img_data, 0);
        // Advisor::vectorize(&img_data, 0);
        // Advisor::block(&img_data, 0);
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