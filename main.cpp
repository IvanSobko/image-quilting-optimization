#include <cstdlib>
#include <string>
#include "image_quilting.h"
#include "png_reader.h"

#include "blocking.h"
#include "src/advance_alg_optimiz.h"
#include "src/benchmarking/timing.h"
#include "src/comp_overlap_optimiz.h"
#include "testing.h"

// input parameters
std::string input_file = "./gallery/input0.png";
std::string output_file = "./gallery/output0.png";
ImgData img_data;

bool runTiming = false;
bool generate = false;
bool test = false;
bool testCorrectnessAndTiming = false;
bool stabilize = false;
bool advisor = false;
bool timingFunctional = false;

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
        } else if (par == "--timingFunctional") {
            timingFunctional = true;
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
        testing.RegisterTestFunction(CompOverlapOptimiz::BasicOpt, "compBasic");
        testing.RegisterTestFunction(CompOverlapOptimiz::AlgOpt, "compBasic+AlgImpr");
        testing.RegisterTestFunction(CompOverlapOptimiz::UnrollOpt, "compBasic+AlgImpr+Unroll");
        testing.RegisterTestFunction(AdvanceAlgOptimiz::StdC_KUnroll_BoundsRefactor,
                                     "Unroll+DividedFunctions");
        testing.RegisterTestFunction(CompOverlapOptimiz::UnrollMaxOpt,
                                     "compBasic+AlgImpr+UnrollTheoreticalMax");
#ifdef __AVX2__
        testing.RegisterTestFunction(CompOverlapOptimiz::VectorizeOpt, "compBasic+AlgImpr+Unroll+Vectorize");
#endif

        std::cout << std::endl;
        testing.TestCorrectnessAndTiming(stabilize);

        testing.RegisterComponentTestFunction(CompOverlapOptimiz::BaseComponent,
                                              CompOverlapOptimiz::BaseComponent, "default");

        testing.RegisterComponentTestFunction(CompOverlapOptimiz::BaseComponent,
                                              CompOverlapOptimiz::BasicOptComponent, "compBasic");

        testing.RegisterComponentTestFunction(CompOverlapOptimiz::BaseComponent,
                                              CompOverlapOptimiz::AlgoOptComponent, "compBasic+AlgImpr");

        testing.RegisterComponentTestFunction(CompOverlapOptimiz::BaseComponent,
                                              CompOverlapOptimiz::UnrollOptComponent,
                                              "compBasic+AlgImpr+Unroll");

        testing.RegisterComponentTestFunction(CompOverlapOptimiz::BaseComponent,
                                              CompOverlapOptimiz::UnrollMaxOptComponent,
                                              "compBasic+AlgImpr+UnrollTheoreticalMax");

#ifdef __AVX2__
        testing.RegisterComponentTestFunction(CompOverlapOptimiz::BaseComponent,
                                              CompOverlapOptimiz::VectorizeOptComponent,
                                              "compBasic+AlgImpr+Unroll+Vectorize");
#endif
        std::cout << std::endl;
        testing.TestComponentsTiming(stabilize);
    }
    // Only runs specific functions to give to Intel Advisor
    else if (advisor) {
        Testing testing = Testing(0);

        testing.RegisterTestingComponentAdvisor(CompOverlapOptimiz::BaseComponent, "default");
        testing.RegisterTestingComponentAdvisor(CompOverlapOptimiz::BasicOptComponent, "compBasic");
        testing.RegisterTestingComponentAdvisor(CompOverlapOptimiz::UnrollOptComponent,
                                                "compBasic+AlgImpr+Unroll");
        testing.RegisterTestingComponentAdvisor(CompOverlapOptimiz::UnrollMaxOptComponent,
                                                "compBasic+AlgImpr+UnrollTheoreticalMax");
#ifdef __AVX2__
        testing.RegisterTestingComponentAdvisor(CompOverlapOptimiz::VectorizeOptComponent,
                                                "compBasic+AlgImpr+Unroll+Vectorize");
#endif
        testing.ComponentsTimingAdvisor();
    }
    // Run the functional timing code
    else if (timingFunctional) {
        // Define the set of image quilting functions
        std::vector<std::pair<std::string, timing::ImageQuiltingFunction>> imageQuiltingFunctions = {
            {"default", Testing::ImageQuiltingFunction},
            {"StdC_Algorithm", CompOverlapOptimiz::AlgOpt},
            {"StdC_Algorithm_ChannelsUnroll", CompOverlapOptimiz::UnrollChnls},
#ifdef __AVX2__
            {"StdC_Algorithm_Vectorize", CompOverlapOptimiz::VectorizeOpt},
#endif
            {"StdC_KUnroll_BoundsRefactor", AdvanceAlgOptimiz::StdC_KUnroll_BoundsRefactor},

            // TODO: do we need to test all possible block sizes? Maybe we know that some of them are not beneficial?
            {"StdC_KUnroll_BoundsRefactor_LoopReorder",
             AdvanceAlgOptimiz::StdC_KUnroll_BoundsRefactor_LoopReorder},
            {"StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking32",
             AdvanceAlgOptimiz::StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking32},
            {"StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking48",
             AdvanceAlgOptimiz::StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking48},
            {"StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking64",
             AdvanceAlgOptimiz::StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking64},
            {"StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking96",
             AdvanceAlgOptimiz::StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking96},
            {"StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking128",
             AdvanceAlgOptimiz::StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking128},
            {"StdC_KSrc2Unroll_BoundsRefactor_LoopReorder_Blocking32",
             AdvanceAlgOptimiz::StdC_KSrc2Unroll_BoundsRefactor_LoopReorder_Blocking32},
            {"StdC_KSrc4Unroll_BoundsRefactor_LoopReorder_Blocking32",
             AdvanceAlgOptimiz::StdC_KSrc4Unroll_BoundsRefactor_LoopReorder_Blocking32},
#ifdef __AVX2__
            {"StdC_KSrc2Unroll_Vector_BoundsRefactor_LoopReorder_Blocking32",
             AdvanceAlgOptimiz::StdC_KSrc2Unroll_Vector_BoundsRefactor_LoopReorder_Blocking32},
            {"StdC_KSrc4Unroll_Vector_BoundsRefactor_LoopReorder_Blocking32",
             AdvanceAlgOptimiz::StdC_KSrc4Unroll_Vector_BoundsRefactor_LoopReorder_Blocking32},
#endif
        };

        // Run on the small - medium inputs; block divisor 4
        for (const auto& pair : imageQuiltingFunctions) {
            timing::run_timing_functional(pair.first, "./timing/input", "input0", "./timing/results",
                                          pair.second, 2, 4, 0);
        }
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