#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#include <utility>
#include <vector>
#include "image_quilting.h"
#include "png_reader.h"

#include "blocking.h"
#include "src/advance_alg_optimiz.h"
#include "src/benchmarking/timing.h"
#include "src/comp_overlap_optimiz.h"
#include "testing.h"

struct Args {
    std::string input_file = "./gallery/input0.png";
    std::string output_file = "./gallery/output0.png";
    uint32_t output_w = 0;
    uint32_t output_h = 0;
    uint32_t block_w = 0;
    uint32_t block_h = 0;
    bool runTiming = false;
    bool generate = false;
    bool test = false;
    bool testCorrectnessAndTiming = false;
    bool stabilize = false;
    bool advisor = false;
    bool timingFunctional = false;
};

Args parse_args(int argc, char* argv[]) {
    Args args;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        auto eq = arg.find('=');
        std::string name = (eq == std::string::npos) ? arg : arg.substr(0, eq);
        std::string value = (eq == std::string::npos) ? std::string() : arg.substr(eq + 1);

        if (name == "--input") {
            if (!value.empty())
                args.input_file = "./gallery/" + value;
        } else if (name == "--output") {
            if (!value.empty())
                args.output_file = "./gallery/" + value;
        } else if (name == "--width" && !value.empty()) {
            args.output_w = static_cast<uint32_t>(std::stoi(value));
        } else if (name == "--height" && !value.empty()) {
            args.output_h = static_cast<uint32_t>(std::stoi(value));
        } else if (name == "--blockW" && !value.empty()) {
            args.block_w = static_cast<uint32_t>(std::stoi(value));
        } else if (name == "--blockH" && !value.empty()) {
            args.block_h = static_cast<uint32_t>(std::stoi(value));
        } else if (name == "--timing") {
            args.runTiming = true;
        } else if (name == "--generate") {
            args.generate = true;
        } else if (name == "--test") {
            args.test = true;
        } else if (name == "--testCorrectnessAndTiming") {
            args.testCorrectnessAndTiming = true;
        } else if (name == "--stabilize") {
            args.stabilize = true;
        } else if (name == "--advisor") {
            args.advisor = true;
        } else if (name == "--timingFunctional") {
            args.timingFunctional = true;
        }
    }
    return args;
}

void test_all_optimizations() {
    Testing testing = Testing(2);
    testing.RegisterTestFunction(Testing::ImageQuiltingFunction, "default");
    testing.RegisterTestFunction(CompOverlapOptimiz::BasicOpt, "basic");
    testing.RegisterTestFunction(CompOverlapOptimiz::AlgOpt, "AlgImpr");
    testing.RegisterTestFunction(CompOverlapOptimiz::UnrollOpt, "Unroll");
    testing.RegisterTestFunction(CompOverlapOptimiz::UnrollMaxOpt, "UnrollMax");
#ifdef __AVX2__
    testing.RegisterTestFunction(CompOverlapOptimiz::VectorizeOpt, "Vectorize");
#endif
    testing.RegisterTestFunction(CompOverlapOptimiz::UnrollChnls, "unroll channels");

    testing.RegisterTestFunction(AdvanceAlgOptimiz::StdC_KUnroll_BoundsRefactor, "unroll_bounds");
    testing.RegisterTestFunction(AdvanceAlgOptimiz::StdC_KUnroll_BoundsRefactor_LoopReorder, "unroll_bounds_loop");
    testing.RegisterTestFunction(AdvanceAlgOptimiz::StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking32,
                                 "unroll_bounds_loop_block32");
    testing.RegisterTestFunction(AdvanceAlgOptimiz::StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking48,
                                 "unroll_bounds_loop_block48");
    testing.RegisterTestFunction(AdvanceAlgOptimiz::StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking64,
                                 "unroll_bounds_loop_block64");
    testing.RegisterTestFunction(AdvanceAlgOptimiz::StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking96,
                                 "unroll_bounds_loop_block96");
    testing.RegisterTestFunction(AdvanceAlgOptimiz::StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking128,
                                 "unroll_bounds_loop_block128");
    testing.RegisterTestFunction(AdvanceAlgOptimiz::StdC_KSrc2Unroll_BoundsRefactor_LoopReorder_Blocking32,
                                 "unroll2_bounds_loop_block32");
    testing.RegisterTestFunction(AdvanceAlgOptimiz::StdC_KSrc4Unroll_BoundsRefactor_LoopReorder_Blocking32,
                                 "unroll4_bounds_loop_block32");
#ifdef __AVX2__
    testing.RegisterTestFunction(AdvanceAlgOptimiz::StdC_KSrc2Unroll_Vector_BoundsRefactor_LoopReorder_Blocking32,
                                 "unroll2_bounds_loop_block32_simd");
    testing.RegisterTestFunction(AdvanceAlgOptimiz::StdC_KSrc4Unroll_Vector_BoundsRefactor_LoopReorder_Blocking32,
                                 "unroll4_bounds_loop_block32_simd");
    testing.RegisterTestFunction(AdvanceAlgOptimiz::StdC_KSrc8Unroll_Vector_BoundsRefactor_LoopReorder_Blocking32,
                                 "unroll8_bounds_loop_block32_simd");
#endif
    testing.TestCorrectness();
}

static void run_benchmarking() {
    // Run timings for 3 block divisors.
    timing::run_timing(2);
    timing::run_timing(4);
    timing::run_timing(8);
}

static void generate_outputs() {
    Testing testing(2);
    testing.GenerateOutputFiles();
}

static void run_base_tests() {
    test_all_optimizations();
}

static void run_tests_and_timing(const Args& args) {
    Testing testing(2);

    testing.RegisterTestFunction(Testing::ImageQuiltingFunction, "default");
    testing.RegisterTestFunction(CompOverlapOptimiz::BasicOpt, "compBasic");
    testing.RegisterTestFunction(CompOverlapOptimiz::AlgOpt, "compBasic+AlgImpr");
    testing.RegisterTestFunction(CompOverlapOptimiz::UnrollOpt, "compBasic+AlgImpr+Unroll");
    testing.RegisterTestFunction(AdvanceAlgOptimiz::StdC_KUnroll_BoundsRefactor, "Unroll+DividedFunctions");
    testing.RegisterTestFunction(CompOverlapOptimiz::UnrollMaxOpt, "compBasic+AlgImpr+UnrollTheoreticalMax");
#ifdef __AVX2__
    testing.RegisterTestFunction(CompOverlapOptimiz::VectorizeOpt, "compBasic+AlgImpr+Unroll+Vectorize");
#endif

    std::cout << std::endl;
    testing.TestCorrectnessAndTiming(args.stabilize);

    // Register component-level tests
    testing.RegisterComponentTestFunction(CompOverlapOptimiz::BaseComponent, CompOverlapOptimiz::BaseComponent, "default");

    testing.RegisterComponentTestFunction(CompOverlapOptimiz::BaseComponent, CompOverlapOptimiz::BasicOptComponent,
                                          "compBasic");

    testing.RegisterComponentTestFunction(CompOverlapOptimiz::BaseComponent, CompOverlapOptimiz::AlgoOptComponent,
                                          "compBasic+AlgImpr");

    testing.RegisterComponentTestFunction(CompOverlapOptimiz::BaseComponent, CompOverlapOptimiz::UnrollOptComponent,
                                          "compBasic+AlgImpr+Unroll");

    testing.RegisterComponentTestFunction(CompOverlapOptimiz::BaseComponent, CompOverlapOptimiz::UnrollMaxOptComponent,
                                          "compBasic+AlgImpr+UnrollTheoreticalMax");

#ifdef __AVX2__
    testing.RegisterComponentTestFunction(CompOverlapOptimiz::BaseComponent, CompOverlapOptimiz::VectorizeOptComponent,
                                          "compBasic+AlgImpr+Unroll+Vectorize");
#endif

    std::cout << std::endl;
    testing.TestComponentsTiming(args.stabilize);
}

static void run_advisor() {
    Testing testing(0);

    testing.RegisterTestingComponentAdvisor(CompOverlapOptimiz::BaseComponent, "default");
    testing.RegisterTestingComponentAdvisor(CompOverlapOptimiz::BasicOptComponent, "compBasic");
    testing.RegisterTestingComponentAdvisor(CompOverlapOptimiz::UnrollOptComponent, "compBasic+AlgImpr+Unroll");
    testing.RegisterTestingComponentAdvisor(CompOverlapOptimiz::UnrollMaxOptComponent,
                                            "compBasic+AlgImpr+UnrollTheoreticalMax");
#ifdef __AVX2__
    testing.RegisterTestingComponentAdvisor(CompOverlapOptimiz::VectorizeOptComponent, "compBasic+AlgImpr+Unroll+Vectorize");
#endif
    testing.ComponentsTimingAdvisor();
}

static void run_timing_functional() {
    // Define the set of image quilting functions to test.
    std::vector<std::pair<std::string, timing::ImageQuiltingFunction>> imageQuiltingFunctions = {
        {"default", Testing::ImageQuiltingFunction},
        {"StdC_Algorithm", CompOverlapOptimiz::AlgOpt},
        {"StdC_Algorithm_ChannelsUnroll", CompOverlapOptimiz::UnrollChnls},
#ifdef __AVX2__
        {"StdC_Algorithm_Vectorize", CompOverlapOptimiz::VectorizeOpt},
#endif
        {"StdC_KUnroll_BoundsRefactor", AdvanceAlgOptimiz::StdC_KUnroll_BoundsRefactor},

        // Different variants / blockings
        {"StdC_KUnroll_BoundsRefactor_LoopReorder", AdvanceAlgOptimiz::StdC_KUnroll_BoundsRefactor_LoopReorder},
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

    auto run_for_divisor = [&](int divisor) {
        for (const auto& pair : imageQuiltingFunctions) {
            timing::run_timing_functional(pair.first, "./timing/input", "input0", "./timing/results", pair.second, 2, divisor,
                                          0);
        }
    };

    // Run for common divisors
    run_for_divisor(2);
    run_for_divisor(4);
    run_for_divisor(8);
}

// Ensure sensible defaults for derived img_data fields.
void set_default(ImgData& img_data) {
    if (img_data.output_w == 0)
        img_data.output_w = img_data.width * 2;
    if (img_data.output_h == 0)
        img_data.output_h = img_data.height * 2;

    uint32_t max_block_w = img_data.width / 3;
    uint32_t max_block_h = img_data.height / 3;
    if (img_data.block_w == 0 || img_data.block_w > max_block_w)
        img_data.block_w = img_data.width / 4;
    if (img_data.block_h == 0 || img_data.block_h > max_block_h)
        img_data.block_h = img_data.height / 4;
}

static void run_default_mode(const Args& args) {
    std::cout << "Running the image quilting algorithm in default mode." << std::endl;
    std::cout << "Input file: " << args.input_file << std::endl;
    std::cout << "Output file: " << args.output_file << std::endl;

    ImgData img_data;

    // Read the input data
    file::read_png_file(args.input_file.c_str(), img_data);
    set_default(img_data);

    // Allocate the output data and run the image quilting algorithm
    img_data.AllocateOutput();
    ImageQuilting imageQuilting(&img_data);
    imageQuilting.Synthesis();

    // Write the output file and free the members of img_data
    file::write_png_file(args.output_file.c_str(), img_data);
    img_data.FreeOutput();
    img_data.FreeInput();
}

int main(int argc, char* argv[]) {
    Args args = parse_args(argc, argv);

    if (args.runTiming) {
        run_benchmarking();
    } else if (args.generate) {
        generate_outputs();
    } else if (args.test) {
        run_base_tests();
    } else if (args.testCorrectnessAndTiming) {
        run_tests_and_timing(args);
    } else if (args.advisor) {
        run_advisor();
    } else if (args.timingFunctional) {
        run_timing_functional();
    } else {
        run_default_mode(args);
    }

    return EXIT_SUCCESS;
}