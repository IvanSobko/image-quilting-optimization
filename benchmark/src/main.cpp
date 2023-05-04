#include <agz-utils/file.h>
#include <agz-utils/image.h>
#include <agz-utils/string.h>

#include <cxxopts.hpp>

#include "../benchmark/tsc_x86.h"

#include "imageQuilting.h"

#define CYCLES_REQUIRED 1e8
#define RDTSC_LATENCY 26
#define REP 2


bool runTiming = true;

struct CmdArgs
{
    std::string inputFilename;
    std::string outputFilename;

    int outputWidth  = 0;
    int outputHeight = 0;

    int blockWidth  = 0;
    int blockHeight = 0;

    int overlapWidth  = 0;
    int overlapHeight = 0;

    bool useMSEBlockSelection = false;
    bool useMinEdgeCut        = false;

    float blockTolerance = 0;
};

std::optional<CmdArgs> parseCmdArgs(int argc, char *argv[])
{
    cxxopts::Options options("ImageQuilting");
    options.add_options()
        ("input",      "input image filename",      cxxopts::value<std::string>())
        ("output",     "output image filename",     cxxopts::value<std::string>())
        ("width",      "output image width",        cxxopts::value<int>())
        ("height",     "output image height",       cxxopts::value<int>())
        ("blockW",     "block width",               cxxopts::value<int>())
        ("blockH",     "block height",              cxxopts::value<int>())
        ("overlapW",   "overlapped area width",     cxxopts::value<int>())
        ("overlapH",   "overlapped area height",    cxxopts::value<int>())
        ("mseBlock",   "use MSE block selection",   cxxopts::value<bool>())
        ("minCutEdge", "use min cost cut edge",     cxxopts::value<bool>())
        ("tolerance",  "block selection tolerance", cxxopts::value<float>())
        ("help",       "help information");

    const auto args = options.parse(argc, argv);

    if(args.count("help"))
    {
        std::cout << options.help({ "" }) << std::endl;
        return {};
    }

    CmdArgs result;

    try
    {
        result.inputFilename  = args["input"] .as<std::string>();
        result.outputFilename = args["output"].as<std::string>();
        result.outputWidth    = args["width"] .as<int>();
        result.outputHeight   = args["height"].as<int>();
        result.blockWidth     = args["blockW"].as<int>();
        result.blockHeight    = args["blockH"].as<int>();

        if(args.count("overlapW"))
            result.overlapWidth = args["overlapW"].as<int>();
        else
            result.overlapWidth = std::max(1, result.blockWidth / 6);

        if(args.count("overlapH"))
            result.overlapHeight = args["overlapH"].as<int>();
        else
            result.overlapHeight = std::max(1, result.blockHeight / 6);

        if(args.count("mseBlock"))
            result.useMSEBlockSelection = args["mseBlock"].as<bool>();
        else
            result.useMSEBlockSelection = true;

        if(args.count("minCutEdge"))
            result.useMinEdgeCut = args["minCutEdge"].as<bool>();
        else
            result.useMinEdgeCut = true;

        if(args.count("tolerance"))
            result.blockTolerance = args["tolerance"].as<float>();
        else
            result.blockTolerance = 0.1f;
    }
    catch(...)
    {
        std::cout << options.help({ "" }) << std::endl;
        return {};
    }

    return result;
}

void run(int argc, char *argv[])
{
    using namespace agz;
    using namespace img;
    using namespace stdstr;

    const auto args = parseCmdArgs(argc, argv);
    if(!args)
        return;

    ImageQuilting imageQuilting;
    imageQuilting.setParams(
        args->blockWidth,   args->blockHeight,
        args->overlapWidth, args->overlapHeight,
        args->blockTolerance);
    imageQuilting.UseMSEBlockSelection(args->useMSEBlockSelection);
    imageQuilting.UseMinEdgeCut(args->useMinEdgeCut);
    
    auto src = Image<Float3>(load_rgb_from_file(
        args->inputFilename).map(
            [](const math::color3b &c)
    {
        return Float3(math::from_color3b<float>(c));
    }));

    auto dst = imageQuilting.generate(
        src, args->outputWidth, args->outputHeight);

    auto dstu8 = dst.map([](const Float3 &c)
    {
        return math::to_color3b<float>(c);
    }).get_data();

    file::create_directory_for_file(args->outputFilename);

    const std::string lowerOutput = to_lower(args->outputFilename);
    if(ends_with(lowerOutput, ".png"))
        save_rgb_to_png_file(args->outputFilename, dstu8);
    else if(ends_with(lowerOutput, ".jpg"))
        save_rgb_to_jpg_file(args->outputFilename, dstu8);
    else if(ends_with(lowerOutput, ".bmp"))
        save_rgb_to_bmp_file(args->outputFilename, dstu8);
    else
        throw std::runtime_error(
            "unsupported output format: " + args->outputFilename);

    std::cout << "Generation complete..." << std::endl;
}

double count_cycles(int argc, char *argv[])
{
    using namespace agz;
    using namespace img;
    using namespace stdstr;

    double cycles = 0;
    int num_runs = 1;
    double multiplier = 1;
    myInt64 start, end;

    const auto args = parseCmdArgs(argc, argv);
    if(!args)
        return 0.0;

    ImageQuilting imageQuilting;
    imageQuilting.setParams(
        args->blockWidth,   args->blockHeight,
        args->overlapWidth, args->overlapHeight,
        args->blockTolerance);
    imageQuilting.UseMSEBlockSelection(args->useMSEBlockSelection);
    imageQuilting.UseMinEdgeCut(args->useMinEdgeCut);
    
    auto src = Image<Float3>(load_rgb_from_file(
        args->inputFilename).map(
            [](const math::color3b &c)
    {
        return Float3(math::from_color3b<float>(c));
    }));

    // Warm-up phase: we determine a number of executions that allows
    // the code to be executed for at least CYCLES_REQUIRED cycles.
    // This helps excluding timing overhead when measuring small runtimes.
    printf("Doing warmup phase...");
    do {
        num_runs = num_runs * multiplier;
        start = start_tsc();
        for (size_t i = 0; i < num_runs; i++) {
            auto dst = imageQuilting.generate(
                src, args->outputWidth, args->outputHeight);
        }
        end = stop_tsc(start);

        cycles = (double)end;
        multiplier = (CYCLES_REQUIRED) / (cycles);

    } while (multiplier > 2);

    // Actual performance measurements repeated REP times.
    // We simply store all results and compute medians during post-processing.
    printf("actually measuring performance (%i times).\n", num_runs * REP);
    double total_cycles = 0;
    for (size_t j = 0; j < REP; j++) {
        start = start_tsc();
        for (size_t i = 0; i < num_runs; ++i) {
            auto dst = imageQuilting.generate(
                src, args->outputWidth, args->outputHeight);
        }
        end = stop_tsc(start) - RDTSC_LATENCY;

        cycles = (double)end / (double)num_runs;
        total_cycles += cycles;
    }

    total_cycles /= REP;
    return total_cycles;
}


int main(int argc, char *argv[])
{
    try
    {
        if(runTiming) {
            for (int i=24; i<=768; i*=2) 
            {
                char* argument[7];
                std::string tmp;

                tmp = "./ImageQuilting";
                argument[0] = tmp.data();
                tmp = "--input=input0_" + std::to_string(i) + "x" + std::to_string(i) + ".png";
                argument[1] = tmp.data();
                tmp = "--output=output.png";
                argument[2] = tmp.data();
                tmp = "--blockW=" + std::to_string(i>>1);
                argument[3] = tmp.data();
                tmp = "--blockH=" + std::to_string(i>>1);
                argument[4] = tmp.data();
                tmp = "--width=" + std::to_string(i<<1);
                argument[5] = tmp.data();
                tmp = "--height=" + std::to_string(i<<1);
                argument[6] = tmp.data();

                double cycles = count_cycles(7, argument);
                std::cout << cycles << std::endl;
            }
        }
        else
            run(argc, argv);
    }
    catch(const std::exception &err)
    {
        std::cerr << err.what() << std::endl;
        return -1;
    }
}
