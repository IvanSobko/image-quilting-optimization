#include <cstdlib>
#include <string>
#include <time.h>
#include "ImageQuilting.h"
#include "PngReader.h"

#include "benchmarking/timing.h"

// input parameters
std::string input_file = "./gallery/input0.png";
std::string output_file = "./gallery/output0.png";
ImgData img_data;

bool benchmark = false;

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
    // TODO: I'm not sure if these are the exact maximum values, but it seems logical that we should bound the block size
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
    parse_args(argc, argv);
    file::read_png_file(input_file.c_str(), img_data);
    set_default();

    // modifies img_data inside and creates output image
    ImageQuilting quilting(img_data);

    if (benchmark) {
        clock_t start = clock();
        double r = timing::rdtsc(&quilting);
        clock_t end = clock();
        double seconds = (double)(end - start) / CLOCKS_PER_SEC;
        printf("RDTSC instruction:\n %lf cycles. Measurement took %.2f sec.\n", r, seconds);
    } else {
        img_data = quilting.synthesis();
        file::write_png_file(output_file.c_str(), img_data);
    }

    return 0;
}