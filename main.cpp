#include "src/pngReader.h"
#include <cstdlib>
#include <string>

// each char * stores one row of an image, for each pixel 4 value {rgba}
unsigned char **img_data = NULL;
uint32_t width, height;

// input parameters
std::string input_file = "../gallery/input0.png";
std::string output_file = "../gallery/output0.png";
uint32_t output_w = 0, output_h = 0;
uint32_t block_w = 0, block_h = 0;

void parse_args(int argc, char *argv[]) {
    std::string delimiter = "=";
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        std::string par = arg.substr(0, arg.find(delimiter));
        std::string val = arg.substr(arg.find(delimiter), std::string::npos);
        if (arg == "--input") {
            input_file = "../gallery/" + val;
        } else if (arg == "--output") {
            output_file = "../gallery/" + val;
        } else if (arg == "--width") {
            output_w = std::stoi(val);
        } else if (arg == "--height") {
            output_h = std::stoi(val);
        } else if (arg == "--blockW") {
            block_w = std::stoi(val);
        } else if (arg == "--blockH") {
            block_h = std::stoi(val);
        }
    }
}

// if invalid or were not set before
void set_default() {
    if (output_w == 0) {
        output_w = width * 2;
    }
    if (output_h == 0) {
        output_h = height * 2;
    }
    // TODO: I'm not sure if these are the exact maximum values, but it seems logical that we should bound the block size
    uint32_t max_block_w = width / 3;
    uint32_t max_block_h = height / 3;
    if (block_w == 0 || block_w > max_block_w) {
        block_w = width / 4;
    }
    if (block_h == 0 || block_h > max_block_h) {
        block_h = height / 4;
    }
}

int main(int argc, char *argv[]) {
    parse_args(argc, argv);
    file::read_png_file(input_file.c_str(), img_data, width, height);
    set_default();
    file::write_png_file(output_file.c_str(), img_data, width, height);
    return 0;
}