#pragma once

#include <cstdint>
#include <cstdlib>
#include <limits>

#define CHANNEL_NUM 4

struct ImgData {
    // Each char * stores one row of an image, for each pixel 4 value {rgba}
    unsigned char** data = NULL;
    unsigned char** output_d = NULL;
    uint32_t width = 0, height = 0;
    uint32_t output_w = 0, output_h = 0;
    uint32_t block_w = 0, block_h = 0;

    // Allocate memory for the input image
    void AllocateInput() {
        data = (unsigned char**)malloc(sizeof(unsigned char*) * height);
        for (int y = 0; y < height; y++) {
            data[y] = (unsigned char*)malloc(width * CHANNEL_NUM);
        }
    }

    // Free memory for the input texture
    void FreeInput() {
        for (int y = 0; y < height; y++) {
            free(data[y]);
        }
        free(data);
        data = NULL;
    }

    // Generate a random number in the range [min, max]
    int GetRandomInt(int min, int max) {
        // https://stackoverflow.com/questions/1190870/i-need-to-generate-random-numbers-in-c
        return rand() % (max - min + 1) + min;
    }

    // Randomize the input image
    void RandomizeInput() {
        const unsigned char unsigned_char_max = std::numeric_limits<unsigned char>::max();
        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                for (int k = 0; k < CHANNEL_NUM - 1; k++)
                    data[i][CHANNEL_NUM * (j) + k] = GetRandomInt(0, unsigned_char_max);
                data[i][CHANNEL_NUM * (j) + 3] = unsigned_char_max;
            }
        }
    }

    // Allocate memory for the output texture and initialize it to black
    void AllocateOutput() {
        // Allocate memory for the output texture
        output_d = (unsigned char**)malloc(sizeof(unsigned char*) * output_h);
        for (int y = 0; y < output_h; y++) {
            output_d[y] = (unsigned char*)malloc(output_w * CHANNEL_NUM);
        }

        // Initialize the output texture to black
        for (int i = 0; i < output_h; i++) {
            for (int j = 0; j < output_w; j++) {
                for (int k = 0; k < CHANNEL_NUM; k++) {
                    output_d[i][CHANNEL_NUM * (j) + k] = 0;
                }
            }
        }
    }

    // Free memory for the output texture
    void FreeOutput() {
        for (int y = 0; y < output_h; y++) {
            free(output_d[y]);
        }
        free(output_d);
        output_d = NULL;
    }
};