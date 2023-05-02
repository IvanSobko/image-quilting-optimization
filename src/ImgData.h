#pragma once

#include <cstdint>
#include <cstdlib>

#define CHANNEL_NUM 4

struct ImgData {
    // Each char * stores one row of an image, for each pixel 4 value {rgba}
    unsigned char** data = NULL;
    unsigned char** output_d = NULL;
    uint32_t width = 0, height = 0;
    uint32_t output_w = 0, output_h = 0;
    uint32_t block_w = 0, block_h = 0;

    // Free memory for the input texture
    void FreeInput(){
        for (int y = 0; y < height; y++) {
            free(data[y]);
        }
        free(data);
        data = NULL;
    }

    // Allocate memory for the output texture and initialize it to black
    void AllocateOutput(){
        // Allocate memory for the output texture
        output_d = (unsigned char **) malloc(sizeof(unsigned char *) * output_h);
        for(int y = 0; y < output_h; y++){
            output_d[y] = (unsigned char *) malloc(output_w * CHANNEL_NUM);
        }

        // Initialize the output texture to black
        for (int i = 0; i < output_h; i++){
            for (int j = 0; j < output_w; j++){
                for (int k = 0; k < CHANNEL_NUM; k++){
                    output_d[i][CHANNEL_NUM * (j) + k] = 0;
                }
            }
        }
    }

    // Free memory for the output texture
    void FreeOutput(){
        for (int y = 0; y < output_h; y++) {
            free(output_d[y]);
        }
        free(output_d);
        output_d = NULL;
    }
};