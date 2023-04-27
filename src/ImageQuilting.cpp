#include "ImageQuilting.h"
#include <cstdlib>
#include <random>
#include <algorithm>

ImgData ImageQuilting::Synthesis() {
    return RandomBlockPlacement();
}

// Synthesize a new texture sample by randomly choosing blocks
ImgData ImageQuilting::RandomBlockPlacement()
{
    // Allocate memory for the synthesized texture
    mData.output_d = (unsigned char **) malloc(sizeof(unsigned char *) * mData.output_h);
    for(int y = 0; y < mData.output_h; y++) {
        mData.output_d[y] = (unsigned char *) malloc(mData.output_w * CHANNEL_NUM);
    }

    // Randomly generate the upper-left corners of blocks
    std::random_device randomDevice;
    std::mt19937 randomNumberGenerator(randomDevice());
    std::uniform_int_distribution<std::mt19937::result_type> randomY(0,mData.height-mData.block_h-1);
    std::uniform_int_distribution<std::mt19937::result_type> randomX(0,mData.width-mData.block_w-1);

    // Compute block parameters
    int numBlocksY = mData.output_h / mData.block_h;
    int numBlocksX = mData.output_w / mData.block_w;

    // Iterate over the block upper-left corners
    for (int blockY = 0; blockY < numBlocksY; blockY++) {
        for (int blockX = 0; blockX < numBlocksX; blockX++) {

            // Top-left corner of the current block
            int y = mData.block_h * blockY;
            int x = CHANNEL_NUM * mData.block_w * blockX;

            // Randomly choose the upper-left corner of a block
            int offsetY = randomY(randomNumberGenerator);
            int offsetX = randomX(randomNumberGenerator);

            // Write the randomly chosen block to the output
            for (int i = 0; i < mData.block_h; i++){
                for (int j = 0; j < mData.block_w; j++){
                    for (int k = 0; k < CHANNEL_NUM; k++){
                        // Write to the CHANNEL_NUM channels
                        mData.output_d[y+i][x+CHANNEL_NUM*j+k] = mData.data[offsetY+i][CHANNEL_NUM*(offsetX+j)+k];
                    }
                }
            }
        }
    }
    return mData;
}