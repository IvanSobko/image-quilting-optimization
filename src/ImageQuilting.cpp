#include "ImageQuilting.h"
#include <cstdlib>
#include <random>
#include <algorithm>

ImgData ImageQuilting::Synthesis() {
    return OverlapConstraints();
}

// Write a block from the source data to the output data specified by the given upper-left corners
void ImageQuilting::WriteBlock(int dstY, int dstX, int srcY, int srcX) {
    for (int i = 0; i < mData.block_h; i++){
        for (int j = 0; j < mData.block_w; j++){
            for (int k = 0; k < CHANNEL_NUM; k++){
                // Write to the CHANNEL_NUM channels
                mData.output_d[dstY+i][dstX+CHANNEL_NUM*j+k] = mData.data[srcY+i][CHANNEL_NUM*(srcX+j)+k];
            }
        }
    }
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

            // Iterate over the block upper-left corners
            for (int blockY = 0; blockY < numBlocksY; blockY++) {
                for (int blockX = 0; blockX < numBlocksX; blockX++) {

                    // Top-left corner of the current block
                    int dstY = mData.block_h * blockY;
                    int dstX = CHANNEL_NUM * mData.block_w * blockX;

                    // Randomly choose the upper-left corner of a block
                    int srcY = randomY(randomNumberGenerator);
                    int srcX = randomX(randomNumberGenerator);

                    // Write the randomly chosen block to the output
                    WriteBlock(dstY, dstX, srcY, srcX);
                }
            }
        }
    }
    return mData;
}

// Compute the left edge overlap between two left and right blocks specified by their upper-left corners
double ImageQuilting::ComputeLeftEdgeOverlap(int leftY, int leftX, int rightY, int rightX){
    // Overlap edge width is 1/6 the size of the block
    int overlapWidth = mData.block_w / 6;

    // Compute the l2 norm of the overlap between the two blocks
    double l2norm = 0;

    return 0;
}

// Synthesize a new texture sample by randomly choosing blocks satisfying overlap constraints
ImgData ImageQuilting::OverlapConstraints()
{
    // Allocate memory for the synthesized texture
    mData.output_d = (unsigned char **) malloc(sizeof(unsigned char *) * mData.output_h);
    for(int y = 0; y < mData.output_h; y++) {
        mData.output_d[y] = (unsigned char *) malloc(mData.output_w * CHANNEL_NUM);
    }

    // Randomly generate the upper-left corners of blocks
    std::random_device randomDevice;
    std::mt19937 randomNumberGenerator(randomDevice());
    std::uniform_int_distribution<std::mt19937::result_type> randomY(0, mData.height-mData.block_h-1);
    std::uniform_int_distribution<std::mt19937::result_type> randomX(0, mData.width-mData.block_w-1);

    // Compute block parameters
    int numBlocksY = mData.output_h / mData.block_h;
    int numBlocksX = mData.output_w / mData.block_w;

    // Iterate over the block upper-left corners
    for (int blockY = 0; blockY < numBlocksY; blockY++) {
        for (int blockX = 0; blockX < numBlocksX; blockX++) {

            // Top-left corner of the current block
            int dstY = mData.block_h * blockY;
            int dstX = CHANNEL_NUM * mData.block_w * blockX;

            // Randomly choose a block and place it
            if (blockY == 0 && blockX == 0) {
                // Randomly choose the upper-left corner of a block
                int srcY = randomY(randomNumberGenerator);
                int srcX = randomX(randomNumberGenerator);

                // Write the randomly chosen block to the output
                WriteBlock(dstY, dstX, srcY, srcX);
            }
        }
    }

    // Return the output
    return mData;
}