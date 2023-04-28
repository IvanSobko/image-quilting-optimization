#include "ImageQuilting.h"
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <random>
#include <float.h>

// Synthesize a new texture
ImgData ImageQuilting::Synthesis(){
    return OverlapConstraints();
}

// Write a block from the source data to the output data given their upper-left corners
void ImageQuilting::WriteBlock(const int dstY, const int dstX, const int srcY, const int srcX){
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
ImgData ImageQuilting::RandomBlockPlacement(){
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
    for (int blockY = 0; blockY < numBlocksY; blockY++){
        for (int blockX = 0; blockX < numBlocksX; blockX++){

            // Top-left corner of the current block
            int y = mData.block_h * blockY;
            int x = CHANNEL_NUM * mData.block_w * blockX;

            // Randomly choose the upper-left corner of a block
            int offsetY = randomY(randomNumberGenerator);
            int offsetX = randomX(randomNumberGenerator);

            // Iterate over the block upper-left corners
            for (int blockY = 0; blockY < numBlocksY; blockY++){
                for (int blockX = 0; blockX < numBlocksX; blockX++){

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

// Compute the vertical edge overlap between block 0 of the output image and block 1 of the input image given their upper-left corners
double ImageQuilting::ComputeVerticalEdgeOverlap(
    const int block0Y, const int block0X, const int block1Y, const int block1X)
{
    // Overlap edge width is 1/6 the size of the block
    int overlapWidth = mData.block_w / 6;
    int block0XOverlapStart = block0X + CHANNEL_NUM * (mData.block_w - overlapWidth);

    // Compute the l2 norm of the overlap between the two blocks
    double l2norm = 0;
    for (int i = 0; i < mData.block_h; i++){
        for (int j = 0; j < overlapWidth; j++){
            for (int k = 0; k < CHANNEL_NUM; k++){
                double x0 = mData.output_d[block0Y+i][block0XOverlapStart+CHANNEL_NUM*j+k];
                double x1 = mData.data[block1Y+i][CHANNEL_NUM*(block1X+j)+k];
                double norm = std::abs(x0 - x1);
                l2norm += norm*norm;
            }
        }
    }
    return std::sqrt(l2norm);
}

// Compute the horizontal edge overlap between block 0 of the output image and block 1 of the input image given their upper-left corners
double ImageQuilting::ComputeHorizontalEdgeOverlap(
    const int block0Y, const int block0X, const int block1Y, const int block1X)
{
    // Overlap edge width is 1/6 the size of the block
    int overlapWidth = mData.block_h / 6;
    int block0YOverlapStart = block0Y + (mData.block_h - overlapWidth);

    // Compute the l2 norm of the overlap between the two blocks
    double l2norm = 0;
    for (int i = 0; i < overlapWidth; i++){
        for (int j = 0; j < mData.block_w; j++){
            for (int k = 0; k < CHANNEL_NUM; k++){
                double x0 = mData.output_d[block0YOverlapStart+i][block0X+CHANNEL_NUM*j+k];
                double x1 = mData.data[block1Y+i][CHANNEL_NUM*(block1X+j)+k];
                double norm = std::abs(x0 - x1);
                l2norm += norm*norm;
            }
        }
    }
    return std::sqrt(l2norm);
}

// BlockValue comparator
int ImageQuilting::BlockValueComparator(const void* blockValue1, const void* blockValue2) {
    return ((BlockValue*)blockValue1)->value < ((BlockValue*)blockValue2)->value;
}

// C style binary search with upper bound
size_t ImageQuilting::bs_upper_bound(
    const void* array, size_t n, const void* upper_bound,
    size_t width, int (*comparator)(const void*, const void*))
{
    // https://stackoverflow.com/questions/6443569/implementation-of-c-lower-bound
    // https://github.com/gcc-mirror/gcc/blob/master/libiberty/bsearch.c
    const char *base = (const char *) array;
    size_t l = 0;
    size_t h = n;
    while (l < h){
        size_t mid = l + (h - l) / 2;
        if (comparator(upper_bound, base + width * mid)) {
            l = mid + 1;
        } else {
            h = mid;
        }
    }
    return l;
}

// Place a vertical edge overlap block with respect to the given block of the output image
void ImageQuilting::PlaceVerticalEdgeOverlapBlock(
    const int blockY, const int blockX, const int maxBlockX, const int maxBlockY, double errorTolerance)
{
    // Compute the value of each block
    int numBlocks = maxBlockY * maxBlockX;
    BlockValue blocks[numBlocks];
    for (int i = 0; i < maxBlockY; i++){
        for (int j = 0; j < maxBlockX; j++){
            int blockIndex = i * maxBlockX + j;
            blocks[blockIndex].y = i;
            blocks[blockIndex].x = j;
            blocks[blockIndex].value = ComputeVerticalEdgeOverlap(blockY, blockX, i, j);
        }
    }
    // Sort the blocks by value
    double minVal = DBL_MAX;
    for (int i=0; i < numBlocks; i++) {
        const auto &block = blocks[i];
        if (block.value < minVal) {
            minVal = block.value;
        }
    }
    std::vector<BlockValue> suitBlocks;
    for (int i=0; i < numBlocks; i++) {
        const auto &block = blocks[i];
        if (block.value < (1.0 + errorTolerance) * minVal) {
            suitBlocks.push_back(block);
        }
    }
//    qsort(blocks, numBlocks, sizeof(BlockValue), BlockValueComparator);
//    BlockValue upperBound;
//    upperBound.value = errorTolerance * blocks[numBlocks-1].value;
//    size_t maxIndex = bs_upper_bound(blocks, numBlocks, &upperBound, sizeof(BlockValue), BlockValueComparator);

    // Sample and place a block
    std::random_device randomDevice;
    std::mt19937 randomNumberGenerator(randomDevice());
    std::uniform_int_distribution<std::mt19937::result_type> randomBlock(0, suitBlocks.size());
    int blockIndex = randomBlock(randomNumberGenerator);
    WriteBlock(blockY, blockX, blocks[blockIndex].y, blocks[blockIndex].x);
}

// Synthesize a new texture sample by randomly choosing blocks satisfying overlap constraints
ImgData ImageQuilting::OverlapConstraints(){
    // Allocate memory for the synthesized texture
    mData.output_d = (unsigned char **) malloc(sizeof(unsigned char *) * mData.output_h);
    for(int y = 0; y < mData.output_h; y++) {
        mData.output_d[y] = (unsigned char *) malloc(mData.output_w * CHANNEL_NUM);
    }

    // Compute block parameters
    int numBlocksY = mData.output_h / mData.block_h;
    int numBlocksX = mData.output_w / mData.block_w;
    int maxBlockY = mData.height - mData.block_h - 1;
    int maxBlockX = mData.width - mData.block_w - 1;

    // Randomly generate the upper-left corners of blocks
    std::random_device randomDevice;
    std::mt19937 randomNumberGenerator(randomDevice());
    std::uniform_int_distribution<std::mt19937::result_type> randomY(0, maxBlockY);
    std::uniform_int_distribution<std::mt19937::result_type> randomX(0, maxBlockX);

    // Iterate over the block upper-left corners
    for (int blockY = 0; blockY < 1; blockY++){
        for (int blockX = 0; blockX < numBlocksX; blockX++){

            // Top-left corner of the current block
            int dstY = mData.block_h * blockY;
            int dstX = CHANNEL_NUM * mData.block_w * blockX;

            // Randomly choose a block and place it
            if (blockY == 0 && blockX == 0){
                // Randomly choose the upper-left corner of a block
                int srcY = randomY(randomNumberGenerator);
                int srcX = randomX(randomNumberGenerator);

                // Write the randomly chosen block to the output
                WriteBlock(dstY, dstX, srcY, srcX);
            }
            // Otherwise place a vertical edge overlap block
            else
                PlaceVerticalEdgeOverlapBlock(dstY, dstX, maxBlockX, maxBlockY, 0.3);
        }
    }

    // Return the output
    return mData;
}