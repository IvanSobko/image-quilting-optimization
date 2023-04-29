#include "ImageQuilting.h"
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <random>
#include <cfloat>

// Synthesize a new texture
ImgData ImageQuilting::Synthesis(){
    return OverlapConstraintsWithMinCut();
}

// Write a block from the source data to the output data given their upper-left corners
void ImageQuilting::WriteBlock(const int dstY, const int dstX, const int srcY, const int srcX){
    for (int i = 0; i < mData.block_h; i++){
        for (int j = 0; j < mData.block_w; j++){
            for (int k = 0; k < CHANNEL_NUM; k++){
                // Write to the CHANNEL_NUM channels
                mData.output_d[dstY + i][CHANNEL_NUM * (dstX + j) + k] = mData.data[srcY + i][
                        CHANNEL_NUM * (srcX + j) + k];
            }
        }
    }
}

// Same as the regular one, but leaves the half of the dst overlapping region untouched
void ImageQuilting::WriteBlockOverlap(int dstY, int dstX, int srcY, int srcX)
{
    int height = mData.block_h - overlapHeight / 2.0;
    int width = mData.block_w - overlapWidth / 2.0;
    for (int i = 0; i < height; i++){
        for (int j = 0; j < width; j++){
            for (int k = 0; k < CHANNEL_NUM; k++){
                // Write to the CHANNEL_NUM channels
                mData.output_d[dstY + i][CHANNEL_NUM * (dstX + j) + k] = mData.data[srcY + i][
                        CHANNEL_NUM * (srcX + j) + k];
            }
        }
    }
}

// Same as the overlapping one, but applies the minimum cut
void ImageQuilting::WriteBlockOverlapWithMinCut(
    int overlapYStart, int overlapXStart, int dstY, int dstX, int srcY, int srcX)
{
    int overlapRegionHeight = dstY - overlapYStart;
    int overlapRegionWidth = dstX - overlapXStart;

    // Vertical min cut
    int verticalPath[mData.block_h];
    if (overlapRegionHeight == 0){

        // Compute the error surface
        double errorSurface[overlapWidth * mData.block_h];
        for (int i = 0; i < mData.block_h; i++){
            for (int j = 0; j < overlapRegionWidth; j++){
                // Compute the per pixel error
                double error = 0;
                for (int k = 0; k < CHANNEL_NUM; k++){
                    double x0 = mData.output_d[dstY+i][CHANNEL_NUM*(overlapXStart+j)+k];
                    double x1 = mData.data[srcY+i][CHANNEL_NUM*(srcX+j)+k];
                    double norm = x0 - x1;
                    error += norm * norm;
                }
                errorSurface[i*overlapWidth+j] = error;
            }
        }

        // Vertical minimum cut using dynamic programming
        double dpTable[overlapWidth * mData.block_h];
        // Fill up the first row with the error surface
        for (int j = 0; j < overlapRegionWidth; j++){
            dpTable[j] = errorSurface[j];
        }
        // DP going forward
        for (int i = 1; i < mData.block_h; i++){
            for (int j = 0; j < overlapRegionWidth; j++){
                // Get the value directly above
                double minError = dpTable[(i-1)*overlapWidth+j];
                // Get the value to the left
                if (j > 0) minError = std::min(minError, dpTable[(i-1)*overlapWidth+(j-1)]);
                // Get the value to the right
                if (j < overlapRegionWidth-1) minError = std::min(minError, dpTable[(i-1)*overlapWidth+(j+1)]);
                dpTable[i*overlapWidth+j] = minError;
            }
        }

        // Find the minimum vertical path
        // Find the minimum of the last row
        double minError = dpTable[(mData.block_h-1)*overlapWidth];
        verticalPath[mData.block_h-1] = 0;
        for (int j = 1; j < overlapRegionWidth; j++){
            double error = dpTable[(mData.block_h-1)*overlapWidth+j];
            if (error < minError){
                minError = error;
                verticalPath[mData.block_h-1] = j;
            }
        }

        // Traverse the dpTable upwards to construct the row
        for (int i = mData.block_h - 2; i >= 0; i--){
            // Get the path from the previous row
            int j = verticalPath[i+1];
            // Get the value directly above
            double localError = dpTable[i*overlapWidth+j];
            verticalPath[i] = j;
            // Get the value to the left
            if (j > 0){
                double leftError = dpTable[i*overlapWidth+j-1];
                if (leftError < localError){
                    localError = leftError;
                    verticalPath[i] = j-1;
                }
            }
            // Get the value to the right
            if (j < overlapRegionWidth-1){
                double rightError = dpTable[i*overlapWidth+j+1];
                if (rightError < localError){
                    localError = rightError;
                    verticalPath[i] = j+1;
                }
            }
        }
    }
    for (int i = 0; i < mData.block_h; i++){
        for (int j = 0; j < mData.block_w; j++){
            if (i < overlapRegionHeight) {
                // horizontal overlap (TODO)
            }
            // >= - general region; another part of the statement - overlap region
            // j > verticalPath[i] - the starting point for the source block
            if (j >= overlapRegionWidth || (j < overlapRegionWidth && j > verticalPath[i])) {
                for (int k = 0; k < CHANNEL_NUM; k++) {
                    // Write to the CHANNEL_NUM channels
//                    mData.output_d[dstY + i][CHANNEL_NUM * (overlapXStart + j) + k] = mData.data[srcY + i][
//                            CHANNEL_NUM * (srcX + j) + k];
                    int val = k == 3 ? 255 : 0;
                    mData.output_d[dstY + i][CHANNEL_NUM * (overlapXStart + j) + k] = val;
                }
            }
        }
    }
}

// Compute the overlap between the current block - block 0 of the output image
// and block 1 of the input image given their upper-left corners
// and the position of the overlap
double ImageQuilting::ComputeOverlap(
        const int overlapYStart, const int overlapXStart,
        const int block0Y, const int block0X,
        const int block1Y, const int block1X)
{
    // will work for all cases: vertical, horizontal, both - will go to 0 if no overlap
    int overlapRegionHeight = block0Y - overlapYStart;
    int overlapRegionWidth = block0X - overlapXStart;

    // Compute the l2 norm of the overlap between the two blocks
    // Compute the horizontal overlap
    double l2norm = 0;
    for (int i = 0; i < overlapRegionHeight; i++){
        for (int j = 0; j < mData.block_w; j++){
            for (int k = 0; k < CHANNEL_NUM; k++){
                double x0 = mData.output_d[overlapYStart+i][CHANNEL_NUM*(block0X+j)+k];
                double x1 = mData.data[block1Y+i][CHANNEL_NUM*(block1X+j)+k];
                double norm = x0 - x1;
                l2norm += norm*norm;
            }
        }
    }

    // Compute the vertical overlap
    for (int i = 0; i < mData.block_h; i++){
        for (int j = 0; j < overlapRegionWidth; j++){
            for (int k = 0; k < CHANNEL_NUM; k++){
                double x0 = mData.output_d[block0Y+i][CHANNEL_NUM*(overlapXStart+j)+k];
                double x1 = mData.data[block1Y+i][CHANNEL_NUM*(block1X+j)+k];
                double norm = x0 - x1;
                l2norm += norm*norm;
            }
        }
    }

    // Compute the corner edge overlap
    for (int i = 0; i < overlapRegionHeight; i++){
        for (int j = 0; j < overlapRegionWidth; j++){
            for (int k = 0; k < CHANNEL_NUM; k++){
                double x0 = mData.output_d[overlapYStart+i][CHANNEL_NUM*(overlapXStart+j)+k];
                double x1 = mData.data[block1Y+i][CHANNEL_NUM*(block1X+j)+k];
                double norm = x0 - x1;
                l2norm += norm*norm;
            }
        }
    }

    return std::sqrt(l2norm);
}

// Place an edge overlap block with respect to the given block of the output image
void ImageQuilting::PlaceEdgeOverlapBlock(
    const int blockY, const int blockX, const int maxBlockX, const int maxBlockY, double errorTolerance)
{
    // calculate an overlap start position and the offset from where to write the block to the output
    int overlapXStart, overlapYStart;
    int drawOffsetX = 0, drawOffsetY = 0;
    if (blockX == 0) {
        // no vertical overline
        overlapXStart = blockX;
    } else {
        overlapXStart = blockX - overlapWidth;
        drawOffsetX = overlapWidth / 2.0;
    }
    if (blockY == 0) {
        // no horizontal overline
        overlapYStart = blockY;
    } else {
        overlapYStart = blockY - overlapHeight;
        drawOffsetY = overlapHeight / 2.0;
    }

    // Compute the value of each block
    int numBlocks = maxBlockY * maxBlockX;
    BlockValue blocks[numBlocks];
    for (int i = 0; i < maxBlockY; i++){
        for (int j = 0; j < maxBlockX; j++){
            int blockIndex = i * maxBlockX + j;
            blocks[blockIndex].y = i;
            blocks[blockIndex].x = j;
            blocks[blockIndex].value = ComputeOverlap(overlapYStart, overlapXStart,
                                                      blockY, blockX,
                                                      i, j);
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
        auto block = blocks[i];
        if (block.value < (1.0 + errorTolerance) * minVal) {
            suitBlocks.push_back(block);
        }
    }

    // Sample and place a block
    std::random_device randomDevice;
    std::mt19937 randomNumberGenerator(randomDevice());
    std::uniform_int_distribution<std::mt19937::result_type> randomBlock(0, suitBlocks.size());
    int blockIndex = randomBlock(randomNumberGenerator);
    WriteBlockOverlap(overlapYStart + drawOffsetY, overlapXStart + drawOffsetX,
                      blocks[blockIndex].y + drawOffsetY, blocks[blockIndex].x + drawOffsetX);
}

// Place an edge overlap block with respect to the given block of the output image
void ImageQuilting::PlaceEdgeOverlapBlockWithMinCut(
    const int blockY, const int blockX, const int maxBlockX, const int maxBlockY, double errorTolerance)
{
    // calculate an overlap start position and the offset from where to write the block to the output
    int overlapXStart, overlapYStart;
    if (blockX == 0) {
        // no vertical overline
        overlapXStart = blockX;
    } else {
        overlapXStart = blockX - overlapWidth;
    }
    if (blockY == 0) {
        // no horizontal overline
        overlapYStart = blockY;
    } else {
        overlapYStart = blockY - overlapHeight;
    }

    // Compute the value of each block
    int numBlocks = maxBlockY * maxBlockX;
    BlockValue blocks[numBlocks];
    for (int i = 0; i < maxBlockY; i++){
        for (int j = 0; j < maxBlockX; j++){
            int blockIndex = i * maxBlockX + j;
            blocks[blockIndex].y = i;
            blocks[blockIndex].x = j;
            blocks[blockIndex].value = ComputeOverlap(overlapYStart, overlapXStart,
                                                      blockY, blockX,
                                                      i, j);
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
        auto block = blocks[i];
        if (block.value < (1.0 + errorTolerance) * minVal) {
            suitBlocks.push_back(block);
        }
    }

    // Sample and place a block
    std::random_device randomDevice;
    std::mt19937 randomNumberGenerator(randomDevice());
    std::uniform_int_distribution<std::mt19937::result_type> randomBlock(0, suitBlocks.size());
    int blockIndex = randomBlock(randomNumberGenerator);
    WriteBlockOverlapWithMinCut(
        overlapYStart, overlapXStart, blockY, blockX, blocks[blockIndex].y, blocks[blockIndex].x);
}


// Synthesize a new texture sample by randomly choosing blocks satisfying overlap constraints
ImgData ImageQuilting::OverlapConstraints(){
    // Allocate memory for the synthesized texture
    mData.output_d = (unsigned char **) malloc(sizeof(unsigned char *) * mData.output_h);
    for(int y = 0; y < mData.output_h; y++) {
        mData.output_d[y] = (unsigned char *) malloc(mData.output_w * CHANNEL_NUM);
    }

    overlapHeight = mData.block_h / 6;
    overlapWidth = mData.block_w / 6;

    int hStep = mData.block_h - overlapHeight;
    int wStep = mData.block_w - overlapWidth;

    // Compute block parameters
    // Two blocks of a full size from each side, all others are blocks with size equal to the Step
    int numBlocksY = (mData.output_h - 2*mData.block_h) / hStep + 2;
    int numBlocksX = (mData.output_w - 2*mData.block_w) / wStep + 2;
    int maxBlockY = mData.height - mData.block_h - 1;
    int maxBlockX = mData.width - mData.block_w - 1;

    // Randomly generate the upper-left corners of blocks
    std::random_device randomDevice;
    std::mt19937 randomNumberGenerator(randomDevice());
    std::uniform_int_distribution<std::mt19937::result_type> randomY(0, maxBlockY);
    std::uniform_int_distribution<std::mt19937::result_type> randomX(0, maxBlockX);

    // Iterate over the block upper-left corners
    for (int blockY = 0; blockY < numBlocksY; blockY++){
        for (int blockX = 0; blockX < numBlocksX; blockX++){

            // Top-left corner of the current block
            int dstY = blockY == 0 ? 0 : mData.block_h + hStep * (blockY - 1);
            int dstX = blockX == 0 ? 0 : mData.block_w + wStep * (blockX - 1);

            // Randomly choose a block and place it
            if (blockY == 0 && blockX == 0){
                // Randomly choose the upper-left corner of a block
                int srcY = randomY(randomNumberGenerator);
                int srcX = randomX(randomNumberGenerator);

                // Write the randomly chosen block to the output
                WriteBlock(dstY, dstX, srcY, srcX);
            } else {
                PlaceEdgeOverlapBlock(dstY, dstX, maxBlockX, maxBlockY, 0.3);
            }
        }
    }

    // Return the output
    return mData;
}


// Synthesize a new texture sample with minimum cut and overlap
ImgData ImageQuilting::OverlapConstraintsWithMinCut(){
    // Allocate memory for the synthesized texture
    mData.output_d = (unsigned char **) malloc(sizeof(unsigned char *) * mData.output_h);
    for(int y = 0; y < mData.output_h; y++) {
        mData.output_d[y] = (unsigned char *) malloc(mData.output_w * CHANNEL_NUM);
    }

    for (int i = 0; i < mData.output_h; i++){
        for (int j = 0; j < mData.output_w; j++){
            for (int k = 0; k < CHANNEL_NUM; k++){
                // Write to the CHANNEL_NUM channels
                int val = k == 3 ? 255: 0;
                mData.output_d[i][CHANNEL_NUM * (j) + k] = val;
            }
        }
    }

    overlapHeight = mData.block_h / 6;
    overlapWidth = mData.block_w / 6;

    int hStep = mData.block_h - overlapHeight;
    int wStep = mData.block_w - overlapWidth;

    // Compute block parameters
    // Two blocks of a full size from each side, all others are blocks with size equal to the Step
    int numBlocksY = (mData.output_h - 2*mData.block_h) / hStep + 2;
    int numBlocksX = (mData.output_w - 2*mData.block_w) / wStep + 2;
    int maxBlockY = mData.height - mData.block_h - 1;
    int maxBlockX = mData.width - mData.block_w - 1;

    // Randomly generate the upper-left corners of blocks
    std::random_device randomDevice;
    std::mt19937 randomNumberGenerator(randomDevice());
    std::uniform_int_distribution<std::mt19937::result_type> randomY(0, maxBlockY);
    std::uniform_int_distribution<std::mt19937::result_type> randomX(0, maxBlockX);

    // Iterate over the block upper-left corners
    for (int blockY = 0; blockY < numBlocksY; blockY++){
        for (int blockX = 0; blockX < numBlocksX; blockX++){

            // Top-left corner of the current block
            int dstY = blockY == 0 ? 0 : mData.block_h + hStep * (blockY - 1);
            int dstX = blockX == 0 ? 0 : mData.block_w + wStep * (blockX - 1);

            // Randomly choose a block and place it
            if (blockY == 0 && blockX == 0){
                // Randomly choose the upper-left corner of a block
                int srcY = randomY(randomNumberGenerator);
                int srcX = randomX(randomNumberGenerator);

                // Write the randomly chosen block to the output
                WriteBlock(dstY, dstX, srcY, srcX);
            } else if (blockY == 0) {
                PlaceEdgeOverlapBlockWithMinCut(dstY, dstX, maxBlockX, maxBlockY, 0.3);
            }
        }
    }

    // Return the output
    return mData;
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
