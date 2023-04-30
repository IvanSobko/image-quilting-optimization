#include "ImageQuilting.h"
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <random>
#include <cfloat>

// Synthesize a new texture
void ImageQuilting::Synthesis(){
    OverlapConstraintsWithMinCut();
}

// Write a block from the source data to the output data given their upper-left corners
void ImageQuilting::WriteBlock(const int dstY, const int dstX, const int srcY, const int srcX){
    for (int i = 0; i < mData->block_h; i++){
        for (int j = 0; j < mData->block_w; j++){
            for (int k = 0; k < CHANNEL_NUM; k++){
                mData->output_d[dstY + i][CHANNEL_NUM * (dstX + j) + k] = mData->data[srcY + i][
                        CHANNEL_NUM * (srcX + j) + k];
            }
        }
    }
}

// Same as WriteBlock but leaves the half of the dst overlapping region untouched
void ImageQuilting::WriteBlockOverlap(int dstY, int dstX, int srcY, int srcX)
{
    // Compute the height and width of the block to write
    int height = mData->block_h - overlapHeight / 2.0;
    int width = mData->block_w - overlapWidth / 2.0;
    // Clamp the height and width to the output image dimensions
    height = std::min(height, std::min(dstY + height, (int)mData->output_h) - dstY);
    width = std::min(width, std::min(dstX + width, (int)mData->output_w) - dstX);
    for (int i = 0; i < height; i++){
        for (int j = 0; j < width; j++){
            for (int k = 0; k < CHANNEL_NUM; k++){
                mData->output_d[dstY + i][CHANNEL_NUM * (dstX + j) + k] = mData->data[srcY + i][
                        CHANNEL_NUM * (srcX + j) + k];
            }
        }
    }
}

// Same as WriteBlockOverlap, but uses a minimum cut to write the new block
void ImageQuilting::WriteBlockOverlapWithMinCut(int overlapType, int dstY, int dstX, int srcY, int srcX)
{
    // Minimum cut paths
    int verticalPath[mData->block_h];
    int horizontalPath[mData->block_w];

    // Compute the correct overlap start pixel coordinates
    int overlapXStart = overlapType != horizontal ? dstX - overlapWidth : dstX;
    int overlapYStart = overlapType != vertical ? dstY - overlapHeight : dstY;

    // Vertical minimum cut
    if (overlapType == vertical /* TODO || overlapType == both */){

        // Compute the error surface
        double errorSurface[overlapWidth * mData->block_h];
        for (int i = 0; i < mData->block_h; i++){
            for (int j = 0; j < overlapWidth; j++){
                // Compute the per pixel error
                double error = 0;
                for (int k = 0; k < CHANNEL_NUM; k++){
                    double x0 = mData->output_d[overlapYStart+i][CHANNEL_NUM*(overlapXStart+j)+k];
                    double x1 = mData->data[srcY+i][CHANNEL_NUM*(srcX+j)+k];
                    double norm = x0 - x1;
                    error += norm * norm;
                }
                errorSurface[i*overlapWidth+j] = error;
            }
        }

        // Vertical minimum cut using dynamic programming
        double dpTable[overlapWidth * mData->block_h];

        // Fill up the first row with the error surface
        for (int j = 0; j < overlapWidth; j++){
            dpTable[j] = errorSurface[j];
        }

        // DP going from the first row to the last row
        for (int i = 1; i < mData->block_h; i++){
            for (int j = 0; j < overlapWidth; j++){
                // Get the value directly above
                double minError = dpTable[(i-1)*overlapWidth+j];
                // Get the value to the left
                if (j > 0) minError = std::min(minError, dpTable[(i-1)*overlapWidth+(j-1)]);
                // Get the value to the right
                if (j < overlapWidth-1) minError = std::min(minError, dpTable[(i-1)*overlapWidth+(j+1)]);
                dpTable[i*overlapWidth+j] = errorSurface[i*overlapWidth+j] + minError;
            }
        }

        // Find the minimum of the last row
        double minError = dpTable[(mData->block_h-1)*overlapWidth];
        verticalPath[mData->block_h-1] = 0;
        for (int j = 1; j < overlapWidth; j++){
            double error = dpTable[(mData->block_h-1)*overlapWidth+j];
            if (error < minError){
                minError = error;
                verticalPath[mData->block_h-1] = j;
            }
        }

        // Traverse the dpTable from the last row to the first row to construct the vertical path
        for (int i = mData->block_h - 2; i >= 0; i--){
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
            if (j < overlapWidth-1){
                double rightError = dpTable[i*overlapWidth+j+1];
                if (rightError < localError){
                    localError = rightError;
                    verticalPath[i] = j+1;
                }
            }
        }
    }

    // Horizontal minimum cut
    if (overlapType == horizontal /* TODO || overlapType == both */){

        // Compute the error surface
        double errorSurface[overlapHeight * mData->block_w];
        for (int i = 0; i < overlapHeight; i++){
            for (int j = 0; j < mData->block_w; j++){
                // Compute the per pixel error
                double error = 0;
                for (int k = 0; k < CHANNEL_NUM; k++){
                    double x0 = mData->output_d[overlapYStart+i][CHANNEL_NUM*(overlapXStart+j)+k];
                    double x1 = mData->data[srcY+i][CHANNEL_NUM*(srcX+j)+k];
                    double norm = x0 - x1;
                    error += norm * norm;
                }
                errorSurface[i*mData->block_w+j] = error;
            }
        }

        // Horizontal minimum cut using dynamic programming
        double dpTable[overlapHeight * mData->block_w];

        // Fill up the first column with the error surface
        for (int i = 0; i < overlapHeight; i++){
            dpTable[i*mData->block_w] = errorSurface[i*mData->block_w];
        }

        // DP going from the first column to the last column
        for (int j = 1; j < mData->block_w; j++){
            for (int i = 0; i < overlapHeight; i++){
                // Get the value directly to the left
                double minError = dpTable[i*mData->block_w+j-1];
                // Get the value to the left and up
                if (i > 0) minError = std::min(minError, dpTable[(i-1)*mData->block_w+(j-1)]);
                // Get the value to the left and below
                if (i < overlapHeight-1) minError = std::min(minError, dpTable[(i+1)*mData->block_w+(j-1)]);
                dpTable[i*mData->block_w+j] = errorSurface[i*mData->block_w+j] + minError;
            }
        }

        // Find the minimum of the last column
        double minError = dpTable[mData->block_w-1];
        horizontalPath[mData->block_w-1] = 0;
        for (int i = 1; i < overlapHeight; i++){
            double error = dpTable[(i+1)*mData->block_w-1];
            if (error < minError){
                minError = error;
                horizontalPath[mData->block_w-1] = i;
            }
        }

        // Traverse the dpTable from the last row to the first row to construct the horizontal path
        for (int j = mData->block_w - 2; j >= 0; j--){
            // Get the path from the right column
            int i = horizontalPath[j+1];
            // Get the value directly on the left
            double localError = dpTable[i*mData->block_w+j];
            horizontalPath[j] = i;
            // Get the value to the left and above
            if (i > 0){
                double leftError = dpTable[(i-1)*mData->block_w+j];
                if (leftError < localError){
                    localError = leftError;
                    horizontalPath[j] = i-1;
                }
            }
            // Get the value to the left and below
            if (i < overlapHeight-1){
                double rightError = dpTable[(i+1)*mData->block_w+j];
                if (rightError < localError){
                    localError = rightError;
                    horizontalPath[j] = i+1;
                }
            }
        }
    }

    // Clamp the height and width to the output image dimensions
    int blockMaxHeight = std::min((int)mData->block_h, std::min(dstY + (int)mData->block_h, (int)mData->output_h) - dstY);
    int blockMaxWidth = std::min((int)mData->block_w, std::min(dstX + (int)mData->block_w, (int)mData->output_w) - dstX);
    int blockMaxOverlapHeight = std::min((int)mData->block_h, std::min(overlapYStart + (int)mData->block_h, (int)mData->output_h) - overlapYStart);
    int blockMaxOverlapWidth = std::min((int)mData->block_w, std::min(overlapXStart + (int)mData->block_w, (int)mData->output_w) - overlapXStart);

    // Vertical overlap
    if (overlapType == vertical) {
        for (int i = 0; i < blockMaxHeight; i++) {
            for (int j = 0; j < blockMaxOverlapWidth; j++) {
                // Write the source block only if we are to the right of it
                if (j > verticalPath[i]) {
                    for (int k = 0; k < CHANNEL_NUM; k++) {
                        mData->output_d[dstY + i][CHANNEL_NUM * (overlapXStart + j) + k] =
                            mData->data[srcY + i][CHANNEL_NUM * (srcX + j) + k];
                    }
                }
            }
        }
    }

    // Horizontal overlap
    if (overlapType == horizontal) {
        for (int i = 0; i < blockMaxOverlapHeight; i++){
            for (int j = 0; j < blockMaxWidth; j++){
                // If the column is to the left of the horizontal path, write the source block
                if (i > horizontalPath[j]) {
                    for (int k = 0; k < CHANNEL_NUM; k++) {
                        mData->output_d[overlapYStart + i][CHANNEL_NUM * (dstX + j) + k] =
                            mData->data[srcY + i][CHANNEL_NUM * (srcX + j) + k];
                    }
                }
            }
        }
    }

    // TODO Corner overlap
    /*
    if (overlapType == both){
        for (int i = 0; i < blockMaxOverlapHeight; i++){
            for (int j = 0; j < blockMaxOverlapWidth; j++){
                if (j > verticalPath[i] && i > horizontalPath[j]) {
                    for (int k = 0; k < CHANNEL_NUM; k++) {
                        mData->output_d[overlapYStart + i][CHANNEL_NUM * (overlapXStart + j) + k] =
                            mData->data[srcY + i][CHANNEL_NUM * (srcX + j) + k];
                    }
                }
            }
        }
    }
    */
}

// Compute the overlap between the current block - block 0 of the output image
// and block 1 of the input image given their upper-left corners
// and the position of the overlap
double ImageQuilting::ComputeOverlap(
        const int overlapType,
        const int block0Y, const int block0X,
        const int block1Y, const int block1X)
{
    int overlapXStart = block0X - overlapWidth;
    int overlapYStart = block0Y - overlapHeight;

    // Clamp the overlapXEnd and overlapYEnd to the output image dimensions
    int overlapYEnd = std::min(overlapHeight, std::min(overlapYStart + overlapHeight, (int)mData->output_h) - overlapYStart);
    int overlapXEnd = std::min(overlapWidth, std::min(overlapXStart + overlapWidth, (int)mData->output_w) - overlapXStart);
    int blockWidthEnd = std::min((int)mData->block_w, std::min(block0X + (int)mData->block_w, (int)mData->output_w) - block0X);
    int blockHeightEnd = std::min((int)mData->block_h, std::min(block0Y + (int)mData->block_h, (int)mData->output_h) - block0Y);

    // Compute the l2 norm of the overlap between the two blocks
    // Compute the horizontal overlap
    double l2norm = 0;
    if (overlapType != vertical) {
        for (int i = 0; i < overlapYEnd; i++){
            for (int j = 0; j < blockWidthEnd; j++){
                for (int k = 0; k < CHANNEL_NUM; k++){
                    double x0 = mData->output_d[overlapYStart+i][CHANNEL_NUM*(block0X+j)+k];
                    double x1 = mData->data[block1Y+i][CHANNEL_NUM*(block1X+j)+k];
                    double norm = x0 - x1;
                    l2norm += norm*norm;
                }
            }
        }
    }

    // Compute the vertical overlap
    if (overlapType != horizontal) {
        for (int i = 0; i < blockHeightEnd; i++){
            for (int j = 0; j < overlapXEnd; j++){
                for (int k = 0; k < CHANNEL_NUM; k++){
                    double x0 = mData->output_d[block0Y+i][CHANNEL_NUM*(overlapXStart+j)+k];
                    double x1 = mData->data[block1Y+i][CHANNEL_NUM*(block1X+j)+k];
                    double norm = x0 - x1;
                    l2norm += norm*norm;
                }
            }
        }
    }

    // Compute the corner edge overlap
    if (overlapType == both) {
        for (int i = 0; i < overlapYEnd; i++) {
            for (int j = 0; j < overlapXEnd; j++) {
                for (int k = 0; k < CHANNEL_NUM; k++) {
                    double x0 = mData->output_d[overlapYStart + i][CHANNEL_NUM * (overlapXStart + j) + k];
                    double x1 = mData->data[block1Y + i][CHANNEL_NUM * (block1X + j) + k];
                    double norm = x0 - x1;
                    l2norm += norm * norm;
                }
            }
        }
    }

    return std::sqrt(l2norm);
}

// Place an edge overlap block with respect to the given block of the output image
void ImageQuilting::PlaceEdgeOverlapBlock(
        const int blockY, const int blockX, const int maxBlockX, const int maxBlockY, double errorTolerance)
{
    // Calculate the overlap start position and the offset from where to write the block to the output
    int overlapXStart = blockX,  overlapYStart = blockY;
    int drawOffsetX = 0, drawOffsetY = 0;

    // Calculate the overlap type
    int overlapType;
    if (blockY == 0){
        overlapType = vertical;
        overlapXStart = blockX - overlapWidth;
        drawOffsetX = overlapWidth / 2;
    }
    else if (blockX == 0){
        overlapType = horizontal;
        overlapYStart = blockY - overlapHeight;
        drawOffsetY = overlapHeight / 2;
    }
    else{
        overlapType = both;
        overlapYStart = blockY - overlapHeight;
        drawOffsetY = overlapHeight / 2;
        overlapXStart = blockX - overlapWidth;
        drawOffsetX = overlapWidth / 2;
    }

    // Compute the value of each block
    int numBlocks = maxBlockY * maxBlockX;
    BlockValue* blocks = (BlockValue*) malloc(sizeof(BlockValue) * numBlocks);
    for (int i = 0; i < maxBlockY; i++){
        for (int j = 0; j < maxBlockX; j++){
            int blockIndex = i * maxBlockX + j;
            blocks[blockIndex].y = i;
            blocks[blockIndex].x = j;
            blocks[blockIndex].value = ComputeOverlap(overlapType, blockY, blockX,
                                                      i, j);
        }
    }

    // Find the minimum block value
    double minVal = DBL_MAX;
    for (int i = 0; i < numBlocks; i++) {
        if (blocks[i].value < minVal) {
            minVal = blocks[i].value;
        }
    }

    // Choose a random block within the tolerance
    double upperBound = (1.0 + errorTolerance) * minVal;
    BlockValue* suitableBlocks = (BlockValue*) malloc(sizeof(BlockValue) * numBlocks);
    int numSuitableBlocks = 0;
    for (int i = 0; i < numBlocks; i++) {
        if (blocks[i].value < upperBound) {
            suitableBlocks[numSuitableBlocks] = blocks[i];
            numSuitableBlocks++;
        }
    }

    // Sample and place a block
    std::random_device randomDevice;
    std::mt19937 randomNumberGenerator(randomDevice());
    std::uniform_int_distribution<std::mt19937::result_type> randomBlock(0, numSuitableBlocks);
    int blockIndex = randomBlock(randomNumberGenerator);
    WriteBlockOverlap(overlapYStart + drawOffsetY, overlapXStart + drawOffsetX,
                      suitableBlocks[blockIndex].y + drawOffsetY, suitableBlocks[blockIndex].x + drawOffsetX);

    // Clean up
    free(blocks);
    free(suitableBlocks);
}

// Place an edge overlap block with respect to the given block of the output image
void ImageQuilting::PlaceEdgeOverlapBlockWithMinCut(
        const int blockY, const int blockX, const int maxBlockX, const int maxBlockY, double errorTolerance)
{
    // Calculate the overlap type
    int overlapType;
    if (blockY == 0)
        overlapType = vertical;
    else if (blockX == 0)
        overlapType = horizontal;
    else
        overlapType = both;

    // Compute the value of each block
    int numBlocks = maxBlockY * maxBlockX;
    BlockValue* blocks = (BlockValue*) malloc(sizeof(BlockValue) * numBlocks);
    for (int i = 0; i < maxBlockY; i++){
        for (int j = 0; j < maxBlockX; j++){
            int blockIndex = i * maxBlockX + j;
            blocks[blockIndex].y = i;
            blocks[blockIndex].x = j;
            blocks[blockIndex].value = ComputeOverlap(overlapType,blockY, blockX,
                                                      i, j);
        }
    }

    // Find the minimum block value
    double minVal = DBL_MAX;
    for (int i = 0; i < numBlocks; i++) {
        if (blocks[i].value < minVal) {
            minVal = blocks[i].value;
        }
    }

    // Choose a random block within the tolerance
    double upperBound = (1.0 + errorTolerance) * minVal;
    BlockValue* suitableBlocks = (BlockValue*) malloc(sizeof(BlockValue) * numBlocks);
    int numSuitableBlocks = 0;
    for (int i = 0; i < numBlocks; i++) {
        if (blocks[i].value < upperBound) {
            suitableBlocks[numSuitableBlocks] = blocks[i];
            numSuitableBlocks++;
        }
    }

    // Sample and place a block
    std::random_device randomDevice;
    std::mt19937 randomNumberGenerator(randomDevice());
    std::uniform_int_distribution<std::mt19937::result_type> randomBlock(0, numSuitableBlocks);
    int blockIndex = randomBlock(randomNumberGenerator);
    WriteBlockOverlapWithMinCut(
        overlapType,
        blockY, blockX,
        suitableBlocks[blockIndex].y, suitableBlocks[blockIndex].x);

    // Clean up
    free(blocks);
    free(suitableBlocks);
}

// Synthesize a new texture by randomly choosing blocks satisfying constraints and applying minimum cuts
void ImageQuilting::OverlapConstraintsWithMinCut(){

    mData->AllocateOutput();

    overlapHeight = mData->block_h / 6;
    overlapWidth = mData->block_w / 6;

    int hStep = mData->block_h - overlapHeight;
    int wStep = mData->block_w - overlapWidth;

    // Compute block parameters
    // The first block is full size; the others are of size step due to overlapping
    int numBlocksY = (mData->output_h - mData->block_h) / hStep + 2;
    int numBlocksX = (mData->output_w - mData->block_w) / wStep + 2;
    int maxBlockY = mData->height - mData->block_h - 1;
    int maxBlockX = mData->width - mData->block_w - 1;

    // Randomly generate the upper-left corners of blocks
    std::random_device randomDevice;
    std::mt19937 randomNumberGenerator(randomDevice());
    std::uniform_int_distribution<std::mt19937::result_type> randomY(0, maxBlockY);
    std::uniform_int_distribution<std::mt19937::result_type> randomX(0, maxBlockX);

    // Iterate over the block upper-left corners
    for (int blockY = 0; blockY < numBlocksY; blockY++){
        for (int blockX = 0; blockX < numBlocksX; blockX++){

            // Top-left corner of the current block
            int dstY = blockY == 0 ? 0 : mData->block_h + hStep * (blockY - 1);
            int dstX = blockX == 0 ? 0 : mData->block_w + wStep * (blockX - 1);

            // Make sure we are inside of the output image
           if (dstY > mData->output_h || dstX > mData->output_w) continue;

            // Randomly choose a block and place it
            if (blockY == 0 && blockX == 0){
                // Randomly choose the upper-left corner of a block
                int srcY = randomY(randomNumberGenerator);
                int srcX = randomX(randomNumberGenerator);

                // Write the randomly chosen block to the output
                WriteBlock(dstY, dstX, srcY, srcX);
            } else {
                PlaceEdgeOverlapBlockWithMinCut(dstY, dstX, maxBlockX, maxBlockY, 0.5);
            }
        }
    }
}

// Synthesize a new texture sample by randomly choosing blocks satisfying overlap constraints
void ImageQuilting::OverlapConstraints(){

    mData->AllocateOutput();

    overlapHeight = mData->block_h / 6;
    overlapWidth = mData->block_w / 6;

    int hStep = mData->block_h - overlapHeight;
    int wStep = mData->block_w - overlapWidth;

    // Compute block parameters
    // The first block is full size; the others are of size step due to overlapping
    int numBlocksY = (mData->output_h - mData->block_h) / hStep + 2;
    int numBlocksX = (mData->output_w - mData->block_w) / wStep + 2;
    int maxBlockY = mData->height - mData->block_h - 1;
    int maxBlockX = mData->width - mData->block_w - 1;

    // Randomly generate the upper-left corners of blocks
    std::random_device randomDevice;
    std::mt19937 randomNumberGenerator(randomDevice());
    std::uniform_int_distribution<std::mt19937::result_type> randomY(0, maxBlockY);
    std::uniform_int_distribution<std::mt19937::result_type> randomX(0, maxBlockX);

    // Iterate over the block upper-left corners
    for (int blockY = 0; blockY < numBlocksY; blockY++){
        for (int blockX = 0; blockX < numBlocksX; blockX++){

            // Top-left corner of the current block
            int dstY = blockY == 0 ? 0 : mData->block_h + hStep * (blockY - 1);
            int dstX = blockX == 0 ? 0 : mData->block_w + wStep * (blockX - 1);

            // Make sure we are inside of the output image
            if (dstY > mData->output_h || dstX > mData->output_w) continue;

            // Randomly choose a block and place it
            if (blockY == 0 && blockX == 0){
                // Randomly choose the upper-left corner of a block
                int srcY = randomY(randomNumberGenerator);
                int srcX = randomX(randomNumberGenerator);

                // Write the randomly chosen block to the output
                WriteBlock(dstY, dstX, srcY, srcX);
            }
            // Otherwise place a block according to the overlap constraints
            else {
                PlaceEdgeOverlapBlock(dstY, dstX, maxBlockX, maxBlockY, 0.3);
            }
        }
    }
}

// Synthesize a new texture sample by randomly choosing blocks
void ImageQuilting::RandomBlockPlacement(){

    mData->AllocateOutput();

    // Randomly generate the upper-left corners of blocks
    std::random_device randomDevice;
    std::mt19937 randomNumberGenerator(randomDevice());
    std::uniform_int_distribution<std::mt19937::result_type> randomY(0,mData->height-mData->block_h-1);
    std::uniform_int_distribution<std::mt19937::result_type> randomX(0,mData->width-mData->block_w-1);

    // Compute block parameters
    int numBlocksY = mData->output_h / mData->block_h;
    int numBlocksX = mData->output_w / mData->block_w;

    // Iterate over the block upper-left corners
    for (int blockY = 0; blockY < numBlocksY; blockY++){
        for (int blockX = 0; blockX < numBlocksX; blockX++){

            // Top-left corner of the current block
            int y = mData->block_h * blockY;
            int x = CHANNEL_NUM * mData->block_w * blockX;

            // Randomly choose the upper-left corner of a block
            int offsetY = randomY(randomNumberGenerator);
            int offsetX = randomX(randomNumberGenerator);

            // Iterate over the block upper-left corners
            for (int blockY = 0; blockY < numBlocksY; blockY++){
                for (int blockX = 0; blockX < numBlocksX; blockX++){

                    // Top-left corner of the current block
                    int dstY = mData->block_h * blockY;
                    int dstX = CHANNEL_NUM * mData->block_w * blockX;

                    // Randomly choose the upper-left corner of a block
                    int srcY = randomY(randomNumberGenerator);
                    int srcX = randomX(randomNumberGenerator);

                    // Write the randomly chosen block to the output
                    WriteBlock(dstY, dstX, srcY, srcX);
                }
            }
        }
    }
}
