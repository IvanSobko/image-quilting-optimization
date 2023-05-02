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
    // Compute the height and width of the block to write
    int height = mData->block_h;
    int width = mData->block_w;
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

// Same as WriteBlock but leaves the half of the dst overlapping region untouched
void ImageQuilting::WriteBlockOverlap(const int overlapType, const int dstY, const int dstX, const int srcY,
                                      const int srcX)
{
    // Compute the height and width of the block to write
    int height = mData->block_h;
    int width = mData->block_w;
    int overlapXStart = dstX;
    int overlapYStart = dstY;
    int srcOverlapXStart = srcX;
    int srcOverlapYStart = srcY;
    if (overlapType == vertical) {
        overlapXStart -= overlapWidth / 2;
        srcOverlapXStart += overlapWidth / 2;
        width -= overlapWidth / 2;
    } else if (overlapType == horizontal) {
        overlapYStart -= overlapHeight / 2;
        srcOverlapYStart += overlapHeight / 2;
        height -= overlapHeight / 2;
    } else {
        overlapXStart -= overlapWidth / 2;
        overlapYStart -= overlapHeight / 2;
        srcOverlapXStart += overlapWidth / 2;
        srcOverlapYStart += overlapHeight / 2;
        width -= overlapWidth / 2;
        height -= overlapHeight / 2;
    }
    // Clamp the height and width to the output image dimensions
    height = std::min(height, std::min(overlapYStart + height, (int)mData->output_h) - overlapYStart);
    width = std::min(width, std::min(overlapXStart + width, (int)mData->output_w) - overlapXStart);
    for (int i = 0; i < height; i++){
        for (int j = 0; j < width; j++){
            for (int k = 0; k < CHANNEL_NUM; k++){
                mData->output_d[overlapYStart + i][CHANNEL_NUM * (overlapXStart + j) + k] =
                        mData->data[srcOverlapYStart + i][CHANNEL_NUM * (srcOverlapXStart + j) + k];
            }
        }
    }
}

// Same as WriteBlockOverlap, but uses a minimum cut to write the new block
void ImageQuilting::WriteBlockOverlapWithMinCut(
    const int overlapType, const int dstY, const int dstX, const int srcY, const int srcX)
{
    // Compute the overlap region that we are working with
    int overlapYStart = overlapType != vertical ? dstY - overlapHeight : dstY;
    int overlapXStart = overlapType != horizontal ? dstX - overlapWidth : dstX;
    int overlapYEnd = std::min(overlapYStart + (int)mData->block_h, (int)mData->output_h);
    int overlapXEnd = std::min(overlapXStart + (int)mData->block_w, (int)mData->output_w);
    int overlapHeightLocal = overlapYEnd - overlapYStart;
    int overlapWidthLocal = overlapXEnd - overlapXStart;
    int numOverlapPixels = overlapHeightLocal * overlapWidthLocal;

    // Minimum cut paths
    int verticalPath[overlapHeightLocal];
    int horizontalPath[overlapWidthLocal];

    // Vertical minimum cut
    if (overlapType == vertical || overlapType == both){

        // Compute the error surface
        double errorSurface[numOverlapPixels];
        for (int i = 0; i < overlapHeightLocal; i++){
            for (int j = 0; j < overlapWidthLocal; j++){
                // Compute the per pixel error
                double error = 0;
                for (int k = 0; k < CHANNEL_NUM; k++){
                    double x0 = mData->output_d[overlapYStart+i][CHANNEL_NUM*(overlapXStart+j)+k];
                    double x1 = mData->data[srcY+i][CHANNEL_NUM*(srcX+j)+k];
                    double norm = x0 - x1;
                    error += norm * norm;
                }
                errorSurface[i*overlapWidthLocal+j] = error;
            }
        }

        // Vertical minimum cut using dynamic programming
        double dpTable[numOverlapPixels];

        // Fill up the first row with the error surface
        for (int j = 0; j < overlapWidthLocal; j++){
            dpTable[j] = errorSurface[j];
        }

        // DP going from the first row to the last row
        for (int i = 1; i < overlapHeightLocal; i++){
            for (int j = 0; j < overlapWidthLocal; j++){
                // Get the value directly above
                double minError = dpTable[(i-1)*overlapWidthLocal+j];
                // Get the value to the left
                if (j > 0) minError = std::min(minError, dpTable[(i-1)*overlapWidthLocal+(j-1)]);
                // Get the value to the right
                if (j < overlapWidthLocal-1) minError = std::min(minError, dpTable[(i-1)*overlapWidthLocal+(j+1)]);
                dpTable[i*overlapWidthLocal+j] = errorSurface[i*overlapWidthLocal+j] + minError;
            }
        }

        // Find the minimum of the last row
        double minError = dpTable[(overlapHeightLocal-1)*overlapWidthLocal];
        verticalPath[overlapHeightLocal-1] = 0;
        for (int j = 1; j < overlapWidthLocal; j++){
            double error = dpTable[(overlapHeightLocal-1)*overlapWidthLocal+j];
            if (error < minError){
                minError = error;
                verticalPath[overlapHeightLocal-1] = j;
            }
        }

        // Traverse the dpTable from the last row to the first row to construct the vertical path
        for (int i = overlapHeightLocal - 2; i >= 0; i--){
            // Get the path from the previous row
            int j = verticalPath[i+1];
            // Get the value directly above
            double localError = dpTable[i*overlapWidthLocal+j];
            verticalPath[i] = j;
            // Get the value to the left
            if (j > 0){
                double leftError = dpTable[i*overlapWidthLocal+j-1];
                if (leftError < localError){
                    localError = leftError;
                    verticalPath[i] = j-1;
                }
            }
            // Get the value to the right
            if (j < overlapWidthLocal-1){
                double rightError = dpTable[i*overlapWidthLocal+j+1];
                if (rightError < localError){
                    localError = rightError;
                    verticalPath[i] = j+1;
                }
            }
        }
    }

    // Horizontal minimum cut
    if (overlapType == horizontal || overlapType == both){

        // Compute the error surface
        double errorSurface[numOverlapPixels];
        for (int i = 0; i < overlapHeightLocal; i++){
            for (int j = 0; j < overlapWidthLocal; j++){
                // Compute the per pixel error
                double error = 0;
                for (int k = 0; k < CHANNEL_NUM; k++){
                    double x0 = mData->output_d[overlapYStart+i][CHANNEL_NUM*(overlapXStart+j)+k];
                    double x1 = mData->data[srcY+i][CHANNEL_NUM*(srcX+j)+k];
                    double norm = x0 - x1;
                    error += norm * norm;
                }
                errorSurface[i*overlapWidthLocal+j] = error;
            }
        }

        // Horizontal minimum cut using dynamic programming
        double dpTable[numOverlapPixels];

        // Fill up the first column with the error surface
        for (int i = 0; i < overlapHeightLocal; i++){
            dpTable[i*overlapWidthLocal] = errorSurface[i*overlapWidthLocal];
        }

        // DP going from the first column to the last column
        for (int j = 1; j < overlapWidthLocal; j++){
            for (int i = 0; i < overlapHeightLocal; i++){
                // Get the value directly to the left
                double minError = dpTable[i*overlapWidthLocal+j-1];
                // Get the value to the left and up
                if (i > 0) minError = std::min(minError, dpTable[(i-1)*overlapWidthLocal+(j-1)]);
                // Get the value to the left and below
                if (i < overlapHeightLocal-1) minError = std::min(minError, dpTable[(i+1)*overlapWidthLocal+(j-1)]);
                dpTable[i*overlapWidthLocal+j] = errorSurface[i*overlapWidthLocal+j] + minError;
            }
        }

        // Find the minimum of the last column
        double minError = dpTable[overlapWidthLocal-1];
        horizontalPath[overlapWidthLocal-1] = 0;
        for (int i = 1; i < overlapHeightLocal; i++){
            double error = dpTable[(i+1)*overlapWidthLocal-1];
            if (error < minError){
                minError = error;
                horizontalPath[overlapWidthLocal-1] = i;
            }
        }

        // Traverse the dpTable from the last row to the first row to construct the horizontal path
        for (int j = overlapWidthLocal - 2; j >= 0; j--){
            // Get the path from the right column
            int i = horizontalPath[j+1];
            // Get the value directly on the left
            double localError = dpTable[i*overlapWidthLocal+j];
            horizontalPath[j] = i;
            // Get the value to the left and above
            if (i > 0){
                double leftError = dpTable[(i-1)*overlapWidthLocal+j];
                if (leftError < localError){
                    localError = leftError;
                    horizontalPath[j] = i-1;
                }
            }
            // Get the value to the left and below
            if (i < overlapHeightLocal-1){
                double rightError = dpTable[(i+1)*overlapWidthLocal+j];
                if (rightError < localError){
                    localError = rightError;
                    horizontalPath[j] = i+1;
                }
            }
        }
    }

    // Write the source block
    for (int i = 0; i < overlapHeightLocal; i++) {
        for (int j = 0; j < overlapWidthLocal; j++) {
            // Use the vertical and horizontal path to determine if we write the given source pixel
            bool write = false;

            // Vertical overlap: write the source block if we are to the right of the path
            if (overlapType == vertical)
                write = j > verticalPath[i];
            // Horizontal overlap: write the source block if we are above the path
            else if (overlapType == horizontal)
                write = i > horizontalPath[j];
            // Corner overlap: write the source block if we are to the right and above the two paths
            else if (overlapType == both)
                write = j > verticalPath[i] && i > horizontalPath[j];

            // Write the given source pixel
            if (write){
                for (int k = 0; k < CHANNEL_NUM; k++) {
                    mData->output_d[overlapYStart + i][CHANNEL_NUM * (overlapXStart + j) + k] =
                        mData->data[srcY + i][CHANNEL_NUM * (srcX + j) + k];
                }
            }
        }
    }
}

// Compute the overlap between the current block - block 0 of the output image
// and block 1 of the input image given their upper-left corners
// and the position of the overlap
double ImageQuilting::ComputeOverlap(const int overlapType, const int dstY, const int dstX, const int srcY, const int srcX)
{
    // Compute the overlap region that we are working with
    int overlapXStart = dstX - overlapWidth;
    int overlapYStart = dstY - overlapHeight;
    int verticalBlockYEnd = std::min(overlapYStart + (int)mData->block_h, (int)mData->output_h);
    int horizontalBlockXEnd = std::min(overlapXStart + (int)mData->block_w, (int)mData->output_w);
    int verticalBlockHeightLocal = verticalBlockYEnd - dstY;
    int horizontalBlockWidthLocal = horizontalBlockXEnd - dstX;

    // Compute the l2 norm of the overlap between the two blocks
    double l2norm = 0;

    // Compute the vertical overlap
    if (overlapType == vertical || overlapType == both) {
        int srcYOffset = overlapType == both ? overlapHeight : 0;
        for (int i = 0; i < verticalBlockHeightLocal; i++){
            for (int j = 0; j < overlapWidth; j++){
                for (int k = 0; k < CHANNEL_NUM; k++){
                    double x0 = mData->output_d[dstY+i][CHANNEL_NUM*(overlapXStart +j)+k];
                    double x1 = mData->data[srcY+srcYOffset+i][CHANNEL_NUM*(srcX + j)+k];
                    double norm = x0 - x1;
                    l2norm += norm*norm;
                }
            }
        }
    }

    // Compute the horizontal overlap
    if (overlapType == horizontal  || overlapType == both) {
        int srcXOffset = overlapType == both ? overlapWidth : 0;
        for (int i = 0; i < overlapHeight; i++){
            for (int j = 0; j < horizontalBlockWidthLocal; j++){
                for (int k = 0; k < CHANNEL_NUM; k++){
                    double x0 = mData->output_d[overlapYStart +i][CHANNEL_NUM*(dstX+j)+k];
                    double x1 = mData->data[srcY + i][CHANNEL_NUM*(srcX+srcXOffset+j)+k];
                    double norm = x0 - x1;
                    l2norm += norm*norm;
                }
            }
        }
    }

    // Compute the corner edge overlap
    if (overlapType == both) {
        for (int i = 0; i < overlapHeight; i++) {
            for (int j = 0; j < overlapWidth; j++) {
                for (int k = 0; k < CHANNEL_NUM; k++) {
                    double x0 = mData->output_d[overlapYStart + i][CHANNEL_NUM * (overlapXStart + j) + k];
                    double x1 = mData->data[srcY + i][CHANNEL_NUM * (srcX + j) + k];
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

    // Calculate the overlap type
    int overlapType;
    if (blockY == 0){
        overlapType = vertical;
        overlapXStart = blockX - overlapWidth;
    }
    else if (blockX == 0){
        overlapType = horizontal;
        overlapYStart = blockY - overlapHeight;
    }
    else{
        overlapType = both;
        overlapYStart = blockY - overlapHeight;
        overlapXStart = blockX - overlapWidth;
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
        if (blocks[i].value <= upperBound) {
            suitableBlocks[numSuitableBlocks] = blocks[i];
            numSuitableBlocks++;
        }
    }

    // Sample and place a block
    std::random_device randomDevice;
    std::mt19937 randomNumberGenerator(randomDevice());
    std::uniform_int_distribution<std::mt19937::result_type> randomBlock(0, numSuitableBlocks-1);
    int blockIndex = randomBlock(randomNumberGenerator);
    WriteBlockOverlap(overlapType, blockY, blockX,
                      suitableBlocks[blockIndex].y, suitableBlocks[blockIndex].x);

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
        if (blocks[i].value <= upperBound) {
            suitableBlocks[numSuitableBlocks] = blocks[i];
            numSuitableBlocks++;
        }
    }

    // Sample and place a block
    std::random_device randomDevice;
    std::mt19937 randomNumberGenerator(randomDevice());
    std::uniform_int_distribution<std::mt19937::result_type> randomBlock(0, numSuitableBlocks-1);
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
    int maxBlockY = mData->height - mData->block_h;
    int maxBlockX = mData->width - mData->block_w;

    // Randomly generate the upper-left corners of blocks
    std::random_device randomDevice;
    std::mt19937 randomNumberGenerator(randomDevice());
    std::uniform_int_distribution<std::mt19937::result_type> randomY(0, maxBlockY-1);
    std::uniform_int_distribution<std::mt19937::result_type> randomX(0, maxBlockX-1);

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
                PlaceEdgeOverlapBlockWithMinCut(dstY, dstX, maxBlockX, maxBlockY, 0.1);
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
    int maxBlockY = mData->height - mData->block_h;
    int maxBlockX = mData->width - mData->block_w;

    // Randomly generate the upper-left corners of blocks
    std::random_device randomDevice;
    std::mt19937 randomNumberGenerator(randomDevice());
    std::uniform_int_distribution<std::mt19937::result_type> randomY(0, maxBlockY-1);
    std::uniform_int_distribution<std::mt19937::result_type> randomX(0, maxBlockX-1);

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
                PlaceEdgeOverlapBlock(dstY, dstX, maxBlockX, maxBlockY, 0);
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
    int numBlocksY = mData->output_h / mData->block_h + 1;
    int numBlocksX = mData->output_w / mData->block_w + 1;

    // Iterate over the block upper-left corners
    for (int blockY = 0; blockY < numBlocksY; blockY++){
        for (int blockX = 0; blockX < numBlocksX; blockX++){

            // Top-left corner of the current block
            int y = mData->block_h * blockY;
            int x = mData->block_w * blockX;

            // Randomly choose the upper-left corner of a block
            int offsetY = randomY(randomNumberGenerator);
            int offsetX = randomX(randomNumberGenerator);

            // Iterate over the block upper-left corners
            for (int blockY = 0; blockY < numBlocksY; blockY++){
                for (int blockX = 0; blockX < numBlocksX; blockX++){

                    // Top-left corner of the current block
                    int dstY = mData->block_h * blockY;
                    int dstX = mData->block_w * blockX;

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