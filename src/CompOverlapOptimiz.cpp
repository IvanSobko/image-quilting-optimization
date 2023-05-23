#include "CompOverlapOptimiz.h"
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <random>
#include <cfloat>

#include <immintrin.h>
// IMPORTANT: make sure it's the same for all functions that you compare
// Added this because we only need to increase the REP # for individual functions
#define IND_FUNC_REP 10000

// Image quilting function wrapper for testing and timing
void CompOverlapOptimiz::BasicOpt(ImgData* imgData, int seed) {
    CompOverlapOptimiz imageQuilting(imgData);
    imageQuilting.Synthesis(seed, opt_indices);
}

void CompOverlapOptimiz::AlgOpt(ImgData* imgData, int seed) {
    CompOverlapOptimiz imageQuilting(imgData);
    imageQuilting.Synthesis(seed, opt_algorithm);
}

void CompOverlapOptimiz::UnrollOpt(ImgData* imgData, int seed) {
    CompOverlapOptimiz imageQuilting(imgData);
    imageQuilting.Synthesis(seed, opt_unroll);
}

void CompOverlapOptimiz::VectorizeOpt(ImgData* imgData, int seed) {
    CompOverlapOptimiz imageQuilting(imgData);
    imageQuilting.Synthesis(seed, opt_vectorize);
}

void CompOverlapOptimiz::GetComponentParameters(ImgData* imgData, int& overlapType, int& dstY, int& dstX,
                                                int& srcY, int& srcX) {
    int maxBlockYSrc = imgData->height - imgData->block_h;
    int maxBlockXSrc = imgData->width - imgData->block_w;
    int maxBlockYDst = imgData->output_h - imgData->block_h;
    int maxBlockXDst = imgData->output_w - imgData->block_w;
    srcY = GetRandomInt(0, maxBlockYSrc-1);
    srcX = GetRandomInt(0, maxBlockXSrc-1);
    dstY = GetRandomInt(imgData->block_h, maxBlockYDst-1);
    dstX = GetRandomInt(imgData->block_w, maxBlockXDst-1);
    overlapType = both;
}

volatile void CompOverlapOptimiz::BaseComponent(ImgData* imgData, int seed) {
    CompOverlapOptimiz imageQuilting(imgData);
    int overlapType, dstY, dstX, srcY, srcX;
    volatile double dummy = 0;
    for (int k = 0; k < IND_FUNC_REP; k++) {
        GetComponentParameters(imgData, overlapType, dstY, dstX, srcY, srcX);
        dummy += imageQuilting.ComputeOverlapBase(overlapType, dstY, dstX, srcY, srcX);
    }
}

volatile void CompOverlapOptimiz::BasicOptComponent(ImgData* imgData, int seed) {
    CompOverlapOptimiz imageQuilting(imgData);
    int overlapType, dstY, dstX, srcY, srcX;
    volatile double dummy = 0;
    for (int k = 0; k < IND_FUNC_REP; k++) {
        GetComponentParameters(imgData, overlapType, dstY, dstX, srcY, srcX);
        dummy += imageQuilting.ComputeOverlapBasicOpt(overlapType, dstY, dstX, srcY, srcX);
    }
}

volatile void CompOverlapOptimiz::AlgoOptComponent(ImgData* imgData, int seed) {
    CompOverlapOptimiz imageQuilting(imgData);
    int overlapType, dstY, dstX, srcY, srcX;
    volatile double dummy = 0;
    for (int k = 0; k < IND_FUNC_REP; k++) {
        GetComponentParameters(imgData, overlapType, dstY, dstX, srcY, srcX);
        dummy += imageQuilting.ComputeOverlapAlgImpr(overlapType, dstY, dstX, srcY, srcX);
    }
}

volatile void CompOverlapOptimiz::UnrollOptComponent(ImgData* imgData, int seed) {
    CompOverlapOptimiz imageQuilting(imgData);
    int overlapType, dstY, dstX, srcY, srcX;
    volatile double dummy = 0;
    for (int k = 0; k < IND_FUNC_REP; k++) {
        GetComponentParameters(imgData, overlapType, dstY, dstX, srcY, srcX);
        dummy += imageQuilting.ComputeOverlapUnroll(overlapType, dstY, dstX, srcY, srcX);
    }
}

void CompOverlapOptimiz::VectorizeOptComponent(ImgData* imgData, int seed) {
    CompOverlapOptimiz imageQuilting(imgData);
    int overlapType, dstY, dstX, srcY, srcX;
    GetComponentParameters(imgData, overlapType, dstY, dstX, srcY, srcX);
    imageQuilting.ComputeOverlapVectorize(overlapType, dstY, dstX, srcY, srcX);
}


// Synthesize a new texture
void CompOverlapOptimiz::Synthesis() {
    flopCount = 0;
    SeedRandomNumberGenerator();
    OverlapConstraintsWithMinCut();
}

// Synthesize a new texture with the given seed
void CompOverlapOptimiz::Synthesis(int seed, int opt) {
    flopCount = 0;
    opt_type = opt;
    SeedRandomNumberGenerator(seed);
    OverlapConstraintsWithMinCut();
}

// Seed the random number generator with the system time
void CompOverlapOptimiz::SeedRandomNumberGenerator() {
    // https://stackoverflow.com/questions/1190870/i-need-to-generate-random-numbers-in-c
    srand(time(0));
}

// Seed the random number generator with a specified seed
void CompOverlapOptimiz::SeedRandomNumberGenerator(int seed) {
    srand(seed);
}

// Generate a random number in the range [min, max]
int CompOverlapOptimiz::GetRandomInt(int min, int max) {
    // https://stackoverflow.com/questions/1190870/i-need-to-generate-random-numbers-in-c
    return rand() % (max - min + 1) + min;
}

// Write a block from the source data to the output data given their upper-left corners
void CompOverlapOptimiz::WriteBlock(const int dstY, const int dstX, const int srcY, const int srcX){
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

// Same as WriteBlockOverlap, but uses a minimum cut to write the new block
void CompOverlapOptimiz::WriteBlockOverlapWithMinCut(
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

    // Declare our error surface and dynamic programming table
    double* errorSurface = (double*)malloc(sizeof(double) * numOverlapPixels);
    double* dpTable = (double*)malloc(sizeof(double) * numOverlapPixels);

    // Vertical minimum cut
    if (overlapType == vertical || overlapType == both) {
        // Compute the error surface
        for (int i = 0; i < overlapHeightLocal; i++) {
            for (int j = 0; j < overlapWidthLocal; j++) {
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
                if (j > 0) {
                    minError = std::min(minError, dpTable[(i - 1) * overlapWidthLocal + (j - 1)]);
                }
                // Get the value to the right
                if (j < overlapWidthLocal - 1) {
                    minError = std::min(minError, dpTable[(i - 1) * overlapWidthLocal + (j + 1)]);
                }
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
                flopCount++;
                if (leftError < localError){
                    localError = leftError;
                    verticalPath[i] = j-1;
                }
            }
            // Get the value to the right
            if (j < overlapWidthLocal-1){
                double rightError = dpTable[i*overlapWidthLocal+j+1];
                flopCount++;
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
                if (i > 0) {
                    minError = std::min(minError, dpTable[(i - 1) * overlapWidthLocal + (j - 1)]);
                }
                // Get the value to the left and below
                if (i < overlapHeightLocal - 1) {
                    minError = std::min(minError, dpTable[(i + 1) * overlapWidthLocal + (j - 1)]);
                }
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
                flopCount++;
                if (leftError < localError){
                    localError = leftError;
                    horizontalPath[j] = i-1;
                }
            }
            // Get the value to the left and below
            if (i < overlapHeightLocal-1){
                flopCount++;
                double rightError = dpTable[(i+1)*overlapWidthLocal+j];
                if (rightError < localError){
                    localError = rightError;
                    horizontalPath[j] = i+1;
                }
            }
        }
    }

    // Clean up
    free(errorSurface);
    free(dpTable);

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
double CompOverlapOptimiz::ComputeOverlapBasicOpt(const int overlapType, const int dstY, const int dstX,
                                                  const int srcY, const int srcX)
{
    // Compute the overlap region that we are working with
    int overlapXStart = overlapType != horizontal ? (dstX - overlapWidth) : dstX;
    int overlapYStart = overlapType != vertical ? (dstY - overlapHeight) : dstY;
    int verticalBlockYEnd = std::min(overlapYStart + (int)mData->block_h, (int)mData->output_h);
    int horizontalBlockXEnd = std::min(overlapXStart + (int)mData->block_w, (int)mData->output_w);
    int verticalBlockHeightLocal = verticalBlockYEnd - dstY;
    int horizontalBlockWidthLocal = horizontalBlockXEnd - dstX;

    // Compute the l2 norm of the overlap between the two blocks
    int l2norm = 0;

    // Compute the vertical overlap
    if (overlapType == vertical || overlapType == both) {
        int srcYOffset = overlapType == both ? overlapHeight : 0;
        int srcYStart = srcY+srcYOffset;
        int dstXStart = CHANNEL_NUM * overlapXStart;
        int srcXStart = CHANNEL_NUM * srcX;
        for (int i = 0; i < verticalBlockHeightLocal; i++){
            unsigned char* outputRow = mData->output_d[dstY+i] + dstXStart;
            unsigned char* srcRow = mData->data[srcYStart+i] + srcXStart;
            for (int j = 0; j < overlapWidth; j++){
                for (int k = 0; k < CHANNEL_NUM; k++){
                    int x0 = *(outputRow++);
                    int x1 = *(srcRow++);
                    int norm = x0 - x1;
                    l2norm += norm*norm;
                }
            }
        }
    }

    // Compute the horizontal overlap
    if (overlapType == horizontal  || overlapType == both) {
        int srcXOffset = overlapType == both ? overlapWidth : 0;
        int dstXStart = CHANNEL_NUM * dstX;
        int srcXStart = CHANNEL_NUM * (srcX+srcXOffset);
        for (int i = 0; i < overlapHeight; i++){
            unsigned char* outputRow = mData->output_d[overlapYStart+i] + dstXStart;
            unsigned char* srcRow = mData->data[srcY+i] + srcXStart;
            for (int j = 0; j < horizontalBlockWidthLocal; j++){
                for (int k = 0; k < CHANNEL_NUM; k++){
                    int x0 = *(outputRow++);
                    int x1 = *(srcRow++);
                    int norm = x0 - x1;
                    l2norm += norm*norm;
                }
            }
        }
    }

    // Compute the corner edge overlap
    if (overlapType == both) {
        int dstXStart = CHANNEL_NUM * overlapXStart;
        int srcXStart = CHANNEL_NUM * srcX;
        for (int i = 0; i < overlapHeight; i++) {
            unsigned char* outputRow = mData->output_d[overlapYStart+i] + dstXStart;
            unsigned char* srcRow = mData->data[srcY+i] + srcXStart;
            for (int j = 0; j < overlapWidth; j++) {
                for (int k = 0; k < CHANNEL_NUM; k++) {
                    int x0 = *(outputRow++);
                    int x1 = *(srcRow++);
                    int norm = x0 - x1;
                    l2norm += norm * norm;
                }
            }
        }
    }

    return std::sqrt(l2norm);
}

// Tried to improve the overlap calculations additionally to ComputeOverlapBasicOpt
double CompOverlapOptimiz::ComputeOverlapAlgImpr(int overlapType, int dstY, int dstX, int srcY, int srcX) {
    // Compute the overlap region that we are working with
    int overlapXStart = overlapType != horizontal ? (dstX - overlapWidth) : dstX;
    int overlapYStart = overlapType != vertical ? (dstY - overlapHeight) : dstY;
    int verticalBlockYEnd = std::min(overlapYStart + (int)mData->block_h, (int)mData->output_h);
    int horizontalBlockXEnd = std::min(overlapXStart + (int)mData->block_w, (int)mData->output_w);
    int verticalBlockHeightLocal = verticalBlockYEnd - dstY;
    int horizontalBlockWidthLocal = horizontalBlockXEnd - overlapXStart;

    // Compute the l2 norm of the overlap between the two blocks
    int l2norm = 0;

    // Compute the horizontal overlap (+corner if needed)
    if (overlapType != vertical) {
        int dstXStart = CHANNEL_NUM * overlapXStart;
        int srcXStart = CHANNEL_NUM * srcX;
        for (int i = 0; i < overlapHeight; i++){
            unsigned char* outputRow = mData->output_d[overlapYStart+i] + dstXStart;
            unsigned char* srcRow = mData->data[srcY+i] + srcXStart;
            for (int j = 0; j < horizontalBlockWidthLocal; j++){
                for (int k = 0; k < CHANNEL_NUM; k++){
                    int x0 = *(outputRow++);
                    int x1 = *(srcRow++);
                    int norm = x0 - x1;
                    l2norm += norm*norm;
                }
            }
        }
    }

    // Compute the vertical overlap
    if (overlapType != horizontal) {
        int srcYOffset = overlapType == both ? overlapHeight : 0;
        int srcYStart = srcY+srcYOffset;
        int dstXStart = CHANNEL_NUM * overlapXStart;
        int srcXStart = CHANNEL_NUM * srcX;
        for (int i = 0; i < verticalBlockHeightLocal; i++){
            unsigned char* outputRow = mData->output_d[dstY+i] + dstXStart;
            unsigned char* srcRow = mData->data[srcYStart+i] + srcXStart;
            for (int j = 0; j < overlapWidth; j++){
                for (int k = 0; k < CHANNEL_NUM; k++){
                    int x0 = *(outputRow++);
                    int x1 = *(srcRow++);
                    int norm = x0 - x1;
                    l2norm += norm*norm;
                }
            }
        }
    }

    return std::sqrt(l2norm);
}

double CompOverlapOptimiz::ComputeOverlapUnroll(int overlapType, int dstY, int dstX, int srcY, int srcX) {
    // Compute the overlap region that we are working with
    int overlapXStart = overlapType != horizontal ? (dstX - overlapWidth) : dstX;
    int overlapYStart = overlapType != vertical ? (dstY - overlapHeight) : dstY;
    int verticalBlockYEnd = std::min(overlapYStart + (int)mData->block_h, (int)mData->output_h);
    int horizontalBlockXEnd = std::min(overlapXStart + (int)mData->block_w, (int)mData->output_w);
    int verticalBlockHeightLocal = verticalBlockYEnd - dstY;
    int horizontalBlockWidthLocal = horizontalBlockXEnd - overlapXStart;

    // Compute the l2 norm of the overlap between the two blocks
    int l2norm = 0;

    // Compute the horizontal overlap (+corner if needed)
    if (overlapType != vertical) {
        int dstXStart = CHANNEL_NUM * overlapXStart;
        int srcXStart = CHANNEL_NUM * srcX;
        for (int i = 0; i < overlapHeight; i++){
            unsigned char* outputRow = mData->output_d[overlapYStart+i] + dstXStart;
            unsigned char* srcRow = mData->data[srcY+i] + srcXStart;
            for (int j = 0; j < horizontalBlockWidthLocal; j+=2){
                int rDst1 = *(outputRow++);
                int rSrc1 = *(srcRow++);
                int gDst1 = *(outputRow++);
                int gSrc1 = *(srcRow++);
                int bDst1 = *(outputRow++);
                int bSrc1 = *(srcRow++);
                int aDst1 = *(outputRow++);
                int aSrc1 = *(srcRow++);

                int rDst2 = *(outputRow++);
                int rSrc2 = *(srcRow++);
                int gDst2 = *(outputRow++);
                int gSrc2 = *(srcRow++);
                int bDst2 = *(outputRow++);
                int bSrc2 = *(srcRow++);
                int aDst2 = *(outputRow++);
                int aSrc2 = *(srcRow++);

                int rDiff1 = rDst1 - rSrc1;
                int gDiff1 = gDst1 - gSrc1;
                int rNorm1 = rDiff1*rDiff1;
                int gNorm1 = gDiff1*gDiff1;

                int bDiff1 = bDst1 - bSrc1;
                int aDiff1 = aDst1 - aSrc1;
                int bNorm1 = bDiff1*bDiff1;
                int aNorm1 = aDiff1*aDiff1;

                int rDiff2 = rDst2 - rSrc2;
                int gDiff2 = gDst2 - gSrc2;
                int rNorm2 = rDiff2*rDiff2;
                int gNorm2 = gDiff2*gDiff2;

                int bDiff2 = bDst2 - bSrc2;
                int aDiff2 = aDst2 - aSrc2;
                int bNorm2 = bDiff2*bDiff2;
                int aNorm2 = aDiff2*aDiff2;

                int sum1 = rNorm1 + gNorm1 + bNorm1 + aNorm1;
                int sum2 = rNorm2 + gNorm2 + bNorm2 + aNorm2;

                l2norm += sum1 + sum2;
            }
        }
    }

    // Compute the vertical overlap
    if (overlapType != horizontal) {
        int srcYOffset = overlapType == both ? overlapHeight : 0;
        int srcYStart = srcY+srcYOffset;
        int dstXStart = CHANNEL_NUM * overlapXStart;
        int srcXStart = CHANNEL_NUM * srcX;
        for (int i = 0; i < verticalBlockHeightLocal; i++){
            unsigned char* outputRow = mData->output_d[dstY+i] + dstXStart;
            unsigned char* srcRow = mData->data[srcYStart+i] + srcXStart;
            for (int j = 0; j < overlapWidth; j+=2){
                int rDst1 = *(outputRow++);
                int rSrc1 = *(srcRow++);
                int gDst1 = *(outputRow++);
                int gSrc1 = *(srcRow++);
                int bDst1 = *(outputRow++);
                int bSrc1 = *(srcRow++);
                int aDst1 = *(outputRow++);
                int aSrc1 = *(srcRow++);

                int rDst2 = *(outputRow++);
                int rSrc2 = *(srcRow++);
                int gDst2 = *(outputRow++);
                int gSrc2 = *(srcRow++);
                int bDst2 = *(outputRow++);
                int bSrc2 = *(srcRow++);
                int aDst2 = *(outputRow++);
                int aSrc2 = *(srcRow++);

                int rDiff1 = rDst1 - rSrc1;
                int gDiff1 = gDst1 - gSrc1;
                int rNorm1 = rDiff1*rDiff1;
                int gNorm1 = gDiff1*gDiff1;

                int bDiff1 = bDst1 - bSrc1;
                int aDiff1 = aDst1 - aSrc1;
                int bNorm1 = bDiff1*bDiff1;
                int aNorm1 = aDiff1*aDiff1;

                int rDiff2 = rDst2 - rSrc2;
                int gDiff2 = gDst2 - gSrc2;
                int rNorm2 = rDiff2*rDiff2;
                int gNorm2 = gDiff2*gDiff2;

                int bDiff2 = bDst2 - bSrc2;
                int aDiff2 = aDst2 - aSrc2;
                int bNorm2 = bDiff2*bDiff2;
                int aNorm2 = aDiff2*aDiff2;

                int sum1 = rNorm1 + gNorm1 + bNorm1 + aNorm1;
                int sum2 = rNorm2 + gNorm2 + bNorm2 + aNorm2;

                l2norm += sum1 + sum2;
            }
        }
    }

    return std::sqrt(l2norm);
}


double CompOverlapOptimiz::ComputeOverlapVectorize(int overlapType, int dstY, int dstX, int srcY, int srcX) {
    // Compute the overlap region that we are working with
    int overlapXStart = overlapType != horizontal ? (dstX - overlapWidth) : dstX;
    int overlapYStart = overlapType != vertical ? (dstY - overlapHeight) : dstY;
    int verticalBlockYEnd = std::min(overlapYStart + (int)mData->block_h, (int)mData->output_h);
    int horizontalBlockXEnd = std::min(overlapXStart + (int)mData->block_w, (int)mData->output_w);
    int verticalBlockHeightLocal = verticalBlockYEnd - dstY;
    int horizontalBlockWidthLocal = horizontalBlockXEnd - overlapXStart;
    // Compute the l2 norm of the overlap between the two blocks
    volatile int l2norm = 0;

    // Compute the horizontal overlap (+corner if needed)
    if (overlapType != vertical) {
        int dstXStart = CHANNEL_NUM * overlapXStart;
        int srcXStart = CHANNEL_NUM * srcX;
        for (int i = 0; i < overlapHeight; i++) {
            unsigned char* outputRow = mData->output_d[overlapYStart + i] + dstXStart;
            unsigned char* srcRow = mData->data[srcY + i] + srcXStart;
            //TODO: check for edge cases (horizontalBlockWidthLocal is not power of 4)
            for (int j = 0; j < horizontalBlockWidthLocal; j += 4) {
                // load 16 8-bit integers(chars)
                __m128i dst = _mm_loadu_si128((__m128i*)(outputRow + j * 4));
                __m128i src = _mm_loadu_si128((__m128i*)(srcRow + j * 4));

                // now each value is int16, we need it to avoid overflow problems
                __m256i dst_16bit = _mm256_cvtepu8_epi16(dst);
                __m256i src_16bit = _mm256_cvtepu8_epi16(src);
                __m256i diff = _mm256_sub_epi16(dst_16bit, src_16bit);

                // convert diff to int32, again: to avoid overflow
                __m256i diff_high = _mm256_cvtepi16_epi32(_mm256_castsi256_si128(diff));
                __m256i diff_low = _mm256_cvtepi16_epi32(_mm256_extractf128_si256(diff, 1));

                __m256i norm_high = _mm256_mullo_epi32(diff_high, diff_high);
                __m256i norm_low = _mm256_mullo_epi32(diff_low, diff_low);

                __m256i norm_sum = _mm256_add_epi32(norm_high, norm_low);

                // The hadd instruction is inefficient, and may be split into two instructions for faster decoding
                __m128i sum =
                    _mm_add_epi32(_mm256_extracti128_si256(norm_sum, 1), _mm256_castsi256_si128(norm_sum));
                sum = _mm_add_epi32(sum, _mm_unpackhi_epi64(sum, sum));
                sum = _mm_add_epi32(sum, _mm_shuffle_epi32(sum, 1));
                l2norm += _mm_cvtsi128_si32(sum);
            }
        }
    }

    // Compute the vertical overlap
    if (overlapType != horizontal) {
        int srcYOffset = overlapType == both ? overlapHeight : 0;
        int srcYStart = srcY + srcYOffset;
        int dstXStart = CHANNEL_NUM * overlapXStart;
        int srcXStart = CHANNEL_NUM * srcX;
        for (int i = 0; i < verticalBlockHeightLocal; i++) {
            unsigned char* outputRow = mData->output_d[dstY + i] + dstXStart;
            unsigned char* srcRow = mData->data[srcYStart + i] + srcXStart;
            for (int j = 0; j < overlapWidth; j += 4) {
                // load 16 8-bit integers(chars)
                __m128i dst = _mm_loadu_si128((__m128i*)(outputRow + j * 4));
                __m128i src = _mm_loadu_si128((__m128i*)(srcRow + j * 4));

                // now each value is int16, we need it to avoid overflow problems
                __m256i dst_16bit = _mm256_cvtepu8_epi16(dst);
                __m256i src_16bit = _mm256_cvtepu8_epi16(src);
                __m256i diff = _mm256_sub_epi16(dst_16bit, src_16bit);

                // convert diff to int32, again: to avoid overflow
                __m256i diff_high = _mm256_cvtepi16_epi32(_mm256_castsi256_si128(diff));
                __m256i diff_low = _mm256_cvtepi16_epi32(_mm256_extractf128_si256(diff, 1));

                __m256i norm_high = _mm256_mullo_epi32(diff_high, diff_high);
                __m256i norm_low = _mm256_mullo_epi32(diff_low, diff_low);

                __m256i norm_sum = _mm256_add_epi32(norm_high, norm_low);

                // The hadd instruction is inefficient, and may be split into two instructions for faster decoding
                __m128i sum =
                    _mm_add_epi32(_mm256_extracti128_si256(norm_sum, 1), _mm256_castsi256_si128(norm_sum));
                sum = _mm_add_epi32(sum, _mm_unpackhi_epi64(sum, sum));
                sum = _mm_add_epi32(sum, _mm_shuffle_epi32(sum, 1));
                l2norm += _mm_cvtsi128_si32(sum);
            }
        }
    }

    return std::sqrt(l2norm);
}

// Base implementation of ComputeOverlap
double CompOverlapOptimiz::ComputeOverlapBase(const int overlapType, const int dstY, const int dstX, const int srcY, const int srcX)
{
    // Compute the overlap region that we are working with
    int overlapXStart = overlapType != horizontal ? (dstX - overlapWidth) : dstX;
    int overlapYStart = overlapType != vertical ? (dstY - overlapHeight) : dstY;
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
void CompOverlapOptimiz::PlaceEdgeOverlapBlockWithMinCut(
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
        for (int j = 0; j < maxBlockX; j++) {
            int blockIndex = i * maxBlockX + j;
            blocks[blockIndex].y = i;
            blocks[blockIndex].x = j;
            if (opt_type == opt_indices) {
                blocks[blockIndex].value = ComputeOverlapBasicOpt(overlapType, blockY, blockX, i, j);
            } else if (opt_type == opt_algorithm) {
                blocks[blockIndex].value = ComputeOverlapAlgImpr(overlapType, blockY, blockX, i, j);
            } else if (opt_type == opt_unroll) {
                blocks[blockIndex].value = ComputeOverlapUnroll(overlapType, blockY, blockX, i, j);
            } else if (opt_type == opt_vectorize) {
                blocks[blockIndex].value = ComputeOverlapVectorize(overlapType, blockY, blockX, i, j);
            }
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
    int blockIndex = GetRandomInt(0, numSuitableBlocks-1);
    WriteBlockOverlapWithMinCut(
            overlapType,
            blockY, blockX,
            suitableBlocks[blockIndex].y, suitableBlocks[blockIndex].x);

    // Clean up
    free(blocks);
    free(suitableBlocks);

    // flops for ComputeOverlap loop
    // Note: approximating verticalBlockHeightLocal and verticalBlockWidthLocal as block_h and block_w
    if (overlapType == vertical) {
        flopCount += numBlocks * (3 * CHANNEL_NUM * overlapWidth * mData->block_h + 1);
    } else if (overlapType == horizontal) {
        flopCount += numBlocks * (3 * CHANNEL_NUM * overlapHeight * mData->block_w + 1);
    } else {
        flopCount += numBlocks * ((3 * CHANNEL_NUM * overlapWidth * mData->block_h)
                                  + (3 * CHANNEL_NUM * overlapHeight * mData->block_w)
                                  + (3 * CHANNEL_NUM * overlapHeight * overlapWidth) + 1);
    }
    // flops for intermediate calculations
    flopCount += 2 * numBlocks + 2;

    // flops for WriteBlockOverlapWithMinCut + some flops are computed in code
    // TODO: overlapHeightLocal and overlapWidthLocal are actually block_h and block_w. Fix flop count
    // Note: approximating overlapHeightLocal and overlapWidthLocal as overlapHeight and overlapWidth
    if (overlapType == vertical) {
        flopCount += 3 * CHANNEL_NUM * overlapWidth * overlapHeight +  3 * overlapWidth * (overlapHeight - 1) +
                     (overlapWidth - 1);
    } else if (overlapType == horizontal) {
        flopCount += 3 * CHANNEL_NUM * overlapWidth * overlapHeight +  3 * overlapHeight * (overlapWidth - 1) +
                     (overlapHeight - 1);
    } else {
        flopCount += (3 * CHANNEL_NUM * overlapWidth * overlapHeight +  3 * overlapWidth * (overlapHeight - 1) +
                      (overlapWidth - 1)) + (3 * CHANNEL_NUM * overlapWidth * overlapHeight +  3 *
                                                                                               overlapHeight * (overlapWidth - 1) + (overlapHeight - 1));
    }
}

// Synthesize a new texture by randomly choosing blocks satisfying constraints and applying minimum cuts
void CompOverlapOptimiz::OverlapConstraintsWithMinCut(){

    // Compute block parameters
    int hStep = mData->block_h - overlapHeight;
    int wStep = mData->block_w - overlapWidth;

    // The first block is full size; the others are of size step due to overlapping
    int numBlocksY = (mData->output_h - mData->block_h) / hStep + 2;
    int numBlocksX = (mData->output_w - mData->block_w) / wStep + 2;
    int maxBlockY = mData->height - mData->block_h;
    int maxBlockX = mData->width - mData->block_w;

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
                int srcY = GetRandomInt(0, maxBlockY-1);
                int srcX = GetRandomInt(0, maxBlockX-1);

                // Write the randomly chosen block to the output
                WriteBlock(dstY, dstX, srcY, srcX);
            } else {
                PlaceEdgeOverlapBlockWithMinCut(dstY, dstX, maxBlockX, maxBlockY, 0.1);
            }
        }
    }
}

int64_t CompOverlapOptimiz::getFlopCount() const {
    return flopCount;
}
