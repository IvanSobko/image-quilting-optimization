#include "AdvanceAlgOptimiz.h"
#include <algorithm>
#include <cfloat>
#include <cstdlib>
#include <iostream>
#include <random>
#include <climits>

#include <immintrin.h>
// IMPORTANT: make sure it's the same for all functions that you compare
// Added this because we only need to increase the REP # for individual functions
#define IND_FUNC_REP 10000

void AdvanceAlgOptimiz::DividedFuncOpt(ImgData* imgData, int seed) {
    AdvanceAlgOptimiz imageQuilting(imgData);
    imageQuilting.Synthesis(seed);
}

// Synthesize a new texture
void AdvanceAlgOptimiz::Synthesis() {
    flopCount = 0;
    SeedRandomNumberGenerator();
    OverlapConstraintsWithMinCut();
}

// Synthesize a new texture with the given seed
void AdvanceAlgOptimiz::Synthesis(int seed) {
    flopCount = 0;
    SeedRandomNumberGenerator(seed);
    OverlapConstraintsWithMinCutBlocking();
}

// Seed the random number generator with the system time
void AdvanceAlgOptimiz::SeedRandomNumberGenerator() {
    // https://stackoverflow.com/questions/1190870/i-need-to-generate-random-numbers-in-c
    srand(time(0));
}

// Seed the random number generator with a specified seed
void AdvanceAlgOptimiz::SeedRandomNumberGenerator(int seed) {
    srand(seed);
}

// Generate a random number in the range [min, max]
int AdvanceAlgOptimiz::GetRandomInt(int min, int max) {
    // https://stackoverflow.com/questions/1190870/i-need-to-generate-random-numbers-in-c
    return rand() % (max - min + 1) + min;
}

// Write a block from the source data to the output data given their upper-left corners
void AdvanceAlgOptimiz::WriteBlock(const int dstY, const int dstX, const int srcY, const int srcX) {
    // Compute the height and width of the block to write
    int height = mData->block_h;
    int width = mData->block_w;
    // Clamp the height and width to the output image dimensions
    height = std::min(height, std::min(dstY + height, (int)mData->output_h) - dstY);
    width = std::min(width, std::min(dstX + width, (int)mData->output_w) - dstX);
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            for (int k = 0; k < CHANNEL_NUM; k++) {
                mData->output_d[dstY + i][CHANNEL_NUM * (dstX + j) + k] =
                        mData->data[srcY + i][CHANNEL_NUM * (srcX + j) + k];
            }
        }
    }
}

// Same as WriteBlockOverlap, but uses a minimum cut to write the new block
void AdvanceAlgOptimiz::WriteBlockOverlapWithMinCut(const int overlapType, const int dstY, const int dstX,
                                                     const int srcY, const int srcX) {
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
                for (int k = 0; k < CHANNEL_NUM; k++) {
                    double x0 = mData->output_d[overlapYStart + i][CHANNEL_NUM * (overlapXStart + j) + k];
                    double x1 = mData->data[srcY + i][CHANNEL_NUM * (srcX + j) + k];
                    double norm = x0 - x1;
                    error += norm * norm;
                }
                errorSurface[i * overlapWidthLocal + j] = error;
            }
        }

        // Vertical minimum cut using dynamic programming

        // Fill up the first row with the error surface
        for (int j = 0; j < overlapWidthLocal; j++) {
            dpTable[j] = errorSurface[j];
        }

        // DP going from the first row to the last row
        for (int i = 1; i < overlapHeightLocal; i++) {
            for (int j = 0; j < overlapWidthLocal; j++) {
                // Get the value directly above
                double minError = dpTable[(i - 1) * overlapWidthLocal + j];
                // Get the value to the left
                if (j > 0) {
                    minError = std::min(minError, dpTable[(i - 1) * overlapWidthLocal + (j - 1)]);
                }
                // Get the value to the right
                if (j < overlapWidthLocal - 1) {
                    minError = std::min(minError, dpTable[(i - 1) * overlapWidthLocal + (j + 1)]);
                }
                dpTable[i * overlapWidthLocal + j] = errorSurface[i * overlapWidthLocal + j] + minError;
            }
        }

        // Find the minimum of the last row
        double minError = dpTable[(overlapHeightLocal - 1) * overlapWidthLocal];
        verticalPath[overlapHeightLocal - 1] = 0;
        for (int j = 1; j < overlapWidthLocal; j++) {
            double error = dpTable[(overlapHeightLocal - 1) * overlapWidthLocal + j];
            if (error < minError) {
                minError = error;
                verticalPath[overlapHeightLocal - 1] = j;
            }
        }

        // Traverse the dpTable from the last row to the first row to construct the vertical path
        for (int i = overlapHeightLocal - 2; i >= 0; i--) {
            // Get the path from the previous row
            int j = verticalPath[i + 1];
            // Get the value directly above
            double localError = dpTable[i * overlapWidthLocal + j];
            verticalPath[i] = j;
            // Get the value to the left
            if (j > 0) {
                double leftError = dpTable[i * overlapWidthLocal + j - 1];
                flopCount++;
                if (leftError < localError) {
                    localError = leftError;
                    verticalPath[i] = j - 1;
                }
            }
            // Get the value to the right
            if (j < overlapWidthLocal - 1) {
                double rightError = dpTable[i * overlapWidthLocal + j + 1];
                flopCount++;
                if (rightError < localError) {
                    localError = rightError;
                    verticalPath[i] = j + 1;
                }
            }
        }
    }

    // Horizontal minimum cut
    if (overlapType == horizontal || overlapType == both) {

        // Compute the error surface
        for (int i = 0; i < overlapHeightLocal; i++) {
            for (int j = 0; j < overlapWidthLocal; j++) {
                // Compute the per pixel error
                double error = 0;
                for (int k = 0; k < CHANNEL_NUM; k++) {
                    double x0 = mData->output_d[overlapYStart + i][CHANNEL_NUM * (overlapXStart + j) + k];
                    double x1 = mData->data[srcY + i][CHANNEL_NUM * (srcX + j) + k];
                    double norm = x0 - x1;
                    error += norm * norm;
                }
                errorSurface[i * overlapWidthLocal + j] = error;
            }
        }

        // Horizontal minimum cut using dynamic programming

        // Fill up the first column with the error surface
        for (int i = 0; i < overlapHeightLocal; i++) {
            dpTable[i * overlapWidthLocal] = errorSurface[i * overlapWidthLocal];
        }

        // DP going from the first column to the last column
        for (int j = 1; j < overlapWidthLocal; j++) {
            for (int i = 0; i < overlapHeightLocal; i++) {
                // Get the value directly to the left
                double minError = dpTable[i * overlapWidthLocal + j - 1];
                // Get the value to the left and up
                if (i > 0) {
                    minError = std::min(minError, dpTable[(i - 1) * overlapWidthLocal + (j - 1)]);
                }
                // Get the value to the left and below
                if (i < overlapHeightLocal - 1) {
                    minError = std::min(minError, dpTable[(i + 1) * overlapWidthLocal + (j - 1)]);
                }
                dpTable[i * overlapWidthLocal + j] = errorSurface[i * overlapWidthLocal + j] + minError;
            }
        }

        // Find the minimum of the last column
        double minError = dpTable[overlapWidthLocal - 1];
        horizontalPath[overlapWidthLocal - 1] = 0;
        for (int i = 1; i < overlapHeightLocal; i++) {
            double error = dpTable[(i + 1) * overlapWidthLocal - 1];
            if (error < minError) {
                minError = error;
                horizontalPath[overlapWidthLocal - 1] = i;
            }
        }

        // Traverse the dpTable from the last row to the first row to construct the horizontal path
        for (int j = overlapWidthLocal - 2; j >= 0; j--) {
            // Get the path from the right column
            int i = horizontalPath[j + 1];
            // Get the value directly on the left
            double localError = dpTable[i * overlapWidthLocal + j];
            horizontalPath[j] = i;
            // Get the value to the left and above
            if (i > 0) {
                double leftError = dpTable[(i - 1) * overlapWidthLocal + j];
                flopCount++;
                if (leftError < localError) {
                    localError = leftError;
                    horizontalPath[j] = i - 1;
                }
            }
            // Get the value to the left and below
            if (i < overlapHeightLocal - 1) {
                flopCount++;
                double rightError = dpTable[(i + 1) * overlapWidthLocal + j];
                if (rightError < localError) {
                    localError = rightError;
                    horizontalPath[j] = i + 1;
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
            if (write) {
                for (int k = 0; k < CHANNEL_NUM; k++) {
                    mData->output_d[overlapYStart + i][CHANNEL_NUM * (overlapXStart + j) + k] =
                            mData->data[srcY + i][CHANNEL_NUM * (srcX + j) + k];
                }
            }
        }
    }
}


// Place an edge overlap block with respect to the given block of the output image
void AdvanceAlgOptimiz::PlaceEdgeOverlapBlockWithMinCut(int overlapType, const int dstY, const int dstX,
                                                         const int maxBlockX, const int maxBlockY,
                                                         double errorTolerance, int bWidth, int bHeight)
{
    int overlapXStart = (overlapType == horizontal) ? dstX : dstX - overlapWidth;
    int overlapYStart = (overlapType == vertical) ? dstY : dstY - overlapHeight;

    int dstXStart = CHANNEL_NUM * overlapXStart;

    // Compute the value of each block
    int numBlocks = maxBlockY * maxBlockX;
    BlockValue* blocks = (BlockValue*)malloc(sizeof(BlockValue) * numBlocks);
    int minVal = INT_MAX;
    for (int srcY = 0; srcY < maxBlockY; srcY++) {
        for (int srcX = 0; srcX < maxBlockX; srcX++) {
            int blockIndex = srcY * maxBlockX + srcX;
            blocks[blockIndex].y = srcY;
            blocks[blockIndex].x = srcX;

            int l2norm = 0;
            int srcXStart = CHANNEL_NUM * srcX;

            if (overlapType == both) {

                for (int i = 0; i < overlapHeight; i++) {
                    unsigned char* outputRow = mData->output_d[overlapYStart + i] + dstXStart;
                    unsigned char* srcRow = mData->data[srcY + i] + srcXStart;
                    for (int j = 0; j < bWidth; j++) {
                        int rDst = outputRow[j*4];
                        int gDst = outputRow[j*4+1];
                        int bDst = outputRow[j*4+2];
                        int aDst = outputRow[j*4+3];
                        int rSrc = srcRow[j*4];
                        int gSrc = srcRow[j*4+1];
                        int bSrc = srcRow[j*4+2];
                        int aSrc = srcRow[j*4+3];

                        int rDiff = rDst - rSrc;
                        int gDiff = gDst - gSrc;
                        int bDiff = bDst - bSrc;
                        int aDiff = aDst - aSrc;

                        int norm1 = rDiff * rDiff;
                        int norm2 = gDiff * gDiff;
                        int norm3 = bDiff * bDiff;
                        int norm4 = aDiff * aDiff;

                        int sum1 = norm1 + norm2;
                        int sum2 = norm3 + norm4;

                        l2norm += sum1 + sum2;
                    }
                }
                // Compute the vertical overlap
                int srcYStart = srcY + overlapHeight;
                for (int i = 0; i < bHeight; i++) {
                    unsigned char* outputRow = mData->output_d[dstY + i] + dstXStart;
                    unsigned char* srcRow = mData->data[srcYStart + i] + srcXStart;
                    for (int j = 0; j < overlapWidth; j++) {
                        int rDst = outputRow[j*4];
                        int gDst = outputRow[j*4+1];
                        int bDst = outputRow[j*4+2];
                        int aDst = outputRow[j*4+3];
                        int rSrc = srcRow[j*4];
                        int gSrc = srcRow[j*4+1];
                        int bSrc = srcRow[j*4+2];
                        int aSrc = srcRow[j*4+3];

                        int rDiff = rDst - rSrc;
                        int gDiff = gDst - gSrc;
                        int bDiff = bDst - bSrc;
                        int aDiff = aDst - aSrc;

                        int norm1 = rDiff * rDiff;
                        int norm2 = gDiff * gDiff;
                        int norm3 = bDiff * bDiff;
                        int norm4 = aDiff * aDiff;

                        int sum1 = norm1 + norm2;
                        int sum2 = norm3 + norm4;

                        l2norm += sum1 + sum2;
                    }
                }

            } else if (overlapType == vertical) {

                for (int i = 0; i < mData->block_h; i++) {
                    unsigned char* outputRow = mData->output_d[dstY + i] + dstXStart;
                    unsigned char* srcRow = mData->data[srcY + i] + srcXStart;
                    for (int j = 0; j < overlapWidth; j++) {
                        int rDst = outputRow[j*4];
                        int gDst = outputRow[j*4+1];
                        int bDst = outputRow[j*4+2];
                        int aDst = outputRow[j*4+3];
                        int rSrc = srcRow[j*4];
                        int gSrc = srcRow[j*4+1];
                        int bSrc = srcRow[j*4+2];
                        int aSrc = srcRow[j*4+3];

                        int rDiff = rDst - rSrc;
                        int gDiff = gDst - gSrc;
                        int bDiff = bDst - bSrc;
                        int aDiff = aDst - aSrc;

                        int norm1 = rDiff * rDiff;
                        int norm2 = gDiff * gDiff;
                        int norm3 = bDiff * bDiff;
                        int norm4 = aDiff * aDiff;

                        int sum1 = norm1 + norm2;
                        int sum2 = norm3 + norm4;

                        l2norm += sum1 + sum2;
                    }
                }

            } else if (overlapType == horizontal) {

                for (int i = 0; i < overlapHeight; i++) {
                    unsigned char* outputRow = mData->output_d[overlapYStart + i] + dstXStart;
                    unsigned char* srcRow = mData->data[srcY + i] + srcXStart;
                    for (int j = 0; j < mData->block_w; j++) {
                        int rDst = outputRow[j*4];
                        int gDst = outputRow[j*4+1];
                        int bDst = outputRow[j*4+2];
                        int aDst = outputRow[j*4+3];
                        int rSrc = srcRow[j*4];
                        int gSrc = srcRow[j*4+1];
                        int bSrc = srcRow[j*4+2];
                        int aSrc = srcRow[j*4+3];

                        int rDiff = rDst - rSrc;
                        int gDiff = gDst - gSrc;
                        int bDiff = bDst - bSrc;
                        int aDiff = aDst - aSrc;

                        int norm1 = rDiff * rDiff;
                        int norm2 = gDiff * gDiff;
                        int norm3 = bDiff * bDiff;
                        int norm4 = aDiff * aDiff;

                        int sum1 = norm1 + norm2;
                        int sum2 = norm3 + norm4;

                        l2norm += sum1 + sum2;
                    }
                }

            }
            blocks[blockIndex].value = l2norm;
            if (l2norm < minVal) {
                minVal = l2norm;
            }
        }
    }

    // Choose a random block within the tolerance
    double upperBound = (1.0 + errorTolerance) * std::sqrt(minVal);
    upperBound = upperBound * upperBound;
    BlockValue* suitableBlocks = (BlockValue*)malloc(sizeof(BlockValue) * numBlocks);
    int numSuitableBlocks = 0;
    for (int i = 0; i < numBlocks; i++) {
        if (blocks[i].value <= upperBound) {
            suitableBlocks[numSuitableBlocks] = blocks[i];
            numSuitableBlocks++;
        }
    }

    // Sample and place a block
    int blockIndex = GetRandomInt(0, numSuitableBlocks - 1);
    WriteBlockOverlapWithMinCut(overlapType, dstY, dstX, suitableBlocks[blockIndex].y,
                                suitableBlocks[blockIndex].x);

    // Clean up
    free(blocks);
    free(suitableBlocks);
}

void AdvanceAlgOptimiz::PlaceEdgeOverlapBlockWithMinCutBlocking(int overlapType, int dstY, int dstX, int maxBlockX,
                                                            int maxBlockY, double errorTolerance, int bWidth, int bHeight) 
{
    int overlapXStart = (overlapType == horizontal) ? dstX : dstX - overlapWidth;
    int overlapYStart = (overlapType == vertical) ? dstY : dstY - overlapHeight;

    int dstXStart = CHANNEL_NUM * overlapXStart;

    // Compute the value of each block
    int numBlocks = maxBlockY * maxBlockX;
    BlockValue* blocks = (BlockValue*)malloc(sizeof(BlockValue) * numBlocks);
    int minVal = INT_MAX;

    for (int srcY = 0; srcY < maxBlockY; srcY++) {
        for (int srcX = 0; srcX < maxBlockX; srcX++) {
            int blockIndex = srcY * maxBlockX + srcX;
            blocks[blockIndex].y = srcY;
            blocks[blockIndex].x = srcX;
            blocks[blockIndex].value = 0;
        }
    }

    if (overlapType == both) {
        for (int i = 0; i < overlapHeight; i++) {
            unsigned char* outputRow = mData->output_d[overlapYStart + i] + dstXStart;
            for (int j = 0; j < bWidth; j++) {
                int rDst = outputRow[j*4];
                int gDst = outputRow[j*4+1];
                int bDst = outputRow[j*4+2];
                int aDst = outputRow[j*4+3];

                for (int srcY = 0; srcY < maxBlockY; srcY++) {
                    unsigned char* srcRow = mData->data[srcY + i];
                    for (int srcX = 0; srcX < maxBlockX; srcX++) {
                        int rSrc = srcRow[j*4];
                        int gSrc = srcRow[j*4+1];
                        int bSrc = srcRow[j*4+2];
                        int aSrc = srcRow[j*4+3];

                        int rDiff = rDst - rSrc;
                        int gDiff = gDst - gSrc;
                        int bDiff = bDst - bSrc;
                        int aDiff = aDst - aSrc;

                        int norm1 = rDiff * rDiff;
                        int norm2 = gDiff * gDiff;
                        int norm3 = bDiff * bDiff;
                        int norm4 = aDiff * aDiff;

                        int sum1 = norm1 + norm2;
                        int sum2 = norm3 + norm4;

                        blocks[srcY * maxBlockX + srcX].value += sum1 + sum2;
                        srcRow+=4;
                    }
                }
            }
        }
        // Compute the vertical overlap
        for (int i = 0; i < bHeight; i++) {
            unsigned char* outputRow = mData->output_d[dstY + i] + dstXStart;
            for (int j = 0; j < overlapWidth; j++) {
                int rDst = outputRow[j*4];
                int gDst = outputRow[j*4+1];
                int bDst = outputRow[j*4+2];
                int aDst = outputRow[j*4+3];
                for (int srcY = 0; srcY < maxBlockY; srcY++) {
                    int srcYStart = srcY + overlapHeight;
                    unsigned char* srcRow = mData->data[srcYStart + i];
                    for (int srcX = 0; srcX < maxBlockX; srcX++) {
                        int rSrc = srcRow[j*4];
                        int gSrc = srcRow[j*4+1];
                        int bSrc = srcRow[j*4+2];
                        int aSrc = srcRow[j*4+3];

                        int rDiff = rDst - rSrc;
                        int gDiff = gDst - gSrc;
                        int bDiff = bDst - bSrc;
                        int aDiff = aDst - aSrc;

                        int norm1 = rDiff * rDiff;
                        int norm2 = gDiff * gDiff;
                        int norm3 = bDiff * bDiff;
                        int norm4 = aDiff * aDiff;

                        int sum1 = norm1 + norm2;
                        int sum2 = norm3 + norm4;

                        blocks[srcY * maxBlockX + srcX].value += sum1 + sum2;
                        srcRow+=4;
                    }
                }
            }
        }
    } else if (overlapType == vertical) {
        for (int i = 0; i < mData->block_h; i++) {
            unsigned char* outputRow = mData->output_d[dstY + i] + dstXStart;

            for (int j = 0; j < overlapWidth; j++) {
                int rDst = outputRow[j*4];
                int gDst = outputRow[j*4+1];
                int bDst = outputRow[j*4+2];
                int aDst = outputRow[j*4+3];
                for (int srcY = 0; srcY < maxBlockY; srcY++) {
                    unsigned char* srcRow = mData->data[srcY+ i];
                    for (int srcX = 0; srcX < maxBlockX; srcX++) {
                        int rSrc = srcRow[j*4];
                        int gSrc = srcRow[j*4+1];
                        int bSrc = srcRow[j*4+2];
                        int aSrc = srcRow[j*4+3];

                        int rDiff = rDst - rSrc;
                        int gDiff = gDst - gSrc;
                        int bDiff = bDst - bSrc;
                        int aDiff = aDst - aSrc;

                        int norm1 = rDiff * rDiff;
                        int norm2 = gDiff * gDiff;
                        int norm3 = bDiff * bDiff;
                        int norm4 = aDiff * aDiff;

                        int sum1 = norm1 + norm2;
                        int sum2 = norm3 + norm4;

                        blocks[srcY * maxBlockX + srcX].value += sum1 + sum2;
                        srcRow+=4;
                    }
                }
            }
        }
    } else if (overlapType == horizontal) {
        for (int i = 0; i < overlapHeight; i++) {
            unsigned char* outputRow = mData->output_d[overlapYStart + i] + dstXStart;
            for (int j = 0; j < mData->block_w; j++) {
                int rDst = outputRow[j*4];
                int gDst = outputRow[j*4+1];
                int bDst = outputRow[j*4+2];
                int aDst = outputRow[j*4+3];

                for (int srcY = 0; srcY < maxBlockY; srcY++) {
                    unsigned char* srcRow = mData->data[srcY+ i];
                    for (int srcX = 0; srcX < maxBlockX; srcX++) {
                        int rSrc = srcRow[j*4];
                        int gSrc = srcRow[j*4+1];
                        int bSrc = srcRow[j*4+2];
                        int aSrc = srcRow[j*4+3];

                        int rDiff = rDst - rSrc;
                        int gDiff = gDst - gSrc;
                        int bDiff = bDst - bSrc;
                        int aDiff = aDst - aSrc;

                        int norm1 = rDiff * rDiff;
                        int norm2 = gDiff * gDiff;
                        int norm3 = bDiff * bDiff;
                        int norm4 = aDiff * aDiff;

                        int sum1 = norm1 + norm2;
                        int sum2 = norm3 + norm4;

                        blocks[srcY * maxBlockX + srcX].value += sum1 + sum2;
                        srcRow+=4;
                    }
                }
            }
        }
    }
    for (int i = 0; i < numBlocks; i++) {
        if (blocks[i].value < minVal) {
            minVal = blocks[i].value;
        }
    }

    // Choose a random block within the tolerance
    double upperBound = (1.0 + errorTolerance) * std::sqrt(minVal);
    upperBound = upperBound * upperBound;
    BlockValue* suitableBlocks = (BlockValue*)malloc(sizeof(BlockValue) * numBlocks);
    int numSuitableBlocks = 0;
    for (int i = 0; i < numBlocks; i++) {
        if (blocks[i].value <= upperBound) {
            suitableBlocks[numSuitableBlocks] = blocks[i];
            numSuitableBlocks++;
        }
    }

    // Sample and place a block
    int blockIndex = GetRandomInt(0, numSuitableBlocks - 1);
    WriteBlockOverlapWithMinCut(overlapType, dstY, dstX, suitableBlocks[blockIndex].y,
                                suitableBlocks[blockIndex].x);

    // Clean up
    free(blocks);
    free(suitableBlocks);
}

void AdvanceAlgOptimiz::PlaceEdgeOverlapBlockWithMinCutBlocking1(int overlapType, int dstY, int dstX, int maxBlockX,
                                                                int maxBlockY, double errorTolerance, int bWidth, int bHeight)
{
    int overlapXStart = (overlapType == horizontal) ? dstX : dstX - overlapWidth;
    int overlapYStart = (overlapType == vertical) ? dstY : dstY - overlapHeight;

    int dstXStart = CHANNEL_NUM * overlapXStart;

    // Compute the value of each block
    int numBlocks = maxBlockY * maxBlockX;
    BlockValue* blocks = (BlockValue*)malloc(sizeof(BlockValue) * numBlocks);
    int minVal = INT_MAX;

    for (int srcY = 0; srcY < maxBlockY; srcY++) {
        for (int srcX = 0; srcX < maxBlockX; srcX++) {
            int blockIndex = srcY * maxBlockX + srcX;
            blocks[blockIndex].y = srcY;
            blocks[blockIndex].x = srcX;
            blocks[blockIndex].value = 0;
        }
    }

    // Blocking parameters
    const int blockSize = 32;

    if (overlapType == both) {
        for (int i = 0; i < overlapHeight; i++) {
            unsigned char* outputRow = mData->output_d[overlapYStart + i] + dstXStart;
            for (int blockYStart = 0; blockYStart < maxBlockY + blockSize; blockYStart += blockSize) {
                for (int blockXStart = 0; blockXStart < maxBlockX + blockSize; blockXStart += blockSize) {
                    for (int j = 0; j < bWidth; j++) {
                        int rDst = outputRow[j * 4];
                        int gDst = outputRow[j * 4 + 1];
                        int bDst = outputRow[j * 4 + 2];
                        int aDst = outputRow[j * 4 + 3];

                        for (int srcY = blockYStart; srcY < std::min(blockYStart + blockSize, maxBlockY);
                             srcY++) {
                            unsigned char* srcRow = mData->data[srcY + i] + CHANNEL_NUM * blockXStart;
                            for (int srcX = blockXStart; srcX < std::min(blockXStart + blockSize, maxBlockX); srcX++) {
                                int rSrc = srcRow[j * 4];
                                int gSrc = srcRow[j * 4 + 1];
                                int bSrc = srcRow[j * 4 + 2];
                                int aSrc = srcRow[j * 4 + 3];

                                int rDiff = rDst - rSrc;
                                int gDiff = gDst - gSrc;
                                int bDiff = bDst - bSrc;
                                int aDiff = aDst - aSrc;

                                int norm1 = rDiff * rDiff;
                                int norm2 = gDiff * gDiff;
                                int norm3 = bDiff * bDiff;
                                int norm4 = aDiff * aDiff;

                                int sum1 = norm1 + norm2;
                                int sum2 = norm3 + norm4;

                                blocks[srcY * maxBlockX + srcX].value += sum1 + sum2;
                                srcRow += 4;
                            }
                        }
                    }
                }
            }
        }
        // Compute the vertical overlap
        for (int i = 0; i < bHeight; i++) {
            unsigned char* outputRow = mData->output_d[dstY + i] + dstXStart;
            for (int blockYStart = 0; blockYStart < maxBlockY + blockSize; blockYStart += blockSize) {
                for (int blockXStart = 0; blockXStart < maxBlockX + blockSize; blockXStart += blockSize) {
                    for (int j = 0; j < overlapWidth; j++) {
                        int rDst = outputRow[j * 4];
                        int gDst = outputRow[j * 4 + 1];
                        int bDst = outputRow[j * 4 + 2];
                        int aDst = outputRow[j * 4 + 3];
                        for (int srcY = blockYStart; srcY < std::min(blockYStart + blockSize, maxBlockY);
                             srcY++) {
                            int srcYStart = srcY + overlapHeight;
                            unsigned char* srcRow = mData->data[srcYStart + i] + CHANNEL_NUM * blockXStart;
                            for (int srcX = blockXStart; srcX < std::min(blockXStart + blockSize, maxBlockX); srcX++) {
                                int rSrc = srcRow[j * 4];
                                int gSrc = srcRow[j * 4 + 1];
                                int bSrc = srcRow[j * 4 + 2];
                                int aSrc = srcRow[j * 4 + 3];

                                int rDiff = rDst - rSrc;
                                int gDiff = gDst - gSrc;
                                int bDiff = bDst - bSrc;
                                int aDiff = aDst - aSrc;

                                int norm1 = rDiff * rDiff;
                                int norm2 = gDiff * gDiff;
                                int norm3 = bDiff * bDiff;
                                int norm4 = aDiff * aDiff;

                                int sum1 = norm1 + norm2;
                                int sum2 = norm3 + norm4;

                                blocks[srcY * maxBlockX + srcX].value += sum1 + sum2;
                                srcRow += 4;
                            }
                        }
                    }
                }
            }
        }
    } else if (overlapType == vertical) {
        for (int i = 0; i < mData->block_h; i++) {
            unsigned char* outputRow = mData->output_d[dstY + i] + dstXStart;
            for (int blockYStart = 0; blockYStart < maxBlockY + blockSize; blockYStart += blockSize) {
                for (int blockXStart = 0; blockXStart < maxBlockX + blockSize; blockXStart += blockSize) {
                    for (int j = 0; j < overlapWidth; j++) {
                        int rDst = outputRow[j * 4];
                        int gDst = outputRow[j * 4 + 1];
                        int bDst = outputRow[j * 4 + 2];
                        int aDst = outputRow[j * 4 + 3];
                        for (int srcY = blockYStart; srcY < std::min(blockYStart + blockSize, maxBlockY);
                             srcY++) {
                            unsigned char* srcRow = mData->data[srcY + i] + CHANNEL_NUM * blockXStart;
                            for (int srcX = blockXStart; srcX < std::min(blockXStart + blockSize, maxBlockX); srcX++) {
                                int rSrc = srcRow[j * 4];
                                int gSrc = srcRow[j * 4 + 1];
                                int bSrc = srcRow[j * 4 + 2];
                                int aSrc = srcRow[j * 4 + 3];

                                int rDiff = rDst - rSrc;
                                int gDiff = gDst - gSrc;
                                int bDiff = bDst - bSrc;
                                int aDiff = aDst - aSrc;

                                int norm1 = rDiff * rDiff;
                                int norm2 = gDiff * gDiff;
                                int norm3 = bDiff * bDiff;
                                int norm4 = aDiff * aDiff;

                                int sum1 = norm1 + norm2;
                                int sum2 = norm3 + norm4;

                                blocks[srcY * maxBlockX + srcX].value += sum1 + sum2;
                                srcRow += 4;
                            }
                        }
                    }
                }
            }
        }
    } else if (overlapType == horizontal) {
        for (int i = 0; i < overlapHeight; i++) {
            unsigned char* outputRow = mData->output_d[overlapYStart + i] + dstXStart;
            for (int blockYStart = 0; blockYStart < maxBlockY + blockSize; blockYStart += blockSize) {
                for (int blockXStart = 0; blockXStart < maxBlockX + blockSize; blockXStart += blockSize) {
                    for (int j = 0; j < mData->block_w; j++) {
                        int rDst = outputRow[j * 4];
                        int gDst = outputRow[j * 4 + 1];
                        int bDst = outputRow[j * 4 + 2];
                        int aDst = outputRow[j * 4 + 3];

                        for (int srcY = blockYStart; srcY < std::min(blockYStart + blockSize, maxBlockY);
                             srcY++) {
                            unsigned char* srcRow = mData->data[srcY + i] + CHANNEL_NUM * blockXStart;
                            for (int srcX = blockXStart; srcX < std::min(blockXStart + blockSize, maxBlockX); srcX++) {
                                int rSrc = srcRow[j * 4];
                                int gSrc = srcRow[j * 4 + 1];
                                int bSrc = srcRow[j * 4 + 2];
                                int aSrc = srcRow[j * 4 + 3];

                                int rDiff = rDst - rSrc;
                                int gDiff = gDst - gSrc;
                                int bDiff = bDst - bSrc;
                                int aDiff = aDst - aSrc;

                                int norm1 = rDiff * rDiff;
                                int norm2 = gDiff * gDiff;
                                int norm3 = bDiff * bDiff;
                                int norm4 = aDiff * aDiff;

                                int sum1 = norm1 + norm2;
                                int sum2 = norm3 + norm4;

                                blocks[srcY * maxBlockX + srcX].value += sum1 + sum2;
                                srcRow += 4;
                            }
                        }
                    }
                }
            }
        }
    }
    for (int i = 0; i < numBlocks; i++) {
        if (blocks[i].value < minVal) {
            minVal = blocks[i].value;
        }
    }

    // Choose a random block within the tolerance
    double upperBound = (1.0 + errorTolerance) * std::sqrt(minVal);
    upperBound = upperBound * upperBound;
    BlockValue* suitableBlocks = (BlockValue*)malloc(sizeof(BlockValue) * numBlocks);
    int numSuitableBlocks = 0;
    for (int i = 0; i < numBlocks; i++) {
        if (blocks[i].value <= upperBound) {
            suitableBlocks[numSuitableBlocks] = blocks[i];
            numSuitableBlocks++;
        }
    }

    // Sample and place a block
    int blockIndex = GetRandomInt(0, numSuitableBlocks - 1);
    WriteBlockOverlapWithMinCut(overlapType, dstY, dstX, suitableBlocks[blockIndex].y,
                                suitableBlocks[blockIndex].x);

    // Clean up
    free(blocks);
    free(suitableBlocks);
}


// Synthesize a new texture by randomly choosing blocks satisfying constraints and applying minimum cuts
void AdvanceAlgOptimiz::OverlapConstraintsWithMinCut() {

    // Compute block parameters
    int hStep = mData->block_h - overlapHeight;
    int wStep = mData->block_w - overlapWidth;

    int maxBlockY = mData->height - mData->block_h;
    int maxBlockX = mData->width - mData->block_w;

    int regBlockW = mData->block_w;
    int regBlockH = mData->block_h - overlapHeight;

    // Randomly choose the upper-left corner of a block
    int srcY = GetRandomInt(0, maxBlockY - 1);
    int srcX = GetRandomInt(0, maxBlockX - 1);

    // Write the randomly chosen block to the output
    WriteBlock(0, 0, srcY, srcX);

    // fill first row
    int dstX = mData->block_w;
    int dstY = mData->block_h;
    for (; dstX < mData->output_w; dstX += wStep) {
        PlaceEdgeOverlapBlockWithMinCut(vertical, 0, dstX, maxBlockX, maxBlockY, 0.1, regBlockW, regBlockH);
    }
    int lastDstX = dstX - wStep;
    int blockWidth = mData->output_w - (lastDstX - overlapWidth);

    // fill all corner cases except borders
    for (; dstY < mData->output_h-hStep; dstY += hStep) {
        // fill first column
        PlaceEdgeOverlapBlockWithMinCut(horizontal, dstY,0, maxBlockX,
                                        maxBlockY, 0.1, regBlockW, regBlockH);
        dstX = mData->block_w;
        for (; dstX < mData->output_w-wStep; dstX += wStep) {
            PlaceEdgeOverlapBlockWithMinCut(both, dstY,dstX, maxBlockX, maxBlockY,
                                            0.1, regBlockW, regBlockH);
        }

        // fill last column
        PlaceEdgeOverlapBlockWithMinCut(both, dstY, dstX, maxBlockX, maxBlockY, 0.1,
                                        blockWidth, regBlockH);
    }

    // fill last row
    int blockHeight = mData->output_h - dstY;
    PlaceEdgeOverlapBlockWithMinCut(horizontal, dstY,0, maxBlockX, maxBlockY, 0.1, regBlockW, regBlockH);
    dstX = mData->block_w;
    for (; dstX < mData->output_w-wStep; dstX += wStep) {
        PlaceEdgeOverlapBlockWithMinCut(both, dstY, dstX, maxBlockX, maxBlockY, 0.1,
                                            mData->block_w, blockHeight);

    }
    // bottom-right corner
    PlaceEdgeOverlapBlockWithMinCut(both, dstY, lastDstX, maxBlockX, maxBlockY, 0.1,
                                        blockWidth, blockHeight);

}


void AdvanceAlgOptimiz::OverlapConstraintsWithMinCutBlocking() {
    // Compute block parameters
    int hStep = mData->block_h - overlapHeight;
    int wStep = mData->block_w - overlapWidth;

    int maxBlockY = mData->height - mData->block_h;
    int maxBlockX = mData->width - mData->block_w;

    int regBlockW = mData->block_w;
    int regBlockH = mData->block_h - overlapHeight;

    // Randomly choose the upper-left corner of a block
    int srcY = GetRandomInt(0, maxBlockY - 1);
    int srcX = GetRandomInt(0, maxBlockX - 1);

    // Write the randomly chosen block to the output
    WriteBlock(0, 0, srcY, srcX);

    // fill first row
    int dstX = mData->block_w;
    int dstY = mData->block_h;
    for (; dstX < mData->output_w; dstX += wStep) {
        PlaceEdgeOverlapBlockWithMinCutBlocking(vertical, 0, dstX, maxBlockX, maxBlockY, 0.1, regBlockW, regBlockH);
    }
    int lastDstX = dstX - wStep;
    int blockWidth = mData->output_w - (lastDstX - overlapWidth);

    // fill all corner cases except borders
    for (; dstY < mData->output_h-hStep; dstY += hStep) {
        // fill first column
        PlaceEdgeOverlapBlockWithMinCutBlocking(horizontal, dstY,0, maxBlockX,
                                        maxBlockY, 0.1, regBlockW, regBlockH);
        dstX = mData->block_w;
        for (; dstX < mData->output_w-wStep; dstX += wStep) {
            PlaceEdgeOverlapBlockWithMinCutBlocking(both, dstY,dstX, maxBlockX, maxBlockY,
                                            0.1, regBlockW, regBlockH);
        }

        // fill last column
        PlaceEdgeOverlapBlockWithMinCutBlocking(both, dstY, dstX, maxBlockX, maxBlockY, 0.1,
                                        blockWidth, regBlockH);
    }

    // fill last row
    int blockHeight = mData->output_h - dstY;
    PlaceEdgeOverlapBlockWithMinCutBlocking(horizontal, dstY,0, maxBlockX, maxBlockY, 0.1, regBlockW, regBlockH);
    dstX = mData->block_w;
    for (; dstX < mData->output_w-wStep; dstX += wStep) {
        PlaceEdgeOverlapBlockWithMinCutBlocking(both, dstY, dstX, maxBlockX, maxBlockY, 0.1,
                                            mData->block_w, blockHeight);

    }
    // bottom-right corner
    PlaceEdgeOverlapBlockWithMinCutBlocking(both, dstY, lastDstX, maxBlockX, maxBlockY, 0.1,
                                        blockWidth, blockHeight);
}

void AdvanceAlgOptimiz::OverlapConstraintsWithMinCutBlocking1() {
    // Compute block parameters
    int hStep = mData->block_h - overlapHeight;
    int wStep = mData->block_w - overlapWidth;

    int maxBlockY = mData->height - mData->block_h;
    int maxBlockX = mData->width - mData->block_w;

    int regBlockW = mData->block_w;
    int regBlockH = mData->block_h - overlapHeight;

    // Randomly choose the upper-left corner of a block
    int srcY = GetRandomInt(0, maxBlockY - 1);
    int srcX = GetRandomInt(0, maxBlockX - 1);

    // Write the randomly chosen block to the output
    WriteBlock(0, 0, srcY, srcX);

    // fill first row
    int dstX = mData->block_w;
    int dstY = mData->block_h;
    for (; dstX < mData->output_w; dstX += wStep) {
        PlaceEdgeOverlapBlockWithMinCutBlocking1(vertical, 0, dstX, maxBlockX, maxBlockY, 0.1, regBlockW, regBlockH);
    }
    int lastDstX = dstX - wStep;
    int blockWidth = mData->output_w - (lastDstX - overlapWidth);

    // fill all corner cases except borders
    for (; dstY < mData->output_h-hStep; dstY += hStep) {
        // fill first column
        PlaceEdgeOverlapBlockWithMinCutBlocking1(horizontal, dstY,0, maxBlockX,
                                                maxBlockY, 0.1, regBlockW, regBlockH);
        dstX = mData->block_w;
        for (; dstX < mData->output_w-wStep; dstX += wStep) {
            PlaceEdgeOverlapBlockWithMinCutBlocking1(both, dstY,dstX, maxBlockX, maxBlockY,
                                                    0.1, regBlockW, regBlockH);
        }

        // fill last column
        PlaceEdgeOverlapBlockWithMinCutBlocking1(both, dstY, dstX, maxBlockX, maxBlockY, 0.1,
                                                blockWidth, regBlockH);
    }

    // fill last row
    int blockHeight = mData->output_h - dstY;
    PlaceEdgeOverlapBlockWithMinCutBlocking1(horizontal, dstY,0, maxBlockX, maxBlockY, 0.1, regBlockW, regBlockH);
    dstX = mData->block_w;
    for (; dstX < mData->output_w-wStep; dstX += wStep) {
        PlaceEdgeOverlapBlockWithMinCutBlocking1(both, dstY, dstX, maxBlockX, maxBlockY, 0.1,
                                                mData->block_w, blockHeight);

    }
    // bottom-right corner
    PlaceEdgeOverlapBlockWithMinCutBlocking1(both, dstY, lastDstX, maxBlockX, maxBlockY, 0.1,
                                            blockWidth, blockHeight);
}

void AdvanceAlgOptimiz::BlockedFuncOpt(ImgData* imgData, int seed) {
    AdvanceAlgOptimiz imageQuilting(imgData);
    SeedRandomNumberGenerator(seed);
    imageQuilting.OverlapConstraintsWithMinCutBlocking1();
}

int64_t AdvanceAlgOptimiz::getFlopCount() const {
    return flopCount;
}
