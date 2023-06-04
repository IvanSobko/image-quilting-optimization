#include <algorithm>
#include <cstdlib>
#include <random>
#include <climits>

#include "AdvanceAlgOptimiz.h"

// Custom parameters for testing and tuning
void AdvanceAlgOptimiz::CustomParameters(ImgData* imgData) {
    imgData->output_h = 2 * imgData->height;
    imgData->output_w = 2 * imgData->width;
    imgData->block_h = imgData->height / 2;
    imgData->block_w = imgData->width / 2;
}

void AdvanceAlgOptimiz::StdC_KUnroll_BoundsRefactor(ImgData* imgData, int seed) {
    AdvanceAlgOptimiz imageQuilting(imgData);
    SeedRandomNumberGenerator(seed);
    imageQuilting.OverlapConstraintsWithMinCut_StdC_KUnroll_BoundsRefactor();
}

void AdvanceAlgOptimiz::StdC_KUnroll_BoundsRefactor_LoopReorder(ImgData* imgData, int seed) {
    AdvanceAlgOptimiz imageQuilting(imgData);
    SeedRandomNumberGenerator(seed);
    imageQuilting.OverlapConstraintsWithMinCut_StdC_KUnroll_BoundsRefactor_LoopReorder();
}

void AdvanceAlgOptimiz::StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking32(ImgData* imgData, int seed) {
    AdvanceAlgOptimiz imageQuilting(imgData);
    SeedRandomNumberGenerator(seed);
    imageQuilting.OverlapConstraintsWithMinCut_StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking32();
}

void AdvanceAlgOptimiz::StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking48(ImgData* imgData, int seed) {
    AdvanceAlgOptimiz imageQuilting(imgData);
    SeedRandomNumberGenerator(seed);
    imageQuilting.OverlapConstraintsWithMinCut_StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking48();
}

void AdvanceAlgOptimiz::StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking64(ImgData* imgData, int seed) {
    AdvanceAlgOptimiz imageQuilting(imgData);
    SeedRandomNumberGenerator(seed);
    imageQuilting.OverlapConstraintsWithMinCut_StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking64();
}

void AdvanceAlgOptimiz::StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking96(ImgData* imgData, int seed) {
    AdvanceAlgOptimiz imageQuilting(imgData);
    SeedRandomNumberGenerator(seed);
    imageQuilting.OverlapConstraintsWithMinCut_StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking96();
}

void AdvanceAlgOptimiz::StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking128(ImgData* imgData, int seed) {
    AdvanceAlgOptimiz imageQuilting(imgData);
    SeedRandomNumberGenerator(seed);
    imageQuilting.OverlapConstraintsWithMinCut_StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking128();
}

void AdvanceAlgOptimiz::StdC_KSrc2Unroll_BoundsRefactor_LoopReorder_Blocking32(ImgData * imgData, int seed) {
    AdvanceAlgOptimiz imageQuilting(imgData);
    SeedRandomNumberGenerator(seed);
    imageQuilting.OverlapConstraintsWithMinCut_StdC_KSrc2Unroll_BoundsRefactor_LoopReorder_Blocking32();
}

void AdvanceAlgOptimiz::StdC_KSrc4Unroll_BoundsRefactor_LoopReorder_Blocking32(ImgData * imgData, int seed) {
    AdvanceAlgOptimiz imageQuilting(imgData);
    SeedRandomNumberGenerator(seed);
    imageQuilting.OverlapConstraintsWithMinCut_StdC_KSrc4Unroll_BoundsRefactor_LoopReorder_Blocking32();
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
                if (leftError < localError) {
                    localError = leftError;
                    verticalPath[i] = j - 1;
                }
            }
            // Get the value to the right
            if (j < overlapWidthLocal - 1) {
                double rightError = dpTable[i * overlapWidthLocal + j + 1];
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
                if (leftError < localError) {
                    localError = leftError;
                    horizontalPath[j] = i - 1;
                }
            }
            // Get the value to the left and below
            if (i < overlapHeightLocal - 1) {
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


// Std C, K unroll, and bounds refactoring optimizations
void AdvanceAlgOptimiz::PlaceEdgeOverlapBlockWithMinCut_StdC_KUnroll_BoundsRefactor(int overlapType, const int dstY, const int dstX,
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
void AdvanceAlgOptimiz::OverlapConstraintsWithMinCut_StdC_KUnroll_BoundsRefactor() {

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
        PlaceEdgeOverlapBlockWithMinCut_StdC_KUnroll_BoundsRefactor(vertical, 0, dstX, maxBlockX, maxBlockY,
                                                                    0.1, regBlockW, regBlockH);
    }
    int lastDstX = dstX - wStep;
    int blockWidth = mData->output_w - (lastDstX - overlapWidth);

    // fill all corner cases except borders
    for (; dstY < mData->output_h-hStep; dstY += hStep) {
        // fill first column
        PlaceEdgeOverlapBlockWithMinCut_StdC_KUnroll_BoundsRefactor(horizontal, dstY, 0, maxBlockX, maxBlockY,
                                                                    0.1, regBlockW, regBlockH);
        dstX = mData->block_w;
        for (; dstX < mData->output_w-wStep; dstX += wStep) {
            PlaceEdgeOverlapBlockWithMinCut_StdC_KUnroll_BoundsRefactor(both, dstY, dstX, maxBlockX,
                                                                        maxBlockY, 0.1, regBlockW, regBlockH);
        }

        // fill last column
        PlaceEdgeOverlapBlockWithMinCut_StdC_KUnroll_BoundsRefactor(both, dstY, dstX, maxBlockX, maxBlockY,
                                                                    0.1, blockWidth, regBlockH);
    }

    // fill last row
    int blockHeight = mData->output_h - dstY;
    PlaceEdgeOverlapBlockWithMinCut_StdC_KUnroll_BoundsRefactor(horizontal, dstY, 0, maxBlockX, maxBlockY,
                                                                0.1, regBlockW, regBlockH);
    dstX = mData->block_w;
    for (; dstX < mData->output_w-wStep; dstX += wStep) {
        PlaceEdgeOverlapBlockWithMinCut_StdC_KUnroll_BoundsRefactor(both, dstY, dstX, maxBlockX, maxBlockY,
                                                                    0.1, mData->block_w, blockHeight);

    }
    // bottom-right corner
    PlaceEdgeOverlapBlockWithMinCut_StdC_KUnroll_BoundsRefactor(both, dstY, lastDstX, maxBlockX, maxBlockY,
                                                                0.1, blockWidth, blockHeight);

}

// Std C, K unroll, bounds refactoring, and loop reorder optimizations
void AdvanceAlgOptimiz::PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KUnroll_BoundsRefactor_LoopReorder(int overlapType, int dstY, int dstX, int maxBlockX,
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
void AdvanceAlgOptimiz::OverlapConstraintsWithMinCut_StdC_KUnroll_BoundsRefactor_LoopReorder() {
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
        PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KUnroll_BoundsRefactor_LoopReorder(
            vertical, 0, dstX, maxBlockX, maxBlockY, 0.1, regBlockW, regBlockH);
    }
    int lastDstX = dstX - wStep;
    int blockWidth = mData->output_w - (lastDstX - overlapWidth);

    // fill all corner cases except borders
    for (; dstY < mData->output_h-hStep; dstY += hStep) {
        // fill first column
        PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KUnroll_BoundsRefactor_LoopReorder(
            horizontal, dstY, 0, maxBlockX, maxBlockY, 0.1, regBlockW, regBlockH);
        dstX = mData->block_w;
        for (; dstX < mData->output_w-wStep; dstX += wStep) {
            PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KUnroll_BoundsRefactor_LoopReorder(
                both, dstY, dstX, maxBlockX, maxBlockY, 0.1, regBlockW, regBlockH);
        }

        // fill last column
        PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KUnroll_BoundsRefactor_LoopReorder(
            both, dstY, dstX, maxBlockX, maxBlockY, 0.1, blockWidth, regBlockH);
    }

    // fill last row
    int blockHeight = mData->output_h - dstY;
    PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KUnroll_BoundsRefactor_LoopReorder(
        horizontal, dstY, 0, maxBlockX, maxBlockY, 0.1, regBlockW, regBlockH);
    dstX = mData->block_w;
    for (; dstX < mData->output_w-wStep; dstX += wStep) {
        PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KUnroll_BoundsRefactor_LoopReorder(
            both, dstY, dstX, maxBlockX, maxBlockY, 0.1, mData->block_w, blockHeight);

    }
    // bottom-right corner
    PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KUnroll_BoundsRefactor_LoopReorder(
        both, dstY, lastDstX, maxBlockX, maxBlockY, 0.1, blockWidth, blockHeight);
}

// Std C, K unroll, bounds refactoring, loop reorder optimizations, and blocking with 32x32 blocks
void AdvanceAlgOptimiz::
    PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking32(int overlapType, int dstY, int dstX, int maxBlockX,
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
void AdvanceAlgOptimiz::OverlapConstraintsWithMinCut_StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking32() {
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
        PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking32(
            vertical, 0, dstX, maxBlockX, maxBlockY, 0.1, regBlockW, regBlockH);
    }
    int lastDstX = dstX - wStep;
    int blockWidth = mData->output_w - (lastDstX - overlapWidth);

    // fill all corner cases except borders
    for (; dstY < mData->output_h-hStep; dstY += hStep) {
        // fill first column
        PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking32(
            horizontal, dstY, 0, maxBlockX, maxBlockY, 0.1, regBlockW, regBlockH);
        dstX = mData->block_w;
        for (; dstX < mData->output_w-wStep; dstX += wStep) {
            PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking32(
                both, dstY, dstX, maxBlockX, maxBlockY, 0.1, regBlockW, regBlockH);
        }

        // fill last column
        PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking32(
            both, dstY, dstX, maxBlockX, maxBlockY, 0.1, blockWidth, regBlockH);
    }

    // fill last row
    int blockHeight = mData->output_h - dstY;
    PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking32(
        horizontal, dstY, 0, maxBlockX, maxBlockY, 0.1, regBlockW, regBlockH);
    dstX = mData->block_w;
    for (; dstX < mData->output_w-wStep; dstX += wStep) {
        PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking32(
            both, dstY, dstX, maxBlockX, maxBlockY, 0.1, mData->block_w, blockHeight);

    }
    // bottom-right corner
    PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking32(
        both, dstY, lastDstX, maxBlockX, maxBlockY, 0.1, blockWidth, blockHeight);
}

// Std C, K unroll, bounds refactoring, loop reorder optimizations, and blocking with 48x48 blocks
void AdvanceAlgOptimiz::
    PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking48(int overlapType, int dstY, int dstX, int maxBlockX,
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
    const int blockSize = 48;

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
void AdvanceAlgOptimiz::OverlapConstraintsWithMinCut_StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking48() {
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
        PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking48(
            vertical, 0, dstX, maxBlockX, maxBlockY, 0.1, regBlockW, regBlockH);
    }
    int lastDstX = dstX - wStep;
    int blockWidth = mData->output_w - (lastDstX - overlapWidth);

    // fill all corner cases except borders
    for (; dstY < mData->output_h-hStep; dstY += hStep) {
        // fill first column
        PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking48(
            horizontal, dstY, 0, maxBlockX, maxBlockY, 0.1, regBlockW, regBlockH);
        dstX = mData->block_w;
        for (; dstX < mData->output_w-wStep; dstX += wStep) {
            PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking48(
                both, dstY, dstX, maxBlockX, maxBlockY, 0.1, regBlockW, regBlockH);
        }

        // fill last column
        PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking48(
            both, dstY, dstX, maxBlockX, maxBlockY, 0.1, blockWidth, regBlockH);
    }

    // fill last row
    int blockHeight = mData->output_h - dstY;
    PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking48(
        horizontal, dstY, 0, maxBlockX, maxBlockY, 0.1, regBlockW, regBlockH);
    dstX = mData->block_w;
    for (; dstX < mData->output_w-wStep; dstX += wStep) {
        PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking48(
            both, dstY, dstX, maxBlockX, maxBlockY, 0.1, mData->block_w, blockHeight);

    }
    // bottom-right corner
    PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking48(
        both, dstY, lastDstX, maxBlockX, maxBlockY, 0.1, blockWidth, blockHeight);
}

// Std C, K unroll, bounds refactoring, loop reorder optimizations, and blocking with 64x64 blocks
void AdvanceAlgOptimiz::
    PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking64(int overlapType, int dstY, int dstX, int maxBlockX,
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
    const int blockSize = 64;

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
void AdvanceAlgOptimiz::OverlapConstraintsWithMinCut_StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking64() {
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
        PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking64(
            vertical, 0, dstX, maxBlockX, maxBlockY, 0.1, regBlockW, regBlockH);
    }
    int lastDstX = dstX - wStep;
    int blockWidth = mData->output_w - (lastDstX - overlapWidth);

    // fill all corner cases except borders
    for (; dstY < mData->output_h-hStep; dstY += hStep) {
        // fill first column
        PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking64(
            horizontal, dstY, 0, maxBlockX, maxBlockY, 0.1, regBlockW, regBlockH);
        dstX = mData->block_w;
        for (; dstX < mData->output_w-wStep; dstX += wStep) {
            PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking64(
                both, dstY, dstX, maxBlockX, maxBlockY, 0.1, regBlockW, regBlockH);
        }

        // fill last column
        PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking64(
            both, dstY, dstX, maxBlockX, maxBlockY, 0.1, blockWidth, regBlockH);
    }

    // fill last row
    int blockHeight = mData->output_h - dstY;
    PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking64(
        horizontal, dstY, 0, maxBlockX, maxBlockY, 0.1, regBlockW, regBlockH);
    dstX = mData->block_w;
    for (; dstX < mData->output_w-wStep; dstX += wStep) {
        PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking64(
            both, dstY, dstX, maxBlockX, maxBlockY, 0.1, mData->block_w, blockHeight);

    }
    // bottom-right corner
    PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking64(
        both, dstY, lastDstX, maxBlockX, maxBlockY, 0.1, blockWidth, blockHeight);
}

// Std C, K unroll, bounds refactoring, loop reorder optimizations, and blocking with 96x96 blocks
void AdvanceAlgOptimiz::
    PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking96(int overlapType, int dstY, int dstX, int maxBlockX,
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
    const int blockSize = 96;

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
void AdvanceAlgOptimiz::OverlapConstraintsWithMinCut_StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking96() {
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
        PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking96(
            vertical, 0, dstX, maxBlockX, maxBlockY, 0.1, regBlockW, regBlockH);
    }
    int lastDstX = dstX - wStep;
    int blockWidth = mData->output_w - (lastDstX - overlapWidth);

    // fill all corner cases except borders
    for (; dstY < mData->output_h-hStep; dstY += hStep) {
        // fill first column
        PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking96(
            horizontal, dstY, 0, maxBlockX, maxBlockY, 0.1, regBlockW, regBlockH);
        dstX = mData->block_w;
        for (; dstX < mData->output_w-wStep; dstX += wStep) {
            PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking96(
                both, dstY, dstX, maxBlockX, maxBlockY, 0.1, regBlockW, regBlockH);
        }

        // fill last column
        PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking96(
            both, dstY, dstX, maxBlockX, maxBlockY, 0.1, blockWidth, regBlockH);
    }

    // fill last row
    int blockHeight = mData->output_h - dstY;
    PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking96(
        horizontal, dstY, 0, maxBlockX, maxBlockY, 0.1, regBlockW, regBlockH);
    dstX = mData->block_w;
    for (; dstX < mData->output_w-wStep; dstX += wStep) {
        PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking96(
            both, dstY, dstX, maxBlockX, maxBlockY, 0.1, mData->block_w, blockHeight);

    }
    // bottom-right corner
    PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking96(
        both, dstY, lastDstX, maxBlockX, maxBlockY, 0.1, blockWidth, blockHeight);
}

// Std C, K unroll, bounds refactoring, loop reorder optimizations, and blocking with 128x128 blocks
void AdvanceAlgOptimiz::
    PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking128(int overlapType, int dstY, int dstX, int maxBlockX,
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
    const int blockSize = 128;

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
void AdvanceAlgOptimiz::OverlapConstraintsWithMinCut_StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking128() {
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
        PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking128(
            vertical, 0, dstX, maxBlockX, maxBlockY, 0.1, regBlockW, regBlockH);
    }
    int lastDstX = dstX - wStep;
    int blockWidth = mData->output_w - (lastDstX - overlapWidth);

    // fill all corner cases except borders
    for (; dstY < mData->output_h-hStep; dstY += hStep) {
        // fill first column
        PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking128(
            horizontal, dstY, 0, maxBlockX, maxBlockY, 0.1, regBlockW, regBlockH);
        dstX = mData->block_w;
        for (; dstX < mData->output_w-wStep; dstX += wStep) {
            PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking128(
                both, dstY, dstX, maxBlockX, maxBlockY, 0.1, regBlockW, regBlockH);
        }

        // fill last column
        PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking128(
            both, dstY, dstX, maxBlockX, maxBlockY, 0.1, blockWidth, regBlockH);
    }

    // fill last row
    int blockHeight = mData->output_h - dstY;
    PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking128(
        horizontal, dstY, 0, maxBlockX, maxBlockY, 0.1, regBlockW, regBlockH);
    dstX = mData->block_w;
    for (; dstX < mData->output_w-wStep; dstX += wStep) {
        PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking128(
            both, dstY, dstX, maxBlockX, maxBlockY, 0.1, mData->block_w, blockHeight);

    }
    // bottom-right corner
    PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking128(
        both, dstY, lastDstX, maxBlockX, maxBlockY, 0.1, blockWidth, blockHeight);
}


// Std C, bounds refactoring, loop reorder, blocking 32x32, and unrolling channels loop and srcX by 2
void AdvanceAlgOptimiz::
    PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KSrc2Unroll_BoundsRefactor_LoopReorder_Blocking32(int overlapType, int dstY, int dstX, int maxBlockX,
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
                        int rDst1 = outputRow[j * 4];
                        int gDst1 = outputRow[j * 4 + 1];
                        int bDst1 = outputRow[j * 4 + 2];
                        int aDst1 = outputRow[j * 4 + 3];

                        for (int srcY = blockYStart; srcY < std::min(blockYStart + blockSize, maxBlockY);
                             srcY++) {
                            unsigned char* srcRow = mData->data[srcY + i] + CHANNEL_NUM * blockXStart;
                            int srcX = blockXStart;
                            int endX = std::min(blockXStart + blockSize, maxBlockX);
                            int blockIndex = srcY * maxBlockX + srcX;
                            for (; srcX < endX-1; srcX+=2) {
                                int rSrc1 = srcRow[j * 4];
                                int gSrc1 = srcRow[j * 4 + 1];
                                int bSrc1 = srcRow[j * 4 + 2];
                                int aSrc1 = srcRow[j * 4 + 3];
                                int rSrc2 = srcRow[j * 4 + 4];
                                int gSrc2 = srcRow[j * 4 + 5];
                                int bSrc2 = srcRow[j * 4 + 6];
                                int aSrc2 = srcRow[j * 4 + 7];

                                int rDiff1 = rDst1 - rSrc1;
                                int gDiff1 = gDst1 - gSrc1;
                                int bDiff1 = bDst1 - bSrc1;
                                int aDiff1 = aDst1 - aSrc1;

                                int norm11 = rDiff1 * rDiff1;
                                int rDiff2 = rDst1 - rSrc2;
                                int gDiff2 = gDst1 - gSrc2;
                                int bDiff2 = bDst1 - bSrc2;

                                int norm12 = gDiff1 * gDiff1;
                                int aDiff2 = aDst1 - aSrc2;
                                int norm13 = bDiff1 * bDiff1;
                                int norm14 = aDiff1 * aDiff1;

                                int norm21 = rDiff2 * rDiff2;
                                int norm22 = gDiff2 * gDiff2;
                                int norm23 = bDiff2 * bDiff2;
                                int norm24 = aDiff2 * aDiff2;

                                int sum1 = norm11 + norm12;
                                int sum2 = norm13 + norm14;
                                int sum3 = norm21 + norm22;
                                int sum4 = norm23 + norm24;

                                blocks[blockIndex++].value += sum1 + sum2;
                                blocks[blockIndex++].value += sum3 + sum4;
                                srcRow += 8;
                            }
                            for (; srcX < endX; srcX++) {
                                int rSrc1 = srcRow[j * 4];
                                int gSrc1 = srcRow[j * 4 + 1];
                                int bSrc1 = srcRow[j * 4 + 2];
                                int aSrc1 = srcRow[j * 4 + 3];

                                int rDiff1 = rDst1 - rSrc1;
                                int gDiff1 = gDst1 - gSrc1;
                                int bDiff1 = bDst1 - bSrc1;
                                int aDiff1 = aDst1 - aSrc1;

                                int norm11 = rDiff1 * rDiff1;
                                int norm12 = gDiff1 * gDiff1;
                                int norm13 = bDiff1 * bDiff1;
                                int norm14 = aDiff1 * aDiff1;

                                int sum1 = norm11 + norm12;
                                int sum2 = norm13 + norm14;

                                blocks[blockIndex++].value += sum1 + sum2;
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
                    int j = 0;
                    for (; j < overlapWidth; j++) {
                        int rDst1 = outputRow[j * 4];
                        int gDst1 = outputRow[j * 4 + 1];
                        int bDst1 = outputRow[j * 4 + 2];
                        int aDst1 = outputRow[j * 4 + 3];
                        for (int srcY = blockYStart; srcY < std::min(blockYStart + blockSize, maxBlockY);
                             srcY++) {
                            int srcYStart = srcY + overlapHeight;
                            unsigned char* srcRow = mData->data[srcYStart + i] + CHANNEL_NUM * blockXStart;
                            int srcX = blockXStart;
                            int endX = std::min(blockXStart + blockSize, maxBlockX);
                            int blockIndex = srcY * maxBlockX + srcX;
                            for (; srcX < endX-1; srcX+=2) {
                                int rSrc1 = srcRow[j * 4];
                                int gSrc1 = srcRow[j * 4 + 1];
                                int bSrc1 = srcRow[j * 4 + 2];
                                int aSrc1 = srcRow[j * 4 + 3];
                                int rSrc2 = srcRow[j * 4 + 4];
                                int gSrc2 = srcRow[j * 4 + 5];
                                int bSrc2 = srcRow[j * 4 + 6];
                                int aSrc2 = srcRow[j * 4 + 7];

                                int rDiff1 = rDst1 - rSrc1;
                                int gDiff1 = gDst1 - gSrc1;
                                int bDiff1 = bDst1 - bSrc1;
                                int aDiff1 = aDst1 - aSrc1;

                                int rDiff2 = rDst1 - rSrc2;
                                int gDiff2 = gDst1 - gSrc2;
                                int bDiff2 = bDst1 - bSrc2;
                                int aDiff2 = aDst1 - aSrc2;

                                int norm11 = rDiff1 * rDiff1;
                                int norm12 = gDiff1 * gDiff1;
                                int norm13 = bDiff1 * bDiff1;
                                int norm14 = aDiff1 * aDiff1;

                                int norm21 = rDiff2 * rDiff2;
                                int norm22 = gDiff2 * gDiff2;
                                int norm23 = bDiff2 * bDiff2;
                                int norm24 = aDiff2 * aDiff2;

                                int sum1 = norm11 + norm12;
                                int sum2 = norm13 + norm14;
                                int sum3 = norm21 + norm22;
                                int sum4 = norm23 + norm24;

                                blocks[blockIndex++].value += sum1 + sum2;
                                blocks[blockIndex++].value += sum3 + sum4;
                                srcRow += 8;
                            }
                            for (; srcX < endX; srcX++) {
                                int rSrc1 = srcRow[j * 4];
                                int gSrc1 = srcRow[j * 4 + 1];
                                int bSrc1 = srcRow[j * 4 + 2];
                                int aSrc1 = srcRow[j * 4 + 3];

                                int rDiff1 = rDst1 - rSrc1;
                                int gDiff1 = gDst1 - gSrc1;
                                int bDiff1 = bDst1 - bSrc1;
                                int aDiff1 = aDst1 - aSrc1;

                                int norm11 = rDiff1 * rDiff1;
                                int norm12 = gDiff1 * gDiff1;
                                int norm13 = bDiff1 * bDiff1;
                                int norm14 = aDiff1 * aDiff1;

                                int sum1 = norm11 + norm12;
                                int sum2 = norm13 + norm14;

                                blocks[blockIndex++].value += sum1 + sum2;
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
                        int rDst1 = outputRow[j * 4];
                        int gDst1 = outputRow[j * 4 + 1];
                        int bDst1 = outputRow[j * 4 + 2];
                        int aDst1 = outputRow[j * 4 + 3];
                        for (int srcY = blockYStart; srcY < std::min(blockYStart + blockSize, maxBlockY);
                             srcY++) {
                            unsigned char* srcRow = mData->data[srcY + i] + CHANNEL_NUM * blockXStart;
                            int srcX = blockXStart;
                            int endX = std::min(blockXStart + blockSize, maxBlockX);
                            int blockIndex = srcY * maxBlockX + srcX;
                            for (; srcX < endX-1; srcX+=2) {
                                int rSrc1 = srcRow[j * 4];
                                int gSrc1 = srcRow[j * 4 + 1];
                                int bSrc1 = srcRow[j * 4 + 2];
                                int aSrc1 = srcRow[j * 4 + 3];
                                int rSrc2 = srcRow[j * 4 + 4];
                                int gSrc2 = srcRow[j * 4 + 5];
                                int bSrc2 = srcRow[j * 4 + 6];
                                int aSrc2 = srcRow[j * 4 + 7];

                                int rDiff1 = rDst1 - rSrc1;
                                int gDiff1 = gDst1 - gSrc1;
                                int bDiff1 = bDst1 - bSrc1;
                                int aDiff1 = aDst1 - aSrc1;

                                int norm11 = rDiff1 * rDiff1;
                                int rDiff2 = rDst1 - rSrc2;
                                int gDiff2 = gDst1 - gSrc2;
                                int bDiff2 = bDst1 - bSrc2;

                                int norm12 = gDiff1 * gDiff1;
                                int aDiff2 = aDst1 - aSrc2;
                                int norm13 = bDiff1 * bDiff1;
                                int norm14 = aDiff1 * aDiff1;

                                int norm21 = rDiff2 * rDiff2;
                                int norm22 = gDiff2 * gDiff2;
                                int norm23 = bDiff2 * bDiff2;
                                int norm24 = aDiff2 * aDiff2;

                                int sum1 = norm11 + norm12;
                                int sum2 = norm13 + norm14;
                                int sum3 = norm21 + norm22;
                                int sum4 = norm23 + norm24;

                                blocks[blockIndex++].value += sum1 + sum2;
                                blocks[blockIndex++].value += sum3 + sum4;
                                srcRow += 8;
                            }
                            for (; srcX < endX; srcX++) {
                                int rSrc1 = srcRow[j * 4];
                                int gSrc1 = srcRow[j * 4 + 1];
                                int bSrc1 = srcRow[j * 4 + 2];
                                int aSrc1 = srcRow[j * 4 + 3];

                                int rDiff1 = rDst1 - rSrc1;
                                int gDiff1 = gDst1 - gSrc1;
                                int bDiff1 = bDst1 - bSrc1;
                                int aDiff1 = aDst1 - aSrc1;

                                int norm11 = rDiff1 * rDiff1;
                                int norm12 = gDiff1 * gDiff1;
                                int norm13 = bDiff1 * bDiff1;
                                int norm14 = aDiff1 * aDiff1;

                                int sum1 = norm11 + norm12;
                                int sum2 = norm13 + norm14;

                                blocks[blockIndex++].value += sum1 + sum2;
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
                        int rDst1 = outputRow[j * 4];
                        int gDst1 = outputRow[j * 4 + 1];
                        int bDst1 = outputRow[j * 4 + 2];
                        int aDst1 = outputRow[j * 4 + 3];

                        for (int srcY = blockYStart; srcY < std::min(blockYStart + blockSize, maxBlockY);
                             srcY++) {
                            unsigned char* srcRow = mData->data[srcY + i] + CHANNEL_NUM * blockXStart;
                            int srcX = blockXStart;
                            int endX = std::min(blockXStart + blockSize, maxBlockX);
                            int blockIndex = srcY * maxBlockX + srcX;
                            for (; srcX < endX-1; srcX+=2) {
                                int rSrc1 = srcRow[j * 4];
                                int gSrc1 = srcRow[j * 4 + 1];
                                int bSrc1 = srcRow[j * 4 + 2];
                                int aSrc1 = srcRow[j * 4 + 3];
                                int rSrc2 = srcRow[j * 4 + 4];
                                int gSrc2 = srcRow[j * 4 + 5];
                                int bSrc2 = srcRow[j * 4 + 6];
                                int aSrc2 = srcRow[j * 4 + 7];

                                int rDiff1 = rDst1 - rSrc1;
                                int gDiff1 = gDst1 - gSrc1;
                                int bDiff1 = bDst1 - bSrc1;
                                int aDiff1 = aDst1 - aSrc1;

                                int norm11 = rDiff1 * rDiff1;
                                int rDiff2 = rDst1 - rSrc2;
                                int gDiff2 = gDst1 - gSrc2;
                                int bDiff2 = bDst1 - bSrc2;

                                int norm12 = gDiff1 * gDiff1;
                                int aDiff2 = aDst1 - aSrc2;
                                int norm13 = bDiff1 * bDiff1;
                                int norm14 = aDiff1 * aDiff1;

                                int norm21 = rDiff2 * rDiff2;
                                int norm22 = gDiff2 * gDiff2;
                                int norm23 = bDiff2 * bDiff2;
                                int norm24 = aDiff2 * aDiff2;

                                int sum1 = norm11 + norm12;
                                int sum2 = norm13 + norm14;
                                int sum3 = norm21 + norm22;
                                int sum4 = norm23 + norm24;

                                blocks[blockIndex++].value += sum1 + sum2;
                                blocks[blockIndex++].value += sum3 + sum4;
                                srcRow += 8;
                            }
                            for (; srcX < endX; srcX++) {
                                int rSrc1 = srcRow[j * 4];
                                int gSrc1 = srcRow[j * 4 + 1];
                                int bSrc1 = srcRow[j * 4 + 2];
                                int aSrc1 = srcRow[j * 4 + 3];

                                int rDiff1 = rDst1 - rSrc1;
                                int gDiff1 = gDst1 - gSrc1;
                                int bDiff1 = bDst1 - bSrc1;
                                int aDiff1 = aDst1 - aSrc1;

                                int norm11 = rDiff1 * rDiff1;
                                int norm12 = gDiff1 * gDiff1;
                                int norm13 = bDiff1 * bDiff1;
                                int norm14 = aDiff1 * aDiff1;

                                int sum1 = norm11 + norm12;
                                int sum2 = norm13 + norm14;

                                blocks[blockIndex++].value += sum1 + sum2;
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

void AdvanceAlgOptimiz::OverlapConstraintsWithMinCut_StdC_KSrc2Unroll_BoundsRefactor_LoopReorder_Blocking32() {
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
        PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KSrc2Unroll_BoundsRefactor_LoopReorder_Blocking32(
                vertical, 0, dstX, maxBlockX, maxBlockY, 0.1, regBlockW, regBlockH);
    }
    int lastDstX = dstX - wStep;
    int blockWidth = mData->output_w - (lastDstX - overlapWidth);

    // fill all corner cases except borders
    for (; dstY < mData->output_h-hStep; dstY += hStep) {
        // fill first column
        PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KSrc2Unroll_BoundsRefactor_LoopReorder_Blocking32(
                horizontal, dstY, 0, maxBlockX, maxBlockY, 0.1, regBlockW, regBlockH);
        dstX = mData->block_w;
        for (; dstX < mData->output_w-wStep; dstX += wStep) {
            PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KSrc2Unroll_BoundsRefactor_LoopReorder_Blocking32(
                    both, dstY, dstX, maxBlockX, maxBlockY, 0.1, regBlockW, regBlockH);
        }

        // fill last column
        PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KSrc2Unroll_BoundsRefactor_LoopReorder_Blocking32(
                both, dstY, dstX, maxBlockX, maxBlockY, 0.1, blockWidth, regBlockH);
    }

    // fill last row
    int blockHeight = mData->output_h - dstY;
    PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KSrc2Unroll_BoundsRefactor_LoopReorder_Blocking32(
            horizontal, dstY, 0, maxBlockX, maxBlockY, 0.1, regBlockW, regBlockH);
    dstX = mData->block_w;
    for (; dstX < mData->output_w-wStep; dstX += wStep) {
        PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KSrc2Unroll_BoundsRefactor_LoopReorder_Blocking32(
                both, dstY, dstX, maxBlockX, maxBlockY, 0.1, mData->block_w, blockHeight);

    }
    // bottom-right corner
    PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KSrc2Unroll_BoundsRefactor_LoopReorder_Blocking32(
            both, dstY, lastDstX, maxBlockX, maxBlockY, 0.1, blockWidth, blockHeight);
}

// Std C, bounds refactoring, loop reorder, blocking 32x32, and unrolling channels loop and srcX by 2
void AdvanceAlgOptimiz::
PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KSrc4Unroll_BoundsRefactor_LoopReorder_Blocking32(int overlapType, int dstY, int dstX, int maxBlockX,
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
                        int j4 = j * 4;
                        int rDst1 = outputRow[j4];
                        int gDst1 = outputRow[j4 + 1];
                        int bDst1 = outputRow[j4 + 2];
                        int aDst1 = outputRow[j4 + 3];

                        for (int srcY = blockYStart; srcY < std::min(blockYStart + blockSize, maxBlockY);
                             srcY++) {
                            unsigned char* srcRow = mData->data[srcY + i] + CHANNEL_NUM * blockXStart;
                            int srcX = blockXStart;
                            int endX = std::min(blockXStart + blockSize, maxBlockX);
                            int blockIndex = srcY * maxBlockX + srcX;
                            for (; srcX < endX-3; srcX+=4) {
                                int rSrc1 = srcRow[j4];
                                int gSrc1 = srcRow[j4 + 1];
                                int bSrc1 = srcRow[j4 + 2];
                                int aSrc1 = srcRow[j4 + 3];
                                int rSrc2 = srcRow[j4 + 4];
                                int gSrc2 = srcRow[j4 + 5];
                                int bSrc2 = srcRow[j4 + 6];
                                int aSrc2 = srcRow[j4 + 7];
                                int rSrc3 = srcRow[j4 + 8];
                                int gSrc3 = srcRow[j4 + 9];
                                int bSrc3 = srcRow[j4 + 10];
                                int aSrc3 = srcRow[j4 + 11];
                                int rSrc4 = srcRow[j4 + 12];
                                int gSrc4 = srcRow[j4 + 13];
                                int bSrc4 = srcRow[j4 + 14];
                                int aSrc4 = srcRow[j4 + 15];

                                int rDiff1 = rDst1 - rSrc1;
                                int gDiff1 = gDst1 - gSrc1;
                                int bDiff1 = bDst1 - bSrc1;
                                int aDiff1 = aDst1 - aSrc1;

                                int norm11 = rDiff1 * rDiff1;
                                int rDiff2 = rDst1 - rSrc2;
                                int gDiff2 = gDst1 - gSrc2;
                                int bDiff2 = bDst1 - bSrc2;

                                int norm12 = gDiff1 * gDiff1;
                                int aDiff2 = aDst1 - aSrc2;
                                int rDiff3 = rDst1 - rSrc3;
                                int gDiff3 = gDst1 - gSrc3;

                                int norm13 = bDiff1 * bDiff1;
                                int bDiff3 = bDst1 - bSrc3;
                                int aDiff3 = aDst1 - aSrc3;
                                int rDiff4 = rDst1 - rSrc4;

                                int norm14 = aDiff1 * aDiff1;
                                int gDiff4 = gDst1 - gSrc4;
                                int bDiff4 = bDst1 - bSrc4;
                                int aDiff4 = aDst1 - aSrc4;

                                int norm21 = rDiff2 * rDiff2;
                                int norm22 = gDiff2 * gDiff2;
                                int norm23 = bDiff2 * bDiff2;
                                int norm24 = aDiff2 * aDiff2;

                                int norm31 = rDiff3 * rDiff3;
                                int norm32 = gDiff3 * gDiff3;
                                int norm33 = bDiff3 * bDiff3;
                                int norm34 = aDiff3 * aDiff3;

                                int norm41 = rDiff4 * rDiff4;
                                int norm42 = gDiff4 * gDiff4;
                                int norm43 = bDiff4 * bDiff4;
                                int norm44 = aDiff4 * aDiff4;

                                int sum1 = norm11 + norm12;
                                int sum2 = norm13 + norm14;
                                int sum3 = norm21 + norm22;
                                int sum4 = norm23 + norm24;
                                int sum5 = norm31 + norm32;
                                int sum6 = norm33 + norm34;
                                int sum7 = norm41 + norm42;
                                int sum8 = norm43 + norm44;

                                blocks[blockIndex++].value += sum1 + sum2;
                                blocks[blockIndex++].value += sum3 + sum4;
                                blocks[blockIndex++].value += sum5 + sum6;
                                blocks[blockIndex++].value += sum7 + sum8;
                                srcRow += 16;
                            }
                            for (; srcX < endX; srcX++) {
                                int rSrc1 = srcRow[j4];
                                int gSrc1 = srcRow[j4 + 1];
                                int bSrc1 = srcRow[j4 + 2];
                                int aSrc1 = srcRow[j4 + 3];

                                int rDiff1 = rDst1 - rSrc1;
                                int gDiff1 = gDst1 - gSrc1;
                                int bDiff1 = bDst1 - bSrc1;
                                int aDiff1 = aDst1 - aSrc1;

                                int norm11 = rDiff1 * rDiff1;
                                int norm12 = gDiff1 * gDiff1;
                                int norm13 = bDiff1 * bDiff1;
                                int norm14 = aDiff1 * aDiff1;

                                int sum1 = norm11 + norm12;
                                int sum2 = norm13 + norm14;

                                blocks[blockIndex++].value += sum1 + sum2;
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
                    int j = 0;
                    for (; j < overlapWidth; j++) {
                        int j4 = j * 4;
                        int rDst1 = outputRow[j4];
                        int gDst1 = outputRow[j4 + 1];
                        int bDst1 = outputRow[j4 + 2];
                        int aDst1 = outputRow[j4 + 3];
                        for (int srcY = blockYStart; srcY < std::min(blockYStart + blockSize, maxBlockY);
                             srcY++) {
                            int srcYStart = srcY + overlapHeight;
                            unsigned char* srcRow = mData->data[srcYStart + i] + CHANNEL_NUM * blockXStart;
                            int srcX = blockXStart;
                            int endX = std::min(blockXStart + blockSize, maxBlockX);
                            int blockIndex = srcY * maxBlockX + srcX;
                            for (; srcX < endX-3; srcX+=4) {
                                int rSrc1 = srcRow[j4];
                                int gSrc1 = srcRow[j4 + 1];
                                int bSrc1 = srcRow[j4 + 2];
                                int aSrc1 = srcRow[j4 + 3];
                                int rSrc2 = srcRow[j4 + 4];
                                int gSrc2 = srcRow[j4 + 5];
                                int bSrc2 = srcRow[j4 + 6];
                                int aSrc2 = srcRow[j4 + 7];
                                int rSrc3 = srcRow[j4 + 8];
                                int gSrc3 = srcRow[j4 + 9];
                                int bSrc3 = srcRow[j4 + 10];
                                int aSrc3 = srcRow[j4 + 11];
                                int rSrc4 = srcRow[j4 + 12];
                                int gSrc4 = srcRow[j4 + 13];
                                int bSrc4 = srcRow[j4 + 14];
                                int aSrc4 = srcRow[j4 + 15];

                                int rDiff1 = rDst1 - rSrc1;
                                int gDiff1 = gDst1 - gSrc1;
                                int bDiff1 = bDst1 - bSrc1;
                                int aDiff1 = aDst1 - aSrc1;

                                int norm11 = rDiff1 * rDiff1;
                                int rDiff2 = rDst1 - rSrc2;
                                int gDiff2 = gDst1 - gSrc2;
                                int bDiff2 = bDst1 - bSrc2;

                                int norm12 = gDiff1 * gDiff1;
                                int aDiff2 = aDst1 - aSrc2;
                                int rDiff3 = rDst1 - rSrc3;
                                int gDiff3 = gDst1 - gSrc3;

                                int norm13 = bDiff1 * bDiff1;
                                int bDiff3 = bDst1 - bSrc3;
                                int aDiff3 = aDst1 - aSrc3;
                                int rDiff4 = rDst1 - rSrc4;

                                int norm14 = aDiff1 * aDiff1;
                                int gDiff4 = gDst1 - gSrc4;
                                int bDiff4 = bDst1 - bSrc4;
                                int aDiff4 = aDst1 - aSrc4;

                                int norm21 = rDiff2 * rDiff2;
                                int norm22 = gDiff2 * gDiff2;
                                int norm23 = bDiff2 * bDiff2;
                                int norm24 = aDiff2 * aDiff2;

                                int norm31 = rDiff3 * rDiff3;
                                int norm32 = gDiff3 * gDiff3;
                                int norm33 = bDiff3 * bDiff3;
                                int norm34 = aDiff3 * aDiff3;

                                int norm41 = rDiff4 * rDiff4;
                                int norm42 = gDiff4 * gDiff4;
                                int norm43 = bDiff4 * bDiff4;
                                int norm44 = aDiff4 * aDiff4;

                                int sum1 = norm11 + norm12;
                                int sum2 = norm13 + norm14;
                                int sum3 = norm21 + norm22;
                                int sum4 = norm23 + norm24;
                                int sum5 = norm31 + norm32;
                                int sum6 = norm33 + norm34;
                                int sum7 = norm41 + norm42;
                                int sum8 = norm43 + norm44;

                                blocks[blockIndex++].value += sum1 + sum2;
                                blocks[blockIndex++].value += sum3 + sum4;
                                blocks[blockIndex++].value += sum5 + sum6;
                                blocks[blockIndex++].value += sum7 + sum8;
                                srcRow += 16;
                            }
                            for (; srcX < endX; srcX++) {
                                int rSrc1 = srcRow[j * 4];
                                int gSrc1 = srcRow[j * 4 + 1];
                                int bSrc1 = srcRow[j * 4 + 2];
                                int aSrc1 = srcRow[j * 4 + 3];

                                int rDiff1 = rDst1 - rSrc1;
                                int gDiff1 = gDst1 - gSrc1;
                                int bDiff1 = bDst1 - bSrc1;
                                int aDiff1 = aDst1 - aSrc1;

                                int norm11 = rDiff1 * rDiff1;
                                int norm12 = gDiff1 * gDiff1;
                                int norm13 = bDiff1 * bDiff1;
                                int norm14 = aDiff1 * aDiff1;

                                int sum1 = norm11 + norm12;
                                int sum2 = norm13 + norm14;

                                blocks[blockIndex++].value += sum1 + sum2;
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
                        int j4 = j * 4;
                        int rDst1 = outputRow[j4];
                        int gDst1 = outputRow[j4 + 1];
                        int bDst1 = outputRow[j4 + 2];
                        int aDst1 = outputRow[j4 + 3];
                        for (int srcY = blockYStart; srcY < std::min(blockYStart + blockSize, maxBlockY);
                             srcY++) {
                            unsigned char* srcRow = mData->data[srcY + i] + CHANNEL_NUM * blockXStart;
                            int srcX = blockXStart;
                            int endX = std::min(blockXStart + blockSize, maxBlockX);
                            int blockIndex = srcY * maxBlockX + srcX;
                            for (; srcX < endX-3; srcX+=4) {
                                int rSrc1 = srcRow[j4];
                                int gSrc1 = srcRow[j4 + 1];
                                int bSrc1 = srcRow[j4 + 2];
                                int aSrc1 = srcRow[j4 + 3];
                                int rSrc2 = srcRow[j4 + 4];
                                int gSrc2 = srcRow[j4 + 5];
                                int bSrc2 = srcRow[j4 + 6];
                                int aSrc2 = srcRow[j4 + 7];
                                int rSrc3 = srcRow[j4 + 8];
                                int gSrc3 = srcRow[j4 + 9];
                                int bSrc3 = srcRow[j4 + 10];
                                int aSrc3 = srcRow[j4 + 11];
                                int rSrc4 = srcRow[j4 + 12];
                                int gSrc4 = srcRow[j4 + 13];
                                int bSrc4 = srcRow[j4 + 14];
                                int aSrc4 = srcRow[j4 + 15];

                                int rDiff1 = rDst1 - rSrc1;
                                int gDiff1 = gDst1 - gSrc1;
                                int bDiff1 = bDst1 - bSrc1;
                                int aDiff1 = aDst1 - aSrc1;

                                int norm11 = rDiff1 * rDiff1;
                                int rDiff2 = rDst1 - rSrc2;
                                int gDiff2 = gDst1 - gSrc2;
                                int bDiff2 = bDst1 - bSrc2;

                                int norm12 = gDiff1 * gDiff1;
                                int aDiff2 = aDst1 - aSrc2;
                                int rDiff3 = rDst1 - rSrc3;
                                int gDiff3 = gDst1 - gSrc3;

                                int norm13 = bDiff1 * bDiff1;
                                int bDiff3 = bDst1 - bSrc3;
                                int aDiff3 = aDst1 - aSrc3;
                                int rDiff4 = rDst1 - rSrc4;

                                int norm14 = aDiff1 * aDiff1;
                                int gDiff4 = gDst1 - gSrc4;
                                int bDiff4 = bDst1 - bSrc4;
                                int aDiff4 = aDst1 - aSrc4;

                                int norm21 = rDiff2 * rDiff2;
                                int norm22 = gDiff2 * gDiff2;
                                int norm23 = bDiff2 * bDiff2;
                                int norm24 = aDiff2 * aDiff2;

                                int norm31 = rDiff3 * rDiff3;
                                int norm32 = gDiff3 * gDiff3;
                                int norm33 = bDiff3 * bDiff3;
                                int norm34 = aDiff3 * aDiff3;

                                int norm41 = rDiff4 * rDiff4;
                                int norm42 = gDiff4 * gDiff4;
                                int norm43 = bDiff4 * bDiff4;
                                int norm44 = aDiff4 * aDiff4;

                                int sum1 = norm11 + norm12;
                                int sum2 = norm13 + norm14;
                                int sum3 = norm21 + norm22;
                                int sum4 = norm23 + norm24;
                                int sum5 = norm31 + norm32;
                                int sum6 = norm33 + norm34;
                                int sum7 = norm41 + norm42;
                                int sum8 = norm43 + norm44;

                                blocks[blockIndex++].value += sum1 + sum2;
                                blocks[blockIndex++].value += sum3 + sum4;
                                blocks[blockIndex++].value += sum5 + sum6;
                                blocks[blockIndex++].value += sum7 + sum8;
                                srcRow += 16;
                            }
                            for (; srcX < endX; srcX++) {
                                int rSrc1 = srcRow[j * 4];
                                int gSrc1 = srcRow[j * 4 + 1];
                                int bSrc1 = srcRow[j * 4 + 2];
                                int aSrc1 = srcRow[j * 4 + 3];

                                int rDiff1 = rDst1 - rSrc1;
                                int gDiff1 = gDst1 - gSrc1;
                                int bDiff1 = bDst1 - bSrc1;
                                int aDiff1 = aDst1 - aSrc1;

                                int norm11 = rDiff1 * rDiff1;
                                int norm12 = gDiff1 * gDiff1;
                                int norm13 = bDiff1 * bDiff1;
                                int norm14 = aDiff1 * aDiff1;

                                int sum1 = norm11 + norm12;
                                int sum2 = norm13 + norm14;

                                blocks[blockIndex++].value += sum1 + sum2;
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
                        int j4 = j * 4;
                        int rDst1 = outputRow[j4];
                        int gDst1 = outputRow[j4 + 1];
                        int bDst1 = outputRow[j4 + 2];
                        int aDst1 = outputRow[j4 + 3];

                        for (int srcY = blockYStart; srcY < std::min(blockYStart + blockSize, maxBlockY);
                             srcY++) {
                            unsigned char* srcRow = mData->data[srcY + i] + CHANNEL_NUM * blockXStart;
                            int srcX = blockXStart;
                            int endX = std::min(blockXStart + blockSize, maxBlockX);
                            int blockIndex = srcY * maxBlockX + srcX;
                            for (; srcX < endX-3; srcX+=4) {
                                int rSrc1 = srcRow[j4];
                                int gSrc1 = srcRow[j4 + 1];
                                int bSrc1 = srcRow[j4 + 2];
                                int aSrc1 = srcRow[j4 + 3];
                                int rSrc2 = srcRow[j4 + 4];
                                int gSrc2 = srcRow[j4 + 5];
                                int bSrc2 = srcRow[j4 + 6];
                                int aSrc2 = srcRow[j4 + 7];
                                int rSrc3 = srcRow[j4 + 8];
                                int gSrc3 = srcRow[j4 + 9];
                                int bSrc3 = srcRow[j4 + 10];
                                int aSrc3 = srcRow[j4 + 11];
                                int rSrc4 = srcRow[j4 + 12];
                                int gSrc4 = srcRow[j4 + 13];
                                int bSrc4 = srcRow[j4 + 14];
                                int aSrc4 = srcRow[j4 + 15];

                                int rDiff1 = rDst1 - rSrc1;
                                int gDiff1 = gDst1 - gSrc1;
                                int bDiff1 = bDst1 - bSrc1;
                                int aDiff1 = aDst1 - aSrc1;

                                int norm11 = rDiff1 * rDiff1;
                                int rDiff2 = rDst1 - rSrc2;
                                int gDiff2 = gDst1 - gSrc2;
                                int bDiff2 = bDst1 - bSrc2;

                                int norm12 = gDiff1 * gDiff1;
                                int aDiff2 = aDst1 - aSrc2;
                                int rDiff3 = rDst1 - rSrc3;
                                int gDiff3 = gDst1 - gSrc3;

                                int norm13 = bDiff1 * bDiff1;
                                int bDiff3 = bDst1 - bSrc3;
                                int aDiff3 = aDst1 - aSrc3;
                                int rDiff4 = rDst1 - rSrc4;

                                int norm14 = aDiff1 * aDiff1;
                                int gDiff4 = gDst1 - gSrc4;
                                int bDiff4 = bDst1 - bSrc4;
                                int aDiff4 = aDst1 - aSrc4;

                                int norm21 = rDiff2 * rDiff2;
                                int norm22 = gDiff2 * gDiff2;
                                int norm23 = bDiff2 * bDiff2;
                                int norm24 = aDiff2 * aDiff2;

                                int norm31 = rDiff3 * rDiff3;
                                int norm32 = gDiff3 * gDiff3;
                                int norm33 = bDiff3 * bDiff3;
                                int norm34 = aDiff3 * aDiff3;

                                int norm41 = rDiff4 * rDiff4;
                                int norm42 = gDiff4 * gDiff4;
                                int norm43 = bDiff4 * bDiff4;
                                int norm44 = aDiff4 * aDiff4;

                                int sum1 = norm11 + norm12;
                                int sum2 = norm13 + norm14;
                                int sum3 = norm21 + norm22;
                                int sum4 = norm23 + norm24;
                                int sum5 = norm31 + norm32;
                                int sum6 = norm33 + norm34;
                                int sum7 = norm41 + norm42;
                                int sum8 = norm43 + norm44;

                                blocks[blockIndex++].value += sum1 + sum2;
                                blocks[blockIndex++].value += sum3 + sum4;
                                blocks[blockIndex++].value += sum5 + sum6;
                                blocks[blockIndex++].value += sum7 + sum8;
                                srcRow += 16;
                            }
                            for (; srcX < endX; srcX++) {
                                int rSrc1 = srcRow[j * 4];
                                int gSrc1 = srcRow[j * 4 + 1];
                                int bSrc1 = srcRow[j * 4 + 2];
                                int aSrc1 = srcRow[j * 4 + 3];

                                int rDiff1 = rDst1 - rSrc1;
                                int gDiff1 = gDst1 - gSrc1;
                                int bDiff1 = bDst1 - bSrc1;
                                int aDiff1 = aDst1 - aSrc1;

                                int norm11 = rDiff1 * rDiff1;
                                int norm12 = gDiff1 * gDiff1;
                                int norm13 = bDiff1 * bDiff1;
                                int norm14 = aDiff1 * aDiff1;

                                int sum1 = norm11 + norm12;
                                int sum2 = norm13 + norm14;

                                blocks[blockIndex++].value += sum1 + sum2;
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

void AdvanceAlgOptimiz::OverlapConstraintsWithMinCut_StdC_KSrc4Unroll_BoundsRefactor_LoopReorder_Blocking32() {
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
        PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KSrc4Unroll_BoundsRefactor_LoopReorder_Blocking32(
                vertical, 0, dstX, maxBlockX, maxBlockY, 0.1, regBlockW, regBlockH);
    }
    int lastDstX = dstX - wStep;
    int blockWidth = mData->output_w - (lastDstX - overlapWidth);

    // fill all corner cases except borders
    for (; dstY < mData->output_h-hStep; dstY += hStep) {
        // fill first column
        PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KSrc4Unroll_BoundsRefactor_LoopReorder_Blocking32(
                horizontal, dstY, 0, maxBlockX, maxBlockY, 0.1, regBlockW, regBlockH);
        dstX = mData->block_w;
        for (; dstX < mData->output_w-wStep; dstX += wStep) {
            PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KSrc4Unroll_BoundsRefactor_LoopReorder_Blocking32(
                    both, dstY, dstX, maxBlockX, maxBlockY, 0.1, regBlockW, regBlockH);
        }

        // fill last column
        PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KSrc4Unroll_BoundsRefactor_LoopReorder_Blocking32(
                both, dstY, dstX, maxBlockX, maxBlockY, 0.1, blockWidth, regBlockH);
    }

    // fill last row
    int blockHeight = mData->output_h - dstY;
    PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KSrc4Unroll_BoundsRefactor_LoopReorder_Blocking32(
            horizontal, dstY, 0, maxBlockX, maxBlockY, 0.1, regBlockW, regBlockH);
    dstX = mData->block_w;
    for (; dstX < mData->output_w-wStep; dstX += wStep) {
        PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KSrc4Unroll_BoundsRefactor_LoopReorder_Blocking32(
                both, dstY, dstX, maxBlockX, maxBlockY, 0.1, mData->block_w, blockHeight);

    }
    // bottom-right corner
    PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KSrc4Unroll_BoundsRefactor_LoopReorder_Blocking32(
            both, dstY, lastDstX, maxBlockX, maxBlockY, 0.1, blockWidth, blockHeight);
}


