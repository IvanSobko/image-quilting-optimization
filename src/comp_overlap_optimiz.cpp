#include "comp_overlap_optimiz.h"
#include <algorithm>
#include <cfloat>
#include <cstdlib>
#include <random>

#include <immintrin.h>
// IMPORTANT: make sure it's the same for all functions that you compare
// Added this because we only need to increase the REP # for individual functions
#define IND_FUNC_REP 10000

// Image quilting function wrapper for testing and timing
double CompOverlapOptimiz::BasicOpt(ImgData* imgData, int seed) {
    CompOverlapOptimiz imageQuilting(imgData);
    imageQuilting.Synthesis(seed, opt_indices);
    return static_cast<double>(imageQuilting.getFlopCount());
}

double CompOverlapOptimiz::AlgOpt(ImgData* imgData, int seed) {
    CompOverlapOptimiz imageQuilting(imgData);
    imageQuilting.Synthesis(seed, opt_algorithm);
    return static_cast<double>(imageQuilting.getFlopCount());
}

double CompOverlapOptimiz::UnrollOpt(ImgData* imgData, int seed) {
    CompOverlapOptimiz imageQuilting(imgData);
    imageQuilting.Synthesis(seed, opt_unroll);
    return static_cast<double>(imageQuilting.getFlopCount());
}

double CompOverlapOptimiz::UnrollMaxOpt(ImgData* imgData, int seed) {
    CompOverlapOptimiz imageQuilting(imgData);
    imageQuilting.Synthesis(seed, opt_unroll_max);
    return static_cast<double>(imageQuilting.getFlopCount());
}

#ifdef __AVX2__
double CompOverlapOptimiz::VectorizeOpt(ImgData* imgData, int seed) {
    CompOverlapOptimiz imageQuilting(imgData);
    imageQuilting.Synthesis(seed, opt_vectorize);
    return static_cast<double>(imageQuilting.getFlopCount());
}
#endif

double CompOverlapOptimiz::UnrollChnls(ImgData* imgData, int seed) {
    CompOverlapOptimiz imageQuilting(imgData);
    imageQuilting.Synthesis(seed, opt_unroll_chnls);
    return static_cast<double>(imageQuilting.getFlopCount());
}

void CompOverlapOptimiz::GetComponentParameters(ImgData* imgData, int& overlapType, int& dstY, int& dstX,
                                                int& srcY, int& srcX) {
    int maxBlockYSrc = imgData->height - imgData->block_h;
    int maxBlockXSrc = imgData->width - imgData->block_w;
    int maxBlockYDst = imgData->output_h - imgData->block_h;
    int maxBlockXDst = imgData->output_w - imgData->block_w;
    srcY = GetRandomInt(0, maxBlockYSrc - 1);
    srcX = GetRandomInt(0, maxBlockXSrc - 1);
    dstY = GetRandomInt(imgData->block_h, maxBlockYDst - 1);
    dstX = GetRandomInt(imgData->block_w, maxBlockXDst - 1);
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

volatile void CompOverlapOptimiz::UnrollMaxOptComponent(ImgData* imgData, int seed) {
    CompOverlapOptimiz imageQuilting(imgData);
    int overlapType, dstY, dstX, srcY, srcX;
    volatile double dummy = 0;
    for (int k = 0; k < IND_FUNC_REP; k++) {
        GetComponentParameters(imgData, overlapType, dstY, dstX, srcY, srcX);
        dummy += imageQuilting.ComputeOverlapUnrollMax(overlapType, dstY, dstX, srcY, srcX);
    }
}

#ifdef __AVX2__
volatile void CompOverlapOptimiz::VectorizeOptComponent(ImgData* imgData, int seed) {
    CompOverlapOptimiz imageQuilting(imgData);
    int overlapType, dstY, dstX, srcY, srcX;
    volatile double dummy = 0;
    for (int k = 0; k < IND_FUNC_REP; k++) {
        GetComponentParameters(imgData, overlapType, dstY, dstX, srcY, srcX);
        dummy += imageQuilting.ComputeOverlapVectorize(overlapType, dstY, dstX, srcY, srcX);
    }
}
#endif

// Synthesize a new texture with the given seed
void CompOverlapOptimiz::Synthesis(int seed, int opt) {
    flopCount = 0;
    opt_type = opt;
    SeedRandomNumberGenerator(seed);
    OverlapConstraintsWithMinCut();
}

// Compute the overlap between the current block - block 0 of the output image
// and block 1 of the input image given their upper-left corners
// and the position of the overlap
double CompOverlapOptimiz::ComputeOverlapBasicOpt(const int overlapType, const int dstY, const int dstX,
                                                  const int srcY, const int srcX) {
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
        int srcYStart = srcY + srcYOffset;
        int dstXStart = CHANNEL_NUM * overlapXStart;
        int srcXStart = CHANNEL_NUM * srcX;
        for (int i = 0; i < verticalBlockHeightLocal; i++) {
            unsigned char* outputRow = mData->output_d[dstY + i] + dstXStart;
            unsigned char* srcRow = mData->data[srcYStart + i] + srcXStart;
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

    // Compute the horizontal overlap
    if (overlapType == horizontal || overlapType == both) {
        int srcXOffset = overlapType == both ? overlapWidth : 0;
        int dstXStart = CHANNEL_NUM * dstX;
        int srcXStart = CHANNEL_NUM * (srcX + srcXOffset);
        for (int i = 0; i < overlapHeight; i++) {
            unsigned char* outputRow = mData->output_d[overlapYStart + i] + dstXStart;
            unsigned char* srcRow = mData->data[srcY + i] + srcXStart;
            for (int j = 0; j < horizontalBlockWidthLocal; j++) {
                for (int k = 0; k < CHANNEL_NUM; k++) {
                    int x0 = *(outputRow++);
                    int x1 = *(srcRow++);
                    int norm = x0 - x1;
                    l2norm += norm * norm;
                }
            }
        }
    }

    // Compute the corner edge overlap
    if (overlapType == both) {
        int dstXStart = CHANNEL_NUM * overlapXStart;
        int srcXStart = CHANNEL_NUM * srcX;
        for (int i = 0; i < overlapHeight; i++) {
            unsigned char* outputRow = mData->output_d[overlapYStart + i] + dstXStart;
            unsigned char* srcRow = mData->data[srcY + i] + srcXStart;
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
        for (int i = 0; i < overlapHeight; i++) {
            unsigned char* outputRow = mData->output_d[overlapYStart + i] + dstXStart;
            unsigned char* srcRow = mData->data[srcY + i] + srcXStart;
            for (int j = 0; j < horizontalBlockWidthLocal; j++) {
                for (int k = 0; k < CHANNEL_NUM; k++) {
                    int x0 = *(outputRow++);
                    int x1 = *(srcRow++);
                    int norm = x0 - x1;
                    l2norm += norm * norm;
                }
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

double CompOverlapOptimiz::ComputeOverlapUnrollChannels(int overlapType, int dstY, int dstX, int srcY,
                                                        int srcX) {
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
        for (int i = 0; i < overlapHeight; i++) {
            unsigned char* outputRow = mData->output_d[overlapYStart + i] + dstXStart;
            unsigned char* srcRow = mData->data[srcY + i] + srcXStart;
            for (int j = 0; j < horizontalBlockWidthLocal; j++) {
                int rDst = outputRow[j * 4];
                int gDst = outputRow[j * 4 + 1];
                int bDst = outputRow[j * 4 + 2];
                int aDst = outputRow[j * 4 + 3];
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

                l2norm += sum1 + sum2;
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
            for (int j = 0; j < overlapWidth; j++) {
                int rDst = outputRow[j * 4];
                int gDst = outputRow[j * 4 + 1];
                int bDst = outputRow[j * 4 + 2];
                int aDst = outputRow[j * 4 + 3];
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

                l2norm += sum1 + sum2;
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
        for (int i = 0; i < overlapHeight; i++) {
            unsigned char* outputRow = mData->output_d[overlapYStart + i] + dstXStart;
            unsigned char* srcRow = mData->data[srcY + i] + srcXStart;
            int j = 0;
            int sum1 = 0;
            int sum2 = 0;
            int sum3 = 0;
            int sum4 = 0;
            for (j = 0; j < horizontalBlockWidthLocal - 1; j += 2) {
                int pos = j * 4;
                int rDst1 = outputRow[pos];
                int gDst1 = outputRow[pos + 1];
                int bDst1 = outputRow[pos + 2];
                int aDst1 = outputRow[pos + 3];
                int rSrc1 = srcRow[pos];
                int gSrc1 = srcRow[pos + 1];
                int bSrc1 = srcRow[pos + 2];
                int aSrc1 = srcRow[pos + 3];

                int rDst2 = outputRow[pos + 4];
                int gDst2 = outputRow[pos + 5];
                int bDst2 = outputRow[pos + 6];
                int aDst2 = outputRow[pos + 7];
                int rSrc2 = srcRow[pos + 4];
                int gSrc2 = srcRow[pos + 5];
                int bSrc2 = srcRow[pos + 6];
                int aSrc2 = srcRow[pos + 7];

                int rDiff1 = rDst1 - rSrc1;
                int gDiff1 = gDst1 - gSrc1;
                int bDiff1 = bDst1 - bSrc1;
                int aDiff1 = aDst1 - aSrc1;

                int rNorm1 = rDiff1 * rDiff1;
                int rDiff2 = rDst2 - rSrc2;
                int gDiff2 = gDst2 - gSrc2;
                int bDiff2 = bDst2 - bSrc2;

                int gNorm1 = gDiff1 * gDiff1;
                int aDiff2 = aDst2 - aSrc2;
                int bNorm1 = bDiff1 * bDiff1;
                int aNorm1 = aDiff1 * aDiff1;
                int rNorm2 = rDiff2 * rDiff2;
                sum1 += rNorm1 + gNorm1;

                int gNorm2 = gDiff2 * gDiff2;
                int bNorm2 = bDiff2 * bDiff2;
                sum2 += bNorm1 + aNorm1;
                int aNorm2 = aDiff2 * aDiff2;

                sum3 += rNorm2 + gNorm2;
                sum4 += bNorm2 + aNorm2;
            }
            l2norm += sum1 + sum2 + sum3 + sum4;
            for (; j < horizontalBlockWidthLocal; j++) {
                int rDst = outputRow[j * 4];
                int gDst = outputRow[j * 4 + 1];
                int bDst = outputRow[j * 4 + 2];
                int aDst = outputRow[j * 4 + 3];
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

                l2norm += sum1 + sum2;
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
            int j = 0;
            int sum1 = 0;
            int sum2 = 0;
            int sum3 = 0;
            int sum4 = 0;
            for (j = 0; j < overlapWidth - 1; j += 2) {
                int pos = j * 4;
                int rDst1 = outputRow[pos];
                int gDst1 = outputRow[pos + 1];
                int bDst1 = outputRow[pos + 2];
                int aDst1 = outputRow[pos + 3];
                int rSrc1 = srcRow[pos];
                int gSrc1 = srcRow[pos + 1];
                int bSrc1 = srcRow[pos + 2];
                int aSrc1 = srcRow[pos + 3];

                int rDst2 = outputRow[pos + 4];
                int gDst2 = outputRow[pos + 5];
                int bDst2 = outputRow[pos + 6];
                int aDst2 = outputRow[pos + 7];
                int rSrc2 = srcRow[pos + 4];
                int gSrc2 = srcRow[pos + 5];
                int bSrc2 = srcRow[pos + 6];
                int aSrc2 = srcRow[pos + 7];

                int rDiff1 = rDst1 - rSrc1;
                int gDiff1 = gDst1 - gSrc1;
                int bDiff1 = bDst1 - bSrc1;
                int aDiff1 = aDst1 - aSrc1;

                int rNorm1 = rDiff1 * rDiff1;
                int rDiff2 = rDst2 - rSrc2;
                int gDiff2 = gDst2 - gSrc2;
                int bDiff2 = bDst2 - bSrc2;

                int gNorm1 = gDiff1 * gDiff1;
                int aDiff2 = aDst2 - aSrc2;
                int bNorm1 = bDiff1 * bDiff1;
                int aNorm1 = aDiff1 * aDiff1;
                int rNorm2 = rDiff2 * rDiff2;
                sum1 += rNorm1 + gNorm1;

                int gNorm2 = gDiff2 * gDiff2;
                int bNorm2 = bDiff2 * bDiff2;
                sum2 += bNorm1 + aNorm1;
                int aNorm2 = aDiff2 * aDiff2;

                sum3 += rNorm2 + gNorm2;
                sum4 += bNorm2 + aNorm2;
            }
            l2norm += sum1 + sum2 + sum3 + sum4;
            for (; j < overlapWidth; j++) {
                int rDst = outputRow[j * 4];
                int gDst = outputRow[j * 4 + 1];
                int bDst = outputRow[j * 4 + 2];
                int aDst = outputRow[j * 4 + 3];
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

                l2norm += sum1 + sum2;
            }
        }
    }

    return std::sqrt(l2norm);
}

// Forced compiler to parallelize computations
double CompOverlapOptimiz::ComputeOverlapUnrollMax(int overlapType, int dstY, int dstX, int srcY, int srcX) {
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
        for (int i = 0; i < overlapHeight; i++) {
            unsigned char* outputRow = mData->output_d[overlapYStart + i] + dstXStart;
            unsigned char* srcRow = mData->data[srcY + i] + srcXStart;
            int norm0 = 0;
            int norm1 = 0;
            int norm2 = 0;
            int norm3 = 0;
            int j;
            for (j = 0; j < horizontalBlockWidthLocal - 3; j += 4) {
                int rDiff1 = outputRow[0] - srcRow[0];
                int gDiff1 = outputRow[1] - srcRow[1];
                int bDiff1 = outputRow[2] - srcRow[2];
                int aDiff1 = outputRow[3] - srcRow[3];

                rDiff1 = rDiff1 * rDiff1;
                int rDiff2 = outputRow[4] - srcRow[4];
                int gDiff2 = outputRow[5] - srcRow[5];
                int bDiff2 = outputRow[6] - srcRow[6];

                gDiff1 = gDiff1 * gDiff1;
                int aDiff2 = outputRow[7] - srcRow[7];
                int rDiff3 = outputRow[8] - srcRow[8];
                int gDiff3 = outputRow[9] - srcRow[9];

                bDiff1 = bDiff1 * bDiff1;
                int bDiff3 = outputRow[10] - srcRow[10];
                int aDiff3 = outputRow[11] - srcRow[11];
                int rDiff4 = outputRow[12] - srcRow[12];

                aDiff1 = aDiff1 * aDiff1;
                int gDiff4 = outputRow[13] - srcRow[13];
                int bDiff4 = outputRow[14] - srcRow[14];
                int aDiff4 = outputRow[15] - srcRow[15];

                rDiff2 = rDiff2 * rDiff2;
                gDiff2 = gDiff2 * gDiff2;
                bDiff2 = bDiff2 * bDiff2;
                aDiff2 = aDiff2 * aDiff2;
                l2norm += rDiff1 + gDiff1 + bDiff1 + aDiff1;

                rDiff3 = rDiff3 * rDiff3;
                gDiff3 = gDiff3 * gDiff3;
                bDiff3 = bDiff3 * bDiff3;
                aDiff3 = aDiff3 * aDiff3;
                l2norm += rDiff2 + gDiff2 + bDiff2 + aDiff2;

                rDiff4 = rDiff4 * rDiff4;
                gDiff4 = gDiff4 * gDiff4;
                bDiff4 = bDiff4 * bDiff4;
                aDiff4 = aDiff4 * aDiff4;

                l2norm += rDiff3 + gDiff3 + bDiff3 + aDiff3;
                l2norm += rDiff4 + gDiff4 + bDiff4 + aDiff4;

                outputRow += 16;
                srcRow += 16;
            }
            for (; j < horizontalBlockWidthLocal; j++) {
                int rDst = *(outputRow++);
                int rSrc = *(srcRow++);
                int gDst = *(outputRow++);
                int gSrc = *(srcRow++);
                int bDst = *(outputRow++);
                int bSrc = *(srcRow++);
                int aDst = *(outputRow++);
                int aSrc = *(srcRow++);

                int rDiff = rDst - rSrc;
                int gDiff = gDst - gSrc;
                int bDiff = bDst - bSrc;
                int aDiff = aDst - aSrc;

                l2norm += rDiff * rDiff + gDiff * gDiff + bDiff * bDiff + aDiff * aDiff;
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
            int j;
            for (j = 0; j < overlapWidth - 3; j += 4) {
                int rDiff1 = outputRow[0] - srcRow[0];
                int gDiff1 = outputRow[1] - srcRow[1];
                int bDiff1 = outputRow[2] - srcRow[2];
                int aDiff1 = outputRow[3] - srcRow[3];

                rDiff1 = rDiff1 * rDiff1;
                int rDiff2 = outputRow[4] - srcRow[4];
                int gDiff2 = outputRow[5] - srcRow[5];
                int bDiff2 = outputRow[6] - srcRow[6];

                gDiff1 = gDiff1 * gDiff1;
                int aDiff2 = outputRow[7] - srcRow[7];
                int rDiff3 = outputRow[8] - srcRow[8];
                int gDiff3 = outputRow[9] - srcRow[9];

                bDiff1 = bDiff1 * bDiff1;
                int bDiff3 = outputRow[10] - srcRow[10];
                int aDiff3 = outputRow[11] - srcRow[11];
                int rDiff4 = outputRow[12] - srcRow[12];

                aDiff1 = aDiff1 * aDiff1;
                int gDiff4 = outputRow[13] - srcRow[13];
                int bDiff4 = outputRow[14] - srcRow[14];
                int aDiff4 = outputRow[15] - srcRow[15];

                rDiff2 = rDiff2 * rDiff2;
                gDiff2 = gDiff2 * gDiff2;
                bDiff2 = bDiff2 * bDiff2;
                aDiff2 = aDiff2 * aDiff2;
                l2norm += rDiff1 + gDiff1 + bDiff1 + aDiff1;

                rDiff3 = rDiff3 * rDiff3;
                gDiff3 = gDiff3 * gDiff3;
                bDiff3 = bDiff3 * bDiff3;
                aDiff3 = aDiff3 * aDiff3;
                l2norm += rDiff2 + gDiff2 + bDiff2 + aDiff2;

                rDiff4 = rDiff4 * rDiff4;
                gDiff4 = gDiff4 * gDiff4;
                bDiff4 = bDiff4 * bDiff4;
                aDiff4 = aDiff4 * aDiff4;

                l2norm += rDiff3 + gDiff3 + bDiff3 + aDiff3;
                l2norm += rDiff4 + gDiff4 + bDiff4 + aDiff4;

                outputRow += 16;
                srcRow += 16;
            }
            for (; j < overlapWidth; j++) {
                int rDst = *(outputRow++);
                int rSrc = *(srcRow++);
                int gDst = *(outputRow++);
                int gSrc = *(srcRow++);
                int bDst = *(outputRow++);
                int bSrc = *(srcRow++);
                int aDst = *(outputRow++);
                int aSrc = *(srcRow++);

                int rDiff = rDst - rSrc;
                int gDiff = gDst - gSrc;
                int bDiff = bDst - bSrc;
                int aDiff = aDst - aSrc;

                l2norm += rDiff * rDiff + gDiff * gDiff + bDiff * bDiff + aDiff * aDiff;
            }
        }
    }

    return std::sqrt(l2norm);
}

#ifdef __AVX2__
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

    int range = horizontalBlockWidthLocal - (horizontalBlockWidthLocal % 4);
    // Compute the horizontal overlap (+corner if needed)
    if (overlapType != vertical) {
        int dstXStart = CHANNEL_NUM * overlapXStart;
        int srcXStart = CHANNEL_NUM * srcX;
        for (int i = 0; i < overlapHeight; i++) {
            unsigned char* outputRow = mData->output_d[overlapYStart + i] + dstXStart;
            unsigned char* srcRow = mData->data[srcY + i] + srcXStart;
            int j;
            for (j = 0; j < range; j += 4) {
                // load 16 8-bit integers(chars) and convert it to 16 16-bit integers
                // cycle 1
                __m128i dst = _mm_load_si128((__m128i*)(outputRow + j * 4));  // lat:6, tp: 0.5
                __m128i src = _mm_load_si128((__m128i*)(srcRow + j * 4));     // lat:6, tp: 0.5
                //cycle 6
                __m256i dst_16bit = _mm256_cvtepu8_epi16(dst);  // lat:3, tp: 1
                __m256i src_16bit = _mm256_cvtepu8_epi16(src);  // lat:3, tp: 1

                //cycle 10
                __m256i diff = _mm256_sub_epi16(dst_16bit, src_16bit);  //lat: 1, tp:0.3

                //cycle 11
                __m256i norm = _mm256_madd_epi16(diff, diff);  //lat 5, tp: 0.5

                // Start calculation of the sum of all norm elements. Overflow will wrap around
                // cycle 16
                __m128i norm_low = _mm256_extracti128_si256(norm, 1);  //lat 3, tp:1
                __m128i norm_high = _mm256_castsi256_si128(norm);      //lat 0

                //cycle 19
                __m128i norm_sum1 = _mm_add_epi32(norm_low, norm_high);  //lat 1, tp: 0.3
                //cycle 20
                __m128i norm_sum2 = _mm_add_epi32(
                    norm_sum1, _mm_unpackhi_epi64(norm_sum1, norm_sum1));  //lat:1 + 1, tp: 1 + 0.3
                //cycle 22
                __m128i norm_sum3 =
                    _mm_add_epi32(norm_sum2, _mm_shuffle_epi32(norm_sum2, 1));  //lat:1 + 1, tp: 1 + 0.3

                l2norm += _mm_cvtsi128_si32(norm_sum3);
            }

            for (; j < horizontalBlockWidthLocal; j++) {
                for (int k = 0; k < CHANNEL_NUM; k++) {
                    int x0 = *(outputRow + 4 * j + k);
                    int x1 = *(srcRow + 4 * j + k);
                    int norm = x0 - x1;
                    l2norm += norm * norm;
                }
            }
        }
    }

    range = overlapWidth - (overlapWidth % 4);
    // Compute the vertical overlap
    if (overlapType != horizontal) {
        int srcYOffset = overlapType == both ? overlapHeight : 0;
        int srcYStart = srcY + srcYOffset;
        int dstXStart = CHANNEL_NUM * overlapXStart;
        int srcXStart = CHANNEL_NUM * srcX;
        for (int i = 0; i < verticalBlockHeightLocal; i++) {
            unsigned char* outputRow = mData->output_d[dstY + i] + dstXStart;
            unsigned char* srcRow = mData->data[srcYStart + i] + srcXStart;
            int j;
            for (j = 0; j < range; j += 4) {
                // load 16 8-bit integers(chars) and convert it to 16 16-bit integers
                //cycle 1
                __m128i dst = _mm_load_si128((__m128i*)(outputRow + j * 4));  // lat:6, tp: 0.5
                __m128i src = _mm_load_si128((__m128i*)(srcRow + j * 4));     // lat:6, tp: 0.5
                //cycle 6
                __m256i dst_16bit = _mm256_cvtepu8_epi16(dst);  // lat:3, tp: 1
                __m256i src_16bit = _mm256_cvtepu8_epi16(src);  // lat:3, tp: 1
                //cycle 10
                __m256i diff = _mm256_sub_epi16(dst_16bit, src_16bit);  //lat: 1, tp:0.3
                //cycle 11
                __m256i norm = _mm256_madd_epi16(diff, diff);  //lat 5, tp: 0.5

                // Start calculation of the sum of all norm elements. Overflow will wrap around
                // cycle 16
                __m128i norm_low = _mm256_extracti128_si256(norm, 1);  //lat 3, tp:1
                __m128i norm_high = _mm256_castsi256_si128(norm);      //lat 0

                //cycle 19
                __m128i norm_sum1 = _mm_add_epi32(norm_low, norm_high);  //lat 1, tp: 0.3
                //cycle 20
                __m128i norm_sum2 = _mm_add_epi32(
                    norm_sum1, _mm_unpackhi_epi64(norm_sum1, norm_sum1));  //lat:1 + 1, tp: 1 + 0.3
                //cycle 22
                __m128i norm_sum3 =
                    _mm_add_epi32(norm_sum2, _mm_shuffle_epi32(norm_sum2, 1));  //lat:1 + 1, tp: 1 + 0.3

                int norm_sum = _mm_cvtsi128_si32(norm_sum3);

                l2norm += _mm_cvtsi128_si32(norm_sum3);  // + norm2_sum;
            }

            for (; j < overlapWidth; j++) {
                for (int k = 0; k < CHANNEL_NUM; k++) {
                    int x0 = *(outputRow + 4 * j + k);
                    int x1 = *(srcRow + 4 * j + k);
                    int norm = x0 - x1;
                    l2norm += norm * norm;
                }
            }
        }
    }

    return std::sqrt(l2norm);
}
#endif
// Base implementation of ComputeOverlap
double CompOverlapOptimiz::ComputeOverlapBase(const int overlapType, const int dstY, const int dstX,
                                              const int srcY, const int srcX) {
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
        for (int i = 0; i < verticalBlockHeightLocal; i++) {
            for (int j = 0; j < overlapWidth; j++) {
                for (int k = 0; k < CHANNEL_NUM; k++) {
                    double x0 = mData->output_d[dstY + i][CHANNEL_NUM * (overlapXStart + j) + k];
                    double x1 = mData->data[srcY + srcYOffset + i][CHANNEL_NUM * (srcX + j) + k];
                    double norm = x0 - x1;
                    l2norm += norm * norm;
                }
            }
        }
    }

    // Compute the horizontal overlap
    if (overlapType == horizontal || overlapType == both) {
        int srcXOffset = overlapType == both ? overlapWidth : 0;
        for (int i = 0; i < overlapHeight; i++) {
            for (int j = 0; j < horizontalBlockWidthLocal; j++) {
                for (int k = 0; k < CHANNEL_NUM; k++) {
                    double x0 = mData->output_d[overlapYStart + i][CHANNEL_NUM * (dstX + j) + k];
                    double x1 = mData->data[srcY + i][CHANNEL_NUM * (srcX + srcXOffset + j) + k];
                    double norm = x0 - x1;
                    l2norm += norm * norm;
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
void CompOverlapOptimiz::PlaceEdgeOverlapBlockWithMinCut(const int blockY, const int blockX,
                                                         const int maxBlockX, const int maxBlockY,
                                                         double errorTolerance) {
    // Calculate the overlap type
    int overlapType;
    if (blockY == 0) {
        overlapType = vertical;
    } else if (blockX == 0) {
        overlapType = horizontal;
    } else {
        overlapType = both;
    }

    // Compute the value of each block
    int numBlocks = maxBlockY * maxBlockX;
    BlockValue* blocks = (BlockValue*)malloc(sizeof(BlockValue) * numBlocks);
    for (int i = 0; i < maxBlockY; i++) {
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
            } else if (opt_type == opt_unroll_max) {
                blocks[blockIndex].value = ComputeOverlapUnrollMax(overlapType, blockY, blockX, i, j);
            } else if (opt_type == opt_vectorize) {
#ifdef __AVX2__
                blocks[blockIndex].value = ComputeOverlapVectorize(overlapType, blockY, blockX, i, j);
#endif
            } else if (opt_type == opt_unroll_chnls) {
                blocks[blockIndex].value = ComputeOverlapUnrollChannels(overlapType, blockY, blockX, i, j);
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
    WriteBlockOverlapWithMinCut(overlapType, blockY, blockX, suitableBlocks[blockIndex].y,
                                suitableBlocks[blockIndex].x);

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
        //NOTE: it's not correct for opt_indices, since there we have 3 loops for overlap 'both'
        flopCount += numBlocks * ((3 * CHANNEL_NUM * overlapWidth * mData->block_h) +
                                  (3 * CHANNEL_NUM * overlapHeight * mData->block_w) + 1);
    }
    // flops for intermediate calculations
    flopCount += 2 * numBlocks + 2;

    // flops for WriteBlockOverlapWithMinCut + some flops are computed in code
    // Note: approximating overlapHeightLocal and overlapWidthLocal as block_h and block_w
    if (overlapType == vertical) {
        flopCount += 3 * CHANNEL_NUM * mData->block_w * mData->block_h +
                     3 * mData->block_w * (mData->block_h - 1) + (mData->block_w - 1);
    } else if (overlapType == horizontal) {
        flopCount += 3 * CHANNEL_NUM * mData->block_w * mData->block_h +
                     3 * mData->block_h * (mData->block_w - 1) + (mData->block_h - 1);
    } else {
        flopCount += (3 * CHANNEL_NUM * mData->block_w * mData->block_h +
                      3 * mData->block_w * (mData->block_h - 1) + (mData->block_w - 1)) +
                     (3 * CHANNEL_NUM * mData->block_w * mData->block_h +
                      3 * mData->block_h * (mData->block_w - 1) + (mData->block_h - 1));
    }
}

// Synthesize a new texture by randomly choosing blocks satisfying constraints and applying minimum cuts
void CompOverlapOptimiz::OverlapConstraintsWithMinCut() {
    int hStep = mData->block_h - overlapHeight;
    int wStep = mData->block_w - overlapWidth;

    // The first block is full size; the others are of size step due to overlapping
    int numBlocksY = (mData->output_h - mData->block_h) / hStep + 2;
    int numBlocksX = (mData->output_w - mData->block_w) / wStep + 2;
    int maxBlockY = mData->height - mData->block_h;
    int maxBlockX = mData->width - mData->block_w;

    // Iterate over the block upper-left corners
    for (int blockY = 0; blockY < numBlocksY; blockY++) {
        for (int blockX = 0; blockX < numBlocksX; blockX++) {

            // Top-left corner of the current block
            int dstY = blockY == 0 ? 0 : mData->block_h + hStep * (blockY - 1);
            int dstX = blockX == 0 ? 0 : mData->block_w + wStep * (blockX - 1);

            // Make sure we are inside of the output image
            if (dstY > mData->output_h || dstX > mData->output_w)
                continue;

            // Randomly choose a block and place it
            if (blockY == 0 && blockX == 0) {
                // Randomly choose the upper-left corner of a block
                int srcY = GetRandomInt(0, maxBlockY - 1);
                int srcX = GetRandomInt(0, maxBlockX - 1);

                // Write the randomly chosen block to the output
                WriteBlock(dstY, dstX, srcY, srcX);
            } else {
                PlaceEdgeOverlapBlockWithMinCut(dstY, dstX, maxBlockX, maxBlockY, 0.1);
            }
        }
    }
}
