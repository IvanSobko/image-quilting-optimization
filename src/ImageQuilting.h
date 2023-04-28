#pragma once

#include "ImgData.h"

class ImageQuilting {
   public:
    ImageQuilting() = delete;
    ImageQuilting(const ImgData& data) { mData = data; }

    // Synthesize a new texture
    ImgData Synthesis();
    void Transfer() {}

   private:
    ImgData mData;

    // Write a block from the source data to the output data given their upper-left corners
    void WriteBlock(int dstY, int dstX, int srcY, int srcX);

    // Synthesize a new texture sample by randomly choosing blocks
    ImgData RandomBlockPlacement();

    // Choose what type of overlapping to estimate
    double ComputeEdgeOverlap(int block0Y, int block0X, int block1Y, int block1X);

    // Compute the vertical edge overlap between block 0 of the output image and block 1 of the input image given their upper-left corners
    double ComputeVerticalEdgeOverlap(int block0Y, int block0X, int block1Y, int block1X);

    // Compute the horizontal edge overlap between block 0 of the output image and block 1 of the input image given their upper-left corners
    double ComputeHorizontalEdgeOverlap(int block0Y, int block0X, int block1Y, int block1X);

    // Compute the corner edge overlap between block 0 of the output image and block 1 of the input image given their upper-left corners
    double ComputeCornerEdgeOverlap(int block0Y, int block0X, int block1Y, int block1X);

    // Struct to sort blocks by their l2 norm
    struct BlockValue{
        int y, x;
        double value;
    };

    // BlockValue comparator
    static int BlockValueComparator(const void* blockValue1, const void* blockValue2);

    // C style binary search with upper bound
    static size_t bs_upper_bound(
        const void* array, size_t n, const void* upper_bound,
        size_t width, int (*comparator)(const void*, const void*));

    // Place an edge overlap block with respect to the given block of the output image
    void PlaceEdgeOverlapBlock(int blockY, int blockX, int maxBlockX, int maxBlockY, double errorTolerance);

    // Synthesize a new texture sample by randomly choosing blocks satisfying overlap constraints
    ImgData OverlapConstraints();
};
