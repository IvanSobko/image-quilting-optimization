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

    // Same as WriteBlock but leaves the half of the dst overlapping region untouched
    void WriteBlockOverlap(int dstY, int dstX, int srcY, int srcX);

    // Same as WriteBlockOverlap, but uses a minimum cut to write the new block
    void WriteBlockOverlapWithMinCut(int overlapType, int dstY, int dstX, int srcY, int srcX);

    // Compute the overlap between the current block - block 0 of the output image
    // and block 1 of the input image given their upper-left corners
    // and the position of the overlap
    double ComputeOverlap(const int overlapType,
                          const int block0Y, const int block0X,
                          const int block1Y, const int block1X);

    // Struct to sort blocks by their l2 norm
    struct BlockValue{
        int y, x;
        double value;
    };

    // Enum representing the type of overlap between blocks
    enum OverlapType {
        vertical = 0,
        horizontal = 1,
        both = 2
    };

    // Place an edge overlap block with respect to the given block of the output image
    void PlaceEdgeOverlapBlock(int blockY, int blockX, int maxBlockX, int maxBlockY, double errorTolerance);

    // Place an edge overlap block with respect to the given block of the output image
    void PlaceEdgeOverlapBlockWithMinCut(int blockY, int blockX, int maxBlockX, int maxBlockY, double errorTolerance);

    // Synthesize a new texture by randomly choosing blocks satisfying constraints and applying minimum cuts
    ImgData OverlapConstraintsWithMinCut();

    // Synthesize a new texture sample by randomly choosing blocks satisfying overlap constraints
    ImgData OverlapConstraints();

    // Synthesize a new texture sample by randomly choosing blocks
    ImgData RandomBlockPlacement();

    int overlapWidth = 0;
    int overlapHeight = 0;
};
