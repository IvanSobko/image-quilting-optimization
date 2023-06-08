#pragma once

#include "img_data.h"

class ImageQuilting {
public:
    ImageQuilting() = delete;
    ImageQuilting(ImgData* data) {
        mData = data;
        overlapHeight = mData->block_h / 6;
        overlapWidth = mData->block_w / 6;
    }

    // Synthesize a new texture with the given seed
    void Synthesis(int seed = -1);
    int64_t getFlopCount() const;

protected:
    // Enum representing the type of overlap between blocks
    enum OverlapType { vertical = 0, horizontal = 1, both = 2 };

    // Struct to sort blocks by their l2 norm
    struct BlockValue {
        int y, x;
        double value;
    };

    // Seed the random number generator with a specified seed
    static void SeedRandomNumberGenerator(int seed);
    // Generate a random number in the range [min, max]
    static int GetRandomInt(int min, int max);

    // Write a block from the source data to the output data given their upper-left corners
    void WriteBlock(int dstY, int dstX, int srcY, int srcX);
    // Same as WriteBlockOverlap, but uses a minimum cut to write the new block
    void WriteBlockOverlapWithMinCut(int overlapType, int dstY, int dstX, int srcY, int srcX);

    int64_t flopCount = 0;
    static constexpr double errorTolerance = 0.1;
    int overlapWidth = 0;
    int overlapHeight = 0;
    // Keep a pointer to the input image data
    ImgData* mData = nullptr;

private:
    // Same as WriteBlock but leaves the half of the dst overlapping region untouched
    void WriteBlockOverlap(int overlapType, int dstY, int dstX, int srcY, int srcX);

    // Compute the overlap between the current block, dst, of the output image
    // and src, of the input image, given their upper-left corners
    // and the position of the overlap
    double ComputeOverlap(int overlapType, int dstY, int dstX, int srcY, int srcX);

    // Place an edge overlap block with respect to the given block of the output image
    void PlaceEdgeOverlapBlock(int blockY, int blockX, int maxBlockX, int maxBlockY);

    // Place an edge overlap block with respect to the given block of the output image
    void PlaceEdgeOverlapBlockWithMinCut(int blockY, int blockX, int maxBlockX, int maxBlockY);
    // Method 1: Synthesize a new texture by randomly choosing blocks satisfying constraints and applying
    // minimum cuts
    void OverlapConstraintsWithMinCut();
    // Method 2: Synthesize a new texture sample by randomly choosing blocks satisfying overlap constraints
    void OverlapConstraints();
    // Method 3: Synthesize a new texture sample by randomly choosing blocks
    void RandomBlockPlacement();
};
