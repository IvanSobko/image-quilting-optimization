#pragma once

#include "ImgData.h"

class AdvanceAlgOptimiz {
public:
    // Algorithm test functions
    static void DividedFuncOpt(ImgData* imgData, int seed);

    AdvanceAlgOptimiz() = delete;
    AdvanceAlgOptimiz(ImgData* data) {
        mData = data;
        overlapHeight = mData->block_h / 6;
        overlapWidth = mData->block_w / 6;
    }

    // Synthesize a new texture
    void Synthesis();
    // Synthesize a new texture with the given seed
    void Synthesis(int seed);

    // Seed the random number generator with the system time
    static void SeedRandomNumberGenerator();
    // Seed the random number generator with a specified seed
    static void SeedRandomNumberGenerator(int seed);
    // Generate a random number in the range [min, max]
    static int GetRandomInt(int min, int max);

    int64_t getFlopCount() const;

private:
    // Keep a pointer to the input image data
    ImgData* mData;

    // Write a block from the source data to the output data given their upper-left corners
    void WriteBlock(int dstY, int dstX, int srcY, int srcX);

    // Same as WriteBlockOverlap, but uses a minimum cut to write the new block
    void WriteBlockOverlapWithMinCut(int overlapType, int dstY, int dstX, int srcY, int srcX);

    // Struct to sort blocks by their l2 norm
    struct BlockValue {
        int y, x;
        double value;
    };

    // Enum representing the type of overlap between blocks
    enum OverlapType { vertical = 0, horizontal = 1, both = 2 };

    // Place an edge overlap block with respect to the given block of the output image
    void PlaceEdgeOverlapBlockWithMinCut(int overlapType, int dstY, int dstX, int maxBlockX, int maxBlockY,
                                         double errorTolerance, int bWidth, int bHeight);

    // Synthesize a new texture by randomly choosing blocks satisfying constraints and applying minimum cuts
    void OverlapConstraintsWithMinCut();

    int overlapWidth = 0;
    int overlapHeight = 0;

    int64_t flopCount = 0;
    int opt_type = 0;
};
