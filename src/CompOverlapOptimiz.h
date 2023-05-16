#pragma once

#include "ImgData.h"

class CompOverlapOptimiz {
public:
    static void BasicOpt(ImgData* imgData, int seed);

    CompOverlapOptimiz() = delete;
    CompOverlapOptimiz(ImgData* data) { mData = data; }

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

    // Same as WriteBlock but leaves the half of the dst overlapping region untouched
    void WriteBlockOverlap(int overlapType, int dstY, int dstX, int srcY, int srcX);

    // Same as WriteBlockOverlap, but uses a minimum cut to write the new block
    void WriteBlockOverlapWithMinCut(int overlapType, int dstY, int dstX, int srcY, int srcX);

    // Compute the overlap between the current block, dst, of the output image
    // and src, of the input image, given their upper-left corners
    // and the position of the overlap
    double ComputeOverlapBasicOpt(int overlapType, int dstY, int dstX, int srcY, int srcX);

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
    void PlaceEdgeOverlapBlockWithMinCut(int blockY, int blockX, int maxBlockX, int maxBlockY, double errorTolerance);

    // Synthesize a new texture by randomly choosing blocks satisfying constraints and applying minimum cuts
    void OverlapConstraintsWithMinCut();

    int overlapWidth = 0;
    int overlapHeight = 0;

    int64_t flopCount = 0;
};
