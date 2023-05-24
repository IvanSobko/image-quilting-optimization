#pragma once

#include "ImgData.h"

class Blocking {
   public:

    // Struct to keep track of blocks and an associated value
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

    // Enum representing the different optimization variations
    enum OptType {
        none = 0,
        refactor = 1,
        blocking = 2
    };
    OptType optType = none;

    // Set custom image quilting parameters
    static void SetParameters(ImgData * imgData);

    // Get the parameters required to call a component test function
    static void GetComponentParameters(
        ImgData* imgData,
        int & dstY, int & dstX, int & maxBlockY, int & maxBlockX, int & overlapType,
        BlockValue** blockValues);

    // Refactor computing the block values into its own function
    static void Base(ImgData* imgData, int seed);
    static void BaseComponent(ImgData* imgData, int seed);

    // Refactor computing the potential block errors to iterate over the output overlap region exactly once
    void ComputeBlockValuesRefactor(
        int dstY, int dstX, int maxBlockY, int maxBlockX, int overlapType,
        BlockValue* blockValues);
    static void Refactor(ImgData* imgData, int seed);
    static void RefactorComponent(ImgData* imgData, int seed);

    // Block computing the block values
    void ComputeBlockValuesBlocked(
        int dstY, int dstX, int maxBlockY, int maxBlockX, int overlapType,
        BlockValue* blockValues);
    static void Blocked(ImgData* imgData, int seed);
    static void BlockedComponent(ImgData* imgData, int seed);

    // Compute the block values for a given destination block
    void ComputeBlockValues(
        int dstY, int dstX, int maxBlockY, int maxBlockX, int overlapType,
        BlockValue* blockValues);

    Blocking() = delete;
    Blocking(ImgData* data);

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
    int overlapWidth = 0;
    int overlapHeight = 0;
    int64_t flopCount = 0;

    // Write a block from the source data to the output data given their upper-left corners
    void WriteBlock(int dstY, int dstX, int srcY, int srcX);

    // Same as WriteBlockOverlap, but uses a minimum cut to write the new block
    void WriteBlockOverlapWithMinCut(int overlapType, int dstY, int dstX, int srcY, int srcX);

    // Compute the overlap between the current block, dst, of the output image
    // and src, of the input image, given their upper-left corners
    // and the position of the overlap
    double ComputeOverlap(int overlapType, int dstY, int dstX, int srcY, int srcX);

    // Place an edge overlap block with respect to the given block of the output image
    void PlaceEdgeOverlapBlockWithMinCut(int blockY, int blockX, int maxBlockX, int maxBlockY, double errorTolerance);

    // Synthesize a new texture by randomly choosing blocks satisfying constraints and applying minimum cuts
    void OverlapConstraintsWithMinCut();
};