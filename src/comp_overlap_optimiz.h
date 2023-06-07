#pragma once

#include "img_data.h"

class CompOverlapOptimiz {
public:
    CompOverlapOptimiz() = delete;
    CompOverlapOptimiz(ImgData* data) {
        mData = data;
        overlapHeight = mData->block_h / 6;
        overlapWidth = mData->block_w / 6;
    }

    // Algorithm test functions
    static double BasicOpt(ImgData* imgData, int seed);
    static double AlgOpt(ImgData* imgData, int seed);
    static double UnrollOpt(ImgData* imgData, int seed);
    static double UnrollMaxOpt(ImgData* imgData, int seed);
#ifdef __AVX2__
    static double VectorizeOpt(ImgData* imgData, int seed);
#endif
    static double UnrollChnls(ImgData* imgData, int seed);

    // Component test functions
    static void GetComponentParameters(ImgData* imgData, int& overlapType, int& dstY, int& dstX, int& srcY,
                                       int& srcX);
    volatile static void BaseComponent(ImgData* imgData, int seed);
    volatile static void BasicOptComponent(ImgData* imgData, int seed);
    volatile static void AlgoOptComponent(ImgData* imgData, int seed);
    volatile static void UnrollOptComponent(ImgData* imgData, int seed);
    volatile static void UnrollMaxOptComponent(ImgData* imgData, int seed);
#ifdef __AVX2__
    volatile static void VectorizeOptComponent(ImgData* imgData, int seed);
#endif

    // Synthesize a new texture with the given seed
    void Synthesis(int seed, int opt);

    // Seed the random number generator with the system time
    static void SeedRandomNumberGenerator();
    // Seed the random number generator with a specified seed
    static void SeedRandomNumberGenerator(int seed);
    // Generate a random number in the range [min, max]
    static int GetRandomInt(int min, int max);

    int64_t getFlopCount() const;

private:
    enum OptType {
        opt_indices = 0,
        opt_algorithm = 1,
        opt_unroll = 2,
        opt_unroll_max = 3,
        opt_vectorize = 4,
        opt_unroll_chnls = 5,
    };

    // Struct to sort blocks by their l2 norm
    struct BlockValue {
        int y, x;
        double value;
    };

    // Enum representing the type of overlap between blocks
    enum OverlapType { vertical = 0, horizontal = 1, both = 2 };

    // Write a block from the source data to the output data given their upper-left corners
    void WriteBlock(int dstY, int dstX, int srcY, int srcX);

    // Same as WriteBlock but leaves the half of the dst overlapping region untouched
    void WriteBlockOverlap(int overlapType, int dstY, int dstX, int srcY, int srcX);

    // Same as WriteBlockOverlap, but uses a minimum cut to write the new block
    void WriteBlockOverlapWithMinCut(int overlapType, int dstY, int dstX, int srcY, int srcX);

    // Base implementation of ComputeOverlap
    double ComputeOverlapBase(int overlapType, int dstY, int dstX, int srcY, int srcX);

    // Compute the overlap with indices optimization and optimized data access
    double ComputeOverlapBasicOpt(int overlapType, int dstY, int dstX, int srcY, int srcX);

    // Tried to improve the overlap calculations additionally to ComputeOverlapBasicOpt
    double ComputeOverlapAlgImpr(int overlapType, int dstY, int dstX, int srcY, int srcX);

    // Added unroll and ILP to ComputeOverlapAlgImpr
    double ComputeOverlapUnroll(int overlapType, int dstY, int dstX, int srcY, int srcX);
    double ComputeOverlapUnrollChannels(int overlapType, int dstY, int dstX, int srcY, int srcX);

    // Forced compiler to parallelize computations
    double ComputeOverlapUnrollMax(int overlapType, int dstY, int dstX, int srcY, int srcX);

#ifdef __AVX2__
    double ComputeOverlapVectorize(int overlapType, int dstY, int dstX, int srcY, int srcX);
#endif

    // Place an edge overlap block with respect to the given block of the output image
    void PlaceEdgeOverlapBlockWithMinCut(int blockY, int blockX, int maxBlockX, int maxBlockY,
                                         double errorTolerance);

    // Synthesize a new texture by randomly choosing blocks satisfying constraints and applying minimum cuts
    void OverlapConstraintsWithMinCut();

    int overlapWidth = 0;
    int overlapHeight = 0;

    int64_t flopCount = 0;
    int opt_type = 0;

    // Keep a pointer to the input image data
    ImgData* mData;
};
