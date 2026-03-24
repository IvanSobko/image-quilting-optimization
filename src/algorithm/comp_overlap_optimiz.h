#pragma once

#include "image_quilting.h"
#include "img_data.h"

class CompOverlapOptimiz : public ImageQuilting {
public:
    CompOverlapOptimiz() = delete;
    CompOverlapOptimiz(ImgData* data) : ImageQuilting(data) {}

    // Synthesize a new texture with the given seed
    void Synthesis(int seed, int opt);

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

private:
    enum OptType {
        opt_indices = 0,
        opt_algorithm = 1,
        opt_unroll = 2,
        opt_unroll_max = 3,
        opt_vectorize = 4,
        opt_unroll_chnls = 5,
    };

    // Synthesize a new texture by randomly choosing blocks satisfying constraints and applying minimum cuts
    void OverlapConstraintsWithMinCut();

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

    int opt_type = 0;
};
