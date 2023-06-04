#pragma once

#include "ImgData.h"

class AdvanceAlgOptimiz {
public:
    AdvanceAlgOptimiz() = delete;
    AdvanceAlgOptimiz(ImgData* data) {
        mData = data;
        overlapHeight = mData->block_h / 6;
        overlapWidth = mData->block_w / 6;
    }

    // Seed the random number generator with the system time
    static void SeedRandomNumberGenerator();
    // Seed the random number generator with a specified seed
    static void SeedRandomNumberGenerator(int seed);
    // Generate a random number in the range [min, max]
    static int GetRandomInt(int min, int max);

    static void CustomParameters(ImgData* imgData);
    static void StdC_KUnroll_BoundsRefactor(ImgData * imgData, int seed);
    static void StdC_KUnroll_BoundsRefactor_LoopReorder(ImgData * imgData, int seed);
    static void StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking32(ImgData * imgData, int seed);
    static void StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking48(ImgData * imgData, int seed);
    static void StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking64(ImgData * imgData, int seed);
    static void StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking96(ImgData * imgData, int seed);
    static void StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking128(ImgData * imgData, int seed);

    // Std C, bounds refactoring, loop reorder, blocking 32x32, and unrolling channels loop and srcX by 2
    static void StdC_KSrc2Unroll_BoundsRefactor_LoopReorder_Blocking32(ImgData * imgData, int seed);

    // Std C, bounds refactoring, loop reorder, blocking 32x32, and unrolling channels loop, srcX by 4
    static void StdC_KSrc4Unroll_BoundsRefactor_LoopReorder_Blocking32(ImgData * imgData, int seed);

private:
    // Keep a pointer to the input image data
    ImgData* mData;
    int overlapWidth = 0;
    int overlapHeight = 0;

    // Write a block from the source data to the output data given their upper-left corners
    void WriteBlock(int dstY, int dstX, int srcY, int srcX);

    // Same as WriteBlockOverlap, but uses a minimum cut to write the new block
    void WriteBlockOverlapWithMinCut(int overlapType, int dstY, int dstX, int srcY, int srcX);

    // Struct to sort blocks by their l2 norm
    struct BlockValue {
        int y, x;
        int value = 0;
    };

    // Enum representing the type of overlap between blocks
    enum OverlapType { vertical = 0, horizontal = 1, both = 2 };

    // Std C, K unroll, and bounds refactoring optimizations
    void PlaceEdgeOverlapBlockWithMinCut_StdC_KUnroll_BoundsRefactor(int overlapType, int dstY, int dstX, int maxBlockX, int maxBlockY,
                                         double errorTolerance, int bWidth, int bHeight);
    void OverlapConstraintsWithMinCut_StdC_KUnroll_BoundsRefactor();

    // Std C, K unroll, bounds refactoring, and loop reorder optimizations
    void PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KUnroll_BoundsRefactor_LoopReorder(int overlapType, int dstY, int dstX, int maxBlockX, int maxBlockY,
                                                double errorTolerance, int bWidth, int bHeight);
    void OverlapConstraintsWithMinCut_StdC_KUnroll_BoundsRefactor_LoopReorder();

    // Std C, K unroll, bounds refactoring, loop reorder optimizations, and blocking with 32x32 blocks
    void PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking32(int overlapType, int dstY, int dstX, int maxBlockX, int maxBlockY,
                                                                                                    double errorTolerance, int bWidth, int bHeight);
    void OverlapConstraintsWithMinCut_StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking32();

    // Std C, K unroll, bounds refactoring, loop reorder optimizations, and blocking with 48x48 blocks
    void PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking48(int overlapType, int dstY, int dstX, int maxBlockX, int maxBlockY,
                                                                                                    double errorTolerance, int bWidth, int bHeight);
    void OverlapConstraintsWithMinCut_StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking48();

    // Std C, K unroll, bounds refactoring, loop reorder optimizations, and blocking with 64x64 blocks
    void PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking64(int overlapType, int dstY, int dstX, int maxBlockX, int maxBlockY,
                                                                                                    double errorTolerance, int bWidth, int bHeight);
    void OverlapConstraintsWithMinCut_StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking64();

    // Std C, K unroll, bounds refactoring, loop reorder optimizations, and blocking with 96x96 blocks
    void PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking96(int overlapType, int dstY, int dstX, int maxBlockX, int maxBlockY,
                                                                                                    double errorTolerance, int bWidth, int bHeight);
    void OverlapConstraintsWithMinCut_StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking96();

    // Std C, K unroll, bounds refactoring, loop reorder optimizations, and blocking with 128x128 blocks
    void PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking128(int overlapType, int dstY, int dstX, int maxBlockX, int maxBlockY,
                                                                                                    double errorTolerance, int bWidth, int bHeight);
    void OverlapConstraintsWithMinCut_StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking128();

    // Std C, bounds refactoring, loop reorder, blocking 32x32, and unrolling channels loop and srcX by 2
    void PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KSrc2Unroll_BoundsRefactor_LoopReorder_Blocking32(int overlapType, int dstY, int dstX, int maxBlockX, int maxBlockY,
                                                                                                        double errorTolerance, int bWidth, int bHeight);
    void OverlapConstraintsWithMinCut_StdC_KSrc2Unroll_BoundsRefactor_LoopReorder_Blocking32();

    // Std C, bounds refactoring, loop reorder, blocking 32x32, and unrolling channels loop, srcX by 4
    void PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KSrc4Unroll_BoundsRefactor_LoopReorder_Blocking32(int overlapType, int dstY, int dstX, int maxBlockX, int maxBlockY,
                                                                                                        double errorTolerance, int bWidth, int bHeight);
    void OverlapConstraintsWithMinCut_StdC_KSrc4Unroll_BoundsRefactor_LoopReorder_Blocking32();
};
