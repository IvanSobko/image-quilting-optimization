#pragma once

#include "img_data.h"

class AdvanceAlgOptimiz {
public:
    AdvanceAlgOptimiz() = delete;
    AdvanceAlgOptimiz(ImgData* data) {
        mData = data;
        overlapHeight = mData->block_h / 6;
        overlapWidth = mData->block_w / 6;
    }

    int64_t calcFlops();

    static void CustomParameters(ImgData* imgData);
    static double StdC_KUnroll_BoundsRefactor(ImgData* imgData, int seed);
    static double StdC_KUnroll_BoundsRefactor_LoopReorder(ImgData* imgData, int seed);
    static double StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking32(ImgData* imgData, int seed);
    static double StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking48(ImgData* imgData, int seed);
    static double StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking64(ImgData* imgData, int seed);
    static double StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking96(ImgData* imgData, int seed);
    static double StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking128(ImgData* imgData, int seed);

    // Std C, bounds refactoring, loop reorder, blocking 32x32, and unrolling channels loop and srcX by 2
    static double StdC_KSrc2Unroll_BoundsRefactor_LoopReorder_Blocking32(ImgData* imgData, int seed);

    // Std C, bounds refactoring, loop reorder, blocking 32x32, and unrolling channels loop, srcX by 4
    static double StdC_KSrc4Unroll_BoundsRefactor_LoopReorder_Blocking32(ImgData* imgData, int seed);

#ifdef __AVX2__
    // Vectorization of KUnroll_BoundsRefactor_LoopReorder_Blocking32 with unroll of 2,4,8
    static double StdC_KSrc2Unroll_Vector_BoundsRefactor_LoopReorder_Blocking32(ImgData* imgData, int seed);
    static double StdC_KSrc4Unroll_Vector_BoundsRefactor_LoopReorder_Blocking32(ImgData* imgData, int seed);
    static double StdC_KSrc8Unroll_Vector_BoundsRefactor_LoopReorder_Blocking32(ImgData* imgData, int seed);
#endif
private:
    // Enum representing the type of overlap between blocks
    enum OverlapType { vertical = 0, horizontal = 1, both = 2 };

    enum OptType {
        opt_unroll_bounds,
        opt_unroll_bounds_loop,
        opt_unroll_bounds_loop_blocking32,
        opt_unroll_bounds_loop_blocking48,
        opt_unroll_bounds_loop_blocking64,
        opt_unroll_bounds_loop_blocking96,
        opt_unroll_bounds_loop_blocking128,
        opt_unroll2_bounds_loop_blocking32,
        opt_unroll4_bounds_loop_blocking32,
        opt_unroll2_simd_bounds_loop_blocking32,
        opt_unroll4_simd_bounds_loop_blocking32,
        opt_unroll8_simd_bounds_loop_blocking32,
    };

    // Struct to sort blocks by their l2 norm
    struct BlockValue {
        int y, x;
        int value = 0;
    };

    // Seed the random number generator with a specified seed
    //TODO: should not be static
    static void SeedRandomNumberGenerator(int seed);
    // Generate a random number in the range [min, max]
    int GetRandomInt(int min, int max);

    // Write a block from the source data to the output data given their upper-left corners
    void WriteBlock(int dstY, int dstX, int srcY, int srcX);

    // Same as WriteBlockOverlap, but uses a minimum cut to write the new block
    void WriteBlockOverlapWithMinCut(int overlapType, int dstY, int dstX, int srcY, int srcX);

    void OverlapConstraintsWithMinCut(int opt_type);

    // Std C, K unroll, and bounds refactoring optimizations
    void PlaceEdgeOverlapBlockWithMinCut_StdC_KUnroll_BoundsRefactor(int overlapType, int dstY, int dstX,
                                                                     int maxBlockX, int maxBlockY, int bWidth,
                                                                     int bHeight);

    // Std C, K unroll, bounds refactoring, and loop reorder optimizations
    void PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KUnroll_BoundsRefactor_LoopReorder(
        int overlapType, int dstY, int dstX, int maxBlockX, int maxBlockY, int bWidth, int bHeight);

    // Std C, K unroll, bounds refactoring, loop reorder optimizations, and blocking with 32x32 blocks
    void PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking32(
        int overlapType, int dstY, int dstX, int maxBlockX, int maxBlockY, int bWidth, int bHeight);

    // Std C, K unroll, bounds refactoring, loop reorder optimizations, and blocking with 48x48 blocks
    void PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking48(
        int overlapType, int dstY, int dstX, int maxBlockX, int maxBlockY, int bWidth, int bHeight);

    // Std C, K unroll, bounds refactoring, loop reorder optimizations, and blocking with 64x64 blocks
    void PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking64(
        int overlapType, int dstY, int dstX, int maxBlockX, int maxBlockY, int bWidth, int bHeight);

    // Std C, K unroll, bounds refactoring, loop reorder optimizations, and blocking with 96x96 blocks
    void PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking96(
        int overlapType, int dstY, int dstX, int maxBlockX, int maxBlockY, int bWidth, int bHeight);

    // Std C, K unroll, bounds refactoring, loop reorder optimizations, and blocking with 128x128 blocks
    void PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking128(
        int overlapType, int dstY, int dstX, int maxBlockX, int maxBlockY, int bWidth, int bHeight);

    // Std C, bounds refactoring, loop reorder, blocking 32x32, and unrolling channels loop and srcX by 2
    void PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KSrc2Unroll_BoundsRefactor_LoopReorder_Blocking32(
        int overlapType, int dstY, int dstX, int maxBlockX, int maxBlockY, int bWidth, int bHeight);

    // Std C, bounds refactoring, loop reorder, blocking 32x32, and unrolling channels loop, srcX by 4
    void PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KSrc4Unroll_BoundsRefactor_LoopReorder_Blocking32(
        int overlapType, int dstY, int dstX, int maxBlockX, int maxBlockY, int bWidth, int bHeight);

#ifdef __AVX2__
    // Vectorizing unroll 2
    void
    PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KSrc2Unroll_Vector_BoundsRefactor_LoopReorder_Blocking32(
        int overlapType, int dstY, int dstX, int maxBlockX, int maxBlockY, int bWidth, int bHeight);

    // Vectorizing unroll 4
    void
    PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KSrc4Unroll_Vector_BoundsRefactor_LoopReorder_Blocking32(
        int overlapType, int dstY, int dstX, int maxBlockX, int maxBlockY, int bWidth, int bHeight);

    // Vectorizing unroll 8
    void
    PlaceEdgeOverlapBlockWithMinCutBlocking_StdC_KSrc8Unroll_Vector_BoundsRefactor_LoopReorder_Blocking32(
        int overlapType, int dstY, int dstX, int maxBlockX, int maxBlockY, int bWidth, int bHeight);
#endif
    static constexpr double errorTolerance = 0.1;

    // Keep a pointer to the input image data
    ImgData* mData;
    int overlapWidth = 0;
    int overlapHeight = 0;
    int64_t flopCount = 0;
};
