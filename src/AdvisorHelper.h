#pragma once

#include "ImgData.h"

// final versions of optimisations
namespace Advisor {
    void baseline(ImgData* imgData, int seed);
    void basicOpt(ImgData* imgData, int seed);
    void unrollMemory(ImgData* imgData, int seed);
    void unrollChnls(ImgData* imgData, int seed);
    void unroll(ImgData* imgData, int seed);
    void vectorize(ImgData* imgData, int seed);
    void block(ImgData* imgData, int seed);
};
