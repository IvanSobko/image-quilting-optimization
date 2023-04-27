#include "ImageQuilting.h"
#include <algorithm>
#include <cstdlib>
#include <random>

ImgData ImageQuilting::synthesis() {
    mData.output_d = (unsigned char**)malloc(sizeof(unsigned char*) * mData.output_h);
    for (int y = 0; y < mData.output_h; y++) {
        mData.output_d[y] = (unsigned char*)malloc(mData.output_w * 4);
    }
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<std::mt19937::result_type> x_rand(0, mData.width - mData.block_w - 1);
    std::uniform_int_distribution<std::mt19937::result_type> y_rand(0, mData.height - mData.block_h - 1);
    for (int y = 0; y < mData.output_h; y += mData.block_h) {
        for (int x = 0; x < mData.output_w; x += mData.block_w) {
            int x_start = x_rand(rng);
            int y_start = y_rand(rng);

            for (int yo = y; yo < std::min(y + mData.block_h, mData.output_h); yo++) {
                for (int xo = x; xo < std::min(x + mData.block_w, mData.output_w); xo++) {
                    int yb = y_start + yo - y;
                    int xb = (x_start + xo - x) * 4;
                    mData.output_d[yo][xo * 4] = mData.data[yb][xb];
                    mData.output_d[yo][xo * 4 + 1] = mData.data[yb][xb + 1];
                    mData.output_d[yo][xo * 4 + 2] = mData.data[yb][xb + 2];
                    mData.output_d[yo][xo * 4 + 3] = mData.data[yb][xb + 3];
                }
            }
        }
    }
    return mData;
}