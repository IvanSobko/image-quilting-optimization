#ifndef TEAM19_IMGDATA_H
#define TEAM19_IMGDATA_H

#include <cstdint>
#include <cstdlib>

struct ImgData {
    // each char * stores one row of an image, for each pixel 4 value {rgba}
    unsigned char** data = NULL;
    unsigned char** output_d = NULL;
    uint32_t width = 0, height = 0;
    uint32_t output_w = 0, output_h = 0;
    uint32_t block_w = 0, block_h = 0;
};

#endif  //TEAM19_IMGDATA_H
