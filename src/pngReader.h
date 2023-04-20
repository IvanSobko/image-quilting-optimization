#ifndef TEAM19_PNGREADER_H
#define TEAM19_PNGREADER_H

#include <cstdlib>
#include <png.h>

namespace file {
    void read_png_file(char const *filename, png_bytepp &data, uint32_t &width, uint32_t &height);
    void write_png_file(char const *filename, png_bytepp &data, uint32_t width, uint32_t height);
}

#endif //TEAM19_PNGREADER_H
