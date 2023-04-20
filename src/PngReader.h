#ifndef TEAM19_PNGREADER_H
#define TEAM19_PNGREADER_H

struct ImgData;

namespace file {
    void read_png_file(char const *filename, ImgData &data);
    void write_png_file(char const *filename, ImgData &data);
}

#endif //TEAM19_PNGREADER_H
