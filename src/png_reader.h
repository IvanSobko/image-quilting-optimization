#pragma once

struct ImgData;

namespace file {
// Read a png file into data.data
void read_png_file(char const* filename, ImgData& data);

// Write the png file and frees both the data.data and data.output_d
void write_png_file(char const* filename, ImgData& data);

}  // namespace file
