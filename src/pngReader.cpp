#include "pngReader.h"
#include <cstdio>

// modified code from: https://gist.github.com/niw/5963798
void file::read_png_file(char const *filename, png_bytepp &data, uint32_t &width, uint32_t &height) {
    FILE *fp = fopen(filename, "rb");
    if (!fp) {
        printf("ERROR: can't open the file %s\n", filename);
        return;
    }

    if (data) {
        printf("WARNING: already has values - will clean\n");
        for (int y = 0; y < height; y++) {
            free(data[y]);
        }
        free(data);
    }

    png_structp png = png_create_read_struct(PNG_LIBPNG_VER_STRING,
                                             NULL, NULL, NULL);
    if (!png) {
        printf("ERROR: failed the structure creation\n");
        return;
    }

    png_infop info = png_create_info_struct(png);
    if (!info) {
        printf("ERROR: failed info creation\n");
        abort();
    }

    png_init_io(png, fp);
    png_read_info(png, info);

    width      = png_get_image_width(png, info);
    height     = png_get_image_height(png, info);
    png_byte color_type = png_get_color_type(png, info);
    png_byte bit_depth  = png_get_bit_depth(png, info);

    // Read any color_type into 8bit depth, RGBA format.
    // See http://www.libpng.org/pub/png/libpng-manual.txt

    if (bit_depth == 16) {
        png_set_strip_16(png);
    }

    if (color_type == PNG_COLOR_TYPE_PALETTE) {
        png_set_palette_to_rgb(png);
    }

    // PNG_COLOR_TYPE_GRAY_ALPHA is always 8 or 16bit depth.
    if (color_type == PNG_COLOR_TYPE_GRAY && bit_depth < 8) {
        png_set_expand_gray_1_2_4_to_8(png);
    }

    if (png_get_valid(png, info, PNG_INFO_tRNS)) {
        png_set_tRNS_to_alpha(png);
    }

    // These color_type don't have an alpha channel then fill it with 0xff.
    if(color_type == PNG_COLOR_TYPE_RGB || color_type == PNG_COLOR_TYPE_GRAY ||
       color_type == PNG_COLOR_TYPE_PALETTE) {
        png_set_filler(png, 0xFF, PNG_FILLER_AFTER);
    }

    if (color_type == PNG_COLOR_TYPE_GRAY ||  color_type == PNG_COLOR_TYPE_GRAY_ALPHA) {
        png_set_gray_to_rgb(png);
    }

    png_read_update_info(png, info);
    data = (png_bytep*)malloc(sizeof(png_bytep) * height);
    for(int y = 0; y < height; y++) {
        data[y] = (png_byte*)malloc(png_get_rowbytes(png,info));
    }

    png_read_image(png, data);
    fclose(fp);
    png_destroy_read_struct(&png, &info, NULL);
}

void file::write_png_file(char const *filename, png_bytepp &data, uint32_t width, uint32_t height) {
    FILE *fp = fopen(filename, "wb");
    if (!fp) {
        printf("ERROR: can't open the file %s\n", filename);
        return;
    }
    if (!data) {
        printf("ERROR: nothing to save\n");
        return;
    }
    png_structp png = png_create_write_struct(PNG_LIBPNG_VER_STRING,
                                              NULL, NULL, NULL);
    if (!png) {
        printf("ERROR: failed the structure creation\n");
        return;
    }

    png_infop info = png_create_info_struct(png);
    if (!info) {
        printf("ERROR: failed info creation\n");
        abort();
    }

    png_init_io(png, fp);

    // Output is 8bit depth, RGBA format.
    png_set_IHDR(
            png,
            info,
            width, height,
            8,
            PNG_COLOR_TYPE_RGBA,
            PNG_INTERLACE_NONE,
            PNG_COMPRESSION_TYPE_DEFAULT,
            PNG_FILTER_TYPE_DEFAULT
    );
    png_write_info(png, info);

    // To remove the alpha channel for PNG_COLOR_TYPE_RGB format,
    // Use png_set_filler().
    //png_set_filler(png, 0, PNG_FILLER_AFTER);

    png_write_image(png, data);
    png_write_end(png, NULL);

    for(int y = 0; y < height; y++) {
        free(data[y]);
    }
    free(data);
    fclose(fp);
    png_destroy_write_struct(&png, &info);
}