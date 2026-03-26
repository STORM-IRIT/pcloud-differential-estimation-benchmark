//
// Created by leo on 07/02/25.
//

#ifndef PNGHEADER_H
#define PNGHEADER_H

#include <png.h>
#include <vector>
#include <string>
#include <iostream>

class PNGHandler {
public:
    struct RGBA {
        unsigned char r, g, b, a;
    };

    struct Image {
        std::vector<RGBA> pixels;
        uint32_t width;
        uint32_t height;
    };

    struct BoundingBox {
        uint32_t minX, minY;
        uint32_t maxX, maxY;
    };

    static bool readPNG(const std::string& filename, Image& img) {
        FILE* fp = fopen(filename.c_str(), "rb");
        if (!fp) {
            std::cerr << "Can't open the file " << filename << std::endl;
            return false;
        }

        unsigned char header[8];
        fread(header, 1, 8, fp);
        if (png_sig_cmp(header, 0, 8)) {
            std::cerr << "File " << filename << " isn't a valid PNG " << std::endl;
            fclose(fp);
            return false;
        }

        png_structp png = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
        if (!png) {
            fclose(fp);
            return false;
        }

        png_infop info = png_create_info_struct(png);
        if (!info) {
            png_destroy_read_struct(&png, nullptr, nullptr);
            fclose(fp);
            return false;
        }

        if (setjmp(png_jmpbuf(png))) {
            png_destroy_read_struct(&png, &info, nullptr);
            fclose(fp);
            return false;
        }

        png_init_io(png, fp);
        png_set_sig_bytes(png, 8);
        png_read_info(png, info);

        img.width = png_get_image_width(png, info);
        img.height = png_get_image_height(png, info);
        png_byte color_type = png_get_color_type(png, info);
        png_byte bit_depth = png_get_bit_depth(png, info);

        if (bit_depth == 16)
            png_set_strip_16(png);

        if (color_type == PNG_COLOR_TYPE_PALETTE)
            png_set_palette_to_rgb(png);

        if (color_type == PNG_COLOR_TYPE_GRAY && bit_depth < 8)
            png_set_expand_gray_1_2_4_to_8(png);

        if (png_get_valid(png, info, PNG_INFO_tRNS))
            png_set_tRNS_to_alpha(png);

        if (color_type == PNG_COLOR_TYPE_RGB ||
            color_type == PNG_COLOR_TYPE_GRAY ||
            color_type == PNG_COLOR_TYPE_PALETTE)
            png_set_filler(png, 0xFF, PNG_FILLER_AFTER);

        if (color_type == PNG_COLOR_TYPE_GRAY ||
            color_type == PNG_COLOR_TYPE_GRAY_ALPHA)
            png_set_gray_to_rgb(png);

        png_read_update_info(png, info);

        img.pixels.resize(img.width * img.height);
        std::vector<png_bytep> row_pointers(img.height);

        for (size_t y = 0; y < img.height; y++) {
            row_pointers[y] = reinterpret_cast<png_byte*>(&img.pixels[y * img.width]);
        }

        png_read_image(png, row_pointers.data());

        png_destroy_read_struct(&png, &info, nullptr);
        fclose(fp);
        return true;
    }

    static bool writePNG(const std::string& filename, const Image& img) {
        FILE* fp = fopen(filename.c_str(), "wb");
        if (!fp) {
            std::cerr << "Can't create file " << filename << std::endl;
            return false;
        }

        png_structp png = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
        if (!png) {
            fclose(fp);
            return false;
        }

        png_infop info = png_create_info_struct(png);
        if (!info) {
            png_destroy_write_struct(&png, nullptr);
            fclose(fp);
            return false;
        }

        if (setjmp(png_jmpbuf(png))) {
            png_destroy_write_struct(&png, &info);
            fclose(fp);
            return false;
        }

        png_init_io(png, fp);

        png_set_IHDR(
            png,
            info,
            img.width, img.height,
            8,
            PNG_COLOR_TYPE_RGBA,
            PNG_INTERLACE_NONE,
            PNG_COMPRESSION_TYPE_DEFAULT,
            PNG_FILTER_TYPE_DEFAULT
        );
        png_write_info(png, info);

        std::vector<png_bytep> row_pointers(img.height);
        for (size_t y = 0; y < img.height; y++) {
            row_pointers[y] = reinterpret_cast<png_byte*>(const_cast<RGBA*>(&img.pixels[y * img.width]));
        }

        png_write_image(png, row_pointers.data());
        png_write_end(png, nullptr);

        png_destroy_write_struct(&png, &info);
        fclose(fp);
        return true;
    }

    static BoundingBox findBoundingBoxSquared(const Image& img, const int offset = 25) {
        BoundingBox box = {img.width, img.height, 0, 0};

        for(uint32_t y = 0; y < img.height; y++) {
            for(uint32_t x = 0; x < img.width; x++) {
                const RGBA& pixel = img.pixels[y * img.width + x];
                if(pixel.a > 0) {
                    box.minX = std::min(box.minX, x);
                    box.minY = std::min(box.minY, y);
                    box.maxX = std::max(box.maxX, x);
                    box.maxY = std::max(box.maxY, y);
                }
            }
        }

        int offsetMinX = std::max(0, static_cast<int>(box.minX) - offset);
        int offsetMinY = std::max(0, static_cast<int>(box.minY) - offset);
        int offsetMaxX = std::min(static_cast<int>(img.width), static_cast<int>(box.maxX) + offset);
        int offsetMaxY = std::min(static_cast<int>(img.height), static_cast<int>(box.maxY) + offset);

        int width = offsetMaxX - offsetMinX;
        int height = offsetMaxY - offsetMinY;
        int size = std::max(width, height);

        box.minX = offsetMinX - (size - width) / 2;
        box.minY = offsetMinY - (size - height) / 2;
        box.maxX = box.minX + size;
        box.maxY = box.minY + size;

        return box;
    }

    static BoundingBox findBoundingBox(const Image& img, const int offset = 25) {
        BoundingBox box = {img.width, img.height, 0, 0};

        for(uint32_t y = 0; y < img.height; y++) {
            for(uint32_t x = 0; x < img.width; x++) {
                const RGBA& pixel = img.pixels[y * img.width + x];
                if(pixel.a > 0) {
                    box.minX = std::min(box.minX, x);
                    box.minY = std::min(box.minY, y);
                    box.maxX = std::max(box.maxX, x);
                    box.maxY = std::max(box.maxY, y);
                }
            }
        }

        int minX = std::max(0, static_cast<int>(box.minX) - offset);
        int minY = std::max(0, static_cast<int>(box.minY) - offset);
        int maxX = std::min(static_cast<int>(img.width), static_cast<int>(box.maxX) + offset);
        int maxY = std::min(static_cast<int>(img.height), static_cast<int>(box.maxY) + offset);

        box.minX = minX;
        box.minY = minY;
        box.maxX = maxX;
        box.maxY = maxY;

        return box;
    }

    static bool cropToBox(const Image& input, Image& output, const BoundingBox& box) {
        uint32_t cropWidth  = box.maxX - box.minX;
        uint32_t cropHeight = box.maxY - box.minY;

        output.width = cropWidth;
        output.height = cropHeight;
        output.pixels.resize(cropWidth * cropHeight, RGBA{0, 0, 0, 0});

        for (uint32_t y = 0; y < cropHeight; ++y) {
            for (uint32_t x = 0; x < cropWidth; ++x) {
                uint32_t srcX = box.minX + x;
                uint32_t srcY = box.minY + y;

                if (srcX < input.width && srcY < input.height) {
                    output.pixels[y * cropWidth + x] = input.pixels[srcY * input.width + srcX];
                }
            }
        }

        return true;
    }

    static bool mergeDiagonalSquared(const Image& img1, const Image& img2, Image& result, const int spacing = 5) {
        if(img1.width != img2.width || img1.height != img2.height) {
            std::cerr << "Images must have the same size." << std::endl;
            return false;
        }

        result.width = img1.width;
        result.height = img1.height;
        result.pixels.resize(result.width * result.height, RGBA{0, 0, 0, 0});

        for(uint32_t y = 0; y < result.height; y++) {
            for(uint32_t x = 0; x < result.width; x++) {
                const int diagonalDist = (x + y) - result.width;
                RGBA& pixel = result.pixels[y * result.width + x];

                if(diagonalDist < -spacing) {
                    pixel = img1.pixels[y * img1.width + x];
                }
                else if(diagonalDist > spacing) {
                    pixel = img2.pixels[y * img2.width + x];
                }
            }
        }

        return true;
    }

  static bool cropDiagonal(const Image& img, const bool& top, Image& result, const float spacing = 5) {

      const int w = img.width;
      const int h = img.height;

      result.width = w;
      result.height = h;
      result.pixels.resize(w * h, RGBA{0, 0, 0, 0});

      const float m = -static_cast<float>(h) / w;

      for (int y = 0; y < h; ++y) {
        for (int x = 0; x < w; ++x) {
          float diagonalY = m * x + (h - 1);
          float dist = y - diagonalY;

          RGBA& pixel = result.pixels[y * w + x];

          if (dist < -spacing && top) {
            pixel = img.pixels[y * w + x];
          } else if (dist > spacing && !top) {
            pixel = img.pixels[y * w + x];
          }
        }
      }

      return true;
    }

    static bool mergeDiagonal(const Image& img1, const Image& img2, Image& result, const int spacing = 5) {
        if (img1.width != img2.width || img1.height != img2.height) {
            std::cerr << "Images must have the same size." << std::endl;
            return false;
        }

        const int w = img1.width;
        const int h = img1.height;

        result.width = w;
        result.height = h;
        result.pixels.resize(w * h, RGBA{0, 0, 0, 0});

        const float m = -static_cast<float>(h) / w;

        for (int y = 0; y < h; ++y) {
            for (int x = 0; x < w; ++x) {
                float diagonalY = m * x + (h - 1);
                float dist = y - diagonalY;

                RGBA& pixel = result.pixels[y * w + x];

                if (dist < -spacing) {
                    pixel = img1.pixels[y * w + x];
                } else if (dist > spacing) {
                    pixel = img2.pixels[y * w + x];
                }
            }
        }

        return true;
    }

    static bool mergeTripleDiagonal(const Image& img1, const Image& img2, const Image& img3, Image& result, const int spacing = 5) {
        if (img1.width != img2.width || img1.width != img3.width ||
            img1.height != img2.height || img1.height != img3.height) {
            std::cerr << "All images must have the same size." << std::endl;
            return false;
        }

        // const auto box1 = findBoundingBox(img1);
        // const auto box2 = findBoundingBox(img2);
        // const auto box3 = findBoundingBox(img3);
        //
        // Image cropped1, cropped2, cropped3;
        // cropToBox(img1, cropped1, box1);
        // cropToBox(img2, cropped2, box2);
        // cropToBox(img3, cropped3, box3);

        result.width = img1.width;
        result.height = img1.height;
        result.pixels.resize(result.width * result.height, RGBA{0, 0, 0, 0});

        const int w = result.width;
        const int h = result.height;

        float d1x_top = (5.0f / 12.0f) * w;
        float d1x_bot = (3.0f / 12.0f) * w;

        float d2x_top = (9.0f / 12.0f) * w;
        float d2x_bot = (7.0f / 12.0f) * w;

        for (int y = 0; y < h; ++y) {
            for (int x = 0; x < w; ++x) {
                RGBA& pixel = result.pixels[y * w + x];

                float d1_x = d1x_top + (d1x_bot - d1x_top) * (float(y) / h);
                float d2_x = d2x_top + (d2x_bot - d2x_top) * (float(y) / h);

                float d1_min = d1_x - spacing;
                float d1_max = d1_x + spacing;
                float d2_min = d2_x - spacing;
                float d2_max = d2_x + spacing;

                if (x < d1_min) {
                    pixel = img1.pixels[y * w + x];
                }
                else if (x > d1_max && x < d2_min) {
                    pixel = img2.pixels[y * w + x];
                }
                else if (x > d2_max) {
                    pixel = img3.pixels[y * w + x];
                }
                else {
                    pixel = RGBA{0, 0, 0, 0};
                }
            }
        }

        return true;
    }

    template <typename VectorType>
    static bool isLight(const Image& envMap, const VectorType& _dir, const VectorType& _normal) {
        const VectorType normal = _normal.normalized();
        const VectorType dir = _dir.normalized();

        const VectorType reflected = dir - 2 * dir.dot(normal) * normal;

        const float phi = std::atan2(reflected.y(), reflected.x());
        const float theta = std::acos(reflected.z());

        const int u = static_cast<int>((phi / M_PI + 1.0) * 0.5 * envMap.width);
        const int v = static_cast<int>((theta / M_PI) * envMap.height);

        if (u < 0 || u >= envMap.width || v < 0 || v >= envMap.height) {
            return false;
        }

        const RGBA& pixel = envMap.pixels[v * envMap.width + u];
        const float intensity = (pixel.r + pixel.g + pixel.b) / (3.0 * 255.0);
        return intensity > 0.5;
    }

};

#endif //PNGHEADER_H
