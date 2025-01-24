#ifndef IMAGE_TEXTURE_H
#define IMAGE_TEXTURE_H

#include "texture.h"
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>

class ImageTexture : public Texture {
public:
    ImageTexture() {}
    ImageTexture(const std::vector<std::vector<vec3>>& img, int w, int h)
        : image(img), width(w), height(h) {}

    virtual vec3 value(double u, double v, const vec3& p) const override {
        (void)p; // 'p' is unused
        
        // Handle out-of-bounds UV coordinates by wrapping
        u = u - std::floor(u);
        v = v - std::floor(v);

        // Convert UV to floating-point pixel coordinates
        double x = u * (width - 1);
        double y = (1.0 - v) * (height - 1); // Flip v to match image coordinates

        int i = static_cast<int>(x);
        int j = static_cast<int>(y);

        // Clamp indices to image dimensions
        i = std::min(std::max(i, 0), width - 1);
        j = std::min(std::max(j, 0), height - 1);

        return image[j][i];
    }

private:
    std::vector<std::vector<vec3>> image;
    int width;
    int height;
};

#endif 
