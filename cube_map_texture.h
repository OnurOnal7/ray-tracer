#ifndef CUBE_MAP_TEXTURE_H
#define CUBE_MAP_TEXTURE_H

#include "texture.h"
#include "vec3.h"
#include <vector>
#include <cmath>
#include <string>

class CubeMapTexture : public Texture {
public:
    CubeMapTexture(const std::vector<std::vector<vec3>>& data, int w, int h) : tex_data(data), tex_width(w), tex_height(h) {}

    virtual vec3 value(double u, double v, const vec3& p) const override {
        (void)u; (void)v; (void)p; 
        return vec3(0,0,0);
    }

    // Given a direction, return color from the environment map
    vec3 sample_direction(const vec3& dir) const {
        // Normalize direction
        vec3 d = unit_vector(dir);

        // Spherical mapping
        double theta = acos(-d.y());
        double phi = atan2(-d.z(), d.x()) + M_PI; 
        double u = phi / (2*M_PI);
        double v = theta / M_PI;

        // Convert u,v to pixel coords
        u = u - floor(u);
        v = v - floor(v);
        double x = u * (tex_width - 1);
        double y = v * (tex_height - 1);

        int x0 = (int)x;
        int y0 = (int)y;
        int x1 = std::min(x0+1, tex_width-1);
        int y1 = std::min(y0+1, tex_height-1);

        double dx = x - x0;
        double dy = y - y0;

        vec3 c00 = tex_data[y0][x0];
        vec3 c10 = tex_data[y0][x1];
        vec3 c01 = tex_data[y1][x0];
        vec3 c11 = tex_data[y1][x1];

        vec3 c0 = c00*(1-dx)+c10*dx;
        vec3 c1 = c01*(1-dx)+c11*dx;
        vec3 color = c0*(1-dy)+c1*dy;

        return color;
    }

private:
    const std::vector<std::vector<vec3>>& tex_data;
    int tex_width, tex_height;
};

#endif
