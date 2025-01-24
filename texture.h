#ifndef TEXTURE_H
#define TEXTURE_H

#include "vec3.h"

class Texture {
public:
    virtual ~Texture() = default;

    // Returns the color at given UV coordinates and position
    virtual vec3 value(double u, double v, const vec3& p) const = 0;
};

#endif
