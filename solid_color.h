#ifndef SOLID_COLOR_H
#define SOLID_COLOR_H

#include "texture.h"

class SolidColor : public Texture {
public:
    SolidColor() {}
    SolidColor(vec3 c) : color(c) {}

    virtual vec3 value(double u, double v, const vec3& p) const override {
        (void)u;
        (void)v;
        (void)p;
        return color;
    }

private:
    vec3 color;
};

#endif
