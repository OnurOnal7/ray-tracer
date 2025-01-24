#ifndef BLEND_TEXTURE_H
#define BLEND_TEXTURE_H

#include "texture.h"
#include <memory>

class BlendTexture : public Texture {
public:
    BlendTexture(std::shared_ptr<Texture> base, std::shared_ptr<Texture> overlay, double blend_factor = 0.5)
        : base(base), overlay(overlay), blend(blend_factor) {}

    virtual vec3 value(double u, double v, const vec3& p) const override {
        vec3 base_color = base->value(u, v, p);
        vec3 overlay_color = overlay->value(u, v, p);
        return (1.0 - blend) * base_color + blend * overlay_color;
    }

private:
    std::shared_ptr<Texture> base;
    std::shared_ptr<Texture> overlay;
    double blend;
};

#endif 
