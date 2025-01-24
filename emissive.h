#ifndef EMISSIVE_H
#define EMISSIVE_H

#include "material.h"
#include "vec3.h"
#include "hit_record.h"
#include "ray.h"
#include <memory>

// Emissive Material
class Emissive : public Material {
public:
    Emissive(const vec3& color) : emit_color(color) {}

    virtual bool scatter(const Ray&, const HitRecord&, vec3&, Ray&) const override {
        return false; // No scattering for emissive materials
    }

    virtual bool sample(const Ray&, const HitRecord&, vec3&, Ray&) const override {
        return false;
    }

    virtual double scatter_pdf(const Ray&, const HitRecord&, const Ray&) const override {
        return 0.0;
    }

    virtual vec3 emit(double u, double v, const vec3& p) const override {
        (void)u; (void)v; (void)p;
        return emit_color;
    }

private:
    vec3 emit_color; // Emission color
};

#endif 
