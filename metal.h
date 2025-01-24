#ifndef METAL_H
#define METAL_H

#include "material.h"
#include "vec3.h"
#include "hit_record.h"
#include "ray.h"
#include "texture.h"        
#include "solid_color.h"   
#include "image_texture.h"  
#include <memory>
#include <algorithm>

// Metal (Specular) Material
class Metal : public Material {
public:
    // Constructor for solid color
    Metal(const vec3& a, double f)
        : albedo(std::make_shared<SolidColor>(a)), fuzz(f < 1 ? f : 1) {}

    // Constructor for texture
    Metal(std::shared_ptr<Texture> a, double f)
        : albedo(a), fuzz(f < 1 ? f : 1) {}

    virtual bool scatter(const Ray& r_in, const HitRecord& rec, vec3& attenuation, Ray& scattered) const override {
        (void)r_in; 
        vec3 reflected = reflect(unit_vector(r_in.direction()), rec.normal);
        scattered = Ray(rec.p, reflected + fuzz * random_in_unit_sphere());
        attenuation = albedo->value(rec.u, rec.v, rec.p);
        return (dot(scattered.direction(), rec.normal) > 0);
    }

    virtual vec3 emit(double u, double v, const vec3& p) const override {
        (void)u; (void)v; (void)p;
        return vec3(0.0, 0.0, 0.0); // No emission for Metal materials
    }

private:
    std::shared_ptr<Texture> albedo; // Reflective color
    double fuzz; // Fuzziness factor
};

#endif 
