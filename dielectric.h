#ifndef DIELECTRIC_H
#define DIELECTRIC_H

#include "material.h"
#include "vec3.h"
#include "hit_record.h"
#include "ray.h"
#include <memory>
#include <cmath>
#include <cstdlib>

// Dielectric (Transparent) Material
class Dielectric : public Material {
public:
    Dielectric(double index_of_refraction) : ir(index_of_refraction) {}

    virtual bool scatter(const Ray& r_in, const HitRecord& rec, vec3& attenuation, Ray& scattered) const override {
        (void)r_in; 
        attenuation = vec3(1.0, 1.0, 1.0); // No attenuation for dielectric
        double refraction_ratio = rec.front_face ? (1.0 / ir) : ir;

        vec3 unit_direction = unit_vector(r_in.direction());
        double cos_theta = fmin(dot(-unit_direction, rec.normal), 1.0);
        double sin_theta = sqrt(1.0 - cos_theta*cos_theta);

        bool cannot_refract = refraction_ratio * sin_theta > 1.0;
        vec3 direction;
        if (cannot_refract || reflectance(cos_theta, refraction_ratio) > random_double()) {
            direction = reflect(unit_direction, rec.normal);
        } 
        else {
            direction = refract(unit_direction, rec.normal, refraction_ratio);
        }

        scattered = Ray(rec.p, direction);
        return true;
    }

    virtual vec3 emit(double u, double v, const vec3& p) const override {
        (void)u; (void)v; (void)p;
        return vec3(0.0, 0.0, 0.0); // No emission for dielectric
    }

private:
    double ir; // Index of Refraction

    // Schlick's approximation for reflectance
    static double reflectance(double cosine, double ref_idx) {
        // Use Schlick's approximation for reflectance
        double r0 = (1 - ref_idx) / (1 + ref_idx);
        r0 = r0 * r0;
        return r0 + (1 - r0) * pow((1 - cosine), 5);
    }
};

#endif
