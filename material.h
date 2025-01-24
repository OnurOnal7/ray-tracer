#ifndef MATERIAL_H
#define MATERIAL_H

#include "ray.h"
#include "hit_record.h"
#include "vec3.h"
#include <memory>

class Material {
public:
    virtual ~Material() = default;

    virtual bool scatter(const Ray& r_in, const HitRecord& rec, vec3& attenuation, Ray& scattered) const = 0;
    
    virtual vec3 emit(double u, double v, const vec3& p) const {
        (void)u; (void)v; (void)p;
        return vec3(0, 0, 0);
    }

    virtual double scatter_pdf(const Ray& r_in, const HitRecord& rec, const Ray& scattered) const {
        (void)r_in; (void)rec; (void)scattered;
        return 0.0; 
    }

    virtual bool sample(const Ray& r_in, const HitRecord& rec, vec3& attenuation, Ray& scattered) const {
        return scatter(r_in, rec, attenuation, scattered);
    }


};

#endif
