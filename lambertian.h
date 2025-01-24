#ifndef LAMBERTIAN_H
#define LAMBERTIAN_H

#include "material.h"
#include "vec3.h"
#include "hit_record.h"
#include "ray.h"
#include "texture.h"         
#include "solid_color.h"      
#include "image_texture.h"   
#include <memory>

inline vec3 random_cosine_direction() {
    double r1 = random_double();
    double r2 = random_double();
    double z = sqrt(1 - r2);
    double phi = 2*M_PI*r1;
    double x = cos(phi)*sqrt(r2);
    double y = sin(phi)*sqrt(r2);
    return vec3(x,y,z);
}

struct onb {
    vec3 u,v,w;
    onb() {}
    onb(const vec3& n) { build_from_w(n); }
    void build_from_w(const vec3& n) {
        w = unit_vector(n);
        vec3 a = (fabs(w.x())>0.9)? vec3(0,1,0): vec3(1,0,0);
        v = unit_vector(cross(w,a));
        u = cross(w,v);
    }
    vec3 local(double a, double b, double c) const {
        return a*u + b*v + c*w;
    }
    vec3 local(const vec3& a) const {
        return a.x()*u + a.y()*v + a.z()*w;
    }
};

// Lambertian (Diffuse) Material
class Lambertian : public Material {
public:
    // Constructor for solid color
    Lambertian(const vec3& a) : albedo(std::make_shared<SolidColor>(a)) {}

    // Constructor for texture
    Lambertian(std::shared_ptr<Texture> a) : albedo(a) {}

    virtual bool scatter(const Ray& r_in, const HitRecord& rec, vec3& attenuation, Ray& scattered) const override {
        (void)r_in; 
        vec3 scatter_direction = rec.normal + random_unit_vector();
        if (scatter_direction.near_zero()) scatter_direction = rec.normal;
        scattered = Ray(rec.p, scatter_direction);
        attenuation = albedo->value(rec.u, rec.v, rec.p);
        return true;
    }

    virtual bool sample(const Ray& r_in, const HitRecord& rec, vec3& attenuation, Ray& scattered) const override {
        (void)r_in;
        onb uvw(rec.normal);
        vec3 direction = uvw.local(random_cosine_direction());
        scattered = Ray(rec.p, unit_vector(direction));
        attenuation = albedo->value(rec.u, rec.v, rec.p);
        return true;
    }

    virtual double scatter_pdf(const Ray& r_in, const HitRecord& rec, const Ray& scattered) const override {
        (void)r_in;
        double cosine = dot(rec.normal, unit_vector(scattered.direction()));
        return (cosine <= 0) ? 0 : cosine/M_PI;
    }

    virtual vec3 emit(double u, double v, const vec3& p) const override {
        (void)u;(void)v;(void)p;
        return vec3(0.0, 0.0, 0.0); // No emission for Lambertian materials
    }

private:
    std::shared_ptr<Texture> albedo; // Diffuse reflectivity
};

#endif 
