#ifndef SPHERE_H
#define SPHERE_H

#include "hittable.h"
#include "vec3.h"
#include "ray.h"
#include "hit_record.h"
#include "material.h"
#include <memory>
#include <cmath>

// Static Sphere Class
class Sphere : public Hittable {
public:
    // Constructors
    Sphere() {}
    Sphere(vec3 cen, double r, std::shared_ptr<Material> m) 
        : center(cen), radius(r), mat_ptr(m) {}
    
    // Override hit function
    virtual bool hit(const Ray& r, double t_min, double t_max, HitRecord& rec) const override;
    virtual bool bounding_box(AABB& output_box) const override;

    // Getters
    vec3 get_center() const { return center; }
    double get_radius() const { return radius; }
    std::shared_ptr<Material> get_material() const { return mat_ptr; }

protected:
    vec3 center;
    double radius;
    std::shared_ptr<Material> mat_ptr;

    // Function to calculate UV coordinates for a sphere
    void get_sphere_uv(const vec3& p, double& u, double& v) const;
};

// Moving Sphere Class
class MovingSphere : public Hittable {
public:
    // Constructors
    MovingSphere() {}
    MovingSphere(vec3 cen0, vec3 cen1, double time0, double time1, double r, std::shared_ptr<Material> m)
        : center0(cen0), center1(cen1), time0(time0), time1(time1), radius(r), mat_ptr(m) {}
    
    // Override hit function using the ray's time to determine sphere position
    virtual bool hit(const Ray& r, double t_min, double t_max, HitRecord& rec) const override;
    virtual bool bounding_box(AABB& output_box) const override;

    // Getters
    vec3 get_center(double time) const { 
        return center0 + ((time - time0) / (time1 - time0)) * (center1 - center0);
    }

private:
    vec3 center0, center1; // Start and end centers
    double time0, time1;   // Start and end times
    double radius;
    std::shared_ptr<Material> mat_ptr;

    // Function to calculate UV coordinates for a sphere
    void get_sphere_uv(const vec3& p, double& u, double& v) const;
};

/*
* Implementations
*/

// Sphere hit function
bool Sphere::hit(const Ray& r, double t_min, double t_max, HitRecord& rec) const {
    vec3 oc = r.origin() - center;
    double a = dot(r.direction(), r.direction());
    double half_b = dot(oc, r.direction());
    double c = dot(oc, oc) - radius * radius;
    double discriminant = half_b * half_b - a * c;
    if (discriminant < 0) return false;
    double sqrt_d = sqrt(discriminant);

    // Find the nearest root that lies in the acceptable range.
    double root = (-half_b - sqrt_d) / a;
    if (root < t_min || root > t_max) {
        root = (-half_b + sqrt_d) / a;
        if (root < t_min || root > t_max)
            return false;
    }

    rec.t = root;
    rec.p = r.at(rec.t);
    vec3 outward_normal = (rec.p - center) / radius;
    rec.set_face_normal(r, outward_normal);
    rec.mat_ptr = mat_ptr;

    get_sphere_uv(outward_normal, rec.u, rec.v);

    return true;
}

// Sphere bounding box
bool Sphere::bounding_box(AABB& output_box) const {
    output_box = AABB(
        center - vec3(radius, radius, radius),
        center + vec3(radius, radius, radius)
    );
    return true;
}

// Sphere UV mapping
void Sphere::get_sphere_uv(const vec3& p, double& u, double& v) const {
    double phi = atan2(p.z(), p.x());
    double theta = asin(p.y());
    u = 1 - (phi + M_PI) / (2 * M_PI);
    v = (theta + M_PI / 2) / M_PI;
}

// MovingSphere hit function
bool MovingSphere::hit(const Ray& r, double t_min, double t_max, HitRecord& rec) const {
    vec3 center = get_center(r.time());
    vec3 oc = r.origin() - center;
    double a = dot(r.direction(), r.direction());
    double half_b = dot(oc, r.direction());
    double c = dot(oc, oc) - radius * radius;
    double discriminant = half_b * half_b - a * c;
    if (discriminant < 0) return false;
    double sqrt_d = sqrt(discriminant);

    // Find the nearest root that lies in the acceptable range.
    double root = (-half_b - sqrt_d) / a;
    if (root < t_min || root > t_max) {
        root = (-half_b + sqrt_d) / a;
        if (root < t_min || root > t_max)
            return false;
    }

    rec.t = root;
    rec.p = r.at(rec.t);
    vec3 outward_normal = (rec.p - center) / radius;
    rec.set_face_normal(r, outward_normal);
    rec.mat_ptr = mat_ptr;

    get_sphere_uv(outward_normal, rec.u, rec.v);

    return true;
}

// MovingSphere bounding box
bool MovingSphere::bounding_box(AABB& output_box) const {
    AABB box0(
        get_center(time0) - vec3(radius, radius, radius),
        get_center(time0) + vec3(radius, radius, radius)
    );
    AABB box1(
        get_center(time1) - vec3(radius, radius, radius),
        get_center(time1) + vec3(radius, radius, radius)
    );
    output_box = AABB::surrounding_box(box0, box1);
    return true;
}

// MovingSphere UV mapping
void MovingSphere::get_sphere_uv(const vec3& p, double& u, double& v) const {
    double phi = atan2(p.z(), p.x());
    double theta = asin(p.y());
    u = 1 - (phi + M_PI) / (2 * M_PI);
    v = (theta + M_PI / 2) / M_PI;
}

#endif 
