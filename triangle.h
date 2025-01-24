#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "hittable.h"
#include "vec3.h"
#include "ray.h"
#include "hit_record.h"
#include "material.h"
#include <memory>
#include <cmath>

// Triangle Class Definition and Implementation
class Triangle : public Hittable {
public:
    // Constructors
    Triangle() {}
    Triangle(const vec3& v0, const vec3& v1, const vec3& v2,
             const vec3& uv0, const vec3& uv1, const vec3& uv2,
             const vec3& n0, const vec3& n1, const vec3& n2,
             std::shared_ptr<Material> m)
        : v0(v0), v1(v1), v2(v2),
          uv0(uv0), uv1(uv1), uv2(uv2),
          n0(n0), n1(n1), n2(n2),
          mat_ptr(m) {}

    // Override hit function using the Möller–Trumbore intersection algorithm
    virtual bool hit(const Ray& r, double t_min, double t_max, HitRecord& rec) const override {
        const double EPSILON = 1e-8;
        vec3 edge1 = v1 - v0;
        vec3 edge2 = v2 - v0;
        vec3 h = cross(r.direction(), edge2);
        
        double a = dot(edge1, h);
        if (a > -EPSILON && a < EPSILON)
            return false;

        double f = 1.0 / a;
        vec3 s = r.origin() - v0;
        double u = f * dot(s, h);
        if (u < 0.0 || u > 1.0)
            return false;

        vec3 q = cross(s, edge1);
        double v = f * dot(r.direction(), q);
        if (v < 0.0 || u + v > 1.0)
            return false;

        double t = f * dot(edge2, q);
        if (t < t_min || t > t_max)
            return false;

        rec.t = t;
        rec.p = r.at(t);

        // Interpolate normals
        double w = 1.0 - u - v;
        vec3 interpolated_normal = unit_vector(n0 * w + n1 * u + n2 * v);

        rec.set_face_normal(r, interpolated_normal);
        rec.mat_ptr = mat_ptr;

        // Interpolate UV using barycentric coords
        double U = uv0.x()*w + uv1.x()*u + uv2.x()*v;
        double V = uv0.y()*w + uv1.y()*u + uv2.y()*v;
        rec.u = U;
        rec.v = V;

        return true;
    }

    // Override bounding_box function
    virtual bool bounding_box(AABB& output_box) const override {
        vec3 min_pt(
            std::min({v0.x(), v1.x(), v2.x()}),
            std::min({v0.y(), v1.y(), v2.y()}),
            std::min({v0.z(), v1.z(), v2.z()})
        );
        vec3 max_pt(
            std::max({v0.x(), v1.x(), v2.x()}),
            std::max({v0.y(), v1.y(), v2.y()}),
            std::max({v0.z(), v1.z(), v2.z()})
        );

        double epsilon = 1e-4;
        min_pt -= vec3(epsilon, epsilon, epsilon);
        max_pt += vec3(epsilon, epsilon, epsilon);

        output_box = AABB(min_pt, max_pt);
        return true;
    }

private:
    vec3 v0, v1, v2;
    vec3 uv0, uv1, uv2;
    vec3 n0, n1, n2; 
    std::shared_ptr<Material> mat_ptr;
};

#endif 
