#ifndef QUAD_H
#define QUAD_H

#include "hittable.h"
#include "vec3.h"
#include "ray.h"
#include "hit_record.h"
#include "material.h"
#include <memory>

class Quad : public Hittable {
public:
    Quad() {}
    Quad(const vec3& a, const vec3& b, const vec3& c, const vec3& d,
         std::shared_ptr<Material> m)
        : v0(a), v1(b), v2(c), v3(d), mat_ptr(m) {}

    virtual bool hit(const Ray& r, double t_min, double t_max, HitRecord& rec) const override {
        // Check hit with triangle1 (v0,v1,v2)
        if (hit_triangle(r, t_min, t_max, rec, v0,v1,v2)) return true;
        // Check triangle2 (v0,v2,v3)
        if (hit_triangle(r, t_min, t_max, rec, v0,v2,v3)) return true;
        return false;
    }

    virtual bool bounding_box(AABB& output_box) const override {
        double minx = std::min({v0.x(), v1.x(), v2.x(), v3.x()});
        double miny = std::min({v0.y(), v1.y(), v2.y(), v3.y()});
        double minz = std::min({v0.z(), v1.z(), v2.z(), v3.z()});
        double maxx = std::max({v0.x(), v1.x(), v2.x(), v3.x()});
        double maxy = std::max({v0.y(), v1.y(), v2.y(), v3.y()});
        double maxz = std::max({v0.z(), v1.z(), v2.z(), v3.z()});

        vec3 minimum(minx,miny,minz);
        vec3 maximum(maxx,maxy,maxz);
        double eps = 1e-4;
        minimum -= vec3(eps,eps,eps);
        maximum += vec3(eps,eps,eps);

        output_box = AABB(minimum, maximum);
        return true;
    }

private:
    bool hit_triangle(const Ray& r, double t_min, double t_max, HitRecord& rec,
                      const vec3& A, const vec3& B, const vec3& C) const {
        const double EPSILON = 1e-8;
        vec3 edge1 = B - A;
        vec3 edge2 = C - A;
        vec3 h = cross(r.direction(), edge2);
        double a = dot(edge1, h);
        if (fabs(a) < EPSILON) return false;

        double f = 1.0/a;
        vec3 s = r.origin() - A;
        double u = f * dot(s, h);
        if (u < 0.0 || u > 1.0) return false;

        vec3 q = cross(s, edge1);
        double v = f * dot(r.direction(), q);
        if (v < 0.0 || u + v > 1.0) return false;

        double t = f * dot(edge2, q);
        if (t < t_min || t > t_max) return false;

        // Hit
        rec.t = t;
        rec.p = r.at(t);
        vec3 outward_normal = unit_vector(cross(edge1, edge2));
        rec.set_face_normal(r, outward_normal);
        rec.mat_ptr = mat_ptr;

        // Project onto the plane formed by the quad
        double uu = u; 
        double vv = v;
        rec.u = uu; 
        rec.v = vv;
        return true;
    }

    vec3 v0,v1,v2,v3;
    std::shared_ptr<Material> mat_ptr;
};

#endif
