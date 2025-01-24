#ifndef HITTABLE_H
#define HITTABLE_H

#include "ray.h"
#include "hit_record.h"
#include "aabb.h"  
#include <memory>

class Hittable {
public:
    virtual ~Hittable() = default;

    // Pure virtual methods
    virtual bool hit(const Ray& r, double t_min, double t_max, HitRecord& rec) const = 0;
    virtual bool bounding_box(AABB& output_box) const = 0;
};

#endif 
