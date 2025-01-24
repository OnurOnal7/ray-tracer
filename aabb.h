#ifndef AABB_H
#define AABB_H

#include "vec3.h"
#include "ray.h"
#include <algorithm>
#include <iostream>

class AABB {
public:
    vec3 minimum;
    vec3 maximum;

    // Constructors
    AABB();
    AABB(const vec3& a, const vec3& b);

    // Hit function
    bool hit(const Ray& r, double t_min, double t_max) const;

    // Static method to compute surrounding box
    static AABB surrounding_box(const AABB& box0, const AABB& box1);
};

#endif
