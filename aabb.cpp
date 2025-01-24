#include "aabb.h"

// Default Constructor
AABB::AABB() 
    : minimum(vec3(INFINITY, INFINITY, INFINITY)), 
      maximum(vec3(-INFINITY, -INFINITY, -INFINITY)) {}

// Parameterized Constructor
AABB::AABB(const vec3& a, const vec3& b) 
    : minimum(a), maximum(b) {}

// Hit Function Implementation
bool AABB::hit(const Ray& r, double t_min, double t_max) const {
    for (int a = 0; a < 3; a++) {
        double invD = 1.0 / r.direction()[a];
        double t0 = (minimum[a] - r.origin()[a]) * invD;
        double t1 = (maximum[a] - r.origin()[a]) * invD;
        if (invD < 0.0)
            std::swap(t0, t1);
        t_min = t0 > t_min ? t0 : t_min;
        t_max = t1 < t_max ? t1 : t_max;
        if (t_max <= t_min)
            return false;
    }
    return true;
}

// Surrounding Box Function Implementation
AABB AABB::surrounding_box(const AABB& box0, const AABB& box1) {
    vec3 small(
        std::min(box0.minimum.x(), box1.minimum.x()),
        std::min(box0.minimum.y(), box1.minimum.y()),
        std::min(box0.minimum.z(), box1.minimum.z())
    );

    vec3 big(
        std::max(box0.maximum.x(), box1.maximum.x()),
        std::max(box0.maximum.y(), box1.maximum.y()),
        std::max(box0.maximum.z(), box1.maximum.z())
    );

    return AABB(small, big);
}
