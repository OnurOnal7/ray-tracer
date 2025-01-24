#ifndef BVH_H
#define BVH_H

#include "hittable.h"
#include "aabb.h"  
#include "vec3.h"
#include "ray.h"
#include "hit_record.h"
#include <memory>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <iostream>

// Comparator Functions Declaration
bool box_x_compare(const std::shared_ptr<Hittable> a, const std::shared_ptr<Hittable> b);
bool box_y_compare(const std::shared_ptr<Hittable> a, const std::shared_ptr<Hittable> b);
bool box_z_compare(const std::shared_ptr<Hittable> a, const std::shared_ptr<Hittable> b);

// BVHNode Class Declaration
class BVHNode : public Hittable {
public:
    std::shared_ptr<Hittable> left;
    std::shared_ptr<Hittable> right;
    AABB box;

    BVHNode();
    BVHNode(std::vector<std::shared_ptr<Hittable>>& objects, size_t start, size_t end);

    virtual bool hit(const Ray& r, double t_min, double t_max, HitRecord& rec) const override;
    virtual bool bounding_box(AABB& output_box) const override;
    
    // Add a virtual destructor to suppress warnings
    virtual ~BVHNode() {}
};

#endif 
