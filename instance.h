#ifndef INSTANCE_H
#define INSTANCE_H

#include "hittable.h"
#include "vec3.h"
#include <memory>

// Translates a hittable object by a given offset.
class InstanceTranslate : public Hittable {
public:
    InstanceTranslate(std::shared_ptr<Hittable> ptr, const vec3& displacement) : ptr(ptr), offset(displacement) {}

    virtual bool hit(const Ray& r, double t_min, double t_max, HitRecord& rec) const override {
        // Translate ray into local space
        Ray moved_r(r.origin() - offset, r.direction(), r.time());
        if (!ptr->hit(moved_r, t_min, t_max, rec)) return false;

        // Translate hit record back to world space
        rec.p += offset;
        rec.set_face_normal(moved_r, rec.normal);
        return true;
    }

    virtual bool bounding_box(AABB& output_box) const override {
        if (!ptr->bounding_box(output_box)) return false;
        output_box = AABB(output_box.minimum + offset, output_box.maximum + offset);
        return true;
    }

private:
    std::shared_ptr<Hittable> ptr;
    vec3 offset;
};

#endif
