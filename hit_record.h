#ifndef HIT_RECORD_H
#define HIT_RECORD_H

#include "ray.h"
#include "vec3.h"
#include <memory>

// Forward declaration of Material
class Material;

// HitRecord structure
struct HitRecord {
    vec3 p;            
    vec3 normal;          
    std::shared_ptr<Material> mat_ptr; 
    double t;                   
    double u;                      
    double v;           
    bool front_face;                       

    // Set the face normal based on ray direction
    void set_face_normal(const Ray& r, const vec3& outward_normal) {
        front_face = dot(r.direction(), outward_normal) < 0;
        normal = front_face ? outward_normal : -outward_normal;
    }
};

#endif 
