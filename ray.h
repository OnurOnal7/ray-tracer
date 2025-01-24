#ifndef RAY_H
#define RAY_H

#include "vec3.h"

class Ray {
public:
    // Constructors
    Ray() {}
    Ray(const vec3& origin, const vec3& direction, double time = 0.0) 
        : orig(origin), dir(direction), tm(time) {}
    
    // Accessors
    const vec3& origin() const { return orig; }
    const vec3& direction() const { return dir; }
    double time() const { return tm; }
    
    // Calculate the point at parameter t along the ray
    vec3 at(double t) const { return orig + t * dir; }

private:
    vec3 orig;
    vec3 dir;
    double tm;
};

#endif 
