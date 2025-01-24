#ifndef CAMERA_H
#define CAMERA_H

#include "vec3.h"
#include "ray.h"
#include <cmath>

class Camera {
public:
    // Constructor with shutter times
    Camera(vec3 lookFrom, vec3 lookAt, vec3 vup, double vfov, double aspect, double aperture = 0.0, double focus_dist = 1.0, double time0 = 0.0, double time1 = 1.0) {
        double theta = degrees_to_radians(vfov);
        double h = tan(theta / 2.0);
        double viewport_height = 2.0 * h;
        double viewport_width = aspect * viewport_height;

        // Camera coordinate system (u, v, w)
        w = unit_vector(lookFrom - lookAt);
        u = unit_vector(cross(vup, w));
        v = cross(w, u);

        // Aperture and lens radius
        lens_radius = aperture / 2.0;

        // Camera viewport
        origin = lookFrom;
        horizontal = viewport_width * u * focus_dist;
        vertical = viewport_height * v * focus_dist;
        lower_left_corner = origin - horizontal / 2 - vertical / 2 - w * focus_dist;

        // Shutter times
        shutter_open = time0;
        shutter_close = time1;
    }

    // Generate a ray from the camera through the viewport at (s, t)
    Ray get_ray(double s, double t) const {
        vec3 rd = lens_radius * random_in_unit_disk();
        vec3 offset = u * rd.x() + v * rd.y();     

        double time = random_double(shutter_open, shutter_close); 

        return Ray(origin + offset, lower_left_corner + s * horizontal + t * vertical - origin - offset, time);
    }

private:
    vec3 origin;
    vec3 lower_left_corner;
    vec3 horizontal;
    vec3 vertical;
    vec3 u, v, w;
    double lens_radius;

    double shutter_open, shutter_close; // Shutter open and close times

    // Utility function to convert degrees to radians
    static double degrees_to_radians(double degrees) {
        return degrees * M_PI / 180.0;
    }
};

#endif 
