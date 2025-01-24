#ifndef VEC3_H
#define VEC3_H

#include <cmath>
#include <iostream>
#include <random>

// Vector class
class vec3 {
public:
    double e[3];  // x, y, z components

    // Constructors
    vec3() { e[0] = e[1] = e[2] = 0.0; }
    vec3(double e0, double e1, double e2) {
        e[0] = e0;
        e[1] = e1;
        e[2] = e2;
    }

    // Accessors
    double x() const { return e[0]; }
    double y() const { return e[1]; }
    double z() const { return e[2]; }

    // Negation operator
    vec3 operator-() const {
        return vec3(-e[0], -e[1], -e[2]);
    }

    // Indexing operator
    double operator[](int i) const { return e[i]; }
    double& operator[](int i) { return e[i]; }

    // Vector addition
    vec3& operator+=(const vec3& v) {
        e[0] += v.e[0];
        e[1] += v.e[1];
        e[2] += v.e[2];
        return *this;
    }

    // Vector subtraction
    vec3& operator-=(const vec3& v) {
        e[0] -= v.e[0];
        e[1] -= v.e[1];
        e[2] -= v.e[2];
        return *this;
    }

    // Scalar multiplication
    vec3& operator*=(double t) {
        e[0] *= t;
        e[1] *= t;
        e[2] *= t;
        return *this;
    }

    // Scalar division
    vec3& operator/=(double t) {
        return *this *= 1 / t;
    }

    // Vector length (magnitude)
    double length() const {
        return std::sqrt(e[0] * e[0] + e[1] * e[1] + e[2] * e[2]);
    }

    // Vector length squared
    double length_squared() const {
        return e[0] * e[0] + e[1] * e[1] + e[2] * e[2];
    }

    // Normalize the vector
    vec3 normalize() const {
        double len = length();
        return vec3(e[0] / len, e[1] / len, e[2] / len);
    }

    // Check if vector is near zero in all dimensions
    bool near_zero() const {
        const double s = 1e-8;
        return (fabs(e[0]) < s) && (fabs(e[1]) < s) && (fabs(e[2]) < s);
    }
};

/* Operator Overloads */

// Vector addition
inline vec3 operator+(const vec3& u, const vec3& v) {
    return vec3(u.e[0] + v.e[0], u.e[1] + v.e[1], u.e[2] + v.e[2]);
}

// Vector subtraction
inline vec3 operator-(const vec3& u, const vec3& v) {
    return vec3(u.e[0] - v.e[0], u.e[1] - v.e[1], u.e[2] - v.e[2]);
}

// Scalar multiplication
inline vec3 operator*(double t, const vec3& v) {
    return vec3(t * v.e[0], t * v.e[1], t * v.e[2]);
}

// Scalar multiplication (vec3 * double)
inline vec3 operator*(const vec3& v, double t) {
    return vec3(v.e[0] * t, v.e[1] * t, v.e[2] * t);
}

// Element-wise vector multiplication
inline vec3 operator*(const vec3& u, const vec3& v) {
    return vec3(u.e[0] * v.e[0], u.e[1] * v.e[1], u.e[2] * v.e[2]);
}

// Scalar division
inline vec3 operator/(const vec3& v, double t) {
    return vec3(v.e[0] / t, v.e[1] / t, v.e[2] / t);
}

// Element-wise vector division
inline vec3 operator/(const vec3& u, const vec3& v) {
    return vec3(u.e[0] / v.e[0], u.e[1] / v.e[1], u.e[2] / v.e[2]);
}

/* Vector Utility Functions */

// Dot product
inline double dot(const vec3& u, const vec3& v) {
    return u.e[0] * v.e[0] + u.e[1] * v.e[1] + u.e[2] * v.e[2];
}

// Cross product
inline vec3 cross(const vec3& u, const vec3& v) {
    return vec3(
        u.e[1] * v.e[2] - u.e[2] * v.e[1],
        u.e[2] * v.e[0] - u.e[0] * v.e[2],
        u.e[0] * v.e[1] - u.e[1] * v.e[0]
    );
}

// Normalize a vector
inline vec3 unit_vector(const vec3& v) {
    return v / v.length();
}

// Generate a random double in [0,1)
inline double random_double() {
    static std::uniform_real_distribution<double> distribution(0.0, 1.0);
    static std::mt19937 generator;
    return distribution(generator);
}

// Utility function to generate a random double in a given range
inline double random_double(double min, double max) {
    return min + (max - min) * random_double();
}

// Generate a random point inside a unit sphere
inline vec3 random_in_unit_sphere() {
    while (true) {
        vec3 p(random_double() * 2.0 - 1.0,
               random_double() * 2.0 - 1.0,
               random_double() * 2.0 - 1.0);
        if (p.length_squared() >= 1.0) continue;
        return p;
    }
}

// Generate a random point inside a unit disk
inline vec3 random_in_unit_disk() {
    while (true) {
        vec3 p = vec3(random_double(-1,1), random_double(-1,1), 0);
        if (p.length_squared() >= 1) continue;
        return p;
    }
}

// Generate a random unit vector
inline vec3 random_unit_vector() {
    return unit_vector(random_in_unit_sphere());
}

// Reflect the vector around the normal
inline vec3 reflect(const vec3& v, const vec3& n) {
    return v - 2 * dot(v, n) * n;
}

// Refract the vector based on the normal and ratio of indices
inline vec3 refract(const vec3& uv, const vec3& n, double etai_over_etat) {
    double cos_theta = fmin(dot(-uv, n), 1.0);
    vec3 r_out_perp = etai_over_etat * (uv + cos_theta * n);
    vec3 r_out_parallel = -sqrt(fabs(1.0 - r_out_perp.length_squared())) * n;
    return r_out_perp + r_out_parallel;
}

#endif
