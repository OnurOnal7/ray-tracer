#ifndef PERLIN_NOISE_TEXTURE_H
#define PERLIN_NOISE_TEXTURE_H

#include "texture.h"
#include "vec3.h"
#include <vector>
#include <cmath>
#include <cstdlib>
#include <algorithm>

// Perlin Noise Texture with extended features: fractal noise, turbulence, and gradient mapping
class PerlinNoiseTexture : public Texture {
public:
    // Constructor with optional scale parameter
    PerlinNoiseTexture(double scale = 1.0) : scale(scale) {
        generate_permutation();
    }

    virtual ~PerlinNoiseTexture() = default;

    // Returns the color at given UV coordinates and position using simple perlin noise (grayscale)
    virtual vec3 value(double u, double v, const vec3& p) const override {
        (void)u;
        (void)v;
        double noise_val = noise(p * scale);
        return vec3(noise_val, noise_val, noise_val);
    }

    // Fractal noise (octave-based Perlin): sums multiple frequencies of noise
    double fractal_noise(const vec3& p, int octaves = 4, double lacunarity = 2.0, double persistence = 0.5) const {
        double total = 0.0;
        double amplitude = 1.0;
        double frequency = 1.0;
        for (int i = 0; i < octaves; i++) {
            total += noise(p * (scale * frequency)) * amplitude;
            amplitude *= persistence;
            frequency *= lacunarity;
        }
        return total;
    }

    // Turbulence: absolute value of fractal noise for a "marbled" look
    double turbulence(const vec3& p, int octaves = 4) const {
        double accum = 0.0;
        double weight = 1.0;
        vec3 temp_p = p * scale;
        for (int i = 0; i < octaves; i++) {
            accum += fabs(noise(temp_p)) * weight;
            temp_p *= 2.0;
            weight *= 0.5;
        }
        return accum;
    }

    double noise(const vec3& pos) const {
        int X = static_cast<int>(floor(pos.x())) & 255;
        int Y = static_cast<int>(floor(pos.y())) & 255;
        int Z = static_cast<int>(floor(pos.z())) & 255;

        double x = pos.x() - floor(pos.x());
        double y = pos.y() - floor(pos.y());
        double z = pos.z() - floor(pos.z());

        double u = fade(x);
        double v = fade(y);
        double w = fade(z);

        int A = perm[X] + Y;
        int AA = perm[A] + Z;
        int AB = perm[A + 1] + Z;
        int B = perm[X + 1] + Y;
        int BA = perm[B] + Z;
        int BB = perm[B + 1] + Z;

        // Add blended results from 8 corners of cube
        double res = lerp(w, 
            lerp(v, 
                lerp(u, grad(perm[AA], x, y, z),
                         grad(perm[BA], x - 1, y, z)),
                lerp(u, grad(perm[AB], x, y - 1, z),
                         grad(perm[BB], x - 1, y - 1, z))),
            lerp(v, 
                lerp(u, grad(perm[AA + 1], x, y, z - 1),
                         grad(perm[BA + 1], x - 1, y, z - 1)),
                lerp(u, grad(perm[AB + 1], x, y - 1, z - 1),
                         grad(perm[BB + 1], x - 1, y - 1, z - 1)))
        );
        return (res + 1.0) / 2.0; 
    }

private:
    double scale;
    std::vector<int> perm; // Permutation vector

    // Generates a random permutation vector
    void generate_permutation() {
        perm.resize(256);
        for (int i = 0; i < 256; ++i) {
            perm[i] = i;
        }
        std::random_device rd;
        std::mt19937 g(rd());
        std::shuffle(perm.begin(), perm.end(), g);
        perm.insert(perm.end(), perm.begin(), perm.end());
    }

    double fade(double t) const {
        return t * t * t * (t * (t * 6 - 15) + 10);
    }

    // Linear interpolation
    double lerp(double t, double a, double b) const {
        return a + t * (b - a);
    }

    double grad(int hash, double x, double y, double z) const {
        int h = hash & 15;
        double u = h < 8 ? x : y;
        double v = h < 4 ? y : (h == 12 || h == 14 ? x : z);
        return ((h & 1) ? -u : u) + ((h & 2) ? -v : v);
    }
};

#endif 
