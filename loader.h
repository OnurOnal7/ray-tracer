#ifndef LOADER_H
#define LOADER_H

#include "vec3.h"
#include "triangle.h"
#include <string>
#include <vector>

// Represents a single vertex in 3D space
struct Vertex {
    vec3 position;
    vec3 uv;
};

// Represents a face of the mesh
struct Face {
    size_t v0, v1, v2;
    int uv0, uv1, uv2;
    int vn0, vn1, vn2;
};

// Loads a 3D model from an OBJ file
bool load_obj(const std::string& filename, 
              std::vector<Vertex>& vertices, 
              std::vector<vec3>& uvs, 
              std::vector<Face>& faces,
              std::vector<vec3>& vertex_normals);

// Loads a texture from a PPM file
std::vector<std::vector<vec3>> load_texture(const std::string& filename, int& tex_width, int& tex_height);

// Gets the color from the texture using UV coordinates
vec3 sample_texture(const std::vector<std::vector<vec3>>& texture, int tex_width, int tex_height, double u, double v);

#endif 
