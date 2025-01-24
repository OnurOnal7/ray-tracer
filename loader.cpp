#include "loader.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>

// Loads a 3D model from an OBJ file
#include "loader.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>

bool load_obj(const std::string& filename, 
              std::vector<Vertex>& vertices, 
              std::vector<vec3>& uvs, 
              std::vector<Face>& faces,
              std::vector<vec3>& vertex_normals) {
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Error: Could not open OBJ file " << filename << std::endl;
        return false;
    }

    vertices.clear();
    uvs.clear();
    faces.clear();
    vertex_normals.clear();

    std::vector<vec3> normals; // from vn lines

    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;

        std::istringstream iss(line);
        std::string prefix;
        iss >> prefix;

        if (prefix == "v") {
            double x, y, z;
            iss >> x >> y >> z;
            vertices.emplace_back(Vertex{vec3(x,y,z), vec3()});
        } 
        else if (prefix == "vt") {
            double u,v;
            iss >> u >> v;
            uvs.emplace_back(vec3(u,v,0));
        } 
        else if (prefix == "vn") {
            double nx,ny,nz;
            iss >> nx >> ny >> nz;
            normals.emplace_back(vec3(nx,ny,nz));
        } 
        else if (prefix == "f") {
            std::vector<std::string> tokens;
            std::string token;
            while (iss >> token) tokens.push_back(token);
            if (tokens.size()<3) {
                std::cerr<<"Error: Face with less than 3 vertices."<<std::endl;
                return false;
            }

            Face face;
            face.uv0=face.uv1=face.uv2=-1;
            face.vn0=face.vn1=face.vn2=-1;

            for (int i=0; i<3; i++) {
                std::string vert = tokens[i];
                int v_i=-1, vt_i=-1, vn_i=-1;

                size_t first_slash = vert.find('/');
                if (first_slash==std::string::npos) {
                    // no texture no normal
                    v_i=std::stoi(vert);
                } 
                else {
                    std::string v_part = vert.substr(0,first_slash);
                    v_i=std::stoi(v_part);

                    size_t second_slash = vert.find('/', first_slash+1);
                    if (second_slash==std::string::npos) {
                        // v/vt
                        std::string vt_part = vert.substr(first_slash+1);
                        if(!vt_part.empty()) vt_i=std::stoi(vt_part);
                    } 
                    else {
                        // v/vt/vn or v//vn
                        std::string vt_part = vert.substr(first_slash+1, second_slash - (first_slash+1));
                        std::string vn_part = vert.substr(second_slash+1);

                        if(!vt_part.empty()) vt_i=std::stoi(vt_part);
                        if(!vn_part.empty()) vn_i=std::stoi(vn_part);
                    }
                }

                v_i -=1;
                if(vt_i>=1) vt_i-=1;
                if(vn_i>=1) vn_i-=1;

                if (i==0) {face.v0=v_i;face.uv0=vt_i;face.vn0=vn_i;}
                else if(i==1){face.v1=v_i;face.uv1=vt_i;face.vn1=vn_i;}
                else {face.v2=v_i;face.uv2=vt_i;face.vn2=vn_i;}
            }

            faces.push_back(face);
        }
    }

    std::cout << "Loaded " << vertices.size() << " vertices, "
              << uvs.size() << " UVs, and "
              << faces.size() << " faces." << std::endl;

    bool all_vn = true;
    for (auto &f : faces) {
        if (f.vn0<0 || f.vn1<0 || f.vn2<0) {all_vn=false;break;}
    }

    vertex_normals.resize(vertices.size(),vec3(0,0,0));
    std::vector<int> face_count(vertices.size(),0);

    if(all_vn && !normals.empty()) {
        for (auto &f:faces) {
            vec3 n0=normals[f.vn0];
            vec3 n1=normals[f.vn1];
            vec3 n2=normals[f.vn2];

            vertex_normals[f.v0]+=n0; face_count[f.v0]++;
            vertex_normals[f.v1]+=n1; face_count[f.v1]++;
            vertex_normals[f.v2]+=n2; face_count[f.v2]++;
        }

        for(size_t i=0; i<vertices.size(); i++){
            if(face_count[i]>0)
                vertex_normals[i]=unit_vector(vertex_normals[i]/(double)face_count[i]);
            else
                vertex_normals[i]=vec3(0,1,0);
        }

    } 
    else {
        for (auto& f : faces) {
            vec3 v0 = vertices[f.v0].position;
            vec3 v1 = vertices[f.v1].position;
            vec3 v2 = vertices[f.v2].position;
            vec3 edge1 = v1 - v0;
            vec3 edge2 = v2 - v0;
            vec3 face_normal = unit_vector(cross(edge1, edge2));

            vertex_normals[f.v0] += face_normal; face_count[f.v0]++;
            vertex_normals[f.v1] += face_normal; face_count[f.v1]++;
            vertex_normals[f.v2] += face_normal; face_count[f.v2]++;
        }

        for (size_t i = 0; i < vertices.size(); i++) {
            if (face_count[i] > 0) {
                vertex_normals[i] = unit_vector(vertex_normals[i] / (double)face_count[i]);
            } 
            else {
                vertex_normals[i] = vec3(0,1,0);
            }
        }
    }

    return true;
}

// Loads a texture from a PPM file
std::vector<std::vector<vec3>> load_texture(const std::string& filename, int& tex_width, int& tex_height) {
    std::ifstream ppm_file(filename, std::ios::binary);
    if (!ppm_file) {
        std::cerr << "Error: Could not open texture file " << filename << std::endl;
        exit(1);
    }

    // Read PPM header (P6 format)
    std::string header;
    ppm_file >> header;
    if (header != "P6") {
        std::cerr << "Error: Unsupported texture format." << std::endl;
        exit(1);
    }

    // Read texture dimensions and color range
    ppm_file >> tex_width >> tex_height;
    int max_color;
    ppm_file >> max_color;
    ppm_file.ignore(256, '\n'); 

    // Read texture data
    std::vector<std::vector<vec3>> texture(tex_height, std::vector<vec3>(tex_width));
    for (int y = 0; y < tex_height; ++y) {
        for (int x = 0; x < tex_width; ++x) {
            unsigned char color[3];
            ppm_file.read(reinterpret_cast<char*>(color), 3);
            texture[y][x] = vec3(color[0] / 255.0, color[1] / 255.0, color[2] / 255.0);
        }
    }

    ppm_file.close();
    return texture;
}

// Gets the color from the texture using UV coordinates
vec3 sample_texture(const std::vector<std::vector<vec3>>& texture, int tex_width, int tex_height, double u, double v) {
    // Wrap UV coordinates
    u = u - std::floor(u);
    v = v - std::floor(v);

    // Convert UV to floating-point pixel coordinates
    double x = u * (tex_width - 1);
    double y = (1 - v) * (tex_height - 1); 

    int x0 = static_cast<int>(x);
    int y0 = static_cast<int>(y);
    int x1 = std::min(x0 + 1, tex_width - 1);
    int y1 = std::min(y0 + 1, tex_height - 1);

    double dx = x - x0;
    double dy = y - y0;

    // Bilinear interpolation
    vec3 c00 = texture[y0][x0];
    vec3 c10 = texture[y0][x1];
    vec3 c01 = texture[y1][x0];
    vec3 c11 = texture[y1][x1];

    vec3 c0 = c00 * (1 - dx) + c10 * dx;
    vec3 c1 = c01 * (1 - dx) + c11 * dx;
    vec3 color = c0 * (1 - dy) + c1 * dy;

    return color;
}
