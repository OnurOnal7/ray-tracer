#include "camera.h"
#include "ray.h"
#include "sphere.h"
#include "triangle.h"
#include "quad.h"
#include "vec3.h"
#include "hittable.h"
#include "aabb.h"   
#include "bvh.h"
#include "loader.h"
#include "lambertian.h"
#include "metal.h"
#include "dielectric.h"
#include "emissive.h"
#include "solid_color.h"     
#include "image_texture.h" 
#include "perlin_noise_texture.h"
#include "cube_map_texture.h"
#include "blend_texture.h"
#include "instance.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <memory>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <future>  
#include <thread>
#include <random>
#include <atomic>

// Global PerlinNoiseTexture instance
static PerlinNoiseTexture global_perlin(2.0);

inline double star_brightness_modulation(const vec3& star_pos) {
    double t = global_perlin.turbulence(star_pos * 0.5, 2);
    return 0.8 + 0.4 * t;
}
inline double shooting_star_variation(const vec3& pos) {
    double val = global_perlin.noise(pos * 0.5);
    return 1.0 + 0.5 * val; 
}
inline double trail_enhancement(const vec3& pos) {
    double detail = global_perlin.fractal_noise(pos, 3, 2.0, 0.5);
    return 1.0 + 0.3 * detail;
}
inline double clamp(double x, double mn, double mx) {
    if (x < mn) return mn;
    if (x > mx) return mx;
    return x;
}

static std::shared_ptr<CubeMapTexture> environment_map = nullptr;

vec3 ray_color(const Ray& r, const std::shared_ptr<Hittable>& world, int depth) {
    if (depth <= 0) return vec3(0, 0, 0);

    HitRecord rec;
    if (!world->hit(r, 0.001, INFINITY, rec)) {
        if (environment_map) {
            return environment_map->sample_direction(r.direction());
        } 
        else {
            return vec3(0.02, 0.03, 0.10);
        }
    }

    vec3 emitted=rec.mat_ptr->emit(rec.u,rec.v,rec.p);
    vec3 attenuation;
    Ray scattered;

    if (!rec.mat_ptr->sample(r, rec, attenuation, scattered)){
        return emitted;
    }

    double pdf = rec.mat_ptr->scatter_pdf(r, rec, scattered);
    if (pdf < 1e-8){
        vec3 old_atten;
        Ray old_scattered;
        if (!rec.mat_ptr->scatter(r, rec, old_atten, old_scattered))
            return emitted;
        return emitted + old_atten * ray_color(old_scattered, world, depth-1);
    }

    double cosine = dot(rec.normal, unit_vector(scattered.direction()));
    if (cosine < 0) cosine = 0;

    return emitted + attenuation * (cosine / pdf) * ray_color(scattered, world, depth - 1);
}

vec3 compute_pixel_color(int i, int j, int image_width, int image_height, const Camera& cam, const std::shared_ptr<Hittable>& world, int max_depth, int max_samples, double variance_threshold){
    vec3 sum(0, 0, 0), sum_sq(0, 0, 0);
    int samples = 0;
    int min_samples = 150, check_interval = 10;

    for (int s = 1; s <= max_samples; s++) {
        double u = (i + (double)rand() / RAND_MAX) / (image_width - 1);
        double v = (j + (double)rand() / RAND_MAX) / (image_height - 1);
        Ray r = cam.get_ray(u, v);
        vec3 c = ray_color(r, world, max_depth);
        sum += c; sum_sq += c * c; samples++;

        if (s >= min_samples && s % check_interval == 0) {
            vec3 mean = sum / (double)samples;
            vec3 mean_sq = sum_sq / (double)samples;
            vec3 var = mean_sq - mean * mean;
            double var_mag = (var.x() + var.y() + var.z()) / 3.0;
            double stddev = sqrt(var_mag);
            if (stddev < variance_threshold)
                break;
        }
    }
    return sum / (double)samples;
}

int main(){
    srand(7);

    // Load environment map
    int env_w, env_h;
    auto env_data = load_texture("space.ppm", env_w,env_h);
    environment_map = std::make_shared<CubeMapTexture>(env_data, env_w, env_h);

    // Ground: blend moon.ppm + perlin
    int moon_w, moon_h;
    auto moon_data = load_texture("moon.ppm", moon_w,moon_h);
    auto moon_tex = std::make_shared<ImageTexture>(moon_data, moon_w,moon_h);
    auto perlin_tex = std::make_shared<PerlinNoiseTexture>(5.0);
    auto blended_ground = std::make_shared<BlendTexture>(moon_tex, perlin_tex, 0.3);
    auto ground_mat = std::make_shared<Lambertian>(blended_ground);

    std::vector<std::shared_ptr<Hittable>> objects;
    objects.emplace_back(std::make_shared<Sphere>(vec3(0, -1002.5, -1), 1000.0, ground_mat));

    // Load satellite
    std::vector<Vertex> vertices;
    std::vector<vec3> uvs;
    std::vector<Face> faces;
    std::vector<vec3> vertex_normals;
    if (!load_obj("satellite.obj", vertices, uvs, faces, vertex_normals)){
        std::cerr << "Failed to load satellite.obj\n";
        return 1;
    }

    int ship_w, ship_h;
    auto ship_data = load_texture("spaceship.ppm", ship_w,ship_h);
    auto ship_tex = std::make_shared<ImageTexture>(ship_data, ship_w, ship_h);
    auto satellite_mat = std::make_shared<Metal>(ship_tex, 0.1);

    std::vector<std::shared_ptr<Hittable>> sat_tris;
    for (auto &f : faces) {
        vec3 v0 = vertices[f.v0].position * 0.1;
        vec3 v1 = vertices[f.v1].position * 0.1;
        vec3 v2 = vertices[f.v2].position * 0.1;

        vec3 uv0(0, 0, 0), uv1(0, 0, 0), uv2(0, 0, 0);
        if (f.uv0 >= 0 && (size_t)f.uv0 < uvs.size()) uv0 = uvs[f.uv0];
        if (f.uv1 >= 0 && (size_t)f.uv1 < uvs.size()) uv1 = uvs[f.uv1];
        if (f.uv2 >= 0 && (size_t)f.uv2 < uvs.size()) uv2 = uvs[f.uv2];

        vec3 n0, n1, n2;
        if (f.vn0 >= 0 && f.vn1 >= 0 && f.vn2 >= 0) {
            n0 = vertex_normals[f.v0];
            n1 = vertex_normals[f.v1];
            n2 = vertex_normals[f.v2];
        } 
        else {
            // fallback: face normal
            vec3 e1 = v1 - v0; vec3 e2 = v2 - v0;
            vec3 fn = unit_vector(cross(e1, e2));
            n0 = n1 = n2 = fn;
        }

        sat_tris.emplace_back(std::make_shared<Triangle>(v0, v1, v2, uv0, uv1, uv2, n0, n1, n2, satellite_mat));
    }

    auto satellite_bvh = std::make_shared<BVHNode>(sat_tris, 0, sat_tris.size());
    objects.emplace_back(std::make_shared<InstanceTranslate>(satellite_bvh, vec3(1.0, -2.45, -4.0)));    

    // Simple white glass texture
    auto glass_texture = std::make_shared<SolidColor>(vec3(1.0, 1.0, 1.0));

    // Blue frost color
    auto blue_frost = std::make_shared<SolidColor>(vec3(0.4, 0.8, 1.0));

    // Blend glass with blue frost using a fixed blend factor
    auto frosted_glass = std::make_shared<BlendTexture>(glass_texture, blue_frost, 0.5);

    auto frosted_glass_material = std::make_shared<Lambertian>(frosted_glass);

    // Define box parameters
    double half = 1.5;
    double bx = -4.5, by = -2.45, bz = -8.0; 
    vec3 bmin(bx - half, by - half, bz - half);
    vec3 bmax(bx + half, by + half, bz + half);

    objects.emplace_back(std::make_shared<Quad>(vec3(bmin.x(), bmin.y(), bmax.z()), vec3(bmax.x(), bmin.y(), bmax.z()), vec3(bmax.x(), bmax.y(), bmax.z()), vec3(bmin.x(), bmax.y(), bmax.z()), frosted_glass_material));
    objects.emplace_back(std::make_shared<Quad>(vec3(bmin.x(), bmin.y(), bmin.z()), vec3(bmin.x(), bmax.y(), bmin.z()), vec3(bmax.x(), bmax.y(), bmin.z()), vec3(bmax.x(), bmin.y(), bmin.z()), frosted_glass_material));
    objects.emplace_back(std::make_shared<Quad>(vec3(bmin.x(), bmin.y(), bmin.z()), vec3(bmin.x(), bmin.y(), bmax.z()), vec3(bmin.x(), bmax.y(), bmax.z()), vec3(bmin.x(), bmax.y(), bmin.z()), frosted_glass_material));
    objects.emplace_back(std::make_shared<Quad>(vec3(bmax.x(), bmin.y(), bmin.z()), vec3(bmax.x(), bmax.y(), bmin.z()), vec3(bmax.x(), bmax.y(), bmax.z()), vec3(bmax.x(), bmin.y(), bmax.z()), frosted_glass_material));
    objects.emplace_back(std::make_shared<Quad>(vec3(bmin.x(), bmin.y(), bmin.z()), vec3(bmax.x(), bmin.y(), bmin.z()), vec3(bmax.x(), bmin.y(), bmax.z()), vec3(bmin.x(), bmin.y(), bmax.z()), frosted_glass_material));
    objects.emplace_back(std::make_shared<Quad>(vec3(bmin.x(), bmax.y(), bmin.z()), vec3(bmin.x(), bmax.y(), bmax.z()), vec3(bmax.x(), bmax.y(), bmax.z()), vec3(bmax.x(), bmax.y(), bmin.z()), frosted_glass_material));

    // Sun
    int sun_w, sun_h;
    auto sun_data = load_texture("sun.ppm", sun_w,sun_h);
    auto sun_tex = std::make_shared<ImageTexture>(sun_data, sun_w,sun_h);
    class TexturedEmissive : public Material {
    public:
        TexturedEmissive(std::shared_ptr<Texture> t) : tex(t){}
        virtual bool scatter(const Ray&, const HitRecord&, vec3&, Ray&) const override {return false;}
        virtual vec3 emit(double u, double v, const vec3&p) const override {
            vec3 c = tex->value(u, v, p);
            return c * 10.0;
        }
    private:
        std::shared_ptr<Texture> tex;
    };
    auto sun_mat = std::make_shared<TexturedEmissive>(sun_tex);
    objects.emplace_back(std::make_shared<Sphere>(vec3(8, 3.5, -10), 1.0, sun_mat));

    // Planets 
    struct PData {
        const char* name;
        vec3 pos;
        double r;
    } planets[]={
        {"mercury", vec3(5.2, 1.5, -8), 0.4},
        {"venus",   vec3(4.5, 2.5, -9.5), 0.5},
        {"earth",   vec3(3.5, 3.5, -9), 0.55},
        {"mars",    vec3(2, 4.5, -10), 0.45},
        {"jupiter", vec3(0.5, 4.25, -11), 0.7},
        {"saturn",  vec3(-2, 4, -12), 0.65},
        {"uranus",  vec3(-3.5, 2, -13), 0.5},
        {"neptune", vec3(-6, 1, -14), 0.6}
    };

    for (auto &pd : planets) {
        int pw, ph;
        std::string fname = std::string(pd.name) + ".ppm";
        auto p_data = load_texture(fname, pw, ph);
        auto p_tex = std::make_shared<ImageTexture>(p_data, pw, ph);
        auto p_mat = std::make_shared<Lambertian>(p_tex);
        objects.emplace_back(std::make_shared<Sphere>(pd.pos, pd.r, p_mat));
    }

    // Emissive spheres for shooting stars
    auto shooting_star_beam_red = std::make_shared<Emissive>(vec3(5.0, 0, 0));   
    auto shooting_star_beam_cyan = std::make_shared<Emissive>(vec3(0, 5, 5));  
    auto shooting_star_beam_purple = std::make_shared<Emissive>(vec3(5, 0, 5));
    auto shooting_star_meteor_red = std::make_shared<Emissive>(vec3(50, 0, 0));    
    auto shooting_star_meteor_cyan = std::make_shared<Emissive>(vec3(0, 50, 50));  
    auto shooting_star_meteor_purple = std::make_shared<Emissive>(vec3(50, 0, 50));   

    objects.emplace_back(std::make_shared<MovingSphere>(vec3(-5, 2, -3), vec3(-3, 2, -5), 0, 1, 0.02, shooting_star_beam_red));
    objects.emplace_back(std::make_shared<MovingSphere>(vec3(-3, 2, -5), vec3(-3, 2, -5.1), 0, 1, 0.03, shooting_star_meteor_red));
    objects.emplace_back(std::make_shared<MovingSphere>(vec3(5.5, 4, -10), vec3(8.5, 4.5, -15), 0, 1, 0.02, shooting_star_beam_cyan));
    objects.emplace_back(std::make_shared<MovingSphere>(vec3(8.5, 4.5, -15), vec3(8.5, 4.5, -15.1), 0, 1, 0.03, shooting_star_meteor_cyan));
    objects.emplace_back(std::make_shared<MovingSphere>(vec3(-4, 1.5, -7), vec3(-2, 2.5, -9), 0, 1, 0.02, shooting_star_beam_purple));
    objects.emplace_back(std::make_shared<MovingSphere>(vec3(-4, 1.5, -7), vec3(-4, 1.5, -6.9), 0, 1, 0.03, shooting_star_meteor_purple));

    // White star material and random stars
    auto star_material = std::make_shared<Emissive>(vec3(1, 1, 1));
    std::vector<double> star_radii = {0.005, 0.01, 0.015, 0.02};
    std::vector<std::shared_ptr<Hittable>> base_stars;
    for (auto r : star_radii) {
        base_stars.push_back(std::make_shared<Sphere>(vec3(0, 0, 0), r, star_material));
    }

    int num_stars = 1000;
    for (int i = 0; i < num_stars; i++) {
        double x = random_double(-15, 15);
        double y = random_double(0, 10);
        double z = random_double(-15, -10);
        int idx = rand() % 4;
        objects.emplace_back(std::make_shared<InstanceTranslate>(base_stars[idx], vec3(x, y, z)));
    }

    std::shared_ptr<Hittable> world = std::make_shared<BVHNode>(objects, 0, objects.size());

    vec3 lookFrom(0, 0, 2);  
    vec3 lookAt(0, 0, -1);   
    vec3 vup(0, 1, 0);     
    double vfov = 50.0;   
    double aspect_ratio = 16.0 / 9.0; 
    double aperture = 0.005;   
    double focus_dist = (lookFrom - lookAt).length();
    double shutter_open = 0.0;
    double shutter_close = 1.0;
    Camera cam(lookFrom, lookAt, vup, vfov, aspect_ratio, aperture, focus_dist, shutter_open, shutter_close);

    int image_width_final = 7680;
    int image_height_final = (int)(image_width_final / aspect_ratio);
    int max_samples = 3000;
    int max_depth = 150;
    double variance_threshold = 0.0003;

    std::vector<float> framebuffer(image_width_final * image_height_final * 3);
    std::atomic<int> scanlines_remaining(image_height_final);

    auto render_task = [&](int start_row, int end_row) {
        for (int j = start_row; j < end_row; j++) {
            for (int i = 0; i < image_width_final; i++) {
                vec3 pixel_color = compute_pixel_color(i, j, image_width_final, image_height_final, cam, world, max_depth, max_samples, variance_threshold);
                
                // Clamp pixel_color to [0,1]
                pixel_color = vec3(
                    clamp(pixel_color.x(), 0.0, 1.0),
                    clamp(pixel_color.y(), 0.0, 1.0),
                    clamp(pixel_color.z(), 0.0, 1.0)
                );
                
                int idx = (j * image_width_final + i) * 3;
                framebuffer[idx] = pixel_color.x();
                framebuffer[idx+1] = pixel_color.y();
                framebuffer[idx+2] = pixel_color.z();
            }
            --scanlines_remaining;
            std::clog << "\rScanlines remaining: " << scanlines_remaining.load() << ' ' << std::flush;
        }
    };

    unsigned int num_threads = std::thread::hardware_concurrency();
    if (num_threads == 0) num_threads = 4;
    int rows_per_task = (image_height_final + (int)num_threads - 1) / (int)num_threads;
    std::vector<std::future<void>> futures;    

    for (unsigned int t = 0; t < num_threads; t++) {
        int start_row = t * rows_per_task;
        if (start_row >= image_height_final) break;
        int end_row = std::min(start_row + rows_per_task, image_height_final);
        futures.push_back(std::async(std::launch::async, render_task, start_row, end_row));
    }

    for (auto &f : futures) f.get();

    std::ofstream pfm_file("output.pfm", std::ios::out|std::ios::binary);
    pfm_file << "PF\n" << image_width_final << " " << image_height_final << "\n-1.0\n";
    pfm_file.write(reinterpret_cast<char*>(framebuffer.data()), framebuffer.size() * sizeof(float));
    pfm_file.close();
    std::cout << "\nPFM image written to output.pfm\n";

    return 0;
}
