#include <limits>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include "geometry.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

const int WIDTH = 1024, HEIGHT = 1024;

bool scene_intersect(const Vector &orig, const Vector &dir, const std::vector<Sphere> &spheres, const std::vector<Triangle> triangles,
    const Dodekaedr & dodekaedr, Vector &hit, Vector &N, Material &material) 
{
    float dist_i;
    float dist = std::numeric_limits<float>::max();
    for (size_t i=0; i < spheres.size(); i++) {
        if (spheres[i].ray_intersect(orig, dir, dist_i) && dist_i < dist) {
            dist = dist_i;
            hit = orig + dir*dist_i;
            N = (hit - spheres[i].center).normalize();
            material = spheres[i].material;
        }
    }
    float triangle_dist = std::numeric_limits<float>::max(); 

    for (size_t i=0; i < triangles.size(); i++) {
        if (triangles[i].ray_intersect(orig, dir, dist_i) && dist_i < dist) {
            dist = dist_i;
            hit = orig + dir*dist_i;
            N = triangles[i].N().normalize();
            material = triangles[i].material;
        }
    }

    Material material_buf;
    Vector N_buf;
    Vector hit_buf;

    if(dodekaedr.ray_intersect(orig, dir, dist_i, N_buf, hit_buf, material_buf))
    {
        if (dist_i<dist)
        {
            dist = dist_i;
            //std::cout<<dist<<std::endl;
            hit = hit_buf;
            N = N_buf;
            material = material_buf;
        }
    }

    float checkerboard_dist = std::numeric_limits<float>::max();
    if (fabs(dir.y)>1e-3)  {
        float d = -(orig.y+4)/dir.y; // the checkerboard plane has equation y = -4, d - коэффициент
        Vector pt = orig + dir*d;
        if (d>0 && fabs(pt.x)<10 && pt.z<-10 && pt.z>-30 && d<dist) {
            dist = d;
            hit = pt;
            N = Vector(0,1,0);
            material.diffuse_color = (int(.5*hit.x+1000) + int(.5*hit.z)) & 1 ? Vector(1,1,1) : Vector(1, .7, .3);
            material.diffuse_color = material.diffuse_color*.3;
        }
    }
    return dist<1000;
}

Vector reflect(const Vector &I, const Vector &N) {
    return I - N*2.f*(I*N);
}

Vector refract(const Vector &I, const Vector &N, float refractive_index) //normalize
{
    float cosi = I*N*(-1);
    float etai = 1, etat = refractive_index;
    Vector n = N;
    if (cosi < 0) 
    { // if the ray is inside the object, swap the indices and invert the normal to get the correct result
        cosi = -cosi;
        std::swap(etai, etat); 
        n = N*(-1);
    }
    float eta = etai / etat;
    float k = 1 - eta*eta*(1 - cosi*cosi);//cos угла преломления
    return I*eta + n*(eta * cosi - sqrtf(k));
}

Vector cast_ray(const Vector &orig, const Vector &dir, const std::vector<Sphere> &spheres, std::vector<Light> lights, const std::vector<Triangle> triangles, 
    const Dodekaedr & dodekaedr, size_t depth=0) 
{
    Vector point, N;
    Material material;
    if (depth>4 || !scene_intersect(orig, dir, spheres, triangles, dodekaedr,
     point, N, material)) {
        return Vector(0.2, 0.7, 0.8); // background color
    }

    Vector reflect_dir = reflect(dir, N).normalize();
    Vector refract_dir = refract(dir, N, material.refractive_index).normalize();

    Vector refract_orig = refract_dir*N < 0 ? point - N*1e-3 : point + N*1e-3;
    Vector reflect_orig = reflect_dir*N < 0 ? point - N*1e-3 : point + N*1e-3; // offset the original point to avoid occlusion by the object itself
    
    Vector reflect_color = cast_ray(reflect_orig, reflect_dir, spheres, lights, triangles, dodekaedr, depth + 1);
    Vector refract_color = cast_ray(refract_orig, refract_dir, spheres, lights, triangles, dodekaedr, depth + 1);

    Vector point1, N1;
    Material material1;
    float diffuse_light_intensity = 0.5, specular_light_intensity = 0;
    for (size_t i=0; i<lights.size(); i++) {
        Vector light_dir      = (lights[i].position - point).normalize();

        float light_distance = (lights[i].position - point).norm();
        
        if (scene_intersect(lights[i].position, light_dir*(-1.0), spheres, triangles, dodekaedr, point1, N1, material1)
            && fabs((point1 - lights[i].position).norm() - light_distance)>0.001)
            continue;

        diffuse_light_intensity  += lights[i].intensity * std::max(0.f, light_dir*N);
        specular_light_intensity += powf(std::max(0.f, reflect(light_dir, N)*dir), material.specular_exponent)*lights[i].intensity;
    }
    return material.diffuse_color * diffuse_light_intensity * material.albedo.x + 
        Vector(1., 1., 1.)*specular_light_intensity * material.albedo.y + reflect_color*material.albedo.z 
        + refract_color*material.albedo.d;
}

//начало координат в верхнем левом угле
void render(const std::vector<Sphere> & spheres, const std::vector<Light> lights, const std::vector<Triangle> triangles, const Dodekaedr & dodekaedr) 
{
    const int width    = WIDTH;
    const int height   = HEIGHT;
    const int fov = M_PI/2.;
    std::vector<Vector> framebuffer(width*height);

    Pixel * image = new Pixel[width*height];

    #pragma omp parallel for
    for (size_t j = 0; j<height; j++) {
        for (size_t i = 0; i<width; i++) {
            float x =  (2*(i + 0.5)/(float)width  - 1)*tan(fov/2.);
            float y = -(2*(j + 0.5)/(float)height - 1)*tan(fov/2.);
            Vector dir = Vector(x, y, -1).normalize();
            Vector buf = cast_ray(Vector(0,0,0), dir, spheres, lights, triangles, dodekaedr);
            Vector &c = buf;
            float max = std::max(buf.x, std::max(buf.y, buf.z));
            if (max>1)
            {
                c = c*(1./max);
            }
            framebuffer[i+j*width] = buf;
        }
    }

    #pragma omp parallel for
    for (size_t j = 0; j<height; j++) {
        for (size_t i = 0; i<width; i++) {
            image[i+j*width].r = (uint8_t)(framebuffer[i+j*width].x*255);
            image[i+j*width].g = (uint8_t)(framebuffer[i+j*width].y*255);
            image[i+j*width].b = (uint8_t)(framebuffer[i+j*width].z*255);
        }
    }

    stbi_write_png("./out.png", width, height, 3, image, width*3);
}

int main() {
    Material      ivory(1.0, Vector4D(0.6,  0.3, 0.1, 0.0), Vector(0.4, 0.4, 0.3),   256.);
    Material      glass(1.5, Vector4D(0.0,  0.5, 0.1, 0.8), Vector(0.6, 0.7, 0.8),  125.);
    Material red_rubber(1.0, Vector4D(0.9,  0.1, 0.0, 0.0), Vector(0.3, 0.1, 0.1),   10.);
    Material     mirror(1.0, Vector4D(0.0, 10.0, 0.8, 0.0), Vector(1.0, 1.0, 1.0), 1425.);

    std::vector<Sphere> spheres;
    //spheres.push_back(Sphere(Vector(-5,    0,   -8), 2,      ivory));
    //spheres.push_back(Sphere(Vector(-2.0, -2, -12), 1,      glass));
    //spheres.push_back(Sphere(Vector( 1.5, -0.5, -18), 3, red_rubber));
    //spheres.push_back(Sphere(Vector( 7,    5,   -18), 4,     mirror));

    std::vector<Triangle> triangles;
    //triangles.push_back(Triangle(red_rubber, Vector(-2, -3, -500), Vector(2, -3, -500), Vector(0, 0, -500)));

    std::vector<Pentagon> pentagons;
    //pentagons.push_back(Pentagon(red_rubber, Vector(-1, -1.63, -5), Vector(1.36, -1.63, -5), 
        //Vector(2.09, 0.62, -5), Vector(0.18, 2.01, -5), Vector(-1.73, 0.62, -5)));

    Vector * vert_list = new Vector[20]{
    Vector(0.5, -1.30901694297790528, 0),
    Vector(0.80901700258255008, -0.80901700258255008, 0.80901700258255008),
    Vector(-0.5, -1.30901694297790528, 0),
    Vector(0, -0.5, 1.30901694297790528),
    Vector(-0.80901700258255008, -0.80901700258255008, 0.80901700258255008),
    Vector(-0.80901700258255008, -0.80901700258255008, -0.80901700258255008),
    Vector(0, -0.5, -1.30901694297790528),
    Vector(0.80901700258255008, -0.80901700258255008, -0.80901700258255008),
    Vector(1.30901694297790528, 0, -0.5),
    Vector(1.30901694297790528, 0, 0.5),
    Vector(0, 0.5, -1.30901694297790528),
    Vector(0.80901700258255008, 0.80901700258255008, -0.80901700258255008),
    Vector(0, 0.5, 1.30901694297790528),
    Vector(0.80901700258255008, 0.80901700258255008, 0.80901700258255008),
    Vector(-0.5, 1.30901694297790528, 0),
    Vector(0.5, 1.30901694297790528, 0),
    Vector(-0.80901700258255008, 0.80901700258255008, -0.80901700258255008),
    Vector(-1.30901694297790528, 0, -0.5),
    Vector(-1.30901694297790528, 0, 0.5),
    Vector(-0.80901700258255008, 0.80901700258255008, 0.80901700258255008)
    };

    for (int i=0;i<20;i++)
        vert_list[i].z -= 10;

    for (int i=0;i<20;i++)
        vert_list[i].y -= 2;

    int * triag_list = new int[108]{
    1, 2, 3,
    2, 4, 3,
    4, 5, 3,
    3, 6, 7,
    7, 8, 1,
    1, 3, 7,
    8, 9, 1,
    9, 10, 1,
    10, 2, 1,
    9, 8, 7,
    7, 11, 9,
    11, 12, 9,
    13, 4, 2,
    2, 10, 13,
    10, 14, 13,
    15, 16, 12,
    12, 11, 15,
    11, 17, 15,
    18, 19, 15,
    18, 15, 17,
    19, 20, 15,
    11, 7, 6,
    6, 18, 11,
    18, 17, 11,
    3, 5, 19,
    19, 18, 6,
    6, 3, 19,
    12, 16, 10,
    12, 10, 9,
    16, 14, 10,
    13, 14, 16,
    16, 15, 13,
    15, 20, 13,
    19, 5, 4,
    4, 13, 19,
    13, 20, 19
    };

    Dodekaedr dodekaedr(vert_list, triag_list, Vector(0,-1,-10), ivory);

    std::vector<Light>  lights;
    lights.push_back(Light(Vector(-20, 20,  20), 1.5));
    lights.push_back(Light(Vector( 30, 50, -25), 1.8));
    lights.push_back(Light(Vector( 30, 20,  30), 1.7));
    //lights.push_back(Light(Vector( -20, -5,  30), 1.7));

    render(spheres, lights, triangles, dodekaedr);
    return 0;
}
