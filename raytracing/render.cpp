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

bool scene_intersect(const Vector &orig, const Vector &dir, const std::vector<Sphere> &spheres,
    const std::vector<Object> & figures, const Cylinder & cylinder, Vector &hit, Vector &N, Material &material) 
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

    Material material_buf;
    Vector N_buf;
    Vector hit_buf;

    if (cylinder.ray_intersect(orig, dir, dist_i, hit_buf, N_buf) && dist_i < dist)
    {
        //std::cout<<"HERE\n";
        if (dist_i < dist)
        {
        dist = dist_i;
        hit = hit_buf;
        N = N_buf;
        material = cylinder.material;
        }
    }

    for (size_t i=0; i < figures.size(); i++) 
    {
        if(figures[i].ray_intersect(orig, dir, dist_i, N_buf, hit_buf, material_buf))
        {
            if (dist_i<dist)
            {
                //std::cout<<"HERE\n";
                dist = dist_i;
                hit = hit_buf;
                N = N_buf;
                material = material_buf;
            }
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

Vector multiply(const Vector & vec1, const Vector & vec2)
{
    return Vector(vec1.x*vec2.x, vec1.y*vec2.y, vec1.z*vec2.z);
}

Vector cast_ray(const Vector &orig, const Vector &dir, const std::vector<Sphere> &spheres, std::vector<Light> lights,
    const std::vector<Object> & figures, const Cylinder & cylinder,  size_t depth=0) 
{
    Vector point, N;
    Material material;
    if (depth>4 || !scene_intersect(orig, dir, spheres, figures, cylinder, point, N, material)) 
    {
        //return Vector(0, 0, 0);
        return Vector(0.2, 0.7, 0.8); // background color
    }

    Vector reflect_dir = reflect(dir, N).normalize();
    Vector refract_dir = refract(dir, N, material.refractive_index).normalize();

    Vector refract_orig = refract_dir*N < 0 ? point - N*1e-3 : point + N*1e-3;
    Vector reflect_orig = reflect_dir*N < 0 ? point - N*1e-3 : point + N*1e-3; // offset the original point to avoid occlusion by the object itself
    
    Vector reflect_color = cast_ray(reflect_orig, reflect_dir, spheres, lights, figures, cylinder, depth + 1);
    Vector refract_color = cast_ray(refract_orig, refract_dir, spheres, lights, figures, cylinder, depth + 1);

    Vector point1, N1;
    Material material1;
    Vector diffuse(0,0,0), specular(0,0,0);
    Vector ambient = Light::ambient;
    float attenuation;
    for (size_t i=0; i<lights.size(); i++) 
    {
        Vector light_dir      = (lights[i].position - point).normalize();

        float light_distance = (lights[i].position - point).norm();
        
        if (scene_intersect(lights[i].position, light_dir*(-1.0), spheres, figures, cylinder, point1, N1, material1)
            && fabs((point1 - lights[i].position).norm() - light_distance)>0.001)
            continue;

        attenuation = 1.0 / (lights[i].constant + lights[i].linear*light_distance + lights[i].quadratic*light_distance*light_distance);
        diffuse  = diffuse + multiply(lights[i].diffuse, material.diffuse * std::max(0.f, light_dir*N)) * attenuation;
        specular = specular + multiply(lights[i].specular, 
            material.specular * powf(std::max(0.f, reflect(light_dir, N)*dir), material.shininess*128))*attenuation;
        
    }

    return multiply(ambient, material.ambient) + diffuse + specular + 
        reflect_color*material.reflect + refract_color*material.refract;
}

//начало координат в верхнем левом угле
void render(const std::vector<Sphere> & spheres, const std::vector<Light> lights, const std::vector<Triangle> triangles, 
    const std::vector<Object> & figures, const Cylinder & cylinder) 
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
            Vector buf = cast_ray(Vector(0,0,0), dir, spheres, lights, figures, cylinder);
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
    Material gold(1.0, Vector(0.24275, 0.1995, 0.0745), Vector(0.75164, 0.60648, 0.22648), Vector(0.628281, 0.555802, 0.366065), 0.4,
        0.3, 0.0);

    Material glass(1.5, Vector(0, 0, 0), Vector(0.0, 0.0, 0.0), Vector(1.0, 1.0, 1.0), 0.75, 0.1, 0.8);
    Material brown(1.0, Vector(1.0f, 0.5f, 0.31f), Vector(1.0f, 0.5f, 0.31f),   Vector(0.5f, 0.5f, 0.5f), 0.25, 0.3, 0);
    /*Material      glass(1.5, Vector4D(0.0,  0.5, 0.1, 0.8), Vector(0.6, 0.7, 0.8),  125.);
    Material red_rubber(1.0, Vector4D(0.9,  0.1, 0.0, 0.0), Vector(0.3, 0.1, 0.1),   10.);
    Material     mirror(1.0, Vector4D(0.0, 10.0, 0.8, 0.0), Vector(1.0, 1.0, 1.0), 1425.);*/

    std::vector<Sphere> spheres;
    spheres.push_back(Sphere(Vector(5,    0,   -15), 2, brown));
    spheres.push_back(Sphere(Vector( -1.5, 3, -12), 1, gold));
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

    float MatrixY[3][3];
    float alpha = M_PI;
    MatrixY[0][0] = cos(alpha);
    MatrixY[0][1] = 0;
    MatrixY[0][2] = sin(alpha);
    MatrixY[1][0] = 0;
    MatrixY[1][1] = 1;
    MatrixY[1][2] = 0;
    MatrixY[2][0] = -sin(alpha);
    MatrixY[2][1] = 0;
    MatrixY[2][2] = cos(alpha);

    float MatrixX[3][3];
    float beta = M_PI/6;
    MatrixX[0][0] = 1;
    MatrixX[0][1] = 0;
    MatrixX[0][2] = 0;
    MatrixX[1][0] = 0;
    MatrixX[1][1] = cos(beta);
    MatrixX[1][2] = -sin(beta);
    MatrixX[2][0] = 0;
    MatrixX[2][1] = sin(beta);
    MatrixX[2][2] = cos(beta);

    for (int i=0;i<20;i++)
    {
        //std::cout<<vert_list[i].x<<" "<<vert_list[i].y<<" "<<vert_list[i].z<<std::endl;
        vert_list[i] = vert_list[i].matrix(MatrixX).matrix(MatrixY);
        //std::cout<<vert_list[i].x<<" "<<vert_list[i].y<<" "<<vert_list[i].z<<std::endl;
    }

    for (int i=0;i<20;i++)
    {
        vert_list[i].z -= 5;
        vert_list[i].y -= 0.5;
    }
   // Object dodekaedr(vert_list, triag_list, 108, Vector(0,0,-5), gold);
    

   // Dodekaedr dodekaedr(vert_list, triag_list, Vector(0,0,-5), gold);

    Vector * cube_vert = new Vector[8]{
        Vector(1.0, 0.0, 0.0),
        Vector(1.0, 1.0, 0.0),
        Vector(0.0, 1.0, 0.0),
        Vector(0.0, 0.0, 0.0),
        Vector(1.0, 0.0, 1.0),
        Vector(1.0, 1.0, 1.0),
        Vector(0.0, 1.0, 1.0),
        Vector(0.0, 0.0, 1.0)
    };

    for (int i=0;i<8;i++)
    {
        cube_vert[i].z -= 0.5;
        cube_vert[i].x -= 0.5;
        cube_vert[i].y -= 0.5;
        std::cout<<cube_vert[i]<<std::endl;
    }

    int * my_triag_list = new int[36]{
        1,2,3,
        1,3,4,
        2,3,6,
        3,6,7,
        1,2,6,
        1,5,6,
        1,4,5,
        4,5,8,
        3,4,7,
        4,7,8,
        5,6,7,
        5,7,8
    };

    for (int i=0;i<8;i++)
    {
        cube_vert[i].z -= 5;
        cube_vert[i].y -= 1.5;
    }

    Cylinder cylinder(Vector(0, -1.75, -5.0), -2.0, -1.5, 1, 1, gold);

    std::vector<Object> figures;
    figures.push_back(Object(vert_list, triag_list, 108, Vector(0,0,-5), glass));
    //figures.push_back(Object(cube_vert, my_triag_list, 36, Vector(0,0,-2), gold));


    std::vector<Light>  lights;
    lights.push_back(Light(Vector(-20, 20,  -5), Vector(0.5, 0.5, 0.5), Vector(1.0, 1.0, 1.0)));
    lights[0].constant = 1.0f;
    lights[0].linear = 0;
    lights[0].quadratic = 0;
    lights.push_back(Light(Vector(0 , 10,  5), Vector(0.5, 0.5, 0.5), Vector(1.0, 1.0, 1.0)));
    lights[0].constant = 1.0f;
    lights[0].linear = 0;
    lights[0].quadratic = 0;
    lights.push_back(Light(Vector(30, 30, 1), Vector(0.5, 0.5, 0.5), Vector(1.0, 1.0, 1.0)));
    lights[1].constant = 1.0f;
    lights[1].linear = 0.022f;
    lights[1].quadratic = 0.0019f;
    //lights.push_back(Light(Vector( 30, 20,  30), Vector(0.5, 0.5, 0.5), Vector(1.0, 1.0, 1.0)));
    //lights.push_back(Light(Vector( -20, -5,  30), 1.7));

    render(spheres, lights, triangles, figures, cylinder);
    return 0;
}
