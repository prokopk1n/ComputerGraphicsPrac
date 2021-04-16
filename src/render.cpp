#include <limits>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <omp.h>
#include "../lib/geometry.h"

#define STB_IMAGE_IMPLEMENTATION
#include "../lib/stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../lib/stb_image_write.h"

const float Scale = 2.0;
const double POWER = 3.0;
const double Bailout = 2.0;
const int Iterations = 15;
const int MAX_DEPTH = 2;
const float REFLECT_COEFF = 0.8;
double MinimumDistance = 0.003;

int WIDTH, HEIGHT;

bool scene_intersect(const Vector &orig, const Vector &dir, const std::vector<Sphere> &spheres,
    const std::vector<Object> & figures, const Cylinder & cylinder, Vector &hit, Vector &N, Material &material, float & final_dist) 
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
                dist = dist_i;
                hit = hit_buf;
                N = N_buf;
                material = material_buf;
            }
        }
    }

    /*if (fabs(dir.y)>1e-3)  
    {
        dist_i = -(orig.y+4)/dir.y; // the checkerboard plane has equation y = -4
        Vector pt = orig + dir*dist_i;
        if (dist_i>0 && fabs(pt.x)<10 && pt.z<5 && pt.z>-30 && dist_i<dist) {
            dist = dist_i;
            hit = pt;
            N = Vector(0,1,0);
            if ((int(.5*hit.x+1000) + int(.5*hit.z)) & 1)
            {
                material.diffuse = Vector(0.55, 0.55, 0.55);
                material.ambient = Vector(0,0,0);
                material.specular = Vector(0.70, 0.70, 0.70);
                material.refractive_index = 1;
                material.shininess = 0.25;
                material.type = kDiffuse;
            }
            else 
            {
                Material cyan_plastic(1.0, Vector(0.0, 0.50980392, 0.50980392), Vector(0.0, 0.50980392, 0.50980392), Vector(0.50196078, 0.50196078, 0.50196078), 0.25, kDiffuse);
                material = cyan_plastic;
            }
        }
    }*/

    final_dist = dist;

    return dist<1000;
}

Vector reflect(const Vector &I, const Vector &N) // отраженный луч
{
    return I - N*2.f*(I*N);
}

Vector refract(const Vector &I, const Vector &N, float refractive_index) //normalize преломленный луч
{
    float cosi = I*N;
    float etai = 1, etat = refractive_index;
    Vector n = N;
    if (cosi < 0) 
    { // if the ray is inside the object, swap the indices and invert the normal to get the correct result
        cosi = -cosi;
    }
    else 
    {
        std::swap(etai, etat); 
        n = N*(-1);
    }
    float eta = etai / etat;
    float k = 1 - eta*eta*(1 - cosi*cosi);//cos угла преломления
    return k<0 ? 0 : I*eta + n*(eta * cosi - sqrtf(k));
}

float fresnel(const Vector &_I, const Vector &_N, float refractive_index) 
{ 
    float kr;
    Vector I = _I.normalize();
    Vector N = _N.normalize();
    float cosi = I * N; 
    float etai = 1, etat = refractive_index; 
    if (cosi > 0) 
    { 
        std::swap(etai, etat); 
    } 
    // Compute sini using Snell's law
    float sint = etai / etat * sqrtf(std::max(0.f, 1 - cosi * cosi)); 
    // Total internal reflection
    if (sint >= 1) { 
        kr = 1; 
    } 
    else { 
        float cost = sqrtf(std::max(0.f, 1 - sint * sint)); 
        cosi = fabsf(cosi); 
        float Rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost)); 
        float Rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost)); 
        kr = (Rs * Rs + Rp * Rp) / 2; 
    } 

    return kr;
    // As a consequence of the conservation of energy, transmittance is given by:
    // kt = 1 - kr;
} 

float DE1(Vector z)
{
    /*Vector a2 = Vector(1,1,1);
	Vector a1 = Vector(-1,-1,1);
	Vector a3 = Vector(1,-1,-1);
	Vector a4 = Vector(-1,1,-1);*/
    /*Vector a1 = Vector(1,1,1);
	Vector a2 = Vector(-1,1,-1);
	Vector a3 = Vector(-1,-1,1);
	Vector a4 = Vector(1,1,-1);*/
    Vector a1 = Vector(0,-0.25,1);
	Vector a2 = Vector(1,-0.25,0);
	Vector a3 = Vector(0,0.75,0);
	Vector a4 = Vector(-1,-0.25,0);

	Vector c;
	int n = 0;
	float dist, d;
	while (n < Iterations) {
		c = a1; 
        dist = (z-a1).norm();
	    d = (z-a2).norm(); 
        if (d < dist) 
        { 
            c = a2; 
            dist=d; 
        }
		d = (z-a3).norm(); 
        if (d < dist) 
        { 
            c = a3; 
            dist=d; 
        }
		d = (z-a4).norm(); 
        if (d < dist) 
        { 
            c = a4; 
            dist=d; 
        }
		z = z*Scale-c*(Scale-1.0);
		n++;
	}

	return z.norm() * pow(Scale, float(-n));
}

Vector getNormal(const Vector & p, float eps = EPS) 
{
	Vector n;
	n.x = DE1(p + Vector(eps, 0.0, 0.0)) - DE1(p - Vector(eps, 0.0, 0.0));
	n.y = DE1(p + Vector(0.0, eps, 0.0)) - DE1(p - Vector(0.0, eps, 0.0));
	n.z = DE1(p + Vector(0.0, 0.0, eps)) - DE1(p - Vector(0.0, 0.0, eps));
	return n.normalize();
}

bool ray_marching(const Vector & orig, const Vector & dir, Vector & point, Vector & N, Material & material, float & dist)
{
    int MaximumRaySteps = 100;
    float totalDistance = 0.0;
	int steps;
	for (steps=0; steps < MaximumRaySteps; steps++) 
    {
		Vector p = orig + dir * totalDistance;
		//float distance = DE_sphere(p, center, R);
        float distance = DE1(p);
        //std::cout<<distance<<std::endl;
		totalDistance += distance;
        if (totalDistance >= 1000)
        {
            //std::cout<<"HERE\n";
            return false;
        }
		if (distance < MinimumDistance) 
            break;
	}

    //std::cout<<"HERE\n";
    if (steps == MaximumRaySteps)
        false;
    
    material = Material(1.0, Vector(0.75164, 0.60648, 0.22648), Vector(0.5, 0, 0), Vector(0.7, 0.6, 0.6), 0.25);
    //Material(1.0, Vector(0.19225, 0.19225,0.19225), Vector(0.50754, 0.50754, 0.50754), Vector(0.508273, 0.508273, 0.508273), 0.4, kDiffuse); 
    point = orig + dir * totalDistance;
    dist = totalDistance;
    N = getNormal(point);
    return totalDistance < 1000;

}

Vector multiply(const Vector & vec1, const Vector & vec2)
{
    return Vector(vec1.x*vec2.x, vec1.y*vec2.y, vec1.z*vec2.z);
}

Vector cast_ray(const Vector &orig, const Vector &dir, const std::vector<Sphere> &spheres, std::vector<Light> lights,
    const std::vector<Object> & figures, const Cylinder & cylinder, size_t depth=0) 
{
    if (depth > MAX_DEPTH)
        //return Vector(0.2, 0.7, 0.8);
        return Vector(0,0,0);
    Vector point, N, point2, N2;
    Material material, material2;
    float dist1, dist2;
    bool result = false;
    dist2 = std::numeric_limits<float>::max();
    bool result_std = scene_intersect(orig, dir, spheres, figures, cylinder, point, N, material, dist1);
    bool result_ray_marching = ray_marching(orig, dir, point2, N2, material2, dist2);
    if (!result_std && !result_ray_marching)
        //return Vector(0.2, 0.7, 0.8);
        return Vector(0,0,0);

    if (dist2 < dist1)
    {
        material = material2;
        N = N2;
        point = point2;
        result = true;
    }

    Vector reflect_dir = reflect(dir, N).normalize();
    Vector reflect_orig = reflect_dir*N < 0 ? point - N*1e-3 : point + N*1e-3;

    Vector refract_dir = refract(dir, N, material.refractive_index).normalize();
    Vector refract_orig = refract_dir*N < 0 ? point - N*1e-3 : point + N*1e-3;
     // offset the original point to avoid occlusion by the object itself
    
    Vector reflect_color(0,0,0);
    float reflect_coeff = 0;
    Vector refract_color(0,0,0);
    float refract_coeff = 0;

    switch (material.type)
    {
        case kReflection:
            reflect_color = cast_ray(reflect_orig, reflect_dir, spheres, lights, figures, cylinder,  depth + 1) * REFLECT_COEFF;
        break;
        case kReflectionAndRefraction:
            reflect_coeff = fresnel(dir, N, material.refractive_index);
            reflect_color = cast_ray(reflect_orig, reflect_dir, spheres, lights, figures, cylinder,  depth + 1) * reflect_coeff;
            if (reflect_coeff < 1)
            {
                refract_color = cast_ray(refract_orig, refract_dir, spheres, lights, figures,cylinder, depth + 1)*(1-reflect_coeff);
            }
        break;
    }

    Vector point1, N1;
    Material material1;
    Vector diffuse(0,0,0), specular(0,0,0);
    Vector ambient = Light::ambient;
    /*if (result)
        ambient = Vector(1,1,1); */
    float attenuation;
    for (size_t i=0; i<lights.size(); i++) 
    {
        Vector light_dir      = (lights[i].position - point).normalize();

        float light_distance = (lights[i].position - point).norm();
        
        if (scene_intersect(lights[i].position, light_dir*(-1.0), spheres, figures, cylinder, point1, N1, material1, dist1)
            && light_distance - (point1 - lights[i].position).norm() >0.001)
                continue;

        if (ray_marching(lights[i].position, light_dir*(-1.0), point1, N1, material1, dist1) 
            && light_distance - (point1 - lights[i].position).norm() >0.001)
                continue;

        attenuation = 1.0 / (lights[i].constant + lights[i].linear*light_distance + lights[i].quadratic*light_distance*light_distance);
        diffuse  = diffuse + multiply(lights[i].diffuse, material.diffuse * std::max(0.f, light_dir*N)) * attenuation;
        specular = specular + multiply(lights[i].specular, 
            material.specular * powf(std::max(0.f, reflect(light_dir, N)*dir), material.shininess*128))*attenuation;
        
    }

    return multiply(ambient, material.ambient) + diffuse + specular + 
        reflect_color + refract_color;
}

Vector cast_ray1(const Vector &orig, const Vector &dir);

//начало координат в верхнем левом угле
void render(const std::vector<Sphere> & spheres, const std::vector<Light> lights, 
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
            Vector buf = cast_ray(Vector(0,0,8), dir, spheres, lights, figures, cylinder);
            //Vector buf = cast_ray1(Vector(0,0,0), dir);
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

    stbi_write_png("./328_prokopkin_v4v8.png", width, height, 3, image, width*3);
}

void prepare(std::vector<Sphere> & spheres, std::vector<Light> & lights, std::vector<Object> & figures, Cylinder & cylinder);

int main(int argc, char ** argv) {

    assert(argc==1 || argc==3);
    if (argc == 3)
    {
        assert(!strcmp(argv[1], "-w"));
        WIDTH = atoi(argv[2]);
        HEIGHT = atoi(argv[2]);
    }
    else 
    {
        WIDTH = 512;
        HEIGHT = 512;
    }

    std::cout<<"Practical task on Computer graphics \"Witch ball\""<<std::endl;
    std::cout<<"Author: Prokopkin Sergei 328"<<std::endl;
    std::cout<<"Variant: Dodecahedron with fractal inside\n";

    std::vector<Sphere> spheres;
    Cylinder cylinder;
    std::vector<Object> figures;
    std::vector<Light>  lights;

    prepare(spheres, lights, figures, cylinder);
    render (spheres, lights, figures, cylinder);
    return 0;
}

void prepare(std::vector<Sphere> & spheres, std::vector<Light> & lights, std::vector<Object> & figures, Cylinder & cylinder)
{
    Material gold(1.0, Vector(0.24275, 0.1995, 0.0745), Vector(0.75164, 0.60648, 0.22648), Vector(0.628281, 0.555802, 0.366065), 0.4, kReflection);
    Material cyan_plastic(1.0, Vector(0.0, 0.1, 0.06), Vector(0.0, 0.50980392, 0.50980392), Vector(0.50196078, 0.50196078, 0.50196078), 0.25, kReflection);
    Material glass(1.5, Vector(0.1, 0.1, 0.1), Vector(0.1, 0.1, 0.1), Vector(1.0, 1.0, 1.0), 0.75, kReflectionAndRefraction);
    Material brown(1.0, Vector(1.0f, 0.5f, 0.31f), Vector(1.0f, 0.5f, 0.31f),   Vector(0.5f, 0.5f, 0.5f), 0.25, kDiffuse);

    spheres.push_back(Sphere(Vector(5,    5,   -5), 1, gold));


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
    //float alpha = -M_PI/3.2;
    float alpha = -M_PI/2.5;
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
    float beta = M_PI/5.67;
    MatrixX[0][0] = 1;
    MatrixX[0][1] = 0;
    MatrixX[0][2] = 0;
    MatrixX[1][0] = 0;
    MatrixX[1][1] = cos(beta);
    MatrixX[1][2] = -sin(beta);
    MatrixX[2][0] = 0;
    MatrixX[2][1] = sin(beta);
    MatrixX[2][2] = cos(beta);

    float size = 2;
    float MatrixSize[3][3];
    MatrixSize[0][0] = size;
    MatrixSize[0][1] = 0;
    MatrixSize[0][2] = 0;
    MatrixSize[1][0] = 0;
    MatrixSize[1][1] = size;
    MatrixSize[1][2] = 0;
    MatrixSize[2][0] = 0;
    MatrixSize[2][1] = 0;
    MatrixSize[2][2] = size;


    for (int i=0;i<20;i++)
    {
        vert_list[i] = vert_list[i].matrix(MatrixX).matrix(MatrixY).matrix(MatrixSize);
    }

    for (int i=0;i<20;i++)
    {
        vert_list[i].z -= 0;
        vert_list[i].y -= 0.5;
    }

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
    }

    alpha = M_PI/4;
    MatrixY[0][0] = cos(alpha);
    MatrixY[0][1] = 0;
    MatrixY[0][2] = sin(alpha);
    MatrixY[1][0] = 0;
    MatrixY[1][1] = 1;
    MatrixY[1][2] = 0;
    MatrixY[2][0] = -sin(alpha);
    MatrixY[2][1] = 0;
    MatrixY[2][2] = cos(alpha);

    for (int i=0;i<8;i++)
    {
        cube_vert[i] = cube_vert[i].matrix(MatrixSize);
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
        cube_vert[i].y += 5;
        cube_vert[i].x -= 5;
    }

    float min = vert_list[0].y;
    for (int i=1;i<20;i++)
    {
        if (vert_list[i].y < min)
            min = vert_list[i].y;
    }

    cylinder = Cylinder(Vector(0, min-0.5, 0), min-1, min, 2, 2, gold);

    figures.push_back(Object(vert_list, triag_list, 108, Vector(0,0,0), glass));
    figures.push_back(Object(cube_vert, my_triag_list, 36, Vector(-5,5,-5), cyan_plastic));


    lights.push_back(Light(Vector(-9, -7,  3), Vector(0.5, 0.5, 0.5), Vector(1.0, 1.0, 1.0)));
    lights[0].constant = 1.0f;
    lights[0].linear = 0.014;
    lights[0].quadratic = 0.0007;
    lights.push_back(Light(Vector(8 , -7,  2), Vector(0.5, 0.5, 0.5), Vector(1.0, 1.0, 1.0)));
    lights[1].constant = 1.0f;
    lights[1].linear = 0.014;
    lights[1].quadratic = 0.0007;
    lights.push_back(Light(Vector(0, 10, 5), Vector(0.5, 0.5, 0.5), Vector(1.0, 1.0, 1.0)));
    lights[2].constant = 1.0f;
    lights[2].linear = 0.022f;
    lights[2].quadratic = 0.0019f;
    //lights.push_back(Light(Vector( 30, 20,  30), Vector(0.5, 0.5, 0.5), Vector(1.0, 1.0, 1.0)));
    //lights.push_back(Light(Vector( -20, -5,  30), 1.7));

}
