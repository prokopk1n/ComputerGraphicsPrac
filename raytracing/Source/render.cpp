#include <limits>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include "../Lib/geometry.h"

#define STB_IMAGE_IMPLEMENTATION
#include "../Lib/stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../Lib/stb_image_write.h"

const double POWER = 3.0;
const double Bailout = 2.0;
const int Iterations = 15;

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
                dist = dist_i;
                hit = hit_buf;
                N = N_buf;
                material = material_buf;
            }
        }
    }

    if (fabs(dir.y)>1e-3)  
    {
        dist_i = -(orig.y+4)/dir.y; // the checkerboard plane has equation y = -4
        Vector pt = orig + dir*dist_i;
        if (dist_i>0 && fabs(pt.x)<10 && pt.z<0 && pt.z>-30 && dist_i<dist) {
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
            }
            else 
            {
                Material cyan_plastic(1.0, Vector(0.0, 0.1, 0.06), Vector(0.0, 0.50980392, 0.50980392), Vector(0.50196078, 0.50196078, 0.50196078), 0.25);
                material = cyan_plastic;
                /*material.diffuse = Vector(0.01, 0.01, 0.01);
                material.ambient = Vector(0,0,0);
                material.specular = Vector(0.50, 0.50, 0.50);
                material.refractive_index = 1;
                material.shininess = 0.25;*/
            }
        }
    }

    return dist<1000;
}

Vector reflect(const Vector &I, const Vector &N) // отраженный луч
{
    return I - N*2.f*(I*N);
}

Vector refract(const Vector &I, const Vector &N, float refractive_index) //normalize преломленный луч
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

float fresnel(const Vector &I, const Vector &N, float refractive_index) 
{ 
    float kr;
    float cosi = I * N; 
    float etai = 1, etat = refractive_index; 
    if (cosi > 0) { std::swap(etai, etat); } 
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

Vector multiply(const Vector & vec1, const Vector & vec2)
{
    return Vector(vec1.x*vec2.x, vec1.y*vec2.y, vec1.z*vec2.z);
}

Vector cast_ray(const Vector &orig, const Vector &dir, const std::vector<Sphere> &spheres, std::vector<Light> lights,
    const std::vector<Object> & figures, const Cylinder & cylinder, size_t depth=0) 
{
    Vector point, N;
    Material material;
    if (depth>4 || !scene_intersect(orig, dir, spheres, figures, cylinder, point, N, material)) 
    {
        //return Vector(0, 0, 0);
        return Vector(0.2, 0.7, 0.8); // background color
    }

    Vector reflect_dir = reflect(dir, N).normalize();
    Vector reflect_orig = reflect_dir*N < 0 ? point - N*1e-3 : point + N*1e-3;

    Vector refract_dir = refract(dir, N, material.refractive_index).normalize();
    Vector refract_orig = refract_dir*N < 0 ? point - N*1e-3 : point + N*1e-3;
     // offset the original point to avoid occlusion by the object itself
    
    Vector reflect_color(0,0,0);
    float reflect_coeff = fresnel(dir, N, material.refractive_index);
    if (reflect_coeff >= EPS)
        reflect_color = cast_ray(reflect_orig, reflect_dir, spheres, lights, figures, cylinder,  depth + 1);
    Vector refract_color(0,0,0);
    float refract_coeff = 1 - reflect_coeff;
    if (reflect_coeff>=EPS)
        refract_color = cast_ray(refract_orig, refract_dir, spheres, lights, figures,cylinder, depth + 1);

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
        reflect_color*reflect_coeff + refract_color*refract_coeff;
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
            //Vector buf = cast_ray(Vector(0,0,0), dir, spheres, lights, figures, cylinder);
            Vector buf = cast_ray1(Vector(0,0,0), dir);
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

float DE1(Vector z)
{
	Vector a1 = Vector(1,1,1);
	Vector a2 = Vector(-1,-1,1);
	Vector a3 = Vector(1,-1,-1);
	Vector a4 = Vector(-1,1,-1);
	Vector c;
	int n = 0;
	float dist, d;
	while (n < Iterations) {
		c = a1; 
        dist = (z-a1).norm();
	    d = (z-a2).norm(); 
        if (d < dist) 
        { 
            c = a2; dist=d; 
        }
		d = (z-a3).norm(); 
        if (d < dist) 
        { 
            c = a3; dist=d; 
        }
		d = (z-a4).norm(); 
        if (d < dist) 
        { 
            c = a4; dist=d; 
        }
		z = z*POWER-c*(POWER-1.0);
		n++;
	}

	return z.norm() * pow(POWER, float(-n));
}

float DE(Vector pos) {
	Vector z = pos;
	float dr = 1.0;
	float r = 0.0;
	for (int i = 0; i < Iterations ; i++) {
		r = z.norm();
		if (r>Bailout) break;
		
		// convert to polar coordinates
		float theta = acos(z.z/r);
		float phi = atan2(z.y,z.x);
		dr =  pow( r, POWER-1.0)*POWER*dr + 1.0;
		
		// scale and rotate the point
		float zr = pow( r,POWER);
		theta = theta*POWER;
		phi = phi*POWER;
		
		// convert back to cartesian coordinates
		z = Vector(sin(theta)*cos(phi), sin(phi)*sin(theta), cos(theta))*zr;
		z = z + pos;
	}
	return 0.5*log(r)*r/dr;
}

Vector getNormal(const Vector & ray_pos) 
{
    Vector surf_normal;
    Vector epsilon_x = Vector(EPS, 0, 0);
    Vector epsilon_y = Vector(0, EPS, 0);
    Vector epsilon_z = Vector(0, 0, EPS);

    Vector ray_perb_x1 = ray_pos + epsilon_x;
    Vector ray_perb_y1 = ray_pos + epsilon_y;
    Vector ray_perb_z1 = ray_pos + epsilon_z;

    Vector ray_perb_x2 = ray_pos - epsilon_x ;
    Vector ray_perb_y2 = ray_pos - epsilon_y;
    Vector ray_perb_z2 = ray_pos - epsilon_z;

    surf_normal.x = DE(ray_perb_x1) - DE(ray_perb_x2);
    surf_normal.y = DE(ray_perb_y1) - DE(ray_perb_y2);
    surf_normal.z = DE(ray_perb_z1) - DE(ray_perb_z2);
    return surf_normal.normalize();
}

double MinimumDistance = 0.001;

Vector cast_ray1(const Vector &orig, const Vector &dir)
{
    int MaximumRaySteps = 20;
    float totalDistance = 0.0;
	int steps;
	for (steps=0; steps < MaximumRaySteps; steps++) 
    {
		Vector p = orig + dir * totalDistance;
		float distance = DE1(p);
		totalDistance += distance;
		if (distance < MinimumDistance) break;
	}
    float result = 1.0-float(steps)/float(MaximumRaySteps);
	return Vector(result, result, result);
}


void prepare(std::vector<Sphere> & spheres, std::vector<Light> & lights, std::vector<Object> & figures, Cylinder & cylinder);

int main() {
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
    Material gold(1.0, Vector(0.24275, 0.1995, 0.0745), Vector(0.75164, 0.60648, 0.22648), Vector(0.628281, 0.555802, 0.366065), 0.4);
    Material cyan_plastic(1.0, Vector(0.0, 0.1, 0.06), Vector(0.0, 0.50980392, 0.50980392), Vector(0.50196078, 0.50196078, 0.50196078), 0.25);
    Material glass(1.5, Vector(0, 0, 0), Vector(0.0, 0.0, 0.0), Vector(1.0, 1.0, 1.0), 0.75);
    Material brown(1.0, Vector(1.0f, 0.5f, 0.31f), Vector(1.0f, 0.5f, 0.31f),   Vector(0.5f, 0.5f, 0.5f), 0.25);
    /*Material      glass(1.5, Vector4D(0.0,  0.5, 0.1, 0.8), Vector(0.6, 0.7, 0.8),  125.);
    Material red_rubber(1.0, Vector4D(0.9,  0.1, 0.0, 0.0), Vector(0.3, 0.1, 0.1),   10.);
    Material     mirror(1.0, Vector4D(0.0, 10.0, 0.8, 0.0), Vector(1.0, 1.0, 1.0), 1425.);*/

    //spheres.push_back(Sphere(Vector(0,    0,   -15), 2, brown));
    spheres.push_back(Sphere(Vector( 0, 0, -4), 0.5, glass));
    //spheres.push_back(Sphere(Vector( 1.5, -0.5, -18), 3, red_rubber));
    //spheres.push_back(Sphere(Vector( 7,    5,   -18), 4,     mirror));

    //triangles.push_back(Triangle(red_rubber, Vector(-2, -3, -500), Vector(2, -3, -500), Vector(0, 0, -500)));


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
        vert_list[i] = vert_list[i].matrix(MatrixX).matrix(MatrixY);
    }

    for (int i=0;i<20;i++)
    {
        vert_list[i].z -= 50;
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
        cube_vert[i].z -= 15;
        cube_vert[i].y += 3;
        cube_vert[i].x -= 3;
    }

    cylinder = Cylinder(Vector(0, -1.75, -40.0), -2.0, -1.5, 1, 1, brown);

    figures.push_back(Object(vert_list, triag_list, 108, Vector(0,-0.75,-50.0), glass));
    figures.push_back(Object(cube_vert, my_triag_list, 36, Vector(-3,3,-10), gold));


    lights.push_back(Light(Vector(-10, 10,  -5), Vector(0.5, 0.5, 0.5), Vector(1.0, 1.0, 1.0)));
    lights[0].constant = 1.0f;
    lights[0].linear = 0.014;
    lights[0].quadratic = 0.0007;
    lights.push_back(Light(Vector(0 , 10,  5), Vector(0.5, 0.5, 0.5), Vector(1.0, 1.0, 1.0)));
    lights[1].constant = 1.0f;
    lights[1].linear = 0.014;
    lights[1].quadratic = 0.0007;
    lights.push_back(Light(Vector(30, 30, 1), Vector(0.5, 0.5, 0.5), Vector(1.0, 1.0, 1.0)));
    lights[2].constant = 1.0f;
    lights[2].linear = 0.022f;
    lights[2].quadratic = 0.0019f;
    //lights.push_back(Light(Vector( 30, 20,  30), Vector(0.5, 0.5, 0.5), Vector(1.0, 1.0, 1.0)));
    //lights.push_back(Light(Vector( -20, -5,  30), 1.7));

}
