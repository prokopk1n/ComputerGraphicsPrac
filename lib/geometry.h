#ifndef __GEOMETRY_H__
#define __GEOMETRY_H__
#include <limits>
#include <cmath>
#include <vector>
#include <cassert>
#include <iostream>

const float EPS = 0.0001;

struct Vector
{
    float x, y, z;
    Vector(float x1 = 0,float y1 = 0, float z1 = 0): x(x1), y(y1), z(z1){}
    Vector normalize() const
    {
        float abs = sqrtf(x*x + y*y + z*z);
        return Vector(x/abs, y/abs, z/abs);
    }

    float norm() const
    {
        return sqrtf(x*x+y*y+z*z);
    }

    Vector operator+(const Vector & object) const
    {
        return Vector(x+object.x, y+object.y, z+object.z);
    }

    Vector operator-(const Vector & object) const
    {
        return Vector(x-object.x, y-object.y, z-object.z);
    }

    float operator*(const Vector & object) const
    {
        return x*object.x + y*object.y + z*object.z;
    }

    Vector operator*(float object) const
    {
        return Vector(x*object, y*object, z*object);
    }

    bool operator==(const Vector & object) const
    {
        return fabs(x-object.x)<EPS && fabs(y-object.y)<EPS && fabs(z-object.z)<EPS;
    }

    Vector matrix(float matrix[3][3]) const
    {
        return Vector(x*matrix[0][0]+y*matrix[0][1]+z*matrix[0][2],x*matrix[1][0]+y*matrix[1][1]+z*matrix[1][2],
            x*matrix[2][0]+y*matrix[2][1]+z*matrix[2][2]);
    }

    friend std::ostream & operator<<(std::ostream & out,const Vector & vec)
    {
        out<<vec.x<<" "<<vec.y<<" "<<vec.z;
        return out;
    }
};

struct Pixel
{
  uint8_t r;
  uint8_t g;
  uint8_t b;
};

struct Light
{
    Vector position;
    Vector diffuse;
    Vector specular;

    static Vector ambient;

    float constant;
    float linear;
    float quadratic;

    Light(const Vector &p, const Vector & _diffuse, const Vector & _specular):
        position(p), diffuse(_diffuse), specular(_specular){}

};

Vector Light::ambient = Vector(0.2,0.2,0.2);

Vector normal(const Vector & vec1, const Vector & vec2)
{
    return Vector(vec1.y*vec2.z-vec2.y*vec1.z, vec2.x*vec1.z - vec1.x*vec2.z, vec1.x*vec2.y - vec2.x*vec1.y);
}

enum MaterialType
{
    kDiffuse,
    kReflection,
    kReflectionAndRefraction
};

struct Material {
    Material(float _index, const Vector & _ambient, const Vector & _diffuse, const Vector & _specular, float _shininess, MaterialType _type = kDiffuse):
        ambient(_ambient), diffuse(_diffuse), specular(_specular), shininess(_shininess), refractive_index(_index), type(_type){}
    Material(){}

    MaterialType type;
    float refractive_index; // показатель преломления
    Vector ambient;//фоновое освещение
    Vector diffuse;
    Vector specular;//блики
    float shininess; // константа блеска
};

struct Sphere
{
    Vector center;
    float radius;
    Material material;

    Sphere(const Vector & center1, float radius1, const Material & material1): center(center1), radius(radius1), material(material1) {}

    // orig - точка
    // dir - направление
    // t0 - расстояние
    bool ray_intersect(const Vector & orig, const Vector & dir, float &t0) const
    {
        Vector L = center - orig;
        float tca = L*dir/dir.norm();
        float d2 = L*L - tca*tca;
        if (d2 > radius*radius) return false;

        float t1 = sqrtf(radius*radius-d2);
        t0 = tca - t1;

        float t2 = tca + t1;
        if (t0<0) t0 = t2; // если мы внутри шара
        if (t0<0) return false;
        return true;
    }
};

struct Triangle
{
    Vector x, y, z;
    Material material;
    int coeff = 1;

    Triangle(const Material & material1, const Vector & x1, const Vector & y1, const Vector & z1, int coeff1=1): 
        x(x1), y(y1), z(z1), material(material1), coeff(coeff1){}

    Vector N() const
    {
        Vector vec1 = y - x;
        Vector vec2 = z - x;
        Vector normal_vec = normal(vec1, vec2);
        return normal_vec;
    }

    bool ray_intersect(const Vector & orig, const Vector & dir, float &t0) const
    {
        Vector vec1 = y - z;
        Vector vec2 = x - z;
        Vector normal_vec = normal(vec1, vec2);
        float A = normal_vec.x;
        float B = normal_vec.y;
        float C = normal_vec.z;
        float D = -(A*x.x+B*x.y+C*x.z);

        //std::cout<<A<<","<<B<<","<<C<<","<<D<<std::endl;
        float t = - (D+A*orig.x+B*orig.y+C*orig.z)/(A*dir.x+B*dir.y+C*dir.z);

        if (t<0)
            return false;

        //std::cout<<t<<std::endl;

        float x1 = orig.x + dir.x*t;
        float y1 = orig.y + dir.y*t; 
        float z1 = orig.z + dir.z*t; 

        Vector point(x1,y1,z1);
        Vector point1 = orig + dir*t;
        //std::cout<<x1<<","<<y1<<","<<z1<<std::endl;
        //std::cout<<point1.x<<","<<point1.y<<","<<point1.z<<std::endl;

        if(fabs(normal(point - x, point - y).norm()/2.0+normal(point-y, point - z).norm()/2.0+normal(point - z, point - x).norm()/2.0 
            - normal_vec.norm()/2.0)<EPS)
        {
            t0 = (point - orig).norm();
            return true;
        }

        return false;

    }
};

struct Cylinder
{
    float y1, y2;
    Vector center;
    float a, b;
    Material material;

    Cylinder(const Vector & _center, float _y1, float _y2, float _a, float _b, Material _material): 
        center(_center),y1(_y1), y2(_y2), a(_a), b(_b), material(_material){}

    Cylinder(){}

    Vector N(const Vector & point) const
    {
        return (Vector(point.x, center.y, point.z) - center).normalize();
    }

    bool ray_intersect(const Vector & orig, const Vector & dir, float &t0, Vector & point, Vector & Norm) const
    {
        double _b = 2.0/(a*a) * (orig.x - center.x)  * dir.x + 2.0/(b*b) * (orig.z - center.z) * dir.z;
        double _a = 1.0/(a*a) * dir.x * dir.x + 1.0 / (b * b) * dir.z * dir.z;
        double _c = 1.0/(a*a)*(orig.x - center.x)*(orig.x - center.x) + 1/(b*b)*(orig.z - center.z)*(orig.z - center.z) - 1;
        double D = _b*_b - 4*_a*_c;
        if (fabs(D)<EPS)
            D=0;
        if (D<0)
        {
            return false;
        }
        
        double t1 = (-_b + sqrtf(D))/(2.0*_a);
        double t2 = (-_b - sqrtf(D))/(2.0*_a);

        if (std::max(t1,t2)<0)
            return false;

        double t = t2;
        Vector point1 = orig + dir*t1;
        Vector point2 = orig + dir*t2;
        point = orig + dir*t2;

        if(point1.y - y2 > EPS && point2.y - y2 > EPS || y1 - point1.y > EPS && y1 - point2.y > EPS)
            return false;
        
        if (t2>0 && point.y<y2+EPS && point.y>y1-EPS)//стандартный случай
        {
            Norm = N(point);
            t0 = t2;
            return true;
        }
        
        if (t2<0)
        {
            if (orig.y<y2+EPS && orig.y>y1-EPS)
            {
                point = orig+dir*t1;
                if (point.y<y2+EPS && point.y>y1-EPS)
                {
                    Norm = N(point);
                    t0 = t1;
                    point = orig + dir*t0;
                    return true;
                }
            }

            float t_buf1 = (y2 - orig.y)/dir.y;
            float t_buf2 = (y1 - orig.y)/dir.y;
            
            t0 = std::min(std::max(t_buf1,0.0f), std::max(t_buf2,0.0f));
            if (fabs(t0)<EPS)
                return false;
            point = orig + dir*t0;
            Norm = fabs(t0-t_buf1)<EPS ? Vector(0,1,0) : Vector(0,-1,0);
            return true;
        }

        if (t2>=0)
        {
            if (point.y>y2)
            {
                float t_buf = (y2 - orig.y)/dir.y;
                t0 = t_buf;
                point = orig + dir*t0;
                Norm = Vector(0,1,0);
                return true;
            }
            else
            {
                float t_buf = (y1 - orig.y)/dir.y;
                t0 = t_buf;
                point = orig + dir*t0;
                Norm = Vector(0,-1,0);
                return true;
            }
        }

        return false;




        /*if ((point2.y - y2)>EPS && (point1.y - y2)<EPS)
        {
            Norm = Vector(0,1,0);
            double t_buf;
            Vector vec_buf = point1 - point2;
            t_buf = (y2 - point2.y)/vec_buf.y;
            point = point2 + vec_buf*t_buf;
            t = (point - orig).norm()/dir.norm();
            if (t<0)
            {
                t_buf = (y1 - point2.y)/vec_buf.y;
                point = point2 + vec_buf*t_buf;
                t = (point - orig).norm()/dir.norm();
            }
        }



        if (t2<0)
        {
            point = orig+dir*t1;
            if (point.y<y2+EPS && point.y>y1-EPS)
            {
                Norm = N(point);
                t0 = t1;
                return true;
            }

            if (point.y>y2+EPS)
            {
                Vector vec_buf = point1 - point2;

            }
        }

        if ((point2.y - y2)>EPS && (point1.y - y2)<EPS)
        {
            Norm = Vector(0,1,0);
            double t_buf;
            Vector vec_buf = point1 - point2;
            t_buf = (y2 - point2.y)/vec_buf.y;
            point = point2 + vec_buf*t_buf;
            t = (point - orig).norm()/dir.norm();
            if (t<0)
            {
                t_buf = (y1 - point2.y)/vec_buf.y;
                point = point2 + vec_buf*t_buf;
                t = (point - orig).norm()/dir.norm();
            }
        }
        else if (((point2.y - y1)<EPS && (point1.y - y1)>EPS))
        {
            Norm = Vector(0,-1,0);
            double t_buf;
            Vector vec_buf = point1 - point2;
            t_buf = (y1 - point2.y)/vec_buf.y;
            point = point2 + vec_buf*t_buf;
            t = (point - orig).norm()/dir.norm();
            if (t<0)
            {
                t_buf = (y2 - point2.y)/vec_buf.y;
                point = point2 + vec_buf*t_buf;
                t = (point - orig).norm()/dir.norm();
            }
        }
        else if (point.y < y1-EPS || point.y > y2+EPS)
        {
            return false;
        }
        else 
        {
            Norm = N(point);
        }
        t0 = t;
        return true;*/
    }

};

struct Object
{
    std::vector<Triangle> triangles;
    Vector psevdo_center;

    Object(){}

    Object(Vector * vert_list, int * triag_list, int list_size, const Vector & psevdo_center1, Material material):psevdo_center(psevdo_center1)
    {
        for (int i=0;i<list_size;i+=3)
        {
            triangles.push_back(Triangle(material, vert_list[triag_list[i]-1], 
                vert_list[triag_list[i+1]-1], vert_list[triag_list[i+2]-1]));
            //std::cout<<i<<std::endl;
        }
    }

    bool ray_intersect(const Vector & orig, const Vector & dir, float &t0, Vector & N, Vector & hit, Material & material) const
    {
        float dist = std::numeric_limits<float>::max();
        float dist_i;
        for (size_t i=0; i < triangles.size(); i++) {
            if (triangles[i].ray_intersect(orig, dir, dist_i) && dist_i < dist) {
                dist = dist_i;
                hit = orig + dir*dist_i;
                N = triangles[i].N().normalize();
                if ((hit-psevdo_center)*N<0)
                    N = N*(-1);
                material = triangles[i].material;
            }
        }

        t0 = dist;
        return dist < 1000;
    }    


};

struct Pentagon
{
    Vector point1, point2, point3, point4, point5;
    Material material;
    int coeff = 1;

    Pentagon(const Material & material1, const Vector & x1, const Vector & y1, const Vector & z1, const Vector & e1, const Vector & d1, int coeff1 = 1):
         point1(x1), point2(y1), point3(z1), point4(e1), point5(d1), material(material1), coeff(coeff1){}

    Vector N() const
    {
        Vector vec1 = point2 - point1;
        Vector vec2 = point3 - point1;
        Vector normal_vec = normal(vec1, vec2);
        return normal_vec;
    }

    float S() const
    {
        static float S = 0;
        if (S == 0)
        {
            float a = (point1 - point2).norm();
            S = a*a/4 * sqrtf(25.0+10.0*sqrtf(5.0));
        }
        return S;
    }

    bool ray_intersect(const Vector & orig, const Vector & dir, float &t0) const
    {
        Vector normal_vec = N();
        float A = normal_vec.x;
        float B = normal_vec.y;
        float C = normal_vec.z;
        float D = -(A*point1.x+B*point1.y+C*point1.z);

        //std::cout<<A<<","<<B<<","<<C<<","<<D<<std::endl;
        float t = - (D+A*orig.x+B*orig.y+C*orig.z)/(A*dir.x+B*dir.y+C*dir.z);
        if (t<=0)
            return false;

        //std::cout<<t<<std::endl;

        float x1 = orig.x + dir.x*t;
        float y1 = orig.y + dir.y*t; 
        float z1 = orig.z + dir.z*t; 

        Vector point(x1,y1,z1);
        //Vector point1 = orig + dir*t;
        //std::cout<<x1<<","<<y1<<","<<z1<<std::endl;
        //std::cout<<point1.x<<","<<point1.y<<","<<point1.z<<std::endl;

        if(fabs(normal(point - point1, point - point2).norm() + normal(point-point2, point - point3).norm()
         + normal(point - point3, point - point4).norm() + normal(point - point4, point - point5).norm()
         + normal(point - point5, point - point1).norm() - S()*2.0)<EPS*100)
        {
            t0 = (point - orig).norm();
            return true;
        }
        /*std::cout<<fabs(normal(point - point1, point - point2).norm() + normal(point-point2, point - point3).norm()
         + normal(point - point3, point - point4).norm()+ normal(point - point4, point - point5).norm() 
         + normal(point - point5, point - point1).norm())<<std::endl;*/

        return false;

    }

};


#endif