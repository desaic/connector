#include <cmath>
#ifndef MATH_HPP
#define MATH_HPP
#ifndef real_t
#define real_t double
#endif
#include <cstddef>//for size_t
#include <vector>
struct Vec3{
  Vec3(real_t _x=0.0f, real_t  _y=0.0f, real_t _z=0.0f)
  {x[0]=_x;x[1]=_y;x[2]=_z;}
  Vec3(const Vec3& _a){
    x[0]=_a.x[0];x[1]=_a.x[1];x[2]=_a.x[2];
  }
  real_t get(int i)const{
    return x[i];
  }
  Vec3 operator-(Vec3 &a){
    return Vec3(x[0]-a.x[0],x[1]-a.x[1],x[2]-a.x[2]);
  }
  Vec3 operator-(){
    return Vec3(-x[0],-x[1],-x[2]);
  }
  Vec3 operator+(const Vec3 &a)const{
    return Vec3(x[0]+a.x[0],x[1]+a.x[1],x[2]+a.x[2]);
  }
  Vec3 &operator-=(const Vec3 &a){
    x[0]-=a.x[0];
    x[1]-=a.x[1];
    x[2]-=a.x[2];
    return *this;
  }
  Vec3 &operator+=(const Vec3 &a){
    x[0]+=a.x[0];
    x[1]+=a.x[1];
    x[2]+=a.x[2];
    return *this;
  }

  Vec3 cross(Vec3 & a){
    return Vec3(x[1]*a.x[2]-x[2]*a.x[1],
		-x[0]*a.x[2]+x[2]*a.x[0],
		x[0]*a.x[1]-x[1]*a.x[0]);
  }

  real_t dot(const Vec3 & v)const {
    return x[0]*v.x[0]+ x[1]*v.x[1]+ x[2]*v.x[2];
  }
  real_t L1n(){
    return std::abs(x[0])+std::abs(x[1])+std::abs(x[2]);
  }
  /**@brief
	L2 Norm of a vector
  */
  double norm(){
    return std::sqrt((double)(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]));
  }
  Vec3 operator/(real_t a){
    return Vec3(x[0]/a,x[1]/a,x[2]/a);
  }
  Vec3& operator/=(real_t a){
    x[0]/=a;
    x[1]/=a;
    x[2]/=a;
    return *this;
  }
  Vec3 operator*(real_t a)const{
    return Vec3(x[0]*a,x[1]*a,x[2]*a);
  }
  Vec3& operator*=(real_t a){
    x[0]*=a;
    x[1]*=a;
    x[2]*=a;
    return *this;
  }

  real_t & operator[](const int i) {return x[i];}
  real_t x[3];
};
Vec3 operator*(real_t c, const Vec3 & v);


#endif
