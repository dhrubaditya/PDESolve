#ifndef FILE_cvec3_SEEN
#define FILE_cvec3_SEEN
/*---------------------------------------*/
#include <iostream>
#include <stdio.h>
#include<math.h>
#include "complex.h"
#include "3vec.h"
using namespace std;
/*-------------------------*/
class cvec3{
public:
  complex x,y,z;
  cvec3();
  cvec3(complex,complex,complex);
  // definining operators 
  cvec3 operator+(cvec3);
  cvec3 operator-(cvec3);
  cvec3 operator*(double);
  cvec3 operator*(complex);
  cvec3 operator/(double);
};

cvec3::cvec3(){
  x = complex(0.0,0.0);
  y = complex(0.0,0.0);
  z = complex(0.0,0.0);
}
cvec3::cvec3(complex a, complex b, complex c){
  x =  a;
  y =  b;
  z =  c;
}
cvec3 cvec3::operator+(cvec3 param){
  cvec3 temp;
  temp.x = x+param.x;
  temp.y = y+param.y;
  temp.z = z+param.z;
  return(temp);
}
cvec3 cvec3::operator-(cvec3 param){
  cvec3 temp;
  temp.x = x-param.x;
  temp.y = y-param.y;
  temp.z = z-param.z;
  return(temp);
}
cvec3 cvec3::operator*(double param){
  cvec3 temp;
  temp.x=x*param;
  temp.y=y*param;
  temp.z=z*param;
  return(temp);
}
cvec3 cvec3::operator*(complex param){
  cvec3 temp;
  temp.x=x*param;
  temp.y=y*param;
  temp.z=z*param;
  return(temp);
}
cvec3 cvec3::operator/(double param){
  cvec3 temp;
  temp.x=x/param;
  temp.y=y/param;
  temp.z=z/param;
  return(temp);
}

complex dot(cvec3 a, cvec3 b){
  complex temp;
  temp = a.x*b.x+a.y*b.y+a.z*b.z;
  return(temp);
}
cvec3 cross(vec3 A, cvec3 b){
  cvec3 temp;
  temp.x = b.z*A.y-b.y*A.z;
  temp.y = b.x*A.z-b.z*A.x;
  temp.z = b.y*A.x-b.x*A.y;
  return(temp);
}
cvec3 cross(cvec3 a, vec3 B){
  cvec3 temp;
  temp.x = a.y*B.z-a.z*B.y;
  temp.y = a.z*B.x-a.x*B.z;
  temp.z = a.x*B.y-a.y*B.x;
  return(temp);
}
cvec3 cross(cvec3 a, cvec3 b){
  cvec3 temp;
  temp.x = a.y*b.z-a.z*b.y;
  temp.y = b.x*a.z-a.x*b.z;
  temp.z = a.x*b.y-a.y*b.x;
  return(temp);
}
cvec3 conjg(cvec3 a){
  cvec3 temp;
  temp.x=conjg(a.x);
  temp.y=conjg(a.y);
  temp.z=conjg(a.z);
}
/*---------------------------------------*/

/*---------------------------------------*/
#endif /* !FILE_cvec3_SEEN */
