#ifndef FILE_complex_SEEN
#define FILE_complex_SEEN
/*---------------------------------------*/
#include <iostream>
#include<math.h>
using namespace std;
/*-------------------------*/
class complex{
public:
  double RE,IM;
  complex();
  complex(int,int);
  complex(float,float);
  complex(double,double);
    // definining operators
  complex operator+(complex);
  complex operator-(complex);
  complex operator*(double);
  complex operator*(complex);
  complex operator/(double);
};

complex::complex(){
  RE = 0.0;
  IM = 0.0;
}
complex::complex(int a, int b){
  RE = double(a);
  IM = double(b);
}
complex::complex(float a, float b){
  RE = double(a);
  IM = double(b);
}
complex::complex(double a, double b){
  RE =  a;
  IM =  b;
}

complex complex::operator+(complex param){
  complex temp;
  temp.RE = RE+param.RE;
  temp.IM = IM+param.IM;
  return(temp);
}
complex complex::operator-(complex param){
  complex temp;
  temp.RE = RE-param.RE;
  temp.IM = IM-param.IM;
  return(temp);
}
complex complex::operator*(double param){
  complex temp;
  temp.RE=param*RE;
  temp.IM=param*IM;
  return(temp);
}
complex complex::operator/(double param){
  complex temp;
  temp.RE=RE/param;
  temp.IM=IM/param;
  return(temp);
}
complex complex::operator*(complex param){
  complex temp;
  temp.RE=param.RE*RE-param.IM*IM;
  temp.IM=param.RE*IM+param.IM*RE;
  return(temp);
}
complex conjg(complex a){
  complex tmp;
  tmp.RE=a.RE;
  tmp.IM=-a.IM;
  return(tmp);
}
double norm(complex a){
  return( sqrt( a.RE*a.RE+a.IM*a.IM) );}
/*---------------------------------------*/
void Pcomplex(complex a){
  cout<<a.RE<<"\t"<<a.IM<<"\n";
}
/*---------------------------------------*/
#endif /* !FILE_complex_SEEN */

