#ifndef FILE_hydro_SEEN
#define FILE_hydro_SEEN
/*---------------------------------------*/
#include <iostream>
#include<math.h>
#include "cvec3.h"
using namespace std;
double nu=1.;
/*-------------------------*/
void dOdt1(double Psi[], const vec3 kk, const int ijk, const int ldiag);
double dOdt2(double Psi[], double DPsi[], double time, const vec3 kk, const int ijk, const int ldiag);
void dOdt_final(double DPsi[], const vec3 kk, const int ijk);
void dOdt_IFFT(double Psi[]);
void dOdt_FFT(double DPsi[]);
/*-------------------------*/
#endif /* !FILE_hydro_SEEN */
