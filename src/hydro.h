#ifndef FILE_hydro_SEEN
#define FILE_hydro_SEEN
/*---------------------------------------*/
#include <iostream>
#include<math.h>
#include "cvec3.h"
#include "FFT.h"
#include "PDE.h"
#include "grid.h"
#include "pencils.h"
using namespace std;
double nu=1.;
/*-------------------------*/
int IsHydro(){return 1;}
void register_hydro_var(int *nvar, int *psi_index, int *dpsi_index);
void register_hydro_aux(int *naux, int *psi_index, int *dpsi_index);
void initialize_FFT_hydro(double Psi[], double DPsi[]);
void allocate_hydro_pencils();
void dOdt1(double Psi[], const int ijk, const int ldiag);
double dOdt2(double Psi[], double DPsi[], double time, const vec3 kk, const int ijk, const int ldiag);
void dOdt_final(double DPsi[], const int ijk);
void dOdt_IFFT(double Psi[]);
void dOdt_FFT(double DPsi[]);
/*-------------------------*/
#endif /* !FILE_hydro_SEEN */
