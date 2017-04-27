#ifndef FILE_PDE_SEEN
#define FILE_PDE_SEEN
/*---------------------------------------*/
#include <iostream>
#include<math.h>
#include "cvec3.h"
using namespace std;
/*-------------------------*/
int NN;
int lhydro;
int ioo=0,iuu=0,idodt=0;
int iuux=0,iuuy=0,iuuz=0;
int ioox=0,iooy=0,iooz=0;
int iooh=0,iuuh=0,idodth=0;
int lmagnetic;
int iaa=0;
int lpscalar;
int ipscalar=0;
int lparticles;
int lIBM;
double GetDPsi(double Psi[], double DPsi[], double time, const int ldiag, double *pdt);
cvec3 Psi2cvec(double Psi[], const int ijk, const int iqq);
void cvec2Psi(cvec3 qq, double Psi[], const int ijk, const int iqq);
vec3 Psi2vec(double Psi[], const int ijk, const int iqq);
void vec2Psi(vec3 qq, double Psi[], const int ijk, const int iqq);
double GetGamma(const int ijk,const int ivar);
#endif /* !FILE_PDE_SEEN */
