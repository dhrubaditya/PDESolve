#ifndef FILE_PDE_SEEN
#define FILE_PDE_SEEN
/*---------------------------------------*/
#include <iostream>
#include<math.h>
#include "cvec3.h"
using namespace std;
/*-------------------------*/
int NN;
int nvar, naux, psi_index, dpsi_index;
int lhydro;
int ioo=0,iuu=0,idodt=0;
int iux=0,iuy=0,iuz=0;
int iuxh=0,iuyh=0,iuzh=0;
int iox=0,ioy=0,ioz=0;
int ioxh=0,ioyh=0,iozh=0;
int iooh=0,iuuh=0,idodth=0;
int idodtx=0,idodty=0,idodtz=0;
int idodtxh=0,idodtyh=0,idodtzh=0;
int lmagnetic;
int iaa=0;
int lpscalar;
int ipscalar=0;
int lparticles;
int lIBM;
double GetDPsi(double Psi[], double DPsi[], double time, const int ldiag, double *pdt);
void Psi2cvec(cvec3 qqh[], double Psi[], const int ij, const int iqqh);
void cvec2Psi(double Psi[], const int ij, const int iqqh, cvec3 qqh[]);
void Psi2vec(vec3 qq[], double Psi[], const int ij, const int iqq);
void vec2Psi(double Psi[], const int ijk, const int iqq,vec3 qq[]);
double GetGamma(const int ij,const int ivar);
#endif /* !FILE_PDE_SEEN */
