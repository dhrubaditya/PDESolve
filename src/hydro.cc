#include <math.h>
#include "PDE.h"
#include "hydro.h"
#include "FFT.h"
using namespace std;
const complex II=complex(0.,1.);
/*-----------------------------*/
/* The flow is integrated in vorticity-streamfunction formalism. 
    in Fourier space:
   \ddt{\omega(k,t)} = -i k \times FT(\omega(r,t)\times u(r,t) + f(r,t))
  where FT denotes Fourier Transform and f(r,t) is the external force. 
/*-----------------------------*/
void dOdt1(double Psi[], const vec3 kk, const int ijk, const int ldiag){
  cvec3 oohat=Psi2cvec(Psi, ijk, iooh);
  double ksqr=sqnorm(kk);
  cvec3 Lambda=oohat/ksqr;
  cvec3 uuhat=cross(kk,Lambda)*II;
  cvec2Psi(uuhat, Psi, ijk, iuuh);
  /* Calculate k space diagnostics here */
  if (ldiag==1){
  }
}
/*-----------------------------*/
void dOdt_IFFT(double Psi[]){
  IFFT3(Psi,iooh,ioo);
  IFFT3(Psi,iuuh,iuu);
}
/*-----------------------------*/
double dOdt2(double Psi[], double DPsi[], double time, const vec3 kk, const int ijk, const int ldiag){
  double dt_hydro=0.;
  vec3 oo=Psi2vec(Psi, ijk, ioo);
  vec3 uu=Psi2vec(Psi, ijk, iuu);
  /* Calculate real space diagnostics here */
  if (ldiag==1){
    /* ------also set dt_hydro -----*/
    dt_hydro=0.001;
  }else{
    dt_hydro=0.;
  }
  vec3 dodt=cross(uu,oo);
  vec2Psi(dodt, DPsi, ijk, idodt);
  return dt_hydro;
}
/*-------------------------------------*/
void dOdt_FFT(double DPsi[]){
  FFT3(DPsi,idodt,idodth);
}
/* ------------------------------------*/
void dOdt_final(double DPsi[], const vec3 kk, const int ijk){
  cvec3 dodth=Psi2cvec(DPsi, ijk, idodth);
  dodth=cross(kk,dodth)*II;
  cvec2Psi(dodth, DPsi, ijk, idodth);
}
