#include <math.h>
#include "PDE.h"
#include "grid.h"
#include "hydro.h"
using namespace std;
#define setDT 0.001
double get_dt(double dt_hydro);
/*-----------------------------*/
double GetDPsi(double Psi[], double DPsi[], double time, const int ldiag, double *pdt){
  double dt_hydro,dt;
  /* First step of evolution */
  for(int ijk=0;ijk<NN;ijk++){
    vec3 kk = get_kk(ijk);
    if (lhydro == 1){dOdt1(Psi, kk, ijk, ldiag);}
  }
  /* Inverse Fourier Transforms */
  if (lhydro == 1){
    dOdt_IFFT(Psi);}
  /* Second step of evolution */
  for(int ijk=0;ijk<NN;ijk++){
    vec3 kk = get_kk(ijk);
    if (lhydro == 1){
      double dt_hydro=dOdt2(Psi, DPsi, time, kk, ijk, ldiag);
    }
  }
  /* Decide on dt */
  dt=get_dt(dt_hydro);
  if (dt!=0.){
    *pdt=dt;}
  /* Direct Fourier transforms */
  if (lhydro == 1) {
    dOdt_FFT(DPsi);}
  /* Final Step */
  for(int ijk=0;ijk<NN;ijk++){
    vec3 kk = get_kk(ijk);
    if (lhydro == 1){
      dOdt_final(DPsi, kk, ijk);}
  }
}
  /* ------------------------------ */
double get_dt(double dt_hydro){
  return dt_hydro;
}
