#include <math.h>
#include "PDE.h"
#include "grid.h"
#include "hydro.h"
using namespace std;
#define setDT 0.001
double get_dt(double dt_hydro);
void register_var();
void register_aux();
/*-----------------------------*/
void registrate(){
  lhydro=IsHydro();
  /* -- Register all the variables but not the auxiliaries */
  register_var();
  /* Now all the dynamical variables are allocated. */
  /* do a consistency check here. */
  register_aux();
}
/*-----------------------------*/
void register_var(){
  if (lhydro==1){
    register_hydro_var(&nvar, &psi_index, &dpsi_index);
  }
}
/*-----------------------------*/
void register_aux(){
  if (lhydro==1){
    register_hydro_aux(&nvar, &psi_index, &dpsi_index);
  }
}
/*-----------------------------*/
void allocate_pencils(){
  if (lhydro==1){allocate_hydro_pencils();}
}
/*-----------------------------*/
void Psi2vec(vec3 qq[], double Psi[], const int ij, const int iqq){
  for(int i3=0;i3<N3;i3++){
    double q1 = Psi[iqq+ij+i3];
    double q2 = Psi[iqq+ij+Ncube+i3];
    double q3 = Psi[iqq+ij+2*Ncube+i3];
    qq[i3]=vec3(q1,q2,q3);
  }
}
/*-----------------------------*/
void vec2Psi(double Psi[], const int ij, const int iqq, vec3 qq[]){
  for(int i3=0;i3<N3;i3++){
    Psi[iqq+ij+i3]=qq[i3].x;
    Psi[iqq+ij+Ncube+i3]=qq[i3].y;
    Psi[iqq+ij+2*Ncube+i3]=qq[i3].z;
  }
}
/*-----------------------------*/
void Psi2cvec(cvec3 qqh[], double Psi[], const int ij, const int iqqh){
  for(int i3=0;i3<(1+N3/2);i3++){
    int ieven=2*i3;
    int iodd=2*i3+1;
    complex q1=complex(Psi[iqqh+ij+ieven],Psi[iqqh+ij+iodd]);
    complex q2=complex(Psi[iqqh+ij+Ncube+ieven],Psi[iqqh+ij+Ncube+iodd]);
    complex q3=complex(Psi[iqqh+ij+2*Ncube+ieven],Psi[iqqh+ij+2*Ncube+iodd]);
    qqh[i3]=cvec3(q1,q2,q3);
  }
}
/*-----------------------------*/
void cvec2Psi(double Psi[], const int ij, const int iqqh, cvec3 qqh[]){
  for(int i3=0;i3<(1+N3/2);i3++){
    int ieven=2*i3;
    int iodd=2*i3+1;
    Psi[iqqh+ij+ieven]=qqh[i3].x.RE;
    Psi[iqqh+ij+iodd]=qqh[i3].x.IM;
    Psi[iqqh+ij+Ncube+ieven]=qqh[i3].y.RE;
    Psi[iqqh+ij+Ncube+iodd]=qqh[i3].y.IM;
    Psi[iqqh+ij+2*Ncube+ieven]=qqh[i3].z.RE;
    Psi[iqqh+ij+2*Ncube+iodd]=qqh[i3].z.IM;
  }
}
/*-----------------------------*/
double GetDPsi(double Psi[], double DPsi[], double time, const int ldiag, double *pdt){
  double dt_hydro,dt;
  /* First step of evolution */
  for(int isqr=0;isqr<Nsqr;isqr++){
    int ij=isqr*(N3+2);
    get_kk(kk,ij);
    if (lhydro == 1){dOdt1(Psi, ij, ldiag);}
  }
  /* Inverse Fourier Transforms */
  if (lhydro == 1){dOdt_IFFT(Psi);}
  /* Second step of evolution; we are in real space now */
  for(int isqr=0;isqr<Nsqr;isqr++){
    int ij=isqr*(N3+2);
    if (lhydro == 1){
      double dt_hydro=dOdt2(Psi, DPsi, time, ij, ldiag);
      dt_hydro=max(dt_hydro,setDT);
    }
  }
  /* Decide on dt */
  *pdt=dt_hydro;
  /* Direct Fourier transforms */
  if (lhydro == 1) {
    dOdt_FFT(DPsi);}
  /* Final Step; we are back in fourier space */
  for(int isqr=0;isqr<Nsqr;isqr++){
    int ij=isqr*(N3+2);
    get_kk(kk,ij);
    if (lhydro == 1){ dOdt_final(DPsi, ij);}
  }
}
  /* ------------------------------ */

