#include "hydro.h"
using namespace std;
const complex II=complex(0.,1.);
static fftw_plan pohat2oox,puhat2uux;
static fftw_plan pohat2ooy,puhat2uuy;
static fftw_plan pohat2ooz,puhat2uuz;
static fftw_plan pou2ouxh,pou2ouyh,pou2ouzh;
/*-----------------------------*/
/* The flow is integrated in vorticity-streamfunction formalism. 
    in Fourier space:
   \ddt{\omega(k,t)} = -i k \times FT(\omega(r,t)\times u(r,t) + f(r,t))
   where FT denotes Fourier Transform and f(r,t) is the external force.*/ 
/*-----------------------------*/
void register_hydro_var(int *nvar,  int *psi_index, int *dpsi_index){
  /* storage for omega in Fourier space */
  iooh=*psi_index;
  ioxh=*psi_index; ioyh=ioxh+Ncube; iozh=ioyh+Ncube;
  *psi_index = iozh+Ncube;
  /* storage for domega is Fourier space */
  idodth=*dpsi_index;
  idodtxh=idodth; idodtyh=idodtxh+Ncube; idodtzh=idodtyh+Ncube;
  *dpsi_index=idodth+Ncube;
  /* omega is real space is stored in the same location as dodt */
  ioo=idodth;
  iox=idodtxh;
  ioy=idodtyh;
  ioz=idodtzh;
  *nvar=*nvar+3; // Three components of omega.
}
/*-------------------------------------*/
void register_hydro_aux(int *naux, int *psi_index, int *dpsi_index){
  /* First auxiliary variable is uu */
  iuuh= *psi_index;
  iuxh=iuuh; iuyh=iux+Ncube; iuzh=iuy+Ncube;
  /* move the point for the next variable */
  *psi_index=iuzh+Ncube;
  /* uu and uuh are stored in the same place */
  iux=iuxh; iuy=iuyh; iuz=iuzh;
  *naux=*naux+3;//Three components of velocity.
}
/*-----------------------------*/
void allocate_hydro_pencils(){
  cvec3 *uuh = (cvec3 *)malloc((1+N3/2)*sizeof(cvec3));
  cvec3 *ooh = (cvec3 *)malloc((1+N3/2)*sizeof(cvec3));
  cvec3 *dodth = (cvec3 *)malloc((1+N3/2)*sizeof(cvec3));
  vec3 *kk = (vec3 *)malloc((1+N3/2)*sizeof(vec3));
  vec3 *uu = (vec3 *)malloc(N3*sizeof(vec3));
  vec3 *oo = (vec3 *)malloc(N3*sizeof(vec3));
  vec3 *dodt = (vec3 *)malloc(N3*sizeof(vec3));
}
/*-----------------------------*/
void initialize_FFT_hydro(double Psi[], double DPsi[]){
  /* plans for ooh-> oo */
  pohat2oox = plan_c2r_FFT3d_out_of_place( Psi, iox, DPsi, ioxh, N1, N2, N3);
  pohat2ooy = plan_c2r_FFT3d_out_of_place(Psi, ioy, DPsi, ioyh, N1, N2, N3);
  pohat2ooz = plan_c2r_FFT3d_out_of_place(Psi, ioz, DPsi, iozh, N1, N2, N3);
  /*plans for uuh-> uu */
  puhat2uux = plan_c2r_FFT3d( Psi, iuxh, iux, N1, N2, N3);
  puhat2uuy = plan_c2r_FFT3d( Psi, iuyh, iuy, N1, N2, N3);
  puhat2uuz = plan_c2r_FFT3d( Psi, iuzh, iuz, N1, N2, N3);
  /* plans for oXu -> oXuh */
  pou2ouxh = plan_r2c_FFT3d( DPsi, idodtx, idodtxh, N1, N2, N3);
  pou2ouyh = plan_r2c_FFT3d( DPsi, idodty, idodtyh, N1, N2, N3);
  pou2ouzh = plan_r2c_FFT3d( DPsi, idodtz, idodtzh, N1, N2, N3);
}
/*-----------------------------*/
void dOdt1(double Psi[], const int ij, const int ldiag){
  double ksqr;
  cvec3 Lambda;
  Psi2cvec(ooh, Psi, ij, iooh);
  for (int i3=0;i3<N3;i3++){
    ksqr=sqnorm(kk[i3]);
    Lambda=ooh[i3]/ksqr;
    uuh[i3]=cross(kk[i3],Lambda)*II;
    /* Calculate k space diagnostics here */
    if (ldiag==1){
    }
  }
  cvec2Psi(Psi, ij, iuuh,uuh);
}
/*-----------------------------*/
void dOdt_IFFT(){
  /* Inverse FFT omegah -> omega */
  IFFT3(pohat2oox); IFFT3(pohat2ooy); IFFT3(pohat2ooz);
  /* Inverse FFT uuh -> uu */
  IFFT3(puhat2uux);  IFFT3(puhat2uuy);  IFFT3(puhat2uuz);
}
/*-----------------------------*/
double dOdt2(double Psi[], double DPsi[], double time, const int ij, const int ldiag){
  double dt_hydro=0.;
  Psi2vec(oo, DPsi, ij, ioo); //Careful, this needs DPsi.
  Psi2vec(uu, Psi, ij, iuu);
  /* Calculate real space diagnostics here */
  for (int i3=0;i3<N3;i3++){
    if (ldiag==1){
      /* ------also set dt_hydro -----*/
      dt_hydro=0.001;
    }else{
      dt_hydro=0.;
    }
    dodt[i3]=cross(uu[i3],oo[i3]);
  }
  vec2Psi(DPsi, ij, idodt, dodt);
  return dt_hydro;
}
/*-------------------------------------*/
void dOdt_FFT(){
  FFT3(pou2ouxh); FFT3(pou2ouyh); FFT3(pou2ouzh);
}
/* ------------------------------------*/
void dOdt_final(double DPsi[], const int ij){
  Psi2cvec(dodth, DPsi, ij, idodth);
  for(int i3=0;i3<N3;i3++){
  dodth[i3]=cross(kk[i3],dodth[i3])*II;
  }
  cvec2Psi(DPsi, ij, idodth, dodth);
}
