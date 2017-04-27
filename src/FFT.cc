#include <math.h>
#include "FFT.h"
using namespace std;
/*-----------------------------*/
fftw_plan plan_r2c_FFT3d(double Psi[], const int iqq, const int iqqh, const int N1, const int N2, const int N3){
  double* in= (double *)&Psi[iqq];
  fftw_complex* out=(fftw_complex*)&Psi[iqqh];
  fftw_plan plan_r2c_3d= fftw_plan_dft_r2c_3d(N1, N2, N3, in, out, FFTW_ESTIMATE);
  return plan_r2c_3d;
}
/*-------------------------*/
fftw_plan plan_c2r_FFT3d(double Psi[], const int iqqh, const int iqq, const int N1, const int N2, const int N3){
  fftw_complex* in= (fftw_complex *)&Psi[iqqh];
  double* out=(double *)&Psi[iqq];
  fftw_plan plan_c2r_3d= fftw_plan_dft_c2r_3d(N1, N2, N3, in, out, FFTW_ESTIMATE);
  return plan_c2r_3d;
}
/*----------------------*/
void  IFFT3(double Psi[], const int iqqh, const int iqq){


}
void  FFT3(double Psi[], const int iqq, const int iqqh){
}
fftw_plan plan_FFT1d(double Psi[], const int iqq, const int iqqh, const int N1){
  double* in = (double*)&Psi[iqq];
  fftw_complex* out = (fftw_complex*)&Psi[iqqh];
  fftw_plan p1d_r2c=fftw_plan_dft_r2c_1d(N1, in, out, FFTW_ESTIMATE);
  return p1d_r2c;
}
/*--------------------------------*/
fftw_plan plan_IFFT1d(double Psi[], const int iqqh, const int iqq, const int N1){
  fftw_complex* in = (fftw_complex*)&Psi[iqqh];
  double* out = (double*)&Psi[iqqh];
  fftw_plan p1d_c2r=fftw_plan_dft_c2r_1d(N1, in, out, FFTW_ESTIMATE);
  return p1d_c2r;
}
/*--------------------------------*/
void FFT1d(fftw_plan plan1d){
  fftw_execute(plan1d);
}
/*--------------------------------*/
