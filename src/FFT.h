#ifndef FILE_FFT_SEEN
#define FILE_FFT_SEEN
/*---------------------------------------*/
#include <iostream>
#include<math.h>
#include<fftw3.h>
using namespace std;
/*-------------------------*/
void  IFFT3(double Psi[], const int iqqh, const int iqq);
void  FFT3(double Psi[], const int iqq, const int iqqh);
fftw_plan plan_FFT1d(double Psi[], const int iqq, const int iqqh, const int N1);
fftw_plan plan_IFFT1d(double Psi[], const int iqqh, const int iqq, const int N1);
void FFT1d(fftw_plan p1d_r2c);
#endif /* !FILE_FFT_SEEN */
