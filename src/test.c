#include <iostream>
#define N 16
#include<fftw3.h>
#include"FFT.h"
int main(){
  double* psi = (double *)fftw_malloc(sizeof(double)*(N+2)); 
  for(int i=0;i<N;i++){
    psi[i]=double(i);
  }
  //psi[0]=1.;
  cout<<"Input\n";
  for(int i=0;i<N;i++){
    cout<<psi[i]<<"\n";
  }
  fftw_plan p1d_r2c = plan_FFT1d(psi, 0,  0, N);
  cout <<"Output\n";
  FFT1d(p1d_r2c);
  for(int i=0;i<N+2;i++){
    cout<<psi[i]<<"\n";
  }
  cout<<"Now Back\n";
  fftw_plan p1d_c2r = plan_IFFT1d(psi, 0,  0, N);
  FFT1d(p1d_c2r);
  for(int i=0;i<N;i++){
    psi[i]=psi[i]/double(N);
    cout<<psi[i]<<"\n";
    }
}
