#include <math.h>
#include"PDE.h"
using namespace std;
#define EXPRK2 1
#define RK2      2
void ETD2RK(double time, unsigned long int NN, double Psi[], double DPsi[], double DPsiOld[], int Nvar);
void RungeKutta2(double time, unsigned long int NN, double Psi[], double DPsi[], double DPsiOld[], int Nvar);
void evolve(int scheme, double time, unsigned long int NN, double Psi[], double DPsi[], double DPsiOld[], int Nvar);
/* ------------------------ */
void evolve(int scheme, double time, unsigned long int NN, double Psi[], double DPsi[], double DPsiOld[], int Nvar){
  switch(scheme){
  case EXPRK2:
    ETD2RK(time, NN, Psi, DPsi, DPsiOld, Nvar);
    break;
  case RK2:
    RungeKutta2(time, NN, Psi, DPsi, DPsiOld, Nvar);
    break;
  }
}
/* ------------------------ */
void ETD2RK(double time, unsigned long int NN, double Psi[], double DPsi[], double DPsiOld[], int Nvar){
  /* The equation we solve is:
     Psi_t = Gamma*Psi + DPsi(\Psi,t). */
  double dt;
  GetDPsi(Psi, DPsi, time, 1, &dt);
  DPsiOld=DPsi;
  for(int ivar=0;ivar<Nvar;ivar++){
    for(int ijk=0+NN*ivar;ijk<NN*(ivar+1);ijk++){
      double Gamma= GetGamma(ijk,ivar);
      double eGdt = exp(Gamma*dt);
      Psi[ijk] = Psi[ijk]*eGdt + DPsi[ijk]*(eGdt-1)/Gamma;
    }
  }
  GetDPsi(Psi, DPsi, time+dt, 0, NULL);
  for(int ivar=0;ivar<Nvar;ivar++){
    for(int ijk=0+NN*ivar;ijk<NN*(ivar+1);ijk++){
      double Gamma= GetGamma(ijk,ivar);
      double eGdt = exp(Gamma*dt);
      Psi[ijk] = Psi[ijk]+(DPsi[ijk]-DPsiOld[ijk])*(eGdt-1-Gamma*dt)/(dt*Gamma*Gamma);
    }
  }
}
/* ------------------------ */
void RungeKutta2(double time, unsigned long int NN, double Psi[], double DPsi[], double DPsiOld[], int Nvar){
  /* The equation we solve is:
     Psi_t = Gamma*Psi + DPsi(\Psi,t). 
     This uses standard RungeKutta 2nd order not 
     exponential algorithm.  */
  double dt;
  GetDPsi(Psi, DPsi, time, 1, &dt);
  DPsiOld=DPsi;
  for(int ivar=0;ivar<Nvar;ivar++){
    for(int ijk=0+NN*ivar;ijk<NN*(ivar+1);ijk++){
      double Gamma= GetGamma(ijk,ivar);
      double eGdt = exp(Gamma*dt);
      Psi[ijk] = Psi[ijk]*eGdt + DPsi[ijk]*(eGdt-1)/Gamma;
    }
  }
  GetDPsi(Psi, DPsi, time+dt, 0, NULL);
  for(int ivar=0;ivar<Nvar;ivar++){
    for(int ijk=0+NN*ivar;ijk<NN*(ivar+1);ijk++){
      double Gamma= GetGamma(ijk,ivar);
      double eGdt = exp(Gamma*dt);
      Psi[ijk] = Psi[ijk]+(DPsi[ijk]-DPsiOld[ijk])*(eGdt-1-Gamma*dt)/(dt*Gamma*Gamma);
    }
  }
}
/* ------------------------ */
/* void rnkt2(double yy[], double tt,double deltat, int lindex){
  double temp[pdim],k1[pdim];
  eval_rhs(k1,tt,yy,lindex);
  for(int idim=0;idim<pdim;idim++){
    temp[idim]=yy[lindex+idim]+k1[idim]*deltat/2.;
  }
  eval_rhs(k1,tt+(deltat/2.),temp,0);
  for(int idim=0;idim<pdim;idim++){
    yy[lindex+idim]=yy[lindex+idim]+deltat*k1[idim];
  }
  }*/
/********************************
__device__ void rnkt4(double yy[], double tt,double deltat, int lindex){
  double  temp[pdim],k1[pdim],k2[pdim],k3[pdim],k4[pdim];
  eval_rhs(k1,tt,yy,lindex);
  for(int idim=0;idim<pdim;idim++){
    temp[idim]=yy[lindex+idim]+k1[idim]*deltat/2.;
  }
  eval_rhs(k2,tt+(deltat/2.),temp,0);
  for(int idim=0;idim<pdim;idim++){
    temp[idim]=yy[lindex+idim]+k2[idim]*deltat/2.;
  }
  eval_rhs(k3,tt+(deltat/2.),temp,0);
  for(int idim=0;idim<pdim;idim++){
    temp[idim]=yy[lindex+idim]+k3[idim]*deltat;
  }
  eval_rhs(k4,tt+deltat,temp,0);
  for(int idim=0;idim<pdim;idim++){
    yy[lindex+idim]=yy[lindex+idim]+deltat*(  (k1[idim]/6.) + (k2[idim]/3.) + (k3[idim]/3.) + (k4[idim]/6.) );
  }
}*/
/*********************************/
