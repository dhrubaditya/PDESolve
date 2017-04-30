#ifndef FILE_grid_SEEN
#define FILE_grid_SEEN
/*---------------------------------------*/
#include <iostream>
#include<math.h>
#include "3vec.h"
using namespace std;
/*-------------------------*/
int N1;
int N2;
int N3;
long int Nsqr=N1*N2; 
long int Ncube=N1*N2*N3;
void get_kk(vec3 kk[], const int ij);
#endif /* !FILE_grid_SEEN */
