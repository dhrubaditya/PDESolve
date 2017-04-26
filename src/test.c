#include <iostream>
#define N 16
void increment(double *px){
  *px=*px+1.;
}
int main(){
  double y=0.;
  for(int i=0;i<N;i++){
    increment(&y);
    printf("%f\n",y);
  }
}
