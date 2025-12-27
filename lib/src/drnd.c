#include<stdio.h>
#include<stdlib.h>
#include<time.h>

void drnd(double *p,int dim)
{
time_t g;
int i;

g=time(NULL);
srand(g);

  for(i=0;i<dim;i++){
    *p=(double)rand()/RAND_MAX;
    p++;
  }
}
