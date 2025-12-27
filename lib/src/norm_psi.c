void norm_psi(double *psi,int L,int dim)
{
  void normalizza(double *,int);

  double *aux;
  int i;

  aux=psi;
  for(i=0;i<L;i++){
    normalizza(aux,dim);
    aux=aux+dim;
  }
}

#include<math.h>
#include"lao.h"

void normalizza(double *v,int dim)
{
  double norma;

  norma=dps2(v,dim);
  if(norma>0.0){
    molt(v,dim,1.0/sqrt(norma));
  }
  else{
    printf("\nWarning: null vector found!\n");
  }
}
