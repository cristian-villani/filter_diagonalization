#include<stdio.h>
#include<unistd.h>

void ntpsi(double *psi,double *buf,int L,int dim)
{
  void norm_psi(double *,int,int);

  double *aux;
  int i;
  FILE *tmp;

  norm_psi(psi,L,dim);

  tmp=fopen("fil_aux","rw");

  for(aux=psi;aux<psi+dim;aux++){
    for(i=0;i<L;i++){
      buf[i]=*(aux+i*dim);
    }
    cwrite8_(tmp,L,buf);
  }

  aux=psi;
  for(i=0;i<dim;i++){
    cread8_(tmp,L,aux);
    aux=aux+L;
  }
  fclose(tmp);
}
