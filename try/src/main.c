#include<stdio.h>

void main()
{
  void reai_(double *,int *,double *,double *,double *,double *);

  double e[100],emin,emax,hm,dh;
  int L;

  L=10;
  emin=20.396;
  emax=20.543;
  hm=(emax-emin)/2.0;
  dh=(emax+emin)/2.0;

  reai_(e,&L,&emin,&emax,&hm,&dh);
}
