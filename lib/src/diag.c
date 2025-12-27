#include<stdio.h>
#include<math.h>

#define SOGLIA 0.000000001

void diag(double *a,double *q,double *z,double *e,double *v,double hm,
          double dh,int N)
{
double cmm(double *,double *,double *,int,int,int,double);
double ortonormalizza(double *,double *,int,int);
void copia(double *,double *,int);
double qr,qr_prec,err,n2z,nqz;
int i,j,k,ncicli;

for(k=0;k<N;k++)
  {
  for(i=0;i<N;i++)
    z[i]=0;
  z[k]=1;
  ortonormalizza(z,v,k,N);
  copia(q,z,N);
/*for(i=0;i<N;i++)
    printf("%f\n",q[i]);*/

  ncicli=0;
  err=1;
  while(err>=SOGLIA)
    {
    ncicli++;
    cmm(a,q,z,N,N,1,1);
    nqz=cmm(q,z,&err,1,N,1,1);
    n2z=ortonormalizza(z,v,k,N);
    qr=nqz/n2z;
    err=fabs(qr_prec-qr);
    qr_prec=qr;
    copia(q,z,N);
    }
  e[k]=nqz;
  printf("%f\n",nqz*dh+hm);
  fflush(stdout);
  copia(&v[N*k],q,N);
  }
}

double ortonormalizza(double *z,double *v,int k,int N)
  {
  double cmm(double *,double *,double *,int,int,int,double);
  int i,m;
  double n2z,zv,fact,err;
  for(m=0;m<k;m++)
    {
    zv=cmm(z,&v[N*m],&v[N*k],1,N,1,1);
    for(i=0;i<N;i++) z[i]=z[i]-zv*v[N*m+i];
    }
  n2z=cmm(z,z,&err,1,N,1,1);
  fact=sqrt(1.0/n2z);
  for(i=0;i<N;i++)
    z[i]=z[i]*fact;
  return n2z;
  }
