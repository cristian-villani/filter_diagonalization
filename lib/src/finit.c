#include"hfile.h"
#include"const.h"

finit_(int *pndim,int *pL,double *psigma,double *pemin,double *pemax,
      double *phmin,double *phmax,int *core,int *random,int *ioen,
      int *restart,double *pprec1,int *pnr8)
{
void ceai(double *,int,double,double,double,double);
void reai_(double *,int *,double *,double *,double *,double *);
double cf0aip(double,int,double);
double cfkaip(int,double,int,double);
void diag(double *,double *,double *,double *,double *,double,double,int);
void chkyn(int *,int,int *);
void chkfstop(int *);
void drnd(double *,int);
void chkinput(int,int,double *,double,double,double,double);
void accaqu_(double *,double *,double *);
void scrivi_hpsi(double **,double *,double *,int,int);
void scrivi_rinfo(double **,double *,double *,double *,int,int,int);
void scrivit_psi(double **,double *,int,int);
void leggi_rinfo(double **,double *,double *,double *,int *,int,int);
void leggit_psi(double *,int,int);
int errore(char *);

double *Smat,*v,*Hmat,*evalH,*buf,*psi;
int i,j,k,n,L,M,J,K,ndim,nel,mem,nr8,packdim,kstart;
int yn[MAXL];
double sigma,emin,emax,hmin,hmax,hm,dh,hsd,usd,prec1;
double fk,qq,q1,sq;
double *e,*p,*z0,*z1,*z2;
double *PS[MAXL];
int file;

pre();

ndim=*pndim;
L=*pL;
K=100000;
sigma=*psigma;
emin=*pemin;
emax=*pemax;
hmin=*phmin;
hmax=*phmax;
nr8=*pnr8;
prec1=*pprec1;

printf("\n\n%s%d%s\n","Memory requested for the matrix construction: ",
       nr8*sizeof(double)," bytes");

v=(double *)malloc(nr8*sizeof(double));
if(v<=(double *)0) { errore("Memory allocation");}

chkinput(ndim,L,&sigma,emin,emax,hmin,hmax);

J=K+1;
/* controlli vari sull'input */
/* calcolo della memoria     */
mem=((L+3)*ndim+3*L*L)*sizeof(double);

nel=ndim*ndim;
printf("\n%s%d%s\n\n","the filter diagonalizer will need at most ",
       mem," bytes");
fflush(stdout);

hm=(hmin+hmax)/2;
dh=(hmax-hmin)/2;

e=(double *)malloc(L*sizeof(double));
if(e<=(double *)0) { errore("Memory allocation");}

if(*ioen!=1){
  ceai(e,L,emin,emax,hm,dh);
}
else{
  reai_(e,&L,&emin,&emax,&hm,&dh); /* reading e vector from input */
}

printf("\n\nletto input\n\n");

z0=(double *)malloc(ndim*sizeof(double));
if(z0<=(double *)0) { errore("Memory allocation");}
z1=(double *)malloc(ndim*sizeof(double));
if(z1<=(double *)0) { errore("Memory allocation");}
z2=(double *)malloc(ndim*sizeof(double));
if(z2<=(double *)0) { errore("Memory allocation");}

for(i=0;i<L;i++)
 {
 PS[i]=(double *)malloc(ndim*sizeof(double));
 if(PS[i]<=(double *)0) { errore("Memory allocation");}
 }

usd=1.0/dh;
hsd=hm/dh;

for(i=0;i<L;i++)
 yn[i]=1;

if(*restart!=1){
  if(*random==0){
  /* carthesian vector usage */
    for(n=0;n<ndim;n++)
      z0[n]=0.0;
    z0[*core]=1.0;
  }
  else{
  /* random vector usage */
    drnd(z0,ndim);
    qq=1.0/sqrt(dps2(z0,ndim));
    molt(z0,ndim,qq);
  }
 
  for(i=0;i<L;i++){
    fk=cf0aip(e[i],J,sigma);
    for(n=0;n<ndim;n++)
      PS[i][n]=fk*z0[n];
  }

  accaqu_(z0,z1,v);
  molt(z1,ndim,usd);
  somma(z1,z0,1.0,-hsd,ndim);
 
  for(i=0;i<L;i++){
    fk=cfkaip(1,e[i],J,sigma);
    for(n=0;n<ndim;n++)
      PS[i][n]=PS[i][n]+fk*z1[n];
  }
  kstart=2;
  printf("Starting Chebishev expansion...\n\n");
  printf("Iteration n. 1 completed\n");
  fflush(stdout);
}
else{
  leggi_rinfo(PS,z0,z1,z2,&kstart,L,ndim);
  printf("Restarting Chebishev expansion from iteration %d\n\n",kstart);
  fflush(stdout);
}

unlink("blkfil_vecs");   /* deletes the old  vectors file  to avoid
                            final re-reading of the vectors in case 
                            of no convergence */
for(k=kstart;k<=K;k++)
 {
 printf("Iteration n. %d completed\n",k);
 fflush(stdout);

 accaqu_(z1,z2,v);
 molt(z2,ndim,usd);
 somma(z2,z1,1.0,-hsd,ndim);

 molt(z2,ndim,2.0);
 sottrai(z2,z0,ndim);
 for(i=0;i<L;i++)
  {
  fk=cfkaip(k,e[i],J,sigma);
  if((fk<SOGLIA2)&&(fk>-SOGLIA2)) yn[i]=0;
  for(n=0;n<ndim;n++)
   PS[i][n]=PS[i][n]+fk*z2[n];
  }
 copia(z0,z1,ndim);
 copia(z1,z2,ndim);
 chkyn(yn,L,&K);
 chkfstop(&K);
 }

if(K==J-1){
  printf("\nSorry, no convergence after %d iterations\n",K);
  fflush(stdout);
  free(e);
  free(z0);
  free(z1);
  for(i=0;i<L;i++){
    free(PS[i]);
  }
}
else if(K==-1){
  scrivi_rinfo(PS,z0,z1,z2,k,L,ndim);
  free(e);
  free(z0);
  free(z1);
  for(i=0;i<L;i++){
    free(PS[i]);
  }
}
else{
  buf=(double*)malloc(L*sizeof(double));
  if(buf<=(double *)0){ sbagliato("Error allocating memory (buf)");}

/*========================================================
 * Calculating Hpsi vectors and writing them to file
 *=======================================================*/
  scrivi_hpsi(PS,z2,v,L,ndim);

/*========================================================
 * writing psi vectors to file and then reading them
 * transposed for further fortran usage
 *=======================================================*/
  scrivit_psi(PS,buf,L,ndim);
  for(i=0;i<L;i++){
    free(PS[i]);
  }

  psi=(double *)malloc(L*ndim*sizeof(double));
  if(psi<=(double *)0){ sbagliato("Error allocating memory (psi)");}

  leggit_psi(psi,L,ndim);
/*========================================================*/

  packdim=L*(L+1)/2;

  free(e);
  free(z0);
  free(z1);

  Hmat=(double*)malloc(packdim*sizeof(double));
  if(Hmat<=(double *)0){ sbagliato("Error allocating memory (Hmat)");}
  Smat=(double*)malloc(packdim*sizeof(double));
  if(Smat<=(double *)0){ sbagliato("Error allocating memory (Smat)");}
  evalH=(double*)malloc(L*sizeof(double));
  if(evalH<=(double *)0){ sbagliato("Error allocating memory (evalH)");}

  printf("\nFilter diagonalization completed.");
  printf("\n----------------------------");
  printf("\nCalling analyze...\n\n");
  fflush(stdout);

  set_robin_input_(&ndim,&L,&emin,&emax,&hmin,&hmax,&prec1);
  analyze_(psi,z2,Hmat,Smat,evalH,buf,v);

  free(psi);
  free(Hmat);
  free(Smat);
  free(evalH);
  free(buf);
}

free(v);
free(z2);
}

void ceai(double e[],int L,double emin,double emax,double hm,double dh)
{
int i;

for(i=1;i<=L;i++)
 e[i-1]=((emin+(i+0.5)*(emax-emin)/(L+1))-hm)/dh;
}

void chkyn(int yn[],int n,int *N)
 {
 int i;

 for(i=0;i<n;i++)
  if(yn[i]==1) return;
 *N=0;
 }

void chkfstop(int *N)
{
  if(unlink("fil_stop")==0){
    printf("\n\nProgram halt requested.\n\n");
    fflush(stdout);
    *N=-1;
  }
}
