#include<stdio.h>

void scrivi_rinfo(double *PS[],double *z0,double *z1,double *z2,int kstart,
                  int L,int ndim)
{
  int i,j,nbytes,file;

  unlink("fil_psi");
  file=creat("fil_psi",0644);
  if(file<=0){
    printf("Error in creating file fil_psi\n");
    exit(1);
  }

  printf("\nWriting psi vectors on fil_psi file...\n");
  fflush(stdout);
  nbytes=ndim*sizeof(double);

  if(write(file,z0,nbytes) != nbytes){
    printf("Error in writing vector z0 to file fil_psi\n");
    exit(1);
  }
  if(write(file,z1,nbytes) != nbytes){
    printf("Error in writing vector z1 to file fil_psi\n");
    exit(1);
  }
  if(write(file,z2,nbytes) != nbytes){
    printf("Error in writing vector z2 to file fil_psi\n");
    exit(1);
  }
  if(write(file,&kstart,sizeof(int))!= sizeof(int)){
    printf("Error in writing kstart to file fil_psi\n");
    exit(1);
  }
  for(i=0;i<L;i++){
    if(write(file,PS[i],nbytes) != nbytes){
      printf("Error in writing PS[%d] vectors to file fil_psi\n",i);
      exit(1);
    }
  }
  close(file);
}

void scrivit_psi(double *PS[],double *buf,int L,int ndim)
{
  int i,j,nbytes,file;

  unlink("fil_psi");
  file=creat("fil_psi",0644);
  if(file<=0){
    printf("Error in creating file fil_psi\n");
    exit(1);
  }

  printf("\nWriting psi transposed vectors on fil_psi file...\n");
  fflush(stdout);
  nbytes=L*sizeof(double);
  for(j=0;j<ndim;j++){
    for(i=0;i<L;i++){
      buf[i]=PS[i][j];
    }
    if(write(file,buf,nbytes) != nbytes){
      printf("Error in writing file fil_psi\n");
      exit(1);
    }
  }
  close(file);
}

void scrivi_hpsi(double *PS[],double *z2,double *v,int L,int ndim)
{
  void normalizza(double *,int);
  void accaqu_(double *,double *,double *);

  int i,j,nbytes,file;

  unlink("blkfil_hpsimat");
  file=creat("blkfil_hpsimat",0644);
  if(file<=0){
    printf("Error in creating file blkfil_hpsimat\n");
    exit(1);
  }

  printf("\nWriting Hpsi vectors on blkfil_hpsimat file...\n");
  fflush(stdout);
  nbytes=ndim*sizeof(double);
  for(i=0;i<L;i++){
    normalizza(PS[i],ndim);
    accaqu_(PS[i],z2,v);
    if(write(file,z2,nbytes) != nbytes){
      printf("Error in writing file blkfil_hpsimat\n");
      exit(1);
    }
  }
close(file);
}

void leggi_rinfo(double *PS[],double *z0,double *z1,double *z2,int *pkstart,
                 int L,int ndim)
{
  int i,nbytes,file;

  file=open("fil_psi",0);
  if(file<=0){
    printf("Error in opening file fil_psi\n");
    exit(1);
  }

  printf("\nReading psi vectors from fil_psi file...\n");
  fflush(stdout);
  nbytes=ndim*sizeof(double);
  if(read(file,z0,nbytes)!=nbytes){
    printf("Error in reading vector z0 from file fil_psi\n");
    exit(1);
  }
  if(read(file,z1,nbytes)!=nbytes){
    printf("Error in reading vector z1 from file fil_psi\n");
    exit(1);
  }
  if(read(file,z2,nbytes)!=nbytes){
    printf("Error in reading vector z2 from file fil_psi\n");
    exit(1);
  }
  if(read(file,pkstart,sizeof(int))!=sizeof(int)){
    printf("Error in reading kstart from file fil_psi\n");
    exit(1);
  }
  for(i=0;i<L;i++){
    if(read(file,PS[i],nbytes)!=nbytes){
      printf("Error in reading PS[%d] vector from file fil_psi\n",i);
      exit(1);
    }
  }
  close(file);
}

void leggit_psi(double *psi,int L,int ndim)
{
  int i,nbytes,file;
  double *aux;

  file=open("fil_psi",0);
  if(file<=0){
    printf("Error in opening file fil_psi\n");
    exit(1);
  }

  printf("\nReading psi transposed vectors from fil_psi file...\n");
  fflush(stdout);
  nbytes=L*sizeof(double);
  aux=psi;
  for(i=0;i<ndim;i++){
    if(read(file,aux,nbytes)!=nbytes){
      printf("Error in reading from file fil_psi\n");
      exit(1);
    }
    aux=aux+L;
  }
  close(file);
}
