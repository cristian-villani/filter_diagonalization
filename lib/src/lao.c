void copia(double *x,double *y,int num)
{
int i;
for(i=0;i<num;i++)
  *(x+i)=*(y+i);
}

void molt(double *a,int n,double d)
{
int i;

for(i=0;i<n;i++)
 {
 *a=*a*d;
 a++;
 }
}

double dps(double *a,double *b,int n)
{
int i;
double ps;

ps=0.0;
for(i=0;i<n;i++)
 {
 ps=ps+*a**b;
 a++;
 b++;
 }
return ps;
}

double dps2(double *a,int n)
{
int i;
double ps;

ps=0.0;
for(i=0;i<n;i++)
 {
 ps=ps+*a**a;
 a++;
 }
return ps;
}

double cmm(double *indx, double *indy, double *indz, int mx,
        int nx, int ny, double alfa)
{
double *ind, *ind1, *ind2, *indz0, var;
int i,j,k;

indz0=indz;
for(j=1;j<=ny;j++)
 {
 ind=indx;
 ind1=indy;
 ind2=indz;
 for(i=1;i<=mx;i++)
  {
  var=0;
  for(k=0;k<nx;k++)
   {
   var+=*(ind)*(*(ind1+k*ny));
   ind++;
   }
  *ind2=alfa*var;
  ind2+=ny;
  }
 indz++;
 indy++;
 }
return(*indz0);
}

void sottrai(double *v1,double *v2,int k)
{
int i;
for(i=0;i<k;i++)
 {
 *v1=*v1-*v2;
 v1++;
 v2++;
 }
}

void somma(double *v1,double *v2,double alfa,double beta,int k)
{
int i;
for(i=0;i<k;i++)
 {
 *v1=*v1*alfa+*v2*beta;
 v1++;
 v2++;
 }
}
