#include"hfile.h"

double cf0aip(double ai,int J,double sigma)
{
int j;
double f0,pi;

pi=acos(-1.0);
f0=0.0;
for(j=1;j<=J;j++)
 f0=f0+exp(-(pow((cos(pi*(j-0.5)/J)-ai)/sigma,2.0)));
f0=f0/J;

return f0;
}

double cfkaip(int k,double ai,int J,double sigma)
{
int j;
double fk,pi;

pi=acos(-1.0);
fk=0;
for(j=1;j<=J;j++)
 fk=fk+exp(-(pow((cos(pi*(j-0.5)/J)-ai)/sigma,2.0)))*cos(k*pi*(j-0.5)/J);
fk=fk*2/J;

return fk;
}
