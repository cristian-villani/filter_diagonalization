#include"hfile.h"
#include"const.h"

void chkinput(int ndim,int L,double *sigma,double emin,double emax,
              double hmin,double hmax)
{
if(ndim<=0){ errore("matrix dimension should be positive integer"); }
if(L<=0){ errore("number of basis should be positive integer"); }
if(L>ndim){ errore("number of basis is greater than matrix dimension"); }
if(L>MAXL){ errore("sorry, the number of basis is too big"); }
if(emin>emax){ errore("emin > emax, wrong energy interval"); }
if(emin==emax){ errore("emin = emax, worng energy interval"); }
if(hmin>hmax){ errore("hmin > hmax, wrong spectral boundaries"); }
if(hmin==hmax){ errore("hmin = hmax, wrong spectral boundaries"); }
if(*sigma<=0){ *sigma=2.0*(emax-emin)/(L*(hmax-hmin)); }

printf("\n\n------------------------------------------------------\n");
printf(" Matrix dimension            :   %d\n",ndim);
printf(" Spectral boundaries         :   [%.3lf : %.3lf]\n",hmin,hmax);
printf(" Number of basis vectors     :   %d\n",L);
printf(" Width of the filter function:   %lf\n",*sigma);
printf(" Energy interval             :   [%.3lf : %.3lf]\n",emin,emax);
printf("------------------------------------------------------\n");
}

int errore(char *messaggio)
{
printf("\n\n######################################################\n");
printf("Stop. Abnormal termination of program.  Error:\n%s\n",messaggio);
printf("######################################################\n");
exit(1);
}
