#include <unistd.h>

void vstopchk(int *k)
{
  if(unlink("fil_stop")==0)
    {
      printf("\n\nProgram halt requested.\n\n");
      exit(0);
    }
}
