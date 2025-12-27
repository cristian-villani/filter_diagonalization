#include <stdio.h>

/* Stub the missing error function */
void sbagliato() {
    fprintf(stderr, "sbagliato called!\n");
}

/* Fortran routines from libfilter.a */
void reai_(double *energy, int *L, double *emin, double *emax, double *hm, double *dh);

int main() {
    int L = 5;                       // Number of energies
    double energy[100];

    // Initialize energies (dummy values)
    for(int i = 0; i < L; i++) {
        energy[i] = 2.0 + i;         // e.g., 2.0, 3.0, 4.0...
    }

    double emin = 1.0, emax = 10.0;
    double hm = (emax - emin)/2.0;
    double dh = (emax + emin)/2.0;

    // Call Fortran subroutine
    reai_(energy, &L, &emin, &emax, &hm, &dh);

    printf("\nRescaled energies in C:\n");
    for(int i = 0; i < L; i++)
        printf("Energy %d: %f\n", i+1, energy[i]);

    return 0;
}
