#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Stub missing functions from library */
void sbagliato() {}
void accaqu_(double *out,double *vec,int *size) {}

/* Fortran routines */
void reai_(double *energy,int *L,double *emin,double *emax,double *hm,double *dh);
void diag(double *H,double *eval,double *vec,double *buf,double **PS,
          double sigma,double hm,double dh,int Ldiag);

int main() {
    int L = 2; // Number of energies / size of matrix
    double energy[100];
    double H[4]; // 2x2 symmetric matrix
    double eval[2]; // eigenvalues
    double vec[4];  // eigenvectors (2x2)
    double buf[4];  // buffer
    double *PS[2];  // placeholder pointers for PS
    double emin = 1.0, emax = 10.0;
    double hm = (emax - emin)/2.0;
    double dh = (emax + emin)/2.0;
    double sigma = 0.0;

    // Initialize symmetric matrix H
    H[0] = 2.0; H[1] = 1.0;
    H[2] = 1.0; H[3] = 3.0;

    // Initialize PS pointers (just pointing to vectors)
    PS[0] = vec; 
    PS[1] = vec + 2;

    // Initialize energies (can be guesses)
    energy[0] = 2.5; 
    energy[1] = 2.5;

    // Call reai_ to process energies
    reai_(energy, &L, &emin, &emax, &hm, &dh);

    // Call diag: compute eigenvalues/vectors
    diag(H, eval, vec, buf, PS, sigma, hm, dh, L);

    printf("\nEigenvalues computed by diag:\n");
    for(int i=0;i<L;i++)
        printf("eval[%d] = %f\n", i, eval[i]);

    printf("\nEigenvectors (column-wise):\n");
    for(int i=0;i<L;i++) {
        for(int j=0;j<L;j++)
            printf("%f ", vec[i*L + j]);
        printf("\n");
    }

    return 0;
}
