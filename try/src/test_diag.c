#include <stdio.h>
#include <stdlib.h>

/* Stub missing library functions if necessary */
void sbagliato() {}
void accaqu_(double *out,double *vec,int *size) {}

/* Fortran routines */
void reai_(double *energy,int *L,double *emin,double *emax,double *hm,double *dh);
void diag(double *H,double *eval,double *vec,double *buf,double **PS,
          double sigma,double hm,double dh,int Ldiag);

int main() {
    int L = 5;
    double energy[100];
    double H[25];       // 5x5 symmetric matrix
    double eval[5];     // eigenvalues
    double vec[25];     // eigenvectors (column-major)
    double buf[25];     // buffer
    double *PS[5];      // placeholder pointers for PS
    double emin = 1.0, emax = 5.0;
    double hm = (emax - emin)/2.0;
    double dh = (emax + emin)/2.0;
    double sigma = 0.0;

    // Initialize a small symmetric 5x5 matrix (for demo)
    for(int i=0;i<L;i++) {
        for(int j=0;j<L;j++) {
            if(i==j) H[i*L+j] = 2.0 + i;  // diagonal
            else H[i*L+j] = 0.5;          // off-diagonal
        }
    }

    // Initialize PS pointers
    for(int i=0;i<L;i++)
        PS[i] = vec + i*L;

    // Initialize energy guesses (any values)
    for(int i=0;i<L;i++)
        energy[i] = 2.0 + 0.1*i;

    // Call reai_ to process energies from input
    reai_(energy, &L, &emin, &emax, &hm, &dh);

    // Call diag to compute eigenvalues/eigenvectors
    diag(H, eval, vec, buf, PS, sigma, hm, dh, L);

    // Print results
    printf("\nEigenvalues computed by diag:\n");
    for(int i=0;i<L;i++)
        printf("eval[%d] = %f\n", i, eval[i]);

    printf("\nEigenvectors (column-wise):\n");
    for(int i=0;i<L;i++) {
        for(int j=0;j<L;j++)
            printf("%f ", vec[i*L+j]);
        printf("\n");
    }

    return 0;
}
