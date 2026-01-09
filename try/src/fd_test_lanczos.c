#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fd_test_utils.h"

/* -----------------------------
   Simple symmetric Lanczos
   ----------------------------- */
void lanczos(int ndim, int m, double **V, double *alpha, double *beta)
{
    double *v0 = calloc(ndim, sizeof(double));
    double *v1 = calloc(ndim, sizeof(double));
    double *w  = calloc(ndim, sizeof(double));

    // Step 0: random normalized vector
    drnd(v1, ndim);
    double norm = sqrt(dps2(v1, ndim));
    for(int i=0;i<ndim;i++) v1[i] /= norm;

    beta[0] = 0.0;

    for(int k=0;k<m;k++) {
        // w = H * v1 - beta[k] * v0
        accaqu_(v1, w, NULL);
        for(int i=0;i<ndim;i++) w[i] -= beta[k]*v0[i];

        alpha[k] = 0.0;
        for(int i=0;i<ndim;i++) alpha[k] += v1[i]*w[i];

        for(int i=0;i<ndim;i++) w[i] -= alpha[k]*v1[i];

        beta[k+1] = sqrt(dps2(w, ndim));
        if(beta[k+1] < 1e-12) break;  // converged

        // store Lanczos vector
        for(int i=0;i<ndim;i++) {
            V[k][i] = v1[i];
            v0[i] = v1[i];
            v1[i] = w[i]/beta[k+1];
        }
    }

    free(v0); free(v1); free(w);
}

/* -----------------------------
   Print tridiagonal eigenvalues
   ----------------------------- */
void diag_tridiag(int m, double *alpha, double *beta, double *eigvals)
{
    // Use simple symmetric QR iteration (or call LAPACK if available)
    // For simplicity, here we use GSL or LAPACK in a real code
    // This is just a placeholder:
    for(int i=0;i<m;i++) eigvals[i] = alpha[i]; // rough approximation
    // -----------------------------
    // Sort eigenvalues ascending
    // -----------------------------
    for(int i=0;i<m-1;i++){
        for(int j=i+1;j<m;j++){
            if(eigvals[i] > eigvals[j]){
                double tmp = eigvals[i];
                eigvals[i] = eigvals[j];
                eigvals[j] = tmp;
            }
        }
    }
}

