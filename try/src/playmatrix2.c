#include <stdio.h>
#include <stdlib.h>

/* Forward declarations of Fortran routines */
void reai_(double *e, int *L, double *emin, double *emax, double *hm, double *dh);
void diag(double *matrix, double *eigenvalues, double *eigenvectors, int n);

/* Create a symmetric 5x5 test matrix */
void create_test_matrix(double *matrix, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            double value = (i == j) ? (n - i) : 0.1*(i+j);
            matrix[i*n + j] = value;
            matrix[j*n + i] = value; // symmetric
        }
    }
}

int main() {
    int n = 5;
    double matrix[25];
    double eigenvalues[5];
    double eigenvectors[25];

    /* Create symmetric matrix */
    create_test_matrix(matrix, n);

    /* Energy parameters */
    int L = 5;
    double target_energies[5] = {0}; // will be read from stdin by reai_
    double emin = 2.0;
    double emax = 4.0;
    double hm = (emax - emin)/2.0;
    double dh = (emax + emin)/2.0;

    printf("Calling reai_...\n");

    /* Call reai_, read energies from stdin (energies.nml) */
    reai_(target_energies, &L, &emin, &emax, &hm, &dh);

    printf("\nRescaled energies:\n");
    for (int i = 0; i < L; i++) {
        printf("energy[%d] = %f\n", i, target_energies[i]);
    }

    /* Call diag to compute eigenvalues/eigenvectors */
    diag(matrix, eigenvalues, eigenvectors, n);

    printf("\nEigenvalues:\n");
    for (int i = 0; i < n; i++)
        printf("%d: %f\n", i, eigenvalues[i]);

    printf("\nEigenvectors (column-wise):\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            printf("%f ", eigenvectors[i*n + j]);
        printf("\n");
    }

    return 0;
}
