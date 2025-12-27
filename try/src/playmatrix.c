#include <stdio.h>
#include <stdlib.h>

/* Forward declarations of Fortran routines */
void reai_(double *e, int *L, double *emin, double *emax, double *hm, double *dh);
void diag(double *matrix, double *eigenvalues, double *eigenvectors, int n);

/* Create a simple symmetric 5x5 matrix */
void create_test_matrix(double *matrix, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            double value = (i == j) ? (n-i) : 0.1*(i+j);
            matrix[i*n + j] = value;
            matrix[j*n + i] = value; // symmetric
        }
    }
}

/* Write a namelist input file for reai_ */
void write_namelist(double energy) {
    FILE *f = fopen("energies.nml", "w");
    if (!f) {
        perror("Cannot open energies.nml for writing");
        exit(1);
    }
    fprintf(f, "&energies\n");
    fprintf(f, "e = %f\n", energy);
    fprintf(f, "/\n");
    fclose(f);
}

int main() {
    int n = 5;
    double matrix[25];
    double eigenvalues[5];
    double eigenvectors[25];

    create_test_matrix(matrix, n);

    double energy;
    printf("Enter the target energy: ");
    scanf("%lf", &energy);

    /* Write the Fortran namelist file */
    write_namelist(energy);

    /* Define a small window around the target */
    double emin = energy - 0.5;
    double emax = energy + 0.5;
    double hm = (emax - emin) / 2.0;
    double dh = (emax + emin) / 2.0;

    double target_energies[1] = {0}; // will be read from the file by reai_
    int L = 1;

    /* Call reai_ (will read energies from energies.nml) */
    reai_(target_energies, &L, &emin, &emax, &hm, &dh);

    /* Compute eigenvalues/eigenvectors */
    diag(matrix, eigenvalues, eigenvectors, n);

    printf("\nEigenvalues found:\n");
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
