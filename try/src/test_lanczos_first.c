#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fd_test_utils.h"
#include "fd_test_lanczos.h"

int main(void)
{
    /* -------------------------------
       Test parameters
       ------------------------------- */
    int ndim, nr8;
    int L = 6;
    int core = 0;
    int random = 1;
    int ioen = 0;
    int restart = 0;
    int L_default;

    char buf[128];

    double sigma, sigma_window, sigma_min, sigma_max, sigma_default;
    double emin = -5.5;
    double emax = 0.0;
    double hmin, hmax;
    double window;
    double prec1 = 1e-10;

    printf("Enter matrix dimension (ndim): ");
    if (scanf("%d%*c", &ndim) != 1 || ndim <= 0) {
       fprintf(stderr, "Invalid matrix size\n");
       exit(1);
    }
    nr8 = ndim;

    /* -------------------------------
       Generate test matrix
       ------------------------------- */
    generate_random_symmetric_matrix(ndim, 1.0);

    // int L_LANCZOS = 20;
    int L_LANCZOS = ndim; // getting the whole spectrum already
    double **V = malloc(L_LANCZOS * sizeof(double *));
    for(int i=0;i<L_LANCZOS;i++) V[i] = malloc(ndim * sizeof(double));

    double alpha[L_LANCZOS], beta[L_LANCZOS+1];
    lanczos(ndim, L_LANCZOS, V, alpha, beta);

    double eigvals[L_LANCZOS];
    diag_tridiag(L_LANCZOS, alpha, beta, eigvals);

    hmin = eigvals[0];                 // smallest Lanczos eigenvalue
    hmax = eigvals[L_LANCZOS-1];       // largest Lanczos eigenvalue
    double margin = 0.1 * (hmax - hmin);
    hmin -= margin;
    hmax += margin;

    printf("Lanczos approximate eigenvalues:\n");
    for(int i=0;i<L_LANCZOS;i++) printf("% .6f\n", eigvals[i]);

    printf("Running filter diagonalization unit test\n");
    printf("Matrix size: %d\n", ndim);


    // -----------------------------
    // Ask user for L
    // -----------------------------
    L_default = ndim / 10;
    if(L_default < 1) L_default = 1;
    printf("Enter subspace dimension L [default %d]: ", L_default);
    if(fgets(buf, sizeof(buf), stdin) && buf[0] != '\n'){
        if(sscanf(buf, "%d", &L) != 1 || L <= 0 || L > ndim){
            printf("Invalid L, using default %d\n", L_default);
            L = L_default;
        }
    } else {
        L = L_default;
    }

    // -----------------------------
    // Ask user for energy window
    // -----------------------------
    printf("Enter target energy window (emin emax) inside(default) %.6f %.6f: ",
           hmin, hmax);
    if(fgets(buf, sizeof(buf), stdin) && buf[0] != '\n'){
        if(sscanf(buf, "%lf %lf%*c", &emin, &emax) != 2 || emin >= emax){
            printf("Invalid energy window, using full spectral bounds\n");
            emin = hmin;
            emax = hmax;
        }
    } else {
        emin = hmin;
        emax = hmax;
    }

   // -----------------------------
    // Compute suggested sigma value
    // -----------------------------
    window = emax - emin;

    // Suggested defaults
    sigma_window   = 0.5 * window;            // half window width
    sigma_min      = 0.02 * window;           // 2% of window width
    sigma_max      = 0.3  * window;           // 30% of window width

    sigma = 0.3 * sigma_window;

    // Clamp sigma inside reasonable bounds
    if(sigma < sigma_min) sigma = sigma_min;
    if(sigma > sigma_max) sigma = sigma_max;

    // -----------------------------
    // Allow user to change sigma
    // -----------------------------
    sigma_default = sigma;
    printf("Enter filter width sigma [suggested %.6f]: ", sigma_default);
    if (scanf("%lf%*c", &sigma) != 1 || sigma <= 0.0) {
       printf("Using default sigma=%.6f\n", sigma_default);
       sigma = sigma_default;
    }

    printf("Using L=%d, sigma=%.6f, emin=%.6f, emax=%.6f\n",
           L, sigma, emin, emax);

    /* -------------------------------
       Call library
       ------------------------------- */
    finit_(&ndim, &L, &sigma, &emin, &emax, &hmin, &hmax, &core,
           &random, &ioen, &restart, &prec1, &nr8);

    // print_results(L, ndim, emin, emax, hmin, hmax, sigma);

    free(G_H);

    printf("\ntest finished with exit code:0\n");
    return 0;
}
