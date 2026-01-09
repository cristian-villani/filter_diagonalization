#include <stdio.h>
#include <stdlib.h>
#include "fd_test_utils.h"

int main(void)
{
    int ndim, L;
    int L_default;
    int core = 0, random = 1, ioen = 0, restart = 0;
    int nr8;
    double sigma;
    double sigma_default;
    double prec1 = 1e-10;
    double emin, emax;
    double hmin, hmax;


    printf("Enter matrix dimension (ndim): ");
    if (scanf("%d%*c", &ndim) != 1 || ndim <= 0) {
       fprintf(stderr, "Invalid matrix size\n");
       exit(1);
    }

    L_default = ndim/100;
    if(L_default == 0) L_default = ndim/10;
    if(L_default == 0) L_default = ndim;

    printf("Enter subspace dimension L [default %d]: ", L_default);
    if (scanf("%d%*c", &L) != 1 || L <= 0 || L > ndim) {
       printf("Using default L=%d\n", L_default);
       L = L_default;
    }

    /* -------------------------------
       Generate matrix and estimate spectrum
       ------------------------------- */
    generate_random_symmetric_matrix(ndim, 1.0);
    estimate_hmin_hmax(ndim, &hmin, &hmax);

    printf("Estimated spectral bounds: hmin=%.6f hmax=%.6f\n", hmin, hmax);
    printf("Enter target energy window (emin emax): ");
    if (scanf("%lf %lf%*c", &emin, &emax) != 2 || emin >= emax) {
        fprintf(stderr, "Invalid energy window\n");
        exit(1);
    }

    sigma_default = 0.08/(hmax - hmin);
    printf("Enter filter width sigma [default %.6f]: ", sigma_default);
    if (scanf("%lf%*c", &sigma) != 1 || sigma <= 0.0) {
       printf("Using default sigma=%.6f\n", sigma_default);
       sigma = sigma_default;
    }

    /* -------------------------------
       Call filter diagonalization
       ------------------------------- */
    finit_(&ndim, &L, &sigma, &emin, &emax, &hmin, &hmax,
           &core, &random, &ioen, &restart, &prec1, &nr8);

    print_results(L, ndim, emin, emax, hmin, hmax, sigma);
    free(G_H);

    printf("\ntest finished with exit code:0\n");
    return 0;
}

