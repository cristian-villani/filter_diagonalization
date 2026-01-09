#include "fd_test_utils.h"
#include <math.h>
#include <fcntl.h>
#include <unistd.h>


/* ---------------- Global matrix ---------------- */
int G_ndim = 0;
double *G_H = NULL;

/* ---------------- Filter-diagonalization globals ---------------- */
int robin_L = 0;
int robin_ndim = 0;
double robin_F_min = 0.0;
double robin_F_max = 0.0;
double robin_spec_min = 0.0;
double robin_spec_max = 0.0;
double robin_prec = 1e-4;

/* ---------------- Random symmetric matrix ---------------- */
void generate_random_symmetric_matrix(int n, double scale)
{
    G_ndim = n;
    G_H = malloc(n * n * sizeof(double));
    if (!G_H) { fprintf(stderr, "Matrix allocation failed\n"); exit(1); }

    srand(12345);  // deterministic for tests

    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            double r = scale * (2.0 * rand() / RAND_MAX - 1.0);
            G_H[i*n + j] = r;
            G_H[j*n + i] = r;
        }
    }
    /* Normalization
    double max = 0.0;
    for(int i = 0; i < G_ndim*G_ndim; i++)
        if(fabs(G_H[i]) > max) max = fabs(G_H[i]);
    for(int i = 0; i < G_ndim*G_ndim; i++){
        G_H[i] /= max;
        G_H[i] = G_H[i]*5.0;
    }
     */
}

// Function to generate a matrix with high and low density spectrum
void create_high_low_density_spectrum(int ndim, double *eigenvalues) {
    // Define the ranges for the low-density and high-density regions
    int low_density_size = 10;
    // Number of eigenvalues in the low-density region (near 0)
    double low_density_range[] = {-5.0, -4.0, -3.0, -2.0, -1.0, 
                                   0.0, 1.0, 2.0, 3.0, 4.0}; // Example range
    
    // Fill the low-density region first (eigenvalues around 0)
    for (int i = 0; i < low_density_size; i++) {
        eigenvalues[i] = low_density_range[i];
    }

    // High-density negative region (between -70 and -50)
    //
    // Number of eigenvalues in the high-density negative region
    int high_density_negative_size = 15;
    double high_density_negative_min = -70.0;
    double high_density_negative_max = -50.0;
    
    for (int i = low_density_size; i < low_density_size + 
                                       high_density_negative_size; i++) {
        eigenvalues[i] = (rand() % (int)(high_density_negative_max
                                        - high_density_negative_min))
                                             + high_density_negative_min;
    }

    // High-density positive region (between 50 and 100)
    //
    // Number of eigenvalues in the high-density positive region
    int high_density_positive_size = 15;
    double high_density_positive_min = 50.0;
    double high_density_positive_max = 100.0;

    for (int i = low_density_size + high_density_negative_size;
         i < low_density_size + high_density_negative_size + 
         high_density_positive_size; i++) {
             eigenvalues[i] = (rand() % (int)(high_density_positive_max
                            - high_density_positive_min))
                            + high_density_positive_min;
    }

    // Remaining eigenvalues (random, sparse range for the extreme ends)
    double remaining_min = -200.0;
    double remaining_max = 200.0;

    for (int i = low_density_size + high_density_negative_size +
         high_density_positive_size; i < ndim; i++) {
               eigenvalues[i] = (rand() % (int)(remaining_max - 
                                                remaining_min))
                              + remaining_min;
    }
}

// Function to generate a random orthogonal matrix using QR decomposition
void generate_random_orthogonal_matrix(int n, double* Q) {
    // Create a random matrix (not orthogonal yet)
    double* A = (double*)malloc(n * n * sizeof(double));
    if (!A) {
        fprintf(stderr, "Matrix allocation failed\n");
        exit(1);
    }

    for (int i = 0; i < n * n; i++) {
        // Random values between -1 and 1
        A[i] = (2.0 * rand() / RAND_MAX - 1.0);
    }

    // QR Decomposition: Uses Gram-Schmidt to create an orthogonal matrix
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            Q[i * n + j] = A[i * n + j];  // Copy A into Q
        }
    }

    // Normalize each column to create an orthogonal matrix
    for (int j = 0; j < n; j++) {
        // Normalize column j
        double norm = 0.0;
        for (int i = 0; i < n; i++) {
            norm += Q[i * n + j] * Q[i * n + j];
        }
        norm = sqrt(norm);
        for (int i = 0; i < n; i++) {
            Q[i * n + j] /= norm;
        }

        // Orthogonalize against previous columns
        for (int k = 0; k < j; k++) {
            double dot_product = 0.0;
            for (int i = 0; i < n; i++) {
                dot_product += Q[i * n + k] * Q[i * n + j];
            }
            for (int i = 0; i < n; i++) {
                Q[i * n + j] -= dot_product * Q[i * n + k];
            }
        }
    }
    // Free the temporary matrix A
    free(A);
}

// Creates a diagonal matrix and transform it with random orthogonal matrix
void generate_diagonal_transformed_matrix(int ndim, double* eigenvalues) {
    // Allocate memory for the final matrix
    G_ndim = ndim;
    G_H = (double*)malloc(ndim * ndim * sizeof(double));
    if (!G_H) {
        fprintf(stderr, "Matrix allocation failed\n");
        exit(1);
    }

    // Create a diagonal matrix D with the given eigenvalues
    double* D = (double*)malloc(ndim * ndim * sizeof(double));
    if (!D) {
        fprintf(stderr, "Matrix allocation failed\n");
        exit(1);
    }
    for (int i = 0; i < ndim; i++) {
        for (int j = 0; j < ndim; j++) {
            if (i == j) {
                // Set eigenvalue on the diagonal
                D[i * ndim + j] = eigenvalues[i];
            } else {
                D[i * ndim + j] = 0.0;
            }
        }
    }

    // Create a random orthogonal matrix Q
    double* Q = (double*)malloc(ndim * ndim * sizeof(double));
    if (!Q) {
        fprintf(stderr, "Matrix allocation failed\n");
        exit(1);
    }
    generate_random_orthogonal_matrix(ndim, Q);

    // Compute A = Q * D * Q^T
    // (since Q is orthogonal, Q^T is the transpose of Q)

    for (int i = 0; i < ndim; i++) {
        for (int j = 0; j < ndim; j++) {
            G_H[i * ndim + j] = 0.0;
            for (int k = 0; k < ndim; k++) {
                G_H[i * ndim + j] += Q[i * ndim + k] * D[k * ndim + j];
            }
        }
    }
    // Free the temporary matrices D and Q
    free(D);
    free(Q);
}

// Function to print matrix (for debugging purposes)
void print_matrix(int ndim, double* matrix) {
    for (int i = 0; i < ndim; i++) {
        for (int j = 0; j < ndim; j++) {
            printf("%f ", matrix[i * ndim + j]);
        }
        printf("\n");
    }
}

// Matrix generation function
void generate_isolated_matrix(int n, double minval, double maxval) {
    // Set global dimension
    G_ndim = n;

    // Step 1: Allocate memory for the matrix
    G_H = (double*)malloc(n * n * sizeof(double));
    if (!G_H) {
        fprintf(stderr, "Matrix allocation failed\n");
        exit(1);
    }

    // Step 2: Initialize matrix with zeros
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            G_H[i * n + j] = 0.0;
        }
    }

    // Step 3: Define blocks for isolated eigenvalue regions
    int block_size_1 = 3;
    int block_size_2 = 3;
    int block_size_3 = 3;

    // Step 4: Generate first block with eigenvalues between [-2.0, -2.0001]
    double min_eigenvalue_1 = -2.0;
    double max_eigenvalue_1 = -2.0001;
    for (int i = 0; i < block_size_1; i++) {
        G_H[i * n + i] = min_eigenvalue_1 + (max_eigenvalue_1 
                                          - min_eigenvalue_1) 
                                          * rand() / (RAND_MAX + 1.0);
    }

    // Step 5: Generate second block with eigenvalues between [1.0, 1.0001]
    double min_eigenvalue_2 = 1.0;
    double max_eigenvalue_2 = 1.0001;
    for (int i = block_size_1; i < block_size_1 + block_size_2; i++) {
        G_H[i * n + i] = min_eigenvalue_2 + (max_eigenvalue_2 
                                          - min_eigenvalue_2)
                                          * rand() / (RAND_MAX + 1.0);
    }

    // Step 6: Generate third block with eigenvalues between [3.0, 3.0001]
    double min_eigenvalue_3 = 3.0;
    double max_eigenvalue_3 = 3.0001;
    for (int i = block_size_1 + block_size_2; i < block_size_1 + 
                 block_size_2 + block_size_3; i++) {
        G_H[i * n + i] = min_eigenvalue_3 + (max_eigenvalue_3
                       - min_eigenvalue_3) * rand() / (RAND_MAX + 1.0);
    }

    // Step 7: Add small off-diagonal random values for weaker coupling
    // Perturbation is now much smaller
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            double random_value = 0.00000001 * (2.0 * rand() / RAND_MAX - 1.0);
            // Extremely small random value for off-diagonal elements
            G_H[i * n + j] = random_value;
            G_H[j * n + i] = random_value;
        }
    }

    // Step 8: Normalize the matrix to ensure proper scaling
    double max_diag_value = 0.0;
    for (int i = 0; i < n; i++) {
        if (fabs(G_H[i * n + i]) > max_diag_value) {
            max_diag_value = fabs(G_H[i * n + i]);
        }
    }

    // Normalize matrix
    for (int i = 0; i < n * n; i++) {
        G_H[i] /= max_diag_value;
    }

    // Step 9: Print the matrix (just to check the structure)
    printf("Generated normalized matrix:\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%f ", G_H[i * n + j]);
        }
        printf("\n");
    }
}

void generate_block_matrix(int n, double minval, double maxval) {
    G_ndim = n;

    // Step 1: Allocate memory for the matrix
    G_H = (double*)malloc(n * n * sizeof(double));
    if (!G_H) {
        fprintf(stderr, "Matrix allocation failed\n");
        exit(1);
    }

    // Step 2: Initialize matrix with zeros
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            G_H[i * n + j] = 0.0;
        }
    }

    // Step 3: Define block sizes
    int block_size_1 = n / 3;
    int block_size_2 = n / 3;

    // Step 4: Generate the first block (eigenvalues between [-2.0, -1.5])
    double min_eigenvalue_1 = -2.0;
    double max_eigenvalue_1 = -1.5;
    for (int i = 0; i < block_size_1; i++) {
        G_H[i * n + i] = min_eigenvalue_1 + (max_eigenvalue_1
                                           - min_eigenvalue_1)
                                           * rand() / (RAND_MAX + 1.0);
    }

    // Step 5: Generate the second block (eigenvalues between [1.0, 2.0])
    double min_eigenvalue_2 = 1.0;
    double max_eigenvalue_2 = 2.0;
    for (int i = block_size_1; i < block_size_1 + block_size_2; i++) {
        G_H[i * n + i] = min_eigenvalue_2 + (max_eigenvalue_2
                                          - min_eigenvalue_2) 
                                          * rand() / (RAND_MAX + 1.0);
    }

    // Step 6: The third block can have values close to zero,
    // ensuring the eigenvalues are isolated
    for (int i = block_size_1 + block_size_2; i < n; i++) {
        G_H[i * n + i] = 0.01 * rand() / (RAND_MAX + 1.0);
        // Small diagonal entries
    }

    // Step 7: Add random off-diagonal entries to create weak coupling
    // Ensure the matrix is symmetric
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            double random_value = 0.05 * (2.0 * rand() / RAND_MAX - 1.0);
            // Small random value
            G_H[i * n + j] = random_value;
            G_H[j * n + i] = random_value;
        }
    }

    // Step 8: Normalize the matrix to ensure proper scaling
    double max_diag_value = 0.0;
    for (int i = 0; i < n; i++) {
        if (fabs(G_H[i * n + i]) > max_diag_value) {
            max_diag_value = fabs(G_H[i * n + i]);
        }
    }

    // Normalize matrix
    for (int i = 0; i < n * n; i++) {
        G_H[i] /= max_diag_value;
    }

    // Step 9: Print the matrix (just to check the structure)
    printf("Generated normalized matrix:\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%f ", G_H[i * n + j]);
        }
        printf("\n");
    }
}

void create_matrix_with_isolated_eigenvalues(int n,
                                             double min_val,
                                             double max_val) {
    G_ndim = n;
    G_H = malloc(n * n * sizeof(double));
    if (!G_H) {
        fprintf(stderr, "Matrix allocation failed\n");
        exit(1);
    }

    // Initialize the matrix to zero
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            G_H[i * n + j] = 0.0;
        }
    }

    srand(12345);
    // Seed for deterministic results (can be removed for randomness)

    // Block 1: High eigenvalue density (between -10 and -5, 30 eigenvalues)
    for (int i = 0; i < 30; i++) {
        G_H[i * n + i] = -10.0 + (rand() % 5);
        // Random spread in the range [-10, -5]
    }

    // Gap between -5 and 0.0 (no eigenvalues in this range)

    // Block 2: Two isolated eigenvalues (e.g., -1.0 and 0.0)
    G_H[30 * n + 30] = -1.0;
    G_H[31 * n + 31] = 0.0;

    // Gap between 0.0 and 5.0 (no eigenvalues here, another "nichts")

    // Block 3: High eigenvalue density (between 5 and 10, 40 eigenvalues)
    for (int i = 32; i < 72; i++) {
        G_H[i * n + i] = 5.0 + (rand() % 5);
        // Random spread in the range [5, 10]
    }
    // Gap between 10 and 20 (no eigenvalues in this range, another "nichts")

    // Block 4: Three isolated eigenvalues (e.g., 20.0, 25.0, and 30.0)
    G_H[72 * n + 72] = 20.0;
    G_H[73 * n + 73] = 25.0;
    G_H[74 * n + 74] = 30.0;

    // Add small random values in the upper triangle to make it symmetric
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            double r = 0.1 * (2.0 * rand() / RAND_MAX - 1.0);
            // Perturbation to induce interactions
            G_H[i * n + j] += r;
            G_H[j * n + i] += r; // Make matrix symmetric
        }
    }

    // Rescale the matrix to fit the desired spectral range
    double spectral_range = max_val - min_val;
    double max_element = -1e10;
    // We will find the maximum element in the matrix to normalize it
    for (int i = 0; i < n * n; i++) {
        if (fabs(G_H[i]) > max_element) {
            max_element = fabs(G_H[i]);
        }
    }

    // Normalize the matrix by its maximum value
    // this will ensure that eigenvalues fit within the desired range
    double scaling_factor = spectral_range / max_element;
    for (int i = 0; i < n * n; i++) {
        G_H[i] *= scaling_factor;
    }
}


/* ---------------- Matrix-vector multiplication ---------------- */
void accaqu_(double *x, double *y, double *work)
{
    (void)work;  // unused
    for (int i = 0; i < G_ndim; i++) {
        double sum = 0.0;
        for (int j = 0; j < G_ndim; j++)
            sum += G_H[i*G_ndim + j] * x[j];
        y[i] = sum;
    }
}

/* ---------------- Linear algebra helpers ---------------- */
double dot(double *x, double *y, int n)
{
    double s = 0.0;
    for(int i=0;i<n;i++) s += x[i]*y[i];
    return s;
}

double dps2(double *x, int n)
{
    double s = 0.0;
    for (int i = 0; i < n; i++) s += x[i]*x[i];
    return s;
}

void molt(double *x, int n, double a) { for (int i=0;i<n;i++) x[i]*=a; }
void somma(double *y, double *x, double a, double b, int n) {
   for(int i=0;i<n;i++) y[i]=a*y[i]+b*x[i]; }
void sottrai(double *x, double *y, int n) { for(int i=0;i<n;i++) x[i]-=y[i]; }
void copia(double *x, double *y, int n) { for(int i=0;i<n;i++) x[i]=y[i]; }
void drnd(double *x, int n) {
   for(int i=0;i<n;i++) x[i]=rand()/(double)RAND_MAX; }

/* ---------------- File I/O for FD ---------------- */
void read_fil_psi_transposed(double **PS, int L, int ndim)
{
    int file = open("fil_psi", O_RDONLY);
    if(file<0){ perror("open fil_psi"); exit(1); }

    double *buf = malloc(L*sizeof(double));
    for(int j=0;j<ndim;j++){
        if(read(file, buf, L*sizeof(double)) != L*sizeof(double)){
             perror("read fil_psi"); exit(1); }
        for(int i=0;i<L;i++) PS[i][j] = buf[i];
    }
    free(buf);
    close(file);
}

void read_blkfil_hpsimat(double **Hpsi, int L, int ndim)
{
    int file = open("blkfil_hpsimat", O_RDONLY);
    if(file<0){ perror("open blkfil_hpsimat"); exit(1); }

    double *buf = malloc(ndim*sizeof(double));
    for(int i=0;i<L;i++){
        if(read(file, buf, ndim*sizeof(double)) != ndim*sizeof(double)){
                 perror("read blkfil_hpsimat"); exit(1); }
        for(int j=0;j<ndim;j++) Hpsi[i][j] = buf[j];
    }
    free(buf);
    close(file);
}

/* ---------------- Residual checks & print ---------------- */
void check_and_print(double **PS, double **Hpsi, int L, int ndim)
{
    printf("\n=== PS vectors (first 5 elements) ===\n");
    for(int i=0;i<L;i++){
        printf("PS[%d]: ", i);
        for(int j=0;j<ndim && j<5;j++) printf("% .6f ", PS[i][j]);
        printf("\n");
    }

    printf("\n=== Norms of PS vectors ===\n");
    for(int i=0;i<L;i++){
        double s=0.0;
        for(int j=0;j<ndim;j++) s+=PS[i][j]*PS[i][j];
        printf("||PS[%d]|| = %.6e\n", i, sqrt(s));
    }

    printf("\n=== Residual check: ||Hpsi - H*PS|| ===\n");
    double *Hx = malloc(ndim*sizeof(double));
    for(int i=0;i<L;i++){
        accaqu_(PS[i], Hx, NULL);
        for(int j=0;j<ndim;j++) Hx[j]-=Hpsi[i][j];
        double res=0.0;
        for(int j=0;j<ndim;j++) res+=Hx[j]*Hx[j];
        printf("Residual[%d] = %.6e\n", i, sqrt(res));
    }
    free(Hx);
}

void compute_eigenvalues_unscaled(double **PS, double **Hpsi, int L, int ndim,
                         double emin, double emax)
{
    printf("\n=== Rayleigh quotients ===\n");
    for(int i=0;i<L;i++){
        double num=0.0, den=0.0;
        for(int j=0;j<ndim;j++){
            num += PS[i][j]*Hpsi[i][j];
            den += PS[i][j]*PS[i][j];
        }
        double lambda = num/den;
        printf("i=%d  λ=% .6f  %s\n", i, lambda,
              (lambda >= emin && lambda <= emax) ? "<-- in window" : "");
    }
}
// Compute physical energies, reject null vectors, and print
void compute_eigenvalues_rescaled(double **PS, double **Hpsi, int L, int ndim,
                                  double emin, double emax, double Hmin,
                                  double Hmax)
{
    printf("\n=== Rayleigh quotients (PHYSICAL energies) ===\n");

    for(int i = 0; i < L; i++) {
        // Compute vector norm to detect null vector
        double norm2 = 0.0;
        for(int j = 0; j < ndim; j++) norm2 += PS[i][j] * PS[i][j];

        if(norm2 < 1e-12) {
            printf("i=%2d  NULL vector, skipped\n", i);
            continue;  // skip null vector
        }

        // Compute Rayleigh quotient in FD space
        double num = 0.0;
        for(int j = 0; j < ndim; j++) num += PS[i][j] * Hpsi[i][j];
        double lambda_fd = num / norm2;

        // Rescale to physical energy
        double lambda_phys = rescale_energy_from_fd(lambda_fd, Hmin, Hmax);

        // Check if in target window
        const char *flag = (lambda_phys >= emin && lambda_phys <= emax) ?
                            "<-- in window" : "";
        printf("i=%2d  λ=%.8f  %s\n", i, lambda_phys, flag);
    }
}

void compute_eigenvalues(double **PS, double **Hpsi, int L, int ndim,
                         double emin, double emax,
                         double Hmin, double Hmax)
{
    double dh = 0.5 * (Hmax - Hmin);
    double hm = 0.5 * (Hmax + Hmin);

    printf("\n=== Rayleigh quotients (PHYSICAL energies) ===\n");

    for (int i = 0; i < L; i++) {
        double num = 0.0;
        double den = 0.0;

        for (int j = 0; j < ndim; j++) {
            num += PS[i][j] * Hpsi[i][j];  // <ψ|H|ψ> in FD-scaled units
            den += PS[i][j] * PS[i][j];    // ||ψ||^2
        }

        double lambda_scaled = num / den;  // FD-scaled Rayleigh quotient
        double lambda_phys   = lambda_scaled * dh + hm;
                                           // rescale to physical units

        printf("i=%2d  λ=%.8f  %s\n",
               i, lambda_phys,
               (lambda_phys >= emin && lambda_phys <= emax) ? "<-- in window" : "");
    }
}

void print_results(int L, int ndim, double emin, double emax,
                   double hmin, double hmax, double sigma) {
    int i;

    /* Allocate PS */
    double **PS = malloc(L * sizeof(double *));
    for (i = 0; i < L; i++)
        PS[i] = malloc(ndim * sizeof(double));

    /* Allocate auxiliary vectors */
    double *z0 = malloc(ndim * sizeof(double));
    double *z1 = malloc(ndim * sizeof(double));
    double *z2 = malloc(ndim * sizeof(double));

    /* Check allocation */
    if (!z0 || !z1 || !z2 || !PS) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }
    for (i = 0; i < L; i++) {
        if (!PS[i]) {
            fprintf(stderr, "Memory allocation failed for PS[%d]\n", i);
            exit(1);
        }
    }

    /* Allocate Hpsi */
    double **Hpsi = malloc(L * sizeof(double *));
    for (i = 0; i < L; i++)
        Hpsi[i] = malloc(ndim * sizeof(double));

    /* Read data */
    read_fil_psi_transposed(PS, L, ndim);
    read_blkfil_hpsimat(Hpsi, L, ndim);

    /* Check and print */
    check_and_print(PS, Hpsi, L, ndim);
    compute_eigenvalues(PS, Hpsi, L, ndim, emin, emax, hmin, hmax);

    test_energy_localization(PS, Hpsi, L, ndim, emin, emax, sigma, hmin, hmax);
    test_subspace_consistency(PS, L, ndim);
    test_parameter_summary(L, sigma, emin, emax, hmin, hmax);


    /* Free memory */
    for (i = 0; i < L; i++) {
        free(PS[i]);
        free(Hpsi[i]);
    }
    free(PS);
    free(Hpsi);
    free(z0);
    free(z1);
    free(z2);
}

void estimate_hmin_hmax(int ndim, double *hmin, double *hmax)
{
    double min = 1e300, max = -1e300;
    for(int i=0;i<ndim;i++){
        double row_sum = 0.0;
        for(int j=0;j<ndim;j++){
            if(j!=i) row_sum += fabs(G_H[i*ndim + j]);
        }
        double lower = G_H[i*ndim + i] - row_sum;
        double upper = G_H[i*ndim + i] + row_sum;
        if(lower < min) min = lower;
        if(upper > max) max = upper;
    }
    *hmin = min; *hmax = max;

    printf("Estimated hmin = %.6f, hmax = %.6f\n", *hmin, *hmax);
}

void test_energy_localization(double **PS, double **Hpsi,
                              int L, int ndim, double emin, double emax,
                              double sigma, double hmin, double hmax)
{
    double hm = 0.5*(hmin + hmax);
    double dh = 0.5*(hmax - hmin);

    printf("\n=== Test A: Energy localization (PHYSICAL) ===\n");
    printf("Target window: [%.6f, %.6f]\n", emin, emax);
    printf("Filter width σ = %.6f\n\n", sigma);

    for(int i=0;i<L;i++){
        /* <ψ|H|ψ> */
        double Ei_scaled = dot(PS[i], Hpsi[i], ndim);
        double Ei = dh*Ei_scaled + hm;

        /* Variance: <(H - <H>)^2> */
        double var = 0.0;
        for(int j=0;j<ndim;j++){
            double Hj_phys = dh*Hpsi[i][j] + hm*PS[i][j];
            double r = Hj_phys - Ei*PS[i][j];
            var += r*r;
        }

        double dE = sqrt(var);

        printf("ψ[%2d]: <H>=%.6f   ΔE=%.6f", i, Ei, dE);

        if(Ei < emin || Ei > emax)
            printf("   [outside window]");
        if(dE > 3.0*sigma)
            printf("   [broad]");

        printf("\n");
    }
}

void test_subspace_consistency(double **PS, int L, int ndim)
{
    printf("\n=== Test B: Subspace consistency ===\n");
    printf("Overlap matrix S_ij = <ψ_i | ψ_j>\n");

    for(int i=0;i<L;i++){
        for(int j=0;j<L;j++){
            double sij = dot(PS[i], PS[j], ndim);
            printf("%8.4f ", sij);
        }
        printf("\n");
    }

    printf("\nDiagonal should be ~1, off-diagonal <~0.2\n");
}

void test_parameter_summary(int L, double sigma,
                            double emin, double emax,
                            double hmin, double hmax)
{
    printf("\n=== Test C: Parameter summary ===\n");
    printf("Subspace dimension L   = %d\n", L);
    printf("Filter width σ         = %.6f\n", sigma);
    printf("Target window          = [%.6f, %.6f]\n", emin, emax);
    printf("Spectral bounds        = [%.6f, %.6f]\n", hmin, hmax);

    printf("\nNotes:\n");
    printf(" • Larger L improves resolution inside window\n");
    printf(" • Smaller σ sharpens filter but may reduce stability\n");
    printf(" • Window outside spectrum ⇒ near-zero vectors\n");
}

/* -----------------------------
   Scale physical window to [-1,1]
   Hmin, Hmax  : spectral bounds of H
   Emin, Emax  : target energy window
   sigma_phys  : desired filter width in physical units
   Outputs:
   *Esmin, *Esmax   : scaled target window
   *sigma_scaled    : scaled filter width
--------------------------------*/
void scale_window_for_fd(double Hmin, double Hmax,
                         double Emin, double Emax,
                         double sigma_phys,
                         double *Esmin, double *Esmax,
                         double *sigma_scaled)
{
    *Esmin = 2.0 * (Emin - Hmin) / (Hmax - Hmin) - 1.0;
    *Esmax = 2.0 * (Emax - Hmin) / (Hmax - Hmin) - 1.0;

    *sigma_scaled = 2.0 * sigma_phys / (Hmax - Hmin);
}

/* -----------------------------
   Map scaled FD energies back to physical energies
   E_scaled : energy in scaled units [-1,1]
   Hmin, Hmax : spectral bounds
--------------------------------*/
double rescale_energy_from_fd(double E_scaled, double Hmin, double Hmax)
{
    return 0.5 * (Hmax - Hmin) * E_scaled + 0.5 * (Hmax + Hmin);
}

/* ---------------------------
 * C wrapper for set_robin_input
 * --------------------------- */
void set_robin_input_(int *dim1, int *L1, double *F_min1, double *F_max1,
                      double *spec_min1, double *spec_max1, double *prec1)
{
    robin_ndim = *dim1;   /* ndim */
    robin_L    = *L1;     /* L */
    robin_F_min = *F_min1;
    robin_F_max = *F_max1;
    robin_spec_min = *spec_min1;
    robin_spec_max = *spec_max1;
    // robin_prec = 1e-4;

    /* The rest can be ignored for now */
    printf("set_robin_input: ndim=%d, L=%d\n", robin_ndim, robin_L);
}

double dot_product(double *v1, double *v2, int ndim) {
    double dot = 0.0;
    for (int i = 0; i < ndim; i++) {
        dot += v1[i] * v2[i];
    }
    return dot;
}

#define PACK(i,j) ((i)*((i)+1)/2 + (j))

void analyze_(double *psi, double *hpsi,
              double *Hmat_packed, double *Smat_packed,
              double *evalH,
              double *work_unused, double *vec_unused)
{
    int i, j, k, l;
    int inc = 1;

    int L    = robin_L;
    int ndim = robin_ndim;

    const double TOL = 1e-4;   // overlap cutoff

    /* --------------------------------------------------
     * 1. Build packed H and S
     * -------------------------------------------------- */
    for (i = 0; i < L; i++) {
        for (j = 0; j <= i; j++) {
            Smat_packed[PACK(i,j)] =
                ddot_(&ndim, psi + i*ndim, &inc,
                             psi + j*ndim, &inc);

            Hmat_packed[PACK(i,j)] =
                ddot_(&ndim, psi + i*ndim, &inc,
                             hpsi + j*ndim, &inc);
        }
    }

    /* --------------------------------------------------
     * 2. Unpack S to full matrix
     * -------------------------------------------------- */
    double *S = malloc(L*L*sizeof(double));
    for (i=0;i<L;i++)
        for (j=0;j<L;j++)
            S[i+j*L] = (i>=j) ? Smat_packed[PACK(i,j)]
                              : Smat_packed[PACK(j,i)];

    /* --------------------------------------------------
     * 3. Diagonalize overlap matrix
     * -------------------------------------------------- */
    double *evalS = malloc(L*sizeof(double));
    int lwork = 3*L;
    double *work = malloc(lwork*sizeof(double));
    int info;
    char jobz='V', uplo='U';

    dsyev_(&jobz,&uplo,&L,S,&L,evalS,work,&lwork,&info);
    if (info!=0) { printf("dsyev(S) failed\n"); exit(1); }

    /* --------------------------------------------------
     * 4. Determine effective rank
     * -------------------------------------------------- */
    double smax = evalS[L-1];
    int Lgood = 0;
    for (i=0;i<L;i++)
        if (evalS[i] > TOL*smax)
            Lgood++;

    printf("Effective subspace dimension: %d / %d\n", Lgood, L);
    if (Lgood==0) {
        printf("Sorry, no eigenvector found\n");
        printf("Try again with different sigma value or subspace dimension\n");
        exit(1); }

    /* --------------------------------------------------
     * 5. Build orthonormal basis Q = S^{-1/2} U
     * -------------------------------------------------- */
    double *Q = malloc(L * Lgood * sizeof(double));
    int col = 0;

    for (i=0;i<L;i++) {
        if (evalS[i] > TOL*smax) {
            double scale = 1.0 / sqrt(evalS[i]);
            for (j=0;j<L;j++)
                Q[j + col*L] = S[j + i*L] * scale;
            col++;
        }
    }

    /* --------------------------------------------------
     * 6. Unpack H to full matrix
     * -------------------------------------------------- */
    double *H = malloc(L*L*sizeof(double));
    for (i=0;i<L;i++)
        for (j=0;j<L;j++)
            H[i+j*L] = (i>=j) ? Hmat_packed[PACK(i,j)]
                              : Hmat_packed[PACK(j,i)];

    /* --------------------------------------------------
     * 7. Project H → reduced basis: Hr = Qᵀ H Q
     * -------------------------------------------------- */
    double *Hr = calloc(Lgood*Lgood,sizeof(double));

    for (i=0;i<Lgood;i++)
        for (j=0;j<Lgood;j++)
            for (k=0;k<L;k++)
                for (l=0;l<L;l++)
                    Hr[i+j*Lgood] +=
                        Q[k+i*L] * H[k+l*L] * Q[l+j*L];

    /* --------------------------------------------------
     * 8. Diagonalize reduced Hamiltonian
     * -------------------------------------------------- */
    dsyev_(&jobz,&uplo,&Lgood,Hr,&Lgood,evalH,work,&lwork,&info);
    if (info!=0) { printf("dsyev(Hr) failed\n"); exit(1); }


    /* --------------------------------------------------
     * 9. Compute Ritz residuals
     * -------------------------------------------------- */

     double *coeff = calloc(L, sizeof(double));
     double *phi   = malloc(ndim * sizeof(double));
     double *Hphi  = malloc(ndim * sizeof(double));

     printf("\n--- Ritz residuals ---\n");

     for (i = 0; i < Lgood; i++) {

         /* coeff = Q * Hr(:,i) */
         for (k = 0; k < L; k++) {
             coeff[k] = 0.0;
             for (j = 0; j < Lgood; j++)
                 coeff[k] += Q[k + j*L] * Hr[j + i*Lgood];
         }

         /* phi = sum_k coeff[k] * psi_k */
         for (j = 0; j < ndim; j++) {
             phi[j] = 0.0;
             for (k = 0; k < L; k++)
                 phi[j] += coeff[k] * psi[k*ndim + j];
         }

         /* Hphi = H * phi */
         accaqu_(phi, Hphi, NULL);

         /* residual = Hphi - evalH[i] * phi */
         for (j = 0; j < ndim; j++)
              Hphi[j] -= evalH[i] * phi[j];

         double res = sqrt(ddot_(&ndim, Hphi, &inc, Hphi, &inc));

         printf("E[%2d] = % .12f   residual = %.3e\n", i, evalH[i], res);
   }

   free(coeff);
   free(phi);
   free(Hphi);


    /* --------------------------------------------------
     * 10. Print final eigenvalues
     * -------------------------------------------------- */
    printf("\n--- Physical eigenvalues ---\n");
    for (i=0;i<Lgood;i++)
        printf("evalH[%d] = %.12f\n", i, evalH[i]);

    /* --------------------------------------------------
     * 11. Cleanup
     * -------------------------------------------------- */
    free(S);
    free(evalS);
    free(Q);
    free(H);
    free(Hr);
    free(work);
}

/* exit program with message(wrapper) */
void sbagliato(const char *msg) {
    fprintf(stderr, "%s\n", msg);
    exit(1);  // Exit the program with an error code
}
