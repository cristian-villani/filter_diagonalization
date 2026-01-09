#ifndef FD_TEST_UTILS_H
#define FD_TEST_UTILS_H

/* ---------------------------
 * Global variables set by set_robin_input
 * --------------------------- */
extern int robin_L;
extern int robin_ndim;

#include <stdio.h>
#include <stdlib.h>

void finit_( int *, int *, double *, double *, double *, double *, double *,
    int *, int *, int *, int *, double *, int *);
/* ---------------- Matrix generation ---------------- */
void generate_random_symmetric_matrix(int, double);
void generate_block_matrix(int, double, double);
void generate_isolated_matrix(int, double, double);
void generate_diagonal_transformed_matrix(int, double *);
void generate_random_orthogonal_matrix(int, double *);
void create_matrix_with_isolated_eigenvalues(int, double, double);
void create_high_low_density_spectrum(int, double *);
void estimate_hmin_hmax(int, double *, double *);

/* ---------------- Matrix-vector multiplication ---------------- */
void accaqu_(double *, double *, double *);

/* ---------------- Linear algebra helpers ---------------- */
double dot(double *, double *, int);
double dps2(double *, int);
void molt(double *, int, double);
void somma(double *, double *, double, double, int);
void sottrai(double *, double *, int);
void copia(double *, double *, int);
void drnd(double *, int);
double dot_product(double *, double *, int);
/* BLAS */
extern double ddot_(int*, double*, int*, double*, int*);
/* LAPACK (Fortran) routine used by analyze_ */
extern void dspgv_(int *, char *, char *, int *, double *, double *,
                   double *, double *, int *, double *, int *);
extern void dsygv_(int *, char *, char *, int *, double *, int *,
                   double *, int *, double *, double *, int *, int *);
extern void dsyev_(char *, char *, int *, double *, int *,
                   double *, double *, int *, int *);
extern void dspev_(char*, char*, int*, double*, double*, double*, int*,
                   double*, int*);

/* FD test utilities (C implementations of missing Fortran) */
extern void set_robin_input_(int *, int *, double *, double *, double *,
                             double *, double *);
void analyze_(double *, double *, double *, double *, double *, double *,
              double *);
void filter_nearly_dependent_vectors(double **, int *, int, double);


/* ---------------- File I/O for FD ---------------- */
void read_fil_psi_transposed(double **, int, int);
void read_blkfil_hpsimat(double **, int, int);

/* ---------------- Residual checks & print ---------------- */
void check_and_print(double **, double **, int, int);
void compute_eigenvalues(double **, double **, int, int, double, double,
                         double, double);
void compute_eigenvalues_rescaled(double **, double **, int, int, double, double,
                                  double, double);
void compute_eigenvalues_unscaled(double **, double **, int, int, double, double);
void print_results(int, int, double, double, double, double, double);
void test_energy_localization(double **, double **, int, int, double,
                              double, double, double, double);
void test_subspace_consistency(double **, int, int);
void test_parameter_summary(int, double, double, double, double, double);

/* ---------------- Scaling functions ---------------- */
void scale_window_for_fd(double, double, double, double,
                         double, double *, double *, double *);
double rescale_energy_from_fd(double, double, double);

/* ---------------- Global matrix storage ---------------- */
extern int G_ndim;
extern double *G_H;

#endif
