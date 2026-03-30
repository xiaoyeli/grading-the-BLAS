#include "init_matrices.hpp"
#include "utils.hpp"
#include <vector>
#include <cmath>

#ifdef CUDA
#include "cublas_dgemm_wrapper.hpp"
#endif

#ifndef DRYRUN
#include <H5Cpp.h>
#endif

#ifdef STOCHASTIC
#include <ctime>
#endif

#ifdef VERBOSE
#include <stdio.h>
#endif

#define HDF5_FILENAME "grading_dgemm.h5"

#ifndef DGEMM_IMPL
#define DGEMM_IMPL "unknown"
#endif

#ifndef INPUT_TYPE
#define INPUT_TYPE "none"
#endif

#ifndef MAX_DIM
#define MAX_DIM 1024
#endif

#ifndef MIN_DIM
#define MIN_DIM 32
#endif

using namespace std;

extern "C" {
    void dgemm_(const char* TRANSA, const char* TRANSB,
                const int* M, const int* N, const int* K,
                const double* ALPHA,
                const double* A, const int* LDA,
                const double* B, const int* LDB,
                const double* BETA,
                double* C, const int* LDC);

#ifdef EXTENDED_PRECISION_REFERENCE
    void ldgemm_(const char* TRANSA, const char* TRANSB,
                const int* M, const int* N, const int* K,
                const long double* ALPHA,
                const long double* A, const int* LDA,
                const long double* B, const int* LDB,
                const long double* BETA,
                long double* C, const int* LDC);
#endif
}

int main(){

#ifndef DRYRUN
    H5::H5File h5_file = init_hdf5(HDF5_FILENAME);
#endif

    char transpose = 'N';
    double alpha = 1.0;
    double beta = 0.0;

    int max_m, max_k, max_n, lda, ldb, ldc;
    max_m = max_k = max_n = lda = ldb = ldc = MAX_DIM;

    vector<double> A(max_m*max_k);
    vector<double> B(max_k*max_n);
    vector<double> C(max_m*max_n);
#ifdef EXTENDED_PRECISION_REFERENCE
    long double alpha_ext;
    long double beta_ext;
    vector<long double> A_ext(max_m*max_k);
    vector<long double> B_ext(max_k*max_n);
    vector<long double> C_ext(max_m*max_n);
#endif

    for ( int dim = MIN_DIM; dim <= MAX_DIM; dim = 2*dim ){
        
        int m, n, k;
        m = n = k = dim;

        // initialize A and B matrices
#ifdef STOCHASTIC
        int seed = time(nullptr);
#else
        int seed = 0;
#endif
        init_matrices(INPUT_TYPE, seed, m, n, k, &A[0], lda, &B[0], ldb);
#ifndef DRYRUN
        write_matrix_to_hdf5(h5_file, &A[0], lda, "A", INPUT_TYPE, SCALING_FACTOR, m, k, seed);
        write_matrix_to_hdf5(h5_file, &B[0], ldb, "B", INPUT_TYPE, SCALING_FACTOR, k, n, seed);
#endif
#ifdef VERBOSE
        printf("A matrix:\n");
        print_matrix_to_stdout_sci(m, k, &A[0], lda);
        printf("B matrix:\n");
        print_matrix_to_stdout_sci(k, n, &B[0], ldb);
#endif

#ifdef EXTENDED_PRECISION_REFERENCE 
        alpha_ext = alpha;
        beta_ext = beta;
        copy_matrix(m, k, &A_ext[0], lda, &A[0], lda);
        copy_matrix(m, k, &B_ext[0], ldb, &B[0], ldb);
        init_zero_matrix(m, n, &C_ext[0], ldc);
        ldgemm_(&transpose, &transpose, &m, &n, &k, &alpha_ext, &A_ext[0], &lda, &B_ext[0], &ldb, &beta_ext, &C_ext[0], &ldc);
#ifndef DRYRUN
        write_matrix_to_hdf5(h5_file, &C_ext[0], ldc, DGEMM_IMPL, INPUT_TYPE, SCALING_FACTOR, m, n, seed);
#endif
        abs_matrix(m, k, &A_ext[0], lda);
        abs_matrix(m, k, &B_ext[0], ldb);
        init_zero_matrix(m, n, &C_ext[0], ldc);
        ldgemm_(&transpose, &transpose, &m, &n, &k, &alpha_ext, &A_ext[0], &lda, &B_ext[0], &ldb, &beta_ext, &C_ext[0], &ldc);
#ifndef DRYRUN
        string new_name = "ABS_A_AND_B_";
        new_name = new_name + DGEMM_IMPL;
        write_matrix_to_hdf5(h5_file, &C_ext[0], ldc, new_name, INPUT_TYPE, SCALING_FACTOR, m, n, seed);
#endif
#else
        init_zero_matrix(m, n, &C[0], ldc);
#ifdef CUDA
        cublas_dgemm_wrapper(m, n, k, &A[0], lda, &B[0], ldb, &C[0], ldc);
#else
        dgemm_(&transpose, &transpose, &m, &n, &k, &alpha, &A[0], &lda, &B[0], &ldb, &beta, &C[0], &ldc);
#endif
#ifdef VERBOSE
        printf("C matrix:\n");
        print_matrix_to_stdout_sci(m, n, &C[0], ldc);
#endif
#ifndef DRYRUN
        write_matrix_to_hdf5(h5_file, &C[0], ldc, DGEMM_IMPL, INPUT_TYPE, SCALING_FACTOR, m, n, seed);
#endif
#endif
    }
#ifndef DRYRUN
    H5Fclose(h5_file.getId());
#endif
    return 0;
}