#include "init_matrices.hpp"
#include "utils.hpp"
#include <vector>
#include <random>

#ifdef CUDA
#include "cublas_gemm_wrapper.hpp"
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

#define MAX_DIM 1024
#define MIN_DIM 32

using namespace std;
using namespace H5;

extern "C" {
    void dgemm_(const char* TRANSA, const char* TRANSB,
                const int* M, const int* N, const int* K,
                const double* ALPHA,
                const double* A, const int* LDA,
                const double* B, const int* LDB,
                const double* BETA,
                double* C, const int* LDC);

    void ldgemm_(const char* TRANSA, const char* TRANSB,
                const int* M, const int* N, const int* K,
                const long double* ALPHA,
                const long double* A, const int* LDA,
                const long double* B, const int* LDB,
                const long double* BETA,
                long double* C, const int* LDC);
}

int main(){

#ifndef DRYRUN
    H5File h5_file = init_hdf5(HDF5_FILENAME);
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
    vector<string> matrix_types;

    // matrix_types.push_back("U(0,1)");
    // matrix_types.push_back("U(-1,1)");
    
    matrix_types.push_back("test1a");
    matrix_types.push_back("test1b");
    matrix_types.push_back("test1c");
    
    matrix_types.push_back("test2a");
    matrix_types.push_back("test2b");
    matrix_types.push_back("test2c");
    
    matrix_types.push_back("U(1,2)");
    matrix_types.push_back("test3");

    for ( const string& matrix_type : matrix_types ){
        for ( int dim = MIN_DIM; dim <= MAX_DIM; dim = 2*dim ){
            
            int m, n, k;
            m = n = k = dim;

            // initialize A and B matrices
#ifdef STOCHASTIC
            int seed = time(nullptr);
#else
            int seed = 0;
#endif
            init_matrices(matrix_type, seed, m, n, k, &A[0], lda, &B[0], ldb);
#ifndef DRYRUN
            write_matrix_to_hdf5(h5_file, &A[0], lda, "A", matrix_type, m, k, seed);
            write_matrix_to_hdf5(h5_file, &B[0], ldb, "B", matrix_type, k, n, seed);
#endif
#ifdef VERBOSE
            printf("A matrix:\n");
            print_matrix_to_stdout_sci(m, k, &A[0], lda);
            printf("B matrix:\n");
            print_matrix_to_stdout_sci(k, n, &B[0], ldb);
#endif
            // call DGEMM
#ifdef EXTENDED_PRECISION_REFERENCE
            alpha_ext = alpha;
            beta_ext = beta;
            copy_matrix(m, k, &A_ext[0], lda, &A[0], lda);
            copy_matrix(m, k, &B_ext[0], ldb, &B[0], ldb);
            init_zero_matrix(m, n, &C_ext[0], ldc);
            ldgemm_(&transpose, &transpose, &m, &n, &k, &alpha_ext, &A_ext[0], &lda, &B_ext[0], &ldb, &beta_ext, &C_ext[0], &ldc);
#ifndef DRYRUN
            write_matrix_to_hdf5(h5_file, &C_ext[0], ldc, DGEMM_IMPL, matrix_type, m, n, seed);
#endif
            abs_matrix(m, k, &A_ext[0], lda);
            abs_matrix(m, k, &B_ext[0], ldb);
            init_zero_matrix(m, n, &C_ext[0], ldc);
            ldgemm_(&transpose, &transpose, &m, &n, &k, &alpha_ext, &A_ext[0], &lda, &B_ext[0], &ldb, &beta_ext, &C_ext[0], &ldc);
#ifndef DRYRUN
            string new_name = "ABS_A_AND_B_";
            new_name = new_name + DGEMM_IMPL;
            write_matrix_to_hdf5(h5_file, &C_ext[0], ldc, new_name, matrix_type, m, n, seed);
#endif
#else
            init_zero_matrix(m, n, &C[0], ldc);
#ifdef CUDA
            cublas_gemm_wrapper(m, n, k, &A[0], lda, &B[0], ldb, &C[0], ldc);
#else
            dgemm_(&transpose, &transpose, &m, &n, &k, &alpha, &A[0], &lda, &B[0], &ldb, &beta, &C[0], &ldc);
#endif
#ifdef VERBOSE
            printf("C matrix:\n");
            print_matrix_to_stdout_sci(m, n, &C[0], ldc);
#endif
#ifndef DRYRUN
            write_matrix_to_hdf5(h5_file, &C[0], ldc, DGEMM_IMPL, matrix_type, m, n, seed);
#endif
#endif
        }
    }
#ifndef DRYRUN
    H5Fclose(h5_file.getId());
#endif
    return 0;
}