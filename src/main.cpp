#include "init_matrices.hpp"
#include "utils.hpp"
#include "cublas_gemm_wrapper.hpp"
#include <stdio.h>
#include <stdint.h>
#include <fstream>
#include <vector>
#include <ctime>
#include <iomanip>
#include <random>
#include <H5Cpp.h>

using namespace std;
using namespace H5;

#define HDF5_FILENAME "grading_dgemm.h5"

extern "C" {
    void dgemm_(const char* TRANSA, const char* TRANSB,
                const int* M, const int* N, const int* K,
                const double* ALPHA,
                const double* A, const int* LDA,
                const double* B, const int* LDB,
                const double* BETA,
                double* C, const int* LDC);
}

#ifndef DGEMM_IMPL
#define DGEMM_IMPL "unknown"
#endif

int main(){

    H5File h5_file = init_hdf5(HDF5_FILENAME);

    char transpose = 'N';
    double alpha = 1.0;
    double beta = 0.0;

    int32_t max_m, max_k, max_n, lda, ldb, ldc;
    int32_t min_dim, max_dim;

#ifdef SANITY
    min_dim = max_dim = 4;
#else
    max_dim = 2048;
    min_dim = 32;
#endif
    max_m = max_k = max_n = lda = ldb = ldc = max_dim;

    vector<double> A(max_m*max_k);
    vector<double> B(max_k*max_n);
    vector<double> C(max_m*max_n);
    vector<string> matrix_types;
#ifdef SANITY
#else
    // matrix_types.push_back("sanity");
    matrix_types.push_back("test1");
    matrix_types.push_back("test3");
    matrix_types.push_back("U_0_1");
    matrix_types.push_back("U_-1_1");
#endif
    for ( const string& matrix_type : matrix_types ){
        for ( int32_t dim = min_dim; dim <= max_dim; dim = 2*dim ){
            
            int32_t m, n, k;
            m = n = k = dim;

            // initialize A and B matrices
#ifdef STOCHASTIC
            int32_t seed = time(nullptr);
#else
            int32_t seed = 0;
#endif
            init_matrices(matrix_type, seed, m, n, k, &A[0], lda, &B[0], ldb);
            write_matrix_to_hdf5(h5_file, &A[0], lda, "A", matrix_type, m, k, seed);
            write_matrix_to_hdf5(h5_file, &B[0], ldb, "B", matrix_type, k, n, seed);
#ifdef VERBOSE
            printf("A matrix:\n");
            print_matrix_to_stdout_dec(m, k, &A[0], lda);
            printf("B matrix:\n");
            print_matrix_to_stdout_dec(k, n, &B[0], ldb);
#endif
            // call DGEMM
            init_zero_matrix(m, n, &C[0], ldc);
#ifdef CUDA
            cublas_gemm_wrapper(m, n, k, &A[0], lda, &B[0], ldb, &C[0], ldc);
#else
            dgemm_(&transpose, &transpose, &m, &n, &k, &alpha, &A[0], &lda, &B[0], &ldb, &beta, &C[0], &ldc);
#endif
#ifdef VERBOSE
            printf("C matrix:\n");
            print_matrix_to_stdout_dec(m, n, &C[0], ldc);
#endif
            write_matrix_to_hdf5(h5_file, &C[0], ldc, DGEMM_IMPL, matrix_type, m, n, seed);
        }
    }

    H5Fclose(h5_file.getId());

    return 0;
}