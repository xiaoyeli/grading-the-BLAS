#include "init_matrices.hpp"
#include "utils.hpp"
#include "reference_gemm.hpp"
#include "cublas_gemm_wrapper.hpp"
#include <stdio.h>
#include <stdint.h>
#include <fstream>
#include <vector>

extern "C" {
    void dgemm_(const char* TRANSA, const char* TRANSB,
                const int* M, const int* N, const int* K,
                const double* ALPHA,
                const double* A, const int* LDA,
                const double* B, const int* LDB,
                const double* BETA,
                double* C, const int* LDC);
}

#ifdef USE_NETLIBBLAS
#define LIB_NAME "netlibblas"
#endif

#ifdef USE_BLIS
#include <blis.h>
#define LIB_NAME "blis"
#endif

#ifdef USE_OPENBLAS
#define LIB_NAME "openblas"
#endif

#ifdef USE_STRASSEN
#define LIB_NAME "strassen"
#endif

#ifdef USE_CUBLAS
#define LIB_NAME "cublas"
#endif

#ifdef USE_OZIMMU
#define LIB_NAME "ozimmu"
#endif

int main(){
    int32_t m, k, n, lda, ldb, ldc;
    int32_t dim, min_dim, max_dim;
    int32_t n_zeros;

    char transpose = 'N';
    double alpha = 1.0;
    double beta = 0.0;

#ifdef SANITY
    min_dim = max_dim = 4;
#else
    max_dim = 8192;
    min_dim = 32;
#endif
    m = n = k = lda = ldb = ldc = max_dim;

    std::vector<double> A(m*k);
    std::vector<double> B(k*n);
    std::vector<double> C(m*n);

    double max_componentwise_abs_error;
    std::vector<double> C_ref(m*n);
    std::vector<double> errors(m*n);

    std::vector<std::string> matrix_types;
#ifdef SANITY
    matrix_types.push_back("sanity");
#else
    matrix_types.push_back("test1");
    matrix_types.push_back("test3");
    matrix_types.push_back("U_0_1");
    matrix_types.push_back("U_-1_1");
#endif
    for ( const auto& matrix_type : matrix_types ){

#ifndef LIB_NAME
        std::ofstream out_bin("results/" + matrix_type + "/reference_values.bin", std::ios::binary);
#else
        std::ifstream in_bin("results/" + matrix_type + "/reference_values.bin", std::ios::binary);
        std::ofstream out_bin("results/" + matrix_type + "/all_componentwise_errors." + LIB_NAME + ".bin", std::ios::binary);
        std::ofstream out_txt("results/" + matrix_type + "/max_componentwise_abs_errors." + LIB_NAME + ".txt");
        out_txt << "m,n,k,max_componentwise_abs_error,n_zeros" << std::endl;
#endif
        for ( dim = min_dim; dim <= max_dim; dim = 2*dim ){

            // initialize A and B matrices
            init_matrices(matrix_type, &A[0], k, &B[0], lda, dim, dim, dim);
#ifdef VERBOSE
            printf("A matrix:\n");
            print_matrix_to_stdout(&A[0], lda, dim, dim);
            printf("B matrix:\n");
            print_matrix_to_stdout(&B[0], ldb, dim, dim);
#endif
            // get reference values, either calculating them or reading them from file
#ifndef LIB_NAME
            init_zero_matrix(&C_ref[0], n, dim, dim);
            reference_gemm(&A[0], lda, &B[0], ldb, dim, dim, dim, &C_ref[0], ldc);
            write_matrix_bytes_to_file(out_bin, &C_ref[0], ldc, dim, dim);
#else
            read_matrix_bytes_from_file(in_bin, &C_ref[0], ldc, dim, dim);
#ifdef VERBOSE
            printf("C reference matrix:\n");
            print_matrix_to_stdout(&C_ref[0], ldc, dim, dim);
#endif
            // call DGEMM for the library under test
            init_zero_matrix(&C[0], n, dim, dim);
#ifdef CUDA
            cublas_gemm_wrapper(&A[0], lda, &B[0], ldb, dim, dim, dim, &C[0], ldc);
#else
            dgemm_(&transpose, &transpose, &dim, &dim, &dim, &alpha, &A[0], &lda, &B[0], &ldb, &beta, &C[0], &ldc);
#endif
#ifdef VERBOSE
            printf("C matrix:\n");
            print_matrix_to_stdout(&C[0], ldc, dim, dim);
#endif
            // calculate and save errors
            postprocess(&C_ref[0], ldc, &C[0], ldc, dim, dim, &errors[0], ldc, &max_componentwise_abs_error, &n_zeros);
            write_matrix_bytes_to_file(out_bin, &errors[0], ldc, dim, dim);
            out_txt << dim << "," << dim << "," << dim << "," << max_componentwise_abs_error << "," << n_zeros << std::endl;
#endif
        }
        out_bin.close();
#ifdef LIB_NAME
        in_bin.close();
        out_txt.close();
#endif
    }

    return 0;
}