#include "init_matrices.hpp"
#include "utils.hpp"
#include "reference_gemm.hpp"
#include "strassen_gemm.hpp"
#include <stdio.h>
#include <stdint.h>
#include <quadmath.h>
#include <fstream>
#include <vector>
#include "blis.h"

int main(){
    int32_t m, k, n, lda, ldb, ldc;
    int32_t i, j, p;
    int32_t min_dim, max_dim;

#ifndef SANITY
    max_dim = 512;
    min_dim = 64;
#else
    min_dim = 2;
    max_dim = 4;
#endif

    m = n = k = lda = ldb = ldc = max_dim;

    std::vector<double> A(m*k);
    std::vector<double> B(k*n);
    std::vector<double> C(m*n);
    std::vector<double> low_prec_buffer(m*n);
#ifndef SANITY
    std::vector<__float128> high_prec_buffer(m*n);
    std::vector<__float128> Aq(m*k);
    std::vector<__float128> Bq(k*n);
    std::vector<__float128> Cq(m*n);
    std::vector<__float128> errors(m*n);

    std::ofstream outfile1("reference_errors.bin", std::ios::binary);
    std::ofstream outfile2("blis_errors.bin", std::ios::binary);
    std::ofstream outfile3("strassen_errors.bin", std::ios::binary);
#endif
    for ( i = min_dim; i <= max_dim; i = 2*i ){
        for ( j = min_dim; j <= max_dim; j = 2*j ){
            for ( p = min_dim; p <= max_dim; p = 2*p ){    

                // initialize A and B matrices
#ifndef SANITY
                init_test1_matrices(&A[0], k, &B[0], lda, i, p, j);
                copy_matrix(&Aq[0], lda, &A[0], lda, i, p);
                copy_matrix(&Bq[0], ldb, &B[0], ldb, p, j);

                // run reference gemm in high precision
                init_zero_matrix(&Cq[0], n, i, j);
                reference_gemm(&Aq[0], lda, &Bq[0], ldb, i, p, j, &Cq[0], ldc);
#else
                init_sanity_matrices(&A[0], k, &B[0], lda, i, p, j);
                printf("A:\n");
                print_matrix_to_stdout(&A[0], lda, i, p);
                printf("B:\n");
                print_matrix_to_stdout(&B[0], ldb, p, j);
#endif
                // run reference gemm
                init_zero_matrix(&C[0], n, i, j);
                reference_gemm(&A[0], lda, &B[0], ldb, i, p, j, &C[0], ldc);
#ifndef SANITY
                // calculate and write errors to file
                copy_matrix(&high_prec_buffer[0], ldc, &C[0], ldc, i, j);
                get_elementwise_errors(&Cq[0], ldc, &high_prec_buffer[0], ldc, i, j, &errors[0], ldc);
                copy_matrix(&low_prec_buffer[0], ldc, &errors[0], ldc, i, j);
                write_matrix_bytes_to_file(outfile1, &low_prec_buffer[0], ldc, i, j);
#else
                printf("C (reference):\n");
                print_matrix_to_stdout(&C[0], ldc, i, j);
#endif
                // run blis gemm and calculate errors
                char transpose = 'N';
                double alpha = 1.0;
                double beta = 0.0;
                init_zero_matrix(&C[0], n, i, j);
                dgemm_(&transpose, &transpose, &i, &j, &p, &alpha, &A[0], &lda, &B[0], &ldb, &beta, &C[0], &ldc);
#ifndef SANITY
                // calculate and write errors to file
                copy_matrix(&high_prec_buffer[0], ldc, &C[0], ldc, i, j);
                get_elementwise_errors(&Cq[0], ldc, &high_prec_buffer[0], ldc, i, j, &errors[0], ldc);
                copy_matrix(&low_prec_buffer[0], ldc, &errors[0], ldc, i, j);
                write_matrix_bytes_to_file(outfile2, &low_prec_buffer[0], ldc, i, j);
#else
                printf("C (blis):\n");
                print_matrix_to_stdout(&C[0], ldc, i, j);
#endif
                // run strassen gemm
                init_zero_matrix(&C[0], n, i, j);
                strassen_gemm(&A[0], lda, &B[0], ldb, i, p, j, &C[0], ldc);
#ifndef SANITY
                // calculate and write errors to file
                copy_matrix(&high_prec_buffer[0], ldc, &C[0], ldc, i, j);
                get_elementwise_errors(&Cq[0], ldc, &high_prec_buffer[0], ldc, i, j, &errors[0], ldc);
                copy_matrix(&low_prec_buffer[0], ldc, &errors[0], ldc, i, j);
                write_matrix_bytes_to_file(outfile3, &low_prec_buffer[0], ldc, i, j);
#else
                printf("C (strassen):\n");
                print_matrix_to_stdout(&C[0], ldc, i, j);
#endif
            }    
        }    
    }

#ifndef SANITY
    outfile1.close();
    outfile2.close();
    outfile3.close();
#endif

    return 0;
}