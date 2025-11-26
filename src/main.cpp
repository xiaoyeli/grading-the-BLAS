#include "init_matrices.hpp"
#include "utils.hpp"
#include "reference_gemm.hpp"
#include "strassen_gemm.hpp"
#include <stdio.h>
#include <stdint.h>
#include <fstream>
#include <vector>
#include "blis.h"

int main(){
    int32_t m, k, n, lda, ldb, ldc;
    int32_t dim, min_dim, max_dim;

#ifndef SANITY
    max_dim = 256;
    min_dim = 32;
#else
    min_dim = 2;
    max_dim = 4;
#endif

    m = n = k = lda = ldb = ldc = max_dim;

    std::vector<double> A(m*k);
    std::vector<double> B(k*n);
    std::vector<double> C(m*n);

#ifndef SANITY
    double max_componentwise_error;
    std::vector<double> C_ref(m*n);
    std::vector<double> errors(m*n);

    std::vector<std::string> matrix_types = {"U_0_1", "U_-1_1", "test1"};

    for ( const auto& matrix_type : matrix_types ){

        std::ofstream blis_out_bin("results/" + matrix_type + "/all_componentwise_errors.blis.bin", std::ios::binary);
        std::ofstream blis_out_txt("results/" + matrix_type + "/max_componentwise_errors.blis.txt");
        std::ofstream strassen_out_bin("results/" + matrix_type + "/all_componentwise_errors.strassen.bin", std::ios::binary);
        std::ofstream strassen_out_txt("results/" + matrix_type + "/max_componentwise_errors.strassen.txt");
        blis_out_txt << "m,n,k,max_componentwise_error" << std::endl;
        strassen_out_txt << "m,n,k,max_componentwise_error" << std::endl;

        for ( dim = min_dim; dim <= max_dim; dim = 2*dim ){

            // initialize A and B matrices
            init_matrices(matrix_type, &A[0], k, &B[0], lda, dim, dim, dim);

            // run reference gemm
            init_zero_matrix(&C_ref[0], n, dim, dim);
            reference_gemm(&A[0], lda, &B[0], ldb, dim, dim, dim, &C_ref[0], ldc);

            // run blis gemm and calculate errors
            char transpose = 'N';
            double alpha = 1.0;
            double beta = 0.0;
            init_zero_matrix(&C[0], n, dim, dim);
            dgemm_(&transpose, &transpose, &dim, &dim, &dim, &alpha, &A[0], &lda, &B[0], &ldb, &beta, &C[0], &ldc);

            // calculate and write errors to file
            max_componentwise_error = get_componentwise_errors(&C_ref[0], ldc, &C[0], ldc, dim, dim, &errors[0], ldc);
            write_matrix_bytes_to_file(blis_out_bin, &errors[0], ldc, dim, dim);
            blis_out_txt << dim << "," << dim << "," << dim << "," << max_componentwise_error << std::endl;

            // run strassen gemm
            init_zero_matrix(&C[0], n, dim, dim);
            strassen_gemm(&A[0], lda, &B[0], ldb, dim, dim, dim, &C[0], ldc);

            // calculate and write errors to file
            max_componentwise_error = get_componentwise_errors(&C_ref[0], ldc, &C[0], ldc, dim, dim, &errors[0], ldc);
            write_matrix_bytes_to_file(strassen_out_bin, &errors[0], ldc, dim, dim);
            strassen_out_txt << dim << "," << dim << "," << dim << "," << max_componentwise_error << std::endl;
        }
        blis_out_bin.close();
        blis_out_txt.close();
        strassen_out_bin.close();
        strassen_out_txt.close();
    }

#else
    for ( dim = min_dim; dim <= max_dim; dim = 2*dim ){

        // initialize A and B matrices
        init_sanity_matrices(&A[0], k, &B[0], lda, dim, dim, dim);
        printf("A:\n");
        print_matrix_to_stdout(&A[0], lda, dim, dim);
        printf("B:\n");
        print_matrix_to_stdout(&B[0], ldb, dim, dim);

        // run reference gemm
        init_zero_matrix(&C[0], n, dim, dim);
        reference_gemm(&A[0], lda, &B[0], ldb, dim, dim, dim, &C[0], ldc);
        printf("C (reference):\n");
        print_matrix_to_stdout(&C[0], ldc, dim, dim);

        // run blis gemm
        char transpose = 'N';
        double alpha = 1.0;
        double beta = 0.0;
        init_zero_matrix(&C[0], n, dim, dim);
        dgemm_(&transpose, &transpose, &dim, &dim, &dim, &alpha, &A[0], &lda, &B[0], &ldb, &beta, &C[0], &ldc);
        printf("C (blis):\n");
        print_matrix_to_stdout(&C[0], ldc, dim, dim);

        // run strassen gemm
        init_zero_matrix(&C[0], n, dim, dim);
        strassen_gemm(&A[0], lda, &B[0], ldb, dim, dim, dim, &C[0], ldc);
        printf("C (strassen):\n");
        print_matrix_to_stdout(&C[0], ldc, dim, dim);
    }

#endif

    return 0;
}