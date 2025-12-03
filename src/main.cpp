#include "init_matrices.hpp"
#include "utils.hpp"
#include "reference_gemm.hpp"
#include "strassen_gemm.hpp"
#include "cublas_gemm_wrapper.hpp"
#include <stdio.h>
#include <stdint.h>
#include <fstream>
#include <vector>
#include <dlfcn.h>

typedef void (*dgemm_t)(char*, char*, int32_t*, int32_t*, int32_t*,
                        double*, double*, int32_t*, double*, int32_t*, double*, double*, int32_t*);

int main(){
    int32_t m, k, n, lda, ldb, ldc;
    int32_t dim, min_dim, max_dim;

    char transpose = 'N';
    double alpha = 1.0;
    double beta = 0.0;
    
    void *handle_netlibblas, *handle_blis, *handle_openblas;
    dgemm_t dgemm_netlibblas, dgemm_blis, dgemm_openblas;

    handle_netlibblas = dlopen("librefblas.so", RTLD_LAZY);
    dgemm_netlibblas = (dgemm_t)dlsym(handle_netlibblas, "dgemm_");

    handle_blis = dlopen("libblis.so", RTLD_LAZY);
    dgemm_blis = (dgemm_t)dlsym(handle_blis, "dgemm_");

    handle_openblas = dlopen("libopenblas.so", RTLD_LAZY);
    dgemm_openblas = (dgemm_t)dlsym(handle_openblas, "dgemm_");

#ifndef SANITY
    max_dim = 64;
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
    double max_componentwise_abs_error;
    std::vector<double> C_ref(m*n);
    std::vector<double> errors(m*n);

    std::vector<std::string> matrix_types;
    matrix_types.push_back("test1");
    // matrix_types.push_back("test3");
    // matrix_types.push_back("U_0_1");
    // matrix_types.push_back("U_-1_1");

    for ( const auto& matrix_type : matrix_types ){

        std::ofstream netlibblas_out_bin("results/" + matrix_type + "/all_componentwise_errors.netlibblas.bin", std::ios::binary);
        std::ofstream netlibblas_out_txt("results/" + matrix_type + "/max_componentwise_abs_errors.netlibblas.txt");
        netlibblas_out_txt << "m,n,k,max_componentwise_abs_error" << std::endl;

        std::ofstream blis_out_bin("results/" + matrix_type + "/all_componentwise_errors.blis.bin", std::ios::binary);
        std::ofstream blis_out_txt("results/" + matrix_type + "/max_componentwise_abs_errors.blis.txt");
        blis_out_txt << "m,n,k,max_componentwise_abs_error" << std::endl;

        std::ofstream openblas_out_bin("results/" + matrix_type + "/all_componentwise_errors.openblas.bin", std::ios::binary);
        std::ofstream openblas_out_txt("results/" + matrix_type + "/max_componentwise_abs_errors.openblas.txt");
        openblas_out_txt << "m,n,k,max_componentwise_abs_error" << std::endl;

        std::ofstream strassen_out_bin("results/" + matrix_type + "/all_componentwise_errors.strassen.bin", std::ios::binary);
        std::ofstream strassen_out_txt("results/" + matrix_type + "/max_componentwise_abs_errors.strassen.txt");
        strassen_out_txt << "m,n,k,max_componentwise_abs_error" << std::endl;

        std::ofstream cublas_out_bin("results/" + matrix_type + "/all_componentwise_errors.cublas.bin", std::ios::binary);
        std::ofstream cublas_out_txt("results/" + matrix_type + "/max_componentwise_abs_errors.cublas.txt");
        cublas_out_txt << "m,n,k,max_componentwise_abs_error" << std::endl;

#ifdef USE_OZIMMU
        std::ofstream ozimmu_out_bin("results/" + matrix_type + "/all_componentwise_errors.ozimmu.bin", std::ios::binary);
        std::ofstream ozimmu_out_txt("results/" + matrix_type + "/max_componentwise_abs_errors.ozimmu.txt");
        ozimmu_out_txt << "m,n,k,max_componentwise_abs_error" << std::endl;
#endif
        for ( dim = min_dim; dim <= max_dim; dim = 2*dim ){

            // initialize A and B matrices
            init_matrices(matrix_type, &A[0], k, &B[0], lda, dim, dim, dim);

            // run reference gemm
            init_zero_matrix(&C_ref[0], n, dim, dim);
            reference_gemm(&A[0], lda, &B[0], ldb, dim, dim, dim, &C_ref[0], ldc);

            // run netlibblas gemm and calculate errors
            init_zero_matrix(&C[0], n, dim, dim);
            dgemm_netlibblas(&transpose, &transpose, &dim, &dim, &dim, &alpha, &A[0], &lda, &B[0], &ldb, &beta, &C[0], &ldc);
            max_componentwise_abs_error = get_componentwise_errors(&C_ref[0], ldc, &C[0], ldc, dim, dim, &errors[0], ldc);
            write_matrix_bytes_to_file(netlibblas_out_bin, &errors[0], ldc, dim, dim);
            netlibblas_out_txt << dim << "," << dim << "," << dim << "," << max_componentwise_abs_error << std::endl;

            // run blis gemm and calculate errors
            init_zero_matrix(&C[0], n, dim, dim);
            dgemm_blis(&transpose, &transpose, &dim, &dim, &dim, &alpha, &A[0], &lda, &B[0], &ldb, &beta, &C[0], &ldc);
            max_componentwise_abs_error = get_componentwise_errors(&C_ref[0], ldc, &C[0], ldc, dim, dim, &errors[0], ldc);
            write_matrix_bytes_to_file(blis_out_bin, &errors[0], ldc, dim, dim);
            blis_out_txt << dim << "," << dim << "," << dim << "," << max_componentwise_abs_error << std::endl;

            // run openblas gemm and calculate errors
            init_zero_matrix(&C[0], n, dim, dim);
            dgemm_openblas(&transpose, &transpose, &dim, &dim, &dim, &alpha, &A[0], &lda, &B[0], &ldb, &beta, &C[0], &ldc);
            max_componentwise_abs_error = get_componentwise_errors(&C_ref[0], ldc, &C[0], ldc, dim, dim, &errors[0], ldc);
            write_matrix_bytes_to_file(openblas_out_bin, &errors[0], ldc, dim, dim);
            openblas_out_txt << dim << "," << dim << "," << dim << "," << max_componentwise_abs_error << std::endl;

            // run strassen gemm and calculate errors
            init_zero_matrix(&C[0], n, dim, dim);
            strassen_gemm(&A[0], lda, &B[0], ldb, dim, dim, dim, &C[0], ldc);
            max_componentwise_abs_error = get_componentwise_errors(&C_ref[0], ldc, &C[0], ldc, dim, dim, &errors[0], ldc);
            write_matrix_bytes_to_file(strassen_out_bin, &errors[0], ldc, dim, dim);
            strassen_out_txt << dim << "," << dim << "," << dim << "," << max_componentwise_abs_error << std::endl;

            // run cublas gemm and calculate errors
            init_zero_matrix(&C[0], n, dim, dim);
            cublas_gemm_wrapper(&A[0], lda, &B[0], ldb, dim, dim, dim, &C[0], ldc);
            max_componentwise_abs_error = get_componentwise_errors(&C_ref[0], ldc, &C[0], ldc, dim, dim, &errors[0], ldc);
            write_matrix_bytes_to_file(cublas_out_bin, &errors[0], ldc, dim, dim);
            cublas_out_txt << dim << "," << dim << "," << dim << "," << max_componentwise_abs_error << std::endl;

#ifdef USE_OZIMMU
            // run ozimmu gemm and calculate errors
            init_zero_matrix(&C[0], n, dim, dim);
            cublas_ozimmu_wrapper(&A[0], lda, &B[0], ldb, dim, dim, dim, &C[0], ldc);
            max_componentwise_abs_error = get_componentwise_errors(&C_ref[0], ldc, &C[0], ldc, dim, dim, &errors[0], ldc);
            write_matrix_bytes_to_file(ozimmu_out_bin, &errors[0], ldc, dim, dim);
            ozimmu_out_txt << dim << "," << dim << "," << dim << "," << max_componentwise_abs_error << std::endl;
#endif
        }
        netlibblas_out_bin.close();
        netlibblas_out_txt.close();
        blis_out_bin.close();
        blis_out_txt.close();
        openblas_out_bin.close();
        openblas_out_txt.close();
        strassen_out_bin.close();
        strassen_out_txt.close();
        cublas_out_bin.close();
        cublas_out_txt.close();
#ifdef USE_OZIMMU
        ozimmu_out_bin.close();
        ozimmu_out_txt.close();
#endif
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

        // run netlibblas gemm
        init_zero_matrix(&C[0], n, dim, dim);
        dgemm_netlibblas(&transpose, &transpose, &dim, &dim, &dim, &alpha, &A[0], &lda, &B[0], &ldb, &beta, &C[0], &ldc);
        printf("C (netlibblas):\n");
        print_matrix_to_stdout(&C[0], ldc, dim, dim);

        // run blis gemm
        init_zero_matrix(&C[0], n, dim, dim);
        dgemm_blis(&transpose, &transpose, &dim, &dim, &dim, &alpha, &A[0], &lda, &B[0], &ldb, &beta, &C[0], &ldc);
        printf("C (blis):\n");
        print_matrix_to_stdout(&C[0], ldc, dim, dim);

        // run openblas gemm
        init_zero_matrix(&C[0], n, dim, dim);
        dgemm_openblas(&transpose, &transpose, &dim, &dim, &dim, &alpha, &A[0], &lda, &B[0], &ldb, &beta, &C[0], &ldc);
        printf("C (openblas):\n");
        print_matrix_to_stdout(&C[0], ldc, dim, dim);

        // run strassen gemm
        init_zero_matrix(&C[0], n, dim, dim);
        strassen_gemm(&A[0], lda, &B[0], ldb, dim, dim, dim, &C[0], ldc);
        printf("C (strassen):\n");
        print_matrix_to_stdout(&C[0], ldc, dim, dim);

        // run cublas gemm
        init_zero_matrix(&C[0], n, dim, dim);
        cublas_gemm_wrapper(&A[0], lda, &B[0], ldb, dim, dim, dim, &C[0], ldc);
        printf("C (cublas):\n");
        print_matrix_to_stdout(&C[0], ldc, dim, dim);

#ifdef USE_OZIMMU
        // run ozimmu gemm
        init_zero_matrix(&C[0], n, dim, dim);
        cublas_ozimmu_wrapper(&A[0], lda, &B[0], ldb, dim, dim, dim, &C[0], ldc);
        printf("C (ozIMMU):\n");
        print_matrix_to_stdout(&C[0], ldc, dim, dim);
#endif
    }

#endif

    return 0;
}