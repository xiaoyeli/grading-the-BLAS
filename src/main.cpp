#include "dgemm_wrapper.hpp"
#include "init_matrices.hpp"
#include "utils.hpp"
#include <vector>
#include <cmath>

#define OUTPUT_DIR "results/raw_data/"

#if defined(VERBOSE)
#include <stdio.h>
#endif

#if REFERENCE_GEMM > 0
#include "reference_gemm.hpp"
#endif

int main(){

    std::string dgemm_impl_name = DGEMM_IMPL_NAME;
#if defined(OUTER_SCALING)
    dgemm_impl_name += "_w_outer_scaling";
#endif
#if defined(INNER_SCALING)
    dgemm_impl_name += "_w_inner_scaling";
#endif

    char transa = 'N';
    char transb = 'N';
    double alpha = 1.0;
    double beta = 0.0;
    
    int max_m, max_k, max_n, lda, ldb, ldc;
    max_m = max_k = max_n = lda = ldb = ldc = 1 << LOG2_MAX_DIM;
    
    std::vector<double> A(max_m*max_k);
    std::vector<double> B(max_k*max_n);
    std::vector<double> C(max_m*max_n);
    std::vector<double> rel_errors(max_m*max_n);

#if defined(INNER_PERMUTE) || defined(OUTER_PERMUTE)
    std::vector<double> C_ref(max_m*max_n);    
#elif REFERENCE_GEMM == 1
    std::vector<double> C_ref(max_m*max_n);
    std::vector<double> C_ref_abs(max_m*max_n);
    std::vector<double> f_of_n(max_m*max_n);
    #elif REFERENCE_GEMM == 2
    std::vector<long double> C_ref(max_m*max_n);
    std::vector<long double> C_ref_abs(max_m*max_n);
    std::vector<double> f_of_n(max_m*max_n);
#endif

    for ( int d = LOG2_MIN_DIM; d <= LOG2_MAX_DIM; ++d ){
        
        int m, n, k;
        m = n = k = 1 << d;

        for ( int s = LOG2_MIN_SCALING_FACTOR; s <= LOG2_MAX_SCALING_FACTOR; s += LOG2_INC_SCALING_FACTOR ){
        
            int seed = 0;
            init_matrices(INPUT_TYPE, s, seed, m, n, k, &A[0], lda, &B[0], ldb);

#if defined(VERBOSE)
            printf("A matrix:\n");
            print_matrix_to_stdout_sci(m, k, &A[0], lda);
            printf("B matrix:\n");
            print_matrix_to_stdout_sci(k, n, &B[0], ldb);
#endif

#if defined(INNER_PERMUTE) || defined(OUTER_PERMUTE)
            init_zero_matrix(m, n, &C_ref[0], ldc);
            init_zero_matrix(m, n, &rel_errors[0], ldc);
            reference_gemm(&transa, &transb, &m, &n, &k, &alpha, &A[0], &lda, &B[0], &ldb, &beta, &C_ref[0], &ldc);
#endif

#if defined(OUTER_PERMUTE)
            int perm_seed = 0;
            for ( int i = 0; i < OUTER_PERMUTE; ++i ){
                
                init_matrices(INPUT_TYPE, s, seed, m, n, k, &A[0], lda, &B[0], ldb);
                
                perm_seed++;
                permute_rows(perm_seed, m, k, &A[0], lda);
                permute_cols(perm_seed + OUTER_PERMUTE, k, n, &B[0], ldb);
#elif defined(INNER_PERMUTE)
            int perm_seed = 0;
            for ( int i = 0; i < INNER_PERMUTE; ++i ){
                
                init_matrices(INPUT_TYPE, s, seed, m, n, k, &A[0], lda, &B[0], ldb);
                
                perm_seed++;
                permute_cols(perm_seed, m, k, &A[0], lda);
                permute_rows(perm_seed, k, n, &B[0], ldb);
#endif

            init_zero_matrix(m, n, &C[0], ldc);
            dgemm_wrapper(transa, transb, m, n, k, alpha, &A[0], lda, &B[0], ldb, beta, &C[0], ldc);

#if defined(OUTER_PERMUTE)
            permute_rows_inv(seed, m, k, &C[0], ldc);
            permute_cols_inv(seed + OUTER_PERMUTE, k, n, &C[0], ldc);
#endif

#if defined(VERBOSE)
            printf("C matrix:\n");
            print_matrix_to_stdout_sci(m, n, &C[0], ldc);
#endif

#if defined(OUTER_PERMUTE) || defined(INNER_PERMUTE)
                calculate_relative_errors(m, n, &C[0], ldc, &C_ref[0], ldc, &rel_errors[0], ldc);
            }

            std::string input_type = INPUT_TYPE;
#if defined(OUTER_PERMUTE)
            input_type += "_w_outer_permute";
#endif
#if defined(INNER_PERMUTE)
            input_type += "_w_inner_permute";
#endif
            std::string path_string = OUTPUT_DIR + construct_metadata_path(dgemm_impl_name, input_type, m, n, k, s, seed, &rel_errors[0]);
            dump_matrix(m, n, &rel_errors[0], ldc, path_string + "/rel_errors.bin");
            display_distribution_stats(m, n, &rel_errors[0], ldc, path_string + "/rel_errors.bin");

#else
            std::string path_string = OUTPUT_DIR + construct_metadata_path("C_ref", INPUT_TYPE, m, n, k, s, seed, &C_ref[0]);
            if ( ! load_matrix(path_string + "/values.bin", m, n, &C_ref[0], ldc ) ){
                init_zero_matrix(m, n, &C_ref[0], ldc);
                reference_gemm(&transa, &transb, &m, &n, &k, &alpha, &A[0], &lda, &B[0], &ldb, &beta, &C_ref[0], &ldc);
                dump_matrix(m, n, &C_ref[0], ldc, path_string + "/values.bin");
            }

            path_string = OUTPUT_DIR + construct_metadata_path("C_ref_abs", INPUT_TYPE, m, n, k, s, seed, &C_ref[0]);
            if ( ! load_matrix(path_string + "/values.bin", m, n, &C_ref_abs[0], ldc ) ){
                abs_matrix(m, k, &A[0], lda);
                abs_matrix(m, k, &B[0], ldb);
                init_zero_matrix(m, n, &C_ref_abs[0], ldc);
                reference_gemm(&transa, &transb, &m, &n, &k, &alpha, &A[0], &lda, &B[0], &ldb, &beta, &C_ref_abs[0], &ldc);
                dump_matrix(m, n, &C_ref_abs[0], ldc, path_string + "/values.bin");
            }

            init_zero_matrix(m, n, &rel_errors[0], ldc);
            init_zero_matrix(m, n, &f_of_n[0], ldc);
            calculate_relative_errors(m, n, &C[0], ldc, &C_ref[0], ldc, &rel_errors[0], ldc);
            calculate_f_of_n(m, n, &C[0], ldc, &C_ref[0], ldc, &C_ref_abs[0], ldc, &f_of_n[0], ldc);
            path_string = OUTPUT_DIR + construct_metadata_path(dgemm_impl_name, INPUT_TYPE, m, n, k, s, seed, &C_ref[0]);
            dump_matrix(m, n, &rel_errors[0], ldc, path_string + "/rel_errors.bin");
            dump_matrix(m, n, &f_of_n[0], ldc, path_string + "/f_of_n.bin");
            display_distribution_stats(m, n, &rel_errors[0], ldc, path_string + "/rel_errors.bin");
#endif
        }
    }
    return 0;
}