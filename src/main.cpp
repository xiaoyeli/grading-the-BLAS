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
    std::string path_string;

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
#if defined(DUMP_A_B)
            // dump A and B; note that dump_matrix overwrites every time
            path_string = OUTPUT_DIR + construct_metadata_path("A", INPUT_TYPE, m, n, k, s);
            dump_matrix(m, k, &A[0], lda, path_string + "/values.fp64");
            path_string = OUTPUT_DIR + construct_metadata_path("B", INPUT_TYPE, m, n, k, s);
            dump_matrix(k, n, &B[0], ldb, path_string + "/values.fp64");
#endif

#if defined(INNER_PERMUTE) || defined(OUTER_PERMUTE)
            // calculate ground truth before permuting A and B using the DGEMM implementation under test
            init_zero_matrix(m, n, &C_ref[0], ldc);
            init_zero_matrix(m, n, &rel_errors[0], ldc);
            dgemm_wrapper(transa, transb, m, n, k, alpha, &A[0], lda, &B[0], ldb, beta, &C_ref[0], ldc);
#endif

#if defined(OUTER_PERMUTE)
            // set up OUTER_PERMUTE iterations
            int perm_seed = 0;
            for ( int i = 0; i < OUTER_PERMUTE; ++i ){
                
                init_matrices(INPUT_TYPE, s, seed, m, n, k, &A[0], lda, &B[0], ldb);
                
                perm_seed++;
                permute_rows(perm_seed, m, k, &A[0], lda);
                permute_cols(perm_seed + OUTER_PERMUTE, k, n, &B[0], ldb);
#elif defined(INNER_PERMUTE)
            // set up INNER_PERMUTE iterations
            int perm_seed = 0;
            for ( int i = 0; i < INNER_PERMUTE; ++i ){
                
                init_matrices(INPUT_TYPE, s, seed, m, n, k, &A[0], lda, &B[0], ldb);
                
                perm_seed++;
                permute_cols(perm_seed, m, k, &A[0], lda);
                permute_rows(perm_seed, k, n, &B[0], ldb);
#endif
                // call the DGEMM implementation under test
                init_zero_matrix(m, n, &C[0], ldc);
                dgemm_wrapper(transa, transb, m, n, k, alpha, &A[0], lda, &B[0], ldb, beta, &C[0], ldc);

#if defined(OUTER_PERMUTE)
                // for OUTER_PERMUTE, unpermute after DGEMM
                permute_rows_inv(seed, m, k, &C[0], ldc);
                permute_cols_inv(seed + OUTER_PERMUTE, k, n, &C[0], ldc);
#endif

#if defined(VERBOSE)
                printf("C matrix:\n");
                print_matrix_to_stdout_sci(m, n, &C[0], ldc);
#endif

#if defined(OUTER_PERMUTE) || defined(INNER_PERMUTE)
                // calculate componentwise error and keeping the maximum over all permutations tested
                calculate_relative_errors(m, n, &C[0], ldc, &C_ref[0], ldc, &rel_errors[0], ldc);

#if defined(DUMP_C)
                // dump ALL of the C matrices computed over all permutations tested
                std::string full_input_type_string = INPUT_TYPE;
#if defined(OUTER_PERMUTE)
                full_input_type_string += "_w_outer_permute_" + std::to_string(perm_seed);
#elif defined(INNER_PERMUTE)
                full_input_type_string += "_w_inner_permute_" + std::to_string(perm_seed);
#endif
                path_string = OUTPUT_DIR + construct_metadata_path(dgemm_impl_name, full_input_type_string, m, n, k, s);
                dump_matrix(m, n, &C[0], ldc, path_string + "/values.fp64");
#endif
            }

            // dump the maximum componentwise error with respect to A*B using the implementation under test over all permutations tested
            std::string full_input_type_string = INPUT_TYPE;
#if defined(OUTER_PERMUTE)
            full_input_type_string += "_w_outer_permute_" + std::to_string(perm_seed);
#elif defined(INNER_PERMUTE)
            full_input_type_string += "_w_inner_permute_" + std::to_string(perm_seed);
#endif
            path_string = OUTPUT_DIR + construct_metadata_path(dgemm_impl_name, full_input_type_string, m, n, k, s);
            dump_matrix(m, n, &rel_errors[0], ldc, path_string + "/rel_errors.fp64");
            display_distribution_stats(m, n, &rel_errors[0], ldc, path_string + "/rel_errors.fp64");

#else

#if defined(DUMP_C)
            path_string = OUTPUT_DIR + construct_metadata_path(dgemm_impl_name, INPUT_TYPE, m, n, k, s);
            dump_matrix(m, n, &C[0], ldc, path_string + "/values.fp64");
#endif
            // either load a cached C_ref or compute from scratch and cache for later
            path_string = OUTPUT_DIR + construct_metadata_path("C_ref", INPUT_TYPE, m, n, k, s);
            if ( ! load_matrix(path_string + "/values.fp" + std::to_string(sizeof(C_ref[0])*8), m, n, &C_ref[0], ldc ) ){
                init_zero_matrix(m, n, &C_ref[0], ldc);
                reference_gemm(&transa, &transb, &m, &n, &k, &alpha, &A[0], &lda, &B[0], &ldb, &beta, &C_ref[0], &ldc);
                dump_matrix(m, n, &C_ref[0], ldc, path_string + "/values.fp" + std::to_string(sizeof(C_ref[0])*8));
            }

            // either load a cached C_ref_abs or compute from scratch and cache for later
            path_string = OUTPUT_DIR + construct_metadata_path("C_ref_abs", INPUT_TYPE, m, n, k, s);
            if ( ! load_matrix(path_string + "/values.fp" + std::to_string(sizeof(C_ref_abs[0])*8), m, n, &C_ref_abs[0], ldc ) ){
                abs_matrix(m, k, &A[0], lda);
                abs_matrix(m, k, &B[0], ldb);
                init_zero_matrix(m, n, &C_ref_abs[0], ldc);
                reference_gemm(&transa, &transb, &m, &n, &k, &alpha, &A[0], &lda, &B[0], &ldb, &beta, &C_ref_abs[0], &ldc);
                dump_matrix(m, n, &C_ref_abs[0], ldc, path_string + "/values.fp" + std::to_string(sizeof(C_ref_abs[0])*8));
            }

            // compute and dump componentwise relative error and componentwise f(n)
            init_zero_matrix(m, n, &rel_errors[0], ldc);
            init_zero_matrix(m, n, &f_of_n[0], ldc);
            calculate_relative_errors(m, n, &C[0], ldc, &C_ref[0], ldc, &rel_errors[0], ldc);
            calculate_f_of_n(m, n, &C[0], ldc, &C_ref[0], ldc, &C_ref_abs[0], ldc, &f_of_n[0], ldc);
            path_string = OUTPUT_DIR + construct_metadata_path(dgemm_impl_name, INPUT_TYPE, m, n, k, s);
            dump_matrix(m, n, &rel_errors[0], ldc, path_string + "/rel_errors.fp64");
            dump_matrix(m, n, &f_of_n[0], ldc, path_string + "/f_of_n.fp64");
            display_distribution_stats(m, n, &rel_errors[0], ldc, path_string + "/rel_errors.fp64");
#endif
        }
    }
    return 0;
}