#include "dgemm_wrapper.hpp"
#include "init_matrices.hpp"
#include "utils.hpp"
#include <vector>
#include <cstdlib>
#include "reference_gemm.hpp"

int main(){

    const char* env_dir = std::getenv("RUN_OUTPUT_DIR");
    if ( !env_dir ){
        fprintf(stderr,
            "ERROR: RUN_OUTPUT_DIR is not set.\n"
            "       Run this binary via 'make run', which sets the output\n"
            "       directory automatically and co-locates the config files.\n");
        return 1;
    }
    std::string output_dir = std::string(env_dir) + "/";
    std::string path_string;

    char transa = 'N';
    char transb = 'N';
    double alpha = 1.0;
    double beta = 0.0;
    int seed = 0;
    
    int max_m, max_k, max_n, lda, ldb, ldc;
    max_m = max_k = max_n = lda = ldb = ldc = 1 << LOG2_MAX_DIM;
    
    std::vector<double> A(max_m*max_k);
    std::vector<double> B(max_k*max_n);
    std::vector<double> C(max_m*max_n);
    std::vector<double> rel_errors(max_m*max_n);
    std::vector<double> C_ref(max_m*max_n);

    for ( int d = LOG2_MIN_DIM; d <= LOG2_MAX_DIM; ++d ){
        
        int m, n, k;
        m = n = k = 1 << d;

        for ( int s = LOG2_MIN_SCALING_FACTOR; s <= LOG2_MAX_SCALING_FACTOR; s += LOG2_INC_SCALING_FACTOR ){
        
            // initialize A and B
            init_matrices(INPUT_TYPE, s, seed, m, n, k, &A[0], lda, &B[0], ldb, output_dir);

            // calculate ground truth before permuting A and B using the DGEMM implementation under test
            init_zero_matrix(m, n, &C_ref[0], ldc);
            init_zero_matrix(m, n, &rel_errors[0], ldc);
            dgemm_wrapper(transa, transb, m, n, k, alpha, &A[0], lda, &B[0], ldb, beta, &C_ref[0], ldc);

#if defined(OUTER_PERMUTE)
            std::string full_input_type_string = INPUT_TYPE;
            full_input_type_string += "_w_outer_permute";

            // set up OUTER_PERMUTE iterations
            int perm_seed = 0;
            for ( int i = 0; i < OUTER_PERMUTE; ++i ){
                
                init_matrices(INPUT_TYPE, s, seed, m, n, k, &A[0], lda, &B[0], ldb);
                
                perm_seed++;
                permute_rows(perm_seed, m, k, &A[0], lda);
                permute_cols(perm_seed + OUTER_PERMUTE, k, n, &B[0], ldb);

#elif defined(INNER_PERMUTE)
            std::string full_input_type_string = INPUT_TYPE;
            full_input_type_string += "_w_inner_permute";

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
                // calculate componentwise error, keeping the maximum over all permutations tested
                calculate_relative_errors(m, n, &C[0], ldc, &C_ref[0], ldc, &rel_errors[0], ldc);
            }

            // write to disk the maximum componentwise relative error (with respect to A*B using the implementation under test) over all permutations tested
            path_string = output_dir + construct_metadata_path(DGEMM_IMPL_NAME, full_input_type_string, m, n, k, s);
            dump_matrix(m, n, &rel_errors[0], ldc, path_string + "/rel_errors.fp64");
        }
    }
    return 0;
}