#include "dgemm_wrapper.hpp"
#include "init_matrices.hpp"
#include "utils.hpp"
#include <vector>
#include <cstdlib>
#include "reference_gemm.hpp"

#define RAW_DATA_DIR "results/raw_data/"

int main(){

    std::string raw_data_dir = RAW_DATA_DIR;
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
    std::vector<double> C_ref(max_m*max_n);

    for ( int d = LOG2_MIN_DIM; d <= LOG2_MAX_DIM; ++d ){
        
        int m, n, k;
        m = n = k = 1 << d;

        for ( int s = LOG2_MIN_SCALING_FACTOR; s <= LOG2_MAX_SCALING_FACTOR; s += LOG2_INC_SCALING_FACTOR ){
        
            // initialize A and B
            init_matrices(INPUT_TYPE, s, seed, m, n, k, &A[0], lda, &B[0], ldb, raw_data_dir);

            // calculate ground truth before permuting A and B using the DGEMM implementation under test
            init_zero_matrix(m, n, &C_ref[0], ldc);
            double max_rel_error = 0.0;
            dgemm_wrapper(transa, transb, m, n, k, alpha, &A[0], lda, &B[0], ldb, beta, &C_ref[0], ldc);

#if defined(OUTER_PERMUTE)
            std::string full_input_type_string = INPUT_TYPE;
            full_input_type_string += "_w_outer_permute";

            // set up OUTER_PERMUTE iterations
            int perm_seed = 0;
            for ( int i = 0; i < OUTER_PERMUTE; ++i ){
                
                init_matrices(INPUT_TYPE, s, seed, m, n, k, &A[0], lda, &B[0], ldb, raw_data_dir);
                
                perm_seed++;
                permute_rows(perm_seed, m, k, &A[0], lda);
                permute_cols(perm_seed + OUTER_PERMUTE, k, n, &B[0], ldb);

#elif defined(INNER_PERMUTE)
            std::string full_input_type_string = INPUT_TYPE;
            full_input_type_string += "_w_inner_permute";

            // set up INNER_PERMUTE iterations
            int perm_seed = 0;
            for ( int i = 0; i < INNER_PERMUTE; ++i ){
                
                init_matrices(INPUT_TYPE, s, seed, m, n, k, &A[0], lda, &B[0], ldb, raw_data_dir);
                
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
                // accumulate max relative error over all permutations tested
                double perm_max = calculate_max_relative_error(m, n, &C[0], ldc, &C_ref[0], ldc);
                if ( std::isnan(perm_max) || perm_max > max_rel_error ) max_rel_error = perm_max;
            }

            // write scalar max relative error to disk
            path_string = raw_data_dir + construct_metadata_path(DGEMM_IMPL_NAME, full_input_type_string, m, n, k, s);
            dump_scalar(max_rel_error, path_string + "/max_rel_error.fp64");
        }
    }
    return 0;
}