#include "dgemm_wrapper.hpp"
#include "init_matrices.hpp"
#include "utils.hpp"
#include <vector>
#include <cstdlib>

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

    for ( int d = LOG2_MIN_DIM; d <= LOG2_MAX_DIM; ++d ){
        
        int m, n, k;
        m = n = k = 1 << d;

        for ( int s = LOG2_MIN_SCALING_FACTOR; s <= LOG2_MAX_SCALING_FACTOR; s += LOG2_INC_SCALING_FACTOR ){
        
            // initialize A and B
            init_matrices(INPUT_TYPE, s, seed, m, n, k, &A[0], lda, &B[0], ldb, output_dir);
            
            // call the DGEMM implementation under test to compute C
            init_zero_matrix(m, n, &C[0], ldc);
            dgemm_wrapper(transa, transb, m, n, k, alpha, &A[0], lda, &B[0], ldb, beta, &C[0], ldc);

            // write C to disk
            path_string = output_dir + construct_metadata_path(DGEMM_IMPL_NAME, INPUT_TYPE, m, n, k, s);
            dump_matrix(m, n, &C[0], ldc, path_string + "/values.fp64");
        }
    }
    
    return 0;
}