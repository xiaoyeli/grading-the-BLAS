#include "dgemm_wrapper.hpp"
#include "init_matrices.hpp"
#include "utils.hpp"
#include <vector>
#include <cstdlib>
#include <cstdint>
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
#if REFERENCE_GEMM == 1
    std::vector<double> C_ref(max_m*max_n);
    std::vector<double> absA_absB_product(max_m*max_n);
#elif REFERENCE_GEMM == 2
    std::vector<long double> C_ref(max_m*max_n);
    std::vector<long double> absA_absB_product(max_m*max_n);
#endif

    for ( int d = LOG2_MIN_DIM; d <= LOG2_MAX_DIM; ++d ){
        
        int m, n, k;
        m = n = k = 1 << d;

        for ( int s = LOG2_MIN_SCALING_FACTOR; s <= LOG2_MAX_SCALING_FACTOR; s += LOG2_INC_SCALING_FACTOR ){
        
            // initialize A and B
            init_matrices(INPUT_TYPE, s, seed, m, n, k, &A[0], lda, &B[0], ldb, raw_data_dir);

            // Build inclusion mask: entry (i,j) is 0 if any summand A[i,t]*B[t,j] < 0,
            // indicating that we are not to compute f_A(n) at C[i,j]
            std::vector<uint8_t> include_mask(m * n, 1);
            for ( int i = 0; i < m; ++i ){
                for ( int j = 0; j < n; ++j ){
                    for ( int t = 0; t < k; ++t ){
                        if ( A[i + t*lda] * B[t + j*ldb] < 0.0 ){
                            include_mask[i + j*m] = 0;
                            break;
                        }
                    }
                }
            }

            // call the DGEMM implementation under test
            init_zero_matrix(m, n, &C[0], ldc);
            dgemm_wrapper(transa, transb, m, n, k, alpha, &A[0], lda, &B[0], ldb, beta, &C[0], ldc);

            // write C to disk
            path_string = raw_data_dir + construct_metadata_path(DGEMM_IMPL_NAME, INPUT_TYPE, m, n, k, s);
            dump_matrix(m, n, &C[0], ldc, path_string + "/values.fp64");

            // either load a cached C_ref or compute from scratch and cache for later
            path_string = raw_data_dir + construct_metadata_path("C_ref", INPUT_TYPE, m, n, k, s);
            if ( ! load_matrix(path_string + "/values.fp" + std::to_string(sizeof(C_ref[0])*8), m, n, &C_ref[0], ldc ) ){
                init_zero_matrix(m, n, &C_ref[0], ldc);
                reference_gemm(&transa, &transb, &m, &n, &k, &alpha, &A[0], &lda, &B[0], &ldb, &beta, &C_ref[0], &ldc);
                dump_matrix(m, n, &C_ref[0], ldc, path_string + "/values.fp" + std::to_string(sizeof(C_ref[0])*8));
            }

            // either load a cached absA_absB_product or compute from scratch and cache for later
            path_string = raw_data_dir + construct_metadata_path("absA_absB_product", INPUT_TYPE, m, n, k, s);
            if ( ! load_matrix(path_string + "/values.fp" + std::to_string(sizeof(absA_absB_product[0])*8), m, n, &absA_absB_product[0], ldc ) ){
                abs_matrix(m, k, &A[0], lda);
                abs_matrix(k, n, &B[0], ldb);
                init_zero_matrix(m, n, &absA_absB_product[0], ldc);
                reference_gemm(&transa, &transb, &m, &n, &k, &alpha, &A[0], &lda, &B[0], &ldb, &beta, &absA_absB_product[0], &ldc);
                dump_matrix(m, n, &absA_absB_product[0], ldc, path_string + "/values.fp" + std::to_string(sizeof(absA_absB_product[0])*8));
            }

            // compute max f_A(n) and write scalar to disk
            double max_f = calculate_max_f_of_n(m, n, &C[0], ldc, &C_ref[0], ldc, &absA_absB_product[0], ldc, include_mask.data(), m);
            path_string = raw_data_dir + construct_metadata_path(DGEMM_IMPL_NAME, INPUT_TYPE, m, n, k, s);
            dump_scalar(max_f, path_string + "/max_f_of_n.fp64");        }
    }
    return 0;
}