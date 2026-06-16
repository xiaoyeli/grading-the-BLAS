#ifndef INIT_MATRICES_HPP
#define INIT_MATRICES_HPP

#include <random>
#include <algorithm>
#include <vector>
#include <cassert>
#include <cmath>
#include <string>
#include "utils.hpp"

#define TEST1A_RATIO 0.5        // ratio of zero rows / cols in A and B
#define TEST1B_RATIO_1 0.5      // ratio of targeted rows / cols in A and B
#define TEST1B_RATIO_2 0.5      // ratio of zero entries in targeted rows of A; the ratio of zero entries in targeted cols of B is 1 - TEST1B_RATIO_2
#define TEST1C_RATIO_1 0.5      // ratio of targeted rows / cols in A and B
#define TEST1C_RATIO_2 0.5      // ratio of zero entries in targeted rows of A; the ratio of zero entries in targeted cols of B is 1 - TEST1C_RATIO_2
#define TEST1C_RATIO_3 0.5      // ratio of targeted rows / cols in which to add a nonzero
#define TEST2C_RATIO_1 0.25     // ratio of targeted rows / cols in A and B
#define TEST2C_RATIO_2 0.25     // ratio of entries to negate in targeted rows / cols of A and B

#ifndef SCALING_FACTOR
#define SCALING_FACTOR 0.0 / 0.0
#endif

template <typename T>
void init_zero_matrix( const int m, const int n, T *A, const int lda ){

    for ( int i = 0; i < m; ++i ){
        for ( int j = 0; j < n; ++j ){
            A[i + j*lda] = 0.0;
        }
    }
    return;
}

template <typename T>
void init_sanity_matrices( const int m, const int n, const int k, T *A, const int lda, T *B, const int ldb ){

    // initialize matrices
    for ( int i = 0; i < m; ++i ){
        for ( int j = 0; j < k; ++j ){
            A[i + j*lda] = i + j*lda;
        }
    }

    for ( int i = 0; i < k; ++i ){
        for ( int j = 0; j < n; ++j ){
            B[i + j*ldb] = 1 + i + j*ldb;
        }
    }  

    return;
}

template <typename T>
void init_uniform_random_matrices(  const T lower_bound, const T upper_bound, int seed, const int m, const int n, const int k, T *A, const int lda, T *B, const int ldb ){
    
    std::mt19937 gen(seed);
    std::uniform_real_distribution<T> generator(lower_bound, upper_bound);

    // initialize matrices
    for ( int i = 0; i < m; ++i ){
        for ( int j = 0; j < k; ++j ){
            A[i + j*lda] = generator(gen);
        }
    }

    for ( int i = 0; i < k; ++i ){
        for ( int j = 0; j < n; ++j ){
            B[i + j*ldb] = generator(gen);
        }
    }  

    return;
}

template <typename T>
void init_guassian_random_matrices(  const T mean, const T std, int seed, const int m, const int n, const int k, T *A, const int lda, T *B, const int ldb ){
    
    std::mt19937 gen(seed);
    std::normal_distribution<T> generator(mean, std);

    // initialize matrices
    for ( int i = 0; i < m; ++i ){
        for ( int j = 0; j < k; ++j ){
            A[i + j*lda] = generator(gen);
        }
    }

    for ( int i = 0; i < k; ++i ){
        for ( int j = 0; j < n; ++j ){
            B[i + j*ldb] = generator(gen);
        }
    }  

    return;
}

template <typename T>
void init_test1a_matrices( int seed, const int dim, T *A, const int lda, T *B, const int ldb,
                            double *C_mask = nullptr, const int ldcm = 0, const std::string& mask_path = "" ){

    init_guassian_random_matrices( static_cast<T>(0.0), static_cast<T>(1.0), seed, dim, dim, dim, A, lda, B, ldb );

    // figure out which rows/columns to zero
    std::mt19937 gen(seed);
    std::vector<bool> zero_rows_mask_A(dim);
    std::vector<bool> zero_cols_mask_B(dim);
    int n_zero = std::max(1, (int)(TEST1A_RATIO*dim));

    for ( int i = 0; i < n_zero; ++i ){
        zero_rows_mask_A[i] = true;
        zero_cols_mask_B[i] = true;
    }
    for ( int i = n_zero; i < dim; ++i ){
        zero_rows_mask_A[i] = false;
        zero_cols_mask_B[i] = false;
    }

    std::shuffle(zero_rows_mask_A.begin(), zero_rows_mask_A.end(), gen);
    std::shuffle(zero_cols_mask_B.begin(), zero_cols_mask_B.end(), gen);

    // zero selected rows/columns; C_mask[i,j] = 1 wherever a full zero row or col forces C_ij = 0
    for ( int i = 0; i < dim; ++i ){
        for ( int j = 0; j < dim; ++j ){
            if ( zero_rows_mask_A[i] ){
                A[i + j*lda] = 0.0;
            }
            if ( zero_cols_mask_B[j] ){
                B[i + j*ldb] = 0.0;
            }
            if ( C_mask ){
                C_mask[i + j*ldcm] = (zero_rows_mask_A[i] || zero_cols_mask_B[j]) ? 1.0 : 0.0;
            }
        }
    }

    if ( C_mask && !mask_path.empty() )
        dump_matrix(dim, dim, C_mask, ldcm, mask_path);

    return;
}

template <typename T>
void init_test1b_matrices( int seed, const int dim, T *A, const int lda, T *B, const int ldb,
                            double *C_mask = nullptr, const int ldcm = 0, const std::string& mask_path = "" ){

    init_guassian_random_matrices( static_cast<T>(0.0), static_cast<T>(1.0), seed, dim, dim, dim, A, lda, B, ldb );

    std::vector<bool> target_rows_mask_A(dim);
    std::vector<bool> target_cols_mask_B(dim);
    std::vector<int> target_rows_idx_A;
    std::vector<int> target_cols_idx_B;
    std::vector<bool> sparsity_mask_A(dim);
    int n_targeted_rows_and_cols = std::max(1, (int)(TEST1B_RATIO_1*dim));
    int n_zeros_in_A_rows = std::max(1, (int)(TEST1B_RATIO_2*dim));

    for ( int i = 0; i < dim; ++i ){
        target_rows_mask_A[i] = false;
        target_cols_mask_B[i] = false;
        sparsity_mask_A[i] = false;
    }

    for ( int i = 0; i < n_targeted_rows_and_cols; ++i ){
        target_rows_mask_A[i] = true;
        target_cols_mask_B[i] = true;
    }
    for ( int i = 0; i < n_zeros_in_A_rows; ++i ){
        sparsity_mask_A[i] = true;
    }

    std::mt19937 gen(seed);
    std::shuffle(target_rows_mask_A.begin(), target_rows_mask_A.end(), gen);
    std::shuffle(target_cols_mask_B.begin(), target_cols_mask_B.end(), gen);

    for ( int i = 0; i < dim; ++i ){
        if ( target_rows_mask_A[i] ){
            target_rows_idx_A.push_back(i);
        }
        if ( target_cols_mask_B[i] ){
            target_cols_idx_B.push_back(i);
        }
    }

    std::shuffle(target_rows_idx_A.begin(), target_rows_idx_A.end(), gen);
    std::shuffle(target_cols_idx_B.begin(), target_cols_idx_B.end(), gen);

    if ( C_mask )
        for ( int i = 0; i < dim; ++i )
            for ( int j = 0; j < dim; ++j )
                C_mask[i + j*ldcm] = 0.0;

    for ( int s = 0; s < target_rows_idx_A.size(); ++s ){
        // every t contributes either A[row,t]=0 or B[t,col]=0 → C[row,col] = 0 exactly
        if ( C_mask )
            C_mask[target_rows_idx_A[s] + target_cols_idx_B[s]*ldcm] = 1.0;
        std::shuffle(sparsity_mask_A.begin(), sparsity_mask_A.end(), gen);
        for ( int t = 0; t < dim; ++t ){
            if ( sparsity_mask_A[t] ){
                A[target_rows_idx_A[s] + t*lda] = 0.0;
            }
            else{
                B[t + target_cols_idx_B[s]*ldb] = 0.0;
            }
        }
    }

    if ( C_mask && !mask_path.empty() )
        dump_matrix(dim, dim, C_mask, ldcm, mask_path);

    return;
}

template <typename T>
void init_test1c_matrices( int seed, const int dim, T *A, const int lda, T *B, const int ldb, int log2_max_scl_factor,
                            double *C_mask = nullptr, const int ldcm = 0, const std::string& mask_path = "" ){

    init_guassian_random_matrices( static_cast<T>(0.0), static_cast<T>(1.0), seed, dim, dim, dim, A, lda, B, ldb );

    std::vector<bool> target_rows_mask_A(dim);
    std::vector<bool> target_cols_mask_B(dim);
    std::vector<int> target_rows_idx_A;
    std::vector<int> target_cols_idx_B;
    std::vector<bool> sparsity_mask_A(dim);
    std::vector<bool> nonzero_mask(dim);
    int n_targeted_rows_and_cols = std::max(1, (int)(TEST1C_RATIO_1*dim));
    std::vector<bool> target_rows_cols_with_nonzero_mask(n_targeted_rows_and_cols);
    int n_zeros_in_A_rows = std::max(1, (int)(TEST1C_RATIO_2*dim));
    int n_targeted_rows_and_cols_w_a_nonzero = std::max(1, (int)(TEST1C_RATIO_3*n_targeted_rows_and_cols));

    for ( int i = 0; i < dim; ++i ){
        target_rows_mask_A[i] = false;
        target_cols_mask_B[i] = false;
        sparsity_mask_A[i] = false;
        nonzero_mask[i] = false;
    }
    for ( int i = 0; i < n_targeted_rows_and_cols; ++i ){
        target_rows_mask_A[i] = true;
        target_cols_mask_B[i] = true;
        target_rows_cols_with_nonzero_mask[i] = false;
    }
    for ( int i = 0; i < n_zeros_in_A_rows; ++i ){
        sparsity_mask_A[i] = true;
    }
    for ( int i = 0; i < n_targeted_rows_and_cols_w_a_nonzero; ++i ){
        target_rows_cols_with_nonzero_mask[i] = true;
    }
    nonzero_mask[0] = true;

    std::mt19937 gen(seed);
    std::shuffle(target_rows_mask_A.begin(), target_rows_mask_A.end(), gen);
    std::shuffle(target_cols_mask_B.begin(), target_cols_mask_B.end(), gen);
    std::shuffle(target_rows_cols_with_nonzero_mask.begin(), target_rows_cols_with_nonzero_mask.end(), gen);

    for ( int i = 0; i < dim; ++i ){
        if ( target_rows_mask_A[i] ){
            target_rows_idx_A.push_back(i);
        }
        if ( target_cols_mask_B[i] ){
            target_cols_idx_B.push_back(i);
        }
    }

    std::shuffle(target_rows_idx_A.begin(), target_rows_idx_A.end(), gen);
    std::shuffle(target_cols_idx_B.begin(), target_cols_idx_B.end(), gen);

    if ( C_mask )
        for ( int i = 0; i < dim; ++i )
            for ( int j = 0; j < dim; ++j )
                C_mask[i + j*ldcm] = 0.0;

    for ( int s = 0; s < target_rows_idx_A.size(); ++s ){
        // 2 = one tiny product survives (magnitude ~ 2^{-2*log2_max_scl_factor})
        // 1 = every A[row,t]*B[t,col] term is zero → C[row,col] = 0 exactly
        if ( C_mask )
            C_mask[target_rows_idx_A[s] + target_cols_idx_B[s]*ldcm] =
                target_rows_cols_with_nonzero_mask[s] ? 2.0 : 1.0;
        std::shuffle(sparsity_mask_A.begin(), sparsity_mask_A.end(), gen);
        std::shuffle(nonzero_mask.begin(), nonzero_mask.end(), gen);
        for ( int t = 0; t < dim; ++t ){
            if ( target_rows_cols_with_nonzero_mask[s] && nonzero_mask[t] ){
                A[target_rows_idx_A[s] + t*lda] = std::ldexp(A[target_rows_idx_A[s] + t*lda], -1.0 * log2_max_scl_factor);
                B[t + target_cols_idx_B[s]*ldb] = std::ldexp(B[t + target_cols_idx_B[s]*ldb], -1.0 * log2_max_scl_factor);
            }
            else if ( sparsity_mask_A[t] ){
                A[target_rows_idx_A[s] + t*lda] = 0.0;
                B[t + target_cols_idx_B[s]*ldb] = std::ldexp(B[t + target_cols_idx_B[s]*ldb], log2_max_scl_factor);
            }
            else{
                B[t + target_cols_idx_B[s]*ldb] = 0.0;
                A[target_rows_idx_A[s] + t*lda] = std::ldexp(A[target_rows_idx_A[s] + t*lda], log2_max_scl_factor);
            }
        }
    }

    if ( C_mask && !mask_path.empty() )
        dump_matrix(dim, dim, C_mask, ldcm, mask_path);

    return;
}

void init_test2a_matrices( int seed, const int dim, double *A, const int lda, double *B, const int ldb, int log2_max_scl_factor ){

    std::vector<double> x(dim);

    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> generator(1.0, 2.0);
    std::vector<double> diag_D(dim);
    std::vector<double> diag_D_inv(dim);
    std::vector<int> exponents(dim);

    for ( int i = 0; i < dim; ++i ){
        x[i] = generator(gen);
        exponents[i] = ((-1 * dim / 2) + 1 + i) * ((double) log2_max_scl_factor / ( dim / 2 ) );
    }

    std::shuffle(exponents.begin(), exponents.end(), gen);

    for ( int i = 0; i < dim; ++i ){
        diag_D[i] = std::ldexp(1.0, exponents[i]);
        diag_D_inv[i] = std::ldexp(1.0, -1 * exponents[i]);
    }
    
    for ( int i = 0; i < dim; ++i ){
        for ( int j = 0; j < dim; ++j ){
            A[i + j*lda] = x[j] * diag_D[j];
        }
    }

    for ( int i = 0; i < dim; ++i ){
        for ( int j = 0; j < dim; ++j ){
            B[i + j*ldb] = x[i] * diag_D_inv[i];
        }
    }
}

void init_test2b_matrices( int seed, const int dim, double *A, const int lda, double *B, const int ldb, int log2_max_scl_factor ){

    std::vector<double> unpermuted_A(dim * dim);
    std::vector<double> unpermuted_B(dim * dim);

    init_test2a_matrices(seed, dim, &unpermuted_A[0], dim, &unpermuted_B[0], dim, log2_max_scl_factor);

    // circular shifting
    for ( int i = 1; i < dim; ++i ){
        for ( int j = 0; j < dim; ++j ){
            unpermuted_A[i + j*dim] = unpermuted_A[0 + ((dim + j - i)%dim)*dim];
        }
    }

    for ( int j = 1; j < dim; ++j ){
        for ( int i = 0; i < dim; ++i ){
            unpermuted_B[i + j*dim] = unpermuted_B[(dim + i - j)%dim + 0*dim];
        }
    }

    // permuting
    std::vector<int> p1(dim);
    std::vector<int> p2(dim);
    std::vector<int> p3(dim);
    
    for ( int i = 0; i < dim; ++i ){
        p1[i] = i;
        p2[i] = i;
        p3[i] = i;
    }

    std::mt19937 gen(seed);
    std::shuffle(p1.begin(), p1.end(), gen);
    std::shuffle(p2.begin(), p2.end(), gen);
    std::shuffle(p3.begin(), p3.end(), gen);

    for ( int i = 0; i < dim; ++i ){
        for ( int j = 0; j < dim; ++j ){
            A[i + j*lda] = unpermuted_A[p1[i] + p3[j]*dim];
        }
    }

    for ( int j = 0; j < dim; ++j ){
        for ( int i = 0; i < dim; ++i ){
            B[i + j*ldb] = unpermuted_B[p3[i] + p2[j]*dim];
        }
    }
}

void init_test2c_matrices( int seed, const int dim, double *A, const int lda, double *B, const int ldb, int log2_max_scl_factor ){

    init_test2b_matrices(seed, dim, A, lda, B, ldb, log2_max_scl_factor);

    std::vector<bool> target_vector_mask(dim);
    std::vector<bool> negation_mask(dim);
    int n_targeted_vectors = std::max(1, (int)(TEST2C_RATIO_1*dim));
    int n_negations = std::max(1, (int)(TEST2C_RATIO_2*dim));

    for ( int i = 0; i < dim; ++i ){
        target_vector_mask[i] = false;
        negation_mask[i] = false;
    }

    for ( int i = 0; i < n_targeted_vectors; ++i ){
        target_vector_mask[i] = true;
    }

    for ( int i = 0; i < n_negations; ++i ){
        negation_mask[i] = true;
    }

    std::mt19937 gen(seed);

    
    // pick a small random subset of rows of A and negate some of their entries
    std::shuffle(target_vector_mask.begin(), target_vector_mask.end(), gen);
    for ( int i = 0; i < dim; ++i ){
        if ( target_vector_mask[i] ){
            std::shuffle(negation_mask.begin(), negation_mask.end(), gen);
            for ( int j = 0; j < dim; ++j ){
                if ( negation_mask[j] ){
                    A[i + j*lda] *= -1.0;
                }
            }
        }
    }

    // pick a small random subset of columns of B and negate some of their entries
    std::shuffle(target_vector_mask.begin(), target_vector_mask.end(), gen);
    for ( int j = 0; j < dim; ++j ){
        if ( target_vector_mask[j] ){
            std::shuffle(negation_mask.begin(), negation_mask.end(), gen);
            for ( int i = 0; i < dim; ++i ){
                if ( negation_mask[i] ){
                    B[i + j*ldb] *= -1.0;
                }
            }
        }
    }

    // Let S be a diagonal matrix with random +/- 1 diagonal entries and let A = A*S and B = S*B
    std::vector<double> S(dim);
    std::uniform_int_distribution<int> sign_dist(0, 1);
    for ( int i = 0; i < dim; ++i ){
        S[i] = (sign_dist(gen) == 0) ? -1.0 : 1.0;
    }

    for ( int i = 0; i < dim; ++i ){
        for ( int j = 0; j < dim; ++j ){
            A[i + j*lda] *= S[j];
            B[i + j*ldb] *= S[i];
        }
    }
}

void init_test3_matrices( int seed, const int dim, double *A, const int lda, double *B, const int ldb, int log2_max_scl_factor ){
        
    init_uniform_random_matrices( static_cast<double>(1.0), static_cast<double>(2.0), seed, dim, dim, dim, A, lda, B, ldb );
    
    std::mt19937 gen(seed+1);
    std::uniform_real_distribution<double> generator(1.0, 2.0);
    std::vector<double> diag_D1(dim);
    std::vector<double> diag_D2(dim);

    for ( int i = 0; i < dim; ++i ){
        diag_D1[i] = std::ldexp(generator(gen), ((-1 * dim / 2) + 1 + i) * ((double) log2_max_scl_factor / ( dim / 2 ) ));
        diag_D2[i] = std::ldexp(generator(gen), ((-1 * dim / 2) + 1 + i) * ((double) log2_max_scl_factor / ( dim / 2 ) ));
    }

    for ( int i = 0; i < dim/2; i = i + 2 ){
        std::swap(diag_D1[i], diag_D1[dim - 1 - i]);
        std::swap(diag_D2[i], diag_D2[dim - 1 - i]);
    }

    // std::shuffle(diag_D1.begin(), diag_D1.end(), gen);
    // std::shuffle(diag_D2.begin(), diag_D2.end(), gen);

    for ( int i = 0; i < dim; ++i ){
        for ( int j = 0; j < dim; ++j ){
            A[i + j*lda] = A[i + j*lda] * diag_D1[i];
            B[i + j*ldb] = B[i + j*ldb] * diag_D2[j];
        }
    }
}

void init_rubric_baseline( int seed, const int dim, double *A, const int lda, double *B, const int ldb, int log2_max_scl_factor ){

    std::vector<double> x(dim);

    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> generator(1.0, 2.0);
    std::vector<double> diag_D(dim);
    std::vector<double> diag_D_inv(dim);
    std::vector<int> exponents(dim);

    for ( int i = 0; i < dim; ++i ){
        x[i] = generator(gen);
        exponents[i] = ((-1 * dim / 2) + 1 + i) * ((double) log2_max_scl_factor / ( dim / 2 ) );
    }

    std::shuffle(exponents.begin(), exponents.end(), gen);

    for ( int i = 0; i < dim; ++i ){
        diag_D[i] = std::ldexp(1.0, exponents[i]);
        diag_D_inv[i] = std::ldexp(1.0, -1 * exponents[i]);
    }
    
    for ( int i = 0; i < dim; ++i ){
        for ( int j = 0; j < dim; ++j ){
            A[i + j*lda] = x[j] * diag_D[j];
        }
    }

    for ( int i = 0; i < dim; ++i ){
        for ( int j = 0; j < dim; ++j ){
            B[i + j*ldb] = x[i] * diag_D_inv[i];
        }
    }
}

template <typename T>
void init_matrices( std::string type, const int log2_max_scl_factor, const int seed, const int m, const int n, const int k, T *A, const int lda, T *B, const int ldb, std::string output_dir, 
                    double *C_mask = nullptr, const int ldcm = 0, const std::string& mask_path = "" ){
    if ( type == "none" ){
        assert(false && "No input type specified." );
    }
    if ( type == "sanity" ){
        init_sanity_matrices(  m, n, k, A, lda, B, ldb );
    }
    else if ( type == "U(0,1)" ){
        init_uniform_random_matrices( static_cast<T>(0.0), static_cast<T>(1.0), seed, m, n, k, A, lda, B, ldb );
    }
    else if ( type == "U(-1,1)" ){
        init_uniform_random_matrices( static_cast<T>(-1.0), static_cast<T>(1.0), seed, m, n, k, A, lda, B, ldb );
    }
    else if ( type == "U(1,2)" ){
        init_uniform_random_matrices( static_cast<T>(1.0), static_cast<T>(2.0), seed, m, n, k, A, lda, B, ldb );
    }
    else{
        assert( m == n && n == k && "test matrices must be square.");
        assert( seed == 0 && "test matrices must be deterministically generated" );

        if ( type == "test1a" ){
            init_test1a_matrices( seed, m, A, lda, B, ldb, C_mask, ldcm, mask_path );
        }
        else if ( type == "test1b" ){
            init_test1b_matrices( seed, m, A, lda, B, ldb, C_mask, ldcm, mask_path );
        }
        else if ( type == "test1c" ){
            init_test1c_matrices( seed, m, A, lda, B, ldb, log2_max_scl_factor, C_mask, ldcm, mask_path );
        }
        else if ( type == "test2a" ){
            init_test2a_matrices( seed, m, A, lda, B, ldb, log2_max_scl_factor);
        }
        else if ( type == "test2b" ){
            init_test2b_matrices( seed, m, A, lda, B, ldb, log2_max_scl_factor);
        }
        else if ( type == "test2c" ){
            init_test2c_matrices( seed, m, A, lda, B, ldb, log2_max_scl_factor);
        }
        else if ( type == "test3" ){
            init_test3_matrices( seed, m, A, lda, B, ldb, log2_max_scl_factor );
        }
        else if ( type == "rubric_baseline" ){
            int dim = m;

            init_guassian_random_matrices( static_cast<T>(0.0), static_cast<T>(1.0), seed, dim, dim, dim, A, lda, B, ldb );
            
            std::mt19937 gen(seed);
            std::vector<double> diag_D(dim);
            std::vector<double> diag_D_inv(dim);
            std::vector<int> exponents(dim);

            for ( int i = 0; i < dim; ++i ){
                exponents[i] = ((-1 * dim / 2) + 1 + i) * ((double) log2_max_scl_factor / ( dim / 2 ) );
            }

            std::shuffle(exponents.begin(), exponents.end(), gen);

            for ( int i = 0; i < dim; ++i ){
                diag_D[i] = std::ldexp(1.0, exponents[i]);
                diag_D_inv[i] = std::ldexp(1.0, -1 * exponents[i]);
            }
            
            for ( int i = 0; i < dim; ++i ){
                for ( int j = 0; j < dim; ++j ){
                    A[i + j*lda] = A[i + j*lda] * diag_D[j];
                }
            }

            for ( int i = 0; i < dim; ++i ){
                for ( int j = 0; j < dim; ++j ){
                    B[i + j*ldb] = B[i + j*ldb] * diag_D_inv[i];
                }
            }
        }
        else{
            assert(false && "Unknown matrix type specified in init_matrices." );
        }
    }

#if defined(VERBOSE)
    printf("A matrix:\n");
    print_matrix_to_stdout_sci(m, k, &A[0], lda);
    printf("B matrix:\n");
    print_matrix_to_stdout_sci(k, n, &B[0], ldb);
#endif
#if defined(DUMP_A_B)
    // dump A and B; note that dump_matrix overwrites every time
    path_string = output_dir + construct_metadata_path("A", type, m, n, k, log2_max_scl_factor);
    dump_matrix(m, k, &A[0], lda, path_string + "/values.fp64");
    path_string = output_dir + construct_metadata_path("B", type, m, n, k, log2_max_scl_factor);
    dump_matrix(k, n, &B[0], ldb, path_string + "/values.fp64");
#endif

    return;
}

#endif