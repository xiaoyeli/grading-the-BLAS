#ifndef INIT_MATRICES_HPP
#define INIT_MATRICES_HPP

#include <stdint.h>
#include <random>
#include <algorithm>
#include <vector>

#define Z_RATIO 0.25

template <typename T>
void init_zero_matrix( T *A, const int32_t lda, const int32_t m, const int32_t n ){

    int32_t i, j;

    for ( i = 0; i < m; ++i ){
        for ( j = 0; j < n; ++j ){
            A[i + j*lda] = 0.0;
        }
    }
    return;
}

template <typename T>
void init_sanity_matrices( T *A, const int32_t lda, T *B, const int32_t ldb, const int32_t m, const int32_t k, const int32_t n ){

    int32_t i, j;

    // initialize matrices
    for ( i = 0; i < m; ++i ){
        for ( j = 0; j < k; ++j ){
            A[i + j*lda] = i + j*lda;
        }
    }

    for ( i = 0; i < k; ++i ){
        for ( j = 0; j < n; ++j ){
            B[i + j*ldb] = 1 + i + j*ldb;
        }
    }  

    return;
}

template <typename T>
void init_test1_matrices( T *A, const int32_t lda, T *B, const int32_t ldb, const int32_t m, const int32_t k, const int32_t n ){

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<T> generator(0.0d, 1.0d);

    int32_t i, j;

    // figure out which rows/columns to zero
    std::vector<bool> zero_rows_A(m);
    std::vector<bool> zero_cols_B(n);
    int32_t n_zero;

    n_zero = std::max(1, (int)(Z_RATIO*m));
    for ( i = 0; i < n_zero; ++i ){
        zero_rows_A[i] = true;
    }
    for ( i = n_zero; i < m; ++i ){
        zero_rows_A[i] = false;
    }

    n_zero = std::max(1, (int)(Z_RATIO*n));
    for ( i = 0; i < n_zero; ++i ){
        zero_cols_B[i] = true;
    }
    for ( i = n_zero; i < n; ++i ){
        zero_cols_B[i] = false;
    }

    std::shuffle(zero_rows_A.begin(), zero_rows_A.end(), rd);
    std::shuffle(zero_cols_B.begin(), zero_cols_B.end(), rd);

    // initialize matrices
    for ( i = 0; i < m; ++i ){
        for ( j = 0; j < k; ++j ){
            if ( zero_rows_A[i] ){
                A[i + j*lda] = 0.0;
            }
            else{
                A[i + j*lda] = generator(rd);
            }
        }
    }

    for ( i = 0; i < k; ++i ){
        for ( j = 0; j < n; ++j ){
            if ( zero_cols_B[j] ){
                B[i + j*ldb] = 0.0;
            }
            else{
                B[i + j*ldb] = generator(rd);
            }
        }
    }  

    return;
}

#endif