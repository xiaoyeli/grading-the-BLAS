#ifndef INIT_MATRICES_H
#define INIT_MATRICES_H

#include <stdint.h>
#include <random>
#include <algorithm>
#include <vector>

#define Z_RATIO 0.33

template <typename T>
void init_zero_matrix( T *A, const uint32_t lda, const uint32_t m, const uint32_t n ){

    uint32_t i, j;

    for ( i = 0; i < m; ++i ){
        for ( j = 0; j < n; ++j ){
            A[i*lda + j] = 0.0;
        }
    }
    return;
}

template <typename T_to, typename T_from>
void copy_matrix( T_to *A, const uint32_t lda, T_from *B, const uint32_t ldb, const uint32_t m, uint32_t n ){
    
    uint32_t i, j;

    for ( i = 0; i < m; ++i ){
        for ( j = 0; j < n; ++j ){
            A[i*lda + j] = B[i*ldb + j];
        }
    }
    return;
}

template <typename T>
void init_test1_matrices( T *A, const uint32_t lda, T *B, const uint32_t ldb, const uint32_t m, const uint32_t k, const uint32_t n){

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<T> generator(0.0d, 1.0d);

    uint32_t i, j;

    // figure out which rows/columns to zero
    std::vector<bool> zero_rows_A(m);
    std::vector<bool> zero_cols_B(n);
    uint32_t n_zero;

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
    for ( i = n_zero; i < m; ++i ){
        zero_cols_B[i] = false;
    }

    std::shuffle(zero_rows_A.begin(), zero_rows_A.end(), rd);
    std::shuffle(zero_cols_B.begin(), zero_cols_B.end(), rd);

    // initialize matrices
    for ( i = 0; i < m; ++i ){
        for ( j = 0; j < k; ++j ){
            if ( zero_rows_A[i] ){
                A[i*lda + j] = 0.0;
            }
            else{
                A[i*lda + j] = generator(rd);
            }
        }
    }

    for ( i = 0; i < k; ++i ){
        for ( j = 0; j < n; ++j ){
            if ( zero_cols_B[j] ){
                B[i*ldb + j] = 0.0;
            }
            else{
                B[i*ldb + j] = generator(rd);
            }
        }
    }  

    return;
}

#endif