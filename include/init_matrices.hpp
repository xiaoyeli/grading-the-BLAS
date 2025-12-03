#ifndef INIT_MATRICES_HPP
#define INIT_MATRICES_HPP

#include <stdint.h>
#include <random>
#include <algorithm>
#include <vector>
#include <cassert>

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
void init_random_matrices( T *A, const int32_t lda, T *B, const int32_t ldb, const int32_t m, const int32_t k, const int32_t n, const T lower_bound, const T upper_bound ){
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<T> generator(lower_bound, upper_bound);
    int32_t i, j;

    // initialize matrices
    for ( i = 0; i < m; ++i ){
        for ( j = 0; j < k; ++j ){
            A[i + j*lda] = generator(rd);
        }
    }

    for ( i = 0; i < k; ++i ){
        for ( j = 0; j < n; ++j ){
            B[i + j*ldb] = generator(rd);
        }
    }  

    return;
}

template <typename T>
void init_test1_matrices( T *A, const int32_t lda, T *B, const int32_t ldb, const int32_t m, const int32_t k, const int32_t n ){

    init_random_matrices( A, lda, B, ldb, m, k, n, static_cast<T>(0.0), static_cast<T>(1.0) );
    int32_t i, j;

    // figure out which rows/columns to zero
    std::random_device rd;
    std::mt19937 gen(rd());
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

    // zero selected rows/columns
    for ( i = 0; i < m; ++i ){
        for ( j = 0; j < k; ++j ){
            if ( zero_rows_A[i] ){
                A[i + j*lda] = 0.0;
            }
        }
    }

    for ( i = 0; i < k; ++i ){
        for ( j = 0; j < n; ++j ){
            if ( zero_cols_B[j] ){
                B[i + j*ldb] = 0.0;
            }
        }
    }  

    return;
}

void init_test3_matrices_fp64( double *A, const int32_t lda, double *B, const int32_t ldb, const int32_t dim ){

    int32_t i, j;
    int32_t max_exponent = 1000;
    int32_t min_exponent = -1000;
    std::vector<double> x(dim);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> generator(1.0, 2.0);
    std::vector<double> diag_D(dim);
    std::vector<double> diag_D_inv(dim);

    for ( i = 0; i < dim; ++i ){
        x[i] = generator(rd);
    }

    for ( i = 0; i < dim; ++i ){
        diag_D[i] = std::ldexp(1.0, min_exponent + (i * (max_exponent - min_exponent)) / (dim - 1));
        diag_D_inv[i] = std::ldexp(1.0, -min_exponent - (i * (max_exponent - min_exponent)) / (dim - 1));
    }

    for ( i = 0; i < dim; ++i ){
        for ( j = 0; j < dim; ++j ){
            A[i + j*lda] = x[j] * diag_D[j];
        }
    }

    for ( i = 0; i < dim; ++i ){
        for ( j = 0; j < dim; ++j ){
            B[i + j*ldb] = x[i] * diag_D_inv[i];
        }
    }
}

template <typename T>
void init_matrices( std::string type, T *A, const int32_t lda, T *B, const int32_t ldb, const int32_t m, const int32_t k, const int32_t n ){
    if ( type == "U_0_1" ){
        init_random_matrices( A, lda, B, ldb, m, k, n, static_cast<T>(0.0), static_cast<T>(1.0) );
    }
    if ( type == "U_-1_1" ){
        init_random_matrices( A, lda, B, ldb, m, k, n, static_cast<T>(-1.0), static_cast<T>(1.0) );
    }
    else if ( type == "test1" ){
        init_test1_matrices( A, lda, B, ldb, m, k, n );
    }
    else if ( type == "test3" ){
        static_assert( std::is_same<T, double>::value, "test3 matrices are only implemented for double precision." );
        assert( m == n && n == k && "test3 matrices are only implemented for square matrices.");
        init_test3_matrices_fp64( A, lda, B, ldb, m );
    }
    else{
        // default to U(0,1)
        init_random_matrices( A, lda, B, ldb, m, k, n, static_cast<T>(0.0), static_cast<T>(1.0) );
    }

    return;
}

#endif