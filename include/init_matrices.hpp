#ifndef INIT_MATRICES_HPP
#define INIT_MATRICES_HPP

#include <stdint.h>
#include <random>
#include <algorithm>
#include <vector>
#include <cassert>

#define Z_RATIO 0.25

template <typename T>
void init_zero_matrix( const int32_t m, const int32_t n, T *A, const int32_t lda ){

    int32_t i, j;

    for ( i = 0; i < m; ++i ){
        for ( j = 0; j < n; ++j ){
            A[i + j*lda] = 0.0;
        }
    }
    return;
}

template <typename T>
void init_sanity_matrices( const int32_t m, const int32_t n, const int32_t k, T *A, const int32_t lda, T *B, const int32_t ldb ){

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
void init_random_matrices(  const T lower_bound, const T upper_bound, int32_t seed, const int32_t m, const int32_t n, const int32_t k, T *A, const int32_t lda, T *B, const int32_t ldb ){
    
    std::mt19937 gen(seed);
    std::uniform_real_distribution<T> generator(lower_bound, upper_bound);
    int32_t i, j;

    // initialize matrices
    for ( i = 0; i < m; ++i ){
        for ( j = 0; j < k; ++j ){
            A[i + j*lda] = generator(gen);
        }
    }

    for ( i = 0; i < k; ++i ){
        for ( j = 0; j < n; ++j ){
            B[i + j*ldb] = generator(gen);
        }
    }  

    return;
}

template <typename T>
void init_test1_matrices( int32_t seed, const int32_t m, const int32_t n, const int32_t k, T *A, const int32_t lda, T *B, const int32_t ldb ){

    init_random_matrices( static_cast<T>(0.0), static_cast<T>(1.0), seed, m, n, k, A, lda, B, ldb );
    int32_t i, j;

    // figure out which rows/columns to zero
    std::mt19937 gen(seed);
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

    std::shuffle(zero_rows_A.begin(), zero_rows_A.end(), gen);
    std::shuffle(zero_cols_B.begin(), zero_cols_B.end(), gen);

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

void init_test3_matrices_fp64( int32_t seed, const int32_t dim, double *A, const int32_t lda, double *B, const int32_t ldb ){

    int32_t i, j;
    int32_t max_exponent = 1000;
    int32_t min_exponent = -1000;
    std::vector<double> x(dim);

    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> generator(1.0, 2.0);
    std::vector<double> diag_D(dim);
    std::vector<double> diag_D_inv(dim);

    for ( i = 0; i < dim; ++i ){
        x[i] = generator(gen);
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

void init_jim_test3a_matrices_fp64( int32_t seed, const int32_t dim, double *A, const int32_t lda, double *B, const int32_t ldb ){

    int32_t i, j;
    int32_t max_exponent = 40;
    std::vector<double> x(dim);

    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> generator(1.0, 2.0);
    std::vector<double> diag_D(dim);
    std::vector<double> diag_D_inv(dim);

    for ( i = 0; i < dim; ++i ){
        x[i] = generator(gen);
    }

    for ( i = 0; i < dim; ++i ){
        diag_D[i] = std::ldexp(1.0, ((-1 * dim / 2) + 1 + i) * (max_exponent / ( dim / 2 ) ));
        diag_D_inv[i] = std::ldexp(1.0, -1 * (((-1 * dim / 2) + 1 + i) * (max_exponent / ( dim / 2 ) )));
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

void init_jim_test3b_matrices_fp64( int32_t seed, const int32_t dim, double *A, const int32_t lda, double *B, const int32_t ldb ){

    init_jim_test3a_matrices_fp64(seed, dim, A, lda, B, ldb);

    for ( int i = 1; i < dim; ++i ){
        for ( int j = 0; j < dim; ++j ){
            A[i + j*lda] = A[0 + ((dim + j - i)%dim)*lda];
        }
    }

    for ( int j = 1; j < dim; ++j ){
        for ( int i = 0; i < dim; ++i ){
            B[i + j*ldb] = B[(dim + i - j)%dim];
        }
    }
}

void init_jim_test3c_matrices_fp64( int32_t seed, const int32_t dim, double *A, const int32_t lda, double *B, const int32_t ldb ){

    std::vector<double> unpermuted_A(dim * dim);
    std::vector<double> unpermuted_B(dim * dim);

    init_jim_test3b_matrices_fp64(seed, dim, &unpermuted_A[0], dim, &unpermuted_B[0], dim);

    std::vector<uint32_t> p(dim);
    
    for ( int i = 0; i < dim; ++i ){
        p[i] = i;
    }

    std::mt19937 gen(seed);
    std::shuffle(p.begin(), p.end(), gen);

    for ( int i = 0; i < dim; ++i ){
        for ( int j = 0; j < dim; ++j ){
            A[i + j*lda] = unpermuted_A[p[i] + j*dim];
        }
    }

    std::shuffle(p.begin(), p.end(), gen);

    for ( int j = 0; j < dim; ++j ){
        for ( int i = 0; i < dim; ++i ){
            B[i + j*lda] = unpermuted_B[i + p[j]*dim];
        }
    }
}

template <typename T>
void init_matrices( std::string type, int32_t seed, const int32_t m, const int32_t n, const int32_t k, T *A, const int32_t lda, T *B, const int32_t ldb ){
    if ( type == "U_0_1" ){
        init_random_matrices( static_cast<T>(0.0), static_cast<T>(1.0), seed, m, n, k, A, lda, B, ldb );
    }
    else if ( type == "U_-1_1" ){
        init_random_matrices( static_cast<T>(-1.0), static_cast<T>(1.0), seed, m, n, k, A, lda, B, ldb );
    }
    else if ( type == "test1" ){
        init_test1_matrices( seed, m, n, k, A, lda, B, ldb );
    }
    else if ( type == "test3" ){
        static_assert( std::is_same<T, double>::value, "test3 matrices are only implemented for double precision." );
        assert( m == n && n == k && "test3 matrices are only implemented for square matrices.");
        init_test3_matrices_fp64( seed, m, A, lda, B, ldb );
    }
    else if ( type == "jim_test3a" ){
        assert( m == n && n == k && "test3 matrices are only implemented for square matrices.");
        assert( seed == 0 && "jim's test3 matrices are deterministically generated" );
        init_jim_test3a_matrices_fp64( seed, m, A, lda, B, ldb);
    }
    else if ( type == "jim_test3b" ){
        assert( m == n && n == k && "test3 matrices are only implemented for square matrices.");
        assert( seed == 0 && "jim's test3 matrices are deterministically generated" );
        init_jim_test3b_matrices_fp64( seed, m, A, lda, B, ldb);
    }
    else if ( type == "jim_test3c" ){
        assert( m == n && n == k && "test3 matrices are only implemented for square matrices.");
        assert( seed == 0 && "jim's test3 matrices are deterministically generated" );
        init_jim_test3c_matrices_fp64( seed, m, A, lda, B, ldb);
    }
    else if ( type == "sanity" ){
        init_sanity_matrices(  m, n, k, A, lda, B, ldb );
    }
    else{
        assert(false && "Unknown matrix type specified in init_matrices." );
    }

    return;
}

#endif