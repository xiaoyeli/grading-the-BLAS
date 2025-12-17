#ifndef UTILS_HPP
#define UTILS_HPP

#include <stdint.h>
#include <fstream>

template <typename T>
void write_matrix_bytes_to_file( std::ofstream& outfile, T *A, const int32_t lda, const int32_t m, const int32_t n ){
    int32_t j;

    // write out column by column
    for ( j = 0; j < n; ++j ){
        outfile.write( reinterpret_cast<char*>(&A[j*lda]), sizeof(T) * m );
    }
    return;
}

template <typename T>
void read_matrix_bytes_from_file( std::ifstream& infile, T *A, const int32_t lda, const int32_t m, const int32_t n ){
    int32_t j;

    // read in column by column
    for ( j = 0; j < n; ++j ){
        infile.read( reinterpret_cast<char*>(&A[j*lda]), sizeof(T) * m );
    }
    return;
}

template <typename T>
void print_matrix_bytes_to_stdout( const T *A, const int32_t lda, const int32_t m, const int32_t n ){
    int32_t i, j;

    for ( i = 0; i < m; ++i ){
        for ( j = 0; j < n; ++j ){
            printf("%a\n", A[i+j*lda]);
        }
    }
    
    return;
}

template <typename T>
void print_matrix_to_stdout( const T *A, const int32_t lda, const int32_t m, const int32_t n ){
    int32_t i, j;

    for ( i = 0; i < m; ++i ){
        for ( j = 0; j < n; ++j ){
            printf("%f ", A[i+j*lda]);
        }
        printf("\n");
    }
    printf("\n");
    
    return;
}

template <typename T>
void postprocess( const T *A_ref, const int32_t lda_ref, const T *A, const int32_t lda, const int32_t m, const int32_t n, T *signed_absolute_errors, const int32_t lde, double * max_componentwise_unsigned_relative_error, double * max_componentwise_unsigned_absolute_error, int32_t *n_zeros ){
    int32_t i, j;

    double unsigned_relative_error;
    double unsigned_absolute_error;

    *max_componentwise_unsigned_relative_error = 0.0;
    *max_componentwise_unsigned_absolute_error = 0.0;
    *n_zeros = 0;
    for ( i = 0; i < m; ++i ){
        for ( j = 0; j < n; ++j ){
            signed_absolute_errors[i + j*lde] = A_ref[i + j*lda_ref] - A[i + j*lda];
            unsigned_absolute_error = fabs(signed_absolute_errors[i + j*lde]);
            unsigned_relative_error = unsigned_absolute_error / fabs(A_ref[i + j*lda_ref]);
            if ( unsigned_absolute_error > *max_componentwise_unsigned_absolute_error ) {
                *max_componentwise_unsigned_absolute_error = unsigned_absolute_error;
            }
            if ( unsigned_relative_error > *max_componentwise_unsigned_relative_error ){
                *max_componentwise_unsigned_relative_error = unsigned_relative_error;
            }
            if ( A[i + j*lda] == 0.0 ){
                (*n_zeros)++;
            }
        }
    }
    return;
}

template <typename T_to, typename T_from>
void copy_matrix( T_to *A, const int32_t lda, const T_from *B, const int32_t ldb, const int32_t m, const int32_t n ){
    
    int32_t i, j;
    for ( i = 0; i < m; ++i ){
        for ( j = 0; j < n; ++j ){
            A[i + j*lda] = B[i + j*ldb];
        }
    }

    return;
}

#endif