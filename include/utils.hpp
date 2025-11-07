#ifndef UTILS_H
#define UTILS_H

#include <stdint.h>
#include <quadmath.h>
#include <fstream>

template <typename T>
void write_matrix_bytes_to_file( std::ofstream& outfile, T *A, const uint32_t lda, const uint32_t m, const uint32_t n ){
    uint32_t i;
    for ( i = 0; i < m; ++i ){
        outfile.write(reinterpret_cast<char*>(&A[i*lda]), sizeof(T) * n);
    }
    return;
}

template <typename T>
void read_matrix_bytes_from_file( std::string filename, T *A, const uint32_t m, const uint32_t n ){
    std::ifstream infile(filename, std::ios::binary);
    infile.read(reinterpret_cast<char*>(A), sizeof(T) * m * n);
    infile.close();
    return;
}

template <typename T>
void print_matrix_bytes_to_stdout( T *A, const uint32_t lda, const uint32_t m, const uint32_t n ){
    uint32_t i, j;

    for ( i = 0; i < m; ++i ){
        for ( j = 0; j < n; ++j ){
            printf("%a\n", A[i*lda+j]);
        }
    }
    
    return;
}

template <typename T>
void print_matrix_to_stdout( T *A, const uint32_t lda, const uint32_t m, const uint32_t n ){
    uint32_t i, j;

    for ( i = 0; i < m; ++i ){
        for ( j = 0; j < n; ++j ){
            printf("%f ", A[i*lda+j]);
        }
        printf("\n");
    }
    printf("\n");
    
    return;
}

void print_matrix_128( __float128 *A, const uint32_t lda, const uint32_t m, const uint32_t n ){
    uint32_t i, j;
    char buf[10];
    int width = 46;

    for ( i = 0; i < m; ++i ){
        for ( j = 0; j < n; ++j ){
            quadmath_snprintf(buf, sizeof(buf), "%+-#*.2Qe", width, A[i*lda+j]);
            printf("%s ", buf);
        }
        printf("\n");
    }
    printf("\n");
    
    return;
}

template <typename T>
void get_elementwise_abs_errors( T *A, const uint32_t lda, T *B, const uint32_t ldb, const uint32_t m, const uint32_t n, T *C, const uint32_t ldc ){
    uint32_t i, j;

    for ( i = 0; i < m; ++i ){
        for ( j = 0; j < n; ++j ){
            C[i*ldc + j] = fabs(B[i*ldb + j] - A[i*lda + j]);
        }
    }

    return;
}

#endif