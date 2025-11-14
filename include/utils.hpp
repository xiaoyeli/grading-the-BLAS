#ifndef UTILS_HPP
#define UTILS_HPP

#include <stdint.h>
#include <quadmath.h>
#include <fstream>

template <typename T>
void write_matrix_bytes_to_file( std::ofstream& outfile, T *A, const int32_t lda, const int32_t m, const int32_t n ){
    int32_t i;
    for ( i = 0; i < n; ++i ){
        outfile.write(reinterpret_cast<char*>(&A[i*lda]), sizeof(T) * m);
    }
    return;
}

template <typename T>
void read_matrix_bytes_from_file( const std::string filename, const T *A, const int32_t m, const int32_t n ){
    std::ifstream infile(filename, std::ios::binary);
    infile.read(reinterpret_cast<char*>(A), sizeof(T) * m * n);
    infile.close();
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

void print_matrix_128( const __float128 *A, const int32_t lda, const int32_t m, const int32_t n ){
    int32_t i, j;
    char buf[10];
    int width = 46;

    for ( i = 0; i < m; ++i ){
        for ( j = 0; j < n; ++j ){
            quadmath_snprintf(buf, sizeof(buf), "%+-#*.2Qe", width, A[i+j*lda]);
            printf("%s ", buf);
        }
        printf("\n");
    }
    printf("\n");
    
    return;
}

template <typename T>
void get_elementwise_abs_errors( const T *A, const int32_t lda, const T *B, const int32_t ldb, const int32_t m, const int32_t n, T *C, const int32_t ldc ){
    int32_t i, j;

    for ( i = 0; i < m; ++i ){
        for ( j = 0; j < n; ++j ){
            C[i + j*ldc] = fabs(B[i + j*ldb] - A[i + j*lda]);
        }
    }

    return;
}

template <typename T>
void get_elementwise_errors( const T *A, const int32_t lda, const T *B, const int32_t ldb, const int32_t m, const int32_t n, T *C, const int32_t ldc ){
    int32_t i, j;

    for ( i = 0; i < m; ++i ){
        for ( j = 0; j < n; ++j ){
            C[i + j*ldc] = B[i + j*ldb] - A[i + j*lda];
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