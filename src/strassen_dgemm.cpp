
#include <stdint.h>
#include <vector>
#include "reference_gemm.hpp"

extern "C" {
    void dgemm_(const char* TRANSA, const char* TRANSB,
                const int* M, const int* N, const int* K,
                const double* ALPHA,
                const double* A, const int* LDA,
                const double* B, const int* LDB,
                const double* BETA,
                double* C, const int* LDC);
}

void add_matrix( const double *A, const int32_t lda, const double *B, const int32_t ldb, double *C, const int32_t ldc, const int32_t m, const int32_t n ){
    int32_t i, j;
    for ( i = 0; i < m; ++i ){
        for ( j = 0; j < n; ++j ){
            C[i + j*ldc] = A[i + j*lda] + B[i + j*ldb];
        }
    }

    return;
}

void sub_matrix( const double *A, const int32_t lda, const double *B, const int32_t ldb, double *C, const int32_t ldc, const int32_t m, const int32_t n ){
    int32_t i, j;
    for ( i = 0; i < m; ++i ){
        for ( j = 0; j < n; ++j ){
            C[i + j*ldc] = A[i + j*lda] - B[i + j*ldb];
        }
    }

    return;
}

void dgemm_(const char* transa, const char* transb,
                const int* m, const int* n, const int* k,
                const double* alpha,
                const double* A, const int* lda,
                const double* B, const int* ldb,
                const double* beta,
                double* C, const int* ldc){
    
    // do a naive matmul if the dimensions of the operands don't match or are not evenly divisible by 2
    if ( (*m != *k) || (*k != *n) || (*m != *n) || ( *m % 2 == 1 ) ){
        reference_gemm( A, *lda, B, *ldb, *m, *k, *n, C, *ldc );
        return;
    }

    int32_t half_dim = *m / 2;
    std::vector<double> buffer1( half_dim * half_dim );
    std::vector<double> buffer2( half_dim * half_dim );
    std::vector<double> M1( half_dim * half_dim );
    std::vector<double> M2( half_dim * half_dim );
    std::vector<double> M3( half_dim * half_dim );
    std::vector<double> M4( half_dim * half_dim );
    std::vector<double> M5( half_dim * half_dim );
    std::vector<double> M6( half_dim * half_dim );
    std::vector<double> M7( half_dim * half_dim );
    
    add_matrix( &A[0], *lda, &A[half_dim + half_dim**lda], *lda, &buffer1[0], half_dim, half_dim, half_dim ); // buffer1 = A11 + A22
    add_matrix( &B[0], *ldb, &B[half_dim + half_dim**ldb], *ldb, &buffer2[0], half_dim, half_dim, half_dim ); // buffer2 = B11 + B22
    reference_gemm( &buffer1[0], half_dim, &buffer2[0], half_dim, half_dim, half_dim, half_dim, &M1[0], half_dim ); // M1 = buffer1 * buffer2
    
    add_matrix( &A[half_dim], *lda, &A[half_dim + half_dim**lda], *lda, &buffer1[0], half_dim, half_dim, half_dim ); // buffer1 = A21 + A22
    reference_gemm( &buffer1[0], half_dim, &B[0], *ldb, half_dim, half_dim, half_dim, &M2[0], half_dim ); // M2 = buffer1 * B11

    sub_matrix( &B[half_dim**ldb], *ldb, &B[half_dim + half_dim**ldb], *ldb, &buffer2[0], half_dim, half_dim, half_dim ); // buffer2 = B12 - B22
    reference_gemm( &A[0], *lda, &buffer2[0], half_dim, half_dim, half_dim, half_dim, &M3[0], half_dim ); // M3 = A11 * buffer2

    sub_matrix( &B[half_dim], *ldb, &B[0], *ldb, &buffer2[0], half_dim, half_dim, half_dim ); // buffer2 = B21 - B11
    reference_gemm( &A[half_dim + half_dim**lda], *lda, &buffer2[0], half_dim, half_dim, half_dim, half_dim, &M4[0], half_dim ); // M4 = A22 * buffer2
    
    add_matrix( &A[0], *lda, &A[half_dim**lda], *lda, &buffer1[0], half_dim, half_dim, half_dim ); // buffer1 = A11 + A12
    reference_gemm( &buffer1[0], half_dim, &B[half_dim + half_dim**ldb], *ldb, half_dim, half_dim, half_dim, &M5[0], half_dim ); // M5 = buffer1 * B22

    sub_matrix( &A[half_dim], *lda, &A[0], *lda, &buffer1[0], half_dim, half_dim, half_dim ); // buffer1 = A21 - A11
    add_matrix( &B[0], *ldb, &B[half_dim**ldb], *ldb, &buffer2[0], half_dim, half_dim, half_dim ); // buffer2 = B11 + B12
    reference_gemm( &buffer1[0], half_dim, &buffer2[0], half_dim, half_dim, half_dim, half_dim, &M6[0], half_dim ); // M6 = buffer1 * buffer2

    sub_matrix( &A[half_dim**lda], *lda, &A[half_dim + half_dim**lda], *lda, &buffer1[0], half_dim, half_dim, half_dim ); // buffer1 = A12 - A22
    add_matrix( &B[half_dim], *ldb, &B[half_dim + half_dim**ldb], *ldb, &buffer2[0], half_dim, half_dim, half_dim ); // buffer2 = B21 + B22
    reference_gemm( &buffer1[0], half_dim, &buffer2[0], half_dim, half_dim, half_dim, half_dim, &M7[0], half_dim ); // M7 = buffer1 * buffer2

    // C11
    add_matrix( &M1[0], half_dim, &M4[0], half_dim, &buffer1[0], half_dim, half_dim, half_dim ); // buffer1 = M1 + M4
    sub_matrix( &buffer1[0], half_dim, &M5[0], half_dim, &buffer2[0], half_dim, half_dim, half_dim ); // buffer2 = buffer1 - M5
    add_matrix( &buffer2[0], half_dim, &M7[0], half_dim, &C[0], *ldc, half_dim, half_dim ); // C11 = buffer2 + M7
    // C12
    add_matrix( &M3[0], half_dim, &M5[0], half_dim, &C[half_dim**ldc], *ldc, half_dim, half_dim ); // C12 = M3 + M5
    // C21
    add_matrix( &M2[0], half_dim, &M4[0], half_dim, &C[half_dim], *ldc, half_dim, half_dim ); // C21 = M2 + M4
    // C22
    sub_matrix( &M1[0], half_dim, &M2[0], half_dim, &buffer1[0], half_dim, half_dim, half_dim ); // buffer1 = M1 - M2
    add_matrix( &buffer1[0], half_dim, &M3[0], half_dim, &buffer2[0], half_dim, half_dim, half_dim ); // buffer2 = buffer1 + M3
    add_matrix( &buffer2[0], half_dim, &M6[0], half_dim, &C[half_dim + half_dim**ldc], *ldc, half_dim, half_dim ); // C22 = buffer2 + M6

    return;
}