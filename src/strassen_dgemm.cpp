#include <vector>
#include <cassert>

#define BASECASE_SIZE 2

extern "C" {
    void dgemm_(const char* TRANSA, const char* TRANSB,
                const int* M, const int* N, const int* K,
                const double* ALPHA,
                const double* A, const int* LDA,
                const double* B, const int* LDB,
                const double* BETA,
                double* C, const int* LDC);
}

void add_matrix( const double *A, const int lda, const double *B, const int ldb, double *C, const int ldc, const int m, const int n ){

    for ( int i = 0; i < m; ++i ){
        for ( int j = 0; j < n; ++j ){
            C[i + j*ldc] = A[i + j*lda] + B[i + j*ldb];
        }
    }

    return;
}

void sub_matrix( const double *A, const int lda, const double *B, const int ldb, double *C, const int ldc, const int m, const int n ){

    for ( int i = 0; i < m; ++i ){
        for ( int j = 0; j < n; ++j ){
            C[i + j*ldc] = A[i + j*lda] - B[i + j*ldb];
        }
    }

    return;
}

void naive_gemm( const double *A, const int lda, const double *B, const int ldb, const int m, const int k, const int n, double *C, const int ldc ){

    for ( int i = 0; i < m; ++i ){
        for ( int j = 0; j < n; ++j ){
            for ( int p = 0; p < k; ++p ){
                C[i + j*ldc] = A[i + p*lda] * B[p + j*ldb] + C[i + j*ldc];
            }
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
    
    assert((*m == *k) && (*k == *n) && (*m == *n) && ( *m % 2 == 0 ) && "Strassen DGEMM only implemented for square matrices whose dimensions are a power of 2");

    // base case
    if ( *m <= BASECASE_SIZE ){
        naive_gemm( A, *lda, B, *ldb, *m, *k, *n, C, *ldc );
        return;
    }

    const int half_dim = *m / 2;
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
    dgemm_( transa, transb, &half_dim, &half_dim, &half_dim, alpha, &buffer1[0], &half_dim, &buffer2[0], &half_dim, beta, &M1[0], &half_dim ); // M1 = buffer1 * buffer2
    
    add_matrix( &A[half_dim], *lda, &A[half_dim + half_dim**lda], *lda, &buffer1[0], half_dim, half_dim, half_dim ); // buffer1 = A21 + A22
    dgemm_( transa, transb, &half_dim, &half_dim, &half_dim, alpha, &buffer1[0], &half_dim, &B[0], ldb, beta, &M2[0], &half_dim ); // M2 = buffer1 * B11

    sub_matrix( &B[half_dim**ldb], *ldb, &B[half_dim + half_dim**ldb], *ldb, &buffer2[0], half_dim, half_dim, half_dim ); // buffer2 = B12 - B22
    dgemm_( transa, transb, &half_dim, &half_dim, &half_dim, alpha, &A[0], lda, &buffer2[0], &half_dim, beta, &M3[0], &half_dim ); // M3 = A11 * buffer2

    sub_matrix( &B[half_dim], *ldb, &B[0], *ldb, &buffer2[0], half_dim, half_dim, half_dim ); // buffer2 = B21 - B11
    dgemm_( transa, transb, &half_dim, &half_dim, &half_dim, alpha, &A[half_dim + half_dim**lda], lda, &buffer2[0], &half_dim, beta, &M4[0], &half_dim ); // M4 = A22 * buffer2
    
    add_matrix( &A[0], *lda, &A[half_dim**lda], *lda, &buffer1[0], half_dim, half_dim, half_dim ); // buffer1 = A11 + A12
    dgemm_( transa, transb, &half_dim, &half_dim, &half_dim, alpha, &buffer1[0], &half_dim, &B[half_dim + half_dim**ldb], ldb, beta, &M5[0], &half_dim ); // M5 = buffer1 * B22

    sub_matrix( &A[half_dim], *lda, &A[0], *lda, &buffer1[0], half_dim, half_dim, half_dim ); // buffer1 = A21 - A11
    add_matrix( &B[0], *ldb, &B[half_dim**ldb], *ldb, &buffer2[0], half_dim, half_dim, half_dim ); // buffer2 = B11 + B12
    dgemm_( transa, transb, &half_dim, &half_dim, &half_dim, alpha, &buffer1[0], &half_dim, &buffer2[0], &half_dim, beta, &M6[0], &half_dim ); // M6 = buffer1 * buffer2

    sub_matrix( &A[half_dim**lda], *lda, &A[half_dim + half_dim**lda], *lda, &buffer1[0], half_dim, half_dim, half_dim ); // buffer1 = A12 - A22
    add_matrix( &B[half_dim], *ldb, &B[half_dim + half_dim**ldb], *ldb, &buffer2[0], half_dim, half_dim, half_dim ); // buffer2 = B21 + B22
    dgemm_( transa, transb, &half_dim, &half_dim, &half_dim, alpha, &buffer1[0], &half_dim, &buffer2[0], &half_dim, beta, &M7[0], &half_dim ); // M7 = buffer1 * buffer2

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