#include <vector>
#include <cassert>

#if defined(DABD_SCALING) || defined(ADDB_SCALING)
#include <cmath>
#endif

#ifdef ZEROING
#define GAMEFACTOR 1000
#endif

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

void naive_dgemm( const double *A, const int lda, const double *B, const int ldb, const int m, const int k, const int n, double *C, const int ldc ){

    for ( int i = 0; i < m; ++i ){
        for ( int j = 0; j < n; ++j ){
            for ( int p = 0; p < k; ++p ){
                C[i + j*ldc] = A[i + p*lda] * B[p + j*ldb] + C[i + j*ldc];
            }
        }
    }

    return;
}

void strassen_dgemm(const char* transa, const char* transb,
                const int* m, const int* n, const int* k,
                const double* alpha,
                const double* A, const int* lda,
                const double* B, const int* ldb,
                const double* beta,
                double* C, const int* ldc){
    
    // base case
    if ( *m <= BASECASE_SIZE ){
        naive_dgemm( A, *lda, B, *ldb, *m, *k, *n, C, *ldc );
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
    strassen_dgemm( transa, transb, &half_dim, &half_dim, &half_dim, alpha, &buffer1[0], &half_dim, &buffer2[0], &half_dim, beta, &M1[0], &half_dim ); // M1 = buffer1 * buffer2
    
    add_matrix( &A[half_dim], *lda, &A[half_dim + half_dim**lda], *lda, &buffer1[0], half_dim, half_dim, half_dim ); // buffer1 = A21 + A22
    strassen_dgemm( transa, transb, &half_dim, &half_dim, &half_dim, alpha, &buffer1[0], &half_dim, &B[0], ldb, beta, &M2[0], &half_dim ); // M2 = buffer1 * B11

    sub_matrix( &B[half_dim**ldb], *ldb, &B[half_dim + half_dim**ldb], *ldb, &buffer2[0], half_dim, half_dim, half_dim ); // buffer2 = B12 - B22
    strassen_dgemm( transa, transb, &half_dim, &half_dim, &half_dim, alpha, &A[0], lda, &buffer2[0], &half_dim, beta, &M3[0], &half_dim ); // M3 = A11 * buffer2

    sub_matrix( &B[half_dim], *ldb, &B[0], *ldb, &buffer2[0], half_dim, half_dim, half_dim ); // buffer2 = B21 - B11
    strassen_dgemm( transa, transb, &half_dim, &half_dim, &half_dim, alpha, &A[half_dim + half_dim**lda], lda, &buffer2[0], &half_dim, beta, &M4[0], &half_dim ); // M4 = A22 * buffer2
    
    add_matrix( &A[0], *lda, &A[half_dim**lda], *lda, &buffer1[0], half_dim, half_dim, half_dim ); // buffer1 = A11 + A12
    strassen_dgemm( transa, transb, &half_dim, &half_dim, &half_dim, alpha, &buffer1[0], &half_dim, &B[half_dim + half_dim**ldb], ldb, beta, &M5[0], &half_dim ); // M5 = buffer1 * B22

    sub_matrix( &A[half_dim], *lda, &A[0], *lda, &buffer1[0], half_dim, half_dim, half_dim ); // buffer1 = A21 - A11
    add_matrix( &B[0], *ldb, &B[half_dim**ldb], *ldb, &buffer2[0], half_dim, half_dim, half_dim ); // buffer2 = B11 + B12
    strassen_dgemm( transa, transb, &half_dim, &half_dim, &half_dim, alpha, &buffer1[0], &half_dim, &buffer2[0], &half_dim, beta, &M6[0], &half_dim ); // M6 = buffer1 * buffer2

    sub_matrix( &A[half_dim**lda], *lda, &A[half_dim + half_dim**lda], *lda, &buffer1[0], half_dim, half_dim, half_dim ); // buffer1 = A12 - A22
    add_matrix( &B[half_dim], *ldb, &B[half_dim + half_dim**ldb], *ldb, &buffer2[0], half_dim, half_dim, half_dim ); // buffer2 = B21 + B22
    strassen_dgemm( transa, transb, &half_dim, &half_dim, &half_dim, alpha, &buffer1[0], &half_dim, &buffer2[0], &half_dim, beta, &M7[0], &half_dim ); // M7 = buffer1 * buffer2

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


void dgemm_(const char* transa, const char* transb,
                const int* m, const int* n, const int* k,
                const double* alpha,
                const double* A, const int* lda,
                const double* B, const int* ldb,
                const double* beta,
                double* C, const int* ldc){

    assert((*m == *k) && (*k == *n) && (*m == *n) && ( *m % 2 == 0 ) && "Strassen DGEMM only implemented for square matrices whose dimensions are a power of 2");

#ifdef DABD_SCALING
    std::vector<double> A_scaled(*m * *k);
    std::vector<double> B_scaled(*k * *n);

    std::vector<double> D1_diag(*m);
    std::vector<double> D2_diag(*m);

    // init diags
    for ( int i = 0; i < *m; ++i ){
        D1_diag[i] = std::abs(A[i]);
        D2_diag[i] = std::abs(B[i * (*ldb)]);
    }

    // inf norms of rows/cols of A and B respectively
    for ( int i = 0; i < *m; ++i ){
        for ( int j = 1; j < *m; ++j ){
            D1_diag[i] = std::abs(A[i + j*(*lda)]) > D1_diag[i] ? std::abs(A[i + j*(*lda)]) : D1_diag[i];
            D2_diag[i] = std::abs(B[j + i*(*ldb)]) > D2_diag[i] ? std::abs(B[j + i*(*ldb)]) : D2_diag[i];
        }
    }

    // scale to yield unit norms
    for ( int i = 0; i < *m; ++i ){
        for ( int j = 0; j < *m; ++j ){
            A_scaled[i + j*(*m)] = D1_diag[i] != 0.0 ? A[i + j*(*lda)] / D1_diag[i] : 0.0;
            B_scaled[j + i*(*k)] = D2_diag[i] != 0.0 ? B[j + i*(*ldb)] / D2_diag[i] : 0.0;
        }
    }
#endif
#ifdef ADDB_SCALING
    std::vector<double> A_scaled(*m * *k);
    std::vector<double> B_scaled(*k * *n);

    std::vector<double> A_col_norms(*m);
    std::vector<double> B_row_norms(*m);

    // find inf norms of cols/rows of A and B respectively
    for ( int i = 0; i < *m; ++i ){
        A_col_norms[i] = std::abs(A[0 + i * (*lda)]);
        B_row_norms[i] = std::abs(B[i + 0 * (*ldb)]);
        for ( int j = 1; j < *m; ++j ){
            A_col_norms[i] = std::abs(A[j + i*(*lda)]) > A_col_norms[i] ? std::abs(A[j + i*(*lda)]) : A_col_norms[i];
            B_row_norms[i] = std::abs(B[i + j*(*ldb)]) > B_row_norms[i] ? std::abs(B[i + j*(*ldb)]) : B_row_norms[i];
        }
    }

    // scaling
    for ( int i = 0; i < *m; ++i ){
        double scaling_factor = A_col_norms[i] * B_row_norms[i] != 0 ? std::ldexp(1.0, std::round((std::log2(B_row_norms[i]) - std::log2(A_col_norms[i])) / 2.0)) : 1.0;        
        for ( int j = 0; j < *m; ++j ){
            A_scaled[j + i*(*m)] = A[j + i*(*lda)] * scaling_factor;
            B_scaled[i + j*(*k)] = B[i + j*(*ldb)] / scaling_factor;
        }
    }

#endif

#if defined(ADDB_SCALING) || defined(DABD_SCALING)
    strassen_dgemm(transa, transb, m, n, k, alpha, &A_scaled[0], m, &B_scaled[0], k, beta, C, ldc);
#else
    strassen_dgemm(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
#endif

#ifdef DABD_SCALING
#ifdef ZEROING
    double C_abs_max = 0.0;
#endif
    for ( int i = 0; i < *m; ++i ){
        for ( int j = 0; j < *m; ++j ){
            C[i + j*(*ldc)] = D1_diag[i] * C[i + j*(*ldc)] * D2_diag[j];
#ifdef ZEROING
            C_abs_max = std::abs(C[i + j*(*ldc)]) > C_abs_max ? std::abs(C[i + j*(*ldc)]) : C_abs_max;
#endif
        }
    }
#ifdef ZEROING
    for ( int i = 0; i < *m; ++i ){
        for ( int j = 0; j < *m; ++j ){
            C[i + j*(*ldc)] = std::abs(C[i + j*(*ldc)]) < GAMEFACTOR * std::numeric_limits<double>::epsilon() * C_abs_max ? 0.0 : C[i + j*(*ldc)];
        }
    }
#endif
#endif
}