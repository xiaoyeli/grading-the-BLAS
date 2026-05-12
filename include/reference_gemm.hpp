#ifndef REFERENCE_GEMM_HPP
#define REFERENCE_GEMM_HPP

template <typename T>
void reference_gemm( const char* transa, const char* transb,
                const int* m, const int* n, const int* k,
                const double* alpha,
                const double* A, const int* lda,
                const double* B, const int* ldb,
                const double* beta,
                T* C, const int* ldc ){

    for ( int i = 0; i < *m; ++i ){
        for ( int j = 0; j < *n; ++j ){
            for ( int p = 0; p < *k; ++p ){
                C[i + j*(*ldc)] = A[i + p*(*lda)] * B[p + j*(*ldb)] + C[i + j*(*ldc)];
            }
        }
    }
}

#endif