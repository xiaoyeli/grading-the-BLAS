extern "C" {
    void dgemm_(const char* TRANSA, const char* TRANSB,
                const int* M, const int* N, const int* K,
                const double* ALPHA,
                const double* A, const int* LDA,
                const double* B, const int* LDB,
                const double* BETA,
                double* C, const int* LDC);
}

void dgemm_(const char* transa, const char* transb,
                const int* m, const int* n, const int* k,
                const double* alpha,
                const double* A, const int* lda,
                const double* B, const int* ldb,
                const double* beta,
                double* C, const int* ldc){

    for ( int i = 0; i < *m; ++i ){
        for ( int j = 0; j < *n; ++j ){
            for ( int p = 0; p < *k; ++p ){
                C[i + j*(*ldc)] = A[i + p*(*lda)] * B[p + j*(*ldb)] + C[i + j*(*ldc)];
            }
        }
    }
}