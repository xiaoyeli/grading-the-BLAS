#ifndef REFERENCE_GEMM_HPP
#define REFERENCE_GEMM_HPP

#include <stdint.h>

template <typename T>
void reference_gemm( const T *A, const int32_t lda, const T *B, const int32_t ldb, const int32_t m, const int32_t k, const int32_t n, T *C, const int32_t ldc ){
    
    int32_t i, j, p;

    for ( i = 0; i < m; ++i ){
        for ( j = 0; j < n; ++j ){
            for ( p = 0; p < k; ++p ){
                C[i + j*ldc] = A[i + p*lda] * B[p + j*ldb] + C[i + j*ldc];
            }
        }
    }

    return;
}

#endif