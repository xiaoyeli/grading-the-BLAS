#ifndef REFERENCE_GEMM_HPP
#define REFERENCE_GEMM_HPP

#include <stdint.h>

template <typename T>
void reference_gemm( T *A, const uint32_t lda, T *B, const uint32_t ldb, uint32_t m, uint32_t k, uint32_t n, T *C, const uint32_t ldc ){
    
    uint32_t i, j, p;

    for ( i = 0; i < m; ++i ){
        for ( p = 0; p < k; ++p ){
            for ( j = 0; j < n; ++j ){
                C[i*ldc + j] = A[i*lda + p] * B[p*ldb + j] + C[i*ldc + j];
            }
        }
    }

    return;
}

#endif