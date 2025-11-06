#ifndef CLASSICAL_GEMM_H
#define CLASSICAL_GEMM_H

#include <stdint.h>

template <typename T>
void classical_gemm( T *A, T *B, uint32_t m, uint32_t k, uint32_t n, T *C ){
    
    uint32_t i, j, p;

    for ( i = 0; i < m; ++i ){
        for ( p = 0; p < k; ++p ){
            for ( j = 0; j < n; ++j ){
                C[i*n + j] = A[i*n + p] * B[p*n + j] + C[i*n + j];
            }
        }
    }

    return;
}

#endif