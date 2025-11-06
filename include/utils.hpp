#ifndef UTILS_H
#define UTILS_H

#include<stdint.h>
#include <quadmath.h>

template <typename T>
void print_matrix( T *A, const uint32_t m, const uint32_t n ){
    uint32_t i, j;

    for ( i = 0; i < m; ++i ){
        for ( j = 0; j < n; ++j ){
            printf("%f ", A[i*n+j]);
        }
        printf("\n");
    }
    printf("\n");
    
    return;
}

void print_matrix_128( __float128 *A, const uint32_t m, const uint32_t n ){
    uint32_t i, j;
    char buf[10];
    int width = 46;

    for ( i = 0; i < m; ++i ){
        for ( j = 0; j < n; ++j ){
            quadmath_snprintf(buf, sizeof(buf), "%+-#*.2Qe", width, A[i*n+j]);
            printf("%s ", buf);
        }
        printf("\n");
    }
    printf("\n");
    
    return;
}

template <typename T>
void get_elementwise_abs_errors( T *A, T *B, const uint32_t m, const uint32_t n, T *C ){
    uint32_t i, j;

    for ( i = 0; i < m; ++i ){
        for ( j = 0; j < n; ++j ){
            C[i*n + j] = fabs(B[i*n + j] - A[i*n + j]);
        }
    }

    return;
}

#endif