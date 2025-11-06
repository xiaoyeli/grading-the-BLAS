#include "init_matrices.hpp"
#include "utils.hpp"
#include "classical_gemm.hpp"
#include <stdio.h>
#include <stdint.h>
#include <quadmath.h>

int main(){
    uint32_t m, k, n;

    // problem dimensions
    m = 16;
    k = 16;
    n = 16;

    double A[m*k];
    double B[k*n];
    double C[m*n];
    __float128 high_prec_buffer[m*n];
    __float128 Aq[m*k];
    __float128 Bq[k*n];
    __float128 Cq[m*n];
    __float128 errors[m*n];
    
    // initialize matrices
    init_test1_matrices(A, B, m, k, n);
    init_zero_matrix(C, m, n);
    copy_matrix(Aq, A, m, k);
    copy_matrix(Bq, B, k, n);
    init_zero_matrix(Cq, m, n);

    // run classical gemm in quad precision
    classical_gemm(Aq, Bq, m, k, n, Cq);
    
    // run classical gemm in double precision and cast
    // to quad precision for error calculation
    classical_gemm(A, B, m, k, n, C);
    copy_matrix(high_prec_buffer, C, m, n);

    // calculate errors
    get_elementwise_abs_errors( Cq, high_prec_buffer, m, n, errors );

    // print errors
    print_matrix_128(errors, m, n);

    return 0;
}