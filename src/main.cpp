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
    double low_prec_buffer[m*n];
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

    // run classical gemm and cast output
    // to high precision for error calculation
    classical_gemm(A, B, m, k, n, C);
    copy_matrix(high_prec_buffer, C, m, n);

    // run classical gemm in high precision
    classical_gemm(Aq, Bq, m, k, n, Cq);
    
    // calculate errors
    get_elementwise_abs_errors( Cq, high_prec_buffer, m, n, errors );

    // print errors
    copy_matrix(low_prec_buffer, errors, m, n);
    print_matrix_bytes_to_stdout(low_prec_buffer, m, n);

    return 0;
}