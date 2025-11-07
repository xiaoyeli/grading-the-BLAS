`#include "init_matrices.hpp"
#include "utils.hpp"
#include "reference_gemm.hpp"
#include <stdio.h>
#include <stdint.h>
#include <quadmath.h>
#include <fstream>

int main(){
    uint32_t m, k, n;
    uint32_t i, j, p;
    
    m = 1024;
    k = 1024;
    n = 1024;

    double A[m*k];
    double B[k*n];
    double C[m*n];
    double low_prec_buffer[m*n];
    __float128 high_prec_buffer[m*n];
    __float128 Aq[m*k];
    __float128 Bq[k*n];
    __float128 Cq[m*n];
    __float128 errors[m*n];
    
    std::ofstream outfile1("reference_errors.bin", std::ios::binary);
    for ( i = 128; i <= m; i = 2*i ){
        for ( j = 128; j <= n; j = 2*j ){
            for ( p = 128; p <= k; p = 2*p ){    

                // initialize matrices
                init_test1_matrices(A, k, B, n, i, p, j);
                init_zero_matrix(C, n, i, j);
                copy_matrix(Aq, k, A, k, i, p);
                copy_matrix(Bq, n, B, n, p, j);
                init_zero_matrix(Cq, n, i, j);
            
                // run reference gemm in high precision
                reference_gemm(Aq, k, Bq, n, i, p, j, Cq, n);
                
                // run reference gemm and calculate errors
                reference_gemm(A, k, B, n, i, p, j, C, n);
                copy_matrix(high_prec_buffer, n, C, n, i, j);
                get_elementwise_abs_errors(Cq, n, high_prec_buffer, n, i, j, errors, n);
                copy_matrix(low_prec_buffer, n, errors, n, i, j);
                write_matrix_bytes_to_file(outfile1, low_prec_buffer, n, i, j);
            }    
        }    
    }    
    outfile1.close();
    
    return 0;
}