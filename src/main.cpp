#include "init_matrices.hpp"
#include "utils.hpp"
#include "reference_gemm.hpp"
#include <stdio.h>
#include <stdint.h>
#include <quadmath.h>
#include <fstream>
#include <vector>

int main(){
    uint32_t m, k, n;
    uint32_t i, j, p;
    
#ifdef SANITY
    m = n = k = i = j = p = 2;
#else
    m = n = k = 1024;
    i = j = p = 128;
#endif

    std::vector<double> A(m*k);
    std::vector<double> B(k*n);
    std::vector<double> C(m*n);
    std::vector<double> low_prec_buffer(m*n);
    std::vector<__float128> high_prec_buffer(m*n);
    std::vector<__float128> Aq(m*k);
    std::vector<__float128> Bq(k*n);
    std::vector<__float128> Cq(m*n);
    std::vector<__float128> errors(m*n);
    
    std::ofstream outfile1("reference_errors.bin", std::ios::binary);
    for ( ; i <= m; i = 2*i ){
        for ( ; j <= n; j = 2*j ){
            for ( ; p <= k; p = 2*p ){    

                // initialize matrices
#ifdef SANITY
                init_sanity_matrices(&A[0], k, &B[0], n, i, p, j);
#else
                init_test1_matrices(&A[0], k, &B[0], n, i, p, j);
#endif
                init_zero_matrix(&C[0], n, i, j);
                copy_matrix(&Aq[0], k, &A[0], k, i, p);
                copy_matrix(&Bq[0], n, &B[0], n, p, j);
                init_zero_matrix(&Cq[0], n, i, j);
            
#ifdef SANITY
                printf("A:\n");
                print_matrix_to_stdout(&A[0], k, i, p);
                printf("B:\n");
                print_matrix_to_stdout(&B[0], n, p, j);
#endif

                // run reference gemm in high precision
                reference_gemm(&Aq[0], k, &Bq[0], n, i, p, j, &Cq[0], n);

                // run reference gemm and calculate errors
                reference_gemm(&A[0], k, &B[0], n, i, p, j, &C[0], n);
                copy_matrix(&high_prec_buffer[0], n, &C[0], n, i, j);
                get_elementwise_abs_errors(&Cq[0], n, &high_prec_buffer[0], n, i, j, &errors[0], n);
                copy_matrix(&low_prec_buffer[0], n, &errors[0], n, i, j);
                write_matrix_bytes_to_file(outfile1, &low_prec_buffer[0], n, i, j);

#ifdef SANITY
                printf("C:\n");
                print_matrix_to_stdout(&C[0], n, i, j);
#endif
            }    
        }    
    }    
    outfile1.close();
    
    return 0;
}