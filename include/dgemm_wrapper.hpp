#ifndef DGEMM_WRAPPER_HPP
#define DGEMM_WRAPPER_HPP

void dgemm_wrapper( const char transa, const char transb,
                    const int m, const int n, const int k,
                    const double alpha,
                    double *A, const int lda,
                    double *B, const int ldb,
                    const double beta,
                    double *C, const int ldc );
#endif