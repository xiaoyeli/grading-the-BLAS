#include <cuda_runtime.h>
#include <cublas_v2.h>

template <typename T>
void cublas_dgemm_wrapper( const int32_t m, const int32_t n, const int32_t k, const T *A, const int32_t lda, const T *B, const int32_t ldb, T *C, const int32_t ldc ){

    cublasHandle_t handle;
    cublasCreate(&handle);
    
    T alpha = 1.0;
    T beta = 0.0;
    
    cublasOperation_t opA = CUBLAS_OP_N;
    cublasOperation_t opB = CUBLAS_OP_N;

    T *d_A, *d_B, *d_C;

#ifdef ADDB_SCALING
    std::vector<double> A_scaled(m * k);
    std::vector<double> B_scaled(k * n);

    std::vector<double> A_col_norms(m);
    std::vector<double> B_row_norms(m);

    // find inf norms of cols/rows of A and B respectively
    for ( int i = 0; i < m; ++i ){
        A_col_norms[i] = std::abs(A[0 + i * lda]);
        B_row_norms[i] = std::abs(B[i + 0 * ldb]);
        for ( int j = 1; j < m; ++j ){
            A_col_norms[i] = std::abs(A[j + i*lda]) > A_col_norms[i] ? std::abs(A[j + i*lda]) : A_col_norms[i];
            B_row_norms[i] = std::abs(B[i + j*ldb]) > B_row_norms[i] ? std::abs(B[i + j*ldb]) : B_row_norms[i];
        }
    }

    // scaling
    for ( int i = 0; i < m; ++i ){
        double scaling_factor = A_col_norms[i] * B_row_norms[i] != 0 ? std::ldexp(1.0, std::round((std::log2(B_row_norms[i]) - std::log2(A_col_norms[i])) / 2.0)) : 1.0;        
        for ( int j = 0; j < m; ++j ){
            A_scaled[j + i*m] = A[j + i*lda] * scaling_factor;
            B_scaled[i + j*k] = B[i + j*ldb] / scaling_factor;
        }
    }
    size_t size_A = m * k * sizeof(T);
    size_t size_B = k * n * sizeof(T);
    size_t size_C = ldc * n * sizeof(T);

    cudaMalloc((void**)&d_A, size_A);
    cudaMalloc((void**)&d_B, size_B);
    cudaMalloc((void**)&d_C, size_C);
    
    cudaMemcpy(d_A, &A_scaled[0], size_A, cudaMemcpyHostToDevice);
    cudaMemcpy(d_B, &B_scaled[0], size_B, cudaMemcpyHostToDevice);
    cudaMemcpy(d_C, C, size_C, cudaMemcpyHostToDevice);

    cublasDgemm(handle, opA, opB, m, n, k, &alpha, d_A, m, d_B, k, &beta, d_C, ldc);
#elif defined(DABD_SCALING)
    std::vector<double> A_scaled(m * k);
    std::vector<double> B_scaled(k * n);

    std::vector<double> D1_diag(m);
    std::vector<double> D2_diag(m);

    // init diags
    for ( int i = 0; i < m; ++i ){
        D1_diag[i] = std::abs(A[i]);
        D2_diag[i] = std::abs(B[i * ldb]);
    }

    // inf norms of rows/cols of A and B respectively
    for ( int i = 0; i < m; ++i ){
        for ( int j = 1; j < m; ++j ){
            D1_diag[i] = std::abs(A[i + j*lda]) > D1_diag[i] ? std::abs(A[i + j*lda]) : D1_diag[i];
            D2_diag[i] = std::abs(B[j + i*ldb]) > D2_diag[i] ? std::abs(B[j + i*ldb]) : D2_diag[i];
        }
    }

    // scale to yield unit norms
    for ( int i = 0; i < m; ++i ){
        for ( int j = 0; j < m; ++j ){
            A_scaled[i + j*m] = D1_diag[i] != 0.0 ? A[i + j*lda] / D1_diag[i] : 0.0;
            B_scaled[j + i*k] = D2_diag[i] != 0.0 ? B[j + i*ldb] / D2_diag[i] : 0.0;
        }
    }
    size_t size_A = m * k * sizeof(T);
    size_t size_B = k * n * sizeof(T);
    size_t size_C = ldc * n * sizeof(T);

    cudaMalloc((void**)&d_A, size_A);
    cudaMalloc((void**)&d_B, size_B);
    cudaMalloc((void**)&d_C, size_C);
    
    cudaMemcpy(d_A, &A_scaled[0], size_A, cudaMemcpyHostToDevice);
    cudaMemcpy(d_B, &B_scaled[0], size_B, cudaMemcpyHostToDevice);
    cudaMemcpy(d_C, C, size_C, cudaMemcpyHostToDevice);

    cublasDgemm(handle, opA, opB, m, n, k, &alpha, d_A, m, d_B, k, &beta, d_C, ldc);
#else
    size_t size_A = lda * k * sizeof(T);
    size_t size_B = ldb * n * sizeof(T);
    size_t size_C = ldc * n * sizeof(T);

    cudaMalloc((void**)&d_A, size_A);
    cudaMalloc((void**)&d_B, size_B);
    cudaMalloc((void**)&d_C, size_C);
    
    cudaMemcpy(d_A, A, size_A, cudaMemcpyHostToDevice);
    cudaMemcpy(d_B, B, size_B, cudaMemcpyHostToDevice);
    cudaMemcpy(d_C, C, size_C, cudaMemcpyHostToDevice);

    cublasDgemm(handle, opA, opB, m, n, k, &alpha, d_A, lda, d_B, ldb, &beta, d_C, ldc);
#endif

    cudaMemcpy(C, d_C, size_C, cudaMemcpyDeviceToHost);

    cudaFree(d_A);
    cudaFree(d_B);
    cudaFree(d_C);

    cublasDestroy(handle);

#ifdef FP64_DIAGONAL
    for ( int i = 0; i < m; ++i ){
        double C_ii = 0.0;
        for ( int j = 0; j < n; ++j ){
            C_ii = A[i + j*lda] * B[j + i*ldb] + C_ii;
        }
        C[i + i*ldc] = C_ii;
    }
#elif defined(DABD_SCALING)
    for ( int i = 0; i < m; ++i ){
        for ( int j = 0; j < m; ++j ){
            C[i + j*ldc] = D1_diag[i] * C[i + j*ldc] * D2_diag[j];
        }
    }
#endif

    return;
}