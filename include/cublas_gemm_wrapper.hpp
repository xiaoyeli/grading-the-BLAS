#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <dlfcn.h>

template <typename T>
void cublas_gemm_wrapper( const T *A, const int32_t lda, const T *B, const int32_t ldb, const int32_t m, const int32_t k, const int32_t n, T *C, const int32_t ldc ){
    cublasHandle_t handle;
    cublasCreate(&handle);

    T *d_A, *d_B, *d_C;
    size_t size_A = lda * k * sizeof(T);
    size_t size_B = ldb * n * sizeof(T);
    size_t size_C = ldc * n * sizeof(T);

    cudaMalloc((void**)&d_A, size_A);
    cudaMalloc((void**)&d_B, size_B);
    cudaMalloc((void**)&d_C, size_C);

    cudaMemcpy(d_A, A, size_A, cudaMemcpyHostToDevice);
    cudaMemcpy(d_B, B, size_B, cudaMemcpyHostToDevice);
    cudaMemcpy(d_C, C, size_C, cudaMemcpyHostToDevice);

    T alpha = 1.0;
    T beta = 0.0;

    cublasOperation_t opA = CUBLAS_OP_N;
    cublasOperation_t opB = CUBLAS_OP_N;

    cublasDgemm(handle, opA, opB, m, n, k, &alpha, d_A, lda, d_B, ldb, &beta, d_C, ldc);

    cudaMemcpy(C, d_C, size_C, cudaMemcpyDeviceToHost);

    cudaFree(d_A);
    cudaFree(d_B);
    cudaFree(d_C);

    cublasDestroy(handle);

    return;
}

typedef void (*cublas_dgemm_t)(
    cublasHandle_t,
    cublasOperation_t,
    cublasOperation_t,
    int,
    int,
    int,
    double*,
    double*,
    int,
    double*,
    int,
    double*,
    double*,
    int
);

template <typename T>
void cublas_ozimmu_wrapper( const T *A, const int32_t lda, const T *B, const int32_t ldb, const int32_t m, const int32_t k, const int32_t n, T *C, const int32_t ldc ){

    static_assert( std::is_same<T, double>::value, "ozimmu only implemented for double precision." );

    void *handle_ozimmu;
    handle_ozimmu = dlopen("libozimmu.so", RTLD_LAZY);
    cublas_dgemm_t ozimmu = (cublas_dgemm_t)dlsym(handle_ozimmu, "cublasDgemm_v2");
    
    cublasHandle_t handle;
    cublasCreate(&handle);

    T *d_A, *d_B, *d_C;
    size_t size_A = lda * k * sizeof(T);
    size_t size_B = ldb * n * sizeof(T);
    size_t size_C = ldc * n * sizeof(T);

    cudaMalloc((void**)&d_A, size_A);
    cudaMalloc((void**)&d_B, size_B);
    cudaMalloc((void**)&d_C, size_C);

    cudaMemcpy(d_A, A, size_A, cudaMemcpyHostToDevice);
    cudaMemcpy(d_B, B, size_B, cudaMemcpyHostToDevice);
    cudaMemcpy(d_C, C, size_C, cudaMemcpyHostToDevice);

    T alpha = 1.0;
    T beta = 0.0;

    cublasOperation_t opA = CUBLAS_OP_N;
    cublasOperation_t opB = CUBLAS_OP_N;

    ozimmu(handle, opA, opB, m, n, k, &alpha, d_A, lda, d_B, ldb, &beta, d_C, ldc);

    cudaMemcpy(C, d_C, size_C, cudaMemcpyDeviceToHost);

    cudaFree(d_A);
    cudaFree(d_B);
    cudaFree(d_C);

    cublasDestroy(handle);

    return;
}