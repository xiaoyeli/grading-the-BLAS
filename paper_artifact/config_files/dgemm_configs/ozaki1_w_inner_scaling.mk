# === target DGEMM implementation ===
DGEMM_IMPL_NAME := OZIMMU_08eea92_GNU_13.2_NVCC_12.9_DEFAULT_NVCC8.0_MAX_SLICES
LIB_PATH := $(CUBLAS_SO_PATH)
LDFLAGS := -L$(dir $(LIB_PATH)) -lcublas
RUN_PREFIX := LD_PRELOAD=$(OZIMMU_SO_PATH) OZIMMU_COMPUTE_MODE=fp64_int8_18 OZIMMU_INTERCEPT_THRESHOLD_M=2 OZIMMU_INTERCEPT_THRESHOLD_N=2 OZIMMU_INTERCEPT_THRESHOLD_K=2               # optional; string prepended to execution command, i.e., "${RUN_PREFIX} ./run_test"

# === compilation tweaks ====
CXX := nvcc
USE_CUDA := 1	# if USE_CUDA=1, set CXX=nvcc
CXXFLAGS :=

# === advanced options ===
OUTER_SCALING := 0  # when OUTER_SCALING=1, calculate D1_inv * ( (D1 * A) * (B * D2) ) * D2_inv
                    # where, respectively, D1 and D2 yield rows and cols of A and B with unit norms

INNER_SCALING := 1  # when INNER_SCALING=1, calculate (A * D) * (D_inv * B)
                    # where D yields cols and rows of A and B respectively with roughly equal norms
