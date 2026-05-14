# === target DGEMM implementation ===
DGEMM_IMPL_NAME := NETLIBBLAS_3.12.1_GNU_13.2_DEFAULT_x86_64_zen3
LIB_PATH := $(NETLIBBLAS_SO_PATH)
LDFLAGS := -L$(dir $(LIB_PATH)) -lrefblas
RUN_PREFIX :=                # optional; string prepended to execution command, i.e., "${RUN_PREFIX} ./run_test"

# === compilation tweaks ====
CXX := g++
USE_CUDA := 0	# if USE_CUDA=1, set CXX=nvcc
CXXFLAGS := -Wl,-rpath,$(dir $(LIB_PATH))

# === advanced options ===
OUTER_SCALING := 0  # when OUTER_SCALING=1, calculate D1_inv * ( (D1 * A) * (B * D2) ) * D2_inv
                    # where, respectively, D1 and D2 yield rows and cols of A and B with unit norms

INNER_SCALING := 0  # when INNER_SCALING=1, calculate (A * D) * (D_inv * B)
                    # where D yields cols and rows of A and B respectively with roughly equal norms
