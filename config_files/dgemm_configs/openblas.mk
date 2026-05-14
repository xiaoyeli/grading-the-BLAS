# === target DGEMM implementation ===
DGEMM_IMPL_NAME := OPENBLAS_0.3.30 # OPENBLAS_0.3.30_GNU_13.2_DEFAULT_x86_64_zen3
SO_PATH := $(OPENBLAS_SO_PATH)
LDFLAGS := -L$(dir $(SO_PATH)) -lopenblas
RUN_PREFIX :=                # optional; string prepended to execution command, i.e., "${RUN_PREFIX} ./run_test"

# === compilation tweaks ====
CXX := g++
USE_CUDA := 0	# if USE_CUDA=1, set CXX=nvcc
CXXFLAGS := 	# optional

# === advanced options ===
OUTER_SCALING := 0  # when OUTER_SCALING=1, calculate D1_inv * ( (D1 * A) * (B * D2) ) * D2_inv
                    # where, respectively, D1 and D2 yield rows and cols of A and B with unit norms

INNER_SCALING := 0  # when INNER_SCALING=1, calculate (A * D) * (D_inv * B)
                    # where D yields cols and rows of A and B respectively with roughly equal norms
