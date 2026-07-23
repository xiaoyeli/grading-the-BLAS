# === target DGEMM implementation ===
DGEMM_IMPL_NAME := oneMKL_2026.0_OpenMP # e.g.,
# MKLROOT is set by sourcing  setvars.sh script in oneMKL installation, or you can explicitly set the path here
MY_MKLROOT := $(MKLROOT)

# This is oneMKL linkline for dynamic linking, LP64, gcc compiler, OpenMP multithreaded version.
LDFLAGS := -m64 -L$(MY_MKLROOT)/lib -Wl,-rpath,$(MY_MKLROOT)/lib -Wl,--no-as-needed -lmkl_intel_lp64 \
-lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl
RUN_PREFIX :=                # optional; string prepended to execution command, i.e., "${RUN_PREFIX} ./run_test"

# === compilation tweaks ====
CXX := g++
USE_CUDA := 0	# if USE_CUDA=1, set CXX=nvcc
CXXFLAGS :=

# === advanced options ===
OUTER_SCALING := 0  # when OUTER_SCALING=1, calculate D1_inv * ( (D1 * A) * (B * D2) ) * D2_inv
                    # where, respectively, D1 and D2 yield rows and cols of A and B with unit norms

INNER_SCALING := 0  # when INNER_SCALING=1, calculate (A * D) * (D_inv * B)
                    # where D yields cols and rows of A and B respectively with roughly equal norms
