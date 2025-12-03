#############################
# User Options
#############################

DGEMM_IMPLEMENTATIONS := REFERENCE NETLIBBLAS BLIS OPENBLAS STRASSEN CUBLAS OZIMMU
DGEMM_IMPL ?= REFERENCE

# Directories
SRC_DIR := src
INC_DIR := include
LIB_DIR := lib
RESULTS_DIRS := results/U_0_1 results/U_-1_1 results/test1 results/test3 results/sanity

# Flags
CXXFLAGS := -O0 -I$(INC_DIR)

# Compiler
ifneq ($(filter $(DGEMM_IMPL), OZIMMU CUBLAS),)
	CXX := nvcc
else
	CXX := g++
endif

# Sanity/Verbosity
ifeq ($(SANITY),1)
	CXXFLAGS += -DSANITY -DVERBOSE
else ifeq ($(VERBOSE),1)
	CXXFLAGS += -DVERBOSE
endif

# Outputs
MAIN_EXE     := main
STRASSEN_SO  := libstrassen_dgemm.so

########################################
ifeq ($(DGEMM_IMPL),REFERENCE)
	LDFLAGS := 
	RUN_PREFIX :=
endif

ifeq ($(DGEMM_IMPL),NETLIBBLAS)
    LDFLAGS += -lrefblas
	CXXFLAGS += -DUSE_NETLIBBLAS
    RUN_PREFIX :=
endif

ifeq ($(DGEMM_IMPL),BLIS)
    LDFLAGS += -lblis
	CXXFLAGS += -DUSE_BLIS
    RUN_PREFIX :=
endif

ifeq ($(DGEMM_IMPL),OPENBLAS)
    LDFLAGS += -lopenblas
	CXXFLAGS += -DUSE_OPENBLAS
    RUN_PREFIX :=
endif

ifeq ($(DGEMM_IMPL),STRASSEN)
    LDFLAGS += -L$(LIB_DIR) -lstrassen_dgemm
	CXXFLAGS += -DUSE_STRASSEN
    RUN_PREFIX := LD_LIBRARY_PATH=$(LD_LIBRARY_PATH):$(LIB_DIR)
endif

ifeq ($(DGEMM_IMPL),CUBLAS)
    LDFLAGS += -lcublas
	CXXFLAGS += -DUSE_CUBLAS -DCUDA
    RUN_PREFIX :=
endif

ifeq ($(DGEMM_IMPL),OZIMMU)
    LDFLAGS := -lcublas
	CXXFLAGS += -DUSE_OZIMMU -DCUDA
    RUN_PREFIX := LD_PRELOAD=$(OZIMMU_SO_PATH)
endif
########################################

all: dirs $(LIB_DIR)/$(STRASSEN_SO) $(MAIN_EXE)

$(LIB_DIR)/$(STRASSEN_SO): $(SRC_DIR)/strassen_dgemm.cpp
	g++ -O0 -fPIC -shared -Iinclude $< -o $@

$(MAIN_EXE): $(SRC_DIR)/main.cpp
	$(CXX) $(CXXFLAGS) $< $(LDFLAGS) -o $@

dirs:
	@mkdir -p $(RESULTS_DIRS) $(LIB_DIR)

run: all
	$(RUN_PREFIX) ./$(MAIN_EXE)

run_all:
	@for implementation in $(DGEMM_IMPLEMENTATIONS); do \
		echo ""; \
		echo "=== Building and running $$implementation ==="; \
		$(MAKE) $(MAKEFLAGS) DGEMM_IMPL=$$implementation run; \
	done

###########################################
# Cleanup
###########################################
clean:
	rm -rf $(MAIN_EXE) $(LIB_DIR)

.PHONY: all dirs clean run run_all_tests $(MAIN_EXE)
