#############################
# User Options
#############################

DGEMM_IMPLEMENTATIONS := REFERENCE_1 NETLIBBLAS BLIS OPENBLAS STRASSEN_1LVL CUBLAS OZIMMU
DGEMM_IMPL ?= REFERENCE_1

# Directories
SRC_DIR := src
INC_DIR := include
LIB_DIR := lib

# Flags
CXXFLAGS := -g -O0 -I$(INC_DIR) -I/opt/cray/pe/hdf5/1.14.3.1/gnu/12.3/include
LDFLAGS := -L/opt/cray/pe/hdf5/1.14.3.1/gnu/12.3/lib -lhdf5_cpp -lhdf5

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

# Determinisim
ifeq ($(DETERMINISTIC),1)
	CXXFLAGS += -DDETERMINISTIC
endif

# Outputs
MAIN_EXE     := main
STRASSEN_1LVL_SO  := libstrassen_1lvl_dgemm.so
REFERENCE_1_SO := libreference_1_dgemm.so

########################################
ifeq ($(DGEMM_IMPL),REFERENCE_1)
	LDFLAGS +=  -L$(LIB_DIR) -lreference_1_dgemm
	CXXFLAGS += -DDGEMM_IMPL=\"REFERENCE_1\"
    RUN_PREFIX := LD_LIBRARY_PATH=$(LD_LIBRARY_PATH):$(LIB_DIR)
endif

ifeq ($(DGEMM_IMPL),NETLIBBLAS)
    LDFLAGS += -lrefblas
	CXXFLAGS += -DDGEMM_IMPL=\"NETLIBBLAS\"
    RUN_PREFIX :=
endif

ifeq ($(DGEMM_IMPL),BLIS)
    LDFLAGS += -lblis
	CXXFLAGS += -DDGEMM_IMPL=\"BLIS\"
    RUN_PREFIX :=
endif

ifeq ($(DGEMM_IMPL),OPENBLAS)
    LDFLAGS += -lopenblas
	CXXFLAGS += -DDGEMM_IMPL=\"OPENBLAS\"
    RUN_PREFIX :=
endif

ifeq ($(DGEMM_IMPL),STRASSEN_1LVL)
    LDFLAGS += -L$(LIB_DIR) -lstrassen_1lvl_dgemm
	CXXFLAGS += -DDGEMM_IMPL=\"STRASSEN_1LVL\"
    RUN_PREFIX := LD_LIBRARY_PATH=$(LD_LIBRARY_PATH):$(LIB_DIR)
endif

ifeq ($(DGEMM_IMPL),CUBLAS)
    LDFLAGS += -lcublas
	CXXFLAGS += -DDGEMM_IMPL=\"CUBLAS\" -DCUDA
    RUN_PREFIX :=
endif

ifeq ($(DGEMM_IMPL),OZIMMU)
    LDFLAGS += -lcublas
	CXXFLAGS += -DDGEMM_IMPL=\"OZIMMU\" -DCUDA
    RUN_PREFIX := LD_PRELOAD=$(OZIMMU_SO_PATH)
endif
########################################

all: dirs $(LIB_DIR)/$(STRASSEN_1LVL_SO) $(LIB_DIR)/$(REFERENCE_1_SO) $(MAIN_EXE)

$(LIB_DIR)/$(STRASSEN_1LVL_SO): $(SRC_DIR)/strassen_1lvl_dgemm.cpp
	g++ -O0 -fPIC -shared -Iinclude $< -o $@

$(LIB_DIR)/$(REFERENCE_1_SO): $(SRC_DIR)/reference_1_dgemm.cpp
	g++ -O0 -fPIC -shared -Iinclude $< -o $@

$(MAIN_EXE): $(SRC_DIR)/main.cpp
	$(CXX) $(CXXFLAGS) $< $(LDFLAGS) -o $@

dirs:
	@mkdir -p $(LIB_DIR)

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
