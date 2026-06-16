SHELL := /bin/bash

BIN_DIR     := bin
SRC_DIR     := src
INCLUDE_DIR := include

include $(DGEMM_CONFIG_MK)
include config_files/test_configs/common.mk

CONFIG_VARS := DGEMM_IMPL_NAME SO_PATH LDFLAGS RUN_PREFIX CXX USE_CUDA CXXFLAGS \
               OUTER_SCALING INNER_SCALING LOG2_MIN_DIM LOG2_MAX_DIM \
               LOG2_MIN_SCALING_FACTOR LOG2_MAX_SCALING_FACTOR LOG2_INC_SCALING_FACTOR \
               REFERENCE_GEMM DUMP_A_B VERBOSE
$(foreach v,$(CONFIG_VARS),$(eval $(v) := $$(strip $$($(v)))))

# test4 permute: default INNER_PERMUTE=10; override on the command line with
# INNER_PERMUTE=k or OUTER_PERMUTE=k (OUTER_PERMUTE takes precedence when non-zero)
OUTER_PERMUTE ?= 0
INNER_PERMUTE ?= 3

# --- common CPPFLAGS (all test binaries) ---
COMMON_CPPFLAGS :=
COMMON_CPPFLAGS += -DDGEMM_IMPL_NAME=\"$(DGEMM_IMPL_NAME)\"
COMMON_CPPFLAGS += -DLOG2_MIN_DIM=$(LOG2_MIN_DIM)
COMMON_CPPFLAGS += -DLOG2_MAX_DIM=$(LOG2_MAX_DIM)
COMMON_CPPFLAGS += -DLOG2_MIN_SCALING_FACTOR=$(LOG2_MIN_SCALING_FACTOR)
COMMON_CPPFLAGS += -DLOG2_MAX_SCALING_FACTOR=$(LOG2_MAX_SCALING_FACTOR)
COMMON_CPPFLAGS += -DLOG2_INC_SCALING_FACTOR=$(LOG2_INC_SCALING_FACTOR)
COMMON_CPPFLAGS += -DREFERENCE_GEMM=$(REFERENCE_GEMM)
ifeq ($(USE_CUDA),1)
COMMON_CPPFLAGS += -DUSE_CUDA
endif
ifeq ($(VERBOSE),1)
COMMON_CPPFLAGS += -DVERBOSE
endif
ifeq ($(OUTER_SCALING),1)
COMMON_CPPFLAGS += -DOUTER_SCALING
endif
ifeq ($(INNER_SCALING),1)
COMMON_CPPFLAGS += -DINNER_SCALING
endif
ifeq ($(DUMP_A_B),1)
COMMON_CPPFLAGS += -DDUMP_A_B
endif

# --- test4 CPPFLAGS: INPUT_TYPE=test2c + one permute mode ---
TEST4_CPPFLAGS := $(COMMON_CPPFLAGS) -DINPUT_TYPE=\"test2c\"
ifneq ($(OUTER_PERMUTE),0)
TEST4_CPPFLAGS += -DOUTER_PERMUTE=$(OUTER_PERMUTE)
else
TEST4_CPPFLAGS += -DINNER_PERMUTE=$(INNER_PERMUTE)
endif

DGEMM_WRAPPER := $(SRC_DIR)/dgemm_wrapper.cpp

.PHONY: run run_test1 run_test2 run_test3 run_test4 $(BIN_DIR)/test1 $(BIN_DIR)/test2 $(BIN_DIR)/test3 $(BIN_DIR)/test4 compute_distribution_stats stats

# --- test binaries ---
$(BIN_DIR)/test1: $(SRC_DIR)/test1.cpp $(DGEMM_WRAPPER)
	@mkdir -p $(BIN_DIR)
	$(CXX) -std=c++17 $(CXXFLAGS) $(COMMON_CPPFLAGS) -DINPUT_TYPE=\"test1c\" -I$(INCLUDE_DIR) $^ $(LDFLAGS) -o $@

$(BIN_DIR)/test2: $(SRC_DIR)/test2.cpp $(DGEMM_WRAPPER)
	@mkdir -p $(BIN_DIR)
	$(CXX) -std=c++17 $(CXXFLAGS) $(COMMON_CPPFLAGS) -DINPUT_TYPE=\"test2c\" -I$(INCLUDE_DIR) $^ $(LDFLAGS) -o $@

$(BIN_DIR)/test3: $(SRC_DIR)/test3.cpp $(DGEMM_WRAPPER)
	@mkdir -p $(BIN_DIR)
	$(CXX) -std=c++17 $(CXXFLAGS) $(COMMON_CPPFLAGS) -DINPUT_TYPE=\"test3\" -I$(INCLUDE_DIR) $^ $(LDFLAGS) -o $@

$(BIN_DIR)/test4: $(SRC_DIR)/test4.cpp $(DGEMM_WRAPPER)
	@mkdir -p $(BIN_DIR)
	$(CXX) -std=c++17 $(CXXFLAGS) $(TEST4_CPPFLAGS) -I$(INCLUDE_DIR) $^ $(LDFLAGS) -o $@

# --- run targets ---
run_test1: $(BIN_DIR)/test1
	@mkdir -p results/raw_data
	$(RUN_PREFIX) ./$(BIN_DIR)/test1

run_test2: $(BIN_DIR)/test2
	@mkdir -p results/raw_data
	$(RUN_PREFIX) ./$(BIN_DIR)/test2

run_test3: $(BIN_DIR)/test3
	@mkdir -p results/raw_data
	$(RUN_PREFIX) ./$(BIN_DIR)/test3

run_test4: $(BIN_DIR)/test4
	@mkdir -p results/raw_data
	$(RUN_PREFIX) ./$(BIN_DIR)/test4

run: run_test1 run_test2 run_test3 run_test4

# --- standalone stats utility ---
$(BIN_DIR)/compute_distribution_stats: $(SRC_DIR)/compute_distribution_stats.cpp
	@mkdir -p $(BIN_DIR)
	$(CXX) -I$(INCLUDE_DIR) $^ -o $@

compute_distribution_stats stats: $(BIN_DIR)/compute_distribution_stats
