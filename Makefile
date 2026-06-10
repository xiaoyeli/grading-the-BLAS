SHELL := /bin/bash

RUN_ID := $(shell date +%Y%m%d_%H%M%S)

MAIN_EXE := run_test
STATS_EXE := compute_distribution_stats
BIN_DIR := bin
SRC_DIR := src
INCLUDE_DIR := include

include $(DGEMM_CONFIG_MK)

include $(TEST_CONFIG_MK)

# strip whitespace from config vars
CONFIG_VARS := DGEMM_IMPL_NAME SO_PATH LDFLAGS RUN_PREFIX CXX USE_CUDA CXXFLAGS \
				OUTER_SCALING INNER_SCALING INPUT_TYPE LOG2_MIN_DIM LOG2_MAX_DIM \
				LOG2_MIN_SCALING_FACTOR LOG2_MAX_SCALING_FACTOR LOG2_INC_SCALING_FACTOR \
				REFERENCE_GEMM DUMP_A_B DUMP_C VERBOSE OUTER_PERMUTE INNER_PERMUTE
$(foreach v,$(CONFIG_VARS),$(eval $(v) := $$(strip $$($(v)))))

# values specified in config.sh
CPPFLAGS += -DDGEMM_IMPL_NAME=\"$(DGEMM_IMPL_NAME)\"
CPPFLAGS += -DLOG2_MIN_DIM=$(LOG2_MIN_DIM)
CPPFLAGS += -DLOG2_MAX_DIM=$(LOG2_MAX_DIM)
CPPFLAGS += -DLOG2_MIN_SCALING_FACTOR=$(LOG2_MIN_SCALING_FACTOR)
CPPFLAGS += -DLOG2_MAX_SCALING_FACTOR=$(LOG2_MAX_SCALING_FACTOR)
CPPFLAGS += -DLOG2_INC_SCALING_FACTOR=$(LOG2_INC_SCALING_FACTOR)
CPPFLAGS += -DINPUT_TYPE=\"$(INPUT_TYPE)\"
CPPFLAGS += -DREFERENCE_GEMM=$(REFERENCE_GEMM)

ifeq ($(USE_CUDA),1)
	CPPFLAGS += -DUSE_CUDA
endif

ifeq ($(VERBOSE),1)
	CPPFLAGS += -DVERBOSE
endif

ifneq ($(OUTER_PERMUTE),0)
	CPPFLAGS += -DOUTER_PERMUTE=$(OUTER_PERMUTE)
endif

ifneq ($(INNER_PERMUTE),0)
	CPPFLAGS += -DINNER_PERMUTE=$(INNER_PERMUTE)
endif

ifeq ($(OUTER_SCALING),1)
	CPPFLAGS += -DOUTER_SCALING
endif

ifeq ($(INNER_SCALING),1)
	CPPFLAGS += -DINNER_SCALING
endif

ifeq ($(DUMP_A_B),1)
	CPPFLAGS += -DDUMP_A_B
endif

ifeq ($(DUMP_C),1)
	CPPFLAGS += -DDUMP_C
endif

MAIN_SRCS := $(SRC_DIR)/main.cpp $(SRC_DIR)/dgemm_wrapper.cpp
STATS_SRCS := $(SRC_DIR)/compute_distribution_stats.cpp

.PHONY: run

$(MAIN_EXE): $(MAIN_SRCS)
	@mkdir -p $(BIN_DIR)
	$(CXX) -std=c++17 $(CXXFLAGS) $(CPPFLAGS) -I$(INCLUDE_DIR) $^ $(LDFLAGS) -o ./$(BIN_DIR)/$@

$(STATS_EXE): $(STATS_SRCS)
	@mkdir -p $(BIN_DIR)
	$(CXX) -I$(INCLUDE_DIR) $^ -o ./$(BIN_DIR)/$@

run: $(MAIN_EXE)
	@mkdir -p results/$(RUN_ID)
	@cp $(DGEMM_CONFIG_MK) results/$(RUN_ID)/
	@cp $(TEST_CONFIG_MK) results/$(RUN_ID)/
	RUN_OUTPUT_DIR="results/$(RUN_ID)/raw_data" $(RUN_PREFIX) ./$(BIN_DIR)/$(MAIN_EXE)
	.venv/bin/python3 scripts/plot_sweep_results.py results/$(RUN_ID)

stats: $(STATS_EXE)