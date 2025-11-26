CXX := g++
CXXFLAGS := -O0 -std=c++17 -Wall -Wextra -Wno-unused-parameter -Iinclude -Isubmodules/blis-2.0/include/sandybridge -MMD -MP
LDFLAGS := -lquadmath -Lsubmodules/blis-2.0/lib/sandybridge/ -lblis
CPPFLAGS ?=

ifeq ($(SANITY),1)
CPPFLAGS += -DSANITY
endif

# directories
SRC_DIR := src
OBJ_DIR := obj
RESULTS_DIRS := results/U_0_1 results/U_-1_1 results/test1

BIN := main

# Source and object files
SRCS := $(wildcard $(SRC_DIR)/*.cpp)
OBJS := $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(SRCS))
DEPS := $(OBJS:.o=.d)

# Default target
all: $(BIN) $(RESULTS_DIRS)

# Linking
$(BIN): $(OBJS)
	$(CXX) $(OBJS) -o $@ $(LDFLAGS)
	mkdir -p $(RESULTS_DIRS)

# Compilation
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

# Ensure obj/ exists
$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

# Ensure results/ and subdirectories exist
$(RESULTS_DIRS):
	mkdir -p $@

# Clean up
clean:
	rm -rf $(OBJ_DIR) $(BIN)

.PHONY: all clean

# Include dependency files so .hpp changes trigger recompilation
-include $(DEPS)