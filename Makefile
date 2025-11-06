CXX := g++
CXXFLAGS := -O0 -std=c++17 -Wall -Wextra -Iinclude -MMD -MP
LDFLAGS := -lquadmath

# directories
SRC_DIR := src
OBJ_DIR := obj
BIN := main

# Source and object files
SRCS := $(wildcard $(SRC_DIR)/*.cpp)
OBJS := $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(SRCS))
DEPS := $(OBJS:.o=.d)

# Default target
all: $(BIN)

# Linking
$(BIN): $(OBJS)
	$(CXX) $(OBJS) -o $@ $(LDFLAGS)

# Compilation
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Ensure obj/ exists
$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

# Clean up
clean:
	rm -rf $(OBJ_DIR) $(BIN)

.PHONY: all clean

# Include dependency files so .hpp changes trigger recompilation
-include $(DEPS)