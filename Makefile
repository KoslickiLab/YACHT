# Compiler and flags
CXX ?= g++
CXXFLAGS = -std=c++17 -Wall -w -O3 -Wsign-compare

# Directories
SRC_DIR = src/cpp
BIN_DIR = src/yacht

# Source files
SRC_FILES = $(SRC_DIR)/yacht_train_core.cpp $(SRC_DIR)/utils.cpp $(SRC_DIR)/compute_similarity.cpp

# Object files
OBJ_FILES = $(SRC_FILES:.cpp=.o)

# Target executable
TARGET = $(BIN_DIR)/run_yacht_train_core $(BIN_DIR)/run_compute_similarity

# Default target
all: $(TARGET)

# Create the bin directory if it doesn't exist
$(BIN_DIR):
	echo "Creating directory: $(BIN_DIR)"
	mkdir -p $(BIN_DIR)

# build the object files
$(OBJ_FILES): %.o: %.cpp
	echo "Compiling: $<"
	$(CXX) $(CXXFLAGS) -c $< -o $@

# build the target run_yacht_train_core
$(BIN_DIR)/run_yacht_train_core: $(OBJ_FILES)
	echo "Building target: $@"
	$(CXX) $(CXXFLAGS) $(SRC_DIR)/yacht_train_core.cpp $(SRC_DIR)/utils.o -o $@

# build the target run_compute_similarity
$(BIN_DIR)/run_compute_similarity: $(OBJ_FILES)
	echo "Building target: $@"
	$(CXX) $(CXXFLAGS) $(SRC_DIR)/compute_similarity.cpp $(SRC_DIR)/utils.o -o $@

# clean up
clean:
	rm -f $(OBJ_FILES) $(TARGET)