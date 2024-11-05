# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++17 -Wall -w -O3 -Wsign-compare

# Directories
SRC_DIR = src/cpp
BIN_DIR = src/yacht

# Source files
SRC_FILES = $(SRC_DIR)/main.cpp

# Object files
OBJ_FILES = $(SRC_FILES:.cpp=.o)

# Target executable
TARGET = $(BIN_DIR)/run_yacht_train_core

# Default target
all: $(TARGET)

# Create the bin directory if it doesn't exist
$(BIN_DIR):
	echo "Creating directory: $(BIN_DIR)"
	mkdir -p $(BIN_DIR)

# build the object files
$(OBJ_FILES): %.o: %.cpp
	echo "Compiling: $<"
	$(CXX) $(CXXFLAGS) -c $< -o $@ -lz

# build the target executable
$(TARGET): $(OBJ_FILES) | $(BIN_DIR)
	echo "Linking to create executable: $(TARGET)"
	$(CXX) $(CXXFLAGS) $(OBJ_FILES) -o $(TARGET) -lz -lpthread

# clean up
clean:
	rm -f $(OBJ_FILES) $(TARGET)