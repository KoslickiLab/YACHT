# Compiler and flags
CXX ?= g++
CXXFLAGS = -std=c++17 -Wall -w -O3 -Wsign-compare

# Directories
SRC_DIR = src/cpp
BIN_DIR = src/yacht

# Source files
SRC_FILES = $(SRC_DIR)/yacht_train_core.cpp $(SRC_DIR)/utils.cpp $(SRC_DIR)/yacht_run_compute_similarity.cpp

# Object files
OBJ_FILES = $(SRC_FILES:.cpp=.o)

# Target executables
TARGET1 = $(BIN_DIR)/yacht_train_core
TARGET2 = $(BIN_DIR)/yacht_run_compute_similarity

# Build rules
all: $(TARGET1) $(TARGET2)

$(TARGET1): $(SRC_DIR)/yacht_train_core.cpp $(SRC_DIR)/utils.cpp
	$(CXX) $(CXXFLAGS) -pthread $(SRC_DIR)/yacht_train_core.cpp $(SRC_DIR)/utils.cpp -o $(TARGET1)

$(TARGET2): $(SRC_DIR)/yacht_run_compute_similarity.cpp $(SRC_DIR)/utils.cpp
	$(CXX) $(CXXFLAGS) -pthread $(SRC_DIR)/yacht_run_compute_similarity.cpp $(SRC_DIR)/utils.cpp -o $(TARGET2)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -pthread -c $< -o $@

clean:
	rm -f $(OBJ_FILES) $(TARGET1) $(TARGET2)

.PHONY: all clean

