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

# target executable
TARGET1 = $(BIN_DIR)/run_yacht_train_core
TARGET2 = $(BIN_DIR)/run_compute_similarity

# Build rules
all: $(TARGET1) $(TARGET2)

$(TARGET1): $(OBJ_FILES)
	$(CXX) $(CXXFLAGS) $(SRC_DIR)/yacht_train_core.cpp $(SRC_DIR)/utils.cpp -o $(TARGET1)

$(TARGET2): $(OBJ_FILES)
	$(CXX) $(CXXFLAGS) $(SRC_DIR)/compute_similarity.cpp $(SRC_DIR)/utils.cpp -o $(TARGET2)

%.o: %.cpp

clean:
	rm -f $(OBJ_FILES) $(TARGET1) $(TARGET2)

.PHONY: all clean

