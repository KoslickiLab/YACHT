#include "MultiSketchIndex.h"

MultiSketchIndex::MultiSketchIndex() {
    // Constructor
}

MultiSketchIndex::~MultiSketchIndex() {
    // Destructor
}


void MultiSketchIndex::add_hash(hash_t hash_value, std::vector<int> sketch_indices) {
    // Add the hash value to the index
    if (multi_sketch_index.find(hash_value) == multi_sketch_index.end()) {
        multi_sketch_index[hash_value] = sketch_indices;
        return;
    }

    for (int i = 0; i < sketch_indices.size(); i++) {
        add_hash(hash_value, sketch_indices[i]);
    }
}


void MultiSketchIndex::add_hash(hash_t hash_value, int sketch_index) {
    // Add the hash value to the index
    if (multi_sketch_index.find(hash_value) == multi_sketch_index.end()) {
        multi_sketch_index[hash_value] = std::vector<int>();
    }
    multi_sketch_index[hash_value].push_back(sketch_index);
}






const std::vector<int>& MultiSketchIndex::get_sketch_indices(hash_t hash_value) {
    // Get the sketch indices for the hash value
    if (multi_sketch_index.find(hash_value) == multi_sketch_index.end()) {
        return std::vector<int>();
    }
    return multi_sketch_index[hash_value];
}



