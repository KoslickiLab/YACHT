#include "MultiSketchIndex.h"

MultiSketchIndex::MultiSketchIndex() {
    // Constructor
    num_of_indices = 4096;
    multiple_sketch_indices = std::vector<std::unordered_map<hash_t, std::vector<int>>>(num_of_indices);
    mutexes = std::vector<std::mutex>(num_of_indices);
}

MultiSketchIndex::~MultiSketchIndex() {
    // Destructor
}


void MultiSketchIndex::add_hash(hash_t hash_value, std::vector<int> sketch_indices) {
    // Add the hash value to the index
    int idx_of_hash = index_of_hash(hash_value);
    mutexes[idx_of_hash].lock();
    multiple_sketch_indices[idx_of_hash][hash_value] = sketch_indices;
    mutexes[idx_of_hash].unlock();
}


void MultiSketchIndex::add_hash(hash_t hash_value, int sketch_index) {
    // Add the hash value to the index
    int idx_of_hash = index_of_hash(hash_value);
    mutexes[idx_of_hash].lock();
    if (!hash_exists(hash_value)) {
        multiple_sketch_indices[idx_of_hash][hash_value] = std::vector<int>();
    }
    multiple_sketch_indices[idx_of_hash][hash_value].push_back(sketch_index);
    mutexes[idx_of_hash].unlock();
}






const std::vector<int>& MultiSketchIndex::get_sketch_indices(hash_t hash_value) {
    // Get the sketch indices for the hash value
    int idx_of_hash = index_of_hash(hash_value);
    if (hash_exists(hash_value)) {
        return multiple_sketch_indices[idx_of_hash][hash_value];
    } else {
        return empty_vector;
    }
}



