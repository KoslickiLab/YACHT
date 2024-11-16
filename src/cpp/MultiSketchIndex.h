#ifndef SKETCH_H
#define SKETCH_H

#include <iostream>
#include <vector>
#include <unordered_map>


#ifndef HASH_T
#define HASH_T
typedef unsigned long long int hash_t;
#endif


/**
 * @brief MultiSketchIndex class, which is used to store the index of many sketches.
 * 
 */
class MultiSketchIndex {
    public:
        MultiSketchIndex();
        ~MultiSketchIndex();

        /**
         * @brief Add a hash value to the index.
         * 
         * @param hash_value The hash value to add.
         * @param sketch_index The index of the sketch in which this hash value appears.
         */
        void add_hash(hash_t hash_value, int sketch_index);
        


        /**
         * @brief Add a hash value to the index.
         * 
         * @param hash_value The hash value to add.
         * @param sketch_indices Indices of the sketches in which this hash value appears.
         */
        void add_hash(hash_t hash_value, std::vector<int> sketch_indices);



        /**
         * @brief Get the sketch indices for a hash value.
         * 
         * @param hash_value The hash value to get the sketch indices for.
         * @return const std::vector<int>& The sketch indices in which the hash value appears.
         */
        const std::vector<int>& get_sketch_indices(hash_t hash_value);


        /**
         * @brief Check if a hash value exists in the index.
         * 
         * @param hash_value The hash value to check.
         * @return true If the hash value exists in the index.
         * @return false If the hash value does not exist in the index.
         */
        bool hash_exists(hash_t hash_value) {
            int idx_of_hash = index_of_hash(hash_value);
            return multiple_sketch_indices[idx_of_hash].find(hash_value) != multiple_sketch_indices[idx_of_hash].end();
        }

        
    private:
        std::vector<std::unordered_map<hash_t, std::vector<int>>> multiple_sketch_indices;
        std::vector<std::mutex>mutexes;
        int num_of_indices;

        int index_of_hash(hash_t hash_value) {
            return hash_value % num_of_indices;
        }

        const std::vector<int> empty_vector;
        
};

#endif