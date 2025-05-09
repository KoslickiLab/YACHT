#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <stdexcept>
#include <thread>
#include <mutex>
#include <algorithm>
#include <utility>
#include <set>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <cmath>
#include <limits>
#include <iomanip>
#include <sstream>
#include <chrono>
#include <random>

#include "json.hpp"

using json = nlohmann::json;

typedef unsigned long long int hash_t;

/**
 * @brief Read the min-hashes from a FMH sketch file
 * 
 * Assumption: the file is a json file, and its not gzipped
 * 
 * @param sketch_path The path to the sketch file
 */
std::vector<hash_t> read_min_hashes(const std::string& sketch_path);






/**
 * @brief Compute the index from the sketches
 * 
 * @param sketches The sketches
 * @param hash_index The reference to the hash index (where the index will be stored)
 */
void compute_index_from_sketches(std::vector<std::vector<hash_t>>& sketches, std::unordered_map<hash_t, std::vector<int>>& hash_index);






/**
 * @brief Get the sketch paths 
 * 
 * @param filelist The file containing the paths of the sketches
 * @param sketch_paths The vector to store the paths
 */
void get_sketch_paths(const std::string& filelist, std::vector<std::string>& sketch_paths);






/**
 * @brief Read the sketches from the sketch paths
 * 
 * @param sketch_paths The paths to the sketches
 * @param sketches The vector to store the sketches
 * @param empty_sketch_ids The vector to store the ids of empty sketches
 * @param num_threads The number of threads to use
 */
void read_sketches(std::vector<std::string>& sketch_paths,
                        std::vector<std::vector<hash_t>>& sketches, 
                        std::vector<int>& empty_sketch_ids, 
                        const uint num_threads);





/**
 * @brief Show the empty sketches
 * 
 * @param empty_sketch_ids The ids of the empty sketches
 */
void show_empty_sketches(const std::vector<int>&);




/**
 * @brief Compute the intersection matrix
 * 
 * @param sketches_query The query sketches
 * @param sketches_ref The reference (target) sketches
 * @param hash_index_ref The index of the reference (target) sketches
 * @param out_dir The output directory to store the results
 * @param similars The vector to store the similar sketches
 * @param containment_threshold The containment threshold
 * @param num_passes The number of passes to use
 * @param num_threads The number of threads to use
 */
void compute_intersection_matrix(const std::vector<std::vector<hash_t>>& sketches_query,
                                const std::vector<std::vector<hash_t>>& sketches_ref, 
                                const std::unordered_map<hash_t, std::vector<int>>& hash_index_ref,
                                const std::string& out_dir, 
                                std::vector<std::vector<int>>& similars,
                                double containment_threshold,
                                const int num_passes, const int num_threads);

void compute_intersection_matrix2(const std::vector<std::vector<hash_t>>& sketches_query,
    const std::vector<std::vector<hash_t>>& sketches_ref, 
    const std::unordered_map<hash_t, std::vector<int>>& hash_index_ref,
    const std::string& out_dir, 
    std::vector<std::vector<int>>& similars,
    double containment_threshold,
    const int num_passes, const int num_threads,
    const std::vector<std::string>& sketch_paths_query,
    const std::vector<std::string>& sketch_paths_target,
    const std::unordered_map<std::string, std::string>& sig_to_name);

#endif