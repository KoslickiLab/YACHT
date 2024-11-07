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

std::vector<hash_t> read_min_hashes(const std::string&);

void compute_index_from_sketches(std::vector<std::vector<hash_t>>&, std::unordered_map<hash_t, std::vector<int>>&);

void get_sketch_names(const std::string&, std::vector<std::string>&);

void read_sketches(std::vector<std::string>& sketch_names,
                        std::vector<std::vector<hash_t>>& sketches, 
                        std::vector<int>& empty_sketch_ids, 
                        const uint num_threads);

void show_empty_sketches(const std::vector<int>&);


void compute_intersection_matrix(const std::vector<std::vector<hash_t>>& sketches_query,
                                const std::vector<std::vector<hash_t>>& sketches_ref, 
                                const std::unordered_map<hash_t, std::vector<int>>& hash_index_ref,
                                const std::string& out_dir, 
                                std::vector<std::vector<int>>& similars,
                                double containment_threshold,
                                const int num_passes, const int num_threads);

#endif