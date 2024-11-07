#include "utils.h"

std::vector<hash_t> read_min_hashes(const std::string& json_filename) {

    // Open the JSON file
    std::ifstream inputFile(json_filename);

    // Check if the file is open
    if (!inputFile.is_open()) {
        std::cerr << "Could not open the file!" << std::endl;
        return {};
    }

    // Parse the JSON data
    json jsonData;
    inputFile >> jsonData;

    // Access and print values
    std::vector<hash_t> min_hashes = jsonData[0]["signatures"][0]["mins"];

    // Close the file
    inputFile.close();

    return min_hashes;
}


void compute_index_from_sketches(std::vector<std::vector<hash_t>>& sketches, std::unordered_map<hash_t, std::vector<int>>& hash_index) {
    // create the index using all the hashes
    for (uint i = 0; i < sketches.size(); i++) {
        for (uint j = 0; j < sketches[i].size(); j++) {
            hash_t hash = sketches[i][j];
            if (hash_index.find(hash) == hash_index.end()) {
                hash_index[hash] = std::vector<int>();
            }
            hash_index[hash].push_back(i);
        }
    }

    size_t num_hashes = hash_index.size();

    // remove the hashes that only appear in one sketch
    std::vector<hash_t> hashes_to_remove;
    for (auto it = hash_index.begin(); it != hash_index.end(); it++) {
        if (it->second.size() == 1) {
            hashes_to_remove.push_back(it->first);
        }
    }
    for (uint i = 0; i < hashes_to_remove.size(); i++) {
        hash_index.erase(hashes_to_remove[i]);
    }

    size_t num_hashes_after_removal = hash_index.size();

    std::cout << "Total number of distinct hashes: " << num_hashes << std::endl;
    std::cout << "Total number of distinct hashes that appear in only one sketch: " << num_hashes - num_hashes_after_removal << std::endl;
    std::cout << "Size of the index: " << num_hashes_after_removal << std::endl;

}



void get_sketch_names(const std::string& filelist, std::vector<std::string>& sketch_names) {
    // the filelist is a file, where each line is a path to a sketch file
    std::ifstream file(filelist);
    if (!file.is_open()) {
        std::cerr << "Could not open the filelist: " << filelist << std::endl;
        return;
    }
    std::string line;
    while (std::getline(file, line)) {
        sketch_names.push_back(line);
    }
}



void read_sketches_one_chunk(int start_index, int end_index, 
                            std::vector<std::string>& sketch_names,
                            std::vector<std::vector<hash_t>>& sketches,
                            std::mutex& mutex_count_empty_sketch,
                            std::vector<int>& empty_sketch_ids) {

    for (int i = start_index; i < end_index; i++) {
        auto min_hashes_genome_name = read_min_hashes(sketch_names[i]);
        sketches[i] = min_hashes_genome_name;
        if (sketches[i].size() == 0) {
            mutex_count_empty_sketch.lock();
            empty_sketch_ids.push_back(i);
            mutex_count_empty_sketch.unlock();
        }
    }
}



void read_sketches(std::vector<std::string>& sketch_names,
                        std::vector<std::vector<hash_t>>& sketches, 
                        std::vector<int>& empty_sketch_ids, 
                        const uint num_threads) {

    uint num_sketches = sketch_names.size();
    for (uint i = 0; i < num_sketches; i++) {
        sketches.push_back( std::vector<hash_t>() );
    }

    std::mutex mutex_count_empty_sketch;
    int chunk_size = num_sketches / num_threads;
    std::vector<std::thread> threads;
    for (int i = 0; i < num_threads; i++) {
        int start_index = i * chunk_size;
        int end_index = (i == num_threads - 1) ? num_sketches : (i + 1) * chunk_size;
        threads.push_back(std::thread(read_sketches_one_chunk, 
                                        start_index, end_index, 
                                        std::ref(sketch_names), std::ref(sketches), 
                                        std::ref(mutex_count_empty_sketch), 
                                        std::ref(empty_sketch_ids)));
    }
    for (int i = 0; i < num_threads; i++) {
        threads[i].join();
    }
    
}




void show_empty_sketches(const std::vector<int>& empty_sketch_ids) {
    int count_empty_sketch = empty_sketch_ids.size();
    std::cout << "Number of empty sketches: " << count_empty_sketch << std::endl;
    if (count_empty_sketch == 0) {
        return;
    }
    std::cout << "Empty sketch ids: ";
    for (int i : empty_sketch_ids) {
        std::cout << i << " ";
    }
    std::cout << std::endl;
}





void compute_intersection_matrix_by_sketches(int query_sketch_start_index, int query_sketch_end_index, 
                                            int thread_id, std::string out_dir, 
                                            int pass_id, int negative_offset,
                                            const std::vector<std::vector<hash_t>>& sketches_query,
                                            const std::vector<std::vector<hash_t>>& sketches_ref,
                                            const std::unordered_map<hash_t, std::vector<int>>& hash_index_ref,
                                            int** intersectionMatrix, 
                                            double containment_threshold,
                                            std::vector<std::vector<int>>& similars) {
    
    const int num_sketches_ref = sketches_ref.size();
    const int num_sketches_query = sketches_query.size();

    // process the sketches in the range [sketch_start_index, sketch_end_index)
    for (uint i = query_sketch_start_index; i < query_sketch_end_index; i++) {
        for (int j = 0; j < sketches_query[i].size(); j++) {
            hash_t hash = sketches_query[i][j];
            if (hash_index_ref.find(hash) != hash_index_ref.end()) {
                std::vector<int> ref_sketch_indices = hash_index_ref.at(hash);
                for (uint k = 0; k < ref_sketch_indices.size(); k++) {
                    intersectionMatrix[i-negative_offset][ref_sketch_indices[k]]++;
                }
            }
        }
    }

    // write the similarity values to file. filename: out_dir/passid_threadid.txt, where id is thread id in 3 digits
    std::string id_in_three_digits_str = std::to_string(thread_id);
    while (id_in_three_digits_str.size() < 3) {
        id_in_three_digits_str = "0" + id_in_three_digits_str;
    }
    std::string pass_id_str = std::to_string(pass_id);
    std::string filename = out_dir + "/" + pass_id_str + "_" + id_in_three_digits_str + ".txt";
    std::ofstream outfile(filename);

    // only write the values if larger than the threshold
    for (int i = query_sketch_start_index; i < query_sketch_end_index; i++) {
        for (uint j = 0; j < num_sketches_ref; j++) {
            // if nothing in the intersection, then skip
            if (intersectionMatrix[i-negative_offset][j] == 0) {
                continue;
            }

            // if either of the sketches is empty, then skip
            if (sketches_query[i].size() == 0 || sketches_ref[j].size() == 0) {
                continue;
            }

            // if the divisor in the jaccard calculation is 0, then skip
            if (sketches_query[i].size() + sketches_ref[j].size() - intersectionMatrix[i-negative_offset][j] == 0) {
                continue;
            }

            double jaccard = 1.0 * intersectionMatrix[i-negative_offset][j] / ( sketches_query[i].size() + sketches_ref[j].size() - intersectionMatrix[i-negative_offset][j] );
            double containment_i_in_j = 1.0 * intersectionMatrix[i-negative_offset][j] / sketches_query[i].size();
            double containment_j_in_i = 1.0 * intersectionMatrix[i-negative_offset][j] / sketches_ref[j].size();
            
            // containment_i_in_j is the containment of query in target, i is the query
            if (containment_i_in_j < containment_threshold) {
                continue;
            }

            outfile << i << "," << j << "," << jaccard << "," << containment_i_in_j << "," << containment_j_in_i << std::endl;
            similars[i].push_back(j);
        }
    }

    outfile.close();

}



void compute_intersection_matrix(const std::vector<std::vector<hash_t>>& sketches_query,
                                const std::vector<std::vector<hash_t>>& sketches_ref, 
                                const std::unordered_map<hash_t, std::vector<int>>& hash_index_ref,
                                const std::string& out_dir, 
                                std::vector<std::vector<int>>& similars,
                                double containment_threshold,
                                const int num_passes, const int num_threads) {
    
    int num_sketches_query = sketches_query.size();
    int num_sketches_ref = sketches_ref.size();

    // allocate memory for the intersection matrix
    int num_query_sketches_each_pass = ceil(1.0 * num_sketches_query / num_passes);
    int** intersectionMatrix = new int*[num_query_sketches_each_pass + 1];
    for (int i = 0; i < num_query_sketches_each_pass + 1; i++) {
        intersectionMatrix[i] = new int[num_sketches_ref];
    }

    // allocate memory for the similars
    for (int i = 0; i < num_sketches_query; i++) {
        similars.push_back(std::vector<int>());
    }

    for (int pass_id = 0; pass_id < num_passes; pass_id++) {
        // set zeros in the intersection matrix
        for (int i = 0; i < num_query_sketches_each_pass+1; i++) {
            for (uint j = 0; j < num_sketches_ref; j++) {
                intersectionMatrix[i][j] = 0;
            }
        }

        // prepare the indices which will be processed in this pass
        int sketch_idx_start_this_pass = pass_id * num_query_sketches_each_pass;
        int sketch_idx_end_this_pass = (pass_id == num_passes - 1) ? num_sketches_query : (pass_id + 1) * num_query_sketches_each_pass;
        int negative_offset = pass_id * num_query_sketches_each_pass;
        int num_query_sketches_this_pass = sketch_idx_end_this_pass - sketch_idx_start_this_pass;
        
        // create threads
        std::vector<std::thread> threads;
        int chunk_size = num_query_sketches_this_pass / num_threads;
        for (int i = 0; i < num_threads; i++) {
            int start_query_index_this_thread = sketch_idx_start_this_pass + i * chunk_size;
            int end_query_index_this_thread = (i == num_threads - 1) ? sketch_idx_end_this_pass : sketch_idx_start_this_pass + (i + 1) * chunk_size;
            threads.push_back(std::thread(compute_intersection_matrix_by_sketches, 
                            start_query_index_this_thread, end_query_index_this_thread, 
                            i, out_dir, pass_id, negative_offset,
                            std::ref(sketches_query), std::ref(sketches_ref), 
                            std::ref(hash_index_ref), intersectionMatrix, 
                            containment_threshold,
                            std::ref(similars)));
        }

        // join threads
        for (int i = 0; i < num_threads; i++) {
            threads[i].join();
        }

        // show progress
        std::cout << "Pass " << pass_id+1 << "/" << num_passes << " done." << std::endl;
    }

    // free the memory allocated for the intersection matrix
    for (int i = 0; i < num_query_sketches_each_pass + 1; i++) {
        delete[] intersectionMatrix[i];
    }
    delete[] intersectionMatrix;
}