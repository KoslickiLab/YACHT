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



void get_sketch_names(const std::string& filelist, std::vector<std::string>& sketch_names, uint& num_sketches) {
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
    num_sketches = sketch_names.size();
}



void read_sketches_one_chunk(int start_index, int end_index, 
                            std::vector<std::string>& sketch_names,
                            std::vector<std::vector<hash_t>>& sketches,
                            std::vector<std::pair<int, int>>& genome_id_size_pairs,
                            std::mutex& mutex_count_empty_sketch,
                            int& count_empty_sketch,
                            std::vector<int>& empty_sketch_ids) {

    for (int i = start_index; i < end_index; i++) {
        auto min_hashes_genome_name = read_min_hashes(sketch_names[i]);
        sketches[i] = min_hashes_genome_name;
        if (sketches[i].size() == 0) {
            mutex_count_empty_sketch.lock();
            count_empty_sketch++;
            empty_sketch_ids.push_back(i);
            mutex_count_empty_sketch.unlock();
        }
        genome_id_size_pairs[i] = {i, sketches[i].size()};
    }
}



void read_sketches(const uint num_sketches, std::vector<std::vector<hash_t>>& sketches, 
                            std::vector<std::pair<int, int>>& genome_id_size_pairs,
                            const uint num_threads, std::vector<std::string>& sketch_names,
                            int& count_empty_sketch, std::vector<int>& empty_sketch_ids, std::mutex& mutex_count_empty_sketch) {

    for (uint i = 0; i < num_sketches; i++) {
        sketches.push_back( std::vector<hash_t>() );
        genome_id_size_pairs.push_back({-1, 0});
    }

    int chunk_size = num_sketches / num_threads;
    std::vector<std::thread> threads;
    for (int i = 0; i < num_threads; i++) {
        int start_index = i * chunk_size;
        int end_index = (i == num_threads - 1) ? num_sketches : (i + 1) * chunk_size;
        threads.push_back(std::thread(read_sketches_one_chunk, start_index, end_index, 
                            std::ref(sketch_names), std::ref(sketches), 
                            std::ref(genome_id_size_pairs), std::ref(mutex_count_empty_sketch), 
                            std::ref(count_empty_sketch), std::ref(empty_sketch_ids)));
    }
    for (int i = 0; i < num_threads; i++) {
        threads[i].join();
    }
    
}