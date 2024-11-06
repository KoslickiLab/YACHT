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