/*
 * Author: Mahmudur Rahman Hera (mahmudhera93@gmail.com)
 * Date: November 1, 2024
 * Description: yacht train core using indexing of sketches to do genome comparison
 */

#include "argparse.hpp"
#include "json.hpp"

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

using namespace std;
using json = nlohmann::json;


struct Arguments {
    string file_list;
    string working_directory;
    string output_filename;
    int number_of_threads;
    int num_of_passes;
    double containment_threshold;
};


typedef Arguments Arguments;
typedef unsigned long long int hash_t;


Arguments arguments;
std::vector<std::string> sketch_names;
uint num_sketches;
vector<vector<hash_t>> sketches;
vector<pair<int, int>> genome_id_size_pairs;
unordered_map<hash_t, vector<int>> hash_index;
int count_empty_sketch = 0;
mutex mutex_count_empty_sketch;
vector<int> empty_sketch_ids;
int ** intersectionMatrix;
vector<vector<int>> similars;



vector<hash_t> read_min_hashes(const string& json_filename) {

    // Open the JSON file
    ifstream inputFile(json_filename);

    // Check if the file is open
    if (!inputFile.is_open()) {
        cerr << "Could not open the file!" << endl;
        return {};
    }

    // Parse the JSON data
    json jsonData;
    inputFile >> jsonData;

    // Access and print values
    vector<hash_t> min_hashes = jsonData[0]["signatures"][0]["mins"];

    // Close the file
    inputFile.close();

    return min_hashes;
}




void read_sketches_one_chunk(int start_index, int end_index) {
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



void read_sketches() {
    for (uint i = 0; i < num_sketches; i++) {
        sketches.push_back( vector<hash_t>() );
        genome_id_size_pairs.push_back({-1, 0});
    }

    int num_threads = arguments.number_of_threads;

    int chunk_size = num_sketches / num_threads;
    vector<thread> threads;
    for (int i = 0; i < num_threads; i++) {
        int start_index = i * chunk_size;
        int end_index = (i == num_threads - 1) ? num_sketches : (i + 1) * chunk_size;
        threads.push_back(thread(read_sketches_one_chunk, start_index, end_index));
    }
    for (int i = 0; i < num_threads; i++) {
        threads[i].join();
    }
    
}


void get_sketch_names(const std::string& filelist) {
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


void parse_arguments(int argc, char *argv[]) {

    argparse::ArgumentParser parser("yacht train using indexing of sketches");
    
    parser.add_argument("file_list")
        .help("file containing list of files to be processed")
        .store_into(arguments.file_list);
    parser.add_argument("working_directory")
        .help("working directory (where temp files are generated)")
        .store_into(arguments.working_directory);
    parser.add_argument("output_filename")
        .help("output filename (where the reduced ref filenames will be written)")
        .store_into(arguments.output_filename);
    parser.add_argument("-t", "--threads")
        .help("number of threads")
        .scan<'i', int>()
        .default_value(1)
        .store_into(arguments.number_of_threads);
    parser.add_argument("-p", "--passes")
        .help("number of passes")
        .scan<'i', int>()
        .default_value(1)
        .store_into(arguments.num_of_passes);
    parser.add_argument("-c", "--containment_threshold")
        .help("containment threshold")
        .scan<'g', double>()
        .default_value(0.9)
        .store_into(arguments.containment_threshold);
    parser.parse_args(argc, argv);

    if (arguments.number_of_threads < 1) {
        throw std::runtime_error("number of threads must be at least 1");
    }

    if (arguments.num_of_passes < 1) {
        throw std::runtime_error("number of passes must be at least 1");
    }

    if (arguments.containment_threshold < 0.0 || arguments.containment_threshold > 1.0) {
        throw std::runtime_error("containment threshold must be between 0.0 and 1.0");
    }

}


void show_arguments() {
    cout << "Working with the following parameters:" << endl;
    cout << "**************************************" << endl;
    cout << "*" << endl;
    cout << "*    file_list: " << arguments.file_list << endl;
    cout << "*    working_directory: " << arguments.working_directory << endl;
    cout << "*    output_filename: " << arguments.output_filename << endl;
    cout << "*    number_of_threads: " << arguments.number_of_threads << endl;
    cout << "*    num_of_passes: " << arguments.num_of_passes << endl;
    cout << "*    containment_threshold: " << arguments.containment_threshold << endl;
    cout << "*" << endl;
    cout << "**************************************" << endl;
}


void show_empty_sketches() {
    cout << "Number of empty sketches: " << count_empty_sketch << endl;
    if (count_empty_sketch == 0) {
        return;
    }
    cout << "Empty sketch ids: ";
    for (int i : empty_sketch_ids) {
        cout << i << " ";
    }
    cout << endl;
}


void compute_index_from_sketches() {
    // create the index using all the hashes
    for (uint i = 0; i < sketches.size(); i++) {
        for (uint j = 0; j < sketches[i].size(); j++) {
            hash_t hash = sketches[i][j];
            if (hash_index.find(hash) == hash_index.end()) {
                hash_index[hash] = vector<int>();
            }
            hash_index[hash].push_back(i);
        }
    }

    size_t num_hashes = hash_index.size();

    // remove the hashes that only appear in one sketch
    vector<hash_t> hashes_to_remove;
    for (auto it = hash_index.begin(); it != hash_index.end(); it++) {
        if (it->second.size() == 1) {
            hashes_to_remove.push_back(it->first);
        }
    }
    for (uint i = 0; i < hashes_to_remove.size(); i++) {
        hash_index.erase(hashes_to_remove[i]);
    }

    size_t num_hashes_after_removal = hash_index.size();

    cout << "Total number of distinct hashes: " << num_hashes << endl;
    cout << "Total number of distinct hashes that appear in only one sketch: " << num_hashes - num_hashes_after_removal << endl;
    cout << "Size of the index: " << num_hashes_after_removal << endl;

}


void compute_intersection_matrix_by_sketches(int sketch_start_index, int sketch_end_index, int thread_id, string out_dir, int pass_id, int negative_offset) {
    
    // process the sketches in the range [sketch_start_index, sketch_end_index)
    for (uint i = sketch_start_index; i < sketch_end_index; i++) {
        for (int j = 0; j < sketches[i].size(); j++) {
            hash_t hash = sketches[i][j];
            if (hash_index.find(hash) != hash_index.end()) {
                vector<int> sketch_indices = hash_index[hash];
                for (uint k = 0; k < sketch_indices.size(); k++) {
                    intersectionMatrix[i-negative_offset][sketch_indices[k]]++;
                }
            }
        }
    }

    // write the similarity values to file. filename: out_dir/passid_threadid.txt, where id is thread id in 3 digits
    string id_in_three_digits_str = to_string(thread_id);
    while (id_in_three_digits_str.size() < 3) {
        id_in_three_digits_str = "0" + id_in_three_digits_str;
    }
    string pass_id_str = to_string(pass_id);
    string filename = out_dir + "/" + pass_id_str + "_" + id_in_three_digits_str + ".txt";
    ofstream outfile(filename);

    // only write the values if larger than the threshold
    for (int i = sketch_start_index; i < sketch_end_index; i++) {
        for (uint j = 0; j < num_sketches; j++) {
            // skip obvious cases
            if ((uint)i == (uint)j) {
                continue;
            }

            // if nothing in the intersection, then skip
            if (intersectionMatrix[i-negative_offset][j] == 0) {
                continue;
            }

            // if either of the sketches is empty, then skip
            if (sketches[i].size() == 0 || sketches[j].size() == 0) {
                continue;
            }

            // if the divisor in the jaccard calculation is 0, then skip
            if (sketches[i].size() + sketches[j].size() - intersectionMatrix[i-negative_offset][j] == 0) {
                continue;
            }

            double jaccard = 1.0 * intersectionMatrix[i-negative_offset][j] / ( sketches[i].size() + sketches[j].size() - intersectionMatrix[i-negative_offset][j] );
            double containment_i_in_j = 1.0 * intersectionMatrix[i-negative_offset][j] / sketches[i].size();
            double containment_j_in_i = 1.0 * intersectionMatrix[i-negative_offset][j] / sketches[j].size();
            
            // containment_i_in_j is the containment of query in target, i is the query
            if (containment_i_in_j < arguments.containment_threshold) {
                continue;
            }

            outfile << i << "," << j << "," << jaccard << "," << containment_i_in_j << "," << containment_j_in_i << endl;
            similars[i].push_back(j);
        }
    }

    outfile.close();

}



void compute_intersection_matrix() {
    // allocate memory for the intersection matrix
    int num_sketches_each_pass = ceil(1.0 * num_sketches / arguments.num_of_passes);
    intersectionMatrix = new int*[num_sketches_each_pass + 1];
    for (int i = 0; i < num_sketches_each_pass + 1; i++) {
        intersectionMatrix[i] = new int[num_sketches];
    }

    // allocate memory for the similars
    for (int i = 0; i < num_sketches; i++) {
        similars.push_back(vector<int>());
    }

    for (int pass_id = 0; pass_id < arguments.num_of_passes; pass_id++) {
        // set zeros in the intersection matrix
        for (int i = 0; i < num_sketches_each_pass+1; i++) {
            for (uint j = 0; j < num_sketches; j++) {
                intersectionMatrix[i][j] = 0;
            }
        }

        // prepare the indices which will be processed in this pass
        int sketch_idx_start_this_pass = pass_id * num_sketches_each_pass;
        int sketch_idx_end_this_pass = (pass_id == arguments.num_of_passes - 1) ? num_sketches : (pass_id + 1) * num_sketches_each_pass;
        int negative_offset = pass_id * num_sketches_each_pass;
        int num_sketches_this_pass = sketch_idx_end_this_pass - sketch_idx_start_this_pass;
        
        // create threads
        vector<thread> threads;
        int chunk_size = num_sketches_this_pass / arguments.number_of_threads;
        for (int i = 0; i < arguments.number_of_threads; i++) {
            int start_index_this_thread = sketch_idx_start_this_pass + i * chunk_size;
            int end_index_this_thread = (i == arguments.number_of_threads - 1) ? sketch_idx_end_this_pass : sketch_idx_start_this_pass + (i + 1) * chunk_size;
            threads.push_back( thread(compute_intersection_matrix_by_sketches, start_index_this_thread, end_index_this_thread, i, arguments.working_directory, pass_id, negative_offset) );
        }

        // join threads
        for (int i = 0; i < arguments.number_of_threads; i++) {
            threads[i].join();
        }

        // show progress
        std::cout << "Pass " << pass_id+1 << "/" << arguments.num_of_passes << " done." << std::endl;
    }

    // free the memory allocated for the intersection matrix
    for (int i = 0; i < num_sketches_each_pass + 1; i++) {
        delete[] intersectionMatrix[i];
    }
    delete[] intersectionMatrix;
}




void do_yacht_train() {

    cout << "Starting yacht train..." << endl;
    
    vector<int> selected_genome_ids;
    vector<bool> genome_id_to_exclude(num_sketches, false);

    // sort the genome ids by size
    sort(genome_id_size_pairs.begin(), genome_id_size_pairs.end(), [](const pair<int, int>& a, const pair<int, int>& b) {
        return a.second < b.second;
    });
    
    for (int i = 0; i < num_sketches; i++) {
        
        cout << "Processing " << i << "..." << '\r';
        int genome_id_this = genome_id_size_pairs[i].first;
        int size_this = genome_id_size_pairs[i].second;
        bool select_this = true;
        
        // show my size
        for (int j = 0; j < similars[genome_id_this].size(); j++) {
            int genome_id_other = similars[genome_id_this][j];
            if (genome_id_to_exclude[genome_id_other]) {
                continue;
            }
            int size_other = sketches[genome_id_other].size();
            if (size_other >= size_this) {
                select_this = false;
                break;
            }
        }
        if (select_this) {
            selected_genome_ids.push_back(genome_id_this);
        } else {
            genome_id_to_exclude[genome_id_this] = true;
        }
    }

    cout << "Writing to output file.." << endl; 

    // write the selected genome ids to file
    ofstream outfile(arguments.output_filename);
    for (int i = 0; i < selected_genome_ids.size(); i++) {
        int genome_id = selected_genome_ids[i];
        string sketch_name = sketch_names[genome_id];
        outfile << sketch_name << endl;
    }
    outfile.close();

}




int main(int argc, char *argv[]) {

    // *********************************************************
    // *****           parse command line arguments       ******
    // *********************************************************
    try {
        parse_arguments(argc, argv);
    } catch (const std::runtime_error &e) {
        std::cerr << e.what() << std::endl;
        cout << "Usage: " << argv[0] << " -h" << endl;
        return 1;
    }

    // show the arguments
    show_arguments();

    // *********************************************************
    // ************     read the input sketches     ************
    // *********************************************************
    auto read_start = chrono::high_resolution_clock::now();
    cout << "Reading all sketches in filelist using all " << arguments.number_of_threads << " threads..." << endl;
    get_sketch_names(arguments.file_list);
    cout << "Total number of sketches to read: " << num_sketches << endl;
    read_sketches();
    auto read_end = chrono::high_resolution_clock::now();
    
    cout << "All sketches read" << endl;
    
    // show empty sketches
    show_empty_sketches();

    // show time taken to read all sketches
    auto read_duration = chrono::duration_cast<chrono::milliseconds>(read_end - read_start);
    cout << "Time taken to read all sketches: " << read_duration.count() << " milliseconds" << endl;




    // ****************************************************************
    // ************* reading complete, now creating index *************
    // ****************************************************************
    auto index_build_start = chrono::high_resolution_clock::now();
    cout << "Building index from sketches..." << endl;
    compute_index_from_sketches();
    auto index_build_end = chrono::high_resolution_clock::now();
    auto index_build_duration = chrono::duration_cast<chrono::milliseconds>(index_build_end - index_build_start);
    cout << "Time taken to build index: " << index_build_duration.count() << " milliseconds" << endl;



    // **********************************************************************
    // *************        compute intersection matrix         *************
    // **********************************************************************
    auto mat_computation_start = chrono::high_resolution_clock::now();
    cout << "Computing intersection matrix..." << endl;
    compute_intersection_matrix();
    auto mat_computation_end = chrono::high_resolution_clock::now();
    auto mat_computation_duration = chrono::duration_cast<chrono::milliseconds>(mat_computation_end - mat_computation_start);
    cout << "Time taken to compute intersection matrix: " << mat_computation_duration.count() << " milliseconds" << endl;



    // **********************************************************************
    // *************                yacht train                 *************
    // **********************************************************************
    auto yacht_train_start = chrono::high_resolution_clock::now();
    cout << "Starting yacht train..." << endl;
    do_yacht_train();
    auto yacht_train_end = chrono::high_resolution_clock::now();
    auto yacht_train_duration = chrono::duration_cast<chrono::milliseconds>(yacht_train_end - yacht_train_start);
    cout << "Time taken to do yacht train: " << yacht_train_duration.count() << " milliseconds" << endl;


    return 0;
}