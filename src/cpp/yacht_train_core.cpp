/*
 * Author: Mahmudur Rahman Hera (mahmudhera93@gmail.com)
 * Date: November 1, 2024
 * Description: yacht train core using indexing of sketches to do genome comparison
 */

#include "argparse.hpp"
#include "json.hpp"
#include "utils.h"
#include "MultiSketchIndex.h"

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



void parse_arguments(int argc, char *argv[], Arguments &arguments) {

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


void show_arguments(Arguments& arguments) {
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




void do_yacht_train(const vector<vector<hash_t>>& sketches, 
                    const vector<vector<int>>& similars, 
                    const vector<string>& sketch_paths,
                    const string& output_filename) {

    cout << "Starting yacht train..." << endl;

    int num_sketches = sketches.size();
    
    vector<int> selected_genome_ids;
    vector<bool> genome_id_to_exclude(num_sketches, false);

    // sort the genome ids by size
    vector<pair<int, int>> genome_id_size_pairs;
    for (int i = 0; i < num_sketches; i++) {
        genome_id_size_pairs.push_back(make_pair(i, sketches[i].size()));
    }
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
            
            if (genome_id_other == genome_id_this) {
                continue;
            }

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

    // write the selected sketch paths to file
    ofstream outfile(output_filename);
    for (int i = 0; i < selected_genome_ids.size(); i++) {
        int genome_id = selected_genome_ids[i];
        string sketch_name = sketch_paths[genome_id];
        outfile << sketch_name << endl;
    }
    outfile.close();

}




int main(int argc, char *argv[]) {

    Arguments arguments;

    std::vector<std::string> sketch_paths;
    vector<vector<hash_t>> sketches;
    MultiSketchIndex ref_sketches_index;
    mutex mutex_count_empty_sketch;
    vector<int> empty_sketch_ids;
    int ** intersectionMatrix;
    vector<vector<int>> similars;

    // *********************************************************
    // *****           parse command line arguments       ******
    // *********************************************************
    try {
        parse_arguments(argc, argv, arguments);
    } catch (const std::runtime_error &e) {
        std::cerr << e.what() << std::endl;
        cout << "Usage: " << argv[0] << " -h" << endl;
        return 1;
    }

    // show the arguments
    show_arguments(arguments);

    // *********************************************************
    // ************     read the input sketches     ************
    // *********************************************************
    auto read_start = chrono::high_resolution_clock::now();
    cout << "Reading all sketches in filelist using all " << arguments.number_of_threads << " threads..." << endl;
    get_sketch_paths(arguments.file_list, sketch_paths);
    cout << "Total number of sketches to read: " << sketch_paths.size() << endl;
    read_sketches(sketch_paths, sketches, empty_sketch_ids, arguments.number_of_threads);
    auto read_end = chrono::high_resolution_clock::now();
    
    cout << "All sketches read" << endl;

    
    // show empty sketches
    show_empty_sketches(empty_sketch_ids);

    // show time taken to read all sketches
    auto read_duration = chrono::duration_cast<chrono::milliseconds>(read_end - read_start);
    cout << "Time taken to read all sketches: " << read_duration.count() << " milliseconds" << endl;




    // ****************************************************************
    // ************* reading complete, now creating index *************
    // ****************************************************************
    auto index_build_start = chrono::high_resolution_clock::now();
    cout << "Building index from sketches..." << endl;
    compute_index_from_sketches(sketches, ref_sketches_index, arguments.number_of_threads);
    auto index_build_end = chrono::high_resolution_clock::now();
    auto index_build_duration = chrono::duration_cast<chrono::milliseconds>(index_build_end - index_build_start);
    cout << "Time taken to build index: " << index_build_duration.count() << " milliseconds" << endl;



    // **********************************************************************
    // *************        compute intersection matrix         *************
    // **********************************************************************
    auto mat_computation_start = chrono::high_resolution_clock::now();
    cout << "Computing intersection matrix..." << endl;
    compute_intersection_matrix(sketches, sketches, ref_sketches_index, 
                                arguments.working_directory, similars, 
                                arguments.containment_threshold, arguments.num_of_passes, 
                                arguments.number_of_threads);
    auto mat_computation_end = chrono::high_resolution_clock::now();
    auto mat_computation_duration = chrono::duration_cast<chrono::milliseconds>(mat_computation_end - mat_computation_start);
    cout << "Time taken to compute intersection matrix: " << mat_computation_duration.count() << " milliseconds" << endl;



    // **********************************************************************
    // *************                yacht train                 *************
    // **********************************************************************
    auto yacht_train_start = chrono::high_resolution_clock::now();
    cout << "Starting yacht train..." << endl;
    do_yacht_train(sketches, similars, sketch_paths, arguments.output_filename);
    auto yacht_train_end = chrono::high_resolution_clock::now();
    auto yacht_train_duration = chrono::duration_cast<chrono::milliseconds>(yacht_train_end - yacht_train_start);
    cout << "Time taken to do yacht train: " << yacht_train_duration.count() << " milliseconds" << endl;


    return 0;
}