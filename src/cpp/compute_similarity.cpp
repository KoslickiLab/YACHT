/*
 * Author: Mahmudur Rahman Hera (mahmudhera93@gmail.com)
 * Date: November 1, 2024
 * Description: yacht train core using indexing of sketches to do genome comparison
 */

#include "argparse.hpp"
#include "json.hpp"
#include "utils.h"

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
    string file_list_query;
    string file_list_target;
    string output_directory;
    int number_of_threads;
    int num_of_passes;
    double containment_threshold;
    bool combine;
    string combined_output_filename;
};


typedef Arguments Arguments;
typedef unsigned long long int hash_t;


void parse_arguments(int argc, char *argv[], Arguments &arguments) {

    argparse::ArgumentParser parser("yacht train using indexing of sketches");
    
    parser.add_argument("file_list_query")
        .help("file containing paths of query sketches")
        .store_into(arguments.file_list_query);
    
    parser.add_argument("file_list_target")
        .help("file containing paths of target sketches")
        .store_into(arguments.file_list_target);
    
    parser.add_argument("output_directory")
        .help("output directory (where similarity values will be written)")
        .store_into(arguments.output_directory);
    
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
    
    // argument: combine files (store true if present)
    parser.add_argument("-C", "--combine")
        .help("combine the output files")
        .default_value(false)
        .implicit_value(true)
        .store_into(arguments.combine);

    // argument: combined output filename
    parser.add_argument("-o", "--output")
        .help("output filename (where the combined output will be written. Not used if -C is not present)")
        .default_value("combined_output.txt")
        .store_into(arguments.combined_output_filename);
    
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


void show_arguments(Arguments& args) {
    cout << "**************************************" << endl;
    cout << "*" << endl;
    cout << "*    file_list_query: " << args.file_list_query << endl;
    cout << "*    file_list_target: " << args.file_list_target << endl;
    cout << "*    output_directory: " << args.output_directory << endl;
    cout << "*    number_of_threads: " << args.number_of_threads << endl;
    cout << "*    num_of_passes: " << args.num_of_passes << endl;
    cout << "*    containment_threshold: " << args.containment_threshold << endl;
    cout << "*    combine: " << bool(args.combine) << endl;
    cout << "*    combined_output_filename: " << args.combined_output_filename << endl;
    cout << "*" << endl;
    cout << "**************************************" << endl;
}



int main(int argc, char** argv) {
    Arguments args;

    // parse the arguments
    try {
        parse_arguments(argc, argv, args);
    } catch (const std::runtime_error &e) {
        std::cerr << e.what() << std::endl;
        cout << "Usage: " << argv[0] << " -h" << endl;
        return 1;
    }

    // show the arguments
    show_arguments(args);

    // read the query sketches
    cout << "Reading query sketches..." << endl;
    vector<string> sketch_paths_query;
    vector<vector<hash_t>> sketches_query;
    vector<int> empty_sketch_ids_query;
    get_sketch_paths(args.file_list_query, sketch_paths_query);
    read_sketches(sketch_paths_query, sketches_query, empty_sketch_ids_query, args.number_of_threads);
    cout << "All query sketches read" << endl;

    // read the target sketches
    cout << "Reading target sketches..." << endl;
    vector<string> sketch_paths_target;
    vector<vector<hash_t>> sketches_target;
    vector<int> empty_sketch_ids_target;
    get_sketch_paths(args.file_list_target, sketch_paths_target);
    read_sketches(sketch_paths_target, sketches_target, empty_sketch_ids_target, args.number_of_threads);

    // show empty sketches
    cout << "Empty sketches in query:" << endl;
    show_empty_sketches(empty_sketch_ids_query);
    cout << "Empty sketches in target:" << endl;
    show_empty_sketches(empty_sketch_ids_target);

    // compute the index from the target sketches
    cout << "Building index from target sketches..." << endl;
    unordered_map<hash_t, vector<int>> hash_index_target;
    compute_index_from_sketches(sketches_target, hash_index_target);

    // compute the similarity matrix
    cout << "Computing similarity matrix..." << endl;
    vector<vector<int>> similars;
    compute_intersection_matrix(sketches_query, sketches_target, 
                                hash_index_target, 
                                args.output_directory, similars, 
                                args.containment_threshold, 
                                args.num_of_passes, 
                                args.number_of_threads);

    cout << "similarity computation completed, results are here: " << args.output_directory << endl;

    if (args.combine) {
        cout << "Combining the output files..." << endl;
        // cat all the files in the output directory
        string command = "cat " + args.output_directory + "/*.txt > " + args.combined_output_filename;
        system(command.c_str());
        cout << "Combined output written to: " << args.combined_output_filename << endl;
    }

    return 0;
}