#include "s_energy_matrix.hh"
#include "hotspot.hh"
#include "SHAPE.hh"
#include "cmdline.hh"
// a simple driver for the HFold
#include <algorithm>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <stdio.h>
#include <string>
#include <sys/stat.h>

bool exists(const std::string path) {
    struct stat buffer;
    return (stat(path.c_str(), &buffer) == 0);
}

void get_input(std::string file, std::string &sequence) {
    if (!exists(file)) {
        std::cout << "Input file does not exist" << std::endl;
        exit(EXIT_FAILURE);
    }
    std::ifstream in(file.c_str());
    std::string str;
    int i = 0;
    while (getline(in, str)) {
        if (str[0] == '>') continue;
        if (i == 0) sequence = str;
        ++i;
    }
    in.close();
}

// check if sequence is valid with regular expression
// check length and if any characters other than GCAUT
void validateSequence(std::string sequence) {

    if (sequence.length() == 0) {
        std::cerr << "Sequence is missing" << std::endl;
        exit(EXIT_FAILURE);
    }
    // return false if any characters other than GCAUT -- future implement check based on type
    for (char c : sequence) {
        if (!(c == 'G' || c == 'C' || c == 'A' || c == 'U' || c == 'T' || c == 'N')) {
            std::cerr << "Sequence contains character " << c << " that is not N,G,C,A,U, or T." << std::endl;
            exit(EXIT_FAILURE);
        }
    }
}

void seqtoRNA(std::string &sequence) {
    for (char &c : sequence) {
        if (c == 'T') c = 'U';
    }
}

int main(int argc, char *argv[]) {
    args_info args_info;

    // get options (call getopt command line parser)
    if (cmdline_parser(argc, argv, &args_info) != 0) {
        exit(1);
    }

    std::string seq;
    if (args_info.inputs_num > 0) {
        seq = args_info.inputs[0];
    } else {
        if (!args_info.input_file_given) std::getline(std::cin, seq);
    }

    if (args_info.input_file_given) {
        if (exists(args_info.input_file_arg)) {
            get_input(args_info.input_file_arg, seq);
        }
        if (seq.length()== 0) {
            std::cerr << "sequence is missing from file" << std::endl;
            exit(EXIT_FAILURE);
        }
    }
    cand_pos_t n = seq.length();
    std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
    if (!args_info.noConv_given) seqtoRNA(seq);
    validateSequence(seq);
    std::string shapefile = args_info.shape_given ? args_info.shape_arg : "";
    SHAPEData ShapeData(shapefile,n);

    std::string file = args_info.paramFile_given ? args_info.paramFile_arg : "params/rna_DirksPierce09.par";
    if (exists(file)) {
        vrna_params_load(file.c_str(), VRNA_PARAMETER_FORMAT_DEFAULT);
    } else if (seq.find('T') != std::string::npos) {
        vrna_params_load_DNA_Mathews2004();
    }

    // Hotspots
    std::vector<Hotspot> hotspot_list;
    s_energy_matrix min_fold(seq,&ShapeData);
    min_fold.params_->model_details.dangles = args_info.dangles_arg;
    min_fold.get_hotspots(hotspot_list,args_info.subopt_arg,args_info.stem_arg,args_info.turn_arg);
    // output to file

    if (args_info.output_file_given) {
        std::ofstream out(args_info.output_file_arg);
        for (size_t i = 0; i < hotspot_list.size(); i++) {
            out << hotspot_list[i].get_left_outer_index() << "\t" << hotspot_list[i].get_left_inner_index() << "\t" << hotspot_list[i].get_right_inner_index() << "\t" << hotspot_list[i].get_right_outer_index() << "\t" << hotspot_list[i].get_energy() << std::endl;
        }
        out.close();

    } else {
        for (size_t i = 0; i < hotspot_list.size(); i++) {
            std::cout << hotspot_list[i].get_left_outer_index() << "\t" << hotspot_list[i].get_left_inner_index() << "\t" << hotspot_list[i].get_right_inner_index() << "\t" << hotspot_list[i].get_right_outer_index() << "\t" << hotspot_list[i].get_energy() << std::endl;
        }
    }
    cmdline_parser_free(&args_info);
    return 0;
}
