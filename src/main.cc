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
        std::cout << "sequence is missing" << std::endl;
        exit(EXIT_FAILURE);
    }
    // return false if any characters other than GCAUT -- future implement check based on type
    for (char c : sequence) {
        if (!(c == 'G' || c == 'C' || c == 'A' || c == 'U' || c == 'T' || c == 'N')) {
            std::cout << "Sequence contains character " << c << " that is not N,G,C,A,U, or T." << std::endl;
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

    std::string fileI = args_info.input_file_given ? args_info.input_file_arg : "";

    std::string fileO = args_info.output_file_given ? args_info.output_file_arg : "";

    int number_of_hotspots = args_info.subopt_given ? args_info.subopt_arg : 100;

    std::string shapeFile = args_info.shape_given ? args_info.shape_arg : "";

    int dangles = args_info.dangles_given ? args_info.dangles_arg : 2;

    int min_stem = args_info.stem_given ? args_info.stem_arg : 4;

    int turn = args_info.turn_given ? args_info.turn_arg : 3;


    if (fileI != "") {

        if (exists(fileI)) {
            get_input(fileI, seq);
        }
        if (seq == "") {
            std::cout << "sequence is missing from file" << std::endl;
        }
    }
    cand_pos_t n = seq.length();
    std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
    if (!args_info.noConv_given) seqtoRNA(seq);
    validateSequence(seq);

    std::string file = args_info.paramFile_given ? args_info.paramFile_arg : "params/rna_DirksPierce09.par";
    if (exists(file)) {
        vrna_params_load(file.c_str(), VRNA_PARAMETER_FORMAT_DEFAULT);
    } else if (seq.find('T') != std::string::npos) {
        vrna_params_load_DNA_Mathews2004();
    }

    SHAPEData ShapeData(shapeFile,n);

    cmdline_parser_free(&args_info);

    // Hotspots
    std::vector<Hotspot> hotspot_list;
    s_energy_matrix min_fold(seq,&ShapeData);
    min_fold.get_hotspots(hotspot_list,number_of_hotspots,min_stem,turn);
    // output to file

    number_of_hotspots = std::min(number_of_hotspots,(int) hotspot_list.size());
    if (fileO != "") {
        std::ofstream out(fileO);
        for (cand_pos_t i = 0; i < number_of_hotspots; i++) {
            out << hotspot_list[i].get_left_outer_index() << "\t" << hotspot_list[i].get_left_inner_index() << "\t" << hotspot_list[i].get_right_inner_index() << "\t" << hotspot_list[i].get_right_outer_index() << "\t" << hotspot_list[i].get_energy() << std::endl;
        }
        out.close();

    } else {
        for (cand_pos_t i = 0; i < number_of_hotspots; i++) {
            std::cout << hotspot_list[i].get_left_outer_index() << "\t" << hotspot_list[i].get_left_inner_index() << "\t" << hotspot_list[i].get_right_inner_index() << "\t" << hotspot_list[i].get_right_outer_index() << "\t" << hotspot_list[i].get_energy() << std::endl;
        }
    }
    return 0;
}
