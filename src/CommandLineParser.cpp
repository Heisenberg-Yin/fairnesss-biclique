#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <sys/resource.h>
#include <fstream>
#include <sstream>
#include <cassert>
#include "CommandLineParser.h"
#include "BipartiteGraph.h"

using namespace std;

CommandLineParser::CommandLineParser(int argc, char **argv) {

    if (argc == 1) {
        append_to_output(usage);
    }
    else if (argc < 6) {
        printf("%d\n",argc);
        string arg = argv[1];

        if ((arg == "-h") or (arg == "--help")) {
            append_to_output(usage);
        }

        else {
            append_to_output(bad_input_str);
            append_to_output(too_less_args_str);
            append_to_output(usage);
        }
    }

    else if (argc == 8) {
        string arg1 = argv[1]; //dir
        string arg2 = argv[2]; //one side or both side,1 means one side and 2 means two sides, one side means V side
        string arg3 = argv[3]; //fair number for each attribute in U alpha
        string arg4 = argv[4]; //fair number for each attribute in V beta
        string arg5 = argv[5]; //fair gap delta 
        string arg6 = argv[6];
        string arg7 = argv[7];
        not_file_str = string("MBEA error: " + arg1 + " is not a file.");
        ui side=atoi(argv[2]);
        ui alpha=atoi(argv[3]);
        ui beta=atoi(argv[4]);
        ui delta=atoi(argv[5]);                        
        ui strategy=atoi(argv[6]); 
        ui algorithm=atoi(argv[7]);
        if(side==1){
            BipartiteGraph finder =BipartiteGraph(arg1.c_str(),alpha,beta,delta,side,strategy);           
            if(algorithm==1){
                finder.naive_enumeration_one_side();
            }else if(algorithm==2){
                finder.one_side_get_num_fair_bicliques();  
            }else if(algorithm==3){
                finder.get_num_bicliques(); 
            }             
        }else if(side==2){
            BipartiteGraph finder =BipartiteGraph(arg1.c_str(),alpha,beta,delta,side,strategy); 
            if(algorithm==1){
                finder.naive_enumeration_two_side();
            }else if(algorithm==2){
                finder.two_side_get_num_fair_bicliques();  
            }  
            //finder.naive_enumeration_two_side();            
            //finder.two_side_get_num_fair_bicliques();
        }
                     
    }
    else {
        append_to_output(bad_input_str);
        append_to_output(too_many_args_str);
        append_to_output(usage);
    }

}

void CommandLineParser::append_to_output(string to_append) {
    output += to_append + "\n";
}


string CommandLineParser::get_output() {
    return output;
}

string CommandLineParser::get_usage() {
    return usage;
}

string CommandLineParser::get_too_many_args_str() {
    return too_many_args_str;
}

string CommandLineParser::get_not_file_str() {
    return not_file_str;
}

string CommandLineParser::get_bad_input_str() {
    return bad_input_str;
}

string CommandLineParser::get_alg_error() {
    return alg_error;
}

string CommandLineParser::get_algorithm() {
    return algorithm;
}

