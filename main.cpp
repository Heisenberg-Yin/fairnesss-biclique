#include <iostream>
#include <string>
#include "BipartiteGraph.h"
#include <vector>
#include "Vertex.h"
#include "BicliqueFinder.h"
#include <algorithm>
#include <sys/resource.h>
#include <fstream>
#include "CommandLineParser.h"


int main(int argc, char* argv[]) {


    CommandLineParser CLP = CommandLineParser(argc, argv);
    std::cout << CLP.get_output() << std::endl;

    return 0;
}