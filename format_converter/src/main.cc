#include "ConstraintFile.H"
#include "SpecificationFile.h"

#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>

using namespace std;

void generate_input_file(const SpecificationFile &specificationFile,
        const ConstraintFile &constraintFile,
        const std::string &acts_inputfile_name) {
  
    std::ofstream acts_infile(acts_inputfile_name);
    if (!acts_infile.is_open()) {
        std::cerr << "open failed" << std::endl;
        exit(0);
    }

    acts_infile << "[System]" << std::endl;
    acts_infile << "Name: " << acts_inputfile_name << std::endl << std::endl;
    acts_infile << "[Parameter]" << std::endl;

    const Options &options = specificationFile.getOptions();
    for (unsigned option = 0; option < options.size(); ++option) {
        acts_infile << 'p' << option << "(int): ";
        acts_infile << 0;
        for (unsigned var = 1; var < options.symbolCount(option); ++var) {
            acts_infile << ',' << var;
        }
        acts_infile << std::endl;
    }
    acts_infile << std::endl << "[Constraint]" << std::endl;
    const Valid::Formula &formula = constraintFile.getFormula();
    for (auto &c : formula) {
        const unsigned option = options.option(c[0].variable());
        acts_infile << 'p' << option;
        c[0].is_negative() ? acts_infile << "!=" : acts_infile << "=";
        acts_infile << c[0].variable() - options.firstSymbol(option);
        for (unsigned i = 1; i < c.size(); ++i) {
            acts_infile << " || ";
            const unsigned option = options.option(c[i].variable());
            acts_infile << 'p' << option;
            c[i].is_negative() ? acts_infile << "!=" : acts_infile << "=";
            acts_infile << c[i].variable() - options.firstSymbol(option);
        }
        acts_infile << ' ' << std::endl;
    }
    acts_infile.close();
}

int main(int argc, char const *argv[]) {

    if (argc != 4) {
        cout << "Error\n";
        return -1;
    }
    
    string modelFile = argv[1];
    string constrFile = argv[2];
    string actsFile = argv[3];

    // cout << modelFile << " " << constrFile << " " << actsFile << "\n";

    SpecificationFile specificationFile(modelFile);
    ConstraintFile constraintFile(constrFile);
    generate_input_file(specificationFile, constraintFile, actsFile);

    
    return 0;
}