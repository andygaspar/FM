#include "fastme.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <fcntl.h>
#include <getopt.h>
#include <ctype.h>
#include <errno.h>
#include <float.h>
#include <libgen.h>
#include <assert.h>
#include <iostream>
#include <fstream>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>


using string = std::string;


double** fill_matrix() {
    double** matrix = new double*[15];
    std::fstream file("mat");
    std::string line;
    std::cout<<" ***\n";
    int row=0;
    while (getline(file, line) and row < 15){
        std::stringstream ss( line );                     
        std::string data;
        int col = 0;
        matrix[row] = new double[15];
        while (getline( ss, data, ' ' ) )           
        {
            matrix[row][col]=stod(data);
            col++;
            std::cout<<stod(data)<<" ";
        }
        std::cout<<"\n";
        row++;
    }
    return matrix;
}

int main (int argc, char **argv) {
    // printf("ok\n");
    // Options *options;
    // options = chooseSettings (argc, argv);
    double** D = fill_matrix();
    printf("ok\n");
    
    // D = loadMattrix(3);
    // D = loadM (options->fpI_data_file, &numSpecies, species);
    // printf("ok\n");
    run(D, 15, argc, argv);
}



