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

int** fill_int_matrix() {
    int** matrix = new int*[28];
    std::fstream file("../init_mat");
    std::string line;
    // std::cout<<" ***\n";
    int row=0;
    while (getline(file, line) and row < 28){
        std::stringstream ss( line );                     
        std::string data;
        int col = 0;
        matrix[row] = new int[28];
        while (getline( ss, data, ' ' ) )           
        {
            matrix[row][col]=stoi(data);
            // std::cout<<matrix[row][col]<<" ";
            col++;
            
        }
        // std::cout<<"\n";
        row++;
    }
    return matrix;
}


double** fill_matrix() {
    double** matrix = new double*[15];
    std::fstream file("../mat");
    std::string line;
    // std::cout<<" ***\n";
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
            // std::cout<<stod(data)<<" ";
        }
        // std::cout<<"\n";
        row++;
    }
    return matrix;
}

int main () {
    // printf("ok\n");
    // Options *options;
    // options = chooseSettings (argc, argv);
    double** D = fill_matrix();
    int** init_adj = fill_int_matrix();
    // printf("ok\n");
    
    // D = loadMattrix(3);
    // D = loadM (options->fpI_data_file, &numSpecies, species);
    // printf("ok\n");
    run(D, 15, init_adj);
}





// 14355450028247832
// 14673830625000001