//
//  Copyright 2002-2007 Rick Desper, Olivier Gascuel
//  Copyright 2007-2014 Olivier Gascuel, Stephane Guindon, Vincent Lefort
//
//  This file is part of FastME.
//
//  FastME is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  FastME is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with FastME.  If not, see <http://www.gnu.org/licenses/>
//


#ifndef FASTME_H_
#define FASTME_H_

#include "interface_options.h"
#include "NNI.h"
#include "bNNI.h"
#include "SPR.h"
#include "random.h"
#include "distance.h"
#include "initialiser.h"

int run(double** d, int n_taxa, int** init_adj, int argc, char **argv);

void printFinalData (Options *options, char **bootTrees, char **matStr);

char **InitMatStrings (int numBoot, int numSpecies, int precision);

// int **rndForBootstraps (Options *options, int len);

// int **bootFilter (int nBoot, int len);

// int **p_bootPositions (int nBoot, int len);

void freeIntMat (int **mat, int size);

void OpenFiles (Options *options);

void printOptions (Options *options);

void printMatrix (double **D, int size, set *nodes, FILE *ofile, int input_type, int precision);

void printMatrixStr (double **D, int size, set *nodes, char *str, int input_type, int precision);

void PrintTimeInfo (time_t t_beg, time_t t_end);

void PrintEstimatedMemorySpace (int nbtax, int nbsites, Options *options);

// int PrintBootstrapInfo (Options *options, int repCounter, int printedRep);

void InitSpeciesAndTrees (Options *options, set *species, char **bootTrees, char *bestTree);
	
tree *ComputeTree (Options *options, double **D, double **A, set *species,
	int numSpecies, int precision);
	
tree *ImproveTree (Options *options, tree *T0, double **D, double **A,
	int *nniCount, int *sprCount, FILE *ofile);

void explainedVariance (double **D, tree *T, int n, int precision,
	int input_type, FILE *fpO);

#endif /*FASTME_H_*/

