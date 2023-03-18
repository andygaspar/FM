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


#include "fastme.h"

boolean isBoostrap;
int verbose;


void make_set(int n_taxa, set *S) {
	node *v = NULL;
	for(int i=0; i< n_taxa; i++) {
		char c= i + '0';
		v = makeNode (&c, -1);
		v->index2 = i;
		S = addToSet (v, S);
	}

}


void recursion_print (edge *e)
{	
	if (NULL!= e) {
	if(e -> head -> middleEdge == NULL) printf("rrrr");

	printf("node %i %i\n" ,e->head->index2, e->head->index);
	recursion_print(e->head->leftEdge);
	recursion_print(e->head->rightEdge);
	}	
}


void printTree (tree *T)
{	
	if(T-> root -> middleEdge == NULL) printf("rrrr");
	printf("size %i \n", T->size);
	printf("node %i %c%c\n" ,T->root->index2, T->root->label[0], T->root->label[1]);
	recursion_print(T->root->leftEdge);
	recursion_print(T->root->rightEdge);

}

int recursion_fill_adj (int from, int internal_idx, long** A, edge *e){
	if (NULL!= e) {

		if(e->head->index2 < 0){
			A[from][internal_idx] = 1;
			A[internal_idx][from] = 1;	
			from = internal_idx;
			internal_idx ++;	
		}
		else{
			A[from][e->head->index2] = 1;
			A[e->head->index2][from] = 1;
			from = e->head->index2;
		}

	
	internal_idx = recursion_fill_adj(from, internal_idx, A, e->head->leftEdge);
	internal_idx = recursion_fill_adj(from, internal_idx, A, e->head->rightEdge);
	
	}	
	return internal_idx;

}

void tree_to_adj_mat(tree *T) {
	long**  A;
	int size = T->size;
    A = (long** ) mCalloc(size, sizeof(long*)); 
    for(int i = 0;i<size;i++) {
        A[i] = (long* ) mCalloc(size, sizeof(long));
    }
	int internal_idx = (T-> size + 2) / 2;
	internal_idx = recursion_fill_adj(T->root->index, internal_idx, A, T->root->leftEdge);
	recursion_fill_adj(T->root->index, internal_idx, A, T->root->rightEdge);

	for(int i = 0;i<size;i++) {
		for(int j = 0; j < size; j++) {
			printf("%i ", A[i][j]);
		}
		printf("\n");
	}
}


// tree adj_to_tree(double** A, int n_taxa) {
// 	tree *T;
// 	node *root;

// 	T = newTree();
// 	char c= 0 + '0';
// 	root = makeNode(c, 0);
// }

/*********************************************************/

int run(double**d, int n_taxa, int argc, char **argv)
{
	Options *options;
	set *species, *species_bk;
	
	species = species_bk = NULL;
	
	time_t t_beg, t_end;

	tree *T = NULL;

#ifdef _OPENMP
	int i = 0;
#endif

	int numSpecies = n_taxa;
	int setCounter = 0;
	
	int nniCount, sprCount, repCounter, repToPrint, printedRep;
	
	int seqLength = -1;
	
	// Random numbers for bootstraps
	int **rnd = NULL;

	double **D, **A;
	D = A = NULL;

	// Strings containing the best tree and the bootstraped trees
	char *bestTree = NULL;
	char **bootTrees = NULL;
	
	// String to record pseudo-matrices
	char **matStr = NULL;
	
	// Input tree file stream position usefull for bootstraps
	fpos_t fpI_tree_pos;

	isBoostrap = FALSE;

	options = chooseSettings (argc, argv);
    printf( "\n . Method used %d", options -> method);
    printf( "\n");

   PrintOptions (options);

#ifdef _OPENMP
	Message ( (char*)"This analysis will run on %d threads.", options->nb_threads);
#else
	printf ("\n");
	fflush (stdout);
	fflush (stderr);
#endif

	time(&t_beg);

	OpenFiles (options);

	// mem alloc for best tree
	bestTree = (char *) mCalloc (MAX_INPUT_SIZE, sizeof (char));
	// if (options->nb_bootstraps > 0)
	// {
	// 	bootTrees = (char **) mCalloc (options->nb_bootstraps, sizeof (char *));
	// }

	// if (PROTEIN == options->input_type && PDIST != options->model)
	// {
	// 	PROTEINdatamodel = InitProtModel (options);
	// }

	sgenrand ( (unsigned long)options->seed);

	while ( setCounter < options->nb_datasets )
	{
		//nniCount = sprCount = tbrCount = repCounter = 0;
		nniCount = sprCount = repCounter = 0;
		setCounter++;
		species = (set *) mCalloc (1, sizeof(set));
		
		printf ("\n#  Analysing dataset %d\n", setCounter);

		if (setCounter == 1)
			printOptions (options);

		fprintf (options->fpO_stat_file,"Dataset %d\n", setCounter);
		if (options->use_NNI && options->use_SPR)
		{
			fprintf (options->fpO_stat_file,"\tTwo tree searches are performed, with NNIs and SPRs respectively.\n");
			fprintf (options->fpO_stat_file,"\tThe resulting tree is the best (shortest) of both.\n\n");
		}

		if (options->nb_datasets > 1)
		{
			// if (options->nb_bootstraps > 0)
			// {
			// 	if (setCounter > 1)
			// 		fprintf (options->fpO_boot_file,"\n");

			// 	fprintf (options->fpO_boot_file,"Dataset %d\n", setCounter);
			// }
			if (options->use_O_mat_file)
				fprintf (options->fpO_mat_file,"Dataset %d\n", setCounter);
		}

/*********************************************************
		GET MATRIX
**********************************************************/
		if (MATRIX == options->input_type)
		{
			// double** G;
			// G = loadM (options->fpI_data_file, &numSpecies, species);
			make_set(numSpecies, species);
			D = d;
			for(int i=0; i<numSpecies; i++) {
				for(int j=0; j<numSpecies; j++) printf("%f ", D[i][j]);
				printf("\n");
			}
			PrintEstimatedMemorySpace (numSpecies, 0, options);
		}

		checkTrgInequality (D, numSpecies, options->trg_ineq);

		species_bk = copySet (species);
		

/*********************************************************
		COMPUTE TREE
**********************************************************/
		if (numSpecies < 4)
			Warning ( (char*)"Cannot compute tree with less than 4 taxa.");
		if (! options->only_mat && numSpecies > 3)
		{
			if (USER == options->method)
				Message ( (char*)"Reading input tree...");
			else
				Message ( (char*)"Computing tree...");

			A = initDoubleMatrix (2*numSpecies-2);
			
			if (USER == options->method)
			{
				if (repCounter == 0)
					fgetpos (options->fpI_tree_file, &fpI_tree_pos);
				else
					fsetpos(options->fpI_tree_file, &fpI_tree_pos);
				T = loadNewickTree (options->fpI_tree_file, numSpecies);
				T = detrifurcate (T);
				compareSets (T, species);
				partitionSizes (T);
			}
			else
			{
				T = ComputeTree (options, D, A, species, numSpecies, options->precision);
			}

			T = ImproveTree (options, T, D, A, &nniCount, &sprCount, options->fpO_stat_file);
			printTree(T);
			tree_to_adj_mat(T);
		
			explainedVariance (D, T, numSpecies, options->precision, options->input_type, options->fpO_stat_file);
				
			if (verbose > 0)
			{
				Message ( (char*)"Performed %d NNI(s) and %d SPR(s) on dataset %d.",
						nniCount, sprCount, setCounter);
			}
			
		} //End if 'compute tree'


		if (! options->only_mat && numSpecies > 3){
			NewickPrintTree (T, options->fpO_tree_file, options->precision);
			}

#ifndef _OPENMP
		freeMatrix (A, 2*numSpecies-2);
#endif

		if (numSpecies > 3)
		{
			if (NULL != T)
				freeTree (T);
			if (! options->only_mat)
				fprintf (options->fpO_tree_file, "\n");
		}
		
		if (NULL != matStr)
			free (matStr);
		
		freeMatrix (D, numSpecies);
		freeSet (species);
		freeSet (species_bk);
		
	} //end datasets loop

	fflush (stdout);
	fflush (stderr);
	
	if (NULL != rnd)		free (rnd);
	if (NULL != bootTrees)	free (bootTrees);
	if (NULL != bestTree)	free (bestTree);

	Free_Input (options);

	time (&t_end);
	PrintTimeInfo (t_beg, t_end);
	
	//exit (EXIT_SUCCESS);

	return 0;
}

/*********************************************************/

void printFinalData (Options *options, char **bootTrees, char **matStr)
{
	int i;

	for (i=0; i<options->nb_bootstraps; i++)
	{
		// print pseudo-trees to file
		// free memory allocated for best tree and pseudo-trees strings
		fprintf (options->fpO_boot_file, "%s", bootTrees[i]);
		if (NULL != bootTrees[i])
			free (bootTrees[i]);

		// print pseudo-matrices to file
		// free memory allocated for pseudo-matrices strings
		if (options->use_O_mat_file && NULL != matStr)
		{
			fprintf (options->fpO_mat_file, "%s", matStr[i]);
			if (NULL != matStr[i])
				free (matStr[i]);
		}
	}
	
	return;
}

/*********************************************************/

char **InitMatStrings (int numBoot, int numSpc, int precision)
{
	int i, len;
	
	// Mem alloc requires:
	// + 8 chars on the first line (indicating the number of taxa)
	// for each taxa:
	// + MAX_NAME_LENGTH for taxa name
	// + the distances
	//   each distance requires a length of:
	//     + precision + 1 for '.'
	//     + 2 for digits before the dot
	//     + 2 for blank spaces
	// + 1 for '\n'
	
	len = 8 + ( numSpc * (MAX_NAME_LENGTH + ( numSpc * (precision+5) ) + 1));
	
	char **matStr = (char **) mCalloc (numBoot, sizeof (char *));
	
	for (i=0; i<numBoot; i++)
	{
		matStr[i] = (char *) mCalloc (len, sizeof (char));
	}
	
	return matStr;
}

/*********************************************************/

void freeIntMat (int **mat, int size)
{
	int i;

	for (i=0; i<size; i++)
		if (NULL != mat[i])
			free(mat[i]);
	
	return;
}

/*********************************************************/

void OpenFiles (Options *options)
{
	options->fpI_data_file = Openfile (options->I_data_file,  (char*)"r");
	options->fpO_stat_file = Openfile (options->O_stat_file, options->open_mode);

	if (! options->only_mat)
		options->fpO_tree_file = Openfile (options->O_tree_file, options->open_mode);
	
	if (options->nb_bootstraps > 0)
		options->fpO_boot_file = Openfile (options->O_boot_file, options->open_mode);

	if (options->use_O_mat_file)
		options->fpO_mat_file = Openfile (options->O_mat_file, options->open_mode);

	if (USER == options->method)
		options->fpI_tree_file = Openfile (options->I_tree_file, (char*)"r");

  return;
}

/*********************************************************/

void printOptions (Options *options)
{
	char *tmp;
	tmp = (char *) mCalloc (MAX_NAME_LENGTH, sizeof (char));

	FILE *f = options->fpO_stat_file;

	if (options->open_mode[0] == 'w')
	{
		fprintf (f, "\n - FastME %s - \n\n");
		if (!options->only_mat && (options->method == NJ || options->method == UNJ || options->method == BIONJ))
			fprintf (f, "\nPapers to be cited:\n");
		else
			fprintf (f, "\nPaper to be cited:\n");
		
		fprintf (f, "\nFastME 2.0 - A comprehensive, accurate and fast distance-based phylogeny inference program.");
		fprintf (f, "\n\tVincent Lefort, Richard Desper and Olivier Gascuel,");
		fprintf (f, "\n\tMolecular Biology and Evolution 32(10), 2798-800, 2015.");
		
		if (!options->only_mat)
		{
			
			if (options->method == NJ)
			{
				fprintf (f, "\nNJ algorithm:");
				fprintf (f, "\n\tSaitou N., Nei M. 1987. The neighbor-joining method: a new method for reconstructing phylogenetic trees.");
				fprintf (f, "\n\tMolecular Biology and Evolution, 4(4):406-25");
			}
		}
		
		fprintf (f, "\n\n-------------------------------------------------------------------------------\n");
		fprintf (f, "Settings for this run:\n\n");

//----------------------------------------------------------------------//

		constantToStr (options->input_type, tmp);
		fprintf (f, "  I "
			"                                     Input data type "
			" %-15s \n", tmp);


//----------------------------------------------------------------------//

		if (! options->only_mat)
		{
			fprintf (f, "  D "
				"                                  Number of datasets "
				" %-15d \n", options->nb_datasets);

//----------------------------------------------------------------------//
			constantToStr (options->method, tmp);
			fprintf (f, "  M "
				"                        Initial tree building method "
				" %-15s \n", tmp);

//----------------------------------------------------------------------//

			if (options->use_NNI)
			{
				constantToStr (options->NNI, tmp);
				fprintf (f, "  N "
					"                                  NNI postprocessing "
					" %-15s \n", tmp);
			}
			else
			{
				fprintf (f, "  N "
					"                                  NNI postprocessing "
					" %-15s \n", "no");
			}

//----------------------------------------------------------------------//

			fprintf (f, "  S "
				"                                  SPR postprocessing "
				" %-15s \n", (options->use_SPR ? "yes" : "no"));

//----------------------------------------------------------------------//

			if (!options->use_NNI && !options->use_SPR && options->branch != NONE)
			{
				constantToStr (options->branch, tmp);
				fprintf (f, "  W "
					"             Branch lengths assigned to the topology "
					" %-15s \n", tmp);
			}

//----------------------------------------------------------------------//
/*
			fprintf (f, "  T "
				"                                  TBR postprocessing "
				" %-15s \n", (options->use_TBR ? "yes" : "no"));*/

//----------------------------------------------------------------------//

			if (options->nb_bootstraps > 0)
			{
				fprintf (f, "\n");
				fprintf (f, "  B "
					"                     Bootstrap: number of replicates "
					" %-15d \n", options->nb_bootstraps);
			}

		}

		fprintf (f, "\n-------------------------------------------------------------------------------\n");
	}

	free (tmp);

	return;
}

/*********************************************************/

void printMatrix (double **D, int size, set *nodes, FILE *ofile, int input_type, int precision)
{
	char *NAstr;
	char format[8];
	int i, val;
	node *v, *w;
	set *S, *T;
	double threshold, d;

	NAstr = (char *) mCalloc (DECIMAL_DIG +2, sizeof(char));
	// Build string of length depending on 'precision' to write
	strncat (NAstr, "   NA", 5);
	// Add [precision - 2] blank spaces
	for (i=0; i<=precision; i++)
		strncat (NAstr, " ", 1);
					
	snprintf (format, 8, "  %%.%df", precision);

	if (PROTEIN == input_type)
		threshold = PROT_DIST_MAX;
	else
		threshold = DNA_DIST_MAX;

	fprintf (ofile, "%d\n", size);
	for (S=nodes; NULL!=S; S=S->secondNode)
	{
		v = S->firstNode;
		fprintf (ofile, "%-10s ", v->label);
		for (T=nodes; NULL!=T; T=T->secondNode)
		{
			w = T->firstNode;
			d = D[v->index2][w->index2];
			if (d > threshold)
			{
				// The following code generates warning because 'NAstr' is not a string literal
				fprintf (ofile, NAstr);
			}
			else
			{
				// Get integer value of the distance
				val = (int) trunc (d);
				// The string containing the distnce has 3 chars before the dot
				if (val <= -10)
				{
					fprintf (ofile, format, d);
				}
				// 2 chars before the dot
				else if (val < 0)
				{
					fprintf (ofile, "%s", " ");
					fprintf (ofile, format, d);
				}
				// 1 char before the dot
				else if (val < 10)
				{
					fprintf (ofile, "%s", "  ");
					fprintf (ofile, format, d);
				}
				// 2 chars before the dot
				else
				{
					fprintf (ofile, "%s", " ");
					fprintf (ofile, format, d);
				}
			}
		}
		fprintf (ofile, "\n");
	}
	fprintf (ofile, "\n");
	free (NAstr);

	return;
}

/*********************************************************/

void printMatrixStr (double **D, int size, set *nodes, char *str, int input_type, int precision)
{
	int i;
	node *v, *w;
	set *S, *T;
	char *tmp;
	char format[8];
	double threshold, d;

	snprintf (format, 8, "%%.%df  ", precision);

	if (PROTEIN == input_type)
		threshold = PROT_DIST_MAX;
	else
		threshold = DNA_DIST_MAX;
	
	tmp = (char *) mCalloc (INPUT_SIZE, sizeof(char));
	if (strlen (tmp))
		strncpy (tmp, "", strlen(tmp));
	snprintf (tmp, INPUT_SIZE, "%d\n", size);
	strncat (str, tmp, strlen (tmp));

	for (S=nodes; NULL!=S; S=S->secondNode)
	{
		v = S->firstNode;
		snprintf (tmp, INPUT_SIZE, "%-10s ", v->label);
		strncat (str, tmp, strlen (tmp));
		for (T=nodes; NULL!=T; T=T->secondNode)
		{
			w = T->firstNode;
			d = D[v->index2][w->index2];
			if (d > threshold)
			{
				//snprintf (tmp, INPUT_SIZE, "NA          ");
				snprintf (tmp, 3, "NA");
				// Add [precision - 2] blank spaces
				for (i=0; i<precision-2; i++)
					strncat (tmp, " ", 1);
			}
			else
			{
				snprintf (tmp, INPUT_SIZE, (const char *)format, d);
			}
			strncat (str, tmp, strlen (tmp));
		}
		strcat (str, "\n");
	}
	strcat (str, "\n");
	free (tmp);

	return;
}

/*********************************************************/

void PrintTimeInfo (time_t t_beg, time_t t_end)
{
	ldiv_t hour, min;
	int sec;

	hour = ldiv (t_end - t_beg, 3600);
	min = ldiv (t_end - t_beg, 60);
	min.quot -= hour.quot * 60;
	sec = (int) (t_end-t_beg) % 60;

	if (min.quot > 9)
	{
		if (sec > 9)
			Message ( (char*)"Time used %dh%dm%ds", hour.quot, min.quot, sec);
		else
			Message ( (char*)"Time used %dh%dm0%ds", hour.quot, min.quot, sec);
	}
	else
	{
		if (sec > 9)
			Message ( (char*)"Time used %dh0%dm%ds", hour.quot, min.quot, sec);
		else
			Message ( (char*)"Time used %dh0%dm0%ds", hour.quot, min.quot, sec);
	}

	return;
}

/*********************************************************/

void PrintEstimatedMemorySpace (int nbtax, int nbsites, Options *options)
{
	unsigned long long int schar, sint, ssint, sdouble, sphydbl, nbT, nbS;
	unsigned long long int mem, memMat, memTree, memNNI, memSPR, Matsize, Treesize;
	lldiv_t ko;
	ldiv_t Mo;
	ldiv_t Go;
	
	memMat = memTree = memNNI = memSPR = 0;
	nbT = (unsigned long long int) nbtax;
	nbS = (unsigned long long int) nbsites;
	schar = (unsigned long long int) sizeof (char);
	sint = (unsigned long long int) sizeof (int);
	ssint = (unsigned long long int) sizeof (short int);
	sdouble = (unsigned long long int) sizeof (double);
	// sphydbl = (unsigned long long int) sizeof (phydbl);

	Matsize = nbT * nbT;
	Treesize = (nbT * 2) - 2;

	/**
	 ** Estimate memory space needed to compute distance matrix
	 */

	if (options->input_type == DNA)
	{
		memMat += (nbT * sdouble);
		memMat += (nbT * ((schar * nbS) + (10 * sint) + (4 * sdouble)));
		memMat += (Matsize * ( (7 * sdouble) + sint) );
		memMat /= 8;
	}
	else if (options->input_type == PROTEIN)
	{
		// Egein struct
		memMat += (4 * 20 * sdouble);
		memMat += (2 * 20 * sint);
		memMat += (4 * 20 * 20 * sdouble);

		// Pij Qmat
		memMat += (3 * 20 * 20 * sdouble);

		// Sequences
		memMat += ( nbT * MAX_NAME_LENGTH * schar );
		memMat += ( nbT * nbS * (schar + ssint) );
		memMat += ( nbS * 3 * ssint );
		memMat += ( nbS * sint );

		// Matrix
		memMat += ( nbT * MAX_NAME_LENGTH * schar );
		memMat += ( Matsize * 4 * sphydbl );
	}

#ifdef _OPENMP
	memMat = memMat * ( (unsigned long long int) options->nb_threads );
#endif

	// String length for trees
	memMat += ( (1 + (unsigned long long int) options->nb_bootstraps) * MAX_INPUT_SIZE * schar );

	/**
	 ** Estimate memory space needed to compute tree
	 */

	if (options->use_NNI)
	{
		// Size of 1 edge ~ 2 double + 2 int
		memNNI = (Treesize + 1) * 2 * 2 * sdouble;
		memNNI += ( (Treesize + 1) * 2 * sint );
	}
	
	if (options->use_SPR)
		memSPR = 2 * Treesize * Treesize * sdouble;
/*
	if (options->use_TBR)
	{
		//memTree += Treesize * Treesize * Treesize * sdouble;
		memTree = Treesize * Treesize * Treesize * sdouble;
	}
*/
	memTree = (memNNI > memSPR) ? memNNI : memSPR;
	mem = (memMat > memTree) ? memMat : memTree;

	ko = lldiv ( (long long) mem, 1024);
	Mo = ldiv (ko.quot, 1024);
	Go = ldiv (Mo.quot, 1024);

	if (Go.quot > 0)
		Warning ( (char*)"This analysis requires at least %d.%d Go of memory space.", Go.quot, (ldiv (Go.rem, 100)).quot );

	else if (Mo.quot > 100)
		Message ( (char*)"This analysis requires at least %d Mo of memory space.", Mo.quot);

	return;
}



// void InitSpeciesAndTrees (Options *options, set *species, char **bootTrees, char *bestTree)
// {
// 	int i = 0;
	
// 	species->firstNode = NULL;
// 	species->secondNode = NULL;

// 	size_t n = sizeof(bestTree);
// 	memset (bestTree, '\0', n);
	
// 	// mem alloc for bootstraped trees strings
// 	for (i=0; i<options->nb_bootstraps; i++)
// 	{
// 		bootTrees[i] = (char *) mCalloc (MAX_INPUT_SIZE, sizeof (char));
// 	}

// 	return;
// }

/*********************************************************/

tree *ComputeTree (Options *options, double **D, double **A, set *species,
	int numSpecies, int precision)
{
	set *slooper;
	node *addNode;
	tree *T = NULL;
	char format[8];
	snprintf (format, 8, "%%.%df", precision);
	
	switch (options->method)
	{
		case USER:
			break;

		case TaxAddBAL:
			for(slooper = species; NULL != slooper; slooper = slooper->secondNode)
			{
				addNode = copyNode (slooper->firstNode);
				T = BMEaddSpecies (T, addNode, D, A);
			}
			assignBMEWeights (T, A);
			break;
		
	}
	
	return T;
}

/*********************************************************/

tree *ImproveTree (Options *options, tree *T0, double **D, double **A,
	int *nniCount, int *sprCount, FILE *ofile)
{
	//T0 = ME
	//T1 = ME + NNI
	//T2 = ME + SPR
	//T3 = ME + NNI + SPR
	tree *T1, *T2, *T3;
	
	T1 = T2 = T3 = NULL;

	char format[8];
	snprintf (format, 8, "%%.%df", options->precision);
	
	// Print tree length
	weighTree (T0);
	if (verbose > 2)
		printf ("\n");
	Message ( (char*)"Tree length is %f.", T0->weight);
	fprintf (options->fpO_stat_file, "\tTree length is ");
	fprintf (options->fpO_stat_file, format, T0->weight);
	fprintf (options->fpO_stat_file, "\n");
	if (T0->weight < 0.)
		Warning ( (char*)"Negative tree length detected. Check your input data.");

	if (!options->use_NNI && !options->use_SPR)
	{
		switch (options->branch)
		{
			case BrBAL:
				if (options->method != TaxAddBAL)
					makeBMEAveragesTable (T0, D, A);

				assignBMEWeights (T0, A);
				break;

			default:
				break;
		}
	}
	
	if (options->use_NNI)
	{
		if (!isBoostrap)
		{
			if (verbose > 2)
				printf ("\n");

			Message ( (char*)"Performing NNI...");
		}
		
		T1 = copyTree (T0);
		
		switch (options->NNI)
		{
			case BALNNI:
				if (options->method != TaxAddBAL)
					makeBMEAveragesTable (T1, D, A);

				bNNI (T1, A, nniCount, options->fpO_stat_file, options->precision);
				assignBMEWeights (T1, A);
				break;

			default:
				break;
		}
		if (!isBoostrap)
			fprintf (ofile, "\tPerformed %d NNI(s).\n\n", *nniCount);
	}


	if (options->use_SPR)
	{
		if (!isBoostrap)
		{
			if (verbose > 2)
				printf ("\n");

			Message ( (char*)"Performing SPR...");
		}
		
		T2 = copyTree (T0);
		
		makeBMEAveragesTable (T2, D, A);
		SPR (T2, D, A, sprCount, options->fpO_stat_file, options->precision);
		assignBMEWeights (T2, A);
		
		if (!isBoostrap)
			fprintf (ofile, "\tPerformed %d SPR(s).\n\n", *sprCount);
	}
	
	if (NULL != T1)
		weighTree (T1);

	if (NULL != T2)
		weighTree (T2);

	if (NULL != T1)
	{
		if (T0->weight > T1->weight)
		{
			freeTree (T0);
			T0 = T1;	//using T0 as the place to store the minimum evolution tree in all cases
			T1 = NULL;
		}
		else
			Message ( (char*)"NNI did not improve tree length.");
	}
	else if (NULL != T1)
		freeTree (T1);

	if (NULL != T2)
	{
		if (T0->weight > T2->weight)
		{
			freeTree (T0);
			T0 = T2;
			T2 = NULL;
		}
		else
			Message ( (char*)"SPR did not improve tree length.");
	}
	else if (NULL != T2)
		freeTree (T2);
	
	return T0;
}

/*********************************************************/

/* Pij : patristic distances computed from output tree
 * Dij : distances from the input matrix
 * Dm  : mean distances from input matrix
 *              ___
 *             \      _             _  2
 *              \    /   P   _  D    \
 *              /    \_   ij     ij _/
 *             /___   
 * V =  1 -  ------------------------------
 *              ___
 *             \      _             _  2
 *              \    /   D   _  D    \
 *              /    \_   ij     m  _/
 *             /___   
 * 
 * D : the input distances matrix
 * T : the tree from which the patristic distances are computed
 * n : number of sequences*/
void explainedVariance (double **D, tree *T, int n, int precision,
	int input_type, FILE *fpO)
{
	int i, j;
	double var, num, dem, Dm, **P;
	char format[7];

	snprintf (format, 7, "%%.%df\n", precision);
	
	// Mean of input distances
	Dm = meanDist (D, n);
	
	// Patristic distances matrix
	P = patristicMatrix (T, n, FALSE);
	
	// Compute explained variance
	num = dem = 0.0;
	for (i=0; i<n; i++)
	{
		for (j=0; j<n; j++)
		{
			num += (P[i][j] - D[i][j]) * (P[i][j] - D[i][j]);
			dem += (D[i][j] - Dm) * (D[i][j] - Dm);
		}
	}
	var = 1.0 - ( num / dem );
	
	freeMatrix (P, n);
	
	fprintf (fpO, "\tExplained variance : "); 
	fprintf (fpO, format, var);
	
	if (input_type == PROTEIN)
	{
		if (var < EXPL_VAR_AA)
			Warning ( (char*)"\n\tExplained variance = %.3f (<%.2f).\n\tCheck your input data.", var, EXPL_VAR_AA);
	}
	else
	{
		if (var < EXPL_VAR_NT)
			Warning ( (char*)"\n\tExplained variance = %.3f (<%.2f).\n\tCheck your input data.", var, EXPL_VAR_NT);
	}
	
	return;
}

