#include <iostream>
#include "fastme.h"
int run(double**d, int n_taxa, int** init_adj)
{
	set *species, *species_bk;
	
	species = species_bk = nullptr;
	
	time_t t_beg, t_end;

	tree *T = nullptr;

	int** sparse_A;
	sparse_A = get_sparse_A(init_adj, n_taxa*2 - 2);

	printf("\n%i size tree\n", n_taxa);
	
	for(int i = 0; i<n_taxa*2 - 2; i++){

		printf("%i   ", i);
		for(int j=0; j<3; j++) printf("%i ", sparse_A[i][j]);
		printf("\n");
	}

	T = adj_to_tree(sparse_A, n_taxa);

	printTree(T);



#ifdef _OPENMP
	int i = 0;
#endif

	int numSpecies = n_taxa;
	int setCounter = 0;
	
	int nniCount, sprCount, repCounter, repToPrint, printedRep;
	
	int seqLength = -1;
	
	// Random numbers for bootstraps
	int **rnd = nullptr;

	double **D, **A;
	D = A = nullptr;

	// Strings containing the best tree and the bootstraped trees
	char *bestTree = nullptr;
	char **bootTrees = nullptr;
	
	// String to record pseudo-matrices
	char **matStr = nullptr;


// #ifdef _OPENMP
// 	Message ( (char*)"This analysis will run on %d threads.", options->nb_threads);
// #else

// #endif


		//nniCount = sprCount = tbrCount = repCounter = 0;
		nniCount = sprCount = repCounter = 0;
		setCounter++;
		species = new set;
		species->firstNode = nullptr;
		species->secondNode = nullptr;

/*********************************************************
		GET MATRIX
**********************************************************/
		make_set(numSpecies, species);

		// checkTrgInequality (D, numSpecies, options->trg_ineq);

		species_bk = copySet (species);
		

/*********************************************************
		COMPUTE TREE
**********************************************************/


			A = initDoubleMatrix (2*numSpecies-2);
			

			// T = ComputeTree (options, D, A, species, numSpecies, options->precision);
			// printf("hhhh");
			// compareSets (T, species);

			

			T = ImproveTree (T, d, A, &nniCount, &sprCount);
			printTree(T);
			printf("********"); // 0.11784069   0.14665112
			int** adj_mat;
			adj_mat = tree_to_adj_mat(T);
			
			
		
			// explainedVariance (D, T, numSpecies, options->precision, options->input_type, options->fpO_stat_file);
		std::cout<<"Performed "<< nniCount<< "NNI(s) and "<<sprCount<<" SPR(s) "<<std::endl;


#ifndef _OPENMP
		freeMatrix (A, 2*numSpecies-2);
#endif
		
		if (nullptr != matStr)
			free (matStr);
		
		freeMatrix (D, numSpecies);
		freeSet (species);
		freeSet (species_bk);
		


	
	if (nullptr != rnd)		free (rnd);
	if (nullptr != bootTrees)	free (bootTrees);
	if (nullptr != bestTree)	free (bestTree);


	//exit (EXIT_SUCCESS);

	return 0;
}


tree *ImproveTree (tree *T0, double **D, double **A,
	int *nniCount, int *sprCount)
{
	//T0 = ME
	//T1 = ME + NNI
	//T2 = ME + SPR
	//T3 = ME + NNI + SPR
	tree *T1, *T2, *T3;
	
	T1 = T2 = T3 = nullptr;
	
	// Print tree length
	weighTree (T0);

	std::cout<<"Tree length is "<< T0->weight<<std::endl;


	std::cout<<"NNI"<<std::endl;
		
    T1 = copyTree (T0);
		
    makeBMEAveragesTable (T1, D, A);

    bNNI (T1, A, nniCount);
    assignBMEWeights (T1, A);
		
    std::cout<<"Performed NNI(s) "<<*nniCount<<std::endl;
	


    std::cout<<"SPR"<<std::endl;
		
    T2 = copyTree (T0);
		
    makeBMEAveragesTable (T2, D, A);
    SPR (T2, D, A, sprCount);
    assignBMEWeights (T2, A);
		
	std::cout<<"Performed SPR(s) "<<*sprCount<<std::endl;

	
	if (nullptr != T1)
		weighTree (T1);

	if (nullptr != T2)
		weighTree (T2);

	if (nullptr != T1)
	{
		if (T0->weight > T1->weight)
		{
			freeTree (T0);
			T0 = T1;	//using T0 as the place to store the minimum evolution tree in all cases
			T1 = nullptr;
		}
	}
	else if (nullptr != T1)
		freeTree (T1);

	if (nullptr != T2)
	{
		if (T0->weight > T2->weight)
		{
			freeTree (T0);
			T0 = T2;
			T2 = nullptr;
		}

	}
	else if (nullptr != T2)
		freeTree (T2);
	
	return T0;
}