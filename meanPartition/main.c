/* MigSelect 2013-2019 Arun Sethuraman, Vitor Sousa, Jody Hey */
/* This code was developed based on IMa2  2009-2010  Jody Hey, Rasmus Nielsen and Sang Chul Choi*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "DPART.h"

// MEAN PARTITION returns the mean partition and coassignment probabilities
// usage is 
// MeanPartition inputfilename numberGroups numberloci numbersamples(AssignmentVectorsSampledMCMC)



int main (int argc, char *argv[])
{

	FILE *outputfile; // output file with the meanpartition and coassignment probabilities
	FILE *inputfile; // input file with the partition
	int *sampledassign; // matrix with the sampled assigned vectors
	int *meanpartition; // array where the meanpartition is saved
	double *coassign; // array to save the coassignment probabilities (size number of pairs of loci)
	double *probcorrectassign; // array to save the probability of a locus being assigned to the cluster identified in the mean partition
	int npairs; // numbers of pairs of loci
	int ngroups; // number of groups
	int nloci; // number of loci
	int nsamples; // sample size of assignment vectors
	int i; 
	char *loadfilebase; // string with the name of the input file



	// READ INPUT ARGS COMMAND
	loadfilebase = argv[1]; // loadfilebase is given as input
	ngroups=atoi(argv[2]); // number of groups
	nloci=atoi(argv[3]); // number of loci
	nsamples=atoi(argv[4]); // number of samples
	
	// allocate memory for sampledassign
	// instead of assigning a matrix, it assigns an array
	sampledassign = (int *) malloc((nloci*nsamples)*sizeof(int));

	// allocate memory for meanpartition
	meanpartition = (int *) malloc(nloci*sizeof(int));

	// allocate memory for the coassignment probabilities
	npairs = nloci*(nloci-1)/2;
	coassign = (double *) malloc(npairs*sizeof(double));

	// allocate memory for the probability of correct assignment
	probcorrectassign = (double *) malloc(nloci*sizeof(double));

	// open the file to read the assignments sampled in the MCMC
	inputfile=fopen(loadfilebase, "r");
	

	// go through the file and save it into 
	for(i=0; i<(nloci*nsamples); i++) {
		fscanf(inputfile, "%i", &sampledassign[i]);
	}

	// Close the inputfile
	fclose(inputfile);

	// Compute the mean partition
	Meanpartdis(sampledassign, nloci, nsamples, meanpartition);

	// Print the mean partition
	printf("Mean Partition:\n");
	for(i=0; i<nloci; i++) {
		printf("%i ", meanpartition[i]);
	}

	// Get the individual probability for each locus
	indprob(meanpartition, sampledassign, nloci, nsamples, probcorrectassign );

	// Print the probability of correct assignmnet
	printf("\nProbability of correct assignments:\n");
	for(i=0; i<nloci; i++) {
		printf("%g ", probcorrectassign[i]);
	}

	// Initialize coassign to zero
	for(i=0; i<npairs; i++) {
		coassign[i]=0;
	}

	// Compute the coassignment probabilities
	Coassignprobs(sampledassign, nloci, nsamples, coassign);

	// Print the coassignment probabilities
	printf("\nCoassignment probabilities:\n");
	for(i=0; i<npairs; i++) {
		printf("%g ", coassign[i]);
	}

	

	// FREE POINTERS
	free(sampledassign);
	free(meanpartition);
	free(coassign);
	free(probcorrectassign);

}
