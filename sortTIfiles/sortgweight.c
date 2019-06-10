/* MigSelect 2013-2019 Arun Sethuraman, Vitor Sousa, Jody Hey */
/* This code was developed based on IMa2  2009-2010  Jody Hey, Rasmus Nielsen and Sang Chul Choi*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

// SORTSAMPLEDASSIGN
// Function to sort the assignment vector
// Returns an array of an indicator variable I that takes the value 0 or 1
// If I[i]=0, the sampled vector z[i] has the 
//	minimum partition distance to the mean assignment.
// If I[1]=1, the sampled vector z[i] needs to be relabeled 
//	to have the minimum partition distance to the mean assignment.
// This variable determines if the gweights need to be swapped or not.
// INPUT:
//	int *index : pointer to the indicator array with values 0 or 1
//	int *sampledz : pointer to matrix with the sampled assignment vectors (vectors have labels 0, 1)
//	int *meanassign : pointer to the array with the mean assignment (mean partiton has labels 1,2)
//	int ni : number of elements in each vector (i.e. in this case it will be the number of loci)
//	int ns : number of sampled vectors
// Function coded by VS
// NOTE: ONLY WORKS WITH 2 CLUSTERS, OR 2 GROUPS
// When there are more than 2 groups, there are many possible relabellings
// and hence, we would need to return the vector with the required relabelling.
void sortsampledassign(int *index, int* sampledz, int* meanassign, int ni, int ns) {

	int *tmpmeanpart; // temporary array where the meanpartition is saved, changing labels to 0,1,...
	int *tmpassign1; // temporary array where the a given sampled assignment vector is saved, applying restriction growth function
					   // in practice we check if the labels in the meanpartition correspond to the labels in the sampled assignment vector
	int *tmpassign2; // temporary array where the a given sampled assignment vector is saved, applying restriction growth function
					 // in practice we check if the labels in the meanpartition correspond to the labels in the sampled assignment vector
					 // with alternative label
	int *countdiff; // array where the number of times a given individual 
						  // in a given assignment has to be removed when 
						  // computing the partition distance is not in the mean partition
	int ss, ii; // auxiliary index for loops
	
	// partition distance
	int partdist1; // partition distance between assignment and label 1
	int partdist2; // partition distance between assignment and label 2

	// Allocate memory for the temporary arrays
	tmpmeanpart = (int *) malloc(ni * sizeof(int));
	tmpassign1 = (int *) malloc(ni * sizeof(int));
	tmpassign2 = (int *) malloc(ni * sizeof(int));
	countdiff = (int *) malloc(ni * sizeof(int));

	// change labels of mean partition from 1,2 to 0,1
	for(ii=0; ii<ni; ii++) {
		tmpmeanpart[ii] = meanassign[ii]-1;
	}

	// initialize the tmpcountdiff to zero
	for(ii=0; ii<ni; ii++) {
		countdiff[ii] = 0;
	}

	// for loop from first sampled to last assignment vector
	// compare each sampled assignment vector with the meanpartition vector
	for ( ss=0; ss<ns; ss++ )
	{
		// Copy the sampled assignment vector to the temporary vector
		for ( ii=0; ii<ni; ii++ ) {
			tmpassign1[ii] = sampledz[ss*ni + ii];
			// Get the assignment array with alternative labelling
			if(tmpassign1[ii]==0) {
				tmpassign2[ii]=1; // labels 0 become 1
			}
			else {
				tmpassign2[ii]=0; // labels 1 become 0
			}
		}

		// Compute the partition distance between the assignment vectors
		// and the two possible labellings
		
		// Set the partition distances to zero
		partdist1 = 0; partdist2 = 0;
		for ( ii=0; ii<ni; ii++ )
		{
			if(tmpassign1[ii]!=tmpmeanpart[ii]) {
				partdist1++;
			}	
			if(tmpassign2[ii]!=tmpmeanpart[ii]) {
				partdist2++;
			}
		}	

		// Check what is the minimum distance
		if(partdist1<partdist2) {
			// If partdist 1 is lower, it means the indicator variable takes the value zero
			index[ss]=0;
		}
		else {
			// If partdist 1 is lower, it means the indicator variable takes the value zero
			index[ss]=1;
		}		
	}

	// Free the temporary arrays
	free(tmpmeanpart);
	free(tmpassign1);
	free(tmpassign2);
	free(countdiff);
}
