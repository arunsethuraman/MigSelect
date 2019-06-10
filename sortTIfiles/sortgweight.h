/* MigSelect 2013-2019 Arun Sethuraman, Vitor Sousa, Jody Hey */
/* This code was developed based on IMa2  2009-2010  Jody Hey, Rasmus Nielsen and Sang Chul Choi*/
#define XFREE(p) do { free((p)); (p) = NULL; } while(0)

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
void sortsampledassign(int *index, int* z, int* meanassign, int ni, int ns);

