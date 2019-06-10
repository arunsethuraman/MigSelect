/* MigSelect 2013-2019 Arun Sethuraman, Vitor Sousa, Jody Hey */
/* This code was developed based on IMa2  2009-2010  Jody Hey, Rasmus Nielsen and Sang Chul Choi*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "DPART.h"

/* calculation of the mean partition (proposed by Huelsenbeck et al. 2007) ------------------------------------------------------------*/

   /*Cost matrix is generated using eq.16 in Konovalov et al. 2005 and solved by the hungarian algorithm.
   Original C++ source code of the hungarian algorithm is provided by 
           http://www.topcoder.com/tc?module=Static&d1=tutorials&d2=hungarianAlgorithm 
   Here it is converted to C */

//#define Nv 50            /* max number of vertices in one part */
//#define INF 100000000   /* just infinity */
//#define whichmax(a, b) ((a) > (b) ? (a) : (b))
//#define whichmin(a, b) ((a) < (b) ? (a) : (b))

int cost[Nv][Nv];          /* cost matrix */
int n, max_match;        /* n workers and n jobs */
int lx[Nv], ly[Nv];        /* labels of X and Y parts */
int xy[Nv];               /* xy[x] - vertex that is matched with x, */
int yx[Nv];               /* yx[y] - vertex that is matched with y */
int Sh[Nv], Th[Nv];        /* sets Sh and Th in algorithm */
int slack[Nv];            /* as in the algorithm description */
int slackx[Nv];           /* slackx[y] such a vertex, that */
                         /* l(slackx[y]) + l(y) - w(slackx[y],y) = slack[y] */
int prev[Nv];             /* array for memorizing alternating paths */

/*----------------------------------------------------------------------------------------------------*/
void add_to_tree(int x, int prevx) 
/* x - current vertex,prevx - vertex from X before x in the alternating path, */
/* so we add edges (prevx, xy[x]), (xy[x], x) */
{
	int y;

    Sh[x] = 1;                      /* add x to Sh */
    prev[x] = prevx;                /* we need this when augmenting */
    for (y = 0; y < n; y++)         /* update slacks, because we add new vertex to Sh */
        if (lx[x] + ly[y] - cost[x][y] < slack[y])
        {
            slack[y] = lx[x] + ly[y] - cost[x][y];
            slackx[y] = x;
        }
}

/*----------------------------------------------------------------------------------------------------*/
void augment(void)                         /* main function of the algorithm */
{

    int x, y, root;                        /* just counters and root vertex */
    int q[Nv], wr = 0, rd = 0;              /* q - queue for bfs, wr,rd - write and read */
	int cx,cy,ty;
	int i, j;
	int delta = INF;

    if (max_match == n) return;            /* check wether matching is already perfect */
                                           /* pos in queue */
    memset(Sh, 0, sizeof(Sh));             /* init set Sh */
    memset(Th, 0, sizeof(Th));             /* init set Th */
    memset(prev, -1, sizeof(prev));        /* init set prev - for the alternating tree */


    for (x = 0; x < n; x++)                /* finding root of the tree */
        if (xy[x] == -1)
        {
            q[wr++] = root = x;
            prev[x] = -2;
            Sh[x] = 1;
            break;
        }

    for (y = 0; y < n; y++)                /* initializing slack array */
    {
        slack[y] = lx[root] + ly[y] - cost[root][y];
        slackx[y] = root;
    }

    /* second part of augment() function */
    while (1)                                                           /* main cycle */
    {
        while (rd < wr)                                                 /* building tree with bfs cycle */
        {
            x = q[rd++];                                                /* current vertex from X part */
            for (y = 0; y < n; y++)                                     /* iterate through all edges in equality graph */
                if (cost[x][y] == lx[x] + ly[y] &&  !Th[y])
                {
                    if (yx[y] == -1) break;                             /* an exposed vertex in Y found, so */
                                                                        /* augmenting path exists! */
                    Th[y] = 1;                                          /* else just add y to Th, */
                    q[wr++] = yx[y];                                    /* add vertex yx[y], which is matched */
                                                                        /* with y, to the queue */
                    add_to_tree(yx[y], x);                              /* add edges (x,y) and (y,yx[y]) to the tree */
                }
            if (y < n) break;                                           /* augmenting path found! */
        }
        if (y < n) break;                                               /* augmenting path found! */

         /* augmenting path not found, so improve labeling. "update_labels" */
		for (j = 0; j < n; j++)                                         /* calculate delta using slack */
        if (!Th[j])
            delta = whichmin(delta, slack[j]);
        for (i = 0; i < n; i++)                                         /* update X labels */
            if (Sh[i]) lx[i] -= delta;
        for (j = 0; j < n; j++)                                         /* update Y labels */
            if (Th[j]) ly[j] += delta; 
        for (j = 0; j < n; j++)                                         /* update slack array */
            if (!Th[j])
                slack[j] -= delta;

        wr = rd = 0;                
        for (y = 0; y < n; y++)        
        /* in this cycle we add edges that were added to the equality graph as a
        result of improving the labeling, we add edge (slackx[y], y) to the tree if
        and only if !Th[y] &&  slack[y] == 0, also with this edge we add another one
        (y, yx[y]) or augment the matching, if y was exposed */
            if (!Th[y] &&  slack[y] == 0)
            {
                if (yx[y] == -1)                                        /* exposed vertex in Y found - augmenting path exists! */
                {
                    x = slackx[y];
                    break;
                }
                else
                {
                    Th[y] = 1;                                          /* else just add y to Th, */
                    if (!Sh[yx[y]])    
                    {
                        q[wr++] = yx[y];                                /* add vertex yx[y], which is matched with */
                                                                        /* y, to the queue */
                        add_to_tree(yx[y], slackx[y]);                  /* and add edges (x,y) and (y, */
                                                                        /* yx[y]) to the tree */
                    }
                }
            }
        if (y < n) break;                                               /* augmenting path found! */
    }

    if (y < n)                                                          /* we found augmenting path! */
    {
        max_match++;                                                    /* increment matching */
        /* in this cycle we inverse edges along augmenting path */

        for (cx = x, cy = y; cx != -2; cx = prev[cx], cy = ty)
        {
            ty = xy[cx];
            yx[cy] = cx;
            xy[cx] = cy;
        }

        augment();                                                      /* recall function, go to step 1 of the algorithm */
    }
}


/*----------------------------------------------------------------------------------------------------*/
int hungarian(void)
{
    int ret = 0;                      /* weight of the optimal matching */
    int x, y;

    max_match = 0;                    /* number of vertices in current matching */
    memset(xy, -1, sizeof(xy));    
    memset(yx, -1, sizeof(yx));

    memset(lx, 0, sizeof(lx));        /* step 0, "init_labels" */
    memset(ly, 0, sizeof(ly));

    for (x = 0; x < n; x++)
        for (y = 0; y < n; y++)
            lx[x] = whichmax(lx[x], cost[x][y]);                    

    augment();                        /* steps 1-3 */
    for (x = 0; x < n; x++)       /* forming answer there */
        ret += cost[x][xy[x]];

    return (ret);
}


/*----Calculate the partition distance between A and B which has Ne elements--------------------------*/
int Partdis ( const int Ne, const int *A, const int *B )
{
	int ii, jj, kk, Pos, sum, Dis;
	int *Abs=NULL, *Bbsbar=NULL;

	/* n represents the number of clusters in A and B */
	for (ii=0, n=0; ii<Ne; ii++)
	{
		if ( n < (int) A[ii] )
			n = (int) A[ii];
		if ( n < (int) B[ii] )
			n = (int) B[ii];
	}
	n += 1; /* because the cluster number in A and B starts with 0 */

	if ( n > Nv )
	{
		printf ( "Number of clusters in sampled partitions is larger than %d\n", Nv);
		return (0);
	}

	Abs = calloc ( Ne * n, sizeof ( int ) );    if(Abs==NULL) alert();
	Bbsbar = calloc ( Ne * n, sizeof ( int ) ); if(Bbsbar==NULL) alert();

	for (ii=0; ii<n; ii++)
		for (jj=0; jj<Ne; jj++ )
		{
			Pos = ii * Ne + jj;
			if ( A[jj] == (int) ii ) Abs [ Pos ] = 1;
			if ( B[jj] != (int) ii ) Bbsbar [ Pos ] = 1;
		}

	for (ii=0; ii<n; ii++ )
	{
		for (jj=0; jj<n; jj++ )
		{
			sum=0;
			for ( kk=0; kk<Ne; kk++ )
				sum += ( Abs [ jj*Ne+kk ] * Bbsbar [ ii*Ne+kk ] );

			cost [ii][jj] = -sum;
		}
	}

	Dis = hungarian ();
	Dis *= -1;

	free(Abs); Abs=NULL; free(Bbsbar);  Bbsbar=NULL;

	return(Dis);

}

/*----Main ( calculate the mean partition distance )--------------------------------------------------*/
// Function that calculates and returns the mean partition
// INPUT:
//	int *sampledz : pointer to array of int with the sampled assignment vectors
//	int ne : ne is the number of elements of vector, i.e. number of loci
//	int ns : ns is the number of sampled assignment vectors
//	int *meanpartition : pointer to array of int where meantpartition is saved

void Meanpartdis ( int *sampledz, int ne, int ns, int *meanpartition)
{
	int ii, jj;
	int *Meanp, *propMeanp, kmeanp=0, kk, min;
	int Distsum, propDistsum=0, ele=0, Move=0, Notmove=0;
	int poplabel, labelupdate;

	// VS
	// ns is the number of sampled assignment vectors
	// Note that in this case the assignment vector is treated has
	// a single array, instead of being seen as a 2D array (i.e. a matrix)
	Distsum = ne * ns;

	// ne is the number of elements
	// meanp has size ne
	Meanp = malloc ( sizeof ( int ) * ne); if(Meanp==NULL) alert();
	// propMeanp also has size ne
	propMeanp = malloc ( sizeof ( int ) * ne); if(propMeanp==NULL) alert();
	
	// is this initializing the meanp with the firt sampledz
	for (ii=0; ii<ne; ii++ ){ Meanp[ii] = sampledz [ii];}

	
	do
	{
		kmeanp = 0;
		for (ii=0; ii<ne; ii++ )
		{
			if ( Meanp[ii] > kmeanp )
				kmeanp = Meanp [ii];
			propMeanp [ii] = Meanp [ii];
		}
		kmeanp += 1;	

		for (kk=0; kk<=kmeanp; kk++ )
		{
			propMeanp[ele] = kk;
			
			propDistsum = 0;
			for (jj=0; jj<ns; jj++ )
			{
				propDistsum += Partdis (ne, propMeanp, sampledz + ( jj * ne ) );

			}
			
			if ( propDistsum < Distsum )
			{
				Distsum = propDistsum;
				for (jj=0; jj<ne; jj++ ) { Meanp[jj] = propMeanp[jj]; }
				Move += 1;
				Notmove = 0;
			}
		}

		ele += 1;

		if ( Move == 0 ) { Notmove += 1;} else { Move = 0; }
		if ( ele == ne )
			ele = 0;

	} while ( Notmove<( ne - 1) );


	min = Meanp [0];
	for ( ii=1; ii<ne; ii++ )
		if ( Meanp [ii] < min )
			min = Meanp [ii];

	jj=0; poplabel=1;
	do
	{
		labelupdate = 0;
		for ( ii=0; ii<ne; ii++ )
			if ( Meanp [ii] == min )
			{
				meanpartition [ii] = poplabel;
				jj ++;
				labelupdate ++;
			}

		min ++;

		if (labelupdate>0)
			poplabel ++;

	} while ( jj<ne );


	
	free(propMeanp); free(Meanp);

}


/* Calculate the co-assignment probs ---------------------------------------------------------------------------------------------------*/
// Function that calculates and returns the mean partition
// INPUT:
//	int *sampledz : pointer to array of int with the sampled assignment vectors
//	int ni : ne is the number of elements of vector, i.e. number of loci
//	int ns : ns is the number of sampled assignment vectors
//	int *coaprob : pointer to array with coassignment probabilities (size - number of possible pair comparisons)
// PRE-REQUISITE: coaprob must be initialized to zero
void Coassignprobs ( int *sampledz, int ni, int ns, double *coaprob )
{
	int ii, jj, ss, rr;

    for ( ss=0; ss<ns; ss++ )
	{
		rr = 0;
		for ( ii=0; ii<(ni-1); ii++ )
		{
			for ( jj=ii+1; jj<ni; jj++ )
			{
				if ( sampledz [ ss * ni + ii ] == sampledz [ ss * ni + jj ] )
					coaprob [rr] += 1.0;
				rr ++;
			}
		}
	}

	for ( ii=0; ii<rr; ii++ )
		coaprob [ii] /= (double) ns;

}

/*----alert when allocation of memory is failed-----------------------------------------------------------------------------------*/
void alert ( void )
{
	puts ( "memory allocation failure. try again" );
}


// INDPROB
// returns the uncertainty probability of each individual being assigned
// to the cluster identified in the meanpartition.
// This is computed as the proportion of times
// the individual is not removed when computing the partition distance
// between a given assignment and the mean partition.
// INPUT:
//	int *meanpart: pointer to the mean partition array (note that meanpart is labelled as 1,2,...)
//	int *sampledz: pointer to the matrix with the sampled assignment vectors (sampled vectors are labbeled as 0,1,...)
//	int ni : ni (number of individuals) is the number of elements of vector, i.e. number of loci
//	int ns : ns (number of samples)is the number of sampled assignment vectors
//	double *returnprob : pointer to the array where the probability is saved
// Function coded by VS
// NOTE: ONLY WORKS WITH 2 CLUSTERS, OR 2 GROUPS
void indprob (int *meanpart, int *sampledz, int ni, int ns, double *returnprob ) {

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
	double probwrong; // probability of wrong assignment, this is the same as countdiff/ns
	
	// array with the partition distance
	int partdist1; // partition distance between assignment and label 1
	int partdist2; // partition distance between assignment and label 2

	// Allocate memory for the temporary arrays
	tmpmeanpart = (int *) malloc(ni * sizeof(int));
	tmpassign1 = (int *) malloc(ni * sizeof(int));
	tmpassign2 = (int *) malloc(ni * sizeof(int));
	countdiff = (int *) malloc(ni * sizeof(int));

	// change labels of mean partition from 1,2 to 0,1
	for(ii=0; ii<ni; ii++) {
		tmpmeanpart[ii] = meanpart[ii]-1;
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
			// Look at the difference between the assignment vector 1
			// and meanpartition
			for ( ii=0; ii<ni; ii++ )
			{
				if(tmpassign1[ii]!=tmpmeanpart[ii]) {
					countdiff[ii]++;
				}	
			}
		}
		else {
			// Look at the difference between the assignment vector 2
			// and meanpartition
			for ( ii=0; ii<ni; ii++ )
			{
				if(tmpassign2[ii]!=tmpmeanpart[ii]) {
					countdiff[ii]++;
				}	
			}
		}		
	}

	// the probability of being assigned to the same cluster as the one
	// in the meanpartition is given by
	// 1 - prob(wrong assignment)
	// the probability of wrong assignment is the frequency of times
	// the individual (or locus) was different from the meanpartition
	// that is countdiff/ns
	for ( ii=0; ii<ni; ii++ )
	{
		probwrong=(double) countdiff[ii]/ (double)ns;
		returnprob[ii] = 1-probwrong;
	}

	// Free the temporary arrays
	free(tmpmeanpart);
	free(tmpassign1);
	free(tmpassign2);
	free(countdiff);
}



// INDPROB_ATTEMPT TO HAVE A GENERAL FUNCTION FOR MORE THAN 2 GROUPS
// I am stuck at the part where I need to list and obtain all possible relabelling
// returns the uncertainty probability of each individual being assigned
// to the cluster identified in the meanpartition.
// This is computed as the proportion of times
// the individual is not removed when computing the partition distance
// between a given assignment and the mean partition.
// INPUT:
//	int *meanpart: pointer to the mean partition array (note that meanpart is labelled as 1,2,...)
//	int *sampledz: pointer to the matrix with the sampled assignment vectors (sampled vectors are labbeled as 0,1,...)
//	int ni : ni (number of individuals) is the number of elements of vector, i.e. number of loci
//	int ns : ns (number of samples)is the number of sampled assignment vectors
//	int ngroups : number of groups
//	double *returnprob : pointer to the array where the probability is saved
// Function coded by VS
// NOTE: ONLY WORKS WITH 2 CLUSTERS, OR 2 GROUPS
/*void indprob (int *meanpart, int *sampledz, int ni, int ns, int ngroups, double *returnprob ) {

	int *tmpmeanpart; // temporary array where the meanpartition is saved, changing labels to 0,1,...
	int **tmpassign; // temporary array where the a given sampled assignment vector is saved, applying restriction growth function
					 // and considering all the possible relabelings
	int *countdiff; // array where the number of times a given individual 
				    // in a given assignment has to be removed when 
					// computing the partition distance is not in the mean partition
	int *partdist; // array with the partition distance of each relabbeling
	int ss, ii, ll; // auxiliary index for loops
	int nrelabel; // number of possible relabelling
	double probwrong; // probability of wrong assignment, this is the same as countdiff/ns
	
	// Get the total number of relabelling
	nrelabel = factorial(ngroups);


	// Allocate memory for the temporary arrays
	tmpmeanpart = (int *) malloc(ni * sizeof(int));
	partdist = (int *) malloc(ngroups * sizeof(int));
	countdiff = (int *) malloc(ni * sizeof(int));

	// tmp assign is a 2D matrix of integers
	// this matrix saves the assignment under different possible relabellings
	tmpassign = (int **) malloc(nrelabel * sizeof(int *));
	for(ii=0; i<nrelabel; i++) {
		tmpassign[ii] = (int *) malloc(ni * sizeof(int));
	}

	// change labels of mean partition from 1,2 to 0,1
	for(ii=0; ii<ni; ii++) {
		tmpmeanpart[ii] = meanpart[ii]-1;
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
			tmpassign[0][ii] = sampledz[ss*ni + ii];
			// Get the assignment array with alternative labelling
			for(ll=0; ll<nrelabel; ll++) {
				
				

				if(tmpassign[0][ii]==0) {
					tmpassign2[ii]=1; // labels 0 become 1
				}
				else {
					tmpassign2[ii]=0; // labels 1 become 0
				}
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
			// Look at the difference between the assignment vector 1
			// and meanpartition
			for ( ii=0; ii<ni; ii++ )
			{
				if(tmpassign1[ii]!=tmpmeanpart[ii]) {
					countdiff[ii]++;
				}	
			}
		}
		else {
			// Look at the difference between the assignment vector 2
			// and meanpartition
			for ( ii=0; ii<ni; ii++ )
			{
				if(tmpassign2[ii]!=tmpmeanpart[ii]) {
					countdiff[ii]++;
				}	
			}
		}		
	}

	// the probability of being assigned to the same cluster as the one
	// in the meanpartition is given by
	// 1 - prob(wrong assignment)
	// the probability of wrong assignment is the frequency of times
	// the individual (or locus) was different from the meanpartition
	// that is countdiff/ns
	for ( ii=0; ii<ni; ii++ )
	{
		probwrong=(double) countdiff[ii]/ (double)ns;
		returnprob[ii] = 1-probwrong;
	}

	// Free the temporary arrays
	free(tmpmeanpart);
	free(tmpassign1);
	free(tmpassign2);
	free(countdiff);
}*/

// FACTORIAL
// returns the factorial of an integer smaller than 13
// INPUT:
//	int n : integer smaller than 13
// PRE-REQUISITES
// n < 13 and n >0
int factorial (int n)
{
  //assert (!(n < 0));
  //assert (n < 13);
  if (n == 0)
    {
      return 1;
    }
  if (n == 1)
    {
      return 1;
    }
  if (n == 2)
    {
      return 2;
    }
  else
    {
      return n * factorial (n - 1);
    }
}
