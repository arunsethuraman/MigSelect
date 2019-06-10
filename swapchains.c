/* MigSelect 2013-2019 Arun Sethuraman, Vitor Sousa, Jody Hey */
/* This code was developed based on IMa2  2009-2010  Jody Hey, Rasmus Nielsen and Sang Chul Choi*/

#undef GLOBVARS
#include "imamp.h"

/*********** LOCAL STUFF **********/
//static unsigned long swapcount[MAXCHAINS][MAXCHAINS];

static double swapweight (int ci, int cj);
//AS: adding a function to calculate only the partial weight on a chain - for swapping
static double calcpartialswapweight(int c);
//AS: swapweight_bwprocesses is a function that attempts/performs swaps between a pair of processors
static double swapweight_bwprocesses(double sumi, double sumj, double betai, double betaj, int n1u, int n1);
static int getnu(int c);

double
swapweight (int ci, int cj)
{
  int li;
  double sumi = 0, sumj = 0, w;
  // VS added these variables to get the prior ratio of assignment
  int i, n1=0, n1u=0;
  double priora=0.0;

  	for(i=0; i<nloci; i++) {
		// assignment a
		n1 += grouploci_mig[ci][i];
		// assignment a'
		n1u += grouploci_mig[cj][i];
	}
	// DPP prior
	// If the update does not change the number of groups k
	/*if((n1>0 && n1<nloci) && (n1u>0 && n1u<nloci)) {
		// Get the prior ratio
		if(n1u>n1) { // if n1'=n1+1
			priora = log(n1)-log(nloci-n1u);
		}
		else {
			priora = log(nloci-n1)-log(n1u);
		}
	}
	// If the update changes the number of groups
	else {
		// this means that the current vector is all zeros
		// and that updated vector n1' has 1 non zero element
		if(n1==0 || n1==nloci) {
			priora = -log((nloci-1));
		}
		// if n1=nloci-1 and n1u=nloci
		if(n1u==nloci || n1u==0) {
			priora = log(nloci-1);
		}
	}*/
	// Equal probability prior (see int updategrouplociassign_mig (int ci) for details)
	if(n1u>n1) { // if n1'=n1+1
		priora = log(n1u)-log(nloci-n1);	
	}
	else {
		priora = log(nloci-n1+1)-log(n1);	
	}



  for (li = 0; li < nloci; li++)
  {
    sumi += C[ci]->G[li].pdg;
    sumj += C[cj]->G[li].pdg;
  }
  sumi += C[ci]->allpcalc.probg;
  sumj += C[cj]->allpcalc.probg;
  // VS 5/18/2012 added the priora
  w = exp ((beta[ci] - beta[cj]) * (sumj - sumi + priora));
  //w = exp ((beta[ci] - beta[cj]) * (sumj - sumi));
  return (w);
}



///AS: Adding a function to only calculate the sums for a particular chain
////this would then be shared with the swapper process/chain
double calcpartialswapweight (int c)
{
	int l;
	double sum = 0, w;
        for (l = 0; l < nloci; l++)
        {
        	 sum += C[c]->G[l].pdg;
        }
       sum += C[c]->allpcalc.probg;
       return (sum);
}
/* End of function calcpartialswapweight() */

int getnu (int c)
{
	int nu = 0;
	int i;
	for (i = 0; i < nloci; i++) {
		nu += grouploci_mig[c][i];
	}
	return nu;
}

double swapweight_bwprocesses(double sumi, double sumj, double betai, double betaj, int n1u, int n1)
{
	
  int i;
  double priora=0.0;
  double w = 0.0;
	
	if(n1u>n1) { // if n1'=n1+1
		priora = log(n1u)-log(nloci-n1);	
	}
	else {
		priora = log(nloci-n1+1)-log(n1);	
	}

	w = exp((betai - betaj) * (sumj - sumi + priora));
	return w;
}
/* End of function swapweight_bwprocesses */

/************ GLOBAL FUNCTIONS ******************/

void
setheat (double hval1, double hval2, int heatmode, int currentid)
{
  int ci;
  int x = 0;
  double h = 0.0;
  int i;
  //AS: changing how temperatures are saved
  //allbetas is a global array that stores all temps
  allbetas = (double *) (malloc ((numchains * numprocesses) * sizeof (double)));
  for ( i = 0; i < numprocesses * numchains; i++) {
	switch (heatmode)
	{
	case HLINEAR:
		allbetas[x] = 1.0 / (1.0 + hval1 + i);
		x++;
		break;
	case HGEOMETRIC:
		allbetas[x] = 1 - (1 - hval2) * (i) * pow (hval1, (double) (numprocesses * numchains - 1 - (i))) / (double) (numprocesses * numchains - 1);
		x++;
		break;
  	}
  } 




  //beta[0] = 1.0;
  i = 0;
  for (ci = currentid * numchains; ci < currentid * numchains + numchains; ci++)
  {
    switch (heatmode)
    {
    case HLINEAR:
      beta[i] = 1.0 / (1.0 + hval1 * ci);
      break;
    case HTWOSTEP:
    case HGEOMETRIC:
      /* sometimes, with hval 1 values > 1 get identical adjacent beta values */
      beta[i] =
        1 - (1 - hval2) * ci * pow (hval1, (double) (numchains * numprocesses - 1 - ci)) / (double) (numchains * numprocesses - 1);
      break;
    }
    //JH 6/11/2010  added this, mostly in case users use h1>1 (which may be ok but can also easily give negative values for beta)
    if (beta[i] <= 0.0 || beta[i] > 1.0)
      IM_err (IMERR_COMMANDLINEHEATINGTERMS, "command line heating terms have caused a heating value out of range. chain %d beta %lf", ci, beta[ci]);
    i++;
  }
}                               /* setheat */

//AS: adding the new swapchains_bwprocesses function
int swapchains_bwprocesses(int currentid, int step, int swaptries, int swapbetasonly, int chainduration, int burnduration, int swapA, int swapB)
{
	double metropolishastingsterm;
	int whichElementA, whichElementB;
	int procIdForA, procIdForB;
	int doISwap, areWeA;
	double sumi = 0.0;
	double sumj = 0.0;
	int n1u, n1;
	int swapvar = 0;
	int aiters = 0;
	int biters = 0;
	double abeta = 0.0;
	double bbeta = 0.0;
	//AS: temporary holders for the assignment values
	//AS: used in swapping, then freed
	//int *assignmigA;
	//int *assignthetaA;
	int p = 0;
	int q = 0;
	int l = 0;
	int sa = 0;
	int sb = 0;
	int x = 0;
	int y = 0;
	int sbmin, sbrange;
	#define SWAPDIST 7
	MPI_Status status;
	//AS: changing how we get the swapping chains, as well as adding swaptries
	for (x = 0; x < swaptries; x++) {
		//printf("Inside swaptries...attempt %d\n",x);
		//printf("numchains * numprocesses = %d\n", numchains * numprocesses);
		swapvar = 0;
		//do {
			sa = rand() % (numprocesses * numchains);
		//} while (sa < 0 || sa >= numchains * numprocesses);
		if (numchains * numprocesses < 2 * SWAPDIST + 3)
		{
			sbmin = 0;
			sbrange = numprocesses * numchains;
		} else {
			sbmin = IMAX(0, sa - SWAPDIST);
			sbrange = IMIN(numchains * numprocesses, sa + SWAPDIST) - sbmin;
		}
		do {
			sb = sbmin +  rand() % sbrange;
		} while (sb == sa /*|| sb < 0 || sb >= numchains * numprocesses*/);
	//printf("sa = %d sb = %d\n", sa, sb);
	abeta = allbetas[sa];
	bbeta = allbetas[sb];
	//printf ("abeta = %f\n", abeta);
	//printf ("bbeta = %f\n", bbeta);
	//AS: recall that sa and sb are in terms of the original allbetas array
	//AS: so need to track back to which beta is sa, and which is sb	
	//AS: now need to find where is abeta and bbeta
	whichElementA = 0;
	whichElementB = 0;
	doISwap = 0;
	areWeA = 0;
	procIdForA = -1;
	procIdForB = -1;
	for (y = 0; y < numchains; y++) {
		if (beta[y] == abeta) {
			procIdForA = currentid;
			doISwap = 1;
			areWeA = 1;
			whichElementA = y;
			if (currentid != 0) {
				//printf("Sending procidforA from %d\n", currentid);
				MPI_Send(&procIdForA, 1, MPI_INT, 0, 12345, MPI_COMM_WORLD);			
			}
			break;
		}
	}
	for (y = 0; y < numchains; y++) {
		if (beta[y] == bbeta) {
			procIdForB = currentid;
			whichElementB = y;
			doISwap = 1;
			if (currentid != 0) {
				//printf("Sending procidforB from %d\n", currentid);
				MPI_Send(&procIdForB, 1, MPI_INT, 0, 34567, MPI_COMM_WORLD);
			}
			break;
		}

	}
	if (currentid == 0 && procIdForA == -1) {
		//printf("Receiving procidforA on head...\n");
		MPI_Recv(&procIdForA, 1, MPI_INT, MPI_ANY_SOURCE, 12345, MPI_COMM_WORLD, &status);
	}
	if (currentid == 0 && procIdForB == -1) {
		//printf("Receiving procidforB on head...\n");
		MPI_Recv(&procIdForB, 1, MPI_INT, MPI_ANY_SOURCE, 34567, MPI_COMM_WORLD, &status);	
	}
	//printf("broadcasting...\n");
	MPI_Bcast(&procIdForA, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&procIdForB, 1, MPI_INT, 0, MPI_COMM_WORLD);
	//printf("pidA = %d\n",procIdForA);
	//printf("pidB = %d\n",procIdForB);

	//whichElementA = 0;
	//whichElementB = 0;
	
	//procIdForA = floor (sa / numchains);
	//procIdForB = floor (sa / numchains);
	
	//whichElementA = sa - procIdForA * numchains;
	//whichElementB = sb - procIdForB * numchains;
	
	//doISwap = areWeA = 0;
	//if (currentid == procIdForA) {
	//	doISwap = 1;
	//	areWeA = 1;
	//} else if (currentid == procIdForB) {
	//	doISwap = 1;
	//}
	if (doISwap == 1) {
		if (procIdForA == procIdForB) {
			swapchains(swaptries, swapbetasonly, currentid);
			//AS: debug only
		//	printf("Trying to swap within same process\n");
		} else {
			if (areWeA == 1) {
				if (procIdForA < procIdForB) {
					swaps_bwprocesses[procIdForB][procIdForA]++;
				} else {
					swaps_bwprocesses[procIdForA][procIdForB]++;
				}
				abeta = beta[whichElementA];
				sumi = calcpartialswapweight(whichElementA);
				n1u = getnu(whichElementA);
				MPI_Send(&abeta, 1, MPI_DOUBLE, procIdForB, 12345, MPI_COMM_WORLD);
				MPI_Recv(&bbeta, 1, MPI_DOUBLE, procIdForB, 23456, MPI_COMM_WORLD, &status);
				MPI_Send(&sumi, 1, MPI_DOUBLE, procIdForB, 4567, MPI_COMM_WORLD);
				MPI_Recv(&sumj, 1, MPI_DOUBLE, procIdForB, 5678, MPI_COMM_WORLD, &status);
				MPI_Send(&n1u, 1, MPI_INT, procIdForB, 6789, MPI_COMM_WORLD);
				MPI_Recv(&n1, 1, MPI_INT, procIdForB, 7890, MPI_COMM_WORLD, &status);
				for (p = 0; p < numprocesses * numchains; p++) {
					if (allbetas[p] == abeta) {
						break;
					}
				}
				for (q = 0; q < numprocesses * numchains; q++) {
					if (allbetas[q] == bbeta) {
						break;
					}
				}
                                //AS: adding tempbased swap counting as on Mon Aug 17 12:34:52 EDT 2015
                                if (p < q) {
                                    tempbasedswapcount[q][p]++;
                                } else {
                                    tempbasedswapcount[p][q]++;
                                }
			}
			if (areWeA != 1) {
				bbeta = beta[whichElementB];
				sumj = calcpartialswapweight(whichElementB);
				n1 = getnu(whichElementB);
				MPI_Recv(&abeta, 1, MPI_DOUBLE, procIdForA, 12345, MPI_COMM_WORLD, &status);
				MPI_Send(&bbeta, 1, MPI_DOUBLE, procIdForA, 23456, MPI_COMM_WORLD);
				MPI_Recv(&sumi, 1, MPI_DOUBLE, procIdForA, 4567, MPI_COMM_WORLD, &status);
				MPI_Send(&sumj, 1, MPI_DOUBLE, procIdForA, 5678, MPI_COMM_WORLD);
				MPI_Recv(&n1u, 1, MPI_INT, procIdForA, 6789, MPI_COMM_WORLD, &status);
				MPI_Send(&n1, 1, MPI_INT, procIdForA, 7890, MPI_COMM_WORLD);
			}
			swapvar = 0;
			if (areWeA == 1) {
				metropolishastingsterm = swapweight_bwprocesses(sumi, sumj, abeta, bbeta, n1u, n1);
				if (metropolishastingsterm >= 1.0 || metropolishastingsterm > uniform()) {
					beta[whichElementA] = bbeta;
					//AS: debug only
					//printf("Swapped between %d and %d\n", procIdForA, procIdForB);
					if (procIdForA < procIdForB) {
						swaps_bwprocesses[procIdForA][procIdForB]++;
					} else {
						swaps_bwprocesses[procIdForB][procIdForA]++;
					}
                                        //AS: adding temp based swapcounting Mon Aug 17 12:36:58 EDT 2015
                                        if (p < q) {
                                            tempbasedswapcount[p][q]++;
                                        } else {
                                            tempbasedswapcount[q][p]++;
                                        }
					swapvar = 1;
					//AS: 10/22/2014 Removed this upon discussion with VS
					/*assignmigA = (int *) malloc (nloci * sizeof(int));
					assignthetaA = (int *) malloc (nloci * sizeof(int));
					for (l = 0; l < nloci; l++) {
						assignmigA[l] = grouploci_mig[whichElementA][l];
						assignthetaA[l] = grouploci_theta[whichElementA][l];
					}
					for (l = 0; l < nloci; l++) {
						MPI_Send(&assignmigA[l], 1, MPI_INT, procIdForB, l, MPI_COMM_WORLD);
						MPI_Send(&assignthetaA[l], 1, MPI_INT, procIdForB, l*2, MPI_COMM_WORLD);
						printf("Sending %d %d \n", assignmigA[l], assignthetaA[l]);
					}
					XFREE(assignmigA);
					XFREE(assignthetaA);*/
				} else {
					swapvar = 0;
				}
				MPI_Send(&swapvar, 1, MPI_INT, procIdForB, 9876, MPI_COMM_WORLD);
			}
			if (areWeA != 1) {
				MPI_Recv(&swapvar, 1, MPI_INT, procIdForA, 9876, MPI_COMM_WORLD, &status);
				if (swapvar != 0) {
					beta[whichElementB] = abeta;
					//AS: 10/22/2014 Removed this upon discussion with VS
					/*assignmigA = (int *) malloc (nloci * sizeof(int));
					assignthetaA = (int *) malloc (nloci * sizeof(int));
					for (l = 0; l < nloci; l++) {
						MPI_Recv(&assignmigA[l], 1, MPI_INT, procIdForA, l, MPI_COMM_WORLD, &status);
						MPI_Recv(&assignthetaA[l], 1, MPI_INT, procIdForA, l*2, MPI_COMM_WORLD, &status);
					}
					for (l = 0; l < nloci; l++) {
						grouploci_mig[whichElementB][l] = assignmigA[l];
						grouploci_theta[whichElementB][l] = assignthetaA[l];
						printf("Receiving %d %d \n", grouploci_mig[whichElementB][l], grouploci_theta[whichElementB][l]);
					}
					XFREE(assignmigA);
					XFREE(assignthetaA);*/
				}
			}
		}
	}
	} //AS: closes swaptries for loop
	return swapvar;
}

/* swaps chains,  adjust heating levels if adaptive heating is invoked.
will do multiple swaps if swaptries > 1.  If a swap was done involving chain0 then swap0ok =1. Sometimes, with
swaptries > 1,  chain 0 will swap with another chain, and then swap back, restoring the current parameter
values and genealogies. This is not detected, so stuff later that checks the return from this function
will think that parameters have changed.  This should not matter */

/* CR 110929.4 get rid of extraneous args in declaration  and remove several
 * variables not being used */
//AS: adding swapbetasonly, and currentid
int
swapchains (int swaptries, int swapbetasonly, int currentid)
{
  int ci, cj, i, swap0ok;
  double metropolishastingsterm;
  void *swapptr;
  // VS pointer to the array of assignment of loci
  int *swapassign;
  int cjmin, cjrange;
  int l; //AS: indexes loci
  int assignvalue; //AS: temporary variable for swapping value, not pointer
  int p, q; //AS: indexes allbetas

#define  MINSWAP  0.1
#define  BETADJUST  1.414
#define  INCADJUST  1.414
#define  MINHEAT  0.0001
#define  MAXHEAT  0.2
#define  MAXINC  100
#define  MININC  0.1
#define  PAUSESWAP 1000
#define SWAPDIST 7
// 5/27/2010  removed HADAPT stuff

//AS: debug only
//printf("Inside swapchains function...\n");

  for (i = 0, swap0ok = 0; i < swaptries; i++)
  {
    do
    {
      ci = (int) (uniform () * numchains);
    } while (ci < 0 || ci >= numchains);

    if (numchains < 2*SWAPDIST + 3)
    {
      cjmin = 0;
      cjrange = numchains;
    }
    else
    {
      cjmin = IMAX(0,ci-SWAPDIST);
      cjrange = IMIN(numchains, ci+SWAPDIST) -cjmin;
    }
    do
    {
      cj = cjmin + (int) (uniform () * cjrange);
    } while (cj == ci || cj < 0 || cj >= numchains);
    //AS: adding temp based swap counts Mon Aug 17 12:52:05 EDT 2015
    for (p = 0 ; p < numprocesses * numchains; p++) {
        if (allbetas[p] == beta[ci]) {
            break;
        }
    }
    for (q = 0 ; q < numprocesses * numchains; q++) {
        if (allbetas[q] == beta[cj]) {
            break;
        }
    }
    if (p < q) {
        tempbasedswapcount[q][p]++;
    } else {
        tempbasedswapcount[p][q]++;
    }


    if (ci < cj)
    {
      swapcount[cj][ci]++;
    }
    else
    {
      swapcount[ci][cj]++;
    }
    metropolishastingsterm = swapweight (ci, cj);
//AS: debug only
//	printf("%f %f\n", grouploci_mig[ci],grouploci_mig[cj]);
//	printf("%f %f\n", grouploci_theta[ci],grouploci_theta[cj]);
    if (metropolishastingsterm >= 1.0
        || metropolishastingsterm > uniform ())

    {
	 if (swapbetasonly) {
	//AS: After discussion with Vitor, I removed swap of the assignment vectors
	//AS: As on 10/22/2014, looks like this works great!
	//AS: Needs further discussion though.
	/*	for (l = 0; l < nloci; l++) {
			assignvalue = grouploci_mig[ci][l];
			grouploci_mig[ci][l] = grouploci_mig[cj][l];
			grouploci_mig[cj][l] = assignvalue;
		}
		for (l = 0; l < nloci; l++) {
			assignvalue = grouploci_theta[ci][l];
			grouploci_theta[ci][l] = grouploci_theta[cj][l];
			grouploci_theta[cj][l] = assignvalue;
		}*/
	} else {


	 // VS - swap the assignment vectors for mig
	  swapassign = grouploci_mig[ci];
	  grouploci_mig[ci] = grouploci_mig[cj];
	  grouploci_mig[cj] = swapassign;

	  // VS - swap the assignment vectors for theta
	  swapassign = grouploci_theta[ci];
	  grouploci_theta[ci] = grouploci_theta[cj];
	  grouploci_theta[cj] = swapassign;
	}
//AS: debug only
//	printf("Swapping within same process here...\n");
//	printf("%f %f\n", grouploci_mig[ci],grouploci_mig[cj]);
//	printf("%f %f\n", grouploci_theta[ci],grouploci_theta[cj]);
	if (swapbetasonly) {
		swapbetas(ci, cj);
	} else {
//AS: debug only      
//	printf("Before swap beta[%d]=%f, beta[%d]=%f\n",ci,beta[ci],cj,beta[cj]);
	swapptr = C[ci];
      C[ci] = C[cj];
      C[cj] = swapptr;
//AS: debug only
//	printf("After swap beta[%d]=%f, beta[%d]=%f\n",ci,beta[ci],cj,beta[cj]);
	}
      if (ci < cj)

      {
        swapcount[ci][cj]++;
      }
      else
      {
        swapcount[cj][ci]++;
      }
    //AS: adding temp based swap counting Mon Aug 17 12:55:03 EDT 2015
      if ( p < q) {
        tempbasedswapcount[p][q]++;
      } else {
        tempbasedswapcount[q][p]++;
      }

      if (ci == 0 || cj == 0)
        swap0ok |= 1;
    }
	if (numprocesses > 1)
		break;
  }
  return swap0ok;
}                               /* swapchains */

//AS: adding function that will only swap values of beta
void swapbetas (int ci, int cj) {
	double btemp = 0.0;
	//AS: debug only
	//printf("Before swap beta[%d]=%f, beta[%d]=%f\n",ci,beta[ci],cj,beta[cj]);
	btemp = beta[ci];
	beta[ci] = beta[cj];
	beta[cj] = btemp;
	//AS: debug only
	//printf("After swap beta[%d]=%f, beta[%d]=%f\n",ci,beta[ci],cj,beta[cj]);
	
	return;
}

//AS: adding currentid
//AS: also adding print options for swaps between processors
void
printchaininfo (FILE * outto, int heatmode, double hval1,
                double hval2, int currentid)
{
  int i, j, x, y, z;
  #ifdef MPI_ENABLED
  MPI_Status status;
  #endif
  if (currentid == 0) {
  fprintf (outto, "\nCHAIN SWAPPING BETWEEN SUCCESSIVE CHAINS: ");
  switch (heatmode)

  {
  case HLINEAR:
    fprintf (outto, " Linear Increment  term: %.4f\n", hval1);
    break;
  case HTWOSTEP:
    fprintf (outto, " Twostep Increment  term1: %.4f term2: %.4f\n", hval1,
             hval2);
    break;
  case HGEOMETRIC:
    fprintf (outto, " Geometric Increment  term1: %.4f term2: %.4f\n",
             hval1, hval2);
    break;
  }
  fprintf (outto,
           "-----------------------------------------------------------------------------\n");
  }

	#ifdef MPI_ENABLED
		//MPI_Barrier(MPI_COMM_WORLD);
		if (numprocesses > 1) {
			for (x = 0; x < numprocesses; x++) {
				for (y = 0; y < numprocesses; y++) {
					MPI_Reduce(&swaps_bwprocesses[x][y], &swaps_rec_bwprocesses[x][y], 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
				}
			}
                        
                        //AS: adding temp based swap counting as on Mon Aug 17 12:09:10 EDT 2015
                        for (x = 0; x < numprocesses * numchains; x++) {
                            for (y = 0 ; y < numprocesses * numchains; y++) {
                                MPI_Reduce(&tempbasedswapcount[x][y], &tempbased_rec_swapcount[x][y], 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
                            }
                        }

			if (currentid == 0) {
				for (x = 0; x < numchains; x++) {
					for (y = 0; y < numchains; y++) {
						swapcount_bwprocesses[x][y] = swapcount[x][y];
					}
				}
			}
			for (x = 1; x < numprocesses; x++) {
				if (currentid == x) {
					for (y = 0; y < numchains; y++) {
						for (z = 0; z < numchains; z++) {
							MPI_Send(&swapcount[y][z], 1, MPI_INT, 0, 1234, MPI_COMM_WORLD);
						}
					}
				}
				if (currentid == 0) {
					for (y = 0; y < numchains; y++) {
						for (z = 0; z < numchains; z++) {
							MPI_Recv(&swapcount_bwprocesses[x * numchains + y][x * numchains + z], 
								1, MPI_INT, x, 1234, MPI_COMM_WORLD, &status);
						}
					}
				}
			}
		}
	#endif
	//AS: debug
	//printf("Have i come outside this collation??\n");
  if (outto == stdout)

  {
    if (currentid == 0) {
    fprintf (outto, "beta terms :");
    for (i = 0; i < numprocesses * numchains; i++)
      fprintf (outto, "|%2d %5.3f", i, allbetas[i]);
    fprintf (outto, "\n");
    fprintf (outto, "Swap rates :");
    for (i = 0; i < numchains * numprocesses - 1; i++)
      //if (swapcount[i + 1][i])
        fprintf (outto, "|%2d %5.3f", i,
                 swapcount_bwprocesses[i][i + 1] / (float) swapcount_bwprocesses[i + 1][i]);
    fprintf (outto, "\n");
    fprintf (outto, "Swap counts:");
    for (i = 0; i < numchains * numprocesses - 1; i++)
      fprintf (outto, "|%2d %5ld", i, swapcount_bwprocesses[i][i + 1]);
    fprintf (outto, "\n\n");
   }
  }
  else
  {
    if (numprocesses > 1 && currentid == 0) {
    fprintf (outto, "Chain   #Swaps  Rate \n");
    for (i = 0; i < numchains * numprocesses - 1; i++)

    {
      //if (swapcount[i + 1][i])
        fprintf (outto, " %3d   %5ld  %7.4f\n", i,
                 swapcount_bwprocesses[i][i + 1],
                 swapcount_bwprocesses[i][i + 1] / (float) swapcount_bwprocesses[i + 1][i]);

      //else
      //  fprintf (outto, " %3d  %7.4f %5ld  \n", i, beta[i],
        //         swapcount[i][i + 1]);
    }
    fprintf (outto, " %3d    na      na\n", i);
    fprintf (outto, "\n");
   }
  }
  if (numprocesses > 1 && currentid == 0) {
	fprintf(outto, "\nCHAIN SWAPPING BETWEEN PROCESSES:\n");
	fprintf(outto, "-----------------------------------\n");
	fprintf(outto, "Processor1  Processor2  #SwapAttempts\n");
	for (i = 0; i < numprocesses; i++) {
		for (j = 0; j < numprocesses; j++) {
			if ( i > j) {
				fprintf(outto, " %3d        %3d        %5ld\n", i, j, swaps_rec_bwprocesses[i][j]);
			}
		}
	}
	fprintf(outto, "\n\n");
	fprintf(outto, "Processor1  Processor2  #Swaps\n");
	for (i = 0; i < numprocesses; i++) {
		for (j = 0; j < numprocesses; j++) {
			if (i < j) {
				fprintf(outto, " %3d       %3d        %5ld\n", i, j, swaps_rec_bwprocesses[i][j]);
			}
		}
	}
	fprintf(outto, "\n\n");
    //AS: adding printing temp based swap counts Mon Aug 17 12:57:23 EDT 2015
    fprintf (outto, "\nCHAIN SWAPPING BETWEEN ADJACENT TEMPERATURES:\n");
    fprintf (outto, "----------------------------------------------\n");
    fprintf (outto, "Temp1     Temp2     #Swaps    #Attempts     Rate\n");
    for (i = 0; i < numprocesses * numchains - 1; i++) {
        if (tempbased_rec_swapcount[i+1][i] > 0) {
            fprintf(outto, " %7.4f    %7.4f    %5ld    %5ld    %7.4f\n", allbetas[i], allbetas[i+1], tempbased_rec_swapcount[i][i+1],
                    tempbased_rec_swapcount[i+1][i], (float) tempbased_rec_swapcount[i][i+1]/(float) tempbased_rec_swapcount[i+1][i]);
            } else if (tempbased_rec_swapcount[i+1][i] == 0) {
                fprintf(outto, " %7.4f    %7.4f   %5ld   %5ld   na\n", allbetas[i], allbetas[i+1], tempbased_rec_swapcount[i][i+1],
                    tempbased_rec_swapcount[i+1][i]);
            }
    }
    fprintf(outto, "\n\n");
                
  }
		
}                               /* printchaininfo */
