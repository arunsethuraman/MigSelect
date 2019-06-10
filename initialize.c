/* MigSelect 2013-2019 Arun Sethuraman, Vitor Sousa, Jody Hey */
/* This code was developed based on IMa2  2009-2010  Jody Hey, Rasmus Nielsen and Sang Chul Choi*/
#undef GLOBVARS
#include "imamp.h"
#include "updateassignment.h"
#include "update_gtree_common.h"

/* 1/8/09 rewrote this using new data structures   */
/* Most initializations  */

/*********** LOCAL STUFF **********/

extern double pi[MAXLOCI][4];
extern int urri[2 * MAXLOCI][2 * MAXLOCI];      // used mostly in update_mc_params.c
extern double urrlow[2 * MAXLOCI][2 * MAXLOCI], urrhi[2 * MAXLOCI][2 * MAXLOCI];        // used mostly in update_mc_params.c

static double **uvals;
static char startpoptreestring[POPTREESTRINGLENGTHMAX];
static double geomeanvar;
static int **numsitesIS;        // maxpops * maxloci   temporarily holds the number of polymorphic sites in loci with infinite sites model
static double uterm;
static double *popsizeprior_fromfile,**mprior_fromfile;
static int maxpossiblemigrateparams;
/* used if modeloptions[ANCESTRALPOPSIZESSETTOLARGESTSAMPLED]==1 
  varrank[] holds the sum of 4Nu estimates across loci for the sampled populations - fill in set_uvals()
  varranktree[i]:
    if i< npops,  then it is just the sampled population number and its corresponding value of varrank[i]
    else  it holds the sampled population number for which the population size parameter for this ancestral population will be the same as

*/
static double varrank[MAXPOPS] = {0.0};
static struct {
  int samppopnum;
  double var; 
}  varranktree[2*MAXPOPS - 1];


// 1/12/09  why is this stuff on struct edgemiginfo here  ??  what does it do??

/* We initialize these four in function [[IMA_initmemory_aftersetpopntree]]. */
/* global variables of update_gtree */
extern struct edgemiginfo oldedgemig;
extern struct edgemiginfo oldsismig;
extern struct edgemiginfo newedgemig;
extern struct edgemiginfo newsismig;

/* prototypes */
static void addpopsizeidtext(int *fpstri, char fpstr[]);
static void fillvarranktree(int nowpop);
static int numvarHKY (int li, int b, int e);
static int checkaresis (int period, int i, int j);
static void setimigarray (int ii);
// VS
static void setimigarray_vs (int ii); // VS function added
static void setup_iparams ();
static void set_nomigrationchecklist ();
static void setuinfo (double summut);
static double set_uvals (void);
static void set_x (struct value_record *v, int isint);
static void init_value_record (struct value_record *v, int isint);
static void init_g_rec (int li);
static void init_a_rec (int li);
static void init_lpgpd_v (void);
static void fillmrec(int j, int i, int pi, int pj);
static void init_migration_counts_times (void);
static void init_mutation_scalar_rec (int li);
static void start_setup_L (char infilename[], int *fpstri, char fpstr[]);
static void add_priorinfo_to_output(char priorfilename[],int *fpstri, char fpstr[]);
// VS modified this function to print all the infor about the groups of loci
static void add_priorinfo_to_output_vs(char priorfilename[],int *fpstri, char fpstr[]); 
static void start_setup_C (void);
static void finish_setup_C ();
static void set_tvalues (int ci);
static void setup_T ();
static void finish_setup_L (void);
static void reportparamcounts(int *fpstri, char fpstr[]);
// VS 5/18/2012 added function to setup the variables related with the assignment of loci
static void setup_assignLoci();

/******** LOCAL FUNCTIONS ************/

/* used if modeloptions[ANCESTRALPOPSIZESSETTOLARGESTSAMPLED]==1  */
void addpopsizeidtext(int *fpstri, char fpstr[])
{
  int i;
  SP"\nAncestral Population Size Parameter Simplification\n");
  SP"--------------------------------------------------\n");
  SP"\nEach ancestral population size set to be the same as the size its largest descendant sampled population\n");
  SP"Sampled Populations and Polymorphism Summaries\n"); 
  SP"\tPop#\tpolysum\n");
  for (i=0;i<npops;i++)
    SP"\t%d\t%.3lf\n",i, varrank[i]);
  SP"\nAncestral Populations Connected to which Sampled Populations\n");
  SP"\tAncPop#\tSampPop#\n");
  for (i=npops;i<numtreepops;i++)
    SP"\t%d\t%d\n",i,varranktree[i].samppopnum);
  SP"\n");
} //addpopsizeidtext(int *fpstri, char fpstr[])

void fillvarranktree(int nowpop)
{
  if (nowpop < npops)
  {
    varranktree[nowpop].samppopnum = nowpop;
    varranktree[nowpop].var = varrank[nowpop];
  }
  else
  {
    fillvarranktree(C[0]->poptree[nowpop].up[0]);
    fillvarranktree(C[0]->poptree[nowpop].up[1]);
    if (varranktree[C[0]->poptree[nowpop].up[0]].var < varranktree[C[0]->poptree[nowpop].up[1]].var)
    {
      varranktree[nowpop] = varranktree[C[0]->poptree[nowpop].up[1]];
    }
    else
    {
      varranktree[nowpop] = varranktree[C[0]->poptree[nowpop].up[0]];
    }
  }

}// fillvarranktree,  recursive 

int
numvarHKY (int li, int b, int e)
{
  int i, j, tot = 0;
  for (i = 0; i < L[li].numsites; i++)

  {
    j = b + 1;
    while (j < e && L[li].seq[b][i] == L[li].seq[j][i])
      j++;
    if (j < L[li].numgenes)
      tot++;
  }
  return tot;
}

int
checkaresis (int period, int i, int j)  // returns 1 if populations i and j are both in period, and they are sister populations
{
  int aresis, k;
  aresis = 0;
  if (ISELEMENT (C[0]->plist[period][i], C[0]->periodset[period])
      && ISELEMENT (C[0]->plist[period][j], C[0]->periodset[period]))
  {
    for (k = period + 1; k < npops; k++)
    {
      if ((C[0]->droppops[k][0] == C[0]->plist[period][i]
           && C[0]->droppops[k][1] == C[0]->plist[period][j])
          || (C[0]->droppops[k][0] == C[0]->plist[period][j]
              && C[0]->droppops[k][1] == C[0]->plist[period][i]))
        aresis = 1;
    }
  }
  return aresis;
}                               //checkaresis 

void
setimigarray (int ii)
{
  int jj;
  int gp; // VS

  if (modeloptions[EXPOMIGRATIONPRIOR])
  {
    if (!calcoptions[USEPRIORFILE])
      imig[ii].pr[gp].mean = mprior; // VS gp
    else
    {
      if (modeloptions[ONEMIGRATIONPARAMETER])
        imig[ii].pr[gp].mean =  max_m_from_priorfile; // VS gp
      else
        imig[ii].pr[gp].mean = imig[ii].pr[gp].max; // VS gp
	  // this was previously set in setup_iparams()
    }
    imig[ii].pr[gp].max = EXPOMIGPLOTSCALE * imig[ii].pr[gp].mean; // VS gp
    imig[ii].pr[gp].min = 0; // VS gp
  }
  else
  {
    if (!calcoptions[USEPRIORFILE])
      imig[ii].pr[gp].max = DMAX (mprior, MPRIORMIN); // VS gp
    else
    {
      if (modeloptions[ONEMIGRATIONPARAMETER])
        imig[ii].pr[gp].max =  max_m_from_priorfile; // VS gp
      //else this was previously set in setup_iparams()
    }
    imig[ii].pr[gp].min = 0; // VS gp
    imig[ii].pr[gp].mean = (imig[ii].pr[gp].max - imig[ii].pr[gp].min) / 2.0; // VS gp
  }
  imig[ii].xy[gp] = calloc (GRIDSIZE, sizeof (struct plotpoint)); // VS gp
  for (jj = 0; jj < GRIDSIZE; jj++)
  {
    imig[ii].xy[gp][jj].x = // VS gp
      imig[ii].pr[gp].min + // VS gp
      ((jj + 0.5) * (imig[ii].pr[gp].max - imig[ii].pr[gp].min)) / GRIDSIZE; // VS gp
  }
  imig[ii].wp.n = 0;

}                               // setimigarray


// SETIMIGARRAY_VS
// Initializes the migration parameter for prior pr values and initializes plotpoints xy.
// This function was based on previous setimigarray and was changed by VS to allow each locus to have its own demographic parameters
// this funciton gets as input the index of the imig array, i.e. imig[ii]
// INPUT:
//	int ii - index of the imig array, e.g. usually imig[0] refers to mig1->2 param, imig[1] to mig2->1 param and so on
// RETURN:
//	initializes imig[ii].pr and imig[ii].xy
void setimigarray_vs (int ii)
{
  // VS
  // where is the memory allocated to imig param related variables?
  // Not completely sure if this function initializes all the imigration params
  // if so, where are the values read from the prior file saved, and where should I initialize the 
  // prpoint pointer of each locus to its group?

  // VS - should I initialize the prpointer variable here?
  // or should the prpointer be initialized when reading the prior input file

  // VS
  // When read from command line nbgroupsloci_mig=1

  // VS added the variable li and gp
  //	int gp - index of the group to which the locus belongs
  int jj, gp;

  // if using exponential priors
  if (modeloptions[EXPOMIGRATIONPRIOR])
  {
    if (!calcoptions[USEPRIORFILE])
      imig[ii].pr[0].mean = mprior; // VS if reading from command line, assume there is only one group
    else
	// if reading from the prior file
    {
		if (modeloptions[ONEMIGRATIONPARAMETER]) {
			// Go through all groups of loci and save max parameter for each
			// and multiply the maximum with the relative mig
			for(gp = 0; gp < nbgroupsloci_mig; gp++) {
				// even if there is only one migration param, different groups may have different relative migration rates
				imig[ii].pr[gp].mean = max_m_from_priorfile*relmig[gp];
				imig[ii].pr[gp].max = EXPOMIGPLOTSCALE * imig[ii].pr[gp].mean;
				imig[ii].pr[gp].min = 0;
			}
		}
		// if there are more than one migration param
		else {
			// Go through all groups of loci and save max parameter for each
			// and multiply the maximum with the relative mig
			for(gp = 0; gp < nbgroupsloci_mig; gp++) {
				imig[ii].pr[gp].mean = imig[ii].pr[gp].max*relmig[gp];
				imig[ii].pr[gp].max = EXPOMIGPLOTSCALE * imig[ii].pr[gp].mean;
				imig[ii].pr[gp].min = 0;
			}
		}
    }
  }
  else 
  // if Uniform priors
  {
    if (!calcoptions[USEPRIORFILE])
      // VS if reading from command line, assume there is only one group
      imig[ii].pr[0].max = DMAX (mprior, MPRIORMIN); 
    else
	// if reading from the prior file
    {
		// if there is only one migration parameter
		if (modeloptions[ONEMIGRATIONPARAMETER]) {
			// Go through all groups of loci and save max parameter for each
			// and multiply it by relmig of each group
			for(gp = 0; gp < nbgroupsloci_mig; gp++) {
				imig[ii].pr[gp].max = max_m_from_priorfile*relmig[gp];
				imig[ii].pr[gp].max = (imig[ii].pr[gp].max - imig[ii].pr[gp].min) / 2.0;
				imig[ii].pr[gp].min = 0;
			}
		}
        //else this was previously set in setup_iparams() - the case for non onemigration param is done in setup_iparams
    } 
  }

  // VS
  // Allocate memory for xy
  // need to allocate memory to each group loci
  // why is this done here and not in setup_iparams?
  for(gp=0; gp<nbgroupsloci_mig; gp++) {
	  imig[ii].xy[gp] = calloc (GRIDSIZE, sizeof (struct plotpoint));
	  for (jj = 0; jj < GRIDSIZE; jj++)
	  {
		imig[ii].xy[gp][jj].x = 
		  imig[ii].pr[gp].min +
		  ((jj + 0.5) * (imig[ii].pr[gp].max - imig[ii].pr[gp].min)) / GRIDSIZE;
	  }
  }
  // set the wp.n - this does not depend on the groups of loci. 
  // It depends on the migration pair to which the ii index refer.
  imig[ii].wp.n = 0;
}                               // setimigarray_vs



// SETUP_IPARAMS
// this function initialises and allocate memory to the iparams arrays
// initialize the splitting rate, population size and migration instances of struct i_param.
// This function is complicated with lots of sections because of many user options.
// Also, set values of position markers for gsampinf
// INPUT:
//  Global variables: 
//		int numpopsizeparams : number of theta params
//		iparam itheta
//		nbgroupsloci_theta
//		lastperiodnumber
//		C
void setup_iparams (void)
{
  int i, j, k, ci = 0, ii, jj, kk, pj, ni, mi, mcheck, m; // VS probably add m
  int gp; // VS - index of group loci
  char temps[PARAMSTRLEN];

  /* set up the population size parameters */ 
  // get the number of theta params depending on the options of the program 
  if (modeloptions[PARAMETERSBYPERIOD])
  {
    numpopsizeparams = (npops * (npops + 1)) / 2;       // every population in every period gets a parameter
  }
  else if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 1)
  {
    numpopsizeparams = npops;
  }
  else if (modeloptions[ANCESTRALPOPSIZESSETTOLARGESTSAMPLED] == 1)
  {
    numpopsizeparams = npops;
    fillvarranktree(numtreepops-1);
  }
  else
  {
    numpopsizeparams = numtreepops;     // every distinct population gets a parameter
  }

  // VS
  // itheta is an array of pointers to and array of i_param structures
  // itheta[i][l] refer to ith theta (i=1,...,2*npop-1) at group locus l (l=1,...,nbgroupsloci)
  itheta = malloc (numpopsizeparams * sizeof (struct i_param));

  // VS
  // Allocate memory to the arrays of size nbloci or nbgrouploci_theta with the structure itheta
  // in order to allow each locus to have its own migration rate
  for(i = 0; i < numpopsizeparams; i++) {
	  itheta[i].pr = (struct priorvalues *) malloc(nbgroupsloci_theta * sizeof(struct priorvalues)); 
	  itheta[i].xy = (struct plotpoint  **) malloc(nbgroupsloci_theta * sizeof(struct plotpoint *));
	  //itheta[i].indexgroup = (int *) malloc(nloci * sizeof(int));
	  // Another option is just to point indexgroup to grouploci_theta (global array of ints initialized in readpriorfile)
	  // grouploci_theta is allocated and initialized in readpriorfile()
	  //itheta[i].indexgroup = grouploci_theta;
  }

  // Initialize variables of itheta
  // setting min, max and xy grid points
  for (i = 0; i < numpopsizeparams; i++)
  {
	// If using prior file need to go though the theta prior groups
	if (calcoptions[USEPRIORFILE]) {
		for(gp=0; gp < nbgroupsloci_theta; gp++) {
			// set max and min for prior limits
			itheta[i].pr[gp].max = popsizeprior_fromfile[i]*reltheta[gp];
			itheta[i].pr[gp].min = 0;
			// allocate memory for the plot xy points
			// initialize the x points to grid of points between  min and max theta
			itheta[i].xy[gp] = calloc (GRIDSIZE, sizeof (struct plotpoint)); 
			for (m = 0; m < GRIDSIZE; m++)
			{
				itheta[i].xy[gp][m].x = itheta[i].pr[gp].min +((m + 0.5) * (itheta[i].pr[gp].max - itheta[i].pr[gp].min)) / GRIDSIZE;
			}

		}
	}
	// If read from command line it is assumed there is only one group
	// and that all loci share the same param prior
	else {
		itheta[i].pr[0].max = thetaprior;
		itheta[i].pr[0].min = 0;
		itheta[i].xy[0] = calloc (GRIDSIZE, sizeof (struct plotpoint));
		for (j = 0; j < GRIDSIZE; j++)
		{
		  itheta[i].xy[0][j].x = itheta[i].pr[0].min + ((j + 0.5) * (itheta[i].pr[0].max - itheta[i].pr[0].min)) / GRIDSIZE;
		}
	}
  }

  // Initialize wp.p, wp.n and wp.r 
  // These variables are the index of a given param
  // and are the same for all groups of loci.
  // Dependening on whether the parameters were speciefies by period
  // or by population affects these index
  if (modeloptions[PARAMETERSBYPERIOD])
  {
	for (i = 0, ii = 0, k = npops; i <= lastperiodnumber; i++, k--)
	{
		for (j = 0; j < k; j++)
		{
			pj = C[0]->plist[i][j];
			sprintf (itheta[ii].str, "q%d,%d", i, pj);
			itheta[ii].b = C[0]->poptree[pj].b;
			itheta[ii].e = C[0]->poptree[pj].e;
			itheta[ii].wp.n = 1;
			itheta[ii].wp.p = malloc (itheta[ii].wp.n * sizeof (int));
			itheta[ii].wp.r = malloc (itheta[ii].wp.n * sizeof (int));
			*itheta[ii].wp.p = i;
			*itheta[ii].wp.r = j;
			ii++;
		}
	}
  }
  // if parameters are specified by population
  else  
  {
	for (ii = 0; ii < numpopsizeparams; ii++)
	{
		sprintf (itheta[ii].str, "q%d", ii);
		itheta[ii].b = C[0]->poptree[ii].b;
		itheta[ii].e = C[0]->poptree[ii].e;
		itheta[ii].wp.n = 0;

		if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 1)
		{
			assert (modeloptions[ANCESTRALPOPSIZESSETTOLARGESTSAMPLED] == 0);
			itheta[ii].wp.n = 1;
			if (itheta[ii].wp.n > 0)
	        {
				itheta[ii].wp.p = malloc (itheta[ii].wp.n * sizeof (int));
				itheta[ii].wp.r = malloc (itheta[ii].wp.n * sizeof (int));
			}
			ni = 0;
			for (i = 0; i < lastperiodnumber; i++)
			{
				for (j = 0; j < npops; j++)
				{
					if (ii == C[0]->plist[i][j])
					{
						itheta[ii].wp.p[ni] = i;
						itheta[ii].wp.r[ni] = j;
						ni++;
					}
				}
				assert (ni == itheta[ii].wp.n);
			}
		}
		else
		{
			if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 0
				&& modeloptions[ANCESTRALPOPSIZESSETTOLARGESTSAMPLED] == 0)
			{
			for (i = 0, k = npops; i <= lastperiodnumber; i++, k--)
				for (j = 0; j < k; j++)
					itheta[ii].wp.n += (ii == C[0]->plist[i][j]);
					if (itheta[ii].wp.n > 0)
					{
						itheta[ii].wp.p = malloc (itheta[ii].wp.n * sizeof (int));
						itheta[ii].wp.r = malloc (itheta[ii].wp.n * sizeof (int));
					}
					ni = 0;
					for (i = 0, k = npops; i <= lastperiodnumber; i++, k--)
						for (j = 0; j < k; j++)
						{
							if (ii == C[0]->plist[i][j])
							{
								itheta[ii].wp.p[ni] = i;
								itheta[ii].wp.r[ni] = j;
								ni++;
							}
						}
					assert (ni == itheta[ii].wp.n);
			}
			else if (modeloptions[ANCESTRALPOPSIZESSETTOLARGESTSAMPLED] == 1)
			{
				/* need to know what sampled population is  */
				for (i = 0, k = npops; i <= lastperiodnumber; i++, k--)
					for (j = 0; j < k; j++)
						itheta[ii].wp.n += (ii == C[0]->plist[i][j] || ii==varranktree[C[0]->plist[i][j]].samppopnum);
						if (itheta[ii].wp.n > 0)
						{
							itheta[ii].wp.p = malloc (itheta[ii].wp.n * sizeof (int));
							itheta[ii].wp.r = malloc (itheta[ii].wp.n * sizeof (int));
						}
					ni = 0;
			 		for (i = 0, k = npops; i <= lastperiodnumber; i++, k--)
						for (j = 0; j < k; j++)
						{
							if (ii == C[0]->plist[i][j] || ii==varranktree[C[0]->plist[i][j]].samppopnum)
							{
								itheta[ii].wp.p[ni] = i;
								itheta[ii].wp.r[ni] = j;
								ni++;
							}
						}
				assert (ni == itheta[ii].wp.n);
			}
		}
	}
  } // end else of if parameters are specified by population and not by period

  // VS CHECK
  // wp refers to the gweight matrix, in order to know what entries in the matrix to sum
  printf("\n\nLOOK INTO ITHETA\n");
  for(i = 0 ; i < numpopsizeparams; i++) {
	// Go through all the relative parameters in the input file
	printf("\nitheta[%i] %s", i, itheta[i].str);
	printf("\nitheta[%i] first period %i", i, itheta[i].b);
	printf("\nitheta[%i] last period %i", i, itheta[i].e);

	for(gp = 0; gp < nbgroupsloci_theta; gp++) {
		printf("\n itheta[%i].pr[%i].min is %g", i, j, itheta[i].pr[gp].min);
		printf("\n itheta[%i].pr[%i].max is %g", i, j, itheta[i].pr[gp].max);

		for(m = 0; m < itheta[i].wp.n; m++) {
			printf("\n itheta[%i].wp.p[%i] is %i", i, m, itheta[i].wp.p[m]);
			printf("\n itheta[%i].wp.r[%i] is %i", i, m, itheta[i].wp.r[m]);
		}
	}
  }

  // INITIALIZE MIGRATION IPARAM STRUCTURE ----------------------------------------
  // Three, somewhat repetitive loops here:
  //   1. count how many migration parameters are needed.
  //   2. build  most of the imig[] structures, including str, b, e; start to build wp part of imig[]
  //   3. finish building wp parts of imig[]
  // 
  //   there is a complex series of checks to ensure that the intended model is being followed
  //   loop through periods
  //    check order:
  //    MIGRATIONBETWEENSAMPLED  - migration only between sampled populations
  //    PARAMETERSBYPERIOD  - every population size and migration parameter applies for only 1 period
  //    NOMIGBETWEENNONSISTERS  - set migration to zero between non-sister populations
  //    SINGLEMIGRATIONBOTHDIRECTIONS - one migration parameter for migration between two populations (symmetric migration)
  //    if ONEMIGRATIONPARAMETER applies  then there is only one migration rate,  whenever there is migration in the model
  if (!modeloptions[PARAMETERSBYPERIOD]) {
    maxpossiblemigrateparams = 2*(npops-1)*(npops-1);
  }
  // count how many migration parameters are needed, use mi
  // set up a standard sequence of checks that get repeated when building imig[]
  mi = 0; // VS mi saves the nb of migration parameters needed
  if (!modeloptions[NOMIGRATION] && npops > 1)
  {
    if (modeloptions[ONEMIGRATIONPARAMETER])
    {
      mi = 1; // there is only one migration pair 
	  // VS
	  // nbgroupsloci_mig = 1; // it is assumed that there is only one group of loci
	  // this is not necessarily the case, because we can have only one migration parameter
	  // but two different groups of loci
      goto outsidefirstloop; // get out of this big condition section
    }
    for (k = 0; k < lastperiodnumber; k++) {
      for (i = 0; i < npops - k - 1; i++) {
        for (j = i + 1; j < npops - k; j++) {
			// tricky condition.  proceed to set mcheck if:
			//	we are considering only sampled populations (k==0)  or
			//	(not considering just sampled populations  AND
			//	EITHER we are considering each parameter to apply only to one period
			//		  OR one of the two populations being considered was not in the previous (k-1) period (meaning it is a new ancestral pop) */   
		    // if mcheck ends up as 1, we have one more parameter to add
			// if it is also the case that migration is not the same in both directions,  then we have two parameters to add
			// mcheck can get set to 0 if we require two populations to be sisters, and they are not
			// or if the prior given for that pair of populations is zero
			if ( k == 0 ||
                  ( modeloptions[MIGRATIONBETWEENSAMPLED]==0 &&
                  (modeloptions[PARAMETERSBYPERIOD]==1 ||
                  ( !ISELEMENT (C[ci]->plist[k][i], C[ci]->periodset[k - 1]) ||
                    !ISELEMENT (C[ci]->plist[k][j], C[ci]->periodset[k - 1])
                  ) ) ) ) 
			{ // if the above tricky condition holds
				if (modeloptions[NOMIGBETWEENNONSISTERS])
				  mcheck = checkaresis (k, i, j); // Check if the two populations are sisters
				else
				{
				  if (calcoptions[USEPRIORFILE]==0 || (calcoptions[USEPRIORFILE] && mprior_fromfile[C[ci]->plist[k][i]][C[ci]->plist[k][j]] > MINPARAMVAL))
					mcheck = 1;
				  else
					mcheck = 0;
				}
				if (mcheck)
				{
				  mi++;
				  mi += modeloptions[SINGLEMIGRATIONBOTHDIRECTIONS]==0;
				}
			} // end of if tricky condition
		} // end of for loop j
	  } // end of for loop i
	} // end of for loop k
	// when getting out of for loop this is where we should end up
	outsidefirstloop:
	assert (!assignmentoptions[POPULATIONASSIGNMENTINFINITE] || mi == npops * (npops - 1));

	//  now do another loop, much like the previous one,
	//  but now start to build the imig[] structures, including str, b, e
	//  start to build the wp part of imig[] and determine wp.n
	nummigrateparams = mi;
	
	// VS 
	// ALLOCATE MEMORY FOR IMIG array - this will have the size of number of migration pairs
	imig = malloc (nummigrateparams * sizeof (struct i_param));
	// Allocate memory the arrays of size nbloci or nbgrouploci_mig with the structure imig
	// in order to allow each locus to have its own migration rate
	for(i = 0; i < nummigrateparams; i++) {
	  imig[i].pr = (struct priorvalues *) malloc(nbgroupsloci_mig * sizeof(struct priorvalues)); 
	  imig[i].xy = (struct plotpoint  **) malloc(nbgroupsloci_mig * sizeof(struct plotpoint *));
	  //imig[i].indexgroup = (int *) malloc(nloci * sizeof(int));
	  //imig[i].indexgroup = grouploci_mig;
	  // indexgroup is pointing to grouploci_mig (global array of ints initialized in readpriorfile)
	  // grouploci_mig is allocated and initialized in readpriorfile()

	}

	// VS 
	// initializing the values of the matrix - this cannot be done here
	// need to check if this should be done inside setmigarray or outside
	// it should be done here. Setmigarray will correct this values for when
	// exponential priors or onemigparam is true
	// i believe this part is return an error if the ONEMIGPARAM option is true
	// hence I think this should be done later
	/*for(gp = 0; gp < nbgroupsloci_mig; gp++)
	{
		imig[mi].pr[gp].max = relmig[gp]*mprior_fromfile[C[ci]->plist[k][i]][C[ci]->plist[k][j]];
		if (imig[mi].pr[gp].max <= MINPARAMVAL)
		  IM_err(IMERR_PRIORFILEVALS,"migration rate set too low from %d to %d",i,j);
	}*/



	// if there is only one migration parameter
	if (modeloptions[ONEMIGRATIONPARAMETER]) {
      mi = 0;
	  // VS 
	  // setimigarray initializes the pr variables and xy plotpoints
	  // setmigarray will allocate memory and initialize variables?
      // setimigarray (mi); 
	  setimigarray_vs (mi); 	
      imig[mi].wp.n = 0;
      imig[mi].b = 0;
      imig[mi].e = lastperiodnumber - 1;
      sprintf (imig[mi].str, "m_all");
	}
	else {
	  mi = -1;  // mi is the nb of migration parameters
	}

	// go throup each period
	// go through populations pairs in each period
	for (k = 0; k < lastperiodnumber; k++) {
		for (i = 0; i < npops - k - 1; i++) {
		  for (j = i + 1; j < npops - k; j++) {
				if ((k == 0) ||
					(!modeloptions[MIGRATIONBETWEENSAMPLED] && (modeloptions[PARAMETERSBYPERIOD] ||
					((!ISELEMENT(C[ci]->plist[k][i], C[ci]->periodset[k - 1]))
					 ||  (!ISELEMENT(C[ci]->plist[k][j], C[ci]->periodset[k - 1]))))))
				{
					if (modeloptions[NOMIGBETWEENNONSISTERS])
						mcheck = checkaresis (k, i, j);
					else
					{
						if (!calcoptions[USEPRIORFILE] || (calcoptions[USEPRIORFILE] && mprior_fromfile[C[ci]->plist[k][i]][C[ci]->plist[k][j]] > MINPARAMVAL))
							mcheck = 1;
						else
							mcheck = 0;
					}
					if (mcheck)
					{
						if (!modeloptions[ONEMIGRATIONPARAMETER])
						// If there are more than one migration parameter
						{
							mi++;
							if (calcoptions[USEPRIORFILE])
							// if reading the prior file
							{
								// VS added the [gp] but need to check were to allocate memory
								// for loop through nbgroups. Mulitply the maxprior by relative prior
								for(gp = 0; gp < nbgroupsloci_mig; gp++)
								{
									imig[mi].pr[gp].max = relmig[gp]*mprior_fromfile[C[ci]->plist[k][i]][C[ci]->plist[k][j]];
									if (imig[mi].pr[gp].max <= MINPARAMVAL)
									  IM_err(IMERR_PRIORFILEVALS,"migration rate set too low from %d to %d",i,j);
									imig[mi].pr[gp].min = 0;
									imig[mi].pr[gp].mean = (imig[mi].pr[gp].max-imig[mi].pr[gp].min)/2;
								}

							}
							// VS
							// Setmigarray function will intialize the values
							// it seems that it is doing what it is done above
							// function setmigarray_vs initializes the variables for each group
							// setimigarray (mi);
							setimigarray_vs (mi);
							imig[mi].b = k;
							imig[mi].e = k;
							imig[mi].wp.n = 1;
						}
						else
							imig[mi].wp.n++;
						
						imig[mi].wp.n += modeloptions[SINGLEMIGRATIONBOTHDIRECTIONS]; // can't also have ONEMIGRATIONPARAMETER
						
						// if parameters are set by period
						if (modeloptions[PARAMETERSBYPERIOD] && !modeloptions[ONEMIGRATIONPARAMETER])
						{
							sprintf (imig[mi].str, "m%d,", k);
							if (modeloptions[SINGLEMIGRATIONBOTHDIRECTIONS])
								sprintf (temps, "%d<>%d", C[ci]->plist[k][i], C[ci]->plist[k][j]);
							else
								sprintf (temps, "%d>%d", C[ci]->plist[k][i], C[ci]->plist[k][j]);
							strcat (imig[mi].str, temps);
						}
						else
						{
							kk = k + 1;
							while (kk < lastperiodnumber && ISELEMENT (C[ci]->plist[k][i], C[ci]->periodset[kk])
								&& ISELEMENT (C[ci]->plist[k][j], C[ci]->periodset[kk]))
							{
								if (!modeloptions[ONEMIGRATIONPARAMETER])
									imig[mi].e = kk;
								kk++;
								imig[mi].wp.n++;
								imig[mi].wp.n += modeloptions[SINGLEMIGRATIONBOTHDIRECTIONS]; // can't also have ONEMIGRATIONPARAMETER
							}
							if (!modeloptions[ONEMIGRATIONPARAMETER])
							{
								if (modeloptions[SINGLEMIGRATIONBOTHDIRECTIONS])  // can't also have ONEMIGRATIONPARAMETER
									sprintf (imig[mi].str, "m%d<>%d", C[ci]->plist[k][i], C[ci]->plist[k][j]);
								else
									sprintf (imig[mi].str, "m%d>%d", C[ci]->plist[k][i], C[ci]->plist[k][j]);
							}
						}
						if (!modeloptions[SINGLEMIGRATIONBOTHDIRECTIONS])
						{
							if (!modeloptions[ONEMIGRATIONPARAMETER])
							{
							  mi++;
							  if (calcoptions[USEPRIORFILE])
							  {
								// VS add a for loop that will initialize the parameters for the different groups of loci
								for(gp = 0; gp < nbgroupsloci_mig; gp++)
								{
									imig[mi].pr[gp].max = relmig[gp]*mprior_fromfile[C[ci]->plist[k][i]][C[ci]->plist[k][j]];
									if (imig[mi].pr[gp].max <= MINPARAMVAL)
									  IM_err(IMERR_PRIORFILEVALS,"migration rate set too low from %d to %d",i,j);
									imig[mi].pr[gp].min = 0;
									imig[mi].pr[gp].mean = (imig[mi].pr[gp].max-imig[mi].pr[gp].min)/2;
								}
							  }
							  // VS
							  // setimigarray (mi);
							  setimigarray_vs (mi);
							  imig[mi].b = k;
							  imig[mi].e = k;
							  imig[mi].wp.n = 1;
							}
							else
							  imig[mi].wp.n++;

							if (modeloptions[PARAMETERSBYPERIOD] && !modeloptions[ONEMIGRATIONPARAMETER])
							{
							  sprintf (imig[mi].str, "m%d,%d>%d", k, C[ci]->plist[k][j], C[ci]->plist[k][i]);
							}
							else
							{
							  kk = k + 1;
							  while (kk < lastperiodnumber && ISELEMENT (C[ci]->plist[k][i], C[ci]->periodset[kk])
									 && ISELEMENT (C[ci]->plist[k][j], C[ci]->periodset[kk])) {
										 if (!modeloptions[ONEMIGRATIONPARAMETER]) {
											imig[mi].e = kk;
										}
									kk++;
									imig[mi].wp.n++;
							  }
							  if (!modeloptions[ONEMIGRATIONPARAMETER])
								sprintf (imig[mi].str, "m%d>%d", C[ci]->plist[k][j], C[ci]->plist[k][i]);
							} // end of else from if(modeloptions[PARAMETERSBYPERIOD] && !modeloptions[ONEMIGRATIONPARAMETER]) 
						} // end if(!modeloptions[SINGLEMIGRATIONBOTHDIRECTIONS])
					} // end if mcheck
				} // end if tricky condition
		  } // end for loop j
		} // end of for loop i
	} // end of for loop k

	/* initialize the parts of wp */
    for (mi = 0; mi < nummigrateparams; mi++)
    {
      if (imig[mi].wp.n > 0)
      {
        imig[mi].wp.p = malloc (imig[mi].wp.n * sizeof (int));
        imig[mi].wp.r = malloc (imig[mi].wp.n * sizeof (int));
        imig[mi].wp.c = malloc (imig[mi].wp.n * sizeof (int));
      }
    }

    /* now do another loop, much like the previous two. this time finish building the wp parts of imig[]
      using wp.n that was determined in the previous loop */
	if (modeloptions[ONEMIGRATIONPARAMETER])
    {
      mi = 0; // Why is mi=0 in this case?
      ni=-1;
    }
    else
      mi = -1;

	// For loops 
    for (k = 0; k < lastperiodnumber; k++){
		for (i = 0; i < npops - k - 1; i++) {
		  for (j = i + 1; j < npops - k; j++) {
			  if ((k == 0) ||
				  (!modeloptions[MIGRATIONBETWEENSAMPLED] &&
				  (modeloptions[PARAMETERSBYPERIOD] ||
				  ((!ISELEMENT(C[ci]->plist[k][i], C[ci]->periodset[k - 1])) ||
				   (!ISELEMENT(C[ci]->plist[k][j], C[ci]->periodset[k - 1]))))))
			  {
				// if no migration between non sister pops
				if (modeloptions[NOMIGBETWEENNONSISTERS])
				  mcheck = checkaresis (k, i, j);
				else
				{
				   if (!calcoptions[USEPRIORFILE] || (calcoptions[USEPRIORFILE] && mprior_fromfile[C[ci]->plist[k][i]][C[ci]->plist[k][j]] > MINPARAMVAL))
					mcheck = 1;
				   else
					 mcheck = 0;
				}
				if (mcheck)
				{
				  mi++;
				  ni = 0;
				  imig[mi].wp.p[ni] = k;
				  imig[mi].wp.r[ni] = i;
				  imig[mi].wp.c[ni] = j;
				  assert (ni < imig[mi].wp.n);
				  if (modeloptions[SINGLEMIGRATIONBOTHDIRECTIONS])  // can't also have ONEMIGRATIONPARAMETER
				  {
					ni++;
					imig[mi].wp.p[ni] = k;
					imig[mi].wp.r[ni] = j;
					imig[mi].wp.c[ni] = i;
					assert (ni < imig[mi].wp.n);
				  }
				  if (!modeloptions[PARAMETERSBYPERIOD])
				  {
					kk = k + 1;
					while (kk < lastperiodnumber &&
						   ISELEMENT (C[ci]->plist[k][i],
									  C[ci]->periodset[kk])
						   && ISELEMENT (C[ci]->plist[k][j],
										 C[ci]->periodset[kk]))
					{
					  ii = 0;
					  while (C[0]->plist[kk][ii] != C[0]->plist[k][i])
						ii++;
					  assert (ii < npops - kk);
					  jj = 0;
					  while (C[0]->plist[kk][jj] != C[0]->plist[k][j])
						jj++;
					  assert (jj < npops - kk);
					  ni++;
					  imig[mi].wp.p[ni] = kk;
					  imig[mi].wp.r[ni] = ii;
					  imig[mi].wp.c[ni] = jj;
					  assert (ni < imig[mi].wp.n);
					  if (modeloptions[SINGLEMIGRATIONBOTHDIRECTIONS])
					  {
						ni++;
						imig[mi].wp.p[ni] = kk;
						imig[mi].wp.r[ni] = jj;
						imig[mi].wp.c[ni] = ii;
						assert (ni < imig[mi].wp.n);
					  }
					  kk++;
					}
				  }
				  if (!modeloptions[SINGLEMIGRATIONBOTHDIRECTIONS])
				  {
					if (!modeloptions[ONEMIGRATIONPARAMETER])
					{
					  mi++;
					  ni = 0;
					}
					else
					  ni++;
					imig[mi].wp.p[ni] = k;
					imig[mi].wp.r[ni] = j;
					imig[mi].wp.c[ni] = i;
					assert (ni < imig[mi].wp.n);
					if (!modeloptions[PARAMETERSBYPERIOD])
					{
					  kk = k + 1;
					  while (kk < lastperiodnumber &&
							 ISELEMENT (C[ci]->plist[k][i],
										C[ci]->periodset[kk])
							 && ISELEMENT (C[ci]->plist[k][j],
										   C[ci]->periodset[kk]))
					  {
						ii = 0;
						while (C[0]->plist[kk][ii] != C[0]->plist[k][i])
						  ii++;
						assert (ii < npops - kk);
						jj = 0;
						while (C[0]->plist[kk][jj] != C[0]->plist[k][j])
						  jj++;
						assert (jj < npops - kk);
						ni++;
						imig[mi].wp.p[ni] = kk;
						imig[mi].wp.r[ni] = jj;
						imig[mi].wp.c[ni] = ii;
						kk++;
						assert (ni < imig[mi].wp.n);
					  }
					}
				  }
				}
			 } // end of tricky condition
		  } // end of for loop j
        } // end of for loop i
    } // end of for loop k                           
  } // initialize wp

	// VS
	// CHECK
	printf("\n\nLOOK INTO IMIG\n");
	for(i = 0 ; i < nummigrateparams; i++) {
		// Go through all the relative parameters in the input file
		printf("\nimig[%i] %s", i, imig[i].str);
		printf("\nimig[%i] first period %i", i, imig[i].b);
		printf("\nimig[%i] last period %i", i, imig[i].e);

		for(gp = 0; gp < nbgroupsloci_mig; gp++) {
			printf("\n imig[%i].pr[%i].min is %g", i, gp, imig[i].pr[gp].min);
			printf("\n imig[%i].pr[%i].max is %g", i, gp, imig[i].pr[gp].max);

			for(m = 0; m < imig[i].wp.n; m++) {
				printf("\n imig[%i].wp.p[%i] is %i", i, m, imig[i].wp.p[m]);
				printf("\n imig[%i].wp.r[%i] is %i", i, m, imig[i].wp.r[m]);
			}
		}
	}

/*
	type     |  # values  |  cumulative total at end
    cc	       numpopsizeparams   numpopsizeparams
	fc	       numpopsizeparams   2*numpopsizeparams
	hcc	       numpopsizeparams   3*numpopsizeparams
	mc         nummigrateparams   3*numpopsizeparams + nummigrateparams
	fm         nummigrateparams   3*numpopsizeparams + 2*nummigrateparams
	qintegrate numpopsizeparams   4*numpopsizeparams + 2*nummigrateparams
	mintegrate nummigrateparams   4*numpopsizeparams + 3*nummigrateparams
	sintegrate      1        4*numpopsizeparams + 3*nummigrateparams +  1
	pdg             1        4*numpopsizeparams + 3*nummigrateparams +  2
	plg             1        4*numpopsizeparams + 3*nummigrateparams +  3
	probg           1        4*numpopsizeparams + 3*nummigrateparams +  4
	t          numsplittimes 4*numpopsizeparams + 3*nummigrateparams +  numsplittimes + 4
*/

  /* set values of position markers for gsampinf */

  if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 1)
  {
    assert (numpopsizeparams == npops);
    if (modeloptions[NOMIGRATION] == 0)
    {
      assert (nummigrateparams == npops * (npops - 1));
    }
  }

  gsamp_ccp = 0;
  gsamp_fcp = gsamp_ccp + numpopsizeparams;
  gsamp_hccp = gsamp_fcp + numpopsizeparams;
  gsamp_mcp = gsamp_hccp + numpopsizeparams;
  gsamp_fmp = gsamp_mcp + nummigrateparams;
  gsamp_qip = gsamp_fmp + nummigrateparams;
  gsamp_mip = gsamp_qip + numpopsizeparams;
  gsamp_pdgp = gsamp_mip + nummigrateparams;
  gsamp_probgp =  gsamp_pdgp + 1;
  gsamp_tp = gsamp_probgp + 1;

  /*// to redo the 48 locus IMa runs of 2006
itheta[0].pr.max = 3.82;
itheta[1].pr.max = 0.87;
itheta[2].pr.max = 1.53;
imig[0].pr.max = 6.0;
imig[1].pr.max = 3.0; */
}											/* setup_iparams */

/* nomigrationchecklist is three arrays (p for period, r for row, and c for column), each of a set length (nomigrationchecklist.n )

nomigrationchecklist is used to check gweight->mc[][][]
gweight->mc[k][i][j] is the number of migration events in the genealogies in period k 
from population C[ci]->plist[k][i] to population C[ci]->plist[k][j]

So the ith element of the arrays in nomigrationchecklist is used to check the value of 
gweight->mc[nomigrationchecklist.p[i]][nomigrationchecklist.r[i]][nomigrationchecklist.c[i]]
(this is done in update_gree_common.c).  If that value if not zero the update gets rejected. 

building nomigrationchecklist requires identifying which pairs of populations, 
C[ci]->plist[k][i] and C[ci]->plist[k][j], in each period k,  should not be
exchanging migrants.

These pairs are those that meet any of the following:
1. modeloptions[MIGRATIONBETWEENSAMPLED]==1 and the populations are not both sampled populations
2. modeloptions[NOMIGBETWEENNONSISTERS]==1  and the populations are not sisters
3. they have had their priors set to zero in a priorfile

These correspond to the following:
1. modeloptions[MIGRATIONBETWEENSAMPLED]==1 && C[ci]->plist[k][i] >= npops && C[ci]->plist[k][j] >= npops
2. modeloptions[NOMIGBETWEENNONSISTERS]==1 && checkaresis (k, i, j)==1
3. calcoptions[USEPRIORFILE]==1 && mprior_fromfile[C[ci]->plist[k][i]][C[ci]->plist[k][j]] <= MINPARAMVAL


To set this up we must go through two loops,  first to cound how many (i.e. set nomigrationchecklist.n)
Then after initializing the arrays in nomigrationchecklist, we must do another similar loop to fill them up. 
*/

/* this should not be needed if modeloptions[NOMIGRATION]==1 */
void
set_nomigrationchecklist ()
{
  int n, i, j, k, ci = 0;

  nomigrationchecklist.n = 0;
  nomigrationchecklist.p = NULL;
  nomigrationchecklist.r = NULL;
  nomigrationchecklist.c = NULL;
  for (k = 0; k < lastperiodnumber; k++)
    for (i = 0; i < npops - k - 1; i++)
      for (j = i + 1; j < npops - k; j++)
      {
        if ((modeloptions[MIGRATIONBETWEENSAMPLED]==1 && (C[ci]->plist[k][i] >= npops || C[ci]->plist[k][j] >= npops)) ||
            (modeloptions[NOMIGBETWEENNONSISTERS]==1 && checkaresis (k, i, j)==0) ||
            (calcoptions[USEPRIORFILE]==1 && mprior_fromfile[C[ci]->plist[k][i]][C[ci]->plist[k][j]] <= MINPARAMVAL))
        {
          nomigrationchecklist.n+=2;
        }
      }
    if (nomigrationchecklist.n > 0)
    {
      nomigrationchecklist.p = malloc (nomigrationchecklist.n * sizeof (int));
      nomigrationchecklist.r = malloc (nomigrationchecklist.n * sizeof (int));
      nomigrationchecklist.c = malloc (nomigrationchecklist.n * sizeof (int));
    }
    for (n = -1, k = 0; k < lastperiodnumber; k++)
      for (i = 0; i < npops - k - 1; i++)
        for (j = i + 1; j < npops - k; j++)
        {
          if ((modeloptions[MIGRATIONBETWEENSAMPLED]==1 && (C[ci]->plist[k][i] >= npops || C[ci]->plist[k][j] >= npops)) ||
            (modeloptions[NOMIGBETWEENNONSISTERS]==1 && checkaresis (k, i, j)==0) ||
            (calcoptions[USEPRIORFILE]==1 && mprior_fromfile[C[ci]->plist[k][i]][C[ci]->plist[k][j]] <= MINPARAMVAL))
          {
              n++;
              nomigrationchecklist.p[n] = k;
              nomigrationchecklist.r[n] = i;
              nomigrationchecklist.c[n] = j;
              n++;
              nomigrationchecklist.p[n] = k;
              nomigrationchecklist.r[n] = j;
              nomigrationchecklist.c[n] = i;
          }
        }
  assert(nomigrationchecklist.n == n+1);
  return;
}                               // set_nomigrationchecklist

#define MAXPRIORSCALETRY 10000000       // max number of times to try getting starting mutation rates compatible with priors
void
setuinfo (double summut)
{
  int i, li, ui, uj, k;
  int numuprior = 0;
  double priorscaletry = 0;
  double U, r;
  int doneu, rcheck, numupair, *upriorpairlist1, *upriorpairlist2;
  double prodcheck, maxr, newr, d;

/*  ul contains a list of locations by locus and within locus of all mutation rate scalars
    uii is the position in that list
	in set_mcparam_values() the priors on mutation rate scalars are set to standard pos (max) and neg (min) values (they are on a log scale)
	   and are stored in the mcinf.pr   (e.g. C[0]->G[0].u[0].mcinf.pr)
	if the input file has prior ranges on mutation yets,  these are read into uperyear.pr  (e.g. C[ci]->G[li].u[ui].uperyear.pr)
	The ratios of these uperyear values among loci must be taken to reset the values of mcinf.pr   
*/
  for (i = 0, li = 0; li < nloci; li++)
    for (ui = 0; ui < L[li].nlinked; ui++)
    {
      ul[i].l = li;
      ul[i].u = ui;
      L[li].uii[ui] = i;
      i++;
    }

  /* reset the mutation rate values to have geometric mean of 1 */
  for (prodcheck = 1, ui = 0; ui < nurates; ui++)
  {
    C[0]->G[ul[ui].l].uvals[ul[ui].u] =
      exp (uvals[ul[ui].l][ul[ui].u] - summut / nurates);
    prodcheck *= C[0]->G[ul[ui].l].uvals[ul[ui].u];
  }
  if (fabs (log (prodcheck)) > 1e-5)
    IM_err(IMERR_MUTSCALEPRODUCTFAIL,"product of mutation scalars not close to or equal to 1: %lf",prodcheck); 
  for (ui = 0; ui < nurates; ui++)
  {
    if (L[ul[ui].l].uperyear_prior[ul[ui].u].min != 0)
      numuprior++;
  }

  if (calcoptions[MUTATIONPRIORRANGE])
  {
    assert (numuprior > 1);
    numupair = numuprior * (numuprior - 1) / 2;
    upriorpairlist1 = calloc ((size_t) numupair, sizeof (int));
    upriorpairlist2 = calloc ((size_t) numupair, sizeof (int));
    k = 0;
    for (ui = 0; ui < nurates - 1; ui++)
    {
      for (uj = ui + 1; uj < nurates; uj++)
      {
        if (L[ul[ui].l].uperyear_prior[ul[ui].u].min != 0
            && L[ul[uj].l].uperyear_prior[ul[uj].u].min != 0)

        {
          upriorpairlist1[k] = ui;
          upriorpairlist2[k] = uj;
          k++;
        }
      }
    }

    /* urri[i][j] has a 0 if neither i nor j has a prior. 1 if i has a prior and j does not,  -1 if i does not have a pior and j does  */
    for (ui = 0; ui < nurates; ui++)
    {
      for (uj = 0; uj < nurates; uj++)
      {
        urri[ui][uj] = 0;
        if (ui != uj
            && L[ul[ui].l].uperyear_prior[ul[ui].u].min != 0
            && L[ul[uj].l].uperyear_prior[ul[uj].u].min != 0)
        {
          urrlow[ui][uj] =
            log (L[ul[ui].l].uperyear_prior[ul[ui].u].min /
                 L[ul[uj].l].uperyear_prior[ul[uj].u].max);
          urrhi[ui][uj] =
            log (L[ul[ui].l].uperyear_prior[ul[ui].u].max /
                 L[ul[uj].l].uperyear_prior[ul[uj].u].min);
          urri[ui][uj] = urri[uj][ui] = 2;
        }
        else
        {
          if (ui != uj
              && L[ul[ui].l].uperyear_prior[ul[ui].u].min != 0
              && L[ul[uj].l].uperyear_prior[ul[uj].u].min == 0)
            urri[ui][uj] = 1;
          if (ui != uj
              && L[ul[uj].l].uperyear_prior[ul[uj].u].min != 0
              && L[ul[ui].l].uperyear_prior[ul[ui].u].min == 0)
            urri[ui][uj] = -1;
        }
      }
    }
    /* need to set all of the uscalers so that their ratios are in the ranges defined by the priors */
    /* it is possible that suiteable sets of scalars will not be able to be found */
    maxr = 3 * L[0].u_rec[0].pr.max;

    do
    {
      doneu = 1;
      for (i = 0; i < numupair; i++)

      {
        ui = upriorpairlist1[i];
        uj = upriorpairlist2[i];
        r =
          log (C[0]->G[ul[ui].l].uvals[ul[ui].u] /
               C[0]->G[ul[uj].l].uvals[ul[uj].u]);
        rcheck = (r >= urrlow[ui][uj] && r <= urrhi[ui][uj]);
        doneu = doneu && rcheck;
        if (!rcheck)
        {
          do
          {
            U = uniform ();
            if (U > 0.5)
              newr = r + maxr * (2.0 * U - 1.0);

            else
              newr = r - maxr * U * 2.0;
            if (newr > maxr)
              newr = 2.0 * maxr - newr;

            else if (newr < -maxr)
              newr = 2.0 * (-maxr) - newr;
            d = exp ((newr - r) / 2);
          } while ((newr <= urrlow[ui][uj] || newr >= urrhi[ui][uj]));
          C[0]->G[ul[ui].l].uvals[ul[ui].u] *= d;
          C[0]->G[ul[uj].l].uvals[ul[uj].u] /= d;
        }
      }
      priorscaletry++;
    } while (!doneu && priorscaletry < MAXPRIORSCALETRY);

    if (priorscaletry >= MAXPRIORSCALETRY)
      IM_err(IMERR_MUTSCALARPRIORRANGEFAIL,"More than %d failed attempts at finding a valid set of mutation scalars, prior ranges probably too restrictive",MAXPRIORSCALETRY);
    for (prodcheck = 1, ui = 0; ui < nurates; ui++)
      prodcheck *= C[0]->G[ul[ui].l].uvals[ul[ui].u];
      //prodcheck *= C[0]->G[ul[ui].l].u[ul[ui].u].uperyear.val;

    if (fabs (log (prodcheck)) > 1e-5)
    {
      IM_err(IMERR_MUTSCALEPRODUCTFAIL,"product of mutation scalars not close to or equal to 1: %lf",prodcheck); 
    }
    XFREE (upriorpairlist1);
    XFREE (upriorpairlist2);
  }
}                               /* setuinfo */

double
set_uvals (void)
  /* set mutation rate scalars */
  /* get a relative mutation rate for a locus by summing up estimates of 4Nu for all the populations
     the geometric mean of these relative values are then used to set the starting mutation rates */
{
  int j, k, li, ui;
  int b, e;
  double w, w2;
  double temp, setutemp, sumq;
  int npopstemp;
  

  if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
  {
    npopstemp = npops + 1;
  }
  else
  {
    npopstemp = npops;
  }


  if (assignmentoptions[POPULATIONASSIGNMENT] == 0)
  {
    /*************************************************************************/
  setutemp = 0;
  for (li = 0; li < nloci; li++)
  {
    if (L[li].model == INFINITESITES
        || L[li].model == JOINT_IS_SW 
        || L[li].model == HKY)
    {
      sumq = 0;
      b = 0;
      e = L[li].samppop[0] - 1;
      for (k = 0; k < npopstemp; k++)
      {
        if (e >= b)
        {
          for (j = 1, w = 0.0; j <= e; j++)
            w += 1 / (double) j;
          if (w == 0)
            w = 1;
          if (L[li].model == HKY)
            temp = numvarHKY (li, b, e) / w;
          else
            temp = numsitesIS[li][k] / w;
        }
        else
        {
          temp = 0;
        }
        if (modeloptions[ANCESTRALPOPSIZESSETTOLARGESTSAMPLED] && k<npops)
        {
          varrank[k] += temp;
        }
        b = e + 1;
        if (k + 1 == npopstemp)
        {
          e = e + L[li].numgenesunknown;
        }
        else
        {
          e = e + L[li].samppop[k + 1];
        }
        sumq += temp;

      }
      sumq = DMAX (0.1, sumq);  // nominal low values in case of zero variation 
      uvals[li][0] = log (sumq);
      setutemp += uvals[li][0];
    }

    if (L[li].model == STEPWISE || L[li].model == JOINT_IS_SW)
    {
      somestepwise = 1;
      if (L[li].model == JOINT_IS_SW)
        ui = 1;
      else
        ui = 0;

      for (; ui < L[li].nlinked; ui++)
      {
        b = 0;
        e = L[li].samppop[0] - 1;
        sumq = 0;
        for (k = 0; k < npopstemp; k++)
        {
          if (e >= b)
          {
            for (j = b, w = 0, w2 = 0; j <= e; j++)
            {
              w += L[li].A[ui][j];
              w2 += SQR ((double) L[li].A[ui][j]);
            }
            if (k == npops)
            {
              temp =
                (double) 2 *(w2 -
                             SQR (w) / L[li].numgenesunknown) /
                (L[li].numgenesunknown);
            }
            else
            {
              temp =
                (double) 2 *(w2 -
                             SQR (w) / L[li].samppop[k]) / (L[li].samppop[k]);
            }
            sumq += temp;
          }
          else
          {
            temp = 0;
          }
          if (modeloptions[ANCESTRALPOPSIZESSETTOLARGESTSAMPLED] && k<npops)
          {
            varrank[k] += temp;
          }
          b = e + 1;
          if (k + 1 == npopstemp)
          {
            e = e + L[li].numgenesunknown;
          }
          else
          {
            e = e + L[li].samppop[k + 1];
          }
        }
        sumq = DMAX (0.1, sumq);        // nominal low values in case of zero variation 
        uvals[li][ui] = log (sumq);
        setutemp += uvals[li][ui];
      }
    }
  }

  // geomeanvar is used to set the prior on thetas, in setup_iparams(), depending on command line options 
  /* sum is the sum of the logs, across loci, of the maximal amount of variation found among the sampled populations */
  /* the scalar of thetaprior is multiplied times the geometric mean of the maximal estimates of variation across loci */
  geomeanvar = exp (setutemp / nurates);

    /*************************************************************************/
  }
  /* SANGCHUL: Sat May 16 16:48:10 EDT 2009
   * For assignment option, starting u values are all 1.
   */
  else if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
  {
    for (li = 0; li < nloci; li++)
    {
      if (L[li].model == INFINITESITES
          || L[li].model == JOINT_IS_SW 
          || L[li].model == HKY)
      {
        uvals[li][0] = 0.0;
      }

      if (L[li].model == STEPWISE 
          || L[li].model == JOINT_IS_SW)
      {
        somestepwise = 1;
        if (L[li].model == JOINT_IS_SW)
        {
          ui = 1;
        }
        else
        {
          ui = 0;
        }
        for (; ui < L[li].nlinked; ui++)
        {
          uvals[li][ui] = 0.0;
        }
      }
    }
    geomeanvar = 1.0;
    setutemp = 0.0;
  }
  return setutemp;
}                               /*set_uvals */

void
set_x (struct value_record *v, int isint)
{
  int i;
  for (i = 0; i < GRIDSIZE; i++)
  {
    if (v->do_logplot)
    {
      v->xy[i].x = exp (v->plotrange.min + (i + 0.5) * v->plotrescale); // the same scale is used for all mutation rate scalars 
    }
    else
    {
      if (isint)
        v->xy[i].x = (double) i;
      else
        v->xy[i].x = v->plotrange.min + ((i + 0.5) * (v->plotrange.max * v->plotrescale - v->plotrange.min)) / GRIDSIZE;
    }
  }
}                               // set_x 

void
init_value_record (struct value_record *v, int isint)
{
  if (v->do_xyplot)
  {
    v->xy = calloc (GRIDSIZE, sizeof (struct plotpoint));
    set_x (v, isint);
  }

  if (v->do_trend)
    v->trend = calloc (TRENDDIM, sizeof (double));

  v->beforemin = v->aftermax = 0;
}                               // init_value_Record

void
init_g_rec (int li)
{
  int num_g_update_types = IM_UPDATE_GENEALOGY_NUMBER;
  L[li].g_rec = malloc (sizeof (struct chainstate_record_updates_and_values));
  sprintf (L[li].g_rec->str, "gtree_%d", li);
  L[li].g_rec->num_uptypes = num_g_update_types;
  L[li].g_rec->upnames = malloc (L[li].g_rec->num_uptypes * sizeof (strnl));
  sprintf (L[li].g_rec->upnames[IM_UPDATE_GENEALOGY_ANY],      "branch     ");
  sprintf (L[li].g_rec->upnames[IM_UPDATE_GENEALOGY_TOPOLOGY], "topology   ");
  sprintf (L[li].g_rec->upnames[IM_UPDATE_GENEALOGY_TMRCA],    "tmrca      ");
  L[li].g_rec->upinf =
    calloc ((size_t) L[li].g_rec->num_uptypes, sizeof (struct update_rate_calc));
  L[li].g_rec->num_vals = 1;
  L[li].g_rec->v =
    malloc (L[li].g_rec->num_vals * sizeof (struct value_record));
  sprintf (L[li].g_rec->v->str, "g%d_tmrca", li);
  sprintf (L[li].g_rec->v->strshort, "tmrca%d", li);
  L[li].g_rec->v->plotrange.min = 0;
  L[li].g_rec->v->do_autoc = 1;
  L[li].g_rec->v->do_xyplot = outputoptions[PRINTTMRCA];
  L[li].g_rec->v->do_trend = 1;
  L[li].g_rec->v->plotrescale = 1.0;
  L[li].g_rec->v->do_logplot = 0;
  

}                               // init_g_rec

void
init_a_rec (int li)
{
  int num_a_update_types = IM_UPDATE_ASSIGNMENT_NUMBER;

  L[li].a_rec = malloc (sizeof (struct chainstate_record_updates_and_values));
  sprintf (L[li].a_rec->str, "asn_%d", li);
  L[li].a_rec->num_uptypes = num_a_update_types;
  L[li].a_rec->upnames = malloc (num_a_update_types * sizeof (strnl));

  sprintf (L[li].a_rec->upnames[IM_UPDATE_ASSIGNMENT_RELABEL], "relabel");
  sprintf (L[li].a_rec->upnames[IM_UPDATE_ASSIGNMENT_BF], "bf");

  L[li].a_rec->upinf = calloc ((size_t) L[li].a_rec->num_uptypes, sizeof (struct update_rate_calc));
  L[li].a_rec->num_vals = 1;
  L[li].a_rec->v = malloc (L[li].a_rec->num_vals * sizeof (struct value_record));
  sprintf (L[li].a_rec->v->str, "a%2d_asn", li);
  strncpy (L[li].a_rec->v->strshort, L[li].a_rec->v->str, PARAMSTRLENSHORT - 1);
  L[li].a_rec->v->plotrange.min = 0;
  L[li].a_rec->v->plotrange.max = 10;
  L[li].a_rec->v->do_autoc = 1;
  L[li].a_rec->v->do_xyplot = 0;
  L[li].a_rec->v->do_trend = 1;
  L[li].a_rec->v->plotrescale = 1.0;
  L[li].a_rec->v->do_logplot = 0;
  init_value_record (L[li].a_rec->v, 0);
  return;
}                               /* init_a_rec */

void
init_a_rec_multichain ()
{
  int ci;
  Cupinf = malloc (numchains * sizeof (im_chainstate_updateonly));

  for (ci = 0; ci < numchains; ci++)
  {
    sprintf (Cupinf[ci].str, "ac_%d", ci);
    Cupinf[ci].num_uptypes = IM_UPDATE_ASSIGNMENT_NUMBER;
    Cupinf[ci].upnames = malloc (IM_UPDATE_ASSIGNMENT_NUMBER * sizeof (strnl));
    Cupinf[ci].upinf = calloc (IM_UPDATE_ASSIGNMENT_NUMBER, sizeof (struct update_rate_calc));
    sprintf (Cupinf[ci].upnames[IM_UPDATE_ASSIGNMENT_RELABEL], "relabel");
    sprintf (Cupinf[ci].upnames[IM_UPDATE_ASSIGNMENT_BF], "bf");
  }
  return;
}


void
init_lpgpd_v (void)
{
  lpgpd_v = malloc (sizeof (struct value_record));
  sprintf (lpgpd_v->str, "Log[P(G)+P(D|G)]");
  sprintf (lpgpd_v->strshort, "Log[P]");
  lpgpd_v->plotrange.min = lpgpd_v->plotrange.max = 0;
  lpgpd_v->do_xyplot = 0;
  lpgpd_v->do_logplot = 0;
  lpgpd_v->do_trend = 1;
  lpgpd_v->do_autoc = 1;
  lpgpd_v->plotrescale = 1.0;
  init_value_record (lpgpd_v, 0);
  return;
}                               // init_lpgpd_rec

void fillmrec(int j, int i, int pi, int pj)
{
  char stemp[6];
  sprintf(stemp,"%d>%d",pi,pj);
  strcpy (migration_counts_times[j][i].str, stemp);
  strcat (migration_counts_times[j][i].str, "_t");
  migration_counts_times[j][i].plotrange = T[numsplittimes - 1].pr; // same range as for splitting times
  strcpy (migration_counts_times[j][i+1].str,stemp);
  strcat (migration_counts_times[j][i+1].str, "_#");
  migration_counts_times[j][i+1].plotrange.max = (double) GRIDSIZE - 1;       // used for counts so the grid position is the count #
  migration_counts_times[j][i+1].plotrange.min = 0;
}

void
init_migration_counts_times (void)
{
  int i, j, k, pi, pj, ji, jj, numhists, mrows, nummigdirs;
  mrows = nloci + (nloci > 1); 
  nummigdirs = 2*(npops-1)*(npops-1);
  numhists = 2 * nummigdirs;
  migration_counts_times = malloc (mrows * sizeof (struct value_record *));
  for (j = 0; j < mrows; j++)
    migration_counts_times[j] =
      malloc (numhists * sizeof (struct value_record));
  for (j = 0; j < mrows; j++)
  {
    i=0;
    for (k= 0; k < lastperiodnumber; k++)
    {
      for (ji=0; ji< npops-k;ji++) 
      {
          pi = C[0]->plist[k][ji];
          if (k==0)
          {
            for (jj=ji+1; jj< npops-k;jj++) 
            {
              pj = C[0]->plist[k][jj];
              if (pi != pj)
              {
                fillmrec(j,i,pi,pj);
                fillmrec(j,i+2,pj,pi);
                i+=4;
              }
            }
          }
          else
          {
            pj = C[0]->addpop[k];
            if (pi != pj)
            {
              fillmrec(j,i,pi,pj);
              fillmrec(j,i+2,pj,pi);
              i+=4;
            }
          }
      }
    }
    assert(i==numhists);
  }
  for (j = 0; j < mrows; j++)
    for (i = 0; i < numhists; i++)
    {
      migration_counts_times[j][i].do_xyplot = 1;
      migration_counts_times[j][i].do_logplot = 0;
      migration_counts_times[j][i].do_trend = 0;
      migration_counts_times[j][i].do_autoc = 0;
      migration_counts_times[j][i].plotrescale = 1.0;
      if (ODD (i))
        init_value_record (&migration_counts_times[j][i], 1);
      else
        init_value_record (&migration_counts_times[j][i], 0);
    }
}                               //init_migration_counts_times



void
init_mutation_scalar_rec (int li)
{
  int i, ui, ai;
  int num_u_update_types = 1;   // number of different types of updates for mutation rate scalars and kappa values
  int num_A_update_types = 1;   // number of different types of updates for A states 
  double uscale;

  // initialize u_rec
  L[li].u_rec =
    malloc (L[li].nlinked *
            sizeof (struct chainstate_record_updates_and_values));
  // set priors and windows
  for (ui = 0; ui < L[li].nlinked; ui++)
  {
    L[li].u_rec[ui].pr.max = log (UMAX);        // 10000
    L[li].u_rec[ui].pr.min = -log (UMAX);
    L[li].u_rec[ui].win = L[li].u_rec[ui].pr.max / nloci;
  }
  uscale = (2 * log (UMAX)) / GRIDSIZE; // 10000
  // name the mutation rate scalars     
  if (L[li].model == STEPWISE)
  {
    for (ai = 0; ai < L[li].nlinked; ai++)
      sprintf (L[li].u_rec[ai].str, "%dSW%d", li, ai);
  }
  else
  {
    sprintf (L[li].u_rec[0].str, "%du ", li);
    if (L[li].model == JOINT_IS_SW)
      for (ai = 1; ai < L[li].nlinked; ai++)
        sprintf (L[li].u_rec[ai].str, "%dSW%d", li, ai);
  }
  for (ui = 0; ui < L[li].nlinked; ui++)
  {
    L[li].u_rec[ui].num_uptypes = num_u_update_types;
    L[li].u_rec[ui].upnames =
      malloc (L[li].u_rec[ui].num_uptypes * sizeof (strnl));
    for (i = 0; i < num_u_update_types; i++)
    {
      sprintf (L[li].u_rec[ui].upnames[i], "scalar update");
    }
    L[li].u_rec[ui].upinf = calloc ((size_t) L[li].u_rec[ui].num_uptypes, sizeof (struct update_rate_calc));
    L[li].u_rec[ui].num_vals = 1;
    L[li].u_rec[ui].v = malloc (L[li].u_rec[ui].num_vals * sizeof (struct value_record));
    for (i = 0; i < L[li].u_rec[ui].num_vals; i++)
    {
      strcpy (L[li].u_rec[ui].v[i].str, L[li].u_rec[ui].str);   // will this ever get used ? 
      L[li].u_rec[ui].v[i].do_xyplot = 1;
      L[li].u_rec[ui].v[i].plotrescale = uscale;
      L[li].u_rec[ui].v[i].do_logplot = 1;
      L[li].u_rec[ui].v[i].do_autoc = 0;
      L[li].u_rec[ui].v[i].do_trend = 1;
      L[li].u_rec[ui].v[i].plotrange = L[li].u_rec[ui].pr;
      init_value_record (&L[li].u_rec[ui].v[i], 0);
    }
  }

  // do kappa_rec 
  if (L[li].model == HKY)
  {
    L[li].kappa_rec = malloc (sizeof (struct chainstate_record_updates_and_values));
    L[li].kappa_rec->pr.max = KAPPAMAX;
    L[li].kappa_rec->pr.min = 0.0;
    L[li].kappa_rec->win = 2.0;
    sprintf (L[li].kappa_rec->str, "%d_Ka", li);
    L[li].kappa_rec->num_uptypes = num_u_update_types;
    L[li].kappa_rec->upnames = malloc (num_u_update_types * sizeof (strnl));
    for (i = 0; i < num_u_update_types; i++)
    {
      sprintf (L[li].kappa_rec->upnames[i], "kappa update");
    }
    L[li].kappa_rec->upinf = calloc ((size_t) num_u_update_types, sizeof (struct update_rate_calc));
    L[li].kappa_rec->num_vals = 1;
    L[li].kappa_rec->v = malloc (L[li].kappa_rec->num_vals * sizeof (struct value_record));
    strcpy (L[li].kappa_rec->v->str, L[li].kappa_rec->str);
    strncpy (L[li].kappa_rec->v->strshort, L[li].kappa_rec->str,
             PARAMSTRLENSHORT - 1);
    L[li].kappa_rec->v->do_xyplot = 1;
    L[li].kappa_rec->v->do_logplot = 0;
    L[li].kappa_rec->v->do_trend = 1;
    L[li].kappa_rec->v->plotrange = L[li].kappa_rec->pr;
    L[li].kappa_rec->v->plotrescale = 1.0;
    L[li].kappa_rec->v->do_autoc = 0;

    init_value_record (L[li].kappa_rec->v, 0);
  }

  // do A_rec 
  if (L[li].model == STEPWISE || L[li].model == JOINT_IS_SW)
  {
    L[li].A_rec = malloc (L[li].nlinked * sizeof (struct chainstate_record_updates_and_values));
    for (ai = 0; ai < L[li].nlinked; ai++)
    {
      if (L[li].umodel[ai] == STEPWISE)
        sprintf (L[li].A_rec[ai].str, "%dSW%d", li, ai);
      else
        sprintf (L[li].A_rec[ai].str, "NOTE STEPWISE");
      L[li].A_rec[ai].num_uptypes = num_A_update_types;
      L[li].A_rec[ai].upnames = malloc (num_A_update_types * sizeof (strnl));
      for (i = 0; i < num_A_update_types; i++)
      {
        sprintf (L[li].A_rec[ai].upnames[i], "STR update");
      }
      L[li].A_rec[ai].upinf =
        calloc ((size_t) num_A_update_types, sizeof (struct update_rate_calc));
      L[li].A_rec[ai].num_vals = 0;
    }
  }
}                               // init_mutation_scalar_rec

void add_priorinfo_to_output(char priorfilename[],int *fpstri, char fpstr[])
{
  int i;
  int gp; // VS gp
  SP "\nParameter Priors\n");
  SP "-----------------\n");
  if (!calcoptions[USEPRIORFILE])
  {
    SP "  Population size parameters maximum value : %.4lf \n", thetaprior);
    if (modeloptions[EXPOMIGRATIONPRIOR])
      SP "  Migration rate parameters exponential distribution mean : %.4lf \n",
        mprior);
    else
      SP "  Migration rate parameters maximum value: %.4lf \n", mprior);
    if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 0)
    {
      SP "  Splitting time : %.4lf\n", tprior);
    }
    SP "\n");
  }
  else
  {
    SP"  Prior distribution terms given in file: %s\n\n",priorfilename);
    SP"  Splitting time maximum values \n");
    SP"\tPeriod\tPrior maximum value\n");
    for (i=0;i<numsplittimes;i++)
      SP"\t%s\t%.3lf\n",T[i].str,T[i].pr.max);
    SP"  Population size parameters maximum values \n");
    SP"\tPopulation\tPrior maximum value\n");
    for (i=0;i<numtreepops;i++)
      SP"\t%s\t%.3lf\n",itheta[i].str,itheta[i].pr[gp].max); // VS gp
    if (modeloptions[EXPOMIGRATIONPRIOR]==0)
    {
      SP"  Migration parameters maximum values \n");
      SP"\tMigration rate\tMaximum value\n");
    }
    else
    {
      SP"  Migration parameter prior means (exponential priors)\n");
      SP"\tMigration rate\tPrior mean value\n");
    }
    for (i=0;i<nummigrateparams;i++)
      SP"\t%s\t%.3lf\n",imig[i].str,imig[i].pr[gp].max); // VS
  }
}

// VS
// ADD_PRIORINFO_TO_OUTPUT_VS
// prints to output information aboput the prior distributions
// modified to print the info about the groups of loci
// INPUT
//	char priorfilename[] : string with the filename of output
//	int *fpstri : array of integers with?
//	char fpstr[] : array of characters with?
void add_priorinfo_to_output_vs(char priorfilename[],int *fpstri, char fpstr[])
{
  int i;
  int gp; // VS gp
  SP "\nParameter Priors\n");
  SP "-----------------\n");
  if (!calcoptions[USEPRIORFILE])
  {
    SP "  Population size parameters maximum value : %.4lf \n", thetaprior);
    if (modeloptions[EXPOMIGRATIONPRIOR])
      SP "  Migration rate parameters exponential distribution mean : %.4lf \n",
        mprior);
    else
      SP "  Migration rate parameters maximum value: %.4lf \n", mprior);
    if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 0)
    {
      SP "  Splitting time : %.4lf\n", tprior);
    }
    SP "\n");
  }
  else
  {
    SP"  Prior distribution terms given in file: %s\n\n",priorfilename);
    SP"  Splitting time maximum values \n");
    SP"\tPeriod\tPrior maximum value\n");
    for (i=0;i<numsplittimes;i++)
      SP"\t%s\t%.3lf\n",T[i].str,T[i].pr.max);
    SP"  Population size parameters maximum values \n");

	// VS
	// here need to loop through the number of groups
	for(gp=0; gp<nbgroupsloci_theta; gp++) {
		SP"group of loci theta %i", gp);
		SP"\tPopulation\tPrior maximum value\n");
		for (i=0;i<numtreepops;i++)
		  SP"\t%s\t%.3lf\n",itheta[i].str,itheta[i].pr[gp].max); // VS gp
	}
	for(gp=0; gp<nbgroupsloci_mig; gp++) {
		SP"group of loci mig %i", gp);
		if (modeloptions[EXPOMIGRATIONPRIOR]==0)
		{
		  SP"  Migration parameters maximum values \n");
		  SP"\tMigration rate\tMaximum value\n");
		}
		else
		{
		  SP"  Migration parameter prior means (exponential priors)\n");
		  SP"\tMigration rate\tPrior mean value\n");
		}
		for (i=0;i<nummigrateparams;i++)
		  SP"\t%s\t%.3lf\n",imig[i].str,imig[i].pr[gp].max); // VS
	}
  }
}




void
start_setup_L (char infilename[], int *fpstri, char fpstr[])
{
  int li, i;

  /* get the number of loci and the number of populations from the top of the datafile */
  read_datafile_top_lines (infilename, fpstri, fpstr, startpoptreestring);
  /* setup a temporary struture to record how much variation there is,  used for picking starting values of mutation scalars */
  numsitesIS = malloc (nloci * sizeof (int *));
  uvals = malloc (nloci * sizeof (double *));   // rows in this matrix are malloced in setup_L
  for (i = 0; i < nloci; i++)
    numsitesIS[i] = calloc ((size_t) npops + 1, sizeof (int));
  readdata (infilename, fpstri, fpstr, numsitesIS);
  for (li = 0; li < nloci; li++)
  {
    init_mutation_scalar_rec (li);
    init_g_rec (li);
    uvals[li] = malloc (L[li].nlinked * sizeof (double));
    if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
    {
      init_a_rec (li);
    }
  }
  if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
  {
    init_a_rec_multichain ();
  }
  uterm = set_uvals ();
  if (modeloptions[ADDGHOSTPOP])
  {
    assert (assignmentoptions[POPULATIONASSIGNMENT] == 0);                          
    if (npops + 1 > MAXPOPS)
    {
      IM_err (IMERR_INPUTFILEINVALID, 
              "ghost population makes number of population (%d) greater than MAXPOPS [%d]", 
              npops, MAXPOPS);
    }
     for (li = 0; li < nloci; li++)
     {
      L[li].samppop[npops] = 0;
     }
    npops++;
    numtreepops += 2;
    lastperiodnumber++;
    numsplittimes++;
  }
  gi_largestngenes = 0;
  gi_largestnumsites = 0;
  for (li = 0; li < nloci; li++)
  {
    if (gi_largestngenes < L[li].numgenes)
    {
      gi_largestngenes = L[li].numgenes;
    }
    if (gi_largestnumsites < L[li].numsites)
    {
      gi_largestnumsites = L[li].numsites;
    }
  }
 return;
}                               //start_setup_L


// START_SETUP_C
// allocates and initializes instances of Chain C structure
// 1. allocate array of C
// 2.1. init_genealogy_weights for C[ci]->allgweight (sum of all gweights)
// 2.2. init_genealogy_weights for C[ci]->G[li].gweight (gweight for each locus)
// 3. setup_poptree
// 4. set_tvalues
void
start_setup_C (void)
{
  int ci, li;
  int gp; // VS gp

  // Allocate memory for chain array
  C = malloc (numchains * sizeof (struct chain *));     //points to an array of chains 
  for (ci = 0; ci < numchains; ci++)
  {
    C[ci] = malloc (sizeof (struct chain));
    sprintf (C[ci]->name, "%d", ci);
    if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
    {
      C[ci]->nasn = malloc (npops * sizeof (int));
    }
  }
  for (ci = 0; ci < numchains; ci++)
  {
    init_genealogy_weights (&C[ci]->allgweight);

    /* CR: 110516.2
     * change malloc to calloc so that memory would be initialized to
     * a known state before use.
     */
    C[ci]->G = calloc ((size_t) nloci, sizeof (struct genealogy));

    for (li = 0; li < nloci; li++)
    {
      init_genealogy_weights (&(C[ci]->G[li].gweight));
    }

	// VS
	// allocate memory for the gweights for groups of loci
	// the groupgweight are pointers to arrays of size nbgroupsloci_theta and nbgroupsloci_mig
	C[ci]->groupgweight_theta = (struct genealogy_weights*) malloc(nbgroupsloci_theta * sizeof(struct genealogy_weights)); 
	C[ci]->groupgweight_mig = (struct genealogy_weights*) malloc(nbgroupsloci_mig * sizeof(struct genealogy_weights)); 
	
	// VS
	// initialize the gweights for the groups of loci for theta and mig
	for(gp=0; gp<nbgroupsloci_theta; gp++) {
		init_genealogy_weights(&(C[ci]->groupgweight_theta[gp]));
	}
	for(gp=0; gp<nbgroupsloci_mig; gp++) {
		init_genealogy_weights(&(C[ci]->groupgweight_mig[gp]));
	}

    setup_poptree (ci, startpoptreestring);
    set_tvalues (ci);

// VS
// CHECK
/*fprintf(fgweights_theta, "\n\n\nCHECK init_genealogy weights for theta");
for(gp=0; gp<nbgroupsloci_theta;gp++) {
	print_gweight_vs(&(C[ci]->groupgweight_theta[gp]), 0, fgweights_theta);
}
fprintf(fgweights_mig, "\n\n\nCHECK init_genealogy weights for mig");
for(gp=0; gp<nbgroupsloci_mig;gp++) {
	print_gweight_vs(&(C[ci]->groupgweight_mig[gp]), 1, fgweights_mig);
}*/


  }
}                               // start_setup_C

// FINISH_SETUP_C
// ends setting up the chain structure
// 1. init_treeweight
// 2. init_gtreecommon
// 3.1. init_probcalc for allpcalc
// 3.2. 
// 4. makeIS - does this proposes a first tree consistent with the data?
// 5. treeweight
// 6. initialize_integrate_tree_prob
void
finish_setup_C ()
{
  int ci, li, i, ai;
  int gp; // VS gp	
  // VS - sum_probg used to compute the P(G) as a sum of P(G_g) for genealogies belonging to each group g
  double sum_probg; 
  int nosimmigration;
  init_treeweight ();
  init_gtreecommon ();
  if (modeloptions[NOMIGRATION]==0)
  {
    set_nomigrationchecklist ();
  }
  for (ci = 0; ci < numchains; ci++)
  {
    for (li = 0; li < nloci; li++)
    {
      C[ci]->G[li].gtree = malloc (L[li].numlines * sizeof (struct edge));
      C[ci]->G[li].uvals = malloc (L[li].nlinked * sizeof (double));
      C[ci]->G[li].pdg_a = malloc (L[li].nlinked * sizeof (double));
    }
  }
  if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
  {
    /* Individual names are set here. */
    imaAsnSet ();
  }

  for (ci = 0; ci < numchains; ci++)
  {
	// probcalc for allpcalc is initialized here
    init_probcalc (&(C[ci]->allpcalc));

	// VS
	// allocate memory for grouppcalc
	C[ci]->grouppcalc_theta = (struct probcalc *) malloc(nbgroupsloci_theta*sizeof(struct probcalc)); 
	C[ci]->grouppcalc_mig = (struct probcalc *) malloc(nbgroupsloci_mig*sizeof(struct probcalc)); 

	// VS
	// the probcalc for each group of loci is initialized here
	for(gp=0; gp<nbgroupsloci_theta; gp++) {
		init_probcalc (&(C[ci]->grouppcalc_theta[gp]));
	}
	for(gp=0; gp<nbgroupsloci_mig; gp++) {
		init_probcalc (&(C[ci]->grouppcalc_mig[gp]));
	}


    for (li = 0; li < nloci; li++)
    {
      for (i = 0; i < L[li].numlines; i++)
      {
        C[ci]->G[li].gtree[i].mig = malloc (MIGINC * sizeof (struct migstruct));
        C[ci]->G[li].gtree[i].cmm = MIGINC;
        C[ci]->G[li].gtree[i].mig[0].mt = -1;
        C[ci]->G[li].gtree[i].up[0] = -1;
        C[ci]->G[li].gtree[i].up[1] = -1;
        C[ci]->G[li].gtree[i].down = -1;
        C[ci]->G[li].gtree[i].time = 0;
        C[ci]->G[li].gtree[i].mut = -1;
        C[ci]->G[li].gtree[i].pop = -1;
        C[ci]->G[li].gtree[i].ei = i;
        C[ci]->G[li].gtree[i].exist = 'T';

        if (L[li].model == STEPWISE || L[li].model == JOINT_IS_SW)
        {
          /* CR: 110516.1
           * change malloc to calloc so that memory would be initialized to
           * a known state before use.
           */
          C[ci]->G[li].gtree[i].A = 
            calloc ((size_t) L[li].nlinked, sizeof (int));
          C[ci]->G[li].gtree[i].dlikeA =
            calloc ((size_t) L[li].nlinked, sizeof (double));
        }

        if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
        {
          if (L[li].model == INFINITESITES)
          {
            C[ci]->G[li].gtree[i].seq = malloc (L[li].numsites * sizeof (int));
            if (i < L[li].numgenes)
              {
                memcpy (C[ci]->G[li].gtree[i].seq, L[li].seq[i], L[li].numsites * sizeof (int));
              }
          }
        }
      }
      if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
      {
        C[ci]->G[li].mut = malloc (L[li].numsites * (sizeof (int)));
      }
    }
    /* we need member ei to be ready for initial assignment */
    imaAsnInitAssign (ci);

    if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 1)
    {
      for (i = 0; i < npops; i++)
      {
        assert (C[ci]->poptree[i].b == 0);
        assert (C[ci]->poptree[i].e == 1);
        C[ci]->poptree[i].time = TIMEMAX;
      }
    }
    else
    {
      for (i = 0; i < 2 * npops - 1; i++)
      {
        if (C[ci]->poptree[i].e == -1)
          C[ci]->poptree[i].time = TIMEMAX;
        else
          C[ci]->poptree[i].time = C[ci]->tvals[C[ci]->poptree[i].e - 1];
      }
    }

    if (ci == 0)                //&& nurates > 1)
    {
      setuinfo (uterm);
    }

    for (li = 0; li < nloci; li++)
    {
      if (ci > 0 && nurates > 1)
      {
        for (ai = 0; ai < L[li].nlinked; ai++)
          C[ci]->G[li].uvals[ai] = C[0]->G[li].uvals[ai];
      }
      if (nurates == 1)
      {
        C[ci]->G[li].uvals[0] = 1.0;
      }

      nosimmigration = !(modeloptions[PARAMETERSBYPERIOD] || nummigrateparams == maxpossiblemigrateparams);
      
      switch (L[li].model)
      {
      case INFINITESITES:
        makeIS (ci, li,nosimmigration);
        treeweight (ci, li);
        C[ci]->G[li].pdg = C[ci]->G[li].pdg_a[0] = likelihoodIS (ci, li, C[ci]->G[li].uvals[0]);
        break;
      case HKY:
        if (assignmentoptions[JCMODEL] == 1)
        {
          makeHKY (ci, li,nosimmigration);
          treeweight (ci, li);
          C[ci]->G[li].pdg_a[0] = likelihoodJC (ci, li, C[ci]->G[li].uvals[0]);
          C[ci]->G[li].pdg = C[ci]->G[li].pdg_a[0];
        }
        else
        {
          for (i = 0; i < 4; i++)
          {
            C[ci]->G[li].pi[i] = pi[li][i];
          }
          C[ci]->G[li].kappaval = 2.0;    // starting kappa value
          makeHKY (ci, li,nosimmigration);
          treeweight (ci, li);
          C[ci]->G[li].pdg_a[0] = likelihoodHKY (ci, li, C[ci]->G[li].uvals[0], 
                                                 C[ci]->G[li].kappaval, 
                                                 -1, -1, -1, -1);
          C[ci]->G[li].pdg = C[ci]->G[li].pdg_a[0];
          copyfraclike (ci, li);
          storescalefactors (ci, li);
        }
        break;
      case STEPWISE:
        somestepwise = 1;

        for (ai = 0; ai < L[li].nlinked; ai++)
          for (i = 0; i < L[li].numgenes; i++)
            C[ci]->G[li].gtree[i].A[ai] = L[li].A[ai][i];
        makeSW (ci, li,nosimmigration);
        treeweight (ci, li);
        C[ci]->G[li].pdg = 0;
        for (ai = 0; ai < L[li].nlinked; ai++)
        {
          C[ci]->G[li].pdg_a[ai] = likelihoodSW (ci, li, ai, C[ci]->G[li].uvals[ai], 1.0);
          C[ci]->G[li].pdg += C[ci]->G[li].pdg_a[ai];
        }

        break;
      case JOINT_IS_SW:
        for (ai = 1; ai < L[li].nlinked; ai++)
          for (i = 0; i < L[li].numgenes; i++)
            C[ci]->G[li].gtree[i].A[ai] = L[li].A[ai][i];
        somestepwise = 1;
        makeJOINT_IS_SW (ci, li,nosimmigration);
        treeweight (ci, li);
        C[ci]->G[li].pdg = C[ci]->G[li].pdg_a[0] =
          likelihoodIS (ci, li, C[ci]->G[li].uvals[0]);
        for (ai = 1; ai < L[li].nlinked; ai++)
        {
          C[ci]->G[li].pdg_a[ai] = likelihoodSW (ci, li, ai, C[ci]->G[li].uvals[ai], 1.0);
          C[ci]->G[li].pdg += C[ci]->G[li].pdg_a[ai];
        }
        break;
      }
     
	  // Save the P(D|G) in allpcalc.pdg
      C[ci]->allpcalc.pdg += C[ci]->G[li].pdg;
	  // Sum the gweights of all loci in allgweight
      sum_treeinfo (&(C[ci]->allgweight), &(C[ci]->G[li].gweight));
	} // end of for from locus 0 to nloci

   	// VS (maybe these two for loops can go to a function called sum_gweights_group)
	// after loop through loci, all the gweight matrix for each locus are initialized
	// at this point we can obtain the sum of gweights for each group of locus
	for(gp=0; gp<nbgroupsloci_theta; gp++) {
		for(li=0; li<nloci; li++) {
			// if the locus belongs to the group
			if(grouploci_theta[ci][li]==gp) {
				 // add the gweights of loci that belong to this group to groupgweight_theta[gp]
				 sum_treeinfo_theta_vs(&(C[ci]->groupgweight_theta[gp]), &(C[ci]->G[li].gweight));
				 // CHECK
				 check_gweight_vs(&(C[ci]->groupgweight_theta[gp]), 1);
			}
		}
	}
	// VS
	// obtain the sum of gweights for each group of locus for mig parameters
	for(gp=0; gp<nbgroupsloci_mig; gp++) {
		for(li=0; li<nloci; li++) {
			// if the locus belongs to the group
			if(grouploci_mig[ci][li]==gp) {
				 // add the gweights of loci that belong to this group to groupgweight_theta[gp]
				 sum_treeinfo_mig_vs(&(C[ci]->groupgweight_mig[gp]), &(C[ci]->G[li].gweight));
				 // CHECK
				 check_gweight_vs(&(C[ci]->groupgweight_mig[gp]), 0);
			}
		}
	}

	// VS - CHANGE THE INITIALIZATION OF P(G) later!!!
	// This is the initial integrate_tree probability
	// This will initialize the allpcalc values for P(G) based on the allgweights
	// when there are groups of loci this is not correct
	// we need to compute the the tree_probability as the product of the integrals
	// obtained for each group of loci
	// this can be done later

	// VS CHECK
	//fprintf(fpcalc_all, "\nAllpcalcFinishSetupBeforeInitializeIntegrateVS ");
	//print_pcalc_groups_vs_file(&(C[ci]->allpcalc), fpcalc_all);

    initialize_integrate_tree_prob_vs (ci, &C[ci]->allgweight, &C[ci]->allpcalc, 0, 2);

	// VS CHECK
	// fprintf(fpcalc_all, "\nAllpcalcFinishSetupAFTERInitializeIntegrateVS ");
	// print_pcalc_groups_vs_file(&(C[ci]->allpcalc), fpcalc_all);


	// VS (maybe these two for loops can go to a function called sum_pcalc_group)
	// Need to initialize the integrate_tree_prob for each group of loci
	// A special functions may be needed here, because the pcalc of each group
	// will be the sum of integrate theta for the theta groups
	// and the sum of integrate mig for the mig groups
	// then the prior probability of the genealogy for all loci
	// will be the sum of the theta.probg and the mig.probg across groups
	// the total prior probability will be saved in allpcalc
	// after loop through loci, all the gweight matrix for each locus are initialized
	// at this point we can obtain the sum of gweights for each group of locus
	for(gp=0; gp<nbgroupsloci_theta; gp++) {
		for(li=0; li<nloci; li++) {
			// if the locus belongs to the group
			if(grouploci_theta[ci][li]==gp) {
				// compute the integrate theta for each theta param on each group of loci
				initialize_integrate_tree_prob_vs (ci, &(C[ci]->groupgweight_theta[gp]), &(C[ci]->grouppcalc_theta[gp]), gp, 1);
				// VS CHECK
				check_pcalc_groups_vs(&(C[ci]->grouppcalc_theta[gp]), 1);
			}
		}
	}
	// VS
	// obtain the sum of gweights for each group of locus for mig parameters
	for(gp=0; gp<nbgroupsloci_mig; gp++) {
		for(li=0; li<nloci; li++) {
			// if the locus belongs to the group
			if(grouploci_mig[ci][li]==gp) {
				// compute the integrate theta for each theta param on each group of loci
				initialize_integrate_tree_prob_vs (ci, &(C[ci]->groupgweight_mig[gp]), &(C[ci]->grouppcalc_mig[gp]), gp, 0);
				// VS CHECK
				check_pcalc_groups_vs(&(C[ci]->grouppcalc_mig[gp]), 0);
			}
		}
	}

	// VS
	// Based on the pcalcs for each group, need to compute the prior probability as the sum of pcalcs
	// save this information in C[ci]->allpcalc
	sum_probg=0;
	for(gp=0; gp<nbgroupsloci_theta; gp++) {
		sum_probg += C[ci]->grouppcalc_theta[gp].probg;
	}
	for(gp=0; gp<nbgroupsloci_mig; gp++) {
		sum_probg += C[ci]->grouppcalc_mig[gp].probg;
	}
	// set the sum_probg as the allpcalc.probg
	C[ci]->allpcalc.probg = sum_probg;
	
		// VS CHECK GWEIGHTS ---------------------------------------
		/*for(li=0; li<nloci;li++) {
			fprintf(fgweights_loci, "\nGweightsFinishSetupC %i ", li);
			print_gweight_vs_file(&(C[ci]->G[li].gweight), 2, fgweights_loci);
		}		
		for(gp=0; gp<nbgroupsloci_theta; gp++) {
			fprintf(fgweights_theta, "\nGroupThetaFinishSetupC Gp=%i ", gp);
			print_gweight_vs_file(&(C[ci]->groupgweight_theta[gp]), 2, fgweights_theta);
		}
		for(gp=0; gp<nbgroupsloci_mig; gp++) {
			fprintf(fgweights_mig, "\nGroupMigFinishSetupC Gp=%i ", gp);
			print_gweight_vs_file(&(C[ci]->groupgweight_mig[gp]), 2, fgweights_mig);
		}

		fprintf(fgweights_all, "\nGweightAllSinishSetupC locus=%i ", li);
		print_gweight_vs_file(&(C[ci]->allgweight), 2, fgweights_all); 
	
		
		// VS CHECK PCALC ---------------------------------------
		
		for(gp=0; gp<nbgroupsloci_theta; gp++) {
			fprintf(fpcalc_theta, "\nGroupThetaFinishSetupC Gp=%i ", gp);
			print_pcalc_groups_vs_file(&(C[ci]->grouppcalc_theta[gp]), fpcalc_theta);
		}
		
		for(gp=0; gp<nbgroupsloci_mig; gp++) {
			fprintf(fpcalc_mig, "\nGroupMigFinishSetupC Gp=%i ", gp);
			print_pcalc_groups_vs_file(&(C[ci]->grouppcalc_mig[gp]), fpcalc_mig);
		}
		fprintf(fpcalc_all, "\nAllpcalcFinishSetup ");
		print_pcalc_groups_vs_file(&(C[ci]->allpcalc), fpcalc_all);
		*/
		
	
  } // end of for loop for each chain
  return;
}                               // finish_setup_C


/* set random times for poptree */
  /* first pick trandom times over interval from 0 to tprior 
     use broken stick model - simple dirichlet distribution
   */
/* to fix values:
in imamp.h  turn on these undefs
#undef DO_RY1UPDATE
#undef DO_RY1UPDATE
and at the end of set_tvalues()
use code like 
 C[ci]->tvals[0]=0.1;C[ci]->tvals[1]=0.2;C[ci]->tvals[2]=0.4;C[ci]->tvals[3]=0.8; */

void
set_tvalues (int ci)
{
  int i;
  double sum;
  double times[MAXPOPS];

  if (npops == 1)
  {
    C[ci]->tvals = malloc (sizeof (double));
    C[ci]->tvals[0] = TIMEMAX;
  }
  else if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 1)
  {
    C[ci]->tvals = malloc (sizeof (double));
    C[ci]->tvals[0] = TIMEMAX;
  }
  else
  {
    C[ci]->tvals = malloc ((lastperiodnumber + 1) * sizeof (double));
    for (sum = 0, i = 0; i < npops; i++)
    {
      times[i] = expo (1.0);
      sum += times[i];
    }

    for (i = 0; i < lastperiodnumber; i++)
    {
      times[i] *= (T[i].pr.max - T[i].pr.min) / sum;
      times[i] += T[i].pr.min;
      if (i > 0)
        times[i] += times[i - 1];
      C[ci]->tvals[i] = times[i];
      assert (C[ci]->tvals[i] < T[i].pr.max && C[ci]->tvals[i] > T[i].pr.min);
    }
    C[ci]->tvals[i] = TIMEMAX;
  }
#ifndef DO_NWUPDATE
#ifndef DO_RY1UPDATE
 if (lastperiodnumber ==2  && T[0].pr.max == 5.0)
 {
   C[ci]->tvals[0]=1.0;C[ci]->tvals[1]=4.0;
 }
 else
 {
  for (i=0;i<lastperiodnumber;i++)
    C[ci]->tvals[i] = (i+1.0)/(lastperiodnumber+1.0) * T[i].pr.max;
 }
 C[ci]->tvals[0] = 5.00; // VS fix tsplit
 //C[ci]->tvals[1] = 0.02;
#endif
#endif  
  return;
}                               /* set_tvalues */

/* JH 2/18/10 inserted this */
#define MAXTRANGE_IF_ISLANDMODEL  100.0  // have to pick something 
void
setup_T ()
{
  int i, li;
  /* JH 2/18/10  moved this here */ 
  if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 1  || npops == 1)
  {
    T = NULL;
    tprior = 0.0;
    for (li=0;li<nloci;li++)
    {
      L[li].g_rec->v->plotrange.max = MAXTRANGE_IF_ISLANDMODEL;
      init_value_record (L[li].g_rec->v, 0);
    }
  }
  else
  {
    T = calloc ((size_t) (npops - 1), sizeof (struct chainstate_record_updates_and_values));
    for (i = 0; i < lastperiodnumber; i++)
    {
      if (!calcoptions[USEPRIORFILE])
        tperiodpriors[i] = tprior;
      T[i].pr.max = tperiodpriors[i];
	  T[i].pr.min = 0;
      //T[i].pr.min = 10; // VS 9/12/2011 put a minimum for the time of split for the analysis of the rabbit data
    }
    for (i = 0; i < lastperiodnumber; i++)
    {
      sprintf (T[i].str, "t%d", i);
      T[i].num_uptypes = IM_UPDATE_TIME_NUMBER;
      T[i].upnames = malloc (T[i].num_uptypes * sizeof (strnl));
#ifdef DO_NWUPDATE
      if (assignmentoptions[POPULATIONASSIGNMENTBF] == 0)
      {
        sprintf (T[i].upnames[IM_UPDATE_TIME_NW], "NielsenWakeley");
      }
      else
      {
        sprintf (T[i].upnames[IM_UPDATE_TIME_NW], "RannalaYang");
      }
#endif
#ifdef DO_RY1UPDATE
      sprintf (T[i].upnames[IM_UPDATE_TIME_RY1], "RannalaYang");
#endif
      T[i].upinf = calloc ((size_t) T[i].num_uptypes, sizeof (struct update_rate_calc));
      T[i].num_vals = 1;
      T[i].v = malloc (T[i].num_vals * sizeof (struct value_record));
      strcpy (T[i].v->str, T[i].str);
      strncpy (T[i].v->strshort, T[i].v->str, PARAMSTRLENSHORT - 1);
      T[i].v->do_xyplot = 1;
      T[i].v->do_trend = 1;
      T[i].v->do_logplot = 0;
      T[i].v->do_autoc = 1;
      T[i].v->plotrange.max = tperiodpriors[i];
      T[i].v->plotrange.min = 0;
      T[i].v->plotrescale = 1.0;
      init_value_record (T[i].v, 0);
      /* 2/18/10 JH  inserted this,  needed for plotting TMRCAs */
      if (i==0)
        for (li=0;li<nloci;li++)
        {
          L[li].g_rec->v->plotrange.max = tperiodpriors[lastperiodnumber-1] * TIMEPRIORMULTIPLIER;
          init_value_record (L[li].g_rec->v, 0);
        }
    }
  }

}                               //setup_T

// VS 5/18/2012 function to initialize the assignment vector variables
// Allocates memory and initializes the variables related with recording and 
// getting statistics about the acceptance rates of a the assignment vectors of loci (not to be confounded with the assignment of individuals)
void setup_assignLoci ()
{
  int i, li;
  // 1 if we are updating both the theta and mig assignment jointly, such that a_m=a_theta
  // 2 if we are independently updating theta and mig assignment
  int tupdassignloci = 2; 
  int numupdates = 1; // there is only one update type for now (this could become and enumerate variable in future)

  assignloci = calloc ((size_t) tupdassignloci, sizeof (struct chainstate_record_updates_and_values));
  
  for (i = 0; i < tupdassignloci; i++)
  {
	  if(i==0) {sprintf (assignloci[i].str, "a_mig");}
	  if(i==1) {sprintf (assignloci[i].str, "a_theta");}
      assignloci[i].num_uptypes = numupdates; // currently numupdates is 1 (there is only one update type)
      assignloci[i].upnames = malloc (assignloci[i].num_uptypes * sizeof (strnl));
      
      sprintf (assignloci[i].upnames[0], "Random Assign"); // there is only one type of assignment, which is to randomly change the label of one locus at a time

	  // upinf has update information about the assignment vectors
      assignloci[i].upinf = calloc ((size_t) assignloci[i].num_uptypes, sizeof (struct update_rate_calc));
      
	  assignloci[i].num_vals = 1; 
	  // NOTE that the number of values to record is the number of loci, i.e. the length of each assignment vector
	  // however, the structure value record is not suited to deal with assignment vector
	  // it is better to create a different structure to deal with the assignment vector variables
      assignloci[i].v = malloc (assignloci[i].num_vals * sizeof (struct value_record));
      strcpy (assignloci[i].v->str, assignloci[i].str);
      strncpy (assignloci[i].v->strshort, assignloci[i].v->str, PARAMSTRLENSHORT - 1);
      assignloci[i].v->do_xyplot = 0;
      assignloci[i].v->do_trend = 0;
      assignloci[i].v->do_logplot = 0;
      assignloci[i].v->do_autoc = 0;
      assignloci[i].v->plotrange.max = nbgroupsloci_mig;
      assignloci[i].v->plotrange.min = 0;
      assignloci[i].v->plotrescale = 1.0;
      init_value_record (assignloci[i].v, 0);
  }

}                               //setup_assignLoci



void
finish_setup_L (void)
{
  int li;
  for (li = 0; li < nloci; li++)
  {
    XFREE (numsitesIS[li]);
    XFREE (uvals[li]);
  }
  XFREE (numsitesIS);
  XFREE (uvals);
}                               // finish_setup_L

void reportparamcounts(int *fpstri, char fpstr[])
{
  SP"\nParameter Counts\n");
  SP"----------------\n");
  SP"   Population sizes : %d\n",numpopsizeparams);
  SP"   Migration rates  : %d\n",nummigrateparams);
  SP"   Parameters in the MCMC simulation\n");
  SP"      Splitting times : %d\n",numsplittimes);
  SP"      Mutation scalars: %d\n",nurates); 
  SP"      HKY Kappa (ti/tv) ratios: %d\n",nkappas); 
} // reportparamcounts(void)

/**********  GLOBAL FUNCTIONS  *******/

void
setup (char infilename[], int *fpstri, char fpstr[], char priorfilename[])
{
  start_setup_L (infilename, fpstri, fpstr);

  if (calcoptions[USEPRIORFILE])
  {
    popsizeprior_fromfile = malloc (numtreepops * sizeof(double));
    mprior_fromfile = alt2d_alloc2Ddouble(numtreepops,numtreepops);
    readpriorfile(priorfilename,popsizeprior_fromfile,mprior_fromfile);
	printf("nbgroupsloci...%d \n\n\n", nbgroupsloci_theta);
  }
  
  setup_T ();

  // VS 5/18/2012 - added the setup of assignment vector update
  setup_assignLoci ();

  start_setup_C (); // VS initilizes the chain and locus structure
  // VS 
  // Note that at this point the prior file was already read, and the number of groups is known
  // this means that the instances of Chain pointing to the genealogy weights of the groups can
  // be allocated inside start_setup_C()
  setup_iparams (); // VS initializes the param structures

  if (modeloptions[ANCESTRALPOPSIZESSETTOLARGESTSAMPLED])
  {
    addpopsizeidtext(fpstri,fpstr);
  }
  // VS
  // finish setup Chain - VS variables allocated and initialized here
  finish_setup_C ();
  // reportparamcounts prints to the screen the number of parameters and related info
  reportparamcounts(fpstri, fpstr);
  finish_setup_L ();            // somethings in L need info from T , free up numsitesIS and uvals
  // finish_setup_L free memory of numsitesIS[li] and uvals[li] these were needed to start_setup_L() and setup_T() but not anymore  
  
  if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 0
      && npops > 1)
  {
    init_t_NW ();
    init_t_RY ();
  }
  init_updategenealogy ();
  init_lpgpd_v ();
  if (outputoptions[MIGRATEHIST])
    init_migration_counts_times ();
 // VS add_priorinfo_to_output(priorfilename,fpstri, fpstr);
  add_priorinfo_to_output_vs(priorfilename,fpstri, fpstr);
  // CHECK
  printf("CHECK fpstr\n%s", fpstr);
  
  init_autoc_pointers ();
  if (calcoptions[USEPRIORFILE])
  {
    XFREE(popsizeprior_fromfile);
    alt2d_free2D(mprior_fromfile);
  }

}                               // setup() 
