/* MigSelect 2013-2019 Arun Sethuraman, Vitor Sousa, Jody Hey */
/* This code was developed based on IMa2  2009-2010  Jody Hey, Rasmus Nielsen and Sang Chul Choi*/

#undef GLOBVARS
#include "imamp.h"

/* 10/3/09  work in progress
  adding code for reading a file that contains priors of parameters */ 


/* format of the prior file 

a title line 
then one or more lines beginning with # 

then popstring
then popstring with popsize priors
then popstring with t priors
then migration prior matrix 


*/

static int ancestralpopnums[2*MAXPOPS - 1]; 
static double *pp;
static char *treestringspot;
static int numpopnodes;
static struct popedge *temppoptree;
static void readprior_parenth0 (void);
static void readprior_parenth (int mode, int tempcurrent, int startparenth);
static void readprior_poptreeread (int mode, char *poptreestring);

#define MAXPRIORTEXTLINE 500

/* readprior_parenth0()
come in on first opening parenthesis
the ith opening parenthesis is associated with a number after its corresponding closing parenthesis

readprior_parenth0() will go through parentheses from left to right,
when it comes to a close parenthesis it records the ancestral node number for that pairing

ancestral node numbers range from npops to 2*npops-2 and proceed from lowest
to highest in order of what time they occured

the ancestral node number is recorded in the array ancestralpopnums[]

the position in the array is the count of which parenthesis pair has just closed, plus npops

in other words if it is the 0'th parenthesis pair (i.e. the first one that opened, meaning it is the
outermost pair),  then the ancestral node number is recorded in ancestralpopnums[npops]

If it is the ith pair that has closed, it is recorded in ancestralpopnums[npops + i]

This seemed to be the only way to get the correct labeling of internal nodes. 

Then when the function readprior_parenth() is called,  the correct times and populatino sizes can be associated 
with these ancestral populations. 
*/

void
readprior_parenth0 (void)
{
  int itemp;
  char *ne;
  int psetlist[MAXPOPS], nextlistspot, popennum;
  nextlistspot = 0;
  popennum = 0;
  ne = treestringspot;
  while (*ne != '\0')
  {
    if (*ne == '(')
    {
      psetlist[nextlistspot] = popennum;
      nextlistspot++;
      popennum++;
      ne++;
    }
    else
    {
      if (*ne == ')')
      {
        ne += 2;
        itemp = strtol (ne, &ne, 10);
        ancestralpopnums[npops + psetlist[nextlistspot - 1]] = itemp;
        nextlistspot--;
      }
      else
      {
        ne++;
      }
    }
  }
}                               /* parenth0 */


void
readprior_parenth (int mode, int tempcurrent, int startparenth)
/* recursive - reenters for every open parentheses */
{
  int itemp, current;
  static int nextnode;
  static int periodi;
  double val;
  char *ne;

  if (startparenth == 1)
  {
    nextnode = -1;
    periodi = 0;
  }
  current = ancestralpopnums[tempcurrent];

  treestringspot++;
  while (isspace (*(treestringspot + 1))) // why  + 1 ? 
    treestringspot++;

  /* next could be:
     - a number  (a simple upnode) read it and get ':' and float number after it
     - an open parenthesis (a complex upnode - call parenth)
     - a comma  skip it
     - a close parenthesis - close the node 
   */
  do
  {
    if (isdigit (*treestringspot))
    {
      /*itemp = atoi (treestringspot);
      treestringspot++; */
      itemp = strtol(treestringspot,&ne,10);							// read the id of population
      treestringspot=ne;
      treestringspot++; /* skip colon */
      val = strtod(treestringspot,&ne);		 // read time or size of population
      treestringspot=ne;
      if (mode ==1)
      {
        //assert(current == itemp);
        temppoptree[itemp].time = val;
        temppoptree[current].up[temppoptree[current].numup] = itemp;
        temppoptree[current].numup++;
        temppoptree[itemp].down = current;
      }
      if (mode==2) 
      {
        pp[itemp] = val;
      }
    }
    if (*treestringspot == ',')
      treestringspot++;
    if (*treestringspot == '(')
    {
      if (nextnode == -1)
        nextnode = npops + 1;
      else
        nextnode++;
      if (mode==1)
      {
        temppoptree[ancestralpopnums[nextnode]].down = current;
        temppoptree[current].up[temppoptree[current].numup] =  ancestralpopnums[nextnode];
        temppoptree[current].numup++;
      }
      readprior_parenth (mode, nextnode, 0);
    }
  } while (*treestringspot != ')');
  treestringspot++;             /* skip parentheses and colon*/
  if (*treestringspot == ':')
  {
    treestringspot++;
    itemp = strtol(treestringspot,&ne,10); // read the id of population
    treestringspot=ne;
    treestringspot++; /* skip colon */
    val = strtod(treestringspot,&ne);		 // read time or size of population
    treestringspot=ne;
    if (itemp < npops)
      IM_err (IMERR_POPTREESTRINGFAIL,
              " wrong number of ancestral populations indicated. string %s, step %d",
              treestringspot, step);
    periodi = itemp - npops;
    if (mode==1)
    {
      assert(current == itemp);
      temppoptree[current].time = val;
      temppoptree[current].b = periodi + 1;
      temppoptree[temppoptree[current].up[0]].e =
      temppoptree[temppoptree[current].up[1]].e = periodi + 1;
    }
    if (mode  == 2)
    {
      pp[itemp] = val;
    }

  }
  else
  { // is it possible to get here ? 
    temppoptree[current].b = periodi + 1;
    temppoptree[temppoptree[current].up[0]].e =
      temppoptree[temppoptree[current].up[1]].e = periodi + 1;
    periodi++;
  }
  if (temppoptree[current].down != -1)
  {
    numpopnodes++;
    current = temppoptree[current].down;
  }
  else
  {
    periodi++;
    temppoptree[current].e = -1;
  }
}                               /* readprior_parenth */

void
readprior_poptreeread (int mode, char *poptreestring)
{
  int i, j;

  /* read in the tree string until enough parentheses are found */
  /* pcount counts parentheses '(' is +1 ')' is -1 repeat until 0 */
  treestringspot = poptreestring;
  if (mode==0)
  {
    numpopnodes = 0;
    for (i = 0; i < npops; i++)
    {
      temppoptree[i].b = 0;
      temppoptree[i].numup = 0;
      temppoptree[i].up = malloc (2 * sizeof (int));
      for (j = 0; j < 2; j++)
        temppoptree[i].up[j] = -1;
      temppoptree[i].down = -1;
    }
    for (; i < numtreepops; i++)
    {
      temppoptree[i].numup = 0;
      temppoptree[i].up = malloc (2 * sizeof (int));
      for (j = 0; j < 2; j++)
        temppoptree[i].up[j] = -1;
      temppoptree[i].down = -1;
    }
    temppoptree[npops].down = -1;
    readprior_parenth0 ();
  }
  else
    readprior_parenth (mode, npops, 1);
  return;
}  /* end readprior_poptreeread */

// VS
// READPRIORFILE
// This function was changed by Vitor Sousa to read extra lines from the command line with the locus
// specific prior information
// In order to allow different groups of loci to have different prior parameters
// we need to add some extra lines to the priorfile
// extra line 1 - nb of groups for theta
// extra line 2 - relative maximum theta for each group (e.g. if 1, it means the priors are the same as defined above in the prior file)
// extra line 3 - index of group of each locus, e.g. 0 1 1 0 0 0, means there are 6 loci and that loci (2,3) belong to group 1 and that loci (1,4,5,6) belong to group 0
// extra line 4 - nb of groups for mig
// extra line 5 - relative maximum mig for each group (e.g. if 0.1, the prior values in matrix defined in priorfile are multiplied by 0.1)
// extra line 6 - index of group of each locus, e.g. 0 1 1 0 0 0, means there are 6 loci and that loci (2,3) belong to group 1 and that loci (1,4,5,6) belong to group 0
// Variables allocated and initialized:
//	nbgroupsloci_theta, nbgroupsloci_mig, 
//  rel_theta, rel_mig, 
//  grouploci_theta, grouploci_mig
void readpriorfile(char priorfilename[],double *popsizepriorvals, double **mpriorvals)
{
  FILE *priorfile; 
  char  *chpt;
  int i,j;
  char *priortextline;
  char *treetext;
  if ((priorfile = fopen (priorfilename, "r")) == NULL)
  {
    IM_err(IMERR_READFILEOPENFAIL,"Error opening file with prior values: %s", priorfilename);
  }
  pp = popsizepriorvals;
  priortextline = malloc(MAXPRIORTEXTLINE*sizeof(char));
  temppoptree = malloc (numtreepops * sizeof (struct popedge));
  // read the first 3 lines of prior file that are not commented with #
  // this lines will have the treephylogeny, thetapriors, tsplit priors
  // this information is saved in the temppoptree variable
  for (i=0;i<3;i++)
  {
    while (fgets(priortextline,MAXPRIORTEXTLINE,priorfile)!= NULL && priortextline[0]=='#') {};
    treetext = &(priortextline[0]);
    readprior_poptreeread(i,treetext); // Read the prior for the tsplit and the prior for the thetas in each pop
									   // These are stored in the variables temppoptree[].time and pp[]
  }

  // initialize the tsplit priors, saved in variable tperiodpriors
  for (i=npops;i<numtreepops;i++)
  {
    tperiodpriors[i-npops]=temppoptree[temppoptree[i].up[0]].time;
    if (temppoptree[temppoptree[i].up[0]].time != temppoptree[temppoptree[i].up[1]].time)
      IM_err(IMERR_PRIORFILEVALS,"two split times specified in tree string in prior file are not equal: %lf, %lf",temppoptree[temppoptree[i].up[0]].time,temppoptree[temppoptree[i].up[1]].time);
  }
  for (i=0;i<numsplittimes-1;i++)
    if (tperiodpriors[i] > tperiodpriors[i+1])
      IM_err(IMERR_PRIORFILEVALS,"earlier max split time greater than max for older split:%lf, %lf", tperiodpriors[i],tperiodpriors[i+1]);
  // read in migration rates
  // get the migration prior values reading the 2d matrix (npop*npop)
  max_m_from_priorfile = -1.0;
  for (i=0;i<numtreepops;i++)									
  {
      while (fgets(priortextline,MAXPRIORTEXTLINE,priorfile)!= NULL && priortextline[0]=='#') {};
	  chpt = &priortextline[0];
	  for (j=0;j<numtreepops;j++)
      {
	    mpriorvals[i][j] = strtod(chpt,&chpt);
        if (mpriorvals[i][j] > max_m_from_priorfile )
          max_m_from_priorfile  = mpriorvals[i][j];
        if (i>j && ( (mpriorvals[i][j]==0.0 && mpriorvals[j][i] > 0.0) ||(mpriorvals[i][j] > 0.0 && mpriorvals[j][i] == 0.0)))
          IM_err(IMERR_PRIORFILEVALS,"reciprocal migration rates both not zero or both not greater than zero: %d->%d %lf; %d->%d, %lf", i,j,mpriorvals[i][j],j,i,mpriorvals[j][i]);
       }
  }

  // VS 
  // read the extra lines of the input file
  // THETA
  // 1st - read number of groups for theta
  // skip the commentary lines
  while (fgets(priortextline,MAXPRIORTEXTLINE,priorfile)!= NULL && priortextline[0]=='#') {};
  chpt = &priortextline[0];
  // read nb of groups for theta
  nbgroupsloci_theta = atoi(chpt);
  // Allocate memory to save the nbgroupsloci_theta
  // This needs to be done after reading the number of groups of loci
  // Not sure if this is the best place to allocate memory
  reltheta = (double *) malloc( nbgroupsloci_theta * sizeof(double));
  grouploci_theta = (int **) malloc( numchains * sizeof(int *));
  for(i=0; i<numchains; i++) {
	grouploci_theta[i] = (int *) malloc( nloci * sizeof(int)); // Each chain has its own array (this is needed when updating the arrays of loci assignment)
  }


  // line 2 - read the relative theta max for each group of loci
  while (fgets(priortextline,MAXPRIORTEXTLINE,priorfile)!= NULL && priortextline[0]=='#') {};
  chpt = &priortextline[0];
  for (i=0; i<nbgroupsloci_theta; i++)									
  {
	  reltheta[i] = strtod(chpt, &chpt);
	  if (reltheta[i]<=0)
	  {
		IM_err(IMERR_PRIORFILEVALS,"The relative theta for the upper limit prior for locus %i has to be positive. Currently it is %g", i, reltheta[i]);
	  }
  }
  // line 3 - read the group index for each locus
  while (fgets(priortextline,MAXPRIORTEXTLINE,priorfile)!= NULL && priortextline[0]=='#') {};
  chpt = &priortextline[0];
  for (i=0; i<nloci;i++) {
	// for each group read the number of loci that belong to the group and read the index
	// Read the prior file into the first chain (index 0), then copy this to all the other chains
	grouploci_theta[0][i] = (int) strtod(chpt, &chpt);
	for(j=1; j<numchains; j++) {
		grouploci_theta[j][i] = grouploci_theta[0][i];
	}
  }
	printf("grouploci_theta has been allocated...\n");

  // MIG
  // 4th - read number of groups for mig
  // skip the commentary lines
  while (fgets(priortextline,MAXPRIORTEXTLINE,priorfile)!= NULL && priortextline[0]=='#') {};
  chpt = &priortextline[0];
  // read nb of groups for mig
  nbgroupsloci_mig = atoi(chpt);
 
  // Allocate memory to save the nbgroupsloci_theta
  // This needs to be done after reading the number of groups of loci
  // Not sure if this is the best place to allocate memory
  relmig = (double *) malloc( nbgroupsloci_mig * sizeof(double));
  grouploci_mig = (int **) malloc( numchains * sizeof(int *)); // Each chain has its own array (this is needed when updating the arrays of loci assignment)
  for(i=0; i<numchains; i++) {
	grouploci_mig[i] = (int *) malloc( nloci * sizeof(int));
  }
	printf("grouloci_mig has been alllocated...\n");

  // line 5 - read the relative mig max for each group of loci
  while (fgets(priortextline,MAXPRIORTEXTLINE,priorfile)!= NULL && priortextline[0]=='#') {};
  chpt = &priortextline[0];
  for (i=0; i<nbgroupsloci_mig; i++)									
  {
	  relmig[i] = strtod(chpt, &chpt);
	 
	  /*if (relmig[i]<=0)
	  {
		IM_err(IMERR_PRIORFILEVALS,"The relative mig for the upper limit prior for locus %i has to be positive. Currently it is %g", i, relmig[i]);
	  }*/
  }
  // line 6 - read the group index for each locus
  while (fgets(priortextline,MAXPRIORTEXTLINE,priorfile)!= NULL && priortextline[0]=='#') {};
  chpt = &priortextline[0];
  for (i=0; i<nloci;i++) {
	// for each group read the number of loci that belong to the group and read the index
    // Read the prior file into the first chain (index 0), then copy this to all the other chains
	grouploci_mig[0][i] = (int) strtod(chpt, &chpt);
	for(j=1; j<numchains; j++) {
		grouploci_mig[j][i] = grouploci_mig[0][i];
	}
  }


  // VS 5/18/2012 Initialize the acceptance counter for mig assignment



// VS
// CHECK
// Check the variables that were read
printf("\n\n\nCHECK PRIORFILE\nnbgroupsloci_theta is %i", nbgroupsloci_theta);
printf("\nreltheta ");
for(i=0; i<nbgroupsloci_theta; i++) {
	printf(" %g", reltheta[i]);
}
printf("\ngrouploci_THETA\n");
for (i=0; i<nloci;i++) {
	// for each group read the number of loci that belong to the group and read the index
	printf("%i ", grouploci_theta[0][i]); // only print the first chain
}
printf("\nnbgroupsloci_mig is %i", nbgroupsloci_mig);
printf("\nrelmig ");
for(i=0; i<nbgroupsloci_mig; i++) {
	printf(" %g", relmig[i]);
}
printf("\ngrouploci_MIG\n");
for (i=0; i<nloci;i++) {
	// for each group read the number of loci that belong to the group and read the index
	printf("%i ", grouploci_mig[0][i]); // only print the first chain
}


  for (i=0; i < numtreepops; i++)
    XFREE(temppoptree[i].up);
  XFREE(temppoptree);
  XFREE(priortextline);
  fclose(priorfile);
} // readpriorfile
   

