/* MigSelect 2013-2019 Arun Sethuraman, Vitor Sousa, Jody Hey */
/* This code was developed based on IMa2  2009-2010  Jody Hey, Rasmus Nielsen and Sang Chul Choi*/
#undef GLOBVARS
#include "imamp.h"

/* JH added 9/17/09 
based on IMadeopt program 
uses differential evolution to find joint posterior density for full and nested models */ 

/*

Differential evolution algorithm
--------------------------------

This is a conventional differential evolution algorithm.  
Many things were tried to see if it could be improved (e.g. multiple populations with migration, 
seeding population with diverse runs using smaller samples of genealogies, 
lots of playing with the weight and recombination terms of the algorithm,
(0.8 and 0.9 seem best, respectively)
lots of playing with the size of the population (100 seems best).
All to no avail.  Don't mess with it. 


Models for joint parameter estimation:
--------------------------------------
for 2 population models all tests are done against the full 5 parameter set
All nested models are tested against this full model. 
A nested model may include both p and m terms. 

For models with more than two populations,  two models are developed
1)The full model for population size parameters includes all of the population sizes
2)The full model  for the migration parameters includes only all of the migration rate parameters. 
A nested model may include either p terms or m terms, but not both. 
Depending on whether it includes p or m terms it is tested against the corresponding full model. 

Define model types as of three types:
0)  both m and p
1)  only p
2)  only m

Formatting of the nested model file 
-----------------------------------
the model file consists of lines each of which begins with one of three words ('model' 'equal' or 'constant'):

the first line is an integer with the number of models to be read in 

model   - indicates that a new model is being specified,  the string 'model ' is followed by the name of the model, a user specified string
  -all text that follows this up to the end of the line is treated as the name of the model, each time this word appears, a new model is started
equal   - a set of parameters that are always the same in value
  - followed by either an 'm' or a 'p'  for migration or population
  - then followed by two or more integers for the parameters that are equal to each other
  - each integer refers to a column position in the ti file whether or not they are population 
  - if 'p' the integers are the ids for the populations
  - if 'm' the integers are the ids for the migration parameters, in the same order as they appear in full IMa2 analysis. 
  - if 's'  the it is a constant term for the split time parameter 
constant -  set of parameters that are equal to a constant
  - followed by either an 'm' or a 'p'  for migration or population
  - then followed by a number, that is the value that the parameters are fixed at
  - then followed by one or more integers that indicate which parameters are fixed at that value
  - if 'p' the integers are the ids for the populations
  - if 'm' the integers are the ids for the migration parameters, in the same order as they appear in full IMa2 analysis. 

If there are two populations,  then having one or more nested models invokes a full model of all 5 parameters as the first model to be run.

If there are three or more populations, then a nested model can have only p terms or m terms but not both. 
If a nested model has m terms then a full model of all migration parameters is also run. 
If a nested model has p terms then a full model of all population size parameters is also run

e.g. for a two population case this is a nested model with both m and p terms in the nested model
1
model  sizes of populations  0, 1 and 2 equal to each other, no migration between 0 and 1 
equal p 0  1  2    // meaning that 1st, 2nd and third populations are set to always have the same parameter values 
constant m  0.0  0 1  // meaning that the first and second migration parameters in the model are set to a constant value of 0.0

understanding xmap
------------------
xmap contains information on how to copy one array for a nested model into a larger array for the full model
the length of xmap is as many parameters as there are in the model.
if xmap[i] is a negative value,  then the absolute value of this is plugged into the full model array at position i as a constant
if xmap[i] is a positive value, then that value is treated as an integer and position i in the full model array is copied over 
  from position xmap[i] in the nested model array 

the differential evolution algorithm works on an array of varying parameters
but if we are using a nested model the length of this list will be shorter than 
the total number of parameters (because some are shared or identical). 

parameters are counted:  first population size,  then migration, and finally a splitting rate if it is in the model. 
total number of parameters in a full model is nparams

Using xmap to generate a full list for function calculation, using a shorter list of varying numbers:
  given a list f of varying numbers,  corresponding to the things that can vary in a nested model
  if we want to generate a full length list x,  that has numbers in all parameter positions: 
	for (i=0;i<nparams;i++)
	{
		if (xmap[i] < 0)
			x[i] = -xmap[i];
		else
			x[i] = f[(int) xmap[i]];
	}
Using xmap to make a shorter list of varying numbers under a nested model,  given a full list of nparam numbers
  the full list is x and shorter list is f:
	for (i = 0,toi = 0;i<nparams;i++)
	{
		{
			if (toi == (int) xmap[i])
			{
				f[toi] = x[i];
				toi++;
			}
		}
	}    

p: 0 1 2  m: 3  4   
e.g. 2 pop model with both migration rates zero

[0,1,2,-MINPARAMVAL,-MINPARAMVAL]

e.g. 2 pops with all population sizes equal
[0,0,0,1,2]
    
e.g. 2 pop model with 3 popsize, 2 migration and one splitrate term

p: 0 1 2  m: 3  4   s:  5

xmap for a nested model with equal population sizes, freely varying migration and splitrate
[0,0,0,1,2,3]

xmap for a nested model with two equal migration rates,  other things vary
[0,1,2,3,3,4]

xmap for a nested model with a constant for population 1 of 2.5 and a constant for splitrate of 0.7
[0,-2.5,1,2,3,-0.7]

*/


#define STARTHI  1E200
#define STARTLOW  -1E200
#define SPREADTOL 1.0e-7 // criteria for stopping search for peak
#define PLOOPTOL  1.0e-6  // criteria for whether or not current peak is the same as best previously found
#define LOOPMATCHCRITERIA 2  // number of times the best peak must be found before stopping
#define MAXRESTART 10     // maximum number of restarts in search for peak 
#define MAXJOINTMODELS 100  // maximum number of models in nested model file
#define MAXPLENGTH 200   // largest number of parameters in the model, should be good for models with up to 10 sampled populations 
#define MAXMODELTEXTLINE 500 // max length of name of model in nestedmodelfile
#define DEFAULTRECRATE 0.9  // recombination rate in differential evolution algorithm.  played wth lots of values this seems best compromise
#define DEFAULTFWEIGHT 0.8  // mutation weight term in differential evolution algorithm.  played wth lots of values this seems best compromise
#define DEFAULTPOPSIZEMULTIPLIER 100  // population size multiplier in differential evolution algorithm.  played wth lots of values this seems best compromise
// VS DEBUG - just fixed deafaultpopsizemultiplier to 1
//#define DEFAULTPOPSIZEMULTIPLIER 1  // population size multiplier in differential evolution algorithm.  played wth lots of values this seems best compromise
#define MAXPTERMS 20  // maximum number of 'constant' and 'equal' lines in a nested model in the nested model file

/* SANGCHUL: Thu Nov 12 22:17:34 EST 2009
 * There are too many global variables. Could these be file static by prepending
 * static? 
 */ 
static double *bestvals, *tempvals;
static double *parameter_lower_bound;
static double *parameter_lower_bound; // VS- these are the prior distributions
static double *parameter_upper_bound;
static double **pop, **trialpop; 
// VS - the notion of population here is in the context of the genetic algorithm used to find the maximum peak
// it is not related with the populations in IM!!!
// we start with a "population" of random values from the priors
// the genetic algorithm then proceeds iteration by iteration (generation by generation)
// changing this "population" until it finds the maximum of the posterior probability
static double *popnest;
static int nparams;
static int popsizemultiplier = DEFAULTPOPSIZEMULTIPLIER;
static int num_nestedmodels;
static int depopsize;
static double holdxmaps[MAXJOINTMODELS][MAXPLENGTH];
static int paramsnotused[MAXJOINTMODELS][MAXPLENGTH]= {{0}};
static int paramsallused[MAXPLENGTH]= {0};
static int numpopsizeparams_nested[MAXJOINTMODELS];
static int nummigrateparams_nested[MAXJOINTMODELS];
static int nparams_nested[MAXJOINTMODELS];
static int nestedmodel_boundary_in_model[MAXJOINTMODELS] ={0};
static int atleastone_boundary_in_model=0;
static char modelnamelines[MAXJOINTMODELS][MAXMODELTEXTLINE];
static double recrate = DEFAULTRECRATE;
static double fweight = DEFAULTFWEIGHT;
static char logpfstr[40];  
static int calculate_ess = 0; // VS - this variable is used to compute how many genealogies are actually being used to compute the max(Likelihood)
static double effective_sample_size;  // VS - this variable is used to compute how many genealogies are actually being used to compute the max(Likelihood)
static double *logx, *divx, *log2diffx;//used by jointp
static double **double_gsamp_fcp;
static double *probgp_popsize, *probgp_migrate;
static int fullmodeltype[MAXJOINTMODELS] = {-1};
static int modeltypecounts[3] = {0}; // 0, both m and p;  1, only p;  2, only m 
static char modelstartstr[3][40] = {"FULL","ALL POPULATION SIZE PARAMETERS", "ALL MIGRATION PARAMETERS"};
static int nparamrange[2]; 
static int nowmodeltype; 
static double parameter_lower_bound_mapped[MAXPLENGTH];
static double parameter_upper_bound_mapped[MAXPLENGTH];

  
// listelement is for a linked list
struct listelement 
{
  double v;
  struct listelement *l;
  struct listelement *r;
}**list, *first, *last, *now; 

/* prototypes */
char *logpstrformat (double pval);
int fillplist(char *c, int *plist, int startparamcount, int notused[]);
int reduce_mappos(int mappos[], int k, int len);
void setup_mapping(char *nestedmodelfilename);
void reversemapvals(double *from, double *to, double *xmap);
void mapvals(double *from, double *to, double *xmap);
void setbounds();
void nextgen(double *lowpd, double *hipd, int modelparams, int modelnum);
void startpop(int modelparams,int modelnum);
double difeloop(int modelparams,int modelnum);
void modelloop(int modelnum);
void copybest(/* int modelparams */);
void printjointpeakvals(FILE *outfile,int notused[], double *printvals);
void startjointpeakouttable(FILE *outfile,char *fname);
double jointp (double *x,int calc_ess, double *effective_n);
void freejointpmem();
void addfirst(double val, int k);
void addlast(double val, int k);
void addright(double val, int k);
void addleft(double val, int k);
void listput(double val, int k);
void setuplist();

/****** LOCAL FUNCTIONS *********/

/* functions for linkedlist */
/* simple linkedlist,  with the biggest value at the right end, where last always points */
void addfirst(double val, int k)
{
  list[k]->v = val;
  list[k]->r = first;
  list[k]->l = NULL;
  first->l = list[k];
  first = list[k];
}

void addlast(double val, int k)
{
  list[k]->v = val;
  list[k]->r = NULL;
  list[k]->l = last;
  last->r = list[k];
  last = list[k];
}

void addright(double val, int k)
{
  list[k]->v = val;
  list[k]->r = now->r;
  list[k]->l = now;
  now->r->l = list[k];
  now->r = list[k];
  now = list[k];
}

void addleft(double val, int k)
{
  list[k]->v = val;
  list[k]->l = now->l;
  list[k]->r = now;
  now->l->r = list[k];
  now->l = list[k];
  now = list[k];
}

void listput(double val, int k)
{
  if (val <= now->v)
  {
    while (val <= now->v && now != first)
      now = now->l; 
    if (now == first && val <= now->v)
      addfirst(val,k);
    else
      addright(val,k);
  }
  else
  {
    while (val > now->v && now != last)
      now = now->r; 
    if (now == last && val > now->v)
      addlast(val,k);
    else
      addleft (val, k);
  }
} // listput

// setuplist
// allocates memory for the list that will be used to save the values of the posterior probability
// INPUT
//	void
// RETURN
//  allocates memory for the list that will be used to save the values of the posterior probability
//	during the optimization algorithm to find the maximum posterior
/* make enough of these listelements at the beginning rather than malloc and free them as we go along */
void setuplist()
{
  int i;
  list = (struct listelement **) malloc(genealogiessaved * sizeof (struct listelement *));
  for (i=0;i<genealogiessaved;i++)
    list[i] = (struct listelement *) malloc(sizeof (struct listelement));
}

/* for formatting a floating point number */ 
char *logpstrformat (double pval)
{
  logpfstr[0] = '\0';
  if (fabs (pval) < 1e-2)
    sprintf (logpfstr, "%.6lf", pval);
  else if (fabs (pval) < 1e-1)
    sprintf (logpfstr, "%.5lf", pval);
  else if (fabs (pval) < 1e-0)
    sprintf (logpfstr, "%.4lf", pval);
  else if (fabs (pval) < 1e1)
    sprintf (logpfstr, "%.3lf", pval);
  else if (fabs (pval) < 1e2)
    sprintf (logpfstr, "%.2lf", pval);
  else if (fabs (pval) < 1e3)
    sprintf (logpfstr, "%.1lf", pval);
  else if (fabs (pval) < 1e4)
    sprintf (logpfstr, "%.0lf", pval);
  return &logpfstr[0];
} 


// fillplist
// this function will fill plist[][] and notused[] and return the nb of elements filled, i.e. number of parameters set "equal" or "constant"
// INPUT:
//	char *c : string with the index of the parameter that is set equal to another one
//	int *plist : pointer to plist that will be filled. plist has the list of the index of parameters that are set equal or to a constant.
//	int startparamcount : begining of the index of param, e.g. for theta=0 for mig=numpopsizeparam
//	int notused[] : array with a list of parameters not used in a given model. Entries with 1 mean that the parameter at that index are set "equal" or "constant"
// RETURN
//  return the number of elements that were filled in plist. This should be the same as the number of parameters that are set equal
int fillplist(char *c, int *plist, int startparamcount, int notused[])
{
	int i;
	i = 0;
	if (isspace(*c))
		c = nextnonspace(c);
	while (c && sscanf (c, "%d",&plist[i])>0)
	{
      plist[i] += startparamcount;
        notused[plist[i]] = 1;
	  c = nextnonspace(c);
	  i++;
	}
	plist[i] = -1;
	return  i;
}  // fill plist



// reduce_mappos
// function that gets the map of parameters for the nested model
// the map is the index of parameters that are set to be equal or constant
// INPUT:
//	int mappos[] : array with the map of full model, e.g. mappos=[0,1,2,3,4] for a model with 2 populations (5 params)
//	int k : index of the param which is set "equal" or "constant"
//	int len : total number of parameters in the nested model??
int reduce_mappos(int mappos[], int k, int len)
{
	int i=0, j;
	while (mappos[i] != k && i < nparams) i++; 
	assert(i < nparams);
	for (j=i;j< nparams;j++)
		mappos[j] = mappos[j+1]; // VS - why is the mappos the same as mappos[j+1]??? Should it be the position of plist??
	return len - 1;
}


// fillxmap
// function that will fill the xmap variable for each nested model
// INPUT
//	double xmap[] : array of doubles to save the mapping
//	int pos[] : ????
//	int pl[][MAXPLENGTH] :
//	double terms[] :
//	int len : number of parameters in the full model
//	int type[] : vector with the type of nested model - constant parameter value or sharing the same parameters
//	int num
// Fillxmap for the groups of loci the solution implemented
// is to use the existing format and just change the variable nparams. 
// Thus, in the input file for the nested model we specify which 
// parameters are the same and that is defined in the index of xmap. 
// An example of the input file would be:
// # models in which all migration rate among groups are equal
// model equal migration rates among groups
// equal m 0 2 
// equal m 1 3
// This is assuming that the migration rate parameters are included as m12_g1 m21_g1 m12_g2 m21_g2
// NOTE: that it seems that the order of the migration parameters is actually m12_g1 m12_g2 m21_g1 m21_g2!!
void fillxmap(double xmap[], int pos[], int pl[][MAXPLENGTH], double terms[], int len, int type[], int num)
{
  int i, j;
  /* put positions into xmap of parameters that are still in the model  */
  for (i=0;i<nparams;i++)
  {
    j = 0;
    /*CR 110912.1 change compound test order to eliminate invalid mem access */
    while ( j < len && pos[j] != i ) j++;
    if (j<len)
      xmap[i] = (double) j;
  }
  /* now put constants and other positions into xmap  (use negative values for consants ) */
  for (i=0;i<num;i++)
  {
    if (type[i])
    {
      j = 0;
      while (pl[i][j] != -1)
      {
        xmap[pl[i][j]] = -terms[i];  // negative of constant is stored 
        j++;
      }
    }
    else
    {
      j = 0;
      while (pl[i][j] != -1)
      {
        xmap[pl[i][j]] = (double) (int) xmap[(int) terms[i]];
        j++;
      }
    }
  }
} //fillxmap


// setup_mapping
// JH reads the model file and sets up the arrays that determine how a full model maps onto a nested model
// JH really horrible code in this function //
// INPUT:
//	char *fname : string with the filename of the input file with the nested model information
// VS - the xmap when dealing with tests with more than 1 group of loci will need to include all the parameters
// e.g. 2 populations sampled, 1 group for theta, 2 groups for migration
// XMAP=[1,2,3,4,5,6,7] for the full model
// corresponding to q1, q2, qA, m12_g1, m12_g2, m21_g1, m21_g2
// where g1 anf g2 refer to the groups of loci for migration parameters.
// The nested model where the two groups share the same parameter would be coded as
// XMAP=[1,2,3,4,4,5,5] for the nested model
// in this case m12_g1=m12_g2 and m21_g1=m21_g2
// NOTE: given that this function reads a file, and it is very complex,
// for the moment I will just fix the xmap for the example above
// this is a temporary solution to see if the LRT works for groups of loci
// if the algorithm is working after tested with this case, then I could change the way
// to code the nested input file and how to change this function accordingly
void setup_mapping(char *fname)// JH fixed signficiant bug in this functino 4/5/2010 
{
  FILE *nestedmodelfile;
  char *modeltextline;
  char keyword[12];
  char *c, m_or_p;
  int  mappos[MAXPLENGTH], mappos_len;
  int plist[MAXPTERMS][MAXPLENGTH];
  double pterms[MAXPTERMS]; // first number given on a constant or equal line,  if constant line it is a constant, else it is the first parameter number
  int ptype[MAXPTERMS];
  int i, n,numm, inmodel, thismodeltype;

  if ((nestedmodelfile = fopen (fname, "r")) == NULL)
  {
    IM_err(IMERR_READFILEOPENFAIL,"Error opening nested model file: %s", fname);
  }
  modeltextline = malloc(MAXMODELTEXTLINE*sizeof(char));
  do
  {
    fgets(modeltextline,MAXMODELTEXTLINE,nestedmodelfile);
  } while (!isdigit(modeltextline[0]));
  sscanf (modeltextline, "%d", &num_nestedmodels);
  inmodel = -1;
  numm = -1;   // JH fixed bug not sure if it caused problems
  while (fgets(modeltextline,MAXMODELTEXTLINE,nestedmodelfile)!= NULL) if (isalpha(modeltextline[0]))
  {
    for (i=0;i< (int) strlen(modeltextline);i++) // make sure text beginning line is lowercase
    {
      if (isspace(modeltextline[i]))   break;
      modeltextline[i] = tolower(modeltextline[i]);
    }
    c = modeltextline;
    sscanf (c, "%s", keyword);
    c = nextnonspace(c);
    if (strcmp(keyword,"model")==0)
    {
      if (inmodel >= 0) // finish up the previous model
      {
        if (numm < 1)
            IM_err (IMERR_NESTEDMODELLSPECIFYLFAIL,"nested model %d does not include any 'constant' or 'equal' specifications",inmodel);
        if (npops ==2 && fullmodeltype[inmodel]!= 0 && vs_test==0) // VS vs_test 11/7/2011 
            IM_err (IMERR_NESTEDMODELLSPECIFYLFAIL,"nested model %d has m or p type but should be both when there are only two sampled populations",inmodel);
        if (npops > 2 && fullmodeltype[inmodel]== 0 && vs_test==1) // VS vs_test 11/7/2011 
            IM_err (IMERR_NESTEDMODELLSPECIFYLFAIL,"nested model %d has type 0 but should be m or p type",inmodel);
        modeltypecounts[fullmodeltype[inmodel]]++;
		// VS fillxmap - this is the function where XMAP is initialized defined
		// VS fillxmap initializes mappos, plist, pterms, mappos_len, ptype?
		// mappos, plist, pterms, mappos_len are defined below, when reading the full model
        fillxmap(holdxmaps[inmodel], mappos,plist,pterms,mappos_len,/*nparams_nested[inmodel],*/ptype,numm);
      }
      inmodel++;
      strcpy(modelnamelines[inmodel],c);
      if (strlen(modelnamelines[inmodel]) && modelnamelines[inmodel][strlen(modelnamelines[inmodel])-1] == '\n') 
      {
        modelnamelines[inmodel][strlen(modelnamelines[inmodel])-1] = 0;
      }
      numm=0;
	  // VS - need to change the numpopsizeparams_nested[inmodel]
	  // VS - need to change the nummigrateparams_nested[inmodel]
      numpopsizeparams_nested[inmodel] = numpopsizeparams*nbgroupsloci_theta;
      nummigrateparams_nested[inmodel] = nummigrateparams*nbgroupsloci_mig;
      mappos_len = nparams;
      for (i=0;i<nparams;i++)
      {
	    mappos[i] = i;
	    holdxmaps[inmodel][i] = -1;
      }
      //if (npops == 2) // VS 10/19/2011 test version where the migration parameters are searched independently even with two populations
	  if(npops == 2 && vs_test==0)
      {
        fullmodeltype[inmodel] = 0;
        thismodeltype = 0;  //JH  fixed a bug 5/10/2010 that popped up under cygwin/linuxv
        nparams_nested[inmodel] = nparams;
      }
      else
      {
        thismodeltype = -1; // VS if with 3 or more populations thismodeltype=-1
      }
    }
    if (strcmp(keyword,"equal")==0 || (strcmp(keyword,"constant")==0)) 
    {
      sscanf (c, "%c",&m_or_p);
      if (thismodeltype == -1)
      {
        if (m_or_p == 'm')
        {
          thismodeltype = 2; //jh fixed bug 5/27/2010
          fullmodeltype[inmodel] = 2; 
          nparams_nested[inmodel] = nummigrateparams*nbgroupsloci_mig;
        }
        if (m_or_p == 'p')
        {
          thismodeltype = 1; //jh fixed bug 5/27/2010
          fullmodeltype[inmodel] = 1;
          nparams_nested[inmodel] = numpopsizeparams*nbgroupsloci_theta;
        }
      }
      else // VS if looking at the 2 populations case?
      {
        if (m_or_p == 'm' && thismodeltype == 1)
            IM_err (IMERR_NESTEDMODELLSPECIFYLFAIL,"nested model %d specified as both p and m types",inmodel);
        if (m_or_p == 'p' && thismodeltype == 2)
            IM_err (IMERR_NESTEDMODELLSPECIFYLFAIL,"nested model %d specified as both p and m types",inmodel);
      }
      c = nextnonspace(c);
      sscanf (c, "%lf",&pterms[numm]); // first #, can be either a constant or a parameter number 
	  // VS pterms stores the indexes of the xmap that are set to be equal or constant
	  // in the case of "equal" the 1st and 2nd elements store the indexes of the equal parameters
	  // in the case of "constant" the first element stores the index and the second the value of the constant
	  if (m_or_p == 'm' && strcmp(keyword,"equal")==0) { //then first element of pterms is a parameter number, else its a constant and not changed
      //  pterms[numm] += numpopsizeparams; // VS - if it is a migration parameter need to add numpopsizeparams
		  pterms[numm] += numpopsizeparams*nbgroupsloci_theta; // VS - if it is a migration parameter need to add numpopsizeparams
	  }
      c = nextwhite(c);
      if (strcmp(keyword,"equal")==0) // found a match
      {
        ptype[numm] = 0; // VS - ptype==0 if parameters are equal
						 // NOTE that this is different from pterms!!
      }
      else  // "constant"
      {
        ptype[numm] = 1;  // VS - ptype==1 if parameters are constant
        if (pterms[numm] <= MINPARAMVAL ) 
        {
	        pterms[numm] = (double) MINPARAMVAL;
            nestedmodel_boundary_in_model[inmodel] =1 ; // VS - What is this??? It seems that this is a variable used to check if the parameters are within the prior limits of the full model
        }
        else
        {
          if (!calcoptions[USEPRIORFILE]) //unlikely situation of setting constant values greater than the prior, but having mig priors read from file
          {
            if (m_or_p == 'm' && pterms[numm] >= mprior) // VS - need to check if mprior is the same as the prior for the groups of mig!!
            {
              pterms[numm] = mprior;
              nestedmodel_boundary_in_model[inmodel] = 1;
            }
            if (m_or_p == 'p' && pterms[numm] >= thetaprior) // VS - need to check if thetaprior is the same as the prior for the groups of mig!!
            {
              pterms[numm] = thetaprior;
              nestedmodel_boundary_in_model[inmodel] = 1;
            }
          }
        }
		// VS note that |= is the Bitwise OR assignment
        atleastone_boundary_in_model |= nestedmodel_boundary_in_model[inmodel];
      }
      if (m_or_p == 'p')
        n = fillplist(c,&plist[numm][0], 0, paramsnotused[inmodel]); // VS plist will store the array of index of parameters that are set equal to other param
					// e.g. in the case that we have q1=q2=q3 the plist would store the index of those
      if (m_or_p == 'm')
        n = fillplist(c,&plist[numm][0],numpopsizeparams*nbgroupsloci_theta,paramsnotused[inmodel]);
      // VS n variable contains the number of paramaters set equal or to constant in model inmodel
	  for (i=0;i<n;i++)
      {
        if (m_or_p == 'p')
          numpopsizeparams_nested[inmodel]--;
        else
          nummigrateparams_nested[inmodel]--;
		// VS reduce_mappos is an important function
        nparams_nested[inmodel] =  reduce_mappos(mappos,plist[numm][i],nparams_nested[inmodel]);
        mappos_len--;
      }
      numm++;
    }
  }
  if (inmodel >= 0) // finish up the last model
  {
    if (numm < 1)
      IM_err (IMERR_NESTEDMODELLSPECIFYLFAIL,"nested model %d does not include any 'constant' or 'equal' specifications",inmodel);
    // VS vs_test 10/19/2011
	// if (npops ==2 && fullmodeltype[inmodel]!= 0)
	if (npops ==2 && vs_test==0 && fullmodeltype[inmodel]!= 0) // VS vs_test
        IM_err (IMERR_NESTEDMODELLSPECIFYLFAIL,"nested model %d has m or p type but should be both",inmodel);
    // VS vs_test 10/19/2011
	//if (npops > 2 && fullmodeltype[inmodel]== 0)
	if ((npops > 2 || vs_test==1) && fullmodeltype[inmodel]== 0) // VS vs_test
        IM_err (IMERR_NESTEDMODELLSPECIFYLFAIL,"nested model %d has type 0 but should be m or p type",inmodel);
    modeltypecounts[fullmodeltype[inmodel]]++;      
    fillxmap(holdxmaps[inmodel], mappos,plist,pterms,mappos_len,/*nparams_nested[inmodel],*/ptype,numm);
  }

  XFREE(modeltextline);
  f_close(nestedmodelfile);
}// setup_mapping;


// reversemapvals
// from a full parameter set down to a nested one
// INPUT:
//	double *from : array with index in the full model of the parameter that becomes parameter of index to in nested model
//  double *to : array with index in the nested model of the parameter
//	double *xmap : xmap array
// In the 2 groups for migration, the example would be to model
//  q1 q2 qA m12a m21a m12b m21b -> XMAP=(0,1,2,3,4,5,6) = FULL MODEL
// In the nested model where a and b have the same parameters
//	q1 q2 qA m12 m21 m12 m21	 -> XMAP=(0,1,2,3,4,3,4) = NESTED MODEL
void reversemapvals(double *from, double *to, double *xmap)
{
  int i, toi;

  /*  */
  for (i = 0,toi = 0;i<nparams;i++)
  {
    if (toi == (int) xmap[i])
    {
      to[toi] = from[i];
      toi++;
    }
  }
} //reversemapvals

//from a nested parameter set to a full one
// VS INPUT:
//  double *from: this is an array with the????
//  double *to: this is an array with the???
//  double *xamp: this is the pointer to the xmap array
// In the 2 groups for migration, the example would be to model
//  q1 q2 qA m12a m21a m12b m21b -> XMAP=(0,1,2,3,4,5,6) = FULL MODEL
// In the nested model where a and b have the same parameters
//	q1 q2 qA m12 m21 m12 m21 -> XMAP=(0,1,2,3,4,3,4)
void mapvals(double *from, double *to, double *xmap)
{
  int i = 0;
  for (i=0;i<nparams;i++)
  {
    if (xmap[i] < 0) // if xmap is negative it means that the parameter in the nested model is equal to a constant
      to[i] = -xmap[i]; // so the value of the parameter at index i is set to the constant
    else
      to[i] = from[(int) xmap[i]]; // if it is positive, then the xmap will have the index of the parameter to which it is the same
  }
} //mapvals */


// setbounds
// function that will set the upper and lower limit of the parameter values according to the prior distributions
// INPUT:
//	void function
// RETURNS:
//  initializes the parameter_lower_bound and parameter_upper_bound
// NOTE: the Likelihood Ratio Tests (LRT) can only be done when all groups have the same prior distribution
// the size of the array paramater_lower_bound is nparams
// when we have groups of loci we need to re-compute the nparams
// e.g. in a case where nbgrouptheta=1 and nbgroupmig=2 with 2 populations
// there will be 7 parameters instead of the 5 parameters in the case where there is only one group
void setbounds()
{
  int i;
  int gp, paramindex; // VS gp
  int nbthetaparams; // VS

  // VS compute some variables that will simplify the code in the for loops and allocation of variables
  nbthetaparams=numpopsizeparams*nbgroupsloci_theta;

  parameter_lower_bound = malloc(nparams * sizeof(double));
  parameter_upper_bound = malloc(nparams * sizeof(double));
  for (i = 0; i < nparams; ++i) // VS - note that nparams=(numpopsizeparams * nbgroupsloci_theta) + (numpopmigparams * nbgroupsloci_mig)
  {
    parameter_lower_bound[i] =  MINPARAMVAL;
    // VS if (i < numpopsizeparams)
	if (i < nbthetaparams) // this is a theta parameter
    {
	  // VS find the group and iparam looking at the global variables g_paramindex and g_gp_ip
      paramindex = g_paramindex[i];
	  gp = g_gp_ip[i];
	
      //parameter_upper_bound[i] = itheta[i].pr[gp].max; // VS gp
	  parameter_upper_bound[i] = itheta[paramindex].pr[gp].max; // VS gp, paramindex
    }
    else // else this is a migration parameter
    {
	  // VS find the group and iparam looking at the global variables g_paramindex and g_gp_ip
      paramindex = g_paramindex[i];
	  gp = g_gp_ip[i];

      // parameter_upper_bound[i] = imig[i-numpopsizeparams].pr[gp].max; // VS gp
	  parameter_upper_bound[i] = imig[paramindex].pr[gp].max; // VS gp, paramindex
    }
  }
} /* setbounds */
        

/* applies the differential evolution algorithm */
void nextgen(double *lowpd, double *hipd, int modelparams,int modelnum)
{
  int A,B,Cvalue, i,j, jmax;
  double temp;
  double u;

  *lowpd = STARTHI;
  *hipd = STARTLOW;

  jmax = nparamrange[0] + modelparams;
  assert(jmax <= nparamrange[1]); 
  /*   CR 110929.2 from JH  9/26/2011 note that for nested models this may 
   *   not use the right bounds, if bounds vary among parameters in 
   *   a class - this is abug 
   */
  for (i=0;i<depopsize;i++)
  {
    A = randposint(depopsize);
    B = randposint(depopsize);
    Cvalue = randposint(depopsize);
    //for (j=nparamrange[0];j<nparamrange[1];j++) 
    for (j=nparamrange[0];j<jmax;j++) 
    {
      u = uniform();
      if (u < recrate)
      {
      temp = pop[Cvalue][j] + fweight* (pop[A][j] - pop[B][j]);
      if (temp < parameter_lower_bound_mapped[j]) // move only part of the way towards the lower bound
        temp = pop[Cvalue][j] -  uniform() * (pop[Cvalue][j] - parameter_lower_bound_mapped[j]);
      assert ( temp >= parameter_lower_bound_mapped[j]);
      if (temp > parameter_upper_bound_mapped[j])// move only part of the way towards the upper bound
        temp = pop[Cvalue][j] +  uniform() * (parameter_upper_bound_mapped[j]- pop[Cvalue][j]);
      assert(temp <= parameter_upper_bound_mapped[j]); 
      trialpop[i][j] = temp;
      }
      else
        trialpop[i][j] = pop[i][j];
    }
    if (modelnum >= 0)
    {
      mapvals(trialpop[i],popnest,holdxmaps[modelnum]);
      trialpop[i][nparams] = jointp(popnest, calculate_ess,&effective_sample_size);
    }
    else
      trialpop[i][nparams] = jointp(trialpop[i], calculate_ess,&effective_sample_size);
  }
  for (i=0;i<depopsize;i++)
  {
    if (trialpop[i][nparams] < pop[i][nparams])
      for (j=0;j<=nparams;j++)
        pop[i][j] = trialpop[i][j];
    if (pop[i][nparams] < *lowpd)
      *lowpd = pop[i][nparams];
    if (pop[i][nparams] > *hipd)
       *hipd = pop[i][nparams];
  }

} //nextgen

// startpop
// initializes the "population" of paramater values to update during algorithm to find maximum
// the notion of population here is in the context of the genetic algorithm used to find the maximum peak
// it is not related with the populations in IM!!!
// we start with a "population" of random values from the priors
// the genetic algorithm then proceeds iteration by iteration (generation by generation)
// changing this "population" until it finds the maximum!!!
// INPUT:
//	int modelparams : number of parameters in the model
//	int modelnum : index of the model
// create a starting population by sampling uniformally over the priors
void startpop(int modelparams,int modelnum)
{
  int i,j, jmax;
  jmax = nparamrange[0] + modelparams; // what is nparamrange storing?
  
  // Depending on whether this is the full model or the nested model
  // we have a different set of parameters. Need to check here these lower and upper bounds
  // and see how this is related with the way the parameter search is done.
  if (modelnum < 0) for (j=0;j<nparams;j++)
  {
    parameter_lower_bound_mapped[j]=parameter_lower_bound[j];
    parameter_upper_bound_mapped[j]=parameter_upper_bound[j];
  }
  else // if we are dealing with a nested model
  {
    // the lower and upper bound are defined by the number of parameters
    // So, the parameter_lower_bound_mapped is an array with the index of the parameters
    reversemapvals(parameter_lower_bound,parameter_lower_bound_mapped,holdxmaps[modelnum]);
    reversemapvals(parameter_upper_bound,parameter_upper_bound_mapped,holdxmaps[modelnum]); 
  }

  // VS - for each element of the "population" of values in the optimization algorithm
  // we call the jointp function
  for (i=0;i<depopsize;i++) 
  {
	// go through the parameters
    for (j=0;j<nparams;j++)
    {
                                /* CR 110929.2 from JH 9/26/2011 this 
                                 * turned off because  parameter vals will be 
                                 * mapped into popnest */
      //if (j>=nparamrange[0] && j<nparamrange[1]) 
      if (j>=nparamrange[0] && j<jmax)   /* CR 110929.2 from JH 9/26/2011 but 
                                          * for nested models this may not 
                                          * use the right bounds, if bounds 
                                          * vary among parameters in a 
                                          * class - this is abug */
        pop[i][j] = parameter_lower_bound_mapped[j] + uniform()*(parameter_upper_bound_mapped[j] - parameter_lower_bound_mapped[j]);
      else
        pop[i][j] = -1;
    }
    if (modelnum >= 0)
    {
      mapvals(pop[i],popnest,holdxmaps[modelnum]);
      pop[i][nparams] = jointp(popnest, calculate_ess,&effective_sample_size);
    }
    else
	  // the jointp function is called by each element of the population
      pop[i][nparams] = jointp(pop[i], calculate_ess,&effective_sample_size);
  }
}// startpop

/* keeps calling nextgen() until a peak has been found */ 
// VS how is a peak defined? 
// VS Is it when all individuals of a pop have the same value?
// VS my understanding is that this happens when all individuals in the pop
// VS have a probability whose range or dispersal is lower than SPREADTOL
double difeloop(int modelparams, int modelnum)
{
  int gen;
  double spread, lowpd, hipd;
  /* static double pd; Never used */

  gen = 0;
  do
  {
	// VS - nextgen will apply the genetic algorithm to sample the new
	// population of values. I believe this is done based on the 
	// posterior probability of each element in this population.
    nextgen(&lowpd, &hipd,modelparams,modelnum);
    spread = hipd - lowpd;
    gen++;
  }
  while (spread > SPREADTOL);
  return  lowpd;
} // difeloop

/* find the best individual in the population and copy it into bestvals[] */
void copybest(/* int modelparams */)
{
  int i,k;
  double temp = STARTHI;
  for (i=0;i<depopsize;i++)
  {
    if (pop[i][nparams] < temp)
    {
      temp = pop[i][nparams];
      k = i;
    }
  }
  for (i=0;i<=nparams;i++)
    bestvals[i] = pop[k][i];
  bestvals[nparams] = temp;
}//copybest


// modelloop(int modelnum)
// find a peak,  record location, then begin a loop that repeats the search
// - when a better peak is found it is saved
// - if the same best peak has been found more than LOOPMATCHCRITERIA times 
// - or if restarts have been done MAXRESTART times, it finishes.
// INPUT
//	int modelnum: index of model
void modelloop(int modelnum)
{
  int i,newstartloop;
  int modelparams;
  int countloop;
  double local_pd;
  double global_pd = STARTHI;
  double jointProb;  /* CR 110929.5  variable name changed */

  if (modelnum < 0)
  {
      switch (nowmodeltype)
      {
        case 0 : modelparams =  nparams; break;
        case 1 : modelparams =  numpopsizeparams*nbgroupsloci_theta;break; // VS added the nbgroupsloci_theta
        case 2 : modelparams =  nummigrateparams*nbgroupsloci_mig;break; // VS added the nbgroupsloci_mig
      }
  }
  else
  {
    modelparams = nparams_nested[modelnum]; // VS - modelparams has the number of parameters of the current model
  }
  depopsize = modelparams * popsizemultiplier;
  // VS - startpop initializes the "population" of values of the parameters
  // that will be used in the genetic algorithm to find the peak
  startpop(modelparams,modelnum); // VS startpop calls the jointp function
  
  // VS - depopsize is the size of elements in the population of parameter values
  // the following loop finds the maximum global posterior probability
  // recall that pop is a matrix with the valures of the parameters in each column
  // the last column saves the posterior probability
  for (i=0; i< depopsize;i++)
  {
    if (pop[i][nparams] < global_pd)
      global_pd = pop[i][nparams];
  }

  newstartloop = 0;
  countloop = 0;
  do
  {
    if (newstartloop > 0)
    {
      startpop(modelparams,modelnum);
      printf("         restart #%d (out of %d maximum)\n",newstartloop,MAXRESTART);
    }
	// VS - the local pd is???
    local_pd = difeloop(modelparams,modelnum); // VS - difeloop calls the nextgen function
    if (fabs(local_pd - global_pd) < PLOOPTOL)
      countloop++;
    else
    {
      if (local_pd < global_pd)
        countloop = 0;
    }
    if (local_pd < global_pd)
    {
      copybest(/* modelparams Does nothing */); // VS - copy the best values to the bestvals array
      global_pd  = local_pd;
    }
    newstartloop ++;
  }
  while (countloop < LOOPMATCHCRITERIA && newstartloop < MAXRESTART); 
  calculate_ess = 1;

  /*  CR 110929.5  although return val from jointP() call is not used 
   *  locally, it may be useful during debuggin, also the return 
   *  variable name was changed 
   */
  if (modelnum >= 0)
  {
    mapvals(bestvals,popnest,holdxmaps[modelnum]); // VS - where do we obtain the bestvals??
     jointProb = jointp(popnest,calculate_ess, &effective_sample_size);
  }
  else
     jointProb = jointp(bestvals, calculate_ess,&effective_sample_size);
  calculate_ess = 0;
} // modelloop


// printjointpeakvals
// prints the results of the joint peak values to the output file
// INPUT
//	FILE *outfile: pointer to output file
//	int notused[]: ???
//	double *printvals: ???
void printjointpeakvals(FILE *outfile,int notused[], double *printvals)
{
  int i;
  for (i=0;i<nparams;i++) // VS - recall that nparams is (numpsizeparams*nbgroupsloci_theta)+(nummigrateparams*nbgroupsloci_mig)
  {
    if (i>= nparamrange[0] && i < nparamrange[1]) // VS nparamrange[0] is usually from 0 to nparams. Maybe this varies when there is a parameter set to equal.
    {
      if (notused[i])
      {
        if (printvals[i] < 0.001)
          FP"\t[%.5lf]",printvals[i]);
        else
          FP"\t[%.4lf]",printvals[i]);
      }
      else
      {
        if (printvals[i] < 0.001)
          FP"\t%.5lf",printvals[i]);
        else
          FP"\t%.4lf",printvals[i]);
      }
    }
    else
      FP"\t-");
  }
  FP"\n");
} //printjointpeakvals



// startjointpeakouttable
// prints information about the parameters and info about the LTR tests done
// INPUT
//	FILE *outfile: pointer to the file 
//	char *fname: string with the file name
void startjointpeakouttable(FILE *outfile,char *fname)
{
  //int gp; // VS gp. Since all migration parameters
  int i, mi, typeloopstart, typeloopstop,nowmodeltype_local,j=0;

  int nbthetaparams, nbparams; // VS includes these variables that simplify reading of the code and for loops
  nbthetaparams = numpopsizeparams*nbgroupsloci_theta;
  nbparams = (numpopsizeparams*nbgroupsloci_theta)+(nummigrateparams*nbgroupsloci_mig);


  FP"Joint Peak Locations and Posterior Probabilities\n");
  FP"================================================\n");
  FP"  estimates based on %d sampled genealogies\n",genealogiessaved);
  if (strlen(fname) > 0)
    FP"  nested model filename:%s\n",fname);
  FP"\nModel#  Model Description\n");
  //if (npops == 2 ) // VS added this such that we analyse the data from 2 pop as when there are more than 2 pops
  // in that case, the joint peak is found only focusing in the joint marginal either for the theta or migration parameters
  if (npops == 2 && vs_test==0)
  {
    typeloopstart = 0;
    typeloopstop = 0;
  }
  else
  {
    typeloopstart = 1;
    typeloopstop = 2;
  }
  for (nowmodeltype_local = typeloopstart; nowmodeltype_local <= typeloopstop;nowmodeltype_local++)
  {
    j++;
    FP"%2d     %s\n",j,modelstartstr[nowmodeltype_local]);
    for (mi = 0;mi<num_nestedmodels;mi++) if (fullmodeltype[mi] == nowmodeltype_local)
    {
      j++;
      FP"%2d   %s\n",j,modelnamelines[mi]);
    }
  }

  FP"\nModel#\tlog(P)\t#terms\tdf\t2LLR\tESS");
  for (i = 0; i < nbthetaparams; i++) { // VS nbgroupsloci_theta
	  // VS add the group information, including an extra for loop
	  FP"\t%s_g%i",itheta[g_paramindex[i]].str, g_gp_ip[i]);  // VS g_paramindex[i] refers to parameters, g_gp_ip[i] refers to group index in iparam structures
  }
  for (i = nbthetaparams; i < nbparams; i++) { // VS for loop between nbthetaparams and nbparams, instead of for loop between 0 and nummigrateparams
	  if (imig[g_paramindex[i]].pr[g_gp_ip[i]].max > MPRIORMIN) // VS g_paramindex[i] and g_gp_ip[i]
			FP"\t%s_g%i",imig[g_paramindex[i]].str, g_gp_ip[i]);
	  // VS - using the global indexes g_paramindex and g_gp_ip, there is no need to include this extra for loop
  }
  FP"\n");
}

/* JH new jointp brought in 9/16/09 */
#define OCUTOFF  10
#define POW10I(a) ((a)+308)
#define PRANGELOG 10
double jointp (double *x, int calc_ess, double *effective_n)
/* calculate the joint likelihhood function  */ 
/* the point at which to calculate the value of the function is in x */ 
/* return the negative of the logarithm of the joint posterior probability at x */
/* if *effective_n == 1.0  calculate the effective number of samples */ 
/* this has been tweaked a lot for speed 
  only saves the best genealogies (uses a linked list)
  precalculates as much stuff as possible that will get reused later */ 
{
  int gi, i; //  VS int gi, i, i1; // VS (style) Unused variable: i1
  int gin, iin, ii;
  double sum, p, q;
  static int ccp, fcp, hccp, mcp, fmp, qip, mip, probgp, pdgp;
  static double log_numtrees;
  static int init = 0;
  int zadj, maxz = -10000000;
  double acumm = 0;
  float *g;
  double *dg;
  double acumm_sqr = 0;
 
  // VS gp
  int gp;
   
  static double pow10[617];
 
  int nbthetaparams, nbmigparams; // VS
  // VS compute some variables that will simplify the code in the for loops and allocation of variables
  nbthetaparams=numpopsizeparams*nbgroupsloci_theta;
  nbmigparams=nummigrateparams*nbgroupsloci_mig;


  if (init == 0)
  {
    init = 1;
    //npop_and_m = (numpopsizeparams*nbgroupsloci_theta) + (nummigrateparams*nbgroupsloci_mig); // VS - nbgroupsloci_theta and nbgroupsloci_mig
    log_numtrees = log ((double) genealogiessaved);
    /********* this section copied from IMa2  initialize.c setup_iparams()  on 8/24/09 */
	// VS - these variables are the column index of the gsampinf matrix
	// VS - they are used to refer to certain quantities in the matrix
	// One possible implementation of the groups of loci is to put all the gsamping matrices into a big matrix
	// This way we would still use these functions exactly in the same way.
	// In the case that there is only 1 group for theta and 1 group for mig, do we still save 2 gsamping matrices?
	// Maybe in that case it would be easy to put these two matrices together into 1 matrix.
	// Another approach is to change completely the settings and consider a gsampinf matrix with a different number of rows and columns.
    /* initialize.c part
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
	*/
	// VS - BUG HERE??? YES!!! as can be seen the plgp is not in the initialize.c
	// VS - maybe it would be better just to use the global variables defined in initialize.c
	ccp = 0;
    fcp = ccp + numpopsizeparams;
    hccp = fcp + numpopsizeparams;
    mcp = hccp + numpopsizeparams;
    fmp = mcp + nummigrateparams;
    qip = fmp + nummigrateparams;
    mip = qip + numpopsizeparams;
    pdgp = mip + nummigrateparams;
    /* CR 110830.1
     * probgp was not correct, the old variable plgp in jointp function 
     * (jointfind.c, line 912) was no longer used. 
     */
    probgp = pdgp + 1;
    /**********/
    /* pow10[],logx, divx, log2diffx and doubling gsampinf[gi][fcp+i] are done to speed things up*/
    for (i=-308;i<=308;i++)
      pow10[POW10I(i)] = pow(10,i);
    logx = (double *) malloc(nparams * sizeof (double)); 
    divx = (double *) malloc(nparams * sizeof (double));
    log2diffx = (double *) malloc(nbthetaparams * sizeof (double)); // VS - nbthetaparams instead of numpopsizeparams
    double_gsamp_fcp = alt2d_alloc2Ddouble(genealogiessaved, nbthetaparams); // VS - nbthetaparams instead of numpopsizeparams
    if (npops > 2 || vs_test==1) // VS 10/11/2011 vs_test added
    {
      probgp_popsize = (double *) malloc(genealogiessaved * sizeof(double));
      probgp_migrate = (double *) malloc(genealogiessaved * sizeof(double));
    }

	// VS - go through each row of gsampinf
	// and through the FC summaries 
	// In the case the there are more than two populations it is adding all the qip and mip
	// i.e. it is computing the prior probability of the coalescent part of the genealogy
	// and computing the prior probability of the migration part of the genealogy
	for (gi = 0; gi < genealogiessaved; gi++)
    {
		// VS - need to initialize the probgp_popsize and probgp_migrate outside the for loop for groups
		// if(npops>=3) {
		if(npops>=3 || vs_test==1) { // VS vs_test added 10/11/2011 
			probgp_popsize[gi] = 0;
			probgp_migrate[gi] = 0;
		}
		// VS - for the full model we need to go through all the gsampinf matrices
		//for(gp=0; gp < nbgroupsloci_theta; gp++) {
			for (i = 0; i < nbthetaparams; i++) { // VS - Here it is multiplyig by 2 each FC value
				// VS this is because (2/theta)*fc=(2*fc)/theta. Here we are computing the 2*fc term.
				double_gsamp_fcp[gi][i] = 2.0 * gsampinf[g_gp_gs[i]][gi][fcp + g_paramindex[i]];  // VS gp - 
				// Note that the indexes here depend on the global variables
				// g_gp_gs[i] is the index of group in gsampinf matrix
				// g_paramindex[i] is the index of parameter
			}
			// VS if there are more than 2 populations, it will only focus in one of the parameters, either theta or migration.
			// if (npops > 2) 
			if (npops > 2 || vs_test==1) // VS 10/11/2011 vs_test
			// VS not possible to find the joint peak due to the large amount of parameters
			{
				//for (i = 0,probgp_popsize[gi] = 0; i < numpopsizeparams; i++) // VS - this would not work because the probgp would be initialized to zero for each group
				for (i = 0; i < nbthetaparams; i++)
					probgp_popsize[gi] += gsampinf[g_gp_gs[i]][gi][qip + g_paramindex[i]]; // VS gp - here it is adding all the Theta Integrals
					// Note that the indexes here depend on the global variables
					// g_gp_gs[i] is the index of group in gsampinf matrix
					// g_paramindex[i] is the index of parameter
				//for (i = 0,probgp_migrate[gi] = 0; i < nummigrateparams; i++) // VS - this would not work because the probgp would be initialized to zero for each group
				for (i = nbthetaparams; i < nparams; i++)
					probgp_migrate[gi] += gsampinf[g_gp_gs[i]][gi][mip + g_paramindex[i]]; // VS gp - here it is adding all the Migration Integrals
					// Note that the indexes here depend on the global variables
					// g_gp_gs[i] is the index of group in gsampinf matrix
					// g_paramindex[i] is the index of parameter
			}
		//} // end of for gp through groups of theta
    } // end of for gi to genealogiessaved
  } // end of IF init==0

  sum = 0;
  /* CR 110929.2 from JH 9/26/2011 fixed a bug  the initialization of 
   * log2diffx was being done even when no population size parameters 
   * are in the model being optimized (i.e. nowmodeltype==2) 
   */
  if (nowmodeltype == 0 || nowmodeltype == 1) {
	  //for (i=0;i<numpopsizeparams;i++)
	  //  log2diffx[i] = LOG2 - log (x[i]);
	  for (i=0;i<nbthetaparams;i++) // VS - in the full model we need to consider all the parameters in all groups
		  // VS - the x is the value at which we want to evaluate the function
		  // the user specifies one x for each parameter
		  log2diffx[i] = LOG2 - log (x[i]); // this is because the log of (2/theta)=log(2)-log(theta)
  }
  /* CR 110929.2 from JH 9/26/2011 fixed a bug  logx and divx were being 
   * set for parameters that are not in the model being optimized 
   * (i.e. nowmodeltype==2)
   */
  //for (i=nparamrange[0];i< nparamrange[1];i++)

  // VS NOTE: the aim is to find the maximum of the posterior f(theta|X)
  // Recall that if we have a sample G ~ f(G|X)
  // we can approximate f(theta|X)= SUM(f(G_i|theta)p(theta)/p(G_i))
  // Given that the prior of theta is uniform, we just need to focus on the ratio
  // f(G_i|theta)/p(G_i). The log of this ratio is 
  // log(f(G_i|theta))-log(p(G_i))
  // This is why in the loops below we have a minus sign for probg
  // and a plus sign for p(G_i|theta)
  //  for (i=0;i<nparams;i++)
  for (i=nparamrange[0];i< nparamrange[1];i++)
  {
    logx[i] = log(x[i]);
    divx[i] = 1.0/x[i];
  }
  for (gi = 0, ii=0; gi < genealogiessaved; gi++)
  {
    gp = 0; // VS - start with the gp=0
    g = gsampinf[gp][gi]; // VS gp - pointer to the gsampinf of group zero!!!
	// NOTE: g is pointing to gsampinf of group zero!!!
    dg = double_gsamp_fcp[gi];
    // if (npops == 2) // VS 10/11/2011 vs_test
	if (npops==2 && vs_test==0) // VS 10/11/2011 vs_test
      p =  - g[probgp]; // VS if there are 2 populations, we look at the P(G) - note that the probgp should be the same for all gsampinf matrices
    else
    {
      if (nowmodeltype == 1)
        p = -probgp_popsize[gi]; // VS if there are more than 2 pops, we look at the sum of QINT for theta
								 // see abobe Note about why this is a minus sign
      else
        p = -probgp_migrate[gi]; // VS if there are more than 2 pops, we look at the sum of MINT for theta
								 // see abobe Note about why this is a minus sign
    }
	
	// VS for loop through the groups of loci to compute P(G_i|Theta)
	// i.e. the remaining gsamping matrices.
	// Note that P(G_i|Theta)=PRODUCT(terms)
	// each term correspond to a parameter of a given group
	// So, for the nested and full model this loop is the same
	// ERROR
	// All the following is wrong because of the index - it is adding i to the gsampinf but it should be g_paramindex[i], at least for the theta parameters
	// for(gp=0; gp < nbgroupsloci_theta+nbgroupsloci_mig; gp++) {

		// g = gsampinf[gp][gi]; // VS points to the correct gsamping matrix, depending on the group

		// VS - these nparamrange are related with the number of parameters in the nested model??
		// it seems that this is where the information about the number of parameters is included.
		// Note that nparamrange will refer to index of all parameters put together
		// i.e. 0<i<(nbgroupsloci_theta*numpopsizeparams)+(nbgroupsloci_mig*nummigrateparams)
		for (i = nparamrange[0]; i < nparamrange[1]; i++) {

			g = gsampinf[g_gp_gs[i]][gi]; // VS points to the correct gsamping matrix, depending on the group

			// VS see abobe Note about why below there is a plus sign
			if (i< nbthetaparams) // VS nbgroupsloci_theta*numopsizeparams
				// VS this is the LOG of p(G|theta). Note that in this case theta is X
				// Log((2/theta)^cc (1/h) Exp(-2/theta fc)) = cc*(log(2)-log(theta))-log(h)-((2*fc)/theta)
				// Note that these log(2)-log(theta) were computed above as well as (2*fc)
				// Note that instead of i index we need to consider the global variable g_paramindex[i]
				p += (g[ccp + g_paramindex[i]]*log2diffx[i]) - g[hccp + g_paramindex[i]] - (dg[i] * divx[i]);
			else {
				// VS this is P(G|mig)
				// for the migration params, this is
				// log( m^mc Exp(-m*fm)) = mc*log(m) - m*fm
				// i1 = i - numpopsizeparams; // VS - note that in this case we use the g_paramindex[i]
				// 
				p += (g[mcp + g_paramindex[i]] * logx[i]) - (g[fmp + g_paramindex[i]] * x[i]);
			}
		} // end loop through params
	//} // end loop through groups
	

    // VS CHECK
	// print the values of x, the index of the genealogy and the corresponding f(G_i|theta) and f(G_i)
	/*printf("%i ", gi);
	for(i=0; i<nparams; i++) {
		printf("%4.5f ", x[i]);
	}
	printf(" %4.5f %4.5f\n", gsampinf[0][gi][probgp], p+gsampinf[0][gi][probgp]);
    */
    /* add the value to the linked list,  but only if it is within PRANGELOG of the biggest value found so far */ 
    if (gi == 0)
    {
      list[0]->v = p; // At this point p=log(f(G_i|theta))-log(p(G_i))
      list[0]->r = list[0]->l = NULL;
      last = first = now = list[0];
    }
    else
    {
      if (last->v - p < PRANGELOG) // only save values that are close to the biggest found so far 
      {
        ii++; // VS ii saves the number of values that have been added to the list
        listput(p, ii);
      }
    }
  } // end loop through genealogies

  // VS - At this point, the list contains the genealogies with the higher posterior probabilities
  // for a given combination of parameter values recorded in x.
  // I guess the reason is that if the posterior (or likelihood) is close to zero,
  // it wont contribute to the sum, and hence we can discard those genealogies.
  // the variable ii saves how many genealogies were recorded in the list
  // The value for the posterior is saved in each element of the list
  // in the double v. The list has the pointers first, last and now.
  // Hence, last->v is the posterior value of the last element in the list, and so on.

  iin = ii; // VS ii has the number of genealogies that have been included in the list with higher probabilities
  gi = 0;
  now = last; 
  // VS Why do we now point to the last element??? Not sure. 
  // VS Maybe it is because the pointer last always keeps the last element added, whereas now is last-1.
  // VS NOTE that now and last are both structures listelements

  do // VS loop through genealogies in the list
  {
	// VS eexp is computing EXP(now->v)
	// VS the result is saved in eexpsum[gi].m and eexpsum[gi].z in scientific notation
	// VS e.g. eexp(x, m, z)=> exp(x)=m*10^z??
    eexp (now->v, &eexpsum[gi].m, &eexpsum[gi].z);

    if (eexpsum[gi].z > maxz)
      maxz = eexpsum[gi].z;
    now = now->l; 
    gi++;
  }
  while (gi < iin && last->v - now->v < PRANGELOG); // VS - go through all genealogies of the list, i.e. 0<gi<iin 
  
  // VS - it seems that at this point we have the sum of the ratios
  // over all the genealogies which have the higher posterior probability.
  // Note that the ratio P(G|Theta)/P(G) in this case is both the posterior and the likelihood function!
  // see details in ima2/br/vs1/code/doc/LRT_Ima_VS.nb

  // VS
  // Go again through the genealogies saved that are in
  // the eexpsum variable
  // recall that the eexpsum variable contains the exp(log(p(G|q))-log(p(G)))
  // is this loop done to obtain the cummulative sum????
  // NOTE that we are working with the scientific notation
  // where x=log(y), and exp(x)=m*10^z
  // hence the sums are done by having a commun z and the we can add the m part
  // then the results is given back in log scale
  // Not sure why we go back and forth from log to exponential scale
  // I think it is because a sum of exponentials is not a sum of logs
  // if we had a product of expoentials we cold treat this as a sum of logs
  // but since we have a sum of exponentials we need to go back to the natural scale
  gin = gi;
  maxz -= OCUTOFF;
  for (gi = 0; gi < gin /*genealogiessaved */; gi++) 
  {
	  zadj = eexpsum[gi].z - maxz;
      if (zadj > -308 && zadj < 308)
	    eexpsum[gi].m *= pow10[POW10I(zadj)]; // VS multiply m=m*10^zadj. Note that after this m=m*10^zadj
      else
      {
        if (zadj <= -308)
          eexpsum[gi].m = 0.0;
        else
          eexpsum[gi].m = DBL_MAX;//trap overflow, shouldn't get here in this program
      }
	  acumm += eexpsum[gi].m; // VS acumm is the cumulative sum of the exponential
      if (calc_ess)
        acumm_sqr += eexpsum[gi].m * eexpsum[gi].m;
  } // end of second for loop to obtain the cummulative exponential sum
  if (calc_ess)
      *effective_n = acumm*acumm/acumm_sqr;

  sum = log (acumm) + maxz * LOG10; // VS why is this how we obtain the posterior. 
  // VS Why are we multiplying by LOG10 and maxz???
  q = log_numtrees - sum;  // negative of log of posterior probability
  // VS log_number of trees has the log of the number of trees.
  // VS this seems to be related with the formula of an average where we need to divide the sum
  // VS by the number of trees. Given that this is the negative of the posterior, and that we have logs
  // VS then it will the the log(nb_trees)-log(sum)
  // VS What I do not understand is what is the maxz and why we use it to multiply by the log(sum)???
    return (double) q;
} /* jointp*/ 

#undef OCUTOFF  

void freejointpmem()
{
  int i;
  alt2d_free2D(trialpop);
  free(bestvals);
  free(tempvals);
  free(parameter_lower_bound);
  free(parameter_upper_bound); 
  free(logx);
  free(divx);
  free(log2diffx); 
  alt2d_free2D(double_gsamp_fcp);
  alt2d_free2D(pop);
  if (npops > 2)
  {
    free(probgp_popsize);
    free(probgp_migrate);
  }
  free(popnest);
  for (i=0;i<genealogiessaved;i++)
    free(list[i]);
  free(list);
}


/***********GLOBAL FUNCTIONS **********/


// findjointpeaks
// finds the set of parameters for a given model with the maximum posterior probability 
// INPUT:
//	FILE *outfile : connection to stream to outputfile
//	char *outfilename : string with the name of the outputfile
//	char *nestfname : string with the name of the file with the nested models
//	int number_of_parameters : number of parameters in the nested model
/* CR 110921.1  Change type declaration of outfile parameter so that this
 * function can pass in the proper parameter type when calling closeopenout()
 */
void findjointpeaks(FILE **outfile,char *outfilename, char *nestfname,int number_of_parameters)
{
  int mi,j=0, typeloopstart, typeloopstop;
  double holdml;
  int holdparams;

  nparams = number_of_parameters;
  // VS Note that when we have groups of loci the number of parameters has to be multiplied by each group
  // nparams = (numpopsizeparams * nbgroupsloci_theta) + (nummigrateparams * nbgroupsloci_mig);

  bestvals = malloc((nparams+1) * sizeof(double)); // VS bestvals is a global variable that saves the parameter values corresponding to the maximum posterior probability
  setbounds(); // VS function that will set the lower and upper limit of the prior distributions of the parameters
  // NOTE: the LTR tests with groups of loci can only be done when all groups share the same prior, i.e. all groups for mig have the same prior
  // all groups for theta have the same prior
  setuplist(); // VS allocate memory for the list of values used in the optimization algorithm to find the maximum posterior probability
  
  depopsize = nparams * popsizemultiplier; // VS size of the "population" in the optimization algorithm (genetic algorithm)
  pop = alt2d_alloc2Ddouble(depopsize, nparams+1); // VS allocate memory for the "population" of parameter values used in the genetic algorithm to find the maximum posterior
  trialpop = alt2d_alloc2Ddouble(depopsize, nparams+1); // VS - allocate memory for a trial pop 
  
  popnest = (double *) malloc(nparams * sizeof(double)); // VS popnest - saves the values of the best parameters for the nested model
  if (strlen(nestfname)> 0) // VS if the nested file name has more than zero characters, read the file
    setup_mapping(nestfname); // VS this initializes the XMAP variable
  else // VS if the nested file name has less than or zero characters, the find the joint peaks for the full model
    num_nestedmodels = 0; 
  /* CR 110921.1  Change type of outfile parameter */
  startjointpeakouttable(*outfile, nestfname);
  // VS if (npops == 2)
  if (npops == 2 && vs_test==0) // VS added this vs_test variable to perform the search for the joint maximum just for the migration parameters
  {
    typeloopstart = 0;
    typeloopstop = 0;
    nparamrange[0] = 0;
    nparamrange[1] = nparams;
  }
  else
  {
    typeloopstart = 1;
    typeloopstop = 2;
  }
  for (nowmodeltype = typeloopstart; nowmodeltype <= typeloopstop; nowmodeltype++)
  {
    j++;
   // if (modeltypecounts[nowmodeltype] > 0)
    {
      printf(" starting model: %s\n", modelstartstr[nowmodeltype]);
      switch (nowmodeltype)
      {
        case 0 : 
              holdparams = nparams; 
              nparamrange[0] = 0;
              nparamrange[1] = nparams;
              break;
        case 1 : 
              holdparams = numpopsizeparams*nbgroupsloci_theta; // VS nbgroupsloci_theta
              nparamrange[0] = 0;
              nparamrange[1] = numpopsizeparams*nbgroupsloci_theta;
          break;
        case 2 : 
              holdparams = nummigrateparams*nbgroupsloci_mig; // VS nbgroupsloci_mig
              nparamrange[0] = numpopsizeparams*nbgroupsloci_theta; // VS nbgroupsloci_theta
              nparamrange[1] = (numpopsizeparams*nbgroupsloci_theta)+(nummigrateparams*nbgroupsloci_mig); // VS nbgroupsloci_theta nbgroupsloci_mig
          break;
      }
      modelloop(-1);
      holdml = bestvals[nparams];
      /* CR 110921.1  Change type of outfile parameter */
      fprintf(*outfile, "%d\t%s\t%d\t-\t-",j,
                         logpstrformat(-bestvals[nparams]), holdparams);
      fprintf(*outfile, "\t%s",logpstrformat(effective_sample_size));
      printjointpeakvals(*outfile,paramsallused, bestvals);
      closeopenout (outfile, outfilename);
      printf ("done model : %s\n", modelstartstr[nowmodeltype]);
      printf("joint density: %f\n",-bestvals[nparams]);
    }

	// VS - go through the nested models
    for (mi = 0;mi<num_nestedmodels;mi++) if (fullmodeltype[mi] == nowmodeltype)
    {
      j++;
      printf(" starting model: %s\n", modelnamelines[mi]);// to stdout
      modelloop(mi);
      mapvals(bestvals,popnest,holdxmaps[mi]);
      /* CR 110921.1  Change type of outfile parameter */
      fprintf(*outfile, "%d\t%s\t%d",j,logpstrformat(-bestvals[nparams]), 
                         nparams_nested[mi]);
      fprintf(*outfile,"\t%d",holdparams-nparams_nested[mi]);
      if (nestedmodel_boundary_in_model[mi])
        fprintf(*outfile,"*");
      fprintf(*outfile, "\t%s", logpstrformat(2*(bestvals[nparams]-holdml)));
      fprintf(*outfile, "\t%s",logpstrformat(effective_sample_size));
      printjointpeakvals(*outfile,paramsnotused[mi], popnest);
      closeopenout (outfile, outfilename);
      // print some info to stdout
      printf ("done model : %s\n", modelnamelines[mi]);
      printf("joint density: %f\n",-bestvals[nparams]);
    }
  }
  if (atleastone_boundary_in_model)
  /* CR 110921.1  Change type of outfile parameter */
    fprintf(*outfile,"    * test distribution of 2LLR is a mixture\n");
  fprintf(*outfile,"\n");
  freejointpmem();
}

