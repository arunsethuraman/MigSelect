/* MigSelect 2013-2019 Arun Sethuraman, Vitor Sousa, Jody Hey */
/* This code was developed based on IMa2  2009-2010  Jody Hey, Rasmus Nielsen and Sang Chul Choi*/

/* funccall.c   functions associated w/ the surface, and that make calls optimization routines */
#undef GLOBVARS
#include "imamp.h"

#define LOG_10_2  0.30102999566398119521
#define OCUTOFF  10
#define LOWVAL 1e-200
#define LOWLOG  -1e200

/*********** LOCAL STUFF **********/


/* function prototypes */
// VS changes these functions, adding the gp index as input
// static double marginp (int param, int firsttree, int lasttree, double x, int dummy);       //calculate marginal probability
// static double marginbis (double (*func) (double, double, int, int), double x1, double x2, double yadust, int pi);       // find a value associated w/ a certain marginal probability
//static double marginp (int param, int gp, int firsttree, int lasttree, double x, int dummy);       //calculate marginal probability
//static double marginp (int param, int paramindex, int gp, int firsttree, int lasttree, double x, int dummy);       //calculate marginal probability
static double marginp (int param, int firsttree, int lasttree, double x, int dummy);       //calculate marginal probability
static double marginbis (double (*func) (double, double, int, int), double x1, double x2, double yadust, int pi);       // find a value associated w/ a certain marginal probability

/****** LOCAL FUNCTIONS *********/

// marginp
// function that computes the marginal posterior density for a given parameter
// INPUT:
//	int param : index of the parameter in the gsampinf matrix
//	int paramindex : index of the parameter in the iparam structures
//	int gp : index of the group of loci for the structure iparam, e.g. gp for theta 0<gpTheta<nbgroupsloci_theta; gp for mig 0<gpMig<nbgroupsloci_mig.
//			NOTE: this is NOT the index for gsampinf structure
//	int firsttree : index of first tree considered
//	int lasttree : index of the last tree considered
//	double x : value at which the marginal distribution is evaluated
//	int dummy ??? is this a dummy variable that indicates something special?
#define OFFSCALEVAL 1
//double
//marginp (int param, int gp, int firsttree, int lasttree, double x, int dummy)
//double
//marginp (int param, int paramindex, int gp, int firsttree, int lasttree, double x, int dummy)
double
marginp (int param, int firsttree, int lasttree, double x, int dummy)
{
  int ei, p;
  double hval, sum, prob, temp, max, min, sumtemp = 0;
  double meani;
  // int paramindex; // created this variable to identify what is the index in the iparam 
  // structure that corresponds to the param, which is the index in the array with size
  // (numpopsizeparams*nbgroupsloci_theta)+(nummigrateparams*nbgroupsloci_mig)
  int gp; // VS gp 
  int nbparams, nbthetaparams, nbmigparams; // VS
  // VS compute some variables that will simplify the code in the for loops and allocation of variables
  nbparams=(numpopsizeparams*nbgroupsloci_theta)+(nummigrateparams*nbgroupsloci_mig); 
  nbthetaparams=numpopsizeparams*nbgroupsloci_theta;
  nbmigparams=nummigrateparams*nbgroupsloci_mig;

  dummy = 0; /* a dummy value is set to remove compiler complaints. */
  sum = 0;
  if (param < nbthetaparams) //
	  // tacking into account the fact that there are several groups of loci
  {
	  max = itheta[g_paramindex[param]].pr[g_gp_ip[param]].max; // VS gp
	  min = itheta[g_paramindex[param]].pr[g_gp_ip[param]].min; // VS gp
  }
  else//param < numpopsizeparams + nummigrateparams
  {
	  max = imig[g_paramindex[param]].pr[g_gp_ip[param]].max; // VS gp
      min = imig[g_paramindex[param]].pr[g_gp_ip[param]].min; // VS gp
	  if (modeloptions[EXPOMIGRATIONPRIOR])
		meani = 1.0 / imig[g_paramindex[param]].pr[g_gp_ip[param]].mean; // VS gp
  }

  if (x < min || x > max)
    return OFFSCALEVAL;
    //return 0;

  // VS before going into a for loop through the genealogies
  // decide what is the index for the group in gsampinf matrix
  // this depends on whether we are looking at the theta or migration parameters
  // VS NOTE about the index gp
  // the index gp given as input for marginp is valid for iparam structures
  // in gsampinf, for theta groups, gp varies between 0<gp<nbgroupsloci_theta-1
  // in gsampinf, for mig groups, gp varies between nbgroupsloci_theta<gp<(nbgroupsloci_theta+nbgroupsloci_mig)-1
  //if (param >= numpopsizeparams*nbgroupsloci_theta) { // if we are in migration params, need to correct gp
  //	gp = gp+nbgroupsloci_theta;
  //}

  // VS - here it is going through each genealogical tree
  // and it is computing P(G_i|theta)/Integral(P(G_i|theta),dtheta)
  // and summing it for all trees
  // Note that param is given as input and refers to the parameter
  // in the array of the best values
  // In this case, what we want to use are the values found above for the paramindex
  // because paramindex will refer to the corresponding index in the gsampinf matrix
  for (ei = firsttree; ei < lasttree; ei++)
  {
    // VS if (param < numpopsizeparams)
	if (param < nbthetaparams)
    {
      // VS p = param;
	  p = g_paramindex[param]; // VS
	  gp = g_gp_gs[param]; // VS

	  // VS 
	  // NOTE about access gsampinf and the meaning of gp
	  // The gp given as input for marginp refers to gp in iparam structure, which is not as in gsampinf
      hval = gsampinf[gp][ei][gsamp_hccp + p]; // VS gp - for theta parameters the index is correct
	  // VS hval cancells out in the marginal probability distribution
	  // see notes in ima2/br/vs1/code/doc/LRT_ImaVS.nb, section "Note about marginal posterior"
	  // However, given that the qint saved in gsampinf also includes hval, we need to consider it here.
      temp =
        -gsampinf[gp][ei][gsamp_qip + p] + gsampinf[gp][ei][gsamp_ccp + p]*(LOG2-log (x)) - // VS gp
        hval - 2 * gsampinf[gp][ei][gsamp_fcp + p] / x;
      sumtemp += exp (temp);
	  //			    p(G|q)*(1/hcc)
	  // temp= SUM(------------------------) 
	  //			INT[P(G|q)]dq*(1/hcc)
    }
    else //param < numpopsizeparams + nummigrateparams
    {
      // VS assert (param < numpopsizeparams + nummigrateparams);
	  assert (param < nbparams); // VS
      // VS p = param - numpopsizeparams;
	  p = g_paramindex[param]; // VS - paramindex is found in the for loop above
	  gp = g_gp_gs[param]; // VS
	  
      if (modeloptions[EXPOMIGRATIONPRIOR])
        temp = -gsampinf[gp][ei][gsamp_mip + p] + log (meani) - x * meani + // VS gp
          INTEGERROUND (gsampinf[gp][ei][gsamp_mcp + p]) * // VS gp
          log (x) - gsampinf[gp][ei][gsamp_fmp + p] * x; // VS gp
      else
        temp = -gsampinf[gp][ei][gsamp_mip + p] + // VS gp
          INTEGERROUND (gsampinf[gp][ei][gsamp_mcp + p]) * // VS gp
          log (x) - gsampinf[gp][ei][gsamp_fmp + p] * x; // VS gp - PROGRAM IS CRASHING HERE
      sumtemp += exp (temp);
    }
  }
  sumtemp /= (lasttree - firsttree + (firsttree == 0)); // sumtemp should be the same as sum if numbers in the range for exp() to work 
  prob = sumtemp;
  return -prob;                 /* negative because a minimization routine is used */
}                               /* marginp */
#undef OFFSCALEVAL

// marginbis
// returns the credible intervals
// INPUT
//	double (*func) (double, double, int, int): function to be evaluated (density) - this is usually margincalc function
//	double x1: ?? 
//	double x2: ??
//	double yadust: scale adjustment for y-axis 
//	int pi: index of parameter in parameter array (with all parameters for all groups putted together)
#define JMAX 40
#define BISTOL  1e-4
double
marginbis (double (*func) (double, double, int, int), double x1, double x2,
           double yadust, int pi)
{
  int j;
  double dx, f, fmid, xmid, rtb;
  f = (*func) (x1, yadust, pi, 1);
  fmid = (*func) (x2, yadust, pi, 1);
  if (f * fmid >= 0.0)
    return DBL_MIN;             /* not found does not appear to be a point corresponding to 95% limit */
  rtb = f < 0.0 ? (dx = x2 - x1, x1) : (dx = x1 - x2, x2);
  for (j = 1; j <= JMAX; j++)
  {
    fmid = (*func) (xmid = rtb + (dx *= 0.5), yadust, pi, 1);
    if (fmid <= 0.0)
      rtb = xmid;
    if (fabs (dx) < BISTOL || fmid == 0.0)
      return rtb;
  }
  return DBL_MAX;               /* Too many bisections in marginbis - cannot find root */
}


/********** GLOBAL FUNCTIONS ***********/

// margincalc
// margincalc does the same thing as marginp, but gets called by different functions
// unlike marginp, margincalc automatically uses all trees 
// can be used for direct calculation
// also can be used for root finding. 
// also by using a nonzero value of yadust can be used to find the value of x that is
// associated with a particular value of y (i.e. the likelihood) 
// if logi==1  return the logarithm 
// INPUT
//	double x: value for parameter value, for which we want to obtain the posterior
//	double yadust: if yadust is nonzero, returns the x value associated with an y value (useful to find the quantiles and credible intervals)
//	int pi: index of parameter in array where all parameters are put together for all groups
//	int logi: 0 - return in natural scale, 1 - return result in log scale
double
margincalc (double x, double yadust, int pi, int logi)
{
  int ei, p;
  double hval, sum, temp,meani;
  int gp; //VS, paramindex; // VS
  
  int nbparams, nbthetaparams, nbmigparams; // VS
  // VS compute some variables that will simplify the code in the for loops and allocation of variables
  nbparams=(numpopsizeparams*nbgroupsloci_theta)+(nummigrateparams*nbgroupsloci_mig); 
  nbthetaparams=numpopsizeparams*nbgroupsloci_theta;
  nbmigparams=nummigrateparams*nbgroupsloci_mig;

  sum = 0;

  if (modeloptions[EXPOMIGRATIONPRIOR] 
		&& (pi >= nbthetaparams) &&	 /* cr 110907.1 pi must be at least equal */ 
                                    /* to numpopsizeparams */
		   (pi < nbparams))  // VS changed the variable names to match the number of parameters when there are groups of loci
  {
	  
	  //meani = 1.0 / imig[pi - numpopsizeparams].pr[gp].mean; 
	  meani = 1.0 / imig[g_paramindex[pi]].pr[g_gp_ip[pi]].mean; // VS g_gp_ip[pi] global variable with index of group in iparam structures corresponding to pi, 
				// VS g_paramindex[pi] is a global variable with the index of parameter in structures iparam
  }

  // VS - decided to put the p and gp assignment outside the for loop.
  // This avoids performing this operation at each step of the for loop
  p = g_paramindex[pi]; // VS g_paramindex gives the index in structures iparam corresponding to pi
  gp = g_gp_gs[pi]; // VS g_gp_gs returns the index of group in the matrix gsampinf 

  // for loop through all genealogies
  for (ei = 0; ei < genealogiessaved; ei++)
  {
    if (pi < nbthetaparams) // VS nbthetaparams=numpopsizeparams*nbgroupsloci_theta
    { /*  for popsize params only  */
	  // p = pi; // VS
      hval = gsampinf[gp][ei][gsamp_hccp + p]; // VS gp
      temp =
        -gsampinf[gp][ei][gsamp_qip + p] + // VS gp
        INTEGERROUND (gsampinf[gp][ei][gsamp_ccp + p]) * (LOG2 - log (x)) - // VS gp
        hval - 2 * gsampinf[gp][ei][gsamp_fcp + p] / x; // VS gp
      sum += exp (temp);
    }
    else //pi < numpopsizeparams + nummigrateparams
    {  /*  for migration params only  */
      // p = pi - numpopsizeparams;
      if (modeloptions[EXPOMIGRATIONPRIOR])
        sum += exp (-gsampinf[gp][ei][gsamp_mip + p] + log (meani) - x * meani + // VS gp
                    INTEGERROUND (gsampinf[gp][ei][gsamp_mcp + p]) * // VS gp
                    log (x) - gsampinf[gp][ei][gsamp_fmp + p] * x); // VS gp
      else
        sum += exp (-gsampinf[gp][ei][gsamp_mip + p] + // VS gp
                    INTEGERROUND (gsampinf[gp][ei][gsamp_mcp + p]) * // VS gp
                    log (x) - gsampinf[gp][ei][gsamp_fmp + p] * x); // VS gp
    }
  }
  sum /= genealogiessaved;

  if (logi)
  {
    if (sum <= 0)
      sum = LOWLOG;
    else
      sum = log (sum);
  }
  sum -= yadust;
  return sum;
}                               /* marginalcalc */





// marginalopt
// function that will find the peaklocations for all the parameters in the model
// returns the mlval and peakloc arrays filled with the maximum likelihood and peak
// location for each marginal posterior of each parameter
// INPUT:
//  int firsttree   : index with the first tree
//	int lasttree    : index of the last tree
//  double *mlval   : returns the array with the maximum likelihood value, which also corresponds to the maximum likelihood value
//  double *peakloc : returns the array with the parameters that correspond to the peak location
// NOTE: in this function the index of the parameters is given considering that all
// the parameters for the different groups are given in an array
// theta1g1 theta2g1 ... theta1g2 theta2g2 ... m12g1 m21g1 ... m12g2 m21g2
#define SEARCHSTARTFRAC  4      // fraction of position in parameter range to start at
void
marginalopt (int firsttree, int lasttree, double *mlval, double *peakloc)
{
  int i;
  double ftol, ax, bx, cx, fa, fb, fc, xmax, ml;
  double axt, bxt, cxt, prior;
  double max0, min0, max1, min1;
  // VS needed to change this function double (*func) (int, int, int, double, int);
  double (*func) (int, int, int, double, int); // VS
//  int gp, paramindex; // VS gp refers to the group index, paramindex refer to the param index

  int nbparams, nbthetaparams, nbmigparams; // VS

  // VS compute some variables that will simplify the code in the for loops and allocation of variables
  nbparams=(numpopsizeparams*nbgroupsloci_theta)+(nummigrateparams*nbgroupsloci_mig); 
  nbthetaparams=numpopsizeparams*nbgroupsloci_theta;
  nbmigparams=nummigrateparams*nbgroupsloci_mig;
    //AS: debug only
    //printf("nbparams = %d, nbthetaparams = %d, nbmigparams = %d\n", nbparams, nbthetaparams, nbmigparams);

  ftol = 1e-7;
  func = marginp; // VS Note that this is the function that will be called by mnbrakmod
  for (i = 0;  i < nbparams; i++)
  {
	if (i < nbthetaparams) // tacking into account the fact that there are several groups of loci
	{
	  prior = itheta[g_paramindex[i]].pr[g_gp_ip[i]].max; // VS g_gp and g_paramindex
        //AS: debug only
          printf("In theta: i = %d, g_paramindex[i] = %d, g_gp_ip[i] = %d\n",i,g_paramindex[i],g_gp_ip[i]); 
	}
	else // param >= nbthetaparams
	{
	  prior = imig[g_paramindex[i]].pr[g_gp_ip[i]].max; // VS g_gp and g_paramindex
          //AS: debug only
          printf("In mig: i = %d, g_paramindex[i] = %d, g_gp_ip[i] = %d\n",i,g_paramindex[i],g_gp_ip[i]); 
	}
        
	ax = prior; 
	bx = prior / 2; 
	// VS need to change the i index given as input
	// this should refer to the parameter, so we can find a way to find what is the parameter to which
	// i refers to by 
	// VS mnbrakmod (i, firsttree, lasttree, &ax, &bx, &cx, &fa, &fb, &fc, func,0);
	mnbrakmod (i, firsttree, lasttree, &ax, &bx, &cx, &fa, &fb, &fc, func,0); 
	axt = ax;
	bxt = bx;
	cxt = cx;
	bx = prior / 2;
	ax = MINPARAMVAL;
	// VS mnbrakmod (i, firsttree, lasttree, &ax, &bx, &cx, &fa, &fb, &fc, func,0);
	mnbrakmod (i, firsttree, lasttree, &ax, &bx, &cx, &fa, &fb, &fc, func,0); 
	if (axt < bxt && axt < cxt)
	  min0 = axt;
	if (bxt < axt && bxt < cxt)
	  min0 = bxt;
	if (cxt < axt && cxt < bxt)
	  min0 = cxt;
	if (axt > bxt && axt > cxt)
	{
	  if (axt > prior) 
		axt = prior; 
	  max0 = axt;
	}
	if (bxt > axt && bxt > cxt)
	{
	  if (bxt > prior) 
		bxt = prior; 
	  max0 = bxt;
	}
	if (cxt > axt && cxt > bxt)
	{
	  if (cxt > prior) 
		cxt = prior; 
	  max0 = cxt;
	}
	if (ax < bx && ax < cx)
	  min1 = ax;
	if (bx < ax && bx < cx)
	  min1 = bx;
	if (cx < ax && cx < bx)
	  min1 = cx;
	if (ax > bx && ax > cx)
	  max1 = ax;
	{
	  if (ax > prior) 
		ax = prior; 
	  max1 = ax;
	}
	if (bx > ax && bx > cx)
	{
	  if (bx > prior)
		bx = prior; 
	  max1 = bx;
	}
	if (cx > ax && cx > bx)
	{
	  if (cx > prior) 
		cx = prior; 
	  max1 = cx;
	}
	if (max0 <= min1 || max1 <= min0)
	{
	  peakloc[i] = -1;
	}
	else
	{
	  //VS ml = -goldenmod (i, firsttree, lasttree, ax, bx, cx, ftol, &xmax, func,0);
	  ml = -goldenmod (i, firsttree, lasttree, ax, bx, cx, ftol, &xmax, func,0); // VS - in the end, since I created the global variables this function remained with the same input
	  mlval[i] = ml;
	  peakloc[i] = xmax;
	}
  } // end of for loop through the paramater values for each group
}                               /* marginalopt */

// margin95
// returns the marginal 95% credible intervals
// INPUT:
//	double mlval[]: array with the maximum likelihood (or posterior) for a set of parameters
//	double peakloc[]: array of size nparam with the set of parameters corresponding to the maximum likelihood value
//	int pi: index of parameter in peakloc
//	int UL: 0 - return lower credible interval
//			1 - return upper credible interval
#define  DOWN95  1.92
double
margin95 (double mlval[], double peakloc[], int pi, int UL)
{
  double x1, x2, x, yadjust;
  double (*func) (double, double, int, int);
  // int gp, paramindex; // VS

  func = margincalc;
  if (UL == 0)                  // lower
  {
    x1 = MINPARAMVAL;
    x2 = peakloc[pi];
  }
  else                          // upper
  {
    x1 = peakloc[pi];
    if (pi < numpopsizeparams*nbgroupsloci_theta) // VS nbgroupsloci_theta
    {
	  x2 = itheta[g_paramindex[pi]].pr[g_gp_ip[pi]].max; // VS g_gp_ip (global variable) is the index for the group for iparam structures given the pi index, 
		// VS g_paramindex (global variable) is the index for the parameter for iparam structures given the pi index
    }
    else //pi < numpopsizeparams + nummigrateparams
    {
	  x2 = imig[g_paramindex[pi]].pr[g_gp_ip[pi]].max; // VS g_gp_ip (global variable) is the index for the group for iparam structures given the pi index, 
		// VS g_paramindex (global variable) is the index for the parameter for iparam structures given the pi index
    }
  }
  yadjust = log (mlval[pi]) - DOWN95;
  x = marginbis (margincalc, x1, x2, yadjust, pi);
  return x;
}                               /* margin95 */

#define NUMTREEINT 2            // consider NUMTREEINT batches of trees for finding peaks,  one way to check for convergence

/* findmarginpeaks() 
	1) find marginal peaks for the main model  - save points in peakloc
	2) find 95% confidence limits on marginal peak locations 
		these calls ultimately go to the function marginp() which determines the marginal function value
    3)  does an LLR test on migration parameters using as a test distribution a distributino that is 50% 0 
    and 50% x^2_1df  This distribution has values
    2.70554  at p=0.05   The ratio of probabilities (as opposed to twice the log ratio) is 3.86813
    5.41189	  at p = 0.01  the ratio of prbabilities is 14.9685
    9.54954	 at p = 0.001  the ratio of probabilities is 118.483

*/

/* findmarginpeaks is a complex function that builds a table of parameter peak locations
estimated from the marginal posterior density
also does LLR tests */ 

void
findmarginpeaks (FILE * outfile, float *holdpeakloc)
{
  int gp; // paramindex; // VS gp, paramindex
  int gp_mig; // VS gp_mig - this is used for finding the marginal of 2NM, because we need to get the posterior for something which is a function of theta and mig params
  int gp_theta; // VS gp_theta - this is used for finding the marginal of 2NM, because we need to get the posterior for something which is a function of theta and mig params
  int i, j, k, ii, ilo, ihi, iihi, iilo, p;
  double temp, prior;
  double **mlval, **peakloc;
  double **popmigmlval, **popmigpeakloc;
  double *migtest, *popmigtest, *temptest;
  int *mpop, *mterm;
  char **popmigstr;
  int nmi, tempmpop, thetai, found, mi;
  int firsttree, lasttree;
  int printerrorfootnote = 0;
  int printsigfootnote = 0;
  int nummigprint;
  char sig[4][4] = { "ns\0", "*\0", "**\0", "***\0" };
  char llrstring[20];
  double maxp, max0p;

  
  int nbparams, nbthetaparams, nbmigparams; // VS

  // VS compute some variables that will simplify the code in the for loops and allocation of variables
  nbparams=(numpopsizeparams*nbgroupsloci_theta)+(nummigrateparams*nbgroupsloci_mig); 
  nbthetaparams=numpopsizeparams*nbgroupsloci_theta;
  nbmigparams=nummigrateparams*nbgroupsloci_mig;

    //AS: debug only
    printf ("nbparams = %d, nbthetaparams = %d, nbmigparams = %d\n", nbparams, nbthetaparams, nbmigparams);

  // VS - changed the number of parameters
  // p = numpopsizeparams+nummigrateparams;
  p = nbparams;
  mlval = orig2d_alloc2Ddouble (NUMTREEINT + 1, p); // this is an array with the density corresponding to the mode for each parameter?
  peakloc = orig2d_alloc2Ddouble (NUMTREEINT + 1, p); // this is an array with the mode of each parameter, obtained by looking at the marginal distribution
  if (outputoptions[POPMIGPARAMHIST])
  {
    // Popmig related variables are used to compute the scaled migration rates 2Nm
    // VS - modified the numpopsizeparams to (numpopsizeparams*nbgroupsloci_theta)
    // VS - modified the nummigrateparams to (nummigrateparams*nbgroupsloci_mig)
    popmigmlval = orig2d_alloc2Ddouble (NUMTREEINT + 1, nbmigparams); // VS
    popmigpeakloc = orig2d_alloc2Ddouble (NUMTREEINT + 1, nbmigparams); // VS
    popmigtest = malloc (nbmigparams* sizeof (double)); // VS
    mpop = malloc (nbmigparams * sizeof (int)); // VS mpop is an array of indexes for theta needed to compute the 2NM
    mterm = malloc (nbmigparams * sizeof (int)); // VS mterm is an array of indexes for the migration parameters to compute 2NM
    popmigstr = malloc (nbmigparams * sizeof (char *)); // VS
    for (i=0;i<nbmigparams;i++) // VS
      popmigstr[i] = malloc(PARAMSTRLEN *sizeof(char));
    nmi = 0;

    if (modeloptions[PARAMETERSBYPERIOD])
    {
      for (k = 0; k < lastperiodnumber; k++)
      {
        for (i = 0; i < npops - k; i++)
        {
          tempmpop = C[0]->plist[k][i]; //AS: leaving this as 0, since the plist should be same across chains Tue Feb 23 15:35:02 EST 2016
          thetai = 0;
          found = 0;
          while (!found && thetai < numpopsizeparams)
          {
            // finding the populations by period
            found = (k == atoi (&itheta[thetai].str[1]) && tempmpop == atoi (&itheta[thetai].str[3]));
            if (!found)
              thetai++;
          }
          assert (thetai < numpopsizeparams); 
          for (mi = 0; mi < nummigrateparams; mi++)
          {
            found = 0;
            if (modeloptions[SINGLEMIGRATIONBOTHDIRECTIONS])
            {
              found = (k == atoi (&imig[mi].str[1])
                       && (tempmpop == atoi (&imig[mi].str[3])
                           || tempmpop == atoi (&imig[mi].str[6])));
            }
            else
            {
              found = (k == atoi (&imig[mi].str[1])
                       && tempmpop == atoi (&imig[mi].str[3]));
            }
            if (found)
            {
              mpop[nmi] = tempmpop; // mpop is the index of theta parameter
              mterm[nmi] = mi; // mterm is the index of the migration parameter
              sprintf (popmigstr[nmi], "%d,2N%d", k, tempmpop);
              strcat (popmigstr[nmi], "M");
              if (modeloptions[SINGLEMIGRATIONBOTHDIRECTIONS])
              {
                strcat (popmigstr[nmi], &imig[mi].str[3]);
              }
              else
              {
                strcat (popmigstr[nmi], &imig[mi].str[1]);
              }
                //AS: debug only
                //printf("%d\t%c\n", nmi, popmigstr[nmi]);
              nmi++;
            }
          }
        }
      }
    } // end of IF parameters are coded by period
    else
    {
      for (i = 0; i < numtreepops - 1; i++)
      {
        thetai = i;
        for (mi = 0; mi < nummigrateparams; mi++) // VS - here it is wrong to go through all the groups of loci. All groups share the same strings for the parameter names.
        {
          if (modeloptions[SINGLEMIGRATIONBOTHDIRECTIONS])
          {
			// VS if there is only a single migration parameter looking for the migration parameters that involve thetai as the source or sink population
            found = ((thetai == atoi (&imig[mi].str[1])) || (thetai == atoi (&imig[mi].str[4])));
          }
          else
          {
            // Vs looking for the migration parameters that have as source population thetai population
            // AS: is this correct??
            found = thetai == atoi (&imig[mi].str[1]);
          }
          if (found)
          {
            mpop[nmi] = thetai; // VS mpop saves the index of source population
            mterm[nmi] = mi; // VS mterm saves the index of the migration parameter
            sprintf (popmigstr[nmi], "2N%d", thetai);
            strcat (popmigstr[nmi], "M");
            strcat (popmigstr[nmi], &imig[mi].str[1]);
            //AS: debug only
            printf ("%d\t%c\n", nmi, popmigstr[nmi]);
            nmi++;
          }
        }
      }
    } // end else of IF parameters are coded by period
  } // end of if(outputoptions[POPMIGPARAMHIST])


  migtest = malloc (nbmigparams * sizeof (double)); // VS nbmigparams
  FP "\nMarginal Peak Locations and Probabilities\n");
  FP "=========================================\n");
  FP "  peak locations are estimated using a peak finding algorithm, which may\n");
  FP "  fail if the curve has multiple peaks. All peaks can also be found, and\n"); 
  FP "  related analyses conducted, by plotting the histograms\n\n");
  printf ("Finding marginal peaks and probabilities \n");
  genealogiessaved = IMIN (MAXGENEALOGIESTOSAVE - 1, genealogiessaved);   // trap cases when too saving too many trees is attempted
  if (genealogiessaved <= 10)
  {
    FP " TOO FEW TREES SAVED - MARGINAL VALUES NOT FOUND \n");
    for (i = 0; i < p; i++)
      holdpeakloc[i] = -1;
  }
  else
  {

	// VS this is getting the peak for 
	// the first half of the genealogies saved
	// the second half of the genealogies saved
	// j here refers to the index of either the first or the second half of the trees
    for (firsttree = 0, lasttree = (int) genealogiessaved / NUMTREEINT, j = 0;
         j < NUMTREEINT; j++)
    {
      // VS marginalopt looks for the marginal distribution of each parameter
	  // it returns the maximum likelihood value in mlval for each parameter, 
	  // and also returns the parameter value that corresponds to the maximum posterior (or likelihood) in peakloc
      marginalopt (firsttree, lasttree, mlval[j], peakloc[j]); // mlval will be initialized and returned by this function
      if (outputoptions[POPMIGPARAMHIST])
        marginalopt_popmig (firsttree, lasttree, popmigmlval[j], popmigpeakloc[j],mpop, mterm);
      firsttree = lasttree + 1;
      lasttree += (int) genealogiessaved / NUMTREEINT;
      if (lasttree > genealogiessaved)
        lasttree = genealogiessaved;
	}
    marginalopt (0, genealogiessaved, mlval[NUMTREEINT], peakloc[NUMTREEINT]);
    
	// VS perform the LRT for the migration rate parameters
	for (i = nbthetaparams; i < nbparams; i++) // VS modified this loop to go through all groups of mig
	{

      maxp = -marginp(i, 0, genealogiessaved, peakloc[NUMTREEINT][i], 0); 
      max0p = -marginp(i, 0, genealogiessaved, MINPARAMVAL, 0); 
      migtest[i-nbthetaparams] = 2 * log(maxp/max0p); // VS nbthetaparams
    } // end for loop from i to total migrate parameters

	// Compute marginal for 2NM when option is on
    if (outputoptions[POPMIGPARAMHIST])
	{
		// Only compute the marginal for 2NM if the number of groups
		// is 1 for theta, 1 for mig
		if(nbgroupsloci_theta==1 && nbgroupsloci_mig==1) {
		
		}
		// is 1 for theta, 2 for mig
		if(nbgroupsloci_theta==1 && nbgroupsloci_mig==2) {
		
		}
		// is 2 for theta, 2 for mig
		if(nbgroupsloci_theta==2 && nbgroupsloci_mig==2) {
		
		}

	  // VS - TO BE DONE LATER - all these functions need to be checked carefully
      marginalopt_popmig (0, genealogiessaved, popmigmlval[NUMTREEINT], popmigpeakloc[NUMTREEINT], mpop, mterm);
     
        //AS: Please note that gp_mig, and gp_theta have not been set here at all. Tue Aug 18 12:13:12 EDT 2015

      // VS modified the calls to marginpopmig
      maxp = -marginpopmig (mterm[i], gp_mig, gp_theta, 0, genealogiessaved, popmigpeakloc[NUMTREEINT][i],mpop[i]); // VS added the group index for gp_theta and gp_mig
      max0p = -marginpopmig (mterm[i], gp_mig, gp_theta, 0, genealogiessaved, MINPARAMVAL,mpop[i]); // VS added the group index for gp_theta and gp_mig
      popmigtest[i] = 2 * log(maxp/max0p);

    } // end if outputoptions[POPMIGPARAMHIST]

	// VS Print the results of the marginal posterior distributions
    k = 0;
    for (; k <= numsplittimes; k++)
	{
      FP "\n");
      FP "Period %d\n", k);
      FP "--------\n");
      if (k < lastperiodnumber)       // set the upper bound on the ii loop,  no migration parameters if k==lastperiodnumber
	  {
        iilo = 1;
        if (outputoptions[POPMIGPARAMHIST]) // VS if printing the 2NM, need to go up to 3 in the ii for loop
          iihi = 3;
        else
          iihi = 2;
      }
      else // VS if it is the last period only print the thetaAncestral
	  {
        iilo = 1;
        iihi = 1;
      }
	  // VS go through a loop to print the different set of parameters
      for (ii = iilo; ii <= iihi; ii++)        // small loop ii==1 for population size parameters, ii==2 for migration parameters, ii=3 for 2NM migration parameters
	  {
        switch (ii)
        {
        case 1:                //population size
          FP "Population Size Parameters\n Param:");
		  for (i = 0; i < nbthetaparams; i++) { // VS nbthetaparams instead of numpopsizeparams
			  if (itheta[g_paramindex[i]].b == k) { // VS index i in iparam structure is given by global variable g_paramindex[i]
				// VS include group of loci
				FP "\t %s_g%i\tP", itheta[g_paramindex[i]].str, g_gp_ip[i]); // VS g_gp_ip[i] (global variable) with the index of group of loci for theta in iparam structure
			  }
		  }
          FP "\n");
          ilo = 0;
          ihi = nbthetaparams; // VS nbthetaparams=numpopsizeparam*nbgroupsloci_theta
          break;
        case 2:                // migration 
          nummigprint = 0;
          for (i = 0; i < nummigrateparams; i++) // VS checking if there are migration parameters that beggin in period k
            nummigprint += imig[i].b == k;
          if (nummigprint) // VS if there are migration rates in period k then print the migration rate posterior information
          {
            FP "Migration Rate Parameters\n Param:");
			for (i = nbthetaparams; i < nbparams; i++) { // VS modified the loop, to go from nbthetaparams to nbparams
			    // VS included the index of group of loci
				if (imig[g_paramindex[i]].pr[g_gp_ip[i]].max > MPRIORMIN) { // VS g_paramindex[i] g_gp_ip[i]
					if (imig[g_paramindex[i]].b == k) { // VS
						FP "\t %s_g%i\tP", imig[g_paramindex[i]].str, g_gp_ip[i]); // VS g_paramindex[i] g_gp_ip[i] 
					}
				}
			}
            FP "\n");
		  }
          ilo = nbthetaparams; // VS nbthetaparams=numpopsizeparams*nbgroupsloci_theta
          ihi = nbparams; // VS nbparams=(nbgroupsloci_theta*numpopsizeparams) + (nbgroupsloci_mig*nummigrateparams)
          temptest = migtest;
          break;
        case 3:                // 2NM 
          nummigprint = 0;
          for (i = 0; i < nummigrateparams; i++) 
            nummigprint += imig[i].b == k;
          if (nummigprint)
		  {
            FP "Population Migration (2NM) Terms\n Term:");
			for (i = 0; i < nummigrateparams; i++) {
			    // VS need to include here a for loop for the groups of loci
				for(gp=0; gp<nbgroupsloci_mig; gp++) {
					if (imig[mterm[i]].pr[gp].max > MPRIORMIN) { // VS gp
						if (imig[mterm[i]].b == k) {
							FP "\t %s_g%i\tP", popmigstr[i], gp);
						}
					}
				}
			}
            FP "\n");
          }
          ilo = nbthetaparams; // VS nbthetaparams=numpopsizeparams*nbgroupsloci_theta
          ihi = nbparams; // VS nbparams=(nbgroupsloci_theta*numpopsizeparams) + (nbgroupsloci_mig*nummigrateparams)
          temptest = popmigtest;
          break;
		}

		// VS the above switch of ii defines the lower (ilo) and upper (ihi) 
		// index to look for the parameters depending on wheter we are dealing
		// with theta parameters, migration m or 2NM
		for (j = 0; j <= NUMTREEINT; j++) {
          if (ii == 0 || (ii == 2 && nummigprint > 0) || (ii == 3 && nummigprint > 0) || ii == 1)
		  {
			if (j < NUMTREEINT) {
				  FP " Set%d", j);
			}
			else {
              FP " All");
			}
			// VS for loop going through all the parameters between ilo and ihi (these depend on the switch ii, i.e. depend on which parameter we are looking at)
			for (i = ilo; i < ihi; i++) {
              if ((ii == 0)
                  || (ii == 1 && itheta[g_paramindex[i]].b == k) // VS g_paramindex[i] instead of i-ilo 
                  || (ii == 2 && imig[g_paramindex[i]].b == k) // VS g_paramindex[i] instead of i-ilo 
                  || (ii == 3 && imig[mterm[g_paramindex[i]]].b == k) // VS g_paramindex[i] instead of i-ilo 
                   )

              {
                if (ii < 3)
				{
                  if (peakloc[j][i] >= 0)
				  {
                    FP "\t%7.3lf\t%7.3lf", peakloc[j][i], mlval[j][i]); // VS - printing the peak location for set of genealogies j and param index i
																		// VS - need to related this i with the groups of loci, and the respective parameters
                  }
                  else
				  {
                    FP "\terror*\t");
                    printerrorfootnote = 1;
                  }
                }
                else
				{
                  if (popmigpeakloc[j][g_paramindex[i]] >= 0) // VS g_paramindex[i] instead of i-ilo 
				  {
                    FP "\t%7.3lf\t%7.3lf", popmigpeakloc[j][i-ilo], popmigmlval[j][i-ilo]); // VS printing the peak values, for index i-ilo. Need to check these and its relation with the param and group index
													// VS g_paramindex[i] instead of i-ilo??? NOT SURE HERE - TO BE DONE LATER 
                  }
                  else
				  {
                    FP "\terror*\t");
                    printerrorfootnote = 1;
                  }
                }
			  } // VS end if 
			} // VS end for loop through i from ilo to ihi (i.e. through the parameter index)
            FP "\n");
		  } // VS end if (ii == 0 || (ii == 2 && nummigprint > 0) || ...
		} // VS end for through j, which refers to the different sets of genealogies used
        
		// VS the 95% credible intervals are only computed for the set with "ALL" genealogies
		if (ii == 0 || ii == 1 || ((ii == 2|| ii==3)&& nummigprint > 0))
		{
          if (ii < 3)
          {
            FP " LR95%%Lo");
			// VS for loop through the parameters
            for (i = ilo; i < ihi; i++)
            {
              if ((ii == 0)
                  || (ii == 1 && itheta[g_paramindex[i]].b == k) // VS g_paramindex[i] instead of i-ilo
                  || (ii == 2 && imig[g_paramindex[i]].b == k)) // VS g_paramindex[i] instead of i-ilo

              {
                if (peakloc[NUMTREEINT][i] >= 0)
                {
                  temp =
                    margin95 (mlval[NUMTREEINT], peakloc[NUMTREEINT], i, 0);
                  if (temp <= 0)
                  {
                    FP "\t<min\t");
                  }
                  else
                  {
                    if (temp >= DBL_MAX || temp <= DBL_MIN)
                      FP "\tna\t");
                    else
                      FP "\t%7.3lf\t", temp);
                  }
                }
                else
                {
                  FP "\t\t");
                }
              }
			}
            FP "\n");
            FP " LR95%%Hi");
			
            for (i = ilo; i < ihi; i++) // VS for loop through the parameters
            {
              if ((ii == 0)
                  || (ii == 1 && itheta[g_paramindex[i]].b == k) // VS g_paramindex[i] (global variable) with the index of the parameter in iparam structures
                  || (ii == 2 && imig[g_paramindex[i]].b == k)) // VS g_paramindex[i] (global variable) with the index of the parameter in iparam structures

              {
                if (ii == 1)
                  prior = itheta[g_paramindex[i]].pr[g_gp_ip[i]].max; // VS g_gp_ip[i] (global variable) is the index of groups in iparam structure corresponding to parameter i
																	  // VS g_paramindex[i] (global variable) with the index of the parameter in iparam structures
                else
                  prior = imig[g_paramindex[i]].pr[g_gp_ip[i]].max; // VS g_gp_ip[i] (global variable) is the index of groups in iparam structure corresponding to parameter i
																	// VS g_paramindex[i] (global variable) with the index of the parameter in iparam structures
                if (peakloc[NUMTREEINT][i] >= 0)
                {
                  temp = margin95 (mlval[NUMTREEINT], peakloc[NUMTREEINT], i, 1);
                  if (temp >= prior || temp <= DBL_MIN)
                  {
                    FP "\t>max\t");
                  }
                  else
                  {
                    if (temp <= DBL_MIN)
                      FP "\tna\t");
                    else
                      FP "\t%7.3lf\t", temp);
                  }
                }
                else
                {
                  FP "\t\t");
                }
			  }
			} // VS end for i from ilo to ihi
            FP "\n");
		  } // VS end if (ii < 3)
          if (ii == 2 || ii == 3)
		  {
            FP " LLRtest ");

			// for loop through the migration parameters for which a LLR tests is done for the migration rates
			for (i = nbthetaparams; i < nbparams; i++) { // VS nbparams=(numpopsizeparams*nbgroupsloci_theta)+(nummigrateparams*nbgroupsloci_mig); 
			  
			  if ((ii==2 && imig[g_paramindex[i]].b == k && imig[g_paramindex[i]].pr[g_gp_ip[i]].max > MPRIORMIN ) || // VS g_paramindex[i], g_gp_ip[i]
                  (ii==3 && imig[mterm[g_paramindex[i]]].b == k && imig[mterm[g_paramindex[i]]].pr[g_gp_ip[i]].max > MPRIORMIN )) // VS g_paramindex[i], g_gp_ip[i]
              {
                if (fabs (temptest[i-nbthetaparams]) > 1e6) // VS - need to decrease nbthetaparams to have the correct index
                  sprintf (llrstring, "bad value\t");
                else
                {
				  printsigfootnote = 1; // VS moved this section here so that the footnote is always printed when there are migration parameters
                  if (temptest[i-nbthetaparams] > 9.54954) // VS - need to decrease nbthetaparams to have the correct index
                  {
                    //printsigfootnote = 1; // VS why isn't the footnote always printed when there are migration events?
				    // VS moved this one level up, after the else statement. 
                    sprintf (llrstring, "%7.3lf%s\t", temptest[i-nbthetaparams], sig[3]); // VS - need to decrease nbthetaparams to have the correct index
                  }
                  else if (temptest[i-nbthetaparams] > 5.41189) // VS - need to decrease nbthetaparams to have the correct index
                    sprintf (llrstring, "%7.3lf%s\t", temptest[i-nbthetaparams], sig[2]); // VS - need to decrease nbthetaparams to have the correct index
                  else if (temptest[i-nbthetaparams] > 2.70554) // VS - need to decrease nbthetaparams to have the correct index
                    sprintf (llrstring, "%7.3lf%s\t", temptest[i-nbthetaparams], sig[1]); // VS - need to decrease nbthetaparams to have the correct index
                  else
                    sprintf (&llrstring[0], "%7.3lf%s\t", temptest[i-nbthetaparams], sig[0]); // VS - need to decrease nbthetaparams to have the correct index
				}
                FP "%s", llrstring);
			  }
			} // end of for loop through parameters
            FP "\n");
		  } // end if (ii == 2 || ii == 3)
          if (ii == 1)
          {
            FP " LastPeriod");
            for (i = 0; i < nbthetaparams; i++) // VS nbthetaparams instead of numpopsizeparams
              if (itheta[g_paramindex[i]].b == k) // VS g_paramindex[i] instead of i
                FP "\t%d\t", itheta[g_paramindex[i]].e); // VS gparamindex[i] instead of i
            FP "\n");
          }
          if (ii == 2 || ii==3)
          {
            FP " LastPeriod");
            for (i = nbthetaparams; i < nbparams; i++) // VS go through nbgroupsloci_mig*nummigrateparams, which is the same as going from nbthetaparams to nbparams
			{
			  if (ii==2 && imig[g_paramindex[i]].b == k && imig[g_paramindex[i]].pr[g_gp_ip[i]].max > MPRIORMIN ) // VS g_paramindex[i], g_gp_ip[i]
                FP "\t%d\t", imig[g_paramindex[i]].e);
              if (ii==3 && imig[mterm[g_paramindex[i]]].b == k && imig[mterm[g_paramindex[i]]].pr[g_gp_ip[i]].max > MPRIORMIN ) // VS gp, paramindex
                FP "\t%d\t", imig[mterm[g_paramindex[i]]].e);
            }
            FP "\n");
          }
		} // end if  if (ii == 0 || ii == 1 || ((ii == 2|| ii==3)&& nummigprint > 0))
	  } // end for loop through ii
	} // end for loop through periods
    for (i = 0; i < p; i++)
      holdpeakloc[i] = (float)
        (peakloc[NUMTREEINT][i] < 0) ? (double) MINPARAMVAL : peakloc[NUMTREEINT][i];
  } // end of if else genealogiessaved <= 10

  if (printerrorfootnote==1)
    FP "*  peak not found possibly due to multiple peaks (check plot of marginal density) \n");
  if (printsigfootnote == 1)
  {
     FP " migration rate likelihood ratio test - see Nielsen and Wakeley (2001)\n migration significance levels :  * p < 0.05;   **  p < 0.01,   *** p < 0.001\n");
  }
  
  FP "\n");


  orig2d_free2D ((void **) mlval, NUMTREEINT + 1);
  orig2d_free2D ((void **) peakloc, NUMTREEINT + 1);
  XFREE (migtest);
  if (outputoptions[POPMIGPARAMHIST])
  {
    orig2d_free2D ((void **) popmigmlval, NUMTREEINT + 1);
    orig2d_free2D ((void **) popmigpeakloc, NUMTREEINT + 1);    
    XFREE(popmigtest);
    XFREE(mpop);
    XFREE(mterm);
    for (i=0;i<nummigrateparams;i++)
      XFREE(popmigstr[i]);
    XFREE(popmigstr);
  }
}                               /* findmarginpeaks */


