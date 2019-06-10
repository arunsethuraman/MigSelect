/* MigSelect 2013-2019 Arun Sethuraman, Vitor Sousa, Jody Hey */
/* This code was developed based on IMa2  2009-2010  Jody Hey, Rasmus Nielsen and Sang Chul Choi*/


#undef GLOBVARS
#include "imamp.h"

/*********** LOCAL STUFF **********/


static double calcx (int ei, int pnum, int mode);

/******* LOCAL FUNCTIONS ***********/

// calcx
// computes the mean or variance of demographic parameters based on the genealogies sampled
// for a given genealogy
// INPUT
//	int ei: genealogy index in gsampinf 
//	int pnum: parameter index (note that this index considers all parameters put together in an array, i.e. 0<pnum<(nbgroupsloci_theta*numpopsizeparams)+(nbgroupsloci_mig*nummigrateparams))
//	int mode: 0 for means, 1 for variances
double
calcx (int ei, int pnum, int mode)
{
  int cc, mc,p;
  double fc, fm, hval, denom, max;
  double tempval;

	int gp; // VS gp
	int nbparams, nbthetaparams, nbmigparams; // VS

	// VS compute some variables that will simplify the code in the for loops and allocation of variables
	nbparams=(numpopsizeparams*nbgroupsloci_theta)+(nummigrateparams*nbgroupsloci_mig); 
	nbthetaparams=numpopsizeparams*nbgroupsloci_theta;
	nbmigparams=nummigrateparams*nbgroupsloci_mig;

  //p = pnum; // VS using p to refer to the index of parameter
  if (pnum < nbthetaparams) // VS nbthetaparams instead of numpopsizeparams
  {
	gp = g_gp_gs[pnum]; // VS - gp is the group index in gsampinf
	p = g_paramindex[pnum]; // VS p is the param index in iparam and gsampinf matrix

    if (itheta[p].pr[g_gp_ip[pnum]].max <= MINPARAMVAL) // VS g_gp_ip[pnum] - Note that here the gp is not the same as gp in gsampinf matrix.
      return -1;
    cc = (int) gsampinf[gp][ei][gsamp_ccp + p]; // VS gp
    fc = gsampinf[gp][ei][gsamp_fcp + p]; // VS gp
    assert(fc >= 0);
    hval = gsampinf[gp][ei][gsamp_hccp + p]; // VS gp
    denom = gsampinf[gp][ei][gsamp_qip + p]; // VS gp
    max = itheta[p].pr[g_gp_ip[pnum]].max; // VS g_gp_ip[pnum]. Note that here the gp is not the same as gp in gsampinf matrix.
    if (mode == 0) // VS computing the mean
    {
      if (cc == 0 && fc == 0)
      {
        tempval = (SQR (max) / 2) / exp (denom);
      }
      else
      {
        if (cc > 1)
        {
          tempval = exp (2 * LOG2 - hval + (2 - cc) * log (fc) + uppergamma (cc - 2, 2 * fc / max) - denom);
        }
        else
        {
          if (cc == 1)
          {
            tempval = exp (LOG2 - hval + log (max * exp (-2 * fc / max) - 2 * fc * exp (uppergamma (0, 2 * fc / max))) - denom);
          }
          else                  //cc==0
          {
            assert (cc == 0);
            tempval = exp (log ((max / 2) * (max - 2 * fc) * exp (-2 * fc / max) + 2 * SQR (fc) * exp (uppergamma (0, 2 * fc / max))) - denom);
          }
        }
      }
    }
    else                        //mode == 1 // VS computing the variance
    {
      if (cc == 0 && fc == 0)
        tempval = (max * SQR (max) / 3) / exp (denom);
      else
      {
        if (cc > 2)
        {
          tempval = exp (uppergamma (cc - 3, 2 * fc / max) + 3 * LOG2 - hval + (3 - cc) * log (fc) - denom);
        }
        else
        {
          if (cc == 2)
          {
            tempval = exp (2 * LOG2 - hval + log (max * exp (-2 * fc / max) - 2 * fc * exp (uppergamma (0, 2 * fc / max))) - denom);
          }
          else
          {
            if (cc == 1)
            {
              tempval = exp (-hval + log (max * (max - 2 * fc) * exp (-2 * fc / max) + 4 * SQR (fc) * exp (uppergamma (0, 2 * fc / max))) - denom);
            }
            else                //cc==0
            {
              assert (cc == 0);
              tempval = exp (-log (3.0) + log (max * (2 * SQR (fc) - fc * max + SQR (max)) * exp (-2 * fc / max) - 4 * pow ((double) fc, 3.0) * exp (uppergamma (0, 2 * fc / max))) - denom);
            }
          }
        }
      }
    }
  }
  else // VS if pnum refers to migration parameters
  {
    // VS p -= numpopsizeparams;
    gp = g_gp_gs[pnum]; // VS - gp is the group index in gsampinf
	p = g_paramindex[pnum]; // VS p is the param index in iparam and gsampinf matrix

    if (imig[p].pr[g_gp_ip[pnum]].max <= MINPARAMVAL) // VS g_gp_ip[pnum]
      return -1;
    mc = (int) gsampinf[gp][ei][gsamp_mcp + p]; // VS gp
    fm = gsampinf[gp][ei][gsamp_fmp + p]; // VS gp
    assert(fm >= 0);
    if (fm < 0.0)
      {
        IM_errloc (AT, "fm cannot be negative.");
      }
    denom = gsampinf[gp][ei][gsamp_mip + p]; // VS gp
    max = imig[p].pr[g_gp_ip[pnum]].max; // VS g_gp_ip[pnum]
    if (mode == 0) // VS means
    {
      if (mc == 0 && fm == 0)
      {
        tempval = (SQR (max) / 2) / exp (denom);
      }
      else
      {
        if (mc > 0)
          tempval = exp (lowergamma ((int) mc + 2, fm * max) - (mc + 2) * log (fm) - denom);
        else
          tempval = (1 - (1 + fm * max) * exp (-fm * max)) / SQR (fm) / exp (denom);
      }
    }
    else                      //mode == 1 // VS variances
    {
      if (mc == 0 && fm == 0)
      {
        tempval = (pow (max, 3.0) / 3) / exp (denom);
      }
      else
      {
        tempval = exp (lowergamma ((int) mc + 3, fm * max) - (mc + 3) * log (fm) - denom);
      }
    }
    
  }
  assert (tempval >= 0);         // it should be greater than or equal to 0   (not sure really,  is it ok to be 0 ? 
  return tempval;
}                               //calcx



/********** GLOBAL functions ********/

void closeopenout (FILE ** p_to_file, char fname[]) // why a pointer to a pointer? so the address can be passed back
{
  f_close (* p_to_file);
  if ((* p_to_file = fopen (fname, "a+")) == NULL)
  {
    IM_err(IMERR_OUTPUTFILECHECK,"Error opening text file for writing");
  }
}


/*
 * This function is used for debugging only.  CR 110825.1
 */
void
checkoutfileclosed (FILE ** outfile, char outfilename[])
{
  if ((*outfile = fopen (outfilename, "w")) == NULL)
  {
    IM_err(IMERR_OUTPUTFILECHECK,"Error  - file for results output is not closed ");
  }
  f_close (*outfile);
}                               // checkoutfileclosed



// printrunbasics
// prints the basic statistics of the run
// INPUT:
//	FILE * outfile:
//	int loadrun:
//	char fpstr[]:
//	int burnsteps:
//	int recordint: 
//	int recordstep:
//	int savegenealogyint:
//	time_t endtime:
//	time_t starttime:
//	double hilike:
//	double hiprob:
void
printrunbasics (FILE * outfile, int loadrun, char fpstr[], int burnsteps,
                int recordint, int recordstep, int savegenealogyint,
                time_t endtime, time_t starttime, double hilike,
                double hiprob/*, int step*/,int z)
{
  double seconds;
  int li;

  FP "%s", fpstr);

  if (!loadrun == 1)
  {
    FP "\n\nMCMC INFORMATION\n");
    FP "===========================\n\n");
    FP "Number of steps in burnin: %10d\n", burnsteps);
    FP "Number of steps in chain following burnin: %10d \n",
      (step - burnsteps));
    FP "Number of steps between recording : %d  Number of record steps: %d \n", recordint, recordstep);
    FP "Number of steps between saving genealogy information: %d  Number of genealogies saved per locus: %d \n", savegenealogyint, genealogiessaved);
    time (&endtime);
    seconds = difftime (endtime, starttime);
    FP "\nTime Elapsed : %d hours, %d minutes, %d seconds \n\n",
      (int) seconds / (int) 3600,
      ((int) seconds / (int) 60) -
      ((int) 60 * ((int) seconds / (int) 3600)),
      (int) seconds - (int) 60 *((int) seconds / (int) 60));
    FP "Highest Sampled Joint P(G) (log) : %10.3f \n", hiprob);
    FP "Highest Joint P(D|G) (log) : %10.3f \n\n", hilike);
    FP "Highest P(D|G) (log) for each Locus \n");
    FP "\tLocus\tP(D|G)\n");
    for (li = 0; li < nloci; li++)
    {
      FP "\t%d\t%.3f\n", li, C[z]->G[li].hilike);
    }
    FP "\n");
  }
}                               /* printrunbasics */

void checkhighs (int ci, int printint, double *hilike, double *hiprob,
                 double *like/*, int step*/) 
{
  int li;
  double temp;

#ifdef _DEBUG
  int cmmmax, cmmli, cmmi, i;

#endif /*  */

/* fill hilocuslike, hilike and hiprob  */
  for (*like = 0, li = 0; li < nloci; li++)
  {
    temp = C[ci]->G[li].pdg;
    *like += temp;
    if (C[ci]->G[li].hilike < temp)
      C[ci]->G[li].hilike = temp;
  }

//      assert(*like == C[ci]->allpcalc.pdg); should be close
  if (*hilike < C[ci]->allpcalc.pdg)
    *hilike = C[ci]->allpcalc.pdg;
  if (*hiprob < C[ci]->allpcalc.probg)
    *hiprob = C[ci]->allpcalc.probg;

#ifdef _DEBUG
  if ((step / (int) printint) * (int) printint == step && step > 0)
  {
    for (cmmmax = 0, li = 0; li < nloci; li++)
    {
      for (i = 0; i < L[li].numlines; i++)
      {
        if (cmmmax < C[0]->G[li].gtree[i].cmm)
        {
          cmmmax = C[0]->G[li].gtree[i].cmm;
          cmmli = li;
          cmmi = i;
        }
      }
    }
    printf ("\n max cmm %d  genealogy %d  edge %d \n", cmmmax, cmmli, cmmi);
  }

#endif /*  */
}                               /* checkhighs */


#define MAXLENX  200
#define MAXLENY  104
#define DEFAULTLENX 150;
#define DEFAULTLENY  40;

/* 4/27/09  JH modified so that the full x axis is used even if there are < TRENDDIM points in the array */ 

void asciitrend (FILE * outfile, struct value_record *v, int trenddoublepoint,int trendspot)
/* structure of graph 
0 name
1 ytop    |
40 ybot   |
41         -------
42        xbot    xtop

trenddoublepoint is the position in the array below which the points have twice the density 
with respect to the points abouve the trenddoublepoint 

To plot the array so that the plotted positions are proportional to time,  we adjust the plotting positions as a function 
of the density difference and the position of trenddoublepoint. 
If a is where trenddoublepoint would be plotted without adjustment and x is the value by which to multiply positions from 0 up to a,  to get 
their adjusted plotposition  then 

Solve[x a + (x/2) (d - a) == d, x]

This is becuas we know that the total distance must be d (i.e. available length of x axis)
then {x -> (2 d)/(a + d)}}

let f(a,d) = (2 d)/(a + d)
Then for a point at position c in the array and that is lower than position a,  the plot position is c*f(a, d)
position a itself is plotted at a * f(a,d)
For a point at position c in the array that is higher than position a,  the plot position is 
 a * f(a,d) + (c - a)*f(a, d)/2
Thus we rescale points up to point trendoublepoint by a factor of x 
and we rescale points above that by a factor of x/2 
*/
{
  char graph[MAXLENY][MAXLENX];
  char tc[20];
  double yt, ymax = -1e200, ymin = 1e200;
  double rescalerbelow, rescalerabove;
  double *y;
  double logscale;
  int i, xspot, xtemp, xmax, yspot, notplot = 0, rescaledoublepoint;
  int xlen, ylen, dim, localdoublepoint;
  y = v->trend;
  if (v->do_logplot)
    logscale = UMAX;
  xlen = DEFAULTLENX;
  ylen = DEFAULTLENY;
  if (trendspot < TRENDDIM)
  {
    dim = trendspot;
    localdoublepoint = trendspot;
  }
  else
  {
    dim = TRENDDIM;
    localdoublepoint = trenddoublepoint;
  }
  for (i = 0; i < dim; i++)
  {
    yt = y[i];
    if (ymax < yt)
      ymax = yt;
    if (ymin > yt)
      ymin = yt;
  }
  if (v->do_logplot)
  {
    if (ymax > (double) logscale)
      ymax = (double) logscale;
    if (ymin < 1 / (double) logscale)
      ymin = 1 / (double) logscale;
  }
  strcpy (graph[0], v->str);
  strcat (graph[0], " trend");
  if (v->do_logplot)
    strcat (graph[0], " - log scale");
  sprintf (tc, "%8.4g", ymax);
  strcpy (graph[1], tc);
  while ((int) strlen (graph[1]) < xlen - 1)
    strcat (graph[1], " ");
  sprintf (tc, "%8.4f", ymin);
  strcpy (graph[ylen - 2], tc);
  while ((int) strlen (graph[ylen - 2]) < xlen - 1)
    strcat (graph[ylen - 2], " ");
  strcpy (graph[ylen - 1], "          ");
  while ((int) strlen (graph[ylen - 1]) < xlen - 1)
    strcat (graph[ylen - 1], "-");
  for (i = 2; i < ylen - 2; i++)
  {
    graph[i][0] = '\0';
    while ((int) strlen (graph[i]) < xlen - 1)
      strcat (graph[i], " ");
  } 
  
  for (i = 1; i < ylen - 1; i++)
    graph[i][10] = '|';

  xmax = (int) (xlen - 12);
  rescaledoublepoint =  INTEGERROUND (((float) (xmax * localdoublepoint) / dim));
  rescalerbelow = (2 * (float) xmax) / (float) (rescaledoublepoint + xmax);
  rescalerabove = rescalerbelow / 2;

  for (i = 0; i < dim; i++)
  {
    if (y[i] >= ymin && y[i] <= ymax)
    {
      if (v->do_logplot)
        yspot = (int) 1 + (ylen - 2) - (int) ((ylen - 2) * (log (y[i]) - log (ymin)) / (log (ymax) - log (ymin)));

      else
        yspot = (int) 1 + (ylen - 2) - (int) ((ylen - 2) * (y[i] - ymin) / (ymax - ymin));
    }
    else
    {
      yspot = -1;
      notplot = 1;
    }
    xtemp = INTEGERROUND (((float) (xmax * i) / dim));
    if (i <= localdoublepoint)
      xspot = (int) (xtemp * rescalerbelow);
    else
      xspot = (int) (rescaledoublepoint * rescalerbelow + (xtemp - rescaledoublepoint) * rescalerabove);
    xspot += (int) 11;
    if (xspot < xlen - 1 && yspot >= 1)
      graph[yspot][xspot] = '*';
  }
  for (i = 0; i < ylen; i++)
    FP "%s \n", graph[i]);
  if (notplot)
    FP "    - not all values plotted, some exceed bounds: %5.4f - %5.0f \n", ymin, ymax);
  FP "\n");
}                               /* asciitrend */

// Constants for asciicurve ACXPLOT and AXYPLOT  determine the area of the plot
#define  ACXLEFTSPACE 10
#define  ACXPLOT  75
#define ACXMAX (ACXLEFTSPACE + ACXPLOT)
#define ACYPLOT 50
#define  ACYBOTSPACE 3
#define ACYMAX (ACYPLOT + ACYBOTSPACE)
#define  ASCIICURVEMINVAL 1e-6

/*for logscale,  any integer not 0  makes the plot use a log scale */
/* does not set the correct scale for t */
/* position of rows of graph 
0 name
1               |
ACYPLOT         |
ACYPLOT + 1     -------
ACYPLOT + 2    
*/

void asciicurve (FILE * outfile, struct plotpoint *a, char *qlabel,
                 int logscale, int recordstep)
{
  char graph[ACYMAX][ACXMAX + 1];
  char tc[13];
  double ymax = -1e10, ymin = 1e10;
  int xmax;
  int i, xspot, yspot;

  for (i = 0; i < ACYMAX; i++)
    graph[i][0] = '\0';
  for (i = 0; i < GRIDSIZE; i++)
  {
    if (ymax < a[i].y)
      ymax = a[i].y;
    if (ymin > a[i].y)
      ymin = a[i].y;
  }
  ymax /= recordstep;
  ymin = 0;                     // don't shift plot on y axis 
  if (logscale)
  {
    xmax = GRIDSIZE - 1;
  }
  else
  {
    xmax = -1;
    i = GRIDSIZE - 1;
    while (xmax == -1)
    {
      if (fabs (a[i].y) > ASCIICURVEMINVAL)
        xmax = i;
      i--;
    }
    if (xmax < 0)
      xmax = GRIDSIZE - 1;
  }
// set up plot name line 
  strcat (graph[0], qlabel);
  strcat (graph[0], " curve");
  while (strlen (graph[0]) < ACXMAX)
    strcat (graph[0], " ");
//set up the upper y axis label
  sprintf (tc, "%8.4g", ymax);
  strcpy (graph[1], tc);
  while (strlen (graph[1]) < ACXMAX)
    strcat (graph[1], " ");
//set up the lower y axis label
  sprintf (tc, "%8.4f", ymin);
  strcpy (graph[ACYPLOT], tc);
  while (strlen (graph[ACYPLOT]) < ACXMAX)
    strcat (graph[ACYPLOT], " ");
//set up the horizontal line on the x axis
  for (i = 0; i < ACXLEFTSPACE; i++)
    strcat (graph[ACYPLOT + 1], " ");
  while (strlen (graph[ACYPLOT + 1]) < ACXMAX)
    strcat (graph[ACYPLOT + 1], "-");
// set up the x axis label line
  for (i = 0; i < ACXLEFTSPACE; i++)
    strcat (graph[ACYPLOT + 2], " ");
  sprintf (tc, "%8.4f", a[0].x);
  strcat (graph[ACYPLOT + 2], tc);
  sprintf (tc, "%8.4f", a[xmax].x);
  if (logscale)
    strcat (graph[ACYPLOT + 2], "           Log Scale");
  while (strlen (graph[ACYPLOT + 2]) < ACXMAX - strlen (tc) - 1)
    strcat (graph[ACYPLOT + 2], " ");
  strcat (graph[ACYPLOT + 2], tc);
  while (strlen (graph[ACYPLOT + 2]) < ACXMAX)
    strcat (graph[ACYPLOT + 2], " ");
// fill the plot with empty space
  for (i = 2; i < ACYPLOT; i++)
    while (strlen (graph[i]) < ACXMAX)
      strcat (graph[i], " ");
  // set up the vertical line on the Y axis
  for (i = 1; i < ACYPLOT + 1; i++)
    graph[i][ACXLEFTSPACE] = '|';
  i = 0;
  while (i <= xmax)
  {
    yspot = (int) 1 + ACYPLOT - (int) (ACYPLOT * (a[i].y / recordstep - ymin) / (ymax - ymin));
    if (logscale)
      xspot = (int) ACXLEFTSPACE + 1 + (int) ((ACXPLOT - 2) * (log (a[i].x) - log (a[0].x)) / (2 * log (a[xmax].x)));
    else
      xspot = (int) ACXLEFTSPACE + 1 + (int) ((ACXPLOT - 2) * (a[i].x - a[0].x) / (a[xmax].x - a[0].x));
    assert (xspot < ACXMAX);
    assert (yspot < ACYMAX);
    if (xspot < ACXMAX && xspot >= (ACXLEFTSPACE + 1) && yspot < ACYMAX && yspot >= 1)
      graph[yspot][xspot] = '*';
    i++;
  }
  for (i = 0; i < ACYMAX; i++)
    FP "%s \n", graph[i]);
  fprintf (outfile, "\n");
}                               /* asciicurve */

// print acceptance rates for each the numrec elements pointed to in the array of pointers rec 
// assume that all elements have the same value for num_uptypes
#undef MAXLENX
#undef MAXLENY
void printacceptancerates (FILE * outto, int numrec,
                           struct chainstate_record_updates_and_values *rec[],
                           const char *printstring)
{

  int i, j;
  char numstr[20];

  fprintf (outto, "\n%s\n", printstring);
  for (i = 0; i < (int) strlen (printstring); i++)
    fprintf (outto, "-");
  fprintf (outto, "\n");

  fprintf (outto, "Update Type:");
  for (i = 0; i < rec[0]->num_uptypes; i++)
    fprintf (outto, "\t%s\t", rec[0]->upnames[i]);
  fprintf (outto, "\n");
  fprintf (outto, "            ");
  for (i = 0; i < rec[0]->num_uptypes; i++)
    fprintf (outto, "\t#Tries\t#Accp\t%%");
  fprintf (outto, "\n");
  for (j = 0; j < numrec; j++)
  {
    fprintf (outto, " %s", rec[j]->str);
    i = (int) strlen (rec[j]->str);
    while (i < 13)
    {
      fprintf (outto, " ");
      i++;
    }
    for (i = 0; i < rec[j]->num_uptypes; i++)
    {
      sprintf (&numstr[0], "%.2e", (float) rec[j]->upinf[i].tries);
      fprintf (outto, "\t%s", shorten_e_num (&numstr[0]));
      sprintf (&numstr[0], "%.2e", (float) rec[j]->upinf[i].accp);
      fprintf (outto, "\t%s", shorten_e_num (&numstr[0]));
      if (rec[j]->upinf[i].tries > 0)
        fprintf (outto, "\t%.2f", (float) 100 * rec[j]->upinf[i].accp / rec[j]->upinf[i].tries);
      else
        fprintf (outto, "\tna");
    }
    fprintf (outto, "\n");
  }
}                               /* printacceptancerates() */

void printacceptancerates_multichain (FILE * outto)
{
  int ci;
  int i;
  char numstr[20];

  fprintf (outto, "\n");
  fprintf (outto, "Update Rates -- Assignment per chain\n");
  fprintf (outto, "------------------------------------\n");

  fprintf (outto, "Update Type:");
  for (i = 0; i < Cupinf[0].num_uptypes; i++)
    {
      fprintf (outto, "\t%s\t\t", Cupinf[0].upnames[i]);
    }
  fprintf (outto, "\n");
  fprintf (outto, "            ");
  for (i = 0; i < Cupinf[0].num_uptypes; i++)
    fprintf (outto, "\t#Tries\t#Accp\t%%");
  fprintf (outto, "\n");

  for (ci = 0; ci < numchains; ci++)
  {
    fprintf (outto, " %s", Cupinf[ci].str);
    i = (int) strlen (Cupinf[ci].str);
    while (i < 13)
    {
      fprintf (outto, " ");
      i++;
    }
    for (i = 0; i < Cupinf[ci].num_uptypes; i++)
    {
      sprintf (&numstr[0], "%.2e", (float) Cupinf[ci].upinf[i].tries);
      fprintf (outto, "\t%s", shorten_e_num (&numstr[0]));
      sprintf (&numstr[0], "%.2e", (float) Cupinf[ci].upinf[i].accp);
      fprintf (outto, "\t%s", shorten_e_num (&numstr[0]));
      if (Cupinf[ci].upinf[i].tries > 0)
        fprintf (outto, "\t%.2f", (float) 100 * Cupinf[ci].upinf[i].accp / Cupinf[ci].upinf[i].tries);
      else
        fprintf (outto, "\tna");
    }
    fprintf (outto, "\n");
  }
}

void printcurrentvals (FILE * outto)
{
  int li, i;
  if (numsplittimes > 0)
  {
    fprintf (outto, "\nCurrent Splitting Time Values\n");
    fprintf (outto, "------------------------------\n");
    if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 0
        && npops > 1)
    {
      fprintf (outto, "Split times");
      for (i = 0; i < lastperiodnumber; i++)
        fprintf (outto, "   t%i: %.3f", i, C[0]->tvals[i]);
      fprintf (outto, "\n\n");
    }
  }
  if (runoptions[PRINTMUTATIONUPDATESTOSCREEN] && nloci > 1)
  {
    fprintf(outto, "Current Locus Likelihoods and Mutation Rate Scalars\n");
    fprintf(outto, "---------------------------------------------------\n");
    fprintf (outto, "Locus#    p(D|G)     u  ");

    fprintf (outto, "\n");
    for (li = 0; li < nloci; li++)
    {
      for (i = 0; i < L[li].nlinked; i++)
      {
        fprintf (outto, "%2i", li);
        if (L[li].nlinked > 1)
          fprintf (outto, "_%-2d", i);
        else
          fprintf (outto, "   ");
        fprintf (outto, " %9.3lg %9.4lg", C[0]->G[li].pdg_a[i],
                 C[0]->G[li].uvals[i]);
        fprintf (outto, "\n");
      }
    }
  }
}                               /* printcurrentvals */



void savegenealogyfile (char genealogyinfosavefilename[], FILE * genealogyinfosavefile,
                   int *lastgenealogysavedvalue, int gsampinflength)
{
  int j, i;
  int gp; // VS gp
  if ((genealogyinfosavefile = fopen (genealogyinfosavefilename, "a")) == NULL)
  {
    //printf ("Error opening treeinfosave file for writing\n");
    IM_err(IMERR_CREATEFILEFAIL,"Error opening treeinfosave file for writing");
  }
  for (j = *lastgenealogysavedvalue + 1; j < genealogiessaved; j++)
  {
    // save everything as a float  - round integers later, when they get used
    for (i = 0; i < gsampinflength; i++)
      fprintf (genealogyinfosavefile, "%.6f\t", (float) gsampinf[gp][j][i]); // VS gp
    fprintf (genealogyinfosavefile, "\n");
  } 
  *lastgenealogysavedvalue = genealogiessaved - 1;
  f_close (genealogyinfosavefile);
  return;
}                               /* savegenealogyfile */


// VS
// SAVEGENEALOGYFILE_VS
// VS void savegenealogyfile_vs gets as input the index of the loci, and allows to print
// the gsampinf for a given locus
// INPUT:
//	int gp : index of grousp of loci loci (note: if 0<gp<nbgrouploci_theta it refers to theta groups;
//			 if nbgrouploci_theta<gp<nbgrouploci_theta+nbgrouploci_mig it refers to mig groups)
//	char genealogyinfosavefilename[] : string with file name
//	FILE * genealogyinfosavefile : stream to the FILE to save the output
//	int *lastgenealogysavedvalue : pointer to integer with the last tree index saved, e.g. 100 (it means that the file already has 100 genealogies saved)
//	int gsampinflength : length of the gsampinf matrix (nb of rows)
void savegenealogyfile_vs (int gp, char genealogyinfosavefilename[], FILE * genealogyinfosavefile,
                   int lastgenealogysavedvalue, int gsampinflength)
{
  int j, i;
 
  // If it is the first genealogy to be saved clear the file by using "w" to write instead of "a" when opening the stream
  if(lastgenealogysavedvalue==-1) {
  	if ((genealogyinfosavefile = fopen (genealogyinfosavefilename, "w")) == NULL)
  	{
    		IM_err(IMERR_CREATEFILEFAIL,"Error opening treeinfosave file for writing");
  	}
  	fprintf (genealogyinfosavefile, "VALUESSTART\n"); // VS added this as the first line so that the LTR tests and marginal distributions work
  }
  else {
  	if ((genealogyinfosavefile = fopen (genealogyinfosavefilename, "a")) == NULL)
  	{
    		IM_err(IMERR_CREATEFILEFAIL,"Error opening treeinfosave file for writing");
  	}
  }

  for (j = lastgenealogysavedvalue + 1; j < genealogiessaved; j++)
  {
    // save everything as a float  - round integers later, when they get used
    for (i = 0; i < gsampinflength; i++)
      fprintf (genealogyinfosavefile, "%.6f\t", (float) gsampinf[gp][j][i]);
    fprintf (genealogyinfosavefile, "\n");
  }
  
  f_close (genealogyinfosavefile);
  return;
}                               /* savegenealogyfile_vs */



// VS
// SAVEASSIGNLOCIFILE
// Saves the sampled assignment vectors into a file
// INPUT:
//	char genealogyinfosavefilename[] : input file name
//	int lastgenealogysavedvalue : number of genealogies already saved (this value is not changed in this function, in contrast to what happens in savegenealogyfile_vs)
// 	int **sampledassign : pointer to the matrix with the saved assigned vectors
void saveassignlocifile(char assignlociinfosavefilename[], int lastgenealogysavedvalue, int **sampledassign) {
	
  FILE * assignsavefile;
  int j, i;

  // If it is the firt assignment vector to be saved, clear the file by using "w" to write instead of "a" when opening the stream
  if(lastgenealogysavedvalue==-1) {
	if ((assignsavefile = fopen (assignlociinfosavefilename, "w")) == NULL) {
    		IM_err(IMERR_CREATEFILEFAIL,"Error opening assignment loci file for writing");
	}
    //fprintf (assignsavefile, "VALUESSTART\n"); // VS added this as the first line so that the LTR tests and marginal distributions work
  }
  else {
  	if ((assignsavefile = fopen (assignlociinfosavefilename, "a")) == NULL)
  	{
    		//printf ("Error opening treeinfosave file for writing\n");
    		IM_err(IMERR_CREATEFILEFAIL,"Error opening treeinfosave file for writing");
  	}
  }

  for (j = lastgenealogysavedvalue+1; j < genealogiessaved; j++)
  {
    // save the assignment loci vectors
    for (i = 0; i < nloci; i++)
      fprintf (assignsavefile, "%i\t", sampledassign[j][i]);
    fprintf (assignsavefile, "\n");
  }
  f_close (assignsavefile);
  return;
}




void print_means_variances_correlations (FILE * outfile)
{
  double *means, *variances, **correlations;
  int i, p, q, np;

  int gp; // VS gp
  int nbparams, nbthetaparams, nbmigparams; // VS

  // VS compute some variables that will simplify the code in the for loops and allocation of variables
  // nbparams=(numpopsizeparams*nbgroupsloci_theta)+(nummigrateparams*nbgroupsloci_mig); 
  nbthetaparams=numpopsizeparams*nbgroupsloci_theta;
  nbmigparams=nummigrateparams*nbgroupsloci_mig;


  np = (numpopsizeparams*nbgroupsloci_theta)+(nummigrateparams*nbgroupsloci_mig); // VS numpopsizeparams + nummigrateparams; // VS - this is the same as nbparams

  means = calloc ((size_t) np, sizeof (double)); // VS nbparams instead of nb
  variances = calloc ((size_t) np, sizeof (double)); // VS nbparams instead of nb
  if (npops > 1)
  {
    correlations = orig2d_alloc2Ddouble (np, np); // VS nbparams instead of np
  }


  for (i = 0; i < genealogiessaved; i++)
  {
    for (p = 0; p < np; p++) 
    {
      means[p] += calcx (i, p, 0);
      variances[p] += calcx (i, p, 1);
    }
  }
  for (p = 0; p < np; p++) 
  {
    if (means[p] >= 0.0)
      means[p] /= genealogiessaved;
    if (variances[p] >= 0.0)
    {
      variances[p] /= genealogiessaved;
      variances[p] -= SQR (means[p]);
    }
  }
  if (npops > 1)
  {
    for (p = 0; p < np; p++)
      for (q = 0; q < np; q++)
        correlations[p][q] = 0;
    for (i = 0; i < genealogiessaved; i++)
    {
      for (p = 0; p < np - 1; p++)
        for (q = p + 1; q < np; q++)
          correlations[p][q] += calcx (i, p, 0) * calcx (i, q, 0);
    }
    for (p = 0; p < np - 1; p++)
      for (q = p + 1; q < np; q++)
      {
        if (correlations[p][q] >= 0.0)
        {
          correlations[p][q] /= genealogiessaved;
          correlations[p][q] -= (means[p] * means[q]);
          correlations[p][q] /= sqrt (variances[p] * variances[q]);
        }
        else
          correlations[p][q] = -1.0;
      }
  }
  FP "\nMEANS, VARIANCES and CORRELATIONS OF PARAMETERS ('$' r > 0.4  '*' r > 0.75)\n");
  FP "=============================================================================\n");
  FP "Param:");
  for (p = 0; p < np; p++)
  {
    if (p < nbthetaparams) // VS numpopsizeparams)
    { // VS print also the group information as _g0, _g1, ...
      FP "\t%s_g%i", itheta[g_paramindex[p]].str, g_gp_ip[p]); // VS g_paramindex[p] instead of p
    }
    else
    { // VS print also the group information as _g0, _g1, ...
      if (imig[g_paramindex[p]].pr[g_gp_ip[p]].max > MPRIORMIN) // VS imig[p - numpopsizeparams].pr[gp].max > MPRIORMIN)// VS gp
          FP "\t%s_g%i", imig[g_paramindex[p]].str, g_gp_ip[p]); // VS FP "\t%s", imig[p - numpopsizeparams].str);
    }
  }
  FP "\nMean:");
  for (p = 0; p < np; p++)
    FP "\t%-.3lf", means[p]);
  FP "\nStdv:");
  for (p = 0; p < np; p++)
    FP "\t%-.3lf", sqrt (variances[p]));
  if (npops > 1)
  {
    FP "\n\nCorrelations\n");
    for (p = 0; p < np; p++)
    {
      if (p < nbthetaparams) // VS (p < numpopsizeparams)
      {
        FP "\t%s_g%i", itheta[g_paramindex[p]].str, g_gp_ip[p]); // VS (p < numpopsizeparams)
      }
      else  //p < nbparams
      {
         if (imig[g_paramindex[p]].pr[g_gp_ip[p]].max > MPRIORMIN) // imig[p - numpopsizeparams].pr[gp].max > MPRIORMIN)  
            FP "\t%s_g%i", imig[g_paramindex[p]].str, g_gp_ip[p]); // VS imig[p - numpopsizeparams].str);
      }
    }
    FP "\n");
    for (p = 0; p < np; p++)
    {
      if (p < nbthetaparams) // VS (p < numpopsizeparams)
        FP "%s_g%i", itheta[g_paramindex[p]].str, g_gp_ip[p]); // VS itheta[p].str
      else //p < numpopsizeparams + nummigrateparams
      {
        FP "%s_g%i", imig[g_paramindex[p]].str, g_gp_ip[p]); // VS FP "%s_g%i", imig[p - numpopsizeparams].str);
      }
      for (q = 0; q < np; q++)
      {
        if (q == p)
          FP "\t  - ");
        else
        {
          if (q > p)
          {
            if (fabs (correlations[p][q]) < 0.4)
            {
              FP "\t%-.3lf", correlations[p][q]);
            }
            else
            {
              if (fabs (correlations[p][q]) < 0.75)
                FP "\t%-.3lf%s", correlations[p][q], "$");
              else
                FP "\t%-.3lf%s", correlations[p][q], "*");
            }
          }
          else
          {
            if (fabs (correlations[q][p]) < 0.4)
            {
              FP "\t%-.3lf", correlations[q][p]);
            }
            else
            {
              if (fabs (correlations[q][p]) < 0.75)
                FP "\t%-.3lf%s", correlations[q][p], "$");
              else
                FP "\t%-.3lf%s", correlations[q][p], "*");
            }
          }

        }
      }
      FP "\n");
    }
  }
  FP "\n");

  XFREE (means);
  XFREE (variances);
  if (npops > 1)
  {
    orig2d_free2D ((void **) correlations, np);
  }
  return;
}                               // print_means_variances_correlations
