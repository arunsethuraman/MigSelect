/* MigSelect 2013-2019 Arun Sethuraman, Vitor Sousa, Jody Hey */
/* This code was developed based on IMa2  2009-2010  Jody Hey, Rasmus Nielsen and Sang Chul Choi*/

/* calculations for assessing if one parameter is greater than another 
print matrices of results */ 

/* JH 4/23/2010  extensively revised this file */ 

#undef GLOBVARS
#include "imamp.h"

static int cci, ccj, wi, wj;
static int treeinc, hitreenum, numtreesused;
static double fci, fcj, hval, denom, qmax, fmi,fmj, mmax;
static double pgt_fcj_gt_0 (double qi);
static double mgt_wj_gt_0 (double mi);
static double  trapzd(double  (*func)(double ), double  a, double  b, int n);
static double  qtrap(double  (*func)(double ), double  a, double  b);
static double gtmig(int mi, int mj);
static double gtpops(int pj, int pi);

#define SWITCH_TO_LOWERGAMMA_CRIT  1e-15   // when the difference between a gamma and uppergamma is smaller than this, use the lowergamma
#define USETREESMAX  20000   // only use this many sampled genealogies for calculations of probabilities that one param is greater than another

double mgt_wj_gt_0 (double mi)
{
  double temp1, temp2, a,b;
  if (mi < MINPARAMVAL)
    return 0.0;
  a = logfact[wj];
  b = uppergamma(wj+1,fmj*mi);
  if (a<=b)
    //IM_err (IMERR_LOGDIFF, " in mgt_wj_gt_0() mi %d a %lf  b %lf",mi,a,b);
    return 0.0;  // sometimes get here for very low mi values, floating point thing 
  if ((a-b) < SWITCH_TO_LOWERGAMMA_CRIT)
  {
    temp1 = lowergamma (wj+1,fmj*mi);
  }
  else
  {
    LogDiff(temp1,a,b);
  }
  temp2 =  wi*log(mi) - fmi*mi -(wj+1)*log(fmj) + temp1;
  temp2 -= denom;
  return exp(temp2);
}

double pgt_fcj_gt_0 (double qi)
{
  double fcj2, fci2,temp1, temp2, temp3, a,b;

  if (qi < MINPARAMVAL)
    return 0.0;
  fcj2 = 2*fcj;
  fci2 = 2*fci;
  if (ccj == 0)
  {
    a = log(qi)-fcj2/qi;
    b = log(fcj2)+ uppergamma(0,fcj2/qi);
    if (a>b)
    {
      LogDiff(temp1,a,b);
      temp2 =  - fci2/qi + cci*log(2/qi) ; 
      temp3  = temp1 + temp2 - hval - denom; 
      return exp(temp3);
    }
    else 
      return 0.0;
  }
  else
  {
    temp1 =  uppergamma(ccj-1,fcj2/qi);
    temp2 = LOG2 + cci*log(2/qi) + (1-ccj)*log(fcj) -fci2/qi; 
    temp3 = temp2 + temp1 - hval - denom;
    return exp(temp3);
  }
}

#define FUNC(x) ((*func)(x))

double  trapzd(double  (*func)(double ), double  a, double  b, int n)
{
   double  x,tnm,sum,del;
   static double  s;
   int it,j;

   if (n == 1) {
      return (s=0.5*(b-a)*(FUNC(a)+FUNC(b)));
   } else {
      for (it=1,j=1;j<n-1;j++) 
        it <<= 1;
      tnm=it;
      del=(b-a)/tnm;
      x=a+0.5*del;
      for (sum=0.0,j=1;j<=it;j++,x+=del) 
        sum += FUNC(x);
      s=0.5*(s+(b-a)*sum/tnm);
      return s;
   }
}
#undef FUNC
/* (C) Copr. 1986-92 Numerical Recipes Software '$&'3$. */

#define EPS 1.0e-4
#define JMAX 20

double  qtrap(double  (*func)(double ), double  a, double  b)
{
   /* void nrerror(char error_text[]); */
   int j;
   double  s,olds;
   olds = -1.0e100;
   for (j=1;j<=JMAX;j++) {
      s=trapzd(func,a,b,j);
      if (j > 5) //Avoid spurious early convergence.
      {
        if (fabs(s-olds) < EPS*fabs(olds) ||(s == 0.0 && olds == 0.0)) 
          return s;
      }
      olds=s;
   }
   return s;
   //nrerror("Too many steps in routine qtrap");
}
#undef EPS
#undef JMAX
/* (C) Copr. 1986-92 Numerical Recipes Software '$&'3$. */

// gtmig
// INPUT
//	int mi: index of migration rate 1
//	int mj: index of migration rate 2
// NOTE: these indexes refer to the migration rate assuming a large array 
//	with all the parameters put together.
// Thus, mi and mj must be larger than nbgroupsloci*numpopsizeparam
double gtmig(int mi, int mj)
{
  int ei;
  //int gp; // VS
  double sum, temp, temp1, temp2, temp3, temp4, a,b,c;
  double  (*func)(double );
  //int nbthetaparams; // VS

  func = mgt_wj_gt_0;

  mmax = imig[g_paramindex[mi]].pr[g_gp_ip[mi]].max; // VS g_paramindex[mi] and g_gp_ip[mi] - TO BE DONE LATER - check how to deal with cases where mmax varies between groups
  sum = 0;

  //nbthetaparams = numpopsizeparams*nbgroupsloci_theta; // VS
  
  for (ei = 0; ei < hitreenum; ei+= treeinc) 
  {
	// VS for accessing gsampinf, the index of group is given by g_gp_gs[mi] or g_gp_gs[mj]
    fmi = gsampinf[g_gp_gs[mi]][ei][gsamp_fmp + g_paramindex[mi]]; // VS gp=g_gp_gs[mi]
    fmj = gsampinf[g_gp_gs[mj]][ei][gsamp_fmp + g_paramindex[mj]]; // VS gp=g_gp_gs[mj]
    wi = (int) gsampinf[g_gp_gs[mi]][ei][gsamp_mcp + g_paramindex[mi]];  // VS gp=g_gp_gs[mi]
    wj = (int) gsampinf[g_gp_gs[mj]][ei][gsamp_mcp + g_paramindex[mj]];  // VS gp=g_gp_gs[mj]
    denom = gsampinf[g_gp_gs[mi]][ei][gsamp_mip + g_paramindex[mi]]+gsampinf[g_gp_gs[mj]][ei][gsamp_mip + g_paramindex[mj]];  // VS gp=g_gp_gs[mi] or gp=g_gp_gs[mj]
    if (wj==0) 
    {
      if (fmj>0.0) //
      {
        if (wi > 0) //fmi > 0.0
        {
          a = logfact[wi];
          b = uppergamma(wi+1,(fmi + fmj)*mmax);
          c = uppergamma(wi+1,fmi*mmax);
          if (a<=b || a <= c)
            temp = 0.0;
            //IM_err (IMERR_LOGDIFF, " pos 1 in gtmig() mi %d  mj %d a %lf  b %lf",mi,mj,a,b);
          else
          {
            if ((a-b) < SWITCH_TO_LOWERGAMMA_CRIT)
            {
              temp1 = lowergamma (wi+1,(fmi + fmj)*mmax);
            }
            else
            {
              LogDiff(temp1,a,b);
            }
            temp1 += -(wi+1)*log(fmi + fmj);
            if ((a-c) < SWITCH_TO_LOWERGAMMA_CRIT)
            {
              temp2 = lowergamma (wi+1,fmi*mmax);
            }
            else
            {
              LogDiff(temp2,a,c);
            }
            temp2 += -(wi+1)*log(fmi);
            if (temp2<=temp1)
              temp = 0.0;
              //IM_err (IMERR_LOGDIFF, " pos 3 in gtmig() mi %d  mj %d temp2 %lf  temp1 %lf",mi,mj,temp2,temp1);
            else
            {
              LogDiff(temp3,temp2,temp1);
              temp4 = temp3 - log(fmj) - denom;
              temp = exp(temp4);
            }
          }
        }
        else
        {
          if (fmi > 0.0)
          {
            temp1 = fmi*(exp(-(fmi+fmj)*mmax)- exp(-fmi*mmax) );
            temp2 = fmj*(1-exp(-fmi*mmax));
            temp3 = (temp1 + temp2) / (fmi*fmj*(fmi+fmj));
            temp4 = log(temp3) - denom;
            temp = exp(temp4);
          }
          else // fmi == 0.0 
          {
            temp1 = (fmj*mmax - 1.0 + exp(-fmj*mmax))/(fmj*fmj);
            temp4 = log(temp1) - denom;
            temp = exp(temp4);
          }
        }
      }
      else  //fmj == 0.0 && wj==0
      {
        if (wi > 0) //fmi > 0.0
        {
          a = logfact[wi+1];
          b = uppergamma(wi+2,fmi*mmax);
          if (a<=b)
            temp = 0.0;
            //IM_err (IMERR_LOGDIFF, " pos 4 in gtmig() mi %d  mj %d a %lf  b %lf",mi,mj,a,b);
          else
          {
            if ((a-b) < SWITCH_TO_LOWERGAMMA_CRIT)
            {
              temp1 = lowergamma (wi+2,fmi*mmax);
            }
            else
            {
              LogDiff(temp1,a,b);
            }
            temp2 = -(wi+2) * log(fmi);
            temp3 = temp2 + temp1 - denom;
            temp = exp(temp3);
          }
        }
        else
        {
          if (fmi > 0.0)
          {
            temp1 = (1.0 - exp(-fmi * mmax)*(fmi*mmax + 1.0))/(fmi*fmi);
            temp2 = log(temp1)-denom;
            temp = exp(temp2);
          }
          else // fmi == 0.0 
          {
            temp1 = log(mmax*mmax/2.0) - denom;
            temp = exp(temp1); 
          } 
        }
      }
    }
    else
    {
      temp = qtrap(func, MINPARAMVAL, mmax);
    }
    temp = DMIN(1.0, temp);  // JH 4/22/2010 values should not be greater than 1 , except for small fp issues
    sum += temp;
  }
  sum /= numtreesused;
  return sum;

} // gtmig

/* return the probability that pram pi >param pj */
double gtpops(int pi, int pj)
{
  int ei;
  //int gp; // VS
  double sum, temp, temp1, temp2, temp3;
  double  (*func)(double );

  func = pgt_fcj_gt_0;

  qmax = itheta[pi].pr[g_gp_ip[pi]].max; // VS g_gp_ip[pi]

  sum = 0;
  for (ei = 0; ei < hitreenum; ei+= treeinc)
  {
    cci = (int) gsampinf[g_gp_gs[pi]][ei][gsamp_ccp + g_paramindex[pi]]; // VS gp=g_gp_gs[pi]
    ccj = (int) gsampinf[g_gp_gs[pj]][ei][gsamp_ccp + g_paramindex[pj]]; // VS gp=g_gp_gs[pj]
    fci = gsampinf[g_gp_gs[pi]][ei][gsamp_fcp + g_paramindex[pi]]; // VS gp=g_gp_gs[pi]
    fcj = gsampinf[g_gp_gs[pj]][ei][gsamp_fcp + g_paramindex[pj]]; // VS gp=g_gp_gs[pj]
    hval = gsampinf[g_gp_gs[pi]][ei][gsamp_hccp + g_paramindex[pi]]+gsampinf[g_gp_gs[pj]][ei][gsamp_hccp + g_paramindex[pj]]; // VS gp=g_gp_gs[pi] or gp=g_gp_gs[pj]
    denom = gsampinf[g_gp_gs[pi]][ei][gsamp_qip + g_paramindex[pi]]+gsampinf[g_gp_gs[pj]][ei][gsamp_qip + g_paramindex[pj]]; // VS gp=g_gp_gs[pi] or gp=g_gp_gs[pj]

    if (ccj==0)
    {
      if (fcj == 0)
      {
        if (fci==0)
        {
          temp = exp(2.0 * log(qmax)- LOG2 - hval - denom);
        }
        else
        {
          if (cci>= 2)
          {
            temp1 = 2*LOG2 + (2-cci)* log(fci);
            temp2 = uppergamma(cci-2,2*fci/qmax);
            temp3 = temp1 + temp2 - hval-denom;
            temp = exp(temp3);
          }
          else
          {
            if (cci == 1)
            {
              temp1 = 4*fci*exp(uppergamma(0,2*fci/qmax));
              temp2 = 2*qmax*exp(-2*fci/qmax) - temp1;
              temp3 = log(temp2) - hval - denom;
              temp=exp(temp3); // JH added 4/22/2010 
            }
            else  // cci==0
            {
              temp1 = exp(uppergamma(0,2*fci/qmax));  // use for ei() of a negative value 
              temp2 = (qmax/2) * (qmax-2*fci) * exp(-2*fci/qmax) + 2*fci*fci*temp1;
              temp3 = log(temp2) - hval -denom;
              temp=exp(temp3);  // JH added 4/22/2010 
            }
          }
        }
      }
      else //ccj=0 fcj > 0 
      {
        temp = qtrap(func, MINPARAMVAL, qmax);

      }
    }
    else
    {
      temp = qtrap(func, MINPARAMVAL, qmax);
    }
    temp = DMIN(1.0, temp);  // JH 4/22/2010 values should not be greater than 1 except for small fp issues
    sum += temp;
  }
  sum /= numtreesused;
  return sum;
} // gtpops 

#define CHECKGREATERTHAN 0.02  //JH 4/22/2010  to see if reciprocal values sum approximately to 1
void print_greater_than_tests (FILE * outfile)
{
  double **gt_popsize, **gt_mig;
  int i, j;
  //int gp; // VS gp
  int printwarning = 0;
 
  int nbparams, nbthetaparams, nbmigparams; // VS

  if (genealogiessaved > USETREESMAX)
  {
    treeinc = (int) genealogiessaved/(int) USETREESMAX;
    numtreesused = USETREESMAX;
    hitreenum = treeinc * USETREESMAX;
  }
  else
  {
    treeinc = 1;
    hitreenum = numtreesused = genealogiessaved;
  }

  // VS compute some variables that will simplify the code in the for loops and allocation of variables
  nbparams=(numpopsizeparams*nbgroupsloci_theta)+(nummigrateparams*nbgroupsloci_mig); 
  nbthetaparams=numpopsizeparams*nbgroupsloci_theta;
  nbmigparams=nummigrateparams*nbgroupsloci_mig;

  gt_popsize= orig2d_alloc2Ddouble(nbthetaparams, nbthetaparams); // VS nbgroupsloci_theta
  for (i=0;i < nbthetaparams; i++) //JH 4/22/2010  calculate full matrix - VS nbthetaparams
    for (j = 0; j < nbthetaparams; j++) // VS nbthetaparams
    {
      if (i != j && itheta[g_paramindex[i]].pr[g_gp_ip[i]].max == itheta[g_paramindex[j]].pr[g_gp_ip[j]].max) // VS g_paramindex, g_gp_ip[i] and g_gp_ip[j]
        gt_popsize[i][j] = gtpops(i,j); 
      else
        gt_popsize[i][j] = -1;
    } 
  if (!modeloptions[EXPOMIGRATIONPRIOR])
  {
    gt_mig = orig2d_alloc2Ddouble (nbmigparams, nbmigparams); // VS nbgroupsloci_mig
    // VS correct the index for migration parameters. The migration indexes start at nbgroupsloci_theta*numposizeparams
	for (i=nbthetaparams; i < nbparams; i++) //JH 4/22/2010  calculate full matrix - VS nbgroupsloci_theta
      for (j=nbthetaparams; j < nbparams; j++)
      {
        if (i!=j && imig[g_paramindex[i]].pr[g_gp_ip[i]].max > MINPARAMVAL && imig[g_paramindex[j]].pr[g_gp_ip[j]].max > MINPARAMVAL && imig[g_paramindex[i]].pr[g_gp_ip[i]].max ==imig[g_paramindex[j]].pr[g_gp_ip[j]].max) // VS - g_paramindex and g_gp_ip
          gt_mig[i-nbthetaparams][j-nbthetaparams] = gtmig(i,j); // VS correcting the index of the parameters for the matrix gt_mig
        else 
          gt_mig[i-nbthetaparams][j-nbthetaparams] = -1.0; // VS correcting the index of the parameters for the matrix gt_mig
      }
  }
  FP "\nPARAMETER COMPARISONS, PROBABILITY THAT ROW PARAMETER IS GREATER THAN COLUMN PARAMETER\n");
  FP "========================================================================================\n");
  if (calcoptions[USEPRIORFILE]) 
    FP"    Comparisons only done for parameters with identical prior distributions \n");
  FP "Population Sizes\n");
  for (i=0;i<nbthetaparams;i++) { // VS nbgroupsloci_theta
    FP "\t%s_g%i", itheta[g_paramindex[i]].str, g_gp_ip[i]); // VS - also print the group
  }
  FP "\n");
  for (i=0;i<nbthetaparams;i++) // VS nbgroupsloci_theta
  {
    FP "%s_g%i", itheta[g_paramindex[i]].str, g_gp_ip[i]); // VS - also print the group of loci
    for (j=0;j<nbthetaparams;j++) // VS nbgroupsloci_theta
    {
      if (i == j)
        FP "\t  - ");
      else
      {
        if (gt_popsize[i][j] < 0.0)  
          FP"\tna");
        else
        {
          if ( fabs(1.0 - (gt_popsize[i][j] + gt_popsize[j][i])) >= CHECKGREATERTHAN) 
          {
            FP "\t%.3lf?", gt_popsize[i][j]);  //JH 4/22/2010  calculate full matrix
            printwarning = 1;
          }
          else
            FP "\t%.3lf", gt_popsize[i][j]);  //JH 4/22/2010  calculate full matrix
        }
      }
    }
    FP "\n");
  }
  orig2d_free2D ((void **) gt_popsize, nbthetaparams); // VS nbgroupsloci_theta
  FP "\nMigration Rates\n");
  if (modeloptions[EXPOMIGRATIONPRIOR])
    FP"  NOT IMPLEMENTED FOR MIGRATION RATES WITH EXPONENTIAL PRIORS \n");
  else
  {
    for (i=nbthetaparams;i<nbparams;i++) // VS - changed these for loops
    {
      FP "\t%s_g%i", imig[g_paramindex[i]].str, g_gp_ip[i]); // VS also print the group
    }
    FP "\n");
    for (i=nbthetaparams;i<nbparams;i++) // VS - changed these for loops for the migration parameters
    {
      FP "%s_g%i", imig[g_paramindex[i]].str, g_gp_ip[i]);
      for (j=nbthetaparams;j<nbparams;j++) // VS - changed these for loops for the migration parameters
      {
        if (i == j)
          FP "\t  - ");
        else
        {
          if (gt_mig[g_paramindex[i]][g_paramindex[j]] < 0.0) // VS - here index i and j need to be 0<i<nbmigparams, as g_paramindex
            FP"\tna");
          else
          {
            if ( fabs(1.0 - (gt_mig[g_paramindex[i]][g_paramindex[j]] + gt_mig[g_paramindex[j]][g_paramindex[i]])) >= CHECKGREATERTHAN) // VS
            {
              FP "\t%.3lf?", gt_mig[g_paramindex[i]][g_paramindex[j]]);  //JH 4/22/2010  calculate full matrix // VS
              printwarning = 1;
            }
            else
              FP "\t%.3lf", gt_mig[g_paramindex[i]][g_paramindex[j]]);  //JH 4/22/2010  calculate full matrix // VS
          }
        }
      }
      FP "\n");
    }
    orig2d_free2D ((void **) gt_mig, nbmigparams); // VS nbnmigparams
  }
  if (printwarning)
    FP"  ""?"" indicates that reciprocal values do not sum to approximately 1, possibly due to a small sample of genealogies\n");
  FP "\n\n");
  return;
}                               // print_greater_than_tests
#undef CHECKGREATERTHAN 

#undef USETREESMAX
#undef SWITCH_TO_LOWERGAMMA_CRIT
