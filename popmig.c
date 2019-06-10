/* MigSelect 2013-2019 Arun Sethuraman, Vitor Sousa, Jody Hey */
/* This code was developed based on IMa2  2009-2010  Jody Hey, Rasmus Nielsen and Sang Chul Choi*/

#undef GLOBVARS
#include "imamp.h"


// calc_probmig
// calculate the probability density of the product 2NM
// return the density for a given value of 2NM=x
// INPUT
//	int thetai: index of theta parameter
//	int mi: index of mig parameter
//	double x: point at which the function is evaluated (i.e. the x axis)
//	int prob_or_like: ???
//	int gp_theta: index of group theta 
//	int gp_mig: index of group mig
// NOTE: the index gp_mig for imig starts at zero, but not for gsampinf structure.
// For gsampinf structure gp_mig must be added to nbgroupsloci_theta to obtain the correct index.
double
calc_popmig (int thetai, int mi, double x, int prob_or_like, int gp_mig, int gp_theta)
{

  //int gp_mig; // VS gp_mig
  //int gp_theta; // VS gp_theta

  int ei, cc, mc;
  int ccp, fcp, hcp, mcp, fmp, qip, mip;
  double sum, temp, temp1, temp2, fc, fm, hc, qintg, mintg, mmax, qmax;
  double a,b;
  ccp = gsamp_ccp + thetai;
  fcp = gsamp_fcp + thetai;
  hcp = gsamp_hccp + thetai;
  mcp = gsamp_mcp + mi;
  fmp = gsamp_fmp + mi;
  qip = gsamp_qip + thetai;
  mip = gsamp_mip + mi;
  mmax = imig[mi].pr[gp_mig].max; // VS gp
  qmax = itheta[thetai].pr[gp_theta].max; // VS gp

  for (sum = 0, ei = 0; ei < genealogiessaved; ei++)

  {
    cc = (int) gsampinf[gp_theta][ei][ccp]; // VS gp
    fc = (double) gsampinf[gp_theta][ei][fcp]; // VS gp
    hc = (double) gsampinf[gp_theta][ei][hcp]; // VS gp
    mc = (int) gsampinf[gp_mig+(nbgroupsloci_theta*numpopsizeparams)][ei][mcp]; // VS gp
    fm = (double) gsampinf[gp_mig+(nbgroupsloci_theta*numpopsizeparams)][ei][fmp]; // VS gp
    qintg = (double) gsampinf[gp_theta][ei][qip]; // VS gp
    mintg = (double) gsampinf[gp_mig+(nbgroupsloci_theta*numpopsizeparams)][ei][mip]; // VS gp
    if (fc == 0 && cc == 0 && fm > 0)
    {
      temp1 = LOG2 - (mc * log (fm)) - hc - qintg - mintg;
      temp2 = log (exp (uppergamma (mc, 2 * fm * x / qmax)) - exp (uppergamma (mc, mmax * fm)));
    }
    else
    {
      if (fm == 0 && mc == 0 && fc > 0)
      {
        temp1 = LOG2 - (cc * log (fc)) - hc - qintg - mintg;
        temp2 = log (exp (uppergamma (cc, 2 * fc / qmax)) - exp (uppergamma (cc, fc * mmax / x)));
      }
      else
      {
        if (fc == 0 && cc == 0 && mc == 0 && fm == 0)
        {
          temp1 = log (2 * log (mmax * qmax / (2 * x))) - hc - qintg - mintg;;
          temp2 = 0;
        }
        else// fc>0, cc>0, mc> 0, fm> 0
        {
          temp1 = LOG2 + (mc * log (x)) - ((cc + mc) * log (fc + fm * x)) - hc - qintg - mintg;
          //temp2 = log ( exp(uppergamma (cc + mc, 2 * (fc + fm * x) / qmax)) - exp(uppergamma (cc + mc, mmax * (fm + fc / x))) );
          a = uppergamma (cc + mc, 2 * (fc + fm * x) / qmax);
          b = uppergamma (cc + mc, mmax * (fm + fc / x));
          /* the difference between two upppergamma()s is the negative of */
          /* the difference between two lowergamma()s.  However numbers   */
          /* for a or b often come out to be very large and equal under   */
          /* under uppergamma  or underlowergramma, in which case it is   */
          /* necessary to trap this and use the other function            */
          if (a==b)
          {
            b = lowergamma (cc + mc, 2 * (fc + fm * x) / qmax);
            a = lowergamma (cc + mc, mmax * (fm + fc / x));
          }
		  if (a>b)
          {
			LogDiff(temp2,a,b);
          }
		  else
          {
            temp1 = temp2 = 0.0;
          }
        }
      }
    }
    if ((temp1 + temp2 < 700 ) && (temp1 + temp2 > -700 )) // skip things that cannot be exped
      sum += exp (temp1 + temp2);
  }
  sum /= genealogiessaved;
  if (prob_or_like)
  {
    temp = 2 * (log (qmax) + log (mmax) - log (2 * x)) / (qmax * mmax);
    sum /= temp;
  }
  return sum;
}                               //calc_popmig


// calc_pop_expomig
// calculate the probability density of the product 2NM when migration has an exponential prior
// returns the density for a given x value for 2NM=x
// INPUT
//	int thetai: index of theta parameter
//	int mi: index of mig parameter
//	double x: point at which the function is evaluated (i.e. the x axis)
//	int prob_or_like: ???
//	int gp_theta: index of group theta 
//	int gp_mig: index of group mig 
// NOTE: the index gp_mig for imig starts at zero, but not for gsampinf structure.
// For gsampinf structure gp_mig must be added to nbgroupsloci_theta to obtain the correct index.
#define OCUTOFF  10
double
calc_pop_expomig (int thetai, int mi, double x, int prob_or_like, int gp_mig, int gp_theta)
{
  //int gp_mig; // VS gp
  //int gp_theta; // VS gp

  int ei, cc, mc;
  int ccp, fcp, hcp, mcp, fmp, qip, mip;
  double sum, temp, temp1, temp2, temp3,fc, fm, hc, qintg, mintg, mmean, qmax;
  int zadj, maxz = -1000000000;
  double acumm = 0;
  ccp = gsamp_ccp + thetai;
  fcp = gsamp_fcp + thetai;
  hcp = gsamp_hccp + thetai;
  mcp = gsamp_mcp + mi;
  fmp = gsamp_fmp + mi;
  qip = gsamp_qip + thetai;
  mip = gsamp_mip + mi;
  mmean = imig[mi].pr[gp_mig].mean; // VS gp
  qmax = itheta[thetai].pr[gp_theta].max; // VS gp

  for (sum = 0, ei = 0; ei < genealogiessaved; ei++)

  {
    cc = (int) gsampinf[gp_theta][ei][ccp]; // VS gp
    fc = (double) gsampinf[gp_theta][ei][fcp]; // VS gp
    hc = (double) gsampinf[gp_theta][ei][hcp]; // VS gp
    mc = (int) gsampinf[gp_mig+(nbgroupsloci_theta*numpopsizeparams)][ei][mcp]; // VS gp
    fm = (double) gsampinf[gp_mig+(nbgroupsloci_theta*numpopsizeparams)][ei][fmp]; // VS gp
    qintg = (double) gsampinf[gp_theta][ei][qip]; // VS gp
    mintg = (double) gsampinf[gp_mig+(nbgroupsloci_theta*numpopsizeparams)][ei][mip]; // VS gp
    temp1 = x + fc *mmean + fm * mmean * x;
    temp2 = LOG2 - hc - log(mmean) - qintg - mintg;
    temp3 = 2 * temp1/(mmean * qmax);

    if (mc == 0 && cc == 0 )
    {
      temp3 = uppergamma(0,temp3);
      temp2 += temp3;
    }
    else
    {
      if (cc == 0) // mc > 0 
      {
        temp3 = uppergamma(mc,temp3);
        temp2 += temp3 - mc* log(fm + 1/mmean);
      }
      else
      {
        if (mc == 0) //cc > 0
        {
          temp3 = uppergamma(cc,temp3);
          temp2 += temp3 + -cc *log(fc + x*(fm + 1/mmean));; 
        }
        else
        {
          temp3 = uppergamma(mc+cc,temp3);
          temp2 += temp3 - cc*log(x) + (cc+mc)*log(mmean*x/temp1); 
        }
      }
    }
    eexp (temp2, &eexpsum[ei].m, &eexpsum[ei].z);
    if (eexpsum[ei].z > maxz)
      maxz = eexpsum[ei].z;
  }
  for (ei = 0; ei < genealogiessaved; ei++)
  {
    zadj = eexpsum[ei].z - (maxz - OCUTOFF);
    eexpsum[ei].m *= pow (10.0, (double) zadj);
    acumm += eexpsum[ei].m;
  }
  sum = log (acumm) + (maxz - OCUTOFF) * LOG10;
  sum -= log((double) genealogiessaved);
  sum = exp(sum);
  if (prob_or_like)
  {
    temp = 2 * exp(uppergamma(0,2*x/(mmean*qmax)))/ (qmax * mmean);
    sum /= temp;
  }
  return sum;
}                               //calc_pop_expomig
#undef OCUTOFF  



// marginpopmig
// do calculations for finding the peak of the marginal density for a 2NM term
// this is similar to calc_popmig().  It is called by marginalopt_popmig()
// INPUT
//	int mi: index of migration parameter in imig
//	int gp_mig: index of group for migration parameter
//	int gp_theta: index of group for theta parameter
//	int firsttree: index of first tree to be used in computations
//	int lasttree:  index of last tree to be used in computations
//	double x: output of the maximum parameter combination???? 
//	int thetai: index of thete parameter in itheta
// NOTE: the index mi and thetai are different from the indexes in the gsampinf matrix
#define OFFSCALEVAL 1.0
double
marginpopmig (int mi, int gp_mig, int gp_theta, int firsttree, int lasttree, double x, int thetai)
{
  int ei, cc, mc;
  int ccp, fcp, hcp, mcp, fmp, qip, mip;
  double sum, temp1, temp2, fc, fm, hc, qintg, mintg, mmax, qmax;
  double a,b;
  double max, min;

  //int gp_theta; // VS gp
  //int gp_mig; // VS gp

  max = itheta[thetai].pr[gp_theta].max*imig[mi].pr[gp_mig].max/2.0; // VS gp
  min = 0;
  if (x < min || x > max)
    return OFFSCALEVAL;
  
  ccp = gsamp_ccp + thetai;
  fcp = gsamp_fcp + thetai;
  hcp = gsamp_hccp + thetai;
  mcp = gsamp_mcp + mi;
  fmp = gsamp_fmp + mi;
  qip = gsamp_qip + thetai;
  mip = gsamp_mip + mi;
  mmax = imig[mi].pr[gp_mig].max; // VS gp
  qmax = itheta[thetai].pr[gp_theta].max; // VS gp

  for (sum = 0,ei = firsttree; ei < lasttree; ei++)

  {
    // VS - note that need to add nbgroupsloci_theta to gp_mig
	// because the index of gp_mig for imig structure is different from the index
	// in gsampinf structure
    cc = (int) gsampinf[gp_theta][ei][ccp]; // VS gp
    fc = (double) gsampinf[gp_theta][ei][fcp]; // VS gp
    hc = (double) gsampinf[gp_theta][ei][hcp]; // VS gp
    mc = (int) gsampinf[gp_mig+(nbgroupsloci_theta*numpopsizeparams)][ei][mcp]; // VS gp
    fm = (double) gsampinf[gp_mig+(nbgroupsloci_theta*numpopsizeparams)][ei][fmp]; // VS gp
    qintg = (double) gsampinf[gp_theta][ei][qip]; // VS gp
    mintg = (double) gsampinf[gp_mig+(nbgroupsloci_theta*numpopsizeparams)][ei][mip]; // VS gp
    if (fc == 0 && cc == 0 && fm > 0)
    {
      temp1 = LOG2 - (mc * log (fm)) - hc - qintg - mintg;
      temp2 = log (exp (uppergamma (mc, 2 * fm * x / qmax)) - exp (uppergamma (mc, mmax * fm)));
    }
    else
    {
      if (fm == 0 && mc == 0 && fc > 0)
      {
        temp1 = LOG2 - (cc * log (fc)) - hc - qintg - mintg;
        temp2 = log (exp (uppergamma (cc, 2 * fc / qmax)) - exp (uppergamma (cc, fc * mmax / x)));
      }
      else
      {
        if (fc == 0 && cc == 0 && mc == 0 && fm == 0)
        {
          temp1 = log (2 * log (mmax * qmax / (2 * x))) - hc - qintg - mintg;;
          temp2 = 0;
        }
        else
        {
          temp1 = LOG2 + (mc * log (x)) - ((cc + mc) * log (fc + fm * x)) - hc - qintg - mintg;
          a = uppergamma (cc + mc, 2 * (fc + fm * x) / qmax);
          b = uppergamma (cc + mc, mmax * (fm + fc / x));
          /* the difference between two upppergamma()s is the negative of */
          /* the difference between two lowergamma()s.  However numbers   */
          /* for a or b often come out to be very large and equal under   */
          /* under uppergamma  or underlowergramma, in which case it is   */
          /* necessary to trap this and use the other function            */
          if (a==b)
          {
            b = lowergamma (cc + mc, 2 * (fc + fm * x) / qmax);
            a = lowergamma (cc + mc, mmax * (fm + fc / x));
          }
		  if (a>b)
          {
			LogDiff(temp2,a,b);
          }
		  else
          {
            temp1 = temp2 = 0.0;
          }

          //temp2 = log (exp (uppergamma (cc + mc, 2 * (fc + fm * x) / qmax)) - exp (uppergamma (cc + mc, mmax * (fm + fc / x))));
        }
      }
    }
    if ((temp1 + temp2 < 700 ) && (temp1 + temp2 > -700 )) // skip things that cannot be exped
      sum += exp (temp1 + temp2);
  }
  sum /= (lasttree - firsttree + (firsttree == 0)); 
  return -sum; 
}                               /* marginpopmig */


// marginpop_expomig
// calculate the probability density of the product 2NM when migration has an exponential prior.
// VS changed this function to include the groups of loci.
// INPUT
//	int mi: index of migration parameter, e.g. mi=0 for m0>1, mi=1 for m0>2 ...  
//	int gp_mig: index of group for migration parameter
//	int gp_theta: index of group for theta parameter
//	int firsttree: index of the first tree that will be used to evaluate this function
//	int lasttree: index of the last tree that will be used to evaluate this function
//	double x: ??
//	int thetai: is this the index of theta parameter
// NOTE: the index gp_mig is different in imig and gsampinf structures.
// This is why in the gsampinf we add gp_mig+(nbgroupsloci_theta*numpopsizeparams) to obtain the correct index
#define OCUTOFF  10
double
marginpop_expomig (int mi, int gp_mig, int gp_theta, int firsttree, int lasttree, double x, int thetai)
{
  int ei, cc, mc;
  int ccp, fcp, hcp, mcp, fmp, qip, mip;
  double sum, temp1, temp2, temp3,fc, fm, hc, qintg, mintg, mmean, qmax;
  int zadj, maxz = -1000000000;
  double acumm = 0;
  double max, min;

  //int gp_mig; // VS gp
  //int gp_theta; // VS these two values should be given as input to the function

  min = 0;
  max = EXPOMIGPLOTSCALE * imig[mi].pr[gp_mig].mean; // VS gp
  if (x < min || x > max)
    return OFFSCALEVAL;

  ccp = gsamp_ccp + thetai;
  fcp = gsamp_fcp + thetai;
  hcp = gsamp_hccp + thetai;
  mcp = gsamp_mcp + mi;
  fmp = gsamp_fmp + mi;
  qip = gsamp_qip + thetai;
  mip = gsamp_mip + mi;
  mmean = imig[mi].pr[gp_mig].mean; // VS gp
  qmax = itheta[thetai].pr[gp_theta].max; // VS gp

  for (sum = 0, ei = firsttree; ei < lasttree; ei++)

  {
    cc = (int) gsampinf[gp_theta][ei][ccp]; // VS gp
    fc = (double) gsampinf[gp_theta][ei][fcp]; // VS gp
    hc = (double) gsampinf[gp_theta][ei][hcp]; // VS gp
    mc = (int) gsampinf[gp_mig+(nbgroupsloci_theta*numpopsizeparams)][ei][mcp]; // VS gp
    fm = (double) gsampinf[gp_mig+(nbgroupsloci_theta*numpopsizeparams)][ei][fmp]; // VS gp
    qintg = (double) gsampinf[gp_theta][ei][qip]; // VS gp
    mintg = (double) gsampinf[gp_mig+(nbgroupsloci_theta*numpopsizeparams)][ei][mip]; // VS gp
    temp1 = x + fc *mmean + fm * mmean * x;
    temp2 = LOG2 - hc - log(mmean) - qintg - mintg;
    temp3 = 2 * temp1/(mmean * qmax);

    if (mc == 0 && cc == 0 )
    {
      temp3 = uppergamma(0,temp3);
      temp2 += temp3;
    }
    else
    {
      if (cc == 0) // mc > 0 
      {
        temp3 = uppergamma(mc,temp3);
        temp2 += temp3 - mc* log(fm + 1/mmean);
      }
      else
      {
        if (mc == 0) //cc > 0
        {
          temp3 = uppergamma(cc,temp3);
          temp2 += temp3 + -cc *log(fc + x*(fm + 1/mmean));; 
        }
        else
        {
          temp3 = uppergamma(mc+cc,temp3);
          temp2 += temp3 - cc*log(x) + (cc+mc)*log(mmean*x/temp1); 
        }
      }
    }
    eexp (temp2, &eexpsum[ei].m, &eexpsum[ei].z);
    if (eexpsum[ei].z > maxz)
      maxz = eexpsum[ei].z;
  }
  for (ei = firsttree; ei < lasttree; ei++)
  {
    zadj = eexpsum[ei].z - (maxz - OCUTOFF);
    eexpsum[ei].m *= pow (10.0, (double) zadj);
    acumm += eexpsum[ei].m;
  }
  sum = log (acumm) + (maxz - OCUTOFF) * LOG10;
  sum -= log(lasttree - firsttree + (firsttree == 0));
  sum = exp(sum);
  return -sum;
}                               //marginpop_expomig
#undef OCUTOFF  
#undef OFFSCALEVAL


// marginalopt_popmig
// find the peak of the marginal density for a 2NM term.  
// This is similar to the marginalopt() function, but this calls marginpopmig()
// VS changed this function to allow for groups of loci
// INPUT
//	int firsttree: index of the first tree that will be used to find the maximum
//	int lasttree: index of the last tree that will be used to find the maximum 
//	double *mlval: maximum likelihood value (OUTPUT)
//	double *peakloc: parameter values that correspond to the maximum posterior probability (OUTPUT)
//	int *mpop: index of theta parameter in itheta structure 
//	int *mterm: index of migration parameter in imig structure 
// NOTE
//	int gp_theta: index of group for theta
//	int gp_mig: index of group for mig - NOTE: gp index does not start at zero for the migration parameters in the gsampinf matrix
//				but it starts at zero for the imig structure!!
// NOT sure if gp_theta and gp_mig should be given as input, or if the for loop should be inside this function
#define SEARCHSTARTFRAC  4      // fraction of position in parameter range to start at
#define POPMIG_MARGINBINS 100
void
marginalopt_popmig (int firsttree, int lasttree, double *mlval, double *peakloc, int *mpop, int *mterm)
{
  int i,j;
  double ftol, ax, bx, cx, fa, xmax, ml;
  double upperbound;
  double maxf;
  int maxj;
  //double (*func) (int, int, int, double, int); // VS changed this function
  double (*func) (int, int, int, int, int, double, int); // VS changed this function to include gp_mig and gp_theta
  
  
  int gp_theta; // VS gp_theta - in this case we are trying to find 2NM we need to give as input the group pairs (gtheta, gmig)
  int gp_mig; // VS gp_mig - in this case we are trying to find 2NM we need to give as input the group pais (gtheta, gmig)

  int nbparams, nbthetaparams, nbmigparams; // VS

  // VS compute some variables that will simplify the code in the for loops and allocation of variables
  nbparams=(numpopsizeparams*nbgroupsloci_theta)+(nummigrateparams*nbgroupsloci_mig); 
  nbthetaparams=numpopsizeparams*nbgroupsloci_theta;
  nbmigparams=nummigrateparams*nbgroupsloci_mig;
 
  ftol = 1e-7;

	
  if (modeloptions[EXPOMIGRATIONPRIOR])
    func = marginpop_expomig;
  else
    func = marginpopmig;

  // If there are two groups of theta and 2 groups of mig
  // the assumption is that the all the loci share the same
  // theta and mig groups. That means that the assignment array
  // for theta and mig are the same.



  // VS for loop through all the mig parameters
  for (i = 0; i < nbmigparams; i++) // VS nbmigparams
  {
	
    gp_mig = g_gp_ip[mterm[i]]; // VS group for migration is the group of index i corresponding to mterm
	gp_theta = g_gp_ip[mpop[i]]; // VS group for theta - not sure about this

    // VS again there is the issue of mapping i into groups
    // Note that in the itheta and imig the index of param goes from 0 to numpopsizeparam for theta
	// and 0 to nummigrateparam for mig. Hence need to get the correct index for mterm and mpop.
	// This is obtained by i-(gp_mig*nbgroupsloci_mig)
    if (modeloptions[EXPOMIGRATIONPRIOR])
      upperbound = EXPOMIGPLOTSCALE * imig[mterm[i]].pr[gp_mig].mean; // VS gp_mig, index i needs to be corrected depending on the group
    else
      upperbound = itheta[mpop[i-(gp_mig*nbgroupsloci_mig)]].pr[gp_theta].max * imig[mterm[i-(gp_mig*nbgroupsloci_mig)]].pr[gp_mig].max/2.0; // VS gp
/* 8/31/09  JH
  it turns out these curves often have multiple peaks, and mnbrak just was not working then.
  replace use of mnbrak with function calls over POPMIG_MARGINBINS to find an interval with 
  a peak.  This is quite crude */
    ax = MINPARAMVAL;
    maxf = 1e100;
    maxj = -1;
    for (j = 0;j< POPMIG_MARGINBINS;j++)
    {
      ax = MINPARAMVAL + j*(upperbound /(double) POPMIG_MARGINBINS);
      fa = func (mterm[i], gp_mig, gp_theta, firsttree, lasttree, ax, mpop[i]); // VS gp
//printf("j  %d  ax  %lf  fa  %lf \n",j, ax, fa);
      if (fa < maxf)
      {
        maxf = fa;
        maxj = j;
      }
    }
//printf("i %d  maxf  %lf   maxj  %d\n",i,maxf, maxj);
    if (maxj == 0)
    {
      ax  =  MINPARAMVAL;
      cx = MINPARAMVAL + upperbound /(double) POPMIG_MARGINBINS;
      bx = (ax + cx)/2.0;
    }
    else
    {
      if (maxj == POPMIG_MARGINBINS -1)
      {
        ax = MINPARAMVAL + ((double) POPMIG_MARGINBINS - 1)*(upperbound /(double) POPMIG_MARGINBINS);
        bx = MINPARAMVAL + ((double) POPMIG_MARGINBINS - 2)*(upperbound /(double) POPMIG_MARGINBINS);
        cx = upperbound;
      }
      else
      {
        bx = MINPARAMVAL + ((double) maxj - 1)*(upperbound /(double) POPMIG_MARGINBINS);
        cx = MINPARAMVAL + ((double) maxj + 1)*(upperbound /(double) POPMIG_MARGINBINS);
        ax = MINPARAMVAL + ((double) maxj)*(upperbound /(double) POPMIG_MARGINBINS);
      }
    }
//printf("ax %lf  bx %lf  cx %lf\n",ax,bx,cx);
    //{
      //VS ml = -goldenmod (mterm[i], firsttree, lasttree, ax, bx, cx, ftol, &xmax, func,mpop[i]);
	  // VS changed this function because we need to give as input the groups for theta and migration
	  ml = -goldenmod_2nm(mterm[i], gp_mig, gp_theta, firsttree, lasttree, ax, bx, cx, ftol, &xmax, func,mpop[i]); // VS gp
      mlval[i] = ml;
      peakloc[i] = xmax;
    //}
//printf("ml  %lf  xmax %lf\n",ml, xmax);
  }
}                               /* marginalopt_popmig */
