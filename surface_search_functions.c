/* MigSelect 2013-2019 Arun Sethuraman, Vitor Sousa, Jody Hey */
/* This code was developed based on IMa2  2009-2010  Jody Hey, Rasmus Nielsen and Sang Chul Choi*/

/* this file is numerical recipes stuff with some small modifications */
#undef GLOBVARS
#include "imamp.h"              // may not really be needed in this particular file

/*********** LOCAL STUFF **********/

#define NRANSI
#define NR_END 1
#define FREE_ARG char*

/* local function prototypes */

static double sign (double a, double b);

/* local functions */


static double
sign (double a, double b)
{
  if (b > 0.0)
    return fabs (a);

  else
    return (-fabs (a));
}


/************ GLOBAL FUNCTIONS ******************/

/* NR stuff for finding minimum using golden section */

// mnbrakmod
// function that will find the minimum of a function using a certain algorithm
// VS modified this function such that it gets as input the group index
// INPUT
//	int ndim: index in the array of parameters (where all the groups are joined together in a contiguous array)
//  int paramindex: index of the parameters in the iparam structure
//	int gp: index of the group for iparam structure, not for gsampinf, i.e. 0<gp<nbgroupsloci_theta or 0<gp<nbgroupsloci_mig
//			NOTE: For gsampinf matrix the index of gp is different. For gsamping we have:
//				gp for theta from 0 to nbgroupsloci_theta-1
//				gp for mig from nbgroupsloci_theta to (nbgroupsloci_theta+nbgroupsloci_mig)-1
//	int firsttree: index of the first tree that will be used to find the maximum
//	int lasttree: index of the last tree that will be used to find the maximum
//	double *ax: ??
//	double *bx: ??
//  double *cx: ??
//	double *fa: ??
//	double *fb: ??
//	double *fc: ??
//	double (*func) (int, int, int, int, double, int): pointer to the function that will be used to find the peak
//	int ifneeded: ??
#define gold            1.618034
#define glimit          100.0
#define tiny            1.0e-20
#define r               0.61803399
/* modified mnbrak() and golden() from NR code.  used only for the one dimensional optimization */
// VS modified this function such that the group index gp is given as input
void
mnbrakmod (int ndim, int firsttree, int lasttree, double *ax, double *bx,
           double *cx, double *fa, double *fb, double *fc,
           double (*func) (int, int, int, double, int), int ifneeded)
{

  /* Programs using routine MNBRAK must supply an EXTERNAL
     FUNCTION func(x:REAL):REAL FOR which a minimum is TO be found */
  /* I converted this to a function and included code to check for
     a failure to find a braket */
  double ulim, u, rr, q, fu, dum;
  *fa = func (ndim, firsttree, lasttree, *ax, ifneeded); 
  *fb = func (ndim, firsttree, lasttree, *bx, ifneeded); 
  if (*fb > *fa)
  {
    dum = *ax;
    *ax = *bx;
    *bx = dum;
    dum = *fb;
    *fb = *fa;
    *fa = dum;
  }
  *cx = fabs (*bx + gold * (*bx - *ax));
  *fc = func (ndim, firsttree, lasttree, *cx, ifneeded); 
_L1:
  /* need to add some excapes for floating point problems */
  if (*fb >= *fc && *fb > -DBL_MAX && !(*fb == 0 && *fc == 0))
  {
    rr = (*bx - *ax) * (*fb - *fc);
    q = (*bx - *cx) * (*fb - *fa);
    u = (double) *bx - ((*bx - *cx) * q - (*bx - *ax) * rr) /
      (2.0 * sign (FMAX (fabs (q - rr), (float) tiny), q - rr));
    ulim = *bx + glimit * (*cx - *bx);
    if ((*bx - u) * (u - *cx) > 0.0)
    {
      fu = func (ndim,firsttree, lasttree, u, ifneeded); 
      if (fu < *fc)
      {
        *ax = *bx;
        *fa = *fb;
        *bx = u;
        *fb = fu;
        goto _L1;
      }
      if (fu > *fb)
      {
        *cx = u;
        *fc = fu;
        goto _L1;
      }
      u = *cx + gold * (*cx - *bx);
      fu = func (ndim, firsttree, lasttree, u, ifneeded);
    }
    else if ((*cx - u) * (u - ulim) > 0.0)
    {
      fu = func (ndim, firsttree, lasttree, u, ifneeded); 
      if (fu < *fc)
      {
        *bx = *cx;
        *cx = u;
        u = *cx + gold * (*cx - *bx);
        *fb = *fc;
        *fc = fu;
        fu = func (ndim, firsttree, lasttree, u, ifneeded);
      }
    }
    else if ((u - ulim) * (ulim - *cx) >= 0.0)
    {
      u = ulim;
      fu = func (ndim, firsttree, lasttree, u, ifneeded); 
    }
    else
    {
      u = *cx + gold * (*cx - *bx);
      fu = func (ndim, firsttree, lasttree, u, ifneeded); 
    }
    *ax = *bx;
    *bx = *cx;
    *cx = u;
    *fa = *fb;
    *fb = *fc;
    *fc = fu;
    goto _L1;
  }
  if (*fb < -DBL_MAX)
    printf (" minus infinity in mnbrakmod \n");
}
#undef gold
#undef glimit
#undef tiny

// goldenmod
// this functions finds the peak, or the mode of the posterior when it is close to the margins??
// VS - modified this function such that the index of the group is given as input
// also modified the number of arguments of function *func
// INPUT:
//	int ndim: index of parameter in array of parameters (where all groups of loci are put together in a contiguous array)
//	int paramindex: index of parameter for iparam structure (in this case, for migrate paramindex varies between 0 and number of migration paramaters)
//	int gp: index of the group of loci for iparam structure
//		NOTE: For gsampinf matrix the index of gp is different. For gsamping we have:
//			  gp for theta from 0 to nbgroupsloci_theta-1
//			  gp for mig from nbgroupsloci_theta to (nbgroupsloci_theta+nbgroupsloci_mig)-1
//	int firsttree: index of the first tree that will be used to find the maximum 
//	int lasttree: index of the last tree
//	double ax:
//	double bx:
//	double cx:
//	double tol_:
//	double *xmin:
//	double (*func) (int, int, int, int, double, int): pointer to the function that will be used to evaluate the values
//	int ifneeded:
//goldenmod (int ndim, int firsttree, int lasttree, double ax, double bx,
//           double cx, double tol_, double *xmin, double (*func) (int, int,
//                                                                 int, double, int), int ifneeded)
double
goldenmod (int ndim, int firsttree, int lasttree, double ax, double bx,
           double cx, double tol_, double *xmin, double (*func) (int, int, int, double, int), int ifneeded)
{

  /* Programs using routine GOLDEN must supply an EXTERNAL
     FUNCTION func(x:REAL):REAL whose minimum is TO be found. */
  double Result, f1, f2, c, x0, x1, x2, x3;
  c = 1.0 - r;
  x0 = ax;
  x3 = cx;
  if (fabs (cx - bx) > fabs (bx - ax))
  {
    x1 = bx;
    x2 = bx + c * (cx - bx);
  }
  else
  {
    x2 = bx;
    x1 = bx - c * (bx - ax);
  }
  f1 = func (ndim, firsttree, lasttree, x1, ifneeded); 
  f2 = func (ndim, firsttree, lasttree, x2, ifneeded); 
  while (fabs (x3 - x0) > tol_ * (fabs (x1) + fabs (x2)) && (x3 > -DBL_MAX))
  {
    if (f2 < f1)
    {
      x0 = x1;
      x1 = x2;
      x2 = r * x1 + c * x3;
      f1 = f2;
      f2 = func (ndim, firsttree, lasttree, x2, ifneeded); 
      continue;
    }
    x3 = x2;
    x2 = x1;
    x1 = r * x2 + c * x0;
    f2 = f1;
    f1 = func (ndim, firsttree, lasttree, x1, ifneeded);
  }
  if (f1 < f2)
  {
    Result = f1;
    *xmin = x1;
  }
  else
  {
    Result = f2;
    *xmin = x2;
  }
  if (x3 < -DBL_MAX)
    printf (" minus infinity in goldenmod \n");
  return Result;
}




// goldenmod_2nm
// this functions finds the peak, or the mode of the posterior when it is close to the margins??
// VS - modified this function such that the index of the group is given as input
// also modified the number of arguments of function *func
// VS needed to create a new function specific to 2nm because in that case
// we need to give as input functions with a different number of arguments
// due to the fact that we need to take into account both the gp_theta and gp_mig.
// The reason is that to compute the 2NM we need information both from theta and mig parameters.
// INPUT:
//	int ndim: index of parameter in array of parameters (where all groups of loci are put together in a contiguous array)
//	int gp_mig: index of the group of loci for mig
//	int gp_theta: index of the group of loci for theta
//	int firsttree: index of the first tree that will be used to find the maximum 
//	int lasttree: index of the last tree
//	double ax:
//	double bx:
//	double cx:
//	double tol_:
//	double *xmin:
//	double (*func) (int, int, int, int, double, int): pointer to the function that will be used to evaluate the values
//	int ifneeded:
// NOTE: the index gp_mig starts at zero for imig structure, but not for gsampinf.
//		this is why when accessing gsampinf we need gp_mig+(nbgroupsloci_theta*numpopsizeparams)
//goldenmod (int ndim, int firsttree, int lasttree, double ax, double bx,
//           double cx, double tol_, double *xmin, double (*func) (int, int,
//                                                                 int, double, int), int ifneeded)
double
goldenmod_2nm (int ndim, int gp_mig, int gp_theta, int firsttree, int lasttree, double ax, double bx,
           double cx, double tol_, double *xmin, double (*func) (int, int, int, int, int, double, int), int ifneeded)
{

  /* Programs using routine GOLDEN must supply an EXTERNAL
     FUNCTION func(x:REAL):REAL whose minimum is TO be found. */
  double Result, f1, f2, c, x0, x1, x2, x3;
  c = 1.0 - r;
  x0 = ax;
  x3 = cx;
  if (fabs (cx - bx) > fabs (bx - ax))
  {
    x1 = bx;
    x2 = bx + c * (cx - bx);
  }
  else
  {
    x2 = bx;
    x1 = bx - c * (bx - ax);
  }
  f1 = func (ndim, gp_mig, gp_theta, firsttree, lasttree, x1, ifneeded); // VS gp
  f2 = func (ndim, gp_mig, gp_theta, firsttree, lasttree, x2, ifneeded); // VS gp
  while (fabs (x3 - x0) > tol_ * (fabs (x1) + fabs (x2)) && (x3 > -DBL_MAX))
  {
    if (f2 < f1)
    {
      x0 = x1;
      x1 = x2;
      x2 = r * x1 + c * x3;
      f1 = f2;
      f2 = func (ndim, gp_mig, gp_theta, firsttree, lasttree, x2, ifneeded); // VS gp
      continue;
    }
    x3 = x2;
    x2 = x1;
    x1 = r * x2 + c * x0;
    f2 = f1;
    f1 = func (ndim, gp_mig, gp_theta, firsttree, lasttree, x1, ifneeded); // VS gp
  }
  if (f1 < f2)
  {
    Result = f1;
    *xmin = x1;
  }
  else
  {
    Result = f2;
    *xmin = x2;
  }
  if (x3 < -DBL_MAX)
    printf (" minus infinity in goldenmod \n");
  return Result;
}
