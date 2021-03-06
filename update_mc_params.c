/* MigSelect 2013-2019 Arun Sethuraman, Vitor Sousa, Jody Hey */
/* This code was developed based on IMa2  2009-2010  Jody Hey, Rasmus Nielsen and Sang Chul Choi*/

#undef GLOBVARS
#include "imamp.h"
#include "updateassignment.h"

/****** global - used mostly here but setup in initialize.c *********/
int urri[2 * MAXLOCI][2 * MAXLOCI];     // double it - some loci might have multiple rates
double urrlow[2 * MAXLOCI][2 * MAXLOCI], urrhi[2 * MAXLOCI][2 * MAXLOCI];

/*********** LOCAL STUFF **********/

/*** GLOBAL FUNCTIONS ***/


#define MAXURCHECK 100          /* not sure how difficult it will be to find random ratios that fit the priors */

/* changeu picks new values on a log scale over a very wide range
the recorded range is set in IMn1.c to be 1/1000 -1000  (i.e. over a
1,000,000 fold range.  However the actual range is set to be much greater
This effectively sets the range to be infinite, and causes mutation rates to be picked without a prior.
By setting the actual range to be very wide, the prior probability within the range of 1/1000 -1000 becomes flat */
int
changeu (int ci, int j, int *k)
{
  /* update the ratio between two mutation rate scalars .  The upate is drawn from a uniform log scale between 1/maxratio and maxratio 
     maxratio is set to be 3 times the maximum value of the individual scalars */
  double U, newr, metropolishastingsterm, olduj, olduk;
  double newpdg[2], newkappa[2], likenewsum;
  double d, r;
  int i, li, lj, lk, ai, aj, ak;
  static int start = 0;
  static double windowsize, maxratio;
  double tempr, temp;
  int urcheck, counturcheck, uindex;

#ifdef _DEBUG
  double prodcheck;
#endif /* _DEBUG */

  if (assignmentoptions[POPULATIONASSIGNMENTCHECKPOINT] == 1)
  {
    assertgenealogy (ci);
  }

  /* windowsize and maxratio are on log scales 
     all mutation rate scalars have the same priors and windowsize */
  if (start == 0)
  {
    maxratio = 3 * L[0].u_rec[0].pr.max;
    if (calcoptions[MUTATIONPRIORRANGE])
      windowsize = L[0].u_rec[0].pr.max;
    else
      windowsize = L[0].u_rec[0].win;
    start = 1;
  }

#ifdef _DEBUG
  for (prodcheck = 1, ai = 0; ai < nurates; ai++)
  {
    prodcheck *= C[ci]->G[ul[ai].l].uvals[ul[ai].u];
  }
  assert (prodcheck < 1.00001 && prodcheck > 0.99999);

#endif /* _DEBUG */

  if (nurates > 2)
  {
    do
    {
      *k = (int) (uniform () * nurates);
    }
    while (*k == j || *k < 0 || *k >= nurates);
  }
  else
  {
    assert (j == 0);
    *k = 1;
  }

  lj = ul[j].l;
  aj = ul[j].u;
  lk = ul[*k].l;
  ak = ul[*k].u;
  assert (*k != j);
  olduj = C[ci]->G[lj].uvals[aj];
  olduk = C[ci]->G[lk].uvals[ak];
  r = log (olduj / olduk);
  if (calcoptions[MUTATIONPRIORRANGE])
  {
    counturcheck = 0;

    do
    {
      if (counturcheck >= MAXURCHECK)
      {
        if (assignmentoptions[POPULATIONASSIGNMENTCHECKPOINT] == 1)
        {
          assertgenealogy (ci);
        }
        return -1;
      }
      U = uniform ();

      if (U > 0.5)
        newr = r + (2.0 * U - 1.0) * windowsize;
      else
        newr = r - windowsize * U * 2.0;

      if (newr > maxratio)
        newr = 2.0 * maxratio - newr;
      else if (newr < -maxratio)
        newr = 2.0 * (-maxratio) - newr;

      d = exp ((newr - r) / 2);

      /* urri[i][j] has a 0 if neither i nor j has a prior. 1 if i has a prior and j does not,  -1 if i does not have a pior and j does  */
      /* must check if the new rates will cause the ratios among scalars to fall outside of the allowed ranges of ratios  */
      /* if they do, then they are rejected */
      urcheck = 0;
      counturcheck++;
      if (urri[j][*k] == 2)     // both scalars have priors 
      {
        urcheck = (newr <= urrlow[j][*k] || newr >= urrhi[j][*k]);
        for (uindex = 0; uindex < nurates; uindex++)
        {
          if (urri[j][uindex] == 2 && uindex != *k)
          {
            tempr =
              log (d * C[ci]->G[lj].uvals[aj] /
                   C[ci]->G[ul[uindex].l].uvals[ul[uindex].u]);
            urcheck = urcheck || (tempr <= urrlow[j][uindex]
                                  || tempr >= urrhi[j][uindex]);
          }
          if (urcheck)
            break;
        }
        for (uindex = 0; uindex < nurates; uindex++)
        {
          if (urri[uindex][*k] == 2 && uindex != j)
          {
            tempr =
              log (C[ci]->G[ul[uindex].l].uvals[ul[uindex].u] /
                   (C[ci]->G[lk].uvals[ak] / d));
            urcheck = urcheck || (tempr <= urrlow[uindex][*k]
                                  || tempr >= urrhi[uindex][*k]);
          }
          if (urcheck)
            break;
        }
      }
      else
      {
        if (urri[j][*k] == 1)   // the uj scalar has a prior, but not the uk scalar
        {
          for (uindex = 0; uindex < nurates; uindex++)
          {
            if (urri[j][uindex] == 2)
            {
              tempr =
                log (d * C[ci]->G[lj].uvals[aj] /
                     C[ci]->G[ul[uindex].l].uvals[ul[uindex].u]);
              urcheck = urcheck || (tempr <= urrlow[j][uindex]
                                    || tempr >= urrhi[j][uindex]);
            }
            if (urcheck)
              break;
          }
        }
        if (urri[j][*k] == -1)  // the uk scalar has a prior, but not the uj scalar
        {
          for (uindex = 0; uindex < nurates; uindex++)
          {
            if (urri[uindex][*k] == 2)
            {
              tempr =
                log (C[ci]->G[ul[uindex].l].uvals[ul[uindex].u] /
                     (C[ci]->G[lk].uvals[ak] / d));
              urcheck = urcheck || (tempr <= urrlow[uindex][*k]
                                    || tempr >= urrhi[uindex][*k]);
            }
            if (urcheck)
              break;
          }
        }
      }
    }
    while (urcheck);
  }
  else
  {
    U = uniform ();
    if (U > 0.5)
      newr = r + (2.0 * U - 1.0) * windowsize;
    else
      newr = r - windowsize * U * 2.0;

    if (newr > maxratio)
      newr = 2.0 * maxratio - newr;
    else if (newr < -maxratio)
      newr = 2.0 * (-maxratio) - newr;

    d = exp ((newr - r) / 2);
  }
  C[ci]->G[lj].uvals[aj] *= d;
  C[ci]->G[lk].uvals[ak] /= d;
  likenewsum = 0;
  for (i = 0; i < 2; i++)
  {
    if (i == 0)
    {
      li = lj;
      ai = aj;
    }
    else
    {
      li = lk;
      ai = ak;
    }
    switch (L[li].umodel[ai])
    {
    case HKY:
      if (assignmentoptions[JCMODEL] == 1)
      {
        newpdg[i] = likelihoodJC (ci, li, C[ci]->G[li].uvals[ai]);
      }
      else
      {
        U = uniform ();
        if (U > 0.5)
        {
          newkappa[i] =
            C[ci]->G[li].kappaval + (2.0 * U - 1.0) * L[li].kappa_rec->win;
          if (newkappa[i] > L[li].kappa_rec->pr.max)
            newkappa[i] = 2.0 * L[li].kappa_rec->pr.max - newkappa[i];
        }
        else
        {
          newkappa[i] = C[ci]->G[li].kappaval - L[li].kappa_rec->win * U * 2.0;
          if (newkappa[i] < 0)
            newkappa[i] = -newkappa[i];
        }
        if (ci == 0)
          L[li].kappa_rec->upinf->tries++;
        newpdg[i] =
          likelihoodHKY (ci, li, C[ci]->G[li].uvals[ai], newkappa[i],
                         -1, -1, -1, -1);
      }
      break;
    case INFINITESITES:
      newpdg[i] = likelihoodIS (ci, li, C[ci]->G[li].uvals[ai]);
      break;
    case STEPWISE:
      newpdg[i] = likelihoodSW (ci, li, ai, C[ci]->G[li].uvals[ai], 1.0);
      break;
    }
    likenewsum += newpdg[i] - C[ci]->G[li].pdg_a[ai];
  }

  metropolishastingsterm = exp (beta[ci] * gbeta * likenewsum);
  U = uniform ();
  if (U < DMIN(1.0, metropolishastingsterm))  //9/13/2010 
  //if (metropolishastingsterm >= 1.0 || metropolishastingsterm > U)
  {
    for (i = 0; i < 2; i++)
    {
      if (i == 0)
      {
        li = lj;
        ai = aj;
      }
      else
      {
        li = lk;
        ai = ak;
      }
      C[ci]->G[li].pdg += newpdg[i] - C[ci]->G[li].pdg_a[ai];
      C[ci]->G[li].pdg_a[ai] = newpdg[i];
      if (L[li].umodel[ai] == HKY)
      {
        assert (ai == 0);
        C[ci]->G[li].kappaval = newkappa[i];
        if (ci == 0)
          L[li].kappa_rec->upinf->accp++;
        copyfraclike (ci, li);
        storescalefactors (ci, li);
      }
    }
    C[ci]->allpcalc.pdg += likenewsum;

#ifdef _DEBUG
    for (prodcheck = 1, ai = 0; ai < nurates; ai++)
    {
      prodcheck *= C[ci]->G[ul[ai].l].uvals[ul[ai].u];
    }
    assert (prodcheck < 1.00001 && prodcheck > 0.99999);

#endif /*  */

    if (assignmentoptions[POPULATIONASSIGNMENTCHECKPOINT] == 1)
    {
      assertgenealogy (ci);
    }
    return 1;
  }
  else
  {

    C[ci]->G[lj].uvals[aj] = olduj;
    C[ci]->G[lk].uvals[ak] = olduk;
    if (L[lj].umodel[aj] == HKY)
      restorescalefactors (ci, lj);
    if (L[lk].umodel[ak] == HKY)
      restorescalefactors (ci, lk);
    if (L[lj].umodel[aj] == STEPWISE)
    {
      temp = likelihoodSW (ci, lj, aj, olduj, 1.0);
      C[ci]->G[lj].pdg += temp - C[ci]->G[lj].pdg_a[aj];
      C[ci]->G[lj].pdg_a[aj] = temp;
    }
    if (L[lk].umodel[ak] == STEPWISE)
    {
      temp = likelihoodSW (ci, lk, ak, olduk, 1.0);
      C[ci]->G[lk].pdg += temp - C[ci]->G[lk].pdg_a[ak];
      C[ci]->G[lk].pdg_a[ak] = temp;
    }

#ifdef _DEBUG
    for (prodcheck = 1, ai = 0; ai < nurates; ai++)
    {
      prodcheck *= C[ci]->G[ul[ai].l].uvals[ul[ai].u];
    }
    assert (prodcheck < 1.00001 && prodcheck > 0.99999);

#endif /*  */
    if (assignmentoptions[POPULATIONASSIGNMENTCHECKPOINT] == 1)
    {
      assertgenealogy (ci);
    }
    return 0;
  }

#ifdef _DEBUG
  for (i = 0, temp = 0; i < nloci; i++)
    for (ai = 0; ai < L[i].nlinked; ai++)
      temp += C[0]->G[i].pdg_a[ai];

#endif /*  */
  assert (fabs (temp - C[0]->allpcalc.pdg) < 1e-10); //AS: not sure if this has to change or not...Tue Feb 23 15:28:19 EST 2016
}                               /* changeu */


/* only use for HKY model with single locus data sets, otherwise kappa is updated in changeu */
int
changekappa (int ci)
{
  int li;
  double U, metropolishastingsterm, newpdg, newkappa;
  li = 0;
  U = uniform ();
  L[li].kappa_rec->upinf->tries++;
  if (U > 0.5)
  {
    newkappa = C[ci]->G[li].kappaval + (2.0 * U - 1.0) * L[li].kappa_rec->win;
    if (newkappa > L[li].kappa_rec->pr.max)
      newkappa = 2.0 * L[li].kappa_rec->pr.max - newkappa;
  }
  else
  {
    newkappa = C[ci]->G[li].kappaval - L[li].kappa_rec->win * U * 2.0;
    if (newkappa < 0)
      newkappa = -newkappa;
  }
  newpdg =
    likelihoodHKY (ci, li, C[ci]->G[li].uvals[0], newkappa, -1, -1, -1, -1);
  metropolishastingsterm = newpdg - C[ci]->G[li].pdg_a[0];
  metropolishastingsterm = exp (beta[ci] * metropolishastingsterm);
  U = uniform ();
  if (metropolishastingsterm >= 1.0 || metropolishastingsterm > U)
  {
    C[ci]->allpcalc.pdg += newpdg - C[ci]->G[li].pdg;
    C[ci]->G[li].pdg = C[ci]->G[li].pdg_a[0] = newpdg;
    C[ci]->G[li].kappaval = newkappa;
    L[li].kappa_rec->upinf->accp++;
    copyfraclike (ci, li);
    storescalefactors (ci, li);
    return (1);
  }
  else
  {
    restorescalefactors (ci, li);
    return (0);
  }
}                               /* changekappa */

#undef MAXPRIORSCALETRY
