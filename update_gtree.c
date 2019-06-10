/* MigSelect 2013-2019 Arun Sethuraman, Vitor Sousa, Jody Hey */
/* This code was developed based on IMa2  2009-2010  Jody Hey, Rasmus Nielsen and Sang Chul Choi*/

#undef GLOBVARS
#include "imamp.h"
#include "update_gtree_common.h"
#include "update_gtree.h"
#include "updateassignment.h"

/*********** LOCAL STUFF **********/

struct edgemiginfo oldedgemig;
struct edgemiginfo oldsismig;
struct edgemiginfo newedgemig;
struct edgemiginfo newsismig;

/* declarded in update_gtree_common.h  not used here 
struct edge *copyedge; 
extern int holddownA[MAXLINKED];
extern int medgedrop;
extern double lmedgedrop;
extern double holdsisdlikeA[MAXLINKED];
extern struct genealogy holdgtree;
*/

static int mrootdrop;
static double lmrootdrop;
static struct genealogy_weights holdgweight_updategenealogy;
static struct genealogy_weights holdallgweight_updategenealogy;
static struct probcalc holdallpcalc_updategenealogy;

// VS - "hold" variables save the gweights for group to which updated locus belongs to 
// (does not need to be an array of groups, because at each time only one locus is updated)
static struct genealogy_weights holdgweight_theta_updategenealogy;
static struct genealogy_weights holdgweight_mig_updategenealogy;
static struct probcalc holdpcalc_theta_updategenealogy;
static struct probcalc holdpcalc_mig_updategenealogy;

// VS - hold variables to save the gweights for groups of loci
// when updating the assignment vectors for the loci that belong to different groups of loci
// This is not related with the assignment of individuals to populations (assignment.c functions)
static struct genealogy_weights holdgweight_update_assignloci[MAXGROUPS];
static struct genealogy_weights holdgweight_theta_update_assignloci[MAXGROUPS];
static struct probcalc holdpcalc_update_assignloci[MAXGROUPS];
static struct probcalc holdpcalc_theta_update_assignloci[MAXGROUPS];
static int holdgroupassign[MAXGROUPS];




int rootmove;                   /* used in update_gtree.c and update_gtree_covar.c */ // 12/16/09 JH note that update_gree_covar.c no longer used
static double mlist[2 * ABSMIGMAX];     // used in mwork_single_edge() mlist very wasteful of space - could use dynamic memory and checkmig() 
static int mplist[2 * ABSMIGMAX];

/* find the time when two populations join */
double
findjointime (int ci, int slidepop, int sispop, double edgeuptime,
              double sisuptime)
{
  int edgeperiod, sisperiod;
  double jointime;
  edgeperiod = findperiod (ci, edgeuptime);
  sisperiod = findperiod (ci, sisuptime);
  while (edgeperiod < sisperiod)

  {
    edgeperiod++;
    if (slidepop == C[ci]->droppops[edgeperiod][0]
        || slidepop == C[ci]->droppops[edgeperiod][1])
      slidepop = C[ci]->poptree[slidepop].down;
  }

  while (sisperiod < edgeperiod)
  {
    sisperiod++;
    if (sispop == C[ci]->droppops[sisperiod][0]
        || sispop == C[ci]->droppops[sisperiod][1])
      sispop = C[ci]->poptree[sispop].down;
  }

  // at this point edgeperiod == sisperiod 
  while (slidepop != sispop)
  {
    edgeperiod++;
    if (slidepop == C[ci]->droppops[edgeperiod][0]
        || slidepop == C[ci]->droppops[edgeperiod][1])
      slidepop = C[ci]->poptree[slidepop].down;
    if (sispop == C[ci]->droppops[edgeperiod][0]
        || sispop == C[ci]->droppops[edgeperiod][1])
      sispop = C[ci]->poptree[sispop].down;
  }

  if (edgeperiod == 0)
    jointime = 0;
  else
    jointime = C[ci]->tvals[edgeperiod - 1];
  return jointime;
}                               //findjointime

void
slider_nomigration (int ci, int li, int slidingedge, int *sis,
                    double *timepoint, double *slidedist)
/* this is very similar to slider, but with a an extra if/else on the slides up for the case of zero migration */
/* timepoint points at C[ci]->G[li].gtree[*slidingedge].time and is the current position of the sliding point, slidedist is the distance it must move 
with multiple populations,  and no migration, the sliding edge can only be in the population it started in, or an ancestral population
so the same goes for the point on the edge on which the slide is currently at. 
At the beginning of a slide the point is necessarily valid  - All down slides are valid for any distance
the upper limit at any point is the 
maximum of the top of the sliding edge and the time at which the sliding edge and the sister edge are in different populations  
if distance is negative  move up
	determine the upper limit ( MAX(top of edge, split time of sis and edge,  time of node of sis) )
	if upper limit is at a node,  pick left or right  at random
		call slider and continue up
	else reflect  (switch sign on remaining distance)
		call slider and and move down
else  move down
	if reach a node,  pick down or up at random
		if down,  call slider and continue down
		if up,  switch sign on remaining distance
			call slider and move up
kinds of upper limits
 top of sliding edge
 top of sister edge (a node) 
 beginning of period when the population of the edge and the sister branch come together into the same population 
*/
{
  double edgeuptime, sisuptime, popjointime;
  struct edge *gtree = C[ci]->G[li].gtree;
  int slidepop, sispop;
  if (*slidedist < 0)
  {
    /* go up */
    *slidedist = -*slidedist;
    if (slidingedge < L[li].numgenes)
      edgeuptime = 0;
    else
      edgeuptime = gtree[gtree[slidingedge].up[0]].time;

    if (*sis < L[li].numgenes)
      sisuptime = 0;
    else
      sisuptime = gtree[gtree[*sis].up[0]].time;

    assert (*timepoint >= edgeuptime);
    slidepop = gtree[slidingedge].pop;
    sispop = gtree[*sis].pop;
    if (slidepop != sispop)
      popjointime = findjointime (ci, slidepop, sispop, edgeuptime, sisuptime);
    else
      popjointime = 0;

    if (popjointime > edgeuptime && popjointime > sisuptime)
    {
      if (*slidedist < *timepoint - popjointime)
        /* slide up and stop,  sis remains the same */
      {
        *timepoint -= *slidedist;
        *slidedist = 0;
        assert (*timepoint > popjointime);
        return;
      }
      else
      {

        /* slide up and reflect, sis remains the same, leave slidedist positive so slidingedge goes down with next call to slider */
        *slidedist -= *timepoint - popjointime;
        *timepoint = popjointime;
        assert (*slidedist > 0);
        slider_nomigration (ci, li, slidingedge, sis, timepoint, slidedist);
        return;
      }
    }
    else
    {
      if (sisuptime == 0 || edgeuptime >= sisuptime)
      {
        if (*slidedist < *timepoint - edgeuptime)

          /* slide up and stop,  sis remains the same */
        {
          *timepoint -= *slidedist;
          *slidedist = 0;
          assert (*timepoint > edgeuptime);
          return;
        }
        else
        {
          /* slide up and reflect, sis remains the same, leave slidedist positive so slidingedge goes down with next call to slider */
          *slidedist -= *timepoint - edgeuptime;
          *timepoint = edgeuptime;
          assert (*slidedist > 0);
          slider_nomigration (ci, li, slidingedge, sis, timepoint, slidedist);
          return;
        }
      }
      else
      {
        /* edgeuptime is less than sis up time, and thus slidingedge can reach a node */
        if (*slidedist < *timepoint - sisuptime)

          /* slide up and stop,  sis remains the same */
        {
          *timepoint -= *slidedist;
          *slidedist = 0;
          assert (*timepoint > sisuptime);
          return;
        }
        else
        {

          /* slide up and reach a node, pick one side at random and recurse */
          *slidedist -= *timepoint - sisuptime;
          *timepoint = sisuptime;
          if (bitran () /*uniform() < 0.5 */ )
          {
            *sis = gtree[*sis].up[0];
          }
          else
          {
            *sis = gtree[*sis].up[1];
          }

          /* reset slidedist to negative, so slidingedge continues up the gtree in next call to slider */
          *slidedist = -*slidedist;
          assert (*slidedist < 0);
          slider_nomigration (ci, li, slidingedge, sis, timepoint, slidedist);
          return;
        }
      }
    }
  }
  else
  {
    /* go down */
    if (gtree[*sis].down == -1 || *timepoint + *slidedist < gtree[*sis].time)
    {

      /* if sis is the root, or distance is less than to next down node, just slide down */
      *timepoint += *slidedist;
      if (*timepoint >= TIMEMAX)
        *timepoint = TIMEMAX;
      *slidedist = 0;
      return;
    }
    else
    {

      /* a down node is reached */
      *slidedist -= gtree[*sis].time - *timepoint;
      *timepoint = gtree[*sis].time;
      if (bitran ())
      {
        /* begin to slide down the down node */
        *sis = gtree[*sis].down;
        slider_nomigration (ci, li, slidingedge, sis, timepoint, slidedist);
        return;
      }
      else
      {
        /* begin to slide up the sis  */
        if (gtree[gtree[*sis].down].up[0] == *sis)
        {
          *sis = gtree[gtree[*sis].down].up[1];
        }
        else
        {
          *sis = gtree[gtree[*sis].down].up[0];
        }
        *slidedist = -*slidedist;
        slider_nomigration (ci, li, slidingedge, sis, timepoint, slidedist);
        return;
      }
    }
  }
}                               /* slider_nomigration */

void
slider (int ci, int li, int slidingedge, int *sis, double *timepoint,
        double *slidedist)
/* this is not ready for case of multiple populations and zero migration */
/* timepoint points at C[ci]->G[li].gtree[*slidingedge].time and is the current position of the sliding point, 
slidedist is the distance it must move 
do not restructure the gtree. just figure out when and on which branch timepoint ends up on, 
this will be the new sisterbranch
use recursion */
// VS - program is crashing here with microsat data.
// The reason is a stack overflow, which may be caused by an infinite recursion in slider function.
{
  double uplimit;
  struct edge *gtree = C[ci]->G[li].gtree;
  if (*slidedist < 0)
  {
    /* go up */
    *slidedist = -*slidedist;
    if (slidingedge < L[li].numgenes)
      uplimit = 0;
    else
      uplimit = gtree[gtree[slidingedge].up[0]].time;
    assert (*timepoint >= uplimit);

    /* if uplimit >= sis up time  - slidingedge cannot reach a node */
    if (gtree[*sis].up[0] == -1 || uplimit >= gtree[gtree[*sis].up[0]].time)
    {
      if (*slidedist < *timepoint - uplimit)
        /* slide up and stop,  sis remains the same */
      {
        *timepoint -= *slidedist;
        *slidedist = 0;
        assert (*timepoint > uplimit);
        return;
      }
      else
      {
        /* slide up and reflect, sis remains the same, leave slidedist positive so slidingedge goes down with next call to slider */
        *slidedist -= *timepoint - uplimit;
        *timepoint = uplimit;
        assert (*slidedist > 0);
        slider (ci, li, slidingedge, sis, timepoint, slidedist);
        return;
      }
    }
    else
    {
      /* uplimit is less than sis up time, and thus slidingedge can reach a node */
      if (*slidedist < *timepoint - gtree[gtree[*sis].up[0]].time)
        /* slide up and stop,  sis remains the same */
      {
        *timepoint -= *slidedist;
        *slidedist = 0;
        assert (*timepoint > gtree[gtree[*sis].up[0]].time);
        return;
      }
      else
      {
        /* slide up and reach a node, pick one side at random and recurse */
        *slidedist -= *timepoint - gtree[gtree[*sis].up[0]].time;
        *timepoint = gtree[gtree[*sis].up[0]].time;
        if (bitran () /*uniform() < 0.5 */ )
        {
          *sis = gtree[*sis].up[0];
        }
        else
        {
          *sis = gtree[*sis].up[1];
        }

        /* reset slidedist to negative, so slidingedge continues up the gtree in next call to slider */
        *slidedist = -*slidedist;
        assert (*slidedist < 0);
        slider (ci, li, slidingedge, sis, timepoint, slidedist);
        return;
      }
    }
  }
  else
  {
    /* go down */
    if (gtree[*sis].down == -1 || *timepoint + *slidedist < gtree[*sis].time)
    {
      /* if sis is the root, or distance is less than to next down node, just slide down */
      *timepoint += *slidedist;
      if (*timepoint >= TIMEMAX)
        *timepoint = TIMEMAX;
      *slidedist = 0;
      return;
    }
    else
    {
      /* a down node is reached */
      *slidedist -= gtree[*sis].time - *timepoint;
      *timepoint = gtree[*sis].time;
      if (bitran ())
      {
        /* begin to slide down the down node */
        *sis = gtree[*sis].down;
        slider (ci, li, slidingedge, sis, timepoint, slidedist);
        return;
      }
      else
      {
        /* begin to slide up the sis  */
        if (gtree[gtree[*sis].down].up[0] == *sis)
        {
          *sis = gtree[gtree[*sis].down].up[1];
        }
        else
        {
          *sis = gtree[gtree[*sis].down].up[0];
        }
        *slidedist = -*slidedist;
        slider (ci, li, slidingedge, sis, timepoint, slidedist);
        return;
      }
    }
  }
}                               /* slider */

void
joinsisdown (int ci, int li, int sis, int *tmrcachange)
{

  /* extend sis, and XFREE up the down edge */
  int i, j, ai, downdown, down;
  double uptime;
  struct edge *gtree = C[ci]->G[li].gtree;
  down = gtree[sis].down;
  i = 0;
  while (gtree[sis].mig[i].mt > -0.5)
    i++;
  j = -1;

  do
  {
    j++;
    checkmig (i + 1, &(gtree[sis].mig), &(gtree[sis].cmm));
    gtree[sis].mig[i] = gtree[down].mig[j];
    i++;
  } while (gtree[down].mig[j].mt > -0.5);
  gtree[sis].time = gtree[down].time;

  /* set the up to which sis now connects */
  gtree[sis].down = gtree[down].down;
  downdown = gtree[sis].down;
  if (downdown != -1)
  {
    rootmove = 0;
    if (gtree[downdown].up[0] == down)
      gtree[downdown].up[0] = sis;
    else
      gtree[downdown].up[1] = sis;
    mrootdrop = 0;
    lmrootdrop = 0;
  }
  else
  {
    rootmove = 1;
    *tmrcachange += 1;
/* figure out total time and number of migrants being dropped */
    i = -1;
    do
    {
      i++;
    } while (gtree[sis].mig[i].mt > -0.5);

    /* mrootdrop and lmrootdrop do not seem to do anything. We may want to
     * delete two variables? */
    mrootdrop = i; 
    if (sis < L[li].numgenes)
    {
      uptime = 0;
    }
    else
    {
      uptime = gtree[gtree[sis].up[0]].time;
    }

    if (uptime < C[ci]->tvals[lastperiodnumber - 1])
    {
      if (C[ci]->G[li].roottime < C[ci]->tvals[lastperiodnumber - 1])
        lmrootdrop = C[ci]->G[li].roottime - uptime;
      else
        lmrootdrop = C[ci]->tvals[lastperiodnumber - 1] - uptime;
    }
    else
    {
      lmrootdrop = 0;
    }

    C[ci]->G[li].root = sis;
    gtree[sis].down = -1;
    gtree[sis].time = TIMEMAX;
    if (L[li].model == STEPWISE || L[li].model == JOINT_IS_SW)
      for (ai = (L[li].model == JOINT_IS_SW); ai < L[li].nlinked; ai++)
        gtree[sis].dlikeA[ai] = 0;
    gtree[sis].mig[0].mt = -1;
    C[ci]->G[li].roottime = uptime;
  }
}                               /* joinsisdown */

void
splitsisdown (int ci, int li, int slidingedge, int down, int newsis)
{

  /* split newsis into two parts, and make a new down edge out of the lower part */
  int i, j, downdown, nowpop;
  double curt;
  struct edge *gtree = C[ci]->G[li].gtree;
  curt = gtree[slidingedge].time;
  gtree[down].time = gtree[newsis].time;
  gtree[newsis].time = curt;

  /* set the up  of the edge to which down now connects, depends on whether newsis is the root */
  downdown = gtree[newsis].down;
  if (downdown != -1)
  {
    if (gtree[downdown].up[0] == newsis)
      gtree[downdown].up[0] = down;
    else
      gtree[downdown].up[1] = down;
  }
  else
  {
    /* newsis is the current root so the root must move down */
    C[ci]->G[li].root = down;
    C[ci]->G[li].roottime = curt;
    rootmove = 1;
    if (C[ci]->G[li].roottime > TIMEMAX)
      IM_err(IMERR_ROOTTIMEMAXFAIL, "roottime greater than TIMEMAX, chain: %d,locus: %d, roottime %lf, TIMEMAX %lf",ci,li,C[ci]->G[li].roottime,TIMEMAX);
    gtree[down].mig[0].mt = -1;
    if (L[li].model == STEPWISE || L[li].model == JOINT_IS_SW)
      for (i = (L[li].model == JOINT_IS_SW); i < L[li].nlinked; i++)
        gtree[down].dlikeA[i] = 0;
  }
  gtree[down].down = downdown;

  /* divide the migration events along newsis into upper part for newsis and lower part for down */
  /* this might have bugs setting the population of gtree[down] */
  i = 0;
  while (gtree[newsis].mig[i].mt > -0.5 && gtree[newsis].mig[i].mt < curt)
    i++;
  if (i > 0)
    nowpop = gtree[newsis].mig[i - 1].mp;
  else
    nowpop = gtree[newsis].pop;
  j = findperiod (ci, curt);
  /* SANGCHUL: Thu Oct  8 14:14:31 EDT 2009
   * We could use
   * nowpop = saC.popndown[nowpop][j];
   * instead of using the followng while statement. 
   * */
  while (C[ci]->poptree[nowpop].e <= j && C[ci]->poptree[nowpop].e != -1)
    nowpop = C[ci]->poptree[nowpop].down;
  gtree[down].pop = nowpop;
  j = 0;
  if (downdown != -1)
  {
    while (gtree[newsis].mig[j + i].mt > -0.5)
    {
      checkmig (j, &(gtree[down].mig), &(gtree[down].cmm));
      gtree[down].mig[j] = gtree[newsis].mig[j + i];
      assert (nowpop != gtree[down].mig[j].mp);
      nowpop = gtree[down].mig[j].mp;
      j++;
    }
  }
  gtree[newsis].mig[i].mt = -1;
  gtree[down].mig[j].mt = -1;
  gtree[newsis].down = gtree[slidingedge].down = down;
  gtree[down].up[0] = newsis;
  gtree[down].up[1] = slidingedge;
  return;
}                               /* splitsisdown */


/* called by addmigration(),  does most of the migration either by adding events or by calling mwork_single_edge() 
    mwork_single_edge is called for the simpler cases, 
    getm() itself handles the case where both edges need migration and the state of the population at the bottom of the edges needs to be 
    simulated  */

void  getm (int ci,struct edgemiginfo *edgem,struct edgemiginfo *sisem, struct edgemiginfo *oldedgem,struct edgemiginfo *oldsisem)
{
  int lastmigperiod;
//  int i, ii;

/*for (ii = 0;ii<2;ii++)
for (i=0;i<5;i++)
{
  checkptype[ii][i] = -100;
  checkrr[ii][i] = 0;
}  //8_30_10 */

  lastmigperiod = IMIN(edgem->e,lastperiodnumber-1);
  if ( sisem->edgeid == -1  /* sisem->mtall <= 0*/)  // no sister edge, or sister edge not in a period where migration can occur
  {
    assert(sisem->mtimeavail[0] == 0);
    edgem->mpall = mwork_single_edge (ci, edgem,oldedgem, lastmigperiod); 
  }
  else
  {
    assert (edgem->e == sisem->e);
    if (edgem->mtall <= 0)   // edge has no length in a period with migration,  so just do sister edge
    {
      assert(edgem->mtimeavail[0] == 0);
      sisem->mpall = mwork_single_edge (ci, sisem,oldsisem, lastmigperiod);
    }
    else
    {  // both edge and sis have length in periods with migration
      mwork_two_edges(ci, edgem, sisem,oldedgem, oldsisem, lastmigperiod, &edgem->mpall, &sisem->mpall);

    }
  }
}  //getm

/* add migration to edge, and to its sister if edge connects to the root, 
 * return the log of the hastings ratio of update probabilities 
 *
 * add migration events to edge that just slid.  Also if it slid down the 
 * root node, and moved the root, then migration events may need to be 
 * added to the sister branch as well 
 * 
 * CR 111006.1 10/6/2011 JH   removed oldmigcount, oltlength, newmigcount 
 * and newtlength these had once been used for calculating the migration 
 * weight but this is now done using the  edgemiginfo structures passed 
 * to and from getmprob()
 *
 * Additional note: local vars mtime and mcount also removed.
 */
double
addmigration (int ci, int li)
{
  int newsis, edge;
  double weight;
  double mproposenum, mproposedenom, temp;
  struct edge *gtree = C[ci]->G[li].gtree;

  assert (C[ci]->G[li].mignum >= 0 && C[ci]->G[li].tlength > 0);
  IMA_reset_edgemiginfo (&newedgemig);
  IMA_reset_edgemiginfo (&newsismig);
  newedgemig.edgeid = edge = oldedgemig.edgeid;
  newedgemig.li = li;
  if (edge < L[li].numgenes)
    newedgemig.upt = 0;
  else
    newedgemig.upt = gtree[gtree[edge].up[0]].time;
  newedgemig.fpop = gtree[gtree[edge].down].pop;
  newedgemig.pop = newedgemig.temppop = gtree[edge].pop;
  newedgemig.dnt = gtree[edge].time;
  newedgemig.mig[0].mt = -1;
  fillmiginfoperiods (ci, &newedgemig);
  if (gtree[edge].down == C[ci]->G[li].root)    /* simulate migration on the sister branch as well */
  {
    //IMA_reset_edgemiginfo (&newsismig);
    newedgemig.fpop = -1;       //pop unknown, as edge must be determined by migration 
    if (gtree[gtree[edge].down].up[0] == edge)
      newsis = gtree[gtree[edge].down].up[1];
    else
      newsis = gtree[gtree[edge].down].up[0];
    if (newsis < L[li].numgenes)
      newsismig.upt = 0;
    else
      newsismig.upt = gtree[gtree[newsis].up[0]].time;
    newsismig.edgeid = newsis;
    newsismig.li = li;
    newsismig.fpop = -1;        //pop unknown, as edge must be determined by migration 
    newsismig.pop = newsismig.temppop = gtree[newsis].pop;
    newsismig.dnt = gtree[newsis].time;
    newsismig.mig[0].mt = -1;
    fillmiginfoperiods (ci, &newsismig);
  }
  else
  {
    newedgemig.fpop = gtree[gtree[edge].down].pop;
    newsismig.edgeid = -1;
    newsismig.mtall = 0; 
    newsismig.b = newsismig.e = -1;// no second edge to deal with
  }

  assert((newsismig.mtall > 0 && newsismig.edgeid >= 0) || newsismig.mtall == 0);


  getm (ci, &newedgemig, &newsismig, &oldedgemig, &oldsismig);

  if (newedgemig.fpop == -1)
  {
    assert(newsismig.edgeid != -1);
    assert(newedgemig.e == lastperiodnumber);
    newedgemig.fpop  = newsismig.fpop  = C[ci]->G[li].root;

  }
  if (gtree[newedgemig.edgeid].down == C[ci]->G[li].root)
    gtree[C[ci]->G[li].root].pop = newedgemig.fpop;


/*checkpp = 1;  //debug 8_30_10 */
    temp = getmprob (ci, &newedgemig, &newsismig, &oldedgemig, &oldsismig);
/*checkpp = 0;  //debug  8_30_10 */

  mproposedenom = temp;
  assert (temp > -1e200 && temp < 1e200);

  /* calculate probability of reverse update    */
  temp = getmprob (ci, &oldedgemig, &oldsismig, &newedgemig, &newsismig);
  mproposenum = temp;
  assert (temp > -1e200 && temp < 1e200);
  weight = mproposenum - mproposedenom;
  return weight;
}                               /* addmigration */

/********GLOBAL FUNCTIONS *******/

// VS 
// changed to initialize the holdgweights and holdpcalc for groups of loci
void
init_updategenealogy (void)
{
  int gp;

  init_genealogy_weights (&holdallgweight_updategenealogy);
  init_genealogy_weights (&holdgweight_updategenealogy);
  init_probcalc (&holdallpcalc_updategenealogy);

  // VS - initialize the holdgweights and holdpcalc for groups of loci
  init_genealogy_weights (&holdgweight_theta_updategenealogy);
  init_genealogy_weights (&holdgweight_mig_updategenealogy);
  init_probcalc (&holdpcalc_theta_updategenealogy);
  init_probcalc (&holdpcalc_mig_updategenealogy);

  // VS - initialize the holdgweights and holdpcalc for groups of loci
  for(gp=0; gp<nbgroupsloci_theta+nbgroupsloci_mig; gp++) {
	init_genealogy_weights(&holdgweight_update_assignloci[gp]);
	init_probcalc (&holdpcalc_update_assignloci[gp]);
	init_genealogy_weights(&holdgweight_theta_update_assignloci[gp]);
	init_probcalc (&holdpcalc_theta_update_assignloci[gp]);
  }
    
  //init_genealogy_weights (&holdgweight_update_assignloci);
  //init_genealogy_weights (&holdgweight_new_update_assignloci);
  //init_probcalc (&holdpcalc_update_assignloci); 
  //init_probcalc (&holdpcalc_new_update_assignloci);



  
}                               // init_updategenealogy

// VS 
// changed to free holdgweights and holdpcalc for groups of loci
void
free_updategenealogy (void)
{
  int gp;

  free_genealogy_weights (&holdallgweight_updategenealogy);
  free_genealogy_weights (&holdgweight_updategenealogy);
  free_probcalc (&holdallpcalc_updategenealogy);

  // VS - free the holdgweights and holdpcalc for groups of loci
  free_genealogy_weights (&holdgweight_theta_updategenealogy);
  free_genealogy_weights (&holdgweight_mig_updategenealogy);
  free_probcalc (&holdpcalc_theta_updategenealogy);
  free_probcalc (&holdpcalc_mig_updategenealogy);

  // VS - free the holdgweights and holdpcalc for updating the assignment 
  // of loci to different groups of loci
  for(gp=0; gp<nbgroupsloci_theta+nbgroupsloci_mig; gp++) {
	free_genealogy_weights(&holdgweight_update_assignloci[gp]);
	free_probcalc (&holdpcalc_update_assignloci[gp]);
	free_genealogy_weights(&holdgweight_theta_update_assignloci[gp]);
	free_probcalc (&holdpcalc_theta_update_assignloci[gp]);
  }

  //free_genealogy_weights (&holdgweight_update_assignloci);
  //free_genealogy_weights (&holdgweight_new_update_assignloci);
  //free_probcalc (&holdpcalc_update_assignloci);
  //free_probcalc (&holdpcalc_new_update_assignloci);

}                               //free_updategenealogy

#define SLIDESTDVMAX 20  

/* steps in picking a new genealogy 
- pick an edge, the bottom of which will slide
- save all info for that edge, sis and the down edge that will be freed up
- join the sis and down edges, freeing up an edge
	if down was the root, then sis becomes the root
-set aside the number for the down edge - this is the freed up edge and will get used again later
-do the sliding and pick a new sis and a location - but do not change the gtree. 
-split the newsis edge into sis and down edges, and divide the migration events accordingly
-connect the original edge to the new spot
-add migration events to this edge
- calculate the probability of the update, in forwrad and reverse directions
- calculate the total Metropolis Hastings terms for genealogy update,  accept or reject 
*/
/* for gtreeprint calls,  use callingsource = 0 */


/* 9/25/08  updated this
revised genealogy updating so that the current migration rate is based on the current number of migration events and the 
current length of the branch that is being updated 
reasoned that this might work better than using the rate that occurs for the entitre genealogy - e.g. help to avoid promoting
correlations and improve mixing  */

/* 7/15/2010  JH revised updategenealogy() and addmigration() and descendant functions to fix a problem associated with the 
Hastings term calculation of the migration path */ 


//prune this   9/21/2010

/* CR 111006.1 10/6/2011 JH   removed local variables mpart and tlengthpart
 * these had once been used for calculating the migration weight
 * but this is now done using the  edgemiginfo structures passed to
 * and from getmprob() in addmigration()
 */ 
int
updategenealogy (int ci, int li, int *topolchange, int *tmrcachange)
{
  int ai, i;
  int edge, oldsis, newsis, freededge, accp;
  double newpdg, newpdg_a[MAXLINKED];
  double migweight, metropolishastingsterm, U;
  double tpw;
  double Aterm[MAXLINKED], Atermsum;
  double slidedist;
  double slideweight, holdslidedist, slidestdv;
  struct genealogy *G = &(C[ci]->G[li]);
  struct edge *gtree = G->gtree;
  int rejectIS;
  double like;
  double holdt[MAXPERIODS];

  // VS
  // variables related with the genealogy groups
  // identify which is the group to which the locus belongs to
  int gp_theta = grouploci_theta[ci][li];
  int gp_mig = grouploci_mig[ci][li];
  //printf("gp_theta and gp_mig here are %d,%d\n\n", gp_theta, gp_mig);
  
  // VS CHECK - open all check files
  /*fgweights_loci = fopen ("gweights_loci.out", "a");
  fgweights_all = fopen ("gweights_all.out", "a");
  fgweights_theta = fopen ("gweights_theta.out", "a");
  fgweights_mig = fopen ("gweights_mig.out", "a");
  fpcalc_all = fopen ("pcalc_all.out", "a");
  fpcalc_theta = fopen ("pcalc_theta.out", "a");
  fpcalc_mig = fopen ("pcalc_mig.out", "a");*/

  if (assignmentoptions[POPULATIONASSIGNMENTCHECKPOINT] == 1)
  {
    assertgenealogyloc (ci, li);
  }
  
  // initialize and make copies structures that hold quantities for calculating prob of genealogy
  copy_treeinfo (&holdgweight_updategenealogy, &G->gweight);
  copy_treeinfo (&holdallgweight_updategenealogy, &C[ci]->allgweight);
  // VS
  // initilize and make copies of group of loci structures that hold quantities for calculating prob of genealogy
  copy_treeinfo (&holdgweight_theta_updategenealogy, &C[ci]->groupgweight_theta[gp_theta]);
  copy_treeinfo (&holdgweight_mig_updategenealogy, &C[ci]->groupgweight_mig[gp_mig]);

  // VS CHECK - print the tree
  //printf("\nTREE before update, ci=%i, li=%i", ci, li);
  //gtreeprint(ci, li);

  // VS CHECK - print the holdgweights
  /*
  fprintf(fgweights_loci, "\nholdgweight_updategenealogy %i ", li);
  print_gweight_vs_file(&holdgweight_updategenealogy, 2, fgweights_loci);
  
  fprintf(fgweights_all, "\nholdALLgweight_updategenealogy %i ", li);
  print_gweight_vs_file(&holdallgweight_updategenealogy, 2, fgweights_all);

  fprintf(fgweights_theta, "\nholdgweight_theta_updategenealogy %i    ", li);
  print_gweight_vs_file(&holdgweight_theta_updategenealogy, 2, fgweights_theta);
  fprintf(fgweights_mig, "\nholdgweight_mig_updategenealogy %i    ", li);
  print_gweight_vs_file(&holdgweight_mig_updategenealogy, 2, fgweights_mig);
  */

  // store summary stats of the genealogy
  storegenealogystats (ci, li, 0);
  for (i = 0; i < lastperiodnumber; i++)
    holdt[i] = C[ci]->tvals[i];  // JH is this actually necessary ? 
  *tmrcachange = 0;
  *topolchange = 0;
  // Atermsum only used for Stepwise mutation model
  Atermsum = 0;

/* pick an edge, identify freedup edge (the down edge) and the sister edge */
  do
  {
    edge = randposint (L[li].numlines);
  } while (gtree[edge].down == -1);
  freededge = gtree[edge].down;
  if ((oldsis = gtree[freededge].up[0]) == edge)
    oldsis = gtree[freededge].up[1];

  /* copy information on the edge,  and if it connects to the root, then the sister edge as well */
  if (gtree[edge].down == G->root)
  {
    fillmiginfo (ci, li, gtree, edge, oldsis);
  }
  else
  {
    fillmiginfo (ci, li, gtree, edge, -1);
  }

  /* store information on the genealogy before changing it */
  storeoldedges (ci, li, edge, oldsis, freededge);
// remove any migrations  from the slidingedge 
  gtree[edge].mig[0].mt = -1;
/* slide edge, pick a distance and slide it  */
  slidestdv = DMIN (SLIDESTDVMAX, G->roottime/3 );
  holdslidedist = slidedist = normdev (0.0, slidestdv); 

// join the sister and the down branches at the point where edge used to connect, this frees up the down branch 
  joinsisdown (ci, li, oldsis, tmrcachange);

                            // use when debugging with gtreeprint()
  //gtree[edge].down = -1;  // not necessary but makes gtreeprint() output easier to read for intermediate stages
  //gtree[freededge].down = -1; // not necessary but makes gtreeprint() output easier to read for intermediate stages

// do the slide and identify the new sister branch and where new connection point for the edge is 
  newsis = oldsis;
  if (modeloptions[NOMIGRATION])
    slider_nomigration (ci, li, edge, &newsis, &(gtree[edge].time),&slidedist);
  else
    slider (ci, li, edge, &newsis, &(gtree[edge].time), &slidedist);
  *topolchange += (oldsis != newsis);

// now separate the new sister branch into a shorter sis branch and a down branch 
  splitsisdown (ci, li, edge, freededge, newsis);

  if (rootmove)
  {
    slideweight = -log (normprob (0.0, slidestdv, holdslidedist));
    slidestdv = DMIN (SLIDESTDVMAX, G->roottime / 3); 
    slideweight += log (normprob (0.0, slidestdv, holdslidedist)); 
  }
  else
  {
    slideweight = 0;
  }

// add migration events 
  if (modeloptions[NOMIGRATION])
  {
    migweight = 0;
  }
  else
  {
    migweight = addmigration (ci, li);
  }

  // copy the migration info in newedgemig and newsismig  to the genealogy
  copynewmig_to_gtree (ci, li);


  // VS CHECK - print the tree
 // printf("\nTREE after G update, ci=%i, li=%i", ci, li);
  //gtreeprint(ci, li);


// determine all the weights needed for calculating the probability of the genealogy
  setzero_genealogy_weights (&G->gweight);
  // TREEWEIGHT
  treeweight (ci, li);

  // VS CHECK
 /* fprintf(fgweights_loci,"\ngweight_AFTER_treeweight %i    ", li);
  print_gweight_vs_file(&G->gweight, 2, fgweights_loci);

  fprintf(fgweights_all,"\nALLgweight_AFTER_treeweight locus=%i    ", li);
  print_gweight_vs_file(&C[ci]->allgweight, 2, fgweights_all);
  */
  
  // VS commented this part - we do not need allgweight anymore
  //sum_subtract_treeinfo (&C[ci]->allgweight, &G->gweight,
  //                       &holdgweight_updategenealogy);

  // VS
  // given that we have groups of loci for theta and for mig
  // need to sum_subtract_treeinfo_vs for the group_gweights
  // these will be done with special functions because when decreasing theta weights
  // only cc, fc and hc are taken into account, 
  // whereas for mig weights, only mc and fm are taken into account
  // update theta
  sum_subtract_treeinfo_vs(&C[ci]->groupgweight_theta[gp_theta], &G->gweight, &holdgweight_updategenealogy, 1); 
  // update mig
  sum_subtract_treeinfo_vs(&C[ci]->groupgweight_mig[gp_mig], &G->gweight, &holdgweight_updategenealogy, 0); 
  // check that the gweihts are set to zero in the correct places (mc=fm=0 for theta, cc=fc=hc=0 for mig)
  check_gweight_vs(&C[ci]->groupgweight_theta[gp_theta], 1);
  check_gweight_vs(&C[ci]->groupgweight_mig[gp_mig], 0);

  // VS CHECK
  /*fprintf(fgweights_theta, "\nC[ci]->groupweight_theta_afterupdate %i ", li);
  print_gweight_vs_file(&C[ci]->groupgweight_theta[gp_theta], 2, fgweights_theta);
  fprintf(fgweights_mig,"\nC[ci]->groupweight_mig_afterupdate %i ", li);
  print_gweight_vs_file(&C[ci]->groupgweight_mig[gp_mig], 2, fgweights_mig);
  */

/* calculate P(D|G)  for new genealogy */
  rejectIS = 0;                 /* use this to catech when P(D|G) for IS model is zero */
  newpdg = 0;

  switch (L[li].model)
  {
  case HKY:
    if (assignmentoptions[JCMODEL] == 1)
    {
      newpdg_a[0] = likelihoodJC (ci, li, G->uvals[0]);
      newpdg = newpdg_a[0];
    }
    else
    {
      newpdg = newpdg_a[0] =
        likelihoodHKY (ci, li, G->uvals[0], G->kappaval, edge,
                       freededge, oldsis, newsis);
    }
    break;
  case INFINITESITES:
    newpdg = newpdg_a[0] = like = likelihoodIS (ci, li, G->uvals[0]);
    rejectIS = (like == REJECTINFINITESITESCONSTANT);
    break;
  case STEPWISE:
    {
      for (ai = 0, newpdg = 0; ai < L[li].nlinked; ai++)
      {
        newpdg_a[ai] =
          G->pdg_a[ai] + finishSWupdateA (ci, li, ai, edge, freededge,
                                          oldsis, newsis,
                                          G->uvals[ai], &Aterm[ai]);

        newpdg += newpdg_a[ai];
        Atermsum += Aterm[ai];
      }

//            checklikelihoodSW(ci, li,G->u[ai].mcinf.val);  
      break;
    }
  case JOINT_IS_SW:
    newpdg = newpdg_a[0] = likelihoodIS (ci, li, G->uvals[0]);
    rejectIS = (newpdg == REJECTINFINITESITESCONSTANT);
    for (ai = 1; ai < L[li].nlinked; ai++)
    {
      newpdg_a[ai] =
        G->pdg_a[ai] + finishSWupdateA (ci, li, ai, edge, freededge,
                                        oldsis, newsis,
                                        G->uvals[ai], &Aterm[ai]);
      newpdg += newpdg_a[ai];
      Atermsum += Aterm[ai];
    }

    //checklikelihoodSW(ci, li,Q[ci]->us[li]);  
    break;
  }
  accp = 0;

/* final weight calculation */
/* tpw is the ratio of new and old prior probability of the genealogies.  It is actually the ratio of the total across all loci,  but
since only genealogy li is being changed at the present time,  the ratio works out to just be the ratio for genealogy li */
  copy_probcalc (&holdallpcalc_updategenealogy, &C[ci]->allpcalc);
  
  // VS - copy the probcalc of each group of loci
  // theta group
  copy_probcalc (&holdpcalc_theta_updategenealogy, &C[ci]->grouppcalc_theta[gp_theta]);
  // check that mintegrate are zero
  check_pcalc_groups_vs(&C[ci]->grouppcalc_theta[gp_theta], 1);
  check_pcalc_groups_vs(&holdpcalc_theta_updategenealogy, 1);
  // mig group
  copy_probcalc (&holdpcalc_mig_updategenealogy, &C[ci]->grouppcalc_mig[gp_mig]);
  // check that qintegrate are zero
  check_pcalc_groups_vs(&C[ci]->grouppcalc_mig[gp_mig], 0);
  check_pcalc_groups_vs(&holdpcalc_mig_updategenealogy, 0);
  
  // VS CHECK
  /*fprintf(fpcalc_all, "\nholdALLprobcalc %i ", li);
  print_pcalc_groups_vs_file(&holdallpcalc_updategenealogy, fpcalc_all);
  fprintf(fpcalc_theta, "\nholdprobcalctheta %i ", li);
  print_pcalc_groups_vs_file(&holdpcalc_theta_updategenealogy, fpcalc_theta);
  fprintf(fpcalc_mig, "\nholdprobcalcmig %i ", li);
  print_pcalc_groups_vs_file(&holdpcalc_mig_updategenealogy, fpcalc_mig);
  */


  /* Find all internal node sequences and mutations of a full genealogy. */
  if (assignmentoptions[POPULATIONASSIGNMENTCHECKPOINT] == 1)
  {
    if (L[li].model == INFINITESITES)
    {
      accp = IMA_genealogy_findIntSeq (ci, li);
      if (accp == 0)
      {
        rejectIS = 1;
      }
      accp = 0;
    }
  }
  if (rejectIS == 0)
  {
    // the metropolis term includes p(D|G) and p(G),  
	// VS - old version
    //tpw = -C[ci]->allpcalc.probg;
	//integrate_tree_prob (ci, &C[ci]->allgweight,
    //                     &holdallgweight_updategenealogy, &C[ci]->allpcalc,
    //                     &holdallpcalc_updategenealogy, &holdt[0]);

	// VS - new version
    // When updating an edge of a given locus just need to consider the gweights of the groups to which the locus being updated belongs to
    // everything else is cancelled out
    tpw = -holdpcalc_theta_updategenealogy.probg - holdpcalc_mig_updategenealogy.probg;

	// CHECK grouppcalc_theta and grouppcalc_mig
	//printf("\n\nCheck C[ci]->grouppcalc_theta and mig BEFORE UPDATE\n");
	//print_pcalc_groups_vs(&C[ci]->grouppcalc_theta[gp_theta]);
	//print_pcalc_groups_vs(&C[ci]->grouppcalc_mig[gp_mig]);


	// VS - new version
	// get the integral for G_gt' for theta
	
	//printf("Calling integrate_tree_prob_vs...4\n\n");
	integrate_tree_prob_vs(ci, &C[ci]->groupgweight_theta[gp_theta], &holdgweight_theta_updategenealogy, &C[ci]->grouppcalc_theta[gp_theta], &holdpcalc_theta_updategenealogy, &holdt[0], gp_theta, 1);
	// get the integral for G_gm' for mig
	
	//printf("Calling integrate_tree_prob_vs...5\n\n");
	integrate_tree_prob_vs(ci, &C[ci]->groupgweight_mig[gp_mig], &holdgweight_mig_updategenealogy, &C[ci]->grouppcalc_mig[gp_mig], &holdpcalc_mig_updategenealogy, &holdt[0], gp_mig, 0);
	// add the new probability for the integral for theta and mig
	tpw += C[ci]->grouppcalc_theta[gp_theta].probg + C[ci]->grouppcalc_mig[gp_mig].probg;
    //tpw += C[ci]->allpcalc.probg;


	// CHECK grouppcalc_theta and grouppcalc_mig
	/*fprintf(fpcalc_theta, "\nC[ci]->grouppcalc_thetaAFTERupdate %i ", li);
	print_pcalc_groups_vs_file(&C[ci]->grouppcalc_theta[gp_theta], fpcalc_theta);
	fprintf(fpcalc_mig, "\nC[ci]->grouppcalc_migAFTERupdate %i ", li);
	print_pcalc_groups_vs_file(&C[ci]->grouppcalc_mig[gp_mig], fpcalc_mig);
	*/
    
	

    metropolishastingsterm = tpw + gbeta*(newpdg - G->pdg);
    U = uniform ();
    metropolishastingsterm = exp (beta[ci] * metropolishastingsterm + migweight + slideweight + Atermsum);
    if (U < DMIN(1.0, metropolishastingsterm))  //9/13/2010 
    {
      /* accept the update */
      C[ci]->allpcalc.pdg -= G->pdg;
      C[ci]->allpcalc.pdg += newpdg;

      G->pdg = newpdg;
      for (ai = 0; ai < L[li].nlinked; ai++)
        G->pdg_a[ai] = newpdg_a[ai];
      if (L[li].model == HKY)
      {
        copyfraclike (ci, li);
        storescalefactors (ci, li);
      }

	  // VS
	  // update the allpcalc.probg
	  // this is the old Log(P(G))+tpw
	  // in this case tpw is the ratio of p(Gg')/P(Gg), where the index g refers to the updated groups
	  C[ci]->allpcalc.probg = holdallpcalc_updategenealogy.probg+tpw;

      accp = 1;
    }
  }
  /* reject the update */
  if (accp == 0)
  {
    // put the edges back 
    restoreedges (ci, li, edge, oldsis, freededge, newsis);

    // copy summary stats back
    storegenealogystats (ci, li, 1);

    // reset HKY terms
    if (L[li].model == HKY)
      restorescalefactors (ci, li);
    // copy back all the weights and results associated with calculating the probability of the genealogy 
    copy_probcalc (&C[ci]->allpcalc, &holdallpcalc_updategenealogy);
    // VS - this is not really needed anymore
	//copy_treeinfo (&C[ci]->allgweight, &holdallgweight_updategenealogy);
    copy_treeinfo (&G->gweight, &holdgweight_updategenealogy);
    *topolchange = 0;
    *tmrcachange = 0;

	// VS
	// copy back all the group weights and pcalc for the groups of loci
	copy_probcalc (&C[ci]->grouppcalc_theta[gp_theta], &holdpcalc_theta_updategenealogy);
	copy_probcalc (&C[ci]->grouppcalc_mig[gp_mig], &holdpcalc_mig_updategenealogy);
	// check that mintegrate are zero
	check_pcalc_groups_vs(&C[ci]->grouppcalc_theta[gp_theta], 1);
	check_pcalc_groups_vs(&holdpcalc_theta_updategenealogy, 1);
	// check that qintegrate are zero
	check_pcalc_groups_vs(&C[ci]->grouppcalc_mig[gp_mig], 0);
	check_pcalc_groups_vs(&holdpcalc_mig_updategenealogy, 0);

	// VS
	// copy back the gweights for groups of loci
	copy_treeinfo (&C[ci]->groupgweight_theta[gp_theta], &holdgweight_theta_updategenealogy);
	copy_treeinfo (&C[ci]->groupgweight_mig[gp_mig], &holdgweight_mig_updategenealogy);
	// check gweights for groups of loci
	check_gweight_vs (&C[ci]->groupgweight_theta[gp_theta], 1);
	check_gweight_vs (&C[ci]->groupgweight_mig[gp_mig], 0);
  }
 
  // CHECK gweights of updated locus
  //fprintf(fgweights_loci,"\ngweight_locus[%i]    ", li);
  //print_gweight_vs_file(&G->gweight, 2, fgweights_loci);


  // CHECK gweights of theta and mig after update
  //fprintf(fgweights_theta, "\nC[ci]->groupweight_theta_afterupdate[0]");
  //print_gweight_vs_file(&C[ci]->groupgweight_theta[0], 2, fgweights_theta);
  //fprintf(fgweights_mig,"\nC[ci]->groupweight_mig_afterupdate[0]");
  //print_gweight_vs_file(&C[ci]->groupgweight_mig[0], 2, fgweights_mig);
  //fprintf(fgweights_mig,"\nC[ci]->groupweight_mig_afterupdate[1]");
  //print_gweight_vs_file(&C[ci]->groupgweight_mig[1], 2, fgweights_mig);

  // VS CHECK probg allpcalc
  // Look at the all pcalc after the update and if it was accepted or rejected
  //fprintf(fpcalc_all, "\nALLprobcalcAFTERacc %i ", li);
  //print_pcalc_groups_vs_file(&C[ci]->allpcalc, fpcalc_all);

  // VS CHECK the pcalc of groups of loci
  //fprintf(fpcalc_theta, "\nC[ci]->grouppcalc_thetaAFTERupdate[0] ");
  //print_pcalc_groups_vs_file(&C[ci]->grouppcalc_theta[0], fpcalc_theta);
  ////fprintf(fpcalc_theta, "\nC[ci]->grouppcalc_thetaAFTERupdate[1] ");
  ////print_pcalc_groups_vs_file(&C[ci]->grouppcalc_theta[1], fpcalc_theta);
  //fprintf(fpcalc_mig, "\nC[ci]->grouppcalc_migAFTERupdate[0] ");
  //print_pcalc_groups_vs_file(&C[ci]->grouppcalc_mig[0], fpcalc_mig);
  //fprintf(fpcalc_mig, "\nC[ci]->grouppcalc_migAFTERupdate[1] ");
  //print_pcalc_groups_vs_file(&C[ci]->grouppcalc_mig[1], fpcalc_mig);

  /* do updates at nodes for stepwise loci, regardless of whether slide update was accepted.  This could go somewhere else  */

  // VS CHECK - Close all check files
  /*fclose(fpcalc_theta);
  fclose(fpcalc_mig);
  fclose(fpcalc_all);

  fclose(fgweights_loci);
  fclose(fgweights_all);
  fclose(fgweights_theta);
  fclose(fgweights_mig);*/

  return accp;
}   /*updategenealogy */


#undef SLIDESTDVMAX



//// UPDATEGROUPLOCIASSIGN_MIG
//// Function that updates the vector of indexes of groups of loci
//// For the moment it only works for two groups of loci
//// Currently we assume that we only update the migration vector.
//// INPUT:
////	int ci : index of chain
//// NOTE: grouploci_mig is a global variable that is updated during this function
//void updategrouplociassign_mig (int ci) {
//
//	char tempfilename[100]; 
//	FILE *fgweights_locus;
//	FILE *fgweights_mig;
//	FILE *fpcalc_mig;
//	FILE *fassign;
//	double tmhratio; // temp Metropolis-Hastings ratio in log scale
//	int i, aux; // auxiliary variables for loops
//	struct genealogy *G; // pointer to the genealogy of chain ci, locus li
//	int li; // index of locus
//	int gp; // index of current group of migration for locus li
//	int new_gp; // index of new group of migration for locus li
//	double mhratio; // metropolis hastings ratio
//	double U; // variable to save sample from uniform distribution [0,1] to decide to accept or reject
//	int sum=0; // variable to check if the array of groups is either 0,0,0,0 or 1,1,1,1
//	double holdt[MAXPERIODS];
//
//	// Store the current assignment vector
//	for(i=0; i<nloci; i++) {
//		holdgroupassign[i] = grouploci_mig[ci][i];
//	}
//	for (i = 0; i < lastperiodnumber; i++) {
//		holdt[i] = C[ci]->tvals[i];
//	}
//
//	// DEBUG - check hold assign vector
//	fassign = fopen("checkholdassign.out","a");
//	// Store the current assignment vector
//	for(i=0; i<nloci; i++) {
//		fprintf(fassign, "%i ", holdgroupassign[i]);
//	}
//	fprintf(fassign, "\n");
//	fclose(fassign);
//
//	// Pick a random locus
//	li = randposint (nloci);
//
//	// Get the current group of locus li
//	gp = grouploci_mig[ci][li]; // NOTE: grouploci_mig is a global variable
//
//	// Point the genealogy to locus li of chain ci
//	G = &(C[ci]->G[li]);
//    // Note that G->gweight is pointing to the gweight 
//	// of locus li, in chain ci
//	// i.e. cc, mc, fc, fm 
//
//	// CHECk all loci gweights
//	fgweights_mig = fopen("checkgweightslocus0.out","a");
//	G = &(C[ci]->G[0]);
//	print_gweight_vs_file(&G->gweight, 2, fgweights_mig);
//	fprintf(fgweights_mig, "\n");
//	fclose(fgweights_mig);
//
//	fgweights_mig = fopen("checkgweightslocus1.out","a");
//	G = &(C[ci]->G[1]);
//	print_gweight_vs_file(&G->gweight, 2, fgweights_mig);
//	fprintf(fgweights_mig, "\n");
//	fclose(fgweights_mig);
//
//	fgweights_mig = fopen("checkgweightslocus2.out","a");
//	G = &(C[ci]->G[2]);
//	print_gweight_vs_file(&G->gweight, 2, fgweights_mig);
//	fprintf(fgweights_mig, "\n");
//	fclose(fgweights_mig);
//
//	fgweights_mig = fopen("checkgweightslocus3.out","a");
//	G = &(C[ci]->G[3]);
//	print_gweight_vs_file(&G->gweight, 2, fgweights_mig);
//	fprintf(fgweights_mig, "\n");
//	fclose(fgweights_mig);
//	
//
//	// Pick a random group for the locus l that is different from
//	// gp_mig
//    do {
//		new_gp = randposint (nbgroupsloci_mig);
//	} while (new_gp == gp);
//
//	// Modify the assignment vector
//	grouploci_mig[ci][li] = new_gp;
//
//	// DEBUG
//	fassign = fopen("checkassign.out","a");
//	// Store the current assignment vector
//	for(i=0; i<nloci; i++) {
//		fprintf(fassign, "%i ", grouploci_mig[ci][i]);
//	}
//	fprintf(fassign, "\n");
//	fclose(fassign);
//
//	// Set the holdgweights to zero
//	setzero_genealogy_weights (&holdgweight_update_assignloci);
//	setzero_genealogy_weights (&holdgweight_new_update_assignloci);
//
//	// Hold the gweights of the group of the locus li, i.e. z_li
//    copy_treeinfo (&holdgweight_update_assignloci, &C[ci]->groupgweight_mig[gp]);
//	// Hold the gweights of the new group of the locus li, i.e. z*_li
//	copy_treeinfo (&holdgweight_new_update_assignloci, &C[ci]->groupgweight_mig[new_gp]);
//
//	// DEBUG - check the holdgweights
//	sprintf(tempfilename, "holdgweightsgroup%i", gp); // create name of file for holdgweights
//	fgweights_mig = fopen(tempfilename,"a");
//	print_gweight_vs_file(&holdgweight_update_assignloci, 2, fgweights_mig);
//	fprintf(fgweights_mig, "\n");
//	fclose(fgweights_mig);
//
//	sprintf(tempfilename, "holdgweightsgroup%i", new_gp); // create name of file for holdgweights
//	fgweights_mig = fopen(tempfilename,"a");
//	print_gweight_vs_file(&holdgweight_new_update_assignloci, 2, fgweights_mig);
//	fprintf(fgweights_mig, "\n");
//	fclose(fgweights_mig);
//
//	// check the gweights
//	check_gweight_vs(&holdgweight_update_assignloci, 0);
//	check_gweight_vs(&holdgweight_new_update_assignloci, 0);
//
//	// Hold the pcalc of the group of the locus li, i.e. z_li
//	// copy group of loci structures with the solution to the Qintegrate and Mintegrate
//	copy_probcalc (&holdpcalc_update_assignloci, &C[ci]->grouppcalc_mig[gp]);
//	// Hold the pcalc of the new group of locus li, i.e. zi*_li
//	copy_probcalc (&holdpcalc_new_update_assignloci, &C[ci]->grouppcalc_mig[new_gp]);
//	
//	// DEBUG - check the holdpcalc for the two groups
//	sprintf(tempfilename, "holdpcalcgroup%i", gp); // create name of file for holdgweights
//	fgweights_mig = fopen(tempfilename,"a");
//	print_pcalc_groups_vs_file(&holdpcalc_update_assignloci, fgweights_mig);
//	fprintf(fgweights_mig, "\n");
//	fclose(fgweights_mig);
//
//	sprintf(tempfilename, "holdpcalcgroup%i", new_gp); // create name of file for holdgweights
//	fgweights_mig = fopen(tempfilename,"a");
//	print_pcalc_groups_vs_file(&holdpcalc_new_update_assignloci, fgweights_mig);
//	fprintf(fgweights_mig, "\n");
//	fclose(fgweights_mig);
//
//	// check that qintegrate are zero
//	check_pcalc_groups_vs(&holdpcalc_update_assignloci, 0);
//	check_pcalc_groups_vs(&holdpcalc_new_update_assignloci, 0);
//
//	// To obtain the new gweights for gp_mig and new_gp_mig
//	// Need to substract gweight of locus li to gp_mig
//	// And add gweight of locus li to new_gp_mig
//	sum_subtract_treeinfo_assignloci_vs(&C[ci]->groupgweight_mig[new_gp], &C[ci]->groupgweight_mig[gp],&G->gweight, 0);  
//	
//	// If the new array has either all elements of the same group, 
//	// we need to set the gweights of that group to zero
//	for(i=1; i<nloci; i++) {
//		if(grouploci_mig[ci][i]==grouploci_mig[ci][i-1]) {
//			sum++;
//		}
//	}
//
//	// If all elements of the array are the same, sum==nloci-1
//	if(sum==(nloci-1)) {
//		// Set gweights of empty groups (i.e. groups for which there are no loci assigned) to be zero.
//		// the group gp index is the same as the index of firt element of grouploci_mig[ci]
//		aux = grouploci_mig[ci][0]; // aux is the group index that is repeated
//		for(i=0;i<nbgroupsloci_mig;i++) {
//			if(i!=aux) { // if the group is not the same as the one to which all loci belong
//				// set the gweights to zero
//				setzero_genealogy_weights(&C[ci]->groupgweight_mig[i]);
//			}
//		}
//	}
//
//	// DEBUG - check the gweights
//	sprintf(tempfilename, "gweightsgroup%i", gp); // create name of file for holdgweights
//	fgweights_mig = fopen(tempfilename,"a");
//	print_gweight_vs_file(&C[ci]->groupgweight_mig[gp], 2, fgweights_mig);
//	fprintf(fgweights_mig, "\n");
//	fclose(fgweights_mig);
//
//	sprintf(tempfilename, "gweightsgroup%i", new_gp); // create name of file for holdgweights
//	fgweights_mig = fopen(tempfilename,"a");
//	print_gweight_vs_file(&C[ci]->groupgweight_mig[new_gp], 2, fgweights_mig);
//	fprintf(fgweights_mig, "\n");
//	fclose(fgweights_mig);
//
//	// check that the gweihts are set to zero in the correct places (mc=fm=0 for theta, cc=fc=hc=0 for mig)
//	check_gweight_vs(&C[ci]->groupgweight_mig[gp], 0);
//	check_gweight_vs(&C[ci]->groupgweight_mig[new_gp], 0);
//
//
//	// With the new weights we can compute the new prior probability
//	// get the integral for G_gm' for gp_mig (NEED TO CHECK THIS)
//	integrate_tree_prob_vs(ci, &C[ci]->groupgweight_mig[gp], &holdgweight_update_assignloci, &C[ci]->grouppcalc_mig[gp], &holdpcalc_update_assignloci, &holdt[0], gp, 0); // Note that holdt was replaced by 0
//	// get the integral for G_gm' for new_gp_mig
//	integrate_tree_prob_vs(ci, &C[ci]->groupgweight_mig[new_gp], &holdgweight_new_update_assignloci, &C[ci]->grouppcalc_mig[new_gp], &holdpcalc_new_update_assignloci, &holdt[0], new_gp, 0); // Note that holdt was replaced by 0
//
//	// DEBUG - check the pcalc for the two groups
//	sprintf(tempfilename, "pcalcgroup%i", gp); // create name of file for holdgweights
//	fgweights_mig = fopen(tempfilename,"a");
//	print_pcalc_groups_vs_file(&C[ci]->grouppcalc_mig[gp], fgweights_mig);
//	fprintf(fgweights_mig, "\n");
//	fclose(fgweights_mig);
//
//	sprintf(tempfilename, "pcalcgroup%i", new_gp); // create name of file for holdgweights
//	fgweights_mig = fopen(tempfilename,"a");
//	print_pcalc_groups_vs_file(&C[ci]->grouppcalc_mig[new_gp], fgweights_mig);
//	fprintf(fgweights_mig, "\n");
//	fclose(fgweights_mig);
//
//	// At this point the structures
//	// &C[ci]->grouppcalc_mig[gp_mig] and &C[ci]->grouppcalc_mig[new_gp_mig]
//	// have stored the new pcalc
//	// The MH ratio is this case is the product of these two probabilities over the old ones
//	tmhratio = C[ci]->grouppcalc_mig[gp].probg + C[ci]->grouppcalc_mig[new_gp].probg;
//	tmhratio -= (holdpcalc_update_assignloci.probg + holdpcalc_new_update_assignloci.probg);
//
//	// DEBUG - check allpcalc.probg before and after update
//	fgweights_mig = fopen("checkPriorBeforeUpdate","a");
//	print_pcalc_groups_vs_file(&C[ci]->allpcalc, fgweights_mig);
//	fprintf(fgweights_mig, "\n");
//	fclose(fgweights_mig);
//	
//    U = uniform ();
//    mhratio = exp (beta[ci] * tmhratio);
//    if (U < DMIN(1.0, mhratio))  //9/13/2010 
//    {
//        // accept the update
//	    
//	    // Need to update the prior probability of the genealogy
//	    // The new prior probability is given by the previous probg
//		// probg-holdpcalc[gp_mig]-holdpcalc[new_gp_mig]+newprobg[gp_mig]+newprobg[new_gp_mig]
//	    C[ci]->allpcalc.probg += tmhratio;
//    }
//	else {
//
//		// Change the assignment vector
//		for(i=0; i<nloci; i++) {
//			grouploci_mig[ci][i] = holdgroupassign[i];
//		}
//	    
//		// copy back all the pcalc for the groups of loci
//		copy_probcalc (&C[ci]->grouppcalc_mig[gp], &holdpcalc_update_assignloci);
//		copy_probcalc (&C[ci]->grouppcalc_mig[new_gp], &holdpcalc_new_update_assignloci);
//		
//		// check that qintegrate are zero
//		check_pcalc_groups_vs(&C[ci]->grouppcalc_mig[gp], 0);
//		check_pcalc_groups_vs(&C[ci]->grouppcalc_mig[new_gp], 0);
//		//check_pcalc_groups_vs(&holdpcalc_update_assignloci, 0);
//		//check_pcalc_groups_vs(&holdpcalc_new_update_assignloci, 0);
//
//		// copy back the gweights for groups of loci
//		copy_treeinfo (&C[ci]->groupgweight_mig[gp], &holdgweight_update_assignloci);
//		copy_treeinfo (&C[ci]->groupgweight_mig[new_gp], &holdgweight_new_update_assignloci);
//		// check gweights for groups of loci
//		check_gweight_vs (&C[ci]->groupgweight_mig[gp], 0);
//		check_gweight_vs (&C[ci]->groupgweight_mig[new_gp], 0);
//	}
//
//	// DEBUG - check allpcalc.probg before and after update
//	fgweights_mig = fopen("checkPriorAfterUpdate","a");
//	print_pcalc_groups_vs_file(&C[ci]->allpcalc, fgweights_mig);
//	fprintf(fgweights_mig, "\n");
//	fclose(fgweights_mig);
//	
//}
//



// UPDATEGROUPLOCIASSIGN_MIG
// Function that updates the vector of indexes of groups of loci
// For the moment it only works for two groups of loci
// Currently we assume that we only update the migration vector.
// VS 5/17/2012 included in this function the Dirichlet process prior when deciding how to accept the vector.
// INPUT:
//	int ci : index of chain
// NOTE: grouploci_mig is a global variable that is updated during this function
// Currently (Dec 2011) only works for 2 populations. 
// Need to include a for look for periods to extend to multiple populations.
// VS 5/17/2012 - I do not understand why this does not work for multiple pop???
int updategrouplociassign_mig (int ci) {

	//FILE *fassign;
	double tmhratio; // temp Metropolis-Hastings ratio in log scale
	int i; // auxiliary variables for loops
	struct genealogy *G; // pointer to the genealogy of chain ci, locus li
	int li; // index of locus
	int gp; // index of current group of migration for locus li
	int new_gp; // index of new group of migration for locus li
	double mhratio; // metropolis hastings ratio
	double U; // variable to save sample from uniform distribution [0,1] to decide to accept or reject
	int sum=0; // variable to check if the array of groups is either 0,0,0,0 or 1,1,1,1
	double holdt[MAXPERIODS]; // holds the times of split
	double sum_probg_gp; // variable that saves the temp prior probability
	int n1=0; // variable to count number of loci belonging to group 1 in current assignment vector
	int n1u=0; // variable to count number of loci belonging to group 1 in updated assignment vector
	double priora=0; // prior ratio of assignment vectors (implementing Dirichlet Process Prior)
	int result=0; // variable returned by the function: 0 if assignment is not updated, 1 if assignment is updated

	// Store the current assignment vector
	for(i=0; i<nloci; i++) {
		holdgroupassign[i] = grouploci_mig[ci][i];
	}
	// DEBUG - check hold assign vector
	/*fassign = fopen("checkholdassign.out","a");
	// Store the current assignment vector
	for(i=0; i<nloci; i++) {
		fprintf(fassign, "%i ", holdgroupassign[i]);
	}
	fprintf(fassign, "\n");
	fclose(fassign);*/

	// Store the times of split
	for (i = 0; i < lastperiodnumber; i++) {
		holdt[i] = C[ci]->tvals[i];
	}

	// Pick a random locus
	li = randposint (nloci);

	// Get the current group of locus li
	gp = grouploci_mig[ci][li]; // NOTE: grouploci_mig is a global variable

	// Pick a random group for the locus l that is different from
	// gp_mig
    do {
		new_gp = randposint (nbgroupsloci_mig);
	} while (new_gp == gp);

	// Modify the assignment vector
	grouploci_mig[ci][li] = new_gp;

	// Copy the gweights into holdgweights
	for(i=0; i<nbgroupsloci_mig; i++) {
		// 1. Save gweights and pcals
		// 1.1. save gweights for mig
		copy_treeinfo (&holdgweight_update_assignloci[i], &C[ci]->groupgweight_mig[i]);
		// check gweights for mig
		check_gweight_vs (&holdgweight_update_assignloci[i], 0);
		// 1.2. save pcalc for groups of loci for mig
		copy_probcalc (&holdpcalc_update_assignloci[i], &C[ci]->grouppcalc_mig[i]);
		// check pcalc
		check_pcalc_groups_vs(&holdpcalc_update_assignloci[i], 0);
		
		// 2. set gweights to zero
		setzero_genealogy_weights (&C[ci]->groupgweight_mig[i]);
	}

	// Go through all loci and add the gweights for the specific groups
	for (i = 0; i < nloci; i++)
    {
		// point to the correct locus
		G = &(C[ci]->G[i]);
		
		// add the gweights of locus i to the corresponding mig group
		sum_treeinfo_mig_vs(&C[ci]->groupgweight_mig[grouploci_mig[ci][i]], &G->gweight);
    }

	// Compute the mintegrate for each group and obtain the prior probability p(G|a')
	sum_probg_gp=0;
    for(i=0; i<nbgroupsloci_mig; i++) {
		// mintegrate for mig
		integrate_tree_prob_vs(ci, &C[ci]->groupgweight_mig[i], &holdgweight_update_assignloci[i], &C[ci]->grouppcalc_mig[i], &holdpcalc_update_assignloci[i], &holdt[0], i, 0);
		// check pcalc for mig
		check_pcalc_groups_vs( &C[ci]->grouppcalc_mig[i], 0);
		// add the probg of group to total sum_probg_gp
		sum_probg_gp +=  C[ci]->grouppcalc_mig[i].probg;
    }

	// Initialize the tmhratio as the prior probability of the new assignment vector
	tmhratio = sum_probg_gp;

	// Obtain the prior probability of the previous genealogy and assignment p(G|a)
	sum_probg_gp=0;
    for(i=0; i<nbgroupsloci_mig; i++) {
		// add the probg of group to total sum_probg_gp
		sum_probg_gp +=  holdpcalc_update_assignloci[i].probg;
    }
	
	// Decrease the sum_prob_gp from the tmhratio
	tmhratio -= sum_probg_gp;

	
	// VS 5/17/2012 included Dirichlet process prior (DPP) on assignment
	// The Dirichlet process prior is explained cleary in Huelsenback and Andolfatto (2007) STRUCTURAMA paper.
	// Basically, it depends on the paramater alpha.
	// The probability of two individuals (loci) being co-assigned is given by 1/(1+alpha).
	// If a priori we assume that two individuals have 0.5 probability of belonging to the same cluster, then alpha=1.
	// When alpha=1, the DPP is given by, for an assignment a given k groups
	// p(a,k|alpha=1,n)= [PROD_(i=1)^(i=k) (ni-1)!] / n!
	// where n is the total number of loci (or individuals), and ni is the number of loci in group i.
	// With two groups, the prior ratio simplifies to
	// p(a')/p(a)= [(n1'-1)!(n2'-1)!]/[(n1-1)!(n2-1)!]
	// Noting that n1'=n1+1 or n1'=n1-1, given that we update one of the loci to a different group
	// and that n2=n-n1, the ratio simplifies to:
	// (1) if n1'=n1+1 -> p(a')/p(a)=n1/(n-n1')
	// (2) if n1'=n1-1 -> p(a')/p(a)=(n-n1)/(n1-1)
	// This is only true if n1>0,n1'>0,n2>0,n2'>0!
	// We assume that n1 is the number of elements in group 1.
	// Get the number of n1
	/*for(i=0; i<nloci; i++) {
		// assignment a
		n1 += holdgroupassign[i];
		// assignment a'
		n1u += grouploci_mig[ci][i]; 
	}
	// DPP prior
	// If the update does not change the number of groups k
	if((n1>0 && n1<nloci) && (n1u>0 && n1u<nloci)) {
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

	// Equal Probability Prior
	// See notes on "Notespriorassignment.nb" in doc folder for details.
	// We assume that n1 is the number of elements in group 1.
	// Equal probability prior, in which the probability of p(n1) are all equally likely
	// In this case p(a,n1)=p(a|n1)p(n1)
	// We assume that all p(n1) are equally likely
	// p(n1,n2)=(n+1)^(-1) for all n1 and n2 possible combinations
	// Note that n2=n-n1, so the probability only depends on one of them p(a|n1)=p(a|n1,n2)
	// Thus, the ratio of two assignment vectors
	// p(a',n1')/p(a,n1)=p(a'|n1')/p(a|n1)=Choose(n,n1)/Choose(n,n1')=[n1'!(n-n1')!]/[n1!(n-n1)!]
	// In out update, we can only have n1'=n1+1 or n1'=n1-1
	// (1) p(a'|n1'=n1+1)/p(a|n1)=(n1+1)/(n-n1)=n1'/(n-n1)
	// (2) p(a'|n1'=n1-1)/p(a|n1)=(n-n1+1)/n1
	// Get the prior ratio
	// Get the number of n1
	for(i=0; i<nloci; i++) {
		// assignment a
		n1 += holdgroupassign[i];
		// assignment a'
		n1u += grouploci_mig[ci][i]; 
	}
	if(n1u>n1) { // if n1'=n1+1
		priora = log(n1u)-log(nloci-n1);	
	}
	else {
		priora = log(nloci-n1+1)-log(n1);	
	}

	// VS CHECK PRINT
	//printf("\n step=%i, n1=%i, n1u=%i, priora=%g", step, n1, n1u, priora);
	
    U = uniform ();
    mhratio = exp (beta[ci] * (tmhratio+priora)); // VS - version with the DPP prior
	//mhratio = exp (beta[ci] * (tmhratio)); 
    if (U < DMIN(1.0, mhratio))  //9/13/2010 
    {
        // accept the update
		result = 1;

	    // Need to update the prior probability of the genealogy
	    // The new prior probability is given by the previous probg
		// probg-holdpcalc[gp_mig]-holdpcalc[new_gp_mig]+newprobg[gp_mig]+newprobg[new_gp_mig]
	    C[ci]->allpcalc.probg += tmhratio;

    }
	else {

		// reject the update
		result = 0;

		// Change the assignment vector
		for(i=0; i<nloci; i++) {
			grouploci_mig[ci][i] = holdgroupassign[i];
		}
	    
		// copy back all the pcalc for the groups of loci
		// copy back the gweights for groups of loci
		for(gp=0; gp<nbgroupsloci_mig; gp++) {
			copy_probcalc (&C[ci]->grouppcalc_mig[gp], &holdpcalc_update_assignloci[gp]);
			// check that qintegrate are zero
			check_pcalc_groups_vs(&C[ci]->grouppcalc_mig[gp], 0);
			// copy the gweights
			copy_treeinfo (&C[ci]->groupgweight_mig[gp], &holdgweight_update_assignloci[gp]);
			// check gweights for groups of loci
			check_gweight_vs (&C[ci]->groupgweight_mig[gp], 0);
		}
	}

	return result;
}


// UPDATEGROUPLOCIASSIGN_THETA
// Function that updates the vector of indexes of groups of loci
// For the moment it only works for two groups of loci
// Currently we assume that we only update the migration vector.
// INPUT:
//	int ci : index of chain
// NOTE: grouploci_mig is a global variable that is updated during this function
//void updategrouplociassign_theta (int ci) {
//
//	int i, aux; // auxiliary variables for loops
//	struct genealogy *G; // pointer to the genealogy of chain ci, locus li
//	int li; // index of locus
//	int gp; // index of current group of migration for locus li
//	int new_gp; // index of new group of migration for locus li
//	double mhratio; // metropolis hastings ratio (in natural scale)
//	double tmhratio; // metropolis hastings ratio in log scale
//	double U; // variable to save sample from uniform distribution [0,1] to decide to accept or reject
//	int sum=0; // variable to check if the array of groups is either 0,0,0,0 or 1,1,1,1
//
//	// Store the current assignment vector
//	for(i=0; i<nloci; i++) {
//		holdgroupassign[i] = grouploci_theta[ci][i];
//	}
//
//	// Pick a random locus
//	li = randposint (nloci);
//
//	// Get the current group of locus li
//	gp = grouploci_theta[ci][li]; // NOTE: grouploci_theta is a global variable
//
//	// Point the genealogy to locus li of chain ci
//	G = &(C[ci]->G[li]);
//    // Note that G->gweight is pointing to the gweight 
//	// of locus li, in chain ci
//	// i.e. cc, mc, fc, fm 
//
//	// Pick a random group for the locus l that is different from
//	// gp_theta
//    do {
//		new_gp = randposint (nbgroupsloci_theta);
//	} while (new_gp == gp);
//
//	// Modify the assignment vector
//	grouploci_theta[ci][li] = new_gp;
//
//	// Set the holdgweights to zero
//	setzero_genealogy_weights (&holdgweight_update_assignloci);
//	setzero_genealogy_weights (&holdgweight_new_update_assignloci);
//
//	// Hold the gweights of the group of the locus li, i.e. z_li
//    copy_treeinfo (&holdgweight_update_assignloci, &C[ci]->groupgweight_theta[gp]);
//	// Hold the gweights of the new group of the locus li, i.e. z*_li
//	copy_treeinfo (&holdgweight_new_update_assignloci, &C[ci]->groupgweight_theta[new_gp]);
//
//	// check the gweights
//	check_gweight_vs(&holdgweight_update_assignloci, 1);
//	check_gweight_vs(&holdgweight_new_update_assignloci, 1);
//
//	// Hold the pcalc of the group of the locus li, i.e. z_li
//	// copy group of loci structures with the solution to the Qintegrate and Mintegrate
//	copy_probcalc (&holdpcalc_update_assignloci, &C[ci]->grouppcalc_theta[gp]);
//	// Hold the pcalc of the new group of locus li, i.e. zi*_li
//	copy_probcalc (&holdpcalc_new_update_assignloci, &C[ci]->grouppcalc_theta[new_gp]);
//	
//	// check that mintegrate are zero
//	check_pcalc_groups_vs(&holdpcalc_update_assignloci, 1);
//	check_pcalc_groups_vs(&holdpcalc_new_update_assignloci, 1);
//
//	// To obtain the new gweights for gp_theta and new_gp_theta
//	// Need to substract gweight of locus li to gp_theta
//	// And add gweight of locus li to new_gp_theta
//	sum_subtract_treeinfo_assignloci_vs(&C[ci]->groupgweight_theta[new_gp], &C[ci]->groupgweight_theta[gp],&G->gweight, 1);  
//	
//	// If the new array has either all elements of the same group, 
//	// we need to set the gweights of that group to zero
//	for(i=1; i<nloci; i++) {
//		if(grouploci_theta[ci][i]==grouploci_theta[ci][i-1]) {
//			sum++;
//		}
//	}
//
//	// If all elements of the array are the same, sum==nloci-1
//	if(sum==(nloci-1)) {
//		// Set gweights of empty groups (i.e. groups for which there are no loci assigned) to be zero.
//		// the group gp index is the same as the index of firt element of grouploci_mig[ci]
//		aux = grouploci_theta[ci][0]; // aux is the group index that is repeated
//		for(i=0;i<nbgroupsloci_theta;i++) {
//			if(i!=aux) { // if the group is not the same as the one to which all loci belong
//				// set the gweights to zero
//				setzero_genealogy_weights(&C[ci]->groupgweight_theta[i]);
//			}
//		}
//	}
//	
//	// check that the gweihts are set to zero in the correct places (mc=fm=0 for theta, cc=fc=hc=0 for mig)
//	check_gweight_vs(&C[ci]->groupgweight_theta[gp], 1);
//	check_gweight_vs(&C[ci]->groupgweight_theta[new_gp], 1);
//
//
//	// With the new weights we can compute the new prior probability
//	// get the integral for G_gm' for gp_theta (NEED TO CHECK THIS)
//	integrate_tree_prob_vs(ci, &C[ci]->groupgweight_theta[gp], &holdgweight_update_assignloci, &C[ci]->grouppcalc_theta[gp], &holdpcalc_update_assignloci, 0, gp, 1); // Note that holdt was replaced by 0
//	// get the integral for G_gm' for new_gp_theta
//	integrate_tree_prob_vs(ci, &C[ci]->groupgweight_theta[new_gp], &holdgweight_new_update_assignloci, &C[ci]->grouppcalc_theta[new_gp], &holdpcalc_new_update_assignloci, 0, new_gp, 1); // Note that holdt was replaced by 0
//
//	// At this point the structures
//	// &C[ci]->grouppcalc_theta[gp_theta] and &C[ci]->grouppcalc_theta[new_gp_theta]
//	// have stored the new pcalc
//	// The MH ratio is this case is the product of these two probabilities over the old ones
//	tmhratio = C[ci]->grouppcalc_theta[gp].probg + C[ci]->grouppcalc_theta[new_gp].probg;
//	tmhratio -= holdpcalc_update_assignloci.probg + holdpcalc_new_update_assignloci.probg;
//
//    U = uniform ();
//    mhratio = exp (beta[ci] * tmhratio);
//    if (U < DMIN(1.0, mhratio))  //9/13/2010 
//    {
//
//	    // Need to update the prior probability of the genealogy
//	    // The new prior probability is given by the previous probg
//		// probg-holdpcalc[gp_theta]-holdpcalc[new_gp_theta]+newprobg[gp_theta]+newprobg[new_gp_theta]
//	    C[ci]->allpcalc.probg += tmhratio;
//    }
//	else {
//
//		// Change the assignment vector
//		for(i=0; i<nloci; i++) {
//			grouploci_theta[ci][i] = holdgroupassign[i];
//		}
//    
//		// copy back all the weights and results associated with calculating the probability of the genealogy 
//		// copy_probcalc (&C[ci]->allpcalc, &holdallpcalc_update_assignloci);
//	    
//		// copy back all the group weights and pcalc for the groups of loci
//		copy_probcalc (&C[ci]->grouppcalc_theta[gp], &holdpcalc_update_assignloci);
//		copy_probcalc (&C[ci]->grouppcalc_theta[new_gp], &holdpcalc_new_update_assignloci);
//		
//		// check that qintegrate are zero
//		check_pcalc_groups_vs(&C[ci]->grouppcalc_theta[gp], 1);
//		check_pcalc_groups_vs(&C[ci]->grouppcalc_theta[new_gp], 1);
//		//check_pcalc_groups_vs(&holdpcalc_update_assignloci, 1);
//		//check_pcalc_groups_vs(&holdpcalc_new_update_assignloci, 1);
//
//		// copy back the gweights for groups of loci
//		copy_treeinfo (&C[ci]->groupgweight_theta[gp], &holdgweight_update_assignloci);
//		copy_treeinfo (&C[ci]->groupgweight_theta[new_gp], &holdgweight_new_update_assignloci);
//		// check gweights for groups of loci
//		check_gweight_vs (&C[ci]->groupgweight_theta[gp], 1);
//		check_gweight_vs (&C[ci]->groupgweight_theta[new_gp], 1);
//	}
//}



// UPDATEGROUPLOCIASSIGN_THETA
// Function that updates the vector of indexes of groups of loci
// For the moment it only works for two groups of loci
// Currently we assume that we only update the migration vector.
// INPUT:
//	int ci : index of chain
// NOTE: grouploci_mig is a global variable that is updated during this function
// Currently (Dec 2011) only works for 2 populations. 
// Need to include a for look for periods to extend to multiple populations.
int updategrouplociassign_theta (int ci) {

	//FILE *fassign;
	double tmhratio; // temp Metropolis-Hastings ratio in log scale
	int i; // auxiliary variables for loops
	struct genealogy *G; // pointer to the genealogy of chain ci, locus li
	int li; // index of locus
	int gp; // index of current group of migration for locus li
	int new_gp; // index of new group of migration for locus li
	double mhratio; // metropolis hastings ratio
	double U; // variable to save sample from uniform distribution [0,1] to decide to accept or reject
	int sum=0; // variable to check if the array of groups is either 0,0,0,0 or 1,1,1,1
	double holdt[MAXPERIODS]; // holds the times of split
	double sum_probg_gp; // variable that saves the temp prior probability
	int result=0; // variable returned by the function: 0 if assignment is not updated, 1 if assignment is updated

	// Store the current assignment vector
	for(i=0; i<nloci; i++) {
		holdgroupassign[i] = grouploci_theta[ci][i];
	}
	// DEBUG - check hold assign vector
	/*fassign = fopen("checkholdassign.out","a");
	// Store the current assignment vector
	for(i=0; i<nloci; i++) {
		fprintf(fassign, "%i ", holdgroupassign[i]);
	}
	fprintf(fassign, "\n");
	fclose(fassign);*/

	// Store the times of split
	for (i = 0; i < lastperiodnumber; i++) {
		holdt[i] = C[ci]->tvals[i];
	}

	// Pick a random locus
	li = randposint (nloci);

	// Get the current group of locus li
	gp = grouploci_theta[ci][li]; // NOTE: grouploci_mig is a global variable

	// Pick a random group for the locus l that is different from gp
    do {
		new_gp = randposint (nbgroupsloci_theta);
	} while (new_gp == gp);

	// Modify the assignment vector
	grouploci_theta[ci][li] = new_gp;

	// Copy the gweights into holdgweights
	for(i=0; i<nbgroupsloci_theta; i++) {
		// 1.1. save gweights for theta
		copy_treeinfo (&holdgweight_update_assignloci[i], &C[ci]->groupgweight_theta[i]);
		// check gweights for theta
		check_gweight_vs (&holdgweight_update_assignloci[i], 1); // 1 for theta
		// 1.2. save pcalc for groups of loci for theta
		copy_probcalc (&holdpcalc_update_assignloci[i], &C[ci]->grouppcalc_theta[i]);
		// check pcalc
		check_pcalc_groups_vs(&holdpcalc_update_assignloci[i], 1); // 1 for theta
		// 2. set gweights to zero
		setzero_genealogy_weights (&C[ci]->groupgweight_theta[i]);
	}

	// Go through all loci and add the gweights for the specific groups
	for (i = 0; i < nloci; i++)
    {
		// point to the correct locus
		G = &(C[ci]->G[i]);
		
		// add the gweights of locus i to the corresponding theta group
		sum_treeinfo_theta_vs(&C[ci]->groupgweight_theta[grouploci_theta[ci][i]], &G->gweight);
    }

	// Compute the qintegrate for each group and obtain the prior probability
	sum_probg_gp=0;
    for(i=0; i<nbgroupsloci_theta; i++) {
		// qintegrate for theta
		//AS: debug
		//printf("calling integrate_tree_prob_vs...1\n\n");
		integrate_tree_prob_vs(ci, &C[ci]->groupgweight_theta[i], &holdgweight_update_assignloci[i], &C[ci]->grouppcalc_theta[i], &holdpcalc_update_assignloci[i], &holdt[0], i, 1); // 1 at last argument indicates to compute qintegrate for theta parameters
		// check pcalc for theta
		check_pcalc_groups_vs( &C[ci]->grouppcalc_theta[i], 1);
		// add the probg of group to total sum_probg_gp
		sum_probg_gp +=  C[ci]->grouppcalc_theta[i].probg;
    }

	// Initialize the tmhratio as the prior probability of the new assignment vector
	tmhratio = sum_probg_gp;

	// Obtain the prior probability of the previous assignment
	sum_probg_gp=0;
    for(i=0; i<nbgroupsloci_theta; i++) {
		// add the probg of group to total sum_probg_gp
		sum_probg_gp +=  holdpcalc_update_assignloci[i].probg;
    }
	
	// Decrease the sum_prob_gp from the tmhratio
	tmhratio -= sum_probg_gp;
	
    U = uniform ();
    mhratio = exp (beta[ci] * tmhratio);
    if (U < DMIN(1.0, mhratio))  //9/13/2010 
    {
        // accept the update
		result = 1;
	    
	    // Need to update the prior probability of the genealogy
	    // The new prior probability is given by the previous probg
		// probg-holdpcalc[gp_mig]-holdpcalc[new_gp_mig]+newprobg[gp_mig]+newprobg[new_gp_mig]
	    C[ci]->allpcalc.probg += tmhratio;
    }
	else {

        // reject the update
		result = 0;

		// Change the assignment vector
		for(i=0; i<nloci; i++) {
			grouploci_theta[ci][i] = holdgroupassign[i];
		}
	    
		// copy back all the pcalc for the groups of loci
		// copy back the gweights for groups of loci
		for(gp=0; gp<nbgroupsloci_theta; gp++) {
			copy_probcalc (&C[ci]->grouppcalc_theta[gp], &holdpcalc_update_assignloci[gp]);
			// check that qintegrate are zero
			check_pcalc_groups_vs(&C[ci]->grouppcalc_theta[gp], 1);
			// copy the gweights
			copy_treeinfo (&C[ci]->groupgweight_theta[gp], &holdgweight_update_assignloci[gp]);
			// check gweights for groups of loci
			check_gweight_vs (&C[ci]->groupgweight_theta[gp], 1);
		}
	}

	return result;
}






// UPDATEGROUPLOCIASSIGN_THETAMIG
// Function that updates the vector of indexes of groups of loci for theta and mig.
// In this case it is assumed that both vectors are the same.
// INPUT:
//	int ci : index of chain
// OUTPUT: Returns 1 if the assignment vector was updated, zero otherwise.
// NOTE: grouploci_mig amd grouploci_theta are global variables updated during this function
// Currently (Dec 2011) only works for 2 populations. 
// Need to include a for look for periods to extend to multiple populations.
// VS 5/17/2012 - I do not understand why this does not work for multiple pop???
int updategrouplociassign_thetamig (int ci) {

	//FILE *fassign;
	double tmhratio; // temp Metropolis-Hastings ratio in log scale
	int i; // auxiliary variables for loops
	struct genealogy *G; // pointer to the genealogy of chain ci, locus li
	int li; // index of locus
	int gp; // index of current group of migration for locus li
	int new_gp; // index of new group of migration for locus li
	double mhratio; // metropolis hastings ratio
	double U; // variable to save sample from uniform distribution [0,1] to decide to accept or reject
	int sum=0; // variable to check if the array of groups is either 0,0,0,0 or 1,1,1,1
	double holdt[MAXPERIODS]; // holds the times of split
	double sum_probg_gp; // variable that saves the temp prior probability
	double priora; // prior ratio of assignment
	int result=0; // variable returned by the function: 0 if assignment is not updated, 1 if assignment is updated
	int n1=0, n1u=0; // number of elements in group 1 for current assignment and updated assignment

	// Check that the two assignment vectors are the same
	for(i=0; i<nloci;i++) {
		if(grouploci_theta[ci][i]!=grouploci_mig[ci][i]) {
			printf("\nError in updategrouplociassign_thetamig. The two assignment vectors are different.");
			exit(1);
		}
	}

	// Store the times of split
	for (i = 0; i < lastperiodnumber; i++) {
		holdt[i] = C[ci]->tvals[i];
	}

	// Store the current assignment vector (only need to save one of them as the two are equal)
	for(i=0; i<nloci; i++) {
		holdgroupassign[i] = grouploci_theta[ci][i];
	}
	
	// Pick a random locus to update
	li = randposint (nloci);

	// Get the current group of locus li
	gp = grouploci_theta[ci][li]; // NOTE: grouploci_theta is a global variable

	// Pick a random group for the locus l that is different from gp
	// again note that in this case the nbgroupsloci_theta=nbgroupsloci_mig
    do {
		new_gp = randposint (nbgroupsloci_theta);
	} while (new_gp == gp);

	// Modify the assignment vectors for theta and mig
	grouploci_theta[ci][li] = new_gp;
	grouploci_mig[ci][li] = new_gp;

	// Copy the gweights into holdgweights for the theta groups
	for(i=0; i<nbgroupsloci_theta; i++) {
		// 1.1. save gweights for theta
		copy_treeinfo (&holdgweight_theta_update_assignloci[i], &C[ci]->groupgweight_theta[i]);
		// check gweights for theta
		check_gweight_vs (&holdgweight_theta_update_assignloci[i], 1); // 1 for theta
		// 1.2. save pcalc for groups of loci for theta
		copy_probcalc (&holdpcalc_theta_update_assignloci[i], &C[ci]->grouppcalc_theta[i]);
		// check pcalc
		check_pcalc_groups_vs(&holdpcalc_theta_update_assignloci[i], 1); // 1 for theta
		// 2. set gweights to zero
		setzero_genealogy_weights (&C[ci]->groupgweight_theta[i]);
	}

	// Copy the gweights into holdgweights for the migration groups
	for(i=0; i<nbgroupsloci_mig; i++) {
		// 1.1. save gweights for mig
		copy_treeinfo (&holdgweight_update_assignloci[i], &C[ci]->groupgweight_mig[i]);
		// check gweights for mig
		check_gweight_vs (&holdgweight_update_assignloci[i], 0); // 0 for mig
		// 1.2. save pcalc for groups of loci for mig
		copy_probcalc (&holdpcalc_update_assignloci[i], &C[ci]->grouppcalc_mig[i]);
		// check pcalc
		check_pcalc_groups_vs(&holdpcalc_update_assignloci[i], 0); // 0 for mig
		// 2. set gweights to zero
		setzero_genealogy_weights (&C[ci]->groupgweight_mig[i]);
	}

	// Go through all loci and add the gweights for the specific groups
	for (i = 0; i < nloci; i++)
    {
		// point to the correct locus
		G = &(C[ci]->G[i]);
		
		// add the gweights of locus i to the corresponding theta group
		sum_treeinfo_theta_vs(&C[ci]->groupgweight_theta[grouploci_theta[ci][i]], &G->gweight);

		// add the gweights of locus i to the corresponding mig group
		sum_treeinfo_mig_vs(&C[ci]->groupgweight_mig[grouploci_mig[ci][i]], &G->gweight);
    }

	// Compute the mintegrate for each group and obtain the prior probability
	sum_probg_gp=0;
    for(i=0; i<nbgroupsloci_theta; i++) {
		// qintegrate for theta
		//AS: debug
		//printf("Calling integrate_tree_prob_vs...2\n\n");
		integrate_tree_prob_vs(ci, &C[ci]->groupgweight_theta[i], &holdgweight_theta_update_assignloci[i], &C[ci]->grouppcalc_theta[i], &holdpcalc_theta_update_assignloci[i], &holdt[0], i, 1); // 1 at last argument indicates to compute qintegrate for theta parameters
		// check pcalc for theta
		check_pcalc_groups_vs( &C[ci]->grouppcalc_theta[i], 1); // 1 for theta
		// add the probg of group to total sum_probg_gp
		sum_probg_gp +=  C[ci]->grouppcalc_theta[i].probg;

		// mintegrate for mig

		//printf("Calling integrate_tree_prob_vs...3\n\n");
		integrate_tree_prob_vs(ci, &C[ci]->groupgweight_mig[i], &holdgweight_update_assignloci[i], &C[ci]->grouppcalc_mig[i], &holdpcalc_update_assignloci[i], &holdt[0], i, 0); // 0 at last argument indicates to compute mintegrate for mig parameters
		// check pcalc for mig
		check_pcalc_groups_vs( &C[ci]->grouppcalc_mig[i], 0); // 0 for mig
		// add the probg of group to total sum_probg_gp
		sum_probg_gp +=  C[ci]->grouppcalc_mig[i].probg;
    }

	// Initialize the tmhratio as the prior probability of the new assignment vector
	tmhratio = sum_probg_gp;

	// Obtain the prior probability of the previous assignment
	sum_probg_gp=0;
    for(i=0; i<nbgroupsloci_theta; i++) {
		// add the probg of group to total sum_probg_gp
		sum_probg_gp +=  holdpcalc_update_assignloci[i].probg; // variable holding mig pcalc
		sum_probg_gp +=  holdpcalc_theta_update_assignloci[i].probg; // variable holding theta pcalc
    }
	
	// Decrease the sum_prob_gp from the tmhratio
	tmhratio -= sum_probg_gp;
	

	// VS 6/8/2012 included Equal Prior probability
	// See notes on "Notespriorassignment.nb" in doc folder for details.
	// We assume that n1 is the number of elements in group 1.
	// Equal probability prior, in which the probability of p(n1) are all equally likely
	// In this case p(a,n1)=p(a|n1)p(n1)
	// We assume that all p(n1) are equally likely
	// p(n1,n2)=(n+1)^(-1) for all n1 and n2 possible combinations
	// Note that n2=n-n1, so the probability only depends on one of them p(a|n1)=p(a|n1,n2)
	// Thus, the ratio of two assignment vectors
	// p(a',n1')/p(a,n1)=p(a'|n1')/p(a|n1)=Choose(n,n1)/Choose(n,n1')=[n1'!(n-n1')!]/[n1!(n-n1)!]
	// In out update, we can only have n1'=n1+1 or n1'=n1-1
	// (1) p(a'|n1'=n1+1)/p(a|n1)=(n1+1)/(n-n1)=n1'/(n-n1)
	// (2) p(a'|n1'=n1-1)/p(a|n1)=(n-n1+1)/n1
	// Get the prior ratio
	// Get the number of n1
	for(i=0, n1=0, n1u=0; i<nloci; i++) {
		// assignment a
		n1 += holdgroupassign[i];
		// assignment a'
		n1u += grouploci_mig[ci][i]; 
	}
	if(n1u>n1) { // if n1'=n1+1
		priora = log(n1u)-log(nloci-n1);	
	}
	else {
		priora = log(nloci-n1+1)-log(n1);	
	}

	// VS CHECK PRINT
	//printf("\n step=%i, n1=%i, n1u=%i, priora=%g", step, n1, n1u, priora);
	
    U = uniform ();
    mhratio = exp (beta[ci] * (tmhratio+priora));
    if (U < DMIN(1.0, mhratio))  //9/13/2010 
    {
        // accept the update
		result = 1;
	    
	    // Need to update the prior probability of the genealogy
	    // The new prior probability is given by the previous probg
		// probg-holdpcalc[gp_mig]-holdpcalc[new_gp_mig]+newprobg[gp_mig]+newprobg[new_gp_mig]
	    C[ci]->allpcalc.probg += tmhratio;
    }
	else {

        // reject the update
		result = 0;

		// Change the assignment vector
		for(i=0; i<nloci; i++) {
			grouploci_theta[ci][i] = holdgroupassign[i]; // change the theta
			grouploci_mig[ci][i] = holdgroupassign[i]; // change the mig
		}
	    
		// copy back all the pcalc for the groups of loci
		// copy back the gweights for groups of loci
		for(gp=0; gp<nbgroupsloci_theta; gp++) {
			
			// copy back the gweights and pcalc for theta
			copy_probcalc (&C[ci]->grouppcalc_theta[gp], &holdpcalc_theta_update_assignloci[gp]);
			// check that mintegrate are zero
			check_pcalc_groups_vs(&C[ci]->grouppcalc_theta[gp], 1); // 1 for theta
			// copy the gweights
			copy_treeinfo (&C[ci]->groupgweight_theta[gp], &holdgweight_theta_update_assignloci[gp]);
			// check gweights for groups of loci
			check_gweight_vs (&C[ci]->groupgweight_theta[gp], 1); // 1 for theta
			
			// copy back the gweights and pcalc for mig
			copy_probcalc (&C[ci]->grouppcalc_mig[gp], &holdpcalc_update_assignloci[gp]);
			// check that mintegrate are zero
			check_pcalc_groups_vs(&C[ci]->grouppcalc_mig[gp], 0); // 0 for mig
			// copy the gweights
			copy_treeinfo (&C[ci]->groupgweight_mig[gp], &holdgweight_update_assignloci[gp]);
			// check gweights for groups of loci
			check_gweight_vs (&C[ci]->groupgweight_mig[gp], 0); // 0 for mig
		}
	}

	return result;
}
