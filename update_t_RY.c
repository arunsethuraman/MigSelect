/* MigSelect 2013-2019 Arun Sethuraman, Vitor Sousa, Jody Hey */
/* This code was developed based on IMa2  2009-2010  Jody Hey, Rasmus Nielsen and Sang Chul Choi*/

/* JH removed RY2 and RY3  7/21/2010 */

#undef GLOBVARS

#include "imamp.h"
#include "update_gtree_common.h"
#include "updateassignment.h"


/*********** LOCAL STUFF **********/

static struct genealogy_weights holdallgweight_t_RY;
static struct genealogy_weights holdgweight_t_RY[MAXLOCI];
static struct probcalc holdallpcalc_t_RY;
// VS - holdpcalc_t_RY 
static struct genealogy_weights holdgweight_theta_t_RY[MAXGROUPS];
static struct genealogy_weights holdgweight_mig_t_RY[MAXGROUPS];
static struct probcalc holdpcalc_theta_t_RY[MAXGROUPS];
static struct probcalc holdpcalc_mig_t_RY[MAXGROUPS];

static int largestsamp;
static int **skipflag;


static double beforesplit (int tnode, double oldt, double newt, double tau_u, double ptime);
static double aftersplit (int tnode, double oldt, double newt, double tau_d, double ptime);
static void storegenealogystats_all_loci (int ci, int mode);

static double forwardRY3 (double ptime, /* double oldt, */double r, double t_u_prior);
static double backwardRY3 (double ptime, /* double newt, */double r,double t_u_prior);
/********* LOCAL FUNCTIONS **************/

void
storegenealogystats_all_loci (int ci, int mode)
{
  static double holdlength[MAXLOCI], holdtlength[MAXLOCI];
  static double holdroottime[MAXLOCI];
  static int holdroot[MAXLOCI];
  static int holdmig[MAXLOCI];
  int li;
  if (mode == 0)
  {
    for (li = 0; li < nloci; li++)
    {
      holdlength[li] = C[ci]->G[li].length;
      holdtlength[li] = C[ci]->G[li].tlength;
      holdroottime[li] = C[ci]->G[li].roottime;
      holdroot[li] = C[ci]->G[li].root;
      holdmig[li] = C[ci]->G[li].mignum;
    }
  }
  else
  {
    for (li = 0; li < nloci; li++)
    {
      C[ci]->G[li].length = holdlength[li];
      C[ci]->G[li].tlength = holdtlength[li];
      C[ci]->G[li].mignum = holdmig[li];
      C[ci]->G[li].roottime = holdroottime[li];
      C[ci]->G[li].root = holdroot[li];
    }
  }
  return;
}                               // storegenealogystats  

double
aftersplit (int tnode, double oldt, double newt, double tau_d, double ptime)
{
  if (tnode == lastperiodnumber - 1)
  {
     return ptime + newt - oldt;
  }
  else
  {
    return tau_d - (tau_d - newt) * (tau_d - ptime) / (tau_d - oldt);
  }
}

double
beforesplit (int tnode, double oldt, double newt, double tau_u, double ptime)
{
  if (tnode == 0)
  {
    return ptime * newt / oldt;
  }
  else
  {
    return tau_u + (ptime - tau_u) * (newt - tau_u) / (oldt - tau_u);
  }
}


/*************GLOBAL FUNCTIONS ******************/


// VS
// changed to initialize gweights and pcalc for groups of loci
void
init_t_RY (void)
{
  int li, j;
  int gp; // VS
  init_genealogy_weights (&holdallgweight_t_RY);
  for (li = 0; li < nloci; li++)
    init_genealogy_weights (&holdgweight_t_RY[li]);
  init_probcalc (&holdallpcalc_t_RY);
  for (largestsamp = 0, j = 0; j < nloci; j++)
    if (largestsamp < L[j].numlines)
      largestsamp = L[j].numlines;
  skipflag = alloc2Dint (nloci, 2 * largestsamp - 1);

  // VS
  // initialize the gweights and pcalc for groups of loci
  for(gp=0; gp<nbgroupsloci_theta; gp++) {
	init_genealogy_weights(&holdgweight_theta_t_RY[gp]);
	init_probcalc (&holdpcalc_theta_t_RY[gp]);
  }
  for(gp=0; gp<nbgroupsloci_mig; gp++) {
 	init_genealogy_weights(&holdgweight_mig_t_RY[gp]);
	init_probcalc (&holdpcalc_mig_t_RY[gp]);
  }

}                               // init_changet_RY

// VS
// changed to free gweights and pcalc for groups of loci
void
free_t_RY (void)
{
  int li;
  int gp; // VS
  free_genealogy_weights (&holdallgweight_t_RY);
  for (li = 0; li < nloci; li++)
  {
    free_genealogy_weights (&holdgweight_t_RY[li]);
  }
  free_probcalc (&holdallpcalc_t_RY);
  orig2d_free2D ((void **) skipflag, nloci);

  // VS 
  // free holdgweight and holdpcalc for each group of loci
  for(gp=0; gp<nbgroupsloci_theta; gp++) {
	free_genealogy_weights(&holdgweight_theta_t_RY[gp]);
	free_probcalc (&holdpcalc_theta_t_RY[gp]);
  }
  for(gp=0; gp<nbgroupsloci_mig; gp++) {
 	free_genealogy_weights(&holdgweight_mig_t_RY[gp]);
	free_probcalc (&holdpcalc_mig_t_RY[gp]);
  }

}                               // free_changet_RY


/* 
Notes on changet_RY()  implements updating of Rannala and Yang (2003)   

This application is pretty much the same as theirs - changing times on a species tree that contains a gene tree. 
The big difference is that IM includes migration.   This means that we have to count migration events and change 
migration times in the same was as we change coalescent times and count how many get changed. 

R&Y also only change coalescent times in populations affected by a changed t.  But because we have migration
there is more entanglement between populations.  It seems best to change all times within an interval that is 
affected by a changing t. 


in R&Y usage
u (upper) means older, deeper in the genealogy
l (lower)  means younger, more recent in the genealogy

In R&Y 
Tu next older split time
Tl  next more recent split time
T - current split time
T*  - proposed time 
t - time of some event 
t* - time of that event after update

If  Tu> t > T
t* = Tu - (Tu-t) (Tu-t*)/(Tu-T)

If Tl< t<= T
t* = Tl + (t-tl) (T*-Tl)/(T-Tl)

nl = number of nodes with Tl< t < T
ml = number of migration events with Tl< t < T

nu = number of nodes with Tu> t > T
mu = number of migration events with Tu> t > T

MH criteria  
p(X|G*,t*)p(G*,t*) (Tu-T*)^(nu+mu) (T*-Tl)^(nl+ml)
---------------------------------------------------
p(X|G,t)p(G,t)      (Tu-T)^(nu+mu) (T-Tl)^(nl+ml)


but this causes confusion with jhey usage in which 
u means upper - more recent. 

so here we use  u  for upper (meaning more recent)
use d for down  (older)

tau current split time
tau*  new value 
tau_d - older node time (i.e. time of next oldest node - deeper in the genealogy)
tau_u - more recent node time (i.e. time of next youngest node - more recent in time)
tau_d  > tau > tau_u

if tau is the first node,  then tau_u = 0. 
if tau is the last node, then tau_d = infinity

for an event at time t where tau_u < t < tau_d

if t > tau  see aftersplit()  t is older than the current split time 
t* = tau_d - (tau_d - tau*)(tau_d - t)/(tau_d-tau)  (A7)

if t <= tau  see beforesplit()
t* = tau_u + (tau* - tau_u)(t - tau_u)/(tau - tau_u) (A8)

m is the number of events moved using A7,  n is the number moved using A8

then Hastings term for the update is:
 tau_d - tau*      tau* - tau_u 
(------------)^m  (------------)^n
 tau-u - tau        tau - tau_u

 For IM,  we use the same except m and n include both includes migation and coalescent events

For edges where downtime < tau_d || uptime > tau_d  are not involved in the update
For any edge that spends zero time in either splitpop  or the ancestral pop, during the tau_u/tau_d interval
it is  possible to not update the coalescent time or migration times of 

The difficulty is that it is possible that an uptime for one edge gets moved because the times on one of its daughter edges got moved. 
This means that for checking about skipping an edge, because it is not involved in any
population associated with the splittin time update
we need to check entire groups of branches that descend from 
an edge that crosses the tau_u boundary. 

use a scoring system for edges  
-1 to ignore because out of bounds above
-2 to ignore because out of bounds below
0 to  deal with, edge is relevant
1  to ignore because not in splitpops or ancestral pop

set up a recursive function  
for an edge that crosses the tau_d line,  check to see if that edge
and all descendent edges that to not cross tau_l  are in the 
populations affected by the splitting time update
If all of those edges are not relevant then 
they all get a skipflag value of 1
If any one of them does spend any time in any of the 
populations involved int the population split
then all of the edges have their skipflag value 
set to 0  

*/

/* let u refer to the more recent time  and d to the older time  */
int
changet_RY1 (int ci, int timeperiod/*, double tv*/)    // after Rannala and Yang (2003)  - rubberband method
//AS: removed tv after discussion with Jody. Refering to ci instead across C[0] occurrences Tue Feb 23 13:48:59 EST 2016
{
  // VS
  int gp; // VS gp
  int gp_theta; // VS group index for theta
  int gp_mig; // VS group index for mig
  double sum_probg_gp; // VS temporary saves prior P(G) as sum across groups

  double metropolishastingsterm, newt, oldt;
  double pdgnew[MAXLOCI + MAXLINKED], pdgnewsum, pdgoldsum, probgnewsum,
    temppdg;
  double t_u_hterm, t_d_hterm, tpw;
  int li, i, j, ecd, ecu, emd, emu, ai, ui;
  double U;
  struct genealogy *G;
  struct edge *gtree;
  double t_d, t_u, t_u_prior, t_d_prior;
  double holdt[MAXPERIODS];

  // VS CHECK - open all check files
  /*fgweights_loci = fopen ("gweights_loci_RY.out", "a");
  fgweights_all = fopen ("gweights_all_RY.out", "a");
  fgweights_theta = fopen ("gweights_theta_RY.out", "a");
  fgweights_mig = fopen ("gweights_mig_RY.out", "a");
  fpcalc_all = fopen ("pcalc_all_RY.out", "a");
  fpcalc_theta = fopen ("pcalc_theta_RY.out", "a");
  fpcalc_mig = fopen ("pcalc_mig_RY.out", "a");*/


  if (assignmentoptions[POPULATIONASSIGNMENTCHECKPOINT] == 1)
  {
    assertgenealogy (ci);
  }

  t_u = (timeperiod == 0) ? 0 : C[ci]->tvals[timeperiod - 1];
  t_d =
    (timeperiod ==
     (lastperiodnumber - 1)) ? TIMEMAX : C[ci]->tvals[timeperiod + 1];
  t_d_prior = DMIN (T[timeperiod].pr.max, t_d);
  t_u_prior = DMAX (T[timeperiod].pr.min, t_u);
  oldt = C[ci]->tvals[timeperiod];
  newt = getnewt (timeperiod, t_u_prior, t_d_prior, oldt, 1);
  
  t_u_hterm = (newt - t_u) / (oldt - t_u);
  if (timeperiod == lastperiodnumber - 1)
  {
    t_d_hterm = 1;
  }
  else
  {
    t_d_hterm = (t_d - newt) / (t_d - oldt);
  }

  // hold info allgweight, allpcalc and tvals
  copy_treeinfo (&holdallgweight_t_RY, &C[ci]->allgweight);  // try turning this off and forcing all recalculations
  copy_probcalc (&holdallpcalc_t_RY, &C[ci]->allpcalc);
  for (i = 0; i < lastperiodnumber; i++)
    holdt[i] = C[ci]->tvals[i];

  // hold sum(P(D|G)) across loci
  pdgoldsum = C[ci]->allpcalc.pdg;
  // VS the allgweight as become obsolete somehow, because it is not used anymore with groups of loci
  setzero_genealogy_weights (&C[ci]->allgweight);
  ecd = ecu = emd = emu = 0;
  pdgnewsum = 0;
  probgnewsum = 0;
  storegenealogystats_all_loci (ci, 0);
  C[ci]->tvals[timeperiod] = newt;
  for (i = 0; i < nurates; i++)
    pdgnew[i] = 0;

  // VS 
  // 1. save on hold gweights and pcalc for groups of loci
  // 2. set the gweights to zero
  for(gp=0; gp<nbgroupsloci_theta; gp++) {
	// 1.1. save gweights for theta
	copy_treeinfo (&holdgweight_theta_t_RY[gp], &C[ci]->groupgweight_theta[gp]);
	// check gweights for theta
	check_gweight_vs (&holdgweight_theta_t_RY[gp], 1);
	// 1.2. save pcalc for groups of loci for theta
	copy_probcalc (&holdpcalc_theta_t_RY[gp], &C[ci]->grouppcalc_theta[gp]);
	// check pcalc
	check_pcalc_groups_vs(&holdpcalc_theta_t_RY[gp], 1);
	// 2. set gweights to zero
	setzero_genealogy_weights (&C[ci]->groupgweight_theta[gp]);
  }
  for(gp=0; gp<nbgroupsloci_mig; gp++) {
	// 1.1. save gweights for mig
	copy_treeinfo (&holdgweight_mig_t_RY[gp], &C[ci]->groupgweight_mig[gp]);
	// check gweights for mig
	check_gweight_vs (&holdgweight_mig_t_RY[gp], 0);
	// 1.2. save pcalc for groups of loci for mig
	copy_probcalc (&holdpcalc_mig_t_RY[gp], &C[ci]->grouppcalc_mig[gp]);
	// check pcalc
	check_pcalc_groups_vs(&holdpcalc_mig_t_RY[gp], 0);
	// 2. set gweights to zero
	setzero_genealogy_weights (&C[ci]->groupgweight_mig[gp]);
  }

  // VS CHECK
  //printf("\nInside RY update. Check holdgweights and holdpcalc theta");
  /*for(gp=0; gp<nbgroupsloci_theta; gp++) {
	fprintf(fgweights_theta, "\nholdgweight_theta_RY gp=%i    ", gp);
    print_gweight_vs_file(&holdgweight_theta_t_RY[gp], 2, fgweights_theta); 
	fprintf(fpcalc_theta, "\nholdpcalc_theta_RY gp=%i    ", gp);
	print_pcalc_groups_vs_file(&holdpcalc_theta_t_RY[gp], fpcalc_theta);
  }
  //printf("\nInside RY update. Check holdgweights and holdpcalc mig");
  for(gp=0; gp<nbgroupsloci_mig; gp++) {
    fprintf(fgweights_mig, "\nholdgweight_mig_RY gp=%i    ", gp);
    print_gweight_vs_file(&holdgweight_mig_t_RY[gp], 2, fgweights_mig); 
	fprintf(fpcalc_mig, "\nholdpcalc_mig_RY gp=%i    ", gp);
	print_pcalc_groups_vs_file(&holdpcalc_mig_t_RY[gp], fpcalc_mig);
  }*/

  for (li = 0; li < nloci; li++)
  {
    G = &(C[ci]->G[li]);
    gtree = G->gtree;

	// VS
	gp_theta = grouploci_theta[ci][li];
	gp_mig = grouploci_mig[ci][li];
//	printf ("%f %f gptheta, gp_mig\n", gp_theta, gp_mig);

    copy_treeinfo (&holdgweight_t_RY[li], &G->gweight);
    for (i = 0; i < L[li].numlines; i++)
    {
      if (gtree[i].down != -1)
      {
        if (gtree[i].time <= oldt && gtree[i].time > t_u)

        {
          //assert (skipflag[li][i] == 0);turn off 9/19/10
          gtree[i].time =
            beforesplit (timeperiod, oldt, newt, /* t_d, */ t_u, gtree[i].time);
          assert (gtree[i].time != newt);
          ecu++;
        }
        else
        {
          if (gtree[i].time > oldt && gtree[i].time < t_d)
          {
           // assert (skipflag[li][i] == 0); turn off 9/19/10
            gtree[i].time =
              aftersplit (timeperiod, oldt, newt, t_d, /* t_u, */ gtree[i].time);
            assert (gtree[i].time != newt);
            ecd++;
          }
          //else  do not change the time
        }
        j = 0;
        while (gtree[i].mig[j].mt > -0.5)
        {
          assert (gtree[i].mig[j].mt < C[ci]->tvals[lastperiodnumber]);//this assertion - AS - has been fixed Tue Feb 23 13:49:36 EST 2016
          if (gtree[i].mig[j].mt <= oldt && gtree[i].mig[j].mt > t_u)
          {
            gtree[i].mig[j].mt =
              beforesplit (timeperiod, oldt, newt, /* t_d, */ t_u,
                           gtree[i].mig[j].mt);
            emu++;
          }
          else
          {
            assert (oldt < C[ci]->tvals[lastperiodnumber]); //same with this one AS - fixed Tue Feb 23 13:49:47 EST 2016
            if (gtree[i].mig[j].mt > oldt && gtree[i].mig[j].mt < t_d)
            {
              gtree[i].mig[j].mt =
                aftersplit (timeperiod, oldt, newt, t_d, /* t_u, */
                            gtree[i].mig[j].mt);
              emd++;
            }
            // else no need to change the time
          }
          j++;
        }
      }
    }
    if (G->roottime <= oldt && G->roottime > t_u
        /* && skipflag[li][G->root] == 0 turn off 9/19/10*/)
      G->roottime =
        beforesplit (timeperiod, oldt, newt, /* t_d, */ t_u, G->roottime);
    else if (G->roottime > oldt && G->roottime < t_d
            /* && skipflag[li][G->root] == 0 turn off 9/19/10*/)
      G->roottime =
        aftersplit (timeperiod, oldt, newt, t_d, /* t_u, */ G->roottime);
    
	// Set the genealogy weights to zero for the locus 
	setzero_genealogy_weights (&G->gweight);
        
    treeweight (ci, li);

    sum_treeinfo (&C[ci]->allgweight, &G->gweight);

	// VS
	// add the gweights of each locus to the weights of the group of loci to which it belongs
	sum_treeinfo_theta_vs(&C[ci]->groupgweight_theta[gp_theta], &G->gweight);
	sum_treeinfo_mig_vs(&C[ci]->groupgweight_mig[gp_mig], &G->gweight);
	// check gweights
	check_gweight_vs (&C[ci]->groupgweight_theta[gp_theta], 1);
	check_gweight_vs (&C[ci]->groupgweight_mig[gp_mig], 0);

	// VS CHECK
	/*printf("\nUpdating locus %i, of gpTheta=%i, gpMig=%i", li, gp_theta, gp_mig);
	printf("\nGweights of locus %i", li);
	print_gweight_vs(&G->gweight,2);
	printf("\nGweights of group theta %i", gp_theta);
	print_gweight_vs(&C[ci]->groupgweight_theta[gp_theta], 0);
	printf("\nGweights of group mig %i", gp_mig);
	print_gweight_vs(&C[ci]->groupgweight_mig[gp_mig], 1);*/
	
	ai = 0;
    ui = L[li].uii[ai];

    switch (L[li].model)
    {
      assert (pdgnew[ui] == 0);
    case HKY:
      if (assignmentoptions[JCMODEL] == 1)
      {
        temppdg = pdgnew[ui] =
          likelihoodJC (ci, li, G->uvals[0]);
      }
      else
      {
        temppdg = pdgnew[ui] =
          likelihoodHKY (ci, li, G->uvals[0], G->kappaval, -1, -1, -1, -1);
      }
      break;
    case INFINITESITES:
      temppdg = pdgnew[ui] = likelihoodIS (ci, li, G->uvals[0]);
      break;
    case STEPWISE:
      temppdg = 0;
      for (; ai < L[li].nlinked; ai++)
      {
        ui = L[li].uii[ai];
        assert (pdgnew[ui] == 0);
        pdgnew[ui] = likelihoodSW (ci, li, ai, G->uvals[ai], 1.0);
        temppdg += pdgnew[ui];
      }
      break;
    case JOINT_IS_SW:
      temppdg = pdgnew[ui] = likelihoodIS (ci, li, G->uvals[0]);
      for (ai = 1; ai < L[li].nlinked; ai++)
      {
        ui = L[li].uii[ai];
        assert (pdgnew[ui] == 0);
        pdgnew[ui] = likelihoodSW (ci, li, ai, G->uvals[ai], 1.0);
        temppdg += pdgnew[ui];
      }
      break;
    }
    pdgnewsum += temppdg;
  }
  tpw = gbeta * (pdgnewsum - pdgoldsum);
  assert (!ODD (ecd));
  assert (!ODD (ecu));
  ecd /= 2;
  ecu /= 2;



  // VS - old version (not used with groups of loci) 
  //integrate_tree_prob (ci, &C[ci]->allgweight, &holdallgweight_t_RY,
  //                     &C[ci]->allpcalc, &holdallpcalc_t_RY, &holdt[0]);   // try enforcing full cacullation
  //tpw += C[ci]->allpcalc.probg - holdallpcalc_t_RY.probg;


    // VS - NEW VERSION WITH GROUPS OF LOCI
    // In new version with groups of loci
    // we will compute the qintegrate for groups of theta and mintegrate for groups of mig
    // and then obtain the sum of the Log(Int(G, param))
	sum_probg_gp=0.0;
	// 1. get qintegrate for theta
	for(gp=0; gp<nbgroupsloci_theta; gp++) {
		// qintegrate for theta
		//AS debug
		//printf("going into integrate_tree_prob_vs....6\n");
		integrate_tree_prob_vs(ci,  &C[ci]->groupgweight_theta[gp], &holdgweight_theta_t_RY[gp], 
                  &C[ci]->grouppcalc_theta[gp], &holdpcalc_theta_t_RY[gp], &holdt[0], gp, 1);
		// check pcalc for theta
		check_pcalc_groups_vs( &C[ci]->grouppcalc_theta[gp], 1);
		// add the probg of group to total sum_probg_gp
		sum_probg_gp +=  C[ci]->grouppcalc_theta[gp].probg;
	}
	//printf("Integrated theta...going into integrating mig\n");
	// 2. get mintegrate for mig
	for(gp=0; gp<nbgroupsloci_mig; gp++) {
		// qintegrate for mig
		
	//	printf("going into integrate_tree_prob_vs....7\n");
		integrate_tree_prob_vs(ci,  &C[ci]->groupgweight_mig[gp], &holdgweight_mig_t_RY[gp], 
					&C[ci]->grouppcalc_mig[gp], &holdpcalc_mig_t_RY[gp], &holdt[0], gp, 0);
		// check pcalc for mig
		check_pcalc_groups_vs( &C[ci]->grouppcalc_mig[gp], 0);
		// add the probg of group to total sum_probg_gp
		sum_probg_gp +=  C[ci]->grouppcalc_mig[gp].probg;
    } 
	//printf("Integrated mig...\n"); 

	// VS
	// sum of the qintegrate and mintegrate to obtain the probg
	// are saved in sum_probg_gp
	C[ci]->allpcalc.probg = sum_probg_gp;
	assert (sum_probg_gp == C[ci]->allpcalc.probg);

	// add the difference in Log(P(G'))-Log(P(G)) to total ratio
	tpw += C[ci]->allpcalc.probg - holdallpcalc_t_RY.probg;


  metropolishastingsterm = beta[ci] * tpw + (ecd + emd) * log (t_d_hterm) +
    (ecu + emu) * log (t_u_hterm);
  //assert(metropolishastingsterm >= -1e200 && metropolishastingsterm < 1e200);
  U = log (uniform ());
  if (U < DMIN(1.0, metropolishastingsterm))  //9/13/2010 
  //if (metropolishastingsterm >= 0.0 || metropolishastingsterm > U)
  {
    for (li = 0; li < nloci; li++)
    {
      C[ci]->G[li].pdg = 0;
      for (ai = 0; ai < L[li].nlinked; ai++)
      {
        C[ci]->G[li].pdg_a[ai] = pdgnew[L[li].uii[ai]];
        C[ci]->G[li].pdg += C[ci]->G[li].pdg_a[ai];
      }
      if (L[li].model == HKY)
      {
        storescalefactors (ci, li);
        copyfraclike (ci, li);
      }
    }
    C[ci]->allpcalc.pdg = pdgnewsum;
    C[ci]->poptree[C[ci]->droppops[timeperiod + 1][0]].time =
      C[ci]->poptree[C[ci]->droppops[timeperiod + 1][1]].time = newt;

    if (assignmentoptions[POPULATIONASSIGNMENTCHECKPOINT] == 1)
    {
      assertgenealogy (ci);
    }
    return 1;
  }
  else
  {
    copy_treeinfo (&C[ci]->allgweight, &holdallgweight_t_RY);
    copy_probcalc (&C[ci]->allpcalc, &holdallpcalc_t_RY);
    assert (pdgoldsum == C[ci]->allpcalc.pdg);
    C[ci]->tvals[timeperiod] = oldt;
    for (li = 0; li < nloci; li++)
    {
      G = &(C[ci]->G[li]);
      gtree = G->gtree;
      storegenealogystats_all_loci (ci, 1);
      copy_treeinfo (&G->gweight, &holdgweight_t_RY[li]);
      for (i = 0; i < L[li].numlines; i++)
      {
        if (gtree[i].down != -1)
        {
          if (gtree[i].time <= newt && gtree[i].time > t_u)
          {
           // assert (skipflag[li][i] == 0); turned off 9/19/10
            gtree[i].time =
              beforesplit (timeperiod, newt, oldt, /* t_d, */ t_u, gtree[i].time);
            //cecu++;
          }

          else
          {
            if (gtree[i].time > newt && gtree[i].time < t_d)
            {
             //assert (skipflag[li][i] == 0); turned off 9/19/10
              gtree[i].time =
                aftersplit (timeperiod, newt, oldt, t_d, /* t_u, */ gtree[i].time);
              //cecl++;
            }
          }
          j = 0;
          while (gtree[i].mig[j].mt > -0.5)
          {
            if (gtree[i].mig[j].mt <= newt && gtree[i].mig[j].mt > t_u)
            {
              gtree[i].mig[j].mt =
                beforesplit (timeperiod, newt, oldt, /* t_d, */
                             t_u, gtree[i].mig[j].mt);
              //cemu++;
            }
            else if (gtree[i].mig[j].mt > newt && gtree[i].mig[j].mt < t_d)
            {
              gtree[i].mig[j].mt =
                aftersplit (timeperiod, newt, oldt, t_d, /* t_u, */
                            gtree[i].mig[j].mt);
              //ceml++;
            }
            j++;
          }
        }
      }
//        assert(fabs(C[ci]->G[li].gtree[  C[ci]->G[li].gtree[C[ci]->G[li].root].up[0]].time - C[ci]->G[li].roottime) < 1e-8);    
    }
    /*    assert(ecu==cecu/2);
       assert(ecd==cecl/2);
       assert(emu==cemu);
       assert(emd==ceml); */
    for (li = 0; li < nloci; li++)
    {
      if (L[li].model == HKY)
        restorescalefactors (ci, li);
      /* have to reset the dlikeA values in the genealogies for stepwise model */
      if (L[li].model == STEPWISE)
        for (ai = 0; ai < L[li].nlinked; ai++)
          likelihoodSW (ci, li, ai, C[ci]->G[li].uvals[ai], 1.0);
      if (L[li].model == JOINT_IS_SW)
        for (ai = 1; ai < L[li].nlinked; ai++)
          likelihoodSW (ci, li, ai, C[ci]->G[li].uvals[ai], 1.0);
      // assert(fabs(C[ci]->G[li].gtree[  C[ci]->G[li].gtree[C[ci]->G[li].root].up[0]].time - C[ci]->G[li].roottime) < 1e-8);    
    }
    if (assignmentoptions[POPULATIONASSIGNMENTCHECKPOINT] == 1)
    {
      assertgenealogy (ci);
    }
	
	// VS
	// restore the gweights and pcalc for groups of loci
	// 1. save gweights and pcalc for groups of loci
	// 2. set the gweights to zero
	for(gp=0; gp<nbgroupsloci_theta; gp++) {
		// 1.1. save gweights for theta
		copy_treeinfo (&C[ci]->groupgweight_theta[gp], &holdgweight_theta_t_RY[gp]);
		// check gweights for theta
		check_gweight_vs (&C[ci]->groupgweight_theta[gp], 1);
		// 1.2. save pcalc for groups of loci for theta
		copy_probcalc (&C[ci]->grouppcalc_theta[gp], &holdpcalc_theta_t_RY[gp]);
		// check pcalc
		check_pcalc_groups_vs(&C[ci]->grouppcalc_theta[gp], 1);
	}
	for(gp=0; gp<nbgroupsloci_mig; gp++) {
		// 1.1. save gweights for mig
		copy_treeinfo (&C[ci]->groupgweight_mig[gp], &holdgweight_mig_t_RY[gp]);
		// check gweights for mig
		check_gweight_vs (&C[ci]->groupgweight_mig[gp], 0);
		// 1.2. save pcalc for groups of loci for mig
		copy_probcalc (&C[ci]->grouppcalc_mig[gp], &holdpcalc_mig_t_RY[gp]);
		// check pcalc
		check_pcalc_groups_vs(&C[ci]->grouppcalc_mig[gp], 0);
	}

  // VS CHECK - Close all check files
  /*fclose(fpcalc_theta);
  fclose(fpcalc_mig);
  fclose(fpcalc_all);

  fclose(fgweights_loci);
  fclose(fgweights_all);
  fclose(fgweights_theta);
  fclose(fgweights_mig);*/



    return 0;
  }
}                               /* changet_RY1 */

