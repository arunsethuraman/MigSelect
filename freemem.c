/* MigSelect 2013-2019 Arun Sethuraman, Vitor Sousa, Jody Hey */
/* This code was developed based on IMa2  2009-2010  Jody Hey, Rasmus Nielsen and Sang Chul Choi*/
#undef GLOBVARS
#include "imamp.h"


/* misc functions for freeing memory at the end */

extern struct edgemiginfo oldedgemig;
extern struct edgemiginfo oldsismig;
extern struct edgemiginfo newedgemig;
extern struct edgemiginfo newsismig;


/********LOCAL PROTOTYPES ********/
static void free_value_record (struct value_record *v);
static void free_iparam (struct i_param *ip, int n, int m_or_p_or_s);
static void free_locus ();
static void free_T (void);
static void free_lpgpd_v (struct value_record *l);
static void free_migration_counts_times (struct value_record **l);
/********LOCAL FUNCTIONS ********/
void
free_value_record (struct value_record *v)
{
  if (v->do_xyplot)
    XFREE (v->xy);
  if (v->do_trend)
    XFREE (v->trend);
  return;
}

void
free_lpgpd_v (struct value_record *l)
{
  free_value_record (l);
  XFREE (l);
}

void
free_migration_counts_times (struct value_record **l)
{
  int i, j;
  int nummigdirs;
  /* CR 110921.1  set the upper bound for the for loop that free's the
   * histogram storage to the same value used when storage was allocated
   */
  int numhists;   /* number of histograms to free   */
  nummigdirs = 2*(npops-1)*(npops-1);
  numhists = 2 * nummigdirs;  /* 1 for each direction,  locus and loci sum */
  for (j = 0; j < nloci + (nloci > 1); j++)
  {
    for (i = 0; i < numhists; i++)
      free_value_record (&l[j][i]);
    XFREE (l[j]);
  }
  XFREE (l);
}

void
free_T (void)
{
  int i;
  for (i = 0; i < lastperiodnumber; i++)
  {
    XFREE (T[i].upnames);
    XFREE (T[i].upinf);
    free_value_record (T[i].v);
    XFREE (T[i].v);
  }
  XFREE (T);
}                               // free_T



//m_or_p_or_s refers to migration (1) or population size (0) or split (-1)
void
free_iparam (struct i_param *ip, int n, int m_or_p_or_s)
{
  int i, gp; // VS gp
  for (i = 0; i < n; i++) // if itheta-> n=numpopsizeparams, if imig -> n=nummigrateparams
  {
	// VS the array of groups was not being free correctly
	// needed to include a for loop for the number of groups
	if(m_or_p_or_s == 0) { // VS if this is a itheta parameter
		for(gp=0; gp<nbgroupsloci_theta; gp++) {
			XFREE (ip[i].xy[gp]);		
		}
	} 
	if(m_or_p_or_s ==1 ) { // VS if migration parameter
		for(gp=0; gp<nbgroupsloci_mig; gp++) {
			XFREE (ip[i].xy[gp]);		
		}
	}
	XFREE (ip[i].xy); 
	// VS where is the ip[i].pr[gp] free??
	XFREE (ip[i].pr); // VS included the free of the ip[i].pr pointer
	
    if (ip[i].wp.n > 0)
    {
      XFREE (ip[i].wp.p);
      XFREE (ip[i].wp.r);
      if (m_or_p_or_s == 1)
        XFREE (ip[i].wp.c);
    }
  }
  if (m_or_p_or_s >= 0)
    XFREE (ip);
  ip = NULL;
  return;
}                               //free_iparam 

void
free_chainstate_record_updates_and_values (struct
                                           chainstate_record_updates_and_values
                                           *rec, int nrec)
{
  int i, j;
  for (i = 0; i < nrec; i++)
  {
    XFREE ((rec + i)->upnames);
    XFREE ((rec + i)->upinf);
    for (j = 0; j < (rec + i)->num_vals; j++)
      free_value_record (&((rec + i)->v[j]));
    if ((rec + i)->num_vals)
      XFREE ((rec + i)->v);
  }
  XFREE (rec);
}                               //free_chainstate_record_updates_and_values

void
free_a_rec_multichain ()
{
  int ci;
  for (ci = 0; ci < numchains; ci++)
  {
    free (Cupinf[ci].upnames);
    Cupinf[ci].upnames = NULL;
    free (Cupinf[ci].upinf);
    Cupinf[ci].upinf = NULL;
  }
  free (Cupinf);
  Cupinf = NULL;
  return;
}

void
free_locus ()
{
  int li, ai, i;
  for (li = 0; li < nloci; li++)
  {
    if (L[li].model == INFINITESITES || L[li].model == HKY
        || L[li].model == JOINT_IS_SW)
    {
      for (i = 0; i < L[li].numgenes; i++)
        XFREE (L[li].seq[i]);
      XFREE (L[li].seq);
    }
    if (L[li].model == STEPWISE || L[li].model == JOINT_IS_SW)
    {
      if (L[li].model == STEPWISE)
        ai = 0;
      else
        ai = 1;
      for (; ai < L[li].nlinked; ai++)
      {
        XFREE (L[li].A[ai]);
      }

      XFREE (L[li].A);
    }
    if (L[li].model == INFINITESITES || L[li].model == JOINT_IS_SW)
      XFREE (L[li].badsite);
    free_chainstate_record_updates_and_values (L[li].u_rec, L[li].nlinked);
    if (L[li].model == HKY)
      free_chainstate_record_updates_and_values (L[li].kappa_rec, 1);
    if (L[li].model == STEPWISE || L[li].model == JOINT_IS_SW)
      free_chainstate_record_updates_and_values (L[li].A_rec, L[li].nlinked);
    free_chainstate_record_updates_and_values (L[li].g_rec, 1);
    if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
    {
      free_chainstate_record_updates_and_values (L[li].a_rec, 1);

    }
  }
}                               // free_locus


/****** GLOBAL FUNCTIONS  *****/


void
freeanymemory ()
{
  int ci;
  int i;
  int j;
  int li;
  int gp; // VS
  int npnodes;
  int nperiods;

  unsetseeds ();

  if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 1)
  {
    npnodes = npops + 1;
    nperiods = 2;
  }
  else if (npops == 1)
  {
    npnodes = 1;
    nperiods = 1;
  }
  else
  {
    npnodes = 2 * npops - 1;
    nperiods = npops;
    for (i = 0; i < lastperiodnumber; i++)
    {
      free_value_record (T[i].v);
    }
  }

  // VS 5/18/2012
  // free memory of the value record structure for assignment of loci into groups
  if(nbgroupsloci_mig>0 || nbgroupsloci_theta>0) {
	// free the recorded values for the update of migration
	XFREE (assignloci[0].upnames); 
	XFREE (assignloci[0].upinf);
	free_value_record (assignloci[0].v);
	XFREE (assignloci[0].v);
	
	// free the recorded values for the update of theta
	XFREE (assignloci[1].upnames); 
	XFREE (assignloci[1].upinf);
	free_value_record (assignloci[1].v);
	XFREE (assignloci[1].v);

	XFREE (assignloci);
  }


  if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
  {
    free_a_rec_multichain ();
  }

  for (ci = 0; ci < numchains; ci++)
  {
    for (i = 0; i < npnodes; i++)
    {
      free (C[ci]->poptree[i].up);
      C[ci]->poptree[i].up = NULL;
    }
    if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
    {
      XFREE (C[ci]->nasn);
    }

    if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 0)
    {
      XFREE (C[ci]->tvals);
    }

    XFREE (C[ci]->poptree);
    // Free gweights for group theta
    free_genealogy_weights (&(C[ci]->allgweight));
	// VS
	// free memory for groupgweight
	for(gp=0; gp<nbgroupsloci_theta; gp++) {
		free_genealogy_weights(&(C[ci]->groupgweight_theta[gp]));
	}
	XFREE(C[ci]->groupgweight_theta); // VS free pointer for groupgweight_theta
	// Free gweights for group mig
	for(gp=0; gp<nbgroupsloci_mig; gp++) {
		free_genealogy_weights(&(C[ci]->groupgweight_mig[gp]));
	}
	XFREE(C[ci]->groupgweight_mig); // VS free pointer for groupgweight_mig
	
	// Free probcalc
    free_probcalc (&(C[ci]->allpcalc));
	// VS
	// free memory for grouppcalc
	for(gp=0; gp<nbgroupsloci_theta; gp++) {
		free_probcalc(&(C[ci]->grouppcalc_theta[gp]));
	}
	XFREE(C[ci]->grouppcalc_theta);
	for(gp=0; gp<nbgroupsloci_mig; gp++) {
		free_probcalc(&(C[ci]->grouppcalc_mig[gp]));
	}
	XFREE(C[ci]->grouppcalc_mig);
	
    if (C[ci]->plist != NULL)
    {
      for (j = 0; j < nperiods; j++)
      {
        XFREE (C[ci]->plist[j]);
      }
      XFREE (C[ci]->plist);
    }
    for (li = 0; li < nloci; li++)
    {
      free_genealogy_weights (&(C[ci]->G[li].gweight));

	  XFREE (C[ci]->G[li].uvals);
      XFREE (C[ci]->G[li].pdg_a);
      if (L[li].model == HKY)
      {
        XFREE (L[li].mult);
        for (i = L[li].numgenes; i < 2 * L[li].numgenes - 1; i++)
        {
          XFREE (C[ci]->G[li].gtree[i].hkyi.scalefactor);
          XFREE (C[ci]->G[li].gtree[i].hkyi.oldscalefactor);
          for (j = 0; j < L[li].numsites; j++)
          {
            XFREE (C[ci]->G[li].gtree[i].hkyi.frac[j]);
            XFREE (C[ci]->G[li].gtree[i].hkyi.newfrac[j]);
          }
          XFREE (C[ci]->G[li].gtree[i].hkyi.frac);
          XFREE (C[ci]->G[li].gtree[i].hkyi.newfrac);
        }
      }

      for (i = 0; i < L[li].numlines; i++)
      {
        XFREE (C[ci]->G[li].gtree[i].mig);
        if (L[li].model == STEPWISE || L[li].model == JOINT_IS_SW)
        {
          XFREE (C[ci]->G[li].gtree[i].A);
          XFREE (C[ci]->G[li].gtree[i].dlikeA);
        }
        if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
        {
          if (L[li].model == INFINITESITES)
          {
            XFREE (C[ci]->G[li].gtree[i].seq);
          }
        }
      }
      XFREE (C[ci]->G[li].gtree);
      if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
      {
        XFREE (C[ci]->G[li].mut);
      }
    }

    XFREE (C[ci]->G);
    XFREE (C[ci]);
  }
  free_locus ();

/*JH 1/6/09  this will get deleted with these mig values are restored to a large array instead of pointers */
  XFREE (oldedgemig.mtimeavail);
  XFREE (oldedgemig.mp);

  XFREE (oldsismig.mtimeavail);
  XFREE (oldsismig.mp);

  XFREE (newedgemig.mtimeavail);
  XFREE (newedgemig.mp);

  XFREE (newsismig.mtimeavail);
  XFREE (newsismig.mp);

  XFREE (C);
  free_iparam (itheta, numpopsizeparams, 0);
  free_iparam (imig, nummigrateparams, 1);
  itheta = NULL;

  // Add by VS
  XFREE (reltheta);
  XFREE (relmig);

  // Add by VS
  for(i=0; i<numchains; i++) {
	XFREE (grouploci_theta[i]);
	XFREE (grouploci_mig[i]);
  }
  XFREE (grouploci_theta);
  XFREE (grouploci_mig);

  if (npops > 2 && npops <= 5 && outputoptions[PRINTJOINTTEST])
    free_multi_t_arrays ();

  if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 0
      && npops > 1)
  {
    free_t_NW ();
    free_t_RY ();
  }
  free_updategenealogy ();
  // free_updategenealogy_covar (); stopped using 12/15/09 JH
  free_treeweight ();

  free_gtreecommon ();
  free_sumlogk ();
  for (i = 0; i < nomigrationchecklist.n; i++)
  {
    XFREE (nomigrationchecklist.p);
    XFREE (nomigrationchecklist.r);
    XFREE (nomigrationchecklist.c);
  }
  if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 1
      || npops == 1)
  {
    /* no split time */
  }
  else
  {
    free_T ();
  }
  free_lpgpd_v (lpgpd_v);
  if (outputoptions[MIGRATEHIST])
    free_migration_counts_times (migration_counts_times);
  free_autoc_pointers ();

  XFREE (L);
  return;
}                               //freeanymemory
