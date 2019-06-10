/* MigSelect 2013-2019 Arun Sethuraman, Vitor Sousa, Jody Hey */
/* This code was developed based on IMa2  2009-2010  Jody Hey, Rasmus Nielsen and Sang Chul Choi*/
#undef GLOBVARS
#include "imamp.h"

/* functions for basic operations on struct genealogy_weights and struct probcalc data structures */

/*********** LOCAL STUFF **********/

/********** GLOBAL FUNCTIONS ***********/

#define ESTARTESIZE 20

/* For Island model, we have no split time, but two periods that are separated
 * by the imaginary split. We have also imaginary population tree where ancestor
 * population resides in the imaginary last period and descendents in the first
 * period.
 * ---------------------------------------------------------------------------
 * numsplittimes: This tells us number of popoulation splits. If we consider a
 * binary population tree, then it would be one less than the number of
 * populations: npops - 1. If we consider Island model of [[npops]] populations,
 * it would be 0. Function [[read_datafile_top_lines]] sets the value.
 * lastperiodnumber: This would be the same as [[numsplittimes]] if we consider
 * a binary population tree. We have the imaginary split time for Island model.
 * It would 1 for Island model. Function [[read_datafile_top_lines]] sets the 
 * value.
 * numtreepops: The number of nodes of a population tree. For Island model, it
 * would be [[npops]] plus one. For a binary population tree, it would be twice
 * the [[npops]] minus one. Function [[read_datafile_top_lines]] sets the value
 * after reading the number of populations from an input file.
 * numpopsizeparams: The number of nodes of a population tree. For Island model,
 * it would be [[npops]]. For a binary population tree, it would be twice the
 * [[npops]] minus one, which is equal to [[numtreepops]].
 */
void
init_genealogy_weights (struct genealogy_weights *gweight)
{
  int i;
  gweight->cc = malloc ((numsplittimes + 1) * sizeof (int *));
  for (i = 0; i < numsplittimes + 1; i++)
    gweight->cc[i] = malloc ((npops - i) * sizeof (int));
  gweight->hcc = malloc ((numsplittimes + 1) * sizeof (double *));
  for (i = 0; i < numsplittimes + 1; i++)
    gweight->hcc[i] = malloc ((npops - i) * sizeof (double));
  gweight->fc = malloc ((numsplittimes + 1) * sizeof (double *));
  for (i = 0; i < numsplittimes + 1; i++)
    gweight->fc[i] = malloc ((npops - i) * sizeof (double));
  
  
  if (modeloptions[NOMIGRATION] == 0)
  {
    gweight->mc = malloc (lastperiodnumber * sizeof (int **));
    gweight->fm = malloc (lastperiodnumber * sizeof (double **));
    for (i = 0; i < lastperiodnumber; i++)
    {
      gweight->mc[i] = alloc2Dint (npops - i, npops - i);
      gweight->fm[i] = orig2d_alloc2Ddouble (npops - i, npops - i);
    }
  }
  setzero_genealogy_weights (gweight);
}

void
setzero_genealogy_weights (struct genealogy_weights *gweight)
{
  int i, j, k;
  for (i = 0; i < numsplittimes + 1; i++)
  {
    for (j = 0; j < npops - i; j++)
    {
      gweight->cc[i][j] = 0;
      gweight->fc[i][j] = 0;
      gweight->hcc[i][j] = 0;
    }
  }
  if (modeloptions[NOMIGRATION] == 0)
  {
    for (k = 0; k < lastperiodnumber; k++)
    {
      for (i = 0; i < npops - k; i++)
      {
        for (j = 0; j < npops - k; j++)
        {
          gweight->mc[k][i][j] = 0;
          gweight->fm[k][i][j] = 0;
        }
      }
    }
  }
  return;
}



void
free_genealogy_weights (struct genealogy_weights *gweight)
{
  int i;

  for (i = 0; i < numsplittimes + 1; i++)
  {
    XFREE (gweight->cc[i]);
    XFREE (gweight->fc[i]);
    XFREE (gweight->hcc[i]);
  }
  XFREE (gweight->cc);
  XFREE (gweight->fc);
  XFREE (gweight->hcc);

  if (modeloptions[NOMIGRATION] == 0)
  {
    for (i = 0; i < lastperiodnumber; i++)
    {
      orig2d_free2D ((void **) gweight->mc[i], npops - i);     /* orig2d_free2D ((void **) gweight->mc[i], npops - i); */
      orig2d_free2D ((void **) gweight->fm[i], npops - i);     /* orig2d_free2D ((void **) gweight->fm[i], npops - i); */
    }
    XFREE (gweight->mc);
    XFREE (gweight->fm);
  }
  return;
}


void
init_probcalc (struct probcalc *pcalc)
{
  int i;
  pcalc->qintegrate = malloc (numpopsizeparams * sizeof (double));
  if (nummigrateparams > 0)
  {
    pcalc->mintegrate = malloc (nummigrateparams * sizeof (double));
  }

  for (i = 0; i < numpopsizeparams; i++)
    pcalc->qintegrate[i] = 0;
  for (i = 0; i < nummigrateparams; i++)
    pcalc->mintegrate[i] = 0;
  pcalc->pdg = 0;
  pcalc->probg = 0;
}

void
free_probcalc (struct probcalc *pcalc)
{
  XFREE (pcalc->qintegrate);
  if (nummigrateparams > 0)
  {
    XFREE (pcalc->mintegrate);
  }
}

// compared memcpy with loops and memcpy did ok
void
copy_treeinfo (struct genealogy_weights *dest, struct genealogy_weights *srce)
{

  int k, i, j;

  for (i = 0; i < numsplittimes + 1; i++)
  {
    for (j = 0; j < npops - i; j++)
    {
      dest->cc[i][j] = srce->cc[i][j];
      dest->hcc[i][j] = srce->hcc[i][j];
      dest->fc[i][j] = srce->fc[i][j];
    }
  }

  if (modeloptions[NOMIGRATION] == 0)
  {
    for (k = 0; k < lastperiodnumber; k++)
    {
      for (i = 0; i < npops - k; i++)
      {
        memcpy (dest->mc[k][i], srce->mc[k][i], (npops - k) * sizeof (int));
        memcpy (dest->fm[k][i], srce->fm[k][i], (npops - k) * sizeof (double));
      }
    }
  }
  return;
}


// copies over the things in probcalc  
// used for genealogy as well as RY and NW updates 
// (NOTE: does not do pdg, i.e. prob(data|genealogy), 
// because this needs to get handled separately 
// in different ways depending on context in update_genealogy(),  changet()  
// and changeu() 
// compared memcpy with loops and memcpy did ok

void
copy_probcalc (struct probcalc *dest, struct probcalc *srce)
{

  memcpy (dest->qintegrate, srce->qintegrate, numpopsizeparams * sizeof (double));
  if (nummigrateparams > 0)
  {
    memcpy (dest->mintegrate, srce->mintegrate, nummigrateparams * sizeof (double));
  }
  /* Note that pdg is not copied!!! */
  dest->probg = srce->probg;
  return;
}

// SUM_TREEINFO
// sums the addto genealogy_weights to addup genealogy_weights
// INPUT:
//	genealogy_weights *addup : pointer to gweights to which add gweights (gweights will change)
//  genealogy_weights *addto : pointer to gweights to add to addup
// RETURN:
//	sums the gweights of addto to addup
// PRE-REQUISITE:
//	addup and addto must be initialized to some values.
void
sum_treeinfo (struct genealogy_weights *addup,
              struct genealogy_weights *addto)
{

  int i, j, k;
  for (i = 0; i < numsplittimes + 1; i++)
  {
    for (j = 0; j < npops - i; j++)
    {
      addup->cc[i][j] += addto->cc[i][j];
      addup->fc[i][j] += addto->fc[i][j];
      assert (addup->cc[i][j] == 0 || addup->fc[i][j] > 0);
      assert (addup->fc[i][j] < DBL_MAX);
      addup->hcc[i][j] += addto->hcc[i][j];

    }
  }

  if (modeloptions[NOMIGRATION] == 0)
  {
    for (k = 0; k < lastperiodnumber; k++)
    {
      for (i = 0; i < npops - k; i++)
      {
        for (j = 0; j < npops - k; j++)

        {
          addup->mc[k][i][j] += addto->mc[k][i][j];
          addup->fm[k][i][j] += addto->fm[k][i][j];
        }
      }
    }
  }

  return;
}

// VS
// SUM_TREEINFO_THETA_VS
// sums the addto genealogy_weights to addup genealogy_weights (only the cc, hc and fc)
// INPUT:
//	genealogy_weights *addup : pointer to gweights to which add gweights (gweights will change)
//  genealogy_weights *addto : pointer to gweights to add to addup
// RETURN:
//	sums the gweights of addto to addup for cc, fc, hc 
//	mc and fm do not change and should be zero
// PRE-REQUISITE:
//	addup and addto must be initialized to some values.
void
sum_treeinfo_theta_vs (struct genealogy_weights *addup,
						struct genealogy_weights *addto)
{
  int i, j;
  for (i = 0; i < numsplittimes + 1; i++)
  {
    for (j = 0; j < npops - i; j++)
    {
      addup->cc[i][j] += addto->cc[i][j];
      addup->fc[i][j] += addto->fc[i][j];
      assert (addup->cc[i][j] == 0 || addup->fc[i][j] > 0);
      assert (addup->fc[i][j] < DBL_MAX);
      addup->hcc[i][j] += addto->hcc[i][j];

    }
  }
  return;
}

// SUM_TREEINFO_MIG_VS
// sums the addto genealogy_weights to addup genealogy_weights
// INPUT:
//	genealogy_weights *addup : pointer to gweights to which add gweights (gweights will change)
//  genealogy_weights *addto : pointer to gweights to add to addup
// RETURN:
//	sums the gweights of addto to addup (only mc and fm)
//	cc, hc and fc are not summed into addup
// PRE-REQUISITE:
//	addup and addto must be initialized to some values.
void
sum_treeinfo_mig_vs (struct genealogy_weights *addup,
					 struct genealogy_weights *addto)
{
  int i, j, k;
  if (modeloptions[NOMIGRATION] == 0)
  {
    for (k = 0; k < lastperiodnumber; k++)
    {
      for (i = 0; i < npops - k; i++)
      {
        for (j = 0; j < npops - k; j++)

        {
          addup->mc[k][i][j] += addto->mc[k][i][j];
          addup->fm[k][i][j] += addto->fm[k][i][j];
        }
      }
    }
  }

  return;
}

// VS
// CHECK_GWEIGHT_VS
// checks if the gweights for theta groups (fc>=0, cc>=0, hc>=0 and mc=0, fm=0)
// checks if the gweights for mig groups (fc==0, cc==0, hc==0 and mc>=0, fm>=0)
// INPUT:
//	genealogy_weights *check : pointer to gweights to evaluate
//  int theta_or_mig : 1 if theta group, 0 if mig group
// RETURN:
//	error printscreen and exits the program if any of the gweights is not correctly defined
// PRE-REQUISITE:
//	check must be initialized to some values
void
check_gweight_vs (struct genealogy_weights *check, int theta_or_mig)
{
  double sum=0.0;

  int i, j, k;

  if(theta_or_mig==1) { // if theta group
	  if (modeloptions[NOMIGRATION] == 0)
	  {
		for (k = 0; k < lastperiodnumber; k++)
		{
		  for (i = 0; i < npops - k; i++)
		  {
			for (j = 0; j < npops - k; j++)

			{
			  sum += check->mc[k][i][j];
			  sum += check->fm[k][i][j];
			}
		  }
		}
	  }
  }
  else { // if migration group
	  
	  for (i = 0; i < numsplittimes + 1; i++)
	  {
		for (j = 0; j < npops - i; j++)
		{
		  sum += check->cc[i][j];
		  sum += check->fc[i][j];
		  assert (check->cc[i][j] == 0 || check->fc[i][j] > 0);
		  assert (check->fc[i][j] < DBL_MAX);
		  sum += check->hcc[i][j];

		}
	  }
  }
  if(sum>0) {
	  printf("\nERROR: gweights different from zero for theta or mig group. Exiting the program.");
	  exit(1);
  }

  return;
}

// VS
// PRINT_GWEIGHT_VS
// prints gweights
// INPUT:
//	genealogy_weights *check : pointer to gweights to evaluate
//  int theta_or_mig : 1 if theta group, 0 if mig group, 2 if both
//  FILE * : pointer to the file to send the output
// RETURN:
//	prints gweights to screen
// PRE-REQUISITE:
//	check must be initialized to some values
void
print_gweight_vs (struct genealogy_weights *check, int theta_or_mig)
{
  int i, j, k;

  if(theta_or_mig==1 || theta_or_mig==2) { // if theta group
	  if (modeloptions[NOMIGRATION] == 0)
	  {
		for (k = 0; k < lastperiodnumber; k++)
		{
		  printf("\nperiod %i\n", k);
		  for (i = 0; i < npops - k; i++)
		  {
			for (j = 0; j < npops - k; j++)
			{
              printf( "mig %i>%i ", i, j);
			  printf( "mc=%i, ", check->mc[k][i][j]);
			  printf( "fm=%g\n", check->fm[k][i][j]);
			}
		  }
		}
	  }
  }
 if(theta_or_mig==0 || theta_or_mig==2) {// if migration group
	  for (i = 0; i < numsplittimes + 1; i++)
	  {
		printf( "\nperiod %i\n", i);
		for (j = 0; j < npops - i; j++)
		{
		  printf( "pop %i, ", j);
		  printf( "cc=%i, ", check->cc[i][j]);
		  printf( "fc=%g, ", check->fc[i][j]);
		  printf( "hcc=%i\n", check->hcc[i][j]);
		}
	  }
  }

  return;
}					/* print_gweight_vs */


// VS
// PRINT_GWEIGHT_VS_FILE
// prints gweights
// INPUT:
//	genealogy_weights *check : pointer to gweights to evaluate
//  int theta_or_mig : 1 if theta group, 0 if mig group, 2 if both
//  FILE * : pointer to the file to send the output
// RETURN:
//	prints gweights to screen
// PRE-REQUISITE:
//	check must be initialized to some values
void
print_gweight_vs_file (struct genealogy_weights *check, int theta_or_mig, FILE *file)
{
  int i, j, k;

 if(theta_or_mig==1 || theta_or_mig==2) {// if theta group
	  for (i = 0; i < numsplittimes + 1; i++)
	  {
		//fprintf(file, "\nperiod %i\n", i);
		for (j = 0; j < npops - i; j++)
		{
		  //fprintf(file, "pop %i, ", j);
		  fprintf(file, "%i ", check->cc[i][j]);
		  fprintf(file, "%g ", check->fc[i][j]);
		  fprintf(file, "%i ", check->hcc[i][j]);
		}
	  }
  }
 if(theta_or_mig==2 || theta_or_mig==2) { // if mig group
	  if (modeloptions[NOMIGRATION] == 0)
	  {
		for (k = 0; k < lastperiodnumber; k++)
		{
		  //fprintf(file, "\nperiod %i\n", k);
		  for (i = 0; i < npops - k; i++)
		  {
			for (j = 0; j < npops - k; j++)
			{
              //fprintf(file, "mig %i>%i ", i, j);
			  fprintf(file, "%i ", check->mc[k][i][j]);
			  fprintf(file, "%g ", check->fm[k][i][j]);
			}
		  }
		}
	  }
  }

  return;
}





#define  subminval 1e-10        // 0 1e-10
void
sum_subtract_treeinfo (struct genealogy_weights *addup,
                       struct genealogy_weights *addtoplus,
                       struct genealogy_weights *addtominus)
{
  int i, j, k;

  for (i = 0; i < numsplittimes + 1; i++)
  {
    for (j = 0; j < npops - i; j++)
    {
      addup->cc[i][j] += addtoplus->cc[i][j] - addtominus->cc[i][j];
      addup->fc[i][j] -= addtominus->fc[i][j];
      addup->fc[i][j] += addtoplus->fc[i][j];
      addup->fc[i][j] = DMAX(0,addup->fc[i][j]);
      addup->hcc[i][j] -= addtominus->hcc[i][j];
      addup->hcc[i][j] += addtoplus->hcc[i][j];
    }
  }

  if (modeloptions[NOMIGRATION] == 0)
  {
    for (k = 0; k < lastperiodnumber; k++)
    {
      for (i = 0; i < npops - k; i++)
      {
        for (j = 0; j < npops - k; j++)
        {
          addup->mc[k][i][j] += addtoplus->mc[k][i][j] - addtominus->mc[k][i][j];
          addup->fm[k][i][j] -= addtominus->fm[k][i][j];
          addup->fm[k][i][j] += addtoplus->fm[k][i][j];
          addup->fm[k][i][j] = DMAX(0,addup->fm[k][i][j]);
        }
      }
    }
  }
  return;
}                               /*sum_subtract_treeinfo */

// VS
// SUM_SUBSTRACY_TREEINFO_VS
// decrease addtominus gweights from addup and add addtoplus to addup
// useful to get the gweights for each group of loci.
// when updating genealogy, only one locus is updated at a time.
// the new gweights for that locus and computed and we need to recompute the gweights
// taking into account all the loci that belong to a given group.
// Instead of re-computing everything again, we just decrease the old and add the new.
// INPUT:
//	genealogy_weights *addup : pointer to gweight that will be saved
//  genealogy_weights *addtoplus : pointer to gweight that will be added to addup
//	genealogy_weights *addtominus : pointer to gweight that will be decreased from addup
//	int theta_or_mig : indicator variable. 1 if THETA params, 0 if MIG params.
// If THETA only cc, hc and fc are taken into account. No operation is done for mc and fm.
// If MIG only mc and fm are taken into account. No operation is done for cc, fc and hc.
void sum_subtract_treeinfo_vs (struct genealogy_weights *addup,
							struct genealogy_weights *addtoplus,
							struct genealogy_weights *addtominus, int theta_or_mig)
{
  int i, j, k;

  for (i = 0; i < numsplittimes + 1; i++)
  {
	  if(theta_or_mig==1) {
		for (j = 0; j < npops - i; j++)
		{
		  addup->cc[i][j] += addtoplus->cc[i][j] - addtominus->cc[i][j];
		  addup->fc[i][j] -= addtominus->fc[i][j];
		  addup->fc[i][j] += addtoplus->fc[i][j];
		  addup->fc[i][j] = DMAX(0,addup->fc[i][j]);
		  addup->hcc[i][j] -= addtominus->hcc[i][j];
		  addup->hcc[i][j] += addtoplus->hcc[i][j];
		}
	  }
  }

  if (modeloptions[NOMIGRATION] == 0)
  {
	  if(theta_or_mig==0) { 
		for (k = 0; k < lastperiodnumber; k++)
		{
		  for (i = 0; i < npops - k; i++)
		  {
			for (j = 0; j < npops - k; j++)
			{
			 
			  addup->mc[k][i][j] += addtoplus->mc[k][i][j] - addtominus->mc[k][i][j];
			  //if (addup->mc[k][i][j] < subminval) addup->mc[k][i][j] = 0;
			  addup->fm[k][i][j] -= addtominus->fm[k][i][j];
			  addup->fm[k][i][j] += addtoplus->fm[k][i][j];
			  addup->fm[k][i][j] = DMAX(0,addup->fm[k][i][j]);
			}
		  }
		}
	  }
  }
  return;
}                               /*sum_subtract_treeinfo_vs */



// VS
// SUM_SUBSTRACY_TREEINFO_ASSIGNLOCI_VS
// decrease addandsubstract gweights from substractfrom 
// and add addandsubstract gweights to addup.
// when updating the assignment of loci to groups of loci
// for migration or theta parameters, we only need to 
// remove the gweight of that locus from old group,
// and add it to the new group.
// INPUT:
//	genealogy_weights *addup : pointer to gweight to where addandsubstract gweights will be added (new group of loci)
//  genealogy_weights *substractfrom : pointer to gweight to where addandsubstract gweights will be substracted (old group of loci)
//	genealogy_weights *addandsubstract : pointer to gweight of locus updated. This will be added to addup and substracted from substractfrom
//	int theta_or_mig : indicator variable. 1 if THETA params, 0 if MIG params.
// If THETA only cc, hc and fc are taken into account. No operation is done for mc and fm.
// If MIG only mc and fm are taken into account. No operation is done for cc, fc and hc.
void sum_subtract_treeinfo_assignloci_vs (struct genealogy_weights *addup,
									  struct genealogy_weights *substractfrom,
									  struct genealogy_weights *addandsubstract, int theta_or_mig)
{
  int i, j, k;

  for (i = 0; i < numsplittimes + 1; i++)
  {
	  if(theta_or_mig==1) {
		for (j = 0; j < npops - i; j++)
		{
		  // Add to addup
		  addup->cc[i][j] += addandsubstract->cc[i][j];
		  addup->fc[i][j] += addandsubstract->fc[i][j];
		  addup->fc[i][j] = DMAX(0,addup->fc[i][j]);
		  addup->hcc[i][j] += addandsubstract->hcc[i][j];

		  // Substract from substractfrom
		  substractfrom->cc[i][j] -= addandsubstract->cc[i][j];
		  substractfrom->fc[i][j] -= addandsubstract->fc[i][j];
		  substractfrom->hcc[i][j] -= addandsubstract->hcc[i][j];
		}
	  }
  }

  if (modeloptions[NOMIGRATION] == 0)
  {
	  if(theta_or_mig==0) { 
		for (k = 0; k < lastperiodnumber; k++)
		{
		  for (i = 0; i < npops - k; i++)
		  {
			for (j = 0; j < npops - k; j++)
			{
			  // Add to addup
			  addup->mc[k][i][j] += addandsubstract->mc[k][i][j];
			  addup->fm[k][i][j] += addandsubstract->fm[k][i][j];
			  addup->fm[k][i][j] = DMAX(0,addup->fm[k][i][j]);

			  // Substract from substractfrom
			  substractfrom->mc[k][i][j] -= addandsubstract->mc[k][i][j];
			  substractfrom->fm[k][i][j] -= addandsubstract->fm[k][i][j];

			}
		  }
		}
	  }
  }
  return;
}                               /*sum_subtract_treeinfo_assignloci_vs */





/* gsampinf is an array of floats that will hold all of the stuff in allgweight  and pcalc for chain 0
nummigrateparams is determined in setup_iparams()
the sequence in the gsampinf array:
the sequence in this array:

	type     |  # values  |  cumulative total at end
    cc	       numpopsizeparams   numpopsizeparams
	fc	       numpopsizeparams   2*numpopsizeparams
	hcc	       numpopsizeparams   3*numpopsizeparams
	mc         nummigrateparams   3*numpopsizeparams + nummigrateparams
	fm         nummigrateparams   3*numpopsizeparams + 2*nummigrateparams
	qintegrate numpopsizeparams   4*numpopsizeparams + 2*nummigrateparams
	mintegrate nummigrateparams   4*numpopsizeparams + 3*nummigrateparams
	pdg             1        4*numpopsizeparams + 3*nummigrateparams +  1 
	probg           1        4*numpopsizeparams + 3*nummigrateparams +  2
	t          numsplittimes 4*numpopsizeparams + 3*nummigrateparams +  numsplittimes + 2
*/

int
calc_gsampinf_length (void)
{
  int i;
  i = 4 * numpopsizeparams;     //cc, hcc, fc,qintegrate 
  if (!modeloptions[NOMIGRATION])
    i += 3 * nummigrateparams;  //mc, fm, mintegrate
  i += numsplittimes;           // times
  i += 2;        // pdg and probg 
  return i;
}

//AS: adding z which indexes the coldchain
void
savegsampinf (float *g, int z)
{

  int i, j, c;
  float f, hc;
  // positions where these types (mc,fc,mintegrate, qintegrate) begin in array g
  struct genealogy_weights *gweight = &C[z]->allgweight;
  for (i = 0; i < numpopsizeparams; i++)
  {
    c = 0;
    f = hc = 0.0;
    for (j = 0; j < itheta[i].wp.n; j++)
    {
      c += (int) gweight->cc[itheta[i].wp.p[j]][itheta[i].wp.r[j]];
      f += (float) gweight->fc[itheta[i].wp.p[j]][itheta[i].wp.r[j]];
      hc += (float) gweight->hcc[itheta[i].wp.p[j]][itheta[i].wp.r[j]];
    }
    g[gsamp_ccp + i] = (float) c;
    assert(g[gsamp_ccp + i] >= 0);
    g[gsamp_fcp + i] = f;
    assert(g[gsamp_fcp + i] >= 0);
    g[gsamp_hccp + i] = hc;
  }
  if (!modeloptions[NOMIGRATION])
  {
    for (i = 0; i < nummigrateparams; i++)
    {
      c = 0;
      f = 0.0;
      for (j = 0; j < imig[i].wp.n; j++)
      {
        c +=
          (int) gweight->mc[imig[i].wp.p[j]][imig[i].wp.r[j]][imig[i].
                                                              wp.c[j]];
        f += (float)
          gweight->fm[imig[i].wp.p[j]][imig[i].wp.r[j]][imig[i].wp.c[j]];
      }
      g[gsamp_mcp + i] = (float) c;
      assert(g[gsamp_mcp + i] >= 0);
      g[gsamp_fmp + i] = f;
      assert(g[gsamp_fmp + i] >= 0);
    }
  }
  for (i = 0; i < numpopsizeparams; i++)
  {
    g[gsamp_qip + i] = (float) C[z]->allpcalc.qintegrate[i];
  }
  for (i = 0; i < nummigrateparams; i++)
  {
    g[gsamp_mip + i] = (float) C[z]->allpcalc.mintegrate[i];
  }
  g[gsamp_pdgp] = (float) C[z]->allpcalc.pdg;
  g[gsamp_probgp] = (float) C[z]->allpcalc.probg;
  for (i = 0; i < numsplittimes; i++)
    g[gsamp_tp + i] = (float) C[z]->tvals[i];
}                               /* savesampinf */


// VS
// SAVEGSAMPINF_VS
// This function saves the gsampinf matrix for each group of loci
// INPUT:
//	float *g : pointer to array (row of matrix where the gsampinf will be saved)
//	int gp : group index (note 0<gp<nbgrouploci_theta or 0<gp<nbgrouploci_mig)
//	int theta_or_mig : 1 if refers to Theta groups, 0 if refers to mig groups
// NOTE: it is important to note that gsampinf in current implementation is a pointer to 3D array.
//	gsampinf[i] is pointing to a 2D array (matrix). Each element i correspond to one group of loci.
//  hence 0<i<(nbgrouploci_theta+nbgrouploci_mig). 
//	Thus, gsampinf[i] refers to theta groups if 0<i<nbgrouploci_theta
//		  gsamping[i] refers to mig groups if (nbgrouploci_theta-1)<i<(nbgrouploci_theta+nbgrouploci_mig)
// AS: adding z, which indexes the coldchain 
void savegsampinf_vs (float *g, int gp, int theta_or_mig, int z) 
{

  int i, j, c;
  float f, hc;
  struct genealogy_weights *gweight;
  
  // positions where these types (mc,fc,mintegrate, qintegrate, sintegrate) begin in array g
  // VS
  // Instead of allgweights it is looking at gweight for each group, which will depend on whether we look to theta or mig groups
  if(theta_or_mig==1) {
	gweight = &C[z]->groupgweight_theta[gp]; 
  } else {
	gweight = &C[z]->groupgweight_mig[gp]; 
  }
  
  // In original IMa2 function only the allgweight is looked up.
  // When different loci may have different demographic parameters
  // we need to also save gsampinfo for each group of locus. 
  // i.e., we need to save the gweights for each group of loci.
  // Done exactly with the same format as gsampinf, 
  // we can use for the functions to obtain marginal posterior distributions
  // for each parameter as if using allgweights.
  for (i = 0; i < numpopsizeparams; i++)
  {
    c = 0;
    f = hc = 0.0;
    for (j = 0; j < itheta[i].wp.n; j++)
    {
      c += (int) gweight->cc[itheta[i].wp.p[j]][itheta[i].wp.r[j]];
      f += (float) gweight->fc[itheta[i].wp.p[j]][itheta[i].wp.r[j]];
      hc += (float) gweight->hcc[itheta[i].wp.p[j]][itheta[i].wp.r[j]];
    }

	// gsamp_ccp, gsamp_fcp and gsamp_hccp are auxiliary variables that save the index of the array
	// since gsampinfo is a matrix, where each row points to a long array, gsamp_ccp etc indicate where is the beggining
	// of a given index
    g[gsamp_ccp + i] = (float) c;
    g[gsamp_fcp + i] = f;
    g[gsamp_hccp + i] = hc;
  } // end of for i=0...numpopsizeparams

  if (!modeloptions[NOMIGRATION])
  {
    for (i = 0; i < nummigrateparams; i++)
    {
      c = 0;
      f = 0.0;
      for (j = 0; j < imig[i].wp.n; j++)
      {
        c += (int) gweight->mc[imig[i].wp.p[j]][imig[i].wp.r[j]][imig[i].wp.c[j]];
        f += (float) gweight->fm[imig[i].wp.p[j]][imig[i].wp.r[j]][imig[i].wp.c[j]];
      }
      g[gsamp_mcp + i] = (float) c;
      g[gsamp_fmp + i] = f;
    }
  } // end if migration

  // go through theta parameters
  for (i = 0; i < numpopsizeparams; i++)
  {
	  // VS - g[gsamp_qip + i] = (float) C[0]->allpcalc.qintegrate[i]; 
	  if(theta_or_mig==1) {
		g[gsamp_qip + i] = (float) C[z]->grouppcalc_theta[gp].qintegrate[i];
	  }
	  else {
		g[gsamp_qip + i] = (float) C[z]->grouppcalc_mig[gp].qintegrate[i];
		// this values should be zero for mig groups
		// because we do not integrate over the theta params for mig groups
	  }
		 
  } // end of theta parameters


  // go through migration parameters
  for (i = 0; i < nummigrateparams; i++)
  {
      // VS g[gsamp_mip + i] = (float) C[0]->allpcalc.mintegrate[i];
	  if(theta_or_mig==1) {
		g[gsamp_mip + i] = (float) C[z]->grouppcalc_theta[gp].mintegrate[i];
		// this values should be zero for theta groups
		// because we do not integrate over the mig params for theta groups
	  }
	  else {
		g[gsamp_mip + i] = (float) C[z]->grouppcalc_mig[gp].mintegrate[i];	
	  }
  }
  g[gsamp_pdgp] = (float) C[z]->allpcalc.pdg;
  g[gsamp_probgp] = (float) C[z]->allpcalc.probg;
  for (i = 0; i < numsplittimes; i++)
    g[gsamp_tp + i] = (float) C[z]->tvals[i];
}                               /* savesampinf_vs */



// VS
// SAVEASSIGNLOCI
// This function saves the assignment of loci into groups array when this is being treated as a random variable in the MCMC (i.e. assignment vectors are updated)
// The structures where the assignment vectors are saved are the global variables int *assignlocisample_mig, and int *assignlocisample_theta;
// INPUT:
//	int *saveassign : pointer to array where the assignment vector is saved (assignlocisample_mig[i] or assignlocisample_theta[i])
//	int *avector : pointer to array of assignment vector (groupsloci_mig[ci] or groupsloci_theta[ci])
// NOTE: grouploci_mig or groupsloci_theta are global variables given as input to this function
void saveassignloci(int *saveassign, int *avector) {

	int i;

	// set the saveassign equal to avector
	for(i=0; i<nloci; i++) {
		saveassign[i] = avector[i];
	}

}








