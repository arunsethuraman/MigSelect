/* MigSelect 2013-2019 Arun Sethuraman, Vitor Sousa, Jody Hey */
/* This code was developed based on IMa2  2009-2010  Jody Hey, Rasmus Nielsen and Sang Chul Choi*/

#define GLOBVARS
#include "imamp.h"
#include "updateassignment.h"
#ifndef WIN32
/* cr 110524.1  config.h part of Auto make process.  File will not be 
 * included if the NOCONFIG define is set in the makefile.  This change 
 * only affects build on platforms other than Windows.
 */
//#ifndef NOCONFIG
//#include <config.h>
//#endif /* NOCONFIG */
#endif

/*AS: MPI definitions */
#ifdef MPI_ENABLED
#include <mpi.h>
#endif

#include <stdlib.h>


extern struct edgemiginfo oldedgemig;
extern struct edgemiginfo oldsismig;
extern struct edgemiginfo newedgemig;
extern struct edgemiginfo newsismig;

#define BURNTRENDSTARTDELAYDEFAULT  10000

static FILE *outfile;
char *outfilename;
/*AS: adding a dummy file name character array to be broadcast */
char outfile_name[FNSIZE];
char *trueassignment;
int gsampinflength; /* updateassignment.c needs to know this. */

/* This is the number of tries of updating assignment, which is equal to number
 * of individuals with unknown origin */
static int snupdatei;              

static char *loadfilebase;
static time_t starttime;
static time_t endtime;
static time_t chainstarttime;
static time_t timer;
static time_t lasttime;
static time_t remained_starttime;
static time_t remained_endtime;
static long burnduration, chainduration;
static int burndone;
static long int burnsteps;
static int burntrendstartdelay;
static int burndurationmode, cdurationmode;
static int genealogiestosave;
static int memforgenealogiessaved = 0;
static int savegenealogyint;
static int recordint;
static int printint;
static double generationtime;
static double scaleumeaninput = 0;
static int swaptries;
static int heatmode;
static char fpstr[50000];       // probably long enough
static int *fpstri;
static char oldoutfilename[FNSIZE];
static double hilocuslike[MAXLOCI];
static int numgenealogyfiles;
static FILE *genealogyinfosavefile;
static char genealogyinfosavefilename[FNSIZE];
static long int recordstep = 0;
static double hilike = -1e20, hiprob = -1e20;
static FILE *checkdonefile;
static double hval1, hval2;
static int trenddoublepoint;
static int trendspot = 0;
static int maxedoutgenealogysave = 0;
static char priorfilename[FNSIZE];
static char mcfwritefilename[FNSIZE];
static char mcfreadfilename[FNSIZE];
static char command_line[1001];
static char heatingterm_str[50], modeloptions_str[50], calcoptions_str[50], outputoptions_str[50];
// VS added a string to save the assignment loci options
static char assignlocioptions_str[50];
static char *infilename;

// AS: adding an infile name character array to be broadcast
static char infile_name[FNSIZE];

static long seed_for_ran1;
static int **migcount, *migfrom, *migto;
static int nummigdirs;
static char migplotfilename[FNSIZE];
static char migrationnamefilename[FNSIZE];
static int migrationnamefrom,migrationnameto;
static FILE *migrationnamefile;
static FILE *migplotfile;
static char nestedmodelfilename[FNSIZE] = "\0";
static char defaultpriorfilename[14]= "imapriors.txt";

// AS: adding swapper and swappee arrays. I don't think I need this, but double check - 9/25/2014
static int *swapper;
static int *swappee;
// AS: in a parallel context, doesn't make sense to swap pointers - so only swapping temperatures
// Functions for doing this will be further added to swapchains.c
static int swapbetasonly = 1;

//AS: currentid indexes the processor on which the function is being called
//AS: if the program is run in serial, currentid = 0, or the head node is the only node running it

/*Local function prototypes  */
static void init_after_start_IMA ();
static void init_IMA ();
static void scan_commandline (int argc, char *argv[], int currentid); //AS: adding currentid
static void print_outputfile_info_string (void);
static void start (int argc, char *argv[], int currentid); //AS: adding currentid
static void qupdate (int currentid, int swapA, int swapB); //AS: adding currentid, and swapping chain ids
// VS modified savegenalogyinfo to allocate memory to gsampinf
static void savegenealogyinfo (int currentid); //AS: adding currentid
// VS added function saveassignlociinfo to save information about the assignment vectors
static void saveassignlociinfo(int currentid); //AS: adding currentid
static void reset_after_burn (int currentid); //AS: adding currentid
static void output_burntrendfile (int currentid); //AS: adding currentid
static int run (int currentid); //AS: adding currentid
static void inctrend (int m, int t, struct value_record *v, double newval);
static void trend_reset (struct value_record *v, int nv);
static void trendrecord (int loadarrayj, int currentid); //AS: adding currentid
static void recordval (struct value_record *v, double val);
static void record_migrations (int z, int currentid); //AS: adding currentid and z - z indexes the coldchain
static void record (int currentid); //AS: adding currentid
static void loadgenealogyvalues (void);
static void callasciicurves (void);
static void callasciitrend (FILE * outfile);
static void printoutput (int currentid); //AS: adding currentid
static void intervaloutput (FILE * outto, int currentid); //AS: adding currentid
static void free_ima_main_stuff ();
static void callprintacceptancerates (FILE * outto);
static void printsteps (FILE * outto, double like);
static void check_to_record (int currentid); //AS: adding currentid
static void record_migration_names();
static void fillswaparrays (int *swapper, int *swappee); //AS: do I need this function? Double check! 9/25/2014

/* SANGCHUL: Mon Jan 12 20:32:41 EST 2009
 * All global variables can be initialized in function [[init_IMA]].
 * All variables must be finalized in function [[free_IMA]].
 * Sometimes we need to initialze some global variables after we call function
 * [[start]]. Function [[init_after_start_IMA]] does this job.
 * */
void
init_IMA ()
{
  if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
  {
    imaAsnInit ();
  }

  lastgenealogysaved = -1;
  return;
}

void
init_after_start_IMA ()
{
#ifdef DEBUG
  int ci;
#endif /* DEBUG */
  char Structurama[FNSIZE];

  if (assignmentoptions[PRINTSTRUCTURAMA] == 1)
    {
      strcpy (Structurama, outfilename);
      strcat (Structurama, ".in");
      IMA_convert_IM2Structurama (Structurama);
      IMA_output_structurama_bat (outfilename);
    }

  if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
  {
    recordassignment_header (outfilename, 0, 0);
#ifdef DEBUG
    if (numchains > 1)
      {
        for (ci = 0; ci < numchains; ci++)
        {
          recordassignment_header (outfilename, ci, 1);
        }
      }
#endif /* DEBUG */
  }
  /* Function IMA_ninds returns the number of individuals with their label being
   * unknown. */
  /* snupdatei = IMA_nindsunknown (); */
  snupdatei = 1;
  return;
}

// VS
// needs to be modified to free the gsampinf 3D array
// AS: I am not entirely sure how I am going to handle the gsampinf 3D array just as yet
// AS: so as of 9/25/2014 - I am leaving this untouched
void
free_ima_main_stuff ()
{
  int i;
  // VS
  int gp; // VS gp
  if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
  {
    imaAsnFree ();
  }
  if (memforgenealogiessaved > 0)
  {
	// VS - added loop through nbgroups of loci
	for(gp=0; gp<nbgroupsloci_theta+nbgroupsloci_mig; gp++) {
		for (i = 0; i < memforgenealogiessaved; i++)
		{
		  XFREE (gsampinf[gp][i]);
		}
		XFREE(gsampinf[gp]);
	}
    XFREE (gsampinf);

	if (runoptions[LOADRUN]==0) {
		// VS free the vectors of assignment of loci into groups
		for(i=0; i < memforgenealogiessaved; i++) {
			XFREE(assignlocisample_mig[i]);
			XFREE(assignlocisample_theta[i]);		
		}
		XFREE(assignlocisample_mig);
		XFREE(assignlocisample_theta);
	}
  }
  XFREE (fpstri);
  if (modeloptions[EXPOMIGRATIONPRIOR] || (runoptions[LOADRUN] && calcoptions[FINDJOINTPOSTERIOR]))
     XFREE(eexpsum);

  if (outputoptions[MIGRATEHIST])
  {
    for (i = 0; i < nloci + (nloci > 1); i++)
    {
      XFREE (migcount[i]);
    }
    XFREE (migcount);
    XFREE (migfrom);
    XFREE (migto);
  }

  if (trueassignment != NULL)
  {
    XFREE (trueassignment);
  }

  if (infilename != NULL)
  {
    XFREE (infilename);
  }

  if (outfilename != NULL)
  {
    XFREE (outfilename);
  }

  if (loadfilebase != NULL)
  {
    XFREE (loadfilebase);
  }
  //AS: do I need this? Double check - as on 9/25/2014
  XFREE(swapper);
  XFREE(swappee);


  freeanymemory ();

}                               //free_ima_main_stuff


#define GENERATIONTIMEDEFAULT   1
#define DEFAULTNUMGENEALOGIES 10000

void
scan_commandline (int argc, char *argv[], int currentid)
{
  static int i, j, k;
  static char ch, ch1;
  static char pstr[256];
  /* option flags that, although values are assigned to the flags, the flag
   * is never used by the code:  Cp, Ep, Wp, Yp.  These variable flags
   * are being left in for completeness
   */
  int Ap, Bp, Cp, Dp, Ep, Fp, Gp, Hfp, Hnp, Hkp, Hap, Hbp, Ip, Jp, Lp, Mp, Op, Pp,
    Qp, Rp, Sp, Tp, Up, Vp, Wp, Yp, Zp;
  // VS
  int Al; // related with assignment of loci options
  int Xp;
  double tempf;
  char *opt;
  const char *sep = ",:";
  int ioption;

  time (&starttime);
  time (&remained_starttime);


  Ap = 0;                       /* Assignment options */
  Bp = 0;                       /* duration of burnin */
  Cp = 0;                       /* calculation options  - flag not used */
  Dp = 0;                       /* number of steps in between genealogy saves */
  Ep = 0;                       /* prior on splitting rate parameter -
                                        flag not used */
  Fp = 0;                       /* name of mcf file */
  Gp = 0;                       /* name of prior file */
  Hfp = 0;                      /* heating model */
  Hnp = 0;                      /* # of chains */
  Hkp = 0;                      /* # of swap attempts */
  Hap = 0;                      /* heat term1 */
  Hbp = 0;                      /* heat term2 */
  Ip = 0;                       /*input file */
  Jp = 0;                       /* used for programmer options */
  Lp = 0;                       /* duration of chain */
  Mp = 0;                       /* migration rate max */
  Op = 0;                       /* output file */
  Pp = 0;                       /* output options */
  Qp = 0;                       /* Theta max scalar */
  Rp = 0;                       /* run options */
  Sp = 0;                       /* random number seed */
  Tp = 0;                       /* Time maximum */
  Up = 0;                       /* generation time in years */
  Vp = 0;                       /* genealogy load file name base */
  Wp = 0;                       /* name of file with nested models -  
                                            flag not used */
  Yp = 0;                       /* mutation rate scalar for loci with mutation rates given in input file - for use with LOADRUN mode  - flag not used */
  Zp = 0;                       /* screen printout frequency */
  Xp = 0;                       /* True Assignment */

  // VS
  Al = 0;

  trueassignment = NULL;
  if (currentid == 0) {
	  printf ("executing program ...Scanning commandline happens only on the head node\n");
  }
  if ((argc == 2 && (char) toupper (argv[1][1]) == 'H') || argc == 1 && currentid == 0)
  {
    printf ("IMa2 Program - copyright 2010 by Jody Hey, Rasmus Nielsen and Sang Chul Choi\n");
    printf ("Release date: %s\n",RELEASE_DATE);
    printf ("This program is run through a command line interface\n");
    printf ("To execute the program, type the program name followed by the necessary command line flags and options \n");
/*  assignment options not included in output for now 4/27/09  jhey 
    printf ("-a  Population assignment options: \n");
    printf ("    1  Invoke population assignment\n");
    printf ("    2  Invoke DNA barcoding\n");
    printf ("    3  Island model\n"); */
/*
    printf ("    0  Turn on check-points\n");
    printf ("    1  Invoke population assignment\n");
    printf ("    2  Relabel update\n");
    printf ("    3  Beerli and Felsenstein update\n");
    printf ("    4  Print info. for DNA Barcoding\n");
    printf ("    5  Island model\n");
    printf ("    6  Local assignment of genes\n");
    printf ("    7  Print input file for Structurama\n");
    printf ("    e.g, -a12 relabel update of assignment with population tree model\n");
    printf ("         -a13 Beerli and Felsenstein update of assignment with population tree model\n");
    printf ("         -a135 Beerli and Felsenstein update of assignment with island model\n");
    printf ("         -a012: turn on checking genealogy integrity additional to option 1,2\n");
    printf ("         -a124: print assignment proportion of a single unknown gene\n");
    printf ("         -a125: relabel update of assignment with island model\n");
    printf ("         -a127: print out the STRUCTURAMA input file additional to option 1,2\n");
    printf ("         -a126: local assignment of genes are allowed\n");
Island model must be in -j option.
1 Population Assignment, BF99 update, -a7 should be turned on
2 DNA Barcoding, NM05 update, -a4 should be turned on
-j6 is no migration. BF99 can be implemented for no migration case.
NM05 can be implemented for no migration case as well.
-a1 : Assignment, Print structurama output, BF99 or NM05 or both
-a2 : DNA Barcoding, Print info. for DNA Barcoding, BF99 or NM05 or both
-a1 -j6 : Assignment with no-migration
-a2 -j6 : DNA Barcoding with no-migration
-a13 : Assignment with island model
-a23 : DNA Barcoding with island model
*/
    printf ("-b  Duration of burn  (MCMC mode only)\n");
    printf ("    - If integer, the number of burnin steps \n");
    printf ("    - If floating point, the time in hours between writing of burntrend file\n");
    printf ("         run continues until file " "IMburn" " is no longer present\n");
    printf ("         in the directory, or if present, does not begin with 'y'\n");
    printf ("-c  Calculation options: \n");
    printf ("    0 Likelihood of data functions return a constant - posterior should equal prior \n");
    printf ("    1 Include ranges on mutation rates as priors on mutation rate scalars\n");
    printf ("    2 Joint posterior density calculations, for LLR tests of nested models use with -w (LOAD-GENEALOGY mode only)\n");
    printf ("    3 Get prior distribution terms from file (requires filename given with -g )\n");
    printf ("-d  Number of steps between genealogy saving (MCMC mode only) (default 100)\n");
    printf ("-f  Name of file with saved Markov chain state generated in previous run - use with -r3\n");
    printf ("-g  Name of file with parameter priors  (requires -c3) default: '%s'\n",defaultpriorfilename);
    //printf("-jh  jh personal options \n");
    //printf("      0  alt data format for SW data -  one data line for each allele in each pop,  1st # is allele length, 2nd is # copies \n");
    printf ("-h  Heating terms (MCMC mode only): \n");
    printf ("  -hf Heating model: l linear (default); t twostep; g geometric\n");
    printf ("  -hn Number of chains \n");
    printf ("  -hk Number of chain swap attempts per step (default = number of chains)\n");
    printf ("  -ha First heating parameter, effect depends on heating model \n");
    printf ("  -hb Second heating parameter, effect depends on heating model  \n");
    printf ("-i  Input file name (no spaces) \n");
    printf ("-j  Model options: \n");
    printf ("    1  Migration only between sister populations (no migration between non-sister populations)\n");
    printf ("    2  One migration parameter for each pair of populations (do not use with -p5)\n");
    printf ("    3  Migration only between sampled populations (ancestral populations have zero migration)\n");
    printf ("    4  Add a non-sampled ghost population to the model \n");
    printf ("    5  Separate population size and migration parameters in each period (lots of parameters) \n");
    printf ("    6  No migration in the model\n");
    printf ("    7  Migration prior follows exponential distribution with mean given by -m or in parameter prior file \n");
    printf ("    8  Each ancestral population size is assumed to identical to that of their largest descendant population\n");
    printf ("    9  One single migration parameter for all pairs of populations (do not use with -p5)\n");
    printf ("-l  Run duration (default: %d genealogies sampled per locus):\n", DEFAULTNUMGENEALOGIES);
    printf ("     If in MCMC mode (i.e. not loading genealogies from a previous run) \n");
    printf ("       - If integer, the number of genealogies to save\n");
    printf ("         This value times -d value sets the # of steps in chain after burnin) \n");
    printf ("	    - If floating point, the time in hours between outputs. \n");
    printf ("         Run continues until file " "IMrun" " is no longer present\n");
    printf ("           in the directory, or if present, does not begin with 'y'\n");
    printf ("     If in load-genealogy mode (i.e. using -r0 to load genealogies from previous run)\n");
    printf ("       - Integer indicates number of genealogies to load from file(s) named with -r\n");
    printf ("-m  Migration prior value (maximum for uniform,  mean if exponential distribution is used \n");
    printf ("-o  Output file name (no spaces) default is 'outfile.txt' \n");
    printf ("-p  Output options: \n");
    printf ("    0 Turn off trend plots in outfile (default is to print trend plots)\n");
    printf ("    1 Turn off plots of marginal curves in outfile (default is to print marginal density plots)\n");
    printf ("    2 Print TMRCA histogram for each genealogy (MCMC mode only)\n");
    printf ("    3 Print histogram of parameters on demographic scales  (requires mutation rate(s) in data file)\n");
    printf ("    4 Print histogram of splitting times divided by prior (do not use with -j0 or when only 2 sampled populations\n");
    printf ("    5 Print estimates and histograms of population migration rate (2NM)\n");
    printf ("    6 Print pairwise probabilities that one parameter is greater than another \n");
    /* CR:110114.2  message text changed */
    printf ("    7 Print histograms of the number of migration events (MCMC mode only)\n");
    printf ("    8 Print joint estimate for splitting times (MCMC mode only, for models with 3, 4 or 5 populations)\n");
    printf ("-q  Maximum for population size parameters (4Nu) \n");
    printf ("-r  Run options \n");
    printf ("    0 LOAD-GENEALOGY Mode - load genealogies from previous run(s); also requires -v \n");
    printf ("    1 Do not save genealogies to a file (default saves sampled genealogies) \n");
    printf ("    2 Save the state of the Markov chain in a file - named with extension .mcf (MCMC mode only)\n");
    printf ("    3 Start run by loading a previously saved *.mcf file; requires -f (data and priors must be the same) \n");
    printf ("    4 Write all mutation related updates rates to stdout during the run (default is to suppress this)\n");
    printf ("    5 Print burntrend file at end of burnin period; use with -b followed by integer (MCMC mode only)\n");
    printf ("-s  Random number seed (default is taken from current time)\n");
    printf ("-t  Maximum time of population splitting\n");
    /*printf ("-t  Maximum time of population splitting ( do not use with -e) \n"); */
    printf ("-u  Generation time in years - for use with -p3 (default is %d) \n", GENERATIONTIMEDEFAULT);
    printf ("-v  Base name (no extension) of *.ti files with genealogy data  (requires use of -r0) \n");
    printf ("-w  Name of file with nested models to be tested (LOAD-GENEALOGY mode only), invokes -c2\n");
/*     printf ("-x  beta for raising the power to likelihood\n"); */
    printf ("-y  Mutation rate scalar for relevant loci - for use with -p3 \n");
    printf ("-z  Number of steps between screen output (default is %d) (MCMC mode only)\n",PRINTINTDEFAULT);
	// VS - assignment of loci into groups of loci options
	printf ("-al Assignment of Loci into group options \n");
	printf ("	 0 Update a_mig \n");
	printf ("	 1 Update a_theta \n");
	printf ("	 2 Update a=a_mig=a_theta \n");
	printf ("	 3 Perform likelihood ratio tests for two populations looking at theta and mig params independently\n");
        printf ("	 4 Do not update assignment of loci into groups. Run with a fixed assignment. \n");
    exit (0);
  }
  else
  {
/*
command line circumstances:
all flags begin with '-'
-most flags are single letter flags
-some are double letter flags: h
// VS some are double letter flags : al 
-some flags are followed by a string or a character
others by a single number (int or float)
others by a string of integers

it is ok to have spaces between a flag and its values 

All flags are followed by at least something
no flag is followed by nothing 
*/
    strcpy (command_line, "");
    strcpy (heatingterm_str,"");
    strcpy (modeloptions_str,"");
    strcpy (calcoptions_str,"");
    strcpy (outputoptions_str,"");
	// VS 
	strcpy (assignlocioptions_str,"");

    for (i = 1; i < argc; i++)
    {
      strcpy (pstr, argv[i]);
      strcat (command_line, " ");
      strcat (command_line, pstr);

      if (strlen (pstr) < 2)
        IM_err (IMERR_COMMANDLINEFORMAT, " one of the command line strings is too short: %s ",pstr);
      if (pstr[0] != '-')
        IM_err (IMERR_COMMANDLINEFORMAT, "command line flag not preceded by '-' : %s", pstr);
      ch = toupper (pstr[1]);
      // VS - get the assignment of loci options
      if(ch == 'A') {
	strcat (assignlocioptions_str, " ");
	strcat (assignlocioptions_str, pstr);
      }
      if (ch == 'H') {
        strcat (heatingterm_str, " ");
        strcat (heatingterm_str, pstr);
      }
      if (ch == 'J') {
        strcat (modeloptions_str, " ");
        strcat (modeloptions_str, pstr);
      }
      if (ch == 'C') {
        strcat (calcoptions_str, " ");
        strcat (calcoptions_str, pstr);
      }
      if (ch == 'P') {
        strcat (outputoptions_str, " ");
        strcat (outputoptions_str, pstr);
      }
      if (ch == 'H') {
        ch1 = toupper (pstr[2]);
      }
      // VS if we have assignment of loci options
      else if (ch == 'A' && strlen (argv[i]) > 2) {
	ch1 = toupper (pstr[2]);
      }
      else {
        ch1 = ' ';
      }
      if (strlen (argv[i]) == 2 || (i < argc - 1 && isdigit (argv[i + 1][0])))  // space separates flag from its number
      {
        i++;
        strcpy (pstr, argv[i]);
        strcat (command_line, " ");
        strcat (command_line, pstr);
      }
      else
      { 
        if ((ch == 'H')) {
          strdelete (pstr, 1, 3);
        }
	// VS if assignment options
	else if (ch == 'A' && strlen (argv[i]) > 2) {
          strdelete (pstr, 1, 3);
        }
        else {
          strdelete (pstr, 1, 2);
        }
      }
      switch ((char) toupper (ch))
      {
      case 'A':
        // VS - If update assignment of loci we have L in ch1
        if(ch1 == 'L') { // the following code assumes that I can only have one option after AL.
            ioption = atoi (&pstr[0]);
            if (ioption < 0 || ioption > 4)
            {
              IM_err (IMERR_COMMANDLINEFORMAT, "option -al %s", pstr);
            }
            assignmentlocioptions[ioption] = 1; //AS: TODO: Do I have to broadcast this?
            Al = 1; // update Al
	   
	}
	else {	// VS - For population assignment options
            j = (int) (strlen (pstr) - 1);
            if (comma_exists (pstr)) {
              for (opt =  (pstr, sep); opt; opt = strtok (NULL, sep))
              {
                ioption = atoi (opt);
                if (ioption < 0 || ioption >= POPULATIONASSIGNMENT_NUMBER)
                {
                  IM_err (IMERR_COMMANDLINEFORMAT, "option -a %s", pstr);
                }
                assignmentoptions[ioption] = 1; //AS: TODO: Do I have to broadcast this?
              }
              Ap = 1;
		
            }
            else
            {
              while (j >= 0)
              {
                ioption = atoi (&pstr[j]);
                if (ioption < 0 || ioption >= POPULATIONASSIGNMENT_NUMBER)
                {
                  IM_err (IMERR_COMMANDLINEFORMAT, "option -a %s", pstr);
                }
                assignmentoptions[ioption] = 1;
                pstr[j] = '\0';
                j--;
                Ap = 1;
		
	
              }
            }
      }
      break;
      case 'B':
        tempf = atof (&pstr[0]);
        /* check to see if the value is floating point, in which case treat it as being in fractions of an hour  and convert to seconds */
        if (strchr (pstr, '.'))
        {
          burnduration = (int) (3600 * tempf);
          burndurationmode = TIMEINF;
          time (&lasttime);
          runoptions[PRINTBURNTREND] = 1;
        }
        else
        {
          burnduration = (int) tempf;
          burndurationmode = TIMESTEPS;
        }
	#ifdef MPI_ENABLED
		MPI_Bcast(&burnduration, 1, MPI_LONG, 0, MPI_COMM_WORLD);
		MPI_Bcast(&burndurationmode, 1, MPI_INT, 0, MPI_COMM_WORLD);
	#endif
	//AS: TODO: error catching - MPI_Bcast returns an int - so check with that for errors


        Bp = 1;
        break;
      case 'C':
        j = (int) (strlen (pstr) - 1);
        while (j >= 0)
        {
          if (!isdigit (pstr[j]))
            IM_err (IMERR_COMMANDLINEFORMAT, "calculation option flag -c should be followed by a digit: %s",pstr);
          calcoptions[atoi (&pstr[j])] = 1;
          pstr[j] = '\0';
          j--;
        }
        break;
      case 'D':
        savegenealogyint = atoi (&pstr[0]);
	#ifdef MPI_ENABLED
		MPI_Bcast(&savegenealogyint, 1, MPI_INT, 0, MPI_COMM_WORLD);
	#endif //AS: TODO: error catching
        Dp = 1;
        break;
      case 'F':
        strcpy (mcfreadfilename, pstr);
        put_spaces_in_filepaths(mcfreadfilename);
	#ifdef MPI_ENABLED
		MPI_Bcast(&mcfreadfilename, 500*sizeof(char), MPI_CHAR, 0, MPI_COMM_WORLD);
	#endif //AS: TODO: error catching

        Fp = 1;
        break;
      case 'G':
        strcpy (priorfilename, pstr);
        put_spaces_in_filepaths(priorfilename);
	#ifdef MPI_ENABLED
		MPI_Bcast(&priorfilename, 500*sizeof(char), MPI_CHAR, 0, MPI_COMM_WORLD);
	#endif //AS: TODO: error catching

        Gp = 1;
        break;
      case 'H':
        switch ((char) toupper (ch1))
        {
        case 'A':
          hval1 = atof (&pstr[0]);
          Hap = 1;
          break;
        case 'B':
          hval2 = atof (&pstr[0]);
          Hbp = 1;
          break;
        case 'N':
          numchains = atoi (&pstr[0]);
          Hnp = 1;
          break;
        case 'K':
          swaptries = atoi (&pstr[0]);
          Hkp = 1;
          break;
        case 'F':
          Hfp = 1;
          switch ((char) toupper (pstr[0]))
          {
          case 'T':
            heatmode = HTWOSTEP;
            break;
          case 'G':
            heatmode = HGEOMETRIC;
            break;
          default:
            heatmode = HLINEAR;
            break;
          }
          break;
        default:
          IM_err (IMERR_COMMANDLINEHEATINGTERMS, "mistake in use of -h flag : %s", pstr);
        }
	#ifdef MPI_ENABLED
		MPI_Bcast(&hval1, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&hval2, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&numchains, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&swaptries, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&heatmode, 1, MPI_INT, 0, MPI_COMM_WORLD);
	#endif //AS: TODO error catching

        break;
      case 'I':
        infilename = strdup (pstr);
        put_spaces_in_filepaths(infilename);
	//AS: copying this out so I can broadcast it instead
	strncpy(infile_name, infilename, sizeof(infilename));
	#ifdef MPI_ENABLED
		MPI_Bcast(&infile_name, 500*sizeof(char), MPI_CHAR, 0, MPI_COMM_WORLD);
	#endif //AS: TODO: error catching
        Ip = 1;
        break;
      case 'J':
        Jp = 0;
        if (!(toupper(pstr[0]) == 'H'))
        {
          j = (int) (strlen (pstr) - 1);
          while (j >= 0)
          {
            if (!isdigit (pstr[j]))
              IM_err (IMERR_COMMANDLINEFORMAT, "model option flag -j should be followed by a digit: %s", pstr);
            modeloptions[atoi (&pstr[j])] = 1;
            pstr[j] = '\0';
            j--;
          }
        }
        else
        {
          j = (int) (strlen (pstr) - 1);
          if (isdigit(pstr[1])  && pstr[1] == '1')
          {
            migrationnamefrom = atoi (&pstr[2]);
            migrationnameto = atoi (&pstr[3]);
            strcpy(migrationnamefilename, &pstr[4]);
            put_spaces_in_filepaths(migrationnamefilename);
            jheyoptions[WRITEMIGRATIONNAME] = 1;
          }
          else
          {
            while (j >= 0)
            {
              if (isdigit(pstr[j]))
              {
                jheyoptions[atoi (&pstr[j])] = 1;
              }
              pstr[j] = '\0';
              j--;
            }
          }
        };
        break;
      case 'L':
        tempf = atof (&pstr[0]);
        /* check to see if the value is floating point, in which case treat it as being in fractions of an hour  and convert to seconds */
        if (strchr (pstr, '.'))
        {
          chainduration = (int) (3600 * tempf);
          cdurationmode = TIMEINF;
          genealogiestosave = -1;
        }
        else
        {
          genealogiestosave = (int) tempf;
          cdurationmode = TIMESTEPS;
        }
        Lp = 1;
	#ifdef MPI_ENABLED
		MPI_Bcast(&chainduration, 1, MPI_LONG, 0, MPI_COMM_WORLD);
		MPI_Bcast(&cdurationmode, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&genealogiestosave, 1, MPI_INT, 0, MPI_COMM_WORLD);
	#endif //AS: TODO error catching
        break;
      case 'M':
        mprior = (double) atof (&pstr[0]);
        if (mprior == 0)
          modeloptions[NOMIGRATION] = 1;
	#ifdef MPI_ENABLED
		MPI_Bcast(&mprior, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
	#endif //AS: TODO error catching
        Mp = 1;
        break;
      case 'O':
        outfilename = strdup (pstr);
        put_spaces_in_filepaths(outfilename);
        Op = 1;
        break;
      case 'P':
        Pp = 1;
        j = (int) (strlen (pstr) - 1);
        while (j >= 0)
        {
          if (!isdigit (pstr[j]))
          {
            IM_err (IMERR_COMMANDLINEFORMAT, "print option flag -p should be followed by a digit : %s",pstr);
          }
          k = atoi (&pstr[j]);
          outputoptions[k] = 1;
          pstr[j] = '\0';
          j--;
        }
        break;
      case 'Q':
        thetaprior = (double) atof (&pstr[0]);
	#ifdef MPI_ENABLED
		MPI_Bcast(&thetaprior, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
	#endif //AS: TODO error catching

        Qp = 1;
        break;
      case 'R':
        Rp = 1;
        j = (int) (strlen (pstr) - 1);
        while (j >= 0)
        {
          if (!isdigit (pstr[j]))
            IM_err (IMERR_COMMANDLINEFORMAT, "run option flag -r should be followed by a digit : %s ", pstr);
          k = atoi (&pstr[j]);
          runoptions[k] = 1;
          pstr[j] = '\0';
          j--;
        }
        break;

      case 'S':
        seed_for_ran1 = atoi (&pstr[0]);
        if (!seed_for_ran1)
          seed_for_ran1 = 1;
        Sp = 1;
        break;
      case 'T':
        tprior = (double) atof (&pstr[0]);
	#ifdef MPI_ENABLED
		MPI_Bcast(&tprior, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
	#endif //AS: TODO error catching

        Tp = 1;
        break;
      case 'U':
        generationtime = atof (&pstr[0]);
	#ifdef MPI_ENABLED
		MPI_Bcast(&generationtime, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
	#endif //AS: TODO error catching

        Up = 1;
        break;
      case 'V':
        loadfilebase = strdup (pstr);
        put_spaces_in_filepaths(loadfilebase);
        Vp = 1;
        break;
      case 'W':
        strcpy (nestedmodelfilename, pstr);
        put_spaces_in_filepaths(nestedmodelfilename);
        calcoptions[FINDJOINTPOSTERIOR] = 1;
        break;
      case 'Y':
        scaleumeaninput = atof (&pstr[0]);
	#ifdef MPI_ENABLED
		MPI_Bcast(&scaleumeaninput, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
	#endif //AS: TODO add error catching

        break;
      case 'Z':
        printint = atoi (&pstr[0]);
	#ifdef MPI_ENABLED
		MPI_Bcast(&printint, 1, MPI_INT, 0, MPI_COMM_WORLD);
	#endif //AS: TODO add error catching
        Zp = 1;
        break;
      case 'X':
        if (strchr (pstr, '.') == NULL)
        {
          trueassignment = strdup (pstr);
          gbeta = 1.0;
        }
        else
        {
          gbeta = atof (&pstr[0]);
          trueassignment = NULL;
        }
        Xp = 1;
	#ifdef MPI_ENABLED
		MPI_Bcast(&gbeta, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
	#endif //AS: TODO error catching
        break;
      default:
        IM_err (IMERR_COMMANDLINEFORMAT, &ch);
      }
    }
  }

  /* Check if command line options are compatible with single population. */
  if (infilename==0)
    IM_err (IMERR_READFILEOPENFAIL,  "pointer to input file not set,  check -i on command line");
  npops = imaInfileNpops (infilename);
	#ifdef MPI_ENABLED
		MPI_Bcast(&npops, 1, MPI_INT, 0, MPI_COMM_WORLD);
	#endif //AS: TODO error catching

  if (npops < 1 || npops > 10)
  {
    IM_err (IMERR_COMMANDLINEFORMAT, "Number of populations must be nonnegative and less than 10");
  }
  if (npops == 1)
  {
    if (modeloptions[NOMIGRATION] == 0)
      modeloptions[NOMIGRATION] = 1;
    if (Mp == 1)
    {
      IM_err (IMERR_COMMANDLINEFORMAT, "model option [single population] should not go with -m");
    }
    if (Tp == 1)
    {
      IM_err (IMERR_COMMANDLINEFORMAT, "model option [single population] should not go with -t");
    }
  }
  
  // VS - assignment of loci options
  if(Al == 1) {
      // Allowed command lines
      // -al0
      // -al1
      // -al2
      // -al3
      // -al4
      if(assignmentlocioptions[NOTUPDATEASSIGLOCI]==1 
              && (assignmentlocioptions[UPDATEASSIGNMIG]==1 
                  || assignmentlocioptions[UPDATEASSIGNTHETA]==1 
                  || assignmentlocioptions[UPDATEASSIGNMIGTHETA]==1)) 
      {
        IM_errloc (AT, "option -al is being misused.");
      }
  
  }

  if (Ap == 1)
  {
    /* Allowed command line options for assignment work.
     * -a0, -a01, -a02, -a03, -a013, -023,
     * -a1, -a13
     * -a2, -a23
     * -a3,
     */
    if (assignmentoptions[POPULATIONASSIGNMENT] == 1
        || assignmentoptions[POPULATIONASSIGNMENTRELABEL] == 1
        || assignmentoptions[POPULATIONASSIGNMENTBF] == 1
        || assignmentoptions[POPULATIONASSIGNMENTASSIGNED] == 1
        || assignmentoptions[PRINTSTRUCTURAMA] == 1
        || assignmentoptions[JCMODEL] == 1
        || assignmentoptions[POPULATIONTREEWCH] == 1)
    {
      IM_errloc (AT, "option -a is being misused.");
    }

    if (assignmentoptions[POPULATIONASSIGNMENTASN] == 1
        && assignmentoptions[POPULATIONASSIGNMENTBAR] == 0)
    {
      assignmentoptions[POPULATIONASSIGNMENT] = 1;
      assignmentoptions[POPULATIONASSIGNMENTBF] = 1;
      assignmentoptions[PRINTSTRUCTURAMA] = 1;
    }
    else if (assignmentoptions[POPULATIONASSIGNMENTASN] == 0
             && assignmentoptions[POPULATIONASSIGNMENTBAR] == 1)
    {
      assignmentoptions[POPULATIONASSIGNMENT] = 1;
      assignmentoptions[POPULATIONASSIGNMENTRELABEL] = 1;
      assignmentoptions[POPULATIONASSIGNMENTASSIGNED] = 1;
    } 
    else if (assignmentoptions[POPULATIONASSIGNMENTCHECKPOINT] == 1)
    {
      /* MCMC State check is allowed without assignment. */
    }
    else if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 1)
    {
      /* Island model is allowed without assignment. */
    }
    else
    {
      IM_errloc (AT, "option -a  is being misused.");
    }
  }

  if (Xp == 0)
  {
    gbeta = 1.0;
	#ifdef MPI_ENABLED
		MPI_Bcast(&gbeta, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
	#endif //AS: TODO error catching
  }

  if (!Ip)
  {
    IM_err (IMERR_MISSINGCOMMANDINFO, " No data file given on command line");
  }
  if (!Lp && !runoptions[LOADRUN])
  {
    genealogiestosave = (int) DEFAULTNUMGENEALOGIES;
    cdurationmode = TIMESTEPS;
	#ifdef MPI_ENABLED
		MPI_Bcast(&genealogiestosave, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&cdurationmode, 1, MPI_INT, 0, MPI_COMM_WORLD);
	#endif //AS: TODO error catcing
  }
  if (runoptions[LOADRUN] || modeloptions[NOMIGRATION])
    outputoptions[MIGRATEHIST] = 0;
  if (!Op)
    strcpy (outfilename, "outfile.txt");
#ifdef _DEBUG
  checkoutfileclosed (&outfile, outfilename);   // just make sure that outfilename does not name a file that is already opened
#endif
  if (runoptions[LOADRUN] && !Vp)
  {
    IM_err (IMERR_MISSINGCOMMANDINFO, " -r0 invoked without -v information, i.e. no base name for files containing genealogys was given on the command line");
  }
  if (!Hnp || runoptions[LOADRUN])
  {
    numchains = DEFCHAINS;      /* default value */
	#ifdef MPI_ENABLED
		MPI_Bcast(&numchains, 1, MPI_INT, 0, MPI_COMM_WORLD);
	#endif //AS: TODO error catchign
  }
  else
  {
    if ((heatmode == HGEOMETRIC && numchains < 4)
        || (heatmode == HTWOSTEP && numchains < 3) )
    {
      IM_err (IMERR_COMMANDLINEHEATINGTERMS, "too few chains specified in heating model");
    }
    if (!Hfp)
    {
      heatmode = HLINEAR;
      if (!Hap)
        hval1 = 0.05;           /* default value */
	#ifdef MPI_ENABLED
		MPI_Bcast(&hval1, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
	#endif //AS: TODO error catching
    }
    else
    {
      if (heatmode > HLINEAR)
      {
        if (heatmode == HTWOSTEP)
        {
          if (!Hap) {
            hval1 = 0.05;
		#ifdef MPI_ENABLED
			MPI_Bcast(&hval1, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
		#endif //AS: TODO error catchign
		}
		
          if (!Hbp) {
            hval2 = 2;
		#ifdef MPI_ENABLED
			MPI_Bcast(&hval2, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
		#endif //AS: TODO error catching
	}
        }
        if (heatmode == HGEOMETRIC)
        {
          if (!Hap) {
            hval1 = 0.95;  // default value 
		#ifdef MPI_ENABLED
			MPI_Bcast(&hval1, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
		#endif //AS: TODO error catching
	}
          if (!Hbp) {
            hval2 = 0.8;
		#ifdef MPI_ENABLED
			MPI_Bcast(&hval2, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
		#endif //AS: TODO error catchign
	}
          // if (hval1 > 1.0)  //6/11/2010 JH  stopped this,  it turns out numbers slightly higher than 1 can be useful when the are 
            // a large number of chains
           // IM_err (IMERR_COMMANDLINEHEATINGTERMS, "ha commandline term is out of range, should be <= 1");
          if (hval1 > 1.1) //6/11/2010  JH  added this to avoid values much larger than 1
            IM_err (IMERR_COMMANDLINEHEATINGTERMS, "for geometric heating it is not useful to have the ha term be greater than 1.1");
          if (hval1 < 0.9)
            IM_err (IMERR_COMMANDLINEHEATINGTERMS, "for geometric heating it is not useful to have the ha term be less than 0.9");
          if (hval2 >= 1.0|| hval2 <= 0.0)
            IM_err (IMERR_COMMANDLINEHEATINGTERMS, "hb commandline term is out of range, should be < 1 and > 0)");
        }
      }
      else if (!Hap) {
        hval1 = 0.05;           /* default value */
	#ifdef MPI_ENABLED
		MPI_Bcast(&hval1, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
	#endif //AS: TODO error catching
	}
    }
  }
  recordint = RECORDINTDEFAULT;
	#ifdef MPI_ENABLED
		MPI_Bcast(&recordint, 1, MPI_INT, 0, MPI_COMM_WORLD);
	#endif //AS: TODO error catching
  if (!Up) {
    generationtime = GENERATIONTIMEDEFAULT;
	#ifdef MPI_ENABLED
		MPI_Bcast(&generationtime, 1, MPI_INT, 0, MPI_COMM_WORLD);
	#endif //AS: TODO error catching
}
  if (!Dp) {
    savegenealogyint = SAVEGENEALOGYDEFAULT;
	#ifdef MPI_ENABLED
		MPI_Bcast(&savegenealogyint, 1, MPI_INT, 0, MPI_COMM_WORLD);
	#endif //AS: TODO error catching
 }
  if (!Zp) {
    printint = PRINTINTDEFAULT;
	#ifdef MPI_ENABLED
		MPI_Bcast(&printint, 1, MPI_INT, 0, MPI_COMM_WORLD);
	#endif //AS: TODO error catching
 }
  if (!Bp && !runoptions[LOADRUN])
  {
    IM_err (IMERR_MISSINGCOMMANDINFO,
            "No burn duration information given on command line (use -b)");
  }
  if (modeloptions[NOMIGRATION] && !Mp)
  {
    mprior = 0;
	#ifdef MPI_ENABLED
		MPI_Bcast(&mprior, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
	#endif //AS: TODO error catching
  }
  if (modeloptions[NOMIGRATION] && Mp && mprior > 0)
  {
    IM_err (IMERR_COMMANDLINEINCOMPAT, " Conflicting command line arguments, no migration set but migration prior > 0 : %lf",mprior);
  }
  if (Gp && !calcoptions[USEPRIORFILE])
  {
    calcoptions[USEPRIORFILE] = 1;
  }
  if (!Tp && !assignmentoptions[POPULATIONASSIGNMENTINFINITE]
      && npops > 1
      && !calcoptions[USEPRIORFILE])
  {
    IM_err (IMERR_MISSINGCOMMANDINFO,
            " No prior information provided for population splitting times (-t)");
  }

  if (!Mp && modeloptions[NOMIGRATION] != 1 && !calcoptions[USEPRIORFILE] && npops > 1)
  {
    IM_err (IMERR_MISSINGCOMMANDINFO,
            " No information provided for maximum value for migration parameter (-m)");
  }
  if (!Qp && !calcoptions[USEPRIORFILE])
  {
    IM_err (IMERR_MISSINGCOMMANDINFO,
            " No information provided for maximum value of 4Nu parameters");
  }

  if (runoptions[LOADMCSTATE] && !Fp)
  {
    IM_err (IMERR_MISSINGCOMMANDINFO,
            "  -r3 invoked without -f information, i.e. no filename given for markov chain state file,  for loading state of markov chain ");
  }
  if (modeloptions[PARAMETERSBYPERIOD] && modeloptions[ANCESTRALPOPSIZESSETTOLARGESTSAMPLED])
    IM_err (IMERR_COMMANDLINEINCOMPAT,  " model options in conflict: -%d with -%d",PARAMETERSBYPERIOD,ANCESTRALPOPSIZESSETTOLARGESTSAMPLED);
  if (!Hkp) // default is to swap a once step or once for every 10 chains each step whichver is greater
    //swaptries = numchains;
    swaptries = IMAX(1,numprocesses * numchains/10);
  else
  {
    if (numprocesses * numchains > 1 && swaptries > numprocesses * numchains * (numprocesses * numchains - 1) / 2)
      swaptries = numprocesses * numchains * (numprocesses * numchains - 1) / 2;
  }
	#ifdef MPI_ENABLED
		MPI_Bcast(&swaptries, 1, MPI_INT, 0, MPI_COMM_WORLD);
	#endif //AS: TODO catch errors

  if (runoptions[PRINTBURNTREND])
  {
    if (runoptions[LOADMCSTATE])
      burntrendstartdelay = 0;
    else
      burntrendstartdelay = BURNTRENDSTARTDELAYDEFAULT;
  }
	#ifdef MPI_ENABLED
		MPI_Bcast(&burntrendstartdelay, 1, MPI_INT, 0, MPI_COMM_WORLD);
	#endif //AS: TODO error catching
  if (!Sp)
    seed_for_ran1 = (long) time (NULL);
  if (strcmp (infilename, outfilename) == 0)
  {
    IM_err (IMERR_COMMANDLINEINCOMPAT, " Input and output file names are identical");
  }
  if (cdurationmode == TIMESTEPS)
  {
    chainduration = genealogiestosave * savegenealogyint;
	#ifdef MPI_ENABLED
		MPI_Bcast(&chainduration, 1, MPI_LONG, 0, MPI_COMM_WORLD);
	#endif //AS: TODO error catching
  }
  if (jheyoptions[WRITEMIGRATIONNAME])
  {
    migrationnamefile = fopen (migrationnamefilename, "w");
  }
  if (( modeloptions[NOMIGBETWEENNONSISTERS] ||
        modeloptions[SINGLEMIGRATIONBOTHDIRECTIONS] ||
        modeloptions[MIGRATIONBETWEENSAMPLED] || 
        modeloptions[ADDGHOSTPOP] || 
        modeloptions[PARAMETERSBYPERIOD] ||
       /* modeloptions[EXPOMIGRATIONPRIOR] ||  */
        npops == 1 ||
        modeloptions[NOMIGRATION]) && calcoptions[USEPRIORFILE])
  {
    IM_err (IMERR_COMMANDLINEINCOMPAT,"incompatibility between a model option and the use of a file with parameter priors");
  }

  if (( modeloptions[SINGLEMIGRATIONBOTHDIRECTIONS] ||
        npops == 1 ||
        modeloptions[NOMIGRATION]) &&  modeloptions[ONEMIGRATIONPARAMETER])
  {
    IM_err (IMERR_COMMANDLINEINCOMPAT,"incompatibility between a model option and the use of a single migration parameter");
  }
  if (calcoptions[FINDJOINTPOSTERIOR]==1 && runoptions[LOADRUN]==0)
  {
    IM_err (IMERR_COMMANDLINEINCOMPAT," finding joint posterior and tests of nested models requires L mode");
  }

  /* IMa stops its running if there is any conflicted options. */
  if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 1)
  {
    if (modeloptions[NOMIGRATION] == 1)
    {
      IM_err (IMERR_COMMANDLINEINCOMPAT, 
              "Island model must be allowed for migration events: do not use -j%d with -a%d", 
              NOMIGRATION, 
              POPULATIONASSIGNMENTINFINITE);
    }
    if ( modeloptions[NOMIGBETWEENNONSISTERS] == 1
        || modeloptions[PARAMETERSBYPERIOD] == 1)
    {
      IM_err (IMERR_COMMANDLINEINCOMPAT, 
              "Island model has only a single period: do not use -j%d, -j%d, or -j%d with -a%d", 
              NOMIGBETWEENNONSISTERS, PARAMETERSBYPERIOD, 
              POPULATIONASSIGNMENTINFINITE);
    }
    if (outputoptions[THISTDIVIDEBYPRIOR] == 1
        || outputoptions[PRINTJOINTTEST] == 1)
    {
      IM_err (IMERR_COMMANDLINEINCOMPAT, 
              "Island model has no split time: do not use -p%d or -p%d with -a%d", 
              THISTDIVIDEBYPRIOR, 
              PRINTJOINTTEST,
              POPULATIONASSIGNMENTINFINITE);
    }
    if (outputoptions[POPMIGPARAMHIST] && (modeloptions[SINGLEMIGRATIONBOTHDIRECTIONS] || modeloptions[ONEMIGRATIONPARAMETER]  ))
    {
       IM_err (IMERR_COMMANDLINEINCOMPAT, 
              "output option for 2NM, -p%d, cannot be used with model options -j%d or -j%d",POPMIGPARAMHIST,SINGLEMIGRATIONBOTHDIRECTIONS,ONEMIGRATIONPARAMETER);
    }
      
    if (Tp == 1)
    {
      IM_err (IMERR_COMMANDLINEINCOMPAT, 
              "Island model has no split time: do not use -t with -a%d", 
              POPULATIONASSIGNMENTINFINITE);
    }


  }
  if (calcoptions[USEPRIORFILE] && !Gp) {
    strcpy(priorfilename,defaultpriorfilename);
	#ifdef MPI_ENABLED
		MPI_Bcast(&priorfilename, 500, MPI_CHAR, 0, MPI_COMM_WORLD);
	#endif //AS: TODO error catching
 }

  // VS
  // By definition when the command line is read the number of groups of loci is set to 1
  nbgroupsloci_theta = 1;
  nbgroupsloci_mig = 1;
	#ifdef MPI_ENABLED
		MPI_Bcast(&nbgroupsloci_theta, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&nbgroupsloci_mig, 1, MPI_INT, 0, MPI_COMM_WORLD);
	#endif //AS: TODO error catching
  //AS: now broadcast everything that was set in scan_commandline() function above
  int swapbetasonly = 1;
  //AS: check if these are still the same size as in IMa2 first before broadcasting!!!! 9/25/2014
  //AS: checked and look like this should be consistent.
  //AS: This version does not have the marginal likelihood calculations.
  #ifdef MPI_ENABLED
	MPI_Bcast(outputoptions, 10, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(modeloptions, 10, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(calcoptions, 5, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&swapbetasonly, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(assignmentlocioptions, 5, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(runoptions, 6, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(assignmentoptions, 13, MPI_INT, 0, MPI_COMM_WORLD);
  #endif //AS: TODO error catching


  return;
}                               // scan_commandline 


/*  this prints basic info to a string, fpstr, that later gets printed to the output file */
void
print_outputfile_info_string (void)
{
  fpstri = malloc (sizeof (int));
  *fpstri = 0;
  SP "IMa version 2.0 - Isolation with Migration Analysis  -  Jody Hey, Rasmus Nielsen, Sang Chul Choi 2010 \n");
  SP "Release date: %s\n\n",RELEASE_DATE);
  SP "\nINPUT AND STARTING INFORMATION \n");
  SP "================================\n");
  SP "\nCommand line string : %s \n", command_line);
  SP "  Input filename : %s \n", infilename);
  SP "  Output filename: %s \n", outfilename);
  SP "  Random number seed : %li \n", seed_for_ran1);
  SP "  Heating terms on command line : %s \n",heatingterm_str);
  SP "  Calculation options on command line : %s \n",calcoptions_str);
  SP "  Model options on command line : %s \n",modeloptions_str);
  SP "  Output options on command line : %s \n",outputoptions_str);
  if (calcoptions[DONTCALCLIKELIHOODMUTATION])
    SP "**NO DATA ** - Data likelihoods set to constant  posterior should equal prior \n");
  if (!runoptions[LOADRUN])
  {
    SP "- Run Duration - \n");
    switch (burndurationmode)
    {
    case TIMESTEPS:
      SP "     Burn period, # steps: %li \n", burnduration);
      break;
    case TIMEINF:
      SP "     Burn period, # seconds: %li (total burn duration depends on IMburn file status)\n", burnduration);
      break;
    };
    if (runoptions[PRINTBURNTREND])
    {
      SP "          -User option for printing trendline during, or at end of burnin period, invoked\n");
      SP "          -initial burn duration prior to beginning recording burntrend : %d steps\n", burntrendstartdelay);
    }

    switch (cdurationmode)
    {
    case TIMESTEPS:
      SP "     Record period, #saves: %d  #steps each: %li   total #steps: %li \n", genealogiestosave, chainduration / genealogiestosave, chainduration);
      break;
    case TIMEINF:
      SP "      Record period, # seconds per interval: %ld \n",
        chainduration);
      SP "      Run period indefinite in length - determined by user using 'IMrun' file\n");
      break;
    }
    
    

    SP "- Metropolis Coupling -\n");
    if (numchains > 1)
    {
      SP "     Metropolis Coupling implemented using %d chains \n", numchains);
      switch (heatmode)
      {
      case HLINEAR:
        SP "     Linear Increment Model   term: %.3f\n", hval1);
        break;
      case HTWOSTEP:
        SP "     Twostep Increment Model   term1: %.3f  term2: %.3f\n",
          hval1, hval2);
        break;
      case HGEOMETRIC:
        SP "     Geometric Increment Model   term1: %.3f  term2: %.3f\n",
          hval1, hval2);
        break;
      }
    }
    else
      SP "     None \n");
  }
  if (calcoptions[USEPRIORFILE])
    SP "Prior distribution terms loaded from file : %s \n",priorfilename);

  if (runoptions[LOADMCSTATE])
  {
    SP "Initial Markov chain state space loaded from file: %s\n",
      mcfreadfilename);
  }
  if (modeloptions[EXPOMIGRATIONPRIOR]==1)
  {
    SP"Exponential priors used for migration rate parameters\n");
  }
  if (runoptions[SAVEMCSTATEFILE])
  {
    SP "State of Markov chain saved to file : %s\n", mcfwritefilename);
  }
  if (!runoptions[DONTSAVEGENEALOGIES])
  {
    SP "All genealogy information used for surface estimation printed to file: %s\n", genealogyinfosavefilename);
  }
  if (outputoptions[DONTSAVEGENEALOGIES])
  {
    SP "All genealogy information used for surface estimation printed to file: %s\n", genealogyinfosavefilename);
  }
  if (outputoptions[MIGRATEHIST])
  {
    /* CR:110114.2  message text changed */
    SP "Distributions for migration event counts saved in file: %s%s\n",
        outfilename,".mpt");
  }

  if (outputoptions[PRINTTMRCA])
    SP "TMRCA  histograms printed \n");
}                               //  print_outputfile_info_string


// sets up filenames and calls the main initialization functions 
void start (int argc, char *argv[], int currentid)
{
  int i;

  scan_commandline (argc, argv, currentid);

  if (runoptions[SAVEMCSTATEFILE])
  {
    //strcpy (mcfwritefilename, outfilename);
    //strcat (mcfwritefilename, ".mcf");
    sprintf(mcfwritefilename, "%s.mcf.%d", outfilename, currentid);
  }
  if (!runoptions[DONTSAVEGENEALOGIES])
  {
    //strcpy (genealogyinfosavefilename, outfilename);
    //strcat (genealogyinfosavefilename, ".ti");
    sprintf(genealogyinfosavefilename, "%s.ti.%d", outfilename, currentid);
  }
  if (cdurationmode == TIMEINF)
  {
    //strcpy (oldoutfilename, outfilename);
    sprintf(oldoutfilename, "%s.old.%d", outfilename, currentid);
    //strcat (oldoutfilename, ".old");
  }
  for (i = 0; i < MAXLOCI; i++)
    hilocuslike[i] = -1e20;
  print_outputfile_info_string ();
  //AS: Setting seeds based on currentid
  setseeds (seed_for_ran1 * currentid);
  setlogfact ();
  setup (infilename, fpstri, fpstr,priorfilename);
  //AS: to change setheat function
  setheat (hval1, hval2, heatmode, currentid);
  //AS: to change checkautoc function
  checkautoc (1, 0, 0, currentid);
  gsampinflength = calc_gsampinf_length ();
  
  if (outputoptions[MIGRATEHIST])
  {
    //strcpy (migplotfilename, outfilename);
    //strcat (migplotfilename, ".mpt");
    sprintf(migplotfilename, "%s.mpt.%d", outfilename, currentid);
    nummigdirs = 2*(npops-1)*(npops-1);
    migcount = malloc ((nloci + (nloci > 1)) * sizeof (int *));
    for (i = 0; i < nloci + (nloci > 1); i++)
    {
      migcount[i] = malloc (nummigdirs * sizeof (int));
    }
    migfrom = malloc (nummigdirs * sizeof (int));
    migto = malloc (nummigdirs * sizeof (int));
    for (i = 0; i < nummigdirs; i++)
    {
      migfrom[i]= atoi(migration_counts_times[0][2*i].str);
      migto[i] = atoi(strchr(migration_counts_times[0][2*i].str,'>')+1);
    } 
  }
  if (runoptions[LOADMCSTATE])
  {
	sprintf(mcfreadfilename, "%s.mcf.%d", mcfreadfilename, currentid);
    readmcf (mcfreadfilename);
  }
}                               /* start */



#define TUPDATEINC  0           //1            // do a single t parameter every TUPDATEINC steps
#define UUPDATEINC  4//9           // u parameters seem to mix well so do every 5 steps

// VS
// QUPDATE
// one of the most important functions of the algorithm
// this will update the genealogy of a given locus
// AS: adding currentid, swapper and swappee information
void qupdate (int currentid, int swapA, int swapB)
{
  // VS debbug files to print assignment vectors
  // AS: are these files to be printed on all processors? TODODODO: check this!!! 9/25/2014
  FILE *checkassignMig;
  FILE *checkassignTheta;

  int i;
  int j, k, li, ci, ui;
  int changed;
  int qswapped = 0;
  //double tv = 0.0; //AS: adding this -  edit - took it out after discusion with Jody Tue Feb 23 14:02:56 EST 2016
  #ifdef MPI_ENABLED
  MPI_Status status;
  #endif
  int topolchange, tmrcachange;
  int periodpick, tupdatemethodpick;
  // VS auxiliary variable to update assignment of loci into groups
  int upta=0;

  static int tui = 0;
  static int count_t_updatetypes = -1;
  static int NWi = -1, RY1i = -1;
  static int uui = 0;
  int z; //AS: indexes coldchain
  if (count_t_updatetypes == -1)
  {
    i = 0;
#ifdef DO_RY1UPDATE  //jh 8/16/10  swapped with NW  to match an older program for debugging
      i++;
      RY1i = i;
#endif
#ifdef DO_NWUPDATE
    if (modeloptions[NOMIGRATION] == 0)
    {
      i++;
      NWi = i;
    }
#endif
    count_t_updatetypes = i;
  }

  z = whichiscoldchain();
 /* #ifdef MPI_ENABLED
  if ( z >= 0) {
	tv = C[z]->tvals[lastperiodnumber];
  }

   MPI_Bcast(&tv, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  #endif*/ // AS: took this out Tue Feb 23 14:03:25 EST 2016
//  printf ("After broadcast in qupdate...\n");
  /* update genealogies */
  for (ci = 0; ci < numchains; ci++)
  {
    for (li = 0;li<nloci;li++)
    {
      //AS:Genealogy update saving should happen on the cold chain, here indexed by z
      if (ci == z)
      {
        L[li].g_rec->upinf[IM_UPDATE_GENEALOGY_ANY].tries++;
        L[li].g_rec->upinf[IM_UPDATE_GENEALOGY_TOPOLOGY].tries++;
        L[li].g_rec->upinf[IM_UPDATE_GENEALOGY_TMRCA].tries++;
      }
	  // 1. UPDATEGENEALOGY()
	  // update genealogy is a very important function where the genealogy of a given locus is updated
      if (updategenealogy (ci, li, &topolchange, &tmrcachange))
      {
          if (ci == z)
          {
            L[li].g_rec->upinf[IM_UPDATE_GENEALOGY_ANY].accp++;
            L[li].g_rec->upinf[IM_UPDATE_GENEALOGY_TOPOLOGY].accp += (topolchange > 0);
            L[li].g_rec->upinf[IM_UPDATE_GENEALOGY_TMRCA].accp += (tmrcachange > 0);
          }
      }
    }

	// UPDATE ASSIGNMENT OF LOCI (to detect loci that belong to different groups)
	// Update assignments of vectors with groups of loci for mig
	if(nbgroupsloci_mig>1) {
		upta = updategrouplociassign_mig (ci);
		
		// if chain z, then increase the number of tries AS
		if(ci==z) {
			assignloci[0].upinf[0].tries++;
			// if the assigbment vector is updated and chain is 0, increase the number of updates
			if(upta==1)  assignloci[0].upinf[0].accp++; //AS: this has to be kept track of on each processor and then "reduced"?? TODO 9/25/2014
		}
		//printf("updated assignment of loci\n");
	}
	// Update assignments of vectors with groups of loci for theta
	if(nbgroupsloci_theta>1) {
		upta = updategrouplociassign_theta (ci);

		// if chain z(cold), then increase the number of tries AS
		if(ci==z) {
			assignloci[1].upinf[0].tries++;
			// if the assigbment vector is updated and chain is 0, increase the number of updates
			if(upta==1)  assignloci[1].upinf[0].accp++;
		}
		//printf("updated assignment of theta\n");

	} //AS: TODO assignloci.upinf has to be "reduced" at regular intervals? 9/25/2014

//	if(nbgroupsloci_mig==nbgroupsloci_theta) {
//		upta = updategrouplociassign_thetamig(ci);
//		
//		// if chain 0, then increase the number of tries
//		if(ci==0) {
//			assignloci[0].upinf[0].tries++;
//			assignloci[1].upinf[0].tries++;
//			// if the assigbment vector is updated and chain is 0, increase the number of updates
//			if(upta==1)  {
//				assignloci[0].upinf[0].accp++;
//				assignloci[1].upinf[0].accp++;
//			}
//		}
//	}
  }

 /* update population assignment */
  if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
  {
    for (ci = 0; ci < numchains; ci++)
    {
      for (i = 0; i < snupdatei; i++)
      {
        if (assignmentoptions[POPULATIONASSIGNMENTRELABEL] == 1)
        {
          updateassignmentrelabel (ci);
        }
        if (assignmentoptions[POPULATIONASSIGNMENTBF] == 1)
        {
          updateassignmentbf (ci);
        }
      }
    }
    if (assignmentoptions[POPULATIONASSIGNMENTASSIGNED] == 1)
    {
      /* for DNA Barcoding */
      imaBarTick ();
    }
  } 

  if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 0 && npops > 1  && tui == TUPDATEINC)
  {
    for (ci = 0; ci < numchains; ci++)
	{
      periodpick = randposint (numsplittimes);
      tupdatemethodpick = randposint (count_t_updatetypes)+1;    
      if (modeloptions[NOMIGRATION])
        tupdatemethodpick = 0;  // must do ry1 , NW does not work with no migration

      if (tupdatemethodpick == RY1i)
      {
#ifdef DO_RY1UPDATE
	//printf("Going into RY update...\n");
        changed = changet_RY1 (ci, periodpick/*, tv*/); //AS: changed this upon discussion with Jody Tue Feb 23 13:51:15 EST 2016
        if (ci == z)
        {
          T[periodpick].upinf[IM_UPDATE_TIME_RY1].tries++;
          if (changed)
            T[periodpick].upinf[IM_UPDATE_TIME_RY1].accp++;
        }
	//printf("Here done RY update...\n");
#endif //DO_RY1UPDATE
      }
		if (tupdatemethodpick == NWi)
	  {
#ifdef DO_NWUPDATE
        /* SANGCHUL: NW update does not see eye to eye with
         * assignment update at the moment. */
	    if (assignmentoptions[POPULATIONASSIGNMENTBF] == 1)
        {
          changed = changet_RY1 (ci, periodpick/*, tv*/);
        }
        else
        {
          changed = changet_NW (ci, periodpick);
        }
        if (ci == z) //AS: z is coldchain
        {
          T[periodpick].upinf[IM_UPDATE_TIME_NW].tries++;
          if (changed)
           T[periodpick].upinf[IM_UPDATE_TIME_NW].accp++;
        }
#endif /* DO_NWUPDATE */
      } 
	}
    tui = 0;
  }
  else
  {
    tui++;
  }

  

  if (uui == UUPDATEINC)
  {
    if (nurates > 1)
      for (ci = 0; ci < numchains; ci++)
      {
        for (j = 0; j < (nurates - (nurates == 2)); j++)
        {
          ui = changeu (ci, j, &k);
          if (ci == z)

          {
            L[ul[j].l].u_rec[ul[j].u].upinf->tries++;
            L[ul[k].l].u_rec[ul[k].u].upinf->tries++;
            if (ui == 1)

            {
              L[ul[j].l].u_rec[ul[j].u].upinf->accp++;
              L[ul[k].l].u_rec[ul[k].u].upinf->accp++;
            }
          }
        }
      }
    else
    {
      if (nloci == 1 && L[0].model == HKY)      /* if there is just one HKY locus kappa needs updating on its own */
        for (ci = 0; ci < numchains; ci++)
          changekappa (ci);
    }
    uui = 0;
  }
  else
  {
    uui++;
  }
//As: debug only
//printf("numchains = %d, numprocesses = %d\n",numchains,numprocesses);
 //AS: If there is only one process
  if (numchains > 1 && numprocesses == 1)
  {
	//printf("huh?!?!?!\n");
    /* CR 110929.4 get rid of extraneous args in function call
     * although return val from swapchains() call is not used locally,it may
     * be useful during debugging
     */
    swapbetasonly = 1;
    qswapped = swapchains (swaptries, swapbetasonly, currentid); //AS: TODO edit the swapchains function 9/25/2014
  }

  // VS
  // Print the assignment vectors to a test output file
  // VS file to check the assignment vectors of loci
  /*if (burndone) {
	  if((step % 100) == 0) {
		  checkassignMig = fopen ("checkAssignLociMig.out", "a");
		  
		  //fprintf(checkassignMig, "step=%i ", step);

		  for(li=0; li<nloci;li++) {
			fprintf(checkassignMig, "%i ", grouploci_mig[0][li]);
		  }	
		  fprintf(checkassignMig, "\n");
		  fclose(checkassignMig);
		  checkassignTheta = fopen ("checkAssignLociTheta.out", "a");
		  for(li=0; li<nloci;li++) {
			fprintf(checkassignTheta, "%i ", grouploci_theta[0][li]);
		  }	
		  fprintf(checkassignTheta, "\n");
		  fclose(checkassignTheta);
	  }
  }*/
  if (numchains > 1 && numprocesses > 1) {
	qswapped = swapchains_bwprocesses(currentid, step, swaptries, swapbetasonly, chainduration, burnduration, swapA, swapB);
  }
  //AS: now check which is coldchain again, since it could have been swapped, before printing intervaloutput
  //AS: Tue Feb 23 15:40:04 EST 2016
  //z = whichiscoldchain();
  //if (currentid == 0) {
  	intervaloutput (stdout, currentid);
  //}
  #ifdef MPI_ENABLED
	MPI_Barrier(MPI_COMM_WORLD);
  #endif
  if (step >= CHECKAUTOCWAIT)
    checkautoc (0, burndone, burnsteps, currentid); //AS: TODO change checkautoc function as on 9/25/2014

  return;
}                               /* qupdate */



void reset_after_burn (int currentid)
{
  int ci;
  int li, ui, i, j;
  int z; //AS: coldchain 
  time (&chainstarttime);
  burnsteps = step - 1;
  if (cdurationmode == TIMEINF)
  {
    time (&lasttime);
    time (&timer);
  }
  // reset acceptance rate accumulators

  if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 1
      || npops == 1)
  {
    /* no split times */
  }
  else
  {
    for (i = 0; i < lastperiodnumber; i++)
      for (j = 0; j < T[i].num_uptypes; j++)
        T[i].upinf[j].accp = T[i].upinf[j].tries = 0;
  }
  if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
  {
    for (ci = 0; ci < numchains; ci++)
    {
      for (j = 0; j < Cupinf[ci].num_uptypes; j++)
      {
        Cupinf[ci].upinf[j].accp = 0;
        Cupinf[ci].upinf[j].tries = 0;
      }
    }
  }

  for (li = 0; li < nloci; li++)
  {
    if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
    {
      for (j = 0; j < L[li].a_rec->num_uptypes; j++)
      {
        L[li].a_rec->upinf[j].accp = 0;
        L[li].a_rec->upinf[j].tries = 0;
      }
      if (assignmentoptions[POPULATIONASSIGNMENTASSIGNED] == 1)
      {
        imaBarReset ();
      }
    }
    for (ui = 0; ui < L[li].nlinked; ui++)
    {
      for (j = 0; j < L[li].u_rec[ui].num_uptypes; j++)
        L[li].u_rec[ui].upinf[j].accp = L[li].u_rec[ui].upinf[j].tries = 0;
      if (L[li].model == HKY)
        for (j = 0; j < L[li].kappa_rec->num_uptypes; j++)
          L[li].kappa_rec->upinf[j].accp = L[li].kappa_rec->upinf[j].tries =
            0;
      if (L[li].model == STEPWISE || L[li].model == JOINT_IS_SW)
      {
        for (j = 0; j < L[li].A_rec[ui].num_uptypes; j++)
          L[li].A_rec[ui].upinf[j].accp = L[li].A_rec[ui].upinf[j].tries = 0;
      }
      for (j = 0; j < L[li].g_rec->num_uptypes; j++)
        L[li].g_rec->upinf[j].accp = L[li].g_rec->upinf[j].tries = 0;
    }
  }
/*updated the way checkautoc() initializes on 3/25/08 */
  checkautoc (1, burndone, burnsteps, currentid);
  if (!runoptions[DONTSAVEGENEALOGIES])
  {
    if ((genealogyinfosavefile = fopen (genealogyinfosavefilename, "w")) == NULL)
    {
      IM_err (IMERR_CREATEFILEFAIL, "Error creating file for holding genealogy information");
    }
    fprintf (genealogyinfosavefile,
             "-------------------------------------------\n\n");
    fprintf (genealogyinfosavefile, "Header for genealogy file:  %s\n\n",
             genealogyinfosavefilename);
    fprintf (genealogyinfosavefile,
             "-------------------------------------------\n\n");
    fprintf (genealogyinfosavefile, "%s\n", fpstr);
    fprintf (genealogyinfosavefile,
             "-------------------------------------------\n\n");
    fprintf (genealogyinfosavefile, "End of header for genealogy file:  %s\n\n",
             genealogyinfosavefilename);
    fprintf (genealogyinfosavefile,
             "-------------------------------------------\n\n");
    fprintf (genealogyinfosavefile, "VALUESSTART\n");
    f_close (genealogyinfosavefile);
  }
  if (runoptions[SAVEMCSTATEFILE])
  {
    writemcf (mcfwritefilename);
  }
  z = whichiscoldchain(); //AS: TODO add whichiscoldchain function 9/25/2014
  if (z >= 0) {
	  gloglikelihood = C[z]->allpcalc.pdg; //AS: changing this - previously assumed that chain 0 was always coldchain
  }
}                               /* reset_after_burn() */

void output_burntrendfile (int currentid)
{
  FILE *burntrendfile;
  char burntrendfilename[FNSIZE];
  assert (runoptions[PRINTBURNTREND]);
  //strcpy (burntrendfilename, outfilename);
  //strcat (burntrendfilename, ".burntrend.out");
  sprintf(burntrendfilename, "%s.burntrend.out.%d", outfilename, currentid);
  if ((burntrendfile = fopen (burntrendfilename, "w")) == NULL)
  {
    IM_err (IMERR_CREATEFILEFAIL,
            "Error opening burntrend output file for writing");
  }
  printf ("\n\n========================\n");
  printf ("Printing Burn Trend File\n");
  printf ("========================\n\n");
  fprintf (burntrendfile,
           "Plots of Runtime Information and Parameter Trends during Burnin \n");
  fprintf (burntrendfile, "========================================\n");
  intervaloutput (burntrendfile, currentid);
  fprintf (burntrendfile, "========================================\n\n");
  fprintf (burntrendfile, "Current Step #: %d \n\n", step);
 //AS: ASCII plots called only on head node
  if (trendspot > 1 && currentid == 0)
    callasciitrend (burntrendfile);
  else if (trendspot <= 1 && currentid == 0) {
    fprintf(burntrendfile, "burn period too short to plot trends,  trend recording begins at step %d \n",burntrendstartdelay); 
  }
  if (currentid > 0) {
	fprintf(burntrendfile, "Check burn trend file on head processor (0) for trends!\n\n");
  }
  fclose (burntrendfile);
}                               // output_burntrendfile

#define CHECKINTERVALSTEPS 1000

/* run() determines the status of a run in terms of whether it is in the burnin phase or not
  and of whether the run should keep going. 
  
  if the burnin period has just ended,  some work is done

  run()  also  checks to see if it is time to print an output file */

int run (int currentid)
{
  static int checkinterval = 0;
  static int printburntrendstep = 0, burnrecordi = 1;
  int tempburndone;
  char ch;
  int gp; // VS gp - auxiliary index for groups of loci
  char tempgenealogyinfosavefilename[1000]; // VS filename string used to save the TI files 
  int x, y, z; //AS: indices for MPI reduce loops
  #ifdef MPI_ENABLED
	MPI_Status status;
  #endif
 
  if (burndone)
  {
    switch (cdurationmode)
    {
    case TIMESTEPS:
      return (step < (chainduration + burnsteps));
	//printf("chainduration and burnsteps here are: %d %d\n", chainduration, burnsteps);
      break;
    case TIMEINF:
      if (checkinterval < CHECKINTERVALSTEPS)
      {
        checkinterval++;
        return (1);
      }
      else
      {
        checkinterval = 0;
        if (maxedoutgenealogysave)
          return (0);
        time (&timer);
        if ((timer - lasttime) > chainduration)
        {
          if ((checkdonefile = fopen ("IMrun", "r")) == NULL)
          {
            return (0);
          }
          else
          {
            ch = (char) getc (checkdonefile);
            f_close (checkdonefile);
            if ((char) toupper (ch) != 'Y')
            {
              return (0);
            }
            else
            {
	//AS: adding functions to MPI Reduce the swap counts
	//AS: 10/22/2014
	/*#ifdef MPI_ENABLED
		//MPI_Barrier(MPI_COMM_WORLD);
		if (numprocesses > 1) {
			for (x = 0; x < numprocesses; x++) {
				for (y = 0; y < numprocesses; y++) {
					MPI_Reduce(&swaps_bwprocesses[x][y], &swaps_rec_bwprocesses[x][y], 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
				}
			}
			if (currentid == 0) {
				for (x = 0; x < numchains; x++) {
					for (y = 0; y < numchains; y++) {
						swapcount_bwprocesses[x][y] = swapcount[x][y];
					}
				}
			}
			for (x = 1; x < numprocesses; x++) {
				if (currentid == x) {
					for (y = 0; y < numchains; y++) {
						for (z = 0; z < numchains; z++) {
							MPI_Send(&swapcount[y][z], 1, MPI_INT, 0, 1234, MPI_COMM_WORLD);
						}
					}
				}
				if (currentid == 0) {
					for (y = 0; y < numchains; y++) {
						for (z = 0; z < numchains; z++) {
							MPI_Recv(&swapcount_bwprocesses[x * numchains + y][x * numchains + z], 
								1, MPI_INT, x, 1234, MPI_COMM_WORLD, &status);
						}
					}
				}
			}
		}
	#endif*/
		if (currentid == 0) {
		//AS: TODO Swap trials and acceptance information has to be reduced 9/25/2014 - see ima_main_mpi.cpp line 2237
              if (!runoptions[DONTSAVEGENEALOGIES] && genealogiessaved > 0)
              {
                // VS
                //savegenealogyfile (genealogyinfosavefilename, genealogyinfosavefile, &lastgenealogysaved, gsampinflength);
				// VS
				// SAVE GENEALOGIES - save the summaries of the trees into files. Each group of loci is one file.
				
				for(gp=0; gp<nbgroupsloci_theta+nbgroupsloci_mig; gp++) {
					if(gp<nbgroupsloci_theta) { // if the group is for theta
						sprintf(tempgenealogyinfosavefilename, "%s.groupTheta%i", genealogyinfosavefilename, gp); // create name of file for group Theta
					}
					else {
						sprintf(tempgenealogyinfosavefilename, "%s.groupMig%i", genealogyinfosavefilename, gp); // create name of file for group Mig
					}
					// VS modified savegenealogyfile_vs to get as input the group index
					savegenealogyfile_vs(gp, tempgenealogyinfosavefilename, genealogyinfosavefile, lastgenealogysaved, gsampinflength);
				}

				// VS SAVE ASSIGNMENT VECTORS if vectors are being updated
				// get file name for assigment vector theta
				sprintf(tempgenealogyinfosavefilename, "%s.assignlociTheta", genealogyinfosavefilename);
				saveassignlocifile(tempgenealogyinfosavefilename, lastgenealogysaved, assignlocisample_theta);	


				// VS get file name for assigment vector mig
				sprintf(tempgenealogyinfosavefilename, "%s.assignlociMig", genealogyinfosavefilename);
				saveassignlocifile(tempgenealogyinfosavefilename, lastgenealogysaved, assignlocisample_mig);	
				
				// VS Change the number of genealogies and assignment vectors saved
				lastgenealogysaved = genealogiessaved - 1;

              }
		if (currentid == 0) {
              printoutput (currentid);
		}
	     }
		#ifdef MPI_ENABLED
			MPI_Barrier(MPI_COMM_WORLD);
		#endif
		//AS: All other processors will wait until the head node completes I/O operations above
              time (&lasttime); // start the clock again 
            }
          }
        }
        return (1);
      }
      break;
    default:
      return (0);
      break;
    }
  }
  else
  {
    tempburndone = 0;
    switch (burndurationmode)
    {
    case TIMESTEPS:
      tempburndone = (step > burnduration);
      break;
    case TIMEINF:
      if (checkinterval < CHECKINTERVALSTEPS)
      {
        checkinterval++;
      }
      else
      {
        checkinterval = 0;
        time (&timer);
        tempburndone = (timer - lasttime) > burnduration;
      }
      break;
    default:
      return (0);
      break;
    }
    // plot burn trend if called for 
    if (runoptions[PRINTBURNTREND] && (tempburndone || step >= burntrendstartdelay))
    {
      if (burnrecordi == recordint)
      {
        trendrecord (-1, currentid);       // record values of parameters that are in mcmc
	//AS: adding currentid to trendrecord  - TODO as on 9/25/2014
        burnrecordi = 1;
      }
      else
      {
        burnrecordi++;
      }
      if (tempburndone)
      {
	//AS: TODO - MPI reduce all the swapcount information prior to output - as on 9/25/2014
	//AS: DONE as on 10/22/2014
	
	/*#ifdef MPI_ENABLED
		//MPI_Barrier(MPI_COMM_WORLD);
		if (numprocesses > 1) {
			for (x = 0; x < numprocesses; x++) {
				for (y = 0; y < numprocesses; y++) {
					MPI_Reduce(&swaps_bwprocesses[x][y], &swaps_rec_bwprocesses[x][y], 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
				}
			}
			if (currentid == 0) {
				for (x = 0; x < numchains; x++) {
					for (y = 0; y < numchains; y++) {
						swapcount_bwprocesses[x][y] = swapcount[x][y];
					}
				}
			}
			for (x = 1; x < numprocesses; x++) {
				if (currentid == x) {
					for (y = 0; y < numchains; y++) {
						for (z = 0; z < numchains; z++) {
							MPI_Send(&swapcount[y][z], 1, MPI_INT, 0, 1234, MPI_COMM_WORLD);
						}
					}
				}
				if (currentid == 0) {
					for (y = 0; y < numchains; y++) {
						for (z = 0; z < numchains; z++) {
							MPI_Recv(&swapcount_bwprocesses[x * numchains + y][x * numchains + z], 
								1, MPI_INT, x, 1234, MPI_COMM_WORLD, &status);
						}
					}
				}
			}
		}
	#endif*/
        output_burntrendfile (currentid); //AS: TODO add currentid to output_burntrendfile

        printburntrendstep = 1;
        /* now check to see if IMburn file is present */
        if (burndurationmode == TIMEINF)
        {
          if ((checkdonefile = fopen ("IMburn", "r")) != NULL)
          {
            ch = (char) getc (checkdonefile);
            f_close (checkdonefile);
            if ((char) toupper (ch) == 'Y')
            {
              tempburndone = 0;
              time (&lasttime); // start the clock again 
            }
          }
        }
      }
      else
      {
        printburntrendstep++;
      }
    }
    if (tempburndone)
    {
      burndone = 1;
      reset_after_burn (currentid); //AS: TODO - add currentid to reset_after_burn as on 9/25/2014
    }
	//printf("Returning 1\n\n");
    return (1);
  }
}                               /* run */

#define  TRENDLASTPT (TRENDDIM - 1)
void inctrend (int m, int t, struct value_record *v, double newval)
{
  int j;
  for (j = m; j < TRENDLASTPT; j++)
  {
    /* assert (v->trend[j + 1] != 0); This asserts sometimes. */
    v->trend[j] = v->trend[j + 1];
  }
  v->trend[t] = newval;
  //assert(v != 0);
}

void trend_reset (struct value_record *v, int nv)
{
  int i, j;
  for (i = 0; i < nv; i++)
    for (j = 0; j < TRENDDIM; j++)
      v[i].trend[j] = 0;
}                               // trend_reset


/* Using trendrecord()
This function records trend lines
It works on instances of struct value_record 
the value_record must first be initialized (probably in initialize.c) 
Code for a particular value_record or array of value_records 
can be placed in trendrecord at two places:
  in the "if (burndone && reset == 0) " section
  and in the "if (recordinc == recordtrendinc)" section


explanation for how trendline data are accumulated:
 - movespot is the point at which values are deleted from the array, 
 - each new value is added to the end of the array (i.e. at trendspot)
 - all values from one position past movespot up to the end of the array are moved down one position
 - this erases the value at movespot and makes room for a new value at the end. 
 -each time the replacement point (movespot) reaches the end of the array the time 
	period between additions doubles 
values are added more slowly as the run proceeds.  
- this is because the time period doubles when movespot reaches the end 
- the values to the left of movespot have half the density (in time) of those to the right 
*/
//AS: adding currentid
void trendrecord (int loadarrayj, int currentid)
{
  int gp; // VS gp
  int z; //AS: coldchain
  double probrec; //AS: adding this for receiving
  static int /*trendspot = 0,*/ recordtrendinc = 1, recordinc = 0, movespot =
    TRENDLASTPT;
  static int init = 1, reset = 0;
  int j, k, li, ui;
  #ifdef MPI_ENABLED
  MPI_Status status;
  #endif

  if (burndone && reset == 0)   // reset all trend-related values after burnin
  {
    init = 1;
    reset = 1;
    if (lpgpd_v->do_trend)
    {
      trend_reset (lpgpd_v, 1);
    }
    for (li = 0; li < nloci; li++)
    {
      if (runoptions[LOADRUN] == 0)
      {
        if (L[li].g_rec->v->do_trend)
        {
          trend_reset (L[li].g_rec->v, 1);
        }
        for (ui = 0; ui < L[li].nlinked; ui++)
        {
          if (L[li].u_rec[ui].v->do_trend)
          {
            trend_reset (L[li].u_rec[ui].v, 1);
          }
        }
        if (L[li].model == HKY)
        {
          if (L[li].kappa_rec->v->do_trend)
          {
            trend_reset (L[li].kappa_rec->v, 1);
          }
        }
        if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
        {
          if (L[li].a_rec->v->do_trend == 1)
          {
            trend_reset (L[li].a_rec->v, 1);
          }
        }
      }
    }
    if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 1
        || npops == 1)
    {
      /* no split time */
    }
    else
    {
      for (k = 0; k < lastperiodnumber; k++)
        if (T[k].v->do_trend)
          trend_reset (T[k].v, 1);
    }

/* ADD ADDITONAL trend_reset() calls here */
  }
  if (init == 1)
  {
    trendspot = 0;
    recordtrendinc = 1;
    recordinc = 0;
    movespot = TRENDLASTPT;
    init = 0;
  }
  recordinc++;
  if (recordinc == recordtrendinc)
  {
    if (runoptions[LOADRUN])
    {
	//AS: I need to think about how to handle this 3D array for gsampinf - as on 9/25/2014
		// VS added the for loop through groups of loci
		for(gp=0; gp<nbgroupsloci_theta+nbgroupsloci_mig; gp++) 
		{	
			if (lpgpd_v->do_trend) {
				inctrend (movespot, trendspot, lpgpd_v,
						gsampinf[gp][loadarrayj][gsamp_pdgp] + // VS gp
						gsampinf[gp][loadarrayj][gsamp_probgp]); // VS gp
			}
		}

      if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 0
          && npops > 1)
      {
		  for (j = 0; j < lastperiodnumber; j++) {
			  
			  // VS added the for loop through groups of loci
			  for(gp=0; gp<nbgroupsloci_theta+nbgroupsloci_mig; gp++) 
			  {	
				if (T[j].v->do_trend) {
					inctrend (movespot, trendspot, T[j].v,
                      gsampinf[gp][loadarrayj][gsamp_tp + j]); // VS gp
				}
			  }
		  }
      }
    }
    else
    {
      if (lpgpd_v->do_trend)
      {
	z = whichiscoldchain();
	if (z >= 0 && currentid == 0) {
	        inctrend (movespot, trendspot, lpgpd_v,
                  C[z]->allpcalc.probg + C[z]->allpcalc.pdg);
	}
	#ifdef MPI_ENABLED
	if (z < 0 && currentid == 0) {
		probrec = 0.0;
		MPI_Recv(&probrec, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 1313, MPI_COMM_WORLD, &status);
		inctrend (movespot, trendspot, lpgpd_v, probrec);
	}
	if (z >= 0 && currentid != 0) {
		probrec = C[z]->allpcalc.probg + C[z]->allpcalc.pdg;
		MPI_Send(&probrec, 1, MPI_DOUBLE, 0, 1313, MPI_COMM_WORLD);
	}
	#endif
      }

      if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 1
          || npops == 1)
      {
        /* no split time */
      }
      else
      {
        for (j = 0; j < lastperiodnumber; j++)
          if (T[j].v->do_trend) {
		z = whichiscoldchain();
		if (z >= 0 && currentid == 0) {
	            inctrend (movespot, trendspot, T[j].v, C[z]->tvals[j]);
		}
		#ifdef MPI_ENABLED
		if (z < 0 && currentid == 0) {
			probrec = 0.0;
			MPI_Recv(&probrec, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 1717, MPI_COMM_WORLD, &status);
			inctrend(movespot, trendspot, T[j].v, probrec);
		}
		if (z >= 0 && currentid != 0) {
			probrec = C[z]->tvals[j];
			MPI_Send(&probrec, 1, MPI_DOUBLE, 0, 1717, MPI_COMM_WORLD);
		}
		#endif
	}

      }
    }
    for (li = 0; li < nloci; li++)
    {
      if (runoptions[LOADRUN] == 0)
      {
        for (ui = 0; ui < L[li].nlinked; ui++)
          if (L[li].u_rec[ui].v->do_trend) {
		z = whichiscoldchain();
		if (z >= 0 && currentid == 0) {	
            inctrend (movespot, trendspot, L[li].u_rec[ui].v,
                      C[z]->G[li].uvals[ui]);
		}
		#ifdef MPI_ENABLED
		if (z < 0 && currentid == 0) {
			probrec = 0.0;
			MPI_Recv(&probrec, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 1414, MPI_COMM_WORLD, &status);
			inctrend(movespot, trendspot, L[li].u_rec[ui].v, probrec);
		}
		if (z >= 0 && currentid != 0) {
			probrec = C[z]->G[li].uvals[ui];
			MPI_Send(&probrec, 1, MPI_DOUBLE, 0, 1414, MPI_COMM_WORLD);
		}
		#endif
	}
	
        if (L[li].model == HKY)
          if (L[li].model == HKY) {
		z = whichiscoldchain();
		if (z >= 0 && currentid == 0) {
	        	inctrend (movespot, trendspot, L[li].kappa_rec->v,
                      		C[z]->G[li].kappaval);
		}
		#ifdef MPI_ENABLED
		if (z < 0 && currentid == 0) {
			probrec = 0.0;
			MPI_Recv(&probrec, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 1515, MPI_COMM_WORLD, &status);
			inctrend(movespot, trendspot, L[li].kappa_rec->v, probrec);
		}
		if (z >= 0 && currentid != 0) {
			probrec = C[z]->G[li].kappaval;
			MPI_Send(&probrec, 1, MPI_DOUBLE, 0, 1515, MPI_COMM_WORLD);
		}
		#endif
	   }
        if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
        {
          inctrend (movespot, trendspot, L[li].a_rec->v, 
                    imaAsnValue (0, li));
        }
      }
    }
/* ADD ADDITONAL inctrend() calls here */
    if (movespot == TRENDLASTPT && trendspot == TRENDLASTPT)
    {
      movespot = 0;
      recordtrendinc *= 2;
    }
    else
    {
      movespot += (movespot < TRENDLASTPT);
    }
    trendspot += (trendspot < TRENDLASTPT);
    recordinc = 0;
  }
  trenddoublepoint = movespot;
}                               /* trendrecord */

/* calculates the bin number of an xy array of a value_record that a value falls in,  increments the count in that bin */
void recordval (struct value_record *v, double val)
{
  int k;
  double logval;
  if (v->do_logplot)
  {
    assert (!(val < 0.0));
    logval = log (val);
    k =
      (int) (GRIDSIZE * (logval + v->plotrange.max) /
             (2.0 * v->plotrange.max));
  }
  else
  {
    k = (int) (GRIDSIZE * val / v->plotrange.max);
  }
  if (k < 0)
  {
    v->beforemin++;
  }
  else if (k >= GRIDSIZE)
  {
    v->aftermax++;
  }
  else
  {
    assert (!(k < 0));          // FIXME: it's been crashing
    v->xy[k].y++;
  }
  return;
}

/* this is used to record the names of actual loci and gene copies that migrate, from to to
write these names to a file */
void record_migration_names(void)
 {
  int i, j, li;
  int from, to;
  struct edge *gtree;
  int z = whichiscoldchain();
  if (z >= 0) {
  from = migrationnamefrom;
  to = migrationnameto;
  for (li = 0; li < nloci; li++)
  {
    gtree = C[z]->G[li].gtree;
    for (i = 0; i < L[li].numlines; i++)
    {
      j = 0;
      while (gtree[i].mig[j].mt > 0)
      {
        if (from == nowedgepop (0, &gtree[i], gtree[i].mig[j].mt) && to == C[z]->G[li].gtree[i].mig[j].mp)
        {
          fprintf(migrationnamefile,"%s ",L[li].name);
          if (i< L[li].numgenes)
            fprintf(migrationnamefile,"%s ",L[li].gNames[i]);
          else
            fprintf(migrationnamefile,"internal ");
        }
        j++;
      }
    }
  }
  fprintf(migrationnamefile,"\n");
 } // record_migration_names
 else {
	return;
 }
}



/* INSTRUCTIUONS to record a numerical value from the markov chain:
----------------------------------------------------
this works on instances of struct value_record

the value_record must be initiatlized (e.g. in initialize.c,  see e.g. init_g_rec)
this includes a call to init_value_record()

insert line(s) code into record()  below,  to make a call to recordval() 

*/
//AS: adding z and currentid - z indexes the coldchain, currentid is the processor ID
void record_migrations (int z, int currentid)
{
  int i, j, k, li, from, to, foundparam;
  struct edge *gtree;
  for (j = 0; j < nloci + (nloci > 1); j++)
    for (i = 0; i < nummigdirs; i++)
    {
      migcount[j][i] = 0;
    }
  for (li = 0; li < nloci; li++)
  {
    gtree = C[z]->G[li].gtree;
    for (i = 0; i < L[li].numlines; i++)
    {
      j = 0;
      while (gtree[i].mig[j].mt > 0)
      {
        from = nowedgepop (0, &gtree[i], gtree[i].mig[j].mt);
        to = C[z]->G[li].gtree[i].mig[j].mp;
        k = 0;
        foundparam = 0;
        while (k < nummigdirs && foundparam == 0)
        {
          if (from == migfrom[k] && to == migto[k])
            foundparam = 1;
          else
            k++;
        }
        assert(k<nummigdirs);
        if (nloci == 1)
        {
          migcount[0][k]++;
          /* CR:110114.2  1 line removed */
        }
        else
        {
          migcount[li + 1][k]++;
          migcount[0][k]++;
          /* CR:110114.2  2 lines removed */
        }
        j++;
      }
    }
  }
  for (i = 0; i < nloci + (nloci > 1); i++)
    for (k = 0; k < nummigdirs; k++)
    {
      recordval (&migration_counts_times[i][2 * k + 1], (double) migcount[i][k]);
    }
}                               //record_migrations

//AS: adding currentid
void record (int currentid)
{
  int j, li, ui;
  double rtime; //AS: adding rtime for receiving 
  struct genealogy *G;
  double uvals, kval, tval; //AS: adding uvals, kval, tvals for receiving
  int z = whichiscoldchain();
  #ifdef MPI_ENABLED
  MPI_Status status;
  #endif

  if (outputoptions[PRINTTMRCA])
    for (li = 0; li < nloci; li++)
    {
      if (z >= 0 && currentid == 0) {
		recordval (L[li].g_rec->v, C[z]->G[li].roottime);
	}
	#ifdef MPI_ENABLED
	if (z >= 0 && currentid != 0) {
		rtime = C[z]->G[li].roottime;
		MPI_Send(&rtime, 1, MPI_DOUBLE, 0, 2323, MPI_COMM_WORLD);
	}
	if (z < 0 && currentid == 0) {
		rtime = 0.0;
		MPI_Recv(&rtime, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 2323, MPI_COMM_WORLD, &status);
		recordval(L[li].g_rec->v, rtime);
	}
	#endif	
    }
  for (li = 0; li < nloci; li++)
  {
    if (z >= 0) {
	    G = &(C[z]->G[li]);
    }
    for (ui = 0; ui < L[li].nlinked; ui++)
    {
	if ( z >= 0 && currentid == 0) {
	      recordval (L[li].u_rec[ui].v, G->uvals[ui]);
	}
	#ifdef MPI_ENABLED
	if (z >= 0 && currentid != 0) {
		uvals = G->uvals[ui];
		MPI_Send(&uvals, 1, MPI_DOUBLE, 0, 2424, MPI_COMM_WORLD);
	}
	if (z < 0 && currentid == 0) {
		uvals = 0.0;
		MPI_Recv(&uvals, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 2424, MPI_COMM_WORLD, &status);
		recordval(L[li].u_rec[ui].v, uvals);
	}
	#endif	
    }
    if (L[li].model == HKY)
    {
	if ( z >= 0 && currentid == 0) {
	      recordval (L[li].kappa_rec->v, G->kappaval);
	}
	#ifdef MPI_ENABLED
	if (z >= 0 && currentid != 0) {
		kval = G->kappaval;
		MPI_Send(&kval, 1, MPI_DOUBLE, 0, 2345, MPI_COMM_WORLD);
	}
	if (z < 0 && currentid == 0) {
		kval = 0.0;
		MPI_Recv(&kval, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 2345, MPI_COMM_WORLD, &status);
		recordval (L[li].kappa_rec->v, kval);
	}
	#endif
    }
  }

  if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 1
      || npops == 1)
  {
    /* we could have lastperiodnumber be 0 */
  }
  else
  {
    for (j = 0; j < lastperiodnumber; j++)
    {
	if (z >= 0) {
	      assert (C[z]->tvals[j] > T[j].pr.min && C[z]->tvals[j] < T[j].pr.max);
	}
	if (z >= 0 && currentid == 0) {
      	      recordval (T[j].v, C[z]->tvals[j]);
	}
	#ifdef MPI_ENABLED
	if (z >= 0 && currentid != 0) {
		tval = C[z]->tvals[j];
		MPI_Send(&tval, 1, MPI_DOUBLE, 0, 2525, MPI_COMM_WORLD);
	}
	if ( z < 0 && currentid == 0) {
		tval = 0.0;
		MPI_Recv(&tval, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 2525, MPI_COMM_WORLD, &status);
		recordval (T[j].v, tval);
	}
	#endif
    }
  }
  if (outputoptions[MIGRATEHIST])
  {
    record_migrations (z, currentid); //AS: TODO add z and currentid to record_migrations 9/25/2014
  }

  if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 1)
  {
    /* no split time */
  }
  else
  {
    if (npops >= 3 && npops <= 5 && outputoptions[PRINTJOINTTEST])
      setup_multi_t_arrays (z); //AS: TODO add z to setup_multi_t_arrays 9/25/2014
  }
  if (!outputoptions[DONTPRINTASCIITREND])
  {
    trendrecord (-1, currentid); //AS: TODO add currentid to trendrecord 9/25/2014
  }
  return;
}                               /* record */

//AS: Adding this function - essentially returns the chain ID of cold chain, -1 if cold chain doesn't exist on this current processor
int whichiscoldchain(void)
{
	int which = -1;
	int i;
	for (i = 0; i < numchains; i++) {
		if (beta[i] == 1.0) {
			which = i;
		}
	}
	return which;
}



// VS
// modified this function to allocate memory to gsampinf
void savegenealogyinfo (int currentid)        // use floats to save space
{
  int i;
  int z; //AS: adding this - coldchain index
  int v; //AS: for loop while sending/receiving
  float **gsampinflocal; //AS: adding this, as a proxy for a 2D array of gsampinf, which I will send/receive in a loop
  int tempcurrentid;
  //AS: I am yet to figure out if this is the most efficient way of doing this, but it's the most intuitive way
  //AS: come back to this - 9/25/2014
  int gp; // VS gp
  #ifdef MPI_ENABLED
  MPI_Status status;
  #endif
  if (genealogiessaved == 0)
  {
    if (cdurationmode == TIMESTEPS)
    {
		// VS
		//gsampinf = malloc (genealogiestosave * sizeof (float *));
		//for (i = 0; i < genealogiestosave; i++)
		//  gsampinf[i] = malloc (gsampinflength * sizeof (float));
		gsampinf = (float ***) malloc ((nbgroupsloci_theta+nbgroupsloci_mig) * sizeof (float **));
		for(gp=0; gp < (nbgroupsloci_theta+nbgroupsloci_mig); gp++) 
		{
			gsampinf[gp] = (float **) malloc (genealogiestosave * sizeof (float *));
			for (i = 0; i < genealogiestosave; i++) 
			{
				gsampinf[gp][i] = (float *) malloc (gsampinflength * sizeof (float));
			}
		}	

		memforgenealogiessaved = genealogiestosave;
    }
    else                        /* cdurationmode == TIMEINF */
    {
		// VS
		//gsampinf = malloc (MAXGENEALOGIESTOSAVE * sizeof (float *));
		//for (i = 0; i < MAXGENEALOGIESTOSAVE; i++)
		//  gsampinf[i] = malloc (gsampinflength * sizeof (float));
		gsampinf = (float ***) malloc ((nbgroupsloci_theta+nbgroupsloci_mig) * sizeof (float **));	for(gp=0; gp < (nbgroupsloci_theta+nbgroupsloci_mig); gp++) 
		{
			gsampinf[gp] = (float **) malloc (MAXGENEALOGIESTOSAVE * sizeof (float *));
			for (i = 0; i < MAXGENEALOGIESTOSAVE; i++) 
			{
				gsampinf[gp][i] = (float *) malloc (gsampinflength * sizeof (float)); // VS gp
			}
		}

		memforgenealogiessaved = MAXGENEALOGIESTOSAVE;
    }
  }
	//AS: creating a local gsampinf 2D array
	z = whichiscoldchain();
	if (currentid != 0 && z >= 0) {
	
		gsampinflocal = (float **) malloc ((nbgroupsloci_theta+nbgroupsloci_mig) * sizeof (float *));	
		for(gp=0; gp < (nbgroupsloci_theta+nbgroupsloci_mig); gp++) 
		{
			gsampinflocal[gp] = (float *) malloc (gsampinflength * sizeof (float));
		}
	}

  if (genealogiessaved >= (MAXGENEALOGIESTOSAVE - 1) && cdurationmode == TIMEINF)
  {
    printf (" maximum possible genealogies saved \n");
    maxedoutgenealogysave = 1;
   #ifdef MPI_ENABLED
	MPI_Bcast(&maxedoutgenealogysave, 1, MPI_INT, currentid, MPI_COMM_WORLD);
   #endif //AS: TODO error check
  }
  else
  {
	//AS: only cold chain has to be saved!
	//AS: cold chain can exist on any processor! So if the cold chain NOT on the head node, 
	//that genealogy has to be saved to gsampinflocal, then sent to the head node to be collated
	//z = whichiscoldchain();
	if (currentid == 0) {
		if (z >= 0) {
    // VS
    //savegsampinf (gsampinf[genealogiessaved]);
	// This functions saves the gsampinfo for each group of loci 
	// savegsampinf (gsampinf[genealogiessaved]);
			for(gp=0; gp < nbgroupsloci_theta; gp++) {
				savegsampinf_vs (gsampinf[gp][genealogiessaved], gp, 1, z);
			}
			for(gp=0; gp < nbgroupsloci_mig; gp++) {
				savegsampinf_vs (gsampinf[gp+nbgroupsloci_theta][genealogiessaved], gp, 0, z);
			}
			//AS: to add z to savegsampinf_vs
		}
	} else if (currentid != 0) {
		if ( z >= 0) {//Indicates that the cold chain exists on this processor, which is not the head node
			for (gp=0; gp < nbgroupsloci_theta; gp++) {
				savegsampinf_vs (gsampinflocal[gp], gp, 1, z);
			}
			for (gp=0; gp < nbgroupsloci_mig; gp++) {
				savegsampinf_vs (gsampinflocal[gp+nbgroupsloci_theta], gp, 0, z);
			}
		}
	}
	//AS: Now onto writing send and receive functions between the cold chain node, head node
	#ifdef MPI_ENABLED
	if (currentid != 0 && numprocesses > 1 && z >= 0) {
		tempcurrentid = currentid;
		MPI_Send(&tempcurrentid, 1, MPI_INT, 0, 1212, MPI_COMM_WORLD);
		for (gp = 0; gp < nbgroupsloci_theta+nbgroupsloci_mig; gp++) {
			for (v = 0; v < gsampinflength; v++) {
				MPI_Send(&gsampinflocal[gp][v], 1, MPI_FLOAT, 0, v*2, MPI_COMM_WORLD);
			}
		}
	} else if (currentid == 0 && z < 0 && numprocesses > 1) {
		tempcurrentid = 0;
		MPI_Recv(&tempcurrentid, 1, MPI_INT, MPI_ANY_SOURCE, 1212, MPI_COMM_WORLD, &status);
		for (gp = 0; gp < nbgroupsloci_theta + nbgroupsloci_mig; gp++) {
			for (v = 0; v < gsampinflength; v++) {
				MPI_Recv(&gsampinf[gp][genealogiessaved][v], 1, MPI_FLOAT, tempcurrentid, v*2, MPI_COMM_WORLD, &status);
			}
		}
	}
	//Free gsampinflocal on other processors
	if (currentid != 0 && numprocesses > 1 && z >= 0) {
		for (gp = 0; gp < nbgroupsloci_theta + nbgroupsloci_mig; gp++) {
			XFREE(gsampinflocal[gp]);
		}
	}
	#endif	
  }
}                               /* savegenealogyinfo */



//AS: TO BE CONTINUED - 9/25/2014
//AS: Continuing - 9/26/2014

// VS 
// SAVEASSIGNLOCIINFO
// Allocates memory for the assignmene t loci vectors if it is the the first iteration
// and saves the current assignment vector for theta and migration parameters into the array assignlocisample_mig and assignlocisample_theta
// AS: SAVEASSIGNLOCIINFO also has to come under the MPI framework
// AS: assignment vectors on a cold chain that don't exist on the head node will be saved locally, then sent back to the head node

void saveassignlociinfo(int currentid) {
  int i;
  int z; //AS: coldchain id
  int l; //AS: indexing loci
  //AS: local arrays
  int *assignlocisample_miglocal;
  int *assignlocisample_thetalocal;
  #ifdef MPI_ENABLED
  MPI_Status status;
  #endif
  if (genealogiessaved == 0)
  {
    if (cdurationmode == TIMESTEPS)
    {
		assignlocisample_mig = (int **) malloc (genealogiestosave  * sizeof (int *));
		assignlocisample_theta = (int **) malloc (genealogiestosave  * sizeof (int *));
		for (i = 0; i < genealogiestosave; i++) {
			assignlocisample_mig[i] = (int *) malloc (nloci * sizeof (int));
			assignlocisample_theta[i] = (int *) malloc (nloci * sizeof (int));
		}	
    }
    else                        /* cdurationmode == TIMEINF */
    {

		assignlocisample_mig = (int **) malloc (MAXGENEALOGIESTOSAVE   * sizeof (int *));
		assignlocisample_theta = (int **) malloc (MAXGENEALOGIESTOSAVE   * sizeof (int *));
		for (i = 0; i < MAXGENEALOGIESTOSAVE; i++) {
			assignlocisample_mig[i] = (int *) malloc (nloci * sizeof (int));
			assignlocisample_theta[i] = (int *) malloc (nloci * sizeof (int));
		}
    }
  }
  //AS: Allocate space for local assignment vector if it's not the head node
  z = whichiscoldchain();
  #ifdef MPI_ENABLED
  if (currentid != 0 && z >= 0) {
	assignlocisample_miglocal = (int *) malloc (nloci * sizeof (int));
	assignlocisample_thetalocal = (int *) malloc (nloci * sizeof (int));
	saveassignloci (assignlocisample_miglocal, grouploci_mig[z]); //AS: z corresponds to cold chain
	saveassignloci (assignlocisample_thetalocal, grouploci_theta[z]); //AS: z corresponds to cold chain
	for (l = 0; l < nloci; l++) {
		MPI_Send(&assignlocisample_miglocal[l], 1, MPI_INT, 0, 1234, MPI_COMM_WORLD);
		MPI_Send(&assignlocisample_thetalocal[l], 1, MPI_INT, 0, 4567, MPI_COMM_WORLD);
	}
	XFREE(assignlocisample_miglocal);
	XFREE(assignlocisample_thetalocal);
  }
  else if (currentid == 0 && z < 0) {
	for (l = 0; l < nloci; l++) {
		MPI_Recv(&assignlocisample_mig[genealogiessaved][l], 1, MPI_INT, MPI_ANY_SOURCE, 1234, MPI_COMM_WORLD, &status);
		MPI_Recv(&assignlocisample_theta[genealogiessaved][l], 1, MPI_INT, MPI_ANY_SOURCE, 4567, MPI_COMM_WORLD, &status);
	}
  }
  #endif
  // save the assignment vectors 
  if (currentid == 0 && z >= 0) {
	  saveassignloci (assignlocisample_mig[genealogiessaved], grouploci_mig[z]); // grouploci_mig[0] corresponds to the assignment vector of chain 0
	  saveassignloci (assignlocisample_theta[genealogiessaved], grouploci_theta[z]); // grouploci_theta[0] corresponds to the assignment vector of chain 0
  }

}



/* SANGCHUL: This could be simplified considering its function.
 * */
// VS - modified this function to read the TI files from different groups of loci
void loadgenealogyvalues (void)
{
  char filenamewildcard[FNSIZE];
  char *ctp, *textline, *dataline, *c, tempc;
  int charspervalue = 12;
  FILE *sfile;
  int i, j, numgenealogies, totalnumgenealogies;
  /* int filefound, nofile; */
  char *defaultdir; 
  struct dirent *dir_entry;
  DIR *dp;
  int numfiles = 0;
  int numtoload, loadall, loaded, notloaded;
  int len_base;
  int l2;
  int len_defaultdir;
  char *basename;
  char *treefilename;
  int gp, m, aux_proceed, laux; // VS added these variables
  // VS - created an array of strings that save the extension of the name of the files
  // used to read the different TI files created for each group of loci
  char **strnfilename;
  float load_fraction;

  // VS allocate memory to strnfilename
  strnfilename = (char **) malloc(sizeof(char *) * (nbgroupsloci_theta+nbgroupsloci_mig)); // each index of strnfilename points to one string
  for(m=0; m<(nbgroupsloci_theta+nbgroupsloci_mig); m++) {
	strnfilename[m] = (char *) malloc(sizeof(char) * 17); // each string has a maximum size of 20 characters
  }

  if (strlen (loadfilebase))
  {
    strcpy (filenamewildcard, loadfilebase);
  }
  else
  {
    strcpy (filenamewildcard, outfilename);
    strtrunc (filenamewildcard, (char) '.');
    strtrunc (filenamewildcard, (char) '-');
  }
  // VS - changed this to *.ti.groupTheta0
  //strcat (filenamewildcard, "*.ti");
  strcat (filenamewildcard, ".ti.groupTheta0");
  SP "\nLOAD TREES (L) MODE INFORMATION\n");
  SP "============================================================================\n");
  SP "  Base filename for loading files with sampled genealogies: %s\n", filenamewildcard);
  SP "  Files loaded with sampled genealogies:\n");
  numgenealogyfiles = 0;
  textline = malloc (300 * sizeof (char));
  dataline = malloc (gsampinflength * charspervalue * sizeof (char));
  ctp = &textline[0];
  numgenealogies = totalnumgenealogies = 0;

  /* Find the default directory based on loadfilebase.
   * /this/directory/a.out -> defaultdir is /this/directory/
   * /a.out                -> /
   * a.out                 -> ./
   */
  imaDirBase (loadfilebase, &defaultdir); 
  len_defaultdir = strlen (defaultdir);
  
  if ((dp = opendir (defaultdir)) == NULL)
  {
    IM_err (IMERR_TIFILE, "cannot open directory: %s", defaultdir);
  }
  basename = strrchr (loadfilebase, '/');
  if (basename == NULL)
    {
      basename = strrchr (loadfilebase, '\\');
      if (basename == NULL)
        basename = loadfilebase;
      else
        basename++;
    }
  else
    {
      basename++;
    }
  len_base = strlen (basename); // VS len_base contains the number of characters of base length
  while ((dir_entry = readdir(dp)) != NULL)
  {
    if (!strncmp(dir_entry->d_name, basename, len_base)) 
    {
      l2 = strlen (dir_entry->d_name);
	  // VS -changed here the file extension so that it only reads the files
	  // for the grouptheta0. The aim os this loop is to find out how many TI files are there
	  // if we have multiple loci, we know that, we should have an equal number of files for 
	  // theta0, theta1, mig2, and so on. So, since this loop looks into all files and count the number
	  // of genealogies, it is just reading the theta0 group.
	  // One way to check if this is correct is to verify that all the other groups have the same
	  // number of genealoagies saved
      //if (!strcmp (&dir_entry->d_name[l2 - 3], ".ti"))
	  if (!strcmp (&dir_entry->d_name[l2 - 15], ".ti.groupTheta0"))
      {  
        numgenealogyfiles++; 
        /* We found one. */
        treefilename = malloc ((len_defaultdir + l2 + 1) * sizeof (char));
        sprintf (treefilename, "%s%s", defaultdir, dir_entry->d_name);

        /* Count the number of gene genealogies of the found file. */
        if ((sfile = fopen (treefilename, "r")) == NULL) 
        {
          IM_err (IMERR_TIFILE, " cannot open .ti file");
        }
        while (fgets (textline, 300, sfile)
               && strstr (textline, "VALUESSTART") == NULL && !feof (sfile));
        numgenealogies = 0;
        while ((tempc = fgetc (sfile)) != EOF)    // count lines
        {
          numgenealogies += (tempc == '\n');
        }
        if (numgenealogies < 1)
        {
          printf ("  *no genealogies loaded from file %s\n", dir_entry->d_name);
          SP "  *no genealogies loaded from file %s\n", dir_entry->d_name);
        }
        else
        {
          printf ("  loaded %d genealogies from genealogy file  %s\n", numgenealogies,
                  dir_entry->d_name);
          SP "  loaded %d genealogies from genealogy file  %s\n", numgenealogies,
            dir_entry->d_name);
        }
        fclose(sfile);
        totalnumgenealogies += numgenealogies;

        XFREE (treefilename);
      }
    }
  }
  closedir(dp);
  // VS - the above loop will look into all files in the directory with the extension *.TI
  // it will count how many genealogies are there in total
  // This is not necessarily useful for the groups of loci implementation
  // Unless we change the above for each group, instead of looking for *.TI files
  // look for *.TI.THETA0 files. Count how many genealogies that file has and that will be
  // the number of genealogies in all files.
  /* We have counted gene genealogies. */

  if (genealogiestosave > 0)
    numtoload = IMIN (totalnumgenealogies, genealogiestosave);
  else
    numtoload = totalnumgenealogies;
  numtoload = IMIN (numtoload, MAXGENEALOGIESTOSAVE);
  memforgenealogiessaved = numtoload;

  // VS - allocate memory for gsampinf taking into account the groups of loci
  // allocate numtoload genealogies to each index of the gsampinf matrix
  //gsampinf = malloc (numtoload * sizeof (float *));
  gsampinf = (float ***) malloc ((nbgroupsloci_theta+nbgroupsloci_mig) * sizeof (float **));	
  for(gp=0; gp < (nbgroupsloci_theta+nbgroupsloci_mig); gp++) 
  {
	gsampinf[gp] = (float **) malloc (numtoload  * sizeof (float *));
  } // gsamping with groups of loci is seen as a 3D matrix



  loadall = numtoload >= totalnumgenealogies;
  load_fraction = (float) numtoload/ (float) totalnumgenealogies;
  loaded = 0;
  notloaded = 0;
  SP "\n  HEADER INFORMATION FROM FIRST GENEALOGY FILE\n");
  SP "  ---------------------------------------\n");
// now go through again and save the genealogies
  /* closedir (dp); */

// VS - copied this section into inside the next for loop
//  if ((dp = opendir (defaultdir)) == NULL)
//  {
//    IM_err (IMERR_TIFILE, "cannot open directory: %s", defaultdir);
//  }
//  numfiles = 0;

  // VS - here it is going again through all the files 
  // the easiest temporary solution is just to repeat this loop for each group of loci
  // changing the strcmp(&dir_entry->d_name[l2 - 3], ".ti")
  // for the corresponding groups
  // Need to create a for loop that will go through the number of groups
  // for each group it creates a string with the final of the file name

  // Create and array of strings with the file name
  for(m=0; m<nbgroupsloci_theta; m++) {
		sprintf(strnfilename[m], "%s%i", ".ti.groupTheta", m);
  }
  for(; m<(nbgroupsloci_theta+nbgroupsloci_mig); m++) {
		sprintf(strnfilename[m], "%s%i", ".ti.groupMig", m);
  }

  // VS - for loop through the TI files, one TI file for each group of loci
  for(m=0; m<(nbgroupsloci_theta+nbgroupsloci_mig); m++) {
	
	  if ((dp = opendir (defaultdir)) == NULL)
	  {
		IM_err (IMERR_TIFILE, "cannot open directory: %s", defaultdir);
	  }
	  numfiles = 0;

	  aux_proceed = 0;	// VS - reset the proceed
	  loaded = 0; // VS - reset the number of loaded genealogies
	  notloaded = 0; // VS - reset the number of notloaded genealogies
	  // VS given that the fraction is used, if this is not done, then the number of genealogies loaded will be wrong
	  
	  while ((dir_entry = readdir(dp)) != NULL)
	  {
		if (!strncmp(dir_entry->d_name, basename, len_base)) 
		{
		  l2 = strlen (dir_entry->d_name);
		  // VS - compare with the string for the TI file
		  //if (!strcmp (&dir_entry->d_name[l2 - 3], ".ti"))
		  // VS - This if statement needs to be changed if m>10
		  // alse, need to consider if we have a Theta group or Mig group,
		  // because the number of characters to compare will be different
		  // to do this created an auxiliary variable aux_proceed
		  // this takes the value 0 zero if we should not go through the IF statement
		  // and if takes the value 1 if we should go through the IF statement
		  // also, created the variable laux, which saves the size of the strnfilename
		  // for each m index
		  if(m<10) {
			  if(m<nbgroupsloci_theta) {
				laux = 15; // VS need to decrease char length of ".ti.groupTheta0"
				aux_proceed=1; 
			  }
			  else {
				laux = 13; // VS need to decrease char length of ".ti.groupMig0"
				aux_proceed=1;
			  }	
		  } else {
			  if(m<nbgroupsloci_theta) {
				laux = 16; // VS need to decrease char length of ".ti.groupTheta10"
				aux_proceed=1; 
			  }
			  else {
				laux = 14; // VS need to decrease char length of ".ti.groupMig10"
				aux_proceed=1;
			  }	
		  }

		  if((aux_proceed==1) && (!strncmp(&dir_entry->d_name[l2 - laux], strnfilename[m], laux))) 
		  {  
			numfiles++;
			/* We found one. */
			treefilename = malloc ((len_defaultdir + l2 + 1) * sizeof (char));
			sprintf (treefilename, "%s%s", defaultdir, dir_entry->d_name);

			/* Count the number of gene genealogies of the found file. */
			if ((sfile = fopen (treefilename, "r")) == NULL) 
			{
			  IM_err (IMERR_TIFILE, " cannot open .ti file");
			}

			while (fgets (textline, 300, sfile) && strstr (textline, "VALUESSTART") == NULL && !feof (sfile))
			  if (numfiles == 1)
			  {
				SP "  ||%s", textline);
			  }
			while (fgets (dataline, gsampinflength * charspervalue, sfile)
				   && !feof (sfile))
			{
			  if (loadall || (loaded == 0)
				  ||  ((float) loaded/ (float) (loaded + notloaded)) <= load_fraction)
				  //JH 7/24/09  change this to using a proportion((float) loaded / (float) notloaded) <= load_notload_ratio)
			  {

                // VS - Allocate memory for each row of gsampinf
				// Need to change this when dealing with groups of loci
			    // in that case gsampinf is a 3D matrix
				//gsampinf[loaded] = malloc (gsampinflength * sizeof (float));
				gsampinf[m][loaded] = (float *) malloc (gsampinflength * sizeof (float));
				c = dataline;
				for (i = 0; i < gsampinflength; i++)
				{
				  // VS - modified gsampinf for a 3D matrix
				  // sscanf (c, "%f", &gsampinf[loaded][i]);
				  sscanf (c, "%f", &gsampinf[m][loaded][i]);
				  j = allwhitespace (c);
				  if (j ==1 || j== -1)
					IM_err (IMERR_TIFILE, "Problem in .ti file %s,  too few values per genealogy, .ti file may have been generated with a different program",treefilename);
				  c = nextwhite (c);
				}
				j = allwhitespace (c);
				if (j==0)
				   IM_err (IMERR_TIFILE, "Problem in .ti file %s,  too many values per genealogys, .ti file may have been generated with a different program",treefilename);
				loaded++;
			  }
			  else
				notloaded++;
			  /* gcounter++; */
			}
			// VS added this printf to screen
			printf ("  loaded %d genealogies from genealogy file  %s\n", loaded,treefilename);

			fclose (sfile);

			XFREE (treefilename);
		  }
		}
	  }
	  closedir(dp);

	  SP "  END OF HEADER INFORMATION FROM FIRST GENEALOGY FILE\n");
	  SP "  ----------------------------------------------\n\n");
	  
	  SP "  Number of files loaded : %d  total number of genealogies used: %d out of a total of: %d \n", numgenealogyfiles, loaded, totalnumgenealogies);
	  printf("  Number of files loaded : %d  total number of genealogies used: %d out of a total of: %d \n", numgenealogyfiles, loaded, totalnumgenealogies);
	  if (numgenealogies < 1)
		IM_err (IMERR_TIFILE, "  no genealogies loaded from .ti file(s)");

	 

  } // VS - end of loop through TI files for each group of loci

  /* closedir (dp); */
  genealogiessaved = loaded;
  for (j = 0; j < genealogiessaved; j++)
  {
    // use full range of t, ignore t.pr.min > 0
    if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 0
        && npops > 1)
    {
      for (i = 0; i < lastperiodnumber; i++)
      {
        recordval (T[i].v, gsampinf[0][j][gsamp_tp + i]); // VS set gp=0. Note that all gsampinf matrices share the same tsplit values in the last column
      }
    }
    if (!outputoptions[DONTPRINTASCIITREND])
      trendrecord (j, 0); //AS: since this is only called by the head node, i am directly passing 0 (currentid)
  }

  // VS - free the variable with the names of the files for each group of loci
  for(m=0; m<(nbgroupsloci_theta+nbgroupsloci_mig);m++) {
	XFREE(strnfilename[m]);
  }
  XFREE(strnfilename);

  XFREE (textline);
  XFREE (dataline);
  XFREE (defaultdir); 
}                               /* modified VS based on sang chul's loadgenealogyvalues */
//AS: I am not entirely sure that this would need any modification, since I'll just call it on the 
//AS: (cont) head node anyway. Talk to VS about this - 9/26/2014

void printsteps (FILE * outto, double like)
{
  int ci;

  ci = 0;
  if (!burndone)
  {
    fprintf (outto, "=BURNIN-PERIOD===============================\n");
    fprintf (outto, "STEP # %d  p(D|G): %.3lf p(G): %.3lf\n", step,
             like, C[ci]->allpcalc.probg);
  }
  else
  {
    fprintf (outto, "=============================================\n");
    if (genealogiessaved > 0)
      fprintf (outto,
               "STEP # %ld # Genealogies Saved: %d p(D|G): %.1lf p(G): %.1f\n",
               step - burnsteps, genealogiessaved, like, C[ci]->allpcalc.probg);
    else
      fprintf (outto, "STEP # %ld  p(D|G): %.3lf p(G): %.3lf\n",
               step - burnsteps, like, C[ci]->allpcalc.probg);
  }
  return;
}

/* To print acceptance rates:
	reclist[] is an array of pointers to struct chainstate_record_updates_and_values
	set values of reclist[] to those structures for which you want to print acceptance rates
	call printacceptancerates ()
	*/

void callprintacceptancerates (FILE * outto)
{
  int i, j, li;
  // length of this array must be fairly long, although it is technically possible to have MAXLOCI * MAXLINKED records,  but very unlikely
  struct chainstate_record_updates_and_values *reclist[MAXLOCI + MAXLINKED];
// t values 
  for (i = 0; i < numsplittimes; i++)
    reclist[i] = (T + i);

  if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 1
      || npops == 1)
  {
    assert (numsplittimes == 0);
    /* no split time */
  }
  else
  {
    printacceptancerates (outto, numsplittimes, reclist,
                          "Update Rates -- Population Splitting Times");
  }
// genealogy updates
  for (li = 0; li < nloci; li++)
    reclist[li] = L[li].g_rec;
  printacceptancerates (outto, nloci, reclist, "Update Rates -- Genealogies");
// mutation rate scalars
  if (nurates > 1
      && (runoptions[PRINTMUTATIONUPDATESTOSCREEN] || outto != stdout))
  {
    for (i = 0, li = 0; li < nloci; li++)
      for (j = 0; j < L[li].nlinked; j++)
      {
        reclist[i] = &L[li].u_rec[j];
        i++;
      }
// kappa values for HKY model
    printacceptancerates (outto, i, reclist,
                          "Update Rates -- Mutation Rate Scalars");
    for (i = 0, li = 0; li < nloci; li++)
      for (j = 0; j < L[li].nlinked; j++)
      {
        if (L[li].umodel[j] == HKY)
        {
          reclist[i] = L[li].kappa_rec;
          i++;
        }
      }
    if (i > 0)
      printacceptancerates (outto, i, reclist,
                            "Update Rates -- HKY Model Kappa parameter");
// STR ancestral allele states 
    for (i = 0, li = 0; li < nloci; li++)
      for (j = 0; j < L[li].nlinked; j++)
      {
        if (L[li].umodel[j] == STEPWISE)
        {
          reclist[i] = &L[li].A_rec[j];
          i++;
        }
      }
    if (i > 0)
    {
      printacceptancerates (outto, i, reclist,
                            "Update Rates -- STR Genealogy Allele States");
    }
  }
  if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
  {
    for (li = 0; li < nloci; li++)
      reclist[li] = L[li].a_rec;

    if (assignmentoptions[POPULATIONASSIGNMENTLOCAL] == 1)
    {
      printacceptancerates (outto, nloci, reclist,
                            "Update Rates -- Assignment");
    }
    // assignment updating information 
    // set reclest
    // call printacceptancerates 
    //
    printacceptancerates_multichain (outto);
  }

  // Assignment of loci update rates 
  // There are only two assignloci vectors (one for migration, and one for theta groups)
  for (i = 0; i < 2; i++)
    reclist[i] = (assignloci + i);
  printacceptancerates (outto, 2, reclist, "Update Rates -- Assignment of loci");


  return;
}                               //callprintacceptancerates


/* set up arrays pointing to information to put in curve plots,  then call asciicurve
   some things that are plotted are based on struct value_record and others on struct i_param
   this is why we cannot simply call asciicurve() with a pointer to a single type of structure  */
void callasciicurves (void)
{
	// int gp; // VS gp gp not used anymore
  struct plotpoint **curvexy;
  char **curvestr;
  int *curve_do_logplot;
  int numcurve = 0;
  int i, j, li, ui;
  int *nrecstep;
  //_CrtCheckMemory( );
// find out how many curves

  int nbparams, nbthetaparams, nbmigparams; // VS

  // VS compute some variables that will simplify the code in the for loops and allocation of variables
  nbparams=(numpopsizeparams*nbgroupsloci_theta)+(nummigrateparams*nbgroupsloci_mig); 
  nbthetaparams=numpopsizeparams*nbgroupsloci_theta;
  nbmigparams=nummigrateparams*nbgroupsloci_mig;

  numcurve += nbthetaparams; // VS numpopsizeparams
  for (i = nbthetaparams; i <nbparams; i++) // VS nummigrateparams
    if (imig[g_paramindex[i]].pr[g_gp_ip[i]].max > MPRIORMIN) // gp
      numcurve++;
  numcurve += numsplittimes;
  if (runoptions[LOADRUN] == 0 && nurates > 1)
    numcurve += nurates;
  if (runoptions[LOADRUN] == 0)
    for (li = 0; li < nloci; li++)
      if (L[li].model == HKY)
        numcurve++;
// allocate
  curvexy = malloc (numcurve * sizeof (struct plotpoint *));
  curvestr = malloc (numcurve * sizeof (char *));
  nrecstep = malloc (numcurve * sizeof (int));
  curve_do_logplot = malloc (numcurve * sizeof (int));
// assign
  j = 0;
  for (i = 0; i < nbthetaparams; i++) // VS numpopsizeparams
  {
    curvexy[j] = itheta[g_paramindex[i]].xy[g_gp_ip[i]]; // VS g_paramindex[i], g_gp_ip[i]
    curvestr[j] = &itheta[g_paramindex[i]].str[0]; // VS g_paramindex[i]
    curve_do_logplot[j] = 0;
    nrecstep[j] = 1;
    j++;
  }
  for (i = nbthetaparams; i < nbparams; i++) // VS nummigrateparams
    if (imig[g_paramindex[i]].pr[g_gp_ip[i]].max > MPRIORMIN) // VS g_paramindex[i], g_gp_ip[i]
    {
      curvexy[j] = imig[g_paramindex[i]].xy[g_gp_ip[i]]; // VS g_paramindex[i], g_gp_ip[i]
      curvestr[j] = &imig[g_paramindex[i]].str[0]; // VS g_paramindex[i]
      curve_do_logplot[j] = 0;
      nrecstep[j] = 1;
      j++;
    }
  if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 1
      || npops == 1)
  {
    /* no split time */
  }
  else
  {
    for (i = 0; i < lastperiodnumber; i++)
    {
      curvexy[j] = T[i].v->xy;
      curvestr[j] = &T[i].v->str[0];
      curve_do_logplot[j] = T[i].v->do_logplot;
      nrecstep[j] = recordstep;
      j++;
    }
  }
  if (runoptions[LOADRUN] == 0 && nurates > 1)
  {
    for (li = 0; li < nloci; li++)
      for (ui = 0; ui < L[li].nlinked; ui++)
      {
        curvexy[j] = L[li].u_rec[ui].v->xy;
        curvestr[j] = &L[li].u_rec[ui].v->str[0];
        curve_do_logplot[j] = L[li].u_rec[ui].v->do_logplot;
        nrecstep[j] = recordstep;
        j++;
      }
  }
  if (runoptions[LOADRUN] == 0)
    for (li = 0; li < nloci; li++)
      if (L[li].model == HKY)
      {
        curvexy[j] = L[li].kappa_rec->v->xy;
        curvestr[j] = &L[li].kappa_rec->v->str[0];
        curve_do_logplot[j] = L[li].kappa_rec->v->do_logplot;
        nrecstep[j] = recordstep;
        j++;
      }
  assert (numcurve == j);
  for (j = 0; j < numcurve; j++)
    asciicurve (outfile, curvexy[j], curvestr[j], curve_do_logplot[j],
                nrecstep[j]);
  // VS - TO BE DONE LATER maybe here it can have added in the curvestr the group of loci at some point


//free
  XFREE (curvexy);
  XFREE (curvestr);
  XFREE (curve_do_logplot);
  XFREE (nrecstep);
}                               //callasciicurve 

// makes calls to asciitrend
void callasciitrend (FILE * outtofile)
{
  int i, li, ui;
  asciitrend (outtofile, lpgpd_v, trenddoublepoint, trendspot);

  if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 1
      || npops == 1)
  {
    /* no split time */
  }
  else
  {
    for (i = 0; i < lastperiodnumber; i++)
      asciitrend (outtofile, T[i].v, trenddoublepoint, trendspot);
  }
  if (nurates > 1 && runoptions[LOADRUN] == 0)
  {
    for (li = 0; li < nloci; li++)
      for (ui = 0; ui < L[li].nlinked; ui++)
        asciitrend (outtofile, L[li].u_rec[ui].v, trenddoublepoint, trendspot);
  }
  for (li = 0; li < nloci; li++)
    if (L[li].model == HKY && runoptions[LOADRUN] == 0)
      asciitrend (outtofile, L[li].kappa_rec->v, trenddoublepoint, trendspot);

  if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
  {
    for (li = 0; li < nloci; li++)
    {
      asciitrend (outtofile, L[li].a_rec->v, trenddoublepoint, trendspot);
    }
  }
  return;
}                               // callasciitrend 



// printoutput
// prints the results for the screen and for the output file
// calls many important functions to make histograms, tables, and statistical tests
// including the LRT (likelihood ratio tests)
// AS: adding currentid to this
void 
printoutput (int currentid)         // mostly calls functions in output.c
{
  int i;
  int z = 0; //AS: coldchain
  double seconds;
  int p;
  float *holdpeakloc;
  double multitpeak[MAXPOPS - 1];
  int gp, paramindex; // VS
  int nbparams, nbthetaparams, nbmigparams; // VS need to check this (style) Variable 'nbmigparams' is assigned a value that is never used 

  // VS compute some variables that will simplify the code in the for loops and allocation of variables
  nbparams=(numpopsizeparams*nbgroupsloci_theta)+(nummigrateparams*nbgroupsloci_mig); 
  nbthetaparams=numpopsizeparams*nbgroupsloci_theta;
  nbmigparams=nummigrateparams*nbgroupsloci_mig;

  // VS initialize the global variables with the indexes of group and param indexes
  // depending on the index i in the output printing functions.
  // In most output functions, the parameters are put all together in as a list (or vector).
  // Thus, the index i for the parameters varies between 0<i<(nbgroupsloci_theta*numpopsizeparams)+(nbgroupsloci_mig*nummigrateparams).
  // But, the index i is not the same as in other structures.
  // The reference to a given parameter depends on whether we 
  // are in the context of iparam structures, gsampinf matrix or printout functions.
  // The solution was to create global variables which are used as lookup tables.
  // Given the value of i, looking at these global variables we can obtain
  // the corresponding index in the different contexts.
  // Example of a case with 2 groups of theta and 2 groups of mig
  //		i	g_param	g_gp	g_gp_gsampinf
  // q0_g0	0	0		0		0
  // q1_g0	1	1		0		0
  // qA_g0	2	2		0		0
  // q0_g1	3	0		1		1
  // q1_g1	4	1		1		1	
  // qA_g1	5	2		1		1
  // m01_g0	6	0		0		2
  // m10_g0	7	1		0		2
  // m01_g1	8	0		1		3
  // m10_g1	9	1		1		3
if (currentid == 0) {
  for(i=0; i<nbparams; i++) {
	  // if we are refering to theta parameters
	  if(i<nbthetaparams){ 
		for(gp=0; gp<nbgroupsloci_theta; gp++) { // go through the groups of theta
		  if(i >= gp*numpopsizeparams && i < (gp+1)*numpopsizeparams) {
				paramindex = i-(gp*numpopsizeparams);
				break; // as soon as it finds the paramindex breaks, out of for loop. Hence the gp is the correct gp
		  }
		}
		// initialize the global variables
		g_paramindex[i] = paramindex; // index for param in itheta
		g_gp_ip[i] = gp; // index of group in itheta
		g_gp_gs[i] = gp; // index of group in gsampinf
	  }
	  else {
		// if we are referring to migration parameters
		for(gp=0; gp<nbgroupsloci_mig; gp++) {
		  if(i >= (gp*nummigrateparams)+nbthetaparams && i < ((gp+1)*nummigrateparams)+nbthetaparams) {
				paramindex = i-((gp*nummigrateparams)+nbthetaparams);
				break; // as soon as it finds the param index breaks, out of for loop. Hence the gp is the correct gp
		  }
		}
		// initialize the global variables
		g_paramindex[i] = paramindex; // index for param in imig
		g_gp_ip[i] = gp; // index of group in imig
		g_gp_gs[i] = gp+nbgroupsloci_theta; // index of group in gsampinf
	  }
  }

  // VS CHECK print the tables with the indeces and check if they are correct
  printf("\n-----------------------------------\n");
  printf("\nCheck Global Variables with Indeces\n");
  for(i=0; i<nbparams; i++) {
	printf("%i, %i, %i, %i\n", i, g_paramindex[i], g_gp_ip[i], g_gp_gs[i]);
  }
  printf("-----------------------------------\n");




  if (runoptions[LOADRUN] == 0)
  {
    /* savegenealogyfile moved out of printoutput because it does not seem to belong to it. */
    if (cdurationmode == TIMEINF)
    {
      remove (oldoutfilename);
      rename (outfilename, oldoutfilename);
    }
  }
 //AS: closes if currentid == 0
//AS: now onto writing the run basics from the processor that has the cold chain
  //z = whichiscoldchain();
//  if (z >= 0) {
	  if ((outfile = fopen (outfilename, "w")) == NULL)
	  {
	    IM_err (IMERR_CREATEFILEFAIL, "Error opening text file for writing");
	  }
	  printrunbasics (outfile, runoptions[LOADRUN], fpstr, burnsteps, recordint,
                  recordstep, savegenealogyint, endtime, starttime, hilike,
                  hiprob/*, step*/, 0);
  
  if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
  {
    imaAsnPrintNumGenesPopn (outfile);
  }
 } //closes if currentid == 0
  // Bayes factor stuff of Sang Chul's fprintf (outfile, "Average loglikelihood: %lf\n", gloglikelihood);
  if (runoptions[LOADRUN] == 0)
  {
   if (currentid == 0) {
    callprintacceptancerates (outfile);
   }
    if (numchains * numprocesses > 1 && runoptions[LOADRUN] != 1)
      printchaininfo (outfile, heatmode, hval1, hval2, currentid); //AS: adding currentid TODO: add currentid to printchaininfo 9/26/2014
    	if (currentid == 0) {
	callprintautoctable (outfile/*, step*/);
	}
  }
  if (currentid == 0) {
  if (outputoptions[PARAMGREATERTHAN] && currentid == 0)
  {
    print_greater_than_tests (outfile);
  }
  if (!modeloptions[EXPOMIGRATIONPRIOR] && calcoptions[DONTCALCLIKELIHOODMUTATION]==0)        //as of 11/19/09 have not yet done the math for case of migration with exponential prior
	  print_means_variances_correlations (outfile);


  // init_surface_calc (); not used as of 8/24/09
/*  get marginal peaks */
  if (!calcoptions[DONTCALCLIKELIHOODMUTATION])
  {

	vs_test=1; // VS 10/11/2010 vs_test
    // VS modified the number of parameters
	//p = numpopsizeparams + nummigrateparams ;
	p = nbparams; // VS
    holdpeakloc = malloc (p * sizeof (float));
    printf ("surface calculations . . .\n");
    if (modeloptions[EXPOMIGRATIONPRIOR] || (runoptions[LOADRUN] && calcoptions[FINDJOINTPOSTERIOR]))
      eexpsum = (struct extendnum *) malloc ((size_t) ((genealogiessaved + 1) * sizeof (struct extendnum)));
    // VS - the variable eexpsum is used to compute the sum of exponentials, with the notation aEb
    findmarginpeaks (outfile, holdpeakloc);

	// VS 11/8/2011 commented this section to get only the marginal posteriors
    if (runoptions[LOADRUN] && calcoptions[FINDJOINTPOSTERIOR])
    {
      closeopenout (&outfile, outfilename);
      // CR 110921.1  Change outfile parameter type to (FILE **) 
      findjointpeaks(&outfile,outfilename,nestedmodelfilename,p);
    }

  }
  if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
  {
    imaAsnPrintNumGenesPopn (outfile);
  }
  else
  {
/* get joint splittime peak */
    if (npops >=3  && npops <= 5 && !runoptions[LOADRUN] && outputoptions[PRINTJOINTTEST])
    {
      return_joint_t (multitpeak);
      FP "\nEstimated joint splitting time from multi-dimensional histogram\n");
      FP "  number of bins per dimension %d\n", NUMTARRAYBINS);
      FP "  Posterior probability of estimated joint value of splitting time: %7.4lf\n", joint_t_prob (&multitpeak[0]));
      FP "---------------------------------------------------------------\n");
      for (i = 0; i < numsplittimes; i++)
        FP "   %s\t%.3lf\n", T[i].str, multitpeak[i]);
      FP "\n\n");
    }
  }


  // VS - disabled the option to print histograms and callasciicurves
  printhistograms (outfile, recordstep, generationtime, scaleumeaninput);
  if (!outputoptions[DONTPRINTASCIICURVE])
  {
    FP "\n\nASCII Curves - Approximate Posterior Densities \n");
    FP "===================================================\n");
	callasciicurves ();
  }
  if (!outputoptions[DONTPRINTASCIITREND])
  {
    if (trendspot <= 1)
    {
      FP "burn period too short to plot trends,  trend recording begins at step %d \n",burntrendstartdelay); 
    }
    else
    {
      FP "\n\nASCII Plots of Parameter Trends \n");
      FP "===================================================\n");
      FP " - note points to the left of '!' on X axis have twice the density in time relative to points to the right\n\n");
	  callasciitrend (outfile);
    }

  }
  time (&endtime);
  seconds = difftime (endtime, starttime);
  FP "Time Elapsed : %d hours, %d minutes, %d seconds \n\n",
    (int) seconds / (int) 3600,
    ((int) seconds / (int) 60) - ((int) 60 * ((int) seconds / (int) 3600)),
    (int) seconds - (int) 60 *((int) seconds / (int) 60));
  FP "\nEND OF OUTPUT\n");
  f_close (outfile);
  //free_surface_calc ();

  if (outputoptions[MIGRATEHIST] && nummigrateparams > 0)
  {

    if ((migplotfile = fopen (migplotfilename, "w")) == NULL)
    {
      IM_err (IMERR_CREATEFILEFAIL,
              "Error opening file for plotting migration amounts and times");
    }
    printmigrationhistograms (migplotfile, recordstep);
    f_close (migplotfile);
  }


  if (!calcoptions[DONTCALCLIKELIHOODMUTATION])
    XFREE (holdpeakloc);
  if (runoptions[SAVEMCSTATEFILE])
  {
    writemcf (mcfwritefilename);
  }
  if (jheyoptions[WRITEMIGRATIONNAME])
  {
    f_close(migrationnamefile);
    migrationnamefile = fopen (migrationnamefilename, "a");
  }
}
  return;
}                               /* printoutput */


void intervaloutput (FILE * outto, int currentid)
{
  int li;
  int j;
  double like;
  int ci;
//debug
  double multitpeak[MAXPOPS - 1];
  int z = 0; //coldchain

//AS: stuff added for printing Tue Feb 23 16:40:48 EST 2016
  #ifdef MPI_ENABLED
  	MPI_Status status;
  #endif

  z = whichiscoldchain();
  ci = 0;
  if(z >= 0)
	  checkhighs (z, printint, &hilike, &hiprob, &like/*, step*/);
  if (z>=0 && currentid != 0){
	MPI_Send(&hilike, 1, MPI_DOUBLE, 0, 8787, MPI_COMM_WORLD);	
	MPI_Send(&hiprob, 1, MPI_DOUBLE, 0, 6767, MPI_COMM_WORLD);
	MPI_Send(&like, 1, MPI_DOUBLE, 0, 5656, MPI_COMM_WORLD);
  }
  if (z < 0 && currentid == 0) {
	MPI_Recv(&hilike, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 8787, MPI_COMM_WORLD, &status);
	MPI_Recv(&hiprob, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 6767, MPI_COMM_WORLD, &status);
	MPI_Recv(&like, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 5656, MPI_COMM_WORLD, &status);
  }

  if (((step / (int) printint) * (int) printint == step && step > 0)
      || outto != stdout)
  {
    if (currentid == 0) {
    printsteps (outto, like);

    callprintacceptancerates (outto);

    printcurrentvals (outto);
    callprintautoctable (outto /*, step*/);
    }
    if (numchains * numprocesses > 1) {
	//AS: debug only
	//printf("printchaininfo called on processor %d\n", currentid);
      printchaininfo (outto, heatmode, hval1, hval2, currentid); //AS: TODO add currentid to printchaininfo 9/26/2014
	}
    if (currentid == 0) {
  	 if (genealogiessaved > 0)
	    {
	      time (&remained_endtime);
	    }
    }
    /* For ASSIGNMENT */
    if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
    {
      printsteps (outto, like);
      imaAsnPrintNumGenesPopn (outto);
      for (ci = 0; ci < numchains; ci++)
      {
        for (j = 0; j < Cupinf[ci].num_uptypes; j++)
        {
          Cupinf[ci].upinf[j].accp = 0;
          Cupinf[ci].upinf[j].tries = 0;
        }
      }
      for (li = 0; li < nloci; li++)
      {
        for (j = 0; j < L[li].a_rec->num_uptypes; j++)
          L[li].a_rec->upinf[j].accp = L[li].a_rec->upinf[j].tries = 0;
/*
        if (assignmentoptions[POPULATIONASSIGNMENTASSIGNED] == 1)
        {
          IMA_assigned_reset ();
        }
*/
      }
    }

  }

//debug
if (currentid == 0) {
if (burndone)
  return_joint_t (multitpeak);
}
  return;
}                               /* intervaloutput */

// check if it is time to call record() and savegenealogyinfo(), and call if it is
void check_to_record (int currentid) //AS: adding currentid 
{
  static int i;
  static int j;
  static int init = 0;

  if (init == 0)
  {
    i = recordint;
    j = savegenealogyint;
    init = 1;
  }
  if (i == recordint)
  {
    record (currentid);                  // record values of parameters that are in mcmc //AS: adding currentid
    recordstep++;
    i = 1;
  }
  else
  {
    i++;
  }

  if (j == savegenealogyint)
  {
    savegenealogyinfo (currentid);            // record values associated with genealogies //AS: adding currentid
    if (jheyoptions[WRITEMIGRATIONNAME])
      record_migration_names();
    if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
    {
      recordassignmentloci (outfilename, 0);
    }

    // VS record the assignment vectors
    saveassignlociinfo(currentid);     //AS: adding currentid

    genealogiessaved++;
    j = 1;
  }
  else
  {
    j++;
  }
}                               // check_to_record

//AS: adding this function for now - I doubt I need it, since I am picking swapper and swappee on the fly
// as of 9/26/2014
void fillswaparrays(int *swapper, int *swappee)
{
	int i; //AS: index
	for (i = 0; i < 1; i++) {
		swapper[i] = rand() % numprocesses;
		swappee[i] = rand() % numprocesses;
	}
	return;
}

int main (int argc, char *argv[])
{
  // VS
  int gp; // VS gp
  char tempgenealogyinfosavefilename[1000]; // VS filename string
  int swapA; //AS: Swapper chain ID
  int swapB; //AS: Swappee chain ID
  int z; //AS: coldchin ID
  int x, y; //AS: indices for swapcount MPI_Reduce operations
  int currentid; //AS: ID of current processor calling this function
//  int swaps; //AS: tracking number of swap tries per MCMC iteration
  //AS: Adding MPI Initializations
  #ifdef MPI_ENABLED
  MPI_Init(&argc, &argv); //AS: Initialize MPI_COMM_WORLD
  MPI_Comm_size(MPI_COMM_WORLD, &numprocesses); //AS: Getting the total number of processors involved
  MPI_Comm_rank(MPI_COMM_WORLD, &currentid); //AS: Getting the currentid - ID of the current processor calling main()
  MPI_Status status;
  //printf("numprocesses = %d, currentid = %d\n",numprocesses,currentid);  
  #endif	
  init_IMA ();
  start (argc, argv, currentid); //AS: adding currentid
  if (runoptions[LOADRUN])
  {
	if (currentid == 0) {
	loadgenealogyvalues ();
    recordstep = genealogiessaved;
    if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
    {
      IMAgenealogyAlign (outfilename, genealogiessaved);
    }
	}
   printoutput (currentid);
   if (currentid == 0) {
   free_ima_main_stuff ();
  }
   return 0;
  }
	

  printf ("Starting Markov chain on processor %d \n", currentid);
  init_after_start_IMA ();
  //printf("Inited after start...\n");
  step = 0;
  recordstep = 0;
  
  while (run (currentid))
  {
	//AS: Mon Feb  1 13:40:45 EST 2016
	//Upon chatting with Jody, realized that there are persistent mixing problems in the current swapping scheme
	//So decided to adhere to his original scheme of picking chains/temperatures that are at some small distance
	//from each other, to avoid wasted swaps across extreme temperatures. This could potentially improve mixing
	//Also, decided to add in swaptries as the number of attempts per MCMC iteration.

	//AS: debug only
	//printf("inside run...\n");
    //AS: pick swapper and swappee chains prior to update in each iteration
    
    if (numprocesses * numchains > 1) {
	    swapA = rand() % (numprocesses * numchains);
	    swapB = rand() % (numprocesses * numchains);
	    while (swapA == swapB) {
		swapB = rand() % (numprocesses * numchains);
    	}
   }
	//AS: debug 
	//printf("before qupdate...%d %d\n",swapA, swapB);	
    qupdate(currentid, swapA, swapB); //AS: adding currentid, swapper and swappee ID's
    if (burndone)
    {
	//AS: debug
	//printf("burn done!\n");
      check_to_record (currentid); // VS commented this as it seems to be the place where the program crash //AS: adding currentid
#ifdef COUNT_PRIOR_ASSIGNMENT
      if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
      {
        IMA_rgf_tickSaasn (0);
      }
#endif // COUNT_PRIOR_ASSIGNMENT
      z = whichiscoldchain();
      if (z >= 0) {
	gloglikelihood += C[z]->allpcalc.pdg; //AS: changing this to be computed from the coldchain only
      	gloglikelihood /= 2.0; //AS: PS - is this used at all though?? 9/26/2014
      }
    }
    step++;
	//AS: debug only
	//printf("step here is %d\n", step);
  }

//AS: all processes need to finish M mode before starting L mode
#ifdef MPI_ENABLED
MPI_Barrier(MPI_COMM_WORLD);
#endif

//AS: debug only
//printf("Done with M mode! yay!\n");

#ifdef COUNT_PRIOR_ASSIGNMENT
  if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
  {
    IMA_rgf_saveSaasn ();
  }
#endif // COUNT_PRIOR_ASSIGNMENT

  // save genealogy info in *.ti file
  if (!runoptions[DONTSAVEGENEALOGIES] && genealogiessaved > 0 && currentid == 0)
  {
        //savegenealogyfile (genealogyinfosavefilename, genealogyinfosavefile, &lastgenealogysaved, gsampinflength);
	// VS - changed the way we deal with lastgenealogy saved since this is used for saving the genealogies and the assignment vectors of loci into groups 
	// lastgenealogysaved = -1; // lastgenealogy saved is initiallized in init_IMA()

	// VS
    // SAVE GENEALOGIES - save the summaries of the trees into files. Each group of loci corresponds to one file.
	for(gp=0; gp<nbgroupsloci_theta+nbgroupsloci_mig; gp++) {
		if(gp<nbgroupsloci_theta) { // if the group is for theta
			sprintf(tempgenealogyinfosavefilename, "%s.groupTheta%i", genealogyinfosavefilename, gp); // create name of file for group Theta
		}
		else {
			sprintf(tempgenealogyinfosavefilename, "%s.groupMig%i", genealogyinfosavefilename, gp); // create name of file for group Mig
		}
		// VS modified savegenealogyfile_vs to get as input the group index
		savegenealogyfile_vs(gp, tempgenealogyinfosavefilename, genealogyinfosavefile, lastgenealogysaved, gsampinflength);
	}
	
	// VS SAVE ASSIGNMENT VECTORS if vectors are being updated
	// get file name for assigment vector theta
	sprintf(tempgenealogyinfosavefilename, "%s.assignlociTheta", genealogyinfosavefilename);
	saveassignlocifile(tempgenealogyinfosavefilename, lastgenealogysaved, assignlocisample_theta);	


	// VS get file name for assigment vector mig
	sprintf(tempgenealogyinfosavefilename, "%s.assignlociMig", genealogyinfosavefilename);
	saveassignlocifile(tempgenealogyinfosavefilename, lastgenealogysaved, assignlocisample_mig);


	// VS Change the value of lastgenealogysaved
	lastgenealogysaved = genealogiessaved - 1;
	
  }

  if (assignmentoptions[POPULATIONASSIGNMENT] == 1
      && assignmentoptions[POPULATIONASSIGNMENTASSIGNED] == 0)
  {
    IMAgenealogyAlign (outfilename, genealogiessaved);
  }
  //AS: Should combine all swap matrices at this point - TODO! 9/26/2014
  //AS: DONE as on 10/22/2014

	/*#ifdef MPI_ENABLED
		//MPI_Barrier(MPI_COMM_WORLD);
		if (numprocesses > 1) {
			for (x = 0; x < numprocesses; x++) {
				for (y = 0; y < numprocesses; y++) {
					MPI_Reduce(&swaps_bwprocesses[x][y], &swaps_rec_bwprocesses[x][y], 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
				}
			}
			if (currentid == 0) {
				for (x = 0; x < numchains; x++) {
					for (y = 0; y < numchains; y++) {
						swapcount_bwprocesses[x][y] = swapcount[x][y];
					}
				}
			}
			for (x = 1; x < numprocesses; x++) {
				if (currentid == x) {
					for (y = 0; y < numchains; y++) {
						for (z = 0; z < numchains; z++) {
							MPI_Send(&swapcount[y][z], 1, MPI_INT, 0, 1234, MPI_COMM_WORLD);
						}
					}
				}
				if (currentid == 0) {
					for (y = 0; y < numchains; y++) {
						for (z = 0; z < numchains; z++) {
							MPI_Recv(&swapcount_bwprocesses[x * numchains + y][x * numchains + z], 
								1, MPI_INT, x, 1234, MPI_COMM_WORLD, &status);
						}
					}
				}
			}
		}
	#endif*/
	
  //AS: printoutput is called only on the head node - doesn't affect this if serial or parallel
  //AS: changed this - has to be called on every node now - so they can collate swap info
  //if (currentid == 0) {
	  printoutput (currentid); // print to the output file
  //}

  if (jheyoptions[WRITEMIGRATIONNAME])
    f_close(migrationnamefile);

  free_ima_main_stuff ();
  #ifdef MPI_ENABLED
  MPI_Finalize();
  #endif
  return 0;
}                               /* main */
