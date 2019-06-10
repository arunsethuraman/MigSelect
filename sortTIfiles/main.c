/* MigSelect 2013-2019 Arun Sethuraman, Vitor Sousa, Jody Hey */
/* This code was developed based on IMa2  2009-2010  Jody Hey, Rasmus Nielsen and Sang Chul Choi*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "DPART.h"
#include "sortgweight.h"

// SORT TI FILES sorts the gweights in TI files according to the assignment vectors
// with minimum partition distance towards the mean partition

// This function only works with 2 groups of loci!

int main (int argc, char *argv[])
{

	// One problem is that there are several TI files, and these must be sorted simultaneously.
	// That would avoid errors due to performing the sorting independently for each file due to 
	// sorting instability. 
	// Given that we can have several ti files with the Theta Ti files and with the Mig
	// the names of the files depend on the number of groups.
	// We need to give as input the number of groups for Theta, and the number of groups for Mig.

	FILE *outputfile; // output file with the meanpartition and coassignment probabilities
	FILE *inputassignfile; // input file with the sampled assignments
	FILE *inputtifile; // input TI files
	int *sampledassign; // matrix with the sampled assigned vectors
	int *tmpassign; // temporary assignment vector for one vector
	int *meanpartition; // array where the meanpartition is saved
	int ngroups; // number of groups
	int ngroupsmig; // number of groups for mig
	int ngroupstheta; // number of groups for theta
	int nloci; // number of loci
	int nsamples; // sample size of assignment vectors
	int i, j, m; 
	char *loadfilebase; // string with the name of the input file for assignment
	char *tifilename; // string with the name of the TI input file
	char tempfilename[300]; // string with temp names of TI input files
	char textline[1000]; // used to skip lines when reading the files
	int *indicator; // indicator variable with size nsamples
				  // indicat[i]=0 if labels do not need to be changed
				  // indicat[i]=1 if labels have to be swapped (indicating that we need to swap the gweights)
	double **gsampinfg1; // matrix with the genealogy weights for group 1
	double **gsampinfg2; // matrix with the genealogy weights for group 2
	double *tmpgweight;  // temporary gweights corresponding to one row of gsampinf (used to relabel if needed)
	int sizeginf; // length of each row of the gsampginf matrix 
        int indpORmaps; // indicator variable that takes value 0 if theta and mig params have independent assignment
                        // and the value 1 if theta and mig params have the same assignment vector

	// READ INPUT ARGS COMMAND
	loadfilebase = argv[1]; // input file name with sampled assignment vectors (this is just a generic name)
	tifilename = argv[2]; // input file name for TI files (this is just the generic name)
	ngroups=atoi(argv[3]); // number of groups
	nloci=atoi(argv[4]); // number of loci
	nsamples=atoi(argv[5]); // number of samples
	ngroupstheta = atoi(argv[6]); // number of groups theta
	ngroupsmig = atoi(argv[7]); // number of groups mig
	sizeginf = atoi(argv[8]); // length of each row of the gsampginf matrix 
        indpORmaps = atoi(argv[9]); // indicator variable that takes value 0 if theta and mig params have independent assignment
                                    // and the value 1 if theta and mig params have the same assignment vector 

	// allocate memory for sampledassign
	// instead of assigning a matrix, it assigns an array
	// this is because of the functions implemented to get the mean partition treat matrices as arrays
	sampledassign = (int *) malloc((nloci*nsamples)*sizeof(int));

	// allocate memory for tmpassign
	tmpassign = (int *) malloc (nloci * sizeof(int));

	// allocate memory for meanpartition mig
	meanpartition = (int *) malloc(nloci*sizeof(int));        

	// allocate memory for indicator variable
	indicator = (int *) malloc(nsamples*sizeof(int));

	// allocate memory for gsampinf matrix
	gsampinfg1 = (double **) malloc(nsamples*sizeof(double *));
	for(i=0; i<nsamples; i++) {
		gsampinfg1[i] = (double *) malloc(sizeginf*sizeof(double));
	}

	// allocate memory for gsampinf matrix
	gsampinfg2 = (double **) malloc(nsamples*sizeof(double *));
	for(i=0; i<nsamples; i++) {
		gsampinfg2[i] = (double *) malloc(sizeginf*sizeof(double));
	}

	// allocate memory for tmpgweight
	tmpgweight = (double *) malloc(sizeginf*sizeof(double));

	// open the file to read the assignments sampled in the MCMC
        // get file name
	sprintf(tempfilename, "%s%s", loadfilebase, ".Unsorted_Mig.out");        
	inputassignfile=fopen(tempfilename, "r");
	if (inputassignfile==NULL) {                
		perror ("Error opening assignment file for Mig.");
                printf("\nError reading assignment file %s\n", inputassignfile);
		exit(1);
	}
	
	// go through the file and save it into array sampledassign
	for(i=0; i<(nloci*nsamples); i++) {
		if(fscanf(inputassignfile, "%i", &sampledassign[i])==EOF) {
			printf("\nError reading assignment file\n");
			exit(1);
		}
	}
	// Close the inputfile
	fclose(inputassignfile);


	// Compute the mean partition
	Meanpartdis(sampledassign, nloci, nsamples, meanpartition);

	// Print the mean partition
	printf("Mean Partition for mig:\n");
	for(i=0; i<nloci; i++) {
		printf("%i ", meanpartition[i]);
	}

	// Sort the sampledassign matrix by sorting the vectors such that they have
	// the minimum partition distance with the mean assignment
	sortsampledassign(indicator, sampledassign, meanpartition, nloci, nsamples);

	// Print the indicator variable
//	printf("\nIndicator variable:\n");
//	for(i=0; i<nsamples; i++) {
//		printf("%i\n", indicator[i]);
//	}


	// Sort assignment vectors
	for(i=0; i<nsamples; i++) {
		// relabel the assignment vector
		if(indicator[i]==1) {
			for(j=0; j<nloci; j++) {
				// Copy old label into tmpassign
				tmpassign[j] = 	sampledassign[i*nloci + j];
				// relabel
				if(tmpassign[j]==0) {sampledassign[i*nloci + j]=1;}
				if(tmpassign[j]==1) {sampledassign[i*nloci + j]=0;}
			}
		}
	}

	// Save the sorted assignment vector
	// open the file to write the assignments sampled in the MCMC
	sprintf(tempfilename, "%s%s", loadfilebase, ".Sorted_Mig.out");
        outputfile=fopen(tempfilename, "w");
	if (outputfile!=NULL) {
		for(i=0; i<nsamples; i++) {
			// relabel the assignment vector
			for(j=0; j<nloci; j++) {
				fprintf(outputfile, "%i\t", sampledassign[i*nloci + j]);
			}
			fprintf(outputfile, "\n");
		}
	}
	else {
		perror ("Error opening file to save sorted assignment");
		exit(1);
	}
	fclose(outputfile);

	// Sort the TI files for the mig groups
	if(ngroupsmig>1)  {
		for(i=0; i<ngroupsmig; i++) {
			// get file name
			sprintf(tempfilename, "%s%s%i", tifilename, ".ti.groupMig", i+ngroupstheta);
			// check file names
			printf("\nReading %s", tempfilename);

			// Open the file
			inputtifile=fopen(tempfilename, "r");
			if (inputtifile==NULL) {
				perror ("Error opening gsampinf file for mig");
				exit(1);
			}

			// skips the first line (starts with VALUESSTART)
			fgets (textline, 1000, inputtifile);

			// read the ti file and save it into gsampinf matrix
			for(j=0; j<nsamples; j++) {
				for(m=0; m<sizeginf; m++) {
					if(i==0) { // if it is the first group
						if(fscanf(inputtifile, "%lf", &gsampinfg1[j][m])==EOF) {
							printf("\nError reading gweight file  for mig g1\n");
							exit(1);
						}
					}
					if(i==1) { // if it is the first group
						if(fscanf(inputtifile, "%lf", &gsampinfg2[j][m])==EOF) {
							printf("\nError reading gweight file  for mig g2\n");
							exit(1);
						}
					}
				}
			}

			// close the file
			fclose(inputtifile);

		}

		// Check the TI file - print first and last row			
		/*printf("\ngsampinf1 before swap\n");
		for(i=0; i<nsamples; i++) {
			for(m=0; m<sizeginf; m++) {
				printf("%g ", gsampinfg1[i][m]);
			}
			printf("\n");
		}
		
		printf("\ngsampinf2 before swap\n");
		for(i=0; i<nsamples; i++) {
			for(m=0; m<sizeginf; m++) {
				printf("%g ", gsampinfg2[i][m]);
			}
			printf("\n");
		}*/

		// go through each row, and swap the gweights where it is needed
		for(i=0; i<nsamples; i++) {

			//printf("\nindicator[%i]=%i", i, indicator[i]);

			// check if the gweights need to be swapped among groups				
			if(indicator[i]==1) {
				// go through the gsampinf line
				for(j=0; j<sizeginf; j++) {
					// copy the gsampinfg1 into tmpgsampinf
					tmpgweight[j] = gsampinfg1[i][j];
					// set gsampingg1 equal to gsampinfg2
					gsampinfg1[i][j] = gsampinfg2[i][j];
					// set gsampingg2 equal to tmpgsampinf
					gsampinfg2[i][j] = tmpgweight[j];

					// CHECK
					//printf("\n tmpweight[%i]=%g", j, tmpgweight[j]);
				}
			}
		
		}


		// Check the TI file - print first and last row			
		/*printf("\ngsampinf1 after swap\n");
		for(i=0; i<nsamples; i++) {
			for(m=0; m<sizeginf; m++) {
				printf("%g ", gsampinfg1[i][m]);
			}
			printf("\n");
		}
		
		printf("\ngsampinf2 after swap\n");
		for(i=0; i<nsamples; i++) {
			for(m=0; m<sizeginf; m++) {
				printf("%g ", gsampinfg2[i][m]);
			}
			printf("\n");
		}*/

		// Save the sorted ti files
		for(i=0; i<ngroupsmig; i++) {
			// get file name
			sprintf(tempfilename, "%s%s%i", tifilename, ".ti.groupMig", i+ngroupstheta);

			// check file names
			printf("\nTrying to open %s", tempfilename);

			// Open the file
			outputfile=fopen(tempfilename, "w");

			if (outputfile!=NULL) { // If able to open the file
				// check file names
				printf("\nWriting %s", tempfilename);
	
				// write VALUESSTART in the first line
				fprintf(outputfile, "VALUESSTART\n");

				// write the ti file 
				for(j=0; j<nsamples; j++) {
					for(m=0; m<sizeginf; m++) {
						if(i==0) { // if it is the first group
							fprintf(outputfile, "%.6f\t", (float) gsampinfg1[j][m]);
						}
						if(i==1) { // if it is the first group
							fprintf(outputfile, "%.6f\t", (float) gsampinfg2[j][m]);
						}
					}
					fprintf(outputfile, "\n");
				}
				// close the output file
				fclose(outputfile);
			}
			else {
				perror ("Error opening gsampinf file to write sorted mig");
				exit(1);
			}	
		}

	}

        // If independent assignment vector for theta and mig, need to recompute mean partition for theta, and re-run the index
        if(indpORmaps==0) {
            // get file name
            sprintf(tempfilename, "%s%s", loadfilebase, ".Unsorted_Theta.out");        
            inputassignfile=fopen(tempfilename, "r");
            if (inputassignfile==NULL) {
                    perror ("Error opening assignment file for Mig.");
                    exit(1);
            }

            // go through the file and save it into array sampledassign
            for(i=0; i<(nloci*nsamples); i++) {
                    if(fscanf(inputassignfile, "%i", &sampledassign[i])==EOF) {
                            printf("\nError reading assignment file\n");
                            exit(1);
                    }
            }
            // Close the inputfile
            fclose(inputassignfile);


            // Compute the mean partition
            Meanpartdis(sampledassign, nloci, nsamples, meanpartition);

            // Print the mean partition
            printf("Mean Partition for theta:\n");
            for(i=0; i<nloci; i++) {
                    printf("%i ", meanpartition[i]);
            }

            // Sort the sampledassign matrix by sorting the vectors such that they have
            // the minimum partition distance with the mean assignment
            sortsampledassign(indicator, sampledassign, meanpartition, nloci, nsamples);        
            
            // Sort assignment vectors
            for(i=0; i<nsamples; i++) {
                    // relabel the assignment vector
                    if(indicator[i]==1) {
                            for(j=0; j<nloci; j++) {
                                    // Copy old label into tmpassign
                                    tmpassign[j] = 	sampledassign[i*nloci + j];
                                    // relabel
                                    if(tmpassign[j]==0) {sampledassign[i*nloci + j]=1;}
                                    if(tmpassign[j]==1) {sampledassign[i*nloci + j]=0;}
                            }
                    }
            }

            // Save the sorted assignment vector
            // open the file to write the assignments sampled in the MCMC
            // get file name
            sprintf(tempfilename, "%s%s", loadfilebase, ".Sorted_Theta.out");
            outputfile=fopen(tempfilename, "w");
            if (outputfile!=NULL) {
                    for(i=0; i<nsamples; i++) {
                            // relabel the assignment vector
                            for(j=0; j<nloci; j++) {
                                    fprintf(outputfile, "%i\t", sampledassign[i*nloci + j]);
                            }
                            fprintf(outputfile, "\n");
                    }
            }
            else {
                    perror ("Error opening file to save sorted assignment");
                    exit(1);
            }
            fclose(outputfile);            
            
            
        }

	
	// Sort the TI files for the theta groups
	if(ngroupstheta>1)  {
		for(i=0; i<ngroupstheta; i++) {
			// get file name
			sprintf(tempfilename, "%s%s%i", tifilename, ".ti.groupTheta", i);
			// check file names
			printf("\nReading %s", tempfilename);

			// Open the file
			inputtifile=fopen(tempfilename, "r");
			if (inputtifile==NULL) {
				perror ("Error opening gsampinf file");
				exit(1);
			}

			// skips the first line (starts with VALUESSTART)
			fgets (textline, 1000, inputtifile);

			// read the ti file and save it into gsampinf matrix
			for(j=0; j<nsamples; j++) {
				for(m=0; m<sizeginf; m++) {
					if(i==0) { // if it is the first group
						if(fscanf(inputtifile, "%lf", &gsampinfg1[j][m])==EOF) {
							printf("\nError reading gweight file g1\n");
							exit(1);
						}
					}
					if(i==1) { // if it is the first group
						if(fscanf(inputtifile, "%lf", &gsampinfg2[j][m])==EOF) {
							printf("\nError reading gweight file g2\n");
							exit(1);
						}
					}
				}
			}

			// close the file
			fclose(inputtifile);

		}

		// Check the TI file - print first and last row			
		/*printf("\ngsampinf1 before swap\n");
		for(i=0; i<nsamples; i++) {
			for(m=0; m<sizeginf; m++) {
				printf("%g ", gsampinfg1[i][m]);
			}
			printf("\n");
		}
		
		printf("\ngsampinf2 before swap\n");
		for(i=0; i<nsamples; i++) {
			for(m=0; m<sizeginf; m++) {
				printf("%g ", gsampinfg2[i][m]);
			}
			printf("\n");
		}*/

		// go through each row, and swap the gweights where it is needed
		for(i=0; i<nsamples; i++) {

			//printf("\nindicator[%i]=%i", i, indicator[i]);

			// check if the gweights need to be swapped among groups				
			if(indicator[i]==1) {
				// go through the gsampinf line
				for(j=0; j<sizeginf; j++) {
					// copy the gsampinfg1 into tmpgsampinf
					tmpgweight[j] = gsampinfg1[i][j];
					// set gsampingg1 equal to gsampinfg2
					gsampinfg1[i][j] = gsampinfg2[i][j];
					// set gsampingg2 equal to tmpgsampinf
					gsampinfg2[i][j] = tmpgweight[j];

					// CHECK
					//printf("\n tmpweight[%i]=%g", j, tmpgweight[j]);
				}
			}
		
		}


		// Check the TI file - print first and last row			
		/*printf("\ngsampinf1 after swap\n");
		for(i=0; i<nsamples; i++) {
			for(m=0; m<sizeginf; m++) {
				printf("%g ", gsampinfg1[i][m]);
			}
			printf("\n");
		}
		
		printf("\ngsampinf2 after swap\n");
		for(i=0; i<nsamples; i++) {
			for(m=0; m<sizeginf; m++) {
				printf("%g ", gsampinfg2[i][m]);
			}
			printf("\n");
		}*/

		// Save the sorted ti files
		for(i=0; i<ngroupstheta; i++) {
			// get file name
			sprintf(tempfilename, "%s%s%i", tifilename, ".ti.groupTheta", i);

			// check file names
			printf("\nTrying to open %s", tempfilename);

			// Open the file
			outputfile=fopen(tempfilename, "w");

			if (outputfile!=NULL) { // If able to open the file
				// check file names
				printf("\nWriting %s", tempfilename);
	
				// write VALUESSTART in the first line
				fprintf(outputfile, "VALUESSTART\n");

				// write the ti file 
				for(j=0; j<nsamples; j++) {
					for(m=0; m<sizeginf; m++) {
						if(i==0) { // if it is the first group
							fprintf(outputfile, "%.6f\t", (float) gsampinfg1[j][m]);
						}
						if(i==1) { // if it is the first group
							fprintf(outputfile, "%.6f\t", (float) gsampinfg2[j][m]);
						}
					}
					fprintf(outputfile, "\n");
				}
				// close the output file
				fclose(outputfile);
			}
			else {
				perror ("Error opening gsampinf file to write sorted theta");
				exit(1);
			}	
		}

	}


	// FREE POINTERS
	free(sampledassign);
	free(meanpartition);
	free(indicator);
	for(i=0; i<nsamples; i++) {
		XFREE(gsampinfg1[i]);
		XFREE(gsampinfg2[i]);
	}
	free(gsampinfg1);
	free(gsampinfg2);
	free(tmpgweight);

}
