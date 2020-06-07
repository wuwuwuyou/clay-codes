/* *
 * Copyright (c) 2014, James S. Plank and Kevin Greenan
 * All rights reserved.
 *
 * Jerasure - A C/C++ Library for a Variety of Reed-Solomon and RAID-6 Erasure
 * Coding Techniques
 *
 * Revision 2.0: Galois Field backend now links to GF-Complete
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 *  - Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 *  - Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 *  - Neither the name of the University of Tennessee nor the names of its
 *    contributors may be used to endorse or promote products derived
 *    from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
 * OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 * AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY
 * WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

/* Jerasure's authors:

   Revision 2.x - 2014: James S. Plank and Kevin M. Greenan.
   Revision 1.2 - 2008: James S. Plank, Scott Simmerman and Catherine D. Schuman.
   Revision 1.0 - 2007: James S. Plank.
 *///
   
/* 

This program takes as input an inputfile, k, m, a coding 
technique, w, and packetsize.  It creates k+m files from 
the original file so that k of these files are parts of 
the original file and m of the files are encoded based on 
the given coding technique. The format of the created files 
is the file name with "_k#" or "_m#" and then the extension.  
(For example, inputfile test.txt would yield file "test_k1.txt".)
*/

#include <assert.h>
#include <time.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <signal.h>
#include <gf_rand.h>
#include <unistd.h>
#include "jerasure.h"
#include "reed_sol.h"
#include "cauchy.h"
#include "liberation.h"
#include "timing.h"

#define N 10
#define M 128
#define r 2
enum Coding_Technique {Reed_Sol_Van, Reed_Sol_R6_Op, Cauchy_Orig, Cauchy_Good, Liberation, Blaum_Roth, Liber8tion, RDP, EVENODD, No_Coding};

char *Methods[N] = {"reed_sol_van", "reed_sol_r6_op", "cauchy_orig", "cauchy_good", "liberation", "blaum_roth", "liber8tion", "no_coding"};

/* Global variables for signal handler */
int readins, n;
enum Coding_Technique method;

/* Function prototypes */
int is_prime(int w);
void ctrl_bs_handler(int dummy);

int jfread(void *ptr, int size, int nmembers, FILE *stream)
{
  if (stream != NULL) return fread(ptr, size, nmembers, stream);

  MOA_Fill_Random_Region(ptr, size);
  return size;
}

static void print_data_and_coding(int k, int m, int w, int size,
	char **data, char **coding)
{
	int i, j, x;
	int n, sp;
//	long l;

	if (k > m) n = k;
	else n = m;
	sp = size * 2 + size / (w / 8) + 8;

	printf("%-*sCoding\n", sp, "Data");
	for (i = 0; i < n; i++) {
		if (i < k) {
			printf("D%-2d:", i);
			for (j = 0; j < size; j += (w / 8)) {
				printf(" ");
				for (x = 0; x < w / 8; x++) {
					printf("%02x", (unsigned char)data[i][j + x]);
				}
			}
			printf("    ");
		}
		else printf("%*s", sp, "");
		if (i < m) {
			printf("C%-2d:", i);
			for (j = 0; j < size; j += (w / 8)) {
				printf(" ");
				for (x = 0; x < w / 8; x++) {
					printf("%02x", (unsigned char)coding[i][j + x]);
				}
			}
		}
		printf("\n");
	}
	printf("\n");
}

int main (int argc, char **argv) {
	FILE *fp, *fp2;				// file pointers
	char *block;				// padding file
	int size, newsize;			// size of file and temp size 
	struct stat status;			// finding file size

	
	enum Coding_Technique tech;		// coding technique (parameter)
	int k, m, w, packetsize;		// parameters
	int buffersize;					// paramter
	int i,j,i1,j1,i2,j2;
	int i22,iii;						// loop control variables
	int blocksize;					// size of k+m files
	int total;
	int extra3;
	int stripe_size;
	
	/* Jerasure Arguments */
	char **data;				
	char **coding;
	char **fdata;				
	char **fcoding;
	//char **ccoding;
	int *matrix;
	int *bitmatrix;
	int **schedule;
	//char *datacopy1;
	//char *datacopy2;
        char *extra1;
	char *extra2;
       // char *datacopy;
	char **extra;

        char **A;					
        char *A1;
	char *A2;
	char **A3;
	char **AA;
	char **AAA;

	
	/* Creation of file name variables */
	char temp[5];
	char *s1, *s2, *extension;
	char *fname;
	int md;
	char *curdir;
	
	/* Timing variables */
	struct timing t1, t2, t3, t4,q1,q2,q3,q4,q5,q6;

	double tsec;
	double totalsec;
	double encode_time;
	double bit_operation_time;
	double sum_time;
	struct timing start;

	/* Find buffersize */
	int up, down;


	signal(SIGQUIT, ctrl_bs_handler);

	/* Start timing */
	timing_set(&t1);
	totalsec = 0.0;
	matrix = NULL;
	bitmatrix = NULL;
	schedule = NULL;
	
	/* Error check Arguments*/
	if (argc != 8) {
		fprintf(stderr,  "usage: inputfile k m coding_technique w packetsize buffersize\n");
		fprintf(stderr,  "\nChoose one of the following coding techniques: \nreed_sol_van, \nreed_sol_r6_op, \ncauchy_orig, \ncauchy_good, \nliberation, \nblaum_roth, \nliber8tion");
		fprintf(stderr,  "\n\nPacketsize is ignored for the reed_sol's");
		fprintf(stderr,  "\nBuffersize of 0 means the buffersize is chosen automatically.\n");
		fprintf(stderr,  "\nIf you just want to test speed, use an inputfile of \"-number\" where number is the size of the fake file you want to test.\n\n");
		exit(0);
	}
	/* Conversion of parameters and error checking */	
	if (sscanf(argv[2], "%d", &k) == 0 || k <= 0) {
		fprintf(stderr,  "Invalid value for k\n");
		exit(0);
	}
	if (sscanf(argv[3], "%d", &m) == 0 || m < 0) {
		fprintf(stderr,  "Invalid value for m\n");
		exit(0);
	}
	if (sscanf(argv[5],"%d", &w) == 0 || w <= 0) {
		fprintf(stderr,  "Invalid value for w.\n");
		exit(0);
	}
	if (argc == 6) {
		packetsize = 0;
	}
	else {
		if (sscanf(argv[6], "%d", &packetsize) == 0 || packetsize < 0) {
			fprintf(stderr,  "Invalid value for packetsize.\n");
			exit(0);
		}
	}
	if (argc != 8) {
		buffersize = 0;
	}
	else {
		if (sscanf(argv[7], "%d", &buffersize) == 0 || buffersize < 0) {
			fprintf(stderr, "Invalid value for buffersize\n");
			exit(0);
		}
		
	}

	/* Determine proper buffersize by finding the closest valid buffersize to the input value  */
	if (buffersize != 0) {
		if (packetsize != 0 && buffersize%(sizeof(long)*w*k*packetsize) != 0) { 
			up = buffersize;
			down = buffersize;
			while (up%(sizeof(long)*w*k*packetsize) != 0 && (down%(sizeof(long)*w*k*packetsize) != 0)) {
				up++;
				if (down == 0) {
					down--;
				}
			}
			if (up%(sizeof(long)*w*k*packetsize) == 0) {
				buffersize = up;
			}
			else {
				if (down != 0) {
					buffersize = down;
				}
			}
		}
		else if (packetsize == 0 && buffersize%(sizeof(long)*w*k) != 0) {
			up = buffersize;
			down = buffersize;
			while (up%(sizeof(long)*w*k) != 0 && down%(sizeof(long)*w*k) != 0) {
				up++;
				down--;
			}
			if (up%(sizeof(long)*w*k) == 0) {
				buffersize = up;
			}
			else {
				buffersize = down;
			}
		}
	}

	/* Setting of coding technique and error checking */
	
	if (strcmp(argv[4], "no_coding") == 0) {
		tech = No_Coding;
	}
	else if (strcmp(argv[4], "reed_sol_van") == 0) {
		tech = Reed_Sol_Van;
		if (w != 8 && w != 16 && w != 32) {
			fprintf(stderr,  "w must be one of {8, 16, 32}\n");
			exit(0);
		}
	}
	else if (strcmp(argv[4], "reed_sol_r6_op") == 0) {
		if (m != 2) {
			fprintf(stderr,  "m must be equal to 2\n");
			exit(0);
		}
		if (w != 8 && w != 16 && w != 32) {
			fprintf(stderr,  "w must be one of {8, 16, 32}\n");
			exit(0);
		}
		tech = Reed_Sol_R6_Op;
	}
	else if (strcmp(argv[4], "cauchy_orig") == 0) {
		tech = Cauchy_Orig;
		if (packetsize == 0) {
			fprintf(stderr, "Must include packetsize.\n");
			exit(0);
		}
	}
	else if (strcmp(argv[4], "cauchy_good") == 0) {
		tech = Cauchy_Good;
		if (packetsize == 0) {
			fprintf(stderr, "Must include packetsize.\n");
			exit(0);
		}
	}
	else if (strcmp(argv[4], "liberation") == 0) {
		if (k > w) {
			fprintf(stderr,  "k must be less than or equal to w\n");
			exit(0);
		}
		if (w <= 2 || !(w%2) || !is_prime(w)) {
			fprintf(stderr,  "w must be greater than two and w must be prime\n");
			exit(0);
		}
		if (packetsize == 0) {
			fprintf(stderr, "Must include packetsize.\n");
			exit(0);
		}
		if ((packetsize%(sizeof(long))) != 0) {
			fprintf(stderr,  "packetsize must be a multiple of sizeof(long)\n");
			exit(0);
		}
		tech = Liberation;
	}
	else if (strcmp(argv[4], "blaum_roth") == 0) {
		if (k > w) {
			fprintf(stderr,  "k must be less than or equal to w\n");
			exit(0);
		}
		if (w <= 2 || !((w+1)%2) || !is_prime(w+1)) {
			fprintf(stderr,  "w must be greater than two and w+1 must be prime\n");
			exit(0);
		}
		if (packetsize == 0) {
			fprintf(stderr, "Must include packetsize.\n");
			exit(0);
		}
		if ((packetsize%(sizeof(long))) != 0) {
			fprintf(stderr,  "packetsize must be a multiple of sizeof(long)\n");
			exit(0);
		}
		tech = Blaum_Roth;
	}
	else if (strcmp(argv[4], "liber8tion") == 0) {
		if (packetsize == 0) {
			fprintf(stderr, "Must include packetsize\n");
			exit(0);
		}
		if (w != 8) {
			fprintf(stderr, "w must equal 8\n");
			exit(0);
		}
		if (m != 2) {
			fprintf(stderr, "m must equal 2\n");
			exit(0);
		}
		if (k > w) {
			fprintf(stderr, "k must be less than or equal to w\n");
			exit(0);
		}
		tech = Liber8tion;
	}
	else {
		fprintf(stderr,  "Not a valid coding technique. Choose one of the following: reed_sol_van, reed_sol_r6_op, cauchy_orig, cauchy_good, liberation, blaum_roth, liber8tion, no_coding\n");
		exit(0);
	}

	/* Set global variable method for signal handler */
	method = tech;

	/* Get current working directory for construction of file names */
	curdir = (char*)malloc(sizeof(char)*1000);	
	assert(curdir == getcwd(curdir, 1000));

        if (argv[1][0] != '-') {

		/* Open file and error check */
		fp = fopen(argv[1], "rb");
		if (fp == NULL) {
			fprintf(stderr,  "Unable to open file.\n");
			exit(0);
		}
	
		/* Create Coding directory */
		i = mkdir("Coding", S_IRWXU);
		if (i == -1 && errno != EEXIST) {
			fprintf(stderr, "Unable to create Coding directory.\n");
			exit(0);
		}
	
		/* Determine original size of file */
		stat(argv[1], &status);	
		size = status.st_size;
        } else {
        	if (sscanf(argv[1]+1, "%d", &size) != 1 || size <= 0) {
                	fprintf(stderr, "Files starting with '-' should be sizes for randomly created input\n");
			exit(1);
		}
        	fp = NULL;
		MOA_Seed(time(0));
        }

	newsize = size;
	
	/* Find new size by determining next closest multiple */
	if (packetsize != 0) {
		if (size%(k*w*packetsize*sizeof(long)) != 0) {
			while (newsize%(k*w*packetsize*sizeof(long)) != 0) 
				newsize++;
		}
	}
	else {
		if (size%(k*w*sizeof(long)) != 0) {
			while (newsize%(k*w*sizeof(long)) != 0) 
				newsize++;
		}
	}
	
	if (buffersize != 0) {
		while (newsize%buffersize != 0) {
			newsize++;
		}
	}


	/*if ((size % M) != 0) {
	while ((newsize % M) != 0)
		newsize++;
	}
	stripe_size = newsize / M;
	temp_stripe_size = stripe_size;
	if ((stripe_size % k) != 0){
		while ((temp_stripe_size % k) != 0)
			temp_stripe_size++;
       }
	newsize = temp_stripe_size * M;
        stripe_size=newsize/M;
	block_size= temp_stripe_size /k;*/



	/* Determine size of k+m files */
	stripe_size=newsize/M;
	blocksize = stripe_size/k;
	
	
	/* Allow for buffersize and determine number of read-ins */
	if (size > buffersize && buffersize != 0) {
		if (newsize%buffersize != 0) {
			readins = newsize/buffersize;
		}
		else {
			readins = newsize/buffersize;
		}
		block = (char *)malloc(sizeof(char)*buffersize);
		blocksize = buffersize/k;
	}
	else {
		readins = 1;
		buffersize = size;
		block = (char *)malloc(sizeof(char)*newsize);
	}
	printf("buffersize:%d\n", buffersize);
	printf("size:%d\n", size);
        printf("newsize:%d\n",newsize);
	printf("stripe_size:%d\n",stripe_size);	
	printf("blocksize:%d\n", blocksize);
	/* Break inputfile name into the filename and extension */	
	s1 = (char*)malloc(sizeof(char)*(strlen(argv[1])+20));
	s2 = strrchr(argv[1], '/');
	if (s2 != NULL) {
		s2++;
		strcpy(s1, s2);
	}
	else {
		strcpy(s1, argv[1]);
	}
	s2 = strchr(s1, '.');
	if (s2 != NULL) {
          extension = strdup(s2);
          *s2 = '\0';
	} else {
          extension = strdup("");
        }
	
	/* Allocate for full file name */
	fname = (char*)malloc(sizeof(char)*(strlen(argv[1])+strlen(curdir)+20));
	sprintf(temp, "%d", k);
	md = strlen(temp);
	
	/* Allocate data and coding */
	data = (char **)malloc(sizeof(char*)*k);
	coding = (char **)malloc(sizeof(char*)*m);
	for (i = 0; i < m; i++) {
		coding[i] = (char *)malloc(sizeof(char)*blocksize);
                if (coding[i] == NULL) { perror("malloc"); exit(1); }
	}

	fdata = (char **)malloc(sizeof(char*)*M);
	fcoding = (char **)malloc(sizeof(char*)*M);
	//ccoding = (char **)malloc(sizeof(char*)*M*m);
	//datacopy1 =  (char *)malloc(sizeof(char)*blocksize);
	//datacopy2 =  (char *)malloc(sizeof(char)*blocksize);
        extra1 = (char *)malloc(sizeof(char)*blocksize);
	extra2 = (char *)malloc(sizeof(char)*blocksize);
        extra = (char **)malloc(sizeof(char*)*128);
		for(i=0;i<128;i++)
		{extra[i] = (char *)malloc(sizeof(char)*blocksize);}
	//printf("coding[0]:%p\n",&coding[0]);
       
	/*A1 =  (char *)malloc(sizeof(char));
	A2 =  (char *)malloc(sizeof(char));
	A3 = (char**)malloc(sizeof(char*)*2);
     for (i = 0; i < 2; ++i){
        A3[i] = (char*)malloc(sizeof(char)*2);}
 
       A = (char**)malloc(sizeof(char*)*2);
     for (i = 0; i < 2; ++i){
        A[i] = (char*)malloc(sizeof(char)*4);}
    AA = (char **)malloc(sizeof(char*)*2);
	for(j = 0; j < 2; j++)  {
	AA[j] = (char *)malloc(sizeof(char)*6);}
		AAA = (char **)malloc(sizeof(char*)*6);
		for(j = 0; j < 2; j++)  {
		AAA[j] = (char *)malloc(sizeof(char)*2);}*/
  

	/* Create coding matrix or bitmatrix and schedule */
	timing_set(&t3);
	switch(tech) {
		case No_Coding:
			break;
		case Reed_Sol_Van:
			matrix = reed_sol_vandermonde_coding_matrix(k, m, w);
			break;
		case Reed_Sol_R6_Op:
			break;
		case Cauchy_Orig:
			matrix = cauchy_original_coding_matrix(k, m, w);
			bitmatrix = jerasure_matrix_to_bitmatrix(k, m, w, matrix);
			schedule = jerasure_smart_bitmatrix_to_schedule(k, m, w, bitmatrix);
			break;
		case Cauchy_Good:
			matrix = cauchy_good_general_coding_matrix(k, m, w);
			bitmatrix = jerasure_matrix_to_bitmatrix(k, m, w, matrix);
			schedule = jerasure_smart_bitmatrix_to_schedule(k, m, w, bitmatrix);
			break;	
		case Liberation:
			bitmatrix = liberation_coding_bitmatrix(k, w);
			schedule = jerasure_smart_bitmatrix_to_schedule(k, m, w, bitmatrix);
			break;
		case Blaum_Roth:
			bitmatrix = blaum_roth_coding_bitmatrix(k, w);
			schedule = jerasure_smart_bitmatrix_to_schedule(k, m, w, bitmatrix);
			break;
		case Liber8tion:
			bitmatrix = liber8tion_coding_bitmatrix(k);
			schedule = jerasure_smart_bitmatrix_to_schedule(k, m, w, bitmatrix);
			break;
		case RDP:
		case EVENODD:
			assert(0);
	}
	timing_set(&start);
	timing_set(&t4);
	totalsec += timing_delta(&t3, &t4);
 printf("matrix: \n");
	jerasure_print_matrix(matrix,k+m,k,w);
printf("\n");

	/* Read in data until finished */
	n = 1;
	total = 0;

	while (n <= readins) {
		/* Check if padding is needed, if so, add appropriate 
		   number of zeros */
		if (total < size && total+buffersize <= size) {
			total += jfread(block, sizeof(char), buffersize, fp);
		}
		else if (total < size && total+buffersize > size) {
			extra3 = jfread(block, sizeof(char), buffersize, fp);
			for (i = extra3; i < buffersize; i++) {
				block[i] = '0';
			}
		}
		else if (total == size) {
			for (i = 0; i < buffersize; i++) {
				block[i] = '0';
			}
		}

printf("clay-encoding: \n");
timing_set(&t3);
timing_set(&q5);
timing_set(&q1);
		for(j = 0; j < M; j++)
	     {
		fdata[j] = (char *)malloc(sizeof(char)*k*blocksize);
		fcoding[j] = (char *)malloc(sizeof(char)*m*blocksize);
		/* Set pointers to point to file data */
		for (i = 0; i < k; i++) {
		   data[i] = block+((j*k+i)*blocksize);}
		

              
		
		/* Encode according to coding method */
		switch(tech) {	
			case No_Coding:
				break;
			case Reed_Sol_Van:
				jerasure_matrix_encode(k, m, w, matrix, data, coding, blocksize);
				break;
			case Reed_Sol_R6_Op:
				reed_sol_r6_encode(k, w, data, coding, blocksize);
				break;
			case Cauchy_Orig:
				jerasure_schedule_encode(k, m, w, schedule, data, coding, blocksize, packetsize);
				break;
			case Cauchy_Good:
				jerasure_schedule_encode(k, m, w, schedule, data, coding, blocksize, packetsize);
				break;
			case Liberation:
				jerasure_schedule_encode(k, m, w, schedule, data, coding, blocksize, packetsize);
				break;
			case Blaum_Roth:
				jerasure_schedule_encode(k, m, w, schedule, data, coding, blocksize, packetsize);
				break;
			case Liber8tion:
				jerasure_schedule_encode(k, m, w, schedule, data, coding, blocksize, packetsize);
				break;
			case RDP:
			case EVENODD:
				assert(0);
		}
 	        
 
			/* for(i1=0;i1<1;i1++){
	                        for(j1=0;j1<blocksize;j1++)
                                   {
					printf("%d ",fcoding[i1][j1]);
			           }
				printf( " \n");
				}*/
	    
                    //printf("Encoding Complete:\n\n");
	           // print_data_and_coding(k, m, w, sizeof(long), data, coding);
				
				
				//for (i=0;i<m*blocksize;i++)
				//{fcoding[j][i]= *(coding[0]+i);}
				

		/*for(i=0;i<k;i++){
                for(i1=0;i1<blocksize;i1++){
                   fdata[j][i*blocksize+i1]=data[i][i1];}
                }

               for(i=0;i<m;i++){
               for(i1=0;i1<blocksize;i1++){
                   fcoding[j][i*blocksize+i1]=coding[i][i1];}
               }	*/



				
				for(iii=0;iii<m;iii++){
				for (i=0;i<blocksize;i++){
				fcoding[j][i22]= *(coding[iii]+i);
				i22++;}}
				i22=0;
				
				 for(i=0;i<k*blocksize;i++)
				  {
				  fdata[j][i] = *(block+j*k*blocksize+i);}

	   
		
	   }		
timing_set(&q2);

printf("encoder complete \n");		
	/*printf( " original fcoding first row m*blocksize:\n");
			 for(i1=0;i1<1;i1++){
	                        for(j1=0;j1<m*blocksize;j1++)
                                   {
					printf("%d ",fcoding[0][j1]);
			           }
				printf( " \n");}
	printf( " original fdata--------  after encoder:\n");
			
			 for(i1=0;i1<1;i1++){
	                        for(j1=0;j1<k*blocksize;j1++)
                                   {
					printf("%d ",fdata[0][j1]);
			           }
				printf( " \n");}*/
				 
	 /*   printf( " fcoding   strip after encoder--------------------  :\n");
			
	
			 for(i1=0;i1<M;i1++){
	                        for(j1=0;j1<m*blocksize;j1++){
                                  
				printf("%d ",fcoding[i1][j1]);}
     				printf( " \n");	}
			   printf( " \n");		              
printf( " original fdata n=1 :\n");
			 for(i1=0;i1<1;i1++){
	                        for(j1=0;j1<blocksize;j1++)
                                   {
					printf("%d ",fdata[i1][j1]);
			           }
				printf( " \n");
				}
			printf( " \n");
printf( " original fdata second row n=1 :\n");
			 for(i1=1;i1<2;i1++){
	                        for(j1=0;j1<blocksize;j1++)
                                   {
					printf("%d ",fdata[i1][j1]);
			           }
				printf( " \n");
				}
			printf( " \n");
printf( " original fdata 3 row n=1 :\n");
			 for(i1=2;i1<3;i1++){
	                        for(j1=0;j1<blocksize;j1++)
                                   {
					printf("%d ",fdata[i1][j1]);
			           }
				printf( " \n");
				}
			printf( " \n");
printf( " original fdata 4 row n=1 :\n");
			 for(i1=3;i1<4;i1++){
	                        for(j1=0;j1<blocksize;j1++)
                                   {
					printf("%d ",fdata[i1][j1]);
			           }
				printf( " \n");
				}
			printf( " \n");
printf( " original fdata last two row n=1 :\n");
			 for(i1=M-2;i1<M-1;i1++){
	                        for(j1=0;j1<blocksize;j1++)
                                   {
					printf("%d ",fdata[i1][j1]);
			           }
				printf( " \n");
				}
			printf( " \n");
printf( " original fdata last row n=1 :\n");
			 for(i1=M-1;i1<M;i1++){
	                        for(j1=0;j1<blocksize;j1++)
                                   {
					printf("%d ",fdata[i1][j1]);
			           }
				printf( " \n");
				}
			printf( " \n");
 printf( " original fdata n=2 :\n");
			 for(i1=0;i1<1;i1++){
	                        for(j1=0;j1<blocksize;j1++)
                                   {
					printf("%d ",fdata[i1][j1+1*blocksize]);
			           }
				printf( " \n");
				}
			printf( " \n");
printf( " original fdata 2 row n=2 :\n");
			 for(i1=1;i1<2;i1++){
	                        for(j1=0;j1<blocksize;j1++)
                                   {
					printf("%d ",fdata[i1][j1+blocksize]);
			           }
				printf( " \n");
				}
			printf( " \n");
printf( " original fdata 3 row n=2 :\n");
			 for(i1=2;i1<3;i1++){
	                        for(j1=0;j1<blocksize;j1++)
                                   {
					printf("%d ",fdata[i1][j1+blocksize]);
			           }
				printf( " \n");
				}
			printf( " \n");

printf( " original fdata 4 row n=2 :\n");
			 for(i1=3;i1<4;i1++){
	                        for(j1=0;j1<blocksize;j1++)
                                   {
					printf("%d ",fdata[i1][j1+blocksize]);
			           }
				printf( " \n");
				}
			printf( " \n");

printf( " original fdata last two row n=2 :\n");
			 for(i1=M-2;i1<M-1;i1++){
	                        for(j1=0;j1<blocksize;j1++)
                                   {
					printf("%d ",fdata[i1][j1+blocksize]);
			           }
				printf( " \n");
				}
			printf( " \n");

printf( " original fdata last row n=2 :\n");
			 for(i1=M-1;i1<M;i1++){
	                        for(j1=0;j1<blocksize;j1++)
                                   {
					printf("%d ",fdata[i1][j1+blocksize]);
			           }
				printf( " \n");
				}
			printf( " \n");










 printf( " original fdata n=3 :\n");
			 for(i1=0;i1<1;i1++){
	                        for(j1=0;j1<blocksize;j1++)
                                   {
					printf("%d ",fdata[i1][j1+2*blocksize]);
			           }
				printf( " \n");
				}
			printf( " \n");

 printf( " original 2 row fdata n=3 :\n");
			 for(i1=1;i1<2;i1++){
	                        for(j1=0;j1<blocksize;j1++)
                                   {
					printf("%d ",fdata[i1][j1+2*blocksize]);
			           }
				printf( " \n");
				}
			printf( " \n");
 printf( " original 3 row fdata n=3 :\n");
			 for(i1=2;i1<3;i1++){
	                        for(j1=0;j1<blocksize;j1++)
                                   {
					printf("%d ",fdata[i1][j1+2*blocksize]);
			           }
				printf( " \n");
				}
			printf( " \n");
 printf( " original fdata n=4 :\n");
			 for(i1=0;i1<1;i1++){
	                        for(j1=0;j1<blocksize;j1++)
                                   {
					printf("%d ",fdata[i1][j1+3*blocksize]);
			           }
				printf( " \n");
				}
			printf( " \n");
 printf( " original 2 row  fdata n=4 :\n");
			 for(i1=1;i1<2;i1++){
	                        for(j1=0;j1<blocksize;j1++)
                                   {
					printf("%d ",fdata[i1][j1+3*blocksize]);
			           }
				printf( " \n");
				}
			printf( " \n");
 printf( " original 3 row  fdata n=4 :\n");
			 for(i1=2;i1<3;i1++){
	                        for(j1=0;j1<blocksize;j1++)
                                   {
					printf("%d ",fdata[i1][j1+3*blocksize]);
			           }
				printf( " \n");
				}
			printf( " \n");
 printf( " original fdata n=5 :\n");
			 for(i1=0;i1<1;i1++){
	                        for(j1=0;j1<blocksize;j1++)
                                   {
					printf("%d ",fdata[i1][j1+4*blocksize]);
			           }
				printf( " \n");
				}
			printf( " \n");
 printf( " original fdata n=6 :\n");
			 for(i1=0;i1<1;i1++){
	                        for(j1=0;j1<blocksize;j1++)
                                   {
					printf("%d ",fdata[i1][j1+5*blocksize]);
			           }
				printf( " \n");
				}
			printf( " \n");
 printf( " original fdata n=7 :\n");
			 for(i1=0;i1<1;i1++){
	                        for(j1=0;j1<blocksize;j1++)
                                   {
					printf("%d ",fdata[i1][j1+6*blocksize]);
			           }
				printf( " \n");
				}
				printf( " \n"); 
printf( " original fdata n=8 :\n");
			 for(i1=0;i1<1;i1++){
	                        for(j1=0;j1<blocksize;j1++)
                                   {
					printf("%d ",fdata[i1][j1+7*blocksize]);
			           }
				printf( " \n");
				}
			printf( " \n");
 printf( " original fdata n=9 :\n");
			 for(i1=0;i1<1;i1++){
	                        for(j1=0;j1<blocksize;j1++)
                                   {
					printf("%d ",fdata[i1][j1+8*blocksize]);
			           }
				printf( " \n");
				}
			printf( " \n");
 printf( " original fdata n=10 :\n");
			 for(i1=0;i1<1;i1++){
	                        for(j1=0;j1<blocksize;j1++)
                                   {
					printf("%d ",fdata[i1][j1+9*blocksize]);
			           }
				printf( " \n");
				}
			printf( " \n");
                   printf( " read fdata before:\n");
			 for(i1=0;i1<1;i1++){
	                        for(j1=0;j1<k*blocksize;j1++)
                                   {
					printf("%d ",fdata[i1][j1]);
			           }
				printf( " \n");
				}
			printf( " encoder fcoding  n=11 :\n");
			 for(i1=0;i1<1;i1++){
	                        for(j1=0;j1<blocksize;j1++)
                                   {
					printf("%d ",fcoding[i1][j1]);
			           }
				printf( " \n");
				}
				printf( " \n");
		       printf( " befor decoded coding  n=12 :\n");
			 for(i1=0;i1<1;i1++){
	                        for(j1=0;j1<blocksize;j1++)
                                   {
					printf("%d ",fcoding[i1][j1+blocksize]);
			           }
				printf( " \n");
				}
				printf( " \n");
			printf( " befor decoded coding  n=13 :\n");
			 for(i1=0;i1<1;i1++){
	                        for(j1=0;j1<blocksize;j1++)
                                   {
					printf("%d ",fcoding[i1][j1+2*blocksize]);
			           }
				printf( " \n");
				}
				printf( " \n");
			printf( " encoder fcoding  n=14:\n");
			 for(i1=0;i1<1;i1++){
	                        for(j1=0;j1<blocksize;j1++)
                                   {
					printf("%d ",fcoding[i1][j1+3*blocksize]);
			           }
				printf( " \n");
				}
				printf( " \n");
                 printf( " encoder fdata  n=2 :\n");
			 for(i1=0;i1<1;i1++){
	                        for(j1=0;j1<blocksize;j1++)
                                   {
					printf("%d ",fdata[i1][j1+blocksize]);
			           }
				printf( " \n");
				}
			printf( " encoder fdata  n=10:\n");
			 for(i1=0;i1<1;i1++){
	                        for(j1=0;j1<blocksize;j1++)
                                   {
					printf("%d ",fdata[i1][j1+9*blocksize]);
			           }
				printf( " \n");
				}*/


				
printf(" 0 \n\n");
printf(" bit_operation_start:\n");		
 timing_set(&q3);       

	 for(i=0;i<M;i++){	
            if( i%2 == 0){////i =0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,46,48,56,64,72,74,76,78,80,82,84,86,88,90,92,94,96,98,100,102,104,106,108,110,112,114,116,118,120,122,124,126
/*printf( " value=extra1:\n");
for(i1 = 0; i1< blocksize; i1++) {
printf("%d ",*(fdata[0] +blocksize+ i1));}
printf( " \n");	     
printf( " value=extra2:\n");
for(i1 = 0; i1< blocksize; i1++) {
printf("%d ",*(fdata[0+1] + i1));}
printf( " \n");	*/
		for (i1=0;i1<blocksize;i1++){
			extra1[i1]=*(fdata[i] +blocksize+ i1);
			extra2[i1]=*(fdata[i+1] + i1);}
/*printf( " extra1:\n");
for(i1 = 0; i1< blocksize; i1++) {
printf("%d ",extra1[i1]);}
printf( " \n");
printf( " extra2:\n");
for(i1 = 0; i1< blocksize; i1++) {
printf("%d ",extra2[i1]);}
printf( " \n");


                
printf( " 1 before first times1   *r:\n");
for(i1 = 0; i1< blocksize; i1++) {
printf("%d ",*(fdata[0] + blocksize+i1));}
printf( " \n");*/


		galois_w08_region_multiply((fdata[i] + blocksize), r, blocksize, fdata[i+1], 0);
	

/*printf( " 2 after first times1 *r:\n");
for(i1 = 0; i1< blocksize; i1++) {
printf("%d ",*(fdata[0+1] +i1));}
printf( " \n");

printf( " 3--extra2value  before second times2   *r:\n");
for(i1 = 0; i1< blocksize; i1++) {
printf("%d ",*(extra2+i1));}
printf( " \n");*/


		galois_w08_region_multiply(extra2 , r, blocksize, (fdata[i] + blocksize), 0);


/*printf( " 4 afte second times2  * r:\n");
for(i1 = 0; i1< blocksize; i1++) {
printf("%d ",*(fdata[0] + blocksize+i1));}
printf( " \n");	



printf( " 5--extra1value  before XOR1   extra1 XOR 4:\n");
for(i1 = 0; i1< blocksize; i1++) {
printf("%d ",*(extra1+i1));}
printf( " \n");*/

		galois_w8_region_xor(extra1, (fdata[i] + blocksize),blocksize);
/*printf( " 6 after XOR1  compare 4:\n");
for(i1 = 0; i1< blocksize; i1++) {
printf("%d ",*(fdata[0] + blocksize+i1));}
printf( " \n");

printf( " 7 before XOR2:\n");
for(i1 = 0; i1< blocksize; i1++) {
printf("%d ",*(fdata[0+1]+i1));}
printf( " \n");

printf( " 9--extra2value    XOR 7:\n");
for(i1 = 0; i1< blocksize; i1++) {
printf("%d ",*(extra2+i1));}
printf( " \n");*/
		galois_w8_region_xor(extra2, (fdata[i+1]), blocksize);

/*printf( " 8 after XOR2  compare 7:\n");
for(i1 = 0; i1< blocksize; i1++) {
printf("%d ",*(fdata[0+1]+i1));}
printf( " \n");	*/		
}}




       for(i=0;i<M;i++){
	   if( i%4 == 0 ){
               // save original value		
		for(j=0;j<2;j++){
		for (i1=0;i1<blocksize;i1++){
		extra[j][i1]= *(fdata[i+j]+3*blocksize+i1);}}                   
		    for(j=2;j<4;j++){
		    for (i1=0;i1<blocksize;i1++){
                extra[j][i1] = *(fdata[i+j]+2*blocksize+i1);}}

 		    for(j=0;j<2;j++){
		galois_w08_region_multiply((fdata[i+j]+3*blocksize), r, blocksize, (fdata[i+j+2]+2*blocksize), 0);}
		
		    for(j=2;j<4;j++){
	
	        galois_w08_region_multiply(extra[j], r, blocksize, (fdata[i+j-2]+3*blocksize), 0);}

		    for(j=0;j<2;j++){
                 
		galois_w8_region_xor((extra[j]), (fdata[i+j]  + 3*blocksize),blocksize);}
		
		    for(j=2;j<4;j++){
	
		 galois_w8_region_xor(extra[j], (fdata[i+j]  + 2*blocksize),blocksize);}}}
		


printf(" 1 \n\n");
	for(i=0;i<M;i++){
	 if( i%8 == 0 ){//i = 0,8,16,32,40,48,56,64,72,80,88,96,104,112,120
		int t4=8;
     

                for(j=0;j<t4/2;j++){
		for (i1=0;i1<blocksize;i1++){
		extra[j][i1]= *(fdata[i+j]+5*blocksize+i1);}}                   
		    for(j=t4/2;j<t4;j++){
		    for (i1=0;i1<blocksize;i1++){
                extra[j][i1] = *(fdata[i+j]+4*blocksize+i1);}}


 		    for(j=0;j<t4/2;j++){
		galois_w08_region_multiply((fdata[i+j]+5*blocksize), r, blocksize, (fdata[i+j+t4/2]+4*blocksize), 0);}
		
		    for(j=t4/2;j<t4;j++){
		galois_w08_region_multiply(extra[j], r, blocksize, (fdata[i+j-t4/2]+5*blocksize), 0);}

		    for(j=0;j<t4/2;j++){
		
		galois_w8_region_xor(extra[j], (fdata[i+j]  + 5*blocksize),blocksize);}
		    for(j=t4/2;j<t4;j++){
		
		galois_w8_region_xor(extra[j], (fdata[i+j]  + 4*blocksize),blocksize);}}}
printf(" 2 \n\n");	
	for(i=0;i<M;i++){
         if( i%16 == 0 ){//i=0,16,32,48,64,80,96,112
		int t3=16;
             

		for(j=0;j<t3/2;j++){
		for (i1=0;i1<blocksize;i1++){
		extra[j][i1]= *(fdata[i+j]+7*blocksize+i1);}}                   
		    for(j=t3/2;j<t3;j++){
		    for (i1=0;i1<blocksize;i1++){
                extra[j][i1] = *(fdata[i+j]+6*blocksize+i1);}}

 		    for(j=0;j<t3/2;j++){
		galois_w08_region_multiply((fdata[i+j]+7*blocksize), r, blocksize, (fdata[i+j+t3/2]+6*blocksize), 0);}
		
		    for(j=t3/2;j<t3;j++){
		galois_w08_region_multiply(extra[j], r, blocksize, (fdata[i+j-t3/2]+7*blocksize), 0);}

		    for(j=0;j<t3/2;j++){
	
		galois_w8_region_xor(extra[j], (fdata[i+j]  + 7*blocksize),blocksize);}
		    for(j=t3/2;j<t3;j++){
		
		galois_w8_region_xor(extra[j], (fdata[i+j]  + 6*blocksize),blocksize);}}}
printf(" 3\n\n");
	for(i=0;i<M;i++){
         if( i%32 == 0 ){//i = 0,32 ,64,96
		int t2=32;
              

		for(j=0;j<t2/2;j++){
		for (i1=0;i1<blocksize;i1++){
		extra[j][i1]= *(fdata[i+j]+9*blocksize+i1);}}                   
		    for(j=t2/2;j<t2;j++){
		    for (i1=0;i1<blocksize;i1++){
                extra[j][i1] = *(fdata[i+j]+8*blocksize+i1);}}

 		    for(j=0;j<t2/2;j++){
		galois_w08_region_multiply((fdata[i+j]+9*blocksize), r, blocksize, (fdata[i+j+t2/2]+8*blocksize), 0);}
		
		    for(j=t2/2;j<t2;j++){
		galois_w08_region_multiply(extra[j], r, blocksize, (fdata[i+j-t2/2]+9*blocksize), 0);}

		    for(j=0;j<t2/2;j++){
	
		galois_w8_region_xor(extra[j], (fdata[i+j]  + 9*blocksize),blocksize);}
		    for(j=t2/2;j<t2;j++){

		galois_w8_region_xor(extra[j], (fdata[i+j]  + 8*blocksize),blocksize);}}}
printf(" 4 \n\n");
     for(i=0;i<M;i++){
       if( i%64 == 0 ){//i=0 ,64
		int t1=64;
		for(j=0;j<t1/2;j++){
		for (i1=0;i1<blocksize;i1++){
		extra[j][i1]= *(fcoding[i+j]+blocksize+i1);}}                   
		    for(j=t1/2;j<t1;j++){
		    for (i1=0;i1<blocksize;i1++){
                extra[j][i1] = *(fcoding[i+j]+i1);}}


 		    for(j=0;j<t1/2;j++){
		galois_w08_region_multiply((fcoding[i+j]+blocksize), r, blocksize, (fcoding[i+j+t1/2]), 0);}
		
		    for(j=t1/2;j<t1;j++){
		galois_w08_region_multiply(extra[j], r, blocksize, (fcoding[i+j-t1/2]+blocksize), 0);}

		    for(j=0;j<t1/2;j++){
	
		galois_w8_region_xor(extra[j], (fcoding[i+j]  + blocksize),blocksize);}
		    for(j=t1/2;j<t1;j++){
	
		galois_w8_region_xor(extra[j], fcoding[i+j]  ,blocksize);}}}
printf(" 5 \n\n");
  for(i=0;i<M;i++){
    if( i%128 == 0 ){//i=0
		int t=128;
              for(j=0;j<t/2;j++){
		for (i1=0;i1<blocksize;i1++){
		extra[j][i1]= *(fcoding[i+j]+3*blocksize+i1);}}                   
		    for(j=t/2;j<t;j++){
		    for (i1=0;i1<blocksize;i1++){
                extra[j][i1] = *(fcoding[i+j]+2*blocksize+i1);}}


 		    for(j=0;j<t/2;j++){
		galois_w08_region_multiply((fcoding[i+j]+3*blocksize), r, blocksize, (fcoding[i+j+t/2]+2*blocksize), 0);}
		
		    for(j=t/2;j<t;j++){
		galois_w08_region_multiply(extra[j], r, blocksize, (fcoding[i+j-t/2]+3*blocksize), 0);}

		    for(j=0;j<t/2;j++){

		  galois_w8_region_xor((extra[j]), (fcoding[i+j]  + 3*blocksize),blocksize);}
		    for(j=t/2;j<t;j++){
	
		galois_w8_region_xor((extra[j]), (fcoding[i+j]  + 2*blocksize),blocksize);}}}
printf(" 6 \n\n");
     
timing_set(&q4);
timing_set(&q6);
timing_set(&t4);
printf("bit operation ended \n");
						//}
	    /*  printf( " after operation fcoding first strip  :\n");
	  		
	
			 for(i1=0;i1<1;i1++){
	                        for(j1=0;j1<m*blocksize;j1++){
                                  
				printf("%d ",fcoding[i1][j1]);}
     				printf( " \n");	}
			   printf( " \n");
printf( " after operation fdata first strip  :\n");
	  		
	
			 for(i1=0;i1<1;i1++){
	                        for(j1=0;j1<k*blocksize;j1++){
                                  
				printf("%d ",fdata[i1][j1]);}
     				printf( " \n");	}
			   printf( " \n");
 printf( " after operation fdata n=3 3 row  :\n");
	  		
	
			 for(i1=2;i1<3;i1++){
	                        for(j1=0;j1<blocksize;j1++){
                                  
				printf("%d ",fdata[i1][j1+2*blocksize]);}
     				printf( " \n");	}
			   printf( " \n");*/




		/* Write data and encoded data to k+m files */
		for	(i = 0; i < k; i++) {
		 
			if (fp == NULL) {
				bzero(fdata[i], blocksize);
 			} else {
				sprintf(fname, "%s/Coding/%s_k%0*d%s", curdir, s1, md, i, extension);
				if (n == 1) {
					fp2 = fopen(fname, "wb");
				}
				else {
				
				fp2 = fopen(fname, "ab");
				}
				for(j=0;j<M;j++){
				fwrite(&fdata[j][(i)*blocksize], sizeof(char), blocksize, fp2);}
				
				fclose(fp2);
			}
			
		}
             //printf(" outComplete33\n\n");
		for	(i = 0; i < m; i++) {
			if (fp == NULL) {
				bzero(fdata[i], blocksize);
 			} else {
				sprintf(fname, "%s/Coding/%s_m%0*d%s", curdir, s1, md, i, extension);
				//if (n == 1) {
				//	fp2 = fopen(fname, "wb");
				//}
				//else {
					fp2 = fopen(fname, "ab");
				//}
				for(j=0;j<M;j++){
				fwrite(&fcoding[j][(i)*blocksize], sizeof(char), blocksize, fp2);}
				fclose(fp2);
			}
		}
		n++;
		/* Calculate encoding time */

		totalsec += timing_delta(&t3, &t4);
		encode_time= timing_delta(&q1, &q2);
		bit_operation_time= timing_delta(&q3, &q4);
		sum_time= timing_delta(&q5, &q6);
	}

	/* Create metadata file */
        if (fp != NULL) {
		sprintf(fname, "%s/Coding/%s_meta.txt", curdir, s1);
		fp2 = fopen(fname, "wb");
		fprintf(fp2, "%s\n", argv[1]);
		fprintf(fp2, "%d\n", size);
		fprintf(fp2, "%d %d %d %d %d\n", k, m, w, packetsize, buffersize);
		fprintf(fp2, "%s\n", argv[4]);
		fprintf(fp2, "%d\n", tech);
		fprintf(fp2, "%d\n", readins);
		fclose(fp2);
	}


	/* Free allocated memory */
	free(s1);
	free(fname);
	free(block);
	free(curdir);
	
	/* Calculate rate in MB/sec and print */
	timing_set(&t2);
	tsec = timing_delta(&t1, &t2);
	printf("Encoding (MB/sec): %0.10f\n", (((double) size)/1024.0/1024.0)/totalsec);
	printf("En_Total (MB/sec): %0.10f\n", (((double) size)/1024.0/1024.0)/tsec);
	printf("encode_time (sec): %0.10f\n", encode_time);
	printf("bit_operation_time (sec): %0.10f\n", bit_operation_time);
	printf("sum_time (sec): %0.10f\n", sum_time);
	return 0;
}


/* is_prime returns 1 if number if prime, 0 if not prime */
int is_prime(int w) {
	int prime55[] = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,
	    73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,
		    181,191,193,197,199,211,223,227,229,233,239,241,251,257};
	int i;
	for (i = 0; i < 55; i++) {
		if (w%prime55[i] == 0) {
			if (w == prime55[i]) return 1;
			else { return 0; }
		}
	}
	assert(0);
}

/* Handles ctrl-\ event */
void ctrl_bs_handler(int dummy) {
	time_t mytime;
	mytime = time(0);
	fprintf(stderr, "\n%s\n", ctime(&mytime));
	fprintf(stderr, "You just typed ctrl-\\ in encoder.c.\n");
	fprintf(stderr, "Total number of read ins = %d\n", readins);
	fprintf(stderr, "Current read in: %d\n", n);
	fprintf(stderr, "Method: %s\n\n", Methods[method]);	
	signal(SIGQUIT, ctrl_bs_handler);
}



























