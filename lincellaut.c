
/* A simple program for calculating various things about linear
 *  cellular automaton.
 *
 * Mar 27, 2022
 *
 */
 
/* The following resources were used as a reference:
https://docs.microsoft.com/en-us/cpp/c-language/cpp-integer-limits?view=msvc-170
https://stackoverflow.com/questions/3219393
*/
 
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h> //So I can get the maximum integer

#include "headers/helper.h"
#include "headers/bigint.h" //Arbitrary precision

#include "headers/linalg.h" 
#include "headers/cycles.h"  //Allows us to use Floyd's Algorithm
#include "headers/modular.h" //Modular square roots and inverses
#include "headers/fibonacci.h"
#include "headers/factors.h" //For LCM()

#include "headers/algebra.h"

#define FREE(v) free(v); v = NULL

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"


int main(int argc, char* argv[])
{
	//Check to see if MAXBUNCH for BigIntT structs is a power of ten
	//If it isn't, warn the user
	for (int i = 10; i < MAXBUNCH; i *= 10)
		if (MAXBUNCH % i != 0)
		{
			printf("MAXBUNCH is not a power of ten. Some functionality involving the use" \
			" of BigIntT structs may not work correctly. Change the MAXBUNCH constant within" \
			" bigint.c to fix this issue.\n");
			break;
		}
	
	//The maximum length for a string in the .config file
	const int MAXSTRLEN = 101;
	
	//Read .config file to get appropriate data loaded
	FILE* system = fopen("config/system.config", "r");
	
	if (system == NULL)
	{
		fprintf(stderr, "Unable to open config file. Check the config folder.\n");
		return EXIT_FAILURE;
	}
	
	int   iterations;
	int   modulus;
	
	//For holding str representation of mod, used for BigIntT number (if needed)
	char* bigintmodstring = malloc(MAXSTRLEN*sizeof(char));
	
	char* updatefilepath  = malloc(MAXSTRLEN*sizeof(char));
	char* initialfilepath = malloc(MAXSTRLEN*sizeof(char));
	char* iterfilepath    = malloc(MAXSTRLEN*sizeof(char));
	
	//Temporarily holds config data
	char* systemData = malloc(MAXSTRLEN*sizeof(char));
	
	char* tempStr; //Holds temporary info when using string functions
	
	
	while (fscanf(system, "%10s", systemData) == 1)
	{
		//printf("%s\n", systemData);
		if (! strcmp(systemData, "mod"))
		{
			//Getting modulus as string first in case we need to convert to BigIntT
			if (fscanf(system, "%101s", bigintmodstring) != 1)
			{
				fprintf(stderr, "Unable to read modulus from config file.\n");
				return EXIT_FAILURE;
			}
			
			if (strlen(bigintmodstring) < 10)
			{
				modulus = (int)strtol(bigintmodstring, &tempStr, 10);
				if (tempStr[0] != '\0')
				{
					fprintf(stderr, "Invalid modulus provided in config file.\n");
					return EXIT_FAILURE;
				}
			}
		}
		
		else if (! strcmp(systemData, "iterations"))
		{
			if (fscanf(system, "%d", &iterations) != 1)
			{
				fprintf(stderr, "Unable to read number of iterations from config file.\n");
				return EXIT_FAILURE;
			}
		}
		
		else if (! strcmp(systemData, "update"))
		{
			if (fscanf(system, "%s", updatefilepath) != 1)
			{
				fprintf(stderr, "Unable to read update matrix path from config file.\n");
				return EXIT_FAILURE;
			}
		}
		
		else if (! strcmp(systemData, "initial"))
		{
			if (fscanf(system, "%s", initialfilepath) != 1)
			{
				fprintf(stderr, "Unable to read initial vector path from config file.\n");
				return EXIT_FAILURE;
			}
		}
		
		else if (! strcmp(systemData, "itername"))
		{
			if (fscanf(system, "%s", iterfilepath) != 1)
			{
				fprintf(stderr, "Unable to read iteration file path from config file.\n");
				return EXIT_FAILURE;
			}
		}
	}
	
	FREE(systemData);
	if (fclose(system) == EOF)
	{
		fprintf(stderr, "Unable to close config file.\n");
		return EXIT_FAILURE;
	}
	
	//If we have command line arguments
	if (argc > 1)
	{
		//If we want to iterate the given update matrix
		if (! strcmp(argv[1], "iterate"))
		{
			//If the user provided a custom number of iterations
			if (argc > 2)
			{
				iterations = (int)strtol(argv[2], &tempStr, 10);
				
				//If we didn't read any digits for iterations
				if (tempStr == argv[2])
				{
					fprintf(stderr, "Invalid number of iterations provided at command line.\n");
					return EXIT_FAILURE;
				}
			}
			
			IntMatrixTP F;
			IntMatrixTP F_2;
			IntMatrixTP F_result;
			
			if (iterations < 0)
			{
				F_2 = read_IntMatrixT(updatefilepath);
				F   = inverse(F_2, modulus);
				
				if (F == NULL)
					printf("An inverse for the update matrix mod %d does not exist.\n", modulus);
				
				else
				{
					iterations *= -1;
					F_2 = free_IntMatrixT(F_2);
					F_2 = read_IntMatrixT(initialfilepath);
				}
			}
			
			else
			{
				F   = read_IntMatrixT(updatefilepath);
				F_2 = read_IntMatrixT(initialfilepath);
			}
			
			if ((F == NULL) || (F_2 == NULL))
				return EXIT_FAILURE;
			
			//Iterate F a few times
			//The minus 1 is for easier conversion between ORBITVIS results
			printf("Iterations: %d\n", iterations);
			if (iterations > 0)
				F_result = iterate(F, F_2, modulus, iterations-1);
			
			//Prevents the matrix from being printed when
			// an inverse doesn't exist and the iterations
			// provided was negative
			if (iterations > 0)
				printm(F_result);
			
			else if (iterations == 0)
				printf("I\n");

			F        = free_IntMatrixT(F);
			F_2      = free_IntMatrixT(F_2);
			F_result = free_IntMatrixT(F_result);
		}
		
		
		//Find the inverse of the update matrix
		else if (!strcmp(argv[1], "inverse"))
		{
			IntMatrixTP F = read_IntMatrixT(updatefilepath);
			IntMatrixTP Finv;
			
			if (F == NULL)
				return EXIT_FAILURE;
			
			printf("Update matrix:\n");
			printm(F);
			
			//If the user specified a modulus at the command line
			if (argc > 2)
			{
				modulus = (int)strtol(argv[2], &tempStr, 10);
				
				if (tempStr[0] != '\0')
				{
					fprintf(stderr, "Invalid modulus passed at command line.\n");
					return EXIT_FAILURE;
				}
			}
			
			Finv = inverse(F, modulus);
			
			if (Finv == NULL)
				printf("An inverse for the update matrix mod %d does not exist.\n", modulus);
			
			else
			{
				printf("Inverse mod %d:\n", modulus);
				printm(Finv);
			}
			
			F    = free_IntMatrixT(F);
			Finv = free_IntMatrixT(Finv);
		}
		
		
		//If we want to calculate the determinant of our matrix
		else if (!strcmp(argv[1], "det"))
		{
			//If the user specified a custom modulus
			if (argc > 2)
			{
				modulus = (int)strtol(argv[2], &tempStr, 10);
				if (tempStr[0] != '\0')
				{
					fprintf(stderr, "Unable to read modulus from command line.\n");
					return EXIT_FAILURE;
				}
			}
			
			IntMatrixTP A = read_IntMatrixT(updatefilepath);
			if (A == NULL)
			{
				fprintf(stderr, "Unable to read .matrix file at %s.\n", updatefilepath);
				return EXIT_FAILURE;
			}
			
			printf("Matrix:\n");
			printm(A);
			printf("Determinant: %d\n", det(A) % modulus);
			
			A = free_IntMatrixT(A);
		}
		
		
		//Find the characteristic equation of the update matrix
		else if (! strcmp(argv[1], "chara"))
		{
			BigIntTP bigMod;
			BigIntMatrixTP bigMatrix;
			
			BigPolyTP bigEqn;
			BigPolyTP* bigEqnFactors;
			
			//If the user provided a custom modulus
			if (argc > 2)
			{
				if (strtoBIT(argv[2], &bigMod) != 1)
				{
					fprintf(stderr, "Unable to read modulus from command line.\n");
					return EXIT_FAILURE;
				}
			}
			
			//Extract the modulus in the .config file otherwise
			else
			{
				if (strtoBIT(bigintmodstring, &bigMod) != 1)
				{
					fprintf(stderr, "Unable to read modulus from config file.\n");
					return EXIT_FAILURE;
				}
			}
			
			//Now, read the update matrix
			bigMatrix = read_BigIntMatrixT(updatefilepath);
			
			if (bigMatrix == NULL)
			{
				fprintf(stderr, "Unable to read matrix from %s.\n", updatefilepath);
				return EXIT_FAILURE;
			}
			
			//Actually calculate the characteristic equation here
			bigEqn = chara_eqn(bigMatrix, bigMod);
			printf("Matrix:\n");
			printbm(bigMatrix);
			printf("Modulus: ");
			printi(bigMod);
			printf("\nCharacteristic equation: ");
			printp(bigEqn);
			printf("\n\n");
			
			printf("Factored characteristic equation:\n");
			bigEqnFactors = factor_BigPolyT(bigEqn, bigMod);
			
			for (int i = 0; i < degree(bigEqn); i += 1)
			{
				printp(bigEqnFactors[i]);
				printf(" ");
				bigEqnFactors[i] = free_BigPolyT(bigEqnFactors[i]);
			}
			FREE(bigEqnFactors);
			printf("\n");
			
			bigMod = free_BigIntT(bigMod);
			bigMatrix = free_BigIntMatrixT(bigMatrix);
			bigEqn = free_BigPolyT(bigEqn);
		}
		
		
		//If we want to use Floyd's Cycle Detection Algorithm
		else if (!strcmp(argv[1], "floyd"))
		{
			//If the user provided a custom modulus
			if (argc > 2)
			{
				modulus = (int)strtol(argv[2], &tempStr, 10);
				if (tempStr[0] == '\0')
				{
					fprintf(stderr, "Invalid modulus passed at command line.\n");
					return EXIT_FAILURE;
				}
			}
			
			IntMatrixTP initial = read_IntMatrixT(initialfilepath);
			IntMatrixTP update  = read_IntMatrixT(updatefilepath);
			CycleInfoTP coolCycle;
			
			if (initial == NULL)
			{
				fprintf(stderr, "Unable to read matrix in %s.\n", initialfilepath);
				return EXIT_FAILURE;
			}
			
			else if (update == NULL)
			{
				fprintf(stderr, "Unable to read matrix in %s.\n", updatefilepath);
				return EXIT_FAILURE;
			}
			
			//If we actually get all the data we need
			else
			{
				coolCycle = floyd(update, initial, modulus);
				
				if (coolCycle == NULL)
					printf("Update matrix provided is not a square matrix.\n");
				
				else
					printcycle(coolCycle);
			}
			
			initial   = free_IntMatrixT(initial);
			update    = free_IntMatrixT(update);
			coolCycle = free_CycleInfoT(coolCycle);
		}
		
		
		//Same as floyd, but uses big datatypes
		else if (! strcmp(argv[1], "bigfloyd"))
		{
			BigIntTP bigModulus;
			BigIntMatrixTP initial;
			BigIntMatrixTP update;
			CycleInfoTP coolCycle = NULL;
			
			//If the user specified a custom modulus
			if (argc > 2)
			{
				if (strtoBIT(argv[2], &bigModulus) == 0)
				{
					fprintf(stderr, "Unable to read modulus from command line.\n");
					return EXIT_FAILURE;
				}
			}
			
			//Use the modulus provided in the .config file
			else
			{
				if (strtoBIT(bigintmodstring, &bigModulus) == 0)
				{
					fprintf(stderr, "Unable to read modulus from .config file. It's possible the modulus" \
					" is too big, in which case it must be provided at the command line.\n");
					return EXIT_FAILURE;
				}
			}
			
			//Now, get the matrices from the provided files
			initial = read_BigIntMatrixT(initialfilepath);
			update  = read_BigIntMatrixT(updatefilepath);
			
			if ((initial == NULL) || (update == NULL))
			{
				fprintf(stderr, "Unable to read .matrix files.\n");
				return EXIT_FAILURE;
			}
			
			//Get cycle info
			big_floyd(update, initial, bigModulus, &coolCycle);
			
			if (coolCycle == NULL)
				printf("Update matrix is not a square matrix.\n");
			
			else
				printcycle(coolCycle);
			
			
			bigModulus = free_BigIntT(bigModulus);
			initial    = free_BigIntMatrixT(initial);
			update     = free_BigIntMatrixT(update);
			coolCycle  = free_CycleInfoT(coolCycle);
		}
		
		
		//Generate some basic rotation matrices for the given modulus
		else if (! strcmp(argv[1], "rots"))
		{
			bool has30 = TRUE;
			bool has45 = TRUE;
			
			int inv2;
			int sqrt3;
			
			//If user provided modulus on CLI
			if (argc > 2)
			{
				modulus = (int)strtol(argv[2], &tempStr, 10);
				if (tempStr[0] != '\0')
				{
					fprintf(stderr, "Unable to read modulus from command line.\n");
					return EXIT_FAILURE;
				}
			}
			
			if ((num_inverse(2, modulus) == -1))
			{
				printf("45degrot and 30degrot don't exist mod %d.\n", modulus);
				has45 = FALSE;
				has30 = FALSE;
			}
			
			else
			{
				if (square_root(2, modulus) == -1)
				{
					printf("45degrot doesn't exist mod %d.\n", modulus);
					has45 = FALSE;
				}
				
				if (square_root(3, modulus) == -1)
				{
					printf("30degrot does not exist mod %d.\n", modulus);
					has30 = FALSE;
				}
			}
			
			
			if (has45)
			{
				printf("45degrot mod %d:\n", modulus);
				
				printf("%2d %2d\n",
				(num_inverse(2, modulus)*square_root(2, modulus)) % modulus,
				modulus - ((num_inverse(2, modulus)*square_root(2, modulus)) % modulus));
				
				printf("%2d %2d\n",
				(num_inverse(2, modulus)*square_root(2, modulus)) % modulus,
				(num_inverse(2, modulus)*square_root(2, modulus)) % modulus);
			}
			
			if (has30)
			{
				inv2  = (modulus - num_inverse(2, modulus)) % modulus;
				sqrt3 = square_root(3, modulus);
				
				printf("30degrot mod %d:\n", modulus);
				
				printf("%2d %2d\n",
				(inv2*sqrt3) % modulus,
				((modulus - inv2) % modulus));
				
				printf("%2d %2d\n",
				inv2, 
				(inv2*sqrt3) % modulus);
			}
		}
		
		
		//If we want to look for matrices that have column vectors
		// that cycle with a cycle length less than the matrix itself
		else if (! strcmp(argv[1], "cycmatsearch"))
		{
			//If user didn't provide enough arguments for the tool, explain how to use it
			if (argc < 5)
			{
				printf(ANSI_COLOR_YELLOW "cycmatsearch " ANSI_COLOR_CYAN "size maxmod cycles..." ANSI_COLOR_RESET \
				": Searches for a matrix whose column vectors cycle with cycle lengths less than the matrix itself.\n");
				printf(" - " ANSI_COLOR_CYAN "size" ANSI_COLOR_RESET \
				": Tells what size matrix to use.\n");
				printf(" - " ANSI_COLOR_CYAN "maxmod" ANSI_COLOR_RESET \
				": Tells which modulus to stop searching at.\n");
				printf(" - " ANSI_COLOR_CYAN "cycles..." ANSI_COLOR_RESET \
				": A list of cycle lengths for their respective vectors.\n\n");
				
				return EXIT_SUCCESS;
			}
			
			FILE* textOutput; 
			char* textOutputName;
			
			BigIntTP maxMod;  //What's the largest modulus we should search to?
			BigIntTP currMod; //What modulus are we currently checking?
			BigIntTP one;
			BigIntTP zero;
			BigIntTP temp;    //For holding temporary results of calculations
			
			BigIntTP** currMatElements; //Holds matrix numbers so we can set the matrix easily
			
			CycleInfoTP theCycle; //Used to get matrices in a cycle
			
			BigIntMatrixTP zeroMat; //The zero matrix
			BigIntMatrixTP currMat; //Holds the matrix we're currently testing for cyclic vectors
			BigIntMatrixTP tempMat; //Holds results of matrix multiplication
			BigIntMatrixTP tempMat2;
			
			bool checkedAllMatrices; //Says whether we've exhaused all matrices under a given modulus
			bool matrixIsCyclic;     //Says whether the matrix we tested has cyclic column vectors or not
			int i, j; //Used for for-loops so we don't have to keep declaring new variables
			int size; //The size of our matrix to use
			int currIteration;
			int nextIteration; //Used to keep track of which iteration of the matrix to check next
			int cycleLCM = 1;
			
			/*
			Searched so far:
			cycmatsearch 2 20 2 3
			cycmatsearch 3 6 2 3 5
			*/
			
			int oneArr[] = {1};
			int start[] = {6};
			
			printf("Currently, the first modulus checked is not 2 for testing purposes.\n");
			
			int* colVectCycles; //Holds the different cycle lengths 
			
			size = (int)strtol(argv[2], &tempStr, 10);
			
			if (tempStr[0] != '\0')
			{
				fprintf(stderr, "Unable to read matrix size from command line.\n");
				return EXIT_FAILURE;
			}
			
			//If we couldn't read the modulus
			if (strtoBIT(argv[3], &maxMod) == 0)
			{
				fprintf(stderr, "Unable to read modulus from command line.\n");
			}
			
			//If the user didn't provide enough cycles for the given matrix size
			if (argc < 4 + size)
			{
				printf("Too few cycles passed for given matrix size.\n");
				return EXIT_FAILURE;
			}
			
			//Get column vector cycles
			colVectCycles = malloc(size*sizeof(int));
			for (i = 0; i < size; i += 1)
			{
				colVectCycles[i] = (int)strtol(argv[4+i], &tempStr, 10);
				if (tempStr[0] != '\0')
				{
					fprintf(stderr, "Unable to read cycle lengths from command line.\n");
					return EXIT_FAILURE;
				}
			}
			
			//Start looking for matrices at mod 2
			currMod = new_BigIntT(start, 1);
			one     = new_BigIntT(oneArr, 1);
			zero    = empty_BigIntT(1);
			temp    = empty_BigIntT(1);
			
			//Initialise our matrix elements
			currMatElements = malloc(size*sizeof(BigIntTP*));
			for (i = 0; i < size; i += 1)
			{
				currMatElements[i] = malloc(size*sizeof(BigIntTP));
				for (j = 0; j < size; j += 1)
					currMatElements[i][j] = empty_BigIntT(1);
			}
			
			//TESTING A SPECIFIC CASE. DELETE THIS LATER
			/*
			int theModValue[] = {5};
			int theResumeValue[] = {4};
			
			currMatElements = malloc(size*sizeof(BigIntTP*));
			for (i = 0; i < size; i += 1)
			{
				currMatElements[i] = malloc(size*sizeof(BigIntTP));
				for (j = 0; j < size; j += 1)
				{
					if ((i == size-1) && (j == size-1))
						currMatElements[i][j] = new_BigIntT(theResumeValue, 1);
					else
						currMatElements[i][j] = new_BigIntT(theModValue, 1);
				}
			}
			*/
			
			//Find the LCM of our cycles
			for (i = 0; i < size; i += 1)
				cycleLCM = LCM(cycleLCM, colVectCycles[i]);
			
			zeroMat  = new_BigIntMatrixT(size, size);
			currMat  = new_BigIntMatrixT(size, size);
			tempMat  = new_BigIntMatrixT(size, size);
			tempMat2 = new_BigIntMatrixT(size, size);
			
			//Initialising a reusable CycleInfoTP to save time on
			// memory allocations.
			theCycle = new_CycleInfoT();
			
			//Forming name for text file output.
			textOutputName = malloc(MAXSTRLEN*sizeof(char));
			textOutputName[0] = '\0';
			strcat(textOutputName, argv[1]);
			strcat(textOutputName, " ");
			strcat(textOutputName, argv[2]);
			strcat(textOutputName, " ");
			strcat(textOutputName, argv[3]);
			strcat(textOutputName, " ");
			for (i = 0; i < size-1; i += 1)
			{
				strcat(textOutputName, argv[4+i]);
				strcat(textOutputName, " ");
			}
			strcat(textOutputName, argv[size+3]);
			strcat(textOutputName, ".txt");
			//printf("%s\n", textOutputName);
			
			textOutput = fopen(textOutputName, "w");
			if (textOutput == NULL)
			{
				fprintf(stderr, "Unable to create/open %s for writing.\n", textOutputName);
				return EXIT_FAILURE;
			}
			
			//Search all moduli until we get to the specified limit
			while (compare_BigIntT(currMod, maxMod) <= 0)
			{
				printf("Currently searching mod ");
				printi(currMod);
				printf("...\n");

				//Search throuh all matrices under the current modulus
				checkedAllMatrices = FALSE;
				while (!checkedAllMatrices)
				{
					//Getting all the relevant information we need to check for vector cycles
					set_big_matrix(currMat, currMatElements);
					big_floyd(currMat, currMat, currMod, &theCycle);
					copy_BigIntMatrixT(rep(theCycle), tempMat);
					
					currIteration = 0;
					
					//Prevent matrices that don't have the correct cycle length
					// from being printed
					matrixIsCyclic = FALSE; 
				
					//Check to see if the vectors actually cycle how we want them to
					//The matrix's cycle must be the LCM of the cycle lengths if
					// the vectors are to behave how we want them to.
					while ((omega(theCycle) == cycleLCM))
					{
						//Assume the matrix is cyclic until proven otherwise
						matrixIsCyclic = TRUE;
						
						//Get a value of nextIteration that's greater than before,
						// then work from there
						for (i = 0; i < size; i += 1)
							if (colVectCycles[i] > currIteration)
							{
								nextIteration = colVectCycles[i];
								break;
							}
						
						//Looking for the next smallest iteration we need to check
						for (i = 1; i < size; i += 1)
							if ((colVectCycles[i] < nextIteration) && (colVectCycles[i] > currIteration))
								nextIteration = colVectCycles[i];
							
						//If we've checked all the iterations we need to check
						if (nextIteration == currIteration)
							break;
								
						//Iterate matrix to the correct iteration
						for (i = 0; i < nextIteration - currIteration; i += 1)
						{
							big_mat_mul(currMat, tempMat, tempMat2);
							modbm(tempMat2, currMod);
							copy_BigIntMatrixT(tempMat2, tempMat);
						}
						
						//Check for repeated vectors
						for (i = 0; i < size; i += 1)
							if (colVectCycles[i] == nextIteration)
								if (compare_BigIntMatrixT_cols(tempMat, rep(theCycle), i) != 1)
								{
									matrixIsCyclic = FALSE;
									break;
								}
								
							
						currIteration = nextIteration;
						
						//If we've determined the matrix doesn't work
						if (!matrixIsCyclic)
							break;
					}
					
					//If we found a matrix that cycles how we want, print it
					if ((matrixIsCyclic))
					{
						printf("Mod ");
						printi(currMod);
						printf("\n");
						printbm(currMat);
						printf("\n");
						
						fprintf(textOutput, "Mod ");
						fprinti(textOutput, currMod);
						fprintf(textOutput, "\n");
						fprintbm(textOutput, currMat);
						fprintf(textOutput, "\n");
					}
					
					
					//Iterate to the next matrix
					checkedAllMatrices = TRUE;
					for (i = 0; i < size; i += 1)
					{
						for (j = 0; j < size; j += 1)
						{
							//Increment specific element; check if it overflows
							add_BigIntT(one, currMatElements[i][j], temp);
							
							if (compare_BigIntT(temp, currMod) >= 0)
							{
								copy_BigIntT(zero, currMatElements[i][j]);
								
								//Print a progress update to the console
								if ((i == size-1) && (j == size-2))
								{
									printi(currMatElements[i][j+1]);
									printf(" / ");
									printi(currMod);
									printf(" checked...\n");
								}
								
							}
							
							//No overflow
							else
							{
								checkedAllMatrices = FALSE;
								copy_BigIntT(temp, currMatElements[i][j]);
								i = size;
								j = size;
							}
						}
					}
				}
				
				//Increment modulus
				add_BigIntT(one, currMod, temp);
				copy_BigIntT(temp, currMod);
			}
			
			if (fclose(textOutput) == EOF)
				fprintf(stderr, "Unable to close %s.\n", textOutputName);
			
			//Freeing memory
			maxMod  = free_BigIntT(maxMod);
			currMod = free_BigIntT(currMod);
			one     = free_BigIntT(one);
			zero    = free_BigIntT(zero);
			temp    = free_BigIntT(temp);
			
			currMat  = free_BigIntMatrixT(currMat);
			zeroMat  = free_BigIntMatrixT(zeroMat);
			tempMat  = free_BigIntMatrixT(tempMat);
			tempMat2 = free_BigIntMatrixT(tempMat2);
			
			theCycle = free_CycleInfoT(theCycle);
			
			FREE(colVectCycles);
			
			for (i = 0; i < size; i += 1)
			{
				for (j = 0; j < size; j += 1)
					currMatElements[i][j] = free_BigIntT(currMatElements[i][j]);
				
				free(currMatElements[i]);
			}
			FREE(currMatElements);
			FREE(textOutputName);
			
			printf("Finished searching.\n");
		}
		
		
		//If we want to generate a Fibonacci cycle for
		// the given initial vector
		else if (! strcmp(argv[1], "fibcycle"))
		{
			int* initVect = malloc(2*sizeof(int));
			
			FILE* initVectFile = fopen(initialfilepath, "r");
			
			//If the user specified a specific modulus to use
			if (argc > 2)
			{
				modulus = (int)strtol(argv[2], &tempStr, 10);
				
				if (tempStr[0] != '\0')
				{
					fprintf(stderr, "Invalid modulus passed at command line.\n");
					return EXIT_FAILURE;
				}
			}
			
			if (initVectFile == NULL)
			{
				fprintf(stderr, "Unable to read file %s.\n", initialfilepath);
				return EXIT_FAILURE;
			}
			
			//Read dimensions of provided vector
			if (fscanf(initVectFile, "%d %d", &initVect[0], &initVect[1]) != 2)
			{
				fprintf(stderr, "Unable to read data from %s.\n", initialfilepath);
				return EXIT_FAILURE;
			}
			
			//If dimensions are incorrect
			if ((initVect[0] != 2) || (initVect[1] != 1))
			{
				fprintf(stderr, "Vector in %s must be a 2 by 1 vector.\n", initialfilepath);
				return EXIT_FAILURE;
			}
			
			//Get actual vector data from file
			if (fscanf(initVectFile, "%d %d", &initVect[0], &initVect[1]) != 2)
			{
				fprintf(stderr, "Unable to read data from %s.\n", initialfilepath);
				return EXIT_FAILURE;
			}
			
			if (fclose(initVectFile) == EOF)
			{
				fprintf(stderr, "Unable to close file %s.\n", initialfilepath);
				return EXIT_FAILURE;
			}
			
			//Print cycle and get cycle length
			printf("Length of cycle: %d\n", generate_orbit(initVect, modulus, TRUE));
			
			FREE(initVect);
		}
		
		
		//If we want to see all possible cycle lengths for a particular modulus
		else if (! strcmp(argv[1], "fibcyclelens"))
		{
			int justOne[] = {1};
			BigIntTP bigModulus;
			
			int* cycleLengths = malloc(30*sizeof(int));
			int  indexCounter = 0;
			int  counter;
			
			//Fix hardcoding later
			BigIntTP startingX = empty_BigIntT(1);
			BigIntTP startingY = empty_BigIntT(1);
			
			BigIntTP currX = empty_BigIntT(1);
			BigIntTP currY = empty_BigIntT(1);
			
			BigIntTP tempInt = empty_BigIntT(1);
			BigIntTP one     = new_BigIntT(justOne, 1);
			BigIntTP zero    = empty_BigIntT(1);
			
			int ok;
			
			//Use modulus provided on CLI
			if (argc == 3)
				ok = strtoBIT(argv[2], &bigModulus);
			
			//Use config modulus instead
			else
				ok = strtoBIT(bigintmodstring, &bigModulus);
			
			//If we read the modulus alright
			if (ok)
			{
				printf("Modulus: ");
				printi(bigModulus);
				printf("\n");
				
				//Iterate through all possible starting vectors
				while (compare_BigIntT(startingX, bigModulus) < 0)
				{
					while (compare_BigIntT(startingY, bigModulus) < 0)
					{
						//Prepping our starting vector
						copy_BigIntT(startingX, currX);
						copy_BigIntT(startingY, currY);
						counter = 0;
						
						//Increment until we get back to the start
						do
						{
							//Perform the iteration
							add_BigIntT(currX, currY, tempInt);
							mod_BigIntT(tempInt, bigModulus, currY);
							copy_BigIntT(currX, tempInt);
							copy_BigIntT(currY, currX);
							copy_BigIntT(tempInt, currY);
							
							counter += 1;
						}
						while ((compare_BigIntT(currX, startingX) != 0) ||
						       (compare_BigIntT(currY, startingY) != 0));
									 
						//Now, check to see if we've already gotten this cycle length
						ok = 1;
						for (int x = 0; x < indexCounter; x += 1)
						{
							if (cycleLengths[x] == counter)
							{
								ok = 0;
								break;
							}
						}
							
						//If we didn't find it, add the new cycle length
						if (ok)
						{
							cycleLengths[indexCounter] = counter;
							indexCounter += 1;
							printf("Cycle length: %d\n", counter);
							printf("Starting vector: [");
							printi(startingX);
							printf(", ");
							printi(startingY);
							printf("]\n");
						}
						
						//Incrementing Y
						add_BigIntT(startingY, one, tempInt);
						copy_BigIntT(tempInt, startingY);
					}
					//Incrementing X, resetting Y
					copy_BigIntT(zero, startingY);
					add_BigIntT(startingX, one, tempInt);
					copy_BigIntT(tempInt, startingX);
				}
			}
			
			else
			{
				fprintf(stderr, "Unable to use provided modulus.\n");
				return EXIT_FAILURE;
			}
			
			bigModulus = (bigModulus == NULL) ? NULL : free_BigIntT(bigModulus);
			startingX = free_BigIntT(startingX);
			startingY = free_BigIntT(startingY);
			currX     = free_BigIntT(currX);
			currY     = free_BigIntT(currY);
			
			tempInt  = free_BigIntT(tempInt);
			one      = free_BigIntT(one);
			zero     = free_BigIntT(zero);
			
			FREE(cycleLengths);
		}
		
		
		//If we want to check the Fibonacci numbers to see if multiples of numbers
		// appear before multiples of powers of those numbers
		else if (! strcmp(argv[1], "fibmultsearch"))
		{
			int hundred[] = {100};
			int oneArr[]  = {001};
			int zeroArr[] = {000};
			int start[]   = {001}; //Where the program starts counting
			//int debugCounter = 0; //For debugging
			BigIntTP upperbound;
			BigIntTP currNum = new_BigIntT(start, 1);
			BigIntTP counter = new_BigIntT(oneArr, 1);
			
			BigIntTP zero    = new_BigIntT(zeroArr, 1);
			BigIntTP one     = new_BigIntT(oneArr, 1);
			
			BigIntTP fibA    = new_BigIntT(oneArr, 1);
			BigIntTP fibB    = empty_BigIntT(1);
			BigIntTP fibTemp = empty_BigIntT(1);
			
			if (argc > 2)
			{
				if (! strtoBIT(argv[2], &upperbound))
				{
					fprintf(stderr, "Invalid upper bound passed.\n");
					return EXIT_FAILURE;
				}
			}
			else
				upperbound = new_BigIntT(hundred, 1);
			
			printf("pos\tnum\t\tfib\n");
			while (compare_BigIntT(upperbound, currNum) >= 0)
			{	
				//Testing to see if our current Fib number is a multiple of currNum
				mod_BigIntT(fibA, currNum, fibTemp);
				
				if (compare_BigIntT(zero, fibTemp) == 0)
				{
					printi(counter);
					printf("\t");
					
					printi(currNum);
					printf("\tdivides\t");
					printi(fibA);
					
					//Now we check to see if it's a nultiple of some power of currNum
					divide_BigIntT(fibA, currNum, fibTemp);
					copy_BigIntT(fibTemp, fibA);
					mod_BigIntT(fibA, currNum, fibTemp);
					
					if (compare_BigIntT(zero, fibTemp) == 0)
					{
						printf(", along with ");
						printi(currNum);
						printf("^2");
					}
					
					printf("\n");
					
					//Prepare for the next number to check
					add_BigIntT(currNum, one, fibTemp);
					copy_BigIntT(fibTemp, currNum);
					copy_BigIntT(one, fibA);
					copy_BigIntT(zero, fibB);
					copy_BigIntT(one, counter);
				}
				
				else
				{
					//Calculate next Fibonacci number
					add_BigIntT(fibA, fibB, fibTemp);
					copy_BigIntT(fibA, fibB);
					copy_BigIntT(fibTemp, fibA);
					
					add_BigIntT(one, counter, fibTemp);
					copy_BigIntT(fibTemp, counter);
				}
			}
			
			upperbound = free_BigIntT(upperbound);
			currNum    = free_BigIntT(currNum);
			counter    = free_BigIntT(counter); 
			zero       = free_BigIntT(zero);
			one        = free_BigIntT(one);
			fibA       = free_BigIntT(fibA);
			fibB       = free_BigIntT(fibB);
			fibTemp    = free_BigIntT(fibTemp);
		}
		
		
		//If we want to iterate each vector in a space and record their cycle 
		// and transient lengths, then repeat for a higher-powered modulus.
		else if (! strcmp(argv[1], "dynamics"))
		{
			int  lowpower = 1; //Holds the higher power of the modulus to check
			int  highpower = 2;
			int  highmod = 1;   //Holds the modulus raised to the correct power
			int  OGmod ;

			int* possibleCycleLengths;
			int* currVectElements;
			int* currVectElementsMult; //Holds multiples of the current vector's elements
			
			int firstIndex; //Used to help create our lookup table for cycle lengths
			
			//Holds counts for all possible cycle length pairs
			//For instance, if a vector % highmod has a cycle length
			// of 20, while % modulus it has a cycle length of 2,
			// then cycleTable[2][20] += 1;
			int** cycleTable;
			
			bool checkedAllVects = FALSE;
			
			IntMatrixTP A; //Holds our update matrix
			IntMatrixTP currVect;
			
			CycleInfoTP theCycle;
			
			//User-provided lower power
			if (argc > 2)
			{
				lowpower = (int)strtol(argv[2], &tempStr, 10);
				if (tempStr[0] != '\0')
				{
					fprintf(stderr, "Unable to read lower power from command line.\n");
					return EXIT_FAILURE;
				}
			}
			
			//User-provided higher power
			if (argc > 3)
			{
				highpower = (int)strtol(argv[3], &tempStr, 10);
				if (tempStr[0] != '\0')
				{
					fprintf(stderr, "Unable to read higher power from command line.\n");
					return EXIT_FAILURE;
				}
			}
			
			//If user provided a modulus for us
			if (argc > 4)
			{
				modulus = (int)strtol(argv[4], &tempStr, 10);
				if (tempStr[0] != '\0')
				{
					fprintf(stderr, "Unable to read modulus from command line.\n");
					return EXIT_FAILURE;
				}
			}
			
			A = read_IntMatrixT(updatefilepath);
			if (A == NULL)
			{
				fprintf(stderr, "Unable to read .matrix file at %s.\n", updatefilepath);
				return EXIT_FAILURE;
			}
			
			else if (rows(A) != cols(A))
			{
				printf("Update matrix provided is not square.\n");
				return EXIT_SUCCESS;
			}
			
			//Calculate moduli raised to the given power
			//It's the uer's responsibility to prevent these from overflowing
			OGmod = modulus;
			for (int i = 0; i < highpower; i += 1, highmod *= modulus);
			for (int i = 1; i < lowpower; i += 1, modulus *= OGmod);
			
			printf("Moduli: %d, %d\n", modulus, highmod);
			
			
			//Now that we have the matrix, calculate the cycle length
			// and use it to reason all possible cycle lengths
			// for alloting memory
			theCycle = floyd(A, A, highmod);
			printf("Matrix's multiplicative order mod %d: %d\n", highmod, omega(theCycle));
			
			printf("(lowmod, highmod)\n\n");
			
			//First number in the pointer says how many cycle lengths
			// are stored in the vector
			possibleCycleLengths = malloc(2*sizeof(int));
			possibleCycleLengths[0] = 0;
			
			for (int i = 1; i <= omega(theCycle)/2; i += 1)
				if (omega(theCycle) % i == 0)
				{
					possibleCycleLengths[0] += 1;
					possibleCycleLengths[possibleCycleLengths[0]] = omega(theCycle) / i;
					possibleCycleLengths = realloc(possibleCycleLengths, 
					                               (possibleCycleLengths[0]+2)*sizeof(int));
				}
				
			//Add 1 as a cycle length
			possibleCycleLengths[0] += 1;
			possibleCycleLengths[possibleCycleLengths[0]] = 1;
			
			free_CycleInfoT(theCycle);
			
			/*printf("Possible vector cycle lengths: ");
			for (int i = 1; i <= possibleCycleLengths[0]; i += 1)
				printf("%d ", possibleCycleLengths[i]);
			printf("\n\n"); */
			
			//Now, construct our cycleTable to hold cycle counts
			//Only allocating half the table since that's all we need
			cycleTable = malloc(possibleCycleLengths[0]*sizeof(int*));
			for (int i = 0; i < possibleCycleLengths[0]; i += 1)
				cycleTable[i] = calloc(i+1, sizeof(int));
			
			
			currVectElements = calloc(rows(A), sizeof(int));
			currVectElementsMult = malloc(rows(A)*sizeof(int));
			currVect = new_IntMatrixT(rows(A), 1);
			
			//Iterate until we've tested every vector
			while (!checkedAllVects)
			{
				//printf("---\n");
				set_column(currVect, currVectElements);
				theCycle = floyd(A, currVect, modulus);
				
				//Print out vector in an easier-to-look-at way for stdout
				/*printf("<");
				for (int i = 0; i < rows(A)-1; i += 1)
					printf("%d ", currVectElements[i]);
				printf("%d>\t:\t", currVectElements[rows(A)-1]);
				printf("w = %d,\tt = %d\n", omega(theCycle), tau(theCycle));
				printf("---\n");  */
				
				//Finding the first index for our table
				firstIndex = omega(theCycle);
				
				for (int i = 1; i <= possibleCycleLengths[0]; i += 1)
					if (firstIndex == possibleCycleLengths[i])
					{
						firstIndex = i-1;
						break;
					}
				
				free_CycleInfoT(theCycle);
				
				//Now that we've looked at the base vector's stats, let's go to % highmod
				for (int i = 0; i < rows(A); i += 1) 
					currVectElementsMult[i] = currVectElements[i];
				
				//Check all vectors which differ from our current vector by
				// multiples of our modulus
				while (!checkedAllVects)
				{
					set_column(currVect, currVectElementsMult);
					theCycle = floyd(A, currVect, highmod);
					
					//Print out vector in an easier-to-look-at way for stdout
					/*printf("<");
					for (int i = 0; i < rows(A)-1; i += 1)
						printf("%d ", currVectElementsMult[i]);
					printf("%d>\t:\t", currVectElementsMult[rows(A)-1]);
					printf("w = %d,\tt = %d\n", omega(theCycle), tau(theCycle)); */
					
					//Adding a count to our cycle table
					for (int i = 1; i <= possibleCycleLengths[0]; i += 1)
						if (omega(theCycle) == possibleCycleLengths[i])
						{
							if (i-1 > firstIndex)
								printf("This message appearing means a vector's cycle was added to an element outside the printed table.\n");

							cycleTable[firstIndex][i-1] += 1;
							break;
						}
					
					free_CycleInfoT(theCycle);
					
					//Now, iterate to the next vector which reduces to currVectElements
					checkedAllVects = TRUE;
					for (int i = rows(A)-1; i >= 0; i -= 1)
					{
						if (currVectElementsMult[i] != highmod - modulus + currVectElements[i])
						{
							checkedAllVects = FALSE;
							currVectElementsMult[i] += modulus;
							i = -1;
						}
						
						//Carry to next vect component
						else
							currVectElementsMult[i] = currVectElements[i];
					}
				}
				
				//Increment through all possible vectors (% modulus)
				checkedAllVects = TRUE;
				for (int i = rows(A)-1; i >= 0; i -= 1)
				{
					if (currVectElements[i] != modulus - 1)
					{
						checkedAllVects = FALSE;
						currVectElements[i] += 1;
						i = -1;
					}
					
					//Carry to next vector component
					else
						currVectElements[i] = 0;
				}
			}
			
			//Now, we print all the numbers we've collected
			printf("Cycle pairs:\n");
			for (int i = 0; i < possibleCycleLengths[0]; i += 1)
				for (int j = 0; j <= i; j += 1)
					//Only print out info is vectors have that pair of cycle lengths
					if (cycleTable[i][j] != 0)
					{
						printf("(%d, %d): %d\n", 
						possibleCycleLengths[i+1],
						possibleCycleLengths[j+1],
						cycleTable[i][j]);
					}

			A        = free_IntMatrixT(A);
			currVect = free_IntMatrixT(currVect);
			
			//Don't need to free theCycle since it gets freed above
			theCycle = NULL;
			
			for (int i = 0; i < possibleCycleLengths[0]; i += 1)
			{
				FREE(cycleTable[i]);
			}
			FREE(cycleTable);
			FREE(possibleCycleLengths);
			
			FREE(currVectElements);
			FREE(currVectElementsMult);
		}
	}
	
	else
	{
		printf(ANSI_COLOR_GREEN "LINCELLAUT by Zach Strong.\n" ANSI_COLOR_RESET);
		printf("Usage: lincellaut <tool> [options]\n\n");
		printf("Tools:\n");
		
		printf(" - " ANSI_COLOR_YELLOW "iterate " ANSI_COLOR_CYAN "[iterations]" ANSI_COLOR_RESET \
		": Iterate the initial matrix by the update matrix a given number of times.\n\n");
		
		printf(" - " ANSI_COLOR_YELLOW "inverse " ANSI_COLOR_CYAN "[modulus]" ANSI_COLOR_RESET \
		": Find the inverse of an update matrix under some modulus.\n\n");
		
		printf(" - " ANSI_COLOR_YELLOW "det " ANSI_COLOR_CYAN "[modulus]" ANSI_COLOR_RESET \
		": Find the determinant of an update matrix under some modulus.\n\n");
		
		printf(" - " ANSI_COLOR_YELLOW "chara " ANSI_COLOR_CYAN "[modulus]" ANSI_COLOR_RESET \
		": Find the characteristic equation of an update matrix under some modulus.\n\n");
		
		printf(" - " ANSI_COLOR_YELLOW "floyd " ANSI_COLOR_CYAN "[modulus]" ANSI_COLOR_RESET \
		": Use Floyd's Cycle Detection Algorithm to find out specific details about the given LCA.\n\n");
		
		printf(" - " ANSI_COLOR_YELLOW "bigfloyd " ANSI_COLOR_CYAN "[modulus]" ANSI_COLOR_RESET \
		": Same as floyd, but uses BigIntMatrixT and BigIntT structs instead of IntMatrixT and ints.\n\n");
		
		printf(" - " ANSI_COLOR_YELLOW "rots " ANSI_COLOR_CYAN "[modulus]" ANSI_COLOR_RESET \
		": Finds and outputs some basic rotation matrices for the given modulus.\n\n");
		
		printf(" - " ANSI_COLOR_YELLOW "cycmatsearch " ANSI_COLOR_CYAN "size maxmod cycles..." ANSI_COLOR_RESET \
		": Searches for a matrix whose column vectors cycle with cycle lengths less than the matrix itself.\n\n");
		
		printf(" - " ANSI_COLOR_YELLOW "fibcycle" ANSI_COLOR_CYAN " [modulus]" ANSI_COLOR_RESET \
		": Generate the Fibonacci cycle that contains the initial vector.\n\n");
		
		printf(" - " ANSI_COLOR_YELLOW "fibcyclelens" ANSI_COLOR_CYAN " [modulus]" ANSI_COLOR_RESET \
		": Calculate all possible Fibonacci cycle lengths.\n\n");
		
		printf(" - " ANSI_COLOR_YELLOW "fibmultsearch " ANSI_COLOR_CYAN "[bound]" ANSI_COLOR_RESET \
		": Searches the Fibonacci numbers, checking whether a " \
		"multiple of each number up to the bound appears before a multiple of a power of the number.\n\n");
		
		printf(" - " ANSI_COLOR_YELLOW "dynamics " ANSI_COLOR_CYAN "[power1] [power2] [modulus]" ANSI_COLOR_RESET \
		": Iterates every vector in a space, recording their transient lengths and cycle lengths. " \
		"It then computes the same numbers for a higher-powered modulus.\n\n");
		
		printf("For a more complete description of LINCELLAUT's usage, " \
		"refer to the included documentation.\n");
		
		//Testing find_factors()
		/*
		int factorArr[1]  = {0};
		int targetArr[1]  = {1};
		int carryArr[1]   = {1};
		int modulusArr[1] = {49};
		
		BigIntTP target  = new_BigIntT(targetArr, 1);
		BigIntTP factor  = new_BigIntT(factorArr, 1);
		BigIntTP carry   = new_BigIntT(carryArr, 1);
		BigIntTP modulus = new_BigIntT(modulusArr, 1);
		
		BigIntTP* possibilities = malloc(sizeof(BigIntTP));
		possibilities[0] = empty_BigIntT(1);
		int possibilityCount = 0;
		
		possibilityCount = find_factors(target, factor, carry, modulus, &possibilities);
		printf("Possibilities found: %d\n", possibilityCount);
		for (int i = 0; i < possibilityCount; i += 1)
		{
			printi(possibilities[i]);
			printf(" ");
		}
		printf("\n");
		
		target  = free_BigIntT(target);
		factor  = free_BigIntT(factor);
		carry   = free_BigIntT(carry);
		modulus = free_BigIntT(modulus);
		
		for (int i = 0; i < possibilityCount+1; i += 1)
			possibilities[i] = free_BigIntT(possibilities[i]);
		FREE(possibilities);
		*/
	}
	
	//Freeing memory
	FREE(updatefilepath);
	FREE(initialfilepath);
	FREE(iterfilepath);
	FREE(bigintmodstring);
	
	return EXIT_SUCCESS;
}