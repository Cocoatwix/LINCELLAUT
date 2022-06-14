
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

#include "headers/bigint.h" //Arbitrary precision

#include "headers/linalg.h" 
#include "headers/cycles.h"  //Allows us to use Floyd's Algorithm
#include "headers/modular.h" //Modular square roots and inverses
#include "headers/fibonacci.h"
#include "headers/factors.h" //For LCM()

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
			Cycles 2, 3, 5 with 3x3 matrices up to and including mod 5
			*/
			
			int oneArr[] = {1};
			int start[] = {2};
			
			//printf("Currently, the first modulus checked is not 2 for testing purposes.\n");
			
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
			
			/*
			//TESTING A SPECIFIC CASE. DELETE THIS LATER
			currMatElements[0][0] = free_BigIntT(currMatElements[0][0]);
			currMatElements[1][1] = free_BigIntT(currMatElements[1][1]);
			
			int scale1[] = {2};
			int scale2[] = {3};
			
			currMatElements[0][0] = new_BigIntT(scale1, 1);
			currMatElements[1][1] = new_BigIntT(scale2, 1);
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
								copy_BigIntT(zero, currMatElements[i][j]);
							
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
	}
	
	else
	{
		printf(ANSI_COLOR_GREEN "LINCELLAUT by Zach Strong.\n" ANSI_COLOR_RESET);
		printf("Usage: lincellaut <tool> [options]\n\n");
		printf("Tools:\n");
		
		printf(" - " ANSI_COLOR_YELLOW "iterate " ANSI_COLOR_CYAN "[iterations]" ANSI_COLOR_RESET \
		": Iterate the initial matrix by the update matrix a given number of times.\n");
		printf("   - " ANSI_COLOR_CYAN "iterations" ANSI_COLOR_RESET \
		": Overrides the number of iterations provided in the .config file.\n\n");
		
		printf(" - " ANSI_COLOR_YELLOW "inverse " ANSI_COLOR_CYAN "[modulus]" ANSI_COLOR_RESET \
		": Find the inverse of an update matrix under some modulus.\n");
		printf("   - " ANSI_COLOR_CYAN "modulus" ANSI_COLOR_RESET \
		": Overrides the modulus provided in the .config file.\n\n");
		
		printf(" - " ANSI_COLOR_YELLOW "floyd " ANSI_COLOR_CYAN "[modulus]" ANSI_COLOR_RESET \
		": Use Floyd's Cycle Detection Algorithm to find out specific details about the given LCA.\n");
		printf("   - " ANSI_COLOR_CYAN "modulus" ANSI_COLOR_RESET \
		": Overrides the modulus provided in the .config file.\n\n");
		
		printf(" - " ANSI_COLOR_YELLOW "bigfloyd " ANSI_COLOR_CYAN "[modulus]" ANSI_COLOR_RESET \
		": Same as floyd, but uses BigIntMatrixT and BigIntT structs instead of IntMatrixT and ints.\n");
		printf("   - " ANSI_COLOR_CYAN "modulus" ANSI_COLOR_RESET \
		": Overrides the modulus provided in the .config file.\n\n");
		
		printf(" - " ANSI_COLOR_YELLOW "rots " ANSI_COLOR_CYAN "[modulus]" ANSI_COLOR_RESET \
		": Finds and outputs some basic rotation matrices for the given modulus.\n");
		printf("   - " ANSI_COLOR_CYAN "modulus" ANSI_COLOR_RESET \
		": Overrides the modulus provided in the .config file.\n\n");
		
		printf(" - " ANSI_COLOR_YELLOW "cycmatsearch " ANSI_COLOR_CYAN "size maxmod cycles..." ANSI_COLOR_RESET \
		": Searches for a matrix whose column vectors cycle with cycle lengths less than the matrix itself.\n");
		printf("   - " ANSI_COLOR_CYAN "size" ANSI_COLOR_RESET \
		": Tells what size matrix to use.\n");
		printf("   - " ANSI_COLOR_CYAN "maxmod" ANSI_COLOR_RESET \
		": Tells which modulus to stop searching at.\n");
		printf("   - " ANSI_COLOR_CYAN "cycles..." ANSI_COLOR_RESET \
		": A list of cycle lengths for their respective vectors.\n\n");
		
		printf(" - " ANSI_COLOR_YELLOW "fibcycle" ANSI_COLOR_CYAN " [modulus]" ANSI_COLOR_RESET \
		": Generate the Fibonacci cycle that contains the initial vector.\n");
		printf("   - " ANSI_COLOR_CYAN "modulus" ANSI_COLOR_RESET \
		": Overrides the modulus provided in the .config file.\n\n");
		
		printf(" - " ANSI_COLOR_YELLOW "fibcyclelens" ANSI_COLOR_CYAN " [modulus]" ANSI_COLOR_RESET \
		": Calculate all possible Fibonacci cycle lengths.\n");
		printf("   - " ANSI_COLOR_CYAN "modulus" ANSI_COLOR_RESET \
		": Overrides the modulus provided in the .config file.\n\n");
		
		printf(" - " ANSI_COLOR_YELLOW "fibmultsearch " ANSI_COLOR_CYAN "[bound]" ANSI_COLOR_RESET \
		": Searches the Fibonacci numbers, checking whether a " \
		"multiple of each number up to the bound appears before a multiple of a power of the number.\n");
		printf("   - " ANSI_COLOR_CYAN "bound" ANSI_COLOR_RESET \
		": Override the default upper bound of 100.\n\n");
		
		printf("For a more complete description of LINCELLAUT's usage, " \
		"refer to the included documentation.\n");
		
		/*
		printf("\n\n\n");
		
		int a[] = {80162984, 97337589, 78534648};
		int b[] = {47565795, 58508592, 1796146};
		
		BigIntTP A = new_BigIntT(a, 3);
		BigIntTP B = new_BigIntT(b, 3);
		BigIntTP C = empty_BigIntT(1);
		
		multiply_BigIntT(A, B, C);
		printi(A);
		printf(" * ");
		printi(B);
		printf(" = ");
		printi(C);
		printf("\n");
		
		A = free_BigIntT(A);
		B = free_BigIntT(B);
		C = free_BigIntT(C);
		
		
		int a11[] = {25807962, 95735574};
		BigIntTP A11 = new_BigIntT(a11, 2);
		int a12[] = {34356792, 1510014};
		BigIntTP A12 = new_BigIntT(a12, 2);
		int a21[] = {88156509, 25413232};
		BigIntTP A21 = new_BigIntT(a21, 2);
		int a22[] = {87029814, 76108034};
		BigIntTP A22 = new_BigIntT(a22, 2);
		int b11[] = {9224407, 4965146};
		BigIntTP B11 = new_BigIntT(b11, 2);
		int b12[] = {18880572, 54166990};
		BigIntTP B12 = new_BigIntT(b12, 2);
		int b21[] = {75487505, 84273811};
		BigIntTP B21 = new_BigIntT(b21, 2);
		int b22[] = {54462723, 10778472};
		BigIntTP B22 = new_BigIntT(b22, 2);
	
		BigIntTP* Amat1 = malloc(2*sizeof(BigIntTP));
		BigIntTP* Amat2 = malloc(2*sizeof(BigIntTP));
		BigIntTP* Bmat1 = malloc(2*sizeof(BigIntTP));
		BigIntTP* Bmat2 = malloc(2*sizeof(BigIntTP));
		
		Amat1[0] = A11;
		Amat1[1] = A12;
		Amat2[0] = A21;
		Amat2[1] = A22;
		
		Bmat1[0] = B11;
		Bmat1[1] = B12;
		Bmat2[0] = B21;
		Bmat2[1] = B22;
		
		BigIntTP** Amat = malloc(2*sizeof(BigIntTP*));
		BigIntTP** Bmat = malloc(2*sizeof(BigIntTP*));
		
		Amat[0] = Amat1;
		Amat[1] = Amat2;
		Bmat[0] = Bmat1;
		Bmat[1] = Bmat2;

		BigIntMatrixTP A = new_BigIntMatrixT(2, 2);
		BigIntMatrixTP B = new_BigIntMatrixT(2, 2);
		BigIntMatrixTP C = new_BigIntMatrixT(2, 2);
		set_big_matrix(A, Amat);
		set_big_matrix(B, Bmat);
		
		printf("A:\n");
		printbm(A);
		printf("B:\n");
		printbm(B);
		
		printf("AB = \n");
		big_mat_mul(A, B, C);
		printbm(C);
		
		A11 = free_BigIntT(A11);
		A12 = free_BigIntT(A12);
		A21 = free_BigIntT(A21);
		A22 = free_BigIntT(A22);
		B11 = free_BigIntT(B11);
		B12 = free_BigIntT(B12);
		B21 = free_BigIntT(B21);
		B22 = free_BigIntT(B22);
		
		FREE(Amat1);
		FREE(Amat2);
		FREE(Bmat1);
		FREE(Bmat2);
		
		FREE(Amat);
		FREE(Bmat);
		
		A = free_BigIntMatrixT(A);
		B = free_BigIntMatrixT(B);
		C = free_BigIntMatrixT(C); */
	}
	
	/*
	//Now testing our ability to create eigenvectors
	printf("Matrix:\n");
	printm(F, TRUE);
	printf("Using eigenvalue %d:\n", values[1]);
	IntMatrixTP E = eigenvector(F, values[1], modulus);
	//printm(E, TRUE);
	
	FREE(values);
	//E = free_IntMatrixT(E); */
	
	//Freeing memory
	FREE(updatefilepath);
	FREE(initialfilepath);
	FREE(iterfilepath);
	FREE(bigintmodstring);
	
	return EXIT_SUCCESS;
}