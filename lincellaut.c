
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

//Debug macros; can be removed if not needed
#define opp(a, b, c, o) printf("("); printp(a); printf(")"); printf(o); \
											  printf("("); printp(b); printf(")"); printf(" == "); \
											  printp(c); printf("\n")
#define opi(a, b, c, o) printf("("); printi(a); printf(")"); printf(o); \
											  printf("("); printi(b); printf(")"); printf(" == "); \
											  printi(c); printf("\n")

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"

/** Frees all "global" variables we allocated memory for. */
#define FREE_VARIABLES FREE(updatefilepath); \
											 FREE(initialfilepath); \
											 FREE(resumefilepath); \
											 FREE(iterfilepath); \
											 FREE(bigintmodstring); \
											 FREE(resumemodstring)

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
	
	int iterations;
	int modulus;
	
	//For holding str representation of mod, used for BigIntT number (if needed)
	char* bigintmodstring = malloc(MAXSTRLEN*sizeof(char));
	char* resumemodstring = malloc(MAXSTRLEN*sizeof(char));
	
	char* updatefilepath  = malloc(MAXSTRLEN*sizeof(char));
	char* initialfilepath = malloc(MAXSTRLEN*sizeof(char));
	char* resumefilepath  = malloc(MAXSTRLEN*sizeof(char));
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
		
		else if (! strcmp(systemData, "resumeMat"))
		{
			if (fscanf(system, "%s", resumefilepath) != 1)
			{
				fprintf(stderr, "Unable to read resume matrix path from config file.\n");
				return EXIT_FAILURE;
			}
		}
		
		else if (! strcmp(systemData, "resumeMod"))
		{
			if (fscanf(system, "%101s", resumemodstring) != 1)
			{
				fprintf(stderr, "Unable to read resume modulus from config file.\n");
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
		FREE_VARIABLES;
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
					FREE_VARIABLES;
					return EXIT_FAILURE;
				}
			}
			
			IntMatrixTP F;
			IntMatrixTP F_2;
			IntMatrixTP F_result = NULL;
			
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
			
			if (cols(F) != rows(F_2))
			{
				fprintf(stderr, "Given matrices cannot be multiplied.\n");
				return EXIT_FAILURE;
			}
			
			//Iterate F a few times
			//The minus 1 is for easier conversion between ORBITVIS results
			printf("Iterations: %d\n", iterations);
			if (iterations == 1)
			{
				F_result = new_IntMatrixT(rows(F), cols(F_2));
				mat_mul(F, F_2, F_result);
				modm(F_result, modulus);
			}
			
			else if (iterations > 0)
				F_result = iterate(F, F_2, modulus, iterations);
			
			//Prevents the matrix from being printed when
			// an inverse doesn't exist and the iterations
			// provided was negative
			if (iterations > 0)
				printm(F_result);
			
			else if (iterations == 0)
				printm(F_2);

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
			int theDet;
			
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
			theDet = det(A);
			
			if (theDet < 0)
				printf("Determinant: %d\n", ((theDet % modulus) + modulus) % modulus);
			else
				printf("Determinant: %d\n", (theDet % modulus));
			
			A = free_IntMatrixT(A);
		}
		
		
		//Find the characteristic equation of the update matrix
		else if (! strcmp(argv[1], "chara"))
		{
			int oneArr[1] = {1};
			
			BigIntTP bigMod = NULL;
			BigIntTP counter = NULL; //For displaying the factors we get
			BigIntTP one = NULL;
			BigIntTP temp = NULL;
			BigIntTP numberOfTerms = NULL;
			
			int smallCounter = 0;
			
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
			printf("Matrix:\n");
			printbm(bigMatrix);
			printf("\nModulus: ");
			printi(bigMod);
			printf("\n");
			bigEqn = chara_poly(bigMatrix, bigMod);
			printf("\nCharacteristic equation: ");
			printp(bigEqn);
			printf("\n\n");
			
			printf("Factored characteristic equation:\n");
			bigEqnFactors = factor_BigPolyT(bigEqn, bigMod);
			
			one  = new_BigIntT(oneArr, 1);
			temp = empty_BigIntT(1);
			numberOfTerms = empty_BigIntT(1);
			copy_BigIntT(constant(bigEqnFactors[0]), numberOfTerms);
			
			for (counter = empty_BigIntT(1); 
			     compare_BigIntT(numberOfTerms, counter) >= 0; 
					 add_BigIntT(one, counter, temp), copy_BigIntT(temp, counter), smallCounter += 1)
			{
				if (smallCounter != 0)
				{
					printf("(");
					printp(bigEqnFactors[smallCounter]);
					printf(")");
				}
				bigEqnFactors[smallCounter] = free_BigPolyT(bigEqnFactors[smallCounter]);
			}
			FREE(bigEqnFactors);
			printf("\n");
			
			bigMod  = free_BigIntT(bigMod);
			counter = free_BigIntT(counter);
			temp    = free_BigIntT(temp);
			one     = free_BigIntT(one);
			
			numberOfTerms = free_BigIntT(numberOfTerms);
			
			bigMatrix = free_BigIntMatrixT(bigMatrix);
			bigEqn = free_BigPolyT(bigEqn);
		}
		
		
		//If we want to see how many vectors are in an IntMatrix's core
		else if (! strcmp(argv[1], "core"))
		{
			BigIntMatrixTP A;
			BigIntMatrixTP I = NULL;
			BigIntMatrixTP vectToTest = NULL;
			BigIntMatrixTP zeroVector;
			BigIntMatrixTP tempVector;
			
			BigIntTP bigMod;
			BigIntTP one;
			BigIntTP zero;
			BigIntTP temp;
			BigIntTP kernelCount;
			BigIntTP orderOfSpace;
			
			CycleInfoTP theCycle = NULL;
			
			BigIntTP** entriesToUse;
			
			int oneArr[1] = {1};
			
			bool checkedAllVects = FALSE;
			bool nonSquareMatrix = FALSE;
			
			//If the user provided a modulus on the CLI
			if (argc > 2)
			{
				if (strtoBIT(argv[2], &bigMod) == 0)
				{
					fprintf(stderr, "Invalid modulus passed at command line.\n");
					return EXIT_FAILURE;
				}
			}
			
			else
			{
				if (strtoBIT(bigintmodstring, &bigMod) == 0)
				{
					fprintf(stderr, "Invalid modulus passed at command line.\n");
					return EXIT_FAILURE;
				}
			}
			
			A = read_BigIntMatrixT(updatefilepath);
			if (A == NULL)
			{
				fprintf(stderr, "Unable to read matrix file at %s.\n", updatefilepath);
				return EXIT_FAILURE;
			}
			
			entriesToUse = malloc(big_cols(A)*sizeof(BigIntTP*));
			for (int i = 0; i < big_cols(A); i += 1)
			{
				entriesToUse[i] = malloc(sizeof(BigIntTP));
				entriesToUse[i][0] = empty_BigIntT(1);
			}
			
			one  = new_BigIntT(oneArr, 1);
			zero = empty_BigIntT(1);
			temp = empty_BigIntT(1);
			kernelCount = empty_BigIntT(1);
			
			vectToTest = new_BigIntMatrixT(big_cols(A), 1);
			zeroVector = new_BigIntMatrixT(big_rows(A), 1);
			tempVector = new_BigIntMatrixT(big_rows(A), 1);
			
			printf("Matrix:\n");
			printbm(A);
			printf("Modulus: ");
			printi(bigMod);
			printf("\n");
			
			//Only use floyd if we actually can
			if ((big_rows(A)) == big_cols(A))
			{
				I = identity_BigIntMatrixT(big_rows(A));
				big_floyd(A, I, bigMod, &theCycle);
				
				printf("Rep:\n");
				printbm(rep(theCycle));
			}
			
			//NOTE: this tool is NOT a good way to gauge how many unique vectors
			// can be made with a given basis, since duplicates can appear without 
			// zero vectors ever appearing.
			else
			{
				nonSquareMatrix = TRUE;
				printf("ker(A):\n");
			}
			
			//Loop until all vectors have been checked
			//This could be a lot more efficient if I understood how best to
			// manipulate how the kernel vectors divide the plane
			while (((nonSquareMatrix) || (tau(theCycle) != 0)) && (!checkedAllVects))
			{
				set_big_matrix(vectToTest, entriesToUse);
				
				if (nonSquareMatrix)
					big_mat_mul(A, vectToTest, tempVector);
				else
					big_mat_mul(rep(theCycle), vectToTest, tempVector);
				
				modbm(tempVector, bigMod);
				
				//If we found a vector that maps to the origin
				if (compare_BigIntMatrixT(tempVector, zeroVector))
				{
					add_BigIntT(kernelCount, one, temp);
					copy_BigIntT(temp, kernelCount);
					
					if (nonSquareMatrix)
					{
						printbm(vectToTest);
						printf("\n");
					}
				}
				
				//Increment to next vector
				checkedAllVects = TRUE;
				for (int i = big_cols(A)-1; i >= 0; i -= 1)
				{
					add_BigIntT(entriesToUse[i][0], one, temp);
					
					if (compare_BigIntT(temp, bigMod) == 0)
						copy_BigIntT(zero, entriesToUse[i][0]);
					
					else
					{
						copy_BigIntT(temp, entriesToUse[i][0]);
						checkedAllVects = FALSE;
						break;
					}
				}
			}
			
			printf("|ker(A)| = ");
			printi(kernelCount);
			
			//It only makes sense to talk about the core of a square matrix
			if (! nonSquareMatrix)
			{
				printf("\n|core(A)| = ");
				
				//Calculate number of vectors in our space,
				// then compact it
				orderOfSpace = new_BigIntT(oneArr, 1);
				for (int i = 0; i < big_rows(A); i += 1)
				{
					multiply_BigIntT(orderOfSpace, bigMod, temp);
					copy_BigIntT(temp, orderOfSpace);
				}
				
				if (compare_BigIntT(zero, kernelCount) != 0)
					divide_BigIntT(orderOfSpace, kernelCount, temp);
				else
					copy_BigIntT(orderOfSpace, temp);
				
				printi(temp);
			}
			printf("\n");
			
		
			for (int i = 0; i < big_cols(A); i += 1)
			{
				entriesToUse[i][0] = free_BigIntT(entriesToUse[i][0]);
				FREE(entriesToUse[i]);
			}
			FREE(entriesToUse);
			
			A          = free_BigIntMatrixT(A);
			I          = free_BigIntMatrixT(I);
			vectToTest = free_BigIntMatrixT(vectToTest);
			zeroVector = free_BigIntMatrixT(zeroVector);
			tempVector = free_BigIntMatrixT(tempVector);
			
			bigMod       = free_BigIntT(bigMod);
			one          = free_BigIntT(one);
			zero         = free_BigIntT(zero);
			temp         = free_BigIntT(temp);
			kernelCount  = free_BigIntT(kernelCount);
			orderOfSpace = free_BigIntT(orderOfSpace);
			
			theCycle = free_CycleInfoT(theCycle);
		}
		
		
		//If we want to use Floyd's Cycle Detection Algorithm
		else if (! strcmp(argv[1], "floyd"))
		{
			//If the user provided a custom modulus
			if (argc > 2)
			{
				modulus = (int)strtol(argv[2], &tempStr, 10);
				if (tempStr[0] != '\0')
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
			
			//Making sure dimensions of given matrices are correct
			else if (rows(initial) != cols(update) ||
			        (rows(initial) != rows(update)))
			{
				fprintf(stderr, "Given matrices have inappropriate dimensions for iterating.\n");
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
			
			else if ((big_rows(initial) != big_cols(update)) ||
			         (big_rows(initial) != big_rows(update)))
		  {
			  fprintf(stderr, "Given matrices have inappropriate dimensions for iterating.\n");
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
			if (argc < 6)
			{
				printf(ANSI_COLOR_YELLOW "cycmatsearch " ANSI_COLOR_CYAN "resume maxmod cycles..." ANSI_COLOR_RESET \
				": Searches for a matrix whose column vectors cycle with cycle lengths less than the matrix itself.\n");
				printf(" - " ANSI_COLOR_CYAN "resume" ANSI_COLOR_RESET \
				": Says whether to resume computation at a particular matrix.\n");
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
			cycmatsearch 2 28 3 5
			cycmatsearch 2 30 3 4
			cycmatsearch 3 6 2 3 5
			cycmatsearch 3 6 2 2 3
			
			cycmatsearch 4 3 2 2 3 3 . . .
			*/
			
			int oneArr[] = {1};
			int start[] = {2}; //Modulus to start with
			
			int* colVectCycles; //Holds the different cycle lengths 
			
			size = (int)strtol(argv[3], &tempStr, 10);
			
			if (tempStr[0] != '\0')
			{
				fprintf(stderr, "Unable to read matrix size from command line.\n");
				return EXIT_FAILURE;
			}
			
			//Initialise our matrix elements
			currMatElements = malloc(size*sizeof(BigIntTP*));
			for (i = 0; i < size; i += 1)
			{
				currMatElements[i] = malloc(size*sizeof(BigIntTP));
				for (j = 0; j < size; j += 1)
					currMatElements[i][j] = empty_BigIntT(1);
			}
			
			//If we need to resume computation from a specific matrix
			if (! strcmp(argv[2], "TRUE"))
			{
				currMat = read_BigIntMatrixT(resumefilepath);
				if (strtoBIT(bigintmodstring, &currMod) == 0)
				{
					fprintf(stderr, "Unable to read resume modulus from config file.\n");
					FREE_VARIABLES;
					return EXIT_FAILURE;
				}
				
				if (currMat == NULL)
				{
					fprintf(stderr, "Unable to read resume matrix from config file.\n");
					FREE_VARIABLES;
					return EXIT_FAILURE;
				}
				
				//Now, set the values for currMatElements
				for (i = 0; i < size; i += 1)
					for (j = 0; j < size; j += 1)
						copy_BigIntT(big_element(currMat, i, j), currMatElements[i][j]);
			}
			
			else
			{
				//Start looking for matrices at mod 2 from the zero matrix
				currMod = new_BigIntT(start, 1);
				currMat = new_BigIntMatrixT(size, size);
			}
			
			//If we couldn't read the modulus
			if (strtoBIT(argv[4], &maxMod) == 0)
			{
				fprintf(stderr, "Unable to read modulus from command line.\n");
				return EXIT_FAILURE;
			}
			
			//If the user didn't provide enough cycles for the given matrix size
			if (argc < 5 + size)
			{
				printf("Too few cycles passed for given matrix size.\n");
				return EXIT_FAILURE;
			}
			
			//Get column vector cycles
			colVectCycles = malloc(size*sizeof(int));
			for (i = 0; i < size; i += 1)
			{
				colVectCycles[i] = (int)strtol(argv[5+i], &tempStr, 10);
				if (tempStr[0] != '\0')
				{
					fprintf(stderr, "Unable to read cycle lengths from command line.\n");
					return EXIT_FAILURE;
				}
			}
			
			one     = new_BigIntT(oneArr, 1);
			zero    = empty_BigIntT(1);
			temp    = empty_BigIntT(1);
			
			/*
			//INITIALISING SPECIFIC NUMBERS SO THAT WE CAN STOP EXECUTION
			// WHEN WE GET TO A SPECIFIC MATRIX
			BigIntTP three = new_BigIntT(threeArr, 1);
			BigIntTP two   = new_BigIntT(twoArr, 1);
			BigIntTP five  = new_BigIntT(fiveArr, 1);

			BigIntTP stopMat[3][3] = {{one,   three, five},
																{two,   zero,  one},
																{three, three, one}};
																*/
			
			
			//Find the LCM of our cycles
			for (i = 0; i < size; i += 1)
				cycleLCM = LCM(cycleLCM, colVectCycles[i]);
			
			zeroMat  = new_BigIntMatrixT(size, size);
			tempMat  = new_BigIntMatrixT(size, size);
			tempMat2 = new_BigIntMatrixT(size, size);
			
			//Initialising a reusable CycleInfoTP to save time on
			// memory allocations.
			theCycle = new_CycleInfoT();
			
			//Forming name for text file output.
			textOutputName = malloc(MAXSTRLEN*sizeof(char));
			textOutputName[0] = '\0';
			strcat(textOutputName, argv[1]); //cycmatsearch
			strcat(textOutputName, " ");
			strcat(textOutputName, argv[3]); //size
			strcat(textOutputName, " ");
			strcat(textOutputName, argv[4]); //maxmod
			strcat(textOutputName, " ");
			for (i = 0; i < size-1; i += 1)
			{
				strcat(textOutputName, argv[5+i]);
				strcat(textOutputName, " ");
			}
			strcat(textOutputName, argv[size+4]);
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
				
				fprintf(textOutput, "~~~Mod ");
				fprinti(textOutput, currMod);
				fprintf(textOutput, "~~~\n");

				//Search through all matrices under the current modulus
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
						
						fprintbm_nopad(textOutput, currMat);
						fprintf(textOutput, "\n");
						
						if (fclose(textOutput) == EOF)
						{
							fprintf(stderr, "Error saving file %s.\n", textOutputName);
							checkedAllMatrices = FALSE;
							break;
						}
						
						textOutput = fopen(textOutputName, "a");
						if (textOutput == NULL)
						{
							fprintf(stderr, "Error reopening file %s after saving.\n", textOutputName);
							checkedAllMatrices = FALSE;
							break;
						}
					}
					
					
					//Iterate to the next matrix
					checkedAllMatrices = TRUE;
					
					/*
					//CHECK TO SEE IF OUR MATRIX IS THE ONE WE WANT TO STOP AT
					for (int i = 0; i < size; i += 1)
					{
						for (int j = 0; j < size; j += 1)
						{
							if (compare_BigIntT(currMatElements[i][j], stopMat[i][j]) != 0)
							{
								checkedAllMatrices = FALSE;
								i = size;
								j = size;
							}
						}
					}
					*/

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
		
		
		//If we want to compute some "cycle converting matrices"
		else if (!strcmp(argv[1], "cycconvmat"))
		{
			BigIntTP bigMod;
			BigIntTP from;
			BigIntTP to;

			BigIntMatrixTP sumMat;
			BigIntMatrixTP A;
			
			//Check to see if user provided a custom modulus
			if (argc > 4)
			{
				if (strtoBIT(argv[4], &bigMod) != 1)
				{
					fprintf(stderr, "Unable to read modulus from command line.\n");
					FREE_VARIABLES;
					return EXIT_FAILURE;
				}
			}
			
			else
			{
				if (strtoBIT(bigintmodstring, &bigMod) != 1)
				{
					fprintf(stderr, "Unable to read modulus from config file.\n");
					FREE_VARIABLES;
					return EXIT_FAILURE;
				}
			}
			
			if (strtoBIT(argv[2], &from) != 1)
			{
				fprintf(stderr, "Unable to read starting cycle length from command line.\n");
				FREE_VARIABLES;
				return EXIT_FAILURE;
			}
			
			if (strtoBIT(argv[3], &to) != 1)
			{
				fprintf(stderr, "Unable to read conversion cycle length from command line.\n");
				FREE_VARIABLES;
				return EXIT_FAILURE;
			}
			
			A = read_BigIntMatrixT(updatefilepath);
			if (A == NULL)
			{
				fprintf(stderr, "Unable to read update matrix provided in config file.\n");
				FREE_VARIABLES;
				return EXIT_FAILURE;
			}

			sumMat = new_BigIntMatrixT(big_rows(A), big_rows(A));
			ccm(A, sumMat, from, to, bigMod);
			
			//Now, print out our resulting "cycle converting matrix"
			printf("Matrix:\n");
			printbm(A);
			printf("Modulus: ");
			printi(bigMod);
			printf("\n\nCCM from ");
			printi(from);
			printf(" to ");
			printi(to);
			printf(":\n");
			printbm(sumMat);
			
			
			from   = free_BigIntT(from);
			to     = free_BigIntT(to);
			bigMod = free_BigIntT(bigMod);
			
			A      = free_BigIntMatrixT(A);
			sumMat = free_BigIntMatrixT(sumMat);
		}
		
		
		//If we want to look for matrices that have CCMs equal to zero
		else if (! strcmp(argv[1], "ccmzerosearch"))
		{
			BigIntTP   bigMod;
			BigIntTP   zero;
			BigIntTP   one;
			BigIntTP   temp;
			BigIntTP   temp2;
			BigIntTP*  cycleLengthFactors;
			BigIntTP** currMatElements;
			
			BigIntMatrixTP currMat;
			BigIntMatrixTP tempCCM;
			BigIntMatrixTP zeroMat;
			
			CycleInfoTP theCycle = NULL;
			
			int matSize;
			int indexCounter; //For freeing and printing
			int oneArr[1] = {1};
			
			//Holds the cycle length of each matrix we check
			int cycleLengthArray[1] = {0};
			BigIntTP bigOmega;
			
			bool checkedAllMatrices = FALSE;
			
			//Checking to see if the user provided a modulus on the command line
			if (argc > 3)
			{
				if (strtoBIT(argv[3], &bigMod) != 1)
				{
					fprintf(stderr, "Unable to read modulus on command line.\n");
					FREE_VARIABLES;
					return EXIT_FAILURE;
				}
			}
			
			else
			{
				if (strtoBIT(bigintmodstring, &bigMod) != 1)
				{
					fprintf(stderr, "Unable to read modulus from config file.\n");
					FREE_VARIABLES;
					return EXIT_FAILURE;
				}
			}
			
			matSize = (int)strtol(argv[2], &tempStr, 10);
			if (tempStr[0] != '\0')
			{
				fprintf(stderr, "Invalid matrix size provided on command line.\n");
				FREE_VARIABLES;
				return EXIT_FAILURE;
			}
			
			currMatElements = malloc(matSize*sizeof(BigIntTP*));
			for (int row = 0; row < matSize; row += 1)
			{
				currMatElements[row] = malloc(matSize*sizeof(BigIntTP));
				for (int col = 0; col < matSize; col += 1)
					currMatElements[row][col] = empty_BigIntT(1);
			}
			
			zero = empty_BigIntT(1);
			one  = new_BigIntT(oneArr, 1);
			currMat = new_BigIntMatrixT(matSize, matSize);
			zeroMat = new_BigIntMatrixT(matSize, matSize);
			tempCCM = new_BigIntMatrixT(matSize, matSize);
			
			//Loop until we're checked every matrix under the given modulus
			while (!checkedAllMatrices)
			{
				set_big_matrix(currMat, currMatElements);
				big_floyd(currMat, currMat, bigMod, &theCycle);
				
				printf(":)\n");
				
				//This program currently assumes that the cycle length of the
				// matrix will be less than MAXBUNCH
				//free_CycleInfoT(theCycle);
				//cycleLengthArray[0] = omega(theCycle);
				cycleLengthArray[0] = 123;
				bigOmega = new_BigIntT(cycleLengthArray, 1);
				cycleLengthFactors = divisors_of_BigIntT(bigOmega);
				
				temp  = empty_BigIntT(1);
				temp2 = empty_BigIntT(1);
				indexCounter = 1;
				while (compare_BigIntT(temp, cycleLengthFactors[0]) < 0)
				{
					printi(cycleLengthFactors[indexCounter]);
					printf(", ");
					
					indexCounter += 1;
					add_BigIntT(temp, one, temp2);
					copy_BigIntT(temp2, temp);
				}
				printf("\n");
				
				//Testing to ensure factorisation function works properly
				bigOmega = free_BigIntT(bigOmega);
				
				checkedAllMatrices = TRUE;
			}
			
			
			for (int row = 0; row < matSize; row += 1)
			{
				for (int col = 0; col < matSize; col += 1)
					currMatElements[row][col] = free_BigIntT(currMatElements[row][col]);
				FREE(currMatElements[row]);
			}
			FREE(currMatElements);
			
			//Freeing this matrix is a massive pain
			indexCounter = 1;
			copy_BigIntT(zero, temp);
			while (compare_BigIntT(temp, cycleLengthFactors[0]) < 0)
			{
				cycleLengthFactors[indexCounter] = free_BigIntT(cycleLengthFactors[indexCounter]);
				indexCounter += 1;
				add_BigIntT(temp, one, temp2);
				copy_BigIntT(temp2, temp);
			}
			cycleLengthFactors[0] = free_BigIntT(cycleLengthFactors[0]);
			FREE(cycleLengthFactors);
			
			currMat = free_BigIntMatrixT(currMat);
			zeroMat = free_BigIntMatrixT(zeroMat);
			tempCCM = free_BigIntMatrixT(tempCCM);
			
			theCycle = free_CycleInfoT(theCycle);
			
			bigMod = free_BigIntT(bigMod);
			zero   = free_BigIntT(zero);
			one    = free_BigIntT(one);
			temp   = free_BigIntT(temp);
			temp2  = free_BigIntT(temp2);
		}
		
		
		//If we want to step over a matrix space to see what the
		// characteristic polynomials look like over specific
		// step sizes.
		else if (! strcmp(argv[1], "charawalk"))
		{
			int counter;
			int oneArr[1] = {1};
			
			BigIntTP mod  = NULL;
			BigIntTP step = NULL; //How much to step in each direction in the matrix space
			BigIntTP temp;
			BigIntTP temp2;
			BigIntTP zero;
			BigIntTP one;
			
			//Holds the BigIntTPs describing what matrix we're
			// currently testing
			BigIntTP** currentMatrixElements;
			
			BigPolyTP  charaPoly;
			BigPolyTP* factorList = NULL;
			
			BigIntMatrixTP startingMatrix;
			BigIntMatrixTP currentMatrix;
			BigIntMatrixTP identity;
			
			CycleInfoTP theCycle;
			
			startingMatrix = read_BigIntMatrixT(updatefilepath);
			if (startingMatrix == NULL)
			{
				fprintf(stderr, "Unable to read .matrix file at %s.\n", updatefilepath);
				return EXIT_FAILURE;
			}
			
			//Getting user provided step size
			if (strtoBIT(argv[2], &step) == 0)
			{
				fprintf(stderr, "Invalid step size passed on command line.\n");
				return EXIT_FAILURE;
			}
			
			//If user provided a specific modulus on CLI
			if (argc > 3)
			{
				if (strtoBIT(argv[3], &mod) == 0)
				{
					fprintf(stderr, "Invalid modulus passed on command line.\n");
					return EXIT_FAILURE;
				}
			}
			
			else
			{
				if (strtoBIT(bigintmodstring, &mod) == 0)
				{
					fprintf(stderr, "Unable to read modulus from .config file.\n");
					return EXIT_FAILURE;
				}
			}
			
			//Initialise currentMatrix
			currentMatrixElements = malloc(big_rows(startingMatrix)*sizeof(BigIntTP*));
			for (int i = 0; i < big_rows(startingMatrix); i += 1)
			{
				currentMatrixElements[i] = malloc(big_rows(startingMatrix)*sizeof(BigIntTP));
				for (int j = 0; j < big_rows(startingMatrix); j += 1)
				{
					currentMatrixElements[i][j] = empty_BigIntT(1);
					copy_BigIntT(big_element(startingMatrix, i, j), currentMatrixElements[i][j]);
				}
			}
			
			identity = identity_BigIntMatrixT(big_rows(startingMatrix));
			currentMatrix = new_BigIntMatrixT(big_rows(startingMatrix), big_rows(startingMatrix));
			set_big_matrix(currentMatrix, currentMatrixElements);
			
			temp  = empty_BigIntT(1);
			temp2 = empty_BigIntT(1);
			zero  = empty_BigIntT(1);
			one   = new_BigIntT(oneArr, 1);
			
			theCycle = new_CycleInfoT();
			
			//Check matrices at given step sizes apart in all directions
			for (int incRow = 0; incRow < big_rows(startingMatrix); incRow += 1)
			{
				for (int incCol = 0; incCol < big_rows(startingMatrix); incCol += 1)
				{
					//Step in one direction in the matrix space until we look back
					do
					{
						add_BigIntT(currentMatrixElements[incRow][incCol], step, temp);
						mod_BigIntT(temp, mod, currentMatrixElements[incRow][incCol]);
						set_big_matrix(currentMatrix, currentMatrixElements);
						
						big_floyd(currentMatrix, identity, mod, &theCycle);
						
						charaPoly  = chara_poly(currentMatrix, mod);
						factorList = factor_BigPolyT(charaPoly, mod);
						
						printbm(currentMatrix);
						printf("Transient length: %d\n", tau(theCycle));
						
						//It would be real beneficial to make a "clear_CycleInfoT() function"
						theCycle = free_CycleInfoT(theCycle);
						theCycle = new_CycleInfoT();
						
						//Print all found factors
						printp(charaPoly);
						printf("\n");
						copy_BigIntT(zero, temp);
						counter = 1;
						while (compare_BigIntT(temp, constant(factorList[0])) < 0)
						{
							printf("(");
							printp(factorList[counter]);
							printf(")");
							
							counter += 1;
							add_BigIntT(temp, one, temp2);
							copy_BigIntT(temp2, temp);
						}
						printf("\n");
					}
					while (!compare_BigIntMatrixT(startingMatrix, currentMatrix));
					
					printf("------------------------------\n");
				}
			}
			
			
			for (int i = 0; i < big_rows(startingMatrix); i += 1)
			{
				for (int j = 0; j < big_rows(startingMatrix); j += 1)
					currentMatrixElements[i][j] = free_BigIntT(currentMatrixElements[i][j]);
				
				FREE(currentMatrixElements[i]);
			}
			FREE(currentMatrixElements);
			
			mod   = free_BigIntT(mod);
			step  = free_BigIntT(step);
			temp  = free_BigIntT(temp);
			temp2 = free_BigIntT(temp2);
			zero  = free_BigIntT(zero);
			one   = free_BigIntT(one);
			
			charaPoly = free_BigPolyT(charaPoly);
			
			FREE(factorList);
			
			startingMatrix = free_BigIntMatrixT(startingMatrix);
			currentMatrix  = free_BigIntMatrixT(currentMatrix);
			identity       = free_BigIntMatrixT(identity);
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
		// and transient lengths, then repeat for higher-powered moduli.
		else if (! strcmp(argv[1], "dynamics"))
		{
			int  highpower = 2;
			int  highmod = 1;   //Holds the modulus raised to the correct power

			int* possibleCycleLengths = NULL;
			int* currVectElements;
			
			//Each vector's lower-powered cycle lengths,
			// used for finding correct previous tuple
			int* currCycleLengths;
			
			int firstIndex; //Used to help create our lookup table for cycle lengths
			int tempModulus;
			
			//Holds counts for all possible cycle length pairs
			//For instance, if a vector % highmod has a cycle length
			// of 20, while % modulus it has a cycle length of 2,
			// then cycleTable[w=2][w=20] += 1;
			int** cycleTable = NULL;
			
			//Each time we go up a modulus, we have to keep track of the pevious
			// pairs of cycles we've found.
			int** prevCycleTuples = NULL;
			
			bool checkedAllVects = FALSE;
			bool matchedTuple;
			
			IntMatrixTP A; //Holds our update matrix
			IntMatrixTP currVect;
			
			CycleInfoTP theCycle = NULL;
			
			//User-provided upper bound for moduli power
			if (argc > 2)
			{
				highpower = (int)strtol(argv[2], &tempStr, 10);
				if (tempStr[0] != '\0')
				{
					fprintf(stderr, "Unable to read upper bound from command line.\n");
					return EXIT_FAILURE;
				}
			}
			
			//If user provided a modulus for us
			if (argc > 3)
			{
				modulus = (int)strtol(argv[3], &tempStr, 10);
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
			
			printf("Matrix:\n");
			printm(A);
			
			for (int modulusCounter = 1; modulusCounter <= highpower; modulusCounter += 1)
			{
				checkedAllVects = FALSE;
				
				highmod *= modulus;
				printf("Current modulus: %d\n", highmod);
				
				//Now that we have the matrix, calculate the cycle length
				// and use it to reason all possible cycle lengths
				// for alloting memory
				theCycle = floyd(A, A, highmod);
				printf("Matrix's multiplicative order mod %d: %d\n", highmod, omega(theCycle));
				
				//printf("(lowmod, highmod)\n\n");
				
				//First number in the pointer says how many cycle lengths
				// are stored in the vector
				if (possibleCycleLengths != NULL)
				{
					FREE(possibleCycleLengths);
				}
				possibleCycleLengths = malloc(2*sizeof(int));
				possibleCycleLengths[0] = 0;
				
				//Adding all divisors of the matrix's cycle length to our list
				//Reallocating the array each time is slow, but the numbers
				// we're dealing with should be low enough so it isn't a problem
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
				
				theCycle = free_CycleInfoT(theCycle);
				
				/*printf("Possible vector cycle lengths: ");
				for (int i = 1; i <= possibleCycleLengths[0]; i += 1)
					printf("%d ", possibleCycleLengths[i]);
				printf("\n\n"); */
				
				//Now, construct our cycleTable to hold cycle counts
				//Only allocating "half" the table since that's all we need
				//This gets reallocated for each power we go up
				if (modulusCounter == 1)
				{
					//We only need one row on the first pass; we have no tuples yet
					cycleTable    = malloc(sizeof(int*));
					cycleTable[0] = calloc(possibleCycleLengths[0], sizeof(int));
				}
				
				else
				{
					cycleTable = malloc(prevCycleTuples[0][0]*sizeof(int*));
					for (int i = 0; i < prevCycleTuples[0][0]; i += 1)
					{
						//This line can be improved by only allocating enough
						// ints for those equal to or greater than the cycles
						// within the tuple
						cycleTable[i] = calloc(possibleCycleLengths[0], sizeof(int));
					}
				}
				
				currVectElements     = calloc(rows(A), sizeof(int));
				currVect             = new_IntMatrixT(rows(A), 1);
				
				//Iterate until we've tested every vector
				while (!checkedAllVects)
				{
					//printf("---\n");
					set_column(currVect, currVectElements);
					
					//Print out vector in an easier-to-look-at way for stdout
					/*
					printf("<");
					for (int i = 0; i < rows(A)-1; i += 1)
						printf("%d ", currVectElements[i]);
					printf("%d>\n", currVectElements[rows(A)-1]);
					printf("---\n"); 
					*/
					
					//Helps find which previous tuple this vector corresponds to
					//Contains the cycle lengths of our vector under all moduli
					// less than the current one. The last element is the
					// current modulus' cycle length
					currCycleLengths = malloc(modulusCounter*sizeof(int));
					
					//Start with highest modulus, work our way down
					//Allows taking vectors mod tempModulus to be easier
					tempModulus = modulus;
					for (int i = 1; i < modulusCounter; i += 1, tempModulus *= modulus);
					
					for (int i = modulusCounter-1; i >= 0; i -= 1)
					{
						theCycle = floyd(A, currVect, tempModulus);
						currCycleLengths[i] = omega(theCycle);
						
						/*
						//Looking for particular cycle lengths. 
						if (omega(theCycle) == 1)
						{
							printm(currVect);
						}
						*/
						theCycle = free_CycleInfoT(theCycle);
						
						tempModulus /= modulus;
						modm(currVect, tempModulus); //Might cause an error mod 1?
					}
					
					//Now, find the correct firstIndex which matches our cycle lengths
					if (modulusCounter == 1)
						firstIndex = 0;
					
					else
					{
						for (int tuple = 1; tuple <= prevCycleTuples[0][0]; tuple += 1)
						{
							//Now, check the specific elements of each tuple to see if they
							// match with ours
							matchedTuple = TRUE;
							for (int tupleElem = 0; tupleElem < modulusCounter-1; tupleElem += 1)
							{
								if (prevCycleTuples[tuple][tupleElem] != currCycleLengths[tupleElem])
								{
									matchedTuple = FALSE;
									break;
								}
							}
							
							//The for loop will reach this at some point since the vectors must fall into
							// one of the defined cycle reductions (tuples)
							if (matchedTuple)
							{
								firstIndex = tuple-1;
								break;
							}
						}
					}
					
					//Now, we simply add our count to the cycleTable entry
					//Find our secondIndex
					for (int secondIndex = 1; 
					secondIndex <= possibleCycleLengths[0]; 
					secondIndex += 1)
					{
						if (possibleCycleLengths[secondIndex] == currCycleLengths[modulusCounter-1])
						{
							cycleTable[firstIndex][secondIndex-1] += 1;
							break;
						}
					}

					//Increment through all possible vectors (% modulus)
					checkedAllVects = TRUE;
					for (int i = rows(A)-1; i >= 0; i -= 1)
					{
						if (currVectElements[i] != highmod - 1)
						{
							checkedAllVects = FALSE;
							currVectElements[i] += 1;
							break;
						}
						
						//Carry to next vector component
						else
							currVectElements[i] = 0;
					}
				}
				
				//Now that we have a completed cycleTable up to a particular modulus,
				// we need to get our list of tuples from it.
				if (prevCycleTuples == NULL)
				{
					prevCycleTuples    = malloc(sizeof(int*));
					prevCycleTuples[0] = malloc(sizeof(int));
					prevCycleTuples[0][0] = 0;
					
					//Add cycle lengths as tuples
					for (int i = 0; i < possibleCycleLengths[0]; i += 1)
					{
						if (cycleTable[0][i] != 0)
						{
							//Increasing tuple count
							prevCycleTuples[0][0] += 1;
							prevCycleTuples = realloc(prevCycleTuples, (prevCycleTuples[0][0]+1)*sizeof(int*));
							prevCycleTuples[prevCycleTuples[0][0]] = malloc(sizeof(int));
							prevCycleTuples[prevCycleTuples[0][0]][0] = possibleCycleLengths[i+1];
						}
					}
					
					//Prevents a segmentation fault if highpower == 1
					if (modulusCounter != highpower)
					{
						FREE(cycleTable[0]);
						FREE(cycleTable);
					}
					
					if (highpower == 1)
					{
						//If we're only checking one iteration, we don't need to
						// prepare this value for a next iteration. It does, however,
						// need to say how many tuples we had last iteration, which 
						// in this case was 1 (the null tuple)
						prevCycleTuples[0][0] = 1;
					}
				}
				
				else
				{
					//Get cycleTable ready for next iteration
					if (modulusCounter != highpower)
					{
						int prevPrevTupleCount = prevCycleTuples[0][0];
						
						//New tuples get put here, then prevCycleTuples points to this after
						int** newPrevCycleTuples = malloc(sizeof(int*));
						newPrevCycleTuples[0] = calloc(1, sizeof(int));
						
						for (int tuple = 0; tuple < prevPrevTupleCount; tuple += 1)
						{
							for (int newCycle = 0; newCycle < possibleCycleLengths[0]; newCycle += 1)
							{
								if (cycleTable[tuple][newCycle] != 0)
								{
									//Add new tuple
									newPrevCycleTuples[0][0] += 1;
									newPrevCycleTuples = realloc(newPrevCycleTuples, (newPrevCycleTuples[0][0]+1)*sizeof(int*));
									newPrevCycleTuples[newPrevCycleTuples[0][0]] = malloc(modulusCounter*sizeof(int));
									
									//Re-add elements from previous tuple, add new cycle length
									for (int tupleElem = 0; tupleElem < modulusCounter-1; tupleElem += 1)
										newPrevCycleTuples[newPrevCycleTuples[0][0]][tupleElem] = prevCycleTuples[tuple+1][tupleElem];
									newPrevCycleTuples[newPrevCycleTuples[0][0]][modulusCounter-1] = possibleCycleLengths[newCycle+1];
								}
							}
						}

						for (int i = 0; i < prevPrevTupleCount; i += 1)
						{
							FREE(cycleTable[i]);
						}
						FREE(cycleTable);
						
						//Now, deallocate prevCycleLengths, have it point to new tuples
						for (int i = 1; i < prevCycleTuples[0][0]; i += 1)
						{
							FREE(prevCycleTuples[i]);
						}
						FREE(prevCycleTuples[0]);
						FREE(prevCycleTuples);
						prevCycleTuples = newPrevCycleTuples;
					}
				}
			}
			
			//Now, we print all the numbers we've collected
			printf("\nCycle tuples:\n");
			for (int i = 0; i < prevCycleTuples[0][0]; i += 1)
			{
				for (int j = 0; j < possibleCycleLengths[0]; j += 1)
				{
					//Only print out info is vectors have that pair of cycle lengths
					if (cycleTable[i][j] != 0)
					{
						//Print relevant tuple
						printf("(");
						for (int tupleElem = 0; tupleElem < highpower-1; tupleElem += 1)
							printf("%d, ", prevCycleTuples[i+1][tupleElem]);
						printf("%d) = %d\n", possibleCycleLengths[j+1], cycleTable[i][j]);
					}
				}
			}
			
			A        = free_IntMatrixT(A);
			currVect = free_IntMatrixT(currVect);
			
			//Don't need to free theCycle since it gets freed above
			theCycle = NULL;
			
			for (int i = 0; i < prevCycleTuples[0][0]; i += 1)
			{
				FREE(cycleTable[i]);
			}
			FREE(cycleTable);
			FREE(possibleCycleLengths);
			
			FREE(currVectElements);
			
			for (int i = 1; i < prevCycleTuples[0][0]; i += 1)
			{
				FREE(prevCycleTuples[i]);
			}
			FREE(prevCycleTuples[0]);
			FREE(prevCycleTuples);
			
			FREE(currCycleLengths);
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
		
		printf(" - " ANSI_COLOR_YELLOW "core " ANSI_COLOR_CYAN "[modulus]" ANSI_COLOR_RESET \
		": Calculates the number of vectors in a matrix's core.\n\n");
		
		printf(" - " ANSI_COLOR_YELLOW "floyd " ANSI_COLOR_CYAN "[modulus]" ANSI_COLOR_RESET \
		": Use Floyd's Cycle Detection Algorithm to find out specific details about the given LCA.\n\n");
		
		printf(" - " ANSI_COLOR_YELLOW "bigfloyd " ANSI_COLOR_CYAN "[modulus]" ANSI_COLOR_RESET \
		": Same as floyd, but uses BigIntMatrixT and BigIntT structs instead of IntMatrixT and ints.\n\n");
		
		printf(" - " ANSI_COLOR_YELLOW "rots " ANSI_COLOR_CYAN "[modulus]" ANSI_COLOR_RESET \
		": Finds and outputs some basic rotation matrices for the given modulus.\n\n");
		
		printf(" - " ANSI_COLOR_YELLOW "cycmatsearch " ANSI_COLOR_CYAN "resume size maxmod cycles..." ANSI_COLOR_RESET \
		": Searches for a matrix whose column vectors cycle with cycle lengths less than the matrix itself.\n\n");
		
		printf(" - " ANSI_COLOR_YELLOW "cycconvmat " ANSI_COLOR_CYAN "from to [mod]" ANSI_COLOR_RESET \
		": Outputs a \"cycle converting matrix\" for the given update matrix.\n\n");
		
		printf(" - " ANSI_COLOR_YELLOW "charawalk" ANSI_COLOR_CYAN " step [modulus]" ANSI_COLOR_RESET \
		": Steps around a matrix space and computes the characteritic polynomial for each matrix it lands on.\n\n");
		
		printf(" - " ANSI_COLOR_YELLOW "fibcycle" ANSI_COLOR_CYAN " [modulus]" ANSI_COLOR_RESET \
		": Generate the Fibonacci cycle that contains the initial vector.\n\n");
		
		printf(" - " ANSI_COLOR_YELLOW "fibcyclelens" ANSI_COLOR_CYAN " [modulus]" ANSI_COLOR_RESET \
		": Calculate all possible Fibonacci cycle lengths.\n\n");
		
		printf(" - " ANSI_COLOR_YELLOW "fibmultsearch " ANSI_COLOR_CYAN "[bound]" ANSI_COLOR_RESET \
		": Searches the Fibonacci numbers, checking whether a " \
		"multiple of each number up to the bound appears before a multiple of a power of the number.\n\n");
		
		printf(" - " ANSI_COLOR_YELLOW "dynamics " ANSI_COLOR_CYAN "[maxPower] [modulus]" ANSI_COLOR_RESET \
		": Iterates every vector in a space, recording their transient lengths and cycle lengths. " \
		"It then computes the same numbers for higher-powered moduli.\n\n");
		
		printf("For a more complete description of LINCELLAUT's usage, " \
		"refer to the included documentation.\n");
	}
	
	//Freeing memory
	FREE_VARIABLES;
	
	return EXIT_SUCCESS;
}
