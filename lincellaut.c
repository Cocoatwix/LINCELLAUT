
/* A simple program for calculating various things about linear
 *  cellular automaton.
 *
 * Mar 27, 2022
 *
 */
 
/* The following resources were used as a reference:
docs.microsoft.com/en-us/cpp/c-language/cpp-integer-limits?view=msvc-170
stackoverflow.com/questions/3219393

stackoverflow.com/questions/1190870
*/
 
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h> //So I can get the maximum integer

#include <time.h> //For seeding random numbers 

#include "headers/helper.h"
#include "headers/bigint.h" //Arbitrary precision

#include "headers/linalg.h" 
#include "headers/cycles.h"  //Allows us to use Floyd's Algorithm
#include "headers/modular.h" //Modular square roots and inverses
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
											 FREE(resumemodstring); \
											 UPDATEMATRIX  = free_BigIntMatrixT(UPDATEMATRIX); \
											 INITIALMATRIX = free_BigIntMatrixT(INITIALMATRIX)
											 
#define SET_BIG_NUM(str, bigint, msg) \
if (strtoBIT((str), &(bigint)) != 1) \
{ \
	fprintf(stderr, "%s\n", (msg)); \
	FREE_VARIABLES; \
	return EXIT_FAILURE; \
}


int main(int argc, char* argv[])
{
	srand(time(NULL)); //Seeding random number generator
	
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
	
	VectorTypeT vectorType = row;
	
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
	
	BigIntMatrixTP UPDATEMATRIX  = NULL;
	BigIntMatrixTP INITIALMATRIX = NULL;
	
	//Temporarily holds config data
	char* systemData  = malloc(MAXSTRLEN*sizeof(char));
	
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
				FREE_VARIABLES;
				return EXIT_FAILURE;
			}
			
			if (strlen(bigintmodstring) < 10)
			{
				modulus = (int)strtol(bigintmodstring, &tempStr, 10);
				if (tempStr[0] != '\0')
				{
					fprintf(stderr, "Invalid modulus provided in config file.\n");
					FREE_VARIABLES;
					return EXIT_FAILURE;
				}
			}
		}
		
		//If the user wants to change how vectors are displayed to the console
		else if (! strcmp(systemData, "vectorType"))
		{
			if (fscanf(system, "%s", systemData) != 1)
			{
				fprintf(stderr, "Unable to read vector type from config file.\n");
				FREE_VARIABLES;
				return EXIT_FAILURE;
			}
			
			else
				if (!strcmp(systemData, "col"))
					vectorType = col;
		}
		
		else if (! strcmp(systemData, "iterations"))
		{
			if (fscanf(system, "%d", &iterations) != 1)
			{
				fprintf(stderr, "Unable to read number of iterations from config file.\n");
				FREE_VARIABLES;
				return EXIT_FAILURE;
			}
		}
		
		else if (! strcmp(systemData, "update"))
		{
			if (fscanf(system, "%s", updatefilepath) != 1)
			{
				fprintf(stderr, "Unable to read update matrix path from config file.\n");
				FREE_VARIABLES;
				return EXIT_FAILURE;
			}
			
			UPDATEMATRIX = read_BigIntMatrixT(updatefilepath);
			if (UPDATEMATRIX == NULL)
			{
				fprintf(stderr, "Unable to set update matrix from config file.\n");
				FREE_VARIABLES;
				return EXIT_FAILURE;
			}
		}
		
		else if (! strcmp(systemData, "initial"))
		{
			if (fscanf(system, "%s", initialfilepath) != 1)
			{
				fprintf(stderr, "Unable to read initial vector path from config file.\n");
				FREE_VARIABLES;
				return EXIT_FAILURE;
			}
			
			INITIALMATRIX = read_BigIntMatrixT(initialfilepath);
			if (INITIALMATRIX == NULL)
			{
				fprintf(stderr, "Unable to set initial matrix from config file.\n");
				FREE_VARIABLES;
				return EXIT_FAILURE;
			}
		}
		
		else if (! strcmp(systemData, "resumeMat"))
		{
			if (fscanf(system, "%s", resumefilepath) != 1)
			{
				fprintf(stderr, "Unable to read resume matrix path from config file.\n");
				FREE_VARIABLES;
				return EXIT_FAILURE;
			}
		}
		
		else if (! strcmp(systemData, "resumeMod"))
		{
			if (fscanf(system, "%101s", resumemodstring) != 1)
			{
				fprintf(stderr, "Unable to read resume modulus from config file.\n");
				FREE_VARIABLES;
				return EXIT_FAILURE;
			}
		}
		
		else if (! strcmp(systemData, "itername"))
		{
			if (fscanf(system, "%s", iterfilepath) != 1)
			{
				fprintf(stderr, "Unable to read iteration file path from config file.\n");
				FREE_VARIABLES;
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
		// This is written with IntMatrices because I don't yet have a "big_inverse()"
		if (! strcmp(argv[1], "iterate"))
		{
			BigIntMatrixTP F;
			BigIntMatrixTP F_2;
			BigIntMatrixTP F_result = NULL;
			BigIntMatrixTP power = NULL;
			
			int tempArr[1];
			BigIntTP bigIters = NULL;
			
			BigIntTP bigMod;
			SET_BIG_NUM(bigintmodstring, bigMod, "Unable to read modulus from config file.");
		
			//If the user provided a custom number of iterations
			if (argc > 2)
			{
				iterations = (int)strtol(argv[2], &tempStr, 10);
				
				//If we didn't read any digits for iterations
				if (tempStr[0] != '\0')
				{
					fprintf(stderr, "Invalid number of iterations provided at command line.\n");
					bigMod = free_BigIntT(bigMod);
					
					FREE_VARIABLES;
					return EXIT_FAILURE;
				}
			}
			
			printf("Iterations: %d\n", iterations);
			if (iterations < 0)
			{
				F_2 = read_BigIntMatrixT(updatefilepath);
				F   = big_inverse(F_2, bigMod);
				
				if (F == NULL)
				{
					printf("An inverse for the update matrix modulo ");
					printi(bigMod);
					printf(" doesn't exist.\n");
				}
				
				else
				{
					iterations *= -1;
					F_2 = free_BigIntMatrixT(F_2);
					F_2 = read_BigIntMatrixT(initialfilepath);
				}
			}
			
			else
			{
				F   = read_BigIntMatrixT(updatefilepath);
				F_2 = read_BigIntMatrixT(initialfilepath);
			}
			
			if ((F == NULL) || (F_2 == NULL))
			{
				fprintf(stderr, "Unable to read matrices from config file.\n");
				F = free_BigIntMatrixT(F);
				F_2 = free_BigIntMatrixT(F_2);
				bigMod = free_BigIntT(bigMod);
				
				FREE_VARIABLES;
				return EXIT_FAILURE;
			}
			
			if (big_cols(F) != big_rows(F_2))
			{
				fprintf(stderr, "Given matrices cannot be multiplied.\n");
				F = free_BigIntMatrixT(F);
				F_2 = free_BigIntMatrixT(F_2);
				bigMod = free_BigIntT(bigMod);
				
				FREE_VARIABLES;
				return EXIT_FAILURE;
			}
			
			//Iterate F a few times
			//The minus 1 is for easier conversion between ORBITVIS results
			if (iterations == 1)
			{
				F_result = new_BigIntMatrixT(big_rows(F), big_cols(F_2));
				big_mat_mul(F, F_2, F_result);
				modbm(F_result, bigMod);
			}
			
			else if (iterations > 0)
			{
				tempArr[0] = iterations;
				bigIters = new_BigIntT(tempArr, 1);
				F_result = new_BigIntMatrixT(big_rows(F_2), big_cols(F_2));
				
				power = new_BigIntMatrixT(big_rows(F), big_cols(F));
				powbm(F, power, bigIters, bigMod);
				big_mat_mul(power, F_2, F_result);
				modbm(F_result, bigMod);
			}
			
			//Prevents the matrix from being printed when
			// an inverse doesn't exist and the iterations
			// provided was negative
			if (iterations > 0)
			{
				if (vectorType == row)
					printbm_row(F_result);
				else
					printbm(F_result);
			}
			
			else if (iterations == 0)
			{
				if (vectorType == row)
					printbm_row(F_2);
				else
					printbm(F_2);
			}
			
			if ((vectorType == row) && (big_rows(F_2) != big_cols(F_2)))
				printf("\n");

			F        = free_BigIntMatrixT(F);
			F_2      = free_BigIntMatrixT(F_2);
			F_result = free_BigIntMatrixT(F_result);
			power    = free_BigIntMatrixT(power);
			
			bigIters = free_BigIntT(bigIters);
		}
		
		
		//Find the inverse of the update matrix
		else if (!strcmp(argv[1], "inverse"))
		{
			BigIntMatrixTP F = read_BigIntMatrixT(updatefilepath);
			BigIntMatrixTP Finv;
			
			BigIntTP bigMod;
			
			if (F == NULL)
			{
				fprintf(stderr, "Unable to read update matrix from config file.\n");
				FREE_VARIABLES;
				return EXIT_FAILURE;
			}
			
			printf("Update matrix:\n");
			printbm(F);
			
			//If the user specified a modulus at the command line
			if (argc > 2)
			{
				SET_BIG_NUM(argv[2], bigMod, "Unable to read modulus from command line.");
			}
			else
			{
				SET_BIG_NUM(bigintmodstring, bigMod, "Unable to read modulus from command line.");
			}
			
			Finv = big_inverse(F, bigMod);
			
			if (Finv == NULL)
			{
				printf("An inverse for the update matrix modulo ");
				printi(bigMod);
				printf(" doesn't exist.\n");
			}
			
			else
			{
				printf("Inverse modulo ");
				printi(bigMod);
				printf(":\n");
				printbm(Finv);
			}
			
			F    = free_BigIntMatrixT(F);
			Finv = free_BigIntMatrixT(Finv);
			
			bigMod = free_BigIntT(bigMod);
		}
		
		
		//Find the characteristic equation of the update matrix
		else if (! strcmp(argv[1], "chara"))
		{
			BigIntTP bigMod = NULL;
			
			BigIntMatrixTP bigMatrix;
			
			BigPolyTP bigEqn;
			BigFactorsTP bigEqnFactors = NULL;
			BigPolyTP* minEqn = NULL;
			
			//If the user provided a custom modulus
			if (argc > 2)
			{
				SET_BIG_NUM(argv[2], bigMod, "Unable to read modulus from command line.");
			}
			
			//Extract the modulus in the .config file otherwise
			else
			{
				SET_BIG_NUM(bigintmodstring, bigMod, "Unable to read modulus from config file.");
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
			printf("Modulus: ");
			printi(bigMod);
			printf("\n");
			bigEqn = chara_poly(bigMatrix, bigMod);
			printf("\nCharacteristic polynomial: ");
			printp(bigEqn);
			
			printf("\nFactored characteristic polynomial: ");
			bigEqnFactors = factor_BigPolyT(bigEqn, bigMod);
			printpf(bigEqnFactors);
			printf("\nMinimal polynomial: ");
			minEqn = min_poly(bigMatrix, NULL, bigMod, bigEqn);
			old_printpf(minEqn);
			printf("\n\n");
			
			bigMod  = free_BigIntT(bigMod);
			
			bigMatrix     = free_BigIntMatrixT(bigMatrix);
			bigEqn        = free_BigPolyT(bigEqn);
			bigEqnFactors = free_BigFactorsT(bigEqnFactors);
			minEqn        = free_BigPolyT_factors(minEqn);
		}
		
		
		//If the user wants to find all matrices with a particular characteristic polynomial
		else if (! strcmp(argv[1], "allcharas"))
		{
			BigIntTP bigMod;
			BigIntTP negOne;
			
			BigIntTP** currMatElements;
			BigIntMatrixTP currMat;
			
			BigIntTP* coeffs;
			BigPolyTP targetEqn;
			BigPolyTP tempPoly;
			BigPolyTP negTargetEqn;
			BigPolyTP tempEqn;
			BigPolyTP negOnePoly;
			
			bool newCycle;
			int tempCurrCycle;
			int numOfCycleLengths = 0;
			int maxArraySize = 0;
			int* cycleLengths = NULL;
			int* cycleLengthCounts = NULL;
			
			BigIntTP** currVectElements;
			BigIntMatrixTP currVect;
			CycleInfoTP theCycle = NULL;
			
			int oneArr[1] = {1};
			BigIntTP one;
			
			if (argc == 2)
			{
				printf("No characteristic polynomial was given on command line.\n");
				FREE_VARIABLES;
				return EXIT_SUCCESS;
			}
			
			SET_BIG_NUM(bigintmodstring, bigMod, "Unable to read modulus from config file.");
			
			one = new_BigIntT(oneArr, 1);
			negOne = empty_BigIntT(1);
			subtract_BigIntT(bigMod, one, negOne);
			negOnePoly = new_BigPolyT(&negOne, 1);
			
			//Get coefficients for targetEqn
			coeffs = malloc((argc-2)*sizeof(BigIntTP));
			for (int i = 2; i < argc; i += 1)
				strtoBIT(argv[i], &coeffs[i-2]);
			
			targetEqn = new_BigPolyT(coeffs, argc-2);
			reduce_BigPolyT(targetEqn);
			
			negTargetEqn = new_BigPolyT(&one, 1);
			tempPoly = new_BigPolyT(&one, 1);
			
			multiply_BigPolyT(targetEqn, negOnePoly, tempPoly);
			mod_BigPolyT(tempPoly, bigMod, negTargetEqn);
			
			currMatElements = new_BigIntT_array(degree(targetEqn), degree(targetEqn));
			currMat = new_BigIntMatrixT(degree(targetEqn), degree(targetEqn));
			
			currVectElements = new_BigIntT_array(degree(targetEqn), 1);
			currVect = new_BigIntMatrixT(degree(targetEqn), 1);
			
			printf("Target characteristic eqn: ");
			printp(targetEqn);
			printf(" or ");
			printp(negTargetEqn);
			printf("\nModulus: ");
			printi(bigMod);
			printf("\n~~~\n");
			
			//Loop through all matrices, find the ones with the correct characteristic polynomial
			do
			{
				set_big_matrix(currMat, currMatElements);
				tempEqn = chara_poly(currMat, bigMod);
				
				//If we found a matrix with the correct characteristic polynomial
				if ((compare_BigPolyT(tempEqn, targetEqn) == 0) ||
				    (compare_BigPolyT(tempEqn, negTargetEqn) == 0))
				{
					printbm(currMat);
					
					//Now, let's find cycle length counts
					do
					{
						set_big_matrix(currVect, currVectElements);
						big_floyd(currMat, currVect, bigMod, &theCycle);
						newCycle = TRUE;
						
						//Try and see if we've found this cycle length before
						for (int i = 0; i < numOfCycleLengths; i += 1)
						{
							if (cycleLengths[i] == omega(theCycle))
							{
								cycleLengthCounts[i] += 1;
								newCycle = FALSE;
								break;
							}
						}
						
						//If we need to add a new cycle length to our list
						if (newCycle)
						{
							numOfCycleLengths += 1;
							
							if (maxArraySize < numOfCycleLengths)
							{
								cycleLengths = realloc(cycleLengths, numOfCycleLengths*sizeof(int));
								cycleLengthCounts = realloc(cycleLengthCounts, numOfCycleLengths*sizeof(int));
								
								maxArraySize = numOfCycleLengths;
							}
							
							cycleLengths[numOfCycleLengths-1] = omega(theCycle);
							cycleLengthCounts[numOfCycleLengths-1] = 1;
						}
					}
					while (!increment_BigIntT_array(currVectElements, degree(targetEqn), 1, one, bigMod));
					
					//Print out our found cycle lengths IN ORDER
					newCycle = TRUE;
					tempCurrCycle = 1;
					printf("1 (%d)", cycleLengthCounts[0]);
					while (newCycle)
					{
						newCycle = FALSE;
						for (int i = 0; i < numOfCycleLengths; i += 1)
							if (cycleLengths[i] > tempCurrCycle)
							{
								newCycle = TRUE;
								tempCurrCycle = cycleLengths[i];
								
								printf(", %d (%d)", tempCurrCycle, cycleLengthCounts[i]);
							}
					}
					printf("\n\n");

					numOfCycleLengths = 0;
				}
				
				tempEqn = free_BigPolyT(tempEqn);
			}
			while (!increment_BigIntT_array(currMatElements, degree(targetEqn), degree(targetEqn), one, bigMod));
			
			
			bigMod = free_BigIntT(bigMod);
			one    = free_BigIntT(one);
			negOne = free_BigIntT(negOne);
			
			currMatElements  = free_BigIntT_array(currMatElements, big_rows(currMat), big_rows(currMat));
			currVectElements = free_BigIntT_array(currVectElements, big_rows(currMat), 1);
			currMat  = free_BigIntMatrixT(currMat);
			currVect = free_BigIntMatrixT(currVect);
			theCycle = free_CycleInfoT(theCycle);
			FREE(cycleLengths);
			FREE(cycleLengthCounts);
			
			for (int i = 0; i < argc-2; i += 1)
				coeffs[i] = free_BigIntT(coeffs[i]);
			FREE(coeffs);
			
			targetEqn    = free_BigPolyT(targetEqn);
			negTargetEqn = free_BigPolyT(negTargetEqn);
			negOnePoly   = free_BigPolyT(negOnePoly);
			tempPoly     = free_BigPolyT(tempPoly);
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
				SET_BIG_NUM(argv[2], bigMod, "Invalid modulus passed at command line.");
			}
			
			else
			{
				SET_BIG_NUM(bigintmodstring, bigMod, "Invalid modulus passed in config file.");
			}
			
			A = read_BigIntMatrixT(updatefilepath);
			if (A == NULL)
			{
				fprintf(stderr, "Unable to read matrix file at %s.\n", updatefilepath);
				return EXIT_FAILURE;
			}
			
			entriesToUse = new_BigIntT_array(big_cols(A), 1);
			
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
				checkedAllVects = increment_BigIntT_array(entriesToUse, big_cols(A), 1, one, bigMod);
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
			
			
			entriesToUse = free_BigIntT_array(entriesToUse, big_cols(A), 1);
			
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
		
		
		//If the user wants to calculate the specific vectors in a
		// matrix's orbits
		else if (! strcmp(argv[1], "orbits"))
		{
			BigIntTP bigMod;
			BigIntTP overflowCheck;
			int smallMod;
			int totalVects = 1;
			int foundVectors = 0;
			int tempVectCount = 0;
			BigIntTP tempIndex1;
			BigIntTP tempIndex2;
			BigIntTP zero;
			
			BigIntMatrixTP  A;
			BigIntMatrixTP  currVect;
			BigIntTP** currVectElements;
			
			CycleInfoTP theCycle = NULL;
			
			int numOfOrbits = 0;
			BigIntMatrixTP* orbitReps = NULL;
			bool foundNewOrbit = FALSE;
			
			int oneArr[1] = {1};
			BigIntTP one;
			BigIntTP temp;
			BigIntTP temp2;
			BigIntTP temp3;
			
			BigIntMatrixTP tempVect;
			BigIntMatrixTP tempVect2;
			
			FILE* outputFile = NULL;
			char* outputFileName = NULL;
			bool fileoutput = FALSE;
			
			FILE* graphFile = NULL;
			char* graphFileName = NULL;
			
			if (argc > 2)
			{
				SET_BIG_NUM(argv[2], bigMod, "Invalid modulus passed on command line.");
				smallMod = (int)strtol(argv[2], &tempStr, 10);
			}
			else
			{
				SET_BIG_NUM(bigintmodstring, bigMod, "Invalid modulus passed in config file.");
				smallMod = (int)strtol(bigintmodstring, &tempStr, 10);
			}
			
			if (argc > 3)
				if (! strcmp(argv[3], "TRUE"))
					fileoutput = TRUE;
				
			A = read_BigIntMatrixT(updatefilepath);
			if (A == NULL)
			{
				fprintf(stderr, "Unable to read update matrix from config file.\n");
				FREE_VARIABLES;
				return EXIT_FAILURE;
			}
			
			if (big_rows(A) != big_cols(A))
			{
				fprintf(stderr, "Given update matrix isn't square.\n");
				FREE_VARIABLES;
				return EXIT_SUCCESS;
			}
				
			//Creating filename and file
			if (fileoutput)
			{
				outputFileName = malloc(MAXSTRLEN*sizeof(char));
				outputFileName[0] = '\0';
				
				graphFileName = malloc(MAXSTRLEN*sizeof(char));
				graphFileName[0] = '\0';
				
				strcat(outputFileName, "orbits ");
				strcat(graphFileName, "graph ");
				append_BigIntT(outputFileName, bigMod);
				append_BigIntT(graphFileName, bigMod);
				strcat(outputFileName, " F");
				strcat(graphFileName, " F");
				for (int row = 0; row < big_rows(A); row += 1)
					for (int col = 0; col < big_rows(A); col += 1)
					{
						append_BigIntT(outputFileName, big_element(A, row, col));
						append_BigIntT(graphFileName, big_element(A, row, col));
					}
				strcat(outputFileName, ".txt");
				strcat(graphFileName, ".graph");
				
				outputFile = fopen(outputFileName, "w");
				if (outputFile == NULL)
					fprintf(stderr, "Unable to create output file. Continuing without saving...\n");
				
				graphFile = fopen(graphFileName, "w");
				if (graphFile == NULL)
					fprintf(stderr, "Unable to create graph file.\n");
			}
			
			currVectElements = new_BigIntT_array(big_rows(A), 1);
			currVect = new_BigIntMatrixT(big_rows(A), 1);
			
			one = new_BigIntT(oneArr, 1);
			temp = empty_BigIntT(1);
			temp2 = empty_BigIntT(1);
			temp3 = empty_BigIntT(1);
			tempVect  = new_BigIntMatrixT(big_rows(A), 1);
			tempVect2 = new_BigIntMatrixT(big_rows(A), 1);
			
			zero = empty_BigIntT(1);
			tempIndex1 = empty_BigIntT(1);
			tempIndex2 = empty_BigIntT(1);
			
			overflowCheck = new_BigIntT(oneArr, 1);
			for (int a = 0; a < big_rows(A); a += 1)
			{
				multiply_BigIntT(bigMod, overflowCheck, temp);
				copy_BigIntT(temp, overflowCheck);
				
				totalVects *= smallMod;
			}
			
			//Iterate through all vectors, find all orbits
			do
			{
				foundNewOrbit = FALSE;
				
				set_big_matrix(currVect, currVectElements);
				big_floyd(A, currVect, bigMod, &theCycle);
				
				if (numOfOrbits == 0)
				{
					foundNewOrbit = TRUE;
					
					numOfOrbits += 1;
					orbitReps = malloc(sizeof(BigIntMatrixTP));
					orbitReps[0] = new_BigIntMatrixT(big_rows(A), 1);
					
					copy_BigIntMatrixT(rep(theCycle), orbitReps[0]);
				}
				
				//Sift through our currently-found orbits to see if our currVect is in them
				else
				{
					foundNewOrbit = TRUE;
					for (int ourRep = 0; ourRep < numOfOrbits; ourRep += 1)
					{
						copy_BigIntMatrixT(orbitReps[ourRep], tempVect);
						do
						{
							//If we found our rep of currVect's cycle in a cycle
							if (compare_BigIntMatrixT(tempVect, rep(theCycle)))
							{
								foundNewOrbit = FALSE;
								break;
							}
							
							big_mat_mul(A, tempVect, tempVect2);
							modbm(tempVect2, bigMod);
							copy_BigIntMatrixT(tempVect2, tempVect);
						}
						while (!compare_BigIntMatrixT(tempVect, orbitReps[ourRep]));
						
						if (!foundNewOrbit)
							break;
					}
					
					if (foundNewOrbit)
					{
						numOfOrbits += 1;
						orbitReps = realloc(orbitReps, numOfOrbits*sizeof(BigIntMatrixTP));
						orbitReps[numOfOrbits-1] = new_BigIntMatrixT(big_rows(A), 1);
						copy_BigIntMatrixT(rep(theCycle), orbitReps[numOfOrbits-1]);
					}
				}
				
				if (foundNewOrbit)
				{
					printf("Orbit #%d:\n", numOfOrbits-1);
					copy_BigIntMatrixT(orbitReps[numOfOrbits-1], tempVect);
					
					if (outputFile != NULL)
						fprintf(outputFile, "Orbit #%d (ω = %d)\n", numOfOrbits-1, omega(theCycle));
					
					for (int v = 0; v < omega(theCycle); v += 1)
					{
						foundVectors += 1;
						
						printbm_row(tempVect);
						printf("\n");
						
						if (outputFile != NULL)
						{
							fprintf(outputFile, "%d: ", v);
							fprintbm_row(outputFile, tempVect);
							fprintf(outputFile, "\n");
						}
						
						big_mat_mul(A, tempVect, tempVect2);
						modbm(tempVect2, bigMod);
						copy_BigIntMatrixT(tempVect2, tempVect);
					}
					
					printf("Vectors in cycles: %d\n", foundVectors);
					printf("\n");
					if (outputFile != NULL)
						fprintf(outputFile, "\n");
				}
				
				//If we've accounted for all the vectors in the module
				//The first condition prevents this from triggering if 
				if ((size(overflowCheck) == 1) && (foundVectors == totalVects))
					break;
			}
			while (!increment_BigIntT_array(currVectElements, big_rows(A), 1, one, bigMod));
			
			if (outputFile != NULL)
				if (fclose(outputFile) == EOF)
					fprintf(stderr, "Unable to save output file.\n");
				
			if (graphFile != NULL)
			{
				fprintf(graphFile, "~b:False\n");
				fprintf(graphFile, "~n:%d\n~l:above\n", foundVectors);
				for (int i = 0; i < big_rows(A); i += 1)
					copy_BigIntT(zero, currVectElements[i][0]);
				
				do
				{
					set_big_matrix(currVect, currVectElements);
					fprintf(graphFile, "%d:", tempVectCount);
					fprintbm_row(graphFile, currVect);
					fprintf(graphFile, "\n");
					tempVectCount += 1;
				}
				while (!increment_BigIntT_array(currVectElements, big_rows(A), 1, one, bigMod));
				
				fprintf(graphFile, "~c:direction\n");
				
				//Go through every cycle, add the cycles to the file
				for (int r = 0; r < numOfOrbits; r += 1)
				{
					copy_BigIntMatrixT(orbitReps[r], tempVect);
					big_mat_mul(A, tempVect, tempVect2);
					modbm(tempVect2, bigMod);
					
					//Loop around until we've recorded every mapping in the orbit
					do
					{
						//Get numbers for each vector
						copy_BigIntT(zero, tempIndex1);
						copy_BigIntT(zero, tempIndex2);
						copy_BigIntT(one, temp);
						
						for (int i = 0; i < big_rows(A); i += 1)
						{
							//Get correct powers of our modulus
							for (int j = 0; j < i; j += 1)
							{
								multiply_BigIntT(temp, bigMod, temp2);
								copy_BigIntT(temp2, temp);
							}
							
							//temp holds correct power
							multiply_BigIntT(temp, big_element(tempVect, i, 0), temp2);
							add_BigIntT(temp2, tempIndex1, temp3);
							copy_BigIntT(temp3, tempIndex1);
							
							multiply_BigIntT(temp, big_element(tempVect2, i, 0), temp2);
							add_BigIntT(temp2, tempIndex2, temp3);
							copy_BigIntT(temp3, tempIndex2);
						}
						
						//Now our tempIndices hold the correct numbers representing each vector
						fprinti(graphFile, tempIndex1);
						fprintf(graphFile, ",");
						fprinti(graphFile, tempIndex2);
						fprintf(graphFile, "\n");
						
						//Update our vectors to the next mapping in the orbit
						copy_BigIntMatrixT(tempVect2, tempVect);
						big_mat_mul(A, tempVect, tempVect2);
						modbm(tempVect2, bigMod);
					}
					while (!compare_BigIntMatrixT(tempVect, orbitReps[r]));
				}
				
				if (fclose(graphFile) == EOF)
					fprintf(stderr, "Unable to save graph file.\n");
			}
			
			FREE(outputFileName);
			FREE(graphFileName);
			bigMod = free_BigIntT(bigMod);
			
			for (int i = 0; i < numOfOrbits; i += 1)
				orbitReps[i] = free_BigIntMatrixT(orbitReps[i]);
			FREE(orbitReps);
			
			currVectElements = free_BigIntT_array(currVectElements, big_rows(A), 1);
			currVect = free_BigIntMatrixT(currVect);
			one = free_BigIntT(one);
			
			tempVect  = free_BigIntMatrixT(tempVect);
			tempVect2 = free_BigIntMatrixT(tempVect);
			
			tempIndex1 = free_BigIntT(tempIndex1);
			tempIndex2 = free_BigIntT(tempIndex2);
			temp = free_BigIntT(temp);
			temp2 = free_BigIntT(temp);
			overflowCheck = free_BigIntT(overflowCheck);
			
			theCycle = free_CycleInfoT(theCycle);
			
			A = free_BigIntMatrixT(A);
		}
		
		
		//If the user wants to calculate orbits over a splitting field for
		// the given update matrix's minimal polynomial
		else if (! strcmp(argv[1], "splitorbits"))
		{
			//Let's just hope we don't have more extensions than this
			char* extensionNames[5] = {"α", "β", "δ", "ε", "θ"};
			
			BigIntTP bigMod;
			BigIntMatrixTP bigF;
			
			BigPolyTP* minPolyFactors; 
			BigPolyTP* extDefns = NULL; //For extracting coeffs from min poly factors
			BigIntTP*  tempFactor;
			bool isSameDefn;
			
			//How many extensions we need to add to our field
			int  numOfExtensionsNeeded = 0;
			int* zeroPointer = NULL; //For placing coefficients into MultiVarExtTs
			
			//Will hold all our extensions
			MultiVarExtTP repExtension = NULL;
			
			void*** genericFelements = NULL;
			GenericMatrixTP genericF = NULL;
			
			void*** currVectElements = NULL;
			GenericMatrixTP currVect = NULL;
			
			GenericMatrixTP tempVect = NULL; //Facilitates calculations
			
			CycleInfoTP theCycle = NULL;
			bool newCycleLength  = FALSE;
			bool newCycle        = FALSE;
			int foundCycleIndex; //Used for holding a position in the arrays below
			
			int     numOfCycleLengthsFound = 0;
			int*    numOfCyclesFound;  //Holds how many cycles of each length we've found
			int*    foundCycleLengths; //Holds all the cycle lengths we've come across so far
			void*** foundReps;         //Holds representatives for all the cycles we've found thus far
			
			bool moreVectors = TRUE;
			
			int oneArr[1] = {1};
			BigIntTP one;
			BigIntTP temp;
			BigIntTP temp2;
			
			int cycleCounter = 0; //How many cycles have we found?
			BigIntTP vectorCounter = empty_BigIntT(1);
			
			char* fileOutputName = NULL;
			FILE* fileOutput = NULL;
			
			if (argc > 2)
			{
				SET_BIG_NUM(argv[2], bigMod, "Unable to read modulus from command line.");
			}
			else
			{
				SET_BIG_NUM(bigintmodstring, bigMod, "Unable to read modulus from .config file.");
			}
			
			bigF = read_BigIntMatrixT(updatefilepath);
			if (bigF == NULL)
			{
				fprintf(stderr, "Unable to read update matrix from %s.\n", updatefilepath);
				bigMod = free_BigIntT(bigMod);
				FREE_VARIABLES;
				return EXIT_FAILURE;
			}
			
			if (big_rows(bigF) != big_cols(bigF))
			{
				fprintf(stderr, "Given matrix isn't square.\n");
				bigMod = free_BigIntT(bigMod);
				bigF = free_BigIntMatrixT(bigF);
				FREE_VARIABLES;
				return EXIT_SUCCESS;
			}
			
			if (argc > 3)
				if (!strcmp(argv[3], "TRUE"))
				{
					fileOutputName = malloc(MAXSTRLEN*sizeof(char));
					fileOutputName[0] = '\0';
					strcat(fileOutputName, "splitorbits F");
					
					//Get elements of matrix
					for (int row = 0; row < big_rows(bigF); row += 1)
						for (int col = 0; col < big_cols(bigF); col += 1)
							append_BigIntT(fileOutputName, big_element(bigF, row, col));
						
					strcat(fileOutputName, " ");
					append_BigIntT(fileOutputName, bigMod);
					strcat(fileOutputName, ".txt");
					
					fileOutput = fopen(fileOutputName, "a");
					if (fileOutput == NULL)
						fprintf(stderr, "Unable to create output file. Continuing without saving...\n");
				}
				
			one   = new_BigIntT(oneArr, 1);
			temp  = empty_BigIntT(1);
			temp2 = empty_BigIntT(1);
			
			printf("Update matrix:\n");
			printbm(bigF);
			printf("Modulus: ");
			printi(bigMod);
			minPolyFactors = min_poly(bigF, NULL, bigMod, NULL);
			printf("\nMin poly: ");
			old_printpf(minPolyFactors);
			printf("\n");
			
			//Now, we check to see how many extensions we need for
			// our splitting field
			copy_BigIntT(one, temp);
			for (int i = 1; compare_BigIntT(constant(minPolyFactors[0]), temp) >= 0; i += 1)
			{
				isSameDefn = FALSE;
				if (degree(minPolyFactors[i]) >= 2)
				{
					//We need to ensure that this factor is unique.
					for (int extDefn = 0; extDefn < numOfExtensionsNeeded; extDefn += 1)
					{
						isSameDefn = TRUE;
						if (degree(extDefns[extDefn]) == degree(minPolyFactors[i]))
						{
							if (compare_BigPolyT(extDefns[extDefn], minPolyFactors[i]) != 0)
								isSameDefn = FALSE;
						}
						
						else
							isSameDefn = FALSE;
						
						if (isSameDefn)
							break;
					}
					
					//Only add new extension if we don't already have it
					if (!isSameDefn)
					{
						numOfExtensionsNeeded += 1;
						
						extDefns = realloc(extDefns, numOfExtensionsNeeded*sizeof(BigPolyTP));
						extDefns[numOfExtensionsNeeded-1] = empty_BigPolyT();
						copy_BigPolyT(minPolyFactors[i], extDefns[numOfExtensionsNeeded-1]);
					}
				}
				
				add_BigIntT(temp, one, temp2);
				copy_BigIntT(temp2, temp);
			}
			
			if (numOfExtensionsNeeded == 0)
			{
				printf("\nThe splitting field for this system contains no new extensions.\n" \
			         "Use another tool to analyse this system.\n");
			}
			
			else
			{
				//We have the number of extensions we need, now create a 
				// representative MultiVarExtT that has all the extensions
				// we need.
				
				repExtension = new_MultiVarExtT(numOfExtensionsNeeded);
				set_MultiVarExtT_mod(repExtension, bigMod);
				
				zeroPointer = calloc(numOfExtensionsNeeded, sizeof(int));
				
				//Adding the extensions we need
				for (int i = 0; i < numOfExtensionsNeeded; i += 1)
				{
					tempFactor = extract_coefficients(extDefns[i]);
					add_extension(repExtension, tempFactor, degree(extDefns[i])+1, extensionNames[i]);
					
					for (int tempCoeff = 0; tempCoeff < degree(extDefns[i])+1; tempCoeff += 1)
						tempFactor[tempCoeff] = free_BigIntT(tempFactor[tempCoeff]);
					FREE(tempFactor);
				}
				
				//Now, initialise our matrix and vector
				genericFelements = malloc(big_rows(bigF)*sizeof(void**));
				currVectElements = malloc(big_rows(bigF)*sizeof(void**));
				for (int i = 0; i < big_rows(bigF); i += 1)
				{
					genericFelements[i] = malloc(big_cols(bigF)*sizeof(void*));
					currVectElements[i] = malloc(sizeof(void*));
					for (int j = 0; j < big_cols(bigF); j += 1)
					{
						genericFelements[i][j] = new_MultiVarExtT(numOfExtensionsNeeded);
						copy_MultiVarExtT(repExtension, genericFelements[i][j]);
						set_MultiVarExtT_coefficient(genericFelements[i][j], zeroPointer, big_element(bigF, i, j));
					}
					currVectElements[i][0] = new_MultiVarExtT(numOfExtensionsNeeded);
					copy_MultiVarExtT(repExtension, currVectElements[i][0]);
				}
				
				genericF = new_MultiVarExtMatrixT(big_rows(bigF), big_cols(bigF), numOfExtensionsNeeded);
				currVect = new_MultiVarExtMatrixT(big_rows(bigF), 1, numOfExtensionsNeeded);
				tempVect = new_MultiVarExtMatrixT(big_rows(bigF), 1, numOfExtensionsNeeded);
				set_GenericMatrixT(genericF, genericFelements);
				set_GenericMatrixT(currVect, currVectElements);
				
				//Setting tempVect's extensions
				set_GenericMatrixT(tempVect, currVectElements);
				
				foundCycleLengths = calloc(1, sizeof(int));
				numOfCyclesFound  = calloc(1, sizeof(int));
				foundReps = NULL;
				
				//Write MultiVarExt definitions to output file
				if (fileOutput != NULL)
				{
					fprintmve(fileOutput, repExtension);
					fprintf(fileOutput, "\n");
					if (fclose(fileOutput) == EOF)
						fprintf(stderr, "Unable to save output file. Continuing without saving...\n");
					else
					{
						fileOutput = fopen(fileOutputName, "a");
						if (fileOutput == NULL)
							fprintf(stderr, "Unable to save output file. Continuing without saving...\n");
					}
				}
				
				//Iterate through all possible vectors over our splitting field
				while (moreVectors)
				{
					//Now, let's actually do the calculations
					generic_floyd(genericF, currVect, &theCycle);
					
					//Check to see if our found cycle length has been found before
					newCycleLength = TRUE;
					for (int cyclen = 0; (cyclen < numOfCycleLengthsFound) && (newCycleLength); cyclen += 1)
						if (foundCycleLengths[cyclen] == omega(theCycle))
						{
							newCycleLength = FALSE;
							foundCycleIndex = cyclen;
						}
						
					//Record data for new cycle found
					if (newCycleLength)
					{
						foundCycleIndex = numOfCycleLengthsFound;
						
						numOfCycleLengthsFound += 1;
						foundCycleLengths = realloc(foundCycleLengths, numOfCycleLengthsFound*sizeof(int));
						numOfCyclesFound  = realloc(numOfCyclesFound, numOfCycleLengthsFound*sizeof(int));
						
						numOfCyclesFound[foundCycleIndex] = 0;
						foundCycleLengths[foundCycleIndex] = omega(theCycle);
						
						foundReps = realloc(foundReps, numOfCycleLengthsFound*sizeof(void**));
						foundReps[numOfCycleLengthsFound-1] = NULL;
					}
					
					//Search to make sure that the cycle we calculated doesn't already have a rep
					newCycle = TRUE;
					for (int repr = 0; (repr < numOfCyclesFound[foundCycleIndex]) && newCycle; repr += 1)
					{
						//Iterate the current representative, see if our new cycle's rep is somewhere
						// in the cycle.
						copy_sim_GenericMatrixT(foundReps[foundCycleIndex][repr], currVect);
						for (int iter = 0; (iter < foundCycleLengths[foundCycleIndex]-1) && newCycle; iter += 1)
						{
							gen_mat_mul(genericF, currVect, tempVect);
							reduce_GenericMatrixT(tempVect);

							if (compare_GenericMatrixT(tempVect, rep(theCycle)))
								newCycle = FALSE;
							
							copy_sim_GenericMatrixT(tempVect, currVect);
						}
					}
					
					//Now, add the new cycle representative if needed
					if (newCycle)
					{
						numOfCyclesFound[foundCycleIndex] += 1;
						foundReps[foundCycleIndex] = realloc(foundReps[foundCycleIndex], numOfCyclesFound[foundCycleIndex]*sizeof(void*));
						foundReps[foundCycleIndex][numOfCyclesFound[foundCycleIndex]-1] = new_MultiVarExtMatrixT(big_rows(bigF), 
																																																		 1, 
																																																		 numOfExtensionsNeeded);
						//Making sure our new foundReps vector has the correct extension definitions and stuff
						set_GenericMatrixT(foundReps[foundCycleIndex][numOfCyclesFound[foundCycleIndex]-1], currVectElements);
						copy_sim_GenericMatrixT(rep(theCycle), foundReps[foundCycleIndex][numOfCyclesFound[foundCycleIndex]-1]);
						
						//Now, we print out the new cycle
						printf("Cycle #%d:\n", cycleCounter++);
						if (fileOutput != NULL)
						{
							fprintf(fileOutput, "Cycle #%d:\n", cycleCounter-1);
							if (fclose(fileOutput) == EOF)
								fprintf(stderr, "Unable to save output file. Continuing without saving...\n");
							else
							{
								fileOutput = fopen(fileOutputName, "a");
								if (fileOutput == NULL)
									fprintf(stderr, "Unable to save output file. Continuing without saving...\n");
							}
						}
						
						set_GenericMatrixT(currVect, currVectElements);
						
						printgm_row(currVect);
						printf("\n");
						if (fileOutput != NULL)
						{
							fprintgm_row(fileOutput, currVect);
							fprintf(fileOutput, "\n");
							if (fclose(fileOutput) == EOF)
								fprintf(stderr, "Unable to save output file. Continuing without saving...\n");
							else
							{
								fileOutput = fopen(fileOutputName, "a");
								if (fileOutput == NULL)
									fprintf(stderr, "Unable to save output file. Continuing without saving...\n");
							}
						}
						
						add_BigIntT(vectorCounter, one, temp);
						copy_BigIntT(temp, vectorCounter);
							
						copy_sim_GenericMatrixT(rep(theCycle), currVect);
						for (int iter = 0; iter < foundCycleLengths[foundCycleIndex]-1; iter += 1)
						{
							add_BigIntT(vectorCounter, one, temp);
							copy_BigIntT(temp, vectorCounter);
							
							gen_mat_mul(genericF, currVect, tempVect);
							reduce_GenericMatrixT(tempVect);
							
							printgm_row(tempVect);
							printf("\n");
							if (fileOutput != NULL)
							{
								fprintgm_row(fileOutput, tempVect);
								fprintf(fileOutput, "\n");
								if (fclose(fileOutput) == EOF)
									fprintf(stderr, "Unable to save output file. Continuing without saving...\n");
								else
								{
									fileOutput = fopen(fileOutputName, "a");
									if (fileOutput == NULL)
										fprintf(stderr, "Unable to save output file. Continuing without saving...\n");
								}
							}
							
							copy_sim_GenericMatrixT(tempVect, currVect);
						}
						printf("Vectors in cycles: ");
						printi(vectorCounter);
						printf("\n\n");
						if (fileOutput != NULL)
						{
							fprintf(fileOutput, "\n");
							if (fclose(fileOutput) == EOF)
								fprintf(stderr, "Unable to save output file. Continuing without saving...\n");
							else
							{
								fileOutput = fopen(fileOutputName, "a");
								if (fileOutput == NULL)
									fprintf(stderr, "Unable to save output file. Continuing without saving...\n");
							}
						}
					}
					
					//Iterate our vector
					moreVectors = FALSE;
					for (int component = 0; component < big_rows(bigF); component += 1)
					{
						if (! increment_MultiVarExtT(currVectElements[component][0]))
						{
							moreVectors = TRUE;
							set_GenericMatrixT(currVect, currVectElements);
							break;
						}
					}
				}
			}
			
			repExtension = free_MultiVarExtT(repExtension);
			
			copy_BigIntT(one, temp);
			for (int i = 1; compare_BigIntT(constant(minPolyFactors[0]), temp) >= 0; i += 1)
			{
				minPolyFactors[i] = free_BigPolyT(minPolyFactors[i]);
				add_BigIntT(temp, one, temp2);
				copy_BigIntT(temp2, temp);
			}
			minPolyFactors[0] = free_BigPolyT(minPolyFactors[0]);
			FREE(minPolyFactors);
				
			if (fileOutput != NULL)
				if (fclose(fileOutput) == EOF)
					fprintf(stderr, "Unable to save output file.\n");
			
			if (fileOutputName != NULL)
			{
				FREE(fileOutputName);
			}
			
			for (int i = 0; i < numOfExtensionsNeeded; i += 1)
				extDefns[i] = free_BigPolyT(extDefns[i]);

			FREE(extDefns);
			
			if (genericFelements != NULL)
			{
				for (int i = 0; i < big_rows(bigF); i += 1)
				{
					for (int j = 0; j < big_cols(bigF); j += 1)
						genericFelements[i][j] = free_MultiVarExtT(genericFelements[i][j]);
					
					FREE(genericFelements[i]);
				}
				FREE(genericFelements);
			}
			
			if (currVectElements != NULL)
			{
				for (int i = 0; i < big_rows(bigF); i += 1)
				{
					currVectElements[i][0] = free_MultiVarExtT(currVectElements[i][0]);
					FREE(currVectElements[i]);
				}
				FREE(currVectElements);
			}
		
			bigF = free_BigIntMatrixT(bigF);
			genericF = free_GenericMatrixT(genericF);
			currVect = free_GenericMatrixT(currVect);
			tempVect = free_GenericMatrixT(tempVect);
			
			if (zeroPointer != NULL)
			{
				FREE(zeroPointer);
			}
			
			theCycle = free_CycleInfoT(theCycle);
			FREE(foundCycleLengths);
			for (int cyclen = 0; cyclen < numOfCycleLengthsFound; cyclen += 1)
			{
				for (int cyc = 0; cyc < numOfCyclesFound[cyclen]; cyc += 1)
				{
					foundReps[cyclen][cyc] = free_GenericMatrixT(foundReps[cyclen][cyc]);
				}
				FREE(foundReps[cyclen]);
			}
			FREE(foundReps);
			FREE(numOfCyclesFound);
			
			one    = free_BigIntT(one);
			temp   = free_BigIntT(temp);
			temp2  = free_BigIntT(temp2);
			bigMod = free_BigIntT(bigMod);
			
			vectorCounter = free_BigIntT(vectorCounter);
		}
		
		
		//If we want to find a representative from each of the update matrix's orbits
		else if (! strcmp(argv[1], "orbitreps"))
		{
			BigIntTP bigMod;
			BigIntTP one;
			int oneArr[1] = {1};
			
			BigIntMatrixTP A;
			BigIntMatrixTP currVect;
			BigIntMatrixTP tempVect;
			BigIntMatrixTP tempVect2;
			
			BigIntMatrixTP** reps = NULL;
			
			int** foundCycleLengths; //[[cycle length, how many reps have cycle length]]
			BigIntTP** currVectElements;
			
			CycleInfoTP theCycle = NULL;
			
			int index;
			bool isNewCycle;
			bool hasRep;
			
			//Says whether to create the .txt output
			bool fileoutput = TRUE;
			
			char* outputfilename;
			FILE* outputFile;
			
			//If user provided modulus at command line
			if (argc > 2)
			{
				SET_BIG_NUM(argv[2], bigMod, "Unable to read modulus from command line.");
			}
			else
			{
				SET_BIG_NUM(bigintmodstring, bigMod, "Unable to read modulus from config file.");
			}
			
			//If the user wanted to change the default file output behaviour
			if (argc > 3)
				if (!strcmp(argv[3], "FALSE"))
					fileoutput = FALSE;
			
			A = read_BigIntMatrixT(updatefilepath);
			if ((A == NULL))
			{
				fprintf(stderr, "Unable to read update matrix file.\n");
				FREE_VARIABLES;
				return EXIT_FAILURE;
			}
			
			if (big_rows(A) != big_cols(A))
			{
				fprintf(stderr, "Given update matrix is not square.\n");
				FREE_VARIABLES;
				return EXIT_SUCCESS;
			}
			
			//First entry says how many cycle lengths we've found so far
			foundCycleLengths = malloc(sizeof(int*));
			foundCycleLengths[0] = calloc(2, sizeof(int));
			
			currVectElements = new_BigIntT_array(big_rows(A), 1);
			
			//Construct output file name
			outputfilename = malloc(MAXSTRLEN*sizeof(char));
			outputfilename[0] = '\0';
			strcat(outputfilename, argv[1]); //orbitreps
			strcat(outputfilename, " ");
			
			//Modulus
			if (argc > 2)
				strcat(outputfilename, argv[2]);
			else
				strcat(outputfilename, bigintmodstring);
			
			strcat(outputfilename, " F");
			
			//Matrix
			for (int row = 0; row < big_rows(A); row += 1)
				for (int col = 0; col < big_rows(A); col += 1)
					append_BigIntT(outputfilename, big_element(A, row, col));
				
			strcat(outputfilename, ".txt");
			
			currVect  = new_BigIntMatrixT(big_rows(A), 1);
			tempVect  = new_BigIntMatrixT(big_rows(A), 1);
			tempVect2 = new_BigIntMatrixT(big_rows(A), 1);
			
			one = new_BigIntT(oneArr, 1);
			
			//Increment through all possible vectors
			do
			{
				//Getting cycle length, adding to list if needed
				set_big_matrix(currVect, currVectElements);
				big_floyd(A, currVect, bigMod, &theCycle);
				
				isNewCycle = TRUE;
				for (int i = 1; i <= foundCycleLengths[0][0]; i += 1)
					if (foundCycleLengths[i][0] == omega(theCycle))
					{
						isNewCycle = FALSE;
						index = i;
						break;
					}
					
				if (isNewCycle)
				{
					//Allocate space for new cycle length info
					foundCycleLengths[0][0] += 1;
					foundCycleLengths = realloc(foundCycleLengths, (foundCycleLengths[0][0]+1)*sizeof(int*));
					foundCycleLengths[foundCycleLengths[0][0]] = calloc(2, sizeof(int));
					foundCycleLengths[foundCycleLengths[0][0]][0] = omega(theCycle);
					index = foundCycleLengths[0][0];
					
					//Also allocate space for new representatives
					reps = realloc(reps, foundCycleLengths[0][0]*sizeof(BigIntMatrixTP*));
					reps[index-1] = NULL;
				}
				
				//Now that we have our cycle length, we have to check every vector
				// in our currVect's orbit to see if we already have a representative
				if (foundCycleLengths[index][1] != 0)
				{
					copy_BigIntMatrixT(rep(theCycle), tempVect);
					hasRep = FALSE;
					
					//Now, check tempVect's orbit for reps
					for (int i = 0; i < omega(theCycle); i += 1)
					{
						for (int r = 0; r < foundCycleLengths[index][1]; r += 1)
							if (compare_BigIntMatrixT(tempVect, reps[index-1][r]))
							{
								hasRep = TRUE;
								break;
							}
							
						if (hasRep)
							break;
						
						big_mat_mul(A, tempVect, tempVect2);
						modbm(tempVect2, bigMod);
						copy_BigIntMatrixT(tempVect2, tempVect);
					}
					
					//If we found a new representative, add it to the list
					if (!hasRep)
					{
						foundCycleLengths[index][1] += 1;
						reps[index-1] = realloc(reps[index-1], foundCycleLengths[index][1]*sizeof(BigIntMatrixTP));
						
						reps[index-1][foundCycleLengths[index][1]-1] = new_BigIntMatrixT(big_rows(A), 1);
						copy_BigIntMatrixT(rep(theCycle), reps[index-1][foundCycleLengths[index][1]-1]);
					}
				}
			
				//If we don't have any representatives yet
				else
				{
					//Increment number of reps we've found
					foundCycleLengths[index][1] += 1;
					reps[index-1] = realloc(reps[index-1], foundCycleLengths[index][1]*sizeof(BigIntMatrixTP));
					
					reps[index-1][0] = new_BigIntMatrixT(big_rows(A), 1);
					copy_BigIntMatrixT(rep(theCycle), reps[index-1][0]);
				}
			}
			while (! increment_BigIntT_array(currVectElements, big_rows(A), 1, one, bigMod));
			
			//Create save file
			if (fileoutput)
			{
				outputFile = fopen(outputfilename, "w");
				if (outputFile == NULL)
					fprintf(stderr, "Unable to create save file. Results will still print to stdout.\n");
			}
			else
				outputFile = NULL;
			
			//Now, print out the representatives!
			for (int cyclenindex = 1; cyclenindex <= foundCycleLengths[0][0]; cyclenindex += 1)
			{
				printf("Cycle length: %d\n", foundCycleLengths[cyclenindex][0]);
				printf("Representatives:\n");
				
				if (outputFile != NULL)
				{
					fprintf(outputFile, "Cycle length: %d\n", foundCycleLengths[cyclenindex][0]);
					fprintf(outputFile, "Reps: \n");
				}
				
				for (int cycrep = 0; cycrep < foundCycleLengths[cyclenindex][1]; cycrep += 1)
				{
					printbm_row(reps[cyclenindex-1][cycrep]);
					printf("\n");
					
					if (outputFile != NULL)
					{						
						fprintbm_row(outputFile, reps[cyclenindex-1][cycrep]);
						fprintf(outputFile, "\n");
					}
				}
				
				if (outputFile != NULL)
					fprintf(outputFile, "\n");
			}
			
			if (fileoutput)
				if (fclose(outputFile) == EOF)
					fprintf(stderr, "Unable to save data.\n");
			
			
			currVectElements = free_BigIntT_array(currVectElements, big_rows(A), 1);
			
			for (int i = 1; i <= foundCycleLengths[0][0]; i += 1)
			{
				for (int j = 0; j < foundCycleLengths[i][1]; j += 1)
					reps[i-1][j] = free_BigIntMatrixT(reps[i-1][j]);
				
				FREE(reps[i-1]);
				FREE(foundCycleLengths[i]);
			}
			FREE(reps);

			FREE(foundCycleLengths[0]);
			FREE(foundCycleLengths);
			
			bigMod = free_BigIntT(bigMod);
			one    = free_BigIntT(one);
			
			theCycle = free_CycleInfoT(theCycle);
			
			A         = free_BigIntMatrixT(A);
			currVect  = free_BigIntMatrixT(currVect);
			tempVect  = free_BigIntMatrixT(tempVect);
			tempVect2 = free_BigIntMatrixT(tempVect2);
			
			FREE(outputfilename);
		}
		
		
		//If we want to generate transient length representatives
		else if (!strcmp(argv[1], "branchreps"))
		{
			BigIntMatrixTP A;
			BigIntMatrixTP I;
			BigIntMatrixTP currVect  = NULL;
			BigIntMatrixTP tempVect  = NULL;
			BigIntMatrixTP tempVect2 = NULL;
			BigIntMatrixTP tempVect3 = NULL;
			BigIntMatrixTP tempRoot  = NULL;
			BigIntTP**     currVectElements;
			
			CycleInfoTP theCycle = NULL;
			
			/* First vector is the "leaf" of the branch, the vector with more steps to get into a cycle
			 * Second vector is the "root" of the branch, the vector closest to a cycle (and which 
			 *  possibly leads to another branch).
			 * Third vector is the "base root" of the branch, the root you'd get to if you follow the branch down
			 */
			BigIntMatrixTP** branches; //Holds the transient branch representations
			int numOfBranches = 0;
			
			//Holds references to the vectors we're using
			//Helps to save on memory usage since there's a lot of redundancy in our tree representation
			BigIntMatrixTP* vectCatalogue = NULL;
			int sizeOfCatalogue = 0;
			
			BigIntTP bigMod;
			BigIntTP one;
			int oneArr[1] = {1};
			
			bool isNewRoot;
			bool foundGrowth; //If we've found the place where a vector grows from
			bool vectInsideVine;
			bool hasConnections = TRUE;
			bool isTopLeaf;   //For determining which leaves are strictly necessary for transient length representatives
			
			char* fileName = NULL;
			FILE* outputFile = NULL;
			
			if (argc > 2)
			{
				SET_BIG_NUM(argv[2], bigMod, "Unable to read modulus from command line.");
			}
			else
			{
				SET_BIG_NUM(bigintmodstring, bigMod, "Unable to read modulus fom command line.");
			}
			
			A = read_BigIntMatrixT(updatefilepath);
			if (A == NULL)
			{
				fprintf(stderr, "Unable to read matrix in %s.\n", updatefilepath);
				FREE_VARIABLES;
				return EXIT_FAILURE;
			}
			
			//Checking to ensure a transient region exists
			I = identity_BigIntMatrixT(big_rows(A));
			big_floyd(A, I, bigMod, &theCycle);
			if (tau(theCycle) == 0)
			{
				printf("Given matrix is invertible under given modulus, so no transient region exists.\n");
				
				bigMod = free_BigIntT(bigMod);
				A = free_BigIntMatrixT(A);
				I = free_BigIntMatrixT(I);
				theCycle = free_CycleInfoT(theCycle);
				FREE_VARIABLES;
				return EXIT_SUCCESS;
			}
			
			fileName = malloc(MAXSTRLEN*sizeof(char));
			fileName[0] = '\0';
			strcat(fileName, "branchreps ");
			append_BigIntT(fileName, bigMod);
			strcat(fileName, " F");
			for (int row = 0; row < big_rows(A); row += 1)
				for (int col = 0; col < big_cols(A); col += 1)
					append_BigIntT(fileName, big_element(A, row, col));
			strcat(fileName, ".txt");
			
			outputFile = fopen(fileName, "a");
			if (outputFile == NULL)
				fprintf(stderr, "Unable to create output file. Continuing without saving...\n");
			
			one = new_BigIntT(oneArr, 1);
			
			branches    = malloc(sizeof(BigIntMatrixTP*));
			branches[0] = malloc(3*sizeof(BigIntMatrixTP));
			
			currVectElements = new_BigIntT_array(big_rows(A), 1);
			
			currVect  = new_BigIntMatrixT(big_rows(A), 1);
			tempVect  = new_BigIntMatrixT(big_rows(A), 1);
			tempVect2 = new_BigIntMatrixT(big_rows(A), 1);
			tempVect3 = new_BigIntMatrixT(big_rows(A), 1);
			tempRoot  = new_BigIntMatrixT(big_rows(A), 1);
			
			printf("Modulus: ");
			printi(bigMod);
			printf("\nMatrix:\n");
			printbm(A);
			printf("~~~~~~~~~~\n");
			
			//Loop over each vector, check if it's in a cycle or not
			do
			{
				//This is here to prevent the program from counting new branches after we blocked off a single leaf
				// from adding to the branch count
				hasConnections = TRUE; //Says whether the leaf we're looking at has other connections feeding into it
				
				set_big_matrix(currVect, currVectElements);
				big_floyd(A, currVect, bigMod, &theCycle);
				
				//If we found a vector that isn't a part of a cycle
				if (tau(theCycle) != 0)
				{
					//If we found our first branch
					if (numOfBranches == 0)
					{
						//printf("First branch.\n");
						
						//Insert vector into leaf
						branches[0][0] = BigIntMatrixT_catalogue_get(&vectCatalogue, &sizeOfCatalogue, currVect);
						
						//Get root (which we already know since this is the only branch and we have the transient length)
						copy_BigIntMatrixT(currVect, tempVect);
						for (int t = 0; t < tau(theCycle); t += 1)
						{
							big_mat_mul(A, tempVect, tempVect2);
							modbm(tempVect2, bigMod);
							copy_BigIntMatrixT(tempVect2, tempVect);
						}
						branches[0][1] = BigIntMatrixT_catalogue_get(&vectCatalogue, &sizeOfCatalogue, tempVect);
						
						//Get base root
						branches[0][2] = branches[0][1];
					}
					
					else
					{
						//Find the root of our current vector (acting as a leaf)
						copy_BigIntMatrixT(currVect, tempRoot);
						for (int t = 0; t < tau(theCycle); t += 1)
						{
							big_mat_mul(A, tempRoot, tempVect);
							modbm(tempVect, bigMod);
							copy_BigIntMatrixT(tempVect, tempRoot);
						}
						
						//Check to see if we've already found this root
						isNewRoot = TRUE;
						for (int r = 0; r < numOfBranches; r += 1)
							if (compare_BigIntMatrixT(branches[r][2], tempRoot))
							{
								isNewRoot = FALSE;
								break;
							}
							
						//Easy case. We just create a new branch like before
						if (isNewRoot)
						{
							/*printf("Has a new root:\n");
							printbm(tempRoot);
							printf("\n");*/
							
							branches[numOfBranches][0] = BigIntMatrixT_catalogue_get(&vectCatalogue, &sizeOfCatalogue, currVect);
							branches[numOfBranches][1] = BigIntMatrixT_catalogue_get(&vectCatalogue, &sizeOfCatalogue, tempRoot);
							branches[numOfBranches][2] = branches[numOfBranches][1];
						}
						
						//Now, we search through the existing branches to find where currVect grows from
						else
						{
							foundGrowth = FALSE;   //Says whether we found a place to attach currVect onto
							vectInsideVine = TRUE; //Says whether the vect is intrinsically inside another branch, meaning we can ignore it
							
							copy_BigIntMatrixT(currVect, tempVect);
							while (!foundGrowth)
							{
								for (int twig = 0; twig < numOfBranches; twig += 1)
								{
									//If we found a branch with the same base root
									if (compare_BigIntMatrixT(branches[twig][2], tempRoot))
									{
										//Is tempVect the root of our current branch?
										if (compare_BigIntMatrixT(tempVect, branches[twig][1]))
										{
											foundGrowth = TRUE;
											if (!vectInsideVine)
											{
												/*printf("Vector connects to base of branch:\n");
												printbm(branches[twig][1]);
												printf("\n");*/
												
												//Create new branch which feeds into the root of our current one
												branches[numOfBranches][0] = BigIntMatrixT_catalogue_get(&vectCatalogue, &sizeOfCatalogue, currVect);
												branches[numOfBranches][1] = BigIntMatrixT_catalogue_get(&vectCatalogue, &sizeOfCatalogue, branches[twig][1]);
												branches[numOfBranches][2] = BigIntMatrixT_catalogue_get(&vectCatalogue, &sizeOfCatalogue, branches[twig][2]);
											}
											else
											{
												//printf("Vector is inside the branch.\n");
												hasConnections = FALSE;
											}
											
											break; //Don't need to do anything if our vector is inside the branch
										}
										
										//Search through entire branch, checking if our tempVect can be found
										copy_BigIntMatrixT(branches[twig][0], tempVect2);
										while (!compare_BigIntMatrixT(tempVect2, branches[twig][1]))
										{
											//If we found where our vector feeds into on the current branch
											if (compare_BigIntMatrixT(tempVect2, tempVect))
											{
												foundGrowth = TRUE;
												
												//If our currVect connects to the leaf
												if (compare_BigIntMatrixT(tempVect2, branches[twig][0]))
												{
													hasConnections = FALSE;
													
													//Check to see if anything else connects to the leaf
													for (int br = 0; br < numOfBranches; br += 1)
													{
														if (compare_BigIntMatrixT(branches[br][1], tempVect2))
														{
															hasConnections = TRUE;
															break;
														}
													}
													
													//If nothing else connects to the leaf, extend the branch so that
													// currVect becomes the leaf
													if (!hasConnections)
													{
														/*printf("Vector connects to a single leaf:\n");
														printbm(branches[twig][0]);
														printf("\n");*/
														branches[twig][0] = BigIntMatrixT_catalogue_get(&vectCatalogue, &sizeOfCatalogue, currVect);
														break;
													}
												}
												
												//If currVect connects somehwere in the middle of the branch,
												// and isn't INSIDE the branch
												if (!vectInsideVine)
												{
													/*printf("Vector connects to branch at:\n");
													printbm(tempVect2);
													printf("\n");*/
													
													//New branch with currVect
													branches[numOfBranches][0] = BigIntMatrixT_catalogue_get(&vectCatalogue, &sizeOfCatalogue, currVect);
													branches[numOfBranches][1] = BigIntMatrixT_catalogue_get(&vectCatalogue, &sizeOfCatalogue, tempVect);
													branches[numOfBranches][2] = BigIntMatrixT_catalogue_get(&vectCatalogue, &sizeOfCatalogue, branches[twig][2]);
													
													numOfBranches += 1;
													branches = realloc(branches, (numOfBranches+1)*sizeof(BigIntMatrixTP*));
													branches[numOfBranches] = malloc(3*sizeof(BigIntMatrixTP));
													
													//Bottom split of old branch
													branches[numOfBranches][0] = BigIntMatrixT_catalogue_get(&vectCatalogue, &sizeOfCatalogue, branches[numOfBranches-1][1]);
													branches[numOfBranches][1] = BigIntMatrixT_catalogue_get(&vectCatalogue, &sizeOfCatalogue, branches[twig][1]);
													branches[numOfBranches][2] = BigIntMatrixT_catalogue_get(&vectCatalogue, &sizeOfCatalogue, branches[twig][2]);
													
													//Top split
													branches[twig][1] = BigIntMatrixT_catalogue_get(&vectCatalogue, &sizeOfCatalogue, branches[numOfBranches][0]);
													
													//Split old branch into two branches
													//We want our new branch's root to be a leaf for one of the split parts
												}
												else
												{
													//printf("Vector is inside the branch.\n");
													hasConnections = FALSE;
												}
												
												break;
											}
											
											//Iterate the vector in the branch we're checking
											big_mat_mul(A, tempVect2, tempVect3);
											modbm(tempVect3, bigMod);
											copy_BigIntMatrixT(tempVect3, tempVect2);
										}
									}
									
									if (foundGrowth)
										break;
								}
								
								//Iterate our vector, check again
								big_mat_mul(A, tempVect, tempVect2);
								modbm(tempVect2, bigMod);
								copy_BigIntMatrixT(tempVect2, tempVect);
								vectInsideVine = FALSE;
							}
						}
					}
					
					//Increase our branch count
					if (hasConnections)
					{
						numOfBranches += 1;
						branches = realloc(branches, (numOfBranches+1)*sizeof(BigIntMatrixTP*));
						branches[numOfBranches] = malloc(3*sizeof(BigIntMatrixTP));
					}
				}
			}
			while (!increment_BigIntT_array(currVectElements, big_rows(A), 1, one, bigMod));
			
			//Print out our branches
			/*
			printf("\n~~~~~Branches:~~~~~\n");
			for (int br = 0; br < numOfBranches; br += 1)
			{
				printf("Leaf:\n");
				printbm_row(branches[br][0]);
				printf("\nRoot:\n");
				printbm_row(branches[br][1]);
				printf("\nBase root:\n");
				printbm_row(branches[br][2]);
				printf("\n~~~~~~~~~~~~~~~~~~~\n");
				getchar();
			}
			*/
			
			//Now, we simply figure out which leaves are "true leaves",
			// that is, leaves that aren't the roots of any branch
			printf("Leaves of the transient tree:\n");
			for (int maybeLeaf = 0; maybeLeaf < numOfBranches; maybeLeaf += 1)
			{
				isTopLeaf = TRUE;
				for (int rootCheck = 0; rootCheck < numOfBranches; rootCheck += 1)
					if (compare_BigIntMatrixT(branches[maybeLeaf][0], branches[rootCheck][1]))
					{
						isTopLeaf = FALSE;
						break;
					}
					
				if (isTopLeaf)
				{
					if (vectorType == row)
					{
						printbm_row(branches[maybeLeaf][0]);
						if (outputFile != NULL)
							fprintbm_row(outputFile, branches[maybeLeaf][0]);
					}
					else
					{
						printbm(branches[maybeLeaf][0]);
						if (outputFile != NULL)
							fprintbm(outputFile, branches[maybeLeaf][0]);
					}
					
					printf("\n");
					if (outputFile != NULL)
						fprintf(outputFile, "\n");
				}
			}
			
			
			for (int i = 0; i <= numOfBranches; i += 1)
			{				
				FREE(branches[i]);
			}
			FREE(branches);
			
			for (int i = 0; i < sizeOfCatalogue; i += 1)
				vectCatalogue[i] = free_BigIntMatrixT(vectCatalogue[i]);
			FREE(vectCatalogue);
			
			free_BigIntT_array(currVectElements, big_rows(A), 1);
			
			one    = free_BigIntT(one);
			bigMod = free_BigIntT(bigMod);
			
			A         = free_BigIntMatrixT(A);
			I         = free_BigIntMatrixT(I);
			currVect  = free_BigIntMatrixT(currVect);
			tempVect  = free_BigIntMatrixT(tempVect);
			tempVect2 = free_BigIntMatrixT(tempVect);
			tempVect3 = free_BigIntMatrixT(tempVect);
			tempRoot  = free_BigIntMatrixT(tempRoot);
			theCycle  = free_CycleInfoT(theCycle);
			
			FREE(fileName);
			if (outputFile != NULL)
				if (fclose(outputFile) == EOF)
					fprintf(stderr, "Unable to save output file.\n");
		}
		
		
		//If the user wants to calculate each unique cycle space
		// within a given LCA system
		else if (! strcmp(argv[1], "orbitspaces"))
		{
			bool calcMinPoly = FALSE; //Do we calculate min polys for each cyclespace?
			bool newCyclespace;
			
			int numArr[1] = {1};
			BigIntTP one;
			
			BigIntTP bigMod;
			BigIntTP** currVectElements;
			
			BigIntMatrixTP A;
			BigIntMatrixTP currVect;
			BigIntMatrixTP tempVect;
			BigIntMatrixTP tempVect2;
			
			BigIntTP** tempSetCyclespace;
			BigIntMatrixTP tempCyclespace;
			BigIntMatrixTP tempCyclespace2;
			
			int numOfCyclespaces = 0;
			BigIntMatrixTP* cyclespaces = NULL;

			BigPolyTP   charaPoly   = NULL;
			BigPolyTP** minpolys    = NULL;
			
			if (argc > 2)
			{
				SET_BIG_NUM(argv[2], bigMod, "Unable to read modulus from command line.");
			}
			else
			{
				SET_BIG_NUM(bigintmodstring, bigMod, "Unable to read modulus from config file.");
			}
			
			//If the user wants to turn on minimal polynomial calculation
			if (argc > 3)
				if (! strcmp(argv[3], "TRUE"))
					calcMinPoly = TRUE;
			
			A = read_BigIntMatrixT(updatefilepath);
			if (A == NULL)
			{
				fprintf(stderr, "Unable to read update matrix from config file.\n");
				bigMod = free_BigIntT(bigMod);
				FREE_VARIABLES;
				return EXIT_FAILURE;
			}
			
			if (big_rows(A) != big_cols(A))
			{
				fprintf(stderr, "Given update matrix is not square.\n");
				A = free_BigIntMatrixT(A);
				bigMod = free_BigIntT(bigMod);
				FREE_VARIABLES;
				return EXIT_SUCCESS;
			}
			
			if (calcMinPoly)
				charaPoly = chara_poly(A, bigMod);
			
			printf("Matrix:\n");
			printbm(A);
			printf("Modulus: ");
			printi(bigMod);
			printf("\n\n");
			
			one = new_BigIntT(numArr, 1);
			currVect  = new_BigIntMatrixT(big_rows(A), 1);
			tempVect  = new_BigIntMatrixT(big_rows(A), 1);
			tempVect2 = new_BigIntMatrixT(big_rows(A), 1);
			currVectElements = new_BigIntT_array(big_rows(A), 1);
			
			tempCyclespace  = new_BigIntMatrixT(big_rows(A), big_cols(A));
			tempCyclespace2 = new_BigIntMatrixT(big_rows(A), big_cols(A));
			tempSetCyclespace = new_BigIntT_array(big_rows(A), big_cols(A));
			
			do
			{
				set_big_matrix(currVect, currVectElements);
				copy_BigIntMatrixT(currVect, tempVect);
				clear_BigIntT_array(tempSetCyclespace, big_rows(A), big_cols(A));
				
				//Calculate the max number of possible linearly-independent vectors
				// for the cyclespace, and put them in tempSetCyclespace
				for (int entry = 0; entry < big_rows(A); entry += 1)
					copy_BigIntT(currVectElements[entry][0], tempSetCyclespace[0][entry]);
				
				for (int v = 1; v < big_rows(A); v += 1)
				{
					big_mat_mul(A, tempVect, tempVect2);
					copy_BigIntMatrixT(tempVect2, tempVect);
					modbm(tempVect, bigMod);
					
					for (int entry = 0; entry < big_rows(A); entry += 1)
						copy_BigIntT(big_element(tempVect, entry, 0), tempSetCyclespace[v][entry]);
				}
				
				/*printbm_row(currVect);
				printf("\n"); */
				set_big_matrix(tempCyclespace, tempSetCyclespace);
				/*printbm(tempCyclespace);
				printf("~~~\n");*/
				big_row_echelon(tempCyclespace, bigMod, tempCyclespace2, NULL);
				big_reduced_row_echelon(tempCyclespace2, bigMod, tempCyclespace, NULL);
				//printbm(tempCyclespace);
				
				//Now, search to see if we've already accounted for the current cyclespace
				newCyclespace = TRUE;
				for (int m = 0; m < numOfCyclespaces; m += 1)
					if (compare_BigIntMatrixT(tempCyclespace, cyclespaces[m]))
					{
						newCyclespace = FALSE;
						break;
					}
				
				if (newCyclespace)
				{
					numOfCyclespaces += 1;
					cyclespaces = realloc(cyclespaces, numOfCyclespaces*sizeof(BigIntMatrixTP));
					cyclespaces[numOfCyclespaces-1] = new_BigIntMatrixT(big_rows(A), big_cols(A));
					copy_BigIntMatrixT(tempCyclespace, cyclespaces[numOfCyclespaces-1]);
					
					//Calculate min poly if needed
					if (calcMinPoly)
					{
						minpolys = realloc(minpolys, numOfCyclespaces*sizeof(BigPolyTP*));
						minpolys[numOfCyclespaces-1] = min_poly(A, currVect, bigMod, charaPoly);
					}
				}
			}
			while (! increment_BigIntT_array(currVectElements, big_rows(A), 1, one, bigMod));
			
			
			//Now, we print out all the cyclespaces we've found
			printf("Unique cyclespaces found:\n");
			for (int m = 0; m < numOfCyclespaces; m += 1)
			{
				big_rowsp(cyclespaces[m]);
				
				//If needed, also print out the min poly associated with the cyclespace
				if (calcMinPoly)
				{
					printf(" : ");
					old_printpf(minpolys[m]);
					minpolys[m] = free_BigPolyT_factors(minpolys[m]);
				}
				
				printf("\n");
				cyclespaces[m] = free_BigIntMatrixT(cyclespaces[m]);
			}
			FREE(cyclespaces);
			if (minpolys != NULL)
			{
				FREE(minpolys);
			}
			
			one    = free_BigIntT(one);
			bigMod = free_BigIntT(bigMod);
			
			charaPoly   = free_BigPolyT(charaPoly);
			
			currVectElements  = free_BigIntT_array(currVectElements, big_rows(A), 1);
			tempSetCyclespace = free_BigIntT_array(tempSetCyclespace, big_rows(A), big_cols(A));
			
			A         = free_BigIntMatrixT(A);
			currVect  = free_BigIntMatrixT(currVect);
			tempVect  = free_BigIntMatrixT(tempVect);
			tempVect2 = free_BigIntMatrixT(tempVect2);
			
			tempCyclespace  = free_BigIntMatrixT(tempCyclespace);
			tempCyclespace2 = free_BigIntMatrixT(tempCyclespace2);
		}
		
		
		//If the user wants to use Floyd's Cycle Detection Algorithm
		else if (! strcmp(argv[1], "floyd"))
		{
			BigIntTP bigModulus;
			BigIntMatrixTP initial;
			BigIntMatrixTP update;
			CycleInfoTP coolCycle = NULL;
			
			//If the user specified a custom modulus
			if (argc > 2)
			{
				SET_BIG_NUM(argv[2], bigModulus, "Unable to read modulus from command line.");
			}
			
			//Use the modulus provided in the .config file
			else
			{
				SET_BIG_NUM(bigintmodstring, bigModulus, "Unable to read modulus from config file. Try giving the modulus as a CLI argument.");
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
			{
				printcycle(coolCycle, vectorType);
				if ((vectorType == row) && (big_rows(initial) != big_cols(initial)))
					printf("\n");
			}
			
			
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
				printf(ANSI_COLOR_YELLOW "cycmatsearch " ANSI_COLOR_CYAN "resume size maxmod cycles..." ANSI_COLOR_RESET \
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
			
			BigIntTP temp;
			BigIntTP temp2;
			BigIntTP tempModCounter;
			BigIntTP tempPercentCounter;
			
			BigIntTP progressUpdateElement; //For keeping track of progress through a modulus
			int progressRow = 2;
			int progressCol = 1; //Used to decide which element in the matrix to look at for progress updates
			bool printProgress = TRUE;
			
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
			cycmatsearch 2 30 6 10
			cycmatsearch 3 6 2 3 5
			cycmatsearch 3 6 2 2 3
			cycmatsearch 2 30 6 15
			cycmatsearch 3 6 6 10 14
			
			cycmatsearch 3 9/9 6 10 14 (54/81)
			
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
			
			currMatElements = new_BigIntT_array(size, size);
			
			//If we need to resume computation from a specific matrix
			if (! strcmp(argv[2], "TRUE"))
			{
				currMat = read_BigIntMatrixT(resumefilepath);
				SET_BIG_NUM(bigintmodstring, currMod, "Unable to read resume modulus from config file.");
				
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
			SET_BIG_NUM(argv[4], maxMod, "Unable to read modulus from command line.");
			
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
			
			one   = new_BigIntT(oneArr, 1);
			zero  = empty_BigIntT(1);
			temp  = empty_BigIntT(1);
			temp2 = empty_BigIntT(1);
			
			tempModCounter        = empty_BigIntT(1);
			tempPercentCounter    = empty_BigIntT(1);
			progressUpdateElement = empty_BigIntT(1);
			
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
			
			textOutput = fopen(textOutputName, "w");
			if (textOutput == NULL)
			{
				fprintf(stderr, "Unable to create/open %s for writing.\n", textOutputName);
				return EXIT_FAILURE;
			}
			
			//Prepare to monitor whether the progress element in our matrix increased or not
			copy_BigIntT(currMatElements[progressRow][progressCol], progressUpdateElement);

			
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
					
					if ((printProgress) && 
					    (compare_BigIntT(currMatElements[progressRow][progressCol], progressUpdateElement) != 0))
					{
						//Calculate proper fraction to show for percentage
						copy_BigIntT(one, tempModCounter);
						copy_BigIntT(zero, tempPercentCounter);
						for (int r = progressRow; r < size; r += 1)
						{
							for (int c = 0; c < size; c += 1)
							{
								//Making sure we start from the correct element in the matrix
								if ((r == progressRow) && (c == 0))
									c = progressCol;
								
								multiply_BigIntT(currMatElements[r][c], tempModCounter, temp);
								add_BigIntT(tempPercentCounter, temp, temp2);
								copy_BigIntT(temp2, tempPercentCounter);
								
								//Increment modCounter
								multiply_BigIntT(tempModCounter, currMod, temp);
								copy_BigIntT(temp, tempModCounter);
							}
						}
						
						printi(tempPercentCounter);
						printf(" / ");
						printi(tempModCounter);
						printf(" searched...\n");
						
						copy_BigIntT(currMatElements[progressRow][progressCol], progressUpdateElement);
					}
					
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
					checkedAllMatrices = increment_BigIntT_array(currMatElements, size, size, one, currMod);
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
			temp2   = free_BigIntT(temp2);
			
			tempModCounter = free_BigIntT(tempModCounter);
			tempPercentCounter = free_BigIntT(tempPercentCounter);
			progressUpdateElement = free_BigIntT(progressUpdateElement);
			
			currMat  = free_BigIntMatrixT(currMat);
			zeroMat  = free_BigIntMatrixT(zeroMat);
			tempMat  = free_BigIntMatrixT(tempMat);
			tempMat2 = free_BigIntMatrixT(tempMat2);
			
			theCycle = free_CycleInfoT(theCycle);
			
			FREE(colVectCycles);
			FREE(textOutputName);
			currMatElements = free_BigIntT_array(currMatElements, size, size);
			
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
				SET_BIG_NUM(argv[4], bigMod, "Unable to read modulus from command line.");
			}
			
			else
			{
				SET_BIG_NUM(bigintmodstring, bigMod, "Unable to read modulus from config file.");
			}
			
			SET_BIG_NUM(argv[2], from, "Unable to read starting cycle length from command line.");
			SET_BIG_NUM(argv[3], to, "Unable to read conversion cycle length from command line.");
			
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
			BigIntTP   prevLastElement; //Holds the last element in the matrix; for debugging
			BigIntTP*  cycleLengthFactors;
			BigIntTP** currMatElements;
			BigIntTP** currVectElements;
			
			BigIntMatrixTP currMat;
			BigIntMatrixTP tempCCM;
			BigIntMatrixTP zeroMat;
			BigIntMatrixTP resumeMat = NULL; //The matrix to resume at if specified
			
			BigIntMatrixTP currVect;
			
			CycleInfoTP theCycle = NULL;
			
			int matSize;
			int indexCounter; //For freeing and printing
			int oneArr[1] = {1};
			
			//Holds the cycle length of each matrix and vector we check
			int cycleLengthArray[1] = {0};
			BigIntTP bigOmega;
			BigIntTP bigVectOmega;
			
			//Holds the number of vectors with each possible cycle length
			int* cycleLengthCounts;
			
			bool checkedAllMatrices = FALSE;
			bool checkedAllVectors  = FALSE;
			bool hasAllZeros        = FALSE;
			
			bool hasInterestingCycles = FALSE;
			
			bool debug = TRUE;
			int matrixCounter = 0;
			
			char* outputfilename;
			FILE* outputFile;
			BigPolyTP  charaPoly;
			BigPolyTP* charaPolyFactors = NULL;
			
			//Checking to see if the user provided a modulus on the command line
			if (argc > 4)
			{
				SET_BIG_NUM(argv[4], bigMod, "Unable to read modulus on command line.");
			}
			
			else
			{
				SET_BIG_NUM(bigintmodstring, bigMod, "Unable to read modulus from config file.");
			}
			
			matSize = (int)strtol(argv[3], &tempStr, 10);
			if (tempStr[0] != '\0')
			{
				fprintf(stderr, "Invalid matrix size provided on command line.\n");
				FREE_VARIABLES;
				return EXIT_FAILURE;
			}
			
			//If the user wants to use the resume matrix
			if (! strcmp(argv[2], "TRUE"))
			{
				resumeMat = read_BigIntMatrixT(resumefilepath);
				if (resumeMat == NULL)
				{
					fprintf(stderr, "Unable to read resume matrix specified in config file. Computation " \
					"will start from the zero matrix.\n");
				}
			}
			
			//Create output name
			outputfilename = malloc(MAXSTRLEN*sizeof(char));
			outputfilename[0] = '\0';
			
			strcat(outputfilename, argv[1]); //ccmzerosearch
			strcat(outputfilename, " ");
			strcat(outputfilename, argv[3]); //matrix size
			strcat(outputfilename, " ");
			
			//modulus
			if (argc > 4)
				strcat(outputfilename, argv[4]);
			else
				strcat(outputfilename, bigintmodstring);
			
			strcat(outputfilename, ".txt");
			outputFile = fopen(outputfilename, "w");
			if (outputFile == NULL)
				fprintf(stderr, "Unable to save found matrices. Continuing without saving...\n");
			else
				printf("Found matrices will be saved at %s\n", outputfilename);
			
			//Initialising elements for matrices and vectors
			currMatElements  = new_BigIntT_array(matSize, matSize);
			currVectElements = new_BigIntT_array(matSize, 1);
			
			zero = empty_BigIntT(1);
			one  = new_BigIntT(oneArr, 1);
			currMat  = new_BigIntMatrixT(matSize, matSize);
			zeroMat  = new_BigIntMatrixT(matSize, matSize);
			currVect = new_BigIntMatrixT(matSize, 1);
			
			temp  = empty_BigIntT(1);
			temp2 = empty_BigIntT(1);
			prevLastElement = empty_BigIntT(1);
			
			//Allowing user to resume computation
			if (resumeMat != NULL)
				for (int i = 0; i < matSize; i += 1)
					for (int j = 0; j < matSize; j += 1)
						copy_BigIntT(big_element(resumeMat, i, j), currMatElements[i][j]);
			
			//Loop until we're checked every matrix under the given modulus
			while (!checkedAllMatrices)
			{
				set_big_matrix(currMat, currMatElements);
				big_floyd(currMat, currMat, bigMod, &theCycle);
				
				if (debug)
				{
					if (compare_BigIntT(currMatElements[matSize-1][matSize-2], prevLastElement) != 0)
					{
							printf("Current matrix:\n");
							printbm(currMat);
					}
				}
				
				//This program currently assumes that the cycle length of the
				// matrix will be less than MAXBUNCH
				cycleLengthArray[0] = omega(theCycle);
				bigOmega = new_BigIntT(cycleLengthArray, 1);
				cycleLengthFactors = divisors_of_BigIntT(bigOmega);
				
				//Check each factor to see what its CCM looks like
				indexCounter = 1;
				copy_BigIntT(zero, temp);
				hasAllZeros = TRUE;
				while (compare_BigIntT(temp, cycleLengthFactors[0]) < 0)
				{
					tempCCM = new_BigIntMatrixT(matSize, matSize);
					ccm(currMat, tempCCM, bigOmega, cycleLengthFactors[indexCounter], bigMod);
					
					//If the CCM is not the zero matrix, we can stop looking
					if (!compare_BigIntMatrixT(tempCCM, zeroMat))
					{
						tempCCM = free_BigIntMatrixT(tempCCM);
						hasAllZeros = FALSE;
						break;
					}
					
					indexCounter += 1;
					add_BigIntT(temp, one, temp2);
					copy_BigIntT(temp2, temp);
					tempCCM = free_BigIntMatrixT(tempCCM);
				}
				
				//Calculate cycle length counts so we can determine if the CCM
				// zeros are justified or not
				cycleLengthCounts = calloc(indexCounter, sizeof(int));
				
				checkedAllVectors = FALSE;
				while (!checkedAllVectors)
				{
					set_big_matrix(currVect, currVectElements);
					big_floyd(currMat, currVect, bigMod, &theCycle);

					cycleLengthArray[0] = omega(theCycle);
					bigVectOmega = new_BigIntT(cycleLengthArray, 1);
					
					//First, check if vector cycle length is equal to matrix
					//This prevents possible index errors below
					if (compare_BigIntT(bigVectOmega, bigOmega) == 0)
						cycleLengthCounts[indexCounter-1] += 1;
					
					else
					{
						//Find which cycle length group this vector belongs to
						for (int w = 1; w <= indexCounter-1; w += 1)
							if (compare_BigIntT(cycleLengthFactors[w], bigVectOmega) == 0)
							{
								cycleLengthCounts[w-1] += 1;
								break;
							}
					}

					//Iterate to next vector
					checkedAllVectors = increment_BigIntT_array(currVectElements, matSize, 1, one, bigMod);
				}
				
				
				//Now, check to see if any "interesting cycle lengths" exist for the matrix
				// This is, cycle lengths that aren't maximal or zero
				hasInterestingCycles = FALSE;
				for (int c = 0; c < indexCounter-1; c += 1)
				{
					//We only care about fixed points if there's more than just the zero vector
					if (c == 0)
					{
						if (cycleLengthCounts[0] > 1)
						{
							hasInterestingCycles = TRUE;
							break;
						}
					}
					
					else
					{
						if (cycleLengthCounts[c] > 0)
						{
							hasInterestingCycles = TRUE;
							break;
						}
					}
				}
				
				//Print out the matrix we found, as well as any relevant information regarding it
				if ((hasAllZeros) && (compare_BigIntT(one, bigOmega) != 0) && (hasInterestingCycles))
				{
					matrixCounter += 1;
					printbm(currMat);
					printf("Cycle length counts: ");
					
					if (outputFile != NULL)
					{
						fprintbm_nopad(outputFile, currMat);
						fprintf(outputFile, "Cycle length counts: ");
					}
					
					indexCounter = 1;
					copy_BigIntT(zero, temp);
					while (compare_BigIntT(temp, cycleLengthFactors[0]) < 0)
					{
						printi(cycleLengthFactors[indexCounter]);
						printf(" (%d), ", cycleLengthCounts[indexCounter-1]);
						
						if (outputFile != NULL)
						{
							fprinti(outputFile, cycleLengthFactors[indexCounter]);
							fprintf(outputFile, " (%d), ", cycleLengthCounts[indexCounter-1]);
						}
						
						indexCounter += 1;
						add_BigIntT(temp, one, temp2);
						copy_BigIntT(temp2, temp);
					}
					printi(bigOmega);
					printf(" (%d)\n", cycleLengthCounts[indexCounter-1]);
					
					charaPoly = chara_poly(currMat, bigMod);
					charaPolyFactors = old_factor_BigPolyT(charaPoly, bigMod);
					printf("Chara poly: ");
					printp(charaPoly);
					printf("\nFactored chara poly: ");
					old_printpf(charaPolyFactors);
					printf("\n\n");
					
					if (outputFile != NULL)
					{
						fprinti(outputFile, bigOmega);
						fprintf(outputFile, " (%d)\n", cycleLengthCounts[indexCounter-1]);
						
						fprintf(outputFile, "Chara poly: ");
						fprintpf(outputFile, charaPolyFactors);
						fprintf(outputFile, "\n\n");
						
						//Saving file
						if (fclose(outputFile) == EOF)
						{
							fprintf(stderr, "Unable to save data. Continuing without saving...\n");
							outputFile = NULL;
						}
						
						else
						{
							outputFile = fopen(outputfilename, "a");
							if (outputFile == NULL)
								fprintf(stderr, "Unable to save data. Continuing without saving...\n");
						}
					}
					
					charaPoly = free_BigPolyT(charaPoly);
					charaPolyFactors = free_BigPolyT_factors(charaPolyFactors);
				}
				
				bigOmega = free_BigIntT(bigOmega);
				bigVectOmega = free_BigIntT(bigVectOmega);
				FREE(cycleLengthCounts);
				
				//Freeing this array is a massive pain
				//I have to do it each loop since the factorisation matrix allocates new
				// arrays each time.
				//I should probably change this in the future for speed.
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
				
				//Iterate to next matrix
				copy_BigIntT(currMatElements[matSize-1][matSize-2], prevLastElement);
				checkedAllMatrices = increment_BigIntT_array(currMatElements, matSize, matSize, one, bigMod);
			}
			
			if (debug)
				printf("Number of matrices with unjustified CCMs: %d\n", matrixCounter);
			
			
			currMatElements  = free_BigIntT_array(currMatElements, matSize, matSize);
			currVectElements = free_BigIntT_array(currVectElements, matSize, 1);
			
			currMat   = free_BigIntMatrixT(currMat);
			zeroMat   = free_BigIntMatrixT(zeroMat);
			currVect  = free_BigIntMatrixT(currVect);
			resumeMat = free_BigIntMatrixT(resumeMat);
			
			theCycle = free_CycleInfoT(theCycle);
			
			bigMod = free_BigIntT(bigMod);
			zero   = free_BigIntT(zero);
			one    = free_BigIntT(one);
			temp   = free_BigIntT(temp);
			temp2  = free_BigIntT(temp2);
			
			prevLastElement = free_BigIntT(prevLastElement);
			
			FREE(outputfilename);
			if (fclose(outputFile) == EOF)
				fprintf(stderr, "Unable to save data.\n");
		}
		
		
		//If the user wants to calculate every vector's 
		// "minimal" annihilating polynomial.
		else if (!strcmp(argv[1], "vectprops"))
		{
			bool checkedAllVects = FALSE;
			BigIntTP bigMod;
			BigIntTP temp;
			
			int numArr[1] = {1};
			BigIntTP one = new_BigIntT(numArr, 1);
			
			BigIntTP** currVectElements;
			BigIntTP*  vectPropElements; //The coeffs for the min poly we're testing on currVect
			
			BigIntMatrixTP bigF;   //Holds our update matrix
			BigIntMatrixTP exprF;  //Holds bigF evaluated in some polynomial
			BigIntMatrixTP currVect;
			BigIntMatrixTP tempVect;
			BigIntMatrixTP zeroVect;
			
			BigIntMatrixTP resumeVect = NULL;
			
			BigPolyTP annihPoly;
			
			char* fileName   = NULL;
			FILE* outputFile = NULL;
			
			if (argc > 2)
			{
				SET_BIG_NUM(argv[2], bigMod, "Unable to read modulus from command line.");
			}
			else
			{
				SET_BIG_NUM(bigintmodstring, bigMod, "Unable to read modulus from config file.");
			}
			
			bigF = read_BigIntMatrixT(updatefilepath);
			if (bigF == NULL)
			{
				fprintf(stderr, "Unable to read matrix from config file.\n");
				bigMod = free_BigIntT(bigMod);
				one    = free_BigIntT(one);
				FREE_VARIABLES;
				return EXIT_FAILURE;
			}
			
			if (big_rows(bigF) != big_cols(bigF))
			{
				fprintf(stderr, "Given update matrix isn't square.\n");
				bigMod = free_BigIntT(bigMod);
				one    = free_BigIntT(one);
				bigF   = free_BigIntMatrixT(bigF);
				FREE_VARIABLES;
				return EXIT_SUCCESS;
			}
			
			
			currVectElements = new_BigIntT_array(big_rows(bigF), 1);
			vectPropElements = malloc((big_rows(bigF)+1)*sizeof(BigIntTP));
			for (int i = 0; i < big_rows(bigF)+1; i += 1)
				vectPropElements[i] = empty_BigIntT(1);
			
			
			//If the user wants to resume computation from a
			// particular vector
			if (argc > 3)
			{
				if (!strcmp(argv[3], "TRUE"))
				{
					resumeVect = read_BigIntMatrixT(resumefilepath);
					if (resumeVect == NULL)
					{
						fprintf(stderr, "Unable to read resume vector from .config file. ");
						fprintf(stderr, "Beginning computation from the zero vector instead.\n");
					}
					
					else if ((big_rows(resumeVect) != big_rows(bigF)) || 
					         (big_cols(resumeVect) != 1))
					{
						fprintf(stderr, "Resume vector is of incorrect dimensions (must be %d by 1). ", big_rows(bigF));
						fprintf(stderr, "Beginning computation from the zero vector instead.\n");
					}
					
					else //Set currVect to resumeVect
						for (int elem = 0; elem < big_rows(bigF); elem += 1)
							copy_BigIntT(big_element(resumeVect, elem, 0), currVectElements[elem][0]);
						
					resumeVect = free_BigIntMatrixT(resumeVect);
				}
			}
	
			annihPoly = new_BigPolyT(vectPropElements, big_rows(bigF)+1);
			currVect = new_BigIntMatrixT(big_rows(bigF), 1);
			tempVect = new_BigIntMatrixT(big_rows(bigF), 1);
			zeroVect = new_BigIntMatrixT(big_rows(bigF), 1);
			exprF = new_BigIntMatrixT(big_rows(bigF), big_cols(bigF));
			temp = empty_BigIntT(1);
			
			//If user wants file output
			if (argc > 4)
			{
				if (!strcmp(argv[4], "TRUE"))
				{
					fileName = malloc(MAXSTRLEN*sizeof(char));
					fileName[0] = '\0';
					strcat(fileName, "vectprops F");
					for (int x = 0; x < big_rows(bigF); x += 1)
						for (int y = 0; y < big_cols(bigF); y += 1)
							append_BigIntT(fileName, big_element(bigF, x, y));
					strcat(fileName, " ");
					append_BigIntT(fileName, bigMod);
					strcat(fileName, ".txt");
					
					outputFile = fopen(fileName, "a");
					if (outputFile == NULL)
						fprintf(stderr, "Unable to open file for saving output. Continuing without saving...\n");
				}
			}
			
			//Iterate through all vectors in the module
			while (!checkedAllVects)
			{
				set_big_matrix(currVect, currVectElements);
				
				//Clear all elements in our polynomial coeffs list
				for (int i = 0; i < big_rows(bigF)+1; i += 1)
					clear_BigIntT(vectPropElements[i]);
				
				//Iterate through all possible "minimal" annihilating polynomials for currVect
				while (TRUE)
				{
					//Increment polynomial
					for (int i = 0; i < big_rows(bigF)+1; i += 1)
					{
						add_BigIntT(vectPropElements[i], one, temp);
						if (compare_BigIntT(temp, bigMod) >= 0)
						{
							clear_BigIntT(temp);
							copy_BigIntT(temp, vectPropElements[i]);
						}
						else
						{
							copy_BigIntT(temp, vectPropElements[i]);
							break;
						}
					}
					
					//Set polynomial, plug bigF into the polynomial
					set_BigPolyT(annihPoly, vectPropElements);
					eval_BigPolyT(annihPoly, bigF, exprF, bigMod);
					
					big_mat_mul(exprF, currVect, tempVect);
					modbm(tempVect, bigMod);
					
					//If we've found currVect's "minimal" annihilating polynomial
					if (compare_BigIntMatrixT(tempVect, zeroVect))
					{
						printbm_row(currVect);
						printf("'s \"min\" poly: ");
						printp(annihPoly);
						printf("\n");
						
						if (outputFile != NULL)
						{
							fprintbm_row(outputFile, currVect);
							fprintf(outputFile, " : ");
							fprintp(outputFile, annihPoly);
							fprintf(outputFile, "\n");
							
							if (fclose(outputFile) == EOF)
								fprintf(stderr, "Unable to save output file. Continuing without saving...\n");
							else
							{
								outputFile = fopen(fileName, "a");
								if (outputFile == NULL)
									fprintf(stderr, "Unable to save output file. Continuing without saving...\n");
							}
						}
						break;
					}
				}
				
				//Increment currVect
				checkedAllVects = increment_BigIntT_array(currVectElements, big_rows(bigF), 1, one, bigMod);
			}
			
			if (outputFile != NULL)
				if (fclose(outputFile) == EOF)
					fprintf(stderr, "Unable to save output file.\n");
				
			if (fileName != NULL)
			{
				FREE(fileName);
			}
			
			currVectElements = free_BigIntT_array(currVectElements, big_rows(bigF), 1);
			
			for (int i = 0; i < big_rows(bigF)+1; i += 1)
				vectPropElements[i] = free_BigIntT(vectPropElements[i]);
			FREE(vectPropElements);
			
			one    = free_BigIntT(one);
			temp   = free_BigIntT(temp);
			bigMod = free_BigIntT(bigMod);
			
			bigF     = free_BigIntMatrixT(bigF);
			exprF    = free_BigIntMatrixT(exprF);
			currVect = free_BigIntMatrixT(currVect);
			tempVect = free_BigIntMatrixT(tempVect);
			zeroVect = free_BigIntMatrixT(zeroVect);
			
			annihPoly = free_BigPolyT(annihPoly);
		}
		
		
		//If the user wants to calculate every annihilating vector under
		// an update matrix up to a certain degree
		else if (!strcmp(argv[1], "vectpolys"))
		{
			bool foundFirstMonic = FALSE;
			bool moreIdealToCheck; //More annihilating polynomials to check in the ideal?
			
			int deg;
			int numArr[1] = {1};

			BigIntTP bigMod, one, negOne, temp, temp2;
			
			BigIntTP** tempPolyCoeffs = NULL;
			BigPolyTP  tempAnnihPoly = NULL;
			BigPolyTP  OGtempAnnihPoly = NULL;
			
			BigPolyTP  tempPoly  = NULL;
			BigPolyTP  tempPoly2 = NULL;
			BigPolyTP  tempPoly3 = NULL;
			BigPolyTP* idealCoords = NULL; //Used to traverse all possible polynomials in our ideal
			BigPolyTP  zeroPoly;
			BigPolyTP  onePoly = NULL;
			
			int rewriteDegree    = 0; //The lowest power of a polynomial variable we can rewrite
			int numOfVectProps   = 0;
			BigPolyTP* vectProps = NULL; //Keeps track of formulas for higher powers of the polynomial variable
			BigPolyTP  rewriteConstant = NULL; //For getting multiples from vectProps
			
			BigIntTP* tempCoeffs;
			BigIntTP* tempCoeffs2;
			BigIntTP* tempCoeffs3;
			
			int numOfAnnihPolys = 0;
			BigPolyTP* annihPolys = NULL;
			
			BigIntMatrixTP A, v, zeroMat, tempMat, tempVect, zeroVect;
			
			deg = (int)strtol(argv[2], &tempStr, 10);
			if (tempStr[0] != '\0')
			{
				fprintf(stderr, "Unable to read degree from command line.\n");
				FREE_VARIABLES;
				return EXIT_FAILURE;
			}
			
			if (argc > 3)
			{
				SET_BIG_NUM(argv[3], bigMod, "Unable to read modulus from command line.");
			}
			else
			{
				SET_BIG_NUM(bigintmodstring, bigMod, "Unable to read modulus from config file.");
			}
			
			//if (argc > 4)
				
			A = UPDATEMATRIX;
			v = INITIALMATRIX;
			
			if (big_rows(A) != big_cols(A))
			{
				fprintf(stderr, "Given update matrix is not square.\n");
				bigMod = free_BigIntT(bigMod);
				FREE_VARIABLES;
				return EXIT_SUCCESS;
			}
			
			if ((big_rows(v) != big_rows(A)) || (big_cols(v) != 1))
			{
				fprintf(stderr, "Given vector's dimensions are invalid. The number of rows must be equal to the matrix's rows \
and the number of columns must be 1.\n");
				bigMod = free_BigIntT(bigMod);
				FREE_VARIABLES;
				return EXIT_SUCCESS;
			}
			
			one    = new_BigIntT(numArr, 1);
			temp   = empty_BigIntT(1);
			temp2  = empty_BigIntT(1);
			negOne = empty_BigIntT(1);
			subtract_BigIntT(bigMod, one, negOne);
			
			tempPolyCoeffs   = new_BigIntT_array(1, deg+1);
			tempAnnihPoly    = new_BigPolyT(tempPolyCoeffs[0], deg+1);
			OGtempAnnihPoly  = new_BigPolyT(tempPolyCoeffs[0], deg+1);
			
			tempPoly  = empty_BigPolyT();
			tempPoly2 = empty_BigPolyT();
			tempPoly3 = empty_BigPolyT();
			
			zeroPoly = empty_BigPolyT();
			onePoly  = constant_BigPolyT(one);
			
			zeroMat  = new_BigIntMatrixT(big_rows(A), big_cols(A));
			tempMat  = new_BigIntMatrixT(big_rows(A), big_cols(A));
			tempVect = new_BigIntMatrixT(big_rows(A), 1);
			zeroVect = new_BigIntMatrixT(big_rows(A), 1);
			
			printf("Matrix:\n");
			printbm(A);
			printf("Vector: ");
			printbm_row(v);
			printf("\nModulus: ");
			printi(bigMod);
			printf("\n");
			
			//Loop over each possible polynomial, check if it's an annihilating one for v
			while (!increment_BigIntT_array(tempPolyCoeffs, 1, deg+1, one, bigMod))
			{
				//When we find the first monic, calculate all future needed expressions
				// representing higher polynomial variable powers
				if ((foundFirstMonic) && (numOfVectProps == 0))
				{
					tempCoeffs  = extract_coefficients(annihPolys[numOfAnnihPolys-1]);
					tempCoeffs2 = malloc(rewriteDegree*sizeof(BigIntTP));
					
					for (int coeff = 0; coeff < rewriteDegree; coeff += 1)
					{
						tempCoeffs2[coeff] = empty_BigIntT(1);
						multiply_BigIntT(tempCoeffs[coeff], negOne, temp);
						mod_BigIntT(temp, bigMod, tempCoeffs2[coeff]);
					}
					
					//Now, create the actual polynomial representing this first polynomial term power
					numOfVectProps += 1;
					vectProps = realloc(vectProps, numOfVectProps*sizeof(BigPolyTP));
					vectProps[numOfVectProps-1] = new_BigPolyT(tempCoeffs2, rewriteDegree);
					
					for (int i = 0; i < rewriteDegree; i += 1)
					{
						tempCoeffs[i]  = free_BigIntT(tempCoeffs[i]);
						tempCoeffs2[i] = free_BigIntT(tempCoeffs2[i]);
					}
					FREE(tempCoeffs);
					FREE(tempCoeffs2);
					
					//One polynomial for each power
					tempCoeffs = extract_coefficients(vectProps[0]);
					for (int p = rewriteDegree+1; p <= deg; p += 1)
					{
						tempCoeffs3 = extract_coefficients(vectProps[p-rewriteDegree-1]);
						tempCoeffs2 = malloc((deg+1)*sizeof(BigIntTP));
						for (int coeff = 0; coeff < rewriteDegree; coeff += 1)
						{
							tempCoeffs2[coeff] = empty_BigIntT(1);
							multiply_BigIntT(tempCoeffs[coeff], tempCoeffs3[rewriteDegree-1], temp);
							if (coeff != 0)
								add_BigIntT(temp, tempCoeffs3[coeff-1], tempCoeffs2[coeff]);
							else
								copy_BigIntT(temp, tempCoeffs2[coeff]);
							
							//Now, take the modulus of each term
							mod_BigIntT(tempCoeffs2[coeff], bigMod, temp);
							copy_BigIntT(temp, tempCoeffs2[coeff]);
						}
						
						//Create the actual polynomial
						numOfVectProps += 1;
						vectProps = realloc(vectProps, numOfVectProps*sizeof(BigPolyTP));
						vectProps[numOfVectProps-1] = new_BigPolyT(tempCoeffs2, degree(annihPolys[numOfAnnihPolys-1]));
						
						for (int i = 0; i < rewriteDegree; i += 1)
						{
							tempCoeffs3[i]  = free_BigIntT(tempCoeffs3[i]);
							tempCoeffs2[i] = free_BigIntT(tempCoeffs2[i]);
						}
						FREE(tempCoeffs3);
						FREE(tempCoeffs2);
					}
					
					for (int i = 0; i < rewriteDegree; i += 1)
						tempCoeffs[i] = free_BigIntT(tempCoeffs[i]);
					FREE(tempCoeffs);
					
					//Now, let's print out our relations to make sure they're reasonable
					/*
					for (int d = degree(annihPolys[numOfAnnihPolys-1]); d-degree(annihPolys[numOfAnnihPolys-1]) < numOfVectProps; d += 1)
					{
						printf("λ^%d = ", d);
						printp(vectProps[d - degree(annihPolys[numOfAnnihPolys-1])]);
						printf("\n");
					}
					*/
				}
				
				//Keeping track of tempAnnihPoly since it can be reduced below
				set_BigPolyT(OGtempAnnihPoly, tempPolyCoeffs[0]);
				
				//If we can reduce our current polynomial using the already-found
				// annihilating polynomial
				if (rewriteDegree != 0)
				{
					resize_BigPolyT(tempAnnihPoly, rewriteDegree);
					set_BigPolyT(tempAnnihPoly, tempPolyCoeffs[0]);
					
					//Rewrite higher-order terms
					for (int term = rewriteDegree; term <= deg; term += 1)
					{
						rewriteConstant = constant_BigPolyT(tempPolyCoeffs[0][term]);
						multiply_BigPolyT(rewriteConstant, vectProps[term-rewriteDegree], tempPoly);
						add_BigPolyT(tempPoly, tempAnnihPoly, tempPoly2);
						mod_BigPolyT(tempPoly2, bigMod, tempAnnihPoly);
					}
				}
				
				else
				{
					resize_BigPolyT(tempAnnihPoly, deg+1);
					set_BigPolyT(tempAnnihPoly, tempPolyCoeffs[0]);
				}
				
				eval_BigPolyT(tempAnnihPoly, A, tempMat, bigMod);
				big_mat_mul(tempMat, v, tempVect);
				modbm(tempVect, bigMod);
				
				if (compare_BigIntMatrixT(tempVect, zeroVect))
				{
					//If we found an annihilating polynomial, let's check to make sure it isn't 
					// in the ideal generated by our current annihilating polynomials
					//We should also check for whether new annihPolys can be reduced to older ones (via substitution like in matprops)
					for (int i = 0; i < numOfAnnihPolys; i += 1)
						clear_BigPolyT(idealCoords[i]);
					
					moreIdealToCheck = TRUE;
					while (moreIdealToCheck)
					{
						//Create an element from the ideal and
						// see if our "new" annihilating polynomial
						// is actually new
						clear_BigPolyT(tempPoly);
						clear_BigPolyT(tempPoly2);
						clear_BigPolyT(tempPoly3);
						for (int poly = 0; poly < numOfAnnihPolys; poly += 1)
						{
							multiply_BigPolyT(idealCoords[poly], annihPolys[poly], tempPoly2);
							add_BigPolyT(tempPoly2, tempPoly, tempPoly3);
							mod_BigPolyT(tempPoly3, bigMod, tempPoly);
						}
						
						//If tempAnnihPoly is in our ideal
						if (compare_BigPolyT(tempPoly, tempAnnihPoly) == 0)
							break;
						
						//Increment idealCoords
						moreIdealToCheck = FALSE;
						for (int coord = 0; coord < numOfAnnihPolys; coord += 1)
						{
							add_BigPolyT(idealCoords[coord], onePoly, tempPoly);
							if (compare_BigIntT(constant(tempPoly), bigMod) >= 0)
								clear_BigPolyT(idealCoords[coord]);
							else
							{
								copy_BigPolyT(tempPoly, idealCoords[coord]);
								moreIdealToCheck = TRUE;
								break;
							}
						}
					}
					
					//If we found a new annihilating polynomial
					if (!moreIdealToCheck)
					{
						if (! foundFirstMonic)
						{
							reduce_BigPolyT(tempAnnihPoly);
							if (compare_BigIntT(leading_term(tempAnnihPoly), one) == 0)
							{
								rewriteDegree = degree(tempAnnihPoly);
								foundFirstMonic = TRUE;
							}
						}
						
						numOfAnnihPolys += 1;
						annihPolys = realloc(annihPolys, numOfAnnihPolys*sizeof(BigPolyTP));
						annihPolys[numOfAnnihPolys-1] = empty_BigPolyT();
						copy_BigPolyT(OGtempAnnihPoly, annihPolys[numOfAnnihPolys-1]);
						
						idealCoords = realloc(idealCoords, numOfAnnihPolys*sizeof(BigPolyTP));
						idealCoords[numOfAnnihPolys-1] = empty_BigPolyT();
					}
				}
			}
			
			if (numOfAnnihPolys > 0)
			{
				printf("Unique annihilating polynomials:\n");
				for (int p = 0; p < numOfAnnihPolys; p += 1)
				{
					printp(annihPolys[p]);
					printf("\n");
				}
			}
			
			one    = free_BigIntT(one);
			temp   = free_BigIntT(temp);
			temp2  = free_BigIntT(temp2);
			negOne = free_BigIntT(negOne);
			bigMod = free_BigIntT(bigMod);
			
			onePoly          = free_BigPolyT(onePoly);
			zeroPoly         = free_BigPolyT(zeroPoly);
			tempPoly         = free_BigPolyT(tempPoly);
			tempPoly2        = free_BigPolyT(tempPoly2);
			tempPoly3        = free_BigPolyT(tempPoly3);
			tempAnnihPoly    = free_BigPolyT(tempAnnihPoly);
			OGtempAnnihPoly  = free_BigPolyT(OGtempAnnihPoly);
			tempPolyCoeffs   = free_BigIntT_array(tempPolyCoeffs, 1, deg+1);
			
			for (int i = 0; i < numOfVectProps; i += 1)
				vectProps[i] = free_BigPolyT(vectProps[i]);
			FREE(vectProps);
			
			for (int i = 0; i < numOfAnnihPolys; i += 1)
			{
				idealCoords[i] = free_BigPolyT(idealCoords[i]);
				annihPolys[i]  = free_BigPolyT(annihPolys[i]);
			}
			FREE(idealCoords);
			FREE(annihPolys);
			
			zeroMat  = free_BigIntMatrixT(zeroMat);
			tempMat  = free_BigIntMatrixT(tempMat);
			tempVect = free_BigIntMatrixT(tempVect);
			zeroVect = free_BigIntMatrixT(zeroVect);
		}
		
		
		//If the user wants to calculate some useful properties of a given matrix
		// (simplify higher powers)
		else if (! strcmp(argv[1], "matprops"))
		{
			int maxpower;
			
			BigIntTP bigMod;
			BigIntTP temp  = NULL;
			BigIntTP temp2 = NULL;
			BigIntTP temp3 = NULL;
			BigIntTP zero  = NULL;
			BigIntTP one   = NULL;
			int oneArr[1]  = {1};
			
			BigIntMatrixTP A;
			BigPolyTP charaPoly;
			BigIntTP* OGcoeffs = NULL; //Holds the coefficients of charaPoly
			BigIntTP* newcoeffs = NULL; //Holds the coefficients for higher powers of A
			
			if (argc > 2)
			{
				maxpower = (int)strtol(argv[2], &tempStr, 10);
				
				if (tempStr[0] != '\0')
				{
					fprintf(stderr, "Unable to read maxpower from command line.\n");
					FREE_VARIABLES;
					return EXIT_FAILURE;
				}
			}
			
			else
			{
				printf(ANSI_COLOR_YELLOW "matprops " ANSI_COLOR_CYAN "maxpower [modulus]" ANSI_COLOR_RESET ": Calculates" \
				" some useful matrix properties.\n");
				printf(" - " ANSI_COLOR_CYAN "maxpower" ANSI_COLOR_RESET \
				": The highest power of the update matrix to calculate an equivalent expression for.\n");
				printf(" - " ANSI_COLOR_CYAN "[modulus]" ANSI_COLOR_RESET \
				": Overrides the modulus specified in the config file.\n");
				FREE_VARIABLES;
				return EXIT_SUCCESS;
			}
			
			if (argc > 3)
			{
				SET_BIG_NUM(argv[3], bigMod, "Invalid modulus passed on command line.");
			}
			
			else
			{
				SET_BIG_NUM(bigintmodstring, bigMod, "Invalid modulus in config file.");
			}
			
			A = read_BigIntMatrixT(updatefilepath);
			if (A == NULL)
			{
				fprintf(stderr, "Unable to read update matrix from config file.\n");
				FREE_VARIABLES;
				return EXIT_FAILURE;
			}
			
			charaPoly = chara_poly(A, bigMod);
			printp(charaPoly);
			printf("\n");
			
			//Check to make sure there's actually some computation to do
			if (degree(charaPoly) <= maxpower)
			{
				printf(":)\n");
				temp  = empty_BigIntT(1);
				temp2 = empty_BigIntT(1);
				temp3 = empty_BigIntT(1);
				zero  = empty_BigIntT(1);
				one   = new_BigIntT(oneArr, 1);
				OGcoeffs = extract_coefficients(charaPoly);
				
				//Now we need to divide by the leading coefficient, if it exists
				//If not, then don't do anything
				if (compare_BigIntT(OGcoeffs[degree(charaPoly)], one) != 0)
					big_num_inverse(OGcoeffs[degree(charaPoly)], bigMod, temp2);
				else
					copy_BigIntT(one, temp2);
				
				//Okay, we have the coefficients. Now to negate them.
				for (int i = 0; i < degree(charaPoly); i += 1)
				{
					subtract_BigIntT(bigMod, OGcoeffs[i], temp);
					mod_BigIntT(temp, bigMod, OGcoeffs[i]);
					multiply_BigIntT(temp2, OGcoeffs[i], temp);
					mod_BigIntT(temp, bigMod, OGcoeffs[i]);
				}
				
				//WHEN A NUMBER IS ZERO, ITS TERM SHOULDN'T BE PRINTED
				
				//This gives us the first equivalent expression
				printf("A^%d = ", degree(charaPoly));
				printi(OGcoeffs[0]);
				printf(" + ");
				for (int i = 1; i < degree(charaPoly)-1; i += 1)
				{
					printi(OGcoeffs[i]);
					printf("A^%d + ", i);
				}
				printi(OGcoeffs[degree(charaPoly)-1]);
				printf("A^%d\n", degree(charaPoly)-1);
				
				newcoeffs = malloc((degree(charaPoly)-1)*sizeof(BigIntTP));
				for (int i = 0; i < degree(charaPoly); i += 1)
				{
					newcoeffs[i] = empty_BigIntT(1);
					copy_BigIntT(OGcoeffs[i], newcoeffs[i]);
				}
				
				//Now, we compute the rest of the requested powers
				for (int i = degree(charaPoly)+1; i <= maxpower; i += 1)
				{
					//Shift all coefficients up a power
					copy_BigIntT(newcoeffs[degree(charaPoly)-1], temp);
					
					for (int n = degree(charaPoly)-2; n >= 0; n -= 1)
						copy_BigIntT(newcoeffs[n], newcoeffs[n+1]);
					copy_BigIntT(zero, newcoeffs[0]);
					
					//Now, simplify the big term into smaller exponent terms
					for (int n = 0; n < degree(charaPoly); n += 1)
					{
						multiply_BigIntT(temp, OGcoeffs[n], temp2);
						add_BigIntT(temp2, newcoeffs[n], temp3);
						mod_BigIntT(temp3, bigMod, newcoeffs[n]);
					}
					
					//Now, we print out the new expression
					printf("A^%d = ", i);
					printi(newcoeffs[0]);
					printf(" + ");
					for (int n = 1; n < degree(charaPoly)-1; n += 1)
					{
						printi(newcoeffs[n]);
						printf("A^%d + ", n);
					}
					printi(newcoeffs[degree(charaPoly)-1]);
					printf("A^%d\n", degree(charaPoly)-1);
				}
			}
			
			else
				printf("Degree of characteristic polynomial is less than maxpower. No computation to do.\n");
			
			for (int i = 0; i <= degree(charaPoly); i += 1)
			{
				if (i < degree(charaPoly))
					newcoeffs[i] = free_BigIntT(newcoeffs[i]);
				
				OGcoeffs[i] = free_BigIntT(OGcoeffs[i]);
			}
			FREE(OGcoeffs);
			FREE(newcoeffs);
			
			bigMod = free_BigIntT(bigMod);
			temp   = free_BigIntT(temp);
			temp2  = free_BigIntT(temp2);
			temp3  = free_BigIntT(temp3);
			zero   = free_BigIntT(zero);
			one    = free_BigIntT(one);
			
			A = free_BigIntMatrixT(A);
			
			charaPoly = free_BigPolyT(charaPoly);
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
			SET_BIG_NUM(argv[2], step, "Invalid step size passed on command line.");
			
			//If user provided a specific modulus on CLI
			if (argc > 3)
			{
				SET_BIG_NUM(argv[3], mod, "Invalid modulus passed on command line.");
			}
			
			else
			{
				SET_BIG_NUM(bigintmodstring, mod, "Unable to read modulus from config file.");
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
						factorList = old_factor_BigPolyT(charaPoly, mod);
						
						printbm(currentMatrix);
						printf("Cycle length: %d\n", omega(theCycle));
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
			
			
			currentMatrixElements = free_BigIntT_array(currentMatrixElements, 
			                                           big_rows(startingMatrix), 
																								 big_rows(startingMatrix));
			
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
		/*
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
		*/
		
		
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
			{
				SET_BIG_NUM(argv[2], bigModulus, "Unable to read modulus provided on command line.");
			}
			
			//Use config modulus instead
			else
			{
				SET_BIG_NUM(bigintmodstring, bigModulus, "Unable to read modulus from config file.");
			}
			
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
				SET_BIG_NUM(argv[2], upperbound, "Invalid upper bound passed.");
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
			int highpower = 2;
			int highmod;   //Holds the modulus raised to the correct power
			int maxmod;    //Holds the highest mod we obtain; used for iterating through matrix lifts

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
			
			IntMatrixTP A; //Holds a lift to our update matrix
			IntMatrixTP originalA; //Holds our original update matrix
			IntMatrixTP currVect;
			
			int** matrixLiftOffsets; //Used to iterate over all possible lifts of the given matrix A
			int matrixLiftIncrement; //The step size when creating new matrix lifts
			
			CycleInfoTP theCycle = NULL;
			
			char* outputfilename = NULL;
			FILE* outputFile = NULL;
			
			//These variables are for recording all possible dynamic configurations
			// [[[1, 1, 1, 10], [1, 1, 4, 20], ...], [[1, 1, 1, 5], [1, 1, 2, 15], ...], ...]
			int*** allConfigs = NULL;      //Holds cleaner representations of all found dynamics configurations
			int* allConfigsLengths = NULL; //Holds how many cycle tuples each config in allConfigs has
			int numOfConfigs = 0;
			
			//Holds an example matrix for each different dynamic behaviour found
			IntMatrixTP* allConfigsMatrices = NULL;
			
			bool findAllConfigs = FALSE;
			bool fileoutput = FALSE;
			bool isUniqueConfig = TRUE; //Used to determine whether we've found a new configuration once we've computed a new one
			
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
			
			//If user wants to find all possible dynamic configurations
			if (argc > 4)
				if (!strcmp(argv[4], "TRUE"))
					findAllConfigs = TRUE;
				
			//If user wants to store an output file
			if (argc > 5)
				if (!strcmp(argv[5], "TRUE"))
					fileoutput = TRUE;
			
			originalA = read_IntMatrixT(updatefilepath);
			if (originalA == NULL)
			{
				fprintf(stderr, "Unable to read .matrix file at %s.\n", updatefilepath);
				return EXIT_FAILURE;
			}
			
			else if (rows(originalA) != cols(originalA))
			{
				printf("Update matrix provided is not square.\n");
				return EXIT_SUCCESS;
			}
			
			printf("Matrix:\n");
			printm(originalA);
			printf("~~~~~~~\n");
			A = new_IntMatrixT(rows(originalA), rows(originalA));
			copy_IntMatrixT(originalA, A);
			
			//Creating filename
			if ((!findAllConfigs) && (fileoutput))
			{
				outputfilename = malloc(MAXSTRLEN*sizeof(char));
				outputfilename[0] = '\0';
				
				strcat(outputfilename, argv[1]); //dynamics
				strcat(outputfilename, " ");
				strcat(outputfilename, argv[2]); //maxPower
				strcat(outputfilename, " ");
				
				//Add modulus
				append_int(outputfilename, modulus);
				strcat(outputfilename, " ");
				
				//Now, include the update matrix elements
				strcat(outputfilename, "F");
				for (int x = 0; x < rows(A); x += 1)
					for (int y = 0; y < rows(A); y += 1)
						append_int(outputfilename, element(A, x, y));
					
				strcat(outputfilename, ".txt");
				printf("Vectors will be saved to %s\n", outputfilename);
				
				//Now, open the output file
				outputFile = fopen(outputfilename, "w");
				if (outputFile == NULL)
					fprintf(stderr, "Unable to create output file. Continuing without saving...\n");
			}
			
			//Initialise things for matrix lifts
			matrixLiftIncrement = modulus;
			
			matrixLiftOffsets = malloc(rows(A)*sizeof(int*));
			for (int i = 0; i < rows(A); i += 1)
				matrixLiftOffsets[i] = calloc(rows(A), sizeof(int));
			
			maxmod = modulus;
			for (int m = 1; m < highpower; m += 1)
				maxmod *= modulus;
			
			currVectElements = calloc(rows(A), sizeof(int));
			currVect         = new_IntMatrixT(rows(A), 1);
			
			//Dealing with exiting when highpower == 1 is done near the bottom of the tool
			while ((increment_int_array(matrixLiftOffsets, rows(A), rows(A), 1, maxmod/matrixLiftIncrement) != TRUE) ||
			       (highpower == 1))
			{
				highmod = 1;
				for (int modulusCounter = 1; modulusCounter <= highpower; modulusCounter += 1)
				{
					checkedAllVects = FALSE;
					
					highmod *= modulus;
					
					//Now that we have the matrix, calculate the cycle length
					// and use it to reason all possible cycle lengths
					// for alloting memory
					theCycle = floyd(A, A, highmod);
					
					if (!findAllConfigs)
					{
						printf("Current modulus: %d\n", highmod);
						printf("Matrix's multiplicative order mod %d: %d\n", highmod, omega(theCycle));
					}
					
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
					
					//Now, construct our cycleTable to hold cycle counts
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
					
					//Helps find which previous tuple this vector corresponds to
					//Contains the cycle lengths of our vector under all moduli
					// less than the current one. The last element is the
					// current modulus' cycle length
					currCycleLengths = malloc(modulusCounter*sizeof(int));
				
					//Iterate until we've tested every vector
					while (!checkedAllVects)
					{
						set_column(currVect, currVectElements);
						
						//Start with highest modulus, work our way down
						//Allows taking vectors mod tempModulus to be easier
						tempModulus = modulus;
						for (int i = 1; i < modulusCounter; tempModulus *= modulus, i += 1);
						
						for (int i = modulusCounter-1; i >= 0; i -= 1)
						{
							theCycle = floyd(A, currVect, tempModulus);
							currCycleLengths[i] = omega(theCycle);
							
							theCycle = free_CycleInfoT(theCycle);
							
							tempModulus /= modulus;
							modm(currVect, tempModulus);
						}
						
						/*
						printf("currTuple: (");
						for (int cc = 0; cc < modulusCounter; cc += 1)
							printf("%d,", currCycleLengths[cc]);
						printf(")\n");*/

						//Now, output to our file the current vector's group only if:
						// we're on the last modulus to check and if the file exists
						if ((modulusCounter == highpower) && (outputFile != NULL))
						{
							fprintf(outputFile, "(");
							for (int i = 0; i < modulusCounter-1; i += 1)
								fprintf(outputFile, "%d, ", currCycleLengths[i]);
							fprintf(outputFile, "%d) : ", currCycleLengths[modulusCounter-1]);
							
							//Now, output the vector
							set_column(currVect, currVectElements);							
							fprintm_row(outputFile, currVect);
							fprintf(outputFile, "\n");
							
							//Now, save the file
							if (fclose(outputFile) == EOF)
							{
								fprintf(stderr, "Unable to properly save data. Continuing without saving...\n");
								outputFile = NULL;
							}
							
							else
								outputFile = fopen(outputfilename, "a");
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
									/*
									printf("Matched prev tuple: (");
									for (int pc = 0; pc < modulusCounter-1; pc += 1)
										printf("%d,", prevCycleTuples[tuple][pc]);
									printf(")\n");*/
									
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
								/*
								printf("Cycle length found: %d\n", possibleCycleLengths[secondIndex]);*/
								cycleTable[firstIndex][secondIndex-1] += 1;
								break;
							}
						}

						//Increment through all possible vectors (% modulus)
						//When I switch this program to use BigIntMatrixTPs, I can use my
						// new function here
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
					
					FREE(currCycleLengths);
					
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
				
				//Now, we print all the numbers we've collected, or prepare to check for a unique tuple
				if (!findAllConfigs)
					printf("\nCycle tuples:\n");
				else
				{
					//Allocating space to hold possible new config and its length
					numOfConfigs += 1;
					allConfigsLengths = realloc(allConfigsLengths, numOfConfigs*sizeof(int));
					allConfigsLengths[numOfConfigs-1] = 0;
					
					allConfigs = realloc(allConfigs, numOfConfigs*sizeof(int**));
					allConfigs[numOfConfigs-1] = NULL;
				}
				
				//Loop through all our arrays and extract the tuple data
				for (int i = 0; i < prevCycleTuples[0][0]; i += 1)
				{
					for (int j = 0; j < possibleCycleLengths[0]; j += 1)
					{
						//Only print out info if vectors have that tuple of cycle lengths
						if (cycleTable[i][j] != 0)
						{
							if (findAllConfigs)
							{
								//Create proper space for storing current tuple
								allConfigsLengths[numOfConfigs-1] += 1;
								allConfigs[numOfConfigs-1] = realloc(allConfigs[numOfConfigs-1], 
																										 allConfigsLengths[numOfConfigs-1]*sizeof(int*));
								allConfigs[numOfConfigs-1][allConfigsLengths[numOfConfigs-1]-1] = malloc((highpower+1)*sizeof(int));
							}
							
							else
								printf("(");
							
							//Here, we either print the tuples or store them in allConfigs, checking whether they're unique
							//printf("(");
							for (int tupleElem = 0; tupleElem < highpower-1; tupleElem += 1)
							{
								if (!findAllConfigs)
									printf("%d, ", prevCycleTuples[i+1][tupleElem]);
								
								//Store tuple in our record of tuple configurations, even if it's not a new one
								//We'll worry about whether it's a new one later on
								// If it's not new, we'll just reallocate our memory back to how it was before
								else
									allConfigs[numOfConfigs-1][allConfigsLengths[numOfConfigs-1]-1][tupleElem] = prevCycleTuples[i+1][tupleElem];
								
								//printf("%d,", prevCycleTuples[i+1][tupleElem]);
							}

							//printf("%d) = %d\n", possibleCycleLengths[j+1], cycleTable[i][j]);
							
							if (!findAllConfigs)
								printf("%d) = %d\n", possibleCycleLengths[j+1], cycleTable[i][j]);
							
							else
							{
								allConfigs[numOfConfigs-1][allConfigsLengths[numOfConfigs-1]-1][highpower-1] = possibleCycleLengths[j+1];
								allConfigs[numOfConfigs-1][allConfigsLengths[numOfConfigs-1]-1][highpower] = cycleTable[i][j];
							}
						}
					}
				}
				
				//Now, we check to see if we've seen this specific tuple before
				//Essentially, we're trying to prove that it's not unique
				isUniqueConfig = TRUE;
				
				//Print out configs 
				/*
				for (int prevConfig = 0; prevConfig < numOfConfigs-1; prevConfig += 1)
				{
					for (int tuple = 0; tuple < allConfigsLengths[prevConfig]; tuple += 1)
					{
						printf("(");
						for (int entry = 0; entry < highpower+1; entry += 1)
							printf("%d,", allConfigs[prevConfig][tuple][entry]);
						printf(")\n");
					}
					printf("~~~\n");
				}
				*/
				
				//Determining whether configurations are equal or not hinges on the
				// assumptipon that configurations with equivalent tuples will appear
				// in exactly the same order.
				//This SHOULD always be the case due to how configs are copied into
				// allConfigs, but I'm leaving a note here just in case.
				for (int prevConfig = 0; prevConfig < numOfConfigs-1; prevConfig += 1)
				{					
					//If the two configs have the same number of tuples
					if (allConfigsLengths[prevConfig] == allConfigsLengths[numOfConfigs-1])
					{
						isUniqueConfig = FALSE;
						for (int currTuple = 0; currTuple < allConfigsLengths[numOfConfigs-1]; currTuple += 1)
						{
							for (int entry = 0; entry < highpower+1; entry += 1)
							{
								//If at least one number in the tuple is different
								if (allConfigs[prevConfig][currTuple][entry] != allConfigs[numOfConfigs-1][currTuple][entry])
								{
									isUniqueConfig = TRUE;
									break;
								}
							}
								
							if (isUniqueConfig)
								break;
						}
					}
					
					//If we've guaranteed that this configuration isn't unique
					if (!isUniqueConfig)
					{
						//printf("#%d is not unique!\n", numOfConfigs);
						//Deallocate new tuples since they're not actually unique and we don't need them
						//This is definitely inefficient, but I don't want a potential memory leak
						for (int n = 0; n < allConfigsLengths[numOfConfigs-1]; n += 1)
						{
							FREE(allConfigs[numOfConfigs-1][n]);
						}
						FREE(allConfigs[numOfConfigs-1]);
						allConfigsLengths[numOfConfigs-1] = 0;
						numOfConfigs -= 1;
						break;
					}
				}
				
				//If it is unique, add the matrix to the array 
				//The allConfigsMatrices == NULL prevents already stored matrices from being overriden when there's only one unique config
				if ((findAllConfigs) && ((isUniqueConfig) || ((numOfConfigs == 1) && (allConfigsMatrices == NULL))))
				{
					allConfigsMatrices = realloc(allConfigsMatrices, numOfConfigs*sizeof(IntMatrixTP));
					allConfigsMatrices[numOfConfigs-1] = new_IntMatrixT(rows(A), rows(A));
					copy_IntMatrixT(A, allConfigsMatrices[numOfConfigs-1]);
				}
				
				if ((outputFile != NULL) && (fclose(outputFile) == EOF))
					fprintf(stderr, "Unable to save data.\n");
				
				//Don't need to free theCycle since it gets freed above
				theCycle = NULL;
				
				for (int i = 0; i < prevCycleTuples[0][0]; i += 1)
				{
					FREE(cycleTable[i]);
				}
				FREE(cycleTable);
				FREE(possibleCycleLengths);
				
				for (int i = 1; i < prevCycleTuples[0][0]; i += 1)
				{
					FREE(prevCycleTuples[i]);
				}
				FREE(prevCycleTuples[0]);
				FREE(prevCycleTuples);
				
				if (outputfilename != NULL)
				{
					FREE(outputfilename);
				}
				
				//Ensuring we don't have to loop a bunch of times when only computing one dynamics configuration
				if ((!findAllConfigs) || (highpower == 1))
					break;
				
				//Set A to the next lift
				else
				{
					for (int r = 0; r < rows(originalA); r += 1)
						for (int c = 0; c < rows(originalA); c += 1)
						{
							matrixLiftOffsets[r][c] *= matrixLiftIncrement;
							matrixLiftOffsets[r][c] += element(originalA, r, c);
						}
						
					set_matrix(A, matrixLiftOffsets);
					modm(A, maxmod);
					
					for (int r = 0; r < rows(originalA); r += 1)
						for (int c = 0; c < rows(originalA); c += 1)
						{
							matrixLiftOffsets[r][c] -= element(originalA, r, c);
							matrixLiftOffsets[r][c] /= matrixLiftIncrement;
						}
				}
			}
			
			//Now, let's try and print out our record of unique dynamics configurations
			for (int config = 0; config < numOfConfigs; config += 1)
			{
				printf("Config #%d:\n", config);
				for (int tuple = 0; tuple < allConfigsLengths[config]; tuple += 1)
				{
					printf("(");
					for (int num = 0; num < highpower-1; num += 1)
						printf("%d, ", allConfigs[config][tuple][num]);
					printf("%d) -> %d\n", allConfigs[config][tuple][highpower-1], allConfigs[config][tuple][highpower]);
				}
				printf("Representative matrix:\n");
				printm(allConfigsMatrices[config]);
				printf("\n");
			}
			
			if (findAllConfigs)
			{
				for (int config = 0; config < numOfConfigs; config += 1)
				{
					allConfigsMatrices[config] = free_IntMatrixT(allConfigsMatrices[config]);
					for (int tuple = 0; tuple < allConfigsLengths[config]; tuple += 1)
					{
						FREE(allConfigs[config][tuple]);
					}
					FREE(allConfigs[config]);
				}
				FREE(allConfigs);
				
				FREE(allConfigsLengths);
				
				FREE(allConfigsMatrices);
			}
			
			for (int i = 0; i < rows(originalA); i += 1)
			{
				FREE(matrixLiftOffsets[i]);
			}
			FREE(matrixLiftOffsets);
			
			FREE(currVectElements);
			currVect = free_IntMatrixT(currVect);
				
			originalA = free_IntMatrixT(originalA);
			A = free_IntMatrixT(A);
		}
		
		
		//If the user wants to see how orbits from higher moduli map to lower moduli orbits
		else if (! strcmp(argv[1], "orbitmaps"))
		{
			BigIntTP bigMod;
			BigIntTP bigModPower; //Holds bigMod^maxpower
			BigIntTP one;
			BigIntTP temp;
			BigIntTP temp2;
			
			int oneArr[1] = {1};
			
			BigIntMatrixTP A;
			BigIntMatrixTP currVect;
			BigIntMatrixTP tempVect;
			BigIntMatrixTP tempVect2;
			BigIntTP** currVectElements;
			
			CycleInfoTP theCycle = NULL;
			
			BigIntMatrixTP** orbitReps;
			int** orbitLengths;
			int** orbitMaps; //Keeps track of how the orbits map onto each other
			int*  numOfOrbits;
			
			int   groupCount = 0;
			int   currCycleLength;
			int   prevCycleLength; //Used to sort mappings by cycle length
			int*  currOrbitMapMatch = NULL; //For organising output intogroupings of reductions
			int** currOrbitMapGroup = NULL; //Holds the vectors that match our match
			
			int maxpower;     //Holds the max power as a regular int
			int reducedOrbit; //Holds which orbit our current orbit feeds into
			
			bool isNewOrbit; //For checking whether particular orbits are already in our list
			bool matchingMap; //For finding maps that match our grouping for output
			bool printedAllMaps = FALSE;
			bool newCycleLength; //For determining whether we found a new cycle length when printing results
			
			FILE* outputFile = NULL;
			char* outputFileName = NULL;
			bool fileoutput = FALSE;
			
			FILE* graphFile = NULL; //GRAPHVIS output file
			char* graphFileName = NULL;
			int tempCounter = 0;
			
			if (argc > 2)
			{
				maxpower = (int)strtol(argv[2], &tempStr, 10);
				if (tempStr[0] != '\0')
				{
					fprintf(stderr, "Unable to read maxpower from command line.\n");
					FREE_VARIABLES;
					return EXIT_FAILURE;
				}
			}
			else
			{
				printf(ANSI_COLOR_YELLOW "orbitmaps " ANSI_COLOR_CYAN "maxpower [modulus] [fileoutput]" ANSI_COLOR_RESET \
				": Calculates orbit representatives for a higher-powered modulus and sees how they map down to lower-powered moduli.\n");
				printf(" - " ANSI_COLOR_CYAN "maxpower" ANSI_COLOR_RESET \
				": The highest (initial) power to use for the modulus.\n");
				printf(" - " ANSI_COLOR_CYAN "[modulus]" ANSI_COLOR_RESET \
				": Says what (base) modulus to use.\n\n");
				printf(" - " ANSI_COLOR_CYAN "[fileoutput]" ANSI_COLOR_RESET \
				": Says whether to create an output file.\n\n");
				
				FREE_VARIABLES;
				return EXIT_SUCCESS;
			}
			
			if (argc > 3)
			{
				SET_BIG_NUM(argv[3], bigMod, "Unable to read modulus from command line.");
			}
			else
			{
				SET_BIG_NUM(bigintmodstring, bigMod, "Unable to read modulus from config file.");
			}
			
			if (argc > 4)
				if (!strcmp(argv[4], "TRUE"))
					fileoutput = TRUE;
			
			A = read_BigIntMatrixT(updatefilepath);
			if (A == NULL)
			{
				fprintf(stderr, "Unable to read update matrix provided in config file.\n");
				FREE_VARIABLES;
				return EXIT_FAILURE;
			}
			
			//Create filename
			if (fileoutput)
			{
				outputFileName = malloc(MAXSTRLEN*sizeof(char));
				outputFileName[0] = '\0';
				graphFileName = malloc(MAXSTRLEN*sizeof(char));
				graphFileName[0] = '\0';
				
				strcat(outputFileName, "orbitmaps ");
				strcat(graphFileName, "graph ");
				append_int(outputFileName, maxpower);
				append_int(graphFileName, maxpower);
				strcat(outputFileName, " ");
				strcat(graphFileName, " ");
				append_BigIntT(outputFileName, bigMod);
				append_BigIntT(graphFileName, bigMod);
				strcat(outputFileName, " F");
				strcat(graphFileName, " F");
				for (int x = 0; x < big_rows(A); x += 1)
					for (int y = 0; y < big_cols(A); y += 1)
					{
						append_BigIntT(outputFileName, big_element(A, x, y));
						append_BigIntT(graphFileName, big_element(A, x, y));
					}
					
				strcat(outputFileName, ".txt");
				strcat(graphFileName, ".graph");
				
				outputFile = fopen(outputFileName, "w");
				if (outputFile == NULL)
					fprintf(stderr, "Unable to create output file. Continuing without saving...\n");
			}
			
			one   = new_BigIntT(oneArr, 1);
			temp  = empty_BigIntT(1);
			temp2 = empty_BigIntT(1);
			
			bigModPower = new_BigIntT(oneArr, 1);
			
			orbitReps    = calloc(maxpower, sizeof(BigIntMatrixTP*));
			orbitLengths = calloc(maxpower, sizeof(int*));
			numOfOrbits  = calloc(maxpower, sizeof(int));
			
			currVectElements = new_BigIntT_array(big_rows(A), 1);
			currVect  = new_BigIntMatrixT(big_rows(A), 1);
			tempVect  = new_BigIntMatrixT(big_rows(A), 1);
			tempVect2 = new_BigIntMatrixT(big_rows(A), 1);
			
			//Calculate bigMod^maxpower
			for (int i = 0; i < maxpower; i += 1)
			{
				multiply_BigIntT(bigMod, bigModPower, temp);
				copy_BigIntT(temp, bigModPower);
			}
			
			printf("Base modulus: ");
			printi(bigMod);
			printf("\nMax power: %d\n", maxpower);
			printf("Matrix:\n");
			printbm(A);
			printf("~~~~~~~~~\n");
			
			//Iterate over each vector in our highest modulus to find all orbits
			do
			{
				set_big_matrix(currVect, currVectElements);
				big_floyd(A, currVect, bigModPower, &theCycle);
				
				if (tau(theCycle) != 0)
					continue;
				
				//If we don't have any orbits yet
				if (numOfOrbits[maxpower-1] == 0)
				{
					numOfOrbits[maxpower-1] += 1;
					orbitReps[maxpower-1] = realloc(orbitReps[maxpower-1], numOfOrbits[maxpower-1]*sizeof(BigIntMatrixTP));
					orbitLengths[maxpower-1] = realloc(orbitLengths[maxpower-1], numOfOrbits[maxpower-1]*sizeof(int));
					
					//Saving orbit data
					orbitReps[maxpower-1][numOfOrbits[maxpower-1]-1] = new_BigIntMatrixT(big_rows(A), 1);
					copy_BigIntMatrixT(rep(theCycle), orbitReps[maxpower-1][numOfOrbits[maxpower-1]-1]); //Does it implicitly cast?
					orbitLengths[maxpower-1][numOfOrbits[maxpower-1]-1] = omega(theCycle);
				}
				
				//If we already have orbits, we have to check and see if 
				// we've already accounted for our current orbit
				else
				{
					isNewOrbit = TRUE;
					
					//Iterate over all the orbits we already have
					for (int i = 0; i < numOfOrbits[maxpower-1]; i += 1)
					{
						copy_BigIntMatrixT(orbitReps[maxpower-1][i], tempVect);
						
						//Check to see if any vector in the cycle matches currVect
						do
						{
							if (compare_BigIntMatrixT(tempVect, currVect))
							{
								isNewOrbit = FALSE;
								break;
							}
							
							big_mat_mul(A, tempVect, tempVect2);
							modbm(tempVect2, bigModPower);
							copy_BigIntMatrixT(tempVect2, tempVect);
						}
						while (!compare_BigIntMatrixT(tempVect, orbitReps[maxpower-1][i])); //Will this implicitly cast?
						
						if (!isNewOrbit)
							break;
					}
					
					if (isNewOrbit)
					{
						numOfOrbits[maxpower-1] += 1;
						orbitReps[maxpower-1] = realloc(orbitReps[maxpower-1], numOfOrbits[maxpower-1]*sizeof(BigIntMatrixTP));
						orbitLengths[maxpower-1] = realloc(orbitLengths[maxpower-1], numOfOrbits[maxpower-1]*sizeof(int));
						
						orbitLengths[maxpower-1][numOfOrbits[maxpower-1]-1] = omega(theCycle);
						orbitReps[maxpower-1][numOfOrbits[maxpower-1]-1] = new_BigIntMatrixT(big_rows(A), 1);
						copy_BigIntMatrixT(rep(theCycle), orbitReps[maxpower-1][numOfOrbits[maxpower-1]-1]);
					}
				}
			}
			while (!increment_BigIntT_array(currVectElements, big_rows(A), 1, one, bigModPower));
			
			orbitMaps = malloc(numOfOrbits[maxpower-1]*sizeof(int*));
			for (int i = 0; i < numOfOrbits[maxpower-1]; i += 1)
			{
				orbitMaps[i] = calloc(maxpower, sizeof(int));
				orbitMaps[i][maxpower-1] = i;
			}
			
			//Now, we have to find the orbitreps for all the lower moduli
			for (int i = maxpower-2; i >= 0; i -= 1)
			{
				//Calculate new modulus
				divide_BigIntT(bigModPower, bigMod, temp);
				copy_BigIntT(temp, bigModPower);
				
				//Go through each rep from the last modulus we calculated
				for (int cyc = 0; cyc < numOfOrbits[i+1]; cyc += 1)
				{
					copy_BigIntMatrixT(orbitReps[i+1][cyc], currVect);
					modbm(currVect, bigModPower);
					
					if (numOfOrbits[i] == 0)
					{
						reducedOrbit = 0;
						numOfOrbits[i] += 1;
						orbitReps[i] = realloc(orbitReps[i], numOfOrbits[i]*sizeof(BigIntMatrixTP));
						orbitLengths[i] = realloc(orbitLengths[i], numOfOrbits[i]*sizeof(int));
						
						big_floyd(A, currVect, bigModPower, &theCycle);
						orbitLengths[i][numOfOrbits[i]-1] = omega(theCycle);
						orbitReps[i][numOfOrbits[i]-1] = new_BigIntMatrixT(big_rows(A), 1);
						copy_BigIntMatrixT(rep(theCycle), orbitReps[i][numOfOrbits[i]-1]);
					}
					
					//Check to see if current higher mod orbit maps to a 
					// lower mod orbit we already have
					else
					{
						isNewOrbit = TRUE;
						
						//Iterate over all our current reps
						for (int newOrbs = 0; newOrbs < numOfOrbits[i]; newOrbs += 1)
						{
							copy_BigIntMatrixT(orbitReps[i][newOrbs], tempVect);
							do
							{
								//If our rep is already accounted for with a cycle
								if (compare_BigIntMatrixT(tempVect, currVect))
								{
									isNewOrbit = FALSE;
									reducedOrbit = newOrbs;
									break;
								}
								
								big_mat_mul(A, tempVect, tempVect2);
								modbm(tempVect2, bigModPower);
								copy_BigIntMatrixT(tempVect2, tempVect);
							}
							while (!compare_BigIntMatrixT(tempVect, orbitReps[i][newOrbs]));
							
							if (!isNewOrbit)
								break;
						}
						
						//Add new orbit to list if necessary
						if (isNewOrbit)
						{
							reducedOrbit = numOfOrbits[i];
							numOfOrbits[i] += 1;
							orbitReps[i] = realloc(orbitReps[i], numOfOrbits[i]*sizeof(BigIntMatrixTP));
							orbitLengths[i] = realloc(orbitLengths[i], numOfOrbits[i]*sizeof(int));
							
							big_floyd(A, currVect, bigModPower, &theCycle);
							orbitLengths[i][numOfOrbits[i]-1] = omega(theCycle);
							orbitReps[i][numOfOrbits[i]-1] = new_BigIntMatrixT(big_rows(A), 1);
							copy_BigIntMatrixT(rep(theCycle), orbitReps[i][numOfOrbits[i]-1]);
						}
					}
					
					//Make sure to keep track of mappings
					for (int map = 0; map < numOfOrbits[maxpower-1]; map += 1)
						if (orbitMaps[map][i+1] == cyc)
							orbitMaps[map][i] = reducedOrbit;
				}
			}
			
			//Now, let's check to see if I got the orbitreps correctly
			for (int e = maxpower-1; e >= 0; e -= 1)
			{
				printf("For mod ");
				printi(bigMod);
				printf("^%d:\n", e+1);
				
				for (int i = 0; i < numOfOrbits[e]; i += 1)
				{
					printf("Rep: ");
					if (vectorType == row)
					{
						printbm_row(orbitReps[e][i]);
						printf(",\tOrbit length: %d\n", orbitLengths[e][i]);
					}
					else if (vectorType == col)
					{
						printf("\n");
						printbm(orbitReps[e][i]);
						printf("Orbit length: %d\n", orbitLengths[e][i]);
					}
					
				}
				
				printf("\n");
			}
			
			currOrbitMapMatch = calloc((maxpower-1), sizeof(int));
			
			//Go through the mappings, organise them based on common mappings
			while (!printedAllMaps)
			{
				groupCount = 0;
				
				for (int map = 0; map < numOfOrbits[maxpower-1]; map += 1)
				{
					matchingMap = TRUE;
					
					//Checking each individual reduction to see if they match
					for (int modLevel = 0; modLevel < maxpower-1; modLevel += 1)
					{
						//If the map we checked isn't equal somewhere
						if (orbitMaps[map][modLevel] != currOrbitMapMatch[modLevel])
						{
							matchingMap = FALSE;
							break;
						}
					}
					
					//Add matching map to the group
					if (matchingMap)
					{
						groupCount += 1;
						currOrbitMapGroup = realloc(currOrbitMapGroup, groupCount*sizeof(int*));
						currOrbitMapGroup[groupCount-1] = orbitMaps[map];
					}
				}
				
				//Now, we should have every matching mapping, sorted by cycle length
				//Print them out, and save any relevant files
				prevCycleLength = 0;
				while (TRUE)
				{
					newCycleLength = FALSE;
					currCycleLength = -1;
					
					//Go through mappings and find the next lowest cycle length
					for (int map = 0; map < groupCount; map += 1)
					{
						if ((orbitLengths[maxpower-1][currOrbitMapGroup[map][maxpower-1]] > prevCycleLength) &&
						    ((currCycleLength == -1) ||
								(orbitLengths[maxpower-1][currOrbitMapGroup[map][maxpower-1]] < currCycleLength)))
						{
							currCycleLength = orbitLengths[maxpower-1][currOrbitMapGroup[map][maxpower-1]];
							newCycleLength = TRUE;
						}
					}
					
					if (newCycleLength)
					{
						for (int map = 0; map < groupCount; map += 1)
						{
							if (orbitLengths[maxpower-1][currOrbitMapGroup[map][maxpower-1]] == currCycleLength)
							{
								for (int modLevel = maxpower-1; modLevel > 0; modLevel -= 1)
								{
									if (vectorType == row)
									{
										printbm_row(orbitReps[modLevel][currOrbitMapGroup[map][modLevel]]);
										printf(" (ω = %d) -> ", orbitLengths[modLevel][currOrbitMapGroup[map][modLevel]]);
									}
									else if (vectorType == col)
									{
										printbm(orbitReps[modLevel][currOrbitMapGroup[map][modLevel]]);
										printf("↓ (ω = %d)\n", orbitLengths[modLevel][currOrbitMapGroup[map][modLevel]]);
									}
									if (outputFile != NULL)
									{
										fprintbm_row(outputFile, orbitReps[modLevel][currOrbitMapGroup[map][modLevel]]);
										fprintf(outputFile, " (ω = %d) -> ", orbitLengths[modLevel][currOrbitMapGroup[map][modLevel]]);
									}
								}
								if (vectorType == row)
								{
									printbm_row(orbitReps[0][currOrbitMapGroup[map][0]]);
									printf(" (ω = %d)\n", orbitLengths[0][currOrbitMapGroup[map][0]]);
								}
								else if (vectorType == col)
								{
									printbm(orbitReps[0][currOrbitMapGroup[map][0]]);
									printf("(ω = %d)\n\n", orbitLengths[0][currOrbitMapGroup[map][0]]);
								}
								if (outputFile != NULL)
								{
									fprintbm_row(outputFile, orbitReps[0][currOrbitMapGroup[map][0]]);
									fprintf(outputFile, " (ω = %d)\n", orbitLengths[0][currOrbitMapGroup[map][0]]);
								}
							}
						}
						
						prevCycleLength = currCycleLength;
					}
					
					else
						break;
				}
				
				//Preventing extra newlines from cluttering the screen
				if (groupCount > 0)
				{
					printf("\n");
					if (outputFile != NULL)
						fprintf(outputFile, "\n");
				}
				
				//Now, we change which mapping we're looking for and repeat
				for (int i = maxpower-2; i >= 0; i -= 1)
				{
					printedAllMaps = TRUE;
					
					currOrbitMapMatch[i] += 1;
					if (i == 0)
					{
						printf("==========\n\n");
						if (outputFile != NULL)
							fprintf(outputFile, "=========\n\n");
					}
					
					if (currOrbitMapMatch[i] == numOfOrbits[i])
						currOrbitMapMatch[i] = 0;
					else
					{
						printedAllMaps = FALSE;
						break;
					}
				}
			}
			
			if (fileoutput)
			{
				graphFile = fopen(graphFileName, "w");
				if (graphFile == NULL)
					fprintf(stderr, "Unable to save .graph file.\n");
			}
			
			if (graphFile != NULL)
			{
				fprintf(graphFile, "~b:False\n"); //Turning on the bounding box
				
				fprintf(graphFile, "~n:");
				for (int i = 0; i < maxpower; i += 1)
					tempCounter += numOfOrbits[i];
				
				fprintf(graphFile, "%d\n~l:above\n", tempCounter);
				tempCounter = 0;
				
				//Give labels to every vector we've saved
				for (int p = 0; p < maxpower; p += 1)
					for (int v = 0; v < numOfOrbits[p]; v += 1)
					{
						fprintf(graphFile, "%d:", tempCounter);
						fprintbm_row(graphFile, orbitReps[p][v]);
						fprintf(graphFile, "\n");
						tempCounter += 1;
					}
					
				fprintf(graphFile, "~c:direction\n");
				
				//Add our reductions to the .graph file
				for (int map = 0; map < numOfOrbits[maxpower-1]; map += 1)
				{
					for (int pow = maxpower-1; pow > 0; pow -= 1)
					{
						//ID of first vector
						tempCounter = 0;
						for (int prevs = 0; prevs < pow; tempCounter += numOfOrbits[prevs], prevs += 1);
						tempCounter += orbitMaps[map][pow];
						fprintf(graphFile, "%d,", tempCounter);
						
						//ID of second vector
						tempCounter = 0;
						for (int prevs = 0; prevs < pow-1;  tempCounter += numOfOrbits[prevs], prevs += 1);
						tempCounter += orbitMaps[map][pow-1];
						fprintf(graphFile, "%d\n", tempCounter);
					}
				}
			}
			
			if (outputFile != NULL)
				if (fclose(outputFile) == EOF)
					fprintf(stderr, "Unable to close output file. Unable to save data.\n");
				
			if (graphFile != NULL)
				if (fclose(graphFile) == EOF)
					fprintf(stderr, "Unable to close .graph file. Unable to save graph.\n");
			
			bigMod      = free_BigIntT(bigMod);
			bigModPower = free_BigIntT(bigModPower);
			one         = free_BigIntT(one);
			temp        = free_BigIntT(temp);
			temp2       = free_BigIntT(temp2);
			
			for (int i = 0; i < maxpower; i += 1)
			{
				for (int j = 0; j < numOfOrbits[i]; j += 1)
					orbitReps[i][j] = free_BigIntMatrixT(orbitReps[i][j]);
				
				FREE(orbitReps[i]);
				FREE(orbitLengths[i]);
			}
			FREE(orbitReps);
			FREE(orbitLengths);
			
			for (int i = 0; i < numOfOrbits[maxpower-1]; i += 1)
			{
				FREE(orbitMaps[i]);
			}
			FREE(orbitMaps);
			
			FREE(numOfOrbits);
			FREE(currOrbitMapMatch);
			FREE(currOrbitMapGroup);
			
			currVectElements = free_BigIntT_array(currVectElements, big_rows(A), 1);
			currVect  = free_BigIntMatrixT(currVect);
			tempVect  = free_BigIntMatrixT(tempVect);
			tempVect2 = free_BigIntMatrixT(tempVect2);
			theCycle  = free_CycleInfoT(theCycle);
			
			FREE(outputFileName);
			FREE(graphFileName);
			
			A = free_BigIntMatrixT(A);
		}
		
		
		//If the user wants to calculate the orbit maps for ALL lifts of a particular matrix
		else if (! strcmp(argv[1], "orbitmaps2"))
		{
			bool isNewVector;
			
			BigIntTP  baseMod; //The base prime
			BigIntTP* modList; //Holds every power of baseMod we'll need, from baseMod to baseMod^maxPower
			int maxPower;
			int currPower = 1; //For keeping track of which layer of moduli we're on
			
			int numArr[1] = {1};
			BigIntTP one;
			
			const BigIntMatrixTP initialA = UPDATEMATRIX;
			
			BigIntMatrixTP*  matrixLifts; //Holds the matrix lifts we're currently considering
			BigIntMatrixTP** orbitReps;   //Holds representatives for all the orbits, from smallest modulus to biggest
			int*  orbitRepsCount;         //How many vectors are in each layer of orbitReps?
			int** orbitRepsCycleLengths;  //What are the cycle lengths of the orbit reps? 
			
			BigIntTP**  vectorLiftElements; //Holds the raw elements responsible for "lifting" each vector
			BigIntTP*** matrixLiftElements; //Holds the raw elements for each matrix lift at each level
			
			BigIntMatrixTP currMat;
			
			const int numOfTempVects = 3;
			BigIntMatrixTP* tempVects;
			CycleInfoTP theCycle = NULL;
			
			if (big_rows(initialA) != big_cols(initialA))
			{
				fprintf(stderr, "Given update matrix isn't square.\n");
				FREE_VARIABLES;
				return EXIT_SUCCESS;
			}
			
			maxPower = (int)strtol(argv[2], &tempStr, 10);
			if (tempStr[0] != '\0')
			{
				fprintf(stderr, "Unable to read maxPower from command line.\n");
				FREE_VARIABLES;
				return EXIT_FAILURE;
			}
			
			if (argc > 3)
			{
				SET_BIG_NUM(argv[3], baseMod, "Unable to read modulus from command line.");
			}
			else
			{
				SET_BIG_NUM(bigintmodstring, baseMod, "Unable to read modulus from config file.");
			}
			
			matrixLiftElements = malloc((maxPower-1)*sizeof(BigIntTP**));
			for (int i = 0; i < maxPower-1; i += 1)
			{
				matrixLiftElements[i] = new_BigIntT_array(big_rows(initialA), big_rows(initialA));
				for (int row = 0; row < big_rows(initialA); row += 1)
					for (int col = 0; col < big_cols(initialA); col += 1)
						copy_BigIntT(big_element(initialA, row, col), matrixLiftElements[i][row][col]);
			}
			vectorLiftElements = new_BigIntT_array(big_rows(initialA), 1);
			currMat = new_BigIntMatrixT(big_rows(initialA), big_rows(initialA));
			
			orbitReps      = calloc(maxPower, sizeof(BigIntMatrixTP*));
			orbitRepsCount = calloc(maxPower, sizeof(int));
			orbitRepsCycleLengths = calloc(maxPower, sizeof(int*));
			
			one = new_BigIntT(numArr, 1);
			tempVects = malloc(numOfTempVects*sizeof(BigIntMatrixTP));
			for (int i = 0; i < numOfTempVects; i += 1)
				tempVects[i] = new_BigIntMatrixT(big_rows(initialA), 1);
			
			modList = malloc(maxPower*sizeof(BigIntTP));
			modList[0] = empty_BigIntT(1);
			copy_BigIntT(baseMod, modList[0]);
			for (int i = 1; i < maxPower; i += 1)
			{
				modList[i] = empty_BigIntT(1);
				multiply_BigIntT(modList[i-1], baseMod, modList[i]);
			}
			
			//Let's calculate our base orbit reps, since those won't change throughout
			// the course of the program running
			orbitRepsCount[0] = big_orbit_reps(initialA, baseMod, orbitReps, orbitRepsCycleLengths);
			
			/*
			for (int repr = 0; repr < orbitRepsCount[0]; repr += 1)
			{
				printbm_row(orbitReps[0][repr]);
				printf(" %d\n", orbitRepsCycleLengths[0][repr]);
			}
			*/
			
			//Now, let's start generating some orbit maps
			currPower += 1;
			while (currPower > 1)
			{
				currPower = maxPower; //This line is only for testing. Remove it when done
				
				set_big_matrix(currMat, matrixLiftElements[maxPower-2]);
				printbm(currMat);
				printf("\n");
				getchar();
				
				//Increment through our matrix lifts in such a way as to reduce
				// the amount of times we need to recompute lifts
				for (int level = maxPower-2; level >= 0; level -= 1)
				{
					//CLEAR VECTORS IN orbitReps[level+1]
					
					//If we're done looking at all the matrix lifts that keep the
					// lower vector lifts the same
					
					//NEED TO CREATE A NEW FUNCTION WHERE THE "CARRY" VALUE WHEN INCREMENTING
					// IS THE SAME AS THE INCREMENTATION VALUE
					if (increment_BigIntT_array(matrixLiftElements[level], 
					                            big_rows(initialA), 
																	    big_rows(initialA), 
																	    modList[level], 
																	    modList[level+1]))
					{
						currPower -= 1;
					}
					
					else
						break;
				}
				
				//Now, make sure that if we changed a lower-moduli matrix lift, we carry that
				// new lift up to the higher matrix lifts so the reduction process that keeps
				// lift vectors constant works properly.
				if (currPower > 1)
					for (int rewrite = 1; rewrite + currPower <= maxPower; rewrite += 1)
					{
						//printf(":) rewrite = %d, currPower = %d, maxPower = %d\n", rewrite, currPower, maxPower);
						for (int row = 0; row < big_rows(initialA); row += 1)
							for (int col = 0; col < big_rows(initialA); col += 1)
								copy_BigIntT(matrixLiftElements[currPower-2][row][col], 
							               matrixLiftElements[currPower-2+rewrite][row][col]);
					}
			}
			
			
			one     = free_BigIntT(one);
			baseMod = free_BigIntT(baseMod);
			
			for (int i = 0; i < maxPower-1; i += 1)
				matrixLiftElements[i] = free_BigIntT_array(matrixLiftElements[i], big_rows(initialA), big_rows(initialA));
			FREE(matrixLiftElements);
			
			vectorLiftElements = free_BigIntT_array(vectorLiftElements, big_rows(initialA), 1);
			currMat = free_BigIntMatrixT(currMat);
			
			for (int i = 0; i < maxPower; i += 1)
			{
				for (int j = 0; j < orbitRepsCount[i]; j += 1)
					orbitReps[i][j] = free_BigIntMatrixT(orbitReps[i][j]);
				FREE(orbitReps[i]);
				FREE(orbitRepsCycleLengths[i]);
			}
			FREE(orbitReps);
			FREE(orbitRepsCount);
			FREE(orbitRepsCycleLengths);
			
			theCycle  = free_CycleInfoT(theCycle);
			for (int i = 0; i < numOfTempVects; i += 1)
				tempVects[i] = free_BigIntMatrixT(tempVects[i]);
			FREE(tempVects);
		}
		
		
		//If the user wants to find an annihilating polynomial which
		// isn't a multiple of the minimal polynomial
		else if (! strcmp(argv[1], "oddminpolysearch"))
		{
			int matSize;
			int maxPower;
			int currPower;
			int polySize;
			BigIntTP bigMod = NULL;
			BigIntTP currMod = NULL;
			BigIntTP** currMatElements = NULL;
			BigIntMatrixTP currMat = NULL;
			
			bool foundMin;
			bool foundMonic;
			BigIntTP* minPolyCoeffs   = NULL;
			BigPolyTP minPoly         = NULL;
			BigPolyTP minMonicPoly    = NULL; //Non-monic annihilating polys can be smaller than minPoly
			BigIntTP* annilPolyCoeffs = NULL;
			BigPolyTP annilPoly       = NULL;
			
			BigPolyTP remainderPoly = NULL;
			BigPolyTP zeroPoly      = NULL;
			
			int oneArr[1] = {1};
			BigIntTP one;
			BigIntTP temp;
			BigIntTP coeffCompare;
			
			BigIntMatrixTP zeroMat = NULL;
			BigIntMatrixTP tempMat = NULL;
			
			char* outputFileName = NULL;
			FILE* outputFile = NULL;
			
			maxPower = (int)strtol(argv[2], &tempStr, 10);
			if (tempStr[0] != '\0')
			{
				fprintf(stderr, "Unable to read maxPower from command line.\n");
				FREE_VARIABLES;
				return EXIT_FAILURE;
			}
			
			matSize = (int)strtol(argv[3], &tempStr, 10);
			if (tempStr[0] != '\0')
			{
				fprintf(stderr, "Unable to read matrix size from command line.\n");
				FREE_VARIABLES;
				return EXIT_FAILURE;
			}
			
			polySize = (int)strtol(argv[4], &tempStr, 10);
			if (tempStr[0] != '\0')
			{
				fprintf(stderr, "Unable to read polynomial size from command line.\n");
				FREE_VARIABLES;
				return EXIT_FAILURE;
			}
			
			if (argc > 5)
			{
				SET_BIG_NUM(argv[5], bigMod, "Unable to read modulus from command line.");
			}
			else
			{
				SET_BIG_NUM(bigintmodstring, bigMod, "Unable to read modulus from config file.");
			}
			
			//Check to see if the user wants to resume computation
			if ((argc > 6) && (! strcmp(argv[6], "TRUE")))
			{
				currMat = read_BigIntMatrixT(resumefilepath);
				if (currMat == NULL)
				{
					fprintf(stderr, "Unable to read resume matrix from %s.\n", resumefilepath);
					FREE_VARIABLES;
					return EXIT_FAILURE;
				}
				
				if ((big_rows(currMat) != big_cols(currMat)) || (big_rows(currMat) != matSize))
				{
					fprintf(stderr, "Resume matrix given isn't the correct dimensions specified.\n");
					FREE_VARIABLES;
					currMat = free_BigIntMatrixT(currMat);
					bigMod = free_BigIntT(bigMod);
					
					return EXIT_FAILURE;
				}
				
				currMatElements = malloc(matSize*sizeof(BigIntTP*));
				for (int i = 0; i < matSize; i += 1)
				{
					currMatElements[i] = malloc(matSize*sizeof(BigIntTP));
					for (int j = 0; j < matSize; j += 1)
					{
						currMatElements[i][j] = empty_BigIntT(1);
						copy_BigIntT(big_element(currMat, i, j), currMatElements[i][j]);
					}
				}
			}
			
			else
			{
				currMat = new_BigIntMatrixT(matSize, matSize);
				
				currMatElements = malloc(matSize*sizeof(BigIntTP*));
				for (int i = 0; i < matSize; i += 1)
				{
					currMatElements[i] = malloc(matSize*sizeof(BigIntTP));
					for (int j = 0; j < matSize; j += 1)
						currMatElements[i][j] = empty_BigIntT(1);
				}
			}
			
			//If the user wants a file output
			if ((argc > 7) && (!strcmp(argv[7], "TRUE")))
			{
				outputFileName = malloc(MAXSTRLEN*sizeof(char));
				outputFileName[0] = '\0';
				
				//Let's hope there's enough space in the string for all this!
				strcat(outputFileName, "oddminpolysearch ");
				strcat(outputFileName, argv[2]);
				strcat(outputFileName, " ");
				strcat(outputFileName, argv[3]);
				strcat(outputFileName, " ");
				strcat(outputFileName, argv[4]);
				strcat(outputFileName, " ");
				append_BigIntT(outputFileName, bigMod);
				strcat(outputFileName, ".txt");	
				
				outputFile = fopen(outputFileName, "w");
				if (outputFile == NULL)
					fprintf(stderr, "Unable to create file %s. COntinuing without saving...\n", outputFileName);
				
			}
			
			one  = new_BigIntT(oneArr, 1);
			temp = empty_BigIntT(1);
			
			zeroPoly = constant_BigPolyT(temp);
			
			currMod = empty_BigIntT(1);
			copy_BigIntT(bigMod, currMod);
			
			minPolyCoeffs   = malloc(polySize*sizeof(BigIntTP));
			annilPolyCoeffs = malloc(polySize*sizeof(BigIntTP));
			for (int i = 0; i < polySize; i += 1)
			{
				minPolyCoeffs[i]   = empty_BigIntT(1);
				annilPolyCoeffs[i] = empty_BigIntT(1);
			}
			
			minPoly      = empty_BigPolyT();
			minMonicPoly = empty_BigPolyT();
			annilPoly    = empty_BigPolyT();
			
			tempMat = new_BigIntMatrixT(matSize, matSize);
			zeroMat = new_BigIntMatrixT(matSize, matSize);
			
			//Loops over all matrices under current modulus
			for (currPower = 1; currPower <= maxPower; currPower += 1)
			{
				if (outputFile != NULL)
				{
					fprintf(outputFile, "~~Mod ");
					fprinti(outputFile, currMod);
					fprintf(outputFile, "~~\n");
					if (fclose(outputFile) == EOF)
					{
						fprintf(stderr, "Unable to save output file. Continuing without saving...\n");
						outputFile = NULL;
					}
					else
					{
						outputFile = fopen(outputFileName, "a");
						if (outputFile == NULL)
							fprintf(stderr, "Unable to save output file. Continuing without saving...\n");
					}
				}
				
				do
				{
					set_big_matrix(currMat, currMatElements);
					resize_BigPolyT(minPoly, polySize);
					resize_BigPolyT(minMonicPoly, polySize);
					
					//Find minimal polynomial
					
					for (int i = 0; i < polySize; i += 1)
					{
						clear_BigIntT(minPolyCoeffs[i]);
						clear_BigIntT(annilPolyCoeffs[i]);
					}
					
					foundMin   = FALSE;
					foundMonic = FALSE;
					while ((!foundMin) || (!foundMonic))
					{
						for (int i = 0; i < polySize; i += 1)
						{
							add_BigIntT(minPolyCoeffs[i], one, temp);
							if (compare_BigIntT(temp, currMod) == 0)
								clear_BigIntT(minPolyCoeffs[i]);
							else
							{
								copy_BigIntT(temp, minPolyCoeffs[i]);
								break;
							}
						}
						
						if (!foundMin)
						{
							set_BigPolyT(minPoly, minPolyCoeffs);
							set_BigPolyT(minMonicPoly, minPolyCoeffs);
							reduce_BigPolyT(minPoly);
							eval_BigPolyT(minPoly, currMat, tempMat, currMod);
							
							resize_BigPolyT(minPoly, polySize);
						}
						else if (!foundMonic)
						{
							set_BigPolyT(minMonicPoly, minPolyCoeffs);
							reduce_BigPolyT(minMonicPoly);
							eval_BigPolyT(minMonicPoly, currMat, tempMat, currMod);
							
							resize_BigPolyT(minMonicPoly, polySize);
						}
						
						//If we've found our minimum polynomial for currMat
						if (compare_BigIntMatrixT(tempMat, zeroMat))
						{
							if (!foundMin)
							{
								reduce_BigPolyT(minPoly);
								foundMin = TRUE;
							}
							
							//Check to see if polynomial is monic
							if (!foundMonic)
							{
								for (int i = polySize-1; i >= 0; i -= 1)
								{
									if (!is_zero(minPolyCoeffs[i]))
									{
										if (compare_BigIntT(one, minPolyCoeffs[i]) == 0)
										{
											reduce_BigPolyT(minMonicPoly);
											foundMonic = TRUE;
										}
										
										break;
									}
								}
							}
						}
					}
					
					//For checking annihilating polynomials, we'll only check
					// monic polynomials of size polySize
					for (int c = 0; c < polySize-1; c += 1)
						copy_BigIntT(currMod, annilPolyCoeffs[c]);
					
					//We've found our minimal polynomial(s)
					//Now to see if all other annihilating polynomials are
					// multiples of them.
					
					//foundMin and foundMonic are now gonna be used to keep track of
					// whether we've found a polynomial that's annihilating
					foundMin = FALSE;
					foundMonic = FALSE;
					
					remainderPoly = empty_BigPolyT();
					coeffCompare = empty_BigIntT(1);
					
					//First, we should check whether our first monic polynomial divides the first annihilating polynomial
					copy_BigPolyT(minMonicPoly, annilPoly);
					while ((!foundMin) || (!foundMonic))
					{
						//Check to see if annilPoly is an annihilating polynomial
						eval_BigPolyT(annilPoly, currMat, tempMat, currMod);
						if (compare_BigIntMatrixT(tempMat, zeroMat))
						{
							if (!foundMin)
							{
								//Check if annilPoly is a multiple of minPoly
								divide_BigPolyT(annilPoly, minPoly, NULL, remainderPoly, currMod);
								if (compare_BigPolyT(remainderPoly, zeroPoly) != 0)
									foundMin = TRUE;
							}
							
							if (!foundMonic)
							{
								//Check if annilPoly is a multiple of minMonicPoly
								divide_BigPolyT(annilPoly, minMonicPoly, NULL, remainderPoly, currMod);
								if (compare_BigPolyT(remainderPoly, zeroPoly) != 0)
									foundMonic = TRUE;
							}
							
							if ((foundMin) && (foundMonic))
							{
								printf("\ncurrMat:\n");
								printbm(currMat);
								printf("Minimal poly: ");
								printp(minPoly);
								printf("\nMonic poly: ");
								printp(minMonicPoly);
								printf("\n");
								printp(annilPoly);
								printf(" is not a multiple of ");
								printp(minPoly);
								printf("\n");
								printp(annilPoly);
								printf(" is not a multiple of ");
								printp(minMonicPoly);
								printf("\n");
								
								if (outputFile != NULL)
								{
									fprintbm(outputFile, currMat);
									fprintf(outputFile, "Minimal poly: ");
									fprintp(outputFile, minPoly);
									fprintf(outputFile, "\nMonic poly: ");
									fprintp(outputFile, minMonicPoly);
									fprintf(outputFile, "\n");
									fprintp(outputFile, minPoly);
									fprintf(outputFile, " !| ");
									fprintp(outputFile, annilPoly);
									fprintf(outputFile, "\n");
									fprintp(outputFile, minMonicPoly);
									fprintf(outputFile, " !| ");
									fprintp(outputFile, annilPoly);
									fprintf(outputFile, "\n\n");
									
									if (fclose(outputFile) == EOF)
									{
										fprintf(stderr, "Unable to save output file. Continuing without saving...\n");
										outputFile = NULL;
									}
									else
									{
										outputFile = fopen(outputFileName, "a");
										if (outputFile == NULL)
											fprintf(stderr, "Unable to save output file. Continuing without saving...\n");
									}
								}
							}
						}
						
						//Increment annilPolyCoeffs
						for (int i = 0; i < polySize; i += 1)
						{
							add_BigIntT(annilPolyCoeffs[i], one, temp);
							if (i != polySize-1)
								copy_BigIntT(currMod, coeffCompare);
							else
								add_BigIntT(one, one, coeffCompare);
							
							if (compare_BigIntT(temp, coeffCompare) >= 0)
							{
								clear_BigIntT(annilPolyCoeffs[i]);
								if (i == polySize-1)
								{
									foundMin = TRUE;
									foundMonic = TRUE;
								}
							}
							else
							{
								copy_BigIntT(temp, annilPolyCoeffs[i]);
								break;
							}
						}
						
						//Geez, I love spending CPU cycles on this :)
						resize_BigPolyT(annilPoly, polySize);
						set_BigPolyT(annilPoly, annilPolyCoeffs);
						reduce_BigPolyT(annilPoly);
					}
					
					remainderPoly = free_BigPolyT(remainderPoly);
					coeffCompare  = free_BigIntT(coeffCompare);
				}
				while (!increment_BigIntT_array(currMatElements, matSize, matSize, one, currMod));
				
				//Get next power of the modulus
				multiply_BigIntT(currMod, bigMod, temp);
				copy_BigIntT(temp, currMod);
			}
			
			if (outputFile != NULL)
			{
				if (fclose(outputFile) == EOF)
					fprintf(stderr, "Output file may not have saved correctly.\n");
			}
			
			bigMod  = free_BigIntT(bigMod);
			currMatElements = free_BigIntT_array(currMatElements, matSize, matSize);
			currMat = free_BigIntMatrixT(currMat);
			
			one          = free_BigIntT(one);
			temp         = free_BigIntT(temp);
			currMod      = free_BigIntT(currMod);
			//coeffCompare = free_BigIntT(coeffCompare);
			
			tempMat = free_BigIntMatrixT(tempMat);
			zeroMat = free_BigIntMatrixT(zeroMat);
			
			for (int i = 0; i < polySize; i += 1)
			{
				minPolyCoeffs[i]   = free_BigIntT(minPolyCoeffs[i]);
				annilPolyCoeffs[i] = free_BigIntT(annilPolyCoeffs[i]);
			}
			FREE(minPolyCoeffs);
			FREE(annilPolyCoeffs);
			
			minPoly      = free_BigPolyT(minPoly);
			minMonicPoly = free_BigPolyT(minMonicPoly);
			annilPoly    = free_BigPolyT(annilPoly);
			
			//remainderPoly = free_BigPolyT(remainderPoly);
			zeroPoly      = free_BigPolyT(zeroPoly);
			
			if (outputFileName != NULL)
			{
				FREE(outputFileName);
			}
		}
	}

	
	else
	{
		printf(ANSI_COLOR_GREEN "LINCELLAUT by Zach Strong.\n" ANSI_COLOR_RESET);
		printf("Usage: lincellaut <tool> [options]\n\n");
		printf("Tools:\n");
		
		printf(" - " ANSI_COLOR_YELLOW "iterate " ANSI_COLOR_CYAN "[iterations]" ANSI_COLOR_RESET "\n");
		printf(" - " ANSI_COLOR_YELLOW "inverse " ANSI_COLOR_CYAN "[modulus]" ANSI_COLOR_RESET "\n");
		printf(" - " ANSI_COLOR_YELLOW "chara " ANSI_COLOR_CYAN "[modulus]" ANSI_COLOR_RESET "\n");
		printf(" - " ANSI_COLOR_YELLOW "allcharas " ANSI_COLOR_CYAN "coeffs..." ANSI_COLOR_RESET "\n");
		printf(" - " ANSI_COLOR_YELLOW "core " ANSI_COLOR_CYAN "[modulus]" ANSI_COLOR_RESET "\n");
		printf(" - " ANSI_COLOR_YELLOW "orbits " ANSI_COLOR_CYAN "[modulus] [fileoutput]" ANSI_COLOR_RESET "\n");
		printf(" - " ANSI_COLOR_YELLOW "splitorbits " ANSI_COLOR_CYAN "[modulus] [fileoutput] " ANSI_COLOR_RESET "\n");
		printf(" - " ANSI_COLOR_YELLOW "orbitreps " ANSI_COLOR_CYAN "[modulus] [fileoutput]" ANSI_COLOR_RESET "\n");
		printf(" - " ANSI_COLOR_YELLOW "branchreps " ANSI_COLOR_CYAN "[modulus]" ANSI_COLOR_RESET "\n");
		printf(" - " ANSI_COLOR_YELLOW "orbitspaces " ANSI_COLOR_CYAN "[modulus] [minpolys] [fileoutput]" ANSI_COLOR_RED " (UNFINISHED)" ANSI_COLOR_RESET "\n");
		printf(" - " ANSI_COLOR_YELLOW "floyd " ANSI_COLOR_CYAN "[modulus]" ANSI_COLOR_RESET "\n");
		printf(" - " ANSI_COLOR_YELLOW "cycmatsearch " ANSI_COLOR_CYAN "resume size maxmod cycles..." ANSI_COLOR_RESET "\n");
		printf(" - " ANSI_COLOR_YELLOW "cycconvmat " ANSI_COLOR_CYAN "from to [mod]" ANSI_COLOR_RESET "\n");
		printf(" - " ANSI_COLOR_YELLOW "ccmzerosearch " ANSI_COLOR_CYAN "resume size [mod]" ANSI_COLOR_RESET "\n");
		printf(" - " ANSI_COLOR_YELLOW "vectprops " ANSI_COLOR_CYAN "[mod] [resume] [fileoutput]" ANSI_COLOR_RESET "\n");
		printf(" - " ANSI_COLOR_YELLOW "vectpolys " ANSI_COLOR_CYAN "degree [mod] [fileoutput]" ANSI_COLOR_RED " (UNFINISHED)" ANSI_COLOR_RESET "\n");
		printf(" - " ANSI_COLOR_YELLOW "matprops " ANSI_COLOR_CYAN "maxpower [modulus]" ANSI_COLOR_RESET "\n");
		printf(" - " ANSI_COLOR_YELLOW "charawalk" ANSI_COLOR_CYAN " step [modulus]" ANSI_COLOR_RESET "\n");
		printf(" - " ANSI_COLOR_YELLOW "fibmultsearch " ANSI_COLOR_CYAN "[bound]" ANSI_COLOR_RESET "\n");
		printf(" - " ANSI_COLOR_YELLOW "dynamics " ANSI_COLOR_CYAN "[maxPower] [modulus] [allConfigs] [fileoutput]" ANSI_COLOR_RESET "\n");
		printf(" - " ANSI_COLOR_YELLOW "orbitmaps " ANSI_COLOR_CYAN "maxPower [modulus] [fileoutput]" ANSI_COLOR_RESET "\n");
		printf(" - " ANSI_COLOR_YELLOW "orbitmaps2 " ANSI_COLOR_CYAN "maxPower [modulus] [fileoutput]" ANSI_COLOR_RED " (UNFINISHED)" ANSI_COLOR_RESET "\n");
		printf(" - " ANSI_COLOR_YELLOW "oddminpolysearch " ANSI_COLOR_CYAN "maxPower size polysize [modulus] [resume] [fileoutput]" ANSI_COLOR_RESET "\n");
		
		printf("\nFor a more complete description of LINCELLAUT's usage, " \
		"refer to the included documentation.\n");
	}
	
	//Freeing memory
	FREE_VARIABLES;
	return EXIT_SUCCESS;
}
