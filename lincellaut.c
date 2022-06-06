
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

#define FREE(v) free(v); v = NULL

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"


int strtoBIT(char* numStr, BigIntTP* theBig)
/** Takes a numerical string and creates a BigIntT
    struct using it, storing it in theBig. 
		
		This function assumes theBig has been declared.
		
		Returns 1 on success, 0 otherwise. */
{
	int substrstart;      //Holds the start of the substring in the for-loop below
	int substrlen;        //Holds how long each substring should be
	int bunchCounter = 0; //For properly storing bunches
	int numStrLength = strlen(numStr);
	
	int*  bunches; //Holds the BigIntT bunches of our number
	char* tempStr; //Holds info about whether the number was read correctly
	char* substr = malloc(5*sizeof(char)); //Holds the substring for each bunch
	
	bunches = malloc((numStrLength/4 + 1)*sizeof(int));
	
	//Note that the constants used in this loop would
	// need to change if we ever change the bunch size for
	// BigIntT structs. We use 4 because the bunch size is
	// currently 9999, or 4 base 10 digits.
	for (int bunch = numStrLength; 
	bunch >= (numStrLength % 4 == 0) ? 1 : 0; //This prevents an extra bunch from being read
	bunch -= 4)
	{
		//Getting the correct substring for the next bunch
		substrstart = (bunch-4 >= 0) ? bunch - 4 : 0;
		substrlen   = (bunch-4 >= 0) ? 4 : numStrLength % 4;
		
		strncpy(substr, numStr+substrstart, substrlen);
		substr[substrlen] = '\0'; //Adding null byte manually
		
		//Store bunches
		bunches[bunchCounter] = (int)strtol(substr, &tempStr, 10);
		bunchCounter += 1;
		
		//If we read an invalid character
		if (tempStr[0] != '\0')
		{
			FREE(bunches);
			FREE(substr);
			return 0;
		}
	}
	
	//Now, we actually create the BigIntT
	*theBig = new_BigIntT(bunches, bunchCounter);
	
	FREE(bunches);
	FREE(substr);
	
	return 1;
}


int main(int argc, char* argv[])
{
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
					F_2 = identity_IntMatrixT(rows(F));
				}
			}
			
			else
			{
				F   = read_IntMatrixT(updatefilepath);
				F_2 = identity_IntMatrixT(rows(F));
				
				if (F == NULL)
					return EXIT_FAILURE;
			}
			
			//Iterate F a few times
			//The minus 1 is for easier conversion between ORBITVIS results
			printf("Iterations: %d\n", iterations);
			if (iterations > 0)
				F_2 = iterate(F, F, modulus, iterations-1);
			
			//Prevents the matrix from being printed when
			// an inverse doesn't exist and the iterations
			// provided was negative
			if (iterations >= 0)
				printm(F_2, TRUE);

			F   = free_IntMatrixT(F);
			F_2 = free_IntMatrixT(F_2);
		}
		
		
		//Find the inverse of the update matrix
		else if (!strcmp(argv[1], "inverse"))
		{
			IntMatrixTP F = read_IntMatrixT(updatefilepath);
			IntMatrixTP Finv;
			
			if (F == NULL)
				return EXIT_FAILURE;
			
			printf("Update matrix:\n");
			printm(F, TRUE);
			
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
				printm(Finv, TRUE);
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
		": Iterate the update matrix a given number of times.\n");
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
		
		printf(" - " ANSI_COLOR_YELLOW "rots " ANSI_COLOR_CYAN "[modulus]" ANSI_COLOR_RESET \
		": Finds and outputs some basic rotation matrices for the given modulus.\n");
		printf("   - " ANSI_COLOR_CYAN "modulus" ANSI_COLOR_RESET \
		": Overrides the modulus provided in the .config file.\n\n");
		
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
	}

	/*
	IntMatrixTP F_3   = new_IntMatrixT(2, 2);
	IntMatrixTP I     = identity_IntMatrixT(2);
	//IntMatrixTP Fmult = new_IntMatrixT(rows(F), cols(F)); */
	
	//Stores our initial vector
	//IntMatrixTP s_0 = read_IntMatrixT(initialfilepath);
	
	//IntMatrixTP s_f; //Stores our final vector
	
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
	
	return EXIT_SUCCESS;
}