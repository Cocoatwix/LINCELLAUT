
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
			
			IntMatrixTP F   = read_IntMatrixT(updatefilepath);
			IntMatrixTP F_2 = identity_IntMatrixT(rows(F));
			
			//Iterate F a few times
			//The minus 1 is for easier conversion between ORBITVIS results
			if (iterations > 0)
			{
				F_2 = iterate(F, F, modulus, iterations-1);
				printf("Iterations: %d\n", iterations);
			}
			printm(F_2, TRUE);

			F   = free_IntMatrixT(F);
			F_2 = free_IntMatrixT(F_2);
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
			printf("Note: for testing purposes this tool currently sets" \
			" its starting value at a value other than 1.\n");
			
			int hundred[] = {100};
			int oneArr[]  = {001};
			int zeroArr[] = {000};
			int start[] = {862};
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
	
	//Testing arbitrary precision stuff
	/*int test1[] = {2, 8, 1};
	int test2[] = {5, 8};
	BigIntTP tp1 = new_BigIntT(test1, 3);
	BigIntTP tp2 = new_BigIntT(test2, 2);
	BigIntTP tp3 = empty_BigIntT(3);
	subtract_BigIntT(tp1, tp2, tp3);
	printi(tp3);
	printf("\n");
	
	tp1 = free_BigIntT(tp1);
	tp2 = free_BigIntT(tp2);
	tp3 = free_BigIntT(tp3); */
	
	/*int bunch1[]   = {789, 456, 3};
	int bunch2[]   = {4, 9653, 23, 45};
	int modBunch[] = {4237, 44};

	BigIntTP num1     = new_BigIntT(bunch1, 3);
	BigIntTP num2     = new_BigIntT(bunch2, 4);
	BigIntTP tempNum  = empty_BigIntT(4);
	BigIntTP tempNum2 = empty_BigIntT(4);
	BigIntTP modNum   = new_BigIntT(modBunch, 2);
	
	printi(num1);
	printf("\n");
	printi(num2);
	printf("\n");
	printi(modNum);
	printf("\n");
	
	printf("num2 - num1 = ");
	subtract_BigIntT(num2, num1, tempNum);
	printi(tempNum);
	printf("\n");
	
	printf("num2 + num1 = ");
	add_BigIntT(num2, num1, tempNum);
	printi(tempNum);
	printf("\n");
	
	printf("num1 mod modNum = ");
	mod_BigIntT(num1, modNum, tempNum);
	printi(tempNum);
	printf("\n");
	
	printf("num2 mod modNum = ");
	mod_BigIntT(num2, modNum, tempNum2);
	printi(tempNum2);
	printf("\n");
	
	num1 = free_BigIntT(num1);
	num2 = free_BigIntT(num2);
	modNum = free_BigIntT(modNum);
	tempNum = free_BigIntT(tempNum);
	tempNum2 = free_BigIntT(tempNum2); */
	
	/*
	//Seeing if we can find a relationship between the moduli 
	// that have specific rotation matrices
	
	printf("45deg\t30deg\t15deg\n");
	for (int x = 3; x < 1000; x += 2)
	{
		//45degrot
		if (square_root(2, x) != -1)
			printf("%2d\t", x);
		else
			printf("..\t");
		
		//30degrot
		if (square_root(3, x) != -1)
			printf("%2d\t", x);
		else
			printf("..\t");
		
		//15degrot
		if ((square_root(3, x) != -1) &&
		    (square_root(2, x) != -1))
			printf("%2d\t", x);
		else
			printf("..\t");
		
		printf("\n");
	}
	*/

	/*
	IntMatrixTP F_3   = new_IntMatrixT(2, 2);
	IntMatrixTP Finv;
	IntMatrixTP I     = identity_IntMatrixT(2);
	//IntMatrixTP Fmult = new_IntMatrixT(rows(F), cols(F)); */
	
	//Stores our initial vector
	//IntMatrixTP s_0 = read_IntMatrixT(initialfilepath);
	
	//IntMatrixTP s_f; //Stores our final vector
	
	//Iterate until we get to the identity or run out of iterations
	/*Finv = inverse(F, modulus);
	copy_IntMatrixT(I, F_2);
	
	if (Finv != NULL)
	{
		iterations = 0;
		
		do
		{
			mat_mul(F, F_2, F_3);
			modm(F_3, modulus);
			copy_IntMatrixT(F_3, F_2);
			
			iterations += 1;
		}
		while (! compare_IntMatrixT(F_2, I));
	}
	
	else
	{
		for (int i = 0; i < iterations; i += 1)
		{
			mat_mul(F, F_2, F_3);
			modm(F_3, modulus);
			copy_IntMatrixT(F_3, F_2);
		}
	}
	printf("Modulus: %d\n", modulus);
	printf("Matrix:\n");
	printm(F, TRUE);
	printf("Iterated matrix:\n");
	printm(F_2, TRUE);
	printf("Iterations: %d\n", iterations); */
	
	//See which points F visits in its orbit
	//INT_MAX
	//printf("Modulus: %d\n", modulus);
	//visit_points(F, modulus, iterations);
	
	//Testing our ability to find eigenvalues
	/*int* values = eigenvalues(F, modulus);
	if (values == NULL)
		printf("No eigenvalues exist for the given system.\n");
	
	else
	{
		printf("Eigenvalues: ");
		for (int i = 1; i <= values[0]; i += 1)
		{
			if (i == values[0])
				printf("%d\n", values[i]);
			else
				printf("%d, ", values[i]);
		}
	}
	
	//Now testing our ability to create eigenvectors
	printf("Matrix:\n");
	printm(F, TRUE);
	printf("Using eigenvalue %d:\n", values[1]);
	IntMatrixTP E = eigenvector(F, values[1], modulus);
	//printm(E, TRUE);
	
	FREE(values);
	//E = free_IntMatrixT(E); */
	
	//Testing the determinant function
	//printf("Determinant of update matrix: %d\n", det(F));
	
	//Testing the inverse function
	/* printf("F:\n");
	printm(F, TRUE);
	
	Finv = inverse(F, MODULUS); 
	if (Finv != NULL)
	{
		printf("The inverse of F is:\n");
		printm(Finv, TRUE);
		
		//Testing both orders to see if the inverse really is the inverse
		printf("F and Finv multipled together give:\n");
		mat_mul(F, Finv, Fmult); modm(Fmult, MODULUS); printm(Fmult, TRUE);
		mat_mul(Finv, F, Fmult); modm(Fmult, MODULUS); printm(Fmult, TRUE);
	}
	else
		printf("F is not invertible mod %d.\n", MODULUS); */
	
	/*CycleInfoTP theCycle = floyd(F, s_0, modulus);
	printcycle(theCycle);
	theCycle = free_CycleInfoT(theCycle); */
	
	/*
	//Generating numbers for rotation matrices
	int numInv = (modulus - num_inverse(2, modulus)) % modulus;
	int sqrt   = square_root(3, modulus);
	printf("2^-1 = %d\n", num_inverse(2, modulus));
	printf("sqrt(3) = %d\n", square_root(3, modulus));
	printf("(2^-1) * sqrt(3) = %d\n", (num_inverse(2, modulus)*square_root(3, modulus)) % modulus);
	printf("-(2^-1) = %d\n", (modulus - num_inverse(2, modulus)) % modulus);
	
	//30degrot
	printf("%d %d\n", 
	(numInv*sqrt) % modulus,
	((modulus - numInv) % modulus));
	printf("%d %d\n",
	numInv,
	(numInv*sqrt) % modulus);
	
	//45degrot
	printf("%d %d\n", 
	(num_inverse(2, modulus)*square_root(2, modulus)) % modulus,
	modulus - ((num_inverse(2, modulus)*square_root(2, modulus)) % modulus));
	printf("%d %d\n",
	(num_inverse(2, modulus)*square_root(2, modulus)) % modulus,
	(num_inverse(2, modulus)*square_root(2, modulus)) % modulus);
	
	*/	
	//Freeing memory
	FREE(updatefilepath);
	FREE(initialfilepath);
	FREE(iterfilepath);
	
	/*
	F_3 = free_IntMatrixT(F_3);
	Finv = Finv != NULL ? free_IntMatrixT(Finv) : NULL;
	I = free_IntMatrixT(I); */
	//s_0 = free_IntMatrixT(s_0);
	//s_f = free_IntMatrixT(s_f);
	
	return EXIT_SUCCESS;
}