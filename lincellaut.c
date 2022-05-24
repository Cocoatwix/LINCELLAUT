
/* A simple program for calculating various things about linear
 *  cellular automaton.
 *
 * Mar 27, 2022
 *
 */
 
/* The following resources were used as a reference:
https://docs.microsoft.com/en-us/cpp/c-language/cpp-integer-limits?view=msvc-170
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
	char* updatefilepath  = malloc(MAXSTRLEN*sizeof(char));
	char* initialfilepath = malloc(MAXSTRLEN*sizeof(char));
	char* iterfilepath    = malloc(MAXSTRLEN*sizeof(char));
	
	//Temporarily holds config data
	char* systemData = malloc(MAXSTRLEN*sizeof(char));
	
	
	while (fscanf(system, "%10s", systemData) == 1)
	{
		//printf("%s\n", systemData);
		if (! strcmp(systemData, "mod"))
		{
			if (fscanf(system, "%d", &modulus) != 1)
			{
				fprintf(stderr, "Unable to read modulus from config file.\n");
				return EXIT_FAILURE;
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
				char* tempStr;
				
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
	}
	
	else
	{
		printf("Usage: lincellaut <tool> [options]\n\n");
		printf("Tools:\n");
		printf(" - iterate [iterations]: Iterate the update matrix a given number of times.\n");
		printf("   - iterations: Overrides the number of iterations provided in the .config file.\n\n");
		
		printf(" - fibcycle: Generate the Fibonacci cycle that contains the initial vector.\n\n");
		
		printf("For a more complete description of LINCELLAUT's usage, refer to the included documentation.\n");
	}
	
	//Testing arbitrary precision stuff
	int bunch1[] = {789, 456, 123};
	int bunch2[] = {4, 9653, 23, 45};
	BigIntTP num1 = new_BigIntT(bunch1, 3);
	BigIntTP num2 = new_BigIntT(bunch2, 4);
	printi(num1);
	printf("\n");
	printi(num2);
	printf("\n");
	
	num1 = free_BigIntT(num1);
	num2 = free_BigIntT(num2);
	
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