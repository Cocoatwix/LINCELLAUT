
/* A simple program for calculating various things about linear
 *  cellular automaton.
 *
 * Mar 27, 2022
 *
 */
 
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "headers/linalg.h" 
#include "headers/cycles.h" //Allows us to use Floyd's Algorithm

#define FREE(v) free(v); v = NULL

/* Floyd's Cycle Detection Algorithm
 * The algorithm will have a stopping time of:
 *  tau (transient length) if tau == 0 mod omega (cycle length)
 *  tau + (omega - (tau mod omega)) otherwise.
 
 
 tau + omega - (tau mod omega) ~~ -omega
 tau - (tau mod omega)
 n*omega + N - N
 n*omega ~~ /omega
 n

 however
 
 n*omega ~~ -omega
 (n-1)*omega ~~ /omega
 (n-1)
 
 maybe we run floyd's cycle algorithm again with different step sizes (+2 and +3?),
 figure out what the stopping time for the new setup is,
 and creating a system of equations to solve for omega and tau?
 
 *
 * If tau < omega, then stopping time = omega
 * Every option for the stopping time is greater than or equal to omega,
 *  so the stopping time gives an upper bound on omega (omega can be no greater than the stopping time).
 *
 * We can also say that if St > 2tau, then St = omega.
 * 2tau < tau + (omega - (tau mod omega)) and rearrange inequality.
 *  You'll find that tau must be less than omega in this case.
 */
 
/* Once we have x_n in the cycle, if we can't conclude omega,
 *  then we can set y = x_n and keep iterating y until it equals
 *  x_n again. The number of steps this takes will give us omega.
 */
 
/* We also know that if our update rule F doesn't have an inverse,
 *  then our automaton has transient regions. This is because for finite sets
 *  that map to themselves, being one-to-one is the same as onto (think about it!),
 *  and if F isn't one-to-one, then there's some configuration that doesn't get mapped to,
 *  hence it would be a transient configuration.
 *
 */


int main()
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
	
	//int   iterations;
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
	
	#ifdef VERBOSE
	printf("Modulus: %d\n", modulus);
	printf("Update: %s\n", updatefilepath);
	printf("Initial: %s\n", initialfilepath);
	#endif //VERBOSE
	
	//Update rule matrix
	IntMatrixTP F     = read_IntMatrixT(updatefilepath);
	//IntMatrixTP Finv;
	//IntMatrixTP Fmult = new_IntMatrixT(rows(F), cols(F));
	
	//Stores our initial vector
	//IntMatrixTP s_0 = read_IntMatrixT(INITIALFILEPATH);
	
	//IntMatrixTP s_f; //Stores our final vector
	
	/* Iterate s_0 a few times
	s_f = iterate(F, s_0, MODULUS, ITERATIONS);
	printm(s_f); */
	
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
		printf("F is not invertible mod %d.\n", MODULUS);
	
	printcycle(floyd(F, s_0, MODULUS)); */
	
	write_iteration(iterfilepath, F, modulus);
	
	//Freeing memory
	FREE(updatefilepath);
	FREE(initialfilepath);
	FREE(iterfilepath);
	
	F = free_IntMatrixT(F);
	//Finv = Finv != NULL ? free_IntMatrixT(Finv) : NULL;
	//s_0 = free_IntMatrixT(s_0);
	//s_f = free_IntMatrixT(s_f);
	
	return EXIT_SUCCESS;
}