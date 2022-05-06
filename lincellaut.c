
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
	IntMatrixTP F_2;
	//IntMatrixTP Finv;
	//IntMatrixTP Fmult = new_IntMatrixT(rows(F), cols(F));
	
	//Stores our initial vector
	//IntMatrixTP s_0 = read_IntMatrixT(initialfilepath);
	
	//IntMatrixTP s_f; //Stores our final vector
	
	//Iterate s_0 a few times
	F_2 = iterate(F, F, modulus, 155);
	printm(F_2, TRUE);
	
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
	
	//write_iteration(iterfilepath, F, modulus);
	
	//Freeing memory
	FREE(updatefilepath);
	FREE(initialfilepath);
	FREE(iterfilepath);
	
	F = free_IntMatrixT(F);
	F_2 = free_IntMatrixT(F_2);
	//Finv = Finv != NULL ? free_IntMatrixT(Finv) : NULL;
	//s_0 = free_IntMatrixT(s_0);
	//s_f = free_IntMatrixT(s_f);
	
	return EXIT_SUCCESS;
}