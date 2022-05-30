
/*
A companion library for
lincellaut.c to explore a
specific LCA system, mainly

2 2
1 1
1 0

in .matrix format.
May 24, 2022
*/

#include <stdlib.h>
#include <stdio.h>

#include "../headers/linalg.h" //For getting digit counts

int generate_orbit(int* v, int modulus, int toPrint)
/** Generates/prints the orbit of the given vector (not in
    IntMatrixT form) under the given modulus and returns
		the length of the orbit. 
		
		This function assumes the vectors have 2 elements. 

		If toPrint == 0, the orbit isn't printed. 
		Else, it is. */
{
	int orbitLength = 0;
	int modDigits   = num_digits(modulus);
	
	int* initial = malloc(2*sizeof(int));
	initial[0] = v[0];
	initial[1] = v[1];
	
	do
	{
		//Fibonacci iteration
		v[0] = v[0] + v[1];
		v[1] = v[0] - v[1];
		v[0] = v[0] % modulus;
		v[1] = v[1] % modulus;
		orbitLength += 1;
		
		if (toPrint)
		{
			printf("[");
			
			for (int x = 0; x < modDigits-num_digits(v[0]); x += 1)
				printf("0");
			printf("%d, ", v[0]);
			
			for (int x = 0; x < modDigits-num_digits(v[1]); x += 1)
				printf("0");
			printf("%d] ", v[1]);
		}
	}
	while ((v[0] != initial[0]) || (v[1] != initial[1]));
	printf("\n");
	
	free(initial);
	initial = NULL;
	
	return orbitLength;
}
