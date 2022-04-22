
#include <stdlib.h>
#include <stdio.h>

#include "../headers/linalg.h" //Letting us use matrices

//Might need to write a "free" function later on
//May want to store the vector that's in the cycle as well
typedef struct cycleinfo
/** Struct to hold any important info about a structure. */
{
	int omega; //Cycle length
	int tau;   //Transient length
	int stoppingTime;
}
CycleInfoT, *CycleInfoTP;


//Not fully implemented yet
void printcycle(CycleInfoTP c)
/** Prints info about a cycle to the screen. */
{
	printf("Floyd's Algorithm stopping time: %d\n", c->stoppingTime);
	printf("Cycle length: %d\n", c->omega);
	
	//-1 is used when we don't know tau
	if (c->tau != -1)
		printf("Transient region length: %d\n", c->tau);
	else
		printf("Transient region length: ??\n");
}


IntMatrixTP iterate(IntMatrixTP F, IntMatrixTP s_0, int modulus, int iterations)
/** Iterates a given vector under some update rule n times.
    Returns a pointer to the result of the iteration. */
{
	IntMatrixTP s_o = new_IntMatrixT(rows(s_0), cols(s_0));
	IntMatrixTP s_e = new_IntMatrixT(rows(s_0), cols(s_0));
	
	copy_IntMatrixT(s_0, s_o);
	
	//The actual iteration
	for (int i = 0; i < iterations; i += 1)
	{
		if (i % 2 == 0)
		{
			mat_mul(F, s_o, s_e);
			modm(s_e, modulus);
		}
		
		else
		{
			mat_mul(F, s_e, s_o);
			modm(s_o, modulus);
		}
	}
	
	//Depending on which vector our result is in, return the right vector
	if (iterations % 2 == 0)
	{
		s_e = free_IntMatrixT(s_e);
		return s_o;
	}
	
	else
	{
		s_o = free_IntMatrixT(s_o);
		return s_e;
	}
}


int write_orbits(const char* fileName, IntMatrixTP F, int modulus)
/** Writes out a .orbit file containing every single orbit in the given
    system. Returns 1 on success, 0 otherwise. */
{
	//F must be a square matrix
	if (rows(F) != cols(F))
		return 0;
	
	FILE* orbitFile = fopen(fileName, "w");
	
	if (orbitFile == NULL)
		return 0;
	
	IntMatrixTP s = new_IntMatrixT(rows(F), 1);
	IntMatrixTP s_temp = new_IntMatrixT(rows(F), 1);
	
	//This will hold all of our starting values for our vectors
	int* vectors = calloc(rows(F), sizeof(int));
	
	//lineNumbers the line numbers where you can find information about
	// specific vectors in the .orbit file.
	// For example, element lineNumbers[1][1] would hold the line
	// number where you can look up <1, 1>'s info.
	// Of course, this can only really be done when F is 2x2.
	int** lineNumbers;
	
	if (rows(F) == 2)
	{
		lineNumbers = malloc(modulus*sizeof(int*));
		for (int i = 0; i < modulus; i += 1)
			lineNumbers[i] = calloc(modulus, sizeof(int));
	}
	
	bool allTested = FALSE;
	
	
	//Loop until we've found the orbits of all vectors
	while (!allTested)
	{
		//Now, I need a way to be able to tell how many
		// matrix multiplications to do before going to the
		// next vector. How do I know when we've reached a
		// cycle?
		set_column(s, vectors);
		
		//Increment our starting vector
		allTested = TRUE;
		for (int element = rows(F)-1; element >= 0; element -= 1)
		{
			//This logic counts through all possible vectors mod modulus
			if (vectors[element] == modulus-1)
				vectors[element] = 0;
			
			else
			{
				vectors[element] += 1;
				allTested = FALSE;
				break;
			}
		}
	}
	
	
	if (rows(F) == 2)
	{
		for (int i = 0; i < modulus; i += 1)
		{
			free(lineNumbers[i]);
			lineNumbers[i] = NULL;
		}
		free(lineNumbers);
		lineNumbers = NULL;
	}
	
	s = free_IntMatrixT(s);
	s_temp = free_IntMatrixT(s_temp);
	free(vectors);
	vectors = NULL;
	
	if (fclose(orbitFile) == EOF)
		return 0;
	
	return 1;
}


//UNFINISHED
CycleInfoTP floyd(IntMatrixTP F, IntMatrixTP s_0, int modulus)
/** Uses Floyd's Cycle Detection Algorithm to calculate the 
    transient length and cycle length of some given set up.
		s_0 is the starting vector. M is the update rule.
		
		Returns a CycleInfoTP that points to the data calculated
		about the cycle. */
{
	int stoppingTime = 0;
	IntMatrixTP Finv = inverse(F, modulus);
	
	CycleInfoTP info = malloc(sizeof(CycleInfoT));
	info->omega = 0;
	info->tau = -1;
	
	//If Finv exists, then no transient regions can exist
	if (Finv != NULL)
		info->tau = 0;
	
	//Initalising x and y to be s_0
	//We need to switch between the two versions to store data w/o overwriting
	IntMatrixTP x_1 = new_IntMatrixT(rows(s_0), cols(s_0));
	IntMatrixTP x_2 = new_IntMatrixT(rows(s_0), cols(s_0));
	IntMatrixTP y_1 = new_IntMatrixT(rows(s_0), cols(s_0)); 
	IntMatrixTP y_2 = new_IntMatrixT(rows(s_0), cols(s_0));
	
	copy_IntMatrixT(s_0, x_1);
	copy_IntMatrixT(s_0, y_1);
	
	//Iterate until x and y are the same
	do
	{
		mat_mul(F, x_1, x_2);
		modm(x_2, modulus);
		copy_IntMatrixT(x_2, x_1);
		
		//Iterate y twice
		mat_mul(F, y_1, y_2);
		modm(y_2, modulus);
		mat_mul(F, y_2, y_1);
		modm(y_1, modulus);
		
		stoppingTime += 1;
	}
	while (!compare_IntMatrixT(x_1, y_1));
	
	info->stoppingTime = stoppingTime;
	
	//Now we have a vector that's confirmed to be in a cycle
	//Now, we determine the cycle length
	
	//If St > 2*L, then w = St.
	if (stoppingTime > 2*rows(F))
		info->omega = stoppingTime;
	
	//If above inequality doesn't hold, calculate omega manually
	else
	{
		do
		{
			mat_mul(F, x_1, x_2);
			modm(x_2, modulus);
			copy_IntMatrixT(x_2, x_1);
			info->omega += 1;
		}
		while (!compare_IntMatrixT(x_1, y_1));
	}
	
	//I need a way to get tau from St and omega
	//If we know omega == 1, then tau is just St - 1 (unless tau = 1)
	
	//We do know that, given an L by L update matrix and a square-free
	// modulus, the maximum transient length is L.
	
	//pg 29 of LCA paper explains why prime powers are more complicated
	//They believe the bound on prime powered systems should be kL, where
	// k is the power of the prime
	
	//Freeing memory
	x_1 = free_IntMatrixT(x_1);
	x_2 = free_IntMatrixT(x_2);
	y_1 = free_IntMatrixT(y_1);
	y_2 = free_IntMatrixT(y_2);
	
	Finv = Finv != NULL ? free_IntMatrixT(Finv) : NULL;
	
	return info;
}