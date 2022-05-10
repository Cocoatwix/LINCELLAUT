
#include <stdlib.h>
#include <stdio.h>

#include "../headers/linalg.h"  //Letting us use matrices
#include "../headers/factors.h" //Letting us check for square-free moduli
#include "../headers/modular.h" //For using modular arithmetic stuff

#define FREE(v) free(v); v = NULL

//Might need to write a "free" function later on
//May want to store the vector that's in the cycle as well
typedef struct cycleinfo
/** Struct to hold any important info about a structure. */
{
	int omega; //Cycle length
	int tau;   //Transient length
	int stoppingTime;
	
	IntMatrixTP inCycle; //Stores a vector that's in the cycle
}
CycleInfoT, *CycleInfoTP;


CycleInfoTP free_CycleInfoT(CycleInfoTP c)
/** Frees memory used by the given CycleInfo struct.
    Returns NULL. */
{
	free_IntMatrixT(c->inCycle);
	return NULL;
}


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
	
	printf("Representative:\n");
	printm(c->inCycle, TRUE);
}


int omega(CycleInfoTP c)
/** Returns the cycle length of a given CycleInfoT object. */
{
	return c->omega;
}


int tau(CycleInfoTP c)
/** Returns the transient length of a given CycleInfoT object. */
{
	return c->tau;
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
 
 
void visit_points(IntMatrixTP F, int modulus, int iterations)
/** Prints out the different "points" that F travels to.
    For example, the matrix [[2 0][0 2]] goes to the identity
		after 3 iterations for mod 7. The identity could be a point.
		
		If iterations is set to -1, this function iterates until
		F reaches the identity. */
{
	//Making sure F is square
	if (rows(F) == cols(F))
	{
		IntMatrixTP iteratedF = identity_IntMatrixT(rows(F));
		IntMatrixTP I         = identity_IntMatrixT(rows(F));
		IntMatrixTP copyF     = new_IntMatrixT(rows(F), rows(F));
		
		int scale; //Used for determining what kind of scale matrix we have
		int elem;  //Used to hold an element from the matrix
		
		printf("0: I\n");
		
		for (int i = 0; i < iterations; i += 1)
		{
			//Perform an iteration
			mat_mul(F, iteratedF, copyF);
			modm(copyF, modulus);
			copy_IntMatrixT(copyF, iteratedF);
			
			//If we know our matrix is a scale matrix
			//This check will need to be more complicated in the future
			if (is_diagonal(iteratedF))
			{
				elem = element(iteratedF, 0, 0);
				printf("%d: %dI or %dI; ", i+1, elem, elem - modulus);
				scale = num_inverse(elem, modulus);
				if (scale > 0)
					printf("Acts as %d^-1 or %d^-1", scale, scale-modulus);
				printf("\n");
			}
			
			//If we get back to the identity, stop
			if (compare_IntMatrixT(I, iteratedF))
				break;
		}
		
		iteratedF = free_IntMatrixT(iteratedF);
		copyF     = free_IntMatrixT(copyF);
	}
}

//UNFINISHED
CycleInfoTP floyd(IntMatrixTP F, IntMatrixTP s_0, int modulus)
/** Uses Floyd's Cycle Detection Algorithm to calculate the 
    transient length and cycle length of some given set up.
		s_0 is the starting vector. F is the update rule.
		
		Returns a CycleInfoTP that points to the data calculated
		about the cycle. Returns NULL on error. */
{
	if (rows(F) != cols(F))
		return NULL;
	
	int stoppingTime = 0;
	IntMatrixTP Finv = inverse(F, modulus);
	
	CycleInfoTP info = malloc(sizeof(CycleInfoT));
	info->omega = 0;
	info->tau = -1;
	info->inCycle = new_IntMatrixT(rows(F), 1);
	
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
	
	//If, after iterating, we end up at the same vector,
	// that means we started on a stationary point
	if ((compare_IntMatrixT(s_0, x_1)) && (stoppingTime == 1))
		info->tau = 0;
	
	info->stoppingTime = stoppingTime;
	copy_IntMatrixT(x_1, info->inCycle);
	
	//Now we have a vector that's confirmed to be in a cycle
	//Now, we determine the cycle length
	
	if (info->tau == 0)
		info->omega = stoppingTime;
	
	//If St > 2*L, and we have a square free modulus, then w = St.
	else if ((stoppingTime > 2*rows(F)) && is_square_free(modulus))
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
	if ((info->omega == 1) && (info->tau != 0))
		info->tau = stoppingTime;
	
	//Maybe we can do casework here?
	//For instance, if we know omega == 2, we can
	// go through each possible case for the Stopping time
	// (even or odd) and make a conclusion about tau that way?
	// This may then extend to omega == 3, 4, etc.
	
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


int write_iteration(const char* fileName, IntMatrixTP F, int modulus)
/** Writes out a .iteration file containing every vector's result when
    multiplied by F mod modulus. Returns 1 on success, 0 otherwise. */
{
	//F must be a square matrix
	if (rows(F) != cols(F))
		return 0;
	
	FILE* iterFile = fopen(fileName, "w");
	
	if (iterFile == NULL)
		return 0;
	
	IntMatrixTP x_1 = new_IntMatrixT(rows(F), 1); //Current vector
	IntMatrixTP x_2 = new_IntMatrixT(rows(F), 1); //Iterated vector
	
	//This will hold all of our starting values for our vectors
	int* vectors = calloc(rows(F), sizeof(int));
	
	bool allTested = FALSE;
	
	
	//Loop until we've found the orbits of all vectors
	while (!allTested)
	{
		//Iterate x_1 once, write to file
		set_column(x_1, vectors);
		
		mat_mul(F, x_1, x_2);
		modm(x_2, modulus);
		for (int i = 0; i < rows(F); i += 1)
			fprintf(iterFile, "%d ", element(x_2, i, 0));
		
		fprintf(iterFile, "\n");
		
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
	
	x_1 = free_IntMatrixT(x_1);
	x_2 = free_IntMatrixT(x_2);
	
	FREE(vectors);
	
	if (fclose(iterFile) == EOF)
		return 0;

	return 1;
}
