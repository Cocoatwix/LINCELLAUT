
#include <stdlib.h>
#include <stdio.h>
//#include <math.h>

int GCD(int a, int b)
/** Returns the GCD of the given integers.
    If both integers given are zero, returns -1. */
{
	if ((a == 0) && (b == 0))
		return -1;
	
	else if (a == 0)
		return b;
	
	else if (b == 0)
		return a;
	
	int greater = a;
	int lesser  = b;
	int r;
	
	if (b > a)
	{
		greater = b;
		lesser  = a;
	}
	
	r = greater % lesser; //Holds the current remainder
	
	while (r != 0)
	{
		greater = lesser;
		lesser = r;
		r = greater % lesser;
	}
	
	return lesser;
}


int LCM(int a, int b)
/** Returns the least common multiple of the two given integers. */
{
	if ((a == 0) || (b == 0))
		return 0;
	
	int newA = a;
	int newB = b;
	
	while (1)
	{
		if (newA < newB)
		{
			newA += a;
		}
		
		else if (newB < newA)
		{
			newB += b;
		}
		
		else
			return newA;
	}
}


int* prime_factors(int n)
/** Returns a pointer of n's prime factors.
    It's the caller's job to free the pointer
		when they're done with it. 
		
		The first number in the pointer says
		how many prime factors were found. */
{
	int* factors = malloc(2*sizeof(int));
	int factorsCount = 0;
	
	//Yes, there are more efficient ways to do this.
	//Leave me alone.
	for (int i = 2; i <= n/i; i += 1)
	{
		//If we found a factor
		if (n % i == 0)
		{
			//printf("New upper bound: %d\n", matCycle/i);
			factors[factorsCount+1] = i;
			factorsCount += 1;
			n /= i;
			factors = realloc(factors, (factorsCount+2)*sizeof(int));
			i = 1; //Start looking for prime factors at 2 again
		}
	}
	
	//Adding the last prime factor which will always be left out
	factors[factorsCount+1] = n;
	factorsCount += 1;
	factors[0] = factorsCount;
	
	return factors;
}


int is_square_free(int n)
/** Returns 1 if n is a square-free number (each factor
    only appears only once in n's factorisation), 0
		otherwise. This is a rudimentary implementation that'll
		be improved later (maybe). */
{
	for (int i = 2; i < n; i += 1)
		if ((n % i == 0) && ((n/i) % i == 0))
			return 0;
		
	return 1;
}