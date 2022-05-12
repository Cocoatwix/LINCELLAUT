
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