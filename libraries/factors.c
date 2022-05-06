
#include <stdlib.h>
//#include <math.h>

int GCD(int a, int b)
/** Returns the GCD of the given integers.
    This is a rudimentary implementation that'll
		be improved later (maybe). */
{
	int bound;
	int gcd = 1;
	
	if (a > b)
		bound = a;
	else
		bound = b;
	
	for (int i = 2; i <= bound; i += 1)
		if ((a % i == 0) && (b % i == 0))
			gcd = i;
		
	return gcd;
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