
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int GCD(int a, int b)
/** Returns the GCd of the given integers.
    This is a rudimentary implementation that'll
		be improved later. */
{
	int bound;
	int gcd = 1;
	
	if ((int)sqrt((double)a) > (int)sqrt((double)b))
		bound = (int)sqrt((double)a);
	else
		bound = (int)sqrt((double)b);
	
	for (int i = 1; i <= bound; i += 1)
		if ((a % i == 0) && (b % i == 0))
			gcd = i;
		
	return gcd;
}