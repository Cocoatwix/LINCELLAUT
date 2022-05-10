
/*
A simple library for facilitating basic modular
arithmetic operations.

May 10, 2022
*/

int num_inverse(int a, int modulus)
/** Returns the inverse of a given number under
    the given modulus. Returns -1 if the inverse DNE.
		
		This is a basic implementation that I'd like to
		improve in the future, but I don't know how. */
{
	int inv = -1;
	
	for (int i = 1; i < modulus; i += 1)
		if ((a*i) % modulus == 1)
		{
			inv = i;
			break;
		}
		
	return inv;
}