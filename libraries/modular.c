
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


int square_root(int a, int modulus)
/** Returns the square root of a mod modulus.
    Returns -1 if the square root DNE. 
		
		According to Wikipedia, there's a more
		efficient way to write this, but I can't
		be bothered to look it up at the moment. */
{
	if (a == 0)
		return 0;
	
	else if (a == 1)
		return 1;
	
	else
	{
		for (int i = 2; i < modulus; i += 1)
			if ((i*i) % modulus == a)
				return i;
			
		return -1;
	}
}
