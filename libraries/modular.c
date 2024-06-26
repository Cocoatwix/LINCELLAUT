
/*
A simple library for facilitating basic modular
arithmetic operations.

May 10, 2022
*/

#include "../headers/bigint.h" //For big_num_inverse()

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


void big_num_inverse(const BigIntTP findInverse, const BigIntTP modulus, BigIntTP result)
/** Same as num_inverse, but for BigIntTs. Stores result in result (must be initialised).
    Stores 0 on error or if inverse doesn't exist. */
{
	int oneArr[] = {1};
	BigIntTP zero  = empty_BigIntT(1);
	BigIntTP one   = new_BigIntT(oneArr, 1);
	BigIntTP temp  = empty_BigIntT(1);
	BigIntTP temp2 = empty_BigIntT(1);
	BigIntTP temp3 = empty_BigIntT(1);
	
	copy_BigIntT(one, temp);
	copy_BigIntT(zero, result);
	
	//Search for inverse
	//Apparently the extended Euclidian algorithm can be used here
	// Will probably implement that later
	if (compare_BigIntT(zero, findInverse) != 0)
	{
		while (compare_BigIntT(temp, modulus) < 0)
		{
			multiply_BigIntT(temp, findInverse, temp2);
			mod_BigIntT(temp2, modulus, temp3);
			
			//Check to see if this is the correct inverse
			if (compare_BigIntT(one, temp3) == 0)
			{
				copy_BigIntT(temp, result);
				break;
			}
			
			//Increment temp
			add_BigIntT(one, temp, temp3);
			copy_BigIntT(temp3, temp);
		}
	}
	
	zero  = free_BigIntT(zero);
	one   = free_BigIntT(one);
	temp  = free_BigIntT(temp);
	temp2 = free_BigIntT(temp2);
	temp3 = free_BigIntT(temp3);
}


//BigIntTP base_prime(const BigIntTP pp)
/** Given a prime-power, this function returns
    what the original prime was. If a composite is
	  given, this returns the smallest prime factor.
	  Returns NULL on error. */ /*
{
	int numArr[1] = {1};
	BigIntTP one = new_BigIntT(numArr, 1);
	numArr[0] = 2;
	BigIntTP two = new_BigIntT(numArr, 1);
	BigIntTP prime  = empty_BigIntT(1);
	BigIntTP result = empty_BigIntT(1);
	BigIntTP temp   = empty_BigIntT(1);
	
	copy_BigIntT(one, prime);
	//Loop until we find a prime divisor of pp
	do
	{
		//Increment prime
		//I could make this faster by only checking numbers 
		// one more or less than a multiple of 6, but eh
		if ((compare_BigIntT(prime, one) == 0) ||
		    (compare_BigIntT(prime, two) == 0))
			add_BigIntT(prime, one, temp);
		else
			add_BigIntT(prime, two, temp);
		
		copy_BigIntT(temp, prime);
		mod_BigIntT(pp, prime, result);
	}
	while (compare_BigIntT(result, zero) != 0);
	
	one = free_BigIntT(one);
	two = free_BigIntT(two);
	temp = free_BigIntT(temp);
	result = free_BigIntT(result);
	
	return prime;
}
*/


int num_root(int a, int modulus)
/** Returns the smallest nonzero number that, when multiplied by a,
    gives 0 when reduced by modulus. 
		Returns 0 if no nonzero number exists. */
{
	int root = 0;
	
	for (int i = 1; i < modulus; i += 1)
		if ((a*i) % modulus == 0)
		{
			root = i;
			break;
		}
	
	return root;
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
