
#include <stdlib.h>
#include <stdio.h>
//#include <math.h>

#include "../headers/bigint.h"
#include "../headers/algebra.h"

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


BigIntTP big_gcd(BigIntTP const A, BigIntTP const B)
/** Same as GCD(), but for BigIntTs. */
{
	BigIntTP zero = empty_BigIntT(1);
	BigIntTP greater;
	BigIntTP lesser;
	BigIntTP r;
	
	if ((compare_BigIntT(A, zero) == 0) && (compare_BigIntT(B, zero) == 0))
	{
		zero = free_BigIntT(zero);
		return NULL;
	}
	
	else if (compare_BigIntT(A, zero) == 0)
	{
		zero = free_BigIntT(zero);
		return B;
	}
	
	else if (compare_BigIntT(B, zero) == 0)
	{
		zero = free_BigIntT(zero);
		return A;
	}
	
	greater = empty_BigIntT(1);
	lesser = empty_BigIntT(1);
	r = empty_BigIntT(1);
	
	if (compare_BigIntT(B, A) > 0)
	{
		copy_BigIntT(B, greater);
		copy_BigIntT(A, lesser);
	}
	else
	{
		copy_BigIntT(A, greater);
		copy_BigIntT(B, lesser);
	}
	
	mod_BigIntT(greater, lesser, r); //r holds the current remainder
	
	while (compare_BigIntT(r, zero) != 0)
	{
		copy_BigIntT(lesser, greater);
		copy_BigIntT(r, lesser);
		mod_BigIntT(greater, lesser, r);
	}
	
	zero    = free_BigIntT(zero);
	greater = free_BigIntT(greater);
	r       = free_BigIntT(r);

	return lesser;
}


int poly_gcd(BigPolyTP const p, BigPolyTP const q, BigPolyTP gcd, BigIntTP const mod, BigPolyTP s, BigPolyTP t)
/** Calculates the GCD of two given polynomials, stores the
    result in the third BigPolyT given. The fourth and fifth 
		BigPolyTs will store the polynomials s and t guaranteed by Bezout's 
		identity, meaning gcd(a, b) = sa + tb. The BigIntT is the 
		modulus to use.
		Returns 1 on success, 0 otherwise. */
{
	BigPolyTP A, B, R, tempQ, zero;
	
	//These are for simplifying the GCD
	BigIntTP  leadingTermInv;
	BigPolyTP invPoly;
	
	A = empty_BigPolyT();
	B = empty_BigPolyT();
	R = empty_BigPolyT();
	tempQ = empty_BigPolyT();
	zero = empty_BigPolyT();
	
	//Checking to see which polynomial given is greater
	if (compare_BigPolyT(p, q) < 0) //q is bigger
	{
		copy_BigPolyT(q, A);
		copy_BigPolyT(p, B);
	}
	else
	{
		copy_BigPolyT(p, A);
		copy_BigPolyT(q, B);
	}
	
	reduce_BigPolyT(A);
	reduce_BigPolyT(B);
	
	//Euclidian algorithm
	do
	{
		divide_BigPolyT(A, B, tempQ, R, mod);
		copy_BigPolyT(B, A);
		copy_BigPolyT(R, B);
		
		/*printf("A = ");
		printp(A);
		printf("\nB = ");
		printp(B);
		printf("\nR = ");
		printp(R);
		printf("\n\n"); */
	}
	while (compare_BigPolyT(zero, R) != 0);
	
	//Now, we should simplify A so that its leading term has a coefficient of 1
	reduce_BigPolyT(A);
	leadingTermInv = empty_BigIntT(1);
	big_num_inverse(leading_term(A), mod, leadingTermInv);
	
	invPoly = constant_BigPolyT(leadingTermInv);
	multiply_BigPolyT(A, invPoly, tempQ);
	mod_BigPolyT(tempQ, mod, A);
	
	printp(A);
	
	//Extended Euclidian algorithm (finding polynomials guaranteed by Bezout's identity)
	
	A = free_BigPolyT(A);
	B = free_BigPolyT(B);
	R = free_BigPolyT(R);
	
	tempQ = free_BigPolyT(tempQ);
	zero = free_BigPolyT(zero);
	
	leadingTermInv = free_BigIntT(leadingTermInv);
	invPoly        = free_BigPolyT(invPoly);
	
	return 0;
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