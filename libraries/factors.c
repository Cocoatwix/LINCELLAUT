
#include <stdlib.h>
#include <stdio.h>
//#include <math.h>

#include "../headers/helper.h" //bool

#include "../headers/bigint.h"
#include "../headers/algebra.h"

#include "../headers/modular.h" //big_num_inverse()

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


BigIntTP big_gcd(const BigIntTP A, const BigIntTP B)
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
		greater = empty_BigIntT(1);
		copy_BigIntT(B, greater);
		return greater; //return B by value
	}
	
	else if (compare_BigIntT(B, zero) == 0)
	{
		zero = free_BigIntT(zero);
		greater = empty_BigIntT(1);
		copy_BigIntT(A, greater);
		return greater; //return A by value
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


int poly_gcd(const BigPolyTP p, const BigPolyTP q, BigPolyTP gcd, const BigIntTP mod, BigPolyTP s, BigPolyTP t)
/** Calculates the GCD of two given polynomials, stores the
    result in the third BigPolyT given. The fourth and fifth 
		BigPolyTs will store the polynomials s and t guaranteed by Bezout's 
		identity, meaning gcd(p, q) = sp + tq. The BigIntT is the 
		modulus to use.
		Returns 1 on success, 0 otherwise. */
{
	bool hasZero = FALSE; //Says whether either of p or q is a zero
	
	BigPolyTP A, B, R, tempQ, zero;
	BigPolyTP cA, cB; //References to s and t
	
	/* Currently, if cA or cB is NULL, the function prevents the terms
	 *  from Bezout's identity to be copied (to avoid crashing).
	 *  However, the function still calculates these terms all the same.
	 *  In the future, it would be nice to be able to tell the function
	 *  not to do these calculations if it sees cA or cB is NULL to 
	 *  increase speed. For now, I'm too lazy :)
	 */
	
	//These are for simplifying the GCD
	BigIntTP  leadingTermInv = NULL;
	BigPolyTP invPoly        = NULL;
	
	int oneArr[1]    = {1};
	BigIntTP  one    = NULL;
	BigIntTP  negOne = NULL;
	BigPolyTP onePoly;
	BigPolyTP negOnePoly;
	
	//[s or t, A or B, t or s, B or A]
	BigPolyTP*  bezoutArray;           //Holds coefficients and stuff for finding s and t
	
	//[[A, B, Q, R], [A, B, Q, R], ...]
	BigPolyTP** divisionRecord = NULL; //Holds all the polynomials needed to find s and t
	int divisionRecordSize = 0;
	
	A     = empty_BigPolyT();
	B     = empty_BigPolyT();
	R     = empty_BigPolyT();
	tempQ = empty_BigPolyT();
	zero  = empty_BigPolyT();
	
	one = new_BigIntT(oneArr, 1);
	onePoly = constant_BigPolyT(one);
	
	negOne = empty_BigIntT(1);
	subtract_BigIntT(mod, one, negOne);
	negOnePoly = constant_BigPolyT(negOne);
	
	//Checking to see which polynomial given is greater
	if (compare_BigPolyT(p, q) < 0) //q is bigger
	{
		copy_BigPolyT(q, A);
		copy_BigPolyT(p, B);
		
		cA = t;
		cB = s;
	}
	else
	{
		copy_BigPolyT(p, A);
		copy_BigPolyT(q, B);
		
		cA = s;
		cB = t;
	}
	
	reduce_BigPolyT(A);
	reduce_BigPolyT(B);
	
	if (compare_BigPolyT(zero, B) == 0)
	{
		mod_BigPolyT(A, mod, gcd);
		if (cA != NULL)
			copy_BigPolyT(onePoly, cA);
		if (cB != NULL)
			copy_BigPolyT(zero, cB);
		
		hasZero = TRUE;
	}
	
	else
	{
		//Euclidian algorithm
		do
		{
			divisionRecordSize += 1;
			divisionRecord = realloc(divisionRecord, divisionRecordSize*sizeof(BigPolyTP*));
			divisionRecord[divisionRecordSize-1] = malloc(4*sizeof(BigPolyTP));
			for (int i = 0; i < 4; i += 1)
				divisionRecord[divisionRecordSize-1][i] = empty_BigPolyT();
			
			divide_BigPolyT(A, B, tempQ, R, mod);
			
			copy_BigPolyT(A, divisionRecord[divisionRecordSize-1][0]);
			copy_BigPolyT(B, divisionRecord[divisionRecordSize-1][1]);
			copy_BigPolyT(tempQ, divisionRecord[divisionRecordSize-1][2]);
			copy_BigPolyT(R, divisionRecord[divisionRecordSize-1][3]);
			
			copy_BigPolyT(B, A);
			copy_BigPolyT(R, B);
		}
		while (compare_BigPolyT(zero, R) != 0);
		
		//Now, we should simplify A so that its leading term has a coefficient of 1
		reduce_BigPolyT(A);
		leadingTermInv = empty_BigIntT(1);
		big_num_inverse(leading_term(A), mod, leadingTermInv);
		
		invPoly = constant_BigPolyT(leadingTermInv);
		multiply_BigPolyT(A, invPoly, tempQ);
		mod_BigPolyT(tempQ, mod, A);

		bezoutArray = malloc(4*sizeof(BigPolyTP));
		for (int i = 0; i < 4; i += 1)
			bezoutArray[i] = empty_BigPolyT();
		
		//Do some casework for getting s and t correctly
		
		//If B | A, meaning gcd(A, B) = B
		if (divisionRecordSize == 1)
		{
			if (cA != NULL)
				copy_BigPolyT(divisionRecord[0][3], cA);
			if (cB != NULL)
				copy_BigPolyT(onePoly, cB);
		}
		
		else
		{
			//Prepare bezoutArray for extended Euclidian algorithm
			if (divisionRecordSize >= 2)
			{
				copy_BigPolyT(onePoly, bezoutArray[0]);
				copy_BigPolyT(divisionRecord[divisionRecordSize-2][0], bezoutArray[1]);
				
				multiply_BigPolyT(divisionRecord[divisionRecordSize-2][2], negOnePoly, tempQ);
				mod_BigPolyT(tempQ, mod, bezoutArray[2]);
				copy_BigPolyT(divisionRecord[divisionRecordSize-2][1], bezoutArray[3]);
			}
			
			//Extended Euclidian algorithm, if needed
			if (divisionRecordSize > 2)
			{
				for (int step = 3; step <= divisionRecordSize; step += 1)
				{
					//Change right number
					if (step % 2 == 1)
					{
						//bezoutArray[3] is being used as a temp variable in this block
						multiply_BigPolyT(bezoutArray[2], divisionRecord[divisionRecordSize-step][2], tempQ);
						multiply_BigPolyT(tempQ, negOnePoly, bezoutArray[3]);
						add_BigPolyT(bezoutArray[3], bezoutArray[0], tempQ);
						mod_BigPolyT(tempQ, mod, bezoutArray[0]);
						
						copy_BigPolyT(divisionRecord[divisionRecordSize-step][0], bezoutArray[3]);
					}
					
					//Change left number
					else
					{
						//bezoutArray[1] is being used as a temp variable in this block
						multiply_BigPolyT(bezoutArray[0], divisionRecord[divisionRecordSize-step][2], tempQ);
						multiply_BigPolyT(tempQ, negOnePoly, bezoutArray[1]);
						add_BigPolyT(bezoutArray[1], bezoutArray[2], tempQ);
						mod_BigPolyT(tempQ, mod, bezoutArray[2]);
						
						copy_BigPolyT(divisionRecord[divisionRecordSize-step][0], bezoutArray[1]);
					}
				}
			}
			
			//Find s and t, copy them to cA and cB
			if (compare_BigPolyT(divisionRecord[0][0], bezoutArray[1]) == 0)
			{
				//Correcting s and t to match with our MONIC gcd
				if (cA != NULL)
				{
					multiply_BigPolyT(bezoutArray[0], invPoly, tempQ);
					mod_BigPolyT(tempQ, mod, cA);
				}
				if (cB != NULL)
				{
					multiply_BigPolyT(bezoutArray[2], invPoly, tempQ);
					mod_BigPolyT(tempQ, mod, cB);
				}
			}
			
			else
			{
				//Correcting s and t to match with our MONIC gcd
				if (cB != NULL)
				{
					multiply_BigPolyT(bezoutArray[0], invPoly, tempQ);
					mod_BigPolyT(tempQ, mod, cB);
				}
				if (cA != NULL)
				{
					multiply_BigPolyT(bezoutArray[2], invPoly, tempQ);
					mod_BigPolyT(tempQ, mod, cA);
				}
			}
		}
		
		//Now, s and t should be copied into their correct variables
		//We also need to copy the GCD itself
		copy_BigPolyT(A, gcd);
	}
	
	for (int i = 0; i < divisionRecordSize; i += 1)
	{
		for (int j = 0; j < 4; j += 1)
			divisionRecord[i][j] = free_BigPolyT(divisionRecord[i][j]);
		
		free(divisionRecord[i]);
		divisionRecord[i] = NULL;
	}
	free(divisionRecord);
	divisionRecord = NULL;
	
	if (!hasZero)
	{
		for (int i = 0; i < 4; i += 1)
			bezoutArray[i] = free_BigPolyT(bezoutArray[i]);
		free(bezoutArray);
		bezoutArray = NULL;
	}
	
	A = free_BigPolyT(A);
	B = free_BigPolyT(B);
	R = free_BigPolyT(R);
	
	tempQ      = free_BigPolyT(tempQ);
	zero       = free_BigPolyT(zero);
	onePoly    = free_BigPolyT(onePoly);
	negOnePoly = free_BigPolyT(negOnePoly);
	
	one            = free_BigIntT(one);
	negOne         = free_BigIntT(negOne);
	leadingTermInv = free_BigIntT(leadingTermInv);
	invPoly        = free_BigPolyT(invPoly);
	
	return 1;
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


int big_lcm(const BigIntTP a, const BigIntTP b, BigIntTP L)
/** Computes the LCM of the first two BigIntTs, stores the result in 
    the third one. This function assumes all BigIntTs have been 
	  initialised. Returns 1 on success, 0 otherwise. */
{
	if ((is_zero(a)) || (is_zero(b)))
		return 0;
	
	BigIntTP* primeFactorsOfa = prime_factors_of_BigIntT(a);
	BigIntTP* primeFactorsOfb = prime_factors_of_BigIntT(b);
	
	BigIntTP temp = empty_BigIntT(1);

	clear_BigIntT(L);
	set_bunch(L, 0, 1);

	int bPos = 1;
	//Go through both arrays of prime factors,
	// multiply the factors together as needed to
	// compute the LCM.
	for (int aPos = 1; aPos <= extract_bunch(primeFactorsOfa[0], 0); aPos += 1)
	{
		multiply_BigIntT(primeFactorsOfa[aPos], L, temp);
		copy_BigIntT(temp, L);
		
		//Make sure we multiply by all prime factors in b that are smaller
		// than the one we just multipled by in a
		while ((bPos <= extract_bunch(primeFactorsOfb[0], 0)) && 
		       (compare_BigIntT(primeFactorsOfa[aPos], primeFactorsOfb[bPos]) > 0))
		{
			multiply_BigIntT(primeFactorsOfb[bPos], L, temp);
			copy_BigIntT(temp, L);
			
			primeFactorsOfb[bPos] = free_BigIntT(primeFactorsOfb[bPos]);
			bPos += 1;
		}
		
		//If there's a duplicate prime factor between a and b, make sure
		// to only multiply by it once
		if ((bPos <= extract_bunch(primeFactorsOfb[0], 0)) &&
		    (compare_BigIntT(primeFactorsOfa[aPos], primeFactorsOfb[bPos]) == 0))
		{
			primeFactorsOfb[bPos] = free_BigIntT(primeFactorsOfb[bPos]);
			bPos += 1;
		}
			
		primeFactorsOfa[aPos] = free_BigIntT(primeFactorsOfa[aPos]);
	}
	
	//Now, multiply by any remaining prime factors in b
	// that we haven't got yet
	for (int i = bPos; i <= extract_bunch(primeFactorsOfb[0], 0); i += 1)
	{
		multiply_BigIntT(primeFactorsOfb[i], L, temp);
		copy_BigIntT(temp, L);
		
		primeFactorsOfb[i] = free_BigIntT(primeFactorsOfb[i]);
	}
	
	primeFactorsOfa[0] = free_BigIntT(primeFactorsOfa[0]);
	free(primeFactorsOfa);

	primeFactorsOfb[0] = free_BigIntT(primeFactorsOfb[0]);
	free(primeFactorsOfb);
	
	temp = free_BigIntT(temp);
	
	return 1;
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


int num_divisors(int n)
/** Returns the number of divisors the number has. */
{
	int divisors = 0;
	for (int count = 1; n/count >= count; count += 1)
	{
		if (n % count == 0)
		{
			divisors += 1;
			if (n/count != count)
				divisors += 1;
		}
	}
	
	return divisors;
}