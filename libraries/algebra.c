
/*
A simple library for manipulating
more typical algebra things, like
polynomials.

yiff day, 2022
*/

#include <stdlib.h>
#include <stdio.h>

#include "../headers/bigint.h"

typedef struct bigpoly
/** Stores a polynomial using BigIntT coefficients. 
    Holds coefficients in little endian form. */
{
	int size; //Holds how many terms are in the polynomial
	BigIntTP* coeffs; //Holds the coefficients of the polynomial
}
BigPolyT, *BigPolyTP;


BigPolyTP free_BigPolyT(BigPolyTP p)
/** Frees the memory of a given BigPolyT. Returns NULL. */
{
	for (int i = 0; i < p->size; i += 1)
		p->coeffs[i] = free_BigIntT(p->coeffs[i]);
	
	free(p->coeffs);
	p->coeffs = NULL;
	
	return NULL;
}


BigPolyTP new_BigPolyT(BigIntTP* const coefficients, int size)
/** Creates a new BigPolyT with specified coefficients and size,
    returns a pointer to it. */
{
	BigPolyTP p = malloc(sizeof(BigPolyT));
	p->size = size;
	
	p->coeffs = malloc(size*sizeof(BigIntTP));
	for (int i = 0; i < size; i += 1)
	{
		p->coeffs[i] = empty_BigIntT(1);
		copy_BigIntT(coefficients[i], p->coeffs[i]);
	}
	
	return p;
}


BigPolyTP constant_BigPolyT(BigIntTP const constant)
/** Creates a constant BigPolyTP, returns it. */
{
	BigPolyTP p = malloc(sizeof(BigPolyT));
	p->size = 1;
	p->coeffs = malloc(sizeof(BigIntTP));
	p->coeffs[0] = empty_BigIntT(1);
	copy_BigIntT(constant, p->coeffs[0]);
	
	return p;
}


BigPolyTP empty_BigPolyT()
/** Creates an empty BigPolyT struct and returns a pointer
    to it. */
{
	BigPolyTP p = malloc(sizeof(BigPolyT));
	p->size = 1;
	
	p->coeffs = malloc(sizeof(BigIntTP));
	p->coeffs[0] = empty_BigIntT(1);
	
	return p;
}


int reduce_BigPolyT(BigPolyTP p)
/** Ensures that the leading term of the BigPolyT
    is nonzero. 
		Returns 1 on success, 0 otherwise. */
{
	BigIntTP zero = empty_BigIntT(1);
	int newSize = 0;
	
	for (int i = p->size-1; i >= 0; i -= 1)
	{
		//If we find the leading nonzero term
		if (compare_BigIntT(p->coeffs[i], zero) != 0)
		{
			newSize = i+1;
			break;
		}
	}
	
	if (newSize != p->size)
	{
		//Freeing the coefficients we don't need anymore
		for (int i = p->size-1; i >= newSize; i -= 1)
			p->coeffs[i] = free_BigIntT(p->coeffs[i]);
		
		p->size = newSize;
		p->coeffs = realloc(p->coeffs, (p->size)*sizeof(BigIntTP));
	}
	
	zero = free_BigIntT(zero);
	return 1;
}


int degree(BigPolyTP const p)
/** Returns the degree of a BigPolyTP. */
{
	return p->size;
}


void printp(BigPolyTP const p)
/** Outputs a BigPolyTP to stdout. */
{
	BigIntTP zero = empty_BigIntT(1);
	
	for (int i = 0; i < p->size; i += 1)
	{
		//Only print if coefficient isn't zero
		if (compare_BigIntT(zero, p->coeffs[i]) != 0)
		{
			if (i == 0)
				printi(p->coeffs[0]);
			
			else if (i == 1)
			{
				printi(p->coeffs[1]);
				printf("λ");
			}
			
			else
			{
				printi(p->coeffs[i]);
				printf("λ^%d", i);
			}
			
			if (i != p->size-1)
				printf(" + ");
		}
	}
	
	zero = free_BigIntT(zero);
}


int add_BigPolyT(BigPolyTP const A, BigPolyTP const B, BigPolyTP sum)
/** Adds A and B together, stores the sum in sum.
    This function assumes sum has been initialised.
		Returns 1 on success, 0 otherwise. */
{
	int bigSize = (A->size > B->size) ? A->size : B->size;
	int smallSize = (A->size < B->size) ? A->size : B->size;
	
	BigIntTP temp = empty_BigIntT(1);
	
	//Making sure our sum is ready to take new terms
	if (sum->size > bigSize)
		for (int i = bigSize; i < sum->size; i += 1)
			sum->coeffs[i] = free_BigIntT(sum->coeffs[i]);
		
	if (sum->size != bigSize)
		sum->coeffs = realloc(sum->coeffs, bigSize*sizeof(BigIntTP));
	
	if (sum->size < bigSize)
		for (int i = sum->size; i < bigSize; i += 1)
			sum->coeffs[i] = empty_BigIntT(1);
			
	sum->size = bigSize;
	
	//Iterate over both polynomials, add them
	for (int a = 0; a < smallSize; a += 1)
	{
		add_BigIntT(A->coeffs[a], B->coeffs[a], temp);
		copy_BigIntT(temp, sum->coeffs[a]);
	}
	
	//Copy remaining terms 
	if (A->size > B->size)
		for (int b = smallSize; b < bigSize; b += 1)
			copy_BigIntT(A->coeffs[b], sum->coeffs[b]);
	
	else
		for (int b = smallSize; b < bigSize; b += 1)
			copy_BigIntT(B->coeffs[b], sum->coeffs[b]);
	
	temp = free_BigIntT(temp);
	
	reduce_BigPolyT(sum);
	return 1;
}


int multiply_BigPolyT(BigPolyTP const A, BigPolyTP const B, BigPolyTP product)
/** Multiples two polynomials together, stores the product in product.
    This function assumes product has been initialised. 
		Returns 1 on success, 0 otherwise. */
{
	BigIntTP temp = empty_BigIntT(1);
	BigIntTP temp2 = empty_BigIntT(1);
	
	//Preparing product to hold the new polynomial
	if (product->size != A->size + B->size - 1)
	{
		//Super inefficient, but prevents errors
		//Will change later
		for (int i = 0; i < product->size; i += 1)
			product->coeffs[i] = free_BigIntT(product->coeffs[i]);
		
		product->size = A->size + B->size - 1;
		product->coeffs = realloc(product->coeffs, (product->size)*sizeof(BigIntTP));
	}
	
	for (int i = 0; i < product->size; i += 1)
		product->coeffs[i] = empty_BigIntT(1);
	
	//Iterate through both polynomials, multiply them through
	for (int a = 0; a < A->size; a += 1)
	{
		for (int b = 0; b < B->size; b += 1)
		{
			multiply_BigIntT(A->coeffs[a], B->coeffs[b], temp);
			add_BigIntT(temp, product->coeffs[a+b], temp2);
			copy_BigIntT(temp2, product->coeffs[a+b]);
		}
	}
	
	temp = free_BigIntT(temp);
	temp2 = free_BigIntT(temp2);
	
	reduce_BigPolyT(product);
	
	return 1;
}


int mod_BigPolyT(BigPolyTP const A, BigIntTP const mod, BigPolyTP residue)
/** Performs a mod operation on each element of a polynomial,
    stores the result in residue. This function assumes all arguments
		have been properly initialised. 
		Returns 1 on success, 0 otherwise. */
{
	//Preparing residue to hold what it needs to hold
	if (residue->size > A->size)
		for (int i = A->size; i < residue->size; i += 1)
			residue->coeffs[i] = free_BigIntT(residue->coeffs[i]);
		
	if (residue->size != A->size)
		residue->coeffs = realloc(residue->coeffs, (A->size)*sizeof(BigIntTP));
	
	if (residue->size < A->size)
		for (int i = residue->size; i < A->size; i += 1)
			residue->coeffs[i] = empty_BigIntT(1);
	
	residue->size = A->size;
	
	//Now, perform the modulo operation
	for (int i = 0; i < A->size; i += 1)
		mod_BigIntT(A->coeffs[i], mod, residue->coeffs[i]);
	
	reduce_BigPolyT(residue);
	
	return 1;
}
