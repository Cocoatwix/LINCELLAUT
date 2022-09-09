
#ifndef ALGEBRA_H
#define ALGEBRA_H

#include "bigint.h"

/** Holds a polynomial datatype using BigIntT coefficients. */
typedef struct bigpoly *BigPolyTP;

/** Frees the memory used by a BigPolyT. Returns NULL. */
BigPolyTP free_BigPolyT(BigPolyTP);

/** Frees an array of BigPolyT factors. */
BigPolyTP* free_BigPolyT_factors(BigPolyTP*);

/** Creates a new BigPolyT with given coefficients and size,
    returns a pointer to it. */
BigPolyTP new_BigPolyT(BigIntTP* const, int);

/** Creates a constant polynomial with the given BigIntTP
    as the constant term. Returns a pointer to the new
		constant. */
BigPolyTP constant_BigPolyT(BigIntTP const);

/** Creates an empty polynomial, returns a pointer to it. */
BigPolyTP empty_BigPolyT();

/** Returns the degree of the given polynomial. */
int degree(BigPolyTP const);

/** Returns the constant term of a BigPolyTP. */
BigIntTP constant(BigPolyTP const);

/** Returns an array of the given BigPolyT's coefficients,
    sorted in ascending order of the terms' exponents.
	Returns NULL on error. */
BigIntTP* extract_coefficients(BigPolyTP const);

/** Copies a polynomial into another. This function assumes
    all arguments have been properly initialised.
		Returns 1 on success, 0 otherwise. */
int copy_BigPolyT(BigPolyTP const, BigPolyTP);

/** Outputs a BigPolyT to stdout (the console). */
void printp(BigPolyTP const);

/** Prints a factored BigPolyT. */
void printpf(BigPolyTP* const);

/** Same as printp and printpf, but it prints to a file stream. */
void fprintp(FILE*, BigPolyTP const);
void fprintpf(FILE*, BigPolyTP* const);

/** Adds the first two polynomials, stores the sum in the third.
    This function assumes all polynomials passed have been
		initialised.
    Returns 1 on success, 0 otherwise. */
int add_BigPolyT(BigPolyTP const, BigPolyTP const, BigPolyTP);

/** Multiples the first two polynomials given, stores the
    product in the third BigPolyTP. This function assumes all 
		polynomials passed have been initialised.
		Returns 1 on success, 0 otherwise. */
int multiply_BigPolyT(BigPolyTP const, BigPolyTP const, BigPolyTP);

/** Reduces a given polynomial by a modulus, stores the result in
    another passed polynomial. This function assumes all arguments
		have been properly initialised.
		Returns 1 on success, 0 oterwise. */
int mod_BigPolyT(BigPolyTP const, BigIntTP const, BigPolyTP);

/** Factors the given BigPolyT under the given modulus.
    Returns a pointer to all the polynomial's factors.
		The first element is a BigPolyT constant which says
		how many factors are in the pointer. */
BigPolyTP* factor_BigPolyT(BigPolyTP const, BigIntTP);

/*
int find_factors(BigIntTP const,
                 BigIntTP const,
								 BigIntTP const,
								 BigIntTP const,
			  				 BigIntTP*);
*/
								 
#endif //ALGEBRA_H