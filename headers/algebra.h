
#ifndef ALGEBRA_H
#define ALGEBRA_H

#include "bigint.h"

/** Holds a polynomial datatype using BigIntT coefficients. */
typedef struct bigpoly *BigPolyTP;

/** Responsible for facilitating multidimensional BigIntT arrays
    of varying dimension. */
//typedef struct bigintdirection *BigIntDirectionTP;

/** Holds a representation of a multivariate field extension
    polynomial expression. */
typedef struct multivarext *MultiVarExtTP;

/** Frees the memory used by a BigPolyT. Returns NULL. */
BigPolyTP free_BigPolyT(BigPolyTP);

/** Frees an array of BigPolyT factors. */
BigPolyTP* free_BigPolyT_factors(BigPolyTP*);

/** Frees the memory used by a MultiVarExtT. Returns NULL. */
MultiVarExtTP free_MultiVarExtT(MultiVarExtTP);

/** Creates a new BigPolyT with given coefficients and size,
    returns a pointer to it. */
BigPolyTP new_BigPolyT(BigIntTP* const, int);

/** Creates a constant polynomial with the given BigIntTP
    as the constant term. Returns a pointer to the new
		constant. */
BigPolyTP constant_BigPolyT(BigIntTP const);

/** Creates an empty polynomial, returns a pointer to it. */
BigPolyTP empty_BigPolyT();

/** Ensures that the leading term of the BigPolyT
    is nonzero. Returns 1 on success, 0 otherwise. */
int reduce_BigPolyT(BigPolyTP);

/** Allocates space for a new MultiVarExtT, returns a
    pointer to it. Returns NULL on error. */
MultiVarExtTP new_MultiVarExtT(int);

/** Sets the modulus value for the given MultiVarExtT.
    Returns 1 on success, 0 otherwise. */
int set_MultiVarExtT_mod(MultiVarExtTP, BigIntTP const);

/** Adds a new extension to a MultiVarExtT.
    Returns 1 on success, 0 otherwise. */
//                extension      minPoly          size name
int add_extension(MultiVarExtTP, BigIntTP* const, int, char* const);

/** Sets a particular value for the given MultiVarExtT's coefficient.
    The MultiVarExtT must be fully set before this function can be used.
	Returns 1 on success, 0 otherwise. */
int set_MultiVarExtT_coefficient(MultiVarExtTP, int* const, BigIntTP const);

/** Uses the extension definitions to reduce the given
    MultiVarExtT. Essentially, it'll zero any terms
	it can. 
	This function can only be used when all extensions
	have been set for the MultiVarExtT.
	Returns 1 on success, 0 otherwise. */
int reduce_MultiVarExtT(MultiVarExtTP);

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

/** Compares two BigPolyTPs to see if they're equal (have the
    same coefficients). Returns 0 if they are, 1 otherwise. */ 
int compare_BigPolyT(BigPolyTP const, BigPolyTP const);

/** Outputs a BigPolyT to stdout (the console). */
void printp(BigPolyTP const);

/** Prints a factored BigPolyT. */
void printpf(BigPolyTP* const);

/** Prints a MultiVarExtT to stdout. */
void printmve(MultiVarExtTP const);

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
BigPolyTP* factor_BigPolyT(BigPolyTP const, BigIntTP const);

/*
int find_factors(BigIntTP const,
                 BigIntTP const,
				 BigIntTP const,
				 BigIntTP const,
			  	 BigIntTP*);
*/

#endif //ALGEBRA_H