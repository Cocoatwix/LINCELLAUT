
#ifndef ALGEBRA_H
#define ALGEBRA_H

#include "bigint.h"

/** Holds a polynomial datatype using BigIntT coefficients. */
typedef struct bigpoly *BigPolyTP;

/** Holds a factored BigPolyT, essentially. */
typedef struct bigfactors *BigFactorsTP;

/** Responsible for facilitating multidimensional BigIntT arrays
    of varying dimension. */
//typedef struct bigintdirection *BigIntDirectionTP;

/** Holds a representation of a multivariate field extension
    polynomial expression. */
typedef struct multivarext *MultiVarExtTP;

/** Frees the memory used by a BigPolyT. Returns NULL. */
BigPolyTP free_BigPolyT(BigPolyTP);

/** Frees the memory used by a BigFactorsT. Returns NULL. */
BigFactorsTP free_BigFactorsT(BigFactorsTP);

/** Frees an array of BigPolyT factors. */
BigPolyTP* free_BigPolyT_factors(BigPolyTP*);

/** Frees the memory used by a MultiVarExtT. Returns NULL. */
void* free_MultiVarExtT(void*);

/** Creates a new BigPolyT with given coefficients and size,
    returns a pointer to it. */
BigPolyTP new_BigPolyT(BigIntTP* const, int);

/** Creates a constant polynomial with the given BigIntTP
    as the constant term. Returns a pointer to the new
		constant. */
BigPolyTP constant_BigPolyT(BigIntTP const);

/** Creates an empty polynomial, returns a pointer to it. */
BigPolyTP empty_BigPolyT();

/** Returns a pointer to an empty BigFactorsT struct.
    Returns NULL on error. */
BigFactorsTP empty_BigFactorsT();

/** Ensures that the leading term of the BigPolyT
    is nonzero. Returns 1 on success, 0 otherwise. */
int reduce_BigPolyT(BigPolyTP);

/** Resizes a given BigPolyT to the given size.
    Returns 1 on success, 0 otherwise. */
int resize_BigPolyT(BigPolyTP, int);

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

/** Removes the most recent extension from the
    given MultiVarExtT.
	Returns 1 on success, 0 otherwise. */
int remove_extension(MultiVarExtTP);

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

/** Returns the leading term's coefficient of the given BigPolyTP. */
BigIntTP leading_term(BigPolyTP const);

/** Returns an array of the given BigPolyT's coefficients,
    sorted in ascending order of the terms' exponents.
	Returns NULL on error. */
BigIntTP* extract_coefficients(BigPolyTP const);

/** Copies a polynomial into another. This function assumes
    all arguments have been properly initialised.
		Returns 1 on success, 0 otherwise. */
int copy_BigPolyT(BigPolyTP const, BigPolyTP);

/** Copies a MultiVarExtT into another MultiVarExtT.
    Both MultiVarExtTs must have the same number of 
	extensions. The function assumes the second
	MultiVarExtT has been initialised. 
	Returns 1 on success, 0 otherwise. */
int copy_MultiVarExtT(const MultiVarExtTP, MultiVarExtTP);

/** Create a new BigFactorsT struct with given factors,
    exponents, and size. Returns a pointer to the struct. */
BigFactorsTP new_BigFactorsT(BigPolyTP* const, int* const, int);

/** Adds a new factor to the given BigFactorsT struct. 
    Returns 1 on success, 0 otherwise. */
int add_factor(BigFactorsTP, BigPolyTP const, int);

/** Returns a BigPolyTP array containing all the
    factors in the given BigFactorsT. The factor will
	appear as many times as its exponent demands. */
BigPolyTP* extract_factors(BigFactorsTP const);

/** Returns the number of factors in the given BigFactorsTP,
    taking into account repeated roots. */
int count_factors(BigFactorsTP const);

/** Compares two BigPolyTPs to see if they're equal (have the
    same coefficients). 
	Returns 0 if they're equal, a negative if the first BigPolyT 
	is less than the second, a positive if the first BigPolyT is 
	greater than the second. */ 
int compare_BigPolyT(BigPolyTP const, BigPolyTP const);

/** Uses the given list of BigIntTs to set the polynomial's
    coefficients. 
	Returns 1 on success, 0 otherwise. */
int set_BigPolyT(BigPolyTP, BigIntTP* const);

/** Sets all coefficients in the given polynomial
    to zero. Returns 1 on success, 0 otherwise. */
int clear_BigPolyT(BigPolyTP);

/** Outputs a BigPolyT to stdout (the console). */
void printp(BigPolyTP const);

/** Prints a factored BigPolyT. */
void old_printpf(BigPolyTP* const);

/** Prints a BigFactorsT to stdout. */
void printpf(BigFactorsTP const);

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

/** Divides two BigPolyTs, stores the quotient in the third BigPolyT,
    stores the remainder in the forth BigPolyT. The BigIntTP is the 
	modulus used. This function assumes the quotient and remainder
	have been initialised.
	Returns 1 on success, 0 otherwise. */
//For now, this function assumes the modulus given is prime (so we're in a field)
//This function also assumes a and b have been reduced
//                  a                b                a/b        remainder  modulus
int divide_BigPolyT(BigPolyTP const, BigPolyTP const, BigPolyTP, BigPolyTP, BigIntTP const);

/** Calculates the first BigPolyT raised to the given BigIntT.
    Stores the result in the second BigPolyT given (which is
	assumed to have been initialised).
	Currently, this function works for nonnegative powers.
	Returns 1 on success, 0 otherwise. */
int pow_BigPolyT(BigPolyTP const, BigIntTP const, BigPolyTP);

/** Reduces a given polynomial by a modulus, stores the result in
    another passed polynomial. This function assumes all arguments
	have been properly initialised.
	Returns 1 on success, 0 oterwise. */
int mod_BigPolyT(BigPolyTP const, BigIntTP const, BigPolyTP);

/** Differentiates the given BigPolyTP with respect to its
    variable, stores the result in the second given BigPolyTP.
	Returns 1 on success, 0 otherwise. */
int diff_BigPolyT(BigPolyTP const, BigPolyTP);

/** Factors the given BigPolyT under the given modulus.
    Returns a pointer to all the polynomial's factors.
		The first element is a BigPolyT constant which says
		how many factors are in the pointer. */
BigPolyTP* old_factor_BigPolyT(BigPolyTP const, BigIntTP const);

/** Factors polynomials under the given modulus. 
	Returns an array of BigPolyTs, each one representing an
	irreducible factor of the given polynomial. The first BigPolyT 
	in the array is a constant telling how many factors are in the array. */
BigFactorsTP factor_BigPolyT(BigPolyTP const, BigIntTP const);

/** Adds two MultiVarExtTs together, assuming they have the same size
    and the same extensions (in the same order). Stores the result
	in the second MultiVarExtT.
	The MultiVarExtTs also have to be fully set.
	This function DOES NOT concern itself with taking
	any sort of modular reduction during addition. Use
	reduce_MultiVarExtT() for that.
	Returns 1 on success, 0 otherwise. */
int inc_sim_MultiVarExtT(MultiVarExtTP const, MultiVarExtTP);

/** Multiplies two MultiVarExtTs with the same extensions together, stores the
    result in the third MultiVarExtTP given. All MultiVarExtTs must have the
	same extensions in the same order, and must be fully set (i.e. can set
	coefficients). 
	This function assumes the modulus used is a prime.
	Returns 1 on success, 0 otherwise. */
int mult_sim_MultiVarExtT(MultiVarExtTP const, MultiVarExtTP const, MultiVarExtTP);

/*
int find_factors(BigIntTP const,
                 BigIntTP const,
				 BigIntTP const,
				 BigIntTP const,
			  	 BigIntTP*);
*/

#endif //ALGEBRA_H