
#ifndef ALGEBRA_H
#define ALGEBRA_H

#include "bigint.h"

/** Holds a polynomial datatype using BigIntT coefficients. */
typedef struct bigpoly *BigPolyTP;

/** Holds an expression containing field extensions. */
typedef struct fieldexp *FieldExpTP;

/** Frees the memory used by a BigPolyT. Returns NULL. */
BigPolyTP free_BigPolyT(BigPolyTP);

/** Frees the memory used by a FieldExpT. Returns NULL. */
FieldExpTP free_FieldExpT(FieldExpTP);

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

/** Creates a new FieldExpT object with <int> unique elements, 
    returns a pointer to it. */
FieldExpTP new_FieldExpT(int);

/** Adds a new extension into our expression with given minPoly, value, and
    name. Returns 1 on success, 0 otherwise. */
//                fieldexp    minpoly          value            name
int add_extension(FieldExpTP, BigPolyTP const, BigPolyTP const, char* const);

/** For a field extension a, find the minimum exponent c such that
    a^c is in the original field. In this case, the field is the
	integers modulo the BigIntT given. This number must be prime,
	otherwise the behaviour of this function is undefined. The 
	BigPolyT given specifies the minimal polynomial for the field 
	extension. 
	Returns the exponent c on success, -1 otherwise. */
//This function assumes the polynomial has been reduced.
int collapse_field_extension(BigPolyTP const, BigIntTP const);

/** Ensures that the leading term of the BigPolyT
    is nonzero. Returns 1 on success, 0 otherwise. */
int reduce_BigPolyT(BigPolyTP);

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

/** Prints a field expression. */
void printfe(FieldExpTP const);

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