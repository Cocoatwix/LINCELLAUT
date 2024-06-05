
#ifndef ALGEBRA_H
#define ALGEBRA_H

#include "helper.h"
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
BigPolyTP new_BigPolyT(const BigIntTP*, int);

/** Reads a polynomial from a .polynomial file specified 
    by the path given. Returns a pointer to the BigPolyT on
	success, NULL otherwise. */
BigPolyTP read_BigPolyT(const char*);

/** Creates a constant polynomial with the given BigIntTP
    as the constant term. Returns a pointer to the new
	constant. */
BigPolyTP constant_BigPolyT(const BigIntTP);

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
    pointer to it. The integer specifies how many
	extensions to keep track of. Returns NULL on error. */
void* new_MultiVarExtT(int);

/** Sets the modulus value for the given MultiVarExtT.
    Returns 1 on success, 0 otherwise. */
int set_MultiVarExtT_mod(MultiVarExtTP, const BigIntTP);

/** Adds a new extension to a MultiVarExtT.
    Returns 1 on success, 0 otherwise. */
//                extension      minPoly          size name
int add_extension(MultiVarExtTP, const BigIntTP*, int, const char*);

/** Removes the most recent extension from the
    given MultiVarExtT.
	Returns 1 on success, 0 otherwise. */
int remove_extension(MultiVarExtTP);

/** Returns a pointer to the coefficient at the given
    position. Each element of the int pointer specifies what power
	of a particular extension we're looking at on the term we're
	interested in. 
	
	For example, if we wanted the coefficient on the (α^7)(β^2)(γ)
	term, where α is the first extension, β is the second, and γ
	is the third, then our int pointer should be {7, 2, 1}.
	
	Returns NULL on error. */
BigIntTP get_MultiVarExtT_coefficient(const MultiVarExtTP, const int*);

/** Sets a particular value for the given MultiVarExtT's coefficient.
    The MultiVarExtT must be fully set before this function can be used.
	See above comment on get_MultiVarExtT_coefficient() to understand
	how to format the int pointer.
	Returns 1 on success, 0 otherwise. */
int set_MultiVarExtT_coefficient(MultiVarExtTP, const int*, const BigIntTP);

/** Increments a MultiVarExt by 1. Once the constant term
    rolls over the modulus, the first extension is incremented,
	then that extensions powers, then the next extension, etc.
	
	Returns 1 if the ENTIRE extension rolls over after calling this
	function, 0 otherwise.
	
	This function assumes the extension has been fully initialised and set. */
bool increment_MultiVarExtT(MultiVarExtTP);

/** Uses the extension definitions to reduce the given
    MultiVarExtT. Essentially, it'll zero any terms
	it can. 
	This function can only be used when all extensions
	have been set for the MultiVarExtT.
	Returns 1 on success, 0 otherwise. */
int reduce_MultiVarExtT(void*);

/** Returns the degree of the given polynomial, independent
    of the internal size of the polynomial. */
int degree(const BigPolyTP);

/** Returns the constant term of a BigPolyTP by reference. */
BigIntTP constant(const BigPolyTP);

/** Returns the leading term's coefficient of the given BigPolyTP. 
    This function assumes the polynomial has been reduced. 
	This can be done by calling reduce_BigPolyT(). */
BigIntTP leading_term(const BigPolyTP);

/** Returns an array of the given BigPolyT's coefficients,
    sorted in ascending order of the terms' exponents.
	Returns NULL on error. */
BigIntTP* extract_coefficients(const BigPolyTP);

/** Copies a polynomial into another. This function assumes
    all arguments have been properly initialised.
	Returns 1 on success, 0 otherwise. */
int copy_BigPolyT(const BigPolyTP, BigPolyTP);

/** Copies a MultiVarExtT into another MultiVarExtT.
    Both MultiVarExtTs must have the same number of 
	extensions. The function assumes the second
	MultiVarExtT has been initialised. 
	Returns 1 on success, 0 otherwise. */
int copy_MultiVarExtT(const void*, void*);

/** Create a new BigFactorsT struct with given factors,
    exponents, and size. Returns a pointer to the struct. */
BigFactorsTP new_BigFactorsT(const BigPolyTP*, const int*, int);

/** Adds a new factor to the given BigFactorsT struct. 
    Returns 1 on success, 0 otherwise. */
int add_factor(BigFactorsTP, const BigPolyTP, int);

/** Returns a BigPolyTP array containing all the
    factors in the given BigFactorsT. The factor will
	appear as many times as its exponent demands. */
BigPolyTP* extract_factors(const BigFactorsTP);

/** Same as extract_factors(), but repeated factors will be multipled
    together into one factor. */
BigPolyTP* extract_coprime_factors(const BigFactorsTP);

/** Returns the number of factors in the given BigFactorsTP,
    taking into account repeated roots. */
int count_factors(const BigFactorsTP);

/** Counts the number of unique factors in the given
    BigFactorsTP. This function assumes there are no
	duplicate factors within the BigFactorsTP factors
	array. */
int count_unique_factors(const BigFactorsTP);

/** If any repeated factors are treated as different factors
    within the BigFactorsTP factors array, this function will
	collect them into a single factor. Returns 1 on success, 
	0 otherwise. */
int collect_factors(BigFactorsTP);

/** Compares two BigPolyTPs to see if they're equal (have the
    same coefficients). 
	Returns 0 if they're equal, a negative if the first BigPolyT 
	is less than the second, a positive if the first BigPolyT is 
	greater than the second. */ 
int compare_BigPolyT(const BigPolyTP, const BigPolyTP);

/** Compares MultiVarExtTs to check whether their 
    expressions are equal. Returns 1 if they're equal, 
	0 otherwise. 
	
	This function reduces the two MultiVarExtTs before
	comparing. */
int compare_MultiVarExtT(void*, void*);

/** Uses the given list of BigIntTs to set the polynomial's
    coefficients. 
	Returns 1 on success, 0 otherwise. */
int set_BigPolyT(BigPolyTP, const BigIntTP*);

/** Sets all coefficients in the given object
    to zero. Returns 1 on success, 0 otherwise. */
int clear_BigPolyT(BigPolyTP);
int clear_MultiVarExtT(void*);

/** Outputs a BigPolyT to stdout (the console). */
void printp(const BigPolyTP);

/** Prints a factored BigPolyT. */
void old_printpf(const BigPolyTP*);

/** Prints a BigFactorsT to stdout. */
void printpf(const BigFactorsTP);

/** Prints a MultiVarExtT. */
void printmve(const MultiVarExtTP);
void fprintmve(FILE*, const MultiVarExtTP);

/** Prints only a MultiVarExtT's expression,
    not its extension definitions. */
void printmve_row(const void*);
void fprintmve_row(FILE*, const void*);

/** Same as printp, printpf, and old_fprintpf but it prints to a file stream. */
void fprintp(FILE*, const BigPolyTP);
void fprintpf(FILE*, const BigFactorsTP);
void old_fprintpf(FILE*, const BigPolyTP*);

/** Adds the first two polynomials, stores the sum in the third.
    This function assumes all polynomials passed have been
		initialised.
    Returns 1 on success, 0 otherwise. */
int add_BigPolyT(const BigPolyTP, const BigPolyTP, BigPolyTP);

/** Multiples the first two polynomials given, stores the
    product in the third BigPolyTP. This function assumes all 
		polynomials passed have been initialised.
		Returns 1 on success, 0 otherwise. */
int multiply_BigPolyT(const BigPolyTP, const BigPolyTP, BigPolyTP);

/** Divides two BigPolyTs, stores the quotient in the third BigPolyT,
    stores the remainder in the forth BigPolyT. The BigIntTP is the 
	modulus used. This function assumes the quotient and remainder
	have been initialised.
	Returns 1 on success, 0 otherwise. */
//For now, this function assumes the modulus given is prime (so we're in a field)
//This function also assumes a and b have been reduced
//                  a                b                a/b        remainder  modulus
int divide_BigPolyT(const BigPolyTP, const BigPolyTP, BigPolyTP, BigPolyTP, const BigIntTP);

/** Calculates the first BigPolyT raised to the given BigIntT.
    Stores the result in the second BigPolyT given (which is
	assumed to have been initialised).
	Currently, this function works for nonnegative powers.
	Returns 1 on success, 0 otherwise. */
int pow_BigPolyT(const BigPolyTP, const BigIntTP, BigPolyTP);
int modulo_pow_BigPolyT(const BigPolyTP, const BigIntTP, const BigIntTP, BigPolyTP);

/** Reduces a given polynomial by a modulus, stores the result in
    another passed polynomial. This function assumes all arguments
	have been properly initialised.
	Returns 1 on success, 0 oterwise. */
int mod_BigPolyT(const BigPolyTP, const BigIntTP, BigPolyTP);

/** Differentiates the given BigPolyTP with respect to its
    variable, stores the result in the second given BigPolyTP.
	Returns 1 on success, 0 otherwise. */
int diff_BigPolyT(const BigPolyTP, BigPolyTP);

/** Factors the given BigPolyT under the given modulus.
    Returns a pointer to all the polynomial's factors.
		The first element is a BigPolyT constant which says
		how many factors are in the pointer. */
BigPolyTP* old_factor_BigPolyT(const BigPolyTP, const BigIntTP);

/** Factors polynomials under the given modulus (assumed prime). 
	Returns a BigFactorsTP representing the factorisation
	of the given polynomial.*/
BigFactorsTP factor_BigPolyT(const BigPolyTP, const BigIntTP);

/** Adds two MultiVarExtTs together, assuming they have the same size
    and the same extensions (in the same order). Stores the result
	in the second MultiVarExtT.
	The MultiVarExtTs also have to be fully set.
	This function DOES NOT concern itself with taking
	any sort of modular reduction during addition. Use
	reduce_MultiVarExtT() for that.
	Returns 1 on success, 0 otherwise. */
int inc_sim_MultiVarExtT(const void*, void*);

/** Multiplies two MultiVarExtTs with the same extensions together, stores the
    result in the third MultiVarExtTP given. All MultiVarExtTs must have the
	same extensions in the same order, and must be fully set (i.e. can set
	coefficients). 
	This function assumes the modulus used is a prime.
	Returns 1 on success, 0 otherwise. */
int mult_sim_MultiVarExtT(const void*, const void*, void*);

/*
int find_factors(BigIntTP const,
                 BigIntTP const,
				 BigIntTP const,
				 BigIntTP const,
			  	 BigIntTP*);
*/

#endif //ALGEBRA_H
