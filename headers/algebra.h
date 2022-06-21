
#ifndef ALGEBRA_H
#define ALGEBRA_H

#include "bigint.h"

/** Holds a polynomial datatype using BigIntT coefficients. */
typedef struct bigpoly *BigPolyTP;

/** Frees the memory used by a BigPolyT. Returns NULL. */
BigPolyTP free_BigPolyT(BigPolyTP);

/** Creates a new BigPolyT with given coefficients and size,
    returns a pointer to it. */
BigPolyTP new_BigPolyT(BigIntTP* const, int);

/** Creates an empty polynomial, returns a pointer to it. */
BigPolyTP empty_BigPolyT();

/** Returns the degree of the given polynomial. */
int degree(BigPolyTP const);

/** Outputs a BigPolyT to stdout (the console). */
void printp(BigPolyTP const);

/** Multiples the first two polynomials given, stores the
    product in the third BigPolyTP. Returns 1 on success, 0 otherwise. */
int multiply_BigPolyT(BigPolyTP const, BigPolyTP const, BigPolyTP);

#endif //ALGEBRA_H