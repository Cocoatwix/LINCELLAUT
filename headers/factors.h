
#ifndef FACTORS_H
#define FACTORS_H //Header guard

#include "bigint.h"
#include "algebra.h"

/** Returns the greatest common divisor of the two given integers. */
int GCD(int, int);

/** Returns a pointer to a BigIntT representing the GCD of the
    two given arguments. */
BigIntTP big_gcd(const BigIntTP, const BigIntTP);

/** Calculates the GCD of two given polynomials, stores the
    result in the third BigPolyT given. The fourth and fifth 
	BigPolyTs will store the polynomials s and t guaranteed by Bezout's 
	identity, meaning gcd(a, b) = sa + tb. The BigIntT is the 
	modulus to use.
	Returns 1 on success, 0 otherwise. */
int poly_gcd(const BigPolyTP, const BigPolyTP, BigPolyTP, const BigIntTP, BigPolyTP, BigPolyTP);

/** Returns the least common multiple of the two arguments. */
int LCM(int, int);

/** Returns a pointer of the given int's prime factors.
    The first number in the pointer says how many factors were found.
		It's guaranteed that the factors will be listed in ascending order.
    It's the user's job to free this pointer. */
int* prime_factors(int);

/** Returns 1 if the given integer is square-free, 0 otherwise. */
int is_square_free(int);

#endif //FACTORS_H