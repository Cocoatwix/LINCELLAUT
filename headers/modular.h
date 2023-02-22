
#ifndef MODULAR_H
#define MODULAR_H

#include "bigint.h"

/** Returns the inverse of the first number
    under some given modulus (second number). 
		Returns -1 if the inverse DNE. */
int num_inverse(int, int);

/** Same as num_inverse(), but for BigIntTs. Stores
    result in third given BigIntT (must be initialised).
	Stores 0 on error or if inverse doesn't exist. */
//                   findInverse     modulus         result
void big_num_inverse(const BigIntTP, const BigIntTP, BigIntTP);

/** Returns the first nonzero number that, when
    multiplied by the first argument, gives 0 when
	reduced by the second number as a modulus.
	Returns 0 if no nonzero number exists. */
int num_root(int, int);

/** Returns the square root of the first given 
    number mod some modulus. Returns -1 if the
	square root DNE. */
int square_root(int, int);

/** Given a prime-power, this function returns
    what the original prime was. If a composite is
	given, this returns the smallest prime factor.
	Returns NULL on error. */
//BigIntTP base_prime(const BigIntTP);

#endif //MODULAR_H