
#ifndef MODULAR_H
#define MODULAR_H

/** Returns the inverse of the first number
    under some given modulus (second number). 
		Returns -1 if the inverse DNE. */
int num_inverse(int, int);

/** Returns the first nonzero number that, when
    multiplied by the first argument, gives 0 when
		reduced by the second number as a modulus.
		Returns 0 if no nonzero number exists. */
int num_root(int, int);

/** Returns the square root of the first given 
    number mod some modulus. Returns -1 if the
		square root DNE. */
int square_root(int, int);

#endif //MODULAR_H