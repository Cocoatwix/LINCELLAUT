
#ifndef FACTORS_H
#define FACTORS_H //Header guard

/** Returns the greatest common divisor of the two given integers. */
int GCD(int, int);

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