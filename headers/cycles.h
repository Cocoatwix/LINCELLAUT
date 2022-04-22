
#ifndef CYCLES_H
#define CYCLES_H //Header guard

#include "linalg.h" //Allows us to use IntMatrixTP datatypes

/** Holds info about a cycle calculated using
    Floyd's Cycle Detection Algorithm. */
typedef struct cycleinfo *CycleInfoTP;

/** Prints out the information contained in the
    CycleInfoT object pointed to. */
void printcycle(CycleInfoTP);

/** Multiples a vector by a given matrix n times, returns
    the final vector. */
IntMatrixTP iterate(IntMatrixTP, IntMatrixTP, int, int);

/** Creates a .orbit file containing the orbits of
    every vector under the given update matrix and modulus. */
int write_orbits(const char*, IntMatrixTP, int);

/** Performs Floyd's Cycle Detection Algorithm on the
    given setup, returns info about the cycle. */
CycleInfoTP floyd(IntMatrixTP, IntMatrixTP, int);

#endif //CYCLES_H