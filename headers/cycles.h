
#ifndef CYCLES_H
#define CYCLES_H //Header guard

#include "linalg.h" //Allows us to use IntMatrixTP datatypes

/** Holds info about a cycle calculated using
    Floyd's Cycle Detection Algorithm. */
typedef struct cycleinfo *CycleInfoTP;

/** Frees the memory allocated to the given CycleInfo struct.
    Returns NULL. */
CycleInfoTP free_CycleInfoT(CycleInfoTP);

/** Prints out the information contained in the
    CycleInfoT object pointed to. */
void printcycle(CycleInfoTP);

/** Multiples a vector by a given matrix n times, returns
    the final vector. 
		<update>, <initial>, <mod>, <iter>*/
IntMatrixTP iterate(IntMatrixTP, IntMatrixTP, int, int);

/** Performs Floyd's Cycle Detection Algorithm on the
    given setup, returns info about the cycle. 
		Returns NULL on error. */
CycleInfoTP floyd(IntMatrixTP, IntMatrixTP, int);

/** Creates a .iteration file containing a single iteration of
    every vector under the given update matrix and modulus. 
		Returns 1 on success, 0 otherwise. */
int write_iteration(const char*, IntMatrixTP, int);

#endif //CYCLES_H
