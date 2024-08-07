
#ifndef CYCLES_H
#define CYCLES_H //Header guard

#include "linalg.h" //Allows us to use IntMatrixTP datatypes

/** Holds info about a cycle calculated using
    Floyd's Cycle Detection Algorithm. */
typedef struct cycleinfo *CycleInfoTP;

/** Frees the memory allocated to the given CycleInfo struct.
    Returns NULL. */
CycleInfoTP free_CycleInfoT(CycleInfoTP);

/** Returns a pointer to a new CycleInfoT struct. */
CycleInfoTP new_CycleInfoT();

/** Returns the cycle length of a given CycleInfoT object. */
int omega(CycleInfoTP);

/** Returns the transient length of a given CycleInfoT object. */
int tau(CycleInfoTP);

/** Returns a pointer to a cycle's inCycle matrix/vector, if it
    has one. Returns NULL otherwise. */
void* rep(CycleInfoTP);

/** Prints out the information contained in the
    CycleInfoT object pointed to. */
void printcycle(CycleInfoTP, VectorTypeT);

/** Multiples a matrix by a given matrix some times, returns
    the final iterated matrix. Returns NULL if F isn't square. 
		<update>, <initial>, <modulus>, <iterations>*/
IntMatrixTP iterate(IntMatrixTP, IntMatrixTP, int, int);

/** Sees what points F ends up at. Points are basically 
    matrices of interest. */
void visit_points(IntMatrixTP, int, int);

/** Performs Floyd's Cycle Detection Algorithm on the
    given setup, returns info about the cycle. 
		Returns NULL on error. */
CycleInfoTP floyd(const IntMatrixTP, const IntMatrixTP, int);

/** Same as floyd(), except with BigIntMatrixT structs.
    Returns 1 on success, 0 otherwise. */
//First matrix is the update matrix, second is the initial matrix
int big_floyd(const BigIntMatrixTP, 
              const BigIntMatrixTP, 
			  const BigIntTP,
			  CycleInfoTP*);
			  
/** Same as floyd(), but with GenericMatrixTs.
    Returns 1 on success, 0 otherwise. */
int generic_floyd(GenericMatrixTP, GenericMatrixTP, CycleInfoTP*);

/** Stores a list of orbit representatives under the given update
    matrix and modulus within the given BigIntMatrixTP array.
	This function assumes the given array has been initialised to NULL.
	If the int** != NULL, then the cycle lengths for each orbit rep
	will also be recorded. It assumes the int* is initialised to NULL.
	Returns the number of orbit reps found on success, -1 otherwise. */
int big_orbit_reps(const BigIntMatrixTP, const BigIntTP, BigIntMatrixTP**, int**);

/** Creates a .iteration file containing a single iteration of
    every vector under the given update matrix and modulus. 
		Returns 1 on success, 0 otherwise. */
int write_iteration(const char*, IntMatrixTP, int);

#endif //CYCLES_H
