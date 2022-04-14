
#ifndef LINALG_H
#define LINALG_H //Header guard

//Whenever you use this struct, there needs to be an interfacing function.
// C can't see "inside" the struct; it needs a function to interface with it.
typedef struct intmatrix *IntMatrixTP;


/** How many rows are in the given matrix? */
int rows(IntMatrixTP);

/** How mnay columns are in the given matrix? */
int cols(IntMatrixTP);

/** What's the element at the given indices in the matrix? */
int element(IntMatrixTP, int, int);


/** Frees the memory of an IntMatrix. Returns NULL. */
IntMatrixTP free_IntMatrixT(IntMatrixTP);

/** Reads in a matrix stored in a .matrix file and
    stores it in a given IntMatrixT struct. This
		function allocates space for M. Returns 1
		on success, 0 otherwise. */
int read_IntMatrixT(char* const, IntMatrixTP*);

/** Copies the values of the first matrix's matrix to the second. 
    Returns 1 on success, 0 otherwise. */
int copy_IntMatrixT(IntMatrixTP const, IntMatrixTP);

/** Compares two matrices to see if they're equal.
    Returns 1 if they are, 0 otherwise. */
int compare_IntMatrixT(IntMatrixTP const, IntMatrixTP const);

/** Creates a new empty IntMatrix of given size, stores it at given pointer. 
    Returns 1 on success, 0 otherwise. */
int new_IntMatrixT(IntMatrixTP*, int, int);

/** Prints a given matrix of given size to the console. */
void printm(IntMatrixTP);

/** Multiples two matrices together, stores result in a third matrix. 
    Returns 1 on success, 0 otherwise. */
int mat_mul(const IntMatrixTP, const IntMatrixTP, IntMatrixTP);

/** Applies a modulus to every element of a given matrix. 
    Returns 1 on success, 0 otherwise. */
int modm(IntMatrixTP, int);

/** Returns the determinant of a given matrix. 
    Returns zero if the matrix is nonzero. */
int det(IntMatrixTP const);

/** Returns the inverse of a given matrix, if it exists.
    Returns NULL otherwise. */
IntMatrixTP inverse(IntMatrixTP);

#endif //LINALG_H