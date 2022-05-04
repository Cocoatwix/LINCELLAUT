
#ifndef LINALG_H
#define LINALG_H //Header guard

//Adding TRUE or FALSE values
typedef enum boolean {FALSE, TRUE} bool;

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

/** Creates a new empty IntMatrix of given size. 
    Returns a pointer to the matrix on success, NULL otherwise. */
IntMatrixTP new_IntMatrixT(int, int);

/** Creates a new identity matrix of given size and
    returns a pointer to the matrix. Returns NULL on error. */
IntMatrixTP identity_IntMatrixT(int);

/** Sets the values of a column vector with elements specified
    in the given int pointer. Returns 1 on success, 0
    otherwise. */
int set_column(IntMatrixTP, int* const);

/** Sets the values of an IntMatrixT's matrix to the
    values supplied in the 2D array. Returns 1 on success,
		0 otherwise. */
int set_matrix(IntMatrixTP, int** const);

/** Reads in a matrix stored in a .matrix file. 
    Returns a pointer to the matrix on success, NULL otherwise. */
IntMatrixTP read_IntMatrixT(char* const);

/** Copies the values of the first matrix's matrix to the second. 
    Returns 1 on success, 0 otherwise. */
int copy_IntMatrixT(IntMatrixTP const, IntMatrixTP);

/** Compares two matrices to see if they're equal.
    Returns 1 if they are, 0 otherwise. */
int compare_IntMatrixT(IntMatrixTP const, IntMatrixTP const);

/** Prints a given matrix of given size to the console. */
void printm(IntMatrixTP, bool);

/** Returns the number of digits a given
    positive integer has. */
int num_digits(int);

/** Multiples two matrices together, stores result in a third matrix. 
    Returns 1 on success, 0 otherwise. */
int mat_mul(IntMatrixTP const, IntMatrixTP const, IntMatrixTP);

/** Applies a modulus to every element of a given matrix. 
    Returns 1 on success, 0 otherwise. */
int modm(IntMatrixTP, int);

/** Returns the determinant of a given matrix. 
    Returns zero if the matrix is nonzero. */
int det(IntMatrixTP const);

/** Returns the inverse of a given matrix with given modulus, if it exists.
    Returns NULL otherwise. */
IntMatrixTP inverse(IntMatrixTP const, int);

#endif //LINALG_H