
#ifndef LINALG_H
#define LINALG_H //Header guard

#include "bigint.h"  //For using BigIntTP type
#include "algebra.h" //For polynomials
#include "helper.h"  //For booleans

//Whenever you use this struct, there needs to be an interfacing function.
// C can't see "inside" the struct; it needs a function to interface with it.
/** For holding matrices. */
typedef struct intmatrix *IntMatrixTP;

/** IntMatrixT, but for BigIntT numbers. */
typedef struct bigintmatrix *BigIntMatrixTP;

/** Used for determining how to print vectors to the console/files. */
typedef enum vt {row, col} VectorTypeE;


/** How many rows are in the given matrix? */
int rows(IntMatrixTP);
int big_rows(BigIntMatrixTP);

/** How mnay columns are in the given matrix? */
int cols(IntMatrixTP);
int big_cols(BigIntMatrixTP);

/** What's the element at the given indices in the matrix? */
int element(IntMatrixTP const, int, int);
BigIntTP big_element(BigIntMatrixTP const, int, int);

/** Increments an array of ints by a given inc. Used for
    incrementing through all possible matrices or vectors.
	Returns TRUE if array rolls over, FALSE otherwise. */
//                       array  row  col  inc  mod
bool increment_int_array(int**, int, int, int, int);

/** Increments an array of BigIntTs by a given inc. Used for
    incrementing through all possible matrices or vectors.
	Returns TRUE if array rolls over, FALSE otherwise. */
//                                array  row  col       increment         modulus
bool increment_BigIntT_array(BigIntTP**, int, int, BigIntTP const, BigIntTP const);

/** Frees the memory of a matrix. Returns NULL. 
    Also returns NULL if the argument is NULL. */
IntMatrixTP free_IntMatrixT(IntMatrixTP);
BigIntMatrixTP free_BigIntMatrixT(BigIntMatrixTP);

/** Creates a new empty matrix of given size. 
    Returns a pointer to the matrix on success, NULL otherwise. */
IntMatrixTP new_IntMatrixT(int, int);
BigIntMatrixTP new_BigIntMatrixT(int, int);

/** Creates a new identity matrix of given size and
    returns a pointer to the matrix. Returns NULL on error. */
IntMatrixTP identity_IntMatrixT(int);
BigIntMatrixTP identity_BigIntMatrixT(int);

/** Sets the values of a column vector with elements specified
    in the given int pointer. Returns 1 on success, 0
    otherwise. */
int set_column(IntMatrixTP, int* const);

/** Sets the values of a matrix to the
    values supplied in the 2D array. Returns 1 on success,
	0 otherwise. */
int set_matrix(IntMatrixTP, int** const);
int set_big_matrix(BigIntMatrixTP, BigIntTP** const);

/** Reads in a matrix stored in a .matrix file. 
    Returns a pointer to the matrix on success, NULL otherwise. */
IntMatrixTP read_IntMatrixT(char* const);
BigIntMatrixTP read_BigIntMatrixT(char* const);

/** Copies the values of the first matrix's matrix to the second. 
    Returns 1 on success, 0 otherwise. */
int copy_IntMatrixT(IntMatrixTP const, IntMatrixTP);
int copy_BigIntMatrixT(BigIntMatrixTP const, BigIntMatrixTP);

/** Returns 1 if the given matrix is diagonal, 0 otherwise. */
int is_diagonal(IntMatrixTP);

/** Compares two matrices to see if they're equal.
    Returns 1 if they are, 0 otherwise. */
int compare_IntMatrixT(IntMatrixTP const, IntMatrixTP const);
int compare_BigIntMatrixT(BigIntMatrixTP const, BigIntMatrixTP const);

/** Compares the given column of the two matrices, returns 1 if they're 
    equal. Returns 0 otherwise. */
int compare_BigIntMatrixT_cols(BigIntMatrixTP const M1, BigIntMatrixTP const M2, int);

/** Prints a given matrix to the console. */
void printm(IntMatrixTP const);
void printbm(BigIntMatrixTP const);

void printm_row(IntMatrixTP const);
void printbm_row(BigIntMatrixTP const); //<0 0 0>

/** Same as functions above, but they print to a file stream. */
void fprintm(FILE*, IntMatrixTP const);
void fprintbm(FILE*, BigIntMatrixTP const);

void fprintm_row(FILE*, IntMatrixTP const);
void fprintbm_row(FILE*, BigIntMatrixTP const); //<0 0 0>

/** Same as function above, except no zero padding is added. */
void fprintbm_nopad(FILE*, BigIntMatrixTP const);

/** Searches a given BigIntMatrixTP array for a matrix equivalent to the one
    passed in, then returns the refernece to it. If the given matrix isn't
	in the array, it's copied to the array (by copy) and its reference is
	returned. 
	
	This function assumes all the matrices in the catalogue are of the same dimensions. */
//                                         catalogue      length  matrix to find
BigIntMatrixTP BigIntMatrixT_catalogue_get(BigIntMatrixTP**, int*, BigIntMatrixTP const);

/** Adds two matrices together, stores result in third matrix.
    Returns 1 on success, 0 otherwise. */
int big_mat_add(BigIntMatrixTP const, BigIntMatrixTP const, BigIntMatrixTP);

/** Multiples two matrices together, stores result in a third matrix. 
    Returns 1 on success, 0 otherwise. */
int mat_mul(IntMatrixTP const, IntMatrixTP const, IntMatrixTP);
int big_mat_mul(BigIntMatrixTP const, BigIntMatrixTP const, BigIntMatrixTP);

/** Calculates powers of first matrix, stores result in second
    matrix. Currently, they only work for positive powers.
	First int is the power, second is the modulus.
	Returns 1 on success, 0 otherwise. */
//int powm(IntMatrixTP const, IntMatrixTP, int, int);
int powbm(BigIntMatrixTP const, BigIntMatrixTP, BigIntTP const, BigIntTP const);

/** Applies a modulus to every element of a given matrix. 
    Returns 1 on success, 0 otherwise. */
int modm(IntMatrixTP, int);
int modbm(BigIntMatrixTP, BigIntTP const);

/** Returns the determinant of a given matrix. 
    Returns zero if the matrix is nonsquare. */
int det(IntMatrixTP const);

/** Returns the inverse of a given matrix with given modulus, if it exists.
    Returns NULL otherwise. */
IntMatrixTP inverse(IntMatrixTP const, int);

/** Same as inverse(), but for BigIntMatrixTs. */
BigIntMatrixTP big_inverse(BigIntMatrixTP const, BigIntTP const);

/** Reduces the given IntMatrixTP to its reduced row
    echelon form. Returns 1 on success, 0 otherwise. */
int rref(IntMatrixTP, int);

/** Calculates a matrix's characteristic polynomial
    mod some modulus. Returns a BigPolyTP
		representing the characteristic polynomial on success,
		NULL otherwise. */
BigPolyTP chara_poly(BigIntMatrixTP const, BigIntTP const);

/** Stores the specified cycle converting matrix for the first one
    in the second one. This function assumes the second matrix
		has been initialised to the zero matrix.
		Returns 1 on success, 0 otherwise. */
//                     A                 CCM            from              to          modulus     
int ccm(BigIntMatrixTP const, BigIntMatrixTP, BigIntTP const, BigIntTP const, BigIntTP const);



/* These two functions are badly tested and specific in their
    use cases. Only use them if you know what you're doing. */

/** Returns the eigenvalues of a given matrix, mod some
    modulus, as an int pointer. Returns NULL if no eigenvalues
		exist. 
		
		This currently only works for 2x2 matrices. */
//int* eigenvalues(IntMatrixTP const, int);

/** Returns an eigenvector of the given matrix, eigenvalue, and
    modulus. It is assumed the given eigenvalue is valid for the
		given system. 
		
		Currently only works for 2x2 matrices. */
//IntMatrixTP eigenvector(IntMatrixTP const, int, int);

#endif //LINALG_H