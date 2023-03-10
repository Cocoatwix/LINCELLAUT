
#ifndef LINALG_H
#define LINALG_H //Header guard

#include "bigint.h"  //For using BigIntTP type
#include "algebra.h" //For polynomials
#include "helper.h"  //For booleans

//Whenever you use these structs, there needs to be an interfacing function.
// C can't see "inside" the structs; it needs a function to interface with it.
/** For holding matrices. */
typedef struct intmatrix *IntMatrixTP;

/** IntMatrixT, but for BigIntT numbers. */
typedef struct bigintmatrix *BigIntMatrixTP;

/** Matrix of some generic objects. */
typedef struct genericmatrix *GenericMatrixTP;

/** Used for determining how to print vectors to the console/files. */
typedef enum vt {row, col} VectorTypeT;


//Wouldn't it have been nice if these were all the same function?

/** How many rows are in the given matrix? */
int rows(const IntMatrixTP);
int big_rows(const BigIntMatrixTP);
int gen_rows(const GenericMatrixTP);

/** How mnay columns are in the given matrix? */
int cols(const IntMatrixTP);
int big_cols(const BigIntMatrixTP);
int gen_cols(const GenericMatrixTP);

/** What's the element at the given indices in the matrix? */
int element(const IntMatrixTP, int, int);
BigIntTP big_element(const BigIntMatrixTP, int, int);

/** Increments an array of ints by a given inc. Used for
    incrementing through all possible matrices or vectors.
	Returns TRUE if array rolls over, FALSE otherwise. */
//                       array  row  col  inc  mod
bool increment_int_array(int**, int, int, int, int);

/** Increments an array of BigIntTs by a given inc. Used for
    incrementing through all possible matrices or vectors.
	Returns TRUE if array rolls over, FALSE otherwise. */
//                                array  row  col       increment         modulus
bool increment_BigIntT_array(BigIntTP**, int, int, const BigIntTP, const BigIntTP);

/** Specifies the free function to use with the given
    GenericMatrixTP. Returns 1 on success, 0 otherwise. */
int set_GenericMatrixT_freeFunction(GenericMatrixTP, void* (*)(void*));

/** Frees the memory of a matrix. Returns NULL. 
    Also returns NULL if the argument is NULL. */
IntMatrixTP     free_IntMatrixT(IntMatrixTP);
BigIntMatrixTP  free_BigIntMatrixT(BigIntMatrixTP);
GenericMatrixTP free_GenericMatrixT(GenericMatrixTP); 

/** Creates a new empty matrix of given size. 
    Returns a pointer to the matrix on success, NULL otherwise. */
IntMatrixTP     new_IntMatrixT(int, int);
BigIntMatrixTP  new_BigIntMatrixT(int, int);
GenericMatrixTP new_GenericMatrixT(int, int);

/** These functions create specialised GenericMatrixTs
    for use with specific structures. They are freed
	using the usual free_GenericMatrixT() */
//                                     rows cols numOfExtensions
GenericMatrixTP new_MultiVarExtMatrixT(int, int, int);

/** Creates a new identity matrix of given size and
    returns a pointer to the matrix. Returns NULL on error. */
IntMatrixTP identity_IntMatrixT(int);
BigIntMatrixTP identity_BigIntMatrixT(int);

/** Sets all elements in the matrix to zero. 
    Returns 1 on success, 0 otherwise. */
int clear_BigIntMatrixT(BigIntMatrixTP);
int clear_GenericMatrixT(GenericMatrixTP);

/** Sets the values of a column vector with elements specified
    in the given int pointer. Returns 1 on success, 0
    otherwise. */
int set_column(IntMatrixTP, const int*);

/** Sets the initialisation value to use with a
    GenericMatrixT's init function.
	Returns 1 on success, 0 otherwise. */
int set_GenericMatrixT_initValue(GenericMatrixTP, int);

/** Sets the function a GenericMatrixT uses to initialise
    elements of its matrix. 
	Returns 1 on success, 0 otherwise. */
int set_GenericMatrixT_initFunction(GenericMatrixTP, void* (*)(int));

/** Sets the function a GenericMatrixT uses to copy
    elements of its matrix to other elements. 
	Returns 1 on success, 0 otherwise. */
int set_GenericMatrixT_copyFunction(GenericMatrixTP, int (*)(const void*, void*));

/** Sets the print function the GenericMatrixT will use to print out
    elements of its matrix. Returns 1 on success, 0 otherwise. */
int set_GenericMatrixT_printFunction(GenericMatrixTP, void (*)(const void*));

/** Sets the function to use for clearing elements of a GenericMatrixT 
    matrix. Returns 1 on success, 0 otherwise. */
int set_GenericMatrixT_clearFunction(GenericMatrixTP, int (*)(void*));

/** Sets the function to use for printing out elements of the
    GenericMatrixT to a file stream. Returns 1 on success,
	0 otherwise. */
int set_GenericMatrixT_fprintFunction(GenericMatrixTP, void (*)(FILE*, const void*));

/** Sets the function the matrix uses to reduce its terms.
    Returns 1 on success, 0 otherwise. */
int set_GenericMatrixT_reduceFunction(GenericMatrixTP, int (*)(void*));

/** Sets the function the matrix will use to compare elements of
    its matrix. Returns 1 on success, 0 otherwise. */
int set_GenericMatrixT_compareFunction(GenericMatrixTP, int (*)(void*, void*));

/** Sets the increment function to use for adding to an element in
    its matrix. Returns 1 on success, 0 otherwise. */
int set_GenericMatrixT_incFunction(GenericMatrixTP, int (*)(const void*, void*));

/** Sets the multiplication function to use on elements of
    its matrix. Returns 1 on success, 0 otherwise. */
int set_GenericMatrixT_multFunction(GenericMatrixTP, int (*)(const void*, const void*, void*));

/** Uses the supplied initFunction to initialise the values in the
    matrix. Returns 1 on success, 0 otherwise. */
int init_GenericMatrixT(GenericMatrixTP);

/** Sets the values of a matrix to the
    values supplied in the 2D array. 
	These functions assume the given arrays are
	the correct size for the matrix.
	Returns 1 on success, 0 otherwise. */
int set_matrix(IntMatrixTP, int** const);
int set_big_matrix(BigIntMatrixTP, BigIntTP** const);
int set_GenericMatrixT(GenericMatrixTP, void** *const);

/** Reduces the given GenericMatrixT using its reduceFunction.
    Returns 1 on success, 0 otherwise. */
int reduce_GenericMatrixT(GenericMatrixTP);

/** Reads in a matrix stored in a .matrix file. 
    Returns a pointer to the matrix on success, NULL otherwise. */
IntMatrixTP read_IntMatrixT(const char*);
BigIntMatrixTP read_BigIntMatrixT(const char*);

/** Copies the values of the first matrix's matrix to the second. 
    Returns 1 on success, 0 otherwise. */
int copy_IntMatrixT(const IntMatrixTP, IntMatrixTP);
int copy_BigIntMatrixT(const BigIntMatrixTP, BigIntMatrixTP);

/** Copies one GenericMatrixT into another, copying both the
    matrix entries and the functions. The matrices must be the 
	same dimensions.
	Returns 1 on success, 0 otherwise. */
int copy_sim_GenericMatrixT(const GenericMatrixTP, GenericMatrixTP);

/** Returns 1 if the given matrix is diagonal, 0 otherwise. */
int is_diagonal(IntMatrixTP);

/** Compares two matrices to see if they're equal.
    Returns 1 if they are, 0 otherwise. */
int compare_IntMatrixT(const IntMatrixTP, const IntMatrixTP);
int compare_BigIntMatrixT(const BigIntMatrixTP, const BigIntMatrixTP);
int compare_GenericMatrixT(GenericMatrixTP, GenericMatrixTP);

/** Compares the given column of the two matrices, returns 1 if they're 
    equal. Returns 0 otherwise. */
int compare_BigIntMatrixT_cols(const BigIntMatrixTP, const BigIntMatrixTP, int);

/** Prints a given matrix to the console. */
void printm(const IntMatrixTP);
void printbm(const BigIntMatrixTP);
void printgm(const GenericMatrixTP);

void printm_row(const IntMatrixTP);
void printbm_row(const BigIntMatrixTP); //<0 0 0>
void printgm_row(const GenericMatrixTP);

/** Same as functions above, but they print to a file stream. */
void fprintm(FILE*, const IntMatrixTP);
void fprintbm(FILE*, const BigIntMatrixTP);

void fprintm_row(FILE*, const IntMatrixTP);
void fprintbm_row(FILE*, const BigIntMatrixTP); //<0 0 0>
void fprintgm_row(FILE*, const GenericMatrixTP);

/** Same as function above, except no zero padding is added. */
void fprintbm_nopad(FILE*, const BigIntMatrixTP);

/** Prints a representation of a matrix's row space to stdout. */
void big_rowsp(const BigIntMatrixTP);

/** Searches a given BigIntMatrixTP array for a matrix equivalent to the one
    passed in, then returns the reference to it. If the given matrix isn't
	in the array, it's copied to the array (by copy) and its reference is
	returned. 
	
	This function assumes all the matrices in the catalogue are of the same dimensions. */
//                                         catalogue      length  matrix to find
BigIntMatrixTP BigIntMatrixT_catalogue_get(BigIntMatrixTP**, int*, const BigIntMatrixTP);

/** Adds two matrices together, stores result in third matrix.
    Returns 1 on success, 0 otherwise. */
int big_mat_add(const BigIntMatrixTP, const BigIntMatrixTP, BigIntMatrixTP);

/** Multiples two matrices together, stores result in a third matrix. 
    Returns 1 on success, 0 otherwise. */
int mat_mul(const IntMatrixTP, const IntMatrixTP, IntMatrixTP);
int big_mat_mul(const BigIntMatrixTP, const BigIntMatrixTP, BigIntMatrixTP);

/** Matrix multiplication for GenericMatrixTs. Assumes all passed
    GenericMatrixTs are initialised correctly, and that all the elements
	in the matrices are compatiable with each other. 
	Returns 1 on success, 0 otherwise. */
int gen_mat_mul(const GenericMatrixTP, const GenericMatrixTP, GenericMatrixTP);

/** Calculates powers of first matrix, stores result in second
    matrix. Currently, they only work for positive powers.
	First int is the power, second is the modulus.
	Returns 1 on success, 0 otherwise. */
//int powm(IntMatrixTP const, IntMatrixTP, int, int);
int powbm(const BigIntMatrixTP, BigIntMatrixTP, const BigIntTP, const BigIntTP);

/** Applies a modulus to every element of a given matrix. 
    Returns 1 on success, 0 otherwise. */
int modm(IntMatrixTP, int);
int modbm(BigIntMatrixTP, const BigIntTP);

/** Takes the first matrix, plugs it into the polynomial, evaluates it,
    then stores the result in the second matrix. The function assumes
	all relevant matrices have been initialised. Returns 1 on success,
	0 otherwise. */
int eval_BigPolyT(const BigPolyTP, const BigIntMatrixTP, BigIntMatrixTP, const BigIntTP);

/** Same as eval_BigPolyT(), but the polynomial is given in factored form. */
int eval_factored_BigPolyT(const BigPolyTP*, const BigIntMatrixTP, BigIntMatrixTP, const BigIntTP);

/** Returns the determinant of a given matrix. 
    Returns zero if the matrix is nonsquare. */
int det(const IntMatrixTP);

/** Returns the inverse of a given matrix with given modulus, if it exists.
    Returns NULL otherwise. */
IntMatrixTP inverse(const IntMatrixTP, int);

/** Computes an upper-triangular form of the first matrix mod the given modulus,
    then stores the result in the second matrix. The third matrix will be altered using
	the same steps that were used to reduce the matrix given. Set the third matrix to NULL
	if not needed. Returns 1 upon computing an upper-triangular form, 0 otherwise. */
int big_row_echelon(const BigIntMatrixTP, const BigIntTP, BigIntMatrixTP, BigIntMatrixTP);

/** Takes the first matrix, in row echelon form, and tries to reduce it to 
    reduced row echelon form. Stores the result in the second matrix.
	The operations performed on the matrix will also be applied to the third
	matrix. Set the third matrix to NULL if this behaviour isn't needed.
	Returns 1 upon reducing the matrix to the identity (or some slice of it), 
	0 otherwise. */
int big_reduced_row_echelon(const BigIntMatrixTP, const BigIntTP, BigIntMatrixTP, BigIntMatrixTP);

/** Same as inverse(), but for BigIntMatrixTs. */
BigIntMatrixTP big_inverse(const BigIntMatrixTP, const BigIntTP);

/** Reduces the given IntMatrixTP to its reduced row
    echelon form. Returns 1 on success, 0 otherwise. */
int rref(IntMatrixTP, int);

/** Calculates a matrix's characteristic polynomial
    mod some modulus. Returns a BigPolyTP
	representing the characteristic polynomial on success,
	NULL otherwise. */
BigPolyTP chara_poly(const BigIntMatrixTP, const BigIntTP);

/** Calculates a matrix's minimal polynomial mod some
    prime modulus if the second BigIntMatrixTP is NULL. 
	Otherwise, calculates a vector's minimal polynomial
	under the given matrix (using the second matrix as the
	vector) mod some prime modulus. If the BigPolyTP is NULL,
	the function calculates the relavent characteristic polynomial itself.
	Otherwise, it uses the one supplied by the caller.
	Returns a BigPolyTP array representing the minimal polynomial 
	on success, NULL otherwise. 
	Output for non-prime moduli is undefined. */
BigPolyTP* min_poly(const BigIntMatrixTP, 
                    const BigIntMatrixTP, 
					const BigIntTP, 
					const BigPolyTP);

/** Stores the specified cycle converting matrix for the first one
    in the second one. This function assumes the second matrix
	has been initialised to the zero matrix.
	Returns 1 on success, 0 otherwise. */
//                     A                 CCM            from              to          modulus     
int ccm(const BigIntMatrixTP, BigIntMatrixTP, const BigIntTP, const BigIntTP, const BigIntTP);

#endif //LINALG_H
