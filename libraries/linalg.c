
/* Apr 7, 2022
 *
 */
 
#include <stdlib.h>
#include <stdio.h>

#include "../headers/helper.h"

#include "../headers/factors.h"
#include "../headers/modular.h" //Allows for taking modular square roots
#include "../headers/bigint.h"  //For BigIntMatrixT
#include "../headers/algebra.h" //For polynomials

/* The following resources were used as a reference:
https://stackoverflow.com/a/6317375
https://www.c-programming-simple-steps.com/c-extern.html
https://stackoverflow.com/questions/13239369/
*/

//Maybe I could make a generic struct that has a union as its matrix
//That way, if we wanted to switch to doubles, we could
//Or maybe a simple void pointer will do

extern const int MAXBUNCH;
extern const int MAXBUNCH;

typedef struct intmatrix
/** Structure for holding matrix data and metadata (dimensions, e.g.). */
{
	int** matrix;
	int m;
	int n;
}
IntMatrixT, *IntMatrixTP;

typedef struct bigintmatrix
/** IntMatrixT struct, but for BigIntT numbers. */
{
	BigIntTP** matrix;
	int m;
	int n; //Should really make these BigIntTs as well, but oh well
} 
BigIntMatrixT, *BigIntMatrixTP;

typedef struct genericmatrix
/** Matrix of some generic object. */
{
	void*** matrix;
	int m;
	int n;
	
	bool isInitialised; //Says whether its matrix elements are initialised
	
	int initValue; //For feeding into the initialisation function
	
	void* (*freeFunction)(void*);
	void* (*initFunction)(int);
	int   (*copyFunction)(const void*, void*);
	void  (*printFunction)(const void*); 
	int   (*clearFunction)(void*); 
	void  (*fprintFunction)(FILE*, const void*);
	int   (*reduceFunction)(void*);
	int   (*compareFunction)(void*, void*);
	
	int   (*incFunction)(const void*, void*);
	int   (*multFunction)(const void*, const void*, void*);
} GenericMatrixT, *GenericMatrixTP;


int rows(const IntMatrixTP M)
/** Returns the number of rows in a given matrix. */
{
	return M->m;
}


int big_rows(const BigIntMatrixTP M)
/** Returns the number of rows in the given matrix. */
{
	return M->m;
}


int gen_rows(const GenericMatrixTP M)
/** Returns the number of rows in M. */
{
	return M->m;
}


int cols(const IntMatrixTP M)
/** Returns the number of columns in a given matrix. */
{
	return M->n;
}


int big_cols(const BigIntMatrixTP M)
/** Returns the number of columns in the given matrix. */
{
	return M->n;
}


int gen_cols(const GenericMatrixTP M)
/** Returns the number of columns in M. */
{
	return M->n;
}


int element(const IntMatrixTP M, int row, int col)
/** Returns the element at M->matrix[row][col].
    Returns 0 for elements that aren't defined. */
{
	if ((row >= M->m) || (col >= M->n) ||
	    (row < 0) || (col < 0))
		return 0;
		
	return M->matrix[row][col];
}


bool increment_int_array(int** intArr, int sizeRow, int sizeCol, int inc, int mod)
/** Increments the given size by size int array by inc. 
    Used to increment through all possible matrices or vectors 
		under a particular modulus.
		Returns TRUE if the array rolls over,
		FALSE otherwise. */
{
	//I could totally restructure this to be more efficient instead of a copy-paste
	// from the BigIntTP version, but this works well enough
	
	int temp  = 0;
	int temp2 = 0;
	int carry = 0;
	bool onceRolledOver = FALSE;
	bool rolledOver = TRUE;
	
	//Initial incrementation
	temp2 = intArr[0][0] + inc;
	
	//While loop allows us to go through the array
	// multiple times if needed
	while (rolledOver)
	{
		for (int row = 0; row < sizeRow; row += 1)
		{
			for (int col = 0; col < sizeCol; col += 1)
			{
				temp = temp2 + carry;
				
				if (temp >= mod)
				{
					temp2 = temp % mod;
					intArr[row][col] = temp2;
					carry = temp / mod;
					
					//Prepare temp2 for carry
					if (col+1 < sizeCol)
						temp2 = intArr[row][col+1];
					else if (row+1 < sizeRow)
						temp2 = intArr[row+1][0];
					else
					{
						temp2 = intArr[0][0];
						carry -= 1; //Prevents the zero vector from being skipped when rolling over
					}
				}
				
				else
				{
					intArr[row][col] = temp;
					rolledOver = FALSE;
					break;
				}
			}
			
			if (!rolledOver)
				break;
		}
		
		if (rolledOver)
			onceRolledOver = TRUE;
	}
	
	return onceRolledOver;
}


/* private */ bool general_increment_BigIntT_array(BigIntTP** intArr, 
																									 int sizeRow, 
																									 int sizeCol, 
																									 const BigIntTP inc, 
																									 const BigIntTP mod,
																									 const BigIntTP elemCarry)
{
	BigIntTP temp  = empty_BigIntT(1);
	BigIntTP temp2 = empty_BigIntT(1);
	BigIntTP carry = empty_BigIntT(1);
	
	int oneArr[1] = {1};
	BigIntTP one = new_BigIntT(oneArr, 1);
	
	bool onceRolledOver = FALSE;
	bool rolledOver = TRUE;

	//Initial incrementation
	add_BigIntT(intArr[0][0], inc, temp2);
	
	//While loop allows us to go through the array
	// multiple times if needed
	while (rolledOver)
	{
		for (int row = 0; row < sizeRow; row += 1)
		{
			for (int col = 0; col < sizeCol; col += 1)
			{
				add_BigIntT(temp2, carry, temp);
				
				if (compare_BigIntT(temp, mod) >= 0)
				{
					mod_BigIntT(temp, mod, temp2);
					copy_BigIntT(temp2, intArr[row][col]);
					
					if (elemCarry == NULL)
						divide_BigIntT(temp, mod, carry);
					else
						copy_BigIntT(elemCarry, carry);
					
					//Prepare temp2 for carry
					if (col+1 < sizeCol)
						copy_BigIntT(intArr[row][col+1], temp2);
					else if (row+1 < sizeRow)
						copy_BigIntT(intArr[row+1][0], temp2);
					else
					{
						copy_BigIntT(intArr[0][0], temp2);
						
						if (elemCarry == NULL)
						{
							subtract_BigIntT(carry, one, temp); //Prevents the zero vector from being skipped when incrementing
							copy_BigIntT(temp, carry);
						}
						else
							clear_BigIntT(carry);
					}
				}
				
				else
				{
					copy_BigIntT(temp, intArr[row][col]);
					rolledOver = FALSE;
					break;
				}
			}
			
			if (!rolledOver)
				break;
		}
		
		if (rolledOver)
			onceRolledOver = TRUE;
	}
	
	one   = free_BigIntT(one);
	temp  = free_BigIntT(temp);
	temp2 = free_BigIntT(temp2);
	carry = free_BigIntT(carry);
	
	return onceRolledOver;
}


bool increment_BigIntT_array(BigIntTP** intArr, 
                             int sizeRow, 
														 int sizeCol, 
														 const BigIntTP inc, 
														 const BigIntTP mod)
/** Increments the given size by size BigIntTP array by inc. 
    Used to increment through all possible matrices or vectors 
		under a particular modulus.
		Returns TRUE if the array rolls over,
		FALSE otherwise. */
{
	return general_increment_BigIntT_array(intArr, sizeRow, sizeCol, inc, mod, NULL);
}


bool step_BigIntT_array(BigIntTP** intArr, 
                             int sizeRow, 
														 int sizeCol, 
														 const BigIntTP inc, 
														 const BigIntTP mod)
/** Same as increment_BigIntT_array(), but the carry
    between elements in the array is the same as the
		incrementation.
		Returns TRUE if the array rolls over,
		FALSE otherwise. */
{
	return general_increment_BigIntT_array(intArr, sizeRow, sizeCol, inc, mod, inc);
}


BigIntTP big_element(const BigIntMatrixTP M, int row, int col)
/** Returns the BigIntTP contained within M at the 
    specified location. Returns NULL if the specified index
		is out of bounds of the given matrix. */
{
	if ((row >= M->m) || (row < 0) ||
	    (col >= M->n) || (col < 0))
		return NULL;
		
	return M->matrix[row][col];
}


IntMatrixTP free_IntMatrixT(IntMatrixTP M)
/** Frees memory taken by M. Returns NULL. */
{
	if (M == NULL)
		return NULL;
	
	//Freeing each row of the matrix
	for (int i = 0; i < M->m; i += 1)
		free(M->matrix[i]);
	
	//Freeing the matrix itself
	free(M->matrix);
	free(M);
	
	return NULL;
}


BigIntMatrixTP free_BigIntMatrixT(BigIntMatrixTP M)
/** Frees memory used by the given matrix. Returns NULL. */
{
	if (M == NULL)
		return NULL;
	
	for (int x = 0; x < M->m; x += 1)
	{
		for (int y = 0; y < M->n; y += 1)
			M->matrix[x][y] = free_BigIntT(M->matrix[x][y]);
		
		free(M->matrix[x]);
	}
	
	free(M->matrix);
	free(M);

	return NULL;
}


int set_GenericMatrixT_freeFunction(GenericMatrixTP a, void* (*theFunction)(void*))
/** Sets the free function to be used for a specific 
    GenericMatrixT. Returns 1 on success, 0 otherwise. */
{
	a->freeFunction = theFunction;
	return 1;
}


GenericMatrixTP free_GenericMatrixT(GenericMatrixTP a)
/** Frees the memory used by a GenericMatrixT. Returns NULL. */
{
	if (a != NULL)
	{
		if (a->freeFunction == NULL)
			fprintf(stderr, "No free function provided to GenericMatrixT. Memory may not free correctly.\n");
		
		for (int row = 0; row < a->m; row += 1)
		{
			for (int col = 0; col < a->n; col += 1)
			{
				if (a->freeFunction != NULL)
					a->matrix[row][col] = a->freeFunction(a->matrix[row][col]);
			}
			
			free(a->matrix[row]);
			a->matrix[row] = NULL;
		}
		free(a->matrix);
		a->matrix = NULL;
		
		a->freeFunction = NULL;
		a->initFunction = NULL;
		a->copyFunction = NULL;
		a->printFunction = NULL;
		a->clearFunction = NULL;
		a->fprintFunction = NULL;
		a->reduceFunction = NULL;
		a->compareFunction = NULL;
		
		a->incFunction = NULL;
		a->multFunction = NULL;
		
		free(a);
	}
	return NULL;
}


IntMatrixTP new_IntMatrixT(int r, int c)
/** Creates an empty m by n matrix.
    Returns a pointer to the matrix on success, NULL otherwise. */
{
	//If given dimensions don't make sense
	if ((r < 0) || (c < 0))
		return NULL;
	
	
	IntMatrixTP M = malloc(sizeof(IntMatrixT));
	
	M->matrix = malloc(r*sizeof(int*));
	
	for (int i = 0; i < r; i += 1)
		M->matrix[i] = calloc(c, sizeof(int));
	
	//Don't forget to set the matrix's dimensions!
	M->m = r;
	M->n = c;
	
	return M;
}


BigIntMatrixTP new_BigIntMatrixT(int r, int c)
/** Creates an empty r by c matrix. Returns a pointer to
    the newly created matrix. Returns NULL on error. */
{
	if ((r <= 0) || (c <= 0))
		return NULL;
	
	BigIntMatrixTP newMat = malloc(sizeof(BigIntMatrixT));
	newMat->matrix = malloc(r*sizeof(BigIntTP*));
	
	for (int x = 0; x < r; x += 1)
	{
		newMat->matrix[x] = malloc(c*sizeof(BigIntTP));
		for (int y = 0; y < c; y += 1)
			newMat->matrix[x][y] = empty_BigIntT(1);
	}
	
	newMat->m = r;
	newMat->n = c;
	
	return newMat;
}


GenericMatrixTP new_GenericMatrixT(int r, int c)
/** Creates a new GenericMatrixT with dimensions
    r by c. Returns a pointer to this new matrix. */
{
	GenericMatrixTP a = malloc(sizeof(GenericMatrixT));
	a->m = r;
	a->n = c;
	
	a->matrix = malloc(r*sizeof(void**));
	for (int i = 0; i < r; i += 1)
		a->matrix[i] = calloc(c, sizeof(void*));
	
	a->isInitialised = FALSE;
	a->initValue = 0;
	
	a->freeFunction = NULL;
	a->initFunction = NULL;
	a->copyFunction = NULL;
	a->printFunction = NULL;
	a->clearFunction = NULL;
	a->fprintFunction = NULL;
	a->reduceFunction = NULL;
	a->compareFunction = NULL;
	
	a->incFunction = NULL;
	a->multFunction = NULL;
	
	return a;
}


int init_GenericMatrixT(GenericMatrixTP a)
/** Uses the previously supplied init function to initialise
    the values of the matrix. Returns 1 on success, 0 otherwise. */
{
	for (int i = 0; i < a->m; i += 1)
		for (int j = 0; j < a->n; j += 1)
			a->matrix[i][j] = a->initFunction(a->initValue);
		
	a->isInitialised = TRUE;

	return 1;
}


GenericMatrixTP new_MultiVarExtMatrixT(int r, int c, int numOfExtensions)
/** Automates the creation of a GenericMatrixT meant for
    holding MultiVarExtTs.
		Returns a pointer to the matrix on success, NULL otherwise. */
{
	GenericMatrixTP theMatrix = new_GenericMatrixT(r, c);

	theMatrix->freeFunction    = free_MultiVarExtT;
	theMatrix->initValue       = numOfExtensions;
	theMatrix->initFunction    = new_MultiVarExtT;
	theMatrix->copyFunction    = copy_MultiVarExtT;
	theMatrix->printFunction   = printmve_row;
	theMatrix->clearFunction   = clear_MultiVarExtT;
	theMatrix->fprintFunction  = fprintmve_row;
	theMatrix->compareFunction = compare_MultiVarExtT;
	theMatrix->reduceFunction  = reduce_MultiVarExtT;

	theMatrix->incFunction  = inc_sim_MultiVarExtT;
	theMatrix->multFunction = mult_sim_MultiVarExtT;
	
	init_GenericMatrixT(theMatrix);
	
	return theMatrix;
}


IntMatrixTP identity_IntMatrixT(int r)
/** Creates an identity matrix of size r by r and returns
    a pointer to it. Returns NULL on error. */
{
	if (r <= 0)
		return NULL;
	
	IntMatrixTP I = new_IntMatrixT(r, r);
	
	for (int i = 0; i < r; i += 1)
		I->matrix[i][i] = 1;
	
	return I;
}


BigIntMatrixTP identity_BigIntMatrixT(int r)
/** Same as identity_IntMatrixT(), but for
    BigIntMatrixTs. */
{
	if (r <= 0)
		return NULL;
	
	BigIntMatrixTP I = new_BigIntMatrixT(r, r);
	int oneArr[1] = {1};
	BigIntTP one = new_BigIntT(oneArr, 1);
	
	for (int i = 0; i < r; i += 1)
		copy_BigIntT(one, I->matrix[i][i]);
	
	one = free_BigIntT(one);
	
	return I;
}


int clear_BigIntMatrixT(BigIntMatrixTP A)
/** Sets all elements in the matrix to zero.
    Returns 1 on success, 0 otherwise. */
{
	BigIntTP zero = empty_BigIntT(1);
	
	for (int r = 0; r < A->m; r += 1)
		for (int c = 0; c < A->n; c += 1)
			copy_BigIntT(zero, A->matrix[r][c]);
	
	zero = free_BigIntT(zero);
	
	return 1;
}


int clear_GenericMatrixT(GenericMatrixTP A)
/** Clears all the elements in A's matrix.
    Returns 1 on success, 0 otherwise. */
{
	if (A->clearFunction == NULL)
		return 0;
	
	for (int row = 0; row < A->m; row += 1)
		for (int col = 0; col < A->n; col += 1)
			A->clearFunction(A->matrix[row][col]);
		
	return 1;
}


int set_column(IntMatrixTP v, const int* elements)
/** Returns a pointer to a new IntMatrixTP vector that contains
    the elements specified in elements. Returns NULL on error. */
{
	for (int i = 0; i < v->m; i += 1)
		v->matrix[i][0] = elements[i];
	
	return 1;
}


int set_big_row(BigIntMatrixTP A, const BigIntTP* newRow, int rowNum)
/** Sets A->matrix[row] to the values in newRow.
    Returns 1 on success, 0 otherwise. */
{
	if ((A == NULL) || (newRow == NULL))
		return 0;
	
	if ((rowNum < 0) || (rowNum >= A->m))
		return 0;
	
	for (int c = 0; c < A->n; c += 1)
		copy_BigIntT(newRow[c], A->matrix[rowNum][c]);
	
	return 1;
}


int set_GenericMatrixT_initFunction(GenericMatrixTP a, void* (*f)(int))
/** Sets the initialisation function for a GenericMatrixT
    (what it uses to initalise new elements of its matrix).
		Returns 1 on success, 0 otherwise. */
{
	a->initFunction = f;
	return 1;
}


int set_GenericMatrixT_initValue(GenericMatrixTP a, int v)
/** Sets the initialisation value for a GenericMatrixT's
    initialisation function.
		Returns 1 on success, 0 otherwise. */
{
	a->initValue = v;
	return 1;
}


int set_GenericMatrixT_copyFunction(GenericMatrixTP a, int (*f)(const void*, void*))
/** Sets the copying function for a GenericMatrixT
    (what it uses to copy elements to elements of its matrix).
		Returns 1 on success, 0 otherwise. */
{
	a->copyFunction = f;
	return 1;
}


int set_GenericMatrixT_printFunction(GenericMatrixTP a, void (*f)(const void*))
/** Sets the function a will use to print elements of its matrix.
    Returns 1 on success, 0 otherwise. */
{
	a->printFunction = f;
	return 1;
}


int set_GenericMatrixT_clearFunction(GenericMatrixTP a, int (*f)(void*))
/** Sets the function a will use to clear out elements of its matrix.
    Returns 1 on success, 0 otherwise. */
{
	a->clearFunction = f;
	return 1;
}


int set_GenericMatrixT_fprintFunction(GenericMatrixTP a, void (*f)(FILE*, const void*))
/** Sets the function a will use to print elements of its matrix
    to a file stream.
    Returns 1 on success, 0 otherwise. */
{
	a->fprintFunction = f;
	return 1;
}


int set_GenericMatrixT_reduceFunction(GenericMatrixTP a, int (*f)(void*))
/** Sets the function the matrix will use to reduce its elements.
    Returns 1 on success, 0 otherwise. */
{
	a->reduceFunction = f;
	return 1;
}


int set_GenericMatrixT_compareFunction(GenericMatrixTP a, int (*f)(void*, void*))
/** Sets the compare function the matrix will use to compare its
    elements. Returns 1 on success, 0 otherwise. */
{
	a->compareFunction = f;
	return 1;
}


int set_GenericMatrixT_incFunction(GenericMatrixTP a, int (*f)(const void*, void*))
/** Sets the increment function to use for adding to an
    element of a's matrix. Returns 1 on success, 0 otherwise. */
{
	a->incFunction = f;
	return 1;
}


int set_GenericMatrixT_multFunction(GenericMatrixTP a, int (*f)(const void*, const void*, void*))
/** Sets the multiplication function to use for multiplying
    two elements of a's matrix. Returns 1 on success, 0 otherwise. */
{
	a->multFunction = f;
	return 1;
}




int set_matrix(IntMatrixTP A, int** const arr)
/* Using the given array, sets the matrix of a given
   initialised IntMatrixT. Returns 1 on success, 0 otherwise. 
	 
	 It is assumed arr will be the correct dimensions for A. */
{
	for (int row = 0; row < A->m; row += 1)
		for (int col = 0; col < A->n; col += 1)
			A->matrix[row][col] = arr[row][col];
		
	return 1;
}


int set_big_matrix(BigIntMatrixTP A, BigIntTP** const arr)
/** Sets the elements of the matrix using the array of BigIntTP
    values. Returns 1 on success, 0 otherwise. 
		
		It's assumed arr will be the correct dimensions for A. */
{
	for (int x = 0; x < A->m; x += 1)
		for (int y = 0; y < A->n; y += 1)
			copy_BigIntT(arr[x][y], A->matrix[x][y]);
		
	return 1;
}


int set_GenericMatrixT(GenericMatrixTP a, void** *const arr)
/** Sets the matrix of a given GenericMatrixT.
    The matrix must first be initalised.
    Returns 1 on success, 0 otherwise. */
{
	if (a->copyFunction == NULL)
		return 0;
	
	if (a->isInitialised != TRUE)
		return 0;
	
	for (int row = 0; row < a->m; row += 1)
		for (int col = 0; col < a->n; col += 1)
			a->copyFunction(arr[row][col], a->matrix[row][col]);
	
	return 1;
}


int reduce_GenericMatrixT(GenericMatrixTP a)
/** Uses a's reduceFunction to reduce its matrix.
    Returns 1 on success, 0 otherwise. */
{
	if (a->reduceFunction == NULL)
		return 0;
	
	for (int row = 0; row < a->m; row += 1)
		for (int col = 0; col < a->n; col += 1)
			a->reduceFunction(a->matrix[row][col]);
	
	return 1;
}


IntMatrixTP read_IntMatrixT(const char* matFilePath)
/** Initialises the given matrix M from the given .matrix file.
    This method assumes matFile points to the beginning of the file.
		Returns a pointer to the matrix upon success, NULL otherwise. */
{
	FILE* matFile = fopen(matFilePath, "r");
	
	if (matFile == NULL)
		return NULL;
	
	//Allocating memory for M
	IntMatrixTP M = malloc(sizeof(IntMatrixT));
	
	//Getting dimensions of the data
	if (fscanf(matFile, "%d %d", &(M->m), &(M->n)) != 2)
	{
		fclose(matFile);
		free(M);
		return NULL;
	}
	
	//Allocating memory for M
	M->matrix = malloc((M->m)*sizeof(int*));
	for (int i = 0; i < M->m; i += 1)
		M->matrix[i] = malloc((M->n)*sizeof(int));
	
	//Filling F with actual values
	for (int row = 0; row < M->m; row += 1)
		for (int column = 0; column < M->n; column += 1)
			if (fscanf(matFile, "%d", &(M->matrix[row][column])) != 1)
			{
				M = free_IntMatrixT(M);
				return NULL;
			}
			

	if (fclose(matFile) == EOF)
	{
		M = free_IntMatrixT(M);
		return NULL;
	}
		
	//Function read matrix's data properly
	return M;
}


BigIntMatrixTP read_BigIntMatrixT(const char* matFilePath)
/** Same as read_IntMatrixT(), but for BigIntMatrixT structs.
    Returns a pointer to the new matrix on success, NULL
		otherwise. */
{
	BigIntMatrixTP M, I;
	char* tempStr;
	
	FILE* matFile = fopen(matFilePath, "r");
	
	if (matFile == NULL)
		return NULL;
	
	M = malloc(sizeof(BigIntMatrixT));
	
	//Get dimensions of the matrix
	if (fscanf(matFile, "%d %d", &(M->m), &(M->n)) != 2)
	{
		fclose(matFile);
		free(M);
		return NULL;
	}
	
	//Make space for the actual matrix data
	M->matrix = malloc((M->m)*sizeof(BigIntTP*));
	for (int row = 0; row < M->m; row+= 1)
	{
		M->matrix[row] = malloc((M->n)*sizeof(BigIntTP));
		
		//Setting matrix elements to NULL to prevent trying to
		// free nothing if M gets freed
		for (int col = 0; col < M->n; col += 1)
			M->matrix[row][col] = NULL;
	}
	
	tempStr = malloc(100*sizeof(char));
	
	//Now, read the data
	for (int row = 0; row < M->m; row += 1)
	{
		for (int col = 0; col < M->n; col += 1)
		{
			if (fscanf(matFile, "%100s", tempStr) != 1)
			{
				free(tempStr);
				fclose(matFile);
				free_BigIntMatrixT(M);
				return NULL;
			}
			
			//Check to see if the user is creating an identity matrix
			if (tempStr[0] == 'I')
			{
				I = identity_BigIntMatrixT(M->m);
				free(tempStr);
				fclose(matFile);
				free_BigIntMatrixT(M);
				return I;
			}
			
			//Actually store our value in the matrix
			strtoBIT(tempStr, &(M->matrix[row][col]));
		}
	}
	
	if (fclose(matFile) == EOF)
	{
		free(tempStr);
		free_BigIntMatrixT(M);
		return NULL;
	}
	
	free(tempStr);
	return M;
}


int copy_IntMatrixT(const IntMatrixTP toCopy, IntMatrixTP copyTo)
/** Copies toCopy to copyTo. This function assumes copyTo has
    already been initialised. Returns 1 on success, 0 otherwise. */
{
	//If the dimensions of the matrices don't match
	if ((toCopy->m != copyTo->m) || 
	    (toCopy->n != copyTo->n))
		return 0;
		
	//Copy toCopy to copyTo
	for (int row = 0; row < toCopy->m; row += 1)
		for (int col = 0; col < toCopy->n; col += 1)
			copyTo->matrix[row][col] = toCopy->matrix[row][col];
		
	return 1;
}


int copy_BigIntMatrixT(const BigIntMatrixTP toCopy, BigIntMatrixTP copyTo)
/** Same as copy_IntMatrixT(), but for BigIntMatrixT structs.
    Returns 1 on success, 0 otherwise. */
{
	if ((toCopy->m != copyTo->m) ||
	    (toCopy->n != copyTo->n))
		return 0;
		
	for (int row = 0; row < copyTo->m; row += 1)
		for (int col = 0; col < copyTo->n; col += 1)
			copy_BigIntT(toCopy->matrix[row][col], copyTo->matrix[row][col]);
		
	return 1;
}


int copy_GenericMatrixT(const GenericMatrixTP toCopy, GenericMatrixTP copyTo)
/** Same as copy_IntMatrixT(), but for GenericMatrixT structs.
    Returns 1 on success, 0 otherwise. 
		
		This function assumes both GenericMatrixTs have all the relevant
		functions set. */
{
	printf("copy_GenericMatrixT() is not yet implemented. This function will now clear \
	given GenericMatrixTs...\n");
	clear_GenericMatrixT(toCopy);
	clear_GenericMatrixT(copyTo);
	return 0;
}


int copy_sim_GenericMatrixT(const GenericMatrixTP toCopy, GenericMatrixTP copyTo)
/** Copies both toCopy's elements and toCopy's functions into
    copyTo. Both matrices must be the same size for this to
		work.
		Returns 1 on success, 0 otherwise. */
{
	if ((toCopy->m != copyTo->m) || (toCopy->n != copyTo->n))
		return 0;
	
	copyTo->initFunction = toCopy->initFunction;
	copyTo->initValue = toCopy->initValue;
	
	//Making sure we'll be able to copy values
	// by-value into copyTo
	if (!copyTo->isInitialised)
		init_GenericMatrixT(copyTo);
	
	copyTo->copyFunction = toCopy->copyFunction;
	for (int row = 0; row < toCopy->m; row += 1)
		for (int col = 0; col < toCopy->n; col += 1)
			copyTo->copyFunction(toCopy->matrix[row][col], copyTo->matrix[row][col]);
		
	copyTo->freeFunction = toCopy->freeFunction;
	copyTo->printFunction = toCopy->printFunction;
	copyTo->clearFunction = toCopy->clearFunction;
	copyTo->fprintFunction = toCopy->fprintFunction;
	copyTo->reduceFunction = toCopy->reduceFunction;
	copyTo->compareFunction = toCopy->compareFunction;
	
	copyTo->incFunction = toCopy->incFunction;
	copyTo->multFunction = toCopy->multFunction;
	
	return 1;
}


int is_diagonal(IntMatrixTP A)
/** Returns 1 if A is diagonal, 0 otherwise. */
{
	for (int row = 0; row < A->m; row += 1)
		for (int col = 0; col < A->n; col += 1)
			if (col != row)
				if (A->matrix[row][col] != 0)
					return 0;
				
	return 1;
}


int compare_IntMatrixT(const IntMatrixTP M1, const IntMatrixTP M2)
/** Compares two matrices to see if they're equal.
    Returns 1 if the two are equivalent, 0 otherwise. */
{
	//Checking to see if dimensions are equal
	if ((M1->m != M2->m) ||
	    (M1->n != M2->n))
		return 0;
		
	//Check each entry to see if any are not equivalent
	for (int row = 0; row < M1->m; row += 1)
		for (int col = 0; col < M1->n; col += 1)
			if (M1->matrix[row][col] != M2->matrix[row][col])
				return 0;
			
	return 1;
}


int compare_BigIntMatrixT(const BigIntMatrixTP M1, const BigIntMatrixTP M2)
/** Same as compare_IntMatrixT(), but with BigIntMatrixT structs.
    Returns 1 if the two matrices are equal, zero otherwise. */
{
	if ((M1 == NULL) || (M2 == NULL))
		return 0;
	
	if ((M1->m != M2->m) || (M1->n != M2->n))
		return 0;
	
	for (int row = 0; row < M1->m; row += 1)
		for (int col = 0; col < M1->n; col += 1)
			if (compare_BigIntT(M1->matrix[row][col], M2->matrix[row][col]) != 0)
				return 0;
			
	return 1;
}


int compare_GenericMatrixT(GenericMatrixTP a, GenericMatrixTP b)
/** Compares a and b to see if their matrices are equal.
    Returns 1 if they are, 0 otherwise. */
{
	if ((a->m != b->m) || (a->n != b->n))
		return 0;
	
	//If the datatypes within the matrices can't be compared
	if (a->compareFunction != b->compareFunction)
		return 0;
	
	for (int row = 0; row < a->m; row += 1)
		for (int col = 0; col < a->n; col += 1)
			if (! a->compareFunction(a->matrix[row][col], b->matrix[row][col]))
				return 0;

	return 1;
}


int compare_BigIntMatrixT_cols(const BigIntMatrixTP M1, const BigIntMatrixTP M2, int c)
/** Compares the given column of the two matrices to see if they're equal.
    Returns 1 if they're equal, 0 otherwise. */
{
	//If one of the matrices isn't big enough for the given column,
	// or if the row dimensions aren't equal
	if ((M1->n <= c) || (M2->n <= c) ||
	    (M1->m != M2->m))
		return 0;
	
	for (int row = 0; row < M1->m; row += 1)
		if (compare_BigIntT(M1->matrix[row][c], M2->matrix[row][c]) != 0)
			return 0;
		
	return 1;
}


int resize_BigIntMatrixT(BigIntMatrixTP A, int newRows, int newCols)
/** Takes A and resizes it to the new given coords.
    Returns 1 on success, 0 otherwise. */
{
	if ((newRows < 1) || (newCols < 1))
		return 0;
	
	if (A == NULL)
		return 0;
	
	//No resizing to do here
	if ((A->m == newRows) && (A->n == newCols))
		return 1;
	
	int lowerRowBound = (A->m > newRows) ? newRows : A->m;
	
	//Adjust columns
	for (int r = 0; r < lowerRowBound; r += 1)
	{
		//If we need to shrink our number of columns
		if (A->n > newCols)
		{
			for (int c = newCols; c < A->n; c += 1)
				A->matrix[r][c] = free_BigIntT(A->matrix[r][c]);
			
			A->matrix[r] = realloc(A->matrix[r], newCols*sizeof(BigIntTP));
		}
		
		//If our columns are expanding
		else
		{
			A->matrix[r] = realloc(A->matrix[r], newCols*sizeof(BigIntTP));
			
			for (int c = A->n; c < newCols; c += 1)
				A->matrix[r][c] = empty_BigIntT(1);
		}
	}
	
	//Adjust rows
	if (A->m > newRows)
	{
		for (int r = newRows; r < A->m; r += 1)
		{
			for (int c = 0; c < A->n; c += 1)
				A->matrix[r][c] = free_BigIntT(A->matrix[r][c]);
			
			free(A->matrix[r]);
		}
		A->matrix = realloc(A->matrix, newRows*sizeof(BigIntTP*));
	}
	
	else
	{
		A->matrix = realloc(A->matrix, newRows*sizeof(BigIntTP*));
		for (int r = A->m; r < newRows; r += 1)
		{
			A->matrix[r] = malloc(newCols*sizeof(BigIntTP));
			for (int c = 0; c < newCols; c += 1)
				A->matrix[r][c] = empty_BigIntT(1);
		}
	}
	
	A->m = newRows;
	A->n = newCols;
	return 1;
}


/* private */ IntMatrixTP det_subIntMatrixT(const IntMatrixTP M, int x, int y)
/** Stores a submatrix of M in N. x and y specify the
    center of the "crosshairs" where the matrix is sliced. Returns 1
		on success, 0 otherwise.
		Freeing the memory created by this function is the programmer's 
		responsibility.
		See documentation for more details. */
{
	//If the given indices are out of the matrix's bounds
	if ((x < 0) || (y < 0) || (x >= M->m) || (y >= M->n))
		return NULL;
	
	//A row and column will be sliced from M, so N's dimensions need to be smaller
	IntMatrixTP N = new_IntMatrixT(M->m-1, M->n-1);
	
	//Extra indices to ensure we skip over the correct rows/cols in M
	int Nrow = 0;
	int Ncol = 0;
	
	for (int Mrow = 0; Mrow < M->m; Mrow += 1)
	{
		if (Mrow != x) //Skipping over sliced row
		{
			Ncol = 0;
			for (int Mcol = 0; Mcol < M->n; Mcol += 1)
			{
				if (Mcol != y) //Skipping over sliced column
				{
					N->matrix[Nrow][Ncol] = M->matrix[Mrow][Mcol];
					Ncol += 1;
				}
			}
			Nrow += 1;
		}
	}
	
	return N;
}
	

void fprintm(FILE* file, const IntMatrixTP M)
/** Same as printm, except it prints to a file stream. */
{
	int maxDigits = 0;
	
	for (int row = 0; row < M->m; row += 1)
		for (int col = 0; col < M->n; col += 1)
			if (num_digits(M->matrix[row][col]) > maxDigits)
				maxDigits = num_digits(M->matrix[row][col]);
	
	for (int mIndex = 0; mIndex < M->m; mIndex += 1)
	{
		for (int nIndex = 0; nIndex < M->n-1; nIndex += 1)
		{
			for (int p = 0; p < maxDigits - num_digits(M->matrix[mIndex][nIndex]); p += 1)
				fprintf(file, "0");
			
			fprintf(file, "%d, ", M->matrix[mIndex][nIndex]);
		}
		
		for (int p = 0; p < maxDigits - num_digits(M->matrix[mIndex][M->n-1]); p += 1)
				fprintf(file, "0");

		fprintf(file, "%d", M->matrix[mIndex][M->n-1]);
		fprintf(file, "\n");
	}
}


void printm(const IntMatrixTP M)
/** Prints an m by n matrix to stdout. */
{
	fprintm(stdout, M);
}


void fprintm_row(FILE* file, const IntMatrixTP M)
/** Same as printm_row(), except it prints vectors to a given file. */
{
	//Making sure we're actually dealing with a column vector
	if ((M->n == 1) && (M->m >= 1))
	{
		fprintf(file, "<");
		for (int mIndex = 0; mIndex < M->m-1; mIndex += 1)
			fprintf(file, "%d ", M->matrix[mIndex][0]);
		
		fprintf(file, "%d>", M->matrix[M->m-1][0]);
	}
	
	else
		fprintm(file, M);
}


void printm_row(const IntMatrixTP M)
/** Same as printm(), except it prints vectors as row vectors. 
    This function also ditches the zero padding from the original. */
{
	fprintm_row(stdout, M);
}


void fprintbm(FILE* file, BigIntMatrixTP M)
/** Prints out a BigIntMatrixT matrix to the console. */
{
	int maxBunches = 0;
	int maxDigits = 0;
	int tempBunch;
	int tempSize;
	
	//Search for the highest amount of bunches in the matrix
	for (int x = 0; x < M->m; x += 1)
		for (int y = 0; y < M->n; y += 1)
		{
			if (maxBunches < size(M->matrix[x][y]))
				maxBunches = size(M->matrix[x][y]);

			tempSize = size(M->matrix[x][y])-1;
			tempBunch = extract_bunch(M->matrix[x][y], tempSize);
			if (maxDigits < num_digits(tempBunch) + tempSize*MAXBUNCHDIGITS)
				maxDigits = num_digits(tempBunch) + tempSize*MAXBUNCHDIGITS;
		}

	//Actually print out the entries
	for (int row = 0; row < M->m; row += 1)
	{
		for (int col = 0; col < M->n; col += 1)
		{
			//Zero padding
			tempSize = size(M->matrix[row][col])-1;
			for (int i = 0; 
			     i < maxDigits - (tempSize*MAXBUNCHDIGITS + num_digits(extract_bunch(M->matrix[row][col], tempSize))); 
					 i += 1)
				fprintf(file, "0");
				
			fprinti(file, M->matrix[row][col]);
			if (col != M->n-1)
				fprintf(file, " ");
		}
		fprintf(file, "\n");
	}
}


void printbm(const BigIntMatrixTP M)
/** Prints out a BigIntMatrixT matrix to the console. */
{
	fprintbm(stdout, M);
}


void fprintbm_row(FILE* file, const BigIntMatrixTP M)
/** Same as printbm_row(), but prints to a file stream. */
{
	if ((M->m > 1) && (M->n > 1))
		fprintbm(file, M);
	
	else if (M->m == 1)
		fprintbm(file, M);
	
	else if (M->n == 1)
	{
		fprintf(file, "<");
		for (int i = 0; i < M->m-1; i += 1)
		{
			fprinti(file, M->matrix[i][0]);
			fprintf(file, " ");
		}
		fprinti(file, M->matrix[M->m-1][0]);
		fprintf(file, ">");
	}
}


void printbm_row(const BigIntMatrixTP M)
/** Prints a BigIntMatrixTP as a row vector. */
{
	fprintbm_row(stdout, M);
}


void fprintbm_nopad(FILE* file, const BigIntMatrixTP M)
/** Same as fprintbm(), except no padding is added to
    the beginning of each matrix element. */
{	
	//Actually print out the entries
	for (int row = 0; row < M->m; row += 1)
	{
		for (int col = 0; col < M->n; col += 1)
		{				
			fprinti(file, M->matrix[row][col]);
			fprintf(file, " ");
		}
		fprintf(file, "\n");
	}
}


void big_rowsp(const BigIntMatrixTP A)
/** Prints a representation of A's row space
    to stdout. */
{
	bool isZero;
	
	printf("span({");
	for (int row = 0; row < A->m; row += 1)
	{
		isZero = row == 0 ? FALSE : TRUE;
		
		for (int col = 0; col < A->n; col += 1)
		{
			if ((!isZero) && (col == 0))
			{
				if (row != 0)
					printf(", ");
				
				printf("<");
			}
			
			if (!isZero)
			{
				printi(big_element(A, row, col));
				if (col != A->n - 1)
					printf(" ");
				
				else
					printf(">");
			}
			
			//Check to see if the row vector is the zero vector
			if ((!is_zero(big_element(A, row, col))) && (isZero))
			{
				isZero = FALSE;
				col = -1;
			}
		}
		
	}
	printf("})");
}


void printgm(const GenericMatrixTP a)
/** Prints a generic matrix to stdout. */
{
	for (int i = 0; i < a->m; i += 1)
	{
		for (int j = 0; j < a->n; j += 1)
		{
			a->printFunction(a->matrix[i][j]);
			(j == a->n-1) ? printf("\n") : printf(", ");
		}
	}
	printf("\n");
}


void fprintgm_row(FILE* file, const GenericMatrixTP a)
/** Same as printgm_row, except it prints to a given
    file stream. */
{
	fprintf(file, "<");
	for (int i = 0; i < a->m; i += 1)
	{
		a->fprintFunction(file, a->matrix[i][0]);
		(i == a->m-1) ? fprintf(file, ">") : fprintf(file, ", ");
	}
}


void printgm_row(const GenericMatrixTP a)
/** Prints a generic matrix as a row vector. */
{
	fprintgm_row(stdout, a);
}


BigIntMatrixTP BigIntMatrixT_catalogue_get(BigIntMatrixTP** catalogue, int* catLength, const BigIntMatrixTP toFind)
/** Searches a BigIntMatrixTP array to find a reference to a matrix which is
    equivalent to the one passed in. If a suitable matrix isn't found, it's added 
		to the array (by copy) and a reference is returned. 
		
		This function assumes all the matrices in the catalogue are of the same dimensions. */
{
	//Search catalogue first to see if the matrix is in there
	for (int i = 0; i < *catLength; i += 1)
		if (compare_BigIntMatrixT(catalogue[0][i], toFind))
			return catalogue[0][i];

		
	//If we can't find toFind in the catalogue, add it
	*catLength += 1;
	catalogue[0] = realloc(catalogue[0], (*catLength)*sizeof(BigIntMatrixTP));
	catalogue[0][(*catLength)-1] = new_BigIntMatrixT(toFind->m, toFind->n);
	copy_BigIntMatrixT(toFind, catalogue[0][(*catLength)-1]);
	return catalogue[0][*catLength-1];
}


int big_mat_add(const BigIntMatrixTP A, const BigIntMatrixTP B, BigIntMatrixTP sum)
/** Computes A + B, stores sum in sum. This function assumes 
    sum has been initialised.
    Returns 1 on success, 0 otherwise. */
{
	//Ensuring dimensions of all matrices are correct
	if ((A->m != B->m) || (A->n != B->n) ||
	    (A->m != sum->m) || (A->n != sum->n))
		return 0;
		
	BigIntTP tempSum = empty_BigIntT(1);
		
	for (int row = 0; row < A->m; row += 1)
		for (int col = 0; col < A->n; col += 1)
		{
			add_BigIntT(A->matrix[row][col], B->matrix[row][col], tempSum);
			copy_BigIntT(tempSum, sum->matrix[row][col]);
		}
		
	tempSum = free_BigIntT(tempSum);
	
	return 1;
}


int mat_mul(const IntMatrixTP A, const IntMatrixTP B, IntMatrixTP result)
/** Computes AB, stores result in result. 
    This function DOES check to make sure result is the proper dimensions.
		This function assumes all three matrices have been initialised.
    Returns 1 upon success, 0 otherwise. */ 
{	
	//Checking to see if result has the appropriate dimensions
	if ((result->m != A->m) || (result->n != B->n) ||
	    (A->n != B->m))
		return 0;
	
	
	for (int mIndex = 0; mIndex < A->m; mIndex += 1) //For each row of A
		for (int rIndex = 0; rIndex < B->n; rIndex += 1) //For each column of B
		{
			//Clear out entry before we start adding multiplications
			result->matrix[mIndex][rIndex] = 0; 
			
			for (int nIndex = 0; nIndex < A->n; nIndex += 1) //Each multiplication
			{
				/*printf("m: %d, n: %d, r: %d, A[m][n]: %d, B[n][r]: %d\n",
				mIndex, nIndex, rIndex, A->matrix[mIndex][nIndex], B->matrix[nIndex][rIndex]);
				printf("Mini-mult: %d\n", A->matrix[mIndex][nIndex] * B->matrix[nIndex][rIndex]); */
				result->matrix[mIndex][rIndex] += (A->matrix[mIndex][nIndex]) * (B->matrix[nIndex][rIndex]);
			}
		}
				
	return 1;
}


int big_mat_mul(const BigIntMatrixTP A, const BigIntMatrixTP B, BigIntMatrixTP result)
/** Computes AB, stores the result in result.
    Returns 1 on success, 0 otherwise. */
{
	if ((A->n != B->m) || (result->m != A->m) || (result->n != B->n))
		return 0;
	
	BigIntTP temp  = empty_BigIntT(1);
	BigIntTP temp2 = empty_BigIntT(1);
	
	for (int row = 0; row < A->m; row += 1)
	{
		for (int col = 0; col < B->n; col += 1)
		{
			//Clear out each row before adding anything
			clear_BigIntT(result->matrix[row][col]);
			
			for (int elem = 0; elem < A->n; elem += 1)
			{
				//This is gonna be super resource heavy
				multiply_BigIntT(A->matrix[row][elem], B->matrix[elem][col], temp);
				add_BigIntT(temp, result->matrix[row][col], temp2);
				copy_BigIntT(temp2, result->matrix[row][col]);
			}
		}
	}
	
	temp  = free_BigIntT(temp);
	temp2 = free_BigIntT(temp2);
	
	return 1;
}


int gen_mat_mul(const GenericMatrixTP A, const GenericMatrixTP B, GenericMatrixTP result)
/** Computes AB, stores result in result.
    Function assumes all given matrices are properly initialised, and that all the
		elements in the matrix are compatiable.
		Returns 1 on success, 0 otherwise. */
{
	void* tempProduct; //For holding intermediary products
	
	if (clear_GenericMatrixT(result) == 0)
		return 0;
	
	//Check to make sure the dimensions of the matrices are
	// such that they can be multiplied
	if (A->n != B->m)
		return 0;
	
	if ((result->m != A->m) || (result->n != B->n))
		return 0;
	
	//Getting tempProduct ready to hold values
	//For MultiVarExtTs, copying a value into it
	// initialises the coefficients, whereas the
	// regular init function just creates the struct.
	// For other structs, the copyFunction here
	// probably isn't necessary. This is basically a
	// lazy way for me to make sure tempProduct is
	// ready.
	tempProduct = result->initFunction(result->initValue);
	result->copyFunction(A->matrix[0][0], tempProduct);
	
	//The actual matrix multiplication
	for (int row = 0; row < A->m; row += 1)
	{
		for (int col = 0; col < B->n; col += 1)
		{
			for (int elem = 0; elem < A->n; elem += 1)
			{
				result->clearFunction(tempProduct);
				result->multFunction(A->matrix[row][elem], B->matrix[elem][col], tempProduct);
				result->incFunction(tempProduct, result->matrix[row][col]);
			}
		}
	}
	
	tempProduct = result->freeFunction(tempProduct);
	
	return 1;
}


int modm(IntMatrixTP M, int mod)
/** Applies a modulus to every element of a given matrix. 
    Returns 1 on success, 0 otherwise. */
{
	for (int row = 0; row < M->m; row += 1)
		for (int col = 0; col < M->n; col += 1)
			M->matrix[row][col] %= mod;
		
	return 1;
}


int modbm(BigIntMatrixTP M, const BigIntTP modulus)
/** Reduces a BigIntMatrixTP by the given modulus.
    Returns 1 on success, 0 otherwise. */
{
	//May rewrite in the future to require a third argument
	// instead of allocating a new matrix.
	//Dependency injection !!!
	
	BigIntTP temp = empty_BigIntT(1);
	
	for (int row = 0; row < M->m; row += 1)
	{
		for (int col = 0; col < M->n; col += 1)
		{
			mod_BigIntT(M->matrix[row][col], modulus, temp);
			reduce_BigIntT(temp);
			copy_BigIntT(temp, M->matrix[row][col]);
		}
	}
	
	temp = free_BigIntT(temp);
	return 1;
}


int powbm(const BigIntMatrixTP A, 
          BigIntMatrixTP AP, 
					const BigIntTP power, 
					const BigIntTP modulus)
/** Calculates A^power, stores result in AP. Currently, this only
    works for positive powers. 
		This function assumes AP has already been initialised to the
		same dimensions as A (a square matrix).
		Returns 1 on success, 0 otherwise. */
{
	int returnVal = 0;
	
	if ((A->m != A->n) ||
	    (AP->m != A->m) ||
			(AP->n != A->n))
		return returnVal;
		
	int oneArr[1] = {1};
	BigIntTP one = new_BigIntT(oneArr, 1);
	
	BigIntTP powerCounter; //For counting how many times we've iterated the matrix
	BigIntTP tempInt;
	
	BigIntMatrixTP tempMat;
	
	//Negative powers can easily be supported with the inverse function,
	// I just don't need that ability right now
	if (compare_BigIntT(power, one) >= 0)
	{
		returnVal = 1;
		
		powerCounter = empty_BigIntT(1);
		tempInt      = empty_BigIntT(1);
		
		tempMat  = identity_BigIntMatrixT(A->m);
		
		//Iterate matrix the desired number of times
		while (compare_BigIntT(powerCounter, power) < 0)
		{
			big_mat_mul(A, tempMat, AP);
			modbm(AP, modulus);
			copy_BigIntMatrixT(AP, tempMat);
			
			//Increment counter
			add_BigIntT(powerCounter, one, tempInt);
			copy_BigIntT(tempInt, powerCounter);
		}
		
		powerCounter = free_BigIntT(powerCounter);
		tempInt      = free_BigIntT(tempInt);
		tempMat      = free_BigIntMatrixT(tempMat);
	}
	
	one = free_BigIntT(one);
	
	return returnVal;
}


int eval_BigPolyT(const BigPolyTP p, const BigIntMatrixTP A, BigIntMatrixTP result, const BigIntTP mod)
/** Plugs A into P, stores the result in result. This assumes result has been initialised.
    Returns 1 on success, 0 otherwise. */
{
	BigIntMatrixTP powMat;
	BigIntMatrixTP tempMat;
	BigIntMatrixTP I;
	
	BigIntTP* polyCoeffs; //Since we can't directly deref BigPolyTPs
	
	int oneArr[1] = {1};
	BigIntTP powCount;  //For calculating powers of A
	BigIntTP multCount; //For calculating multiples of A^k
	BigIntTP temp;
	BigIntTP one;
	
	I       = identity_BigIntMatrixT(big_rows(A));
	powMat  = new_BigIntMatrixT(big_rows(A), big_cols(A));
	tempMat = new_BigIntMatrixT(big_rows(A), big_cols(A));
	
	copy_BigIntMatrixT(tempMat, result); //Clear result before doing anything
	
	powCount  = empty_BigIntT(1);
	multCount = new_BigIntT(oneArr, 1);
	one       = new_BigIntT(oneArr, 1);
	temp      = empty_BigIntT(1);
	
	polyCoeffs = extract_coefficients(p);
	
	//Iterate over all coefficients, compute the relevant terms,
	// store them in result.
	for (int i = 0; i <= degree(p); i += 1)
	{
		//powbm() doesn't currently support ^0, so I have to improvise
		if (i > 0)
			powbm(A, powMat, powCount, mod);
		else
			copy_BigIntMatrixT(I, powMat);
	
		while (compare_BigIntT(multCount, polyCoeffs[i]) <= 0)
		{
			big_mat_add(result, powMat, tempMat);
			copy_BigIntMatrixT(tempMat, result);
			modbm(result, mod);
			
			add_BigIntT(multCount, one, temp);
			copy_BigIntT(temp, multCount);
		}
		
		//Resetting multCount, but using one instead of zero for convenience
		copy_BigIntT(one, multCount);
		
		add_BigIntT(one, powCount, temp);
		copy_BigIntT(temp, powCount);
		
		//Free coefficients as we go to save time
		polyCoeffs[i] = free_BigIntT(polyCoeffs[i]);
	}
	
	free(polyCoeffs);
	polyCoeffs = NULL;
	
	powMat  = free_BigIntMatrixT(powMat);
	tempMat = free_BigIntMatrixT(tempMat);
	I       = free_BigIntMatrixT(I);
	
	powCount  = free_BigIntT(powCount);
	multCount = free_BigIntT(multCount);
	one       = free_BigIntT(one);
	temp      = free_BigIntT(temp);
	
	return 1;
}


int eval_factored_BigPolyT(const BigPolyTP* factors, const BigIntMatrixTP A, BigIntMatrixTP result, const BigIntTP mod)
/** Same as eval_BigPolyT(), but the polynomial is given in factored form. */
{
	BigIntMatrixTP tempMat;
	BigIntMatrixTP tempMat2;
	
	BigIntTP numOfFactors;
	BigIntTP factorCounter;
	int factorCounterInt = 1;
	
	BigIntTP one;
	BigIntTP temp;
	int oneArr[1] = {1};
	
	numOfFactors = constant(factors[0]);
	factorCounter = new_BigIntT(oneArr, 1);
	temp = empty_BigIntT(1);
	one  = new_BigIntT(oneArr, 1);
	
	tempMat  = new_BigIntMatrixT(big_rows(A), big_rows(A));
	tempMat2 = new_BigIntMatrixT(big_rows(A), big_rows(A));
	
	while (compare_BigIntT(factorCounter, numOfFactors) <= 0)
	{
		clear_BigIntMatrixT(tempMat);
		eval_BigPolyT(factors[factorCounterInt], A, tempMat, mod);
		
		if (factorCounterInt == 1)
			copy_BigIntMatrixT(tempMat, result);
		else
		{
			big_mat_mul(tempMat, result, tempMat2);
			modbm(tempMat2, mod);
			copy_BigIntMatrixT(tempMat2, result);
		}
		
		factorCounterInt += 1;
		add_BigIntT(one, factorCounter, temp);
		copy_BigIntT(temp, factorCounter);
	}
	
	//numOfFactors  = free_BigIntT(numOfFactors);
	factorCounter = free_BigIntT(factorCounter);
	temp = free_BigIntT(temp);
	one  = free_BigIntT(one);
	
	tempMat  = free_BigIntMatrixT(tempMat);
	tempMat2 = free_BigIntMatrixT(tempMat2);
	return 1;
}


//It may be better to implement this non-recursively in the future for speed
//This also hasn't been extensively tested yet.
int det(IntMatrixTP M)
/** Returns the determinant of a matrix. If the
    matrix is nonsquare, returns zero. */
{
	//If matrix is nonsquare, return zero
	if (M->m != M->n)
		return 0;
	
	//If the matrix is 1x1, return the only element
	if (M->m == 1)
		return M->matrix[0][0];
	
	
	//Recursively find the determinant
	else
	{
		int sum = 0;
		IntMatrixTP tempMatrix; //Holds submatrices for determinant calculation
		
		for (int col = 0; col < M->n; col += 1)
		{
			tempMatrix = det_subIntMatrixT(M, 0, col); //Preparing submatrix
			
			//Alternate sign for terms to add
			if (col % 2 == 0)
				sum += M->matrix[0][col] * det(tempMatrix);
			else
				sum += (-1) * M->matrix[0][col] * det(tempMatrix);
			
			free_IntMatrixT(tempMatrix);
		}
		
		tempMatrix = NULL;
		return sum;
	}
}


/* private */ int row_swap(IntMatrixTP M, int x, int y)
/** Swaps the two given rows in the given matrix.
    Returns 1 on success, 0 otherwise. */
{
	//If rows given are out of range for the given matrix
	if ((x < 0) || (y < 0) ||
	    (x >= M->m) || (y >= M->m))
		return 0;
		
	//Swap row pointers
	int* tempRow = M->matrix[x];
	M->matrix[x] = M->matrix[y];
	M->matrix[y] = tempRow;

	return 1;
}


/* private */ int big_row_swap(BigIntMatrixTP M, int x, int y)
/** Swaps the two given rows in the given matrix.
    Returns 1 on success, 0 otherwise. */
{
	if (M == NULL)
		return 0;
	
	//If rows given are out of range for the given matrix
	if ((x < 0) || (y < 0) ||
	    (x >= M->m) || (y >= M->m))
		return 0;
		
	//Swap row pointers
	BigIntTP* tempRow = M->matrix[x];
	M->matrix[x] = M->matrix[y];
	M->matrix[y] = tempRow;

	return 1;
}


/* private */ int row_multiply(IntMatrixTP M, int row, int multiple, int modulus)
/** Multiplies a row by the given multiple, then takes it mod modulus.
    Returns 1 on success, 0 otherwise. */
{
	if ((row < 0) || (row >= M->m))
		return 0;
	
	//Iterate over row, multiply, take mod
	for (int i = 0; i < M->n; i += 1)
		M->matrix[row][i] = (multiple * M->matrix[row][i]) % modulus;
	
	return 1;
}


/* private */ int big_row_multiply(BigIntMatrixTP M, int row, const BigIntTP multiple, const BigIntTP modulus)
/** Multiplies a row by the given multiple, then takes it mod modulus.
    Returns 1 on success, 0 otherwise. */
{
	BigIntTP temp;
	
	if (M == NULL)
		return 0;
	
	if ((row < 0) || (row >= M->m))
		return 0;
	
	temp = empty_BigIntT(1);
	
	//Iterate over row, multiply, take mod
	for (int i = 0; i < M->n; i += 1)
	{
		multiply_BigIntT(M->matrix[row][i], multiple, temp);
		mod_BigIntT(temp, modulus, M->matrix[row][i]);
	}
	
	temp = free_BigIntT(temp);
	
	return 1;
}


/* private */ int row_add(IntMatrixTP M, int addTo, int addFrom, int modulus)
/** Adds row addFrom to row addTo in the given matrix.
    Returns 1 on success, 0 otherwise. */
{
	if ((addTo < 0) || (addFrom < 0) ||
	    (addTo >= M->m) || (addFrom >= M->m))
		return 0;
		
	//Iterate over row, add relevant entries, take modulus
	for (int i = 0; i < M->n; i += 1)
	{
		M->matrix[addTo][i] += (M->matrix[addFrom][i]);
		M->matrix[addTo][i] %= modulus;
	}
	
	return 1;
}	


/* private */ int big_row_add(BigIntMatrixTP M, int addTo, int addFrom, const BigIntTP modulus)
/** Adds row addFrom to row addTo in the given matrix.
    Returns 1 on success, 0 otherwise. */
{
	BigIntTP temp;
	
	if (M == NULL)
		return 0;
	
	if ((addTo < 0) || (addFrom < 0) ||
	    (addTo >= M->m) || (addFrom >= M->m))
		return 0;
		
	temp = empty_BigIntT(1);
		
	//Iterate over row, add relevant entries, take modulus
	for (int i = 0; i < M->n; i += 1)
	{
		add_BigIntT(M->matrix[addTo][i], M->matrix[addFrom][i], temp);
		mod_BigIntT(temp, modulus, M->matrix[addTo][i]);
	}
	
	temp = free_BigIntT(temp);
	
	return 1;
}	


IntMatrixTP inverse(const IntMatrixTP M, int modulus)
/** If the inverse of M exists, a pointer to M's inverse
    is returned. Otherwise, returns NULL. 
		
		This function assumes the given matrix conforms to the 
		given modulus (all entries are less than the modulus). */
{
	//If we're given a non-square matrix
	if ((M->m != M->n))
		return NULL;
	
	bool hasLeadEntry;
	
	int numTimesToAdd; //Tells how many times we add a row to another row
	int tempInverse;
	
	IntMatrixTP toReduce = new_IntMatrixT(M->m, M->m);
	IntMatrixTP inv      = identity_IntMatrixT(M->m);
	copy_IntMatrixT(M, toReduce);
	
	#ifdef VERBOSE
	printf("Matrix to reduce:\n");
	printm(toReduce);
	#endif
	
	//Converting toReduce to upper triangular (row echelon) form
	for (int focusRow = 0; focusRow < M->m; focusRow += 1)
	{
		hasLeadEntry = FALSE;
		
		//Find a row with a non-zero first entry
		for (int nonzero = focusRow; nonzero < M->m; nonzero += 1)
		{
			//Checking to see whether the leading entry is one
			//OR
			//Checking to see whether the leading entry is nonzero and has an inverse
			if ((toReduce->matrix[nonzero][focusRow] == 1) ||
					((toReduce->matrix[nonzero][focusRow] != 0) &&
					 (GCD(toReduce->matrix[nonzero][focusRow], modulus) == 1)))
			{
				hasLeadEntry = TRUE;
				
				//Don't need to swap if the row is already in place
				if (nonzero == focusRow)
				{
					#ifdef VERBOSE
					printf("Nonzero leading entry is already in place; no swap necessary.\n");
					#endif
					break;
				}
				
				row_swap(toReduce, focusRow, nonzero);
				row_swap(inv, focusRow, nonzero);
				
				#ifdef VERBOSE
				printf("Swapped row %d and %d.\n", focusRow, nonzero);
				printm(toReduce);
				#endif

				break;
			}
		}
		
		//If we couldn't find a nonzero entry for our pivot column
		if (!hasLeadEntry)
		{
			#ifdef VERBOSE
			printf("No nonzero, invertible leading entry could be found.\n");
			printf("Attempting to add rows instead.\n");
			#endif
			
			//Check to see if we can add rows to get an invertible entry
			for (int rowToMaybeAdd = focusRow+1; rowToMaybeAdd < M->m; rowToMaybeAdd += 1)
			{
				//If we found a row that we can add to our current row to get an invertible number
				for (int mult = 1; mult < modulus; mult += 1)
				{
					tempInverse = num_inverse(mult*M->matrix[rowToMaybeAdd][focusRow] + M->matrix[focusRow][focusRow], modulus);
					if (tempInverse != -1)
					{
						//Add rows to get an invertible number as a leading entry
						for (int t = 0; t < mult; t += 1)
						{
							row_add(toReduce, focusRow, rowToMaybeAdd, modulus);
							row_add(inv, focusRow, rowToMaybeAdd, modulus);
						}
						
						#ifdef VERBOSE
						printf("Added row %d to row %d a total of %d times.\n", rowToMaybeAdd, focusRow, mult);
						printm(toReduce);
						#endif
						
						//Now, get the leading entry to one
						row_multiply(toReduce, focusRow, tempInverse, modulus);
						row_multiply(inv, focusRow, tempInverse, modulus);
						
						#ifdef VERBOSE
						printf("Multipled row %d by %d.\n", focusRow, tempInverse);
						printm(toReduce);
						#endif
						
						hasLeadEntry = TRUE;
						break;
					}
				}
				
				if (hasLeadEntry)
					break;
			}
		}
		
		//If we couldn't find a suitable row to create an invertible element
		if (!hasLeadEntry)
		{
			toReduce = free_IntMatrixT(toReduce);
			inv = free_IntMatrixT(inv);
			
			#ifdef VERBOSE
			printf("No suitable row was found to create an invertible entry.\n");
			#endif
			
			return NULL;
		}

		//If the leading entry is a 1, we don't need to find an inverse
		if (toReduce->matrix[focusRow][focusRow] != 1)
		{
			tempInverse = num_inverse(toReduce->matrix[focusRow][focusRow], modulus);
			row_multiply(toReduce, focusRow, tempInverse, modulus);
			row_multiply(inv, focusRow, tempInverse, modulus);
			
			#ifdef VERBOSE
			printf("Multiplied row %d by %d.\n", focusRow, tempInverse);
			printm(toReduce);
			#endif
		}
		
		//Now, we clear out all nonzero entries in our current pivot column
		for (int i = focusRow+1; i < M->m; i += 1)
		{
			numTimesToAdd = (modulus - toReduce->matrix[i][focusRow]) % modulus;
			for (int v = 0; v < numTimesToAdd; v += 1)
			{
				row_add(toReduce, i, focusRow, modulus);
				row_add(inv, i, focusRow, modulus);
			}
			
			#ifdef VERBOSE
			if (numTimesToAdd > 0)
			{
				printf("Added row %d to row %d a total of %d times.\n", focusRow, i, numTimesToAdd);
				printm(toReduce);
			}
			#endif
		}
	}
	
	#ifdef VERBOSE
	printf("The matrix should now be in row echelon form.\n\n");
	#endif
	
	//Converting toReduce to reduced row echelon form (the identity)
	for (int i = M->m-1; i >= 0; i -= 1)
	{
		//Each element above our leading entry
		for (int element = i-1; element >= 0; element -= 1)
		{
			numTimesToAdd = (modulus - toReduce->matrix[element][i]) % modulus;
			for (int v = 0; v < numTimesToAdd; v += 1)
			{
				row_add(toReduce, element, i, modulus);
				row_add(inv, element, i, modulus);
			}
			
			#ifdef VERBOSE
			if (numTimesToAdd > 0)
			{
				printf("Added row %d to row %d a total of %d times.\n", i, element, numTimesToAdd);
				printm(toReduce);
			}
			#endif
		}
	}
	
	toReduce = free_IntMatrixT(toReduce);
	return inv;
}


int big_row_echelon(const BigIntMatrixTP M, 
                    const BigIntTP modulus, 
										BigIntMatrixTP result,
										BigIntMatrixTP aux)
/** Finds an upper triangular form of M, if it exists, and
    stores it in result. result must be the same dimensions as M.
		The matrix aux has the same operations done to it as M.
		Returns 1 upon finding an upper triangular form,
		0 otherwise. */
{
	bool hasLeadEntry;
	bool focusHasInverse;
	
	//For when columns of zeros appear, place pivot on correct row
	int pivotColOffset = 0;
	
	//Row with least non-invertible element when reducing rows with non-invertible elements
	int leastNonInvIndex = 0;
	
	//For keeping track of how non invertible elements are
	BigIntTP leastNonInvCount;
	BigIntTP leastNonInvTempCount;
	
	int numArr[1] = {1};
	BigIntTP zero, one, bigMult, temp, temp2, tempInverse;
	BigIntTP tempGCD = NULL;
	
	//Making sure our matrices are of the correct shape
	if ((result->m != M->m) || (result->n != M->n))
		return 0;
	
	if (aux != NULL)
		if ((aux->m != M->m))
			return 0;
	
	copy_BigIntMatrixT(M, result);
	modbm(result, modulus); //Making sure no wacky entries throw our algorithm off
	
	zero = empty_BigIntT(1);
	one  = new_BigIntT(numArr, 1);
	
	temp        = empty_BigIntT(1);
	temp2       = empty_BigIntT(1);
	bigMult     = empty_BigIntT(1);
	tempInverse = empty_BigIntT(1);
	
	//For determining which numbers are least non invertible 
	// when reducing rows with non-invertible elements
	leastNonInvTempCount = empty_BigIntT(1);
	leastNonInvCount     = empty_BigIntT(1);
	
	for (int focusRow = 0; focusRow < M->m; focusRow += 1)
	{
		if (focusRow+pivotColOffset >= M->n)
			break;
		
		#ifdef VERBOSE
		printf("focusRow = %d, pivotColOffset = %d\n", focusRow, pivotColOffset);
		#endif
		
		//Make sure each row has to find an invertible leading entry
		hasLeadEntry = FALSE;
		
		//Find a row with a non-zero first entry
		for (int nonzero = focusRow; nonzero < M->m; nonzero += 1)
		{
			//Checking to see whether the leading entry is one
			//OR
			//Checking to see whether the leading entry is nonzero and has an inverse
			tempGCD = big_gcd(result->matrix[nonzero][focusRow+pivotColOffset], modulus);
			
			if (compare_BigIntT(tempGCD, one) == 0)
			{
				hasLeadEntry = TRUE;
				
				//Don't need to swap if the row is already in place
				if (nonzero == focusRow)
				{
					#ifdef VERBOSE
					printf("Nonzero leading entry is already in place; no swap necessary.\n");
					#endif
					break;
				}
				
				big_row_swap(result, focusRow, nonzero);
				big_row_swap(aux, focusRow, nonzero);
				
				#ifdef VERBOSE
				printf("Swapped row %d and %d.\n", focusRow, nonzero);
				printbm(result);
				#endif

				break;
			}
			tempGCD = free_BigIntT(tempGCD);
		}
		tempGCD = free_BigIntT(tempGCD);
		
		//If we couldn't find a nonzero entry for our pivot column
		if (!hasLeadEntry)
		{
			#ifdef VERBOSE
			printf("No nonzero, invertible leading entry could be found.\n");
			#endif
		}
		
		//If we couldn't find a suitable row to create an invertible element
		if (!hasLeadEntry)
		{
			#ifdef VERBOSE
			printf("Settling for the \"least\" non-invertible element.\n");
			#endif
			
			//Finding an element that can eliminate all elements below it
			//THIS ONLY WORKS WITH PRIME-POWER MODULI
			leastNonInvIndex = 0;
			copy_BigIntT(modulus, leastNonInvCount);
			hasLeadEntry = FALSE; //Have we found a nonzero element?
			for (int leastNonInv = focusRow; leastNonInv < M->m; leastNonInv += 1)
			{
				if (compare_BigIntT(result->matrix[leastNonInv][focusRow+pivotColOffset], zero) != 0)
				{
					hasLeadEntry = TRUE;
					free_BigIntT(leastNonInvTempCount);
					
					leastNonInvTempCount = big_gcd(result->matrix[leastNonInv][focusRow+pivotColOffset], modulus);
					
					//leastNonInvTempCount now holds the largest power of basePrime in result->matrix[leastNonInv][focusRow+pivotColOffset]
					if (compare_BigIntT(leastNonInvTempCount, leastNonInvCount) < 0)
					{
						copy_BigIntT(leastNonInvTempCount, leastNonInvCount);
						leastNonInvIndex = leastNonInv;
					}
				}
			}
			
			if (!hasLeadEntry) //If all elements were zero, shift our focus column over
			{
				focusRow -= 1;
				pivotColOffset += 1;
				
				#ifdef VERBOSE
				printf("Current pivot column is completely zero. Offsetting column by 1.\n");
				#endif
				
				continue;
			}
			
			//Swap our found row to the focus row
			else
			{
				//Get the simplest leading non-invertible element possible
				divide_BigIntT(result->matrix[leastNonInvIndex][focusRow+pivotColOffset], leastNonInvCount, temp);
				divide_BigIntT(modulus, leastNonInvCount, temp2);
				big_num_inverse(temp, temp2, tempInverse);
				
				//Using embedded dynamics to find "inverse"
				//multiply_BigIntT(tempInverse, leastNonInvCount, temp);
				if (compare_BigIntT(temp, one) != 0)
				{
					big_row_multiply(result, leastNonInvIndex, tempInverse, modulus);
					big_row_multiply(aux, leastNonInvIndex, tempInverse, modulus);
					
					#ifdef VERBOSE
					printf("Multiplied row %d by ", leastNonInvIndex);
					printi(temp);
					printf(".\n");
					printbm(result);
					#endif
				}
				
				if (focusRow != leastNonInvIndex)
				{
					big_row_swap(result, focusRow, leastNonInvIndex);
					big_row_swap(aux, focusRow, leastNonInvIndex);
					
					#ifdef VERBOSE
					printf("Swapped row %d with %d.\n", focusRow, leastNonInvIndex);
					printbm(result);
					#endif
				}
				
				//Tell next conditional we don't need to do any work
				focusHasInverse = FALSE;
			}
		}
		
		else //If we're dealing with the normal case where the focus element has an inverse
			focusHasInverse = TRUE;

		//If the leading entry is a 1 or is non invertible, we don't need to find an inverse
		if ((compare_BigIntT(result->matrix[focusRow][focusRow+pivotColOffset], one) != 0) &&
		    (focusHasInverse))
		{
			big_num_inverse(result->matrix[focusRow][focusRow+pivotColOffset], modulus, tempInverse);
			if (compare_BigIntT(tempInverse, one) != 0)
			{
				big_row_multiply(result, focusRow, tempInverse, modulus);
				big_row_multiply(aux, focusRow, tempInverse, modulus);
				
				#ifdef VERBOSE
				printf("Multiplied row %d by ", focusRow);
				printi(tempInverse);
				printf(".\n");
				printbm(result);
				#endif
			}
		}
		
		//Now, we clear out all nonzero entries in our current pivot column
		for (int i = focusRow+1; i < M->m; i += 1)
		{
			subtract_BigIntT(modulus, result->matrix[i][focusRow+pivotColOffset], temp);
			mod_BigIntT(temp, modulus, temp2); //temp2 is acting as "numTimesToAdd" here
			
			//Added so that non invertible elements don't need to add more times than necessary
			if (!focusHasInverse)
			{
				divide_BigIntT(temp2, leastNonInvCount, temp);
				copy_BigIntT(temp, temp2);
			}
			
			if (compare_BigIntT(temp2, zero) != 0)
			{
				#ifdef VERBOSE
					printf("Added row %d to row %d a total of ", focusRow, i);
					printi(temp2);
					compare_BigIntT(one, temp2) == 0 ? printf(" time.\n") : printf(" times.\n");
				#endif
				
				while (compare_BigIntT(temp2, zero) != 0)
				{
					big_row_add(result, i, focusRow, modulus);
					big_row_add(aux, i, focusRow, modulus);
					
					//Decrement temp2
					subtract_BigIntT(temp2, one, temp);
					copy_BigIntT(temp, temp2);
				}
				
				#ifdef VERBOSE
				printbm(result);
				#endif
			}
		}
	}
	
	zero = free_BigIntT(zero);
	one  = free_BigIntT(one);
	
	temp        = free_BigIntT(temp);
	temp2       = free_BigIntT(temp2);
	bigMult     = free_BigIntT(bigMult);
	tempInverse = free_BigIntT(tempInverse);
	
	leastNonInvCount     = free_BigIntT(leastNonInvCount);
	leastNonInvTempCount = free_BigIntT(leastNonInvTempCount);
	
	return 1;
}


int big_reduced_row_echelon(const BigIntMatrixTP M, 
                            const BigIntTP modulus, 
														BigIntMatrixTP result, 
														BigIntMatrixTP aux)
/** Takes M, in row echelon form, and tries to reduce it to 
    reduced row echelon form. Stores the result in result.
		The operations performed on M will also be applied to aux.
		Set aux to NULL if this behaviour isn't needed.
		Returns 1 upon reducing M to the identity (or some slice of it), 
		0 otherwise. */
{
	int perfectReducedRowEchelon = 1;
	int leadingEntryCol = 0;
	BigIntTP temp, temp2, zero, one;
	
	if ((M->m != result->m) ||
	    (M->n != result->n))
		return 0;
		
	if (aux != NULL)
		if ((M->m != aux->m) ||
	      (M->n != aux->n))
			return 0;

	copy_BigIntMatrixT(M, result);
	
	int numArr[1] = {1};
	temp  = empty_BigIntT(1);
	temp2 = empty_BigIntT(1);
	zero  = empty_BigIntT(1);
	one   = new_BigIntT(numArr, 1);

	for (int i = M->m-1; i >= 0; i -= 1)
	{
		//We can't assume the leading entry is on the diagonal
		//Let's find it!
		leadingEntryCol = -1;
		for (int j = 0; j < M->n; j += 1)
			if (compare_BigIntT(result->matrix[i][j], zero) != 0)
			{
				leadingEntryCol = j;
				break;
			}
		
		//Each element above our leading entry
		if (leadingEntryCol != -1)
		{
			for (int element = i-1; element >= 0; element -= 1)
			{
				subtract_BigIntT(modulus, result->matrix[element][leadingEntryCol], temp);
				mod_BigIntT(temp, modulus, temp2); //temp2 is acting as "numTimesToAdd" here
				
				//Adjusting for non-one leading terms
				//This allows terms to get as small as possible
				divide_BigIntT(temp2, result->matrix[i][leadingEntryCol], temp);
				
				//Some manual adjustment for style
				if ((compare_BigIntT(temp2, zero) != 0) && (compare_BigIntT(temp, zero) == 0))
					copy_BigIntT(one, temp2);
				else
					mod_BigIntT(temp, modulus, temp2);
				
				if ((perfectReducedRowEchelon == 1) && 
				    (compare_BigIntT(result->matrix[i][leadingEntryCol], one) != 0))
					perfectReducedRowEchelon = 0;
				
				if (compare_BigIntT(temp2, zero) != 0)
				{
					#ifdef VERBOSE
						printf("Added row %d to row %d a total of ", i, element);
						printi(temp2);
						compare_BigIntT(temp2, one) == 0 ? printf(" time.\n") : printf(" times.\n");
					#endif

					//Why did I do it this way
					while (compare_BigIntT(temp2, zero) != 0)
					{
						big_row_add(result, element, i, modulus);
						big_row_add(aux, element, i, modulus);
						
						//Decrement temp2
						subtract_BigIntT(temp2, one, temp);
						copy_BigIntT(temp, temp2);
					}
					
					#ifdef VERBOSE
					printbm(result);
					#endif
				}
			}
		}
		
		else
			perfectReducedRowEchelon = 0;
	}
	
	#ifdef VERBOSE
	printf("The matrix is now as reduced as possible.\n");
	#endif
	
	one   = free_BigIntT(one);
	temp  = free_BigIntT(temp);
	zero  = free_BigIntT(zero);
	temp2 = free_BigIntT(temp2);
	
	return perfectReducedRowEchelon;
}


int big_eliminate_bottom(BigIntMatrixTP A, const BigIntTP bigMod, BigPolyTP annih)
/** Attempts to eliminate the bottom row of the given matrix.
    This function assumes the matrix is in row echelon form.
		If annih != NULL, then whatever row operations annihilate the bottom
		row will be recorded as a polynomial and stored in annih.
		Returns 1 upon successfully eliminating the row, 0 otherwise. */
{
	if ((A == NULL) || (bigMod == NULL))
		return 0;
	
	int oneArr[1] = {1};
	BigIntTP  zero      = empty_BigIntT(1);
	BigIntTP  one       = NULL;
	BigIntTP  temp      = NULL;
	BigIntTP  GCDentry  = NULL;
	BigIntTP  GCDbottom = NULL;
	
	BigIntTP* annihCoeffs; //Holds a poly representation of the row ops we do
	BigPolyTP annihPoly;
	
	if (annih != NULL)
	{
		annihCoeffs = malloc((big_rows(A)-1)*sizeof(BigIntTP));
		for (int i = 0; i < big_rows(A)-1; i += 1)
			annihCoeffs[i] = empty_BigIntT(1);
		
		one  = new_BigIntT(oneArr, 1);
		temp = empty_BigIntT(1);
	}
	
	//Iterate over each row
	for (int r = 0; r < A->m-1; r += 1)
	{
		//Find the leading entry on each row
		for (int entry = 0; entry < A->n; entry += 1)
		{
			if (compare_BigIntT(zero, A->matrix[r][entry]) != 0)
			{
				//Now, check to see if we can eliminate the corresponding element
				// in the bottom row
				GCDbottom = big_gcd(bigMod, A->matrix[A->m-1][entry]);
				GCDentry  = big_gcd(bigMod, A->matrix[r][entry]);
				
				//Eliminate corresponding entry
				//This only works for prime and prime-power moduli
				if (compare_BigIntT(GCDentry, GCDbottom) <= 0)
				{
					while (compare_BigIntT(zero, A->matrix[A->m-1][entry]) != 0)
					{
						big_row_add(A, A->m-1, r, bigMod);
						
						//Store row operations as polynomial
						if (annih != NULL)
						{
							add_BigIntT(annihCoeffs[r], one, temp);
							mod_BigIntT(temp, bigMod, annihCoeffs[r]);
						}
					}
				}
				
				else
				{
					zero      = free_BigIntT(zero);
					GCDentry  = free_BigIntT(GCDentry);
					GCDbottom = free_BigIntT(GCDbottom);
					temp      = free_BigIntT(temp);
					one       = free_BigIntT(one);
					
					if (annih != NULL)
					{
						for (int i = 0; i < big_rows(A)-1; i += 1)
							annihCoeffs[i] = free_BigIntT(annihCoeffs[i]);
						free(annihCoeffs);
					}
					return 0;
				}
				
				GCDentry  = free_BigIntT(GCDentry);
				GCDbottom = free_BigIntT(GCDbottom);
				break;
			}
		}
	}
	
	//Now, we check to see if the bottom row has been eliminated
	for (int c = 0; c < A->n; c += 1)
	{
		if (compare_BigIntT(A->matrix[A->m-1][c], zero) != 0)
		{
			zero = free_BigIntT(zero);
			temp = free_BigIntT(temp);
			one  = free_BigIntT(one);
			
			if (annih != NULL)
			{
				for (int i = 0; i < big_rows(A)-1; i += 1)
					annihCoeffs[i] = free_BigIntT(annihCoeffs[i]);
				free(annihCoeffs);
			}
			return 0;
		}
	}
	
	//Now, let's create our polynomial representation of the bottom row
	if (annih != NULL)
	{
		annihPoly = new_BigPolyT(annihCoeffs, big_rows(A)-1);
		copy_BigPolyT(annihPoly, annih);
		
		annihPoly = free_BigPolyT(annihPoly);
		for (int i = 0; i < big_rows(A)-1; i += 1)
			annihCoeffs[i] = free_BigIntT(annihCoeffs[i]);
		free(annihCoeffs);
	}
	
	zero = free_BigIntT(zero);
	temp = free_BigIntT(temp);
	one  = free_BigIntT(one);
	return 1;
}


BigIntMatrixTP big_inverse(const BigIntMatrixTP M, const BigIntTP modulus)
/** If the inverse of M exists, a pointer to M's inverse
    is returned. Otherwise, returns NULL. 
		
		This function assumes the given matrix conforms to the 
		given modulus (all entries are less than the modulus). 
		
		This function is the same as inverse(), but for BigIntMatrixTs. */
{
	//If we're given a non-square matrix
	if ((M->m != M->n))
		return NULL;
	
	BigIntMatrixTP toReduce  = new_BigIntMatrixT(M->m, M->m);
	BigIntMatrixTP toReduce2 = new_BigIntMatrixT(M->m, M->m);
	BigIntMatrixTP inv       = identity_BigIntMatrixT(M->m);
	copy_BigIntMatrixT(M, toReduce);
	
	#ifdef VERBOSE
	printf("Matrix to reduce:\n");
	printbm(toReduce);
	#endif
	
	//Converting toReduce to upper triangular (row echelon) form
	big_row_echelon(M, modulus, toReduce, inv);
	
	#ifdef VERBOSE
	printf("The matrix should now be in row echelon form.\n\n");
	#endif
	
	//Converting toReduce to reduced row echelon form (the identity)
	if (big_reduced_row_echelon(toReduce, modulus, toReduce2, inv) == 0)
		inv = free_BigIntMatrixT(inv);
	
	toReduce  = free_BigIntMatrixT(toReduce);
	toReduce2 = free_BigIntMatrixT(toReduce2);
	return inv;
}


int rref(IntMatrixTP M, int mod)
/** Row reduces a matrix to the identity, or as close to
    the identity as possible if the matrix is nonsquare.
		Returns 1 if the reduction is successful, 0 otherwise. */
{
	//I should really do proper error checking in these functions, but oh well
	//One day we will
	
	bool foundLeadEntry;
	int minDim = (M->m < M->n) ? M->m : M->n;
	int addFactor;
	
	//Says which row the lead entry should be put on.
	//This is required as the lead entry isn't necessarily on a diagonal
	int rowToSwap = 0;
	
	//Iterate over each pivot column possible
	for (int pivotCol = 0; pivotCol < minDim; pivotCol += 1)
	{
		//Search for suitable leading entry
		foundLeadEntry = FALSE;
		for (int leadEntry = rowToSwap; leadEntry < M->m; leadEntry += 1)
		{
			if (num_inverse(M->matrix[leadEntry][pivotCol], mod) != -1)
			{
				foundLeadEntry = TRUE;
				row_swap(M, rowToSwap, leadEntry);
				break;
			}
		}
		
		if (foundLeadEntry)
		{
			//Now, we use our leading entry to zero out all other entries below it
			row_multiply(M, rowToSwap, num_inverse(M->matrix[rowToSwap][pivotCol], mod), mod);
			rowToSwap += 1;
			
			for (int row = rowToSwap; row < M->m; row += 1)
			{
				//Prevents the value from changing
				addFactor = (mod - M->matrix[row][pivotCol]) % mod;
				for (int t = 0; t < addFactor; t += 1)
					row_add(M, row, rowToSwap-1, mod);
			}
		}
	}
	
	printf("REF:\n");
	printm(M);
	
	//Now we zero out elements above the leading elements
	for (int row = 0; row < M->m; row += 1)
	{
		//Find the leading entry
		for (int leadEntry = 0; leadEntry < M->n; leadEntry += 1)
		{
			if (M->matrix[row][leadEntry] == 1)
			{				
				//We found the leading entry, now zero out the elements above it
				for (int toZeroRow = row-1; toZeroRow >= 0; toZeroRow -= 1)
				{
					addFactor = (mod - M->matrix[toZeroRow][leadEntry]) % mod;
					for (int t = 0; t < addFactor; t += 1)
						row_add(M, toZeroRow, row, mod);
				}
				break;
			}
		}
	}
	
	return 1;
}


/* private */ int chara_poly_sub_matrix(BigPolyTP** const A, int size, int x, int y, BigPolyTP** new)
/** Calculates the submatrices required for chara_poly_recurse(), stores
    result in new. This function assumes new is the correct size.
		Returns 1 on success, 0 otherwise. */
{
	//These two variables facilitate taking the correct rows and columns
	// out of A
	int xSkip = 0;
	int ySkip;
	
	for (int row = 0; row < size; row += 1)
	{
		//Preventing indexing errors
		ySkip = 0;
		
		if (row == x)
		{
			xSkip = -1;
			continue;
		}
		
		for (int col = 0; col < size; col += 1)
		{
			if (col == y)
			{
				ySkip = -1;
				continue;
			}
			
			//We can get away with just copying addresses here
			new[row + xSkip][col + ySkip] = A[row][col];
		}
	}
	
	return 1;
}


/* private */ BigPolyTP chara_poly_recurse(BigPolyTP** const A, BigIntTP mod, int size)
/** Helper function for chara_poly(). Returns the characteristic
    equation for a sub matrix. */
{
	//It might be nice to rewrite this function to pass
	// along constant data instead of recalculating it each time...
	
	BigPolyTP result = empty_BigPolyT();
	BigPolyTP temp   = empty_BigPolyT();
	BigPolyTP temp2  = empty_BigPolyT();
	
	BigPolyTP subResult = NULL;
	
	BigPolyTP** subArr; //Holds a smaller array for calculating sub determinants
	
	int oneArr[] = {1};
	BigIntTP one = new_BigIntT(oneArr, 1);
	
	BigIntTP minusOne = empty_BigIntT(1);
	subtract_BigIntT(mod, one, minusOne);
	
	BigPolyTP minusOneConst = constant_BigPolyT(minusOne);
	
	//Base case
	if (size == 2)
	{
		//2x2 determinant, basically
		multiply_BigPolyT(A[0][1], A[1][0], temp2);
		multiply_BigPolyT(temp2, minusOneConst, temp);
		multiply_BigPolyT(A[0][0], A[1][1], temp2);
		add_BigPolyT(temp, temp2, result);
	}
	
	//Recursive case
	else
	{
		//Allocate space for our sub array
		subArr = malloc((size-1)*sizeof(BigPolyTP*));
		for (int i = 0; i < size-1; i += 1)
			subArr[i] = malloc((size-1)*sizeof(BigPolyTP));
		
		for (int i = 0; i < size; i += 1)
		{
			chara_poly_sub_matrix(A, size, 0, i, subArr);
			subResult = chara_poly_recurse(subArr, mod, size-1);
			
			multiply_BigPolyT(A[0][i], subResult, temp);
			
			//Managing plus and minus signs accordingly
			if (i % 2 == 0)
			{
				add_BigPolyT(temp, result, temp2);
				copy_BigPolyT(temp2, result);
			}
			
			else
			{
				multiply_BigPolyT(temp, minusOneConst, temp2);
				add_BigPolyT(temp2, result, temp);
				copy_BigPolyT(temp, result);
			}
			
			subResult = free_BigPolyT(subResult);
		}
		
		//We only need to deallocate subArr itself
		//The BigPolyTPs inside are just references
		//No need to free an array that doesn't exist
		for (int j = 0; j < size-1; j += 1)
		{
			free(subArr[j]);
			subArr[j] = NULL;
		}
		free(subArr);
		subArr = NULL;
	}
	
	//Reduce result by modulus
	mod_BigPolyT(result, mod, temp);
	copy_BigPolyT(temp, result);
	
	temp = free_BigPolyT(temp);
	temp2 = free_BigPolyT(temp2);
	subResult = free_BigPolyT(subResult);
	minusOneConst = free_BigPolyT(minusOneConst);
	
	one = free_BigIntT(one);
	minusOne = free_BigIntT(minusOne);
	
	return result;
}


//det_subIntMatrixT(IntMatrixTP const M, int x, int y)
BigPolyTP chara_poly(const BigIntMatrixTP A, const BigIntTP mod)
/** Spits a matrix's characteristic equation to the
    screen mod some modulus. A must be a square matrix. 
		Returns the characteristic equation of the passed
		matrix. Returns NULL on error. */
{
	BigPolyTP chara = NULL;
	BigPolyTP** algMat;
	
	//For creating the BigPolyTP 2D array below
	BigIntTP* term;
	
	BigIntTP one;
	BigIntTP minusOne;
	int oneArr[] = {1};
	
	if (A->m != A->n)
		return NULL;
	
	//Special case, will program later
	if (A->m == 1)
		return NULL;
	
	algMat = malloc((A->m)*sizeof(BigPolyTP*));
	term = malloc(2*sizeof(BigIntTP));
	term[0] = empty_BigIntT(1);
	term[1] = empty_BigIntT(1);
	
	one = new_BigIntT(oneArr, 1);
	minusOne = empty_BigIntT(1);
	
	subtract_BigIntT(mod, one, minusOne);
	
	//Creating a 2D array of polynomials
	//Super, SUPER inefficient but I can't be bothered
	// to write anything better.
	for (int row = 0; row < A->m; row += 1)
	{
		algMat[row] = malloc((A->m)*sizeof(BigPolyTP));
		for (int col = 0; col < A->m; col += 1)
		{
			//Diagonal entry, element minus lambda
			if (row == col)
			{
				//multiply_BigIntT(A->matrix[row][col], minusOne, term[0]);
				copy_BigIntT(A->matrix[row][col], term[0]);
				copy_BigIntT(minusOne, term[1]);
				algMat[row][col] = new_BigPolyT(term, 2);
			}
			
			//Just copy the number from the matrix if not a diagonal
			else
				algMat[row][col] = new_BigPolyT(&(A->matrix[row][col]), 1);
		}
	}
	
	//Ok, so we have our algebraic matrix
	//Now, use it to get our characteristic equation
	chara = chara_poly_recurse(algMat, mod, A->m);
	
	//Now, free the horrible matrix we created
	for (int row = 0; row < A->m; row += 1)
	{
		for (int col = 0; col < A->m; col += 1)
			algMat[row][col] = free_BigPolyT(algMat[row][col]);
		
		free(algMat[row]);
		algMat[row] = NULL;
	}
	free(algMat);
	algMat = NULL;
	
	term[0] = free_BigIntT(term[0]);
	term[1] = free_BigIntT(term[1]);
	free(term);
	term = NULL;
	
	one = free_BigIntT(one);
	minusOne = free_BigIntT(minusOne);
	
	return chara;
}


BigPolyTP* min_poly(const BigIntMatrixTP A, 
                    const BigIntMatrixTP v, 
									  const BigIntTP mod,
									  const BigPolyTP somePoly)
/** Calculates A's minimum polynomial, using a given charaPoly,
		and returns it in its factored form on success, NULL otherwise. 
		This won't be tested for composite moduli. */
{	
	BigPolyTP* tempMinPoly = NULL;
	BigPolyTP* tempPoly = NULL;
	BigPolyTP  tempCalcPoly = NULL; //For performing calculations
	BigFactorsTP factoredCharaPoly = NULL;
	BigPolyTP* unneededFactors = NULL;
	BigPolyTP  charaPoly = NULL;
	//BigFactorsTP factoredMinPoly = NULL;
	int unneededFactorsCount = 0;
	
	BigPolyTP onePoly    = NULL;
	BigPolyTP negOnePoly = NULL;
	
	BigIntTP numOfFactors = NULL;
	BigIntTP tempCounter = NULL;
	BigIntTP one = NULL;
	BigIntTP zero = NULL;
	BigIntTP temp = NULL;
	int numOfFactorsInt = 0;
	
	int oneArr[1] = {1};
	int factorToExclude = 0;
	
	BigIntMatrixTP tempMat = NULL;
	BigIntMatrixTP zeroMat = NULL;
	
	//Used for checking whether a matrix annihilates a vector
	BigIntMatrixTP tempVect = NULL;
	BigIntMatrixTP zeroVect = NULL;
	if (v != NULL)
	{
		zeroVect = new_BigIntMatrixT(big_rows(v), 1);
		tempVect = new_BigIntMatrixT(big_rows(v), 1);
	}
	
	bool allFactorsAreNeeded = FALSE;
	int indexCorrection;
	
	//Calculate charaPoly only if needed
	if (somePoly == NULL)
		charaPoly = chara_poly(A, mod);
	else
		charaPoly = somePoly;

	if (charaPoly != NULL)
	{
		factoredCharaPoly = factor_BigPolyT(charaPoly, mod);
		numOfFactorsInt = count_factors(factoredCharaPoly);
		
		one  = new_BigIntT(oneArr, 1);
		zero = empty_BigIntT(1);
		temp = empty_BigIntT(1);
		
		onePoly = constant_BigPolyT(one);
		
		subtract_BigIntT(mod, one, temp);
		negOnePoly = constant_BigPolyT(temp);
		tempCalcPoly = empty_BigPolyT();
		
		tempCounter = new_BigIntT(oneArr, 1);
		numOfFactors = empty_BigIntT(1);
		
		//Get BigIntT version of numOfFactors
		for (int i = 0; i < numOfFactorsInt; i += 1)
		{
			add_BigIntT(one, numOfFactors, temp);
			copy_BigIntT(temp, numOfFactors);
		}
		
		//Making sure tempMinPoly looks like an old BigPolyT_factors array
		tempMinPoly = extract_factors(factoredCharaPoly);
		tempMinPoly = realloc(tempMinPoly, (numOfFactorsInt+1)*sizeof(BigPolyTP));
		for (int i = numOfFactorsInt-1; i >= 0; i -= 1)
			tempMinPoly[i+1] = tempMinPoly[i];
		tempMinPoly[0] = constant_BigPolyT(numOfFactors);
		
		tempPoly = malloc((numOfFactorsInt)*sizeof(BigPolyTP));
		tempPoly[0] = NULL;
		
		tempMat = new_BigIntMatrixT(big_rows(A), big_rows(A));
		zeroMat = new_BigIntMatrixT(big_rows(A), big_rows(A));
		
		//Now, iterate over each of our factors,
		// finding out which ones aren't needed to 
		// minimise the matrix
		while (! allFactorsAreNeeded)
		{
			allFactorsAreNeeded = TRUE; //Assume we need all factors until proven otherwise
			factorToExclude = 0;
			
			//The zero vector doesn't need any factors
			if ((v != NULL) && (compare_BigIntMatrixT(v, zeroVect)))
				allFactorsAreNeeded = FALSE;
			
			//Create constant polynomial to represent number of factors in our temporary polynomial
			if (tempPoly[0] == NULL)
			{
				subtract_BigIntT(numOfFactors, one, temp);
				tempPoly[0] = constant_BigPolyT(temp);
			}
			else //Prevents a new BigPolyT from needing to be created each loop
			{
				add_BigPolyT(tempPoly[0], negOnePoly, tempCalcPoly);
				mod_BigPolyT(tempCalcPoly, mod, tempPoly[0]);
			}
			
			for (int tempCounter = 1; tempCounter <= numOfFactorsInt; tempCounter += 1)
			{
				factorToExclude += 1;
				
				//Construct temporary min poly to test
				for (int i = 1; i <= numOfFactorsInt; i += 1)
				{
					indexCorrection = i > factorToExclude ? i-1 : i;
					if (i != factorToExclude)
						tempPoly[indexCorrection] = tempMinPoly[i];
				}
				
				//Now, tempPoly should hold all the factors except for the one to exclude
				//Let's evaluate the matrix in this polynomial and see if it's zero
				clear_BigIntMatrixT(tempMat);
				
				//If we end up removing all factors, let the polynomial be the one polynomial
				if (compare_BigIntT(constant(tempPoly[0]), zero) == 0)
					eval_BigPolyT(onePoly, A, tempMat, mod);
				else
					eval_factored_BigPolyT(tempPoly, A, tempMat, mod);
				
				if (v != NULL)
				{
					big_mat_mul(tempMat, v, tempVect);
					modbm(tempVect, mod);
				}
				
				//If the matrix is zero, we don't need the factor we removed
				//We're storing the unneeded factors so that we can free their memory at the end
				// (the user won't be able to free them)
				if ((((v == NULL) && (compare_BigIntMatrixT(tempMat, zeroMat))) ||
				     ((v != NULL) && (compare_BigIntMatrixT(tempVect, zeroVect)))) && 
						(compare_BigIntT(temp, zero) != 0))
				{
					allFactorsAreNeeded = FALSE;
					unneededFactorsCount += 1;
					unneededFactors = realloc(unneededFactors, unneededFactorsCount*sizeof(BigPolyTP));
					unneededFactors[unneededFactorsCount-1] = tempMinPoly[factorToExclude];
				}
				
				if (! allFactorsAreNeeded)
					break;
			}
			
			//Copy tempPoly over to tempMinPoly 
			if (! allFactorsAreNeeded)
			{
				numOfFactorsInt -= 1;
				copy_BigIntT(temp, numOfFactors); //numOfFactors = numOfFactors - 1
				
				tempMinPoly = realloc(tempMinPoly, (numOfFactorsInt+1)*sizeof(BigPolyTP));
				for (int i = 0; i <= numOfFactorsInt; i += 1)
					if (i != 0) //Making sure we don't lose the minpoly's number of factors
						tempMinPoly[i] = tempPoly[i];
				
				//Giving tempMinPoly the currect number of factors
				//This is so that the free statement below doesn't remove the number of factors
				if ((numOfFactorsInt != 0) || ((v != NULL) && (compare_BigIntMatrixT(v, zeroVect))))
				{
					add_BigPolyT(tempMinPoly[0], negOnePoly, tempCalcPoly);
					mod_BigPolyT(tempCalcPoly, mod, tempMinPoly[0]);
				}
				
				tempPoly = realloc(tempPoly, numOfFactorsInt*sizeof(BigPolyTP));
			}
			
			//This case occurs when we're trying to find the minimal
			// polynomial of the zero vector
			if (numOfFactorsInt == 0)
				allFactorsAreNeeded = TRUE;
		}
	}
	
	//Free number of constants since we don't need it anymore
	if (numOfFactorsInt != 0)
	{
		tempPoly[0] = free_BigPolyT(tempPoly[0]);
		free(tempPoly);
	}
	tempPoly = NULL;
	
	//Prevents needed factors from being deallocated
	// when finding a vector's minimal polynomial
	if ((numOfFactorsInt == 0) && ((v != NULL) && (!compare_BigIntMatrixT(v, zeroVect))))
		unneededFactorsCount -= 1;
	
	//Free factors that aren't needed for minimal polynomial
	for (int i = 0; i < unneededFactorsCount; i += 1)
		unneededFactors[i] = free_BigPolyT(unneededFactors[i]);
	free(unneededFactors);
	unneededFactors = NULL;
	
	factoredCharaPoly = free_BigFactorsT(factoredCharaPoly);
	
	numOfFactors = free_BigIntT(numOfFactors);
	tempCounter  = free_BigIntT(tempCounter);
	one  = free_BigIntT(one);
	zero = free_BigIntT(zero);
	temp = free_BigIntT(temp);
	
	if (somePoly == NULL)
		charaPoly = free_BigPolyT(charaPoly);
	
	onePoly      = free_BigPolyT(onePoly);
	negOnePoly   = free_BigPolyT(negOnePoly);
	tempCalcPoly = free_BigPolyT(tempCalcPoly);
	
	tempMat  = free_BigIntMatrixT(tempMat);
	zeroMat  = free_BigIntMatrixT(zeroMat);
	zeroVect = free_BigIntMatrixT(zeroVect);
	tempVect = free_BigIntMatrixT(tempVect);

	return tempMinPoly;
}


int ccm(const BigIntMatrixTP A, 
        BigIntMatrixTP CCM, 
				const BigIntTP from, 
				const BigIntTP to,
				const BigIntTP mod)
/** Calculates a cycle converting matrix for A (C_(from->to)) and
    stores it in CCM. This function assumes CCM has been initialised
		to the zero matrix.
		Returns 1 on success, 0 otherwise. */
{
	//Making sure matrices are the correct dimensions
	if ((A->m != A->n) ||
	    (A->m != CCM->m) || (A->n != CCM->n))
		return 0;
		
	BigIntTP tempPower;
	BigIntTP tempInt;
	
	BigIntMatrixTP tempPowerMat = identity_BigIntMatrixT(big_rows(A));
	BigIntMatrixTP tempMat = new_BigIntMatrixT(big_rows(A), big_rows(A));
		
	tempPower = empty_BigIntT(1);
	tempInt   = empty_BigIntT(1);


	//Add up all the relevant matrix powers
	while (compare_BigIntT(tempPower, from) < 0)
	{
		powbm(A, tempPowerMat, tempPower, mod);
		big_mat_add(tempPowerMat, CCM, tempMat);
		modbm(tempMat, mod);
		copy_BigIntMatrixT(tempMat, CCM);
		
		//Increment counter
		add_BigIntT(tempPower, to, tempInt);
		copy_BigIntT(tempInt, tempPower);
	}
	
	tempPower = free_BigIntT(tempPower);
	tempInt   = free_BigIntT(tempInt);

	tempPowerMat = free_BigIntMatrixT(tempPowerMat);
	tempMat      = free_BigIntMatrixT(tempMat);
	
	return 1;
}
