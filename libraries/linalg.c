
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
	
	int initValue; //For feeding into the initialisation function
	
	void* (*freeFunction)(void*);
	void* (*initFunction)(int);
	void* (*copyFunction)(const void*, void*);
} GenericMatrixT, *GenericMatrixTP;


int rows(IntMatrixTP M)
/** Returns the number of rows in a given matrix. */
{
	return M->m;
}


int big_rows(BigIntMatrixTP M)
/** Returns the number of rows in the given matrix. */
{
	return M->m;
}


int cols(IntMatrixTP M)
/** Returns the number of columns in a given matrix. */
{
	return M->n;
}


int big_cols(BigIntMatrixTP M)
/** Returns the number of columns in the given matrix. */
{
	return M->n;
}


int element(IntMatrixTP const M, int row, int col)
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


bool increment_BigIntT_array(BigIntTP** intArr, 
                             int sizeRow, 
														 int sizeCol, 
														 BigIntTP const inc, 
														 BigIntTP const mod)
/** Increments the given size by size BigIntTP array by inc. 
    Used to increment through all possible matrices or vectors 
		under a particular modulus.
		Returns TRUE if the array rolls over,
		FALSE otherwise. */
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
					
					divide_BigIntT(temp, mod, carry);
					
					//Prepare temp2 for carry
					if (col+1 < sizeCol)
						copy_BigIntT(intArr[row][col+1], temp2);
					else if (row+1 < sizeRow)
						copy_BigIntT(intArr[row+1][0], temp2);
					else
					{
						copy_BigIntT(intArr[0][0], temp2);
						subtract_BigIntT(carry, one, temp); //Prevents the zero vector from being skipped when incrementing
						copy_BigIntT(temp, carry);
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
				if (a->freeFunction != NULL)
					a->matrix[row][col] = a->freeFunction(a->matrix[row][col]);
			
			free(a->matrix[row]);
			a->matrix[row] = NULL;
		}
		free(a->matrix);
		a->matrix = NULL;
		
		a->freeFunction = NULL;
		a->initFunction = NULL;
		a->copyFunction = NULL;
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
		a->matrix[i] = malloc(c*sizeof(void*));
	
	a->initValue = 0;
	
	a->freeFunction = NULL;
	a->initFunction = NULL;
	a->copyFunction = NULL;
	
	return a;
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


int set_column(IntMatrixTP v, int* const elements)
/** Returns a pointer to a new IntMatrixTP vector that contains
    the elements specified in elements. Returns NULL on error. */
{
	for (int i = 0; i < v->m; i += 1)
		v->matrix[i][0] = elements[i];
	
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


int set_GenericMatrixT_copyFunction(GenericMatrixTP a, void* (*f)(const void*, void*))
/** Sets the copying function for a GenericMatrixT
    (what it uses to copy elements to elements of its matrix).
		Returns 1 on success, 0 otherwise. */
{
	a->copyFunction = f;
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


int set_GenericMatrixT(GenericMatrixTP a, const void*** arr)
/** Sets the matrix of a given GenericMatrixT.
    Returns 1 on success, 0 otherwise. */
{
	if ((a->copyFunction == NULL) || (a->initFunction == NULL))
		return 0;
	
	for (int row = 0; row < a->m; row += 1)
	{
		for (int col = 0; col < a->n; col += 1)
		{
			a->matrix[row][col] = a->initFunction(a->initValue);
			a->copyFunction(arr[row][col], a->matrix[row][col]);
		}
	}
	
	return 1;
}


IntMatrixTP read_IntMatrixT(char* const matFilePath)
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


BigIntMatrixTP read_BigIntMatrixT(char* const matFilePath)
/** Same as read_IntMatrixT(), but for BigIntMatrixT structs.
    Returns a pointer to the new matrix on success, NULL
		otherwise. */
{
	BigIntMatrixTP M;
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
		M->matrix[row] = malloc((M->n)*sizeof(BigIntTP));
	
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


int copy_IntMatrixT(IntMatrixTP const toCopy, IntMatrixTP copyTo)
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


int copy_BigIntMatrixT(BigIntMatrixTP const toCopy, BigIntMatrixTP copyTo)
/** Same as copy_IntMatrixT(), but for BigImtMatrixT structs.
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


int compare_IntMatrixT(IntMatrixTP const M1, IntMatrixTP const M2)
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


int compare_BigIntMatrixT(BigIntMatrixTP const M1, BigIntMatrixTP const M2)
/** Same as compare_IntMatrixT(), but with BigIntMatrixT structs.
    Returns 1 if the two matrices are equal, zero otherwise. */
{
	if ((M1->m != M2->m) || (M1->n != M2->n))
		return 0;
	
	for (int row = 0; row < M1->m; row += 1)
		for (int col = 0; col < M1->n; col += 1)
			if (compare_BigIntT(M1->matrix[row][col], M2->matrix[row][col]) != 0)
				return 0;
			
	return 1;
}


int compare_BigIntMatrixT_cols(BigIntMatrixTP const M1, BigIntMatrixTP const M2, int c)
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


/* private */ IntMatrixTP det_subIntMatrixT(IntMatrixTP const M, int x, int y)
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
	

//I'll probably make this more efficient later
void printm(IntMatrixTP const M)
/** Prints an m by n matrix to stdout. 
    If zeroPad == TRUE, numbers will have
		zeros added to the left of them to align them on the console. */
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
				printf("0");
			
			printf("%d, ", M->matrix[mIndex][nIndex]);
		}
		
		for (int p = 0; p < maxDigits - num_digits(M->matrix[mIndex][M->n-1]); p += 1)
				printf("0");

		printf("%d", M->matrix[mIndex][M->n-1]);
		printf("\n");
	}
}


void printm_row(IntMatrixTP const M)
/** Same as printm(), except it prints vectors as row vectors. 
    This function also ditches the zero padding from the original. */
{
	//Making sure we're actually dealing with a column vector
	if ((M->n == 1) && (M->m >= 1))
	{
		printf("<");
		for (int mIndex = 0; mIndex < M->m-1; mIndex += 1)
			printf("%d ", M->matrix[mIndex][0]);
		
		printf("%d>", M->matrix[M->m-1][0]);
	}
	
	else
		printm(M);
}


void fprintm(FILE* file, IntMatrixTP const M)
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


void fprintm_row(FILE* file, IntMatrixTP const M)
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


void printbm(BigIntMatrixTP const M)
/** Prints out a BigIntMatrixT matrix to the console. */
{
	int maxBunches = 0;
	
	//Search for the highest amount of bunches in the matrix
	for (int x = 0; x < M->m; x += 1)
		for (int y = 0; y < M->n; y += 1)
			if (maxBunches < size(M->matrix[x][y]))
				maxBunches = size(M->matrix[x][y]);
			
	//Actually print out the entries
	for (int row = 0; row < M->m; row += 1)
	{
		for (int col = 0; col < M->n; col += 1)
		{
			//Zero padding
			for (int i = 0; i < maxBunches - size(M->matrix[row][col]); i += 1)
				for (int z = 1; z < MAXBUNCH; z *= 10)
					printf("0");
				
			printi_pad(M->matrix[row][col]);
			printf(" ");
		}
		printf("\n");
	}
}


void printbm_row(BigIntMatrixTP const M)
/** Prints a BigIntMatrixTP as a row vector. */
{
	if ((M->m > 1) && (M->n > 1))
		printbm(M);
	
	else if (M->m == 1)
	{
		printf("<");
		for (int i = 0; i < M->n-1; i += 1)
		{
			printi(M->matrix[0][i]);
			printf(" ");
		}
		printi(M->matrix[0][M->n-1]);
		printf(">");
	}
		
	
	else if (M->n == 1)
	{
		printf("<");
		for (int i = 0; i < M->m-1; i += 1)
		{
			printi(M->matrix[i][0]);
			printf(" ");
		}
		printi(M->matrix[M->m-1][0]);
		printf(">");
	}
}


void fprintbm(FILE* file, BigIntMatrixTP M)
/** Prints out a BigIntMatrixT matrix to the console. */
{
	int maxBunches = 0;
	
	//Search for the highest amount of bunches in the matrix
	for (int x = 0; x < M->m; x += 1)
		for (int y = 0; y < M->n; y += 1)
			if (maxBunches < size(M->matrix[x][y]))
				maxBunches = size(M->matrix[x][y]);
			
	//Actually print out the entries
	for (int row = 0; row < M->m; row += 1)
	{
		for (int col = 0; col < M->n; col += 1)
		{
			//Zero padding
			for (int i = 0; i < maxBunches - size(M->matrix[row][col]); i += 1)
				for (int z = 1; z < MAXBUNCH; z *= 10)
					fprintf(file, "0");
				
			fprinti(file, M->matrix[row][col]);
			fprintf(file, " ");
		}
		fprintf(file, "\n");
	}
}


void fprintbm_row(FILE* file, BigIntMatrixTP const M)
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


void fprintbm_nopad(FILE* file, BigIntMatrixTP const M)
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


void printgm(GenericMatrixTP const a)
/** Prints a generic matrix to stdout. */
{
	printf("Not implemented yet.\n");
}


void printgm_row(GenericMatrixTP const a)
/** Prints a generic matrix as a row vector. */
{
	printf("Not implemented yet.\n");
}


BigIntMatrixTP BigIntMatrixT_catalogue_get(BigIntMatrixTP** catalogue, int* catLength, BigIntMatrixTP const toFind)
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


int big_mat_add(BigIntMatrixTP const A, BigIntMatrixTP const B, BigIntMatrixTP sum)
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


int mat_mul(IntMatrixTP const A, IntMatrixTP const B, IntMatrixTP result)
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


int big_mat_mul(BigIntMatrixTP const A, BigIntMatrixTP const B, BigIntMatrixTP result)
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


int modm(IntMatrixTP M, int mod)
/** Applies a modulus to every element of a given matrix. 
    Returns 1 on success, 0 otherwise. */
{
	for (int row = 0; row < M->m; row += 1)
		for (int col = 0; col < M->n; col += 1)
			M->matrix[row][col] %= mod;
		
	return 1;
}


int modbm(BigIntMatrixTP M, BigIntTP const modulus)
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


int powbm(BigIntMatrixTP const A, 
          BigIntMatrixTP AP, 
					BigIntTP const power, 
					BigIntTP const modulus)
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


int eval_BigPolyT(BigPolyTP const p, BigIntMatrixTP const A, BigIntMatrixTP result, BigIntTP const mod)
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
	powMat  = new_BigIntMatrixT(big_rows(A), big_rows(A));
	tempMat = new_BigIntMatrixT(big_rows(A), big_rows(A));
	
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
	
	powCount  = free_BigIntT(multCount);
	multCount = free_BigIntT(powCount);
	one       = free_BigIntT(one);
	temp      = free_BigIntT(temp);
	
	return 1;
}


int eval_factored_BigPolyT(BigPolyTP* const factors, BigIntMatrixTP const A, BigIntMatrixTP result, BigIntTP const mod)
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


/* private */ int big_row_multiply(BigIntMatrixTP M, int row, BigIntTP const multiple, BigIntTP const modulus)
/** Multiplies a row by the given multiple, then takes it mod modulus.
    Returns 1 on success, 0 otherwise. */
{
	BigIntTP temp;
	
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


/* private */ int big_row_add(BigIntMatrixTP M, int addTo, int addFrom, BigIntTP const modulus)
/** Adds row addFrom to row addTo in the given matrix.
    Returns 1 on success, 0 otherwise. */
{
	BigIntTP temp;
	
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


IntMatrixTP inverse(IntMatrixTP const M, int modulus)
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


BigIntMatrixTP big_inverse(BigIntMatrixTP const M, BigIntTP const modulus)
/** If the inverse of M exists, a pointer to M's inverse
    is returned. Otherwise, returns NULL. 
		
		This function assumes the given matrix conforms to the 
		given modulus (all entries are less than the modulus). 
		
		This function is the same as inverse(), but for BigIntMatrixTs. */
{
	//If we're given a non-square matrix
	if ((M->m != M->n))
		return NULL;
	
	bool hasLeadEntry;
	
	BigIntTP tempInverse;
	BigIntTP zero;
	BigIntTP one;
	int oneArr[1] = {1};
	
	BigIntTP bigMult; //Holds the loop variable mult so we can multiply matrix elements with it
	BigIntTP temp;
	BigIntTP temp2;
	
	BigIntMatrixTP toReduce = new_BigIntMatrixT(M->m, M->m);
	BigIntMatrixTP inv      = identity_BigIntMatrixT(M->m);
	copy_BigIntMatrixT(M, toReduce);
	
	#ifdef VERBOSE
	printf("Matrix to reduce:\n");
	printbm(toReduce);
	#endif
	
	zero        = empty_BigIntT(1);
	one         = new_BigIntT(oneArr, 1);
	temp        = empty_BigIntT(1);
	temp2       = empty_BigIntT(1);
	bigMult     = empty_BigIntT(1);
	tempInverse = empty_BigIntT(1);
	
	//Converting toReduce to upper triangular (row echelon) form
	for (int focusRow = 0; focusRow < M->m; focusRow += 1)
	{
		//Make sure each row has to find an invertible leading entry
		hasLeadEntry = FALSE;
		
		//Find a row with a non-zero first entry
		for (int nonzero = focusRow; nonzero < M->m; nonzero += 1)
		{
			//Checking to see whether the leading entry is one
			//OR
			//Checking to see whether the leading entry is nonzero and has an inverse
			if ((compare_BigIntT(toReduce->matrix[nonzero][focusRow], one) == 0) ||
					((compare_BigIntT(toReduce->matrix[nonzero][focusRow], zero) != 0) &&
					 (compare_BigIntT(big_gcd(toReduce->matrix[nonzero][focusRow], modulus), one) == 0)))
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
				
				big_row_swap(toReduce, focusRow, nonzero);
				big_row_swap(inv, focusRow, nonzero);
				
				#ifdef VERBOSE
				printf("Swapped row %d and %d.\n", focusRow, nonzero);
				printbm(toReduce);
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
				//Hopefully finging a row that we can add to our current row to get an invertible number
				copy_BigIntT(one, bigMult);
				while (compare_BigIntT(bigMult, modulus) < 0)
				{
					multiply_BigIntT(bigMult, M->matrix[rowToMaybeAdd][focusRow], temp);
					add_BigIntT(temp, M->matrix[focusRow][focusRow], temp2);
					big_num_inverse(temp2, modulus, tempInverse);
					
					if (compare_BigIntT(tempInverse, zero) != 0)
					{
						//Add rows to get an invertible number as a leading entry
						copy_BigIntT(bigMult, temp);
						while (compare_BigIntT(temp, zero) > 0)
						{
							big_row_add(toReduce, focusRow, rowToMaybeAdd, modulus);
							big_row_add(inv, focusRow, rowToMaybeAdd, modulus);
							
							//Decrement temp
							subtract_BigIntT(temp, one, temp2);
							copy_BigIntT(temp2, temp);
						}
						
						#ifdef VERBOSE
						printf("Added row %d to row %d a total of ", rowToMaybeAdd, focusRow);
						printi(bigMult);
						printf(" times.\n");
						printbm(toReduce);
						#endif
						
						//Now, get the leading entry to one
						big_row_multiply(toReduce, focusRow, tempInverse, modulus);
						big_row_multiply(inv, focusRow, tempInverse, modulus);
						
						#ifdef VERBOSE
						printf("Multipled row %d by ", focusRow);
						printi(tempInverse);
						printf(".\n");
						printbm(toReduce);
						#endif
						
						hasLeadEntry = TRUE;
						break;
					}
					
					//Increment bigMult
					add_BigIntT(bigMult, one, temp);
					copy_BigIntT(temp, bigMult);
				}
				
				if (hasLeadEntry)
					break;
			}
		}
		
		//If we couldn't find a suitable row to create an invertible element
		if (!hasLeadEntry)
		{
			toReduce = free_BigIntMatrixT(toReduce);
			inv = free_BigIntMatrixT(inv);
			
			zero        = free_BigIntT(zero);
			one         = free_BigIntT(one);
			tempInverse = free_BigIntT(tempInverse);
			bigMult     = free_BigIntT(bigMult);
			temp        = free_BigIntT(temp);
			temp2       = free_BigIntT(temp2);
			
			#ifdef VERBOSE
			printf("No suitable row was found to create an invertible entry.\n");
			#endif
			
			return NULL;
		}

		//If the leading entry is a 1, we don't need to find an inverse
		if (compare_BigIntT(toReduce->matrix[focusRow][focusRow], one) != 0)
		{
			big_num_inverse(toReduce->matrix[focusRow][focusRow], modulus, tempInverse);
			big_row_multiply(toReduce, focusRow, tempInverse, modulus);
			big_row_multiply(inv, focusRow, tempInverse, modulus);
			
			#ifdef VERBOSE
			printf("Multiplied row %d by ", focusRow);
			printi(tempInverse);
			printf(".\n");
			printbm(toReduce);
			#endif
		}
		
		//Now, we clear out all nonzero entries in our current pivot column
		for (int i = focusRow+1; i < M->m; i += 1)
		{
			subtract_BigIntT(modulus, toReduce->matrix[i][focusRow], temp);
			mod_BigIntT(temp, modulus, temp2); //temp2 is acting as "numTimesToAdd" here
			
			#ifdef VERBOSE
			printf("Added row %d to row %d a total of ", focusRow, i);
			printi(temp2);
			printf(" times.\n");
			#endif
			
			while (compare_BigIntT(temp2, zero) != 0)
			{
				big_row_add(toReduce, i, focusRow, modulus);
				big_row_add(inv, i, focusRow, modulus);
				
				//Decrement temp2
				subtract_BigIntT(temp2, one, temp);
				copy_BigIntT(temp, temp2);
			}
			
			#ifdef VERBOSE
			printbm(toReduce);
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
			subtract_BigIntT(modulus, toReduce->matrix[element][i], temp);
			mod_BigIntT(temp, modulus, temp2); //temp2 is acting as "numTimesToAdd" here
			
			#ifdef VERBOSE
			printf("Added row %d to row %d a total of ", i, element);
			printi(temp2);
			printf(" times.\n");
			#endif

			while (compare_BigIntT(temp2, zero) != 0)
			{
				big_row_add(toReduce, element, i, modulus);
				big_row_add(inv, element, i, modulus);
				
				//Decrement temp2
				subtract_BigIntT(temp2, one, temp);
				copy_BigIntT(temp, temp2);
			}
			
			#ifdef VERBOSE
			printbm(toReduce);
			#endif
		}
	}
	
	zero        = free_BigIntT(zero);
	one         = free_BigIntT(one);
	tempInverse = free_BigIntT(tempInverse);
	bigMult     = free_BigIntT(bigMult);
	temp        = free_BigIntT(temp);
	temp2       = free_BigIntT(temp2);
	
	toReduce = free_BigIntMatrixT(toReduce);
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
	
	BigPolyTP subResult = empty_BigPolyT();
	
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
BigPolyTP chara_poly(BigIntMatrixTP const A, BigIntTP mod)
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


BigPolyTP* min_poly(BigIntMatrixTP const A, BigIntTP const mod)
/** Calculates A's minimum polynomial, returns it in its factored 
    form on success, NULL otherwise. This won't be tested for 
		composite moduli since minimum polynomials don't really work 
		under composite moduli. */
{
	BigPolyTP* tempMinPoly = NULL;
	BigPolyTP* tempPoly = NULL;
	BigFactorsTP factoredCharaPoly = NULL;
	BigPolyTP* unneededFactors = NULL;
	BigPolyTP  charaPoly;
	int unneededFactorsCount = 0;
	
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
	
	bool allFactorsAreNeeded = FALSE;
	int indexCorrection;
	
	charaPoly = chara_poly(A, mod);
	if (charaPoly != NULL)
	{
		factoredCharaPoly = factor_BigPolyT(charaPoly, mod);
		numOfFactorsInt = count_factors(factoredCharaPoly);
		
		one  = new_BigIntT(oneArr, 1);
		zero = empty_BigIntT(1);
		temp = empty_BigIntT(1);
		
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
		
		tempMat = new_BigIntMatrixT(big_rows(A), big_rows(A));
		zeroMat = new_BigIntMatrixT(big_rows(A), big_rows(A));
		
		//Now, iterate over each of our factors,
		// finding out which ones aren't needed to 
		// minimise the matrix
		while (! allFactorsAreNeeded)
		{
			allFactorsAreNeeded = TRUE; //Assume we need all factors until proven otherwise
			factorToExclude = 0;
			
			//Create constant polynomial to represent number of factors in our temporary polynomial
			subtract_BigIntT(numOfFactors, one, temp);
			tempPoly[0] = constant_BigPolyT(temp);
			
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
				eval_factored_BigPolyT(tempPoly, A, tempMat, mod);
				
				//If the matrix is zero, we don't need the factor we removed
				//We're storing the unneeded factors so that we can free their memory at the end
				// (the user won't be able to free them)
				if ((compare_BigIntMatrixT(tempMat, zeroMat) == 1) && (compare_BigIntT(temp, zero) != 0))
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
					tempMinPoly[i] = tempPoly[i];
				
				tempPoly = realloc(tempPoly, numOfFactorsInt*sizeof(BigPolyTP));
				
				/*
				printf(":)\n");
				printf("tempPoly: ");
				old_printpf(tempPoly);
				printf("\ntempMinPoly: ");
				old_printpf(tempMinPoly);
				printf("\n");
				*/
				
				//Free constant so that we don't need to worry about it later
				tempPoly[0] = free_BigPolyT(tempPoly[0]); 
			}
			
			//Change tempPoly[0] to store the correct value
			else
			{
				tempMinPoly[0] = free_BigPolyT(tempPoly[0]); 
				tempMinPoly[0] = constant_BigPolyT(numOfFactors);
			}
		}
	}
	
	free(tempPoly);
	tempPoly = NULL;
	
	//Free factors that aren't needed for minimum polynomial
	for (int i = 0; i < unneededFactorsCount; i += 1)
		unneededFactors[i] = free_BigPolyT(unneededFactors[i]);
	free(unneededFactors);
	unneededFactors = NULL;
	
	factoredCharaPoly = free_BigFactorsT(factoredCharaPoly);
	charaPoly = free_BigPolyT(charaPoly);
	
	numOfFactors = free_BigIntT(numOfFactors);
	tempCounter  = free_BigIntT(tempCounter);
	one  = free_BigIntT(one);
	zero = free_BigIntT(zero);
	temp = free_BigIntT(temp);
	
	tempMat = free_BigIntMatrixT(tempMat);
	zeroMat = free_BigIntMatrixT(zeroMat);
	
	return tempMinPoly;
}


int ccm(BigIntMatrixTP const A, 
        BigIntMatrixTP CCM, 
				BigIntTP const from, 
				BigIntTP const to,
				BigIntTP const mod)
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


//This function currently only works with 2 by 2 matrices
int* eigenvalues(IntMatrixTP const F, int modulus)
/** Returns the found eigenvalues of the given matrix,
    mod the given modulus. The number of values returned
		with the pointer will always be the first number in
		the pointer. 
		
		Returns NULL if no eigenvalues exist or if the
    given matrix isn't square. */
{
	if ((F->m != 2) || (F->m != F->n))
		return NULL;
	
	//The quadratic formula doesn't really work here due
	// to the system being modulated, so the best we can do
	// is guess and check for solutions. There's probably some
	// more efficient way to do this, but I can't be bothered
	// right now.
	int* values = malloc((modulus + 1)*sizeof(int));
	values[0] = 0;
	
	int coefficient = F->matrix[0][0] + F->matrix[1][1];
	int constant    = F->matrix[0][0]*F->matrix[1][1] - F->matrix[0][1]*F->matrix[1][0];
	
	for (int i = 0; i < modulus; i += 1)
		if ((i*i - coefficient*i + constant) % modulus == 0)
		{
			values[0] += 1;
			values[values[0]] = i;
		}

	return values;
}


//Currently only works with 2x2 matrices
IntMatrixTP eigenvector(IntMatrixTP const F, int eigenvalue, int modulus)
/** Returns an eigenvector of the given matrix and eigenvalue, under the
    given modulus. It is assumed that the given eigenvalue is valid for
		the given matrix and modulus. */
{
	if ((F->m != 2) || (F->m != F->n))
		return NULL;
	
	IntMatrixTP toReduce = new_IntMatrixT(F->m, F->n);
	copy_IntMatrixT(F, toReduce);
	toReduce->matrix[0][0] -= eigenvalue;
	toReduce->matrix[1][1] -= eigenvalue;
	
	//Making sure the entries in our matrix are positive
	if (toReduce->matrix[0][0] < 0)
		toReduce->matrix[0][0] += modulus;
	if (toReduce->matrix[1][1] < 0)
		toReduce->matrix[1][1] += modulus;
	
	int currentRow;         //Holds the current row we're reducing
	int entryInverse = -1;  //Holds a number inverse for reducing
	int numTimesToAdd;      //Holds how many times we add one row to another
	
	//Now, we reduce the matrix
	//We don't need to worry about not getting leading entries 
	// since we're assuming the given parameters will work
	
	for (currentRow = 0; currentRow < F->m; currentRow += 1)
	{
		//Make sure the currentRow has a leading entry
		for (int row = currentRow; row < F->m; row += 1)
		{
			//If we found a row with a leading entry
			if (toReduce->matrix[row][currentRow] != 0)
			{
				//Check to see if the leading entry is invertible
				entryInverse = num_inverse(toReduce->matrix[row][row], modulus);
				if (entryInverse != -1)
				{
					//Swap rows if needed
					if (row != currentRow)
						row_swap(toReduce, currentRow, row);
					
					break;
				}
			}
		}
		
		//Reduce currentRow's leading entry to 1, if needed
		if (entryInverse != -1)
			row_multiply(toReduce, currentRow, entryInverse, modulus);
		
		//Now, get rid of all other entries in the current pivot column
		for (int row = currentRow + 1; row < F->m; row += 1)
		{
			//Checking to see if next entry in pivot column is nonzero
			if (toReduce->matrix[row][currentRow] != 0)
			{
				numTimesToAdd = (modulus - toReduce->matrix[row][currentRow]) % modulus;
				for (int t = 0; t < numTimesToAdd; t += 1)
					row_add(toReduce, currentRow, row, modulus);
			}
		}
		//Pivot column should now be cleared out
	}
	
	//Now, we convert to reduced row echelon form
	for (int currentRow = F->m - 1; currentRow > 0; currentRow -= 1)
	{
		//If the leading entry in our currentRow isn't zero
		if (toReduce->matrix[currentRow][currentRow] != 0)
		{
			for (int row = currentRow - 1; row >= 0; row -= 1)
			{
				numTimesToAdd = (modulus - toReduce->matrix[row][currentRow]) % modulus;
				for (int t = 0; t < numTimesToAdd; t += 1)
					row_add(toReduce, currentRow, row, modulus);
			}
		}
	}
	
	printm(toReduce);
	
	//Now, we need to go through the matrix to see what form
	// the eigenvector needs to take
	
	//Have some form of matrix that keeps track of which variables can
	// be substituted, and go through the reduced matrix and get values.
	
	toReduce = free_IntMatrixT(toReduce);
	return NULL;
}
