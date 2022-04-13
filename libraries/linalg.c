
/* A simple linear algebra library to facilitate 
 *  research with Dr. Mendivil.
 *
 * Apr 7, 2022
 *
 */
 
#include <stdlib.h>
#include <stdio.h>

// https://stackoverflow.com/a/6317375

//Maybe I could make a generic struct that has a union as its matrix
//That way, if we wanted to switch to doubles, we could
//Or maybe a simple void pointer will do

typedef struct intmatrix
/** Structure for holding matrix data and metadata (dimensions, e.g.). */
{
	int** matrix;
	int m;
	int n;
}
IntMatrixT, *IntMatrixTP;


int rows(IntMatrixTP M)
/** Returns the number of rows in a given matrix. */
{
	return M->m;
}


int cols(IntMatrixTP M)
/** Returns the number of columns in a given matrix. */
{
	return M->n;
}


int element(IntMatrixTP M, int row, int col)
/** Returns the element at M->matrix[row][col].
    Returns 0 for elements that aren't defined. */
{
	if ((row >= M->m) || (col >= M->n) ||
	    (row < 0) || (col < 0))
		return 0;
		
	else
		return M->matrix[row][col];
}


IntMatrixTP free_IntMatrixT(IntMatrixTP M)
/** Frees memory taken by M. Returns NULL. */
{
	//Freeing each row of the matrix
	for (int i = 0; i < M->m; i += 1)
		free(M->matrix[i]);
	
	//Freeing the matrix itself
	free(M->matrix);
	free(M);
	
	return NULL;
}


int new_IntMatrixT(IntMatrixTP* M, int r, int c)
/** Creates an empty m by n matrix.
    Returns 1 on success, 0 otherwise. */
{
	//If given dimensions don't make sense
	if ((r < 0) || (c < 0))
		return 0;
	
	
	*M = malloc(sizeof(IntMatrixT));
	
	(*M)->matrix = malloc(r*sizeof(int*));
	
	for (int i = 0; i < r; i += 1)
		(*M)->matrix[i] = malloc(c*sizeof(int));
	
	//Don't forget to set the matrix's dimensions!
	(*M)->m = r;
	(*M)->n = c;
	
	return 1;
}


int read_IntMatrixT(char* const matFilePath, IntMatrixTP* M)
/** Initialises the given matrix M from the given .matrix file.
    This method assumes matFile points to the beginning of the file.
		Returns 1 upon success, 0 otherwise. **/
{
	FILE* matFile = fopen(matFilePath, "r");
	
	if (matFile == NULL)
	{
		fprintf(stderr, "Unable to open matrix file.\n");
		return 0;
	}
	
	//Allocating memory for M
	*M = malloc(sizeof(IntMatrixT));
	
	//Getting dimensions of the data
	if (fscanf(matFile, "%d %d", &((*M)->m), &((*M)->n)) != 2)
	{
		fprintf(stderr, "Unable to read matrix dimensions.\n");
		return 0;
	}
	
	//Allocating memory for M
	(*M)->matrix = malloc(((*M)->m)*sizeof(int*));
	for (int i = 0; i < (*M)->m; i += 1)
		(*M)->matrix[i] = malloc(((*M)->n)*sizeof(int));
	
	//Filling F with actual values
	for (int row = 0; row < (*M)->m; row += 1)
		for (int column = 0; column < (*M)->n; column += 1)
			if (fscanf(matFile, "%d", &((*M)->matrix[row][column])) != 1)
			{
				fprintf(stderr, "Unable to read matrix data.\n");
				return 0;
			}
			

	if (fclose(matFile) == EOF)
	{
		fprintf(stderr, "Unable to close matrix file.\n");
		return 0;
	}
		
	//Function read matrix's data properly
	return 1;
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


//Maybe pad numbers with zeros in the future to make output easier to read?
//This would require dividing each number by 10 numerous times to see
// how many digits comprises it. (or take log base 10)
//Maybe we could have a switch to turn zero padding on and off...
//I'll probably make this more efficient later
void printm(IntMatrixTP M)
/** Prints an m by n matrix to stdout. */
{	
	for (int mIndex = 0; mIndex < M->m; mIndex += 1)
	{
		for (int nIndex = 0; nIndex < M->n-1; nIndex += 1)
			printf("%d, ", M->matrix[mIndex][nIndex]);
		printf("%d", M->matrix[mIndex][M->n-1]);
		printf("\n");
	}
}


int mat_mul(const IntMatrixTP A, const IntMatrixTP B, IntMatrixTP result)
/** Computes AB, stores result in result. 
    This function DOES check to make sure result is the proper dimensions.
		This function assumes all three matrices have been initialised.
    Returns 1 upon success, 0 otherwise. */ 
{	
	//Checking to see if result has the appropriate dimensions
	if ((result->m != A->m) || (result->n != B->n))
		return 0;
	
	
	for (int mIndex = 0; mIndex < A->m; mIndex += 1) //For each row of A
		for (int rIndex = 0; rIndex < B->n; rIndex += 1) //For each column of B
		{
			//Clear out entry before we start adding multiplications
			result->matrix[mIndex][rIndex] = 0; 
			
			for (int nIndex = 0; nIndex < A->n; nIndex += 1) //Each multiplication
			{
				/*printf("m: %d, n: %d, r: %d, A[m][n]: %d, B[n][r]: %d\n",
				mIndex, nIndex, rIndex, A[mIndex][nIndex], B[nIndex][rIndex]); */
				result->matrix[mIndex][rIndex] += (A->matrix[mIndex][nIndex]) * (B->matrix[nIndex][rIndex]);
			}
		}
				
	return 1;
}

//Dr. Mendivil told me the modulus has the potential to be negative.
// I tested it with some positive ints, and it always seems to give
// positive numbers. This may be something to consider in the future,
// but not right now.

int modm(IntMatrixTP M, int mod)
/** Applies a modulus to every element of a given matrix. 
    Returns 1 on success, 0 otherwise. */
{
	for (int row = 0; row < M->m; row += 1)
		for (int col = 0; col < M->n; col += 1)
			M->matrix[row][col] %= mod;
		
	return 1;
}

//UNFINISHED
int det(IntMatrixTP M)
/** Returns the determinant of a matrix. */
{
	return 1;
}