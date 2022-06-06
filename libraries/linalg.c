
/* Apr 7, 2022
 *
 */
 
#include <stdlib.h>
#include <stdio.h>

#include "../headers/factors.h"
#include "../headers/modular.h" //Allows for taking modular square roots

/* The following resources were used as a reference:
https://stackoverflow.com/a/6317375
*/

//Maybe I could make a generic struct that has a union as its matrix
//That way, if we wanted to switch to doubles, we could
//Or maybe a simple void pointer will do

typedef enum boolean {FALSE, TRUE} bool;

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


int set_column(IntMatrixTP v, int* const elements)
/** Returns a pointer to a new IntMatrixTP vector that contains
    the elements specified in elements. Returns NULL on error. */
{
	for (int i = 0; i < v->m; i += 1)
		v->matrix[i][0] = elements[i];
	
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


IntMatrixTP read_IntMatrixT(char* const matFilePath)
/** Initialises the given matrix M from the given .matrix file.
    This method assumes matFile points to the beginning of the file.
		Returns 1 upon success, 0 otherwise. **/
{
	FILE* matFile = fopen(matFilePath, "r");
	
	if (matFile == NULL)
		return NULL;
	
	//Allocating memory for M
	IntMatrixTP M = malloc(sizeof(IntMatrixT));
	
	//Getting dimensions of the data
	if (fscanf(matFile, "%d %d", &(M->m), &(M->n)) != 2)
	{
		M = free_IntMatrixT(M);
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


int num_digits(int num)
/** Calculuates how many digits an integer has and returns
    that number. The function assumes the number given is
		positive. */
{
	int count = 0;
	while (TRUE)
	{
		num /= 10;
		count += 1;
		
		if (num == 0)
			return count;
	}
}
	

//I'll probably make this more efficient later
void printm(IntMatrixTP M, bool zeroPad)
/** Prints an m by n matrix to stdout. 
    If zeroPad == TRUE, numbers will have
		zeros added to the left of them to align them on the console. */
{
	int maxDigits = 0;
	
	if (zeroPad)
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


int mat_mul(IntMatrixTP const A, IntMatrixTP const B, IntMatrixTP result)
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
				mIndex, nIndex, rIndex, A->matrix[mIndex][nIndex], B->matrix[nIndex][rIndex]);
				printf("Mini-mult: %d\n", A->matrix[mIndex][nIndex] * B->matrix[nIndex][rIndex]); */
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
	
	IntMatrixTP toReduce = new_IntMatrixT(M->m, M->m);
	IntMatrixTP inv      = identity_IntMatrixT(M->m);
	copy_IntMatrixT(M, toReduce);
	
	#ifdef VERBOSE
	printf("Matrix to reduce:\n");
	printm(toReduce, TRUE);
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
					break;
				
				row_swap(toReduce, focusRow, nonzero);
				row_swap(inv, focusRow, nonzero);
				
				#ifdef VERBOSE
				printf("Swapped row %d and %d.\n", focusRow, nonzero);
				printm(toReduce, TRUE);
				#endif

				break;
			}
		}
		
		//If we couldn't find a nonzero entry for our pivot column
		if (!hasLeadEntry)
		{
			toReduce = free_IntMatrixT(toReduce);
			inv = free_IntMatrixT(inv);
			
			#ifdef VERBOSE
			printf("No nonzero, invertible leading entry could be found.\n");
			#endif
			
			return NULL;
		}

		//If the leading entry is a 1, we don't need to find an inverse
		if (toReduce->matrix[focusRow][focusRow] != 1)
		{
			//Finding the inverse of our leading entry (i)
			for (int i = 0; i < modulus; i += 1)
			{
				if (toReduce->matrix[focusRow][focusRow]*i % modulus == 1)
				{
					row_multiply(toReduce, focusRow, i, modulus);
					row_multiply(inv, focusRow, i, modulus);
					
					#ifdef VERBOSE
					printf("Multiplied row %d by %d.\n", focusRow, i);
					printm(toReduce, TRUE);
					#endif
					
					break;
				}
			}
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
				printm(toReduce, TRUE);
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
				printm(toReduce, TRUE);
			}
			#endif
		}
	}
	
	toReduce = free_IntMatrixT(toReduce);
	return inv;
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
	
	printm(toReduce, TRUE);
	
	//Now, we need to go through the matrix to see what form
	// the eigenvector needs to take
	
	//Have some form of matrix that keeps track of which variables can
	// be substituted, and go through the reduced matrix and get values.
	
	toReduce = free_IntMatrixT(toReduce);
	return NULL;
}
