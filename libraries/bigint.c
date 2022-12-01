
/*
A super basic arbitrary precision library for
integers. Will be used to compute big examples 
of LCA systems if needed.

May 24, 2022
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h> //For strtoBIT()
#include <math.h>

#include "../headers/linalg.h" //num_digits()
#include "../headers/helper.h" //append_int()

//How big each bunch in a BigIntT can be
//This limit was chosen so that the multiply_by_ten() function
// doesn't cause an overflow.
const int MAXBUNCH       = 100000000;
const int MAXBUNCHDIGITS = 8;

/*
BigIntT structs hold numbers in little endian style,
meaning the smallest-valued number bunches will appear 
first in the int*.
*/

/*
For the time being, BigIntT structs are unsigned, 
since I only need unsigned functionality.
May or may not add negative signifiers in the future.
*/
typedef struct bigint
{
	int* theInt; //A pointer of integers representing the integer
	int    size; //Keeps track of how many numbers are in theInt
}
BigIntT, *BigIntTP;


BigIntTP free_BigIntT(BigIntTP n)
/** Frees the memory used by a given BigIntT and
    returns NULL. */
{
	if (n != NULL)
	{
		free(n->theInt);
		n->theInt = NULL;
		free(n);
	}
	return NULL;
}


BigIntTP new_BigIntT(const int* initNum, int len)
/** Returns a pointer to an initialised
    BigIntT struct. Returns NULL on error. */
{
	BigIntTP newBigInt = malloc(sizeof(BigIntT));
	newBigInt->theInt  = malloc(len*sizeof(int));
	newBigInt->size    = len;
	
	//Set our new BigInt's bunches
	for (int i = 0; i < len; i += 1)
		newBigInt->theInt[i] = initNum[i];
	
	return newBigInt;
}


BigIntTP empty_BigIntT(int zeros)
/** Returns a pointer to a BigIntT struct
    which has been filled with zeros and has
		size zeros. */
{
	BigIntTP newBigInt = malloc(sizeof(BigIntT));
	newBigInt->size    = zeros;
	newBigInt->theInt  = calloc(zeros, sizeof(int));
	
	return newBigInt;
}


int is_zero(const BigIntTP a)
/** Returns 1 if a is zero, 0 otherwise. */
{
	for (int i = 0; i < a->size; i += 1)
		if (a->theInt[i] != 0)
			return 0;
		
	return 1;
}


BigIntTP** new_BigIntT_array(int rows, int cols)
/** Creates a 2D array of BigIntTs with the specified rows
    and cols, and initialised it with zeros. Returns a pointer
    to the array on success, NULL otherwise. */
{
	BigIntTP** arr = malloc(rows*sizeof(BigIntTP*));
	for (int i = 0; i < rows; i += 1)
	{
		arr[i] = malloc(cols*sizeof(BigIntTP));
		for (int j = 0; j < cols; j += 1)
			arr[i][j] = empty_BigIntT(1);
	}
	
	return arr;
}


BigIntTP** free_BigIntT_array(BigIntTP** arr, int rows, int cols)
/** Frees an array of BigIntTs. Returns NULL. */
{
	for (int i = 0; i < rows; i += 1)
	{
		for (int j = 0; j < cols; j += 1)
			arr[i][j] = free_BigIntT(arr[i][j]);
		
		free(arr[i]);
		arr[i] = NULL;
	}
	free(arr);
	return NULL;
}


int strtoBIT(const char* numStr, BigIntTP* theBig)
/** Takes a numerical string and creates a BigIntT
    struct using it, storing it in theBig. 
		
		This function assumes theBig has been declared, and
		that MAXBUNCH is a power of ten.
		
		Returns 1 on success, 0 otherwise. */
{
	const int bunchLength = num_digits(MAXBUNCH) - 1;
	
	int substrstart;      //Holds the start of the substring in the for-loop below
	int substrlen;        //Holds how long each substring should be
	int bunchCounter = 0; //For properly storing bunches
	int numStrLength = strlen(numStr);
	
	int*  bunches; //Holds the BigIntT bunches of our number
	char* tempStr; //Holds info about whether the number was read correctly
	char* substr = malloc((bunchLength+1)*sizeof(char)); //Holds the substring for each bunch
	
	bunches = malloc((numStrLength/bunchLength + 1)*sizeof(int));
	
	//Note that the constants used in this loop would
	// need to change if we ever change the bunch size for
	// BigIntT structs. We use 4 because the bunch size is
	// currently 9999, or 4 base 10 digits.
	for (int bunch = numStrLength; 
	bunch >= (numStrLength % bunchLength == 0) ? 1 : 0; //This prevents an extra bunch from being read
	bunch -= bunchLength)
	{
		//Getting the correct substring for the next bunch
		substrstart = (bunch-bunchLength >= 0) ? bunch - bunchLength : 0;
		substrlen   = (bunch-bunchLength >= 0) ? bunchLength : numStrLength % bunchLength;
		
		strncpy(substr, numStr+substrstart, substrlen);
		substr[substrlen] = '\0'; //Adding null byte manually
		
		//Store bunches
		bunches[bunchCounter] = (int)strtol(substr, &tempStr, 10);
		bunchCounter += 1;
		
		//If we read an invalid character
		if (tempStr[0] != '\0')
		{
			free(bunches);
			free(substr);
			return 0;
		}
	}
	
	//Now, we actually create the BigIntT
	*theBig = new_BigIntT(bunches, bunchCounter);
	
	free(bunches);
	free(substr);
	
	return 1;
}


int size(const BigIntTP n)
/** Returns the size of the BigIntT passed. */
{
	return n->size;
}


int extract_bunch(const BigIntTP n, int b)
/** Returns the integer at the given bunch b.
    Returns -1 on error. */
{
	if ((b < 0) || (b >= n->size))
		return -1;
	
	else
		return n->theInt[b];
}


int set_bunch(BigIntTP n, int pos, int val)
/** Sets the value at a particular bunch.
    This function was made to facilitate laziness. 
		Returns 1 on success, 0 otherwise. */
{
	if ((pos < 0) || (pos >= n->size))
		return 0;
	
	n->theInt[pos] = val;
	return 1;
}


int append_BigIntT(char* dest, const BigIntTP src)
/** Appends a BigIntT to the given string.
    This function assumes there's enough space in the string to
		hold the BigIntT.
    Returns 1 on success, 0 otherwise. */
{
	for (int b = src->size-1; b >= 0 ; b -= 1)
	{
		//Add padding zeros if necessary
		if (b != src->size-1)
			for (int i = 0; i < MAXBUNCHDIGITS - num_digits(src->theInt[b]); i += 1)
				append_int(dest, 0);
			
		append_int(dest, src->theInt[b]);
	}
	
	return 1;
}


void printi(const BigIntTP n)
/** Prints a BigIntT struct to the screen as a normal int would. */
{
	int power = 0;
	int tempMaxBunch = MAXBUNCH;
	
	if (MAXBUNCH % 10 == 0)
	{
		//Getting what power of 10 our MAXBUNCH is
		while (tempMaxBunch > 1)
		{
			tempMaxBunch /= 10;
			power += 1;
		}
		
		for (int i = (n->size)-1; i >= 0; i -= 1)
		{
			if (i != n->size-1)
				for (int d = 0; d < power - num_digits(n->theInt[i]); d += 1)
					printf("0");

			printf("%d", n->theInt[i]);
		}
	}
}


void printi_pad(const BigIntTP n)
/** Same as printi(), except zero padding is added. */
{
	int power = 0;
	int tempMaxBunch = MAXBUNCH;
	
	if (MAXBUNCH % 10 == 0)
	{
		//Getting what power of 10 our MAXBUNCH is
		while (tempMaxBunch > 1)
		{
			tempMaxBunch /= 10;
			power += 1;
		}
		
		for (int i = n->size-1; i >= 0; i -= 1)
		{
			//Zero padding
			for (int d = 0; d < power - num_digits(n->theInt[i]); d += 1)
				printf("0");
			
			printf("%d", n->theInt[i]);
		}
	}
}


void fprinti(FILE* file, BigIntTP n)
/** Same as printi, but outputs to a file stream. */
{
	int power = 0;
	int tempMaxBunch = MAXBUNCH;
	
	if (MAXBUNCH % 10 == 0)
	{
		//Getting what power of 10 our MAXBUNCH is
		while (tempMaxBunch > 1)
		{
			tempMaxBunch /= 10;
			power += 1;
		}
		
		for (int i = n->size-1; i >= 0; i -= 1)
		{
			//Zero padding, only added if we're not on the most significant bunch
			if (i != n->size-1)
				for (int d = 0; d < power - num_digits(n->theInt[i]); d += 1)
					fprintf(file, "0");
			
			fprintf(file, "%d", n->theInt[i]);
		}
	}
}


int reduce_BigIntT(BigIntTP toReduce)
/** Removes any unnecessary bunches from toReduce
    and sets its size accordingly. Returns 1 on success,
		0 otherwise. */
{
	int counter = 0;
	
	//Iterate backwards through the number, seeing
	// if there's any bunches we can remove
	for (int i = (toReduce->size)-1; i >= 0; i -= 1)
	{
		if (toReduce->theInt[i] != 0)
			break;
		
		counter += 1;
	}
	
	//Get rid of extra bunches, update size
	toReduce->size -= counter;
	if (toReduce->size <= 0)
		toReduce->size = 1;
	
	//Only spend the time reallocating if something actually changed
	if (counter > 0)
		toReduce->theInt = realloc(toReduce->theInt, (toReduce->size)*sizeof(int));

	return 1;
}


int clear_BigIntT(BigIntTP toClear)
/** Function to reset a BigIntT number to 0.
    Returns 1 on success, 0 otherwise. */
{
	for (int i = 0; i < toClear->size; i += 1)
		toClear->theInt[i] = 0;
	
	return 1;
}


int resize_BigIntT(BigIntTP toResize, int newSize)
/** Changes the size of the provided BigIntT and
    reallocates its memory appropriately. 
		Returns 1 on success, 0 otherwise. */
{
	//Again, only do stuff if we have to
	if (toResize->size != newSize)
	{
		toResize->size   = newSize;
		toResize->theInt = realloc(toResize->theInt, (toResize->size)*sizeof(int));
	}
	
	return 1;
}


int copy_BigIntT(const BigIntTP setTo, BigIntTP toSet)
/** Copies setTo to toSet and returns 1 on completion.
    Returns 0 on error. */
{
	//Reallocate toSet's memory appropriately
	if (toSet->size != setTo->size)
	{
		toSet->size = setTo->size;
		toSet->theInt = realloc(toSet->theInt, (toSet->size)*sizeof(int));
	}
	
	//Copy contents
	for (int i = 0; i < toSet->size; i += 1)
		toSet->theInt[i] = setTo->theInt[i];
	
	return 1;
}


int compare_BigIntT(const BigIntTP AA, const BigIntTP BB)
/** Returns negative if A < B, 0 if A == B, positive if A > B */
{
	BigIntTP A, B;
	int returnVal = 0; //Assume the numbers are equal until proven otherwise
	
	//Check to see if AA and BB are reduced. If so,
	// no need to allocate new BigIntTPs to reduce!
	
	//The == 1 on each expression prevents zeros from being reallocated
	if ((AA->theInt[AA->size-1] != 0) || (AA->size == 1))
		A = AA;
	
	else
	{
		A = empty_BigIntT(1);
		copy_BigIntT(AA, A);
		reduce_BigIntT(A);
	}
		
	if ((BB->theInt[BB->size-1] != 0) || (BB->size == 1))
		B = BB;
	
	else
	{
		B = empty_BigIntT(1);
		copy_BigIntT(BB, B);
		reduce_BigIntT(B);
	}
	
	//Easy cases where the sizes are different
	if (((A->size) != (B->size)))
		returnVal = A->size - B->size;
	
	//Compare most significant bunches, then next most significant, etc.
	else
	{
		//Iterate through all bunches, checking for equality
		for (int currIndex = A->size-1; currIndex >= 0; currIndex -= 1)
		{
			if ((A->theInt[currIndex]) < (B->theInt[currIndex]))
			{
				returnVal = -1;
				break;
			}
			
			else if ((A->theInt[currIndex]) > (B->theInt[currIndex]))
			{
				returnVal = 1;
				break;
			}
		}
	}
	
	//Freeing created pointers, if needed
	if (A != AA)
		A = free_BigIntT(A);
	if (B != BB)
		B = free_BigIntT(B);
	
	return returnVal;
}


int add_bunches(const BigIntTP n, int numOfBunches, BigIntTP result)
/** Add numOfBunches bunches to n and store the result in result.
    This function assumes result has been initialised.
		
		Returns 1 on success, 0 otherwise. */
{
	//Making sure result is ready to store what we need it to
	result->size   = n->size + numOfBunches;
	result->theInt = realloc(result->theInt, (result->size)*sizeof(int));
		
	//Copy the rest of n into result
	for (int i = 0; i < n->size; i += 1)
		result->theInt[i+numOfBunches] = n->theInt[i];
	
	//Adding extra bunches
	for (int i = 0; i < numOfBunches; i += 1)
		result->theInt[i] = 0;
	
	return 1;
}


int subtract_bunches(const BigIntTP n, int numOfBunches, BigIntTP result)
/** Subtracts numOfBunches bunches to n and stores the result in result.
    This function assumes result has been initialised.
		
		Returns 1 on success, 0 otherwise. */
{
	//Can't subtract more bunches than what's already in n
	if (numOfBunches > n->size)
		return 0;
	
	//Preparing result to store what it needs to
	//clear_BigIntT(result);
	result->size = n->size - numOfBunches;
	if (result->size == 0)
		result->size = 1;
	result->theInt = realloc(result->theInt, (result->size)*sizeof(int));
	
	//Copy bunches from n to result
	for (int bunch = numOfBunches; bunch < n->size; bunch += 1)
		result->theInt[bunch-numOfBunches] = n->theInt[bunch];
	
	return 1;
}


int multiply_by_ten(BigIntTP n)
/** Multiplies n by 10. Returns 1.

    This function only works when MAXBUNCH is a power of ten. */
{
	int digitToCarry = 0;
	int lastDigitToCarry;
	
	//If we need to give n an extra bunch to hold itself multiplied by 10
	if (n->theInt[n->size-1] >= MAXBUNCH/10)
	{
		resize_BigIntT(n, n->size+1);
		n->theInt[n->size-1] = 0; //Initialising new bunch
	}
	
	for (int bunch = 0; bunch < n->size; bunch += 1)
	{
		lastDigitToCarry = digitToCarry;
		digitToCarry = 10*n->theInt[bunch] / MAXBUNCH;
		n->theInt[bunch] %= MAXBUNCH/10;
		n->theInt[bunch] *= 10;
		n->theInt[bunch] += lastDigitToCarry;
	}
	
	return 1;
}


int divide_by_ten(BigIntTP n)
/** Divides the given n by ten. If n < 10, n becomes 0.
    Returns 1. 
		
		This function only works if MAXBUNCH is a power of 10. */
{
	//Digit that gets dragged to the next most significant bunch
	int digitToCarry = 0;
	int lastDigitToCarry;
	
	for (int bunch = n->size-1; bunch >= 0; bunch -= 1)
	{
		lastDigitToCarry = digitToCarry;
		digitToCarry = n->theInt[bunch] % 10;
		n->theInt[bunch] /= 10;
		n->theInt[bunch] += (MAXBUNCH/10) * lastDigitToCarry;
	}
	
	reduce_BigIntT(n);
	
	return 1;
}


int add_BigIntT(const BigIntTP A, const BigIntTP B, BigIntTP sum)
/** Computes A + B and stores the result in sum. This function assumes
    sum has been initialised.
		Returns 1 on success, 0 otherwise. */
{
	//Just pointers, so we don't need to free them
	BigIntTP smol; //Pointer to smaller of the two arguments
	BigIntTP big;  //        ... bigger of the two arguments
	
	if (compare_BigIntT(A, B) < 0)
	{
		smol = A;
		big  = B;
	}
	else
	{
		smol = B;
		big  = A;
	}
	
	//Properly initialise sum
	if (sum->size < (big->size) + 1)
		resize_BigIntT(sum, (big->size) + 1);
	
	clear_BigIntT(sum);
	
	//Perform addition
	for (int i = 0; i < smol->size; i += 1)
	{
		if (A->theInt[i] + B->theInt[i] + sum->theInt[i] > MAXBUNCH-1)
			sum->theInt[i+1] += 1;
		
		sum->theInt[i] += (A->theInt[i] + B->theInt[i]) % MAXBUNCH; 
		
		if ((A->size <= i) || (B->size <= i) || (sum->size <= i+1))
		{
			printf("Invalid read!!!!! :o\n");
			printf("A = ");
			printi(A);
			printf("\nB = ");
			printi(B);
			printf("\nsum = ");
			printi(sum);
			printf("\ni = %d\n", i);
			printf("A->size = %d\n", A->size);
			printf("B->size = %d\n", B->size);
			printf("sum->size = %d\n", sum->size);
			if (A == smol)
				printf("A is smol\n");
			else
				printf("B is smol\n");
			getchar();
		}
	}
	
	//Add extra bunches that weren't added above
	for (int i = smol->size; i < big->size; i += 1)
	{
		//If we get another carry
		if (sum->theInt[i] + big->theInt[i] > MAXBUNCH-1)
			sum->theInt[i+1] += 1;
		
		sum->theInt[i] += big->theInt[i] % MAXBUNCH;
	}
	
	reduce_BigIntT(sum);

	return 1;
}


int subtract_BigIntT(const BigIntTP subFrom, const BigIntTP toSub, BigIntTP difference)
/** Calculates subFrom - toSub, stores result in difference. 
    Returns 1 on success, 0 on error. 
		
		This function assumes difference has been initialised. */
{
	//Only positive BigIntT structs can be handled currently, so
	// if toSub > subFrom, then return 0.
	if (compare_BigIntT(subFrom, toSub) < 0)
		return 0;
	
	//Make sure difference is the correct size to hold result
	if (difference->size < subFrom->size)
		resize_BigIntT(difference, subFrom->size);
	
	clear_BigIntT(difference);
	
	//Start from least significant bunch, work up
	for (int b = 0; b < toSub->size; b += 1)
	{
		//The - difference->theInt[b] accounts for borrows that lead to other borrows 
		if (subFrom->theInt[b] < toSub->theInt[b] - difference->theInt[b])
		{
			difference->theInt[b] += MAXBUNCH + (subFrom->theInt[b] - toSub->theInt[b]);
			difference->theInt[b+1] -= 1;
		}
		
		else
			difference->theInt[b] += subFrom->theInt[b] - toSub->theInt[b];
	}
	
	//Add extra bunches that weren't involved in the subtraction
	for (int b = toSub->size; b < subFrom->size; b += 1)
	{
		//If a borrow leaked past the original small number
		if (difference->theInt[b] + subFrom->theInt[b] < 0)
		{
			difference->theInt[b+1] -= 1;
			difference->theInt[b] += MAXBUNCH;
		}
		
		else
			difference->theInt[b] += subFrom->theInt[b];
	}
	
	reduce_BigIntT(difference);
	
	return 1;
}


int multiply_BigIntT(const BigIntTP A, const BigIntTP B, BigIntTP product)
/** Multiples A and B, stores the product in product.
    Returns 1 on success, 0 otherwise. */
{
	BigIntTP zero    = empty_BigIntT(1);
	BigIntTP tempLot = empty_BigIntT(1); //Used for holding how much we add at a time
	BigIntTP temp    = empty_BigIntT(1);
	
	int bunchMagnitude;
	int bunchValue;
	
	//If either A or B is zero
	if ((compare_BigIntT(A, zero) == 0) || (compare_BigIntT(B, zero) == 0))
	{
		copy_BigIntT(zero, product);
		zero    = free_BigIntT(zero);
		tempLot = free_BigIntT(tempLot);
		temp    = free_BigIntT(temp);
		return 1;
	}
	
	//Prepare our variables for the multiplication
	clear_BigIntT(product);
	
	//This is required to prevent indexing errors in add_BigPolyT() later in this function
	// Without it, a zero with extra bunches can appear, throwing off bunch indexing
	reduce_BigIntT(product);
	
	copy_BigIntT(B, tempLot);

	//Prepare to add B a bunch of times
	if (A->size > 1)
	{
		add_bunches(tempLot, A->size - 1, temp);
		copy_BigIntT(temp, tempLot);
	}

	
	//Multiply by ten a few times to get the size perfect
	for (int i = 0; i < num_digits(A->theInt[A->size-1])-1; i += 1)
		multiply_by_ten(tempLot);
	
	//tempLot should now be the correct magnitude for adding
	//Now, we can actually start the multiplication (repeated addition)
	for (int bunchCounter = A->size - 1; bunchCounter >= 0; bunchCounter -= 1)
	{
		bunchValue = A->theInt[bunchCounter];
		
		//If we need to calculate bunchMagnitude
		if (bunchCounter == A->size - 1)
		{
			bunchMagnitude = 1;
			
			for (int i = 0; i < num_digits(bunchValue) - 1; i += 1)
				bunchMagnitude *= 10;
		}
		
		//Bunches past the first, so we know their size
		else
			bunchMagnitude = MAXBUNCH / 10;
		
		//Now, extract the digits from the bunchValue, use them to add tempLot
		//(bunchMagnitude > 0) ensures tempLot remains the correct magnitude 
		while ((bunchValue != 0) || (bunchMagnitude > 0))
		{
			for (int i = 0; i < bunchValue / bunchMagnitude; i += 1)
			{
				add_BigIntT(tempLot, product, temp);
				copy_BigIntT(temp, product);
			}
			
			bunchValue %= bunchMagnitude;
			bunchMagnitude /= 10;
			divide_by_ten(tempLot);
		}
	}
	
	zero    = free_BigIntT(zero);
	tempLot = free_BigIntT(tempLot);
	temp    = free_BigIntT(temp);
	
	reduce_BigIntT(product);
	
	return 1;
}


int divide_BigIntT(const BigIntTP toDivide, const BigIntTP divideBy, BigIntTP quotient)
/** Divides the first BigIntT by the second, and stores the
    result in the third BigIntT. This function assumes
		all BigIntT structs have been initialised.
		
		This function discards any remainder.
		
		Returns 1 on success, 0 otherwise. */
{
	int comparison; //Holds comparison values so we don't have to keep recomputing them
	int one[]    = {1};
	int tenArr[] = {10};
	
	//Can't divide by zero in this context
	BigIntTP zero = empty_BigIntT(1);
	if (compare_BigIntT(divideBy, zero) == 0)
	{
		zero = free_BigIntT(zero);
		return 0;
	}
	
	BigIntTP ten = new_BigIntT(tenArr, 1);
	
	//Holds a divisor of reasonable size so that we
	// aren't subtracting a baby number all day
	BigIntTP tempDivisor = empty_BigIntT(1);
	copy_BigIntT(divideBy, tempDivisor);
	
	//Holds the size of our tempDivisor, to properly
	// keep track of how many times we've subtracted it
	// from toDivide
	BigIntTP divisorMagnitude = new_BigIntT(one, 1);
	
	//For any miscellaneous operations we need
	BigIntTP temp  = empty_BigIntT(1);
	BigIntTP temp2 = empty_BigIntT(1);
	
	//Don't worry about the size of quotient,
	// it sorts itself out once it gets chucked into a function
	clear_BigIntT(quotient);
	//quotient->size = 1; //I don't think I need this line
	
	//If there's actually some meaningful division to do
	// (if toDivide > divideBy)
	if (compare_BigIntT(toDivide, divideBy) >= 0)
	{
		//Find a reasonable divisor to start with
		comparison = compare_BigIntT(toDivide, tempDivisor);
		if (comparison > 0)
		{
			add_bunches(divisorMagnitude, comparison-1, temp);
			copy_BigIntT(temp, divisorMagnitude); 
			
			add_bunches(tempDivisor, comparison-1, temp);
			copy_BigIntT(temp, tempDivisor);
		}
		
		/*printf("Size of toDivide's MSB: %d\n", num_digits(toDivide->theInt[toDivide->size-1]));
		printf("Size of tempDivisor's MSB: %d\n", num_digits(tempDivisor->theInt[toDivide->size-1]));
		
		printf("toDivide: ");
		printi(toDivide);
		printf("\ntempDivisor: ");
		printi(tempDivisor);
		printf("\n"); */
		
		//Get tempDivisor as close to toDivide as possible
		//Shouldn't use a for loop here since numbers with the name number of
		// digits in their most significant bunch can be different in size
		while (compare_BigIntT(tempDivisor, toDivide) <= 0)
		{
			multiply_by_ten(tempDivisor);
			multiply_by_ten(divisorMagnitude);
		}
		
		divide_by_ten(tempDivisor);
		divide_by_ten(divisorMagnitude);
		
		//Now, we perform the actual division (which is just repeated subtraction)
		copy_BigIntT(toDivide, temp);
		
		while (compare_BigIntT(temp, divideBy) >= 0)
		{
			//The if statement prevents from falsely subtracting a bigger number from a smaller one
			if (subtract_BigIntT(temp, tempDivisor, temp2))
			{
				copy_BigIntT(temp2, temp);
				
				//Keep track of how many times we've subtracted divideBy
				add_BigIntT(quotient, divisorMagnitude, temp2);
				copy_BigIntT(temp2, quotient);
			}
			
			else
			{
				//Reduce tempDivisor by as many bunches as needed
				comparison = compare_BigIntT(tempDivisor, temp);
				if ((comparison > 0) && (compare_BigIntT(zero, temp) != 0))
				{
					//Subtract off bunches
					subtract_bunches(tempDivisor, comparison-1, temp2);
					copy_BigIntT(temp2, tempDivisor);
					
					subtract_bunches(divisorMagnitude, comparison-1, temp2);
					copy_BigIntT(temp2, divisorMagnitude);
					
					//Now, just divide by ten until your tempDivisor is small enough
					while (compare_BigIntT(tempDivisor, temp) > 0)
					{
						divide_by_ten(tempDivisor);
						divide_by_ten(divisorMagnitude);
					}
				}
			}
		}
	}
	
	tempDivisor      = free_BigIntT(tempDivisor);
	divisorMagnitude = free_BigIntT(divisorMagnitude);
	temp             = free_BigIntT(temp);
	temp2            = free_BigIntT(temp2);
	zero             = free_BigIntT(zero);
	ten              = free_BigIntT(ten);
	
	reduce_BigIntT(quotient);
	return 1;
}


int mod_BigIntT(const BigIntTP toMod, const BigIntTP modulus, BigIntTP residue)
/** Calculates toMod % modulus and stores it in residue.
    Returns 1 on success, 0 otherwise. */
{
	//This implementation is terrible. I apologise.
	
	BigIntTP tempModulus = empty_BigIntT(modulus->size);
	BigIntTP temp        = empty_BigIntT(toMod->size); 
	
	//Holds comparison values so we don't need to keep recomputing the value
	int comparison; 
	
	copy_BigIntT(modulus, tempModulus);
	copy_BigIntT(toMod, residue);
	
	//If toMod >> modulus, it'll be very inefficient to
	// subtract singular multiples of modulus.
	// Rather, we can find the greatest value of 
	// modulus*(MAXBUNCH^n)(10^m) that's smaller than toMod and
	// use that. Then, we just find the next smallest one and
	// repeat until we've reduced our number.
	
	//Find a reasonable multiple of the modulus
	// to start with
	comparison = compare_BigIntT(toMod, tempModulus);
	
	if (comparison > 0)
	{
		add_bunches(modulus, comparison-1, tempModulus);
		
		while (compare_BigIntT(tempModulus, toMod) <= 0)
			multiply_by_ten(tempModulus);

		divide_by_ten(tempModulus);
	}
	
	while (compare_BigIntT(modulus, residue) <= 0)
	{
		subtract_BigIntT(residue, tempModulus, temp);
		copy_BigIntT(temp, residue);
		
		//If our residue is less than the modulus,
		// shrink the modulus until it's less than residue again.
		while (compare_BigIntT(tempModulus, residue) > 0)
			divide_by_ten(tempModulus);
	}
	
	temp        = free_BigIntT(temp);
	tempModulus = free_BigIntT(tempModulus);
	reduce_BigIntT(residue);
	
	return 1;
}


BigIntTP* divisors_of_BigIntT(const BigIntTP toFactor)
/** Returns an array of every possible factor of the given BigIntT.
    The first number says how many factors are contained in the array. */
{
	BigIntTP* factorList = malloc(sizeof(BigIntTP));
	
	 //This allows us to cut back on the number of unnecessary factors we check
	BigIntTP half = empty_BigIntT(1);
	BigIntTP tempFactor = empty_BigIntT(1);
	BigIntTP temp = empty_BigIntT(1);
	
	int oneArr[1] = {1};
	int twoArr[1] = {2};
	BigIntTP zero = empty_BigIntT(1);
	BigIntTP one  = new_BigIntT(oneArr, 1);
	BigIntTP two  = new_BigIntT(twoArr, 1);
	
	factorList[0] = empty_BigIntT(1);
	divide_BigIntT(toFactor, two, half);
	copy_BigIntT(one, tempFactor);
	
	int size = 0;
	
	//Search up until we get to half of toFactor
	while (compare_BigIntT(tempFactor, half) <= 0)
	{
		//Check to see if tempFactor divides toFactor
		mod_BigIntT(toFactor, tempFactor, temp);
		if (compare_BigIntT(zero, temp) == 0)
		{
			add_BigIntT(factorList[0], one, temp);
			copy_BigIntT(temp, factorList[0]);
			size += 1;
			
			factorList = realloc(factorList, (size+1)*sizeof(BigIntTP));
			factorList[size] = empty_BigIntT(1);
			copy_BigIntT(tempFactor, factorList[size]);
		}
		
		add_BigIntT(tempFactor, one, temp);
		copy_BigIntT(temp, tempFactor);
	}
	
	half = free_BigIntT(half);
	temp = free_BigIntT(temp);
	zero = free_BigIntT(zero);
	one  = free_BigIntT(one);
	two  = free_BigIntT(two);
	
	tempFactor = free_BigIntT(tempFactor);
	
	return factorList;
}
