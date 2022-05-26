
/*
A super basic arbitrary precision library for
integers. Will be used to compute big examples 
of LCA systems if needed.

May 24, 2022
*/

#include <stdlib.h>
#include <stdio.h>

#include "../headers/linalg.h" //num_digits()

//Using <= 1073741823 allows addition to work + above
//Using <= 46340 allows multiplication to work
//Using == 10000 allows for proper printing + above
const int MAXBUNCH = 10000;

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
	free(n->theInt);
	n->theInt = NULL;
	free(n);
	return NULL;
}


BigIntTP new_BigIntT(int* const initNum, int len)
/** Returns a pointer to an initialised
    BigIntT struct. Returns NULL on error. */
{
	BigIntTP newBigInt = malloc(sizeof(BigIntTP));
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
	BigIntTP newBigInt = malloc(sizeof(BigIntTP));
	newBigInt->size    = zeros;
	newBigInt->theInt  = calloc(zeros, sizeof(int));
	
	return newBigInt;
}


void printi(BigIntTP n)
/** Prints a BigIntT struct to the screen as a normal int would. */
{
	if (MAXBUNCH == 10000)
	{
		for (int i = n->size-1; i >= 0; i -= 1)
		{
			//Zero padding
			if (i != n->size-1)
				for (int d = 0; d < 4 - num_digits(n->theInt[i]); d += 1)
					printf("0");
			
			printf("%d", n->theInt[i]);
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
	toResize->size = newSize;
	toResize->theInt = realloc(toResize->theInt, (toResize->size)*sizeof(int));
	
	return 1;
}


int copy_BigIntT(BigIntTP const setTo, BigIntTP toSet)
/** Copies setTo to toSet and returns 1 on completion.
    Returns 0 on error. */
{
	//Reallocate toSet's memory appropriately
	if (toSet->size != setTo->size)
	{
		toSet->size = setTo->size;
		toSet->theInt = realloc(toSet->theInt, (toSet->size)*sizeof(int));
	}
	
	//Clear out toSet before doing anything
	for (int i = 0; i < toSet->size; i += 1)
		toSet->theInt[i] = 0;
	
	//Copy contents
	for (int i = 0; i < setTo->size; i += 1)
		toSet->theInt[i] = setTo->theInt[i];
	
	return 1;
}


int compare_BigIntT(BigIntTP const AA, BigIntTP const BB)
/** Returns -1 if A < B, 0 if A == B, 1 if A > B */
{
	//Making sure no weirdness occurs with size comparisons
	BigIntTP A = empty_BigIntT(1);
	BigIntTP B = empty_BigIntT(1);
	copy_BigIntT(AA, A);
	copy_BigIntT(BB, B);
	reduce_BigIntT(A);
	reduce_BigIntT(B);
	
	int returnVal;
	
	//Easy cases where the sizes are different
	if ((A->size) < (B->size))
		returnVal = -1;
	
	else if ((A->size) > (B->size))
		returnVal = 1;
	
	//Compare most significant bunches
	else
	{
		if ((A->theInt[A->size-1]) < (B->theInt[B->size-1]))
			returnVal = -1;
		
		else if ((A->theInt[A->size-1]) == (B->theInt[B->size-1]))
			returnVal = 0;
		
		else
			returnVal = 1;
	}
	
	A = free_BigIntT(A);
	B = free_BigIntT(B);
	
	return returnVal;
}


int add_bunches(BigIntTP const n, int numOfBunches, BigIntTP result)
/** Add numOfBunches bunches to n and store the result in result.
    This function assumes result has been initialised.
		
		Returns 1 on success, 0 otherwise. */
{
	//Making sure result is ready to store what we need it to
	clear_BigIntT(result);
	result->size = n->size + numOfBunches;
	result->theInt = realloc(result->theInt, (result->size)*sizeof(int));
	
	//Adding extra bunches
	for (int i = 0; i < numOfBunches; i += 1)
		result->theInt[i] = 0;
	
	//Copy the rest of n into result
	for (int i = 0; i < n->size; i += 1)
		result->theInt[i+numOfBunches] = n->theInt[i];
	
	return 1;
}


int subtract_bunches(BigIntTP const n, int numOfBunches, BigIntTP result)
/** Subtracts numOfBunches bunches to n and stores the result in result.
    This function assumes result has been initialised.
		
		Returns 1 on success, 0 otherwise. */
{
	//Can't subtract more bunches than what's already in n
	if (numOfBunches > n->size)
		return 0;
	
	//Preparing result to store what it needs to
	clear_BigIntT(result);
	result->size = (n->size) - numOfBunches;
	if (result->size == 0)
		result->size = 1;
	result->theInt = realloc(result->theInt, (result->size)*sizeof(int));
	
	//Copy bunches from n to result
	for (int bunch = numOfBunches; bunch < n->size; bunch += 1)
		result->theInt[bunch-numOfBunches] = n->theInt[bunch];
	
	return 1;
}


int add_BigIntT(BigIntTP const A, BigIntTP const B, BigIntTP sum)
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
	
	for (int i = 0; i < smol->size; i += 1)
	{
		if (A->theInt[i] + B->theInt[i] + sum->theInt[i] > MAXBUNCH-1)
			sum->theInt[i+1] += 1;
		
		sum->theInt[i] += (A->theInt[i] + B->theInt[i]) % MAXBUNCH; 
	}
	
	//Add extra bunches that weren't added above
	for (int i = smol->size; i < big->size; i += 1)
	{
		//If we get another carry
		if (sum->theInt[i] + big->theInt[i] > MAXBUNCH-1)
			sum->theInt[i+1] += 1;
		
		sum->theInt[i] += big->theInt[i] % MAXBUNCH;
	}
	
	/* printf("Result of sum: ");
	printi(sum);
	printf("\n"); */
	reduce_BigIntT(sum);
	
	
	return 1;
}


int subtract_BigIntT(BigIntTP const subFrom, BigIntTP const toSub, BigIntTP difference)
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


int divide_BigIntT(BigIntTP const toDivide, BigIntTP const divideBy, BigIntTP quotient)
/** Divides the first BigIntT by the second, and stores the
    result in the third BigIntT. This function assumes
		all BigIntT structs have been initialised.
		
		This function discards any remainder.
		
		Returns 1 on success, 0 otherwise. */
{
	int one[] = {1};
	
	//Can't divide by zero in this context
	BigIntTP zero = empty_BigIntT(1);
	if (compare_BigIntT(divideBy, zero) == 0)
	{
		zero = free_BigIntT(zero);
		return 0;
	}
	
	//Holds a divisor of reasonable size so that we
	// aren't subtracting a baby number all day
	BigIntTP tempDivisor = empty_BigIntT(1);
	copy_BigIntT(divideBy, tempDivisor);
	
	//Holds the size of our tempDivisor, to properly
	// keep track of how many times we've subtracted it
	// from toDivide
	BigIntTP divisorMagnitude = new_BigIntT(one, 1);
	
	//For any miscellaneous operations we need this for
	BigIntTP temp  = empty_BigIntT(1);
	BigIntTP temp2 = empty_BigIntT(1);
	
	//Don;t worry about the size of quotient,
	// it sorts itself out once it gets chucked into a function
	clear_BigIntT(quotient);
	quotient->size = 1;
	
	//If there's actually some meaningful division to do
	// (if toDivide > divideBy)
	if (compare_BigIntT(toDivide, divideBy) >= 0)
	{
		//Find a reasonable divisor to start with
		add_bunches(tempDivisor, 1, temp);
		copy_BigIntT(temp, tempDivisor);
			
		while (compare_BigIntT(toDivide, tempDivisor) > 0)
		{
			add_bunches(divisorMagnitude, 1, temp);
			copy_BigIntT(temp, divisorMagnitude); 
			
			add_bunches(tempDivisor, 1, temp);
			copy_BigIntT(temp, tempDivisor);
		}
		
		//We'll now have a tempDivisor which is one bunch too big
		//Make sure to subtract it!
		
		//If you truly want arbitrary precision here, you can't directly
		// use divisorMagnitude's size, as that'll max out at 2147483647.
		// Rather, you need to manually subtract through divisorMagnitude with
		// a for-loop to add the correct number of bunches here.
		add_bunches(divideBy, (divisorMagnitude->size)-1, tempDivisor);
		
		//Now, we can perform the actual division (which is just repeated
		// subtraction here).
		copy_BigIntT(toDivide, temp);
		
		while (compare_BigIntT(temp, divideBy) >= 0)
		{
			subtract_BigIntT(temp, tempDivisor, temp2);
			copy_BigIntT(temp2, temp);
			
			//Keep track of how many times we've subtracted divideBy
			add_BigIntT(quotient, divisorMagnitude, temp2);
			copy_BigIntT(temp2, quotient);
			
			//Reduce tempDivisor by one bunch
			if (compare_BigIntT(temp, tempDivisor) < 0)
			{
				subtract_bunches(tempDivisor, 1, temp2);
				copy_BigIntT(temp2, tempDivisor);
				
				subtract_bunches(divisorMagnitude, 1, temp2);
				copy_BigIntT(temp2, divisorMagnitude);
			}
		}
	}
	
	tempDivisor      = free_BigIntT(tempDivisor);
	divisorMagnitude = free_BigIntT(divisorMagnitude);
	temp             = free_BigIntT(temp);
	temp2            = free_BigIntT(temp2);
	zero             = free_BigIntT(zero);
	
	reduce_BigIntT(quotient);
	return 1;
}


int mod_BigIntT(BigIntTP const toMod, BigIntTP const modulus, BigIntTP residue)
/** Calculates toMod % modulus and stores it in residue.
    Returns 1 on success, 0 otherwise. */
{
	//This implementation is terrible. I apologise.
	BigIntTP tempModulus = empty_BigIntT(modulus->size);
	
	//This should be rewritten to use a BigIntT here
	int modCounter = 0;
	
	copy_BigIntT(modulus, tempModulus);
	
	copy_BigIntT(toMod, residue);
	
	BigIntTP temp = empty_BigIntT(toMod->size); 
	
	//If toMod >> modulus, it'll be very inefficient to
	// subtract singular multiples of modulus.
	// Rather, we can find the greatest value of 
	// modulus*(10000^n) that's smaller than toMod and
	// use that. Then, we just find the next smallest one and
	// repeat until we've reduced our number.
	
	//Find a reasonable multiple of the modulus
	// to start with
	while (compare_BigIntT(toMod, tempModulus) > 0)
	{
		modCounter += 1;
		add_bunches(modulus, modCounter, tempModulus);
	}
	
	modCounter -= 1;
	if (modCounter >= 0)
		add_bunches(modulus, modCounter, tempModulus);
	else
	{
		copy_BigIntT(modulus, tempModulus);
		modCounter = 0;
	}
	
	
	while (compare_BigIntT(modulus, residue) <= 0)
	{
		subtract_BigIntT(residue, tempModulus, temp);
		copy_BigIntT(temp, residue);
		
		//If our residue is less than the modulus,
		// shrink the modulus to the next power of MAXBUNCH down.
		if (compare_BigIntT(residue, tempModulus) < 0)
		{
			if (modCounter > 0)
			{
				modCounter -= 1;
				add_bunches(modulus, modCounter, tempModulus);
			}
		}
	}
	
	temp        = free_BigIntT(temp);
	tempModulus = free_BigIntT(tempModulus);
	reduce_BigIntT(residue);
	
	return 1;
}