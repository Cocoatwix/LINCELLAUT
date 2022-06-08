
/*
A super basic arbitrary precision library for
integers. Will be used to compute big examples 
of LCA systems if needed.

May 24, 2022
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "../headers/linalg.h" //num_digits()

//How big each bunch in a BigIntT can be
const int MAXBUNCH = 100000000;

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


int size(BigIntTP n)
/** Returns the size of the BigIntT passed. */
{
	return n->size;
}


void printi(BigIntTP n)
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
		
		for (int i = n->size-1; i >= 0; i -= 1)
		{
			//Zero padding
			for (int d = 0; d < power - num_digits(n->theInt[i]); d += 1)
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
	
	//Copy contents
	for (int i = 0; i < toSet->size; i += 1)
		toSet->theInt[i] = setTo->theInt[i];
	
	return 1;
}


int compare_BigIntT(BigIntTP const AA, BigIntTP const BB)
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


int add_bunches(BigIntTP const n, int numOfBunches, BigIntTP result)
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


int subtract_bunches(BigIntTP const n, int numOfBunches, BigIntTP result)
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
	
	//Perform addition
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


int multiply_BigIntT(BigIntTP const A, BigIntTP const B, BigIntTP product)
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
		return 1;
	}
	
	//Prepare our variables for the multiplication
	clear_BigIntT(product);
	copy_BigIntT(B, tempLot);
	
	//Prepare to add B a bunch of times
	if (A->size > 1)
	{
		add_bunches(tempLot, A->size - 1, temp);
		copy_BigIntT(temp, tempLot);
	}
	
	//Multiply by ten a few times to get the size perfect
	for (int i = 0; i < num_digits(A->theInt[A->size-1]) - 1; i += 1)
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
		while (bunchValue != 0)
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
	
	return 1;
}


int divide_BigIntT(BigIntTP const toDivide, BigIntTP const divideBy, BigIntTP quotient)
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
	quotient->size = 1;
	
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
		
		//Get tempDivisor as close to toDivide as possible
		while (compare_BigIntT(tempDivisor, toDivide) <= 0)
		{
			multiply_by_ten(tempDivisor);
			multiply_by_ten(divisorMagnitude);
		}
		
		divide_by_ten(tempDivisor);
		divide_by_ten(divisorMagnitude);
		
		// If you truly want arbitrary precision here, you can't directly
		// use divisorMagnitude's size, as that'll max out at 2147483647.
		// Rather, you need to manually subtract through divisorMagnitude with
		// a for-loop to add the correct number of bunches here.
		// add_bunches(divideBy, (divisorMagnitude->size)-1, tempDivisor);
		
		//Now, we can perform the actual division (which is just repeated
		// subtraction here).
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


int mod_BigIntT(BigIntTP const toMod, BigIntTP const modulus, BigIntTP residue)
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
	// modulus*(MAXBUNCH^n) that's smaller than toMod and
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