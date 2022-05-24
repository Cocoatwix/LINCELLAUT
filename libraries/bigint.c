
/*
A super basic arbitrary precision library for
integers. Will be used to compute big examples 
of LCA systems if needed.

May 24, 2022
*/

#include <stdlib.h>
#include <stdio.h>

//Using <= 46340 allows multiplication to work
//Using <= 1073741823 allows addition to work, but not
// multiplication or correct printing of the numbers
const int MAXBUNCH = 46340;

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
	BigIntTP newBigInt = malloc(sizeof(BigIntTP));;
	newBigInt->size    = zeros;
	newBigInt->theInt  = calloc(zeros, sizeof(int));
	
	return newBigInt;
}


void printi(BigIntTP n)
/** Prints a BigIntT struct to the screen as a normal int would. */
{
	for (int i = n->size-1; i >= 0; i -= 1)
		printf("%d", n->theInt[i]);
}


int reduce_BigIntT(BigIntTP toReduce)
/** Removes any unnecessary bunches from toReduce
    and sets its size accordingly. Returns 1 on success,
		0 otherwise. */
{
	int counter = 0;
	
	//Iterate backwards through the number, seeing
	// if there's any bunches we can remove
	for (int i = toReduce->size-1; i >= 0; i -= 1)
	{
		if (toReduce->theInt[i] != 0)
			break;
		
		counter += 1;
	}
	
	//Get rid of extra bunches, update size
	toReduce->size -= counter;
	toReduce->theInt = realloc(toReduce->theInt, toReduce->size);
	
	return 1;
}


int copy_BigIntT(BigIntTP const setTo, BigIntTP toSet)
/** Copies setTo to toSet and returns 1 on completion.
    Returns 0 on error. */
{
	//Clear out toSet before doing anything
	for (int i = 0; i < toSet->size; i += 1)
		toSet->theInt[i] = 0;
	
	//Reallocate toSet's memory appropriately
	if (toSet->size != setTo->size)
		toSet->theInt = realloc(toSet->theInt, setTo->size);
	
	//Copy contents
	for (int i = 0; i < setTo->size; i += 1)
		toSet->theInt[i] = setTo->theInt[i];
	
	toSet->size = setTo->size;
	
	return 1;
}


int compare_BigIntT(BigIntTP const A, BigIntTP const B)
/** Returns -1 if A < B, 0 if A == B, 1 if A > B */
{
	//Easy cases where the sizes are different
	if (A->size < B->size)
		return -1;
	
	else if (A->size > B->size)
		return 1;
	
	//Compare most significant bunches
	else
	{
		if (A->theInt[A->size] < B->theInt[B->size])
			return -1;
		
		else if (A->theInt[A->size] == B->theInt[B->size])
			return 0;
		
		else
			return 1;
	}
}


int subtract_BigIntT(BigIntTP const subFrom, BigIntTP const toSub, BigIntTP difference)
/** Calculates subFrom - toSub, stores result in difference. 
    Returns 1 on success, 0 on error. 
		
		This function assumes difference has been properly initialised to 0. */
{
	//Only positive BigIntT structs can be handled currently, so
	// if toSub > subFrom, then return 0.
	if (compare_BigIntT(subFrom, toSub) < 0)
		return 0;
	
	//Start from least significant bunch, work up
	for (int b = 0; b < toSub->size; b += 1)
	{
		if (subFrom->theInt[b] < toSub->theInt[b])
		{
			difference->theInt[b] += MAXBUNCH + (subFrom->theInt[b] - toSub->theInt[b]);
			difference->theInt[b+1] -= 1;
		}
		
		else
			difference->theInt[b] += subFrom->theInt[b] - toSub->theInt[b];
	}
	
	reduce_BigIntT(difference);
	
	return 1;
}


int mod_BigIntT(BigIntTP const toMod, BigIntTP const modulus, BigIntTP residue)
/** Calculates toMod % modulus and stores it in residue.
    Returns 1 on success, 0 otherwise. */
{
	//This implementation is terrible. I apologise.
	
	copy_BigIntT(toMod, residue);
	BigIntTP temp;
	
	//Making sure temp is the right size to hold result of subtraction
	if (modulus->size > toMod->size)
		temp = empty_BigIntT(modulus->size);
	else
		temp = empty_BigIntT(toMod->size);
	
	while (compare_BigIntT(modulus, residue) <= 0)
	{
		subtract_BigIntT(residue, modulus, temp);
		copy_BigIntT(temp, residue);
	}
	
	reduce_BigIntT(residue);
	return 1;
}