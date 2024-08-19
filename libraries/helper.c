
/*
A collection of miscellaneous helper functions. 

August 3, 2022
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/* Should've made this long ago. */
//At some point, it may be worthwhile to rewrite this
// as a hashmap instead of an ordered pair.
typedef enum dicttype {STR, STRARR, INT} DictionaryTypeT;

typedef union dictval
{
	char* str;
	char** strArr;
	int intNum;
}
DictionaryValueT, *DictionaryValueTP;

typedef struct dict
{
	//Says what kind of data this dictionary holds
	DictionaryTypeT type;
	
	char* key;
	DictionaryValueTP value;
	DictionaryValueTP auxValue;
}
DictionaryT, *DictionaryTP;


DictionaryTP free_DictionaryT(DictionaryTP d)
{
	free(d->key);
	
	switch (d->type)
	{
		case (STR):
		{
			free(d->value->str);
			break;
		}
		
		case (STRARR):
		{
			/*
			for (int i = 0; i < d->auxValue->intNum; i += 1)
				free(d->value->strArr[i]);
			free(d->value->strArr);
			*/
			//In this case, we're assuming the user will handle freeing the array
			break;
		}
		
		default:
		{
			fprintf(stderr, "Freeing dictionary of type %d is not yet implemented.\n", d->type);
		}
	}
	
	free(d->value);
	free(d->auxValue);
	
	d->key = NULL;
	d->value = NULL;
	d->auxValue = NULL;
	
	return NULL;
}


DictionaryTP new_DictionaryT(const char* key, const void* value, const char* type)
/** Creates a new DictionaryT object, returns a pointer to it. */
{
	DictionaryTP d = malloc(sizeof(DictionaryT));
	d->key = malloc((strlen(key)+1)*sizeof(char));
	d->value = malloc(sizeof(DictionaryValueT));
	
	strcpy(d->key, key);
	
	if (!strcmp(type, "STR"))
	{
		d->type = STR;
		d->value->str = malloc((strlen((char*)value)+1)*sizeof(char));
		strcpy(d->value->str, (char*)value);
	}
	
	else if (!strcmp(type, "STRARR"))
	{
		d->type = STRARR;
		d->value->strArr = (char**)value;
	}
	
	else if (!strcmp(type, "INT"))
	{
		d->type = INT;
		d->value->intNum = *((int*)value);
	}
	
	return d;
}


DictionaryTP search_DictionaryT_array(const DictionaryTP* arr, int arrLen, const char* key)
/** Sifts through a DictionaryTP array to find the first instance
    of a dictionary with the given key. Returns found DictionaryT 
		on success, NULL otherwise. */
{
	if (arr == NULL)
		return NULL;
	
	//Yippee! Linear time complexity!
	for (int i = 0; i < arrLen; i += 1)
		if (!strcmp(arr[i]->key, key))
			return arr[i];
		
	return NULL;
}


void* value_of_DictionaryT(const DictionaryTP d)
/** Returns the value held in the dictionary. */
{
	switch (d->type)
	{
		case (INT):
		{
			return &(d->value->intNum);
		}
		
		case (STR):
		{
			return d->value->str;
		}
		
		case (STRARR):
		{
			return d->value->strArr;
		}
		
		default:
		{
			return NULL;
		}
	}
}


void* aux_of_DictionaryT(const DictionaryTP d)
/** Returns the auxillary value of a dictionary. */
{
	switch (d->type)
	{
		case (STRARR):
		{
			return &(d->auxValue->intNum);
		}
		
		default:
		{
			return NULL;
		}
	}
}


int num_digits(int num)
/** Calculuates how many digits an integer has and returns
    that number. The function assumes the number given is
		positive. */
{
	int count = 0;
	while (1)
	{
		num /= 10;
		count += 1;
		
		if (num == 0)
			return count;
	}
}

void append_int(char* str, int toApp)
/** Appends the given integer to the given string. */
{
	int tempInt;
	int tempPower = 1;
	
	char* tempDigit = malloc(2*sizeof(char));
	tempDigit[1] = '\0';
	
	while (toApp / tempPower >= 10)
		tempPower *= 10;
	
	tempInt = toApp;
	while (tempPower > 0)
	{
		//Converting numbers to characters
		tempDigit[0] = (tempInt / tempPower) + 48;
		strcat(str, tempDigit);
		tempInt %= tempPower;
		tempPower /= 10;
	}
	
	free(tempDigit);
	tempDigit = NULL;
}

void quick_sort(int* list, int size)
/** Sorts list from least to greatest using the
    quick sort algorithm. */
{
	//en.wikipedia.org/wiki/Quicksort#Algorithm
	//Thanks, Wikipedia, for helping me realise
	// I can do this algorithm in-place rather
	// than have to allocate a million different arrays
	
	if (size == 1)
		return;
	
	int pivot = 0;    //Assume the first element is the Sedgewick
	int hold;         //For intermediate calculation
	int bound = size; //For sectioning off elements >= pivot
	
	//Choose the median of the first, middle, and last
	// element as the pivot, as Wikipedia says that's
	// faster.
	if (size > 2)
	{
		//If the middle element is the Sedgewick
		if (((list[size/2] <= list[0]) && (list[size/2] >= list[size-1])) || 
		    ((list[size/2] >= list[0]) && (list[size/2] <= list[size-1])))
		{
			hold = list[size/2];
			list[size/2] = list[0];
			list[0] = hold;
		}
		
		//If the last element is the Sedgewick
		else if (((list[size-1] >= list[0]) && (list[size-1] <= list[size/2])) ||
		         ((list[size-1] <= list[0]) && (list[size-1] >= list[size/2])))
		{
			hold = list[size-1];
			list[size-1] = list[0];
			list[0] = hold;
		}
	}
	
	//Swap elements until pivot divides list into two
	// elements less than it and elements greater than
	// or equal to it
	for (int pSwap = 1; pSwap < bound; pSwap += 1)
	{
		//Put smaller elements to the left of the pivot
		if (list[pivot] > list[pSwap])
		{
			hold = list[pivot];
			list[pivot] = list[pSwap];
			list[pSwap] = hold;
			pivot = pSwap;
		}
		
		//Move greater elements to the far-right of the list
		else
		{
			hold = list[bound-1];
			list[bound-1] = list[pSwap];
			list[pSwap] = hold;
			pSwap -= 1;
			bound -= 1;
		}
	}
	
	//Now, call the same algorithm again (if we need to)
	// to sort the two chunks of the list we've just made
	if (pivot > 0)
		quick_sort(list, pivot+1);
	if (size-pivot-1 > 1)
		quick_sort(list+(pivot+1), size-pivot-1);
}
