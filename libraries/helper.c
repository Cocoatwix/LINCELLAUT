
/*
A collection of miscellaneous helper functions. 

August 3, 2022
*/

#include <stdlib.h>
#include <string.h>

void append_int(char* str, int toApp)
/** Appends the given integer to the given string. */
{
	int tempInt;
	int tempPower = 1;
	
	char* tempDigit = malloc(2*sizeof(char));
	tempDigit[1] = '\0';
	
	while (1)
	{
		if (toApp / tempPower >= 10)
			tempPower *= 10;
		else
			break;
	}
	
	tempInt = toApp;
	while (tempPower > 0)
	{
		//Converting numbers to characters
		tempDigit[0] = (tempInt / tempPower) + 48;
		strcat(str, tempDigit);
		tempInt %= 10;
		tempPower /= 10;
	}
	
	free(tempDigit);
	tempDigit = NULL;
}