
/*
A collection of miscellaneous helper functions. 

August 3, 2022
*/

#include <stdlib.h>
#include <string.h>

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