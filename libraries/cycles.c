
#include <stdlib.h>
#include <stdio.h>

#include <string.h>

#include "../headers/linalg.h" //Letting us use matrices

#define FREE(v) free(v); v = NULL

//Might need to write a "free" function later on
//May want to store the vector that's in the cycle as well
typedef struct cycleinfo
/** Struct to hold any important info about a structure. */
{
	int omega; //Cycle length
	int tau;   //Transient length
	int stoppingTime;
	
	IntMatrixTP inCycle; //Stores a vector that's in the cycle
}
CycleInfoT, *CycleInfoTP;


CycleInfoTP free_CycleInfoT(CycleInfoTP c)
/** Frees memory used by the given CycleInfo struct.
    Returns NULL. */
{
	free_IntMatrixT(c->inCycle);
	return NULL;
}


//Not fully implemented yet
void printcycle(CycleInfoTP c)
/** Prints info about a cycle to the screen. */
{
	printf("Floyd's Algorithm stopping time: %d\n", c->stoppingTime);
	printf("Cycle length: %d\n", c->omega);
	
	//-1 is used when we don't know tau
	if (c->tau != -1)
		printf("Transient region length: %d\n", c->tau);
	else
		printf("Transient region length: ??\n");
}


IntMatrixTP iterate(IntMatrixTP F, IntMatrixTP s_0, int modulus, int iterations)
/** Iterates a given vector under some update rule n times.
    Returns a pointer to the result of the iteration. */
{
	IntMatrixTP s_o = new_IntMatrixT(rows(s_0), cols(s_0));
	IntMatrixTP s_e = new_IntMatrixT(rows(s_0), cols(s_0));
	
	copy_IntMatrixT(s_0, s_o);
	
	//The actual iteration
	for (int i = 0; i < iterations; i += 1)
	{
		if (i % 2 == 0)
		{
			mat_mul(F, s_o, s_e);
			modm(s_e, modulus);
		}
		
		else
		{
			mat_mul(F, s_e, s_o);
			modm(s_o, modulus);
		}
	}
	
	//Depending on which vector our result is in, return the right vector
	if (iterations % 2 == 0)
	{
		s_e = free_IntMatrixT(s_e);
		return s_o;
	}
	
	else
	{
		s_o = free_IntMatrixT(s_o);
		return s_e;
	}
}


//UNFINISHED
CycleInfoTP floyd(IntMatrixTP F, IntMatrixTP s_0, int modulus)
/** Uses Floyd's Cycle Detection Algorithm to calculate the 
    transient length and cycle length of some given set up.
		s_0 is the starting vector. M is the update rule.
		
		Returns a CycleInfoTP that points to the data calculated
		about the cycle. Returns NULL on error. */
{
	if (rows(F) != cols(F))
		return NULL;
	
	int stoppingTime = 0;
	IntMatrixTP Finv = inverse(F, modulus);
	
	CycleInfoTP info = malloc(sizeof(CycleInfoT));
	info->omega = 0;
	info->tau = -1;
	info->inCycle = new_IntMatrixT(rows(F), 1);
	
	//If Finv exists, then no transient regions can exist
	if (Finv != NULL)
		info->tau = 0;
	
	//Initalising x and y to be s_0
	//We need to switch between the two versions to store data w/o overwriting
	IntMatrixTP x_1 = new_IntMatrixT(rows(s_0), cols(s_0));
	IntMatrixTP x_2 = new_IntMatrixT(rows(s_0), cols(s_0));
	IntMatrixTP y_1 = new_IntMatrixT(rows(s_0), cols(s_0)); 
	IntMatrixTP y_2 = new_IntMatrixT(rows(s_0), cols(s_0));
	
	copy_IntMatrixT(s_0, x_1);
	copy_IntMatrixT(s_0, y_1);
	
	//Iterate until x and y are the same
	do
	{
		mat_mul(F, x_1, x_2);
		modm(x_2, modulus);
		copy_IntMatrixT(x_2, x_1);
		
		//Iterate y twice
		mat_mul(F, y_1, y_2);
		modm(y_2, modulus);
		mat_mul(F, y_2, y_1);
		modm(y_1, modulus);
		
		stoppingTime += 1;
	}
	while (!compare_IntMatrixT(x_1, y_1));
	
	info->stoppingTime = stoppingTime;
	copy_IntMatrixT(x_1, info->inCycle);
	
	//Now we have a vector that's confirmed to be in a cycle
	//Now, we determine the cycle length
	
	//If St > 2*L, then w = St.
	if (stoppingTime > 2*rows(F))
		info->omega = stoppingTime;
	
	//If above inequality doesn't hold, calculate omega manually
	else
	{
		do
		{
			mat_mul(F, x_1, x_2);
			modm(x_2, modulus);
			copy_IntMatrixT(x_2, x_1);
			info->omega += 1;
		}
		while (!compare_IntMatrixT(x_1, y_1));
	}
	
	//I need a way to get tau from St and omega
	//If we know omega == 1, then tau is just St - 1 (unless tau = 1)
	
	//We do know that, given an L by L update matrix and a square-free
	// modulus, the maximum transient length is L.
	
	//pg 29 of LCA paper explains why prime powers are more complicated
	//They believe the bound on prime powered systems should be kL, where
	// k is the power of the prime
	
	//Freeing memory
	x_1 = free_IntMatrixT(x_1);
	x_2 = free_IntMatrixT(x_2);
	y_1 = free_IntMatrixT(y_1);
	y_2 = free_IntMatrixT(y_2);
	
	Finv = Finv != NULL ? free_IntMatrixT(Finv) : NULL;
	
	return info;
}


//UNFINISHED
int write_orbits(const char* fileName, IntMatrixTP F, int modulus)
/** Writes out a .orbits file containing every single orbit in the given
    system. Returns 1 on success, 0 otherwise. 
		
		Currently, this function only works with 2 by 2 matrices. */
{
	//F must be a square 2 by 2 matrix
	if ((rows(F) != cols(F)) && (rows(F) != 2))
		return 0;
	
	FILE* orbitsFile = fopen(fileName, "w");
	FILE* linesFile; //Used for holding line numbers for orbits
	
	//Will hold the filename for the line number file
	//The +4 accounts for the new extension and the null character
	char* linesFileName = malloc((strlen(fileName)+4)*sizeof(char));
	
	if (orbitsFile == NULL)
		return 0;
	
	int currentLine = 0; //Where we are in orbitsFile
	int charCount   = 0; //Instead of counting lines, we count chars.
	
	IntMatrixTP s   = new_IntMatrixT(rows(F), 1); //current vector to get orbit
	IntMatrixTP x_1 = new_IntMatrixT(rows(s), 1); //Used for iterating F on s
	IntMatrixTP x_2 = new_IntMatrixT(rows(s), 1);
	
	//Holds info about the current vector's cycle
	CycleInfoTP currentCycle;
	IntMatrixTP repvect;
	
	//This will hold all of our starting values for our vectors
	int* vectors = calloc(rows(F), sizeof(int));
	
	//Holds the vector from the rep's cycle to copy
	int* cycleVect = malloc(2*sizeof(int));
	
	//lineNumbers the line numbers where you can find information about
	// specific vectors in the .orbit file.
	// For example, element lineNumbers[1][1] would hold the line
	// number where you can look up <1, 1>'s info.
	// Of course, this can only really be done when F is 2x2.
	
	//Line numbers start at 0
	int** lineNumbers = malloc(modulus*sizeof(int*));
	for (int i = 0; i < modulus; i += 1)
	{
		lineNumbers[i] = calloc(modulus, sizeof(int));
		for (int j = 0; j < modulus; j += 1)
			lineNumbers[i][j] = -1;
	}
	
	//Holds what character position each vector's orbit
	// starts at.
	int** charNumbers = malloc(modulus*sizeof(int*));
	for (int i = 0; i < modulus; i += 1)
		charNumbers[i] = malloc(modulus*sizeof(int));
	
	bool allTested = FALSE;
	
	
	//Loop until we've found the orbits of all vectors
	while (!allTested)
	{
		set_column(s, vectors);
		
		//Our plan for generating orbits
		// 1. Check to see if the vector's orbit has already been generated
		// 2. Find a rep vector that's in our vector's eventual cycle
		// 3. Generate the rep vector's orbit if needed, update line numbers
		// 4. Generate the original vector's orbit up to the rep vector
		// 5. Put line number array in a different file.
		
		//If this vector's orbit hasn't been generated yet
		if (lineNumbers[element(s, 0, 0)][element(s, 1, 0)] == -1)
		{
			currentCycle = floyd(F, s, modulus);
			repvect = currentCycle->inCycle;
			
			//If we need to generate repvect's orbit
			if (lineNumbers[element(repvect, 0, 0)][element(repvect, 1, 0)] == -1)
			{
				//Iterate until we get back to repvect
				copy_IntMatrixT(repvect, x_1);
				
				//Setting char count so that we can find the orbit in the file later
				charNumbers[element(repvect, 0, 0)][element(repvect, 1, 0)] = charCount;
				lineNumbers[element(repvect, 0, 0)][element(repvect, 1, 0)] = currentLine;
				
				do
				{
					mat_mul(F, x_1, x_2);
					modm(x_2, modulus);
					copy_IntMatrixT(x_2, x_1);
					
					fprintf(orbitsFile, "%d %d\n", 
					element(x_2, 0, 0), element(x_2, 1, 0));
					
					charCount += num_digits(element(x_2, 0, 0)) + 
					             num_digits(element(x_2, 1, 0)) + 2; //+2 for space, \n
					currentLine += 1;
									
					//If we get to a vector we already have the orbit for, we can stop
					if (lineNumbers[element(x_2, 0, 0)][element(x_2, 1, 0)] != -1)
						break;
				}
				while (! compare_IntMatrixT(repvect, x_2));
				
				//Separating different vector orbits
				fprintf(orbitsFile, "-\n");
				charCount += 2;
				currentLine += 1; //+1 accounts for extra space
			}
			
			//At this point, our cycle rep's orbit has been added to the .orbits file
			
			//Making sure that our rep is a different vector
			if (! compare_IntMatrixT(s, repvect))
			{
				lineNumbers[element(s, 0, 0)][element(s, 1, 0)] = currentLine;
				charNumbers[element(s, 0, 0)][element(s, 1, 0)] = charCount;
				
				//Generate orbit for s upto repvect
				copy_IntMatrixT(s, x_1);
				
				while (TRUE)
				{
					mat_mul(F, x_1, x_2);
					modm(x_2, modulus);
					copy_IntMatrixT(x_2, x_1);
					
					fprintf(orbitsFile, "%d %d\n", 
					element(x_2, 0, 0), element(x_2, 1, 0));
					currentLine += 1;
					
					charCount += num_digits(element(x_2, 0, 0)) + 
					             num_digits(element(x_2, 1, 0)) + 2; //+2 for space, \n
							
					//If we get to a vector we already have the orbit for, we can stop
					// It's up to the user to find that previous vector and recreate
					// this vector's orbit faithfully.
					if (lineNumbers[element(x_2, 0, 0)][element(x_2, 1, 0)] != -1)
						break;
				}

				fprintf(orbitsFile, "-\n");
				charCount += 2;
				currentLine += 1;
				
				//Now that we ended the orbit at a vector we already have
				// the orbit for, the user can just jump to that one via
				// the line numbers file
			}
		}
		
		//Increment our starting vector
		allTested = TRUE;
		for (int element = rows(F)-1; element >= 0; element -= 1)
		{
			//This logic counts through all possible vectors mod modulus
			if (vectors[element] == modulus-1)
				vectors[element] = 0;
			
			else
			{
				vectors[element] += 1;
				allTested = FALSE;
				break;
			}
		}
	}
	
	//Now that we've generated all the orbits, we should output
	// a file containing the line numbers for all the vectors' orbits
	strcpy(linesFileName, fileName);
	strcat(linesFileName, "loc");
	linesFile = fopen(linesFileName, "w");
	if (linesFile == NULL)
		return 0;
	
	for (int i = 0; i < modulus; i += 1)
	{
		for (int j = 0; j < modulus; j += 1)
			fprintf(linesFile, "%d ", lineNumbers[i][j]);
		fprintf(linesFile, "\n");
	}
	
	#ifdef VERBOSE
	printf("Line numbers:\n");
	for (int i = 0; i < modulus; i += 1)
	{
		for (int j = 0; j < modulus; j += 1)
			printf("%d ", lineNumbers[i][j]);
		printf("\n");
	}
	#endif //VERBOSE
	
	//Freeing memory
	for (int i = 0; i < modulus; i += 1)
	{
		FREE(lineNumbers[i]);
		FREE(charNumbers[i]);
	}
	FREE(lineNumbers);
	FREE(charNumbers);
	
	FREE(linesFileName);
	
	s   = free_IntMatrixT(s);
	x_1 = free_IntMatrixT(x_1);
	x_2 = free_IntMatrixT(x_2);
	
	//This also frees repvect
	currentCycle = free_CycleInfoT(currentCycle);
	
	FREE(vectors);
	FREE(cycleVect);
	
	if (fclose(orbitsFile) == EOF)
		return 0;
	
	if (fclose(linesFile) == EOF)
		return 0;
	
	return 1;
}
