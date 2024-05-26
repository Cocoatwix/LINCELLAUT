
#include <stdio.h>
#include <stdlib.h>

#include <math.h>

//Basic program to sift through all possible values of m_x, m_y, lambda, and c
// to check which polynomial pairs are possible when creating matrices with
// particular eigenspaces using my patented method

// May 7, 2024

typedef enum boolean {FALSE, TRUE} bool;

//For holding points that represent matrices we've found with
// factorable characteristic polynomials
typedef struct matpoint
{
	int x_1, x_2, y_1, y_2;
}
MatPointT, *MatPointTP;

bool equal_MatPointT(MatPointTP a, MatPointTP b)
/** Checks to see if two MatPointT objects are equal. 
    Returns TRUE if they are, FALSE otherwise. */
{
	if ((a->x_1 != b->x_1) ||
	    (a->y_1 != b->y_1) ||
			(a->x_2 != b->x_2) ||
			(a->y_2 != b->y_2))
		return FALSE;
		
	else
		return TRUE;
}


void print_MatPointT(MatPointTP a)
/** Prints a representation of the MatPoint to stdout. */
{
	printf("(%dx + %dy, %dx + %dy)", a->x_1, a->y_1, a->x_2, a->y_2);
}


int main(int argc, char** argv)
{
	//This list holds all the points representing
	// matrices we've come across so far. It allows us
	// to reject duplicate matrices we've found
	int listOfMatPointsCount        = 0;
	int listOfRegularMatPointsCount = 0;
	MatPointTP listOfMatPoints      = NULL;
	
	//For detecting whether we've seen a particular matrix before
	bool uniqueMatPoint;
	
	int mod = 3;                //The modulus used for the matrices
	int vars[4] = {0, 0, 0, 0}; //m_x, m_y, lambda, c
	int m_x, m_y, lambda, c;
	
	//Use command line argument to choose modulus. If no argument is given, default to 3
	if (argc > 1)
		mod = (int)strtol(argv[1], NULL, 10);
	
	//Holds each point representing a matrix we calculate
	MatPointT tempPoint;

	//First, find pairs of the usual kind (c = 0, 1, 2, ...)
	while (vars[0] < mod)
	{
		m_x    = vars[0];
		m_y    = vars[1];
		lambda = vars[2];
		c      = vars[3];
		
		//Calculate a new matrix
		tempPoint.x_1 = (m_x + lambda) % mod;
		tempPoint.y_1 = (m_x*c) % mod;
		tempPoint.x_2 = (m_y + lambda) % mod;
		tempPoint.y_2 = (m_y*c + lambda*c + lambda) % mod;
		
		//Now, check every matrix point we've calculated so far
		// to see if it's a duplicate
		uniqueMatPoint = TRUE;
		for (int i = 0; i < listOfMatPointsCount; i += 1)
		{
			if (equal_MatPointT(&(listOfMatPoints[i]), &tempPoint))
			{
				uniqueMatPoint = FALSE;
				break;
			}
		}
		
		//If the tempPoint is unique, add it to our list
		if (uniqueMatPoint)
		{
			listOfMatPointsCount += 1;
			listOfMatPoints = realloc(listOfMatPoints, sizeof(MatPointT)*listOfMatPointsCount);
			listOfMatPoints[listOfMatPointsCount - 1] = tempPoint;
		}
		
		//Increment variables
		for (int i = 3; i >= 0; i -= 1)
		{
			vars[i] += 1;
			if ((vars[i] >= mod) && (i != 0))
				vars[i] = 0;
			else
				break;
		}
	}
	
	//Print out what we've found so far
	/*
	for (int i = 0; i < listOfMatPointsCount; i += 1)
	{
		print_MatPointT(&(listOfMatPoints[i]));
		printf("\n");
	}
	printf("---\n");
	*/
	
	//Keeping track of how many "regular" points we found
	listOfRegularMatPointsCount = listOfMatPointsCount;
	
	//Now, find pairs of the odd kind (c = odd)
	vars[0] = 0;
	tempPoint.x_2 = 0; //For all pairs of the odd kind, x_2 = 0
	while (vars[0] < mod)
	{
		m_x    = vars[0];
		m_y    = vars[1];
		lambda = vars[2];
		
		//Calculate our matrix
		tempPoint.x_1 = lambda % mod;
		tempPoint.y_1 = m_x % mod;
		tempPoint.y_2 = m_y % mod;
		
		//Now, check every matrix point we've calculated so far
		// to see if it's a duplicate
		uniqueMatPoint = TRUE;
		for (int i = 0; i < listOfMatPointsCount; i += 1)
		{
			if (equal_MatPointT(&(listOfMatPoints[i]), &tempPoint))
			{
				uniqueMatPoint = FALSE;
				break;
			}
		}
		
		//If the tempPoint is unique, add it to our list
		if (uniqueMatPoint)
		{
			listOfMatPointsCount += 1;
			listOfMatPoints = realloc(listOfMatPoints, sizeof(MatPointT)*listOfMatPointsCount);
			listOfMatPoints[listOfMatPointsCount - 1] = tempPoint;
		}

		//Increment variables
		//There's no c needed for the odd kind, so we only use vars[0] to vars[2]
		for (int i = 2; i >= 0; i -= 1)
		{
			vars[i] += 1;
			if ((vars[i] >= mod) && (i != 0))
				vars[i] = 0;
			else
				break;
		}
	}
	
	//Now, let's print out the remaining "odd" points
	/*
	for (int i = listOfRegularMatPointsCount; i < listOfMatPointsCount; i += 1)
	{
		print_MatPointT(&(listOfMatPoints[i]));
		printf("\n");
	}
	*/
	printf("Number of unique matrices found: %d\n", listOfMatPointsCount);
	printf("Thus, there are %.0lf 2x2 matrices mod %d with irreducible characteristic polynomials.\n",
	pow((double)mod, 4.0) - listOfMatPointsCount, mod);
	
	//Deallocate our list of matrix points
	free(listOfMatPoints);
	listOfMatPoints = NULL;
	
	return EXIT_SUCCESS;
}
