
/*
Shared library created for use with ORBITVIS

May 4, 2022
*/

#include <stdlib.h>

#include "../headers/linalg.h"
#include "../headers/cycles.h"


//This function returns a single integer instead
// of a vector to make Python usage easier
int C_step(int x, int y, int (*arr)[2][2], int modulus, int t)
/** This function iterates the vector <x, y> under an update matrix F
    mod modulus t times. Let <a, b> be the resulting 
		vector. Then, this function returns modulus*a + b. */
{
	//Initialise update matrix and initial vector, then iterate.
	IntMatrixTP F   = new_IntMatrixT(2, 2);
	IntMatrixTP s_0 = new_IntMatrixT(2, 1);
	IntMatrixTP s_f;
	
	int a, b;
	
	//This may be one of the worst things I've ever written
	//This is proof of my incompetency regarding Python's ctypes
	int** betterArr = malloc(2*sizeof(int*));
	betterArr[0] = malloc(2*sizeof(int*));
	betterArr[1] = malloc(2*sizeof(int*));
	betterArr[0][0] = (*arr)[0][0];
	betterArr[0][1] = (*arr)[0][1];
	betterArr[1][0] = (*arr)[1][0];
	betterArr[1][1] = (*arr)[1][1];
	
	set_matrix(F, betterArr);
	int vect[2] = {x, y};
	set_column(s_0, vect);
	
	s_f = iterate(F, s_0, modulus, t);
	a = element(s_f, 0, 0);
	b = element(s_f, 1, 0);
	
	free(betterArr[0]);
	free(betterArr[1]);
	betterArr[0] = NULL;
	betterArr[1] = NULL;
	free(betterArr);
	betterArr = NULL;
	
	F   = free_IntMatrixT(F);
	s_0 = free_IntMatrixT(s_0);
	s_f = free_IntMatrixT(s_f);
	
	return a*modulus + b;
}