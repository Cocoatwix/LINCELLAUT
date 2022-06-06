
/*
Shared library created for use with ORBITVIS

May 4, 2022
*/

#include <stdlib.h>

#include "../headers/linalg.h"
#include "../headers/cycles.h"
#include "../headers/modular.h"


int** array_to_pointer(int (*arr)[2][2])
/** Helper function to convert a 2D array pointer
    to a double int pointer. It returns the double int
		pointer. */
{
	//This may be one of the worst things I've ever written
	//This is proof of my incompetency regarding Python's ctypes
	int** betterArr = malloc(2*sizeof(int*));
	betterArr[0]    = malloc(2*sizeof(int*));
	betterArr[1]    = malloc(2*sizeof(int*));
	betterArr[0][0] = (*arr)[0][0];
	betterArr[0][1] = (*arr)[0][1];
	betterArr[1][0] = (*arr)[1][0];
	betterArr[1][1] = (*arr)[1][1];
	
	return betterArr;
}


void free_double_pointer(int*** arr)
/** Helper function to free memory of a double pointer. */
{
	free((*arr)[0]);
	free((*arr)[1]);
	(*arr)[0] = NULL;
	(*arr)[1] = NULL;
	free(*arr);
	*arr = NULL;
}


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
	
	int** betterArr = array_to_pointer(arr);
	
	set_matrix(F, betterArr);
	int vect[2] = {x, y};
	set_column(s_0, vect);
	
	s_f = iterate(F, s_0, modulus, t);
	a = element(s_f, 0, 0);
	b = element(s_f, 1, 0);
	
	free_double_pointer(&betterArr);
	
	F   = free_IntMatrixT(F);
	s_0 = free_IntMatrixT(s_0);
	s_f = free_IntMatrixT(s_f);
	
	return a*modulus + b;
}


int get_orbit_info(int (*v)[2], int (*arr)[2][2], int modulus)
/** This function gets the cycle length and transient length of
    the given vector under the given system and returns it as a
		single integer: tau + modulus*omega */
{
	//Initialise update matrix and initial vector, then iterate.
	IntMatrixTP F   = new_IntMatrixT(2, 2);
	IntMatrixTP s_0 = new_IntMatrixT(2, 1);
	CycleInfoTP c;
	
	int** betterArr = array_to_pointer(arr);
	int   omegaNum, tauNum;
	
	set_matrix(F, betterArr);
	set_column(s_0, *v);
	c = floyd(F, s_0, modulus);
	
	omegaNum = omega(c);
	tauNum = tau(c);
	
	free_double_pointer(&betterArr);
	F   = free_IntMatrixT(F);
	s_0 = free_IntMatrixT(s_0);
	c   = free_CycleInfoT(c);
	
	//In the future, this should also return tau in some way
	return omegaNum;
}


int get_orbit_info_array(int (*update)[2][2], int modulus)
/** This function gets the cycle length and transient length of
    the given matrix and returns it as a
		single integer: tau + modulus*omega */
{
	//Initialise update matrix then iterate.
	IntMatrixTP F = new_IntMatrixT(2, 2);
	CycleInfoTP c;
	
	int** betterArr = array_to_pointer(update);
	int   omegaNum, tauNum;
	int   k = 2; //Needed for ensuring our encoding method for multiple numbers doesn't break
	
	set_matrix(F, betterArr);
	c = floyd(F, F, modulus);
	
	omegaNum = omega(c);
	tauNum   = tau(c);
	
	free_double_pointer(&betterArr);
	F = free_IntMatrixT(F);
	c = free_CycleInfoT(c);
	
	//In the future, this should also return tau in some way
	return tauNum + modulus*k*omegaNum;
}
