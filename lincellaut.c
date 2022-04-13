
/* A simple program for calculating various things about linear
 *  cellular automaton. Coded for my research position with
 *  Prof Mendivil.
 *
 * Mar 27, 2022
 *
 */
 
#include <stdlib.h>
#include <stdio.h>

#include "headers/linalg.h" //Homemade linear algebra library
#include "headers/cycles.h" //Allows us to use Floyd's Algorithm

/* Floyd's Cycle Detection Algorithm
 *
 * - Set x_0 = y_0
 * - x_(n+1) = F(x_n)
 * - y_(n+1) = F(F(y_n))
 * - Repeat the above process until x_n == y_n
 *  - x_n will then be an element in a cycle.
 *
 * The algorithm will have a stopping time of:
 *  tau (transient length) if tau == 0 mod omega (cycle length)
 *  tau + (omega - (tau mod omega)) otherwise.
 
 
 tau + omega - (tau mod omega) ~~ -omega
 tau - (tau mod omega)
 n*omega + N - N
 n*omega ~~ /omega
 n

 however
 
 n*omega ~~ -omega
 (n-1)*omega ~~ /omega
 (n-1)
 
 maybe we run floyd's cycle algorithm again with different step sizes (+2 and +3?),
 figure out what the stopping time for the new setup is,
 and creating a system of equations to solve for omega and tau?
 
 we could also check if the given matrix is invertible to see if transient regions
 even exist...
 
 *
 * If tau < omega, then stopping time = omega
 * Every option for the stopping time is greater than or equal to omega,
 *  so the stopping time gives an upper bound on omega (omega can be no greater than the stopping time).
 *
 * We can also say that if St > 2tau, then St = omega.
 * 2tau < tau + (omega - (tau mod omega)) and rearrange inequality.
 *  You'll find that tau must be less than omega in this case.
 */
 
/* Once we have x_n in the cycle, if we can't conclude omega,
 *  then we can set y = x_n and keep iterating y until it equals
 *  x_n again. The number of steps this takes will give us omega.
 */
 
/* We also know that if our update rule F doesn't have an inverse,
 *  then our automaton has transient regions. This is because for finite sets
 *  that map to themselves, being one-to-one is the same as onto (think about it!),
 *  and if F isn't one-to-one, then there's some configuration that doesn't get mapped to,
 *  hence it would be a transient configuration.
 *
 */


int main()
{
	//FIle paths to where matrices are stored
	char* const UPDATEFILEPATH = "matrices/update.matrix";
	char* const INITIALFILEPATH = "matrices/initial.matrix";
	
	//const int ITERATIONS = 12345;
	const int MODULUS = 25;
	
	IntMatrixTP   F; //Update rule matrix
	
	IntMatrixTP s_0; //Stores our initial vector
	//IntMatrixTP s_f; //Stores our final vector
	
	//Initalising our matrices/vectors
	read_IntMatrixT(UPDATEFILEPATH, &F);
	read_IntMatrixT(INITIALFILEPATH, &s_0);
	
	/* Iterate s_0 a few times
	s_f = iterate(F, s_0, MODULUS, ITERATIONS);
	printm(s_f); */
	
	printcycle(floyd(F, s_0, MODULUS));
	
	//Freeing memory
	F   = free_IntMatrixT(F);
	s_0 = free_IntMatrixT(s_0);
	//s_f = free_IntMatrixT(s_f);
	
	return EXIT_SUCCESS;
}

/*
0 1 2 3 4 5
0 2 4 0 2 4
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
even transient length w + 2n
even cycle length w

0 1 2 3 4 5 6 7 | x
0 2 4 6 0 2 4 6 | y

y's row will be shifted to the left n times
0 1 2 3 4 5 6 7 | x
2 4 6 0 2 4 6 0 | y

number of steps = t + (w - 2n mod w)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
odd transient length w + 2n + 1
even cycle length w

0 1 2 3 4 5 6 7 | x
1 3 5 7 1 3 5 7 | y

y's row will be shifted to the left n times
number of steps = t + (w - 2n+1 mod w) //The +1 accounts for the different starting position
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
even transient length w + 2n + 1
odd cycle length w

0 1 2 3 4 5 6 7 8 | x
0 2 4 6 8 1 3 5 7 | y

With odd cycle lengths, the starting number for y will always
be 2n + 1 mod w:
0 1 2 3 4 5 6 7 8 | x
1 3 5 7 0 2 4 6 8 | y
number of steps = t + (w - 2n + 1 mod w)

odd transient length w + 2n
odd cycle length w

0 1 2 3 4 5 6 7 8 | x
0 2 4 6 8 1 3 5 7 | y
Starting number will always be 2n mod w:
0 1 2 3 4 5 6 7 8 | x
2 4 6 8 1 3 5 7 0 | y
number of steps = t + (w - 2n mod w)

This proves the claim in the paper.
Just mess around with numbers and see that this works.

*/