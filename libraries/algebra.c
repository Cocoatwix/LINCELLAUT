
/*
A simple library for manipulating
more typical algebra things, like
polynomials.

yiff day, 2022
*/

#include <stdlib.h>
#include <stdio.h>

#include <complex.h>
#include <string.h>

#include "../headers/helper.h" //bool 

#include "../headers/bigint.h"
#include "../headers/linalg.h"
#include "../headers/modular.h"

const double PI = 3.1415926535897932384626;

typedef struct bigpoly
/** Stores a polynomial using BigIntT coefficients. 
    Holds coefficients in little endian form. */
{
	int size; //Holds how many terms are in the polynomial
	BigIntTP* coeffs; //Holds the coefficients of the polynomial
	
	char* variable; //The symbol used for displaying polynomials
}
BigPolyT, *BigPolyTP;


const int MAXVARLEN = 20;
typedef struct fieldexp
/** Holds an algebraic expression containing field extensions. */
{
	//Constant/polynomial expressions representing the current state of the expression
	BigPolyTP* expressions; 
	
	//List of elements contained in our expression
	//Regular elements are represented by constant polynomials
	//Extensions are represented by their minimal polynomial
	BigPolyTP* elements;
	int numOfElements;
	int numOfElementsUsed;
	
	//Names to give field extensions
	char**    extNames;
}
FieldExpT, *FieldExpTP;


BigPolyTP free_BigPolyT(BigPolyTP p)
/** Frees the memory of a given BigPolyT. Returns NULL. */
{
	if (p != NULL)
	{
		for (int i = 0; i < p->size; i += 1)
			p->coeffs[i] = free_BigIntT(p->coeffs[i]);
		
		free(p->coeffs);
		p->coeffs = NULL;
		
		free(p->variable);
		p->variable = NULL;
	}
	
	return NULL;
}


FieldExpTP free_FieldExpT(FieldExpTP e)
/** Frees the memory of a given FieldExpT. Returns NULL. */
{
	if (e == NULL)
		return NULL;
	
	for (int i = 0; i < e->numOfElements; i += 1)
	{
		e->expressions[i] = free_BigPolyT(e->expressions[i]);
		e->elements[i]    = free_BigPolyT(e->elements[i]);
		free(e->extNames[i]);
		e->extNames[i] = NULL;
	}
	free(e->expressions);
	e->expressions = NULL;
	free(e->elements);
	e->elements = NULL;
	free(e->extNames);
	e->extNames = NULL;
	
	return NULL;
}


BigPolyTP* free_BigPolyT_factors(BigPolyTP* factors)
/** Frees an array of BigPolyT factors. */
{
	if (factors != NULL)
	{
		BigIntTP one;
		BigIntTP temp;
		BigIntTP indexCounter;
		
		int oneArr[1] = {1};
		int index = 1;
		
		one = new_BigIntT(oneArr, 1);
		indexCounter = new_BigIntT(oneArr, 1);
		temp = empty_BigIntT(1);
		
		while (compare_BigIntT(indexCounter, factors[0]->coeffs[0]) <= 0)
		{
			factors[index] = free_BigPolyT(factors[index]);
			
			index += 1;
			add_BigIntT(one, indexCounter, temp);
			copy_BigIntT(temp, indexCounter);
		}
		factors[0] = free_BigPolyT(factors[0]);
		free(factors);
		
		one  = free_BigIntT(one);
		temp = free_BigIntT(temp);
		indexCounter = free_BigIntT(indexCounter);
	}
	
	return NULL;
}


BigPolyTP new_BigPolyT(BigIntTP* const coefficients, int size)
/** Creates a new BigPolyT with specified coefficients and size,
    returns a pointer to it. */
{
	BigPolyTP p = malloc(sizeof(BigPolyT));
	p->size = size;
	
	p->coeffs = malloc(size*sizeof(BigIntTP));
	for (int i = 0; i < size; i += 1)
	{
		p->coeffs[i] = empty_BigIntT(1);
		copy_BigIntT(coefficients[i], p->coeffs[i]);
	}
	
	p->variable = malloc((MAXVARLEN+1)*sizeof(char));
	p->variable[0] = '\0';
	strcat(p->variable, "λ");
	
	return p;
}


BigPolyTP constant_BigPolyT(BigIntTP const constant)
/** Creates a constant BigPolyTP, returns it. */
{
	BigPolyTP p = malloc(sizeof(BigPolyT));
	p->size = 1;
	p->coeffs = malloc(sizeof(BigIntTP));
	p->coeffs[0] = empty_BigIntT(1);
	copy_BigIntT(constant, p->coeffs[0]);
	
	p->variable = malloc((MAXVARLEN+1)*sizeof(char));
	p->variable[0] = '\0';
	strcat(p->variable, "λ");
	
	return p;
}


BigPolyTP empty_BigPolyT()
/** Creates an empty BigPolyT struct and returns a pointer
    to it. */
{
	BigPolyTP p = malloc(sizeof(BigPolyT));
	p->size = 1;
	
	p->coeffs = malloc(sizeof(BigIntTP));
	p->coeffs[0] = empty_BigIntT(1);
	
	p->variable = malloc((MAXVARLEN+1)*sizeof(char));
	p->variable[0] = '\0';
	strcat(p->variable, "λ");
	
	return p;
}


int set_BigPolyT_var(BigPolyTP p, char* const name)
/** Set the name of the variable used for a polynomial. 
    Returns 1 on success, 0 otherwise. */
{
	if (strlen(name) > (unsigned int)MAXVARLEN)
	{
		fprintf(stderr, "Given variable name is too large.\n");
		return 0;
	}
	
	//Clear variable string before setting new name
	p->variable[0] = '\0';
	strcat(p->variable, name);
	return 1;
}


FieldExpTP new_FieldExpT(int size)
/** Creates a new FieldExpT object, returns a pointer to it. */
{
	FieldExpTP exp = malloc(sizeof(FieldExpT));
	
	exp->expressions = malloc(size*sizeof(BigPolyTP));
	exp->elements    = malloc(size*sizeof(BigPolyTP));
	exp->extNames    = malloc(size*sizeof(char*));
	
	for (int i = 0; i < size; i += 1)
	{
		exp->expressions[i] = empty_BigPolyT();
		exp->elements[i]    = empty_BigPolyT();
		exp->extNames[i]    = malloc((MAXVARLEN+1)*sizeof(char));
		exp->extNames[i][0] = '\0';
	}
	
	exp->numOfElements = size;
	exp->numOfElementsUsed = 0;
	return exp;
}


int add_extension(FieldExpTP extexp, BigPolyTP const minPoly, BigPolyTP const value, char* const name)
/** Adds a new extension into our expression with given minPoly, value, and
    name. Returns 1 on success, 0 otherwise. */
{
	int index = extexp->numOfElementsUsed;
	
	//Resize FieldExpT if we have to
	if (index == extexp->numOfElements)
	{
		extexp->numOfElements += 1;
		extexp->expressions = realloc(extexp->expressions, (extexp->numOfElements)*sizeof(BigPolyTP));
		extexp->elements    = realloc(extexp->elements, (extexp->numOfElements)*sizeof(BigPolyTP));
		extexp->extNames    = realloc(extexp->extNames, (extexp->numOfElements)*sizeof(char*));
		
		extexp->expressions[extexp->numOfElements - 1] = empty_BigPolyT();
		extexp->elements[extexp->numOfElements - 1]    = empty_BigPolyT();
		extexp->extNames[extexp->numOfElements - 1]    = malloc((MAXVARLEN+1)*sizeof(char));
		extexp->extNames[extexp->numOfElements - 1][0] = '\0';
	}
	
	extexp->numOfElementsUsed += 1;
	
	copy_BigPolyT(value, extexp->expressions[index]);
	copy_BigPolyT(minPoly, extexp->elements[index]);
	
	if (strlen(name) <= (unsigned int)MAXVARLEN) //I think this cast is okay...
	{
		strcat(extexp->extNames[index], name);
		set_BigPolyT_var(extexp->expressions[index], name);
	}
	else
	{
		fprintf(stderr, "Given extension name is too long.\n");
		strcat(extexp->extNames[index], "?");
		set_BigPolyT_var(extexp->expressions[index], "?");
	}
	
	return 1;
}


int collapse_field_extension(BigPolyTP const minPoly, BigIntTP const mod)
/** For a field extension a, find the minimum exponent c such that
    a^c is in the original field. In this case, the field is the
    integers modulo the BigIntT given. The BigPolyT given specifies
    the minimal polynomial for the field extension. 
    Returns the exponent c on success, -1 otherwise. */
{
	int exponent = 1;
	BigIntTP* minPolyCoeffs;
	BigIntTP* rewrittenCoeffs;
	BigIntTP* tempCoeffs;
	BigIntTP* tempCoeffs2;
	
	int oneArr[1] = {1};
	BigIntTP one;
	BigIntTP zero;
	BigIntTP negOne;
	BigIntTP temp;
	
	BigIntTP leadingTermInv;
	
	bool collapsed = FALSE;
	
	//This function assumes minPoly has been reduced.
		
	one            = new_BigIntT(oneArr, 1);
	zero           = empty_BigIntT(1);
	negOne         = empty_BigIntT(1);
	temp           = empty_BigIntT(1);
	leadingTermInv = empty_BigIntT(1);
		
	minPolyCoeffs   = extract_coefficients(minPoly);
	rewrittenCoeffs = malloc((minPoly->size - 1)*sizeof(BigIntTP));
	tempCoeffs      = malloc((minPoly->size - 1)*sizeof(BigIntTP));
	tempCoeffs2     = malloc((minPoly->size - 1)*sizeof(BigIntTP));
	big_num_inverse(minPolyCoeffs[minPoly->size - 1], mod, leadingTermInv);
	subtract_BigIntT(mod, one, negOne);
	
	//Now, go through our coefficients and negate them, then divide by the leading term
	for (int i = 0; i < minPoly->size - 1; i += 1)
	{
		//Move to other side of equal sign
		rewrittenCoeffs[i] = empty_BigIntT(1);
		multiply_BigIntT(minPolyCoeffs[i], negOne, temp);
		mod_BigIntT(temp, mod, rewrittenCoeffs[i]);
		
		//Multiply by inverse of leading term
		multiply_BigIntT(rewrittenCoeffs[i], leadingTermInv, temp);
		mod_BigIntT(temp, mod, rewrittenCoeffs[i]);
	}
	
	//Now our relation is rewritten in a form we can use.
	//Time to iterate repeatedly until we get a constant
	for (int i = 0; i < minPoly->size - 1; i += 1)
	{
		tempCoeffs[i]  = empty_BigIntT(1);
		tempCoeffs2[i] = empty_BigIntT(1);
		copy_BigIntT(rewrittenCoeffs[i], tempCoeffs[i]);
	}
	
	//Check to see if we have a constant
	collapsed = TRUE;
	for (int i = 1; i < minPoly->size - 1; i += 1)
	{
		if (compare_BigIntT(tempCoeffs[i], zero) != 0)
		{
			collapsed = FALSE;
			break;
		}
	}
	
	//Iterate our extension until it collapses to an element in the field
	while (!collapsed)
	{
		exponent += 1;
		
		//Iterate through all coefficients, update them
		for (int i = 0; i < minPoly->size - 1; i += 1)
		{
			multiply_BigIntT(rewrittenCoeffs[i], tempCoeffs[minPoly->size - 2], temp);
			if (i > 0) //If we're not updating the constant term
			{
				add_BigIntT(temp, tempCoeffs[i-1], tempCoeffs2[i]);
				copy_BigIntT(tempCoeffs2[i], temp);
			}
			mod_BigIntT(temp, mod, tempCoeffs2[i]);
		}
		
		//tempCoeffs2 should now hold our updated elements. Copy them to tempCoeffs
		for (int i = 0; i < minPoly->size - 1; i += 1)
		{
			copy_BigIntT(tempCoeffs2[i], tempCoeffs[i]);
			//printi(tempCoeffs[i]);
			//printf(", ");
		}
		//printf("\n");
		
		//Check to see if we have a constant
		collapsed = TRUE;
		for (int i = 1; i < minPoly->size - 1; i += 1)
		{
			if (compare_BigIntT(tempCoeffs[i], zero) != 0)
			{
				collapsed = FALSE;
				break;
			}
		}
	}
	
	
	for (int i = 0; i < minPoly->size; i += 1)
	{
		if (i < minPoly->size - 1)
		{
			rewrittenCoeffs[i] = free_BigIntT(rewrittenCoeffs[i]);
			tempCoeffs[i]      = free_BigIntT(tempCoeffs[i]);
			tempCoeffs2[i]     = free_BigIntT(tempCoeffs2[i]);
		}
		
		minPolyCoeffs[i] = free_BigIntT(minPolyCoeffs[i]);
	}
	free(minPolyCoeffs);
	free(rewrittenCoeffs);
	free(tempCoeffs);
	free(tempCoeffs2);
	minPolyCoeffs   = NULL;
	rewrittenCoeffs = NULL;
	tempCoeffs      = NULL;
	tempCoeffs2     = NULL;
	
	one    = free_BigIntT(one);
	zero   = free_BigIntT(zero);
	temp   = free_BigIntT(temp);
	negOne = free_BigIntT(negOne);
	
	return exponent;
}


int reduce_FieldExpT(FieldExpTP exp, BigIntTP const mod)
/** Reduces a FieldExpT by rewriting the expressions using the
    extension relations. Returns 1 on success, 0 otherwise. */
{
	BigIntTP* coeffList;       //Holds extracted coeffs from exp
	BigIntTP* rewrittenCoeffs; //Holds the relation expression
	BigIntTP* tempCoeffs;      //Holds the relation expression for higher powers
	BigIntTP* tempCoeffs2;     //Used for intermediate computations
	
	int oneArr[1] = {1};
	BigIntTP one;
	BigIntTP temp;
	BigIntTP zero;
	BigIntTP negOne;
	BigIntTP leadingTermInv;
	
	int sizeOfElem;
	
	one            = new_BigIntT(oneArr, 1);
	zero           = empty_BigIntT(1);
	temp           = empty_BigIntT(1);
	negOne         = empty_BigIntT(1);
	leadingTermInv = empty_BigIntT(1);
	subtract_BigIntT(mod, one, negOne);
	
	//Iterate over each expression we have
	for (int elem = 0; elem < exp->numOfElementsUsed; elem += 1)
	{
		sizeOfElem = exp->elements[elem]->size;
		
		//Reduce each expression using minPolys
		reduce_BigPolyT(exp->expressions[elem]);
		if (exp->expressions[elem]->size >= sizeOfElem) //If we actually have stuff to reduce
		{
			coeffList = extract_coefficients(exp->elements[elem]);
			rewrittenCoeffs = malloc((sizeOfElem-1)*sizeof(BigIntTP));
			big_num_inverse(coeffList[sizeOfElem - 1], mod, leadingTermInv);
			
			//Move all but highest term to other side of =, divide by leading term
			for (int i = 0; i < sizeOfElem - 1; i += 1)
			{
				rewrittenCoeffs[i] = empty_BigIntT(1);
				multiply_BigIntT(coeffList[i], negOne, temp);
				multiply_BigIntT(temp, leadingTermInv, rewrittenCoeffs[i]);
				mod_BigIntT(rewrittenCoeffs[i], mod, temp);
				copy_BigIntT(temp, rewrittenCoeffs[i]);
				
				//Print out elements to see what I'm doing
				printi(rewrittenCoeffs[i]);
				printf(", ");
			}
			printf("\n");
			
			//Free memory for next extension
			for (int i = 0; i < sizeOfElem; i += 1)
			{
				if (i != sizeOfElem-1)
					rewrittenCoeffs[i] = free_BigIntT(rewrittenCoeffs[i]);
				
				coeffList[i] = free_BigIntT(coeffList[i]);
			}
			free(coeffList);
			free(rewrittenCoeffs);
			coeffList       = NULL;
			rewrittenCoeffs = NULL;
			
			copy_BigIntT(zero, leadingTermInv);
		}
		
		//Take modulus of each expression
	}
	
	one            = free_BigIntT(one);
	zero           = free_BigIntT(zero);
	temp           = free_BigIntT(zero);
	negOne         = free_BigIntT(negOne);
	leadingTermInv = free_BigIntT(leadingTermInv);
	
	return 1;
}


int reduce_BigPolyT(BigPolyTP p)
/** Ensures that the leading term of the BigPolyT
    is nonzero. 
		Returns 1 on success, 0 otherwise. */
{
	BigIntTP zero = empty_BigIntT(1);
	int newSize = 0;
	
	for (int i = p->size-1; i >= 0; i -= 1)
	{
		//If we find the leading nonzero term
		if (compare_BigIntT(p->coeffs[i], zero) != 0)
		{
			newSize = i+1;
			break;
		}
	}
	
	if (newSize != p->size)
	{
		//Freeing the coefficients we don't need anymore
		for (int i = p->size-1; i >= newSize; i -= 1)
			p->coeffs[i] = free_BigIntT(p->coeffs[i]);
		
		p->size = newSize;
		p->coeffs = realloc(p->coeffs, (p->size)*sizeof(BigIntTP));
	}
	
	zero = free_BigIntT(zero);
	return 1;
}


int degree(BigPolyTP const p)
/** Returns the degree of a BigPolyTP. */
{
	return p->size - 1;
}


BigIntTP constant(BigPolyTP const p)
/** Returns the constant term of a BigPolyTP. */
{
	return p->coeffs[0];
}


BigIntTP* extract_coefficients(BigPolyTP const p)
/** Returns a list of the BigPolyT's coefficients, in
    order of ascending term exponents.
		Returns NULL on error. */
{
	BigIntTP* c = malloc((p->size)*sizeof(BigIntTP));
	for (int i = 0; i < p->size; i += 1)
	{
		c[i] = empty_BigIntT(1);
		copy_BigIntT(p->coeffs[i], c[i]);
	}
	
	return c;
}


int copy_BigPolyT(BigPolyTP const toCopy, BigPolyTP copyTo)
/** Copies a polynomial and stores it in another BigPolyTP.
    Returns 1 on success, 0 otherwise. */
{
	//Deallocate unneeded terms
	if (copyTo->size > toCopy->size)
		for (int i = toCopy->size; i < copyTo->size; i += 1)
			copyTo->coeffs[i] = free_BigIntT(copyTo->coeffs[i]);
		
	if (copyTo->size != toCopy->size)
		copyTo->coeffs = realloc(copyTo->coeffs, (toCopy->size)*sizeof(BigIntTP));
	
	//Initialise terms as needed
	if (copyTo->size < toCopy->size)
		for (int i = copyTo->size; i < toCopy->size; i += 1)
			copyTo->coeffs[i] = empty_BigIntT(1);
		
	copyTo->size = toCopy->size;
	
	//Now, we can safely copy the terms into the polynomial
	for (int i = 0; i < copyTo->size; i += 1)
		copy_BigIntT(toCopy->coeffs[i], copyTo->coeffs[i]);
	
	//Copy variable name used
	strcpy(toCopy->variable, copyTo->variable);
	
	//Just to be safe
	reduce_BigPolyT(copyTo);
	
	return 1;
}


int compare_BigPolyT(BigPolyTP const A, BigPolyTP const B)
/** Compares two polynomials to see if they're the same.
    Returns 0 if they are, 1 otherwise. */
{
	int small = A->size < B->size ? A->size : B->size;
	int big   = A->size > B->size ? A->size : B->size;
	int r = 0;
	
	BigPolyTP theBig = degree(A) > degree(B) ? A : B;
	BigIntTP zero = empty_BigIntT(1);
	
	//Comparing terms to see if they match
	for (int i = 0; i < small; i += 1)
		if (compare_BigIntT(A->coeffs[i], B->coeffs[i]) != 0)
		{
			r = 1;
			break;
		}
		
	//For the bigger polynomial, we ensure all the entries above the smallest
	// on the other polynomial are zero
	if (r == 0)
		for (int i = small; i < big; i += 1)
			if (compare_BigIntT(theBig->coeffs[i], zero) != 0)
			{
				r = 1;
				break;
			}
	
	zero = free_BigIntT(zero);
	return r;
}


void printp(BigPolyTP const p)
/** Outputs a BigPolyTP to stdout. */
{
	BigIntTP zero = empty_BigIntT(1);
	bool printPlus = FALSE;
	
	for (int i = 0; i < p->size; i += 1)
	{
		//Logic surrounding when to print a +
		if (i != 0)
			if (compare_BigIntT(zero, p->coeffs[i-1]) != 0) //If last number wasn't zero
				printPlus = TRUE;
					
		//Only print if coefficient isn't zero
		if (compare_BigIntT(zero, p->coeffs[i]) != 0)
		{
			//Prevents a + when all other terms are zero afterwards
			if (printPlus)
			{
				printf(" + ");
				printPlus = FALSE;
			}
			
			if (i == 0)
				printi(p->coeffs[0]);
			
			else if (i == 1)
			{
				printi(p->coeffs[1]);
				printf("(%s)", p->variable);
			}
			
			else
			{
				printi(p->coeffs[i]);
				printf("(%s)^%d", p->variable, i);
			}
		}
	}
	
	zero = free_BigIntT(zero);
}


/* private */ void printp_term(BigPolyTP const p, int term)
/** Prints a specific term in the polynomial. Used for nicely printing
    FieldExpTs. */
{
	if (term == 0)
		printi(p->coeffs[0]);
	
	else if (term == 1)
	{
		printi(p->coeffs[1]);
		printf("(%s)", p->variable);
	}
	
	else
	{
		printi(p->coeffs[term]);
		printf("(%s)^%d", p->variable, term);
	}
}


void printpf(BigPolyTP* factors)
/** Prints a factorised BigPolyT to stdout. */
{
	BigIntTP indexCounter;
	BigIntTP one;
	BigIntTP temp;
	int oneArr[1] = {1};
	
	one = new_BigIntT(oneArr, 1);
	indexCounter = new_BigIntT(oneArr, 1);
	temp = empty_BigIntT(1);
	int index = 1;
	
	while (compare_BigIntT(indexCounter, factors[0]->coeffs[0]) <= 0)
	{
		printf("(");
		printp(factors[index]);
		printf(")");
		
		index += 1;
		add_BigIntT(indexCounter, one, temp);
		copy_BigIntT(temp, indexCounter);
	}
	
	temp = free_BigIntT(temp);
	one  = free_BigIntT(one);
	indexCounter = free_BigIntT(indexCounter);
}


void printfe(FieldExpTP const fe)
/** Print a FieldExpTP to stdout. */
{
	BigIntTP zero = empty_BigIntT(1);
	
	//Find the maximum number of elements 
	int maxNumOfTerms = 0;
	for (int elem = 0; elem < fe->numOfElementsUsed; elem += 1)
		if (fe->expressions[elem]->size > maxNumOfTerms)
			maxNumOfTerms = fe->expressions[elem]->size;
	
	for (int term = 0; term < maxNumOfTerms; term += 1) //Keep printing terms for each element until we run out
	{
		for (int elem = 0; elem < fe->numOfElementsUsed; elem += 1)
		{
			/*
			printf("term: %d, elem: %d\n", term, elem);
			printf("First bool: %d\n", fe->expressions[elem]->size < term);
			printf("Second bool: %d\n", compare_BigIntT(zero, fe->expressions[elem]->coeffs[term]) != 0);
			*/
			if ((fe->expressions[elem]->size > term) && 
			    (compare_BigIntT(zero, fe->expressions[elem]->coeffs[term]) != 0))
			{
				if ((term != 0) || (elem != 0))
					printf(" + ");
			
				printp_term(fe->expressions[elem], term);
			}
		}
	}
	
	zero = free_BigIntT(zero);
}


void fprintp(FILE* file, BigPolyTP const p)
/** Outputs a BigPolyTP to a given file stream.. */
{
	BigIntTP zero = empty_BigIntT(1);
	bool printPlus = FALSE;
	
	for (int i = 0; i < p->size; i += 1)
	{
		//Logic surrounding when to print a +
		if (i != 0)
			if (compare_BigIntT(zero, p->coeffs[i-1]) != 0) //If last number wasn't zero
				printPlus = TRUE;
					
		//Only print if coefficient isn't zero
		if (compare_BigIntT(zero, p->coeffs[i]) != 0)
		{
			//Prevents a + when all other terms are zero afterwards
			if (printPlus)
			{
				fprintf(file, " + ");
				printPlus = FALSE;
			}
			
			if (i == 0)
				fprinti(file, p->coeffs[0]);
			
			else if (i == 1)
			{
				fprinti(file, p->coeffs[1]);
				fprintf(file, "(%s)", p->variable);
			}
			
			else
			{
				fprinti(file, p->coeffs[i]);
				fprintf(file, "(%s)^%d", p->variable, i);
			}
		}
	}
	
	zero = free_BigIntT(zero);
}


void fprintpf(FILE* file, BigPolyTP* factors)
/** Prints a factorised BigPolyT to stdout. */
{
	BigIntTP indexCounter;
	BigIntTP one;
	BigIntTP temp;
	int oneArr[1] = {1};
	
	one = new_BigIntT(oneArr, 1);
	indexCounter = new_BigIntT(oneArr, 1);
	temp = empty_BigIntT(1);
	int index = 1;
	
	while (compare_BigIntT(indexCounter, factors[0]->coeffs[0]) <= 0)
	{
		fprintf(file, "(");
		fprintp(file, factors[index]);
		fprintf(file, ")");
		
		index += 1;
		add_BigIntT(indexCounter, one, temp);
		copy_BigIntT(temp, indexCounter);
	}
	
	temp = free_BigIntT(temp);
	one  = free_BigIntT(one);
	indexCounter = free_BigIntT(indexCounter);
}


int add_BigPolyT(BigPolyTP const A, BigPolyTP const B, BigPolyTP sum)
/** Adds A and B together, stores the sum in sum.
    This function assumes sum has been initialised.
		Returns 1 on success, 0 otherwise. */
{
	int bigSize = (A->size > B->size) ? A->size : B->size;
	int smallSize = (A->size < B->size) ? A->size : B->size;
	
	BigIntTP temp = empty_BigIntT(1);
	
	//Making sure our sum is ready to take new terms
	if (sum->size > bigSize)
		for (int i = bigSize; i < sum->size; i += 1)
			sum->coeffs[i] = free_BigIntT(sum->coeffs[i]);
		
	if (sum->size != bigSize)
		sum->coeffs = realloc(sum->coeffs, bigSize*sizeof(BigIntTP));
	
	if (sum->size < bigSize)
		for (int i = sum->size; i < bigSize; i += 1)
			sum->coeffs[i] = empty_BigIntT(1);
			
	sum->size = bigSize;
	
	//Iterate over both polynomials, add them
	for (int a = 0; a < smallSize; a += 1)
	{
		add_BigIntT(A->coeffs[a], B->coeffs[a], temp);
		copy_BigIntT(temp, sum->coeffs[a]);
	}
	
	//Copy remaining terms 
	if (A->size > B->size)
		for (int b = smallSize; b < bigSize; b += 1)
			copy_BigIntT(A->coeffs[b], sum->coeffs[b]);
	
	else
		for (int b = smallSize; b < bigSize; b += 1)
			copy_BigIntT(B->coeffs[b], sum->coeffs[b]);
	
	temp = free_BigIntT(temp);
	
	reduce_BigPolyT(sum);
	return 1;
}


int multiply_BigPolyT(BigPolyTP const A, BigPolyTP const B, BigPolyTP product)
/** Multiples two polynomials together, stores the product in product.
    This function assumes product has been initialised. 
		Returns 1 on success, 0 otherwise. */
{
	BigIntTP temp = empty_BigIntT(1);
	BigIntTP temp2 = empty_BigIntT(1);
	
	//Preparing product to hold the new polynomial
	if (product->size != A->size + B->size - 1)
	{
		//Super inefficient, but prevents errors
		//Will change later
		for (int i = 0; i < product->size; i += 1)
			product->coeffs[i] = free_BigIntT(product->coeffs[i]);
		
		product->size = A->size + B->size - 1;
		product->coeffs = realloc(product->coeffs, (product->size)*sizeof(BigIntTP));
		
		for (int i = 0; i < product->size; i += 1)
			product->coeffs[i] = empty_BigIntT(1);
	}
	
	//Zeroing out all coefficients
	else
		for (int i = 0; i < product->size; i += 1)
			copy_BigIntT(temp, product->coeffs[i]);
	
	//Iterate through both polynomials, multiply them through
	for (int a = 0; a < A->size; a += 1)
	{
		for (int b = 0; b < B->size; b += 1)
		{
			multiply_BigIntT(A->coeffs[a], B->coeffs[b], temp);
			add_BigIntT(temp, product->coeffs[a+b], temp2);
			copy_BigIntT(temp2, product->coeffs[a+b]);
		}
	}
	
	temp = free_BigIntT(temp);
	temp2 = free_BigIntT(temp2);
	
	reduce_BigPolyT(product);
	
	return 1;
}


int mod_BigPolyT(BigPolyTP const A, BigIntTP const mod, BigPolyTP residue)
/** Performs a mod operation on each element of a polynomial,
    stores the result in residue. This function assumes all arguments
		have been properly initialised. 
		Returns 1 on success, 0 otherwise. */
{
	//Preparing residue to hold what it needs to hold
	if (residue->size > A->size)
		for (int i = A->size; i < residue->size; i += 1)
			residue->coeffs[i] = free_BigIntT(residue->coeffs[i]);
		
	if (residue->size != A->size)
		residue->coeffs = realloc(residue->coeffs, (A->size)*sizeof(BigIntTP));
	
	if (residue->size < A->size)
		for (int i = residue->size; i < A->size; i += 1)
			residue->coeffs[i] = empty_BigIntT(1);
	
	residue->size = A->size;
	
	//Now, perform the modulo operation
	for (int i = 0; i < A->size; i += 1)
		mod_BigIntT(A->coeffs[i], mod, residue->coeffs[i]);
	
	reduce_BigPolyT(residue);
	
	return 1;
}


/* private */ int find_factors(BigIntTP const target,
                               BigIntTP const factor,
															 BigIntTP const carry,
															 BigIntTP const modulus,
															 BigIntTP** p,
															 BigIntTP const one)
/* Finds all the possible numbers one can multiply factor by (mod modulus),
   then add carry, to get target. Stores all the possibilities within 
	 p[0]. Returns the number of numbers in p[0], 0 if there
	 are none. */
{
	//NOTE: possibilities should always start with a size of 1
	
	int possibilityCount = 0;
	
	BigIntTP count = empty_BigIntT(1);
	BigIntTP temp2 = empty_BigIntT(1);
	BigIntTP temp3 = empty_BigIntT(1);
	
	//Check every possible number we could multiply by
	while (compare_BigIntT(count, modulus) != 0)
	{
		multiply_BigIntT(factor, count, temp2);
		add_BigIntT(temp2, carry, temp3);
		mod_BigIntT(temp3, modulus, temp2);
		
		//If we found a possibile number
		if (compare_BigIntT(temp2, target) == 0)
		{
			copy_BigIntT(count, p[0][possibilityCount]);
			possibilityCount += 1;
			
			//Resize array for next possibility
			p[0] = realloc(p[0], (possibilityCount+1)*sizeof(BigIntTP));
			p[0][possibilityCount] = empty_BigIntT(1);
		}
		
		//Increment our count
		add_BigIntT(count, one, temp2);
		copy_BigIntT(temp2, count);
	}
	
	count = free_BigIntT(count);
	temp2 = free_BigIntT(temp2);
	temp3 = free_BigIntT(temp3);
	
	return possibilityCount;
}


/* private */ bool factor_check_recurse(BigIntTP* const prod,  //The polynomial we're factoring
                                        BigIntTP* const fact,  //The factor we're checking
																	      BigIntTP* other,       //The other factor 
																				BigIntTP* carries,     //The inside of the area model
																				BigIntTP const mod,    //The modulus used
																				int size,              //Size of prod
																				int place,             //How far into the multiplication we are
																				BigIntTP const zero,
																				BigIntTP const one)
/** Recursive function to check all possible "other" factors
    given our current factor to check. Basically, it sees whether our
		factor corresponds to any other possible factor. If it does, then
		it's a valid factor. If not, it's not a factor of our polynomial. 
		Returns TRUE if fact is a factor of prod, FALSE otherwise. */
{
	//The way recursion is handled in this function can be improved a lot
	//Currently, I just want to get the logic right, then I'll optimise and
	// rewrite and stuff
	
	//Is the factor we're given a possible factor?
	bool isPossible = TRUE;
	bool recursiveCall = FALSE;
	int firstNonZeroTerm;
	int lastNonZeroTerm;
	
	int degree; //Degree of prod
	
	/*printf("Current polynomial to make:\n");
	for (int i = 0; i < size; i += 1)
	{
		printi(prod[i]);
		printf(" ");
	}
	printf("\n\n"); */
	
	//Find the first nonzero entry in our polynomial
	for (int i = 0; i < size; i += 1)
		if (compare_BigIntT(zero, fact[i]) != 0)
		{
			firstNonZeroTerm = i;
			break;
		}
		
	//Find the last nonzero term; the degree of our factor
	for (int i = size-1; i >= 0; i -= 1)
		if (compare_BigIntT(zero, fact[i]) != 0)
		{
			lastNonZeroTerm = i;
			break;
		}
		
	for (int i = size-1; i >= 0; i -= 1)
		if (compare_BigIntT(zero, prod[i]) != 0)
		{
			degree = i;
			break;
		}
		
	//The first nonzero element in prod must come at least as late as fact's
	//Otherwise, fact can't possibly be a factor
	for (int i = 0; i < size; i += 1)
		if (compare_BigIntT(zero, prod[i]) != 0)
			if (i < firstNonZeroTerm)
				return FALSE;
	
	/*printf("place: %d\n\n", place);
	printf("Factor to check: ");
	for (int i = 0; i < size; i += 1)
	{
		printi(fact[i]);
		printf(" ");
	}
	printf("\n\n");
	
	printf("Current carries: ");
	for (int i = 0; i < size; i += 1)
	{
		printi(carries[i]);
		printf(" ");
	}
	printf("\n\n"); */
	
	BigIntTP* toMultiplyBy = malloc(sizeof(BigIntTP));
	BigIntTP* newCarries   = malloc(size*sizeof(BigIntTP));
	for (int i = 0; i < size; i += 1)
		newCarries[i] = empty_BigIntT(1);
	
	BigIntTP temp  = empty_BigIntT(1);
	BigIntTP temp2 = empty_BigIntT(1);
	
	//Holds how many possible multiplication possibilities there are
	// from find_factor()
	int possibilities = 0;
	
	//Iterate over the rest of the polynomial, do the multiplication
	for (int i = place; i <= degree-lastNonZeroTerm; i += 1)
	{
		//printf("i: %d\n", i);
		//Preparing for numbers to be copied into this pointer
		toMultiplyBy[0] = empty_BigIntT(1);
		
		//Find what we need to multiply the factor by to get the correct polynomial
		possibilities = find_factors(prod[i+firstNonZeroTerm], 
		                             fact[firstNonZeroTerm], 
																 carries[i+firstNonZeroTerm], 
																 mod, 
																 &toMultiplyBy,
																 one);
		
		/*printf("Possibilities found: ");
		for (int pp = 0; pp < possibilities; pp += 1)
		{
			printi(toMultiplyBy[pp]);
			printf(" ");
		}
		printf("\n\n"); */
		
		//If there's no way our current factor can be a factor
		if (possibilities == 0)
		{
			//No reason to keep looking
			i = size;
			isPossible = FALSE;
		}
		
		//There's one term that gives what we want in the multiplication
		else if (possibilities == 1)
		{
			//Keep track of what our other factor looks like
			copy_BigIntT(toMultiplyBy[0], other[i]);
			
			//Calculate our carries, continue
			for (int c = firstNonZeroTerm; c <= lastNonZeroTerm; c += 1)
			{
				multiply_BigIntT(toMultiplyBy[0], fact[c], temp);
				add_BigIntT(temp, carries[c+i], temp2);
				mod_BigIntT(temp2, mod, carries[c+i]);
			}
		}
		
		//Recursive case
		else
		{
			recursiveCall = TRUE;
			
			//Assume all possibilities won't work until proven otherwise
			isPossible = FALSE;
			
			//We need to check if at least one possibility works
			for (int f = 0; f < possibilities; f += 1)
			{
				//Calculate new carries for specific possibility
				for (int c = firstNonZeroTerm; c < size; c += 1)
					copy_BigIntT(carries[c], newCarries[c]);
				
				for (int c = firstNonZeroTerm+i; c <= lastNonZeroTerm; c += 1)
				{
					multiply_BigIntT(toMultiplyBy[f], fact[c], temp);
					add_BigIntT(temp, newCarries[c+i], temp2);
					mod_BigIntT(temp2, mod, newCarries[c+i]);
				}
				
				//Keep track of what our other factor looks like
				copy_BigIntT(toMultiplyBy[f], other[i]);
				
				if (factor_check_recurse(prod, fact, other, newCarries, mod, size, place+1, zero, one))
				{					
					isPossible = TRUE;
					break;
				}
			}
	
			break; 
		}
		
		//Get toMultiplyBy ready for next term
		for (int j = 0; j < possibilities+1; j += 1)
			toMultiplyBy[j] = free_BigIntT(toMultiplyBy[j]);
		
		toMultiplyBy = realloc(toMultiplyBy, sizeof(BigIntTP));
		
		//Print out other so I can get a sense for what this thing looks like
		/*printf("Other factor: ");
		for (int ii = 0; ii < size; ii += 1)
		{
			printi(other[ii]);
			printf(" ");
		}
		printf("\n\n"); */
	}

	free(toMultiplyBy);
	toMultiplyBy = NULL;
	
	for (int i = 0; i < size; i += 1)
		newCarries[i] = free_BigIntT(newCarries[i]);
	free(newCarries);
	newCarries = NULL;
	
	temp  = free_BigIntT(temp);
	temp2 = free_BigIntT(temp2);
	
	//Gotta check the rest of the carries to see if they match the target polynomial
	//We only need to check this in the non-recursive case
	if ((isPossible) && (!recursiveCall))
	{
		/*printf("Checking carries...\n");
		printf("Current carries: ");
		for (int i = 0; i < size; i += 1)
		{
			printi(carries[i]);
			printf(" ");
		}
		printf("\n\n"); */
		
		//Check to see if the remaining carries finish the polynomial
		//This check only happens on a non-recursive fallout
		for (int i = degree-lastNonZeroTerm; i < size; i += 1)
			if (compare_BigIntT(prod[i], carries[i]) != 0)
			{
				isPossible = FALSE;
				break;
			}
	}
	
	//Print out other so I can get a sense for what this thing looks like
	/*printf("Other factor: ");
	for (int ii = 0; ii < size; ii += 1)
	{
		printi(other[ii]);
		printf(" ");
	}
	printf("\n\n");
		
	printf("%d\n", isPossible);
	if (isPossible)
		getchar();
	printf("~~~~~~~~~~~~~~~~~~~~~~~~~~\n"); */
	return isPossible;
}


BigPolyTP* factor_BigPolyT(BigPolyTP const A, BigIntTP const mod)
/** Factors the given BigPolyTP and returns the factors
    in a pointer. The first element in the factor will
    be a constant BigPolyTP telling how many factors there
		are total. 
		This function assumes the degree of A is at least 2. */
{
	int oneArr[1] = {1};
	
	int degree = A->size; //Holds the degree of our polynomial as we reduce it
	
	BigIntTP zero = empty_BigIntT(1);
	BigIntTP one  = new_BigIntT(oneArr, 1);
	BigIntTP temp = empty_BigIntT(1);
	BigIntTP minusOne = empty_BigIntT(1);
	subtract_BigIntT(mod, one, minusOne);
	
	BigPolyTP* factors = malloc(2*sizeof(BigPolyTP));
	BigPolyTP oneConstant = constant_BigPolyT(one);
	BigPolyTP tempPoly = empty_BigPolyT();
	BigPolyTP numOfFactors = constant_BigPolyT(zero);
	
	int numOfFactorsSmall = 0;
	factors[0] = numOfFactors;
	
	//Abstract representation of poly's coefficients
	//It'll be more efficient to use just numbers
	BigIntTP* currentPoly  = malloc((A->size)*sizeof(BigIntTP));
	BigIntTP* factorToTest = malloc((A->size)*sizeof(BigIntTP));
	
	//When checking a factor, this holds the other factor
	BigIntTP* polyRemainder = malloc((A->size)*sizeof(BigIntTP));
	
	//This holds the extra bits we generate when doing multiplication
	BigIntTP* polyCarries   = malloc((A->size)*sizeof(BigIntTP));
	
	bool allPolyTested;
	
	//Preparing the relevant pointers
	for (int i = 0; i < A->size; i += 1)
	{
		currentPoly[i] = empty_BigIntT(1);
		copy_BigIntT(A->coeffs[i], currentPoly[i]);
		
		//This allows us to skip constants as polynomials
		//(possible constant factors are handled at the end)
		if (i == 1)
			factorToTest[i] = new_BigIntT(oneArr, 1);
		else
			factorToTest[i] = empty_BigIntT(1);
		
		polyRemainder[i] = empty_BigIntT(1);
		polyCarries[i]   = empty_BigIntT(1);
	}
	
	//Iterate our current test factor
	allPolyTested = FALSE;
	while (! allPolyTested)
	{
		//Clearing out polyRemainder and polyCarries for next check
		for (int i = 0; i < A->size; i += 1)
		{
			copy_BigIntT(zero, polyRemainder[i]);
			copy_BigIntT(zero, polyCarries[i]);
		}
		
		//If we found a factor of our polynomial
		if (factor_check_recurse(currentPoly, 
		                         factorToTest,
														 polyRemainder,
														 polyCarries,
														 mod,
														 A->size,
														 0,
														 zero,
														 one))
		{
			add_BigPolyT(numOfFactors, oneConstant, tempPoly);
			numOfFactorsSmall += 1;
			copy_BigPolyT(tempPoly, factors[0]);
			
			factors = realloc(factors, (numOfFactorsSmall+2)*sizeof(BigPolyTP));
			factors[numOfFactorsSmall] = new_BigPolyT(factorToTest, A->size);
			
			/*printf("\nFactor: ");
			//Also check what the factor was that reduced our polynomial down to what it is
			printp(factors[numOfFactorsSmall]);
			printf("\n"); */
			
			//Divide polynomial by our newfound factor, continue
			//Also, start checking the lower factors again since we're
			// essentially checking a new polynomial now that we've divided
			for (int i = 0; i < A->size; i += 1)
			{
				//This allows us to start checking at 1λ again
				if (i == 0)
					copy_BigIntT(minusOne, factorToTest[0]);
				
				else
					copy_BigIntT(zero, factorToTest[i]);
				
				copy_BigIntT(polyRemainder[i], currentPoly[i]);
				//printi(currentPoly[i]);
				//printf(" ");
			}
			//printf("\n");
			
			//Calculating degree of newly reduced polynomial
			for (int i = A->size-1; i >= 0; i -= 1)
				if (compare_BigIntT(zero, currentPoly[i]) != 0)
				{
					degree = i;
					break;
				}
		}
		
		//Check to see if our currentPoly is a constant
		//If it is, just add that constant to our factors and call it a day
		allPolyTested = TRUE;
		for (int i = 0; i < A->size; i += 1)
		{
			if (i != 0)
				if (compare_BigIntT(currentPoly[i], zero) != 0)
				{
					allPolyTested = FALSE;
					break;
				}
		}
		
		//If we only have a constant left
		if (allPolyTested)
		{
			add_BigPolyT(numOfFactors, oneConstant, tempPoly);
			numOfFactorsSmall += 1;
			copy_BigPolyT(tempPoly, factors[0]);
			
			factors = realloc(factors, (numOfFactorsSmall+2)*sizeof(BigPolyTP));
			factors[numOfFactorsSmall] = new_BigPolyT(currentPoly, A->size);
		}
		
		//Else, continue iterating as normal
		else
		{
			allPolyTested = TRUE;
			for (int i = 0; i <= degree; i += 1)
			{
				add_BigIntT(factorToTest[i], one, temp);
				mod_BigIntT(temp, mod, factorToTest[i]);
				
				if (compare_BigIntT(factorToTest[i], zero) != 0)
				{
					allPolyTested = FALSE;
					break;
				}
			}
		}
	}
	
	//Now, reduce each polynomial so they have the correct degree
	for (int i = 1; i <= numOfFactorsSmall; i += 1)
		reduce_BigPolyT(factors[i]);
	
	/*
	Store our given polynomial
	Iterate through all possible polynomials of "size" smaller than current polynomial
	Check to see if its possible to divide that polynomial by the current
		- To do this:
			- find term to multiply smallest term by to get current poly's smallest term
				- This may result in multiple possibilities, or none
					- If none, ditch polynomial or check previous alternate possibility
					- If multiple, keep track of possibilities on a stack, pop them off as you check
						- We'll need a variable to keep track of what level we need to revert to when popping
						
			- keep track of extra terms created in the process
			- Repeat process for all other terms in ascending order, always keeping track of extra terms created
			- At the end, check whether our extra terms match the given polynomial
		
	  - If it is, store factor, reduce current polynomial, return to step 2.
	  - If not, continue looking.
	Give user factors when done.
	*/
	
	zero     = free_BigIntT(zero);
	temp     = free_BigIntT(temp);
	one      = free_BigIntT(one);
	minusOne = free_BigIntT(minusOne);
	
	oneConstant = free_BigPolyT(oneConstant);
	tempPoly    = free_BigPolyT(tempPoly);
	
	for (int i = 0; i < A->size; i += 1)
	{
		currentPoly[i]   = free_BigIntT(currentPoly[i]);
		factorToTest[i]  = free_BigIntT(factorToTest[i]);
		polyRemainder[i] = free_BigIntT(polyRemainder[i]);
		polyCarries[i]   = free_BigIntT(polyCarries[i]);
	}
	free(currentPoly);
	free(factorToTest);
	free(polyRemainder);
	free(polyCarries);
	currentPoly   = NULL;
	factorToTest  = NULL;
	polyRemainder = NULL;
	polyCarries   = NULL;
	
	return factors;
}
