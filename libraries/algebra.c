
/*
A simple library for manipulating
more typical algebra things, like
polynomials.

yiff day, 2022
*/

#include <stdlib.h>
#include <stdio.h>

#include <string.h> //For MultiVarExtT extension names

#include "../headers/helper.h" //bool 

#include "../headers/bigint.h"
#include "../headers/linalg.h"
#include "../headers/modular.h" //big_num_inverse()
#include "../headers/factors.h" //poly_gcd()

const int MAXVARLEN = 20;

typedef struct bigpoly
/** Stores a polynomial using BigIntT coefficients. 
    Holds coefficients in little endian form. */
{
	int size; //Holds how many terms are in the polynomial
	BigIntTP* coeffs; //Holds the coefficients of the polynomial
}
BigPolyT, *BigPolyTP;


typedef struct bigfactors
{
	int        size;      //How many factors we're holding
	int*       exponents; //The exponent on each factor
	BigPolyTP* factors;   //The array of factors
}
BigFactorsT, *BigFactorsTP;


/* private */ typedef struct bigintdirector
/** This is the struct responsible for facilitating
    multidimensional BigIntT arrays with varying dimension. */
{
	BigIntTP* coeffs;
	struct bigintdirector** next;
	
	bool hasNext; //If TRUE, next is active. FALSE = coeffs is active.
	int size;     //How big its array is
}
BigIntDirectorT, *BigIntDirectorTP;


typedef struct multivarext
/** Stores a representation of a multivariate field extension. */
{
	BigIntDirectorTP coeffs;
	BigIntTP** extensions;   //Holds the minPolys for each of our extensions
	int* extensionSizes;     //Holds how many numbers are in each extension's defn
	
	BigIntTP mod;
	
	int numOfExtensions;    //Holds how many extensions this struct represents
	int numOfExtensionsSet; //How many extensions have we defined?
	char** extNames;        //Variable names for printing out the struct
}
MultiVarExtT, *MultiVarExtTP;


BigPolyTP free_BigPolyT(BigPolyTP p)
/** Frees the memory of a given BigPolyT. Returns NULL. */
{
	if (p != NULL)
	{
		for (int i = 0; i < p->size; i += 1)
			p->coeffs[i] = free_BigIntT(p->coeffs[i]);
		
		free(p->coeffs);
		p->coeffs = NULL;
	}
	
	return NULL;
}


BigFactorsTP free_BigFactorsT(BigFactorsTP f)
/** Frees the memory used by a BigFactorsT struct.
    Returns NULL. */
{
	if (f != NULL)
	{
		for (int i = 0; i < f->size; i += 1)
			f->factors[i] = free_BigPolyT(f->factors[i]);
		free(f->factors);
		f->factors = NULL;
		
		free(f->exponents);
		f->exponents = NULL;
	}

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


/* private */ BigIntDirectorTP free_BigIntDirectorT(BigIntDirectorTP d)
/** Frees the memory used by a BigIntDirectorT. Returns NULL. */
{
	if (d == NULL)
		return NULL;
	
	if (d->hasNext)
	{
		for (int i = 0; i < d->size; i += 1)
			d->next[i] = free_BigIntDirectorT(d->next[i]);
		
		free(d->next);
		d->next = NULL;
	}
	
	else //If the director is pointing directly to a BigIntT array
	{
		for (int i = 0; i < d->size; i += 1)
			d->coeffs[i] = free_BigIntT(d->coeffs[i]);
		
		free(d->coeffs);
		d->coeffs = NULL;
	}
	
	free(d);
	return NULL;
}


MultiVarExtTP free_MultiVarExtT(MultiVarExtTP ext)
/** Frees the memory used by a MultiVarExtT. Returns NULL. */
{
	ext->coeffs = free_BigIntDirectorT(ext->coeffs);
	
	for (int i = 0; i < ext->numOfExtensions; i += 1)
	{
		free(ext->extNames[i]);
		ext->extNames[i] = NULL;
		
		if (ext->extensions[i] != NULL) //If we try and free before all extensions have been defined
		{
			for (int elem = 0; elem < ext->extensionSizes[i]; elem += 1)
				ext->extensions[i][elem] = free_BigIntT(ext->extensions[i][elem]);
			
			free(ext->extensions[i]);
			ext->extensions[i] = NULL;
		}
	}
	free(ext->extNames);
	ext->extNames = NULL;
	
	ext->mod = free_BigIntT(ext->mod);
	
	free(ext->extensions);
	ext->extensions = NULL;
	
	free(ext->extensionSizes);
	ext->extensionSizes = NULL;
	
	free(ext);
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
	
	return p;
}


int reduce_BigPolyT(BigPolyTP p)
/** Ensures that the leading term of the BigPolyT
    is nonzero. 
		Returns 1 on success, 0 otherwise. */
{
	BigIntTP zero = empty_BigIntT(1);
	int newSize = 1;
	
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


int resize_BigPolyT(BigPolyTP p, int newSize)
/** Resizes a given BigPolyT to the given size.
    Returns 1 on success, 0 otherwise. */
{
	if (p->size == newSize)
		return 1;
	
	//If we need to dealloc a few BigIntTs
	if (p->size > newSize)
		for (int i = p->size; i < newSize; i += 1)
			p->coeffs[i] = free_BigIntT(p->coeffs[i]);
	
	p->coeffs = realloc(p->coeffs, newSize*sizeof(BigIntTP));
	
	//If we have new coefficients to initialise
	if (p->size < newSize)
		for (int i = p->size; i < newSize; i += 1)
			p->coeffs[i] = empty_BigIntT(1);
	
	p->size = newSize;
	return 1;
}


MultiVarExtTP new_MultiVarExtT(int size)
/** Allocates space for a new MultiVarExtT, returns a
    pointer to the object. Returns NULL on error. */
{
	MultiVarExtTP ext   = malloc(sizeof(MultiVarExtT));
	ext->coeffs         = malloc(sizeof(BigIntDirectorT));
	ext->extensions     = malloc(size*sizeof(BigIntTP*));
	ext->mod            = empty_BigIntT(1);
	ext->extensionSizes = malloc(size*sizeof(int));
	ext->extNames       = malloc(size*sizeof(char*));
	
	for (int i = 0; i < size; i += 1)
	{
		ext->extNames[i]    = malloc((MAXVARLEN+1)*sizeof(char));
		ext->extNames[i][0] = '\0';
		
		ext->extensions[i] = NULL;
	}
	
	//Set the proper mode for our BigIntDirectorT
	//We can only set the other members later when we add extension definitions
	if (size == 1)
		ext->coeffs->hasNext = FALSE;
	else
		ext->coeffs->hasNext = TRUE;
	
	ext->numOfExtensions    = size;
	ext->numOfExtensionsSet = 0;
	
	return ext;
}


int set_MultiVarExtT_mod(MultiVarExtTP ext, BigIntTP const bigMod)
/** Sets the modulus value for the given MultiVarExtT.
    Returns 1 on success, 0 otherwise. */
{
	return copy_BigIntT(bigMod, ext->mod);
}


int add_extension(MultiVarExtTP ext, BigIntTP* const minPoly, int sizeOfMinPoly, char* const name)
/** Adds a new extension to a MultiVarExtT.
    Returns 1 on success, 0 otherwise. */
{
	int index = ext->numOfExtensionsSet;
	int coeffsIndex[index+1]; //Holds where we are in the BigIntDirectorTP mess
	
	BigIntDirectorTP ref;     //As we add more coefficient space, this will hold our current location in the mess
	
	bool moreToAdd = TRUE;
	
	if (index == ext->numOfExtensions)
		return 0;
	
	ext->extensionSizes[index] = sizeOfMinPoly;
	
	ext->extensions[index] = malloc(sizeOfMinPoly*sizeof(BigIntTP));
	for (int i = 0; i < sizeOfMinPoly; i += 1)
	{
		ext->extensions[index][i] = empty_BigIntT(1);
		copy_BigIntT(minPoly[i], ext->extensions[index][i]);
	}
	
	if (strlen(name) <= (long unsigned int)MAXVARLEN)
		strcat(ext->extNames[index], name);
	else
		strcat(ext->extNames[index], "?"); //If the name given was too long
	
	//Now, we need to allocate coefficient space to accommodate the new extension
	for (int i = 0; i < index+1; i += 1)
		coeffsIndex[i] = 0;
	
	while (moreToAdd)
	{
		//Finding the next spot to add a coefficient to
		ref = ext->coeffs;
		for (int i = 0; i < index; i += 1)
		{
			ref = ref->next[coeffsIndex[i]];
		}
		
		if (index == ext->numOfExtensions - 1) //If we're adding BigIntTs finally
		{
			//printf("Adding new BigIntDirectorT that is in the last layer\n");
			if (ref->coeffs == NULL)
			{
				ref->coeffs  = malloc(sizeOfMinPoly*sizeof(BigIntTP));
				ref->hasNext = FALSE;
				ref->size    = sizeOfMinPoly;
			}
			
			ref->coeffs[coeffsIndex[index]] = empty_BigIntT(1);
		}
		
		else //If we're adding more BigIntDirectorTs
		{
			//printf("Adding new BigIntDirectorT that isn't in the last layer\n");
			if (ref->next == NULL)
			{
				ref->next    = malloc(sizeOfMinPoly*sizeof(BigIntDirectorTP));
				ref->coeffs  = NULL;
				ref->hasNext = TRUE;
				ref->size    = sizeOfMinPoly;
			}
			
			ref->next[coeffsIndex[index]] = malloc(sizeof(BigIntDirectorT));
			ref->next[coeffsIndex[index]]->next = NULL;
		}
		
		//Now, we need to increment our coeffsIndex so we add everything we need to
		moreToAdd = FALSE;
		for (int i = 0; i <= index; i += 1)
		{
			coeffsIndex[i] += 1;
			if (coeffsIndex[i] >= ext->extensionSizes[i])
				coeffsIndex[i] = 0;
			else
			{
				moreToAdd = TRUE;
				break;
			}
		}
	}
	
	ext->numOfExtensionsSet += 1;
	return 1;
}


int set_MultiVarExtT_coefficient(MultiVarExtTP ext, int* const coeffPos, BigIntTP const coeff)
/** Sets a particular value for the given MultiVarExtT's coefficient.
    The MultiVarExtT must be fully set before this function can be used.
		Returns 1 on success, 0 otherwise. */
{
	BigIntDirectorTP ref;
	
	if (ext->numOfExtensionsSet != ext->numOfExtensions)
		return 0;
	
	ref = ext->coeffs;
	
	for (int i = 0; i < ext->numOfExtensions-1; i += 1)
	{
		if (coeffPos[i] >= ref->size) //Making sure we don't get any IndexErrors
			return 0;
		
		ref = ref->next[coeffPos[i]];
	}
	
	if (coeffPos[ext->numOfExtensions-1] >= ref->size)
		return 0;
	
	copy_BigIntT(coeff, ref->coeffs[coeffPos[ext->numOfExtensions-1]]);
	
	return 1;
}


int reduce_MultiVarExtT(MultiVarExtTP ext)
/** Uses the extension definitions to reduce the given
    MultiVarExtT. Returns 1 on success, 0 otherwise. */
{
	if (ext->numOfExtensionsSet != ext->numOfExtensions)
		return 0;
	
	int currLoc[ext->numOfExtensions]; //Holds where we are in the coefficient array
	BigIntTP* currCoeffs = NULL;
	
	int oneArr[1] = {1};
	BigIntDirectorTP ref;
	BigIntTP leadingTermInv = empty_BigIntT(1);
	BigIntTP minPolyMult    = empty_BigIntT(1);
	BigIntTP temp           = empty_BigIntT(1);
	BigIntTP temp2          = empty_BigIntT(1);
	BigIntTP one            = new_BigIntT(oneArr, 1);
	
	bool moreToCheck;
	
	for (int i = 0; i < ext->numOfExtensions; i += 1)
		currLoc[i] = 0;
	
	
	for (int focusTerm = 0; focusTerm < ext->numOfExtensions; focusTerm += 1)
	{
		currCoeffs = realloc(currCoeffs, (ext->extensionSizes[focusTerm])*sizeof(BigIntTP));
		moreToCheck = TRUE;
		
		while (moreToCheck) //While there's still more terms to check for focusTerm's minPoly
		{
			do //Loop over focusTerm's full exponent range
			{
				//Get our set of coefficients to check
				ref = ext->coeffs;
				for (int i = 0; i < ext->numOfExtensions-1; i += 1)
					ref = ref->next[currLoc[i]];
				
				//Taking the mod before getting the reference of our coefficient
				mod_BigIntT(ref->coeffs[currLoc[ext->numOfExtensions-1]], ext->mod, temp);
				copy_BigIntT(temp, ref->coeffs[currLoc[ext->numOfExtensions-1]]);
				currCoeffs[currLoc[focusTerm]] = ref->coeffs[currLoc[ext->numOfExtensions-1]];
				
				currLoc[focusTerm] += 1;
				currLoc[focusTerm] %= ext->extensionSizes[focusTerm];
			}
			while (currLoc[focusTerm] != 0);
			
			//Print out currCoeffs so I can see what the hell's happening
			/*printf("focusTerm: %d\n", focusTerm);
			printf("currCoeffs before changing: [");
			for (int i = 0; i < ext->extensionSizes[focusTerm]-1; i += 1)
			{
				printi(currCoeffs[i]);
				printf(", ");
			}
			printi(currCoeffs[ext->extensionSizes[focusTerm]-1]);
			printf("]\n");*/
			
			//Now, currCoeffs should hold the coefficients we want to check for focusTerm's minPoly
			//Because we have a list of references, we can just change the list to change the actual coeffs
			
			//We'll reduce the highest power term
			if (! is_zero(currCoeffs[ext->extensionSizes[focusTerm]-1])) //If there's an actual leading term to reduce
			{
				//printf("Reducing...\n");
				//Find the additive inverse
				subtract_BigIntT(ext->mod, currCoeffs[ext->extensionSizes[focusTerm]-1], leadingTermInv);
				
				//Find what multiple of the minPoly is needed to reduce the leading term
				clear_BigIntT(minPolyMult);
				do
				{
					add_BigIntT(minPolyMult, one, temp2);
					copy_BigIntT(temp2, minPolyMult);
					
					multiply_BigIntT(ext->extensions[focusTerm][ext->extensionSizes[focusTerm]-1], minPolyMult, temp2);
					mod_BigIntT(temp2, ext->mod, temp);
				}
				while (compare_BigIntT(temp, leadingTermInv) != 0);
				
				//minPolyMult should now hold the correct multiple of the minPoly to eliminate the leading term
				//Let's actually reduce the terms now
				for (int i = 0; i < ext->extensionSizes[focusTerm]; i += 1)
				{
					multiply_BigIntT(ext->extensions[focusTerm][i], minPolyMult, temp);
					add_BigIntT(temp, currCoeffs[i], temp2);
					mod_BigIntT(temp2, ext->mod, currCoeffs[i]);
				}
				
				//The leading term should now be gone, meaning the expression is as reduced as it can get
				/*printf("currCoeffs after changing: [");
				for (int i = 0; i < ext->extensionSizes[focusTerm]-1; i += 1)
				{
					printi(currCoeffs[i]);
					printf(", ");
				}
				printi(currCoeffs[ext->extensionSizes[focusTerm]-1]);
				printf("]\n");*/
				
				printf("focusTerm: %d\n", focusTerm);
				printf("currLoc: [");
				for (int i = 0; i < ext->numOfExtensions-1; i += 1)
					printf("%d, ", currLoc[i]);
				printf("%d]\n", currLoc[ext->numOfExtensions-1]);
				printmve(ext);
				printf("\n");
			}
			
			//Increment currLoc
			moreToCheck = FALSE;
			for (int notFocus = 0; notFocus < ext->numOfExtensions; notFocus += 1)
			{
				if (notFocus != focusTerm)
				{
					currLoc[notFocus] += 1;
					currLoc[notFocus] %= ext->extensionSizes[notFocus];
					
					if (currLoc[notFocus] != 0) //Once all entries return to zero, we've checked all terms
					{
						moreToCheck = TRUE;
						break;
					}
				}
			}
		}
	}
	
	free(currCoeffs);
	currCoeffs = NULL;
	
	temp           = free_BigIntT(temp);
	temp2          = free_BigIntT(temp2);
	one            = free_BigIntT(one);
	leadingTermInv = free_BigIntT(leadingTermInv);
	minPolyMult    = free_BigIntT(minPolyMult);
	
	return 0;
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


BigIntTP leading_term(BigPolyTP const p)
/** Returns the coefficient on the leading term of p. 
    This function assumes p has been reduced. */
{	
	return p->coeffs[p->size - 1];
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
	
	//Just to be safe
	reduce_BigPolyT(copyTo);
	
	return 1;
}


int compare_BigPolyT(BigPolyTP const A, BigPolyTP const B)
/** Compares two polynomials to see if they're the same.
    Returns 0 if they're equal, a negative if the first
		poly is smaller than the second, a positive if the
		first poly is greater than the second. */
{
	int small = A->size < B->size ? A->size : B->size;
	int big   = A->size > B->size ? A->size : B->size;
	int r = 0;
	int theBigIsBigger = degree(A) > degree(B) ? 1 : -1;
	
	BigPolyTP theBig = degree(A) > degree(B) ? A : B;
	BigIntTP zero = empty_BigIntT(1);
	
	//For the bigger polynomial, we ensure all the entries above the smallest
	// on the other polynomial are zero
	for (int i = small; i < big; i += 1)
	{
		r = compare_BigIntT(zero, theBig->coeffs[i]);
		if (r != 0)
		{
			r = theBigIsBigger; //O_O
			break;
		}
	}
	
	//Comparing terms to see if they match
	if (r == 0)
	{
		for (int i = small-1; i >= 0; i -= 1)
		{
			r = compare_BigIntT(A->coeffs[i], B->coeffs[i]);
			if (r != 0)
				break;
		}
	}
	
	zero = free_BigIntT(zero);
	return r;
}


int set_BigPolyT(BigPolyTP p, BigIntTP* const coeffList)
/** Uses the given list of BigIntTs to set the polynomial's
    coefficients. 
		The function assumes the list is the sam length as the
		polynomial.
		Returns 1 on success, 0 otherwise. */
{
	for (int i = 0; i < p->size; i += 1)
		copy_BigIntT(coeffList[i], p->coeffs[i]);
	
	return 1;
}


int clear_BigPolyT(BigPolyTP p)
/** Sets all coefficients to zero.
    Returns 1 on success, 0 otherwise. */
{
	BigIntTP zero = empty_BigIntT(1);
	
	for (int i = 0; i < p->size; i += 1)
		copy_BigIntT(zero, p->coeffs[i]);
	
	zero = free_BigIntT(zero);
	
	return 1;
}


void printp(BigPolyTP const p)
/** Outputs a BigPolyTP to stdout. */
{
	BigIntTP zero = empty_BigIntT(1);
	bool printPlus = FALSE;
	bool isZero    = TRUE; //If p == 0, we should print a zero
	
	for (int i = 0; i < p->size; i += 1)
	{	
		//Only print if coefficient isn't zero
		if (compare_BigIntT(zero, p->coeffs[i]) != 0)
		{
			isZero = FALSE;
			
			//Prevents a + when all other terms are zero afterwards
			if (printPlus)
				printf(" + ");
			
			printPlus = TRUE;
			
			if (i == 0)
				printi(p->coeffs[0]);
			
			else if (i == 1)
			{
				printi(p->coeffs[1]);
				printf("λ");
			}
			
			else
			{
				printi(p->coeffs[i]);
				printf("λ^%d", i);
			}
		}
	}
	
	if (isZero)
		printf("0");
	
	zero = free_BigIntT(zero);
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


void printmve(MultiVarExtTP const ext)
/** Prints a MultiVarExtTP to stdout. */
{
	//Will hold which coefficient we're printing out
	int coeffPos[ext->numOfExtensionsSet];
	
	BigIntDirectorTP ref;
	BigIntTP intRef;
	
	bool printPlus;
	
	printf("Extension definitions:\n");
	for (int i = 0; i < ext->numOfExtensionsSet; i += 1)
	{
		printPlus = FALSE;
		for (int coeff = 0; coeff < ext->extensionSizes[i]; coeff += 1)
		{
			if (! is_zero(ext->extensions[i][coeff]))
			{
				if (printPlus)
					printf(" + ");
				
				printi(ext->extensions[i][coeff]);
				
				if (coeff == 1) //linear
					printf("(%s)", ext->extNames[i]);
				
				else if (coeff > 1)
					printf("(%s)^%d", ext->extNames[i], coeff);
				
				printPlus = TRUE;
			}
		}
		
		printf(" = 0\n");
	}
	
	if (ext->numOfExtensionsSet == ext->numOfExtensions) //If we can actually print the expression
	{
		for (int i = 0; i < ext->numOfExtensions; i += 1)
			coeffPos[i] = 0;
		
		printPlus = FALSE;
		
		//while we still have coefficients to iterate through
		while (coeffPos[ext->numOfExtensions-1] < ext->extensionSizes[ext->numOfExtensions-1])
		{
			//Get the next coefficient to print
			ref = ext->coeffs;
			for (int i = 0; i < ext->numOfExtensions-1; i += 1)
			{
				ref = ref->next[coeffPos[i]];
			}
			
			intRef = ref->coeffs[coeffPos[ext->numOfExtensions-1]];
			
			if (! is_zero(intRef))
			{
				if (printPlus)
					printf(" + ");
				
				printi(intRef);
				printPlus = TRUE;
				
				//Print appropriate extension names
				for (int i = 0; i < ext->numOfExtensions; i += 1)
				{
					if (coeffPos[i] == 1)
						printf("(%s)", ext->extNames[i]);
					
					else if (coeffPos[i] > 1)
						printf("(%s)^%d", ext->extNames[i], coeffPos[i]);
				}
			}
			
			//Increment coeffPos to next coefficient
			for (int i = 0; i < ext->numOfExtensions; i += 1)
			{
				coeffPos[i] += 1;
				if ((coeffPos[i] >= ext->extensionSizes[i]) &&
				    (i != ext->numOfExtensions-1))
					coeffPos[i] = 0;
				else
					break;
			}
		}
	}
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
				fprintf(file, "λ");
			}
			
			else
			{
				fprinti(file, p->coeffs[i]);
				fprintf(file, "λ^%d", i);
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


int divide_BigPolyT(BigPolyTP const a, BigPolyTP const b, BigPolyTP quotient, BigPolyTP remainder, BigIntTP const mod)
/** Divides two BigPolyTs, stores the quotient in the third BigPolyT,
    stores the remainder in the forth BigPolyT. The BigIntTP is the 
		modulus used.
		This function assumes a and b are reduced.
		Returns 1 on success, 0 otherwise. */
{
	int oneArr[1] = {1};
	
	BigPolyTP tempPoly;
	BigIntTP  temp;
	BigIntTP  temp2;
	BigIntTP  zero;
	BigIntTP  one;
	BigIntTP  negOne;
	BigIntTP  bLeadingInv; //Holds the inverse of b's leading entry's coefficient
	
	BigIntTP  tempDivCoeff; //Just to simplify expressions
	int       tempDivIndex; //Same as above
	
	BigIntTP* qCoeffs;    //Holds the coefficients for the quotient
	BigIntTP* bCoeffs;    //Holds the coefficients of the polynomial b
	BigIntTP* tempCoeffs; //Holds the working set of coefficients as we do the division
	
	int quotientSize = a->size - b->size + 1;
	
	bool bIsZero = TRUE;
	
	//First, let's check and see if we even need to do any division
	if (a->size < b->size)
	{
		if (remainder != NULL)
			copy_BigPolyT(a, remainder);
		
		zero = empty_BigIntT(1);
		tempPoly = constant_BigPolyT(zero);
		
		if (quotient != NULL)
			copy_BigPolyT(tempPoly, quotient);
		
		tempPoly = free_BigPolyT(tempPoly);
		zero = free_BigIntT(zero);
		
		return 1;
	}
	
	//If b == 0, exit
	temp = empty_BigIntT(1);
	for (int i = 0; i < b->size; i += 1)
	{
		mod_BigIntT(b->coeffs[i], mod, temp);
		
		if (! is_zero(temp))
		{
			bIsZero = FALSE;
			break;
		}
	}
		
	if (bIsZero) //Zero division error
	{
		temp = free_BigIntT(temp);
		return 0;
	}
	
	//Let's check to see if quotient is the correct size
	if ((quotient != NULL) && (quotient->size != quotientSize))
		resize_BigPolyT(quotient, quotientSize);
	
	if (remainder != NULL)
		resize_BigPolyT(remainder, a->size);
	
	qCoeffs = malloc(quotientSize*sizeof(BigIntTP));
	for (int i = 0; i < quotientSize; i += 1)
		qCoeffs[i] = empty_BigIntT(1);
	
	tempCoeffs = extract_coefficients(a);
	bCoeffs    = extract_coefficients(b);
	
	temp2       = empty_BigIntT(1);
	bLeadingInv = empty_BigIntT(1);
	big_num_inverse(bCoeffs[b->size-1], mod, bLeadingInv);
	
	one = new_BigIntT(oneArr, 1);
	negOne = empty_BigIntT(1);
	subtract_BigIntT(mod, one, negOne);
	
	//Now, iterate through tempCoeffs, find multiples of b that
	// annihilate the coefficients in it, store them in qCoeffs
	for (int c = a->size-1; c > a->size-1 - quotientSize; c -= 1) //Check this loop if things go wrong
	{
		//Calculate coefficient needed for next step in long division
		multiply_BigIntT(tempCoeffs[c], bLeadingInv, temp);
		mod_BigIntT(temp, mod, qCoeffs[quotientSize - a->size + c]); //This line might also be wrong
		
		tempDivCoeff = qCoeffs[quotientSize - a->size + c];
		tempDivIndex = quotientSize - a->size + c; //This tells us what power of variable the coeff is attached to
		
		//Now qCoeffs[quotientSize - a->size + c] holds the coefficient for this long division step
		//We need to see how this coefficient will affect the rest of the polynomial, not just the leading term
		//Loop through the rest of tempCoeffs, update as necessary
		for (int u = 0; u < b->size; u += 1)
		{
			multiply_BigIntT(tempDivCoeff, bCoeffs[u], temp);
			multiply_BigIntT(temp, negOne, temp2);
			add_BigIntT(temp2, tempCoeffs[tempDivIndex + u], temp);
			mod_BigIntT(temp, mod, tempCoeffs[tempDivIndex + u]);
		}
		
		//Now, tempCoeffs should be ready for the next division cycle
		//I think I can just let it loop now
	}
	
	//Set quotient and remainder
	if (quotient != NULL)
	{
		set_BigPolyT(quotient, qCoeffs);
		reduce_BigPolyT(quotient);
	}
	
	if (remainder != NULL)
	{
		set_BigPolyT(remainder, tempCoeffs);
		reduce_BigPolyT(remainder);
	}
	
	//It feels really gross to free multiple different arrays
	// like this.
	for (int i = 0; i < quotientSize; i += 1)
		qCoeffs[i] = free_BigIntT(qCoeffs[i]);
	free(qCoeffs);
	qCoeffs = NULL;
	
	for (int i = 0; i < a->size; i += 1)
		tempCoeffs[i] = free_BigIntT(tempCoeffs[i]);
	free(tempCoeffs);
	tempCoeffs = NULL;
	
	for (int i = 0; i < b->size; i += 1)
		bCoeffs[i] = free_BigIntT(bCoeffs[i]);
	free(bCoeffs);
	bCoeffs = NULL;
	
	bLeadingInv = free_BigIntT(bLeadingInv);
	temp        = free_BigIntT(temp);
	one         = free_BigIntT(one);
	negOne      = free_BigIntT(negOne);
	temp2       = free_BigIntT(temp2);
	
	return 1;
}


int pow_BigPolyT(BigPolyTP const p, BigIntTP const pow, BigPolyTP exp)
/** Calculates the first BigPolyT raised to the given BigIntT.
    Stores the result in the second BigPolyT given (which is
		assumed to have been initialised).
		Returns 1 on success, 0 otherwise. */
{
	//This is my attempt at implementing exponentiation by squaring
	int returnVal = 0;
	int numArr[1] = {1};
	
	BigIntTP temp    = empty_BigIntT(1);
	BigIntTP tempPow = empty_BigIntT(1);
	BigIntTP zero    = empty_BigIntT(1);
	BigIntTP one     = new_BigIntT(numArr, 1);
	BigIntTP two;
	
	BigPolyTP onePoly  = constant_BigPolyT(one);
	BigPolyTP tempPoly = empty_BigPolyT();
	
	BigPolyTP* results = NULL; //Holds a lost of polynomials we need to multiply together to get exp
	int numOfResults = 0;
	
	numArr[0] = 2;
	two = new_BigIntT(numArr, 1);
	
	if (compare_BigIntT(pow, zero) == 0)
	{
		copy_BigPolyT(onePoly, exp);
		returnVal = 1;
	}
	
	else if (compare_BigIntT(pow, one) == 0)
	{
		copy_BigPolyT(p, exp);
		returnVal = 1;
	}
	
	else
	{
		copy_BigIntT(pow, tempPow);
		numOfResults += 1;
		results = realloc(results, sizeof(BigPolyTP));
		results[0] = empty_BigPolyT();
		
		copy_BigPolyT(p, results[0]);
	
		while (compare_BigIntT(tempPow, one) > 0) //Compute the power
		{
			mod_BigIntT(tempPow, two, temp);
			
			//If tempPow is even
			if (compare_BigIntT(temp, zero) == 0)
			{
				//x^{2n} = (x^{2})^{n}
				divide_BigIntT(tempPow, two, temp);
				copy_BigIntT(temp, tempPow);
			}
			
			//If tempPow is odd
			else
			{
				//x^{2n+1} = x(x^{2})^{n}
				numOfResults += 1;
				results = realloc(results, numOfResults*sizeof(BigPolyTP));
				results[numOfResults-1] = empty_BigPolyT();
				copy_BigPolyT(results[0], results[numOfResults-1]);
				
				subtract_BigIntT(tempPow, one, temp);
				divide_BigIntT(temp, two, tempPow);
			}
			
			//Square our current result
			multiply_BigPolyT(results[0], results[0], tempPoly);
			copy_BigPolyT(tempPoly, results[0]);
		}
		
		//Now, we simply multiply all the results together to get exp
		copy_BigPolyT(onePoly, exp);
		for (int i = 0; i < numOfResults; i += 1)
		{
			multiply_BigPolyT(exp, results[i], tempPoly);
			copy_BigPolyT(tempPoly, exp);
		}
			
		returnVal = 1;
	}
	
	for (int i = 0; i < numOfResults; i += 1)
		results[i] = free_BigPolyT(results[i]);
	free(results);
	results = NULL;
	
	onePoly  = free_BigPolyT(onePoly);
	tempPoly = free_BigPolyT(tempPoly);
	
	tempPow = free_BigIntT(tempPow);
	one     = free_BigIntT(one);
	two     = free_BigIntT(two);
	zero    = free_BigIntT(zero);
	temp    = free_BigIntT(temp);
	
	return returnVal;
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


int diff_BigPolyT(BigPolyTP const p, BigPolyTP dp)
/** Differentiates p, stores result in dp.
    Returns 1 on success, 0 othrwise. */
{
	BigIntTP* coeffs = extract_coefficients(p);
	BigIntTP  temp   = empty_BigIntT(1);
	BigIntTP  count  = empty_BigIntT(1);
	
	int oneArr[1] = {1};
	BigIntTP one = new_BigIntT(oneArr, 1);
	
	//Derivative calculation
	for (int i = 0; i < p->size; i += 1)
	{
		multiply_BigIntT(coeffs[i], count, temp);
		copy_BigIntT(temp, coeffs[i]);
		
		add_BigIntT(count, one, temp);
		copy_BigIntT(temp, count);
	}
	
	//Constants get sent to zero under a derivative calculation
	coeffs[0] = free_BigIntT(coeffs[0]);
	coeffs += 1;
	
	resize_BigPolyT(dp, p->size-1);
	set_BigPolyT(dp, coeffs);
	
	for (int i = 0; i < p->size-1; i += 1)
		coeffs[i] = free_BigIntT(coeffs[i]);
	coeffs -= 1;
	free(coeffs);
	coeffs = NULL;
	
	temp  = free_BigIntT(temp);
	one   = free_BigIntT(one);
	count = free_BigIntT(count);
	
	return 1;
}


BigPolyTP* old_factor_BigPolyT(BigPolyTP const A, BigIntTP const mod)
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


BigFactorsTP factor_BigPolyT(BigPolyTP const p, BigIntTP const mod)
/** Factors polynomials using our brand-new knowledge gained
    from Wikipedia! 
		Returns a BigFactorsT, representing the factorisation
		of the given polynomial.*/
{
	//en.wikipedia.org/wiki/Factorization_of_polynomials_over_finite_fields#Factoring_algorithms
	
	//The factors in this array are all raised to an
	// exponent equal to its corresponding exponent in squareFreeFactorsExponents.
	BigPolyTP* squareFreeFactors = NULL;
	int numOfSquareFreeFactors = 0;
	
	//These help keep track of the correct exponents for the 
	// factors in squareFreeFactors
	int  previousExponentMultiple = 1;
	int  exponentMultiple = 1;
	int* squareFreeFactorsExponents = NULL;
	
	//Holds polynomials which are all products of same-degree irreducible polynomials
	BigPolyTP* distinctDegreeFactors = NULL;
	int numOfDistinctDegreeFactors = 0;
	bool split; //Says whether, during distinct-degree factorisation, the polynomial ever got factored
	
	BigPolyTP currSFF; //Just a pointer, no need to free
	int* distinctDegreeFactorsExponents = NULL; //The degree given by the distinct-degree factorisation
	int* ddfsffExponents = NULL; //The exponents carried over from the squareFreeFactors
	
	//Holds all irreducible factors of our given polynomial
	BigFactorsTP equalDegreeFactorisation;
	
	BigPolyTP diffP;
	BigPolyTP monicP;
	BigPolyTP repeatedFactors;
	BigPolyTP allFactorsNoMult;
	BigPolyTP factorsWithSpecificMult;
	BigPolyTP isolateFactorsWithSpecificMult;
	
	BigPolyTP GCD;
	
	BigIntTP  leadingTermInv = NULL;
	BigPolyTP invPoly = NULL;
	
	BigIntTP counter;
	int intCounter;
	
	int       tempCoeffsSize = 0;
	BigIntTP  temp;
	BigIntTP* tempCoeffs = NULL;
	BigPolyTP tempPoly   = NULL;
	BigPolyTP tempPoly2  = NULL;
	
	int oneArr[1] = {1};
	//int twoArr[1] = {2};
	BigIntTP one    = new_BigIntT(oneArr, 1);
	BigIntTP negOne = empty_BigIntT(1);
	//BigIntTP two    = new_BigIntT(twoArr, 1);
	BigPolyTP onePoly = constant_BigPolyT(one);
	
	//First, let's make sure p is monic
	monicP = empty_BigPolyT();
	copy_BigPolyT(p, monicP);
	reduce_BigPolyT(monicP);
	
	if (compare_BigIntT(leading_term(monicP), one) != 0)
	{
		leadingTermInv = empty_BigIntT(1);
		big_num_inverse(leading_term(monicP), mod, leadingTermInv);
		invPoly = constant_BigPolyT(leadingTermInv);
		
		multiply_BigPolyT(monicP, invPoly, tempPoly);
		mod_BigPolyT(tempPoly, mod, monicP);
	}
	
	//Initialise terms we need for square-free factorization
	tempPoly  = empty_BigPolyT();
	tempPoly2 = empty_BigPolyT();
	diffP     = empty_BigPolyT();
	
	repeatedFactors  = empty_BigPolyT();
	allFactorsNoMult = empty_BigPolyT();
	
	factorsWithSpecificMult        = empty_BigPolyT();
	isolateFactorsWithSpecificMult = empty_BigPolyT();
	
	GCD = empty_BigPolyT();
	
	temp    = empty_BigIntT(1);
	counter = empty_BigIntT(1);
	
	diff_BigPolyT(monicP, tempPoly);
	mod_BigPolyT(tempPoly, mod, diffP);
	
	//Repeated factors w/ one less multiplicity
	poly_gcd(monicP, diffP, repeatedFactors, mod, NULL, NULL);
	
	//All factors except those with multiplicity a multiple of mod, no multiplicities
	divide_BigPolyT(monicP, repeatedFactors, allFactorsNoMult, NULL, mod);
	
	//This do loop lets us run the algorithm without recursion
	do
	{
		//Now, we isolate all the factors with the same multiplicity
		// in allFactorsNoMult, and shove them in squareFreeFactors
	
		for (int i = 0; compare_BigPolyT(onePoly, allFactorsNoMult) != 0; i += 1)
		{
			poly_gcd(allFactorsNoMult, repeatedFactors, isolateFactorsWithSpecificMult, mod, NULL, NULL);
			divide_BigPolyT(allFactorsNoMult, isolateFactorsWithSpecificMult, factorsWithSpecificMult, NULL, mod);
			
			//Store factors with current multiplicity (i+1)
			if (compare_BigPolyT(factorsWithSpecificMult, onePoly) != 0) //Prevents a 1 from being stored
			{
				numOfSquareFreeFactors += 1;
				
				squareFreeFactors = realloc(squareFreeFactors, numOfSquareFreeFactors*sizeof(BigPolyTP));
				squareFreeFactors[numOfSquareFreeFactors-1] = empty_BigPolyT();
				copy_BigPolyT(factorsWithSpecificMult, squareFreeFactors[numOfSquareFreeFactors-1]);
				
				squareFreeFactorsExponents = realloc(squareFreeFactorsExponents, numOfSquareFreeFactors*sizeof(int));
				squareFreeFactorsExponents[numOfSquareFreeFactors-1] = (i+1)*exponentMultiple;
			}
			
			//Remember that the multiplicity goes up each element in squareFreeFactors
			// We don't actually compute the power here; we don't really need to
			// Just keep track of it!
			
			//Get rid of factors in allFactorsNoMult which we've already accounted for
			copy_BigPolyT(isolateFactorsWithSpecificMult, allFactorsNoMult);
			
			//Switch our attention to the next multiplicity
			divide_BigPolyT(repeatedFactors, isolateFactorsWithSpecificMult, tempPoly, NULL, mod);
			copy_BigPolyT(tempPoly, repeatedFactors);
		}
		
		//Now, repeatedFactors now holds factors with multiplicity a multiple of mod
		//We need to isolate these using a different method
		if (compare_BigPolyT(repeatedFactors, onePoly) != 0)
		{
			//Plug x^{1/p} into repeatedFactors
			tempCoeffsSize = 0;
			clear_BigIntT(counter);
			for (int j = 0; j < repeatedFactors->size; j += 1)
			{
				//Check to see if our term has an exponent a multiple of mod
				mod_BigIntT(counter, mod, temp);
				if (is_zero(temp))
				{
					//Shove coefficient into array
					tempCoeffsSize += 1;
					tempCoeffs = realloc(tempCoeffs, tempCoeffsSize*sizeof(BigIntTP));
					tempCoeffs[tempCoeffsSize-1] = repeatedFactors->coeffs[j];
				}
				
				add_BigIntT(counter, one, temp);
				copy_BigIntT(temp, counter);
			}
			
			//Make new polynomial with array of coeffs
			//Run algorithm again
			resize_BigPolyT(monicP, tempCoeffsSize);
			set_BigPolyT(monicP, tempCoeffs);
			
			//Also, get our new exponent multiple
			copy_BigIntT(one, counter);
			while (compare_BigIntT(counter, mod) != 0)
			{
				exponentMultiple += previousExponentMultiple;
				
				add_BigIntT(counter, one, temp);
				copy_BigIntT(temp, counter);
			}
			previousExponentMultiple = exponentMultiple;
			
			//Get ready for next iteration through the square-free factorisation algorithm
			diff_BigPolyT(monicP, tempPoly);
			mod_BigPolyT(tempPoly, mod, diffP);
			
			//Repeated factors w/ one less multiplicity
			poly_gcd(monicP, diffP, repeatedFactors, mod, NULL, NULL);
			
			//All factors except those with multiplicity a multiple of mod, no multiplicities
			divide_BigPolyT(monicP, repeatedFactors, allFactorsNoMult, NULL, mod);
		}
	}
	while (compare_BigPolyT(onePoly, allFactorsNoMult) != 0);
	
	printf("\nsquareFreeFactors: \n");
	for (int i = 0; i < numOfSquareFreeFactors; i += 1)
	{
		printf("(");
		printp(squareFreeFactors[i]);
		printf(")^%d", squareFreeFactorsExponents[i]);
	}
	printf("\n\n");
	
	subtract_BigIntT(mod, one, temp);
	copy_BigIntT(temp, negOne);
	
	//Now we have a square-free factorisation of our initial polynomial p
	//I guess we have to factor each bit now
	
	//Note that the square-free polynomials need to be monic, but
	// I'm pretty sure they always will be
	for (int sfFactor = 0; sfFactor < numOfSquareFreeFactors; sfFactor += 1)
	{
		split = FALSE;
		currSFF = squareFreeFactors[sfFactor];

		intCounter = 1;
		while (degree(currSFF) >= 2*intCounter)
		{
			//We need to compute x^{p^counter} - x mod currSFF
			
			//Calculate x^p
			if (intCounter == 1)
			{
				//We'll have to create a new polynomial here
				if (size(mod) > 1)
				{
					fprintf(stderr, "Given modulus is too big for the factoring function to handle." \
													" Any following behaviour is undefined.\n");
				}

				resize_BigPolyT(tempPoly2, extract_bunch(mod, 0)+1);
				clear_BigPolyT(tempPoly2);
				copy_BigIntT(one, tempPoly2->coeffs[extract_bunch(mod, 0)]); // x^p
				copy_BigIntT(negOne, tempPoly2->coeffs[1]);                  //-x
				
				//tempPoly2 should now hold x^p - x
			}
			
			else
			{
				add_BigIntT(tempPoly->coeffs[1], one, temp);
				mod_BigIntT(temp, mod, tempPoly->coeffs[1]);
				
				//tempPoly should now be x^{p^{counter-1}}
				//Let's make it x^{p^{counter}} - x
				pow_BigPolyT(tempPoly, mod, tempPoly2);
				add_BigIntT(tempPoly2->coeffs[1], negOne, temp);
				mod_BigIntT(temp, mod, tempPoly2->coeffs[1]);
				
				//tempPoly2 should now be x^{p^{counter}} - x
			}
			
			//Dividing here is used to take a polynomial modulus
			divide_BigPolyT(tempPoly2, currSFF, NULL, tempPoly, mod);
			
			//Now that we have x^{p^{counter}} - x, let's actually do the computation
			poly_gcd(tempPoly, currSFF, GCD, mod, NULL, NULL);
			if (compare_BigPolyT(GCD, onePoly) != 0)
			{
				split = TRUE;
				
				//Keeping track of factors and exponents here
				numOfDistinctDegreeFactors += 1;
				distinctDegreeFactors = realloc(distinctDegreeFactors, numOfDistinctDegreeFactors*sizeof(BigPolyTP));
				distinctDegreeFactors[numOfDistinctDegreeFactors-1] = empty_BigPolyT();
				copy_BigPolyT(GCD, distinctDegreeFactors[numOfDistinctDegreeFactors-1]);
				
				distinctDegreeFactorsExponents = realloc(distinctDegreeFactorsExponents, numOfDistinctDegreeFactors*sizeof(int));
				distinctDegreeFactorsExponents[numOfDistinctDegreeFactors-1] = intCounter;
				
				ddfsffExponents = realloc(ddfsffExponents, numOfDistinctDegreeFactors*sizeof(int));
				ddfsffExponents[numOfDistinctDegreeFactors-1] = squareFreeFactorsExponents[sfFactor];
				
				divide_BigPolyT(currSFF, GCD, tempPoly2, NULL, mod);
				copy_BigPolyT(tempPoly2, currSFF);
				reduce_BigPolyT(currSFF);
			}
			
			intCounter += 1;
		}
		
		//Add the remaining bit of currSFF into our list of factors
		if ((compare_BigPolyT(currSFF, onePoly) != 0) || (!split))
		{
			numOfDistinctDegreeFactors += 1;
			distinctDegreeFactors = realloc(distinctDegreeFactors, numOfDistinctDegreeFactors*sizeof(BigPolyTP));
			distinctDegreeFactors[numOfDistinctDegreeFactors-1] = empty_BigPolyT();
			copy_BigPolyT(currSFF, distinctDegreeFactors[numOfDistinctDegreeFactors-1]);

			distinctDegreeFactorsExponents = realloc(distinctDegreeFactorsExponents, numOfDistinctDegreeFactors*sizeof(int));
			
			if (compare_BigPolyT(currSFF, onePoly) != 0)
				distinctDegreeFactorsExponents[numOfDistinctDegreeFactors-1] = degree(currSFF); //Should be okay since currSFF is always reduced
			else
				distinctDegreeFactorsExponents[numOfDistinctDegreeFactors-1] = 1; 
			
			ddfsffExponents = realloc(ddfsffExponents, numOfDistinctDegreeFactors*sizeof(int));
			ddfsffExponents[numOfDistinctDegreeFactors-1] = squareFreeFactorsExponents[sfFactor];
		}
	}
	
	//Now, our polynomial should be split into distinct-degree factors
	printf("\ndistinctDegreeFactors:\n");
	for (int i = 0; i < numOfDistinctDegreeFactors; i += 1)
	{
		printf("(");
		printp(distinctDegreeFactors[i]);
		printf(")^%d", ddfsffExponents[i]);
	}
	printf("\ndistinctDegreeFactorsExponents:\n");
	for (int i = 0; i < numOfDistinctDegreeFactors; i += 1)
		printf("(%d)", distinctDegreeFactorsExponents[i]);
	printf("\n\n");
	
	//Now, our polynomial should be split into distinct-degree factors.
	//All that's left is to factor the reducible factors
	
	for (int ddFactor = 0; ddFactor < numOfDistinctDegreeFactors; ddFactor += 1)
	{
		//If we find a reducible factor
		if (degree(distinctDegreeFactors[ddFactor]) != distinctDegreeFactorsExponents[ddFactor])
		{
			
		}
	}
	

	free(tempCoeffs);
	tempCoeffs = NULL;
	
	diffP     = free_BigPolyT(diffP);
	monicP    = free_BigPolyT(monicP);
	invPoly   = free_BigPolyT(invPoly);
	tempPoly  = free_BigPolyT(tempPoly);
	tempPoly2 = free_BigPolyT(tempPoly2);
	
	onePoly = free_BigPolyT(onePoly);
	
	repeatedFactors                = free_BigPolyT(repeatedFactors);
	allFactorsNoMult               = free_BigPolyT(allFactorsNoMult);
	factorsWithSpecificMult        = free_BigPolyT(factorsWithSpecificMult);
	isolateFactorsWithSpecificMult = free_BigPolyT(isolateFactorsWithSpecificMult);
	
	GCD = free_BigPolyT(GCD);
	
	one            = free_BigIntT(one);
	negOne         = free_BigIntT(negOne);
	//two            = free_BigIntT(two);
	temp           = free_BigIntT(temp);
	leadingTermInv = free_BigIntT(leadingTermInv);
	counter        = free_BigIntT(counter);
	
	free(squareFreeFactorsExponents);
	squareFreeFactorsExponents = NULL;
	//REMEMBER TO FREE SQUAREFREEFACTORS AT SOME POINT
	
	return NULL;
}
