
#ifndef BIGINT_H //Header guard
#define BIGINT_H

#include <stdio.h> //Prevents compiler from complaining about FILE*

const int MAXBUNCH;
const int MAXBUNCHDIGITS;

/** Struct for holding arbitrary precision number. */
/** BigIntT structs store bunches little endian style. */
typedef struct bigint *BigIntTP;

/** Frees the memory used by a BigIntT and
    returns NULL. */
BigIntTP free_BigIntT(BigIntTP);

/** Returns a pointer to a new BigIntT struct,
    initialised with the bunches specified in the
		int* pointer and with int bunches. */
BigIntTP new_BigIntT(const int*, int);

/** Returns a pointer to a new BigIntT struct that's
    initialised with zeros. The argument specifies the
		size of the struct to create. */
BigIntTP empty_BigIntT(int);

/** Returns 1 if the given BigIntT is zero,
    0 otherwise. */
int is_zero(const BigIntTP);

/** Returns a pointer to a BigIntTP 2D array. Each
    BigIntTP will be initialised to zero.
	Returns NULL on error. */
BigIntTP** new_BigIntT_array(int, int);

/** Frees all BigIntTs in a BigIntTP 2D array, then
    frees the array itself. Returns NULL. */
BigIntTP** free_BigIntT_array(BigIntTP**, int, int);

/** Reads in a string and stores its numerical value
    in the given BigIntTP. The BigIntTP is assumed NOT
		to have been initialised, only declared.
		Returns 1 on success, 0 otherwise. */
int strtoBIT(const char*, BigIntTP*);

/** Returns the size of the BigIntT passed. */
int size(const BigIntTP);

/** Returns the integer stored at the given bunch in the
    given BigIntT. Returns -1 on error. */
int extract_bunch(const BigIntTP, int);

/** Sets the value of a bunch to the given number.
    Returns 1 on success, 0 otherwise. */
//            num       pos  bunch value
int set_bunch(BigIntTP, int, int);

/** Appends a BigIntTP to a string.
    This function assumes there's enough space in the string to
	hold the BigIntT.
    Returns 1 on success, 0 otherwise. */
int append_BigIntT(char*, const BigIntTP);

/** Prints a BigIntT struct as an integer would be printed.
    This function does not print a newline or spaces
		around the number. */
void printi(const BigIntTP);

/** Same as printi, except zero padding is added 
    to the beginning of the number. */
void printi_pad(const BigIntTP);

/** Same as printi, but outputs to a file stream. */
void fprinti(FILE*, BigIntTP);

/** Gets rid of unused bunches in the given BigIntT
    and changes the size accordingly. 
		Returns 1 on success, 0 otherwise. */
int reduce_BigIntT(BigIntTP);

/** Resets a BigIntT number to 0.
    Returns 1 on success, 0 otherwise. */
int clear_BigIntT(BigIntTP);

/** Changes the size of the provided BigIntT and
    reallocates its memory appropriately. 
		Returns 1 on success, 0 otherwise. */
int resize_BigIntT(BigIntTP, int);

/** Copies the const BigIntT to the other BigIntT. 
    This function assumes both BigIntT structs are
	initialised.
    Returns 1 on success, 0 otherwise. */
int copy_BigIntT(const BigIntTP, BigIntTP);

/** Returns negative if the first BigIntT is smaller than the
    second, 0 if they're equal, positive otherwise. This function
		assumes both BigIntT structs are initialised. */
int compare_BigIntT(const BigIntTP, const BigIntTP);

/** Adds the first two BigIntT structs together and
    stores the result in the third BigIntT. This function
		assumes the third BigIntT has been initialised.
		Returns 1 on success, 0 otherwise. */
int add_BigIntT(const BigIntTP, const BigIntTP, BigIntTP);

/** Subtracts the second BigIntT from the first BigIntT,
    stores result in the third BigIntT. Note that
		negative BigIntT structs aren't supported yet, so
		a bigger number can't be subtracted from a smaller
		number.
		This function assumes all BigIntT structs are initialised.
		Returns 1 on success, 0 otherwise. */
int subtract_BigIntT(const BigIntTP, const BigIntTP, BigIntTP);

/** Multiples the first BigIntTP but the second, and stores the
    product in the third. This function assumes all three 
		BigIntTP structs have been initialised.
		Returns 1 on success, 0 otherwise. */
int multiply_BigIntT(const BigIntTP, const BigIntTP, BigIntTP);

/** Divides the first BigIntT by the second, and stores the
    result in the third BigIntT. This function assumes
		all BigIntT structs have been initialised.
		
		Note that this function performs floor division, so any
		remainder after division will be discarded. Use
		mod_BigIntT to obtain the remainder.
		
		Returns 1 on success, 0 otherwise. */
int divide_BigIntT(const BigIntTP, const BigIntTP, BigIntTP);

/** Computes first BigIntT mod second BigIntT, stores
    result in third BigIntT, which assumes the third
	BigInt has been initialised.
	Returns 1 on success, 0 otherwise. */
int mod_BigIntT(const BigIntTP, const BigIntTP, BigIntTP);

/** Returns an array pf all the possible factors of the given BigIntTP. 
    The first number tells you how many factors there are. */
BigIntTP* divisors_of_BigIntT(const BigIntTP);


#endif //BIGINT_H