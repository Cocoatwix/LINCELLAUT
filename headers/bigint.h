
#ifndef BIGINT_H //Header guard
#define BIGINT_H

/** Struct for holding arbitrary precision number. */
/** BigIntT structs store bunches little endian style. */
typedef struct bigint *BigIntTP;

/** Frees the memory used by a BigIntT and
    returns NULL. */
BigIntTP free_BigIntT(BigIntTP);

/** Returns a pointer to a new BigIntT struct,
    initialised with the bunches specified in the
		int* pointer and with int bunches. */
BigIntTP new_BigIntT(int* const, int);

/** Returns a pointer to a new BigIntT struct that's
    initialised with zeros. The argument specifies the
		size of the struct to create. */
BigIntTP empty_BigIntT(int);

/** Prints a BigIntT struct as an integer would be printed.
    This function does not print a newline or spaces
		around the number. */
void printi(BigIntTP);

/** Gets rid of unused bunches in the given BigIntT
    and changes the size accordingly. 
		Returns 1 on success, 0 otherwise. */
int reduce_BigIntT(BigIntTP);

/** Copies the const BigIntT to the other BigIntT. 
    This function assumes both BigIntT structs are
		initialised.
    Returns 1 on success, 0 otherwise. */
int copy_BigIntT(BigIntTP const, BigIntTP);

/** Returns -1 if the first BigIntT is smaller than the
    second, 0 if they're equal, 1 otherwise. This function
		assumes both BigIntT structs are initialised. */
int compare_BigIntT(BigIntTP const, BigIntTP const);

/** Subtracts the second BigIntT from the first BigIntT,
    stores result in the third BigIntT. Note that
		negative BigIntT structs aren't supported yet, so
		a bigger number can't be subtracted from a smaller
		number.
		This function assumes all BigIntT structs are initialised
		appropriately (third BigIntT should be big enough to
		store both other structs).
		Returns 1 on success, 0 otherwise. */
int subtract_BigIntT(BigIntTP const, BigIntTP const, BigIntTP);

/** Computes first BigIntT mod second BigIntT, stores
    result in third BigIntT, which should be big enough
		to store the two other structs if needed.
		Returns 1 on success, 0 otherwise. */
int mod_BigIntT(BigIntTP const, BigIntTP const, BigIntTP);

#endif //BIGINT_H