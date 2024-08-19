
#ifndef HELPER_H
#define HELPER_H

typedef enum boolean {FALSE, TRUE} bool;
typedef struct dict *DictionaryTP;

DictionaryTP free_DictionaryT(DictionaryTP);

/*
 * Valid types for a DictionaryT:
 *  - STR
 *  - INT
 *  - STRARR (not yet implemented)
 */

/** Initialise a new DictionaryT with a given key, value, and type. */
//                                   key        value         type
DictionaryTP new_DictionaryT(const char*, const void*, const char*);

/** Sifts through a DictionaryTP array to find the first instance
    of a dictionary with the given key. Returns found DictionaryT 
	on success, NULL otherwise. */
//                                                         length     key
DictionaryTP search_DictionaryT_array(const DictionaryTP*, int, const char*);

/** Returns the value a DictionaryT is holding. */
void* value_of_DictionaryT(const DictionaryTP);

/** Returns the auxillary value a DictionaryT is holding. 
    Returns NULL if the DictionaryT has no auxillary value. */
void* aux_of_DictionaryT(const DictionaryTP);

/** Returns the number of digits a given
    positive integer has. */
int num_digits(int);

/** Appends an integer to the given string. */
void append_int(char*, int);

/** Sorts the given list of integers from least to greatest
    using the quicksort algorithm. */
//              list  size
void quick_sort(int*, int);

#endif //HELPER_H