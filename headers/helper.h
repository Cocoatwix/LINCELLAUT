
#ifndef HELPER_H
#define HELPER_H

typedef enum boolean {FALSE, TRUE} bool;

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