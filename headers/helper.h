
#ifndef HELPER_H
#define HELPER_H

typedef enum boolean {FALSE, TRUE} bool;

/** Returns the number of digits a given
    positive integer has. */
int num_digits(int);

/** Appends an integer to the given string. */
void append_int(char*, int);

#endif //HELPER_H