/*
 * include/bits.h
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2000-2001
 *
 * This file defines various routines to get / set bits of a IEEE floating
 * point number.  This is used by the library for debugging purposes.
 */

#ifndef _QD_BITS_H
#define _QD_BITS_H

#include <iostream>
#include <qd/qd_config.h>

/* Returns the exponent of the double precision number.
   Returns INT_MIN is x is zero, and INT_MAX if x is INF or NaN. */
int get_double_expn(double x);

/* Prints 
     SIGN  EXPN  MANTISSA
   of the given double.  If x is NaN, INF, or Zero, this
   prints out the strings NaN, +/- INF, and 0.             */
void print_double_info(std::ostream &os, double x);


#endif  /* _QD_BITS_H */

