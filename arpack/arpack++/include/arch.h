/*
  ARPACK++ v1.0 8/1/1997
  c++ interface to ARPACK code.

  MODULE arch.h
  Modified version of arch.h (from LAPACK++ v. 1.1).
  Machine dependent functions and variable types.

  ARPACK Authors
     Richard Lehoucq
     Danny Sorensen
     Chao Yang
     Dept. of Computational & Applied Mathematics
     Rice University
     Houston, Texas
*/


#ifndef ARCH_H
#define ARCH_H

// ARPACK++ arcomplex type definition.
// If you are not using g++ (or CC) and also are not intending 
// use complex variables, comment out the following line.

#include "arcomp.h"

// STL vector class.
// If the Standard Template Library is not available at your system
// and you do not want to install it, comment out the following line.

//#include <vector.h>

// UMFPACK parameters.
// These parameters are used by UMFPACK library functions. Normally
// they are not modified by the user. To use the default value, set
// the parameter to zero. For a complete description of all UMFPACK 
// parameters, see the library documentation. 

#define UICNTL7 0 // icntl(7). Block size for the blas (machine-dependent).
#define UICNTL5 0 // icntl(5). Number of columns to examine during pivot search.
#define UCNTL2  0 // cntl(2). Amalgamation parameter.
#define UKEEP7  0 // keep(7). Absolute number of elements a column must have
                  // to be considered "dense".
#define UKEEP8  0 // keep(8). Relative number of elements a column must have
                  // to be considered "dense". Dense columns have more
                  // than max{0,UMFABDEN,UMFREDEN*sqrt(n)} elements.

// Linkage names between C, C++, and Fortran (platform dependent)

#if  defined(RIOS) && !defined(CLAPACK)
#define F77NAME(x) x
#else
//#include <generic.h> 
#define F77NAME(x) x ## _
//name2(x,_)
#endif

#if defined(SGI) && !defined(SGI_DEC)
#define SGI_DEC

extern "C" {
	void mkidxname() {}
	void mkdatname() {}
}
#endif


// Type conversion.

typedef int integer;
typedef int logical;

#ifdef __SUNPRO_CC

  typedef int bool;
  int true  = 1;
  int false = 0;

#endif


#endif // ARCH_H
