/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE debug.h.
   Interface to ARPACK FORTRAN debugging facilities.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef DEBUG_H
#define DEBUG_H

#include "arch.h"
#include "arpackf.h"


inline void TraceOff()

/*
  This function sets all ARPACK FORTRAN debug variables to zero.
*/

{

  F77NAME(debug).logfil = 6;
  F77NAME(debug).ndigit = 0;
  F77NAME(debug).mgetv0 = 0;
  F77NAME(debug).msaupd = 0;
  F77NAME(debug).msaup2 = 0;
  F77NAME(debug).msaitr = 0;
  F77NAME(debug).mseigt = 0;
  F77NAME(debug).msapps = 0;
  F77NAME(debug).msgets = 0;
  F77NAME(debug).mseupd = 0;
  F77NAME(debug).mnaupd = 0;
  F77NAME(debug).mnaup2 = 0;
  F77NAME(debug).mnaitr = 0;
  F77NAME(debug).mneigt = 0;
  F77NAME(debug).mnapps = 0;
  F77NAME(debug).mngets = 0;
  F77NAME(debug).mneupd = 0;
  F77NAME(debug).mcaupd = 0;
  F77NAME(debug).mcaup2 = 0;
  F77NAME(debug).mcaitr = 0;
  F77NAME(debug).mceigt = 0;
  F77NAME(debug).mcapps = 0;
  F77NAME(debug).mcgets = 0;
  F77NAME(debug).mceupd = 0;

} // TraceOff.


inline void sTraceOn(const int digit, const int getv0, const int aupd, 
                     const int aup2, const int aitr, const int eigt,
                     const int apps, const int gets, const int eupd)

/*
  This function sets the values of all ARPACK FORTRAN debug 
  variables corresponding to real symmetric eigenvalue problems.
*/

{

  F77NAME(debug).logfil = 6;
  F77NAME(debug).ndigit = digit;
  F77NAME(debug).mgetv0 = getv0;
  F77NAME(debug).msaupd = aupd;
  F77NAME(debug).msaup2 = aup2;
  F77NAME(debug).msaitr = aitr;
  F77NAME(debug).mseigt = eigt;
  F77NAME(debug).msapps = apps;
  F77NAME(debug).msgets = gets;
  F77NAME(debug).mseupd = eupd;

} // sTraceOn.


inline void nTraceOn(const int digit, const int getv0, const int aupd, 
                     const int aup2, const int aitr, const int eigt,
                     const int apps, const int gets, const int eupd)

/*
  This function sets the values of all ARPACK FORTRAN debug 
  variables corresponding to real nonsymmetric eigenvalue problems.
*/

{

  F77NAME(debug).logfil = 6;
  F77NAME(debug).ndigit = digit;
  F77NAME(debug).mgetv0 = getv0;
  F77NAME(debug).mnaupd = aupd;
  F77NAME(debug).mnaup2 = aup2;
  F77NAME(debug).mnaitr = aitr;
  F77NAME(debug).mneigt = eigt;
  F77NAME(debug).mnapps = apps;
  F77NAME(debug).mngets = gets;
  F77NAME(debug).mneupd = eupd;

} // nTraceOn.


inline void cTraceOn(const int digit, const int getv0, const int aupd, 
                     const int aup2, const int aitr, const int eigt,
                     const int apps, const int gets, const int eupd)

/*
  This function sets the values of all ARPACK FORTRAN debug 
  variables corresponding to complex eigenvalue problems.
*/

{

  F77NAME(debug).logfil = 6;
  F77NAME(debug).ndigit = digit;
  F77NAME(debug).mgetv0 = getv0;
  F77NAME(debug).mcaupd = aupd;
  F77NAME(debug).mcaup2 = aup2;
  F77NAME(debug).mcaitr = aitr;
  F77NAME(debug).mceigt = eigt;
  F77NAME(debug).mcapps = apps;
  F77NAME(debug).mcgets = gets;
  F77NAME(debug).mceupd = eupd;

} // cTraceOn.


#endif // DEBUG_H
