/*
  ARPACK++ v1.0 8/1/1997
  c++ interface to ARPACK code.

  MODULE arcomp.h
  arcomplex complex type definition.

  ARPACK Authors
     Richard Lehoucq
     Danny Sorensen
     Chao Yang
     Dept. of Computational & Applied Mathematics
     Rice University
     Houston, Texas
*/

#ifndef ARCOMP_H
#define ARCOMP_H

#include <complex>
using namespace std; 
#if defined(__GNUG__)|| defined(__MWERKS__)
  
#define arcomplex complex

#endif

#if defined(__SUNPRO_CC) || defined(__sgi)

  template <class FLOAT>
  class arcomplex: public complex
  {
   public:

    arcomplex(FLOAT x, FLOAT y): complex(x,y) { }
    arcomplex(): complex() { }
    arcomplex(complex x): complex(x) { }

  };

#endif

#endif // ARCOMP_H



