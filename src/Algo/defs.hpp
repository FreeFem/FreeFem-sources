#ifndef DEFS_H
#define DEFS_H

#ifndef GLIB
#define GLIB

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <stddef.h>
//#include <cassert>
#include <cmath>
#include <string.h>
#include <list>
using namespace std;

//#include "mvvtp.h"
//#include "mvblas.h"

#endif

// DEFINE simple constants
#ifndef Pi
#define Pi (3.141592653589793)
#endif
#ifndef TINY // sert dans BrentLS et Matrix
#define TINY 		1.0e-20
#endif
#ifndef True
#define True (1)
#endif
#ifndef False
#define False (0)
#endif
#ifndef Yes
#define Yes (1)
#endif
#ifndef No
#define No (0)
#endif

//Define some simple template functions
#ifndef Sgn
#define Sgn(x) ((x) < 0 ? -1.0 : 1.0)
#endif
#ifndef Abs
#define Abs(x) ((x) < 0 ? -(x) : (x))
#endif
#ifndef Max
#define Max(x,y) ((x) > (y) ? (x) : (y))
#endif
#ifndef Min
#define Min(x,y) ((x) < (y) ? (x) : (y))
#endif

#ifndef inValidSize
#define inValidSize() cerr<<"Warning: incompatible sizes!"<<endl;
#endif

// Opérations sur les listes
//
template<class Type>
void affiche(list<Type> l)
{
  typename list<Type>::iterator il;
  for (il=l.begin();il!=l.end();il++) cout <<(*il)<<" ";
  cout<<endl;
}

template<class Type>
list<Type> normalize(list<Type> l) {
  list<Type> v(l);
  Type scale = l.front();
  typename list<Type>::iterator il;
	
  if (scale == 0) {
	cerr << "Zero first element, cannot be normed! \n";
	return v;
  }
  
  else {
	for (il=v.begin();il!=v.end();il++)
	  *il /= scale;
	return v;
  }
}

#endif


