// -*- Mode : c++ -*-
//
// SUMMARY  :      
// USAGE    :        
// ORG      : 
// AUTHOR   : Frederic Hecht
// E-MAIL   : hecht@ann.jussieu.fr
//

/*
 
 This file is part of Freefem++
 
 Freefem++ is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation; either version 2.1 of the License, or
 (at your option) any later version.
 
 Freefem++  is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.
 
 You should have received a copy of the GNU Lesser General Public License
 along with Freefem++; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */
#ifndef MESHTYPE_H
#define MESHTYPE_H
#include <limits.h>
namespace bamg {
template<class T> inline T Square (const T &a) { return a*a;} 
template<class T> inline T Min (const T &a,const T &b){return a < b ? a : b;}
template<class T> inline T Max (const T &a,const T & b){return a > b ? a : b;}
template<class T> inline T Abs (const T &a){return a <0 ? -a : a;}
template<class T> inline double Norme (const T &a){return sqrt(a*a);}
template<class T> inline void Exchange (T& a,T& b) {T c=a;a=b;b=c;}
// for pb on microsoft compiler 
template<class T> inline T Max3 (const T &a,const T & b,const T & c){return Max(Max(a,b),c);}
template<class T> inline T Min3 (const T &a,const T & b,const T & c){return Min(Min(a,b),c);}

typedef float  Real4;
typedef double Real8;
typedef short  Int1;
typedef short  Int2;
typedef long   Int4;

#if LONG_BIT > 63
// for alpha and silicon 
 typedef int  Icoor1;  
 typedef long   Icoor2;
 const Icoor1 MaxICoor = 1073741823; // 2^30-1
 const  Icoor2 MaxICoor22 = Icoor2(2)*Icoor2(MaxICoor) * Icoor2(MaxICoor) ;

#elif  defined(BAMG_LONG_LONG)
 typedef long  Icoor1;
 typedef long long   Icoor2;
 const Icoor1 MaxICoor =   1073741823;// 2^30-1
// not a const due to a bug in hp compiler
#define  MaxICoor22 2305843004918726658LL
//const  Icoor2 MaxICoor22 = Icoor2(2)*Icoor2(MaxICoor) * Icoor2(MaxICoor) ;
#else
 typedef int  Icoor1;
 typedef double   Icoor2;
 const Icoor1 MaxICoor = 8388608; // 2^23
 const  Icoor2 MaxICoor22 = Icoor2(2)*Icoor2(MaxICoor) * Icoor2(MaxICoor) ;
#endif
 class Triangles;
extern void MeshError(int Err,Triangles *Th=0) ;
}
#endif
