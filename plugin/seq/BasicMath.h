/****************************************************************************/
/* This file is part of FreeFEM.                                            */
/*                                                                          */
/* FreeFEM is free software: you can redistribute it and/or modify          */
/* it under the terms of the GNU Lesser General Public License as           */
/* published by the Free Software Foundation, either version 3 of           */
/* the License, or (at your option) any later version.                      */
/*                                                                          */
/* FreeFEM is distributed in the hope that it will be useful,               */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of           */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            */
/* GNU Lesser General Public License for more details.                      */
/*                                                                          */
/* You should have received a copy of the GNU Lesser General Public License */
/* along with FreeFEM. If not, see <http://www.gnu.org/licenses/>.          */
/****************************************************************************/
// SUMMARY : A few common functions which are absent in C++
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
// AUTHORS : Jean-Marie Mirebeau
// E-MAIL  :
// DATE    : 11/07/11

#ifndef BASIC_MATH
#define BASIC_MATH
#include <fstream>
#include <sstream>
#include "assert.h"
using namespace std;

template< class ring >
inline ring square(ring u) {
  return u * u;
}

template< class ring >
inline ring cube(ring u) {
  return u * u * u;
}

inline double sqrt3(double r) { return r >= 0 ? pow(r, 1 / 3.) : -pow(-r, 1 / 3.); }

template< class ordered_zero >
inline int sign(ordered_zero a) {
  return a > 0 ? 1 : (a == ordered_zero(0) ? 0 : -1);
}

inline int mod(int a, int N) {
  const double u = a % N;
  return u >= 0 ? u : u + N;
}    // math modulo. Beware: % is wrong for negative numbers

inline int smallmod(int a, int N) {
  const int NHalf = N / 2;
  return -NHalf + mod(a + NHalf, N);
}    // math symmetrized modulo, in [-N/2,N/2] or so

// inline int round0(double a){return a>0 ? ceil(a-0.5) : floor(a+0.5);} //rounding towards the
// closest integer to 0.

/*!
 * Some debug macros
 * \brief try_debug helps to locate errors in combination with external libraries
 */
#if DEBUG
#define assert_msg(condition, message) \
  if (!(condition)) {                  \
    std::cerr << message << " : ";     \
    assert(condition);                 \
    assert(false);                     \
  }
#define try_debug_msg(instructions, message) \
  try {                                      \
    instructions;                            \
  } catch (...) {                            \
    assert_msg(false, message);              \
  }
#define try_debug(instructions) \
  try {                         \
    instructions;               \
  } catch (...) {               \
    assert(false);              \
  }
#else
#define assert_msg(condition, message) assert(condition)
#define try_debug_msg(instructions, message) \
  { instructions; }
#define try_debug(instructions) \
  { instructions; }
#endif

/*!
 * Display of mathematical structures in diverse formats
 * \biref typical usage: some_ostream << some_Format_Math << data1 << data2 <<data3 << endl;
 */

enum Format_Math { Standard, Mathematica };    // more formats to come if required

class ostream_math {    // This class only contains a reference to an ostream, and a math format
                        // type
 public:
  Format_Math format;
  ostream &os;
  ostream_math(ostream &OS, Format_Math Format) : os(OS), format(Format) {}
};

inline ostream_math operator<<(ostream &f, Format_Math Format) { return ostream_math(f, Format); }

extern ostream_math coutMath;
// ostream_math coutMath = cout << Mathematica;
// default bahavior: ostream
inline ostream_math operator<<(ostream_math f, ostream &(*func)(ostream &)) {
  func(f.os);
  return f;
}

template< class E >
inline ostream_math operator<<(ostream_math f, const E &n) {
  f.os << n;
  return f;
}

// Display of doubles: 1.5e+34 -> 1.5*10^+34
inline ostream_math operator<<(ostream_math f, double x) {
  if (f.format == Mathematica) {
    ostringstream oss;
    oss << x;
    string pxstr = oss.str( );
    const char *px = pxstr.c_str( );
    if (px[0] == 'N') {
      f << "Indeterminate";
    } else if (px[0] == 'i') {
      f << "Infinity";
    } else if (px[0] == '-' && px[1] == 'i') {
      f << "-Infinity";
    } else {
      for (int i = 0; i < 20 && px[i] > 0; i++) {
        if (px[i] == 'e') {
          char Buffer[20];

          for (int j = 0; j < i; j++) {
            Buffer[j] = px[j];
          }

          Buffer[i] = 0;
          f << Buffer << "*10^" << px + i + 1;
          return f;
        }
      }

      f << px;
    }    // if px[0]
  } else {
    f.os << x;
  }

  return f;
}

template< class ForwardIterator >
void print_array(ostream_math f, ForwardIterator first, ForwardIterator last,
                 bool one_per_line = false) {
  string sep = one_per_line ? ",\n" : ",";

  f << "{";
  if (first != last) {
    f << *first++;
  }

  while (first != last) {
    f << sep << *first++;
  }

  f << "}";
}

/*
 * //Display of arrays
 * template<class E> void print_array(ostream & f, const E * tab, int N, bool one_per_line=false){
 *  if(one_per_line) for(int i=0; i<N; i++) f << tab[i] << endl;
 *  else for(int i=0; i<N; i++) f << tab[i] << " ";
 * }
 *
 * template<class E> void print_array(ostream_math f, const E * tab, int N, bool
 * one_per_line=false){ if(f.format==Mathematica) { if(N==0) {f << "{}"; return;} f << "{"; for(int
 * i=0; i<N; i++) {f << tab[i]; if(i<N-1) f << ",";} f << "}";} else {print_array(f.os, tab, N,
 * one_per_line); return;}
 * }
 *
 * //Sampling function values. Not the best way to do it, but how to avoid it without copy-pasting
 * everything above ? template<class E> void print_array(ostream_math f, E (*func)(int), int N, bool
 * one_per_line=false){ E * tab = new E [N]; for(int i=0; i<N; i++) tab[i]=func(i); print_array(f,
 * tab, N, one_per_line); delete tab;
 * }
 *
 * template<class E> inline void print_array(ostream & f, E (*func)(int), int N, bool
 * one_per_line=false){ print_array(f << Standard, func, N, one_per_line);}
 */

// conversion to string
template< class E >
std::string to_string(const E &e) {
  ostringstream oss;
  oss << e;
  return oss.str( );
}

#endif
