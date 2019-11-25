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
// SUMMARY : ...
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
// AUTHORS : ...
// E-MAIL  : ...

// Example C++ function "add random function ", dynamically loaded into "load.edp"

// *INDENT-OFF* //
// ff-c++-LIBRARY-dep:
// ff-c++-cpp-dep:
// *INDENT-ON* //

// TODO: remove this block as soon as autoconf is removed from FreeFem++
#ifndef CMAKE
#include "config.h"
#endif

#include <iostream>
#include <cfloat>
using namespace std;
#include "error.hpp"
#include "AFunction.hpp"

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <ctime>

// FFCS: random() and srandom() do not seem to be available in MinGW
#ifdef _WIN32
#define random rand
#define srandom srand
#endif

unsigned long good_seed( ) {
  unsigned long random_seed, random_seed_a, random_seed_b;
  std::ifstream file("/dev/random", std::ios::binary);

  if (file.is_open( )) {
    unsigned long memblock[10] = {};
    size_t size = sizeof(int);
    file.read((char *)(void *)memblock, size);
    file.close( );
    random_seed_a = memblock[0];
  }    // end if
  else {
    random_seed_a = 0;
  }

  random_seed_b = std::time(0);
  random_seed = random_seed_a xor random_seed_b;
  if (verbosity > 1) {
    cout << " good_seed =" << random_seed << endl;
  }

  return random_seed;
}    // end good_seed()

template< class R >
class OneOperator_0 : public OneOperator {
  class E_F0_F : public E_F0mps {
   public:
    typedef R (*func)( );
    func f;
    E_F0_F(func ff) : f(ff) {}

    AnyType operator( )(Stack) const { return SetAny< R >(f( )); }

    operator aType( ) const { return atype< R >( ); }
  };

  typedef R (*func)( );
  func f;

 public:
  E_F0 *code(const basicAC_F0 &) const { return new E_F0_F(f); }

  OneOperator_0(func ff) : OneOperator(map_type[typeid(R).name( )]), f(ff) {}
};

// FFCS: rand() wants an int on Win32 and Win64, so we need to have a new function "long ffrandom()"
// to create the OneOperator object which requires longs (and the MinGW compiler (g++ 4.5) refuses
// to compile if OneOperator and the underlying function do not have the same type). FFCS patched
// this source with "long ffrandom(){return random();}" but FF now comes with "genrand_int31()".

#ifdef _WIN32

void init_by_array(unsigned long init_key[], int key_length);
long genrand_int31(void);
void init_genrand(unsigned long);
// hach for window ... F. HEcht
long ffsrandom(long s) {
  init_genrand((unsigned int)s);
  return 0;
}

long ffrandom( ) { return genrand_int31( ); }

long ffsrandomdev( ) {
  init_genrand(good_seed( ));
  return 0;
}

#else

long ffsrandom(long s) {
  srandom((unsigned int)s);
  return 0;
}

long ffrandom( ) { return random( ); }

long ffsrandomdev( ) {
#ifdef HAVE_SRANDOMDEV
  srandomdev( );
#else
  srandom(good_seed( ));
#endif

  return 0;
}

#endif

static void init( ) {
  Global.Add("srandomdev", "(", new OneOperator_0< long >(ffsrandomdev));
  Global.Add("srandom", "(", new OneOperator1< long >(ffsrandom));
  Global.Add("random", "(", new OneOperator_0< long >(ffrandom));
}

LOADFUNC(init);
/* These real versions are due to Isaku Wada, 2002/01/09 added */
