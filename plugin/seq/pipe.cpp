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

/* clang-format off */
//ff-c++-LIBRARY-dep:
//ff-c++-cpp-dep:
/* clang-format on */

// to compile ff-c++ pipe.cpp
// warning do not compile under windows...

#include "ff++.hpp"
#include <cstdio>
#include <unistd.h>

#ifdef _WIN32
#include <Windows.h>

long ffsleep(long s) {
  Sleep(s * 1000);
  return 0;
}

long ffusleep(long s) {
  Sleep(s / 1000);
  return 0;
}

#else
long ffsleep(long s) { return sleep(s); }

long ffusleep(long s) { return usleep(s); }

#endif

// #if __clang_major__ > 4
#ifndef _WIN32
#include "pstream.h"

struct pstream {
  redi::basic_pstream< char > *fb;
  // stdiofilebuf * fb;
  ostream *os;
  istream *is;
  pstream(redi::basic_pstream< char > *ff, std::ios_base::openmode mode) : fb(ff), os(0), is(0) {
    if (verbosity > 10) {
      cout << " mode " << mode << endl;
    }

    basic_iostream< char > *bs = fb;    // ->rdbuf();
    ffassert(bs);
    if (mode & ios_base::in) {
      is = dynamic_cast< istream * >(bs);
      assert(is);
    }    // new istream(fb);

    if (mode & ios_base::out) {
      os = dynamic_cast< ostream * >(bs);
      assert(os);
    }    // new ostream(fb);

    if (verbosity > 10) {
      cout << is << " " << os << " ******* " << endl;
    }
  }

  ~pstream( ) {
    if (fb) {
      delete (fb);
    }

    os = 0;
    is = 0;
    fb = 0;
  }

  void flush( ) {
    if (os) {
      os->flush( );
    }
  }
};
static pstream **pstream_init(pstream **const &p, string *const &a, string *const &b) {
  string mode = b ? *b : "w";

  if (mode.length( ) == 0) {
    mode = "wr";
  }

  std::ios_base::openmode om = ios_base::in | ios_base::out;
  if (mode == "r+") {
    om = ios_base::in | ios_base::out;
  } else if (mode == "w") {
    om = ios_base::out;
  } else if (mode == "r") {
    om = ios_base::in;
  } else {
    ExecError("Invalide mode pstream r,r+,w ");
  }

  if (verbosity > 10) {
    *ffapi::cout( ) << "pstream_init: om " << om << "(" << ios_base::in << ios_base::out
                    << ") mode:" << mode << " '" << *a << "'" << endl;
  }

  redi::basic_pstream< char > *pp = new redi::pstream(a->c_str( ), om);
  ffassert(pp);
  *p = new pstream(pp, om);

  if (!*p || !pp) {
    cerr << " Error opening pipe  " << *a << endl;
    ExecError("Error opening pipe");
  }

  return p;
};
static pstream **pstream_init(pstream **const &p, string *const &a) {
  return pstream_init(p, a, 0);
};

#else
// VERSION GNU
#include <ext/stdio_filebuf.h>

typedef __gnu_cxx::stdio_filebuf< char > stdiofilebuf;
struct pstream {
  FILE *f;
  stdiofilebuf *fb;
  ostream *os;
  istream *is;
  pstream(FILE *ff, std::ios_base::openmode mode)
    : f(ff), fb(new stdiofilebuf(f, mode)), os(0), is(0) {
    if (verbosity > 10) {
      cout << " mode " << mode << endl;
    }

    if (mode & ios_base::in) {
      is = new istream(fb);
    }

    if (mode & ios_base::out) {
      os = new ostream(fb);
    }

    if (verbosity > 10) {
      cout << is << " " << os << " ******* " << endl;
    }
  }

  ~pstream( ) {
    if (f) {
      pclose(f);
    }

    if (os) {
      delete os;
    }

    if (is) {
      delete is;
    }

    if (fb) {
      delete (fb);
    }

    f = 0;
    os = 0;
    is = 0;
    fb = 0;
  }

  void flush( ) {
    if (os) {
      os->flush( );
    }

    if (f) {
      fflush(f);
    }
  }
};
// typedef redi::pstream pstream;
// typedef std::string string;
static pstream **pstream_init(pstream **const &p, string *const &a, string *const &b) {
  string mode = b ? *b : "w";

  if (mode.length( ) == 0) {
    mode = "wr";
  }

  std::ios_base::openmode om = ios_base::in | ios_base::out;
  if (mode == "r+") {
    om = ios_base::in | ios_base::out;
  } else if (mode == "w") {
    om = ios_base::out;
  } else if (mode == "r") {
    om = ios_base::in;
  } else {
    ExecError("Invalide mode pstream r,r+,w ");
  }

  if (verbosity > 10) {
    *ffapi::cout( ) << "pstream_init: om " << om << "(" << ios_base::in << ios_base::out
                    << ") mode:" << mode << " '" << *a << "'" << endl;
  }

#ifdef _WIN32
  FILE *pp = _popen(a->c_str( ), mode.c_str( ));
#else
  FILE *pp = popen(a->c_str( ), mode.c_str( ));
#endif
  *p = new pstream(pp, om);

  if (!*p || !pp) {
    cerr << " Error opening pipe  " << *a << endl;
    ExecError("Error opening pipe");
  }

  return p;
};
static pstream **pstream_init(pstream **const &p, string *const &a) {
  return pstream_init(p, a, 0);
};
#endif

AnyType pstream2o(Stack, const AnyType &a) {
  pstream *p = *PGetAny< pstream * >(a);

  ffassert(p->os);
  return SetAny< ostream * >(p->os);
}

AnyType pstream2i(Stack, const AnyType &a) {
  pstream *p = *PGetAny< pstream * >(a);

  ffassert(p->is);
  return SetAny< istream * >(p->is);
}

class istream_good {
 public:
  istream_good(istream *ff) : f(ff) {}

  istream *f;
  operator bool( ) const { return f->good( ); }
};

long cflush(pstream **ppf) {
  pstream &f = **ppf;

  f.flush( );
  return 0;
};

inline istream_good to_istream_good(pstream **f) {
  ffassert((**f).is);
  return istream_good((**f).is);
}

inline bool get_eof(pstream **p) { return (**p).is ? (**p).is->eof( ) : EOF; }

static void inittt( ) {
  Dcl_TypeandPtr< pstream * >(0, 0, ::InitializePtr< pstream * >, ::DeletePtr< pstream * >);
  atype< istream * >( )->AddCast(new E_F1_funcT< istream *, pstream ** >(pstream2i));
  atype< ostream * >( )->AddCast(new E_F1_funcT< ostream *, pstream ** >(pstream2o));
  TheOperators->Add("<-", new OneOperator2_< pstream **, pstream **, string * >(pstream_init));
  TheOperators->Add("<-",
                    new OneOperator3_< pstream **, pstream **, string *, string * >(pstream_init));
  zzzfff->Add("pstream", atype< pstream ** >( ));
  Add< pstream ** >("good", ".", new OneOperator1< istream_good, pstream ** >(to_istream_good));
  Add< pstream ** >("eof", ".", new OneOperator1< bool, pstream ** >(get_eof));
  Global.Add("flush", "(", new OneOperator1< long, pstream ** >(cflush));
  Global.Add("sleep", "(", new OneOperator1< long, long >(ffsleep));
  Global.Add("usleep", "(", new OneOperator1< long, long >(ffusleep));
#ifdef _WIN32
  Global.New("onWIN32", CConstant< bool >(true));
#else
  Global.New("onWIN32", CConstant< bool >(false));
#endif
}

LOADFUNC(inittt);
