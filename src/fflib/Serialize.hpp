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
// AUTHORS : Frederic Hecht
// E-MAIL  : frederic.hecht@sorbonne-universite.fr

#ifndef SERIALIZE_HPP_
#define SERIALIZE_HPP_

#include <cstring>
#include "endian.hpp"

struct MPIrank;
class Serialize {
  // we store a refcounter in the pointer p a adress p-sizeof(long)
  // so we can use the copy constructor
  protected:
    size_t lg;
    const char *what;
    char *p;
    public:
    Serialize(size_t lgg, const char *wht)
      : lg(lgg), what(wht), p((new char[lg + sizeof(long)]) + sizeof(long))
    { count() = 0; }

    void resize(size_t lgn) { // Add nov 2010 FH of asyncrone recv MPI ...
      if (lgn > lg) {
        char *p0 = new char[lgn+sizeof(long)];
        memcpy(p0, p-sizeof(long), lg+sizeof(long));
        delete [](p-sizeof(long));
        p = p0+sizeof(long);
      }
      lg = lgn;
    }
    ~Serialize() { if(count()-- == 0) delete [](p-sizeof(long)); }
    size_t size() const { return lg; }

    inline int havebordermesh() {
      size_t pp=2*sizeof(int);
      int bordermesh=0;
      get( pp,bordermesh);
      return bordermesh;
    }

    // mpi routine
    void mpisend(const MPIrank &, long tag, const void *comm);
    Serialize(const MPIrank &, const char *wht, long tag, const void *comm);
    // end mpi routine
    operator void *() { return p; }
    operator char *() { return p; }
    bool samewhat(const char *w) const { return strncmp(what, w, 8) == 0; }

    Serialize(const Serialize & s) :
      lg(s.lg),
      what(s.what),
      p(s.p)
      { count()++; }

    template<typename T> inline void get(size_t &k, T &x) const {
      T xx;//= r_endian(x);
      assert(k <= lg+sizeof(T));
      memcpy(&xx, p + k, sizeof(T));
      k += sizeof(T);
      x = r_endian(xx);
    }

    template<typename T> inline void put(size_t & k, const T & x) {
      if (!(k <= lg+sizeof(T))) {
        cout << " assert put " << k << " <=" << lg + sizeof(T) << endl;
        assert((k <= lg+sizeof(T)));
      }
      T xx = w_endian(x);
      memcpy(p + k, &xx, sizeof(T));
      k += sizeof(T);
    }

  protected:
    long &count() const { return *(long*)(void*)(p-sizeof(long)); }
  private:
    void operator = (Serialize &s); // no affectation
};

#endif // SERIALIZE_HPP_
