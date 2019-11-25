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

// *INDENT-OFF* //
// ff-c++-LIBRARY-dep: pthread
// ff-c++-cpp-dep: libff-mmap-semaphore.c
// *INDENT-ON* //

#include "ff++.hpp"
#include "libff-mmap-semaphore.h"

struct ff_pointeur_mmap {
  typedef FF_P_mmap T;
  T *map;
  void init( ) { map = 0; }

  void init(string *nm, long len) {
    map = new T;
    ffmmap_init(map, nm->c_str( ), len);
  }

  void destroy( ) {
    if (map) {
      ffmmap_destroy(map);
    }

    delete map;
    map = 0;
  }
};
struct ff_pointeur_sem {
  typedef FF_P_sem T;
  T *sem;
  void init( ) { sem = 0; }

  void init(string *nm, int crea) {
    sem = new T;
    ffsem_init(sem, nm->c_str( ), crea);
  }

  void destroy( ) {
    if (sem) {
      ffsem_destroy(sem);
    }

    delete sem;
    sem = 0;
  }
};

typedef ff_pointeur_sem FF_p_sem;
typedef ff_pointeur_mmap FF_p_mmap;
typedef ff_pointeur_sem *pFF_p_sem;
typedef ff_pointeur_mmap *pFF_p_mmap;

template< class T >
long Read(pFF_p_mmap const &map, const long &offset, KN< T > *const &data) {
  long n = data->N( ), ln = n * sizeof(T);
  T *pt = *data;

  return ffmmap_read(map->map, (void *)pt, ln, offset);
}

template< class T >
long Write(pFF_p_mmap const &map, const long &offset, KN< T > *const &data) {
  long n = data->N( ), ln = n * sizeof(T);
  T *pt = *data;

  return ffmmap_write(map->map, (void *)pt, ln, offset);
}

template< class T >
long Write(pFF_p_mmap const &map, const long &offset, T *const &pt) {
  return ffmmap_write(map->map, (void *)pt, sizeof(T), offset);
}

template< class T >
long Read(pFF_p_mmap const &map, const long &offset, T *const &pt) {
  return ffmmap_read(map->map, (void *)pt, sizeof(T), offset);
}

long ff_msync(pFF_p_mmap const &map, long const &len, long const &off) {
  return ffmmap_msync(map->map, off, len);
}

long ff_msync(pFF_p_mmap const &map) { return ffmmap_msync(map->map, 0, 0); }

AnyType pmmapinit(Stack, const AnyType &a) {
  pFF_p_mmap p = GetAny< pFF_p_mmap >(a);

  p->init( );
  return p;
}

AnyType pmmadel(Stack, const AnyType &a) {
  pFF_p_mmap p = GetAny< pFF_p_mmap >(a);

  p->destroy( );
  return Nothing;
}

AnyType pseminit(Stack, const AnyType &a) {
  pFF_p_sem p = GetAny< pFF_p_sem >(a);

  p->init( );
  return p;
}

AnyType psemdel(Stack, const AnyType &a) {
  pFF_p_sem p = GetAny< pFF_p_sem >(a);

  p->destroy( );
  return Nothing;
}

long ff_post(pFF_p_sem ps) { return ffsem_post(ps->sem); }

long ff_wait(pFF_p_sem ps) { return ffsem_wait(ps->sem); }

long ff_trywait(pFF_p_sem ps) { return ffsem_trywait(ps->sem); }

pFF_p_sem setpsem(FF_p_sem *const &sem, string *const &nm) {
  sem->init(nm, 0);
  return sem;
}

pFF_p_sem setpsem(FF_p_sem *const &sem, string *const &nm, const bool &creat) {
  sem->init(nm, creat);
  return sem;
}

pFF_p_mmap setpmmap(pFF_p_mmap const &map, string *const &nm, long const &len) {
  map->init(nm, len);
  return map;
}

pFF_p_mmap setpmmap(pFF_p_mmap const &map, string *const &nm) {
  map->init(nm, 0);
  return map;
}

void ff_HandleError(const char *msg, int err) {
  cerr << " Error " << msg << " err= " << err << endl;
  ExecError(msg);
}

static void inittt( ) {
  ff_mmap_sem_verb = verbosity;
  libf_HandleError = ff_HandleError;
  Dcl_TypeandPtr< FF_p_mmap >(0, 0, pmmapinit, pmmadel);
  Dcl_TypeandPtr< FF_p_sem >(0, 0, pseminit, psemdel);

  zzzfff->Add("Pmmap", atype< FF_p_mmap * >( ));
  zzzfff->Add("Psemaphore", atype< FF_p_sem * >( ));
  TheOperators->Add("<-", new OneOperator3_< pFF_p_mmap, pFF_p_mmap, string *, long >(setpmmap));
  TheOperators->Add("<-", new OneOperator2_< pFF_p_mmap, pFF_p_mmap, string * >(setpmmap));

  TheOperators->Add("<-", new OneOperator2_< FF_p_sem *, FF_p_sem *, string * >(setpsem));
  TheOperators->Add("<-", new OneOperator3_< FF_p_sem *, FF_p_sem *, string *, bool >(setpsem));
  Global.Add("Wait", "(", new OneOperator1< long, FF_p_sem * >(ff_wait));
  Global.Add("trywait", "(", new OneOperator1< long, FF_p_sem * >(ff_trywait));
  Global.Add("Post", "(", new OneOperator1< long, FF_p_sem * >(ff_post));

  Global.Add("msync", "(", new OneOperator3_< long, FF_p_mmap *, long, long >(ff_msync));
  Global.Add("msync", "(", new OneOperator1_< long, FF_p_mmap * >(ff_msync));

  Global.Add("Read", "(", new OneOperator3_< long, FF_p_mmap *, long, long * >(Read));
  Global.Add("Write", "(", new OneOperator3_< long, FF_p_mmap *, long, long * >(Write));
  Global.Add("Read", "(", new OneOperator3_< long, FF_p_mmap *, long, KN< long > * >(Read));
  Global.Add("Write", "(", new OneOperator3_< long, FF_p_mmap *, long, KN< long > * >(Write));

  Global.Add("Read", "(", new OneOperator3_< long, FF_p_mmap *, long, double * >(Read));
  Global.Add("Write", "(", new OneOperator3_< long, FF_p_mmap *, long, double * >(Write));
  Global.Add("Read", "(", new OneOperator3_< long, FF_p_mmap *, long, KN< double > * >(Read));
  Global.Add("Write", "(", new OneOperator3_< long, FF_p_mmap *, long, KN< double > * >(Write));
  Global.Add("Read", "(", new OneOperator3_< long, FF_p_mmap *, long, Complex * >(Read));
  Global.Add("Write", "(", new OneOperator3_< long, FF_p_mmap *, long, Complex * >(Write));
  Global.Add("Read", "(", new OneOperator3_< long, FF_p_mmap *, long, KN< Complex > * >(Read));
  Global.Add("Write", "(", new OneOperator3_< long, FF_p_mmap *, long, KN< Complex > * >(Write));
}

LOADFUNC(inittt);
