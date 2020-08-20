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

#include <cstdint>
#include <fstream>
#include <iostream>
#include <cfloat>
#include <cmath>
#include <cstring>
#include <complex>
#include <cstdlib>

using namespace std;
#include "error.hpp"
#include "AFunction.hpp"
// FH: I have move AFunction_ext.hpp to fflib dir.
#include "AFunction_ext.hpp"
// Add J. Morice for AFunction_ext.hpp
#include "ufunction.hpp"
using namespace std;
#include "rgraph.hpp"
#include "RNM.hpp"
// after RNM otherwise
// trouble with index in RNM (I do no understande FH)
#include <set>
#include <vector>
#include <map>

#include "fem.hpp"

#include "FESpacen.hpp"
#include "FESpace.hpp"
#include "HashMatrix.hpp"
#include "SparseLinearSolver.hpp"

#include "MeshPoint.hpp"
#include "Mesh2dn.hpp"
#include "Mesh3dn.hpp"
#include "MeshSn.hpp"
#include "Operator.hpp"
#include "lex.hpp"
#include "libmesh5.h"
#include "lgfem.hpp"
#include "lgmesh3.hpp"
#include "lgsolver.hpp"
#include "problem.hpp"

//FFCS redirection
#include "../fflib/ffapi.hpp"
void ff_atend( void (*atendff)());

#undef MPICH_SKIP_MPICXX
#define MPICH_SKIP_MPICXX
#undef MPICH_IGNORE_CXX_SEEK
#define MPICH_IGNORE_CXX_SEEK
#include <mpi.h>

// Remark: on mipich MPI_Comm, MPI_Resquest, MPI_Group, MPI_Op are int
// => encapsulation

//static long verbosity = 1000;
template<class MPI_type, int DIFF>
struct fMPI {
  MPI_type v;
  operator MPI_type &() { return v; }
  operator MPI_type *() { return &v; }
  operator MPI_type () const { return v; }

  // MPI_type * operator &() { return &v; }
  void operator = (const MPI_type vv) { v=vv; }
  fMPI(const MPI_type vv=0) : v(vv) {}
  bool operator != (MPI_type vv) const { return vv != v; }
  bool operator == (MPI_type vv) const { return vv == v; }
};

// the encapsulation for the for MPI type (int on mpich )
typedef fMPI<MPI_Comm, 1> fMPI_Comm;
typedef fMPI<MPI_Group, 2> fMPI_Group;
typedef fMPI<MPI_Request, 3> fMPI_Request;
typedef fMPI<MPI_Op, 4> fMPI_Op;
// end of encapsulation ..

// to send a sparse matrix we send header, line array, colmun array, value array.
// end afer the fist resquest we need to do allocation.
// so in this cas the communication is done in
// 2 step
//    1 the header,
//       alloc time
//    2 the message
// a couple request, pointer.
// Not use of IPROBE because probelem of wait.
class MPIrank;
class DoOnWaitMPI_Request;

map<MPI_Request*, DoOnWaitMPI_Request *> ToDoOnWaitMPI_Request;

void GetPendingWait() ;

template<class T> struct MPI_TYPE { static MPI_Datatype TYPE(){ return MPI_BYTE; } };
template<> struct MPI_TYPE<long> { static MPI_Datatype TYPE(){ return MPI_LONG; } };
template<> struct MPI_TYPE<int> { static MPI_Datatype TYPE(){ return MPI_INT; } };
template<> struct MPI_TYPE<double> { static MPI_Datatype TYPE(){ return MPI_DOUBLE; } };
template<> struct MPI_TYPE<char> { static MPI_Datatype TYPE(){ return MPI_BYTE; } };

#ifdef HAVE_MPI_DOUBLE_COMPLEX
template<> struct MPI_TYPE<Complex> { static MPI_Datatype TYPE(){ return MPI_DOUBLE_COMPLEX; } };
#endif
template<class T> struct MPI_WHAT {};
template<> struct MPI_WHAT<long> { static const int WHAT=101; };
template<> struct MPI_WHAT<double> { static const int WHAT=102; };
template<> struct MPI_WHAT<Complex> { static const int WHAT=103; };
template<> struct MPI_WHAT<KN<long>* > { static const int WHAT=104; };
template<> struct MPI_WHAT<KN<double>* > { static const int WHAT=105; };
template<> struct MPI_WHAT<KN<Complex>* > { static const int WHAT=106; };
template<> struct MPI_WHAT<KNM<long>* > { static const int WHAT=107; };
template<> struct MPI_WHAT<KNM<double>* > { static const int WHAT=108; };
template<> struct MPI_WHAT<KNM<Complex>* > { static const int WHAT=109; };

template<class T> struct MPI_TAG {};
template<> struct MPI_TAG<long> { static const int TAG=5; };
template<> struct MPI_TAG<double> { static const int TAG=4; };
template<> struct MPI_TAG<Complex > { static const int TAG=6; };
template<> struct MPI_TAG<KN<long>* > { static const int TAG=11; };
template<> struct MPI_TAG<KN<double>* > { static const int TAG=12; };
template<> struct MPI_TAG<KN<Complex>* > { static const int TAG=13; };
template<> struct MPI_TAG<KNM<long>* > { static const int TAG=14; };
template<> struct MPI_TAG<KNM<double>* > { static const int TAG=15; };
template<> struct MPI_TAG<KNM<Complex>* > { static const int TAG=16; };
template<> struct MPI_TAG<Mesh* > { static const int TAG=1000; };
template<> struct MPI_TAG<Mesh3* > { static const int TAG=1010; };
template<> struct MPI_TAG<MeshS* > { static const int TAG=1040; };
template<> struct MPI_TAG<MeshL* > { static const int TAG=1050; };
template<> struct MPI_TAG<const Mesh* > { static const int TAG=1000; };
template<> struct MPI_TAG<const Mesh3* > { static const int TAG=1010; };
template<> struct MPI_TAG<const MeshS* > { static const int TAG=1040; };
template<> struct MPI_TAG<const MeshL* > { static const int TAG=1050; };

template<> struct MPI_TAG<Matrice_Creuse<double> *> { static const int TAG=1020; };
template<> struct MPI_TAG<Matrice_Creuse<Complex> *> { static const int TAG=1030; };

void f_initparallele(int &, char **&);
void f_init_lgparallele();

extern long mpirank ;
extern long mpisize ;

// for syncro communication
MPI_Request *Syncro_block = reinterpret_cast<MPI_Request* > (1);

const size_t sizempibuf = 1024*320;

template<class R>
long WSend (R *v, int l, int who, int tag, MPI_Comm comm, MPI_Request *rq) {
  long ret = 0;
  MPI_Request rq0, *request = &rq0;
  if(verbosity>100)
    cout << mpirank << " send to " << who << " tag " << tag << " " << rq << " " << comm << " syncro " << (rq == Syncro_block) << endl;
  if (rq == Syncro_block || rq == 0)
    ret = MPI_Send((void *)v, l, MPI_TYPE<R>::TYPE(), who, tag, comm);
  else {
    ret = MPI_Isend((void *)v, l, MPI_TYPE<R>::TYPE(), who, tag, comm, request);
    *rq = *request;
  }
  return ret;
}

template<class T>
void CheckContigueKNM (const KNM_<T> &t) {
  if (t.step != 1 && !t.IsVector1()) {
    cout << " step = "<< t.step << " size " << t.N() << " " << &t[0] << " " << &t[1] << endl;
    ExecError("Sorry the array is not contiguous (step != 1) ");
  }
}

template<class T>
void CheckContigueKN (const KN_<T> &t) {
  if (t.step != 1 && t.N() > 1) {
    cout<< " step = " << t.step << " size " << t.N() << " " << &t[0] << " " << &t[1] << endl;
    ExecError("Sorry the array is not contiguous (step != 1) ");
  }
}

template<>
long WSend<Complex> (Complex *v, int n, int who, int tag, MPI_Comm comm, MPI_Request *rq) {
  long ret = 0;
  MPI_Request rq0, *request = &rq0;
  if (verbosity > 100)
    cout << mpirank << " send to " << who << " tag " << tag << " " << rq << " " << comm << " syncro "<< (rq == Syncro_block) << endl;
  if (rq == Syncro_block || rq == 0) {
#ifdef HAVE_MPI_DOUBLE_COMPLEX
    ret = MPI_Send(reinterpret_cast<void*>(v), n, MPI_DOUBLE_COMPLEX, who, tag, comm);
#else
    n *= 2;
    ret = MPI_Send(reinterpret_cast<void*>(v), n, MPI_DOUBLE, who, tag, comm);
#endif
  } else {
#ifdef HAVE_MPI_DOUBLE_COMPLEX
    ret = MPI_Isend(reinterpret_cast<void*>(v), n, MPI_DOUBLE_COMPLEX, who, tag, comm, request);
#else
    n *= 2;
    ret = MPI_Isend(reinterpret_cast<void*>(v), n, MPI_DOUBLE, who, tag, comm, request);
    n /= 2;
#endif
    if(rq) *rq = *request;
    else MPI_Request_free(request);
  }
  return ret;
}

template<class R>
long WRecv (R *v, int n, int who, int tag, MPI_Comm comm, MPI_Request *rq) {
  MPI_Status status;
  if (rq && (rq != Syncro_block))
    return MPI_Irecv(reinterpret_cast<void*>(v), n, MPI_TYPE<R>::TYPE(), who, tag, comm, rq);
  else
    return MPI_Recv(reinterpret_cast<void*>(v), n, MPI_TYPE<R>::TYPE(), who, tag, comm, &status);
}

template<>
long WRecv<Complex> (Complex *v, int n, int who, int tag, MPI_Comm comm, MPI_Request *rq) {
  MPI_Status status;
#ifdef HAVE_MPI_DOUBLE_COMPLEX
  if (rq && (rq != Syncro_block))
    return MPI_Irecv(reinterpret_cast<void*>(v), n, MPI_DOUBLE_COMPLEX, who, tag, comm, rq);
  else
  return MPI_Recv(reinterpret_cast<void*>(v), n, MPI_DOUBLE_COMPLEX, who, tag, comm, &status);
#else
  n *= 2;
  if (rq && (rq != Syncro_block))
    return MPI_Irecv(reinterpret_cast<void*>(v), n, MPI_DOUBLE, who, tag, comm, rq);
  else
    return MPI_Recv(reinterpret_cast<void*>(v), n, MPI_DOUBLE, who, tag, comm, &status);
#endif
}

template<class R>
void WBcast (R *v, int n, int who, MPI_Comm comm) {
  assert(v && n > 0);
  MPI_Bcast(reinterpret_cast<void*>(v), n, MPI_TYPE<R>::TYPE(), who, comm);
}

template<>
void WBcast<Complex> (Complex *v, int n, int who, MPI_Comm comm) {
  assert(v && n > 0);
#ifdef HAVE_MPI_DOUBLE_COMPLEX
  MPI_Bcast(reinterpret_cast<void*>(v), n, MPI_DOUBLE_COMPLEX /*MPI_TYPE<R>::TYPE()*/, who, comm);
#else
  n *= 2;
  MPI_Bcast(reinterpret_cast<void*>(v), n, MPI_DOUBLE, who, comm);
#endif
}

struct MPIrank {
  int who,lmpirank;
  MPI_Comm comm;
  MPI_Request *rq;
  // mutable bool block;

  MPIrank(int i=0, MPI_Comm com=MPI_COMM_WORLD, MPI_Request *rqq=0)
    : who(i), comm(com), rq(rqq) {
    int n;
    MPI_Comm_size(comm, &n);
    MPI_Comm_rank(comm, &lmpirank);
  }

  long Send(double a) const { return WSend(&a, 1, who, MPI_TAG< double >::TAG, comm, rq); }
  long Send(long a) const { return WSend(&a, 1, who, MPI_TAG< long >::TAG, comm, rq); }
  long Send(Complex a) const { return WSend(&a, 1, who, MPI_TAG< Complex >::TAG, comm, rq); }

  long Recv(double & a) const { return WRecv(&a, 1, who, MPI_TAG< double >::TAG ,comm, rq); }
  long Recv(long & a) const { return WRecv(&a, 1, who, MPI_TAG< long >::TAG ,comm, rq); }
  long Recv(Complex & a) const { return WRecv(&a, 1, who, MPI_TAG< Complex >::TAG ,comm, rq); }

  const MPIrank & Bcast(double & a) const { WBcast(&a, 1, who, comm); return *this; }
  const MPIrank & Bcast(long & a) const { WBcast(&a, 1, who, comm); return *this; }
  const MPIrank & Bcast(Complex & a) const { WBcast(&a, 1, who, comm); return *this; }

  template<class R>
  long Recv (KN<R> & a) const {
    assert(&a);
    CheckContigueKN(a);

    if (verbosity > 99)
      cout << " ---- " << who << " >> " << &a << " " << a.N() << " " << a.step << " " << MPI_TAG<KN<R>* >::TAG
           << " from " << mpirank << " " << (R *)a << endl;
    int n = a.N();
    long ll = WRecv((R *)a, n, who, MPI_TAG<KN<R>* >::TAG, comm, rq);
    if (verbosity > 99)
      cout << " ++++ " << who << " >> " << &a << " " << a.N() << " " << MPI_TAG<KN<R>* >::TAG
           << " from " << mpirank << " " << (R *)a << endl;
    ffassert(a.N() == n);
    return ll;
  }

  template<class R>
  long Send (const KN<R> *aa) const {
    const KN<R> & a = *aa;
    ffassert(&a);
    int n = a.N();
    CheckContigueKN(*aa);
    if (verbosity > 99)
      cout << " .... " << who << " >> " << &a << " " << a.N() << " " << a.step << " " << MPI_TAG<KN<R>* >::TAG
           << " from " << mpirank << " " << (R *)a << endl;
    return WSend((R *)a, n, who, MPI_TAG<KN<R>* >::TAG, comm, rq);
  }

  template<class R>
  const MPIrank & Bcast (const KN<R> &a) const {
    assert(&a);
    int n = a.N();
    CheckContigueKN(a);

    WBcast((R *)a, n, who, comm);
    ffassert(a.N() == n);
    return *this;
  }

  // KNM *********************************** Add FH. Nov 2016
  template<class R>
  long Recv (KNM<R> & a) const {
    assert(&a);
    CheckContigueKNM(a);

    if (verbosity > 99)
      cout << " ---- " << who << " >> " << &a << " " << a.N() << "x" << a.M() << " " << MPI_TAG<KNM<R>* >::TAG
           << " from " << mpirank << " " << (R *)a << endl;
    int n = a.N()*a.M();
    long ll = WRecv((R *)a, n, who, MPI_TAG<KNM<R>* >::TAG, comm, rq);
    if (verbosity > 99)
      cout << " ++++ " << who << " >> " << &a << " " << a.N() << "x" << a.M() << " " << MPI_TAG<KNM<R>* >::TAG
           << " from " << mpirank << " " << (R *)a << endl;
    ffassert(a.N()*a.M() == n);
    return ll;
  }

  template<class R>
  long Send (const KNM<R> *aa) const {
    const KNM<R> & a = *aa;
    ffassert(&a);
    int n = a.N()*a.M();
    CheckContigueKNM(*aa);
    if (verbosity > 99)
      cout << " .... " << who << " >> " << &a << " " << a.N() << "x" << a.M() << " " << a.step << " " << MPI_TAG<KNM<R>* >::TAG
           << " from " << mpirank << " " << (R *)a << endl;
    return WSend((R *)a, n, who, MPI_TAG<KNM<R>* >::TAG, comm, rq);
  }

  template<class R>
  const MPIrank & Bcast (const KNM<R> &a) const {
    assert(&a);
    int n = a.N()*a.M();
    CheckContigueKNM(a);

    WBcast((R *)a, n, who, comm);
    ffassert(a.N()*a.M() == n);
    return *this;
  }

  // KNM ***********************************
  const MPIrank & Bcast (Fem2D::Mesh const *&a) const {
    if(verbosity>100)
      cout << " MPI Bcast (mesh *) " << a << endl;
    Serialize *buf = 0;
    long nbsize = 0;
    if (who == lmpirank) {
      buf = new Serialize((*a).serialize());
      nbsize = buf->size();
    }
    WBcast(&nbsize, 1, who, comm);
    if (who != lmpirank)
      buf = new Serialize(nbsize, Fem2D::Mesh::magicmesh);
    assert(nbsize);
    if (verbosity > 200)
      cout << " size to bcast : " << nbsize << " mpirank : " << mpirank << endl;

    WBcast((char *)(*buf), nbsize, who, comm);

    if (who != lmpirank) {
      if (a) (*a).decrement();
      Fem2D::Mesh *pTh = new Fem2D::Mesh(*buf);
      Fem2D::R2 Pn, Px;
      pTh->BoundingBox(Pn, Px);
      pTh->quadtree = new Fem2D::FQuadTree(pTh, Pn, Px, pTh->nv);
      a = pTh;
    }
    delete buf;
    return *this;
  }

  const MPIrank & Bcast (Fem2D::Mesh3 const *&a) const {
    if(verbosity>100)
      cout << " MPI Bcast (const mesh3 *) " << a << endl;
    Serialize *buf = 0;
    long nbsize[2]={0,0};
    if (who == lmpirank) {
        if((*a).meshS) {
            buf = new Serialize((*a).serialize_withBorderMesh());
            nbsize[1]=1;
        }
        else
            buf = new Serialize((*a).serialize());
      nbsize[0] = buf->size();
    }
    WBcast(&nbsize[0], 2, who, comm);
    if (who != lmpirank)
      buf = new Serialize(nbsize[0], Fem2D::GenericMesh_magicmesh);
    assert(nbsize);
    if (verbosity > 200)
      cout << " size to bcast : " << nbsize[0] << " mpirank : " << mpirank << endl;

    WBcast((char *)(*buf), nbsize[0], who, comm);
    Fem2D::Mesh3 * aa;
    if (who != lmpirank) {
      if (a) (*a).decrement();
        if(nbsize[1])
            aa = new Fem2D::Mesh3(*buf, 1);
        else
            aa = new Fem2D::Mesh3(*buf);
        
      aa->BuildGTree();
      a = aa;
    }
    delete buf;
    return *this;
  }

  const MPIrank & Bcast (Fem2D::MeshS const *&a) const {
    if(verbosity>100)
      cout << " MPI Bcast (const meshS *) " << a << endl;
    Serialize *buf = 0;
    long nbsize[2]={0,0};
    if (who == lmpirank) {
        if((*a).meshL) {
            buf = new Serialize((*a).serialize_withBorderMesh());
            nbsize[1]=1;
        }
        else
            buf = new Serialize((*a).serialize());
      nbsize[0] = buf->size();
    }
    WBcast(&nbsize[0], 2, who, comm);
    if (who != lmpirank)
      buf = new Serialize(nbsize[0], Fem2D::GenericMesh_magicmesh);
    assert(nbsize);
    if (verbosity > 200)
      cout << " size to bcast : " << nbsize[0] << " mpirank : " << mpirank << endl;
     WBcast((char *)(*buf), nbsize[0], who, comm);
    Fem2D::MeshS * aa;
    if (who != lmpirank) {
      if (a) (*a).decrement();
        if(nbsize[1])
            aa = new Fem2D::MeshS(*buf, 1);
        else
            aa = new Fem2D::MeshS(*buf);
        
      aa->BuildGTree();
      a = aa;
    }
    delete buf;
    return *this;
  }
    
  const MPIrank & Bcast (Fem2D::MeshL const *&a) const {
      if(verbosity>100)
          cout << " MPI Bcast (const meshL *) " << a << endl;
      Serialize *buf = 0;
      long nbsize = 0;
      if (who == lmpirank) {
          buf = new Serialize((*a).serialize());
          nbsize = buf->size();
      }
      WBcast(&nbsize, 1, who, comm);
      if (who != lmpirank)
          buf = new Serialize(nbsize, Fem2D::GenericMesh_magicmesh);
      assert(nbsize);
      if (verbosity > 200)
          cout << " size to bcast : " << nbsize << " mpirank : " << mpirank << endl;
        
      WBcast((char *)(*buf), nbsize, who, comm);
        
      if (who != lmpirank) {
          if (a) (*a).decrement();
          Fem2D::MeshL * aa = new Fem2D::MeshL(*buf);
          aa->BuildGTree();
          a = aa;
      }
      delete buf;
      return *this;
  }

    
  template<class R>
  const MPIrank & Bcast(Matrice_Creuse<R> &a) const {
    if(verbosity>100)
      cout << mpirank << ": MPI Bcast " << who << " (Matrice_Creuse &) " << &a << " " << a.A << endl;
    MatriceMorse<R> *mA = 0;
    int ldata[4] = {0, 0, 0, 0};
    if (who == lmpirank) {
      if(a.A) {
        mA = a.pHM();
        ldata[0] = mA->n;
        ldata[1] = mA->m;
        ldata[2] = mA->nnz;
        ldata[3] = mA->half;
      }
    }
    int n4 = 4;
    WBcast(ldata, n4, who, comm);
    if (who != lmpirank && ldata[0])
      mA = new MatriceMorse<R>(ldata[0], ldata[1], ldata[2], ldata[3]);
    else
      CheckPtrHashMatrix(mA, "Bcast 2");
    if(ldata[0]) {
      WBcast(mA->i, ldata[2], who, comm);
      WBcast(mA->j, ldata[2], who, comm);
      WBcast(mA->aij, ldata[2], who, comm);

      mA->Increaze(ldata[2], ldata[2]);
    }
    CheckPtrHashMatrix(mA, "Bcast 2");
    if (who != lmpirank)
      a.A.master(mA);
    // else
    //	delete mA;
    return *this;
  }

  // old async or sync Nov 2010 ...
  template<class R> long Send (Matrice_Creuse<R> * const &a) const;
  template<class R> long Recv (Matrice_Creuse<R> &a) const;
  long Send (Fem2D::Mesh const *a) const;
  long Send (Fem2D::Mesh3 const *a) const;
  long Send (Fem2D::MeshS const *a) const;
  long Send (Fem2D::MeshL const *a) const;
  long Recv (Fem2D::Mesh const *&a) const;
  long Recv (Fem2D::Mesh3 const *&a) const;
  long Recv (Fem2D::MeshS const *&a) const;
  long Recv (Fem2D::MeshL const *&a) const;

  operator int () const { return who; }
};

// for MPI_WAIT_resquets (complex MPI async MPI recv request ) ...
class DoOnWaitMPI_Request : public MPIrank {
public:
  bool sync;
  DoOnWaitMPI_Request( MPIrank mpr) : MPIrank(mpr), sync((rq == 0 || rq == Syncro_block)) {}
  virtual bool Do(MPI_Request *rrq) = 0; // false -> end
  bool DoSR() { // do the Send/Recv Op.
    bool ret = false;
    if (verbosity > 100)
      cout << mpirank << "   --- Do Send/Recv : " << " " << rq << " " << sync << endl;
    if(sync) { // wait ...
      bool c = 1;
      if (verbosity > 100)
        cout << mpirank << "   --- Do way : " << c << " " << rq << endl;
      while (c) {
        c = Do(rq);
        if (verbosity > 100)
          cout << mpirank << "   --- Do return : " << c << " " << rq << endl;
      }

      ret = true;// clean
    } else
      ToDoOnWaitMPI_Request[rq] = this; // add request for WAIT ..
    return ret;
  }
  virtual ~DoOnWaitMPI_Request(){}
  private:
  DoOnWaitMPI_Request (const DoOnWaitMPI_Request &);
  DoOnWaitMPI_Request & operator= ( DoOnWaitMPI_Request &);
};

void DoOnWaitMPIRequest (MPI_Request *rq) {
  if (rq) {
    map<MPI_Request*, DoOnWaitMPI_Request *>:: iterator drd = ToDoOnWaitMPI_Request.find(rq);
    if (drd != ToDoOnWaitMPI_Request.end()) {
      if (verbosity > 100)
        cout << " Do on DoOnWaitMPIRequest " << rq << " " << endl;
      if (!drd->second->Do(rq)) {
        delete drd->second;
        ToDoOnWaitMPI_Request.erase(drd); // finish ...
      }
    }
  }
}

void DeSerialize (Serialize * sTh, Fem2D::Mesh const ** ppTh) {
  if (*ppTh) (**ppTh).decrement();
  Fem2D::Mesh *pTh = new Fem2D::Mesh(*sTh);
  *ppTh = pTh;

  Fem2D::R2 Pn,Px;
  pTh->BoundingBox(Pn, Px);
  pTh->quadtree = new Fem2D::FQuadTree(pTh, Pn, Px, pTh->nv);
}

void DeSerialize (Serialize *sTh, const Fem2D::Mesh3 **ppTh) {
  if (*ppTh) (**ppTh).decrement();
    Fem2D::Mesh3 *pTh;
  int havebordermesh = sTh->havebordermesh();
  if( !havebordermesh) pTh = new Fem2D::Mesh3(*sTh);
  else if(havebordermesh==1) {
      if (verbosity>99) cout << " DeSerialize mesh3:meshS " << endl;
      pTh = new Fem2D::Mesh3(*sTh, havebordermesh); // have a meshS
  }
  pTh->BuildGTree();
  *ppTh = pTh;
}

void DeSerialize (Serialize *sTh, const Fem2D::MeshS **ppTh) {
  if (*ppTh) (**ppTh).decrement();
    Fem2D::MeshS *pTh;
  int havebordermesh = sTh->havebordermesh();
  if( !havebordermesh) pTh = new Fem2D::MeshS(*sTh);
  else if(havebordermesh==1) {
      if (verbosity>99) cout << " DeSerialize mesh3:meshS " << endl;
      pTh = new Fem2D::MeshS(*sTh, havebordermesh); // have a meshS
  }
  pTh->BuildGTree();
  *ppTh = pTh;
}

void DeSerialize (Serialize *sTh, const Fem2D::MeshL **ppTh) {
  if (*ppTh) (**ppTh).decrement();
  Fem2D::MeshL *pTh = new Fem2D::MeshL(*sTh);
  pTh->BuildGTree();
  *ppTh = pTh;
}

template<class R>
class RevcWMatd : public DoOnWaitMPI_Request {
public:
  typedef Matrice_Creuse<R> Mat;
  Matrice_Creuse<R> *pmat;
  MatriceMorse<R> *mA;
  int state;
  int ldata[4];
  RevcWMatd(const MPIrank *mpirank, Mat * pm)
    : DoOnWaitMPI_Request(*mpirank), pmat(pm), mA(0), state(0) {
    int tag = MPI_TAG<Matrice_Creuse<R>* >::TAG;
    int ll = WRecv(ldata, 4, who, tag, comm, rq);
    ffassert(ll == MPI_SUCCESS);
  }

  bool  Do(MPI_Request *rrq) {
    state++;
    int tag = MPI_TAG<Mat *>::TAG;
    if (verbosity > 100)
      cout << mpirank << "  ---R: ldata " << ldata[0] << " " << ldata[1] << " " << ldata[2] << " " <<ldata[3] << " " << state << endl;

    int ll = 0;

    switch (state) {
      case 1:
        mA = new MatriceMorse<R>(ldata[0], ldata[1], ldata[2], ldata[3]);
        ll = WRecv(mA->i, ldata[2], who, tag+1, comm, rq);
        break;
      case 2:
        ll = WRecv(mA->j, ldata[2], who, tag+2, comm, rq);
        break;
      case 3:
        ll = WRecv(mA->aij, ldata[2], who, tag+3, comm, rq);
        break;
      default:
        mA->Increaze(ldata[2], ldata[2]);
        CheckPtrHashMatrix(mA, " WRecv ");
        pmat->A.master(mA);
        mA = 0;
        return false;
        break;
    }

    ffassert(ll == MPI_SUCCESS);

    return true; // OK
  }

  ~RevcWMatd() {
    if(mA) delete mA;
  }
};

template<class R>
class SendWMatd : public DoOnWaitMPI_Request {
public:
  typedef Matrice_Creuse<R> Mat;
  Matrice_Creuse<R> *pmat;
  MatriceMorse<R> *mA;
  int state;
  int ldata[4];
  SendWMatd(const MPIrank *mpirank, Mat *pm)
    : DoOnWaitMPI_Request(*mpirank), pmat(pm), mA(0), state(0) {
    mA = pmat->pHM();
    CheckPtrHashMatrix(mA, " SendWMatd ");
    ldata[0] = mA->n;
    ldata[1] = mA->m;
    ldata[2] = mA->nnz;
    ldata[3] = mA->half;
    int tag = MPI_TAG<Matrice_Creuse<R>* >::TAG;
    int ll = WSend(ldata, 4, who, tag, comm, rq);
    ffassert(ll == MPI_SUCCESS);
  }
  bool Do(MPI_Request *rrq) {
    state++;
    int tag = MPI_TAG<Mat *>::TAG;
    if (verbosity > 100)
      cout << mpirank << "  ---S ldata " << ldata[0] << " " << ldata[1] << " " << ldata[2] << " " <<ldata[3] << " "<< state << endl;

    CheckPtrHashMatrix(mA, " SendWMatd 2");
    int ll = 0;
    switch (state) {
      case 1:
        ll = WSend(mA->i, ldata[2], who, tag+1, comm, rq);
        break;
      case 2:
        ll = WSend(mA->j, ldata[2], who, tag+2, comm, rq);
        break;
      case 3:
        ll = WSend(mA->aij, ldata[2], who, tag+3, comm, rq);
        break;
      default:
        return false;
        break;
    }

    ffassert(ll == MPI_SUCCESS);

    return true; // OK
  }

  ~SendWMatd() {
  }
};

template<class Mesh>
class RevcWMeshd : public DoOnWaitMPI_Request, Serialize {
public:
  Mesh const **ppTh;
  int state;

  RevcWMeshd (const MPIrank *mpirank, Mesh const **ppThh)
    : DoOnWaitMPI_Request(*mpirank), Serialize(sizempibuf, Fem2D::Mesh::magicmesh),
    ppTh(ppThh), state(0) {
    // remark: the first data in p is the size in long long

    int tag = MPI_TAG<Mesh *>::TAG;
    if (verbosity > 99)
      cout << " -- RevcWMeshd " << rq << " " << comm << " " << p << endl;
    int ll = WRecv(p, sizempibuf, who, tag, comm, rq); // wait first part Warning async => not wait.
  }

  bool Do (MPI_Request *rrq) {
    long long lsz;
    int tag = MPI_TAG<Mesh *>::TAG;
    ffassert(rq == rrq);
    size_t kk = 0;
    get(kk, lsz);
    if (verbosity > 199)
      cout << mpirank << "     -- lsk = " << lsz << " p= " << p << " p[]= "
           << (int)p[0] << (int)p[1] << (int)p[2] << (int)p[3] << endl;

    long l1 = lsz - sizempibuf;
    if (verbosity > 100)
      cout << mpirank << " Do RevcWMeshd " << lsz << " " << state << " cont : " << (l1 > 0) << " " << rq << " " << comm << endl;

    if (0 == state++ && l1 > 0 ) { // recv first part ..
      if (verbosity > 100)
        cout << mpirank << " + Do RevcWMeshd " << lsz << " " << state << " cont : " << (l1 > sizempibuf) << " " << rq << " " << l1 << endl;
      resize(lsz);
      int ll = WRecv(p+sizempibuf, l1, who, tag+state, comm, rq);
      return true; // continue ..
    } else resize(lsz);
    // we have the all buffer => DeSerialize
    if(lsz)
        DeSerialize(this, ppTh);
    if (verbosity > 100)
      cout << "    " << mpirank << " recived from " << who << " serialized " << what << ", l="
           << lsz << ", tag=" << tag << " rq = " << rq << " " << *ppTh << endl;
    return false; // OK
  }

  ~RevcWMeshd() {}
};

template<class Mesh>
class SendWMeshd : public DoOnWaitMPI_Request, Serialize {
public:
  const Mesh **ppTh;
  MPI_Request rqSecond;

  SendWMeshd(const MPIrank *mpirank, const Mesh ** ppThh)
    : DoOnWaitMPI_Request(*mpirank), Serialize((**ppThh).serialize()),
    ppTh(ppThh) {
    {
      long long lsz;
      size_t kk = 0;
      get(kk, lsz);
      ffassert(lsz == lg); // verif
    }
    int tag = MPI_TAG<Mesh *>::TAG;

    if (verbosity > 100)
      cout << " -- SendWMeshd " << rq << " " << comm << " " << p << " " << lg << " "<< " p[]= "
           << (int)p[0] << (int)p[1] << (int)p[2] << (int)p[3] << endl;

    if (lg <= sizempibuf)
      WSend(p, lg, who, tag, comm, rq);
    else
      WSend(p, sizempibuf, who, tag, comm, rq);
    long l1 = lg-sizempibuf;
    if (l1 > 0) { // send the second part
      if (verbosity > 100)
        cout << mpirank << " Do SendWMeshd " << lg << " cont : " << (l1 > 0) << " " << rq << " " << comm << endl;
      WSend(p+sizempibuf, l1, who, tag+1, comm, &rqSecond);
    }
    else
      rqSecond = MPI_REQUEST_NULL;
  }

    
  SendWMeshd(const MPIrank *mpirank, const Mesh ** ppThh, bool havebordermesh)
  : DoOnWaitMPI_Request(*mpirank), Serialize((**ppThh).serialize_withBorderMesh()),
  ppTh(ppThh) {
      {
          long long lsz;
          size_t kk = 0;
          get(kk, lsz);
          ffassert(lsz == lg); // verif
          kk=2*sizeof(int);
          int bordermesh=0;
          get(kk, bordermesh);
          ffassert(bordermesh == 1);
          
      }
      int tag = MPI_TAG<Mesh *>::TAG;
      
      if (verbosity > 100)
          cout << " -- SendWMeshd with border mesh" << rq << " " << comm << " " << p << " " << lg << " "<< " p[]= "
          << (int)p[0] << (int)p[1] << (int)p[2] << (int)p[3] << endl;
      
      if (lg <= sizempibuf)
          WSend(p, lg, who, tag, comm, rq);
      else
          WSend(p, sizempibuf, who, tag, comm, rq);
      long l1 = lg-sizempibuf;
      if (l1 > 0) { // send the second part
        if (verbosity > 100)
          cout << mpirank << " Do SendWMeshd " << lg << " cont : " << (l1 > 0) << " " << rq << " " << comm << endl;
        WSend(p+sizempibuf, l1, who, tag+1, comm, &rqSecond);
      }
      else
        rqSecond = MPI_REQUEST_NULL;
  }
  
  bool Do(MPI_Request *rrq) {
      if(rqSecond != MPI_REQUEST_NULL)
        MPI_Wait(&rqSecond, MPI_STATUS_IGNORE);
      return true;// Fini
  }

  ~SendWMeshd() {count()=0;}
};


template<class R>
  long MPIrank::Send(Matrice_Creuse<R> * const &  a) const
  {
	SendWMatd<R> *rwm= new SendWMatd<R>(this,a);
	if( rwm->DoSR() ) delete rwm;
	return MPI_SUCCESS;
  }

  template<class R>
  long MPIrank::Recv(Matrice_Creuse<R>  &  a) const
  {
      {
	RevcWMatd<R> *rwm= new RevcWMatd<R>(this,&a);
        if( rwm->DoSR() ) delete rwm;
        return MPI_SUCCESS;
      }
  }


long MPIrank::Send(const Fem2D::Mesh *  a) const {
    if(verbosity>100)
      cout << " MPI << (mesh *) " << a << endl;
    if(a) {
        SendWMeshd<Mesh> *rwm= new SendWMeshd<Mesh>(this,&a);
        if( rwm->DoSR() ) delete rwm;
    }
    else {
        WSend((char*)NULL, 0, who, MPI_TAG<Fem2D::Mesh*>::TAG, comm, rq);
    }
    return MPI_SUCCESS;
  }


long MPIrank::Send (const Fem2D::Mesh3 *  a) const {
    if(verbosity>100)
      cout << " MPI << (mesh3 *) " << a << endl;
    if(a) {
        SendWMeshd<Mesh3> *rwm= new SendWMeshd<Mesh3>(this,&a);
        if( rwm->DoSR() ) delete rwm;
    }
    else {
        WSend((char*)NULL, 0, who, MPI_TAG<Fem2D::Mesh3*>::TAG, comm, rq);
    }
    return MPI_SUCCESS;
  }

long MPIrank::Send (const Fem2D::MeshS *  a) const {
    if(verbosity>100)
      cout << " MPI << (meshS *) " << a << endl;
    ffassert(a);
    SendWMeshd<MeshS> *rwm= new SendWMeshd<MeshS>(this,&a);
    if( rwm->DoSR() ) delete rwm;
    return MPI_SUCCESS;
  }

long MPIrank::Send (const Fem2D::MeshL *  a) const {
  if(verbosity>100)
    cout << " MPI << (meshL *) " << a << endl;
  ffassert(a);
  SendWMeshd<MeshL> *rwm= new SendWMeshd<MeshL>(this,&a);
  if( rwm->DoSR() ) delete rwm;
  return MPI_SUCCESS;
}

// new version asyncrone ...  Now 2010 ...
long MPIrank::Recv(const Fem2D::Mesh *& a) const  {
    if(verbosity>100)
	cout << " MPI >> (mesh *) &" << a << " " << &a << endl;
    RevcWMeshd<Mesh> *rwm= new RevcWMeshd<Mesh>(this,&a);
    if( rwm->DoSR() ) delete rwm;
    if((rq==0 || rq == Syncro_block))
      ffassert( a );
    return MPI_SUCCESS;
}

long MPIrank::Recv(const Fem2D::Mesh3 *& a) const  {
    if(verbosity>100)
      cout << " MPI >> (mesh3 *) &" << a << " " << &a << endl;
    RevcWMeshd<Mesh3> *rwm= new RevcWMeshd<Mesh3>(this,&a);
    if( rwm->DoSR() ) delete rwm;
    if((rq==0 || rq == Syncro_block))
      ffassert( a );
    return MPI_SUCCESS;
}

long MPIrank::Recv(const Fem2D::MeshS *& a) const  {
    if(verbosity>100)
      cout << " MPI >> (meshS *) &" << a << " " << &a << endl;
    RevcWMeshd<MeshS> *rwm= new RevcWMeshd<MeshS>(this,&a);
    if( rwm->DoSR() ) delete rwm;
    if((rq==0 || rq == Syncro_block))
      ffassert( a );
    return MPI_SUCCESS;
}

long MPIrank::Recv(const Fem2D::MeshL *& a) const  {
    if(verbosity>100)
      cout << " MPI >> (meshL *) &" << a << " " << &a << endl;
    RevcWMeshd<MeshL> *rwm= new RevcWMeshd<MeshL>(this,&a);
    if( rwm->DoSR() ) delete rwm;
    if((rq==0 || rq == Syncro_block))
      ffassert( a );
    return MPI_SUCCESS;
}

void Serialize::mpisend(const MPIrank & rank,long tag,const void * vmpirank)
{
  const MPIrank * mpirank=static_cast<const MPIrank *> (vmpirank);
  MPI_Comm comm=mpirank->comm;
  MPI_Request *rq=mpirank->rq;
  ffassert(rq==0 || rq == Syncro_block);
  char * pp = p-sizeof(long);
  long countsave=count(); // save count
  count()=lg; // store length in count
  long l=lg+sizeof(long);
  if(verbosity>100)
    cout << " -- send from  " << mpirank << " to " << rank << " serialized " << what
	 <<   ", l=" << l << ", tag=" << tag << " " << (l < sizempibuf) << endl;
  if (l <=sizempibuf)
    WSend(pp,l, rank, tag,comm,rq);
  else {
    WSend(pp,sizempibuf,  rank, tag,comm,rq);
    WSend(pp+sizempibuf,l-sizempibuf, rank, tag+1,comm,rq);
  }
  if(verbosity>100)
    cout << "    ok send is arrived " << endl;
  count()=countsave; // restore count
}


Serialize::Serialize(const MPIrank & rank,const char * wht,long tag,const void * vmpirank)
  :what(wht)
{
  const MPIrank * mpirank=static_cast<const MPIrank *> (vmpirank);
  MPI_Comm comm=mpirank->comm;
  MPI_Request *rq=mpirank->rq;

  if(verbosity>100)
    cout << " -- waiting " << mpirank << " from  " << rank << " serialized " << what
	 << " tag = " << tag <<  endl;
  if(!(rq==0 || rq == Syncro_block))
    {
      ExecError("Not async recv of complex  objet!  Sorry to hard to code (FH!).");
      ffassert(rq==0 || rq == Syncro_block);
    }


  char * buf= new char [sizempibuf];
  WRecv(buf, sizempibuf,  rank, tag,comm,rq);
  lg = * (long *) (void *) buf;
  long l=lg+sizeof(long);
  char * pp= new char[l]  ;
  if ( l <= sizempibuf)
    memcpy(pp,buf,l);
  else
    {
      memcpy(pp,buf,sizempibuf);
      WRecv(pp+sizempibuf,l-sizempibuf,  rank, tag+1,comm,rq)  ;
    }

  if(verbosity>100)
    cout << "    " << mpirank << " recived from " << rank << " serialized " << what <<   ", l="
	 << l << ", tag=" << tag << endl;
  delete [] buf;
  p=pp+sizeof(long);
  count()=0;

}

template<class A>
struct Op_Readmpi : public binary_function<MPIrank,A*,MPIrank> {
  static MPIrank  f(MPIrank const  & f,A *  const  & a)
  {
    f.Recv(*a);
    return f;
  }
};

template<class A>
struct Op_Recvmpi : public binary_function<MPIrank,A*,long> {
  static MPIrank  f(MPIrank const  & f,A *  const  & a)
  {
    ffassert(f.rq ==0 || f.rq == Syncro_block); // Block
    return f.Recv(*a);

  }
};
template<class A>
struct Op_IRecvmpi : public binary_function<MPIrank,A*,long> {
  static MPIrank  f(MPIrank const  & f,A *  const  & a)
  {
    ffassert(f.rq !=0 || f.rq != Syncro_block); // no Block
    return f.Recv(*a);

  }
};


template<class A>
struct Op_Writempi : public binary_function<MPIrank,A,MPIrank> {
  static MPIrank  f(MPIrank const  & f,A   const  &  a)
  {
    f.Send(a);
    return f;
  }
};


template<class A>
struct Op_Bcastmpi : public binary_function<MPIrank,A*,MPIrank> {
  static MPIrank  f(MPIrank const  & f,A *  const  & a)
  {
    f.Bcast(*a);
    return f;
   }
};

template<class A>
struct Op_ISendmpi : public binary_function<MPIrank,A,long> {
  static MPIrank  f(MPIrank const  & f,A   const  & a)
  {
    ffassert(f.rq != Syncro_block);
    return f.Send(a);
  }
};
template<class A>
struct Op_Sendmpi : public binary_function<MPIrank,A,long> {
  static MPIrank  f(MPIrank const  & f,A  const  & a)
  {
	MPIrank ff(f.who,f.comm,Syncro_block);
	return ff.Send(a);
  }
};


template<class R>
struct Op_All2All : public binary_function<KN_<R>,KN_<R>,long> {
  static long  f( KN_<R>  const  & s, KN_<R>  const  &r)
  {
      CheckContigueKN(s);
      CheckContigueKN(r);

    MPI_Comm comm=MPI_COMM_WORLD;
    int mpisizew;
    MPI_Comm_size(comm, &mpisizew); /* local */
    int chunk = s.N()/mpisizew;
    ffassert(s.N()==mpisizew*chunk && r.N()==s.N());

    return MPI_Alltoall( (void *) (R*) s, chunk, MPI_TYPE<R>::TYPE(),
			 (void *) (R*) r, chunk, MPI_TYPE<R>::TYPE(), comm);
  }
};


template<class R>
struct Op_Allgather1 : public binary_function<R*,KN_<R>,long> {
  static long  f( R*  const  & s, KN_<R>  const  &r)
    {
      MPI_Comm comm=MPI_COMM_WORLD;
      int mpisizew;
	CheckContigueKN(r);

      MPI_Comm_size(comm, &mpisizew); /* local */
      int chunk = 1;
      ffassert(r.N()==mpisizew);

      return MPI_Allgather( (void *) (R*) s, chunk, MPI_TYPE<R>::TYPE(),
			    (void *) (R*) r, chunk, MPI_TYPE<R>::TYPE(), comm);
    }
};

template<class R>
struct Op_Allgather : public binary_function<KN_<R>,KN_<R>,long> {
  static long  f( KN_<R>  const  & s, KN_<R>  const  &r)
    {
	CheckContigueKN(s);
	CheckContigueKN(r);

      MPI_Comm comm=MPI_COMM_WORLD;
      int mpisizew;
      MPI_Comm_size(comm, &mpisizew); /* local */
      int chunk = r.N()/mpisizew;
      ffassert(r.N()==mpisizew*chunk && chunk == s.N());
      return MPI_Allgather( (void *) (R*) s, chunk, MPI_TYPE<R>::TYPE(),
			    (void *) (R*) r, chunk, MPI_TYPE<R>::TYPE(), comm);
    }
};

template<class R>
struct Op_All2All3 : public ternary_function<KN_<R>,KN_<R>,fMPI_Comm,long> {
  static long  f(Stack, KN_<R>  const  & s, KN_<R>  const  &r,fMPI_Comm const & cmm )
    {
	CheckContigueKN(s);
	CheckContigueKN(r);

      MPI_Comm comm=cmm;
      int mpisizew;
      MPI_Comm_size(comm, &mpisizew); /* local */
      int chunk = s.N()/mpisizew;
      ffassert(s.N()==mpisizew*chunk && r.N()==s.N());

      return MPI_Alltoall( (void *) (R*) s, chunk, MPI_TYPE<R>::TYPE(),
			   (void *) (R*) r, chunk, MPI_TYPE<R>::TYPE(), comm);
    }
};

template<class R>
struct Op_Allgather3 : public ternary_function<KN_<R>,KN_<R>,fMPI_Comm,long> {
  static long  f(Stack, KN_<R>  const  & s, KN_<R>  const  &r,fMPI_Comm const & cmm)
  {
      CheckContigueKN(r);
      CheckContigueKN(s);

    MPI_Comm comm=cmm;
    int mpisizew;
    MPI_Comm_size(comm, &mpisizew); /* local */
    int chunk = r.N()/mpisizew;    // bug corrected by J. Morice
    //ffassert(s.N()==mpisizew*chunk && r.N()==s.N());
    ffassert(s.N()==chunk && r.N()==s.N()*mpisizew);

    return MPI_Allgather( (void *) (R*) s, chunk, MPI_TYPE<R>::TYPE(),
			  (void *) (R*) r, chunk, MPI_TYPE<R>::TYPE(), comm);
  }
};


template<class R>
struct Op_Allgather13 : public ternary_function<R*,KN_<R>,fMPI_Comm,long> {
  static long  f(Stack, R*  const  & s, KN_<R>  const  &r,fMPI_Comm const & cmm)
  {

      CheckContigueKN(r);

    MPI_Comm comm=cmm;
    int mpisizew;
    MPI_Comm_size(comm, &mpisizew); /* local */
    int chunk = 1;    // bug corrected by J. Morice
    ffassert( r.N()==mpisizew);

    return MPI_Allgather( (void *) (R*) s, chunk, MPI_TYPE<R>::TYPE(),
			  (void *) (R*) r, chunk, MPI_TYPE<R>::TYPE(), comm);
  }
};
// Add J. Morice

template<class R>
long  Op_All2Allv( KN_<R>  const  & s, KN_<R>  const  &r, KN_<long> const &sendcnts, KN_<long> const &sdispls, KN_<long> const &recvcnts, KN_<long> const &rdispls)
{
    CheckContigueKN(s);
    CheckContigueKN(r);

  MPI_Comm comm=MPI_COMM_WORLD;
  int mpirankv=MPI_UNDEFINED;
  MPI_Comm_rank(comm, &mpirankv);

  int mpisizew;
  MPI_Comm_size(comm, &mpisizew); /* local */
  ffassert( sendcnts.N() == sdispls.N() && sendcnts.N() == recvcnts.N() && sendcnts.N() == rdispls.N() && sendcnts.N() == mpisizew );

  KN<int> INTsendcnts(sendcnts.N());
  KN<int> INTsdispls(sdispls.N());
  KN<int> INTrecvcnts(recvcnts.N());
  KN<int> INTrdispls(rdispls.N());

  for(int ii=0; ii< sendcnts.N(); ii++){
    INTsendcnts[ii] = sendcnts[ii];
    INTsdispls[ii]  = sdispls[ii];
    INTrecvcnts[ii]  = recvcnts[ii];
    INTrdispls[ii]  = rdispls[ii];
  }

  return MPI_Alltoallv( (void *) (R*) s, INTsendcnts, INTsdispls, MPI_TYPE<R>::TYPE(),
			(void *) (R*) r, INTrecvcnts, INTrdispls, MPI_TYPE<R>::TYPE(), comm);
}



template<class R>
struct Op_Allgatherv : public quad_function<KN_<R>,KN_<R>,KN_<long>,KN_<long>,long> {
  static long f( Stack ,KN_<R>  const  & s, KN_<R>  const  &r, KN_<long> const & recvcount, KN_<long> const & displs)
  {
      CheckContigueKN(r);
      CheckContigueKN(s);

    MPI_Comm comm=MPI_COMM_WORLD;
    int mpisizew;
    MPI_Comm_size(comm, &mpisizew);
    ffassert( recvcount.N() == displs.N() && recvcount.N() == mpisizew);
    long sum=0;
    for(int ii=0; ii< recvcount.N(); ii++)
      sum+=recvcount[ii];
    ffassert( sum == r.N() );

    KN<int> INTrecvcount(recvcount.N());
    KN<int> INTdispls(displs.N());
    for(int ii=0; ii< recvcount.N(); ii++){
      INTrecvcount[ii]= recvcount[ii];
      INTdispls[ii]= displs[ii];
    }
    return MPI_Allgatherv( (void *) (R*) s, s.N(), MPI_TYPE<R>::TYPE(),
			   (void *) (R*) r, INTrecvcount, INTdispls,MPI_TYPE<R>::TYPE(), comm);
  }
};

template<class R>
long  Op_All2All3v(KN_<R>  const  & s, KN_<R>  const  &r,fMPI_Comm const & cmm, KN_<long> const &sendcnts, KN_<long> const &sdispls, KN_<long> const &recvcnts, KN_<long> const &rdispls )
{
    CheckContigueKN(r);
    CheckContigueKN(s);

  MPI_Comm comm=cmm;
  int mpirankv=MPI_UNDEFINED;
  MPI_Comm_rank(comm, &mpirankv);

  int mpisizew;
  MPI_Comm_size(comm, &mpisizew); /* local */
  ffassert( sendcnts.N() == sdispls.N() && sendcnts.N() == recvcnts.N() && sendcnts.N() == rdispls.N() && sendcnts.N() == mpisizew );

  KN<int> INTsendcnts(sendcnts.N());
  KN<int> INTsdispls(sdispls.N());
  KN<int> INTrecvcnts(recvcnts.N());
  KN<int> INTrdispls(rdispls.N());

  for(int ii=0; ii< sendcnts.N(); ii++){
    INTsendcnts[ii] = sendcnts[ii];
    INTsdispls[ii]  = sdispls[ii];
    INTrecvcnts[ii]  = recvcnts[ii];
    INTrdispls[ii]  = rdispls[ii];
  }

  return MPI_Alltoallv( (void *) (R*) s, INTsendcnts, INTsdispls, MPI_TYPE<R>::TYPE(),
			(void *) (R*) r, INTrecvcnts, INTrdispls, MPI_TYPE<R>::TYPE(), comm);
}


template<class R>
long Op_Allgatherv3(KN_<R>  const  & s, KN_<R>  const  &r,fMPI_Comm const & cmm, KN_<long> const & recvcount, KN_<long> const & displs)
{
    CheckContigueKN(r);
    CheckContigueKN(s);

  MPI_Comm comm=cmm;
  int mpisizew;
  MPI_Comm_size(comm, &mpisizew);
  ffassert( recvcount.N() == displs.N() && recvcount.N() == mpisizew);
  long sum=0;
  for(int ii=0; ii< recvcount.N(); ii++)
    sum+=recvcount[ii];
  ffassert( sum == r.N() );
  KN<int> INTrecvcount(recvcount.N());
  KN<int> INTdispls(displs.N());
  for(int ii=0; ii< recvcount.N(); ii++){
    INTrecvcount[ii]= recvcount[ii];
    INTdispls[ii]= displs[ii];
  }

  return MPI_Allgatherv( (void *) (R*) s, s.N(), MPI_TYPE<R>::TYPE(),
			 (void *) (R*) r, INTrecvcount, INTdispls,MPI_TYPE<R>::TYPE(), comm);
}


template<class R>
struct Op_Scatter1 : public   ternary_function<KN_<R>, R* ,MPIrank,long> {
  static long  f(Stack, KN_<R>  const  & s, R*  const  &r,  MPIrank const & root)
  {
      CheckContigueKN(s);


    int mpisizew;
    MPI_Comm_size(root.comm, &mpisizew);
    int chunk = 1;

    return MPI_Scatter( (void *) (R*) s, chunk, MPI_TYPE<R>::TYPE(),
			(void *) (R*) r, chunk, MPI_TYPE<R>::TYPE(),root.who,root.comm);
  }
};


// Fin add J. Morice


template<class R>
struct Op_Scatter3 : public   ternary_function<KN_<R>,KN_<R>,MPIrank,long> {
  static long  f(Stack, KN_<R>  const  & s, KN_<R>  const  &r,  MPIrank const & root)
  {

      CheckContigueKN(r);
      CheckContigueKN(s);

    int mpisizew;
    MPI_Comm_size(root.comm, &mpisizew);
    int chunk = r.N(); // FH  correct  jan 2012 ...

    return MPI_Scatter( (void *) (R*) s, chunk, MPI_TYPE<R>::TYPE(),
			(void *) (R*) r, chunk, MPI_TYPE<R>::TYPE(),root.who,root.comm);
  }
};

// Add J. Morice
template<class R>
long Op_Scatterv3( KN_<R>  const  & s, KN_<R>  const  &r,  MPIrank const & root, KN_<long> const &sendcnts, KN_<long> const &displs)
{
    CheckContigueKN(r);
    CheckContigueKN(s);

  int mpirankv=MPI_UNDEFINED;
  if(root.comm != MPI_COMM_NULL)
    MPI_Comm_rank(root.comm, &mpirankv);

  int mpisizew;
  MPI_Comm_size(root.comm, &mpisizew); /* local */
  KN<int> INTsendcnts(mpirankv == root.who ? sendcnts.N() : 0);
  KN<int> INTdispls(mpirankv == root.who ? sendcnts.N() : 0);
  for(int ii=0; ii< INTsendcnts.N(); ii++){
    INTsendcnts[ii]= sendcnts[ii];
    INTdispls[ii]= displs[ii];
  }

  return MPI_Scatterv( (void *) (R*) s, INTsendcnts, INTdispls, MPI_TYPE<R>::TYPE(),
		       (void *) (R*) r, r.N(), MPI_TYPE<R>::TYPE(),root.who,root.comm);
}

// fin J. Morice
template<class R> uint64_t CodeIJ(const MatriceMorse<R> * pa) { return pa->CodeIJ();}

template<class R>
struct Op_ReduceMat  : public   quad_function<Matrice_Creuse<R>*,Matrice_Creuse<R> *,MPIrank,fMPI_Op,long> {
    static long  f(Stack, Matrice_Creuse<R>*  const  & s,Matrice_Creuse<R>*  const  &r,  MPIrank const & root, fMPI_Op const &op)
    {
	ffassert( r && s);
	MatriceCreuse<R> * sA=s->A;
	MatriceCreuse<R> * rA=r->A;
        int mpirankw,mpisizew;
        MPI_Comm_rank(root.comm, &mpirankw);
        MPI_Comm_size(root.comm, &mpisizew);
        int who = root.who;
	ffassert( sA );
	MatriceMorse<R> * sM = s->pHM();
        MatriceMorse<R> * rM=0;
        if(  !rA && (mpirankw==who) ) { // build a zero matric copy of sM on proc root
            MatriceMorse<R> *rm=new MatriceMorse<R>(*sM); //new MatriceMorse<R>(sM.n,sM.m,sM.nnz,sM.half,0,sM.lg,sM.cl);
            *rm=R(); // set the matrix to Zero ..
            r->A.master(rm);
            rA=r->A;

        }
        if(r)
          rM = r->pHM();

        ffassert( sM);
        cout << " Op_ReduceMat " << rM << " " << who << " " <<mpirankw << " " << sA << " " << rA << endl;
        ffassert( rM ||  (mpirankw!=who) );



        uint64_t scode =  CodeIJ(sM) ;
        uint64_t rcode = rM ? CodeIJ(rM): scode;
        int chunk = rM ? (int) rM->nnz : (int) sM->nnz;
        //  moralement pattern de rM commun
       //  verif the code of matrix
        sM->COO(); // sort
        if(rM) rM->COO(); // sort
        R * saij=sM->aij;
        if(rcode != scode)
        { //  build array to send ..
            if( verbosity> 9) cout <<mpirank <<  " MPI_Reduce sparse  matrix not same pattern in send/recv build data to send  "<< endl;
            int err =0;
            *rM=R();
            for (int k=0; k< sM->nnz;++k )
            {

                R * prij = rM->pij(sM->i[k],sM->j[k]);
                if( prij) *prij = sM->aij[k];
                else err++;
            }
            if (err)  {
                cerr << mpirank << " **** MPI_Reduce sparse  matrix: warning missing term in pattern "<< err << endl;
            }
            saij = new R [sM->nnz];
            copy(rM->aij,rM->aij+rM->nnz,saij);
            *rM=R();
        }
        KN<uint64_t>  allcode(mpisizew);
        if(mpisizew>1)
        {
            KN<uint64_t> code(mpisizew);
            MPI_Gather( (void *) & rcode  , 1,  MPI_UNSIGNED_LONG_LONG,
                       (void *)  &code[0] , 1,  MPI_UNSIGNED_LONG_LONG, root.who ,root.comm);
            int ok=1;
            if(mpirankw==root.who) //  verif same patern ...
                for(int i=0; i<mpisizew;++i)
                    ok = ok &&  (rcode== code[i]);
            if( !ok)
            {
                cerr << " Fatal error  MPI_reduce mat: the pattern of the all recv matrix are not the same "<< endl;
                cerr << " set the recv matrix with same patten" <<endl;
                 for(int i=0; i<mpisizew;++i)
                     cout << " proc "<< i << " pattern code: "<< code[i] << " != " << rcode << endl;;
                ffassert(0);
            }

        }

	long ret= MPI_Reduce( (void *)  saij,(void *)  rM->aij, chunk , MPI_TYPE<R>::TYPE(),op,root.who,root.comm);
        if( saij != sM->aij) delete [] saij;
        return ret;
    }
};

template<class R>
struct Op_AllReduceMat  : public   quad_function<Matrice_Creuse<R>*,Matrice_Creuse<R> *,fMPI_Comm,fMPI_Op,long> {
    static long  f(Stack, Matrice_Creuse<R>*  const  & s,Matrice_Creuse<R>*  const  &r,  fMPI_Comm const & comm, fMPI_Op const &op)
    {
        ffassert( r && s);
        MatriceCreuse<R> * sA=s->A;
        MatriceCreuse<R> * rA=r->A;
        ffassert( sA );

        MatriceMorse<R> * sM = s->pHM();
        if( ! rA ) { // build a zero matric copy of sM
            MatriceMorse<R> *rm=new MatriceMorse<R>(*sM); //new MatriceMorse<R>(sM.n,sM.m,sM.nnz,sM.half,0,sM.lg,sM.cl);
            *rm=R(); // set the matrix to Zero ..
            r->A.master(rm);
            rA=r->A;

        }

        MatriceMorse<R> * rM = r->pHM();

	ffassert( sM && rM);
        sM->COO(); // sort
        rM->COO(); // sort
	int nnz = (int) sM->nnz;
        uint64_t rcode = CodeIJ(rM);
        uint64_t scode = ( sM != rM) ? CodeIJ(sM) : rcode;
        R * saij=sM->aij;
        if(rcode != scode)
        { //  build array to send ..
             int err =0;
            *rM=R();
            for (int k=0; k< sM->nnz;++k )
            {

                R * prij = rM->pij(sM->i[k],sM->j[k]);
                if( prij) *prij = sM->aij[k];
                else err++;
            }
            nnz = (int) rM->nnz;
            saij = new R [nnz];
            nnz = (int) rM->nnz;
            if( verbosity> 9)
                cout <<mpirank <<  " ** Warning: MPI_AllReduce sparse  matrix not same pattern in send/recv build data to send, build send:  R nnz "<< nnz << " S nnz " << sM->nnz<<   endl;

            copy(rM->aij,rM->aij+rM->nnz,saij);
            *rM=R();

        }
        int mpirankw,mpisizew;
        MPI_Comm_rank(comm, &mpirankw);
        MPI_Comm_size(comm, &mpisizew);

        KN<uint64_t>  allcode(mpisizew);
        if((mpisizew>1)  )
        {
            int chunk = 1;
            KN<uint64_t> code(mpisizew);
            MPI_Gather( (void *) & rcode  , 1,  MPI_UNSIGNED_LONG_LONG,
                         (void *)  &code[0] , 1,  MPI_UNSIGNED_LONG_LONG, 0 ,comm);
            int ok=1;
            if(mpirankw==0)
                for(int i=1; i<mpisizew;++i)
                    ok |= (rcode== code[i]);
            if(verbosity>99) cout << " MPI_Allreduce mat: revif " <<  mpirankw  << " " << ok << " / " << sM->aij << " -> " << rM->aij << " " << comm <<" " <<  chunk << endl;
            if( !ok)
            {
                cerr << " Fatal error  MPI_Allreduce matrix the pattern of the all recv matrix are not the same "<< endl;
                cerr << " set the recv matrix with same patten" <<endl;

            }
            ffassert(ok);
        }
	long ret= MPI_Allreduce( (void *)  saij,(void *)  rM->aij, nnz , MPI_TYPE<R>::TYPE(),op,comm);
        if( saij != sM->aij) delete [] saij;
        return ret;
    }
};



template<class R>
struct Op_Reduce  : public   quad_function<KN_<R>,KN_<R>,MPIrank,fMPI_Op,long> {
  static long  f(Stack, KN_<R>  const  & s, KN_<R>  const  &r,  MPIrank const & root, fMPI_Op const &op)
  {
      CheckContigueKN(r);
      CheckContigueKN(s);

    int chunk = s.N();
    ffassert(chunk==r.N());

    return MPI_Reduce( (void *) (R*) s,(void *) (R*) r, chunk , MPI_TYPE<R>::TYPE(),op,root.who,root.comm);
  }
};

template<class R>
struct Op_AllReduce  : public   quad_function<KN_<R>,KN_<R>,fMPI_Comm,fMPI_Op,long> {
  static long  f(Stack, KN_<R>  const  & s, KN_<R>  const  &r,  fMPI_Comm const & comm,fMPI_Op const &op)
  {
      CheckContigueKN(r);
      CheckContigueKN(s);

    int chunk = s.N();
    ffassert(chunk==r.N());
    return MPI_Allreduce( (void *) (R*) s,(void *) (R*) r, chunk , MPI_TYPE<R>::TYPE(),op,comm);
  }
};

template<class R>
struct Op_AllReduce1  : public   quad_function<R *,R *,fMPI_Comm,fMPI_Op,long> {
    static long  f(Stack, R *  const  & s, R *  const  &r,  fMPI_Comm const & comm,fMPI_Op const &op)
    {
	int chunk = 1;
	return MPI_Allreduce( (void *) (R*) s,(void *) (R*) r, 1 , MPI_TYPE<R>::TYPE(),op,comm);
    }
};

template<class R>
struct Op_Reduce1  : public   quad_function<R*,R*,MPIrank,fMPI_Op,long> {
  static long  f(Stack, R*  const  & s, R*  const  &r,  MPIrank const & root, fMPI_Op const &op)
  {
    int chunk = 1;
    return MPI_Reduce( (void *) (R*) s,(void *) (R*) r, chunk , MPI_TYPE<R>::TYPE(),op,root.who,root.comm);
  }
};

// Add J. Morice
template<class R>
struct Op_Gather1 : public   ternary_function<R*,KN_<R>,MPIrank,long> {
    static long  f(Stack, R* const  & s, KN_<R>  const  &r,  MPIrank const & root)
  {

    int mpisizew,myrank;
    MPI_Comm_size(root.comm, &mpisizew);
    MPI_Comm_rank( root.comm, &myrank)  ;
    int chunk = 1;

    return MPI_Gather( (void *) (R*) s, chunk, MPI_TYPE<R>::TYPE(),
			   (void *) (R*) r, chunk, MPI_TYPE<R>::TYPE(),root.who,root.comm);
  }
};

// Fin Add J. Morice


template<class R>
struct Op_Gather3 : public   ternary_function<KN_<R>,KN_<R>,MPIrank,long> {
    static long  f(Stack, KN_<R>  const  & s, KN_<R>  const  &r,  MPIrank const & root)
  {
      CheckContigueKN(r);
      CheckContigueKN(s);
      int mpisizew,myrank;
      MPI_Comm_size(root.comm, &mpisizew);
      MPI_Comm_rank(root.comm, &myrank)  ;

    int chunk = s.N();

    return MPI_Gather( (void *) (R*) s, chunk, MPI_TYPE<R>::TYPE(),
			   (void *) (R*) r, chunk, MPI_TYPE<R>::TYPE(),root.who,root.comm);
  }
};

// Add by J. Morice

template<class R>
long  Op_Gatherv3(KN_<R>  const  & s, KN_<R>  const  &r,  MPIrank const & root, KN_<long> const & recvcount, KN_<long> const & displs)
{

    CheckContigueKN(r);
    CheckContigueKN(s);

  int mpirankw;
  MPI_Comm_rank(root.comm, &mpirankw);
  int mpisizew;
  MPI_Comm_size(root.comm, &mpisizew);
  KN<int> INTrecvcount(mpirankw == root.who ? recvcount.N() : 0);
  KN<int> INTdispls(mpirankw == root.who ? recvcount.N() : 0);
  for(int ii=0; ii< INTrecvcount.N(); ii++){
    INTrecvcount[ii]= recvcount[ii];
    INTdispls[ii]= displs[ii];
  }

  return MPI_Gatherv( (void *) (R*) s, s.N(), MPI_TYPE<R>::TYPE(),
		      (void *) (R*) r, INTrecvcount, INTdispls,MPI_TYPE<R>::TYPE(),root.who,root.comm);
}

// Fin Add J. Morice


// Add J. Morice communications entre processeurs complex

template<>
struct Op_All2All<Complex> : public binary_function<KN_<Complex>,KN_<Complex>,long> {
  static long  f( KN_<Complex>  const  & s, KN_<Complex>  const  &r)
  {
      CheckContigueKN(r);
      CheckContigueKN(s);

    MPI_Comm comm=MPI_COMM_WORLD;
    int mpisizew;
    MPI_Comm_size(comm, &mpisizew); /* local */
    int chunk = s.N()/mpisizew;
    ffassert(s.N()==mpisizew*chunk && r.N()==s.N());

#ifdef HAVE_MPI_DOUBLE_COMPLEX
    return MPI_Alltoall( (void *) (Complex*) s, chunk, MPI_DOUBLE_COMPLEX,
			 (void *) (Complex*) r, chunk, MPI_DOUBLE_COMPLEX, comm);
#else
    chunk*=2;
    return MPI_Alltoall( reinterpret_cast<void*> ( (Complex*) s), chunk, MPI_DOUBLE,
			 reinterpret_cast<void*> ( (Complex*) r), chunk, MPI_DOUBLE, comm);
#endif
  }
};

template<>
struct Op_Allgather1<Complex> : public binary_function<Complex *,KN_<Complex>,long> {
  static long  f( Complex *  const  & s, KN_<Complex>  const  &r)
  {
      CheckContigueKN(r);


    MPI_Comm comm=MPI_COMM_WORLD;
    int mpisizew;
    MPI_Comm_size(comm, &mpisizew); /* local */
    int chunk = 1;
    ffassert( r.N()== mpisizew);
#ifdef HAVE_MPI_DOUBLE_COMPLEX
    return MPI_Allgather( (void *) (Complex*) s, chunk, MPI_DOUBLE_COMPLEX,
			  (void *) (Complex*) r, chunk, MPI_DOUBLE_COMPLEX, comm);
#else
    chunk*=2;
    return MPI_Allgather( reinterpret_cast<void*> ( (Complex*) s), chunk, MPI_DOUBLE,
			  reinterpret_cast<void*> ( (Complex*) r), chunk, MPI_DOUBLE, comm);
#endif
  }
};


template<>
struct Op_Allgather<Complex> : public binary_function<KN_<Complex>,KN_<Complex>,long> {
  static long  f( KN_<Complex>  const  & s, KN_<Complex>  const  &r)
    {
	CheckContigueKN(r);
	CheckContigueKN(s);

      MPI_Comm comm=MPI_COMM_WORLD;
      int mpisizew;
      MPI_Comm_size(comm, &mpisizew); /* local */
      int chunk = r.N()/mpisizew;
      ffassert( r.N()==chunk*mpisizew && chunk==s.N() );
#ifdef HAVE_MPI_DOUBLE_COMPLEX
      return MPI_Allgather( (void *) (Complex*) s, chunk, MPI_DOUBLE_COMPLEX,
			    (void *) (Complex*) r, chunk, MPI_DOUBLE_COMPLEX, comm);
#else
      chunk*=2;
      return MPI_Allgather( reinterpret_cast<void*> (  (Complex*) s), chunk, MPI_DOUBLE,
			    reinterpret_cast<void*> ( (Complex*) r), chunk, MPI_DOUBLE, comm);
#endif
    }
};

template<>
struct Op_All2All3<Complex> : public ternary_function<KN_<Complex>,KN_<Complex>,fMPI_Comm,long> {
  static long  f(Stack, KN_<Complex>  const  & s, KN_<Complex>  const  &r,fMPI_Comm const & cmm )
    {
	CheckContigueKN(r);
	CheckContigueKN(s);

      MPI_Comm comm=cmm;
      int mpisizew;
      MPI_Comm_size(comm, &mpisizew); /* local */
      int chunk = s.N()/mpisizew;
      ffassert(s.N()==mpisizew*chunk && r.N()==s.N());
#ifdef HAVE_MPI_DOUBLE_COMPLEX
      return MPI_Alltoall( (void *) (Complex*) s, chunk, MPI_DOUBLE_COMPLEX,
			   (void *) (Complex*) r, chunk, MPI_DOUBLE_COMPLEX, comm);
#else
      chunk*=2;
      return MPI_Alltoall( (void *) (Complex*) s, chunk, MPI_DOUBLE,
			   (void *) (Complex*)  (r), chunk, MPI_DOUBLE, comm);
#endif
    }
};

template<>
struct Op_Allgather3<Complex> : public ternary_function<KN_<Complex>,KN_<Complex>,fMPI_Comm,long> {
  static long  f(Stack, KN_<Complex>  const  & s, KN_<Complex>  const  &r,fMPI_Comm const & cmm)
  {
      CheckContigueKN(r);
      CheckContigueKN(s);

    MPI_Comm comm=cmm;
    int mpisizew;
    MPI_Comm_size(comm, &mpisizew); /* local */
    int chunk = r.N()/mpisizew;    // bug corrected by J. Morice
    ffassert(s.N()==chunk && r.N()==s.N()*mpisizew);
#ifdef HAVE_MPI_DOUBLE_COMPLEX
    return MPI_Allgather( (void *) (Complex*) s, chunk, MPI_DOUBLE_COMPLEX,
			  (void *) (Complex*) r, chunk, MPI_DOUBLE_COMPLEX, comm);
#else
    chunk*=2;
    return MPI_Allgather( (void *) (Complex*) (s), chunk, MPI_DOUBLE,
			 (void *) (Complex*) (r), chunk, MPI_DOUBLE, comm);
#endif
  }
};

template<>
struct Op_Allgather13<Complex> : public ternary_function<Complex *,KN_<Complex>,fMPI_Comm,long> {
  static long  f(Stack, Complex *  const  & s, KN_<Complex>  const  &r,fMPI_Comm const & cmm)
  {
      CheckContigueKN(r);


    MPI_Comm comm=cmm;
    int mpisizew;
    MPI_Comm_size(comm, &mpisizew); /* local */
    int chunk = 1;
    ffassert( r.N()==mpisizew);
#ifdef HAVE_MPI_DOUBLE_COMPLEX
    return MPI_Allgather( (void *) (Complex*) s, chunk, MPI_DOUBLE_COMPLEX,
			  (void *) (Complex*) r, chunk, MPI_DOUBLE_COMPLEX, comm);
#else
    chunk*=2;
    return MPI_Allgather((void *) (Complex*)(s), chunk, MPI_DOUBLE,
			   (void *) (Complex*) (r), chunk, MPI_DOUBLE, comm);
#endif
  }
};


template<>
long  Op_All2Allv<Complex>( KN_<Complex>  const  & s, KN_<Complex>  const  &r, KN_<long> const &sendcnts, KN_<long> const &sdispls, KN_<long> const &recvcnts, KN_<long> const &rdispls)
{
    CheckContigueKN(r);
    CheckContigueKN(s);

  MPI_Comm comm=MPI_COMM_WORLD;
  int mpirankv=MPI_UNDEFINED;
  MPI_Comm_rank(comm, &mpirankv);

  int mpisizew;
  MPI_Comm_size(comm, &mpisizew); /* local */
  ffassert( sendcnts.N() == sdispls.N() && sendcnts.N() == recvcnts.N() && sendcnts.N() == rdispls.N() && sendcnts.N() == mpisizew );

  KN<int> INTsendcnts(sendcnts.N());
  KN<int> INTsdispls(sdispls.N());
  KN<int> INTrecvcnts(recvcnts.N());
  KN<int> INTrdispls(rdispls.N());


#ifdef HAVE_MPI_DOUBLE_COMPLEX
  for(int ii=0; ii< sendcnts.N(); ii++){
    INTsendcnts[ii] = sendcnts[ii];
    INTsdispls[ii]  = sdispls[ii];
    INTrecvcnts[ii]  = recvcnts[ii];
    INTrdispls[ii]  = rdispls[ii];
  }
  return MPI_Alltoallv( (void *) (Complex*) s, INTsendcnts, INTsdispls, MPI_DOUBLE_COMPLEX,
			(void *) (Complex*) r, INTrecvcnts, INTrdispls, MPI_DOUBLE_COMPLEX, comm);
#else
  for(int ii=0; ii< sendcnts.N(); ii++){
    INTsendcnts[ii] = 2*sendcnts[ii];
    INTsdispls[ii]  = 2*sdispls[ii];
    INTrecvcnts[ii] = 2*recvcnts[ii];
    INTrdispls[ii]  = 2*rdispls[ii];
  }
  return MPI_Alltoallv( reinterpret_cast<void*> ( (Complex*) s), INTsendcnts, INTsdispls, MPI_DOUBLE,
			reinterpret_cast<void*> ( (Complex*) r), INTrecvcnts, INTrdispls, MPI_DOUBLE, comm);
#endif
}



template<>
struct Op_Allgatherv<Complex> : public quad_function<KN_<Complex>,KN_<Complex>,KN_<long>,KN_<long>,long> {
  static long f( Stack ,KN_<Complex>  const  & s, KN_<Complex>  const  &r, KN_<long> const & recvcount, KN_<long> const & displs)
  {
      CheckContigueKN(r);
      CheckContigueKN(s);

    MPI_Comm comm=MPI_COMM_WORLD;
    int mpisizew;
    MPI_Comm_size(comm, &mpisizew);
    ffassert( recvcount.N() == displs.N() && recvcount.N() == mpisizew);
    long sum=0;
    for(int ii=0; ii< recvcount.N(); ii++)
      sum+=recvcount[ii];
    ffassert( sum == r.N() );


    KN<int> INTrecvcount(recvcount.N());
    KN<int> INTdispls(displs.N());
#ifdef HAVE_MPI_DOUBLE_COMPLEX
    for(int ii=0; ii< recvcount.N(); ii++){
      INTrecvcount[ii]= recvcount[ii];
      INTdispls[ii]= displs[ii];
    }
    return MPI_Allgatherv( (void *) (Complex*)s, s.N(), MPI_DOUBLE_COMPLEX,
			   (void *) (Complex*)r, INTrecvcount, INTdispls,MPI_DOUBLE_COMPLEX, comm);
#else
    for(int ii=0; ii< recvcount.N(); ii++){
      INTrecvcount[ii]= 2*recvcount[ii];
      INTdispls[ii]= 2*displs[ii];
    }
    return MPI_Allgatherv( reinterpret_cast<void*> ( (Complex*) s), 2*s.N(), MPI_DOUBLE,
			   reinterpret_cast<void*> ( (Complex*) r), INTrecvcount, INTdispls,MPI_DOUBLE, comm);
#endif
  }
};

template<>
long  Op_All2All3v<Complex>(KN_<Complex>  const  & s, KN_<Complex>  const  &r,fMPI_Comm const & cmm, KN_<long> const &sendcnts, KN_<long> const &sdispls, KN_<long> const &recvcnts, KN_<long> const &rdispls )
{
    CheckContigueKN(r);
    CheckContigueKN(s);

  MPI_Comm comm=cmm;
  int mpirankv=MPI_UNDEFINED;
  MPI_Comm_rank(comm, &mpirankv);

  int mpisizew;
  MPI_Comm_size(comm, &mpisizew); /* local */
  ffassert( sendcnts.N() == sdispls.N() && sendcnts.N() == recvcnts.N() && sendcnts.N() == rdispls.N() && sendcnts.N() == mpisizew );

  //ffassert(s.N()==sendcnts[mpirankv] && r.N()==recvbuf[mpirankv]);

  KN<int> INTsendcnts(sendcnts.N());
  KN<int> INTsdispls(sdispls.N());
  KN<int> INTrecvcnts(recvcnts.N());
  KN<int> INTrdispls(rdispls.N());

#ifdef HAVE_MPI_DOUBLE_COMPLEX
  for(int ii=0; ii< sendcnts.N(); ii++){
    INTsendcnts[ii] = sendcnts[ii];
    INTsdispls[ii]  = sdispls[ii];
    INTrecvcnts[ii]  = recvcnts[ii];
    INTrdispls[ii]  = rdispls[ii];
  }

  return MPI_Alltoallv( (void *) (Complex*)s, INTsendcnts, INTsdispls, MPI_DOUBLE_COMPLEX,
			(void *) (Complex*)r, INTrecvcnts, INTrdispls, MPI_DOUBLE_COMPLEX, comm);
#else
  for(int ii=0; ii< sendcnts.N(); ii++){
    INTsendcnts[ii] = 2*sendcnts[ii];
    INTsdispls[ii]  = 2*sdispls[ii];
    INTrecvcnts[ii]  = 2*recvcnts[ii];
    INTrdispls[ii]  = 2*rdispls[ii];
  }

  return MPI_Alltoallv( reinterpret_cast<void*> ( (Complex*) s), INTsendcnts, INTsdispls, MPI_DOUBLE,
			reinterpret_cast<void*> ( (Complex*) r), INTrecvcnts, INTrdispls, MPI_DOUBLE, comm);
#endif
}


template<>
long Op_Allgatherv3<Complex>(KN_<Complex>  const  & s, KN_<Complex>  const  &r,fMPI_Comm const & cmm, KN_<long> const & recvcount, KN_<long> const & displs)
{
    CheckContigueKN(r);
    CheckContigueKN(s);

  MPI_Comm comm=cmm;
  int mpisizew;
  MPI_Comm_size(comm, &mpisizew);
  ffassert( recvcount.N() == displs.N() && recvcount.N() == mpisizew);
  long sum=0;
  for(int ii=0; ii< recvcount.N(); ii++)
    sum+=recvcount[ii];
  ffassert( sum == r.N() );
  KN<int> INTrecvcount(recvcount.N());
  KN<int> INTdispls(displs.N());

#ifdef HAVE_MPI_DOUBLE_COMPLEX
  for(int ii=0; ii< recvcount.N(); ii++){
    INTrecvcount[ii]= recvcount[ii];
    INTdispls[ii]= displs[ii];
  }

  return MPI_Allgatherv( (void *) (Complex*)s, s.N(), MPI_DOUBLE_COMPLEX,
			 (void *) (Complex*)r, INTrecvcount, INTdispls,MPI_DOUBLE_COMPLEX, comm);

#else
  for(int ii=0; ii< recvcount.N(); ii++){
    INTrecvcount[ii]= 2*recvcount[ii];
    INTdispls[ii]= 2*displs[ii];
  }

  return MPI_Allgatherv( reinterpret_cast<void*> ( (Complex*) s), 2*s.N(), MPI_DOUBLE,
			 reinterpret_cast<void*> ( (Complex*) r), INTrecvcount, INTdispls,MPI_DOUBLE, comm);
#endif

}

template<>
struct Op_Scatter1<Complex> : public   ternary_function<KN_<Complex>,Complex *,MPIrank,long> {
  static long  f(Stack, KN_<Complex> const  & s, Complex *  const  &r,  MPIrank const & root)
  {

      CheckContigueKN(s);

    int mpisizew;
    MPI_Comm_size(root.comm, &mpisizew);
    int chunk = 1;
#ifdef HAVE_MPI_DOUBLE_COMPLEX
    return MPI_Scatter( (void *) (Complex*)s, chunk, MPI_DOUBLE_COMPLEX,
			(void *) (Complex*)r, chunk, MPI_DOUBLE_COMPLEX,root.who,root.comm);
#else
    chunk*=2;
    return MPI_Scatter( reinterpret_cast<void*> ( (Complex*) s), chunk, MPI_DOUBLE,
			reinterpret_cast<void*> ( (Complex*) r), chunk, MPI_DOUBLE,root.who,root.comm);
#endif
  }
};


template<>
struct Op_Scatter3<Complex> : public   ternary_function<KN_<Complex>,KN_<Complex>,MPIrank,long> {
  static long  f(Stack, KN_<Complex>  const  & s, KN_<Complex>  const  &r,  MPIrank const & root)
  {
      CheckContigueKN(r);
      CheckContigueKN(s);

    int mpisizew;
    MPI_Comm_size(root.comm, &mpisizew);
    int chunk = r.N();// correct 2012 FH
#ifdef HAVE_MPI_DOUBLE_COMPLEX
    return MPI_Scatter( (void *) (Complex*)s, chunk, MPI_DOUBLE_COMPLEX,
			(void *) (Complex*)r, chunk, MPI_DOUBLE_COMPLEX,root.who,root.comm);
#else
    chunk*=2;
    return MPI_Scatter( reinterpret_cast<void*> ( (Complex*) s), chunk, MPI_DOUBLE,
			reinterpret_cast<void*> ( (Complex*) r), chunk, MPI_DOUBLE,root.who,root.comm);
#endif
  }
};

template<>
long Op_Scatterv3<Complex>( KN_<Complex>  const  & s, KN_<Complex>  const  &r,  MPIrank const & root, KN_<long> const &sendcnts, KN_<long> const &displs)
{

    CheckContigueKN(r);
    CheckContigueKN(s);

  int mpirankv=MPI_UNDEFINED;
  if(root.comm != MPI_COMM_NULL)
    MPI_Comm_rank(root.comm, &mpirankv);

  int mpisizew;
  MPI_Comm_size(root.comm, &mpisizew); /* local */
  long sumsize=0;
  for(int ii=0; ii<sendcnts.N(); ii++){
    sumsize += sendcnts[ii];
  }

  KN<int> INTsendcnts(sendcnts.N());
  KN<int> INTdispls(displs.N());
#ifdef HAVE_MPI_DOUBLE_COMPLEX
  for(int ii=0; ii< sendcnts.N(); ii++){
    INTsendcnts[ii]= sendcnts[ii];
    INTdispls[ii]= displs[ii];
  }

  return MPI_Scatterv( (void *) (Complex*)s, INTsendcnts, INTdispls, MPI_DOUBLE_COMPLEX,
		       (void *) (Complex*)r, r.N(), MPI_DOUBLE_COMPLEX,root.who,root.comm);
#else
  for(int ii=0; ii< sendcnts.N(); ii++){
    INTsendcnts[ii]= 2*sendcnts[ii];
    INTdispls[ii]= 2*displs[ii];
  }

  return MPI_Scatterv( reinterpret_cast<void*> ( (Complex*) s), INTsendcnts, INTdispls, MPI_DOUBLE,
		       reinterpret_cast<void*> ( (Complex*) r), 2*r.N(), MPI_DOUBLE,root.who,root.comm);
#endif

}


template<>
struct Op_Reduce<Complex>  : public   quad_function<KN_<Complex>,KN_<Complex>,MPIrank,fMPI_Op,long> {
  static long  f(Stack, KN_<Complex>  const  & s, KN_<Complex>  const  &r,  MPIrank const & root, fMPI_Op const &op)
  {
      CheckContigueKN(r);
      CheckContigueKN(s);

    int chunk = s.N();
    ffassert(chunk==r.N());
#ifdef HAVE_MPI_DOUBLE_COMPLEX
    return MPI_Reduce( (void *) (Complex*)s,(void *) (Complex*)r, chunk , MPI_DOUBLE_COMPLEX,op,root.who,root.comm);
#else
    chunk*=2;
    return MPI_Reduce( reinterpret_cast<void*> ( (Complex*) s), reinterpret_cast<void*> ( (Complex*) r), chunk , MPI_DOUBLE,op,root.who,root.comm);
#endif
  }
};

template<>
struct Op_AllReduce<Complex>  : public   quad_function<KN_<Complex>,KN_<Complex>,fMPI_Comm,fMPI_Op,long> {
  static long  f(Stack, KN_<Complex>  const  & s, KN_<Complex>  const  &r,  fMPI_Comm const & comm,fMPI_Op const &op)
  {
      CheckContigueKN(r);
      CheckContigueKN(s);

    int chunk = s.N();
    ffassert(chunk==r.N());
#ifdef HAVE_MPI_DOUBLE_COMPLEX
    return MPI_Allreduce( (void *) (Complex*)s,(void *) (Complex*)r, chunk , MPI_DOUBLE_COMPLEX,op,comm);
#else
    chunk *=2;
    return MPI_Allreduce( reinterpret_cast<void*> ( (Complex*) s), reinterpret_cast<void*> ( (Complex*) r), chunk , MPI_DOUBLE,op,comm);
#endif
  }
};

template<>
struct Op_Reduce1<Complex>  : public   quad_function<Complex*,Complex*,MPIrank,fMPI_Op,long> {
  static long  f(Stack, Complex*  const  & s, Complex*  const  &r,  MPIrank const & root, fMPI_Op const &op)
  {
#ifdef HAVE_MPI_DOUBLE_COMPLEX
    int chunk = 1;
    return MPI_Reduce( (void *) s, (void *) r, chunk , MPI_DOUBLE_COMPLEX,op,root.who,root.comm);
#else
    int chunk = 2;
    return MPI_Reduce( reinterpret_cast<void*> ( (Complex*) s), reinterpret_cast<void*> ( (Complex*) r), chunk , MPI_DOUBLE,op,root.who,root.comm);
#endif
  }
};

// Add J. Morice
template<>
struct Op_Gather1<Complex> : public   ternary_function<Complex* ,KN_<Complex>,MPIrank,long> {
    static long  f(Stack, Complex * const  & s, KN_<Complex>  const  &r,  MPIrank const & root)
  {

      CheckContigueKN(r);


      int mpisizew,myrank;
      MPI_Comm_size(root.comm, &mpisizew);
      MPI_Comm_rank( root.comm, &myrank)  ;

    int chunk = 1;
#ifdef HAVE_MPI_DOUBLE_COMPLEX

    return MPI_Gather( (void *) (Complex*) s, chunk, MPI_DOUBLE_COMPLEX,
			   (void *) (Complex*) r, chunk, MPI_DOUBLE_COMPLEX, root.who, root.comm);
#else
    chunk = 2;
    return MPI_Gather( reinterpret_cast<void*> ( (Complex*) s), chunk, MPI_DOUBLE,
		       reinterpret_cast<void*> ( (Complex*) r), chunk, MPI_DOUBLE, root.who, root.comm);
#endif
  }
};

// Fin Add J. Morice


template<>
struct Op_Gather3<Complex> : public   ternary_function<KN_<Complex>,KN_<Complex>,MPIrank,long> {
    static long  f(Stack, KN_<Complex>  const  & s, KN_<Complex>  const  &r,  MPIrank const & root)
  {
      CheckContigueKN(r);
      CheckContigueKN(s);

      int mpisizew,myrank;
      MPI_Comm_size(root.comm, &mpisizew);
      MPI_Comm_rank( root.comm, &myrank)  ;

    int chunk = s.N();
#ifdef HAVE_MPI_DOUBLE_COMPLEX

    return MPI_Gather( (void *) (Complex*)s, chunk, MPI_DOUBLE_COMPLEX,
			   (void *) (Complex*)r, chunk, MPI_DOUBLE_COMPLEX,root.who,root.comm);
#else
    chunk *= 2;
    return MPI_Gather( reinterpret_cast<void*> ( (Complex*) s), chunk, MPI_DOUBLE,
		       reinterpret_cast<void*> ( (Complex*) r), chunk, MPI_DOUBLE,root.who,root.comm);
#endif
  }
};


template<>
long  Op_Gatherv3<Complex>(KN_<Complex>  const  & s, KN_<Complex>  const  &r,  MPIrank const & root, KN_<long> const & recvcount, KN_<long> const & displs)
{
    CheckContigueKN(r);
    CheckContigueKN(s);

  int mpirankw;
  MPI_Comm_rank(root.comm, &mpirankw);
  int mpisizew;
  MPI_Comm_size(root.comm, &mpisizew);

  if( mpirankw == root.who){
    long sum=0;
    for(int ii=0; ii< recvcount.N(); ii++)
      sum+=recvcount[ii];
    ffassert( sum == r.N() );
  }
  KN<int> INTrecvcount(recvcount.N());
  KN<int> INTdispls(displs.N());
#ifdef HAVE_MPI_DOUBLE_COMPLEX
  for(int ii=0; ii< recvcount.N(); ii++){
    INTrecvcount[ii]= recvcount[ii];
    INTdispls[ii]= displs[ii];
  }
  return MPI_Gatherv( (void *) (Complex*)s, s.N(), MPI_DOUBLE_COMPLEX,
		      (void *) (Complex*)r, INTrecvcount, INTdispls,MPI_DOUBLE_COMPLEX,root.who,root.comm);
#else
  for(int ii=0; ii< recvcount.N(); ii++){
    INTrecvcount[ii]= 2*recvcount[ii];
    INTdispls[ii]= 2*displs[ii];
  }
  return MPI_Gatherv( reinterpret_cast<void*> ( (Complex*) s), 2*s.N(), MPI_DOUBLE,
		      reinterpret_cast<void*> ( (Complex*) r), INTrecvcount, INTdispls,MPI_DOUBLE,root.who,root.comm);
#endif
}

// Fin Add J. Morice communication entre complexe

MPIrank mpiwho(long i) { return MPIrank(i);}
MPIrank mpiwho(long i,fMPI_Comm comm) { return MPIrank(i,comm,Syncro_block);}
MPIrank mpiwho(fMPI_Comm comm,long i) { return MPIrank(i,comm,Syncro_block);}
MPIrank mpiwhob(long i) { return MPIrank(i);}
MPIrank mpiwhob(long i,fMPI_Comm comm) { return MPIrank(i,comm,Syncro_block);}

MPIrank mpiwho_(const long &i,const fMPI_Comm &comm,fMPI_Request * const &rq) { return MPIrank(i,comm,*rq);}
MPIrank mpiwho_(const long &i,fMPI_Request * const &rq) { return MPIrank(i, MPI_COMM_WORLD ,*rq);}

long mpiWait(fMPI_Request * frq) {
 MPI_Request * rq= *frq;
 MPI_Status status;
 long res=MPI_SUCCESS;
  while(rq && *rq!=MPI_REQUEST_NULL)
    {
      res = MPI_Wait(rq,&status);
      DoOnWaitMPIRequest(rq);
    }
    return res;

}

long mpiBarrier(fMPI_Comm * comm)
{
  return MPI_Barrier(*comm);
}

long mpiWaitAny(KN<MPI_Request>* rq)
{
  MPI_Status status;
  int index;
  do {
      MPI_Waitany(rq->N(),*rq,&index,&status);
      if(index != MPI_UNDEFINED)
          DoOnWaitMPIRequest(&(*rq)[index]);
  }  while (MPI_UNDEFINED != index && (*rq)[index] != MPI_REQUEST_NULL);
  return index;
}

long mpiWaitAll(KN<MPI_Request>* rq)
{
  MPI_Status* statuses = new MPI_Status[rq->N()]();
  MPI_Waitall(rq->N(),*rq,statuses);
  for(int i = 0; i < rq->N(); ++i) {
      if(statuses[i].MPI_TAG != MPI_ANY_TAG && statuses[i].MPI_SOURCE != MPI_ANY_SOURCE)
          DoOnWaitMPIRequest(&(*rq)[i]);
  }
  MPI_Waitall(rq->N(),*rq,statuses);
  for(int i = 0; i < rq->N(); ++i) {
      if(statuses[i].MPI_TAG != MPI_ANY_TAG && statuses[i].MPI_SOURCE != MPI_ANY_SOURCE)
          DoOnWaitMPIRequest(&(*rq)[i]);
  }
  delete [] statuses;
  return 0L;
}

MPIrank * set_copympi( MPIrank* const & a,const MPIrank & b){ *a=b;return a;}


long mpiSize(fMPI_Comm  cmm) {
    int s=0;
    if(cmm !=  (MPI_Comm) MPI_COMM_NULL)
    MPI_Comm_size(cmm, &s); /* local */
    return s;
}
long mpiRank(fMPI_Comm  cmm) {
    int s=MPI_UNDEFINED;
    if(cmm != (MPI_Comm)MPI_COMM_NULL)
    MPI_Comm_rank(cmm, &s); /* local */
    return s;
}

AnyType InitializeGroup(Stack stack,const AnyType &x){
    MPI_Group *g=*PGetAny<fMPI_Group>(x);
    *g=MPI_GROUP_NULL;
     MPI_Comm_group(MPI_COMM_WORLD, g);
    return  g;
}
AnyType DeleteGroup(Stack stack,const AnyType &x){
    MPI_Group *g=*PGetAny<fMPI_Group>(x);
    if(g && (*g != MPI_GROUP_NULL))MPI_Group_free(g);
    return  Nothing;
}
AnyType InitializeComm(Stack stack,const AnyType &x){
    MPI_Comm *comm= *PGetAny<fMPI_Comm>(x);
    *comm=MPI_COMM_NULL;
    MPI_Comm_dup(MPI_COMM_WORLD, comm);
    return  comm;
}
AnyType DeleteComm(Stack stack,const AnyType &x){
    MPI_Comm *comm= *PGetAny<fMPI_Comm>(x);
    if(comm && (*comm != MPI_COMM_NULL && *comm != MPI_COMM_WORLD))// add MPI_COMM_WORLD FH 11/2010 FH
      MPI_Comm_free(comm);
    return  Nothing;
}
AnyType InitializeRequest(Stack stack,const AnyType &x){
    MPI_Request *comm=*PGetAny<fMPI_Request>(x);
    *comm=MPI_REQUEST_NULL;

    return  comm;
}
AnyType DeleteRequest(Stack stack,const AnyType &x){
    MPI_Request *comm=*PGetAny<fMPI_Request>(x);
    if(comm && ( *comm!=MPI_REQUEST_NULL )) MPI_Request_free(comm);
    return  Nothing;
}
//  Hack to Bypass a bug in freefem FH  ...
template<>
class ForEachType<MPI_Group>:  public basicForEachType{public:// correction july 2009..... FH  Hoooo....  (Il y a un bug DUR DUR FH  ...)
    ForEachType(Function1 iv=0,Function1 id=0,Function1 OOnReturn=0):basicForEachType(typeid(MPI_Group),sizeof(MPI_Group),0,0,iv,id,OOnReturn) { }
};

template<>
class ForEachType<fMPI_Comm>:  public basicForEachType{public:// coorection july 2009..... FH  Hoooo....  (Il y a un bug DUR DUR FH  ...)
    ForEachType(Function1 iv=0,Function1 id=0,Function1 OOnReturn=0):basicForEachType(typeid(fMPI_Comm),sizeof(fMPI_Comm),0,0,iv,id,OOnReturn) {}
};

template<>
class ForEachType<fMPI_Request>:  public basicForEachType{public:// correction july 2009..... FH  Hoooo....  (Il y a un bug DUR DUR FH  ...)
    ForEachType(Function1 iv=0,Function1 id=0,Function1 OOnReturn=0):basicForEachType(typeid(fMPI_Request),sizeof(fMPI_Request),0,0,iv,id,OOnReturn) {}
};
// end Hack  ...


fMPI_Group* def_group( fMPI_Group* const & a,fMPI_Comm * const &comm, KN_<long>  const & b)
{
  MPI_Group group;
  MPI_Comm_group(*comm,& group);
  KN<int> ranks(b);
  MPI_Group_incl(group, ranks.N(),(int *) ranks, *a);
  MPI_Group_free(&group) ;
  // ici def a ..
  return a;}

fMPI_Group* def_group( fMPI_Group* const & a, KN_<long>  const & b)
{
    MPI_Group group;
    MPI_Comm comm=MPI_COMM_WORLD;
    MPI_Comm_group(comm,& group);
    KN<int> ranks(b);
    MPI_Group_incl(group, ranks.N(),(int *) ranks, *a);
    MPI_Group_free(&group) ;
    // ici def a ..
return a;}

fMPI_Group* def_group( fMPI_Group* const & a,fMPI_Comm * const &comm)
{
  MPI_Comm_group(*comm,*a);
  return a;}

fMPI_Comm* def_comm( fMPI_Comm* const & a,fMPI_Group* const & g)
{
    int ok=MPI_Comm_create(MPI_COMM_WORLD,*g,*a);
    return a;
}

fMPI_Comm* def_comm( fMPI_Comm* const & a,fMPI_Comm* const & b,fMPI_Group* const & g)
{
    MPI_Comm_create(*b,*g,*a);
    return a;
}

fMPI_Group* def_group( fMPI_Group* const & a,fMPI_Group  * const & group,KN_<long>  const & b)
{

  KN<int> ranks(b);
  MPI_Group_incl(*group, ranks.N(), (int *) ranks, *a);
  // ici def a ..
  return a;
}

struct Def_def_Commsplit  : public  quad_function<fMPI_Comm*,fMPI_Comm*,long,long, fMPI_Comm*>
{
    static fMPI_Comm * f(Stack,fMPI_Comm* const & a,fMPI_Comm* const & comm,const long &color ,const long &key  )
    {
	MPI_Comm_split(*comm,color, key, *a);
	return a;
    }
};

fMPI_Comm * mpiCommsplit(fMPI_Comm* const & a,const MPIrank &p1,const long &key  )
{
    MPI_Comm_split(p1.comm, p1.who, key, *a);
    return a;
}


struct Def_def_Intercommcreate  : public  quad_function<fMPI_Comm*,MPIrank,MPIrank,long, fMPI_Comm*>
{
static fMPI_Comm * f(Stack,fMPI_Comm* const & a, MPIrank const & p1, MPIrank const & p2, long const & tag )
{
    int err;
    err=MPI_Intercomm_create(p1.comm, p1.who, p2.comm,p2.who,tag, *a);
    return a;
}
};

fMPI_Comm * def_intercommmerge(fMPI_Comm* const & a,fMPI_Comm* const & b, const long & high)
{
    MPI_Intercomm_merge(*b, high, *a);
    return a;
}


template<typename K,typename KK>
AnyType ClearReturnpKK(Stack stack, const AnyType & a)
{
  // a ne faire que pour les variables local au return...
  //  pour l'instant on copie pour fqire mqrche
  // a repense  FH  mqi 2009....
    KK * m = GetAny<KK * >(a);
    m->increment();
    Add2StackOfPtr2FreeRC(stack,m);
    if(verbosity>400)
	cout << "ClearReturnpKK:: increment + Add2StackOfPtr2FreeRC nb ref  " <<  -m->next  << endl;
    return m;
}

template<typename K,typename KK,typename KK_>
AnyType ClearReturnpKK_(Stack stack, const AnyType & a)
{
  // il faut faire un copie du tableau
  KK_ * m = GetAny<KK_ * >(a);
  KK *cm=new KK(*m);

  Add2StackOfPtr2Free(stack,cm);// detruire la copie
  if(verbosity>400)
    cout << "ClearReturnpKK_:: copie  Add2StackOfPtr2Free "  << endl;
  return (KK_ *) cm;
}
template<typename K,typename KK,typename KK_>
AnyType ClearReturnKK_(Stack stack, const AnyType & a)
{
  // il faut faire un copie du tableau
  KK_  m = GetAny<KK_>(a);
  KK *cm=new KK(m);

  Add2StackOfPtr2Free(stack,cm);// detruire la copie
  if(verbosity>400)
    cout << "ClearReturnKK_:: copie  Add2StackOfPtr2Free   "  << endl;
    return SetAny<KK_>(*cm);
}
fMPI_Request * get_elementp_( KN<MPI_Request> * const & a,const long & b){
  if( a==0 || b<0 || a->N() <= b)
    { if(a) cerr << " Out of bound  0 <=" << b << " < "  << a->N() << " KN<MPI_Request> * " << endl;
      ExecError("Out of bound in operator []");}
  return  reinterpret_cast<fMPI_Request *> (&((*a)[b]));}// bofBof ...

KN<MPI_Request> * set_init0( KN<MPI_Request> * const & a,const long & b)
   {
     a->init(b);
     for(int i=0;i<b;++i)
       (*a)[i]=MPI_REQUEST_NULL;
     return a;
   }
bool toBool(fMPI_Comm *comm)
{
  return (comm && (*comm != (MPI_Comm)MPI_COMM_NULL));
}
void * topVoid(fMPI_Comm *comm) {    return comm; }

template<typename T>
class Quad_Op : public E_F0 {
  typedef typename T::result_type R;
  typedef typename T::first_argument_type A;
  typedef typename T::second_argument_type B;
  typedef typename T::third_argument_type C;
  typedef typename T::fourth_argument_type D;

  typedef  typename T::result_type Result;
  Expression a,b,c,d;
public:
  AnyType operator()(Stack s)  const
  {return  SetAny<R>(static_cast<R>(T::f(s, GetAny<A>((*a)(s)) ,
					 GetAny<B>((*b)(s)) ,
					 GetAny<C>((*c)(s)),
					   GetAny<D>((*d)(s))
					 )));}
  Quad_Op(Expression aa,Expression bb,Expression cc,Expression dd) : a(aa),b(bb),c(cc),d(dd) {}
  bool MeshIndependent() const {
    return a->MeshIndependent() && b->MeshIndependent() && c->MeshIndependent()  && d->MeshIndependent();}
};


// Fin add J. Morice

void f_end_parallele()
{
    /// FFCS: MPI_Finalize() needs to be called later than this (in
    /// ffcs/src/server.cpp)
    ffapi::mpi_finalize();
    if(verbosity> 2 || (verbosity>1&&mpirank ==0)) cout << "FreeFem++-mpi finalize correctly .\n" << flush ;
    else if(verbosity>5)  cout << '.' << endl ;
}
void f_initparallele(int &argc, char **& argv)
{
  /// FFCS: MPI_Init() needs to be called earlier (in ffcs/src/server.cpp)
  ffapi::mpi_init(argc,argv);

  int mpirank1,mpisize1;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank1); /* local */
  MPI_Comm_size(MPI_COMM_WORLD, &mpisize1); /* local */

  mpirank = mpirank1;//MPI::COMM_WORLD.Get_rank();
  mpisize =mpisize1;// MPI::COMM_WORLD.Get_size();
  if(verbosity> 2 || (verbosity>1&&mpirank ==0))
  cout << "initparallele rank " <<  mpirank << " on " << mpisize << endl;
  ff_atend(f_end_parallele); // set end MPI //
}

double ffMPI_Wtime() {return MPI_Wtime();}
double ffMPI_Wtick() {return MPI_Wtick();}

class splitComm_Op : public E_F0mps {
public:
    Expression comm;
    Expression p;
    Expression splitComm;
    static const int n_name_param = 2;
    static basicAC_F0::name_and_type name_param[];
    Expression nargs[n_name_param];
    splitComm_Op(const basicAC_F0& args, Expression param1, Expression param2, Expression param3) : comm(param1), p(param2), splitComm(param3) {
        args.SetNameParam(n_name_param, name_param, nargs);
    }
    AnyType operator()(Stack stack) const;
};

basicAC_F0::name_and_type splitComm_Op::name_param[] = {
    {"topology", &typeid(long)},
    {"exclude", &typeid(bool)}
};

class splitComm : public OneOperator {
public:
    splitComm() : OneOperator(atype<long>(), atype<pcommworld>(), atype<long*>(), atype<pcommworld>()) {}
    E_F0* code(const basicAC_F0& args) const
    {
        return new splitComm_Op(args, t[0]->CastTo(args[0]), t[1]->CastTo(args[1]), t[2]->CastTo(args[2]));
    }
};

static inline bool splitCommunicator(const MPI_Comm& in, MPI_Comm& out, const bool& exclude, unsigned short& p, const unsigned short& T) {
    int size, rank;
    MPI_Comm_size(in, &size);
    MPI_Comm_rank(in, &rank);
    if(p > size / 2 && size > 1) {
        p = size / 2;
        if(rank == 0)
            std::cout << "WARNING -- the number of master processes was set to a value greater than MPI_Comm_size, the value has been reset to " << p << std::endl;
    }
    p = std::max(p, static_cast<unsigned short>(1));
    if(exclude) {
        MPI_Group oldGroup, newGroup;
        MPI_Comm_group(in, &oldGroup);
        int* pm = new int[p];
        if(T == 1)
            for(int i=0;i<p; ++i) pm[i]=i;
        else if(T == 2) {
            float area = size * size / (2.0 * p);
            *pm = 0;
            for(unsigned short i = 1; i < p; ++i)
                pm[i] = static_cast<int>(size - std::sqrt(std::max(size * size - 2 * size * pm[i - 1] - 2 * area + pm[i - 1] * pm[i - 1], 1.0f)) + 0.5);
        }
        else
            for(unsigned short i = 0; i < p; ++i)
                pm[i] = i * (size / p);
        bool excluded = std::binary_search(pm, pm + p, rank);
        if(excluded)
            MPI_Group_incl(oldGroup, p, pm, &newGroup);
        else
            MPI_Group_excl(oldGroup, p, pm, &newGroup);
        MPI_Comm_create(in, newGroup, &out);
        MPI_Group_free(&oldGroup);
        MPI_Group_free(&newGroup);
        delete [] pm;
        return excluded;
    }
    else {
        MPI_Comm_dup(in, &out);
        return false;
    }
}


AnyType splitComm_Op::operator()(Stack stack) const {
    bool exclude = nargs[1] ? GetAny<bool>((*nargs[1])(stack)) : false;
    MPI_Comm* orig_comm = (MPI_Comm*)GetAny<pcommworld>((*comm)(stack));
    MPI_Comm* new_comm = (MPI_Comm*)GetAny<pcommworld>((*splitComm)(stack));
    long* pp = GetAny<long*>((*p)(stack));
    unsigned short p = *pp;
    long topology = nargs[0] ? GetAny<long>((*nargs[0])(stack)) : 0;
    bool excluded = splitCommunicator(*orig_comm, *new_comm, exclude, p, topology);
    *pp = p;
    return static_cast<long>(excluded);
}


void f_init_lgparallele()
  {
    if(verbosity && mpirank == 0) cout << "parallelempi ";
    using namespace Fem2D;
    Dcl_TypeandPtr<MPIrank>(0);

    Dcl_TypeandPtr<fMPI_Group>(0,0,InitializeGroup,DeleteGroup);
    Dcl_TypeandPtr<fMPI_Comm>(0,0,InitializeComm,DeleteComm);
    Dcl_Type<fMPI_Op>();
    Dcl_TypeandPtr<fMPI_Request>(0,0,InitializeRequest,DeleteRequest); // bof bof ...
    Dcl_TypeandPtr_<KN_<MPI_Request> ,KN<MPI_Request>*  >
      (0,0,0,::Destroy<KN<MPI_Request> >,
       ::ClearReturnKK_<MPI_Request,KN<MPI_Request>,KN_<MPI_Request> >,
       ::ClearReturnpKK<MPI_Request,KN<MPI_Request> >);

    zzzfff->Add("mpiGroup",atype<fMPI_Group*>());
    zzzfff->Add("mpiComm",atype<fMPI_Comm*>());
    zzzfff->Add("mpiRequest",atype<fMPI_Request*>());

    map_type_of_map[make_pair(atype<long>(),atype<fMPI_Request>())]=atype<KN<MPI_Request>*>(); // vector


    map_type[typeid(MPIrank).name()]->AddCast(new E_F1_funcT<MPIrank,MPIrank*>(UnRef<MPIrank>));
    map_type[typeid(fMPI_Group).name()]->AddCast(new E_F1_funcT<fMPI_Group,fMPI_Group*>(UnRef<fMPI_Group>));
    map_type[typeid(fMPI_Comm).name()]->AddCast(new E_F1_funcT<fMPI_Comm,fMPI_Comm*>(UnRef<fMPI_Comm>));
    map_type[typeid(bool).name()]->AddCast(new OneOperator1<bool,fMPI_Comm*>(toBool));
    map_type[typeid(void*).name()]->AddCast(new OneOperator1<void *,fMPI_Comm*>(topVoid));


    TheOperators->Add("<-",
		      new OneOperator2_<MPIrank*,MPIrank*,MPIrank>(&set_copympi));

    // constructor example ...
    TheOperators->Add("<-",
		      new OneOperator2_<fMPI_Group*,fMPI_Group*,KN_<long> >(&def_group),
		      new OneOperator3_<fMPI_Group*,fMPI_Group*,fMPI_Group*,KN_<long> >(&def_group),
		      new OneOperator3_<fMPI_Group*,fMPI_Group*,fMPI_Comm*,KN_<long> >(&def_group),
		      new OneOperator2_<fMPI_Group*,fMPI_Group*,fMPI_Comm*>(&def_group));


      TheOperators->Add("<-",
			new OneOperator2_<fMPI_Comm*,fMPI_Comm*,fMPI_Group* >(&def_comm),
			new OneOperator3_<fMPI_Comm*,fMPI_Comm*,fMPI_Comm*,fMPI_Group* >(&def_comm),
			new OneOperator3_<fMPI_Comm*,fMPI_Comm*,MPIrank,long >(&mpiCommsplit),
			new OneOperator3_<fMPI_Comm*,fMPI_Comm*,fMPI_Comm*,long >(&def_intercommmerge),
			new OneQuadOperator< Def_def_Intercommcreate, Quad_Op<Def_def_Intercommcreate>  >,
			new OneQuadOperator< Def_def_Commsplit, Quad_Op<Def_def_Commsplit>  >

			);

    Global.Add("processor","(",new OneOperator1<MPIrank,long>(mpiwho));
    Global.Add("processor","(",new OneOperator2<MPIrank,long,fMPI_Comm>(mpiwho));
      Global.Add("processor","(",new OneOperator2<MPIrank,fMPI_Comm,long>(mpiwho));

    Global.Add("processorblock","(",new OneOperator1<MPIrank,long>(mpiwhob));
    Global.Add("processorblock","(",new OneOperator2<MPIrank,long,fMPI_Comm>(mpiwhob));

    TheOperators->Add(">>",
		      new OneBinaryOperator<Op_Readmpi<double> >,
		      new OneBinaryOperator<Op_Readmpi<Complex> >,
		      new OneBinaryOperator<Op_Readmpi<long> > ,
		      new OneBinaryOperator<Op_Readmpi<KN<double> > > ,
		      new OneBinaryOperator<Op_Readmpi<KN<long> > > ,
		      new OneBinaryOperator<Op_Readmpi<KN<Complex> > > ,
		      new OneBinaryOperator<Op_Readmpi<const Mesh *> > ,
		      new OneBinaryOperator<Op_Readmpi<const Mesh3 *> > ,
              new OneBinaryOperator<Op_Readmpi<const MeshS *> > ,
              new OneBinaryOperator<Op_Readmpi<const MeshL *> > ,
              new OneBinaryOperator<Op_Readmpi<Matrice_Creuse<R> > > ,
		      new OneBinaryOperator<Op_Readmpi<Matrice_Creuse<Complex> > >
		      );

   TheOperators->Add(">>",
                     new OneBinaryOperator<Op_Readmpi<KNM<double> > > ,
                     new OneBinaryOperator<Op_Readmpi<KNM<long> > > ,
                     new OneBinaryOperator<Op_Readmpi<KNM<Complex> > > );

     TheOperators->Add("<<",
		       new OneBinaryOperator<Op_Writempi<double> >,
		       new OneBinaryOperator<Op_Writempi<Complex> >,
		       new OneBinaryOperator<Op_Writempi<long> > ,
		       new OneBinaryOperator<Op_Writempi<KN<double> * > > ,
		       new OneBinaryOperator<Op_Writempi<KN<long> * > > ,
		       new OneBinaryOperator<Op_Writempi<KN<Complex> * > > ,
		       new OneBinaryOperator<Op_Writempi<const Mesh *> > ,
		       new OneBinaryOperator<Op_Writempi<const Mesh3 *> > ,
               new OneBinaryOperator<Op_Writempi<const MeshS *> > ,
               new OneBinaryOperator<Op_Writempi<const MeshL *> > ,
		       new OneBinaryOperator<Op_Writempi<Matrice_Creuse<R> * > > ,
		       new OneBinaryOperator<Op_Writempi<Matrice_Creuse<Complex>* > >

		       );

      TheOperators->Add("<<",
                        new OneBinaryOperator<Op_Writempi<KNM<double> * > > ,
                        new OneBinaryOperator<Op_Writempi<KNM<long> * > > ,
                        new OneBinaryOperator<Op_Writempi<KNM<Complex> * > > );

     Global.Add("Send","(", new OneBinaryOperator<Op_Sendmpi<double> >);
     Global.Add("Send","(", new OneBinaryOperator<Op_Sendmpi<Complex> >);
     Global.Add("Send","(", new OneBinaryOperator<Op_Sendmpi<long> >);
     Global.Add("Send","(", new OneBinaryOperator<Op_Sendmpi<KN<double> *> >);
     Global.Add("Send","(", new OneBinaryOperator<Op_Sendmpi<KN<long> *> >);
     Global.Add("Send","(", new OneBinaryOperator<Op_Sendmpi<KN<Complex> *> >);
     Global.Add("Send","(", new OneBinaryOperator<Op_Sendmpi<KNM<double> *> >);
     Global.Add("Send","(", new OneBinaryOperator<Op_Sendmpi<KNM<long> *> >);
     Global.Add("Send","(", new OneBinaryOperator<Op_Sendmpi<KNM<Complex> *> >);
     Global.Add("Send","(", new OneBinaryOperator<Op_Sendmpi<const Mesh *> >);
     Global.Add("Send","(", new OneBinaryOperator<Op_Sendmpi<const Mesh3 *> >);
     Global.Add("Send","(", new OneBinaryOperator<Op_Sendmpi<const MeshS *> >);
     Global.Add("Send","(", new OneBinaryOperator<Op_Sendmpi<const MeshL *> >);
     Global.Add("Send","(", new OneBinaryOperator<Op_Sendmpi<Matrice_Creuse<R> *> >);
     Global.Add("Send","(", new OneBinaryOperator<Op_Sendmpi<Matrice_Creuse<Complex> *> >);

     Global.Add("Isend","(", new OneBinaryOperator<Op_ISendmpi<double> >);
     Global.Add("Isend","(", new OneBinaryOperator<Op_ISendmpi<long> >);
     Global.Add("Isend","(", new OneBinaryOperator<Op_ISendmpi<Complex> >);
     Global.Add("Isend","(", new OneBinaryOperator<Op_ISendmpi<KN<double> *> >);
     Global.Add("Isend","(", new OneBinaryOperator<Op_ISendmpi<KN<long> *> >);
     Global.Add("Isend","(", new OneBinaryOperator<Op_ISendmpi<KN<Complex> *> >);
     Global.Add("Isend","(", new OneBinaryOperator<Op_ISendmpi<KNM<double> *> >);
     Global.Add("Isend","(", new OneBinaryOperator<Op_ISendmpi<KNM<long> *> >);
     Global.Add("Isend","(", new OneBinaryOperator<Op_ISendmpi<KNM<Complex> *> >);
     Global.Add("Isend","(", new OneBinaryOperator<Op_ISendmpi<const Mesh *> >);
     Global.Add("Isend","(", new OneBinaryOperator<Op_ISendmpi<const Mesh3 *> >);
     Global.Add("Isend","(", new OneBinaryOperator<Op_ISendmpi<const MeshS *> >);
     Global.Add("Isend","(", new OneBinaryOperator<Op_ISendmpi<const MeshL *> >);
     Global.Add("Isend","(", new OneBinaryOperator<Op_ISendmpi<Matrice_Creuse<R> *> >);
     Global.Add("Isend","(", new OneBinaryOperator<Op_ISendmpi<Matrice_Creuse<Complex> *> >);

      Global.Add("Recv","(", new OneBinaryOperator<Op_Recvmpi<double> >);
      Global.Add("Recv","(", new OneBinaryOperator<Op_Recvmpi<long> >);
      Global.Add("Recv","(", new OneBinaryOperator<Op_Recvmpi<Complex> >);
      Global.Add("Recv","(", new OneBinaryOperator<Op_Recvmpi<KN<double> > >);
      Global.Add("Recv","(", new OneBinaryOperator<Op_Recvmpi<KN<long> > >);
      Global.Add("Recv","(", new OneBinaryOperator<Op_Recvmpi<KN<Complex> > >);
      Global.Add("Recv","(", new OneBinaryOperator<Op_Recvmpi<KNM<double> > >);
      Global.Add("Recv","(", new OneBinaryOperator<Op_Recvmpi<KNM<long> > >);
      Global.Add("Recv","(", new OneBinaryOperator<Op_Recvmpi<KNM<Complex> > >);
      Global.Add("Recv","(", new OneBinaryOperator<Op_Recvmpi<const Mesh *> >);
      Global.Add("Recv","(", new OneBinaryOperator<Op_Recvmpi<const Mesh3 *> >);
      Global.Add("Recv","(", new OneBinaryOperator<Op_Recvmpi<const MeshS *> >);
      Global.Add("Recv","(", new OneBinaryOperator<Op_Recvmpi<const MeshL *> >);
      Global.Add("Recv","(", new OneBinaryOperator<Op_Recvmpi<Matrice_Creuse<R> > >);
      Global.Add("Recv","(", new OneBinaryOperator<Op_Recvmpi<Matrice_Creuse<Complex> > >);

      Global.Add("Irecv","(", new OneBinaryOperator<Op_IRecvmpi<double> >);
      Global.Add("Irecv","(", new OneBinaryOperator<Op_IRecvmpi<long> >);
      Global.Add("Irecv","(", new OneBinaryOperator<Op_IRecvmpi<Complex> >);
      Global.Add("Irecv","(", new OneBinaryOperator<Op_IRecvmpi<KN<double> > >);
      Global.Add("Irecv","(", new OneBinaryOperator<Op_IRecvmpi<KN<long> > >);
      Global.Add("Irecv","(", new OneBinaryOperator<Op_IRecvmpi<KN<Complex> > >);
      Global.Add("Irecv","(", new OneBinaryOperator<Op_IRecvmpi<KNM<double> > >);
      Global.Add("Irecv","(", new OneBinaryOperator<Op_IRecvmpi<KNM<long> > >);
      Global.Add("Irecv","(", new OneBinaryOperator<Op_IRecvmpi<KNM<Complex> > >);
      Global.Add("Irecv","(", new OneBinaryOperator<Op_IRecvmpi<const Mesh *> >);
      Global.Add("Irecv","(", new OneBinaryOperator<Op_IRecvmpi<const Mesh3 *> >);
      Global.Add("Irecv","(", new OneBinaryOperator<Op_IRecvmpi<const MeshS *> >);
      Global.Add("Irecv","(", new OneBinaryOperator<Op_IRecvmpi<const MeshL *> >);
      Global.Add("Irecv","(", new OneBinaryOperator<Op_IRecvmpi<Matrice_Creuse<R> > >);
      Global.Add("Irecv","(", new OneBinaryOperator<Op_IRecvmpi<Matrice_Creuse<Complex> > >);




      Global.Add("broadcast","(",new OneBinaryOperator<Op_Bcastmpi<double> >);
      Global.Add("broadcast","(",new OneBinaryOperator<Op_Bcastmpi<Complex> >);
      Global.Add("broadcast","(",new OneBinaryOperator<Op_Bcastmpi<long> >);
      Global.Add("broadcast","(",new OneBinaryOperator<Op_Bcastmpi<KN<double> > >);
      Global.Add("broadcast","(",new OneBinaryOperator<Op_Bcastmpi<KN<long> > >);
      Global.Add("broadcast","(",new OneBinaryOperator<Op_Bcastmpi<KN<Complex> > >);
      Global.Add("broadcast","(",new OneBinaryOperator<Op_Bcastmpi<KNM<double> > >);
      Global.Add("broadcast","(",new OneBinaryOperator<Op_Bcastmpi<KNM<long> > >);
      Global.Add("broadcast","(",new OneBinaryOperator<Op_Bcastmpi<KNM<Complex> > >);
      Global.Add("broadcast","(",new OneBinaryOperator<Op_Bcastmpi<const Mesh *> >);
      Global.Add("broadcast","(",new OneBinaryOperator<Op_Bcastmpi<const Mesh3 *> >);
      Global.Add("broadcast","(",new OneBinaryOperator<Op_Bcastmpi<const MeshS *> >);
      Global.Add("broadcast","(",new OneBinaryOperator<Op_Bcastmpi<const MeshL *> >);
      Global.Add("broadcast","(",new OneBinaryOperator<Op_Bcastmpi<Matrice_Creuse<R> > >);
      Global.Add("broadcast","(",new OneBinaryOperator<Op_Bcastmpi<Matrice_Creuse<Complex> > >);

      Global.Add("mpiAlltoall","(",new OneBinaryOperator<Op_All2All< long > >);
      Global.Add("mpiAlltoall","(",new OneBinaryOperator<Op_All2All< double > >);
      Global.Add("mpiAllgather","(",new OneBinaryOperator<Op_Allgather< long > >);
      Global.Add("mpiAllgather","(",new OneBinaryOperator<Op_Allgather< double > >);
      Global.Add("mpiAlltoall","(",new OneTernaryOperator3<Op_All2All3< long > >);
      Global.Add("mpiAlltoall","(",new OneTernaryOperator3<Op_All2All3< double > >);
      Global.Add("mpiAllgather","(",new OneTernaryOperator3<Op_Allgather3< long > >);
      Global.Add("mpiAllgather","(",new OneTernaryOperator3<Op_Allgather3< double > >);

      Global.Add("mpiAllgather","(",new OneBinaryOperator<Op_Allgather1< long > >);   // Add J. Morice
      Global.Add("mpiAllgather","(",new OneBinaryOperator<Op_Allgather1< double > >); // Add J. Morice

      Global.Add("mpiAllgather","(",new OneTernaryOperator3<Op_Allgather13< long > >);  // Add J. Morice
      Global.Add("mpiAllgather","(",new OneTernaryOperator3<Op_Allgather13< double > >);// Add J. Morice

      Global.Add("mpiScatter","(",new OneTernaryOperator3<Op_Scatter1< long > >);    // Add J. Morice
      Global.Add("mpiScatter","(",new OneTernaryOperator3<Op_Scatter1< double > >);  // Add J. Morice
      Global.Add("mpiScatter","(",new OneTernaryOperator3<Op_Scatter3< long > >);
      Global.Add("mpiScatter","(",new OneTernaryOperator3<Op_Scatter3< double > >);

      Global.Add("mpiGather","(",new OneTernaryOperator3<Op_Gather1< long > >);   // Add J. Morice
      Global.Add("mpiGather","(",new OneTernaryOperator3<Op_Gather1< double > >); // Add J. Morice
      Global.Add("mpiGather","(",new OneTernaryOperator3<Op_Gather3< long > >); // correction J. Morice Scatter --> Gather
      Global.Add("mpiGather","(",new OneTernaryOperator3<Op_Gather3< double > >);

      // Add J. Morice communication with vector of different size
      Global.Add("mpiAlltoallv","(",new OneOperator6_<long, KN_<double>, KN_<double>, KN_<long>, KN_<long>, KN_<long>, KN_<long> >( Op_All2Allv<double> ) );
      Global.Add("mpiAlltoallv","(",new OneOperator6_<long, KN_<long>, KN_<long>, KN_<long>, KN_<long>, KN_<long>, KN_<long> >( Op_All2Allv<long> ) );

      Global.Add("mpiAlltoallv","(",new OneOperator7_<long, KN_<long>, KN_<long>, fMPI_Comm, KN_<long>, KN_<long>, KN_<long>, KN_<long> >( Op_All2All3v<long> ) );
      Global.Add("mpiAlltoallv","(",new OneOperator7_<long, KN_<double>, KN_<double>, fMPI_Comm, KN_<long>, KN_<long>, KN_<long>, KN_<long> >( Op_All2All3v<double> ) );

      Global.Add("mpiAllgatherv","(",new OneQuadOperator<Op_Allgatherv< long >, Quad_Op<Op_Allgatherv< long > > > );
      Global.Add("mpiAllgatherv","(",new OneQuadOperator<Op_Allgatherv< double >, Quad_Op<Op_Allgatherv< double > > >);
      Global.Add("mpiAllgatherv","(",new OneOperator5_<long, KN_<long>, KN_<long>, fMPI_Comm, KN_<long>, KN_<long> >(Op_Allgatherv3< long >) );
      Global.Add("mpiAllgatherv","(",new OneOperator5_<long, KN_<double>, KN_<double>, fMPI_Comm, KN_<long>, KN_<long> >(Op_Allgatherv3< double >) );

      Global.Add("mpiScatterv","(",new OneOperator5_<long, KN_<long>, KN_<long>, MPIrank, KN_<long>, KN_<long> >(Op_Scatterv3< long >) );
      Global.Add("mpiScatterv","(",new OneOperator5_<long, KN_<double>, KN_<double>, MPIrank, KN_<long>, KN_<long> >(Op_Scatterv3< double >) );

      Global.Add("mpiGatherv","(",new OneOperator5_<long, KN_<long>, KN_<long>, MPIrank, KN_<long>, KN_<long> >( Op_Gatherv3< long > ) );
      Global.Add("mpiGatherv","(",new OneOperator5_<long, KN_<double>, KN_<double>, MPIrank, KN_<long>, KN_<long> >( Op_Gatherv3< double > ) );
      // Fin Add J. Morice

      Global.Add("mpiReduce","(",new OneQuadOperator<Op_Reduce< double >, Quad_Op<Op_Reduce< double > > >);
      Global.Add("mpiReduce","(",new OneQuadOperator<Op_Reduce1< double >, Quad_Op<Op_Reduce1< double > > >);
      Global.Add("mpiAllReduce","(",new OneQuadOperator<Op_AllReduce< double >, Quad_Op<Op_AllReduce< double > > >);
      Global.Add("mpiAllReduce","(",new OneQuadOperator<Op_AllReduce1< double >, Quad_Op<Op_AllReduce1< double > > >); // add FH jan 2011

      // Add J. Morice
      Global.Add("mpiReduce","(",new OneQuadOperator<Op_Reduce< long >, Quad_Op<Op_Reduce< long > > >);
      Global.Add("mpiReduce","(",new OneQuadOperator<Op_Reduce1< long >, Quad_Op<Op_Reduce1< long > > >);
      Global.Add("mpiAllReduce","(",new OneQuadOperator<Op_AllReduce< long >, Quad_Op<Op_AllReduce< long > > >);
      Global.Add("mpiAllReduce","(",new OneQuadOperator<Op_AllReduce1< long >, Quad_Op<Op_AllReduce1< long > > >); // add FH jan 2011
      //    Global.Add("mpiReduceScatter","(",new OneQuadOperator<Op_Reducescatter< long >, Quad_Op<Op_Reducescatter< long > > >);
      // fin Add J. Morice

      // Add J. Morice :: complex communication between processor
      Global.Add("mpiAlltoall","(",new OneBinaryOperator<Op_All2All< Complex > >);
      Global.Add("mpiAllgather","(",new OneBinaryOperator<Op_Allgather< Complex > >);
      Global.Add("mpiAlltoall","(",new OneTernaryOperator3<Op_All2All3< Complex > >);
      Global.Add("mpiAllgather","(",new OneTernaryOperator3<Op_Allgather3< Complex > >);

      Global.Add("mpiAllgather","(",new OneBinaryOperator<Op_Allgather1< Complex > >);
      Global.Add("mpiAllgather","(",new OneTernaryOperator3<Op_Allgather13< Complex > >);

      Global.Add("mpiScatter","(",new OneTernaryOperator3<Op_Scatter3< Complex > >);
      Global.Add("mpiGather","(",new OneTernaryOperator3<Op_Gather3< Complex > >);

      Global.Add("mpiScatter","(",new OneTernaryOperator3<Op_Scatter1< Complex > >);
      Global.Add("mpiGather","(",new OneTernaryOperator3<Op_Gather1< Complex > >);

      // Add J. Morice communication with vector of different size
      Global.Add("mpiAlltoallv","(",new OneOperator6_<long, KN_<Complex>, KN_<Complex>, KN_<long>, KN_<long>, KN_<long>, KN_<long> >( Op_All2Allv<Complex> ) );
      Global.Add("mpiAlltoallv","(",new OneOperator7_<long, KN_<Complex>, KN_<Complex>, fMPI_Comm, KN_<long>, KN_<long>, KN_<long>, KN_<long> >( Op_All2All3v<Complex> ) );
      Global.Add("mpiAllgatherv","(",new OneQuadOperator<Op_Allgatherv< Complex>, Quad_Op<Op_Allgatherv< Complex > > >);

      Global.Add("mpiAllgatherv","(",new OneOperator5_<long, KN_<Complex>, KN_<Complex>, fMPI_Comm, KN_<long>, KN_<long> >(Op_Allgatherv3< Complex >) );
      Global.Add("mpiScatterv","(",new OneOperator5_<long, KN_<Complex>, KN_<Complex>, MPIrank, KN_<long>, KN_<long> >(Op_Scatterv3< Complex >) );
      Global.Add("mpiGatherv","(",new OneOperator5_<long, KN_<Complex>, KN_<Complex>, MPIrank, KN_<long>, KN_<long> >( Op_Gatherv3< Complex > ) );



      Global.Add("mpiReduce","(",new OneQuadOperator<Op_ReduceMat< Complex >, Quad_Op<Op_ReduceMat< Complex > > >);// add FH april 2011
      Global.Add("mpiReduce","(",new OneQuadOperator<Op_ReduceMat< double >, Quad_Op<Op_ReduceMat< double > > >);// add FH april 2011
      Global.Add("mpiAllReduce","(",new OneQuadOperator<Op_AllReduceMat< Complex >, Quad_Op<Op_AllReduceMat< Complex > > >);// add FH april 2011
      Global.Add("mpiAllReduce","(",new OneQuadOperator<Op_AllReduceMat< double >, Quad_Op<Op_AllReduceMat< double > > >);// add FH april 2011

      Global.Add("mpiReduce","(",new OneQuadOperator<Op_Reduce< Complex >, Quad_Op<Op_Reduce< Complex > > >);
      Global.Add("mpiReduce","(",new OneQuadOperator<Op_Reduce1< Complex >, Quad_Op<Op_Reduce1< Complex > > >);
      Global.Add("mpiAllReduce","(",new OneQuadOperator<Op_AllReduce< Complex >, Quad_Op<Op_AllReduce< Complex > > >);
#ifdef HAVE_MPI_DOUBLE_COMPLEX
    Global.Add("mpiAllReduce","(",new OneQuadOperator<Op_AllReduce1< Complex >, Quad_Op<Op_AllReduce1< Complex > > >);// add FH jan 2011
#endif
      // Fin Add J. Morice :: complex communication between processor

      Global.New("mpirank",CConstant<long>(mpirank));
      Global.New("mpisize",CConstant<long>(mpisize));
     static long mpiUndefined=MPI_UNDEFINED, mpiAnySource =  MPI_ANY_SOURCE,mpiAnyTag=MPI_ANY_TAG ;
     static fMPI_Comm fmpiWorld=MPI_COMM_WORLD;
     static fMPI_Comm fmpiSelf=MPI_COMM_SELF;

     Global.New("mpiUndefined",CConstant<long>(mpiUndefined));
     Global.New("mpiAnySource",CConstant<long>(mpiAnySource));
      Global.New("mpiAnyTag",CConstant<long>(mpiAnyTag));



      Global.New("mpiCommWorld",CConstant<fMPI_Comm*>(&fmpiWorld));
      Global.New("mpiCommSelf",CConstant<fMPI_Comm*>(&fmpiSelf));
      // add FH


      Global.Add("mpiWtime","(",new OneOperator0<double>(ffMPI_Wtime));
      Global.Add("mpiWtick","(",new OneOperator0<double>(ffMPI_Wtick));
      Global.Add("processor","(",new OneOperator3_<MPIrank,long,fMPI_Comm,fMPI_Request*>(mpiwho_));
      Global.Add("processor","(",new OneOperator2_<MPIrank,long,fMPI_Request*>(mpiwho_));
      Global.Add("mpiWait","(",new OneOperator1<long,fMPI_Request*>(mpiWait));
      Global.Add("mpiWaitAny","(",new OneOperator1<long,KN<MPI_Request>*>(mpiWaitAny));
      Global.Add("mpiWaitAll","(",new OneOperator1<long,KN<MPI_Request>*>(mpiWaitAll));
      Global.Add("mpiSize","(",new OneOperator1<long,fMPI_Comm>(mpiSize));
      Global.Add("mpiRank","(",new OneOperator1<long,fMPI_Comm>(mpiRank));
      Global.Add("mpiBarrier","(",new OneOperator1<long,fMPI_Comm*>(mpiBarrier));

      static   fMPI_Op op_max(MPI_MAX);
      static   fMPI_Op op_min(MPI_MIN);
      static   fMPI_Op op_sum(MPI_SUM);
      static   fMPI_Op op_prod(MPI_PROD);
      static   fMPI_Op op_land(MPI_LAND);

      static   fMPI_Op op_lor(MPI_LOR);
      static   fMPI_Op op_lxor(MPI_LXOR);
      static   fMPI_Op op_band(MPI_BAND);
      static   fMPI_Op op_bor(MPI_BOR);
      static   fMPI_Op op_bxor(MPI_BXOR);
      static   fMPI_Op op_maxloc(MPI_MAXLOC);
      static   fMPI_Op op_minloc(MPI_MINLOC);

      Global.New("mpiMAX",CConstant<fMPI_Op>(op_max));
      Global.New("mpiMIN",CConstant<fMPI_Op>(op_min));
      Global.New("mpiSUM",CConstant<fMPI_Op>(op_sum));
      Global.New("mpiPROD",CConstant<fMPI_Op>(op_prod));
      Global.New("mpiLAND",CConstant<fMPI_Op>(op_land));

      Global.New("mpiLOR",CConstant<fMPI_Op>(op_lor));
      Global.New("mpiLXOR",CConstant<fMPI_Op>(op_lxor));
      Global.New("mpiBAND",CConstant<fMPI_Op>(op_band));
      Global.New("mpiBXOR",CConstant<fMPI_Op>(op_bxor));
      //  sur des pair bof bof ...
      Global.New("mpiMAXLOC",CConstant<fMPI_Op>(op_maxloc));
      Global.New("mpiMINLOC",CConstant<fMPI_Op>(op_minloc));


      TheOperators->Add("<-",
			new OneOperator2_<KN<MPI_Request> *,KN<MPI_Request> *,long>(&set_init0)
			);
      atype<KN<MPI_Request>* >()->Add("[","",new OneOperator2_<fMPI_Request*,KN<MPI_Request>*,long >(get_elementp_));

       Global.Add("splitComm", "(", new splitComm);

  }


// set the 3 ptr
extern void (*initparallele)(int &argc, char **& argv) ;
extern void (*init_lgparallele)();


void init_ptr_parallelepmi();
void init_ptr_parallelepmi(){
initparallele=&f_initparallele ;
init_lgparallele=&f_init_lgparallele;
};
