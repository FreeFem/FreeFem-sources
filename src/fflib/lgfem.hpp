// -*- Mode : c++ -*-
//
// SUMMARY  :
// USAGE    :
// ORG      : pfes3_tefk
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
// ATTENTION pfes est la classe qui modelise un pointeur
// sur un FESpace (donc un espace d'element fini
//  mais si la maillage change alors
//  l'espace change
//   cf la fonction set qui reconstruit FESpace

#ifndef lgfem_hpp_
#define lgfem_hpp_
#include "array_resize.hpp"

extern Block *currentblock;

void init_lgmat( );    // initialisation for sparse mat functionnallity

class generic_v_fes;       // J. Morice : Added in 06/22
class vect_generic_v_fes;  // J. Morice : Added in 06/22
class v_fes;
class v_fes3;
class v_fesS;
class v_fesL;

typedef generic_v_fes *pgenericfes;           // J. Morice : Added in 06/22
typedef vect_generic_v_fes *pvectgenericfes;  // J. Morice : Added in 06/22

typedef v_fes *pfes;
typedef v_fes3 *pfes3;
typedef v_fesS *pfesS;
typedef v_fesL *pfesL;

//typedef v_fesVect *pfesVect; // J. Morice

namespace Fem2D {
  class Mesh3;

}
using Fem2D::Mesh;
using Fem2D::Mesh3;
using Fem2D::MeshL;
using Fem2D::MeshS;

typedef const Mesh *pmesh;
typedef const Mesh3 *pmesh3;
typedef const MeshS *pmeshS;
typedef const MeshL *pmeshL;

using Fem2D::FESpace;
using Fem2D::R;
using Fem2D::TypeOfFE;
//  typedef double R;
namespace {
  using namespace Fem2D;
  using Fem2D::Vertex;

  class lgVertex {
   public:
    typedef double R;
    CountPointer< const Mesh > pTh;
    Vertex *v;
    void Check( ) const {
      if (!v || !pTh) {
        ExecError("Too bad! Unset Vertex!");
      }
    }
    void init( ) {
      v = 0;
      pTh.init( );
    }
    lgVertex(const Mesh *Th, long kk) : pTh(Th), v(&(*pTh)(kk)) {}
    lgVertex(const Mesh *Th, Vertex *kk) : pTh(Th), v(kk) {}
    operator int( ) const {
      Check( );
      return (*pTh)(v);
    }
    operator R2 *( ) {
      Check( );
      return v;
    }
    R x( ) const {
      Check( );
      return v->x;
    }
    R y( ) const {
      Check( );
      return v->y;
    }
    //  R z() const {Check() ; return v->z;}
    long lab( ) const {
      Check( );
      return v->lab;
    }
    void destroy( ) { pTh.destroy( ); }
  };

  class lgElement {
   public:
    struct Adj {    //
      const Mesh *pTh;
      Triangle *k;
      Adj(const lgElement &pp) : pTh(pp.pTh), k(pp.k) {}
      lgElement adj(long &e) const {
        int ee;
        ffassert(pTh && k && e >= 0 && e < 3);
        long kk = pTh->ElementAdj((*pTh)(k), ee = e);
        e = ee;
        return lgElement(pTh, kk);
      }
    };
    CountPointer< const Mesh > pTh;
    Triangle *k;

    lgElement( ) : k(0) {}
    void Check( ) const {
      if (!k || !pTh) {
        ExecError("Unset Triangle,Sorry!");
      }
    }
    void init( ) {
      k = 0;
      pTh.init( );
    }
    void destroy( ) { pTh.destroy( ); }
    lgElement(const Mesh *Th, long kk) : pTh(Th), k(&(*Th)[kk]) {}
    lgElement(const Mesh *Th, Triangle *kk) : pTh(Th), k(kk) {}
    operator int( ) const {
      Check( );
      return (*pTh)(k);
    }
    lgVertex operator[](const long &i) const {
      Check( );
      return lgVertex(pTh, &(*k)[i]);
    }
    long lab( ) const {
      Check( );
      return k ? k->lab : LONG_MIN;
    }
    double area( ) const {
      Check( );
      return k->area;
    }
    long n( ) const { return k ? 3 : 0; }
    bool operator==(const lgElement &l) const { return pTh == l.pTh && k == l.k; }
    bool operator!=(const lgElement &l) const { return pTh != l.pTh || k != l.k; }
    bool operator<(const lgElement &l) const { return pTh == l.pTh && k < l.k; }
    bool operator<=(const lgElement &l) const { return pTh == l.pTh && k <= l.k; }
  };
  // add FH  August 2009 ...

  class lgBoundaryEdge {
   public:
    struct BE {
      const Mesh *p;
      BE(const Mesh *pp) : p(pp) {}
      BE(const Mesh **pp) : p(*pp) {}
      operator const Mesh *( ) const { return p; }
    };

    CountPointer< const Mesh > pTh;
    BoundaryEdge *k;

    lgBoundaryEdge( ) : k(0) {}
    void Check( ) const {
      if (!k || !pTh) {
        ExecError("Unset BoundaryEdge,Sorry!");
      }
    }
    void init( ) {
      k = 0;
      pTh.init( );
    }
    void destroy( ) { pTh.destroy( ); }
    lgBoundaryEdge(const Mesh *Th, long kk) : pTh(Th), k(&(*pTh).be(kk)) {}
    lgBoundaryEdge(const Mesh *Th, BoundaryEdge *kk) : pTh(Th), k(kk) {}
    lgBoundaryEdge(const BE &be, long kk) : pTh(be.p), k(&(*pTh).be(kk)) {}
    lgBoundaryEdge(const BE &be, BoundaryEdge *kk) : pTh(be.p), k(kk) {}
    operator int( ) const {
      Check( );
      return (*pTh)(k);
    }
    lgVertex operator[](const long &i) const {
      Check( );
      return lgVertex(pTh, &(*k)[i]);
    }
    long lab( ) const {
      Check( );
      return k ? k->lab : 0;
    }
    double length( ) const {
      Check( );
      return k->length( );
    }
    long n( ) const { return k ? 2 : 0; }
    lgElement Element( ) const {
      Check( );
      int ee;
      return lgElement(pTh, (*pTh).BoundaryElement((*pTh)(k), ee));
    }
    long EdgeElement( ) const {
      Check( );
      int ee;
      (*pTh).BoundaryElement((*pTh)(k), ee);
      return ee;
    }
    R2 NBoundaryElement() const {Check() ;int ee;
          int kk=(*pTh).BoundaryElement((*pTh)(k),ee);
          return (*pTh)[kk].n(ee);}

  };

}    // namespace

void GetPeriodic(const int d, Expression perio, int &nbcperiodic, Expression *&periodic);

bool BuildPeriodic2(int nbcperiodic, Expression *periodic, const Mesh &Th, Stack stack, int &nbdfv,
                   KN< int > &ndfv, int &nbdfe, KN< int > &ndfe);

// <<v_fes>> uses [[file:~/ff/src/femlib/RefCounter.hpp::RefCounter]]

// virtual class for generic v_fes
// Morice :: 06/2022
class generic_v_fes : public RefCounter {
  typedef void* pfes;
  typedef void* FESpace;

  public:
  virtual bool differentMesh(){ ffassert(0); return true; }; 
  
  // declaration avec les virtuelles pure
  // function for vect_generic_v_fes
  virtual void *getpVh()=0;         // get the pointer of FESpace
  virtual void *getpVh()const=0;         // get the pointer of FESpace
  virtual const void *getppTh()=0;  // get the pointer on the pointer of the mesh
  virtual int getN()=0;             // get the number of item of the FESpace
};

// 2d
class v_fes : public generic_v_fes{ //: public RefCounter {
 public:
  typedef ::pfes pfes;
  typedef ::FESpace FESpace;

  static const int dHat = 2, d = 2;
  const int N;
  const pmesh *ppTh;    // adr du maillage
  CountPointer< FESpace > pVh;

  Stack stack;    // the stack is use whith periodique expression

  int nbcperiodic;
  Expression *periodic;

  operator FESpace *( ) {
    throwassert(dHat == 2);
    if (!pVh || *ppTh != &pVh->Th) pVh = CountPointer< FESpace >(update( ), true);
    return pVh;
  }

  FESpace *update( );

  v_fes(int NN, const pmesh *t, Stack s, int n, Expression *p)
    : N(NN), ppTh(t), pVh(0), stack(s), nbcperiodic(n), periodic(p) {}
  v_fes(int NN, const v_fes *f, Stack s, int n, Expression *p)
    : N(NN), ppTh(f->ppTh), pVh(0), stack(s), nbcperiodic(n), periodic(p) {}

  void destroy( ) {
    ppTh = 0;
    pVh = 0;
    delete this;
  }
  virtual ~v_fes( ) {}
  bool buildperiodic(Stack stack, int &nbdfv, KN< int > &ndfv, int &nbdfe, KN< int > &ndfe);
  virtual FESpace *buildupdate(int &nbdfv, KN< int > &ndfv, int &nbdfe, KN< int > &ndfe) = 0;
  virtual FESpace *buildupdate( ) = 0;
  
  // Added to return the pointer to FESpace 
  const void *getppTh(){ return ppTh; }; // Morice : get the pointer on the pointer of the mesh
  void *getpVh(){ return pVh; };         // Morice : get the pointer of FESpace
  void *getpVh()const{ return pVh; };         // Morice : get the pointer of FESpace
  int getN(){ return N;};                // Morice : get the number of item of the FESpace
};

// 3D volume
class v_fes3 : public generic_v_fes { //: public RefCounter {
 public:
  typedef pfes3 pfes;
  typedef FESpace3 FESpace;

  static const int dHat = 3, d = 3;
  const int N;
  const pmesh3 *ppTh;    // adr du maillage
  CountPointer< FESpace3 > pVh;

  Stack stack;    // the stack is use whith periodique expression

  int nbcperiodic;
  Expression *periodic;

  operator FESpace3 *( ) {
    throwassert(dHat == 3);
    if (!pVh || *ppTh != &pVh->Th) pVh = CountPointer< FESpace3 >(update( ), true);
    return pVh;
  }
  FESpace3 *update( );

  v_fes3(int NN, const pmesh3 *t, Stack s, int n, Expression *p)
    : N(NN), ppTh(t), pVh(0), stack(s), nbcperiodic(n), periodic(p) {}
  v_fes3(int NN, const v_fes3 *f, Stack s, int n, Expression *p)
    : N(NN), ppTh(f->ppTh), pVh(0), stack(s), nbcperiodic(n), periodic(p) {}

  // void destroy(){ ppTh=0;pVh=0; delete this;}
  virtual ~v_fes3( ) {}
  bool buildperiodic(Stack stack, KN< int > &ndfe);
  virtual FESpace3 *buildupdate(KN< int > &ndfe) { return 0; }
  virtual FESpace3 *buildupdate( ) { return 0; };

  const void * getppTh(){ return ppTh; }; // Morice : get the pointer on the pointer of the mesh
  void * getpVh(){ return pVh; };         // Morice : get the pointer of FESpace
  void * getpVh()const{ return pVh; };         // Morice : get the pointer of FESpace
  int getN(){ return N;};                // Morice : get the number of item of the FESpace
};

// 3D surface
class v_fesS : public generic_v_fes { //: public RefCounter {
 public:
  typedef pfesS pfes;
  typedef FESpaceS FESpace;

  static const int dHat = 2, d = 3;
  const int N;
  const pmeshS *ppTh;    // adr du maillage

  CountPointer< FESpaceS > pVh;

  Stack stack;    // the stack is use whith periodique expression

  int nbcperiodic;
  Expression *periodic;

  operator FESpaceS *( ) {
    throwassert(dHat == 2);
    if (!pVh || *ppTh != &pVh->Th) pVh = CountPointer< FESpaceS >(update( ), true);
    return pVh;
  }
  FESpaceS *update( );

  v_fesS(int NN, const pmeshS *t, Stack s, int n, Expression *p)    /// TODO
    : N(NN), ppTh(t), pVh(0), stack(s), nbcperiodic(n), periodic(p) {
  }    // take a pmesh3 and use the pmeshS
  v_fesS(int NN, const v_fesS *f, Stack s, int n, Expression *p)
    : N(NN), ppTh(f->ppTh), pVh(0), stack(s), nbcperiodic(n), periodic(p) {}

  // void destroy(){ ppTh=0;pVh=0; delete this;}
  virtual ~v_fesS( ) {}
  bool buildperiodic(Stack stack, KN< int > &ndfe);
  virtual FESpaceS *buildupdate(KN< int > &ndfe) { return 0; }
  virtual FESpaceS *buildupdate( ) { return 0; };

  const void * getppTh(){ return ppTh; }; // Morice : get the pointer on the pointer of the mesh
  void * getpVh(){ return pVh; };         // Morice : get the pointer of FESpace
  void * getpVh()const{ return pVh; };         // Morice : get the pointer of FESpac
  int getN(){ return N;};                // Morice : get the number of item of the FESpace
};

// 3D curve
class v_fesL : public generic_v_fes { //: public RefCounter {
 public:
  typedef pfesL pfes;
  typedef FESpaceL FESpace;

  static const int dHat = 1, d = 3;
  const int N;
  const pmeshL *ppTh;    // adr du maillage

  CountPointer< FESpaceL > pVh;

  Stack stack;    // the stack is use whith periodique expression

  int nbcperiodic;
  Expression *periodic;

  operator FESpaceL *( ) {
    throwassert(dHat == 1);
    if (!pVh || *ppTh != &pVh->Th) pVh = CountPointer< FESpaceL >(update( ), true);
    return pVh;
  }
  FESpaceL *update( );

  v_fesL(int NN, const pmeshL *t, Stack s, int n, Expression *p)    /// TODO
    : N(NN), ppTh(t), pVh(0), stack(s), nbcperiodic(n), periodic(p) {
  }    // take a pmesh3 and use the pmeshS
  v_fesL(int NN, const v_fesL *f, Stack s, int n, Expression *p)
    : N(NN), ppTh(f->ppTh), pVh(0), stack(s), nbcperiodic(n), periodic(p) {}

  // void destroy(){ ppTh=0;pVh=0; delete this;}
  virtual ~v_fesL( ) {}
  bool buildperiodic(Stack stack, KN< int > &ndfe);
  virtual FESpaceL *buildupdate(KN< int > &ndfe) { return 0; }
  virtual FESpaceL *buildupdate( ) { return 0; };

  const void *getppTh(){ return ppTh; }; // Morice : get the pointer on the pointer of the mesh
  void *getpVh(){ return pVh; };         // Morice : get the pointer of FESpace
  void * getpVh()const{ return pVh; };         // Morice : get the pointer of FESpac
  int getN(){ return N;};                // Morice : get the number of item of the FESpace
};

// 2d
class pfes_tef : public v_fes {
 public:
  const TypeOfFE *tef;
  pfes_tef(const pmesh *t, const TypeOfFE *tt, Stack s = NullStack, int n = 0, Expression *p = 0)
    : v_fes(tt->N, t, s, n, p), tef(tt) {
    operator FESpace *( );
  }
  FESpace *buildupdate(int &nbdfv, KN< int > &ndfv, int &nbdfe, KN< int > &ndfe) {
    return *ppTh ? new FESpace(**ppTh, *tef, nbdfv, (int *)ndfv, nbdfe, (int *)ndfe) : 0;
  }
  FESpace *buildupdate( ) { return *ppTh ? new FESpace(**ppTh, *tef) : 0; }
};

class pfes_tefk : public v_fes {
 public:
  const TypeOfFE **tef;
  const int k;
  pfes_tefk(const pmesh *t, const TypeOfFE **tt, int kk, Stack s = NullStack, int n = 0,
            Expression *p = 0)
    : v_fes(sum(tt, &Fem2D::TypeOfFE::N, kk), t, s, n, p), tef(tt), k(kk) {
    //cout << "pfes_tefk const" << tef << " " << this << endl; // Morice
    operator FESpace *( );
  }
  FESpace *buildupdate( ) {
    // cout << "pfes_tefk upd:" << tef << " " << this <<  endl;
    assert(tef);
    return *ppTh ? new FESpace(**ppTh, tef, k) : 0;
  }
  virtual ~pfes_tefk( ) { delete[] tef; }
  FESpace *buildupdate(int &nbdfv, KN< int > &ndfv, int &nbdfe, KN< int > &ndfe) {
    assert(tef);
    return *ppTh ? new FESpace(**ppTh, tef, k, nbdfv, ndfv, nbdfe, ndfe) : 0;
  }
};

// 3D volume
class pfes3_tef : public v_fes3 {
 public:
  const TypeOfFE3 *tef;
  pfes3_tef(const pmesh3 *t, const TypeOfFE3 *tt, Stack s = NullStack, int n = 0, Expression *p = 0)
    : v_fes3(tt->N, t, s, n, p), tef(tt) {
    operator FESpace3 *( );
  }
  FESpace3 *buildupdate(KN< int > &ndfe) {
    return *ppTh ? new FESpace3(**ppTh, *tef, ndfe.size( ) / 2, ndfe) : 0;
  }
  FESpace3 *buildupdate( ) { return *ppTh ? new FESpace3(**ppTh, *tef) : 0; }
};

class pfes3_tefk : public v_fes3 {
 public:
  const TypeOfFE3 **tef;
  const int k;
  KN< GTypeOfFE< Mesh3 > const * > atef;
  GTypeOfFESum< Mesh3 > tefs;

  static int sum(const Fem2D::TypeOfFE3 **l, int const Fem2D::TypeOfFE3::*p, int n) {
    int r = 0;
    for (int i = 0; i < n; i++) r += l[i]->*p;
    return r;
  }

  pfes3_tefk(const pmesh3 *t, const Fem2D::TypeOfFE3 **tt, int kk, Stack s = NullStack, int n = 0,
             Expression *p = 0)
    : v_fes3(sum((const Fem2D::TypeOfFE3 **)tt, &Fem2D::TypeOfFE3::N, kk), t, s, n, p), tef(tt),
      k(kk), atef(kk, tt), tefs(atef)

  {
    // cout << "pfes_tefk const" << tef << " " << this << endl;
    operator FESpace3 *( );
  }
  FESpace3 *buildupdate( ) {
    // cout << "pfes_tefk upd:" << tef << " " << this <<  endl;
    // assert(tef);
    return *ppTh ? new FESpace3(**ppTh, tefs) : 0;
  }
  virtual ~pfes3_tefk( ) { delete[] tef; }
  FESpace3 *buildupdate(KN< int > &ndfe) {
    return *ppTh ? new FESpace3(**ppTh, tefs, ndfe.size( ) / 2, ndfe) : 0;
  }
};

// 3D surface
class pfesS_tef : public v_fesS {
 public:
  const TypeOfFES *tef;
  pfesS_tef(const pmeshS *t, const TypeOfFES *tt, Stack s = NullStack, int n = 0, Expression *p = 0)
    : v_fesS(tt->N, t, s, n, p), tef(tt) {
    operator FESpaceS *( );
  }
  FESpaceS *buildupdate(KN< int > &ndfe) {
    return *ppTh ? new FESpaceS(**ppTh, *tef, ndfe.size( ) / 2, ndfe) : 0;
  }
  FESpaceS *buildupdate( ) { return *ppTh ? new FESpaceS(**ppTh, *tef) : 0; }
};

class pfesS_tefk : public v_fesS {
 public:
  const TypeOfFES **tef;
  const int k;
  KN< GTypeOfFE< MeshS > const * > atef;
  GTypeOfFESum< MeshS > tefs;

  static int sum(const Fem2D::TypeOfFES **l, int const Fem2D::TypeOfFES::*p, int n) {
    int r = 0;
    for (int i = 0; i < n; i++) r += l[i]->*p;
    return r;
  }

  pfesS_tefk(const pmeshS *t, const Fem2D::TypeOfFES **tt, int kk, Stack s = NullStack, int n = 0,
             Expression *p = 0)
    : v_fesS(sum((const Fem2D::TypeOfFES **)tt, &Fem2D::TypeOfFES::N, kk), t, s, n, p), tef(tt),
      k(kk), atef(kk, tt), tefs(atef)

  {
    // cout << "pfes_tefk const" << tef << " " << this << endl;
    operator FESpaceS *( );
  }
  FESpaceS *buildupdate( ) {
    // cout << "pfes_tefk upd:" << tef << " " << this <<  endl;
    // assert(tef);
    return *ppTh ? new FESpaceS(**ppTh, tefs) : 0;
  }
  virtual ~pfesS_tefk( ) { delete[] tef; }
  FESpaceS *buildupdate(KN< int > &ndfe) {
    return *ppTh ? new FESpaceS(**ppTh, tefs, ndfe.size( ) / 2, ndfe) : 0;
  }
};

// 3D curve
class pfesL_tef : public v_fesL {
 public:
  const TypeOfFEL *tef;
  pfesL_tef(const pmeshL *t, const TypeOfFEL *tt, Stack s = NullStack, int n = 0, Expression *p = 0)
    : v_fesL(tt->N, t, s, n, p), tef(tt) {
    operator FESpaceL *( );
  }
  FESpaceL *buildupdate(KN< int > &ndfe) {
    return *ppTh ? new FESpaceL(**ppTh, *tef, ndfe.size( ) / 2, ndfe) : 0;
  }
  FESpaceL *buildupdate( ) { return *ppTh ? new FESpaceL(**ppTh, *tef) : 0; }
};

class pfesL_tefk : public v_fesL {
 public:
  const TypeOfFEL **tef;
  const int k;
  KN< GTypeOfFE< MeshL > const * > atef;
  GTypeOfFESum< MeshL > tefs;

  static int sum(const Fem2D::TypeOfFEL **l, int const Fem2D::TypeOfFEL::*p, int n) {
    int r = 0;
    for (int i = 0; i < n; i++) r += l[i]->*p;
    return r;
  }

  pfesL_tefk(const pmeshL *t, const Fem2D::TypeOfFEL **tt, int kk, Stack s = NullStack, int n = 0,
             Expression *p = 0)
    : v_fesL(sum((const Fem2D::TypeOfFEL **)tt, &Fem2D::TypeOfFEL::N, kk), t, s, n, p), tef(tt),
      k(kk), atef(kk, tt), tefs(atef)

  {
    // cout << "pfes_tefk const" << tef << " " << this << endl;
    operator FESpaceL *( );
  }
  FESpaceL *buildupdate( ) {
    // cout << "pfes_tefk upd:" << tef << " " << this <<  endl;
    // assert(tef);
    return *ppTh ? new FESpaceL(**ppTh, tefs) : 0;
  }
  virtual ~pfesL_tefk( ) { delete[] tef; }
  FESpaceL *buildupdate(KN< int > &ndfe) {
    return *ppTh ? new FESpaceL(**ppTh, tefs, ndfe.size( ) / 2, ndfe) : 0;
  }
};
/*
struct infoBlockCompositeFESpace{

  // struct useful for composite FESpace
  KN<int> beginBlock;
  KN<int> indexBlock;  
  KN<int> localIndexInTheBlock;

  infoBlockCompositeFESpace( const KN<int> & b, const KN<int> & i, const KN<int> & li ): beginBlock(b), indexBlock(i), localIndexInTheBlock(li) {};

};
*/

class vect_generic_v_fes {  //: public RefCounter{
  public:

  typedef pgenericfes A;

  public:
  int N;            // size of the vector ( the number of FESpace )
  A *vect;          // the vector containing the adress of generic_v_fes 
  KN<int> typeFE;   // ===== type of FE ====
                 // The value of a type correspond to the second parameter of << forEachTypePtrfspace >>
                 // in for pfes, pfes3, pfesS, pfesL in [[ fflib/lgfem.cpp ]].

  //CountPointer< FESpace > pVh;  // A voir avec Frederic pour savoir si il faut rajouter un 
                                  // un CountPointer

  
  Stack stack; // the stack is use whith periodique expression

  
  // constructor
  vect_generic_v_fes(const E_FEarray &args, Stack s): N( args.size() ), typeFE( args.size() ), stack(s){ 
  
 
  // constructor
  //vect_generic_v_fes(const E_FEarray &args, Stack s=NullStack): N( args.size() ), typeFE( args.size() ){ 
    // == PAC(e) == //
    // Morice :: attention <<pfes>> est redéfinis pour cette classe.
    //           Il faut faire appel au type de <<pfes>> associé à v_fes:
    //                 ::pfes == typedef v_fes *pfes;
    //
    typedef ::pfes pfesTrue; 

    aType atfs[] = {atype< ::pfes * >( ), atype< pfes3 * >( ), atype< pfesS * >( ), atype< pfesL * >( )};
    vect = new A[N];
    
    if(verbosity>3){
      for(int i=0; i<args.size(); i++){
        cout << "info FESpace["<<i<<"] === " << endl;
        cout << "nbitem= " << args[i].nbitem() << endl;
      }
    }
    
    for( int i=0; i<args.size(); i++ ){
        //cout << "i=" << i << endl;
        if( atype< ::pfes * >() == args[i].left( ) ){
          //a = to< atfs[iii] >(args[i]);
          pfesTrue * tmp = GetAny< pfesTrue *> ( args[i].eval(stack)) ;
          vect[i] = *tmp;
          //vect[i] = args[i].f;
          typeFE[i] = args[i].TYPEOFID();
        }
        else if( atype< pfes3 * >() == args[i].left( ) ){
          //a = to< atfs[iii] >(args[i]);
          pfes3 * tmp = GetAny< pfes3*> ( args[i].eval(stack)) ;
          vect[i] = *tmp;
          //vect[i] = args[i].RightValue();
          typeFE[i] = args[i].TYPEOFID();
        }
        else if( atype< pfesS * >() == args[i].left( ) ){
          //a = to< atfs[iii] >(args[i]);
          pfesS * tmp = GetAny< pfesS*> ( args[i].eval(stack)) ;
          vect[i] = *tmp;
          //vect[i] = args[i].RightValue();
          typeFE[i] = args[i].TYPEOFID();
        }
        else if( atype< pfesL * >() == args[i].left( ) ){
          //a = to< atfs[iii] >(args[i]);
          pfesL *tmp = GetAny< pfesL* > ( args[i].eval(stack)) ;
          vect[i] = *tmp;
          //vect[i] = args[i].RightValue();
          typeFE[i] = args[i].TYPEOFID();
        }
        else{
          cerr << "Error in the constructor of <<vect_generic_v_fes>> " << endl;
          cerr << "The aType of the FESpace " << args[i].left( ) << " is not considered." << endl;
          cerr << "The aType considered for FESpace is " << atfs[0]  << "." << endl;
          cerr << "The aType considered for FESpace is " << atfs[1]  << "." << endl;
          cerr << "The aType considered for FESpace is " << atfs[2]  << "." << endl;
          cerr << "The aType considered for FESpace is " << atfs[3]  << "." << endl;
          cerr << "type pgenericfes= " << atype< pgenericfes *>() << endl;
          ffassert(0);
        }
    }
    /*
    for( int i=0; i<args.size(); i++ ){
      cout << "getpVh(), i= " << i << "; " << vect[i]->getpVh() << endl;
    }
    */
    update();
  }

  // destructor
  void destroy( ) {
    // pVh = 0;
    delete [] vect;  // A demander à PH pour la memoire
    //vect = 0;
    delete this;
  }

  ~vect_generic_v_fes(){
  }
  /**
    * @brief update the composant of the composite FESpace if a mesh have changed
    * 
    */
  void update( ){
    // loop over the FESpace
    for(int i=0; i<N; i++){
      int tt = typeFE[i];
      if(tt == 2){
        pfes pFES = dynamic_cast<pfes>( vect[i] );
        ffassert(pFES);
        // cout << "est ce que je vais  :"<< (!mac->pVh || *(mac->ppTh) != &mac->pVh->Th) << endl;
        pmesh* ppTh = (pmesh *) vect[i]->getppTh();
        FESpace * pVh = (FESpace *) vect[i]->getpVh(); 
        // pfes pFES = dynamic_cast<pfes>( vect[i] );
        // cout << &pVh->Th << " "<< &mac->pVh->Th << endl;
        // cout << "est ce que je vais 2: "<< (!pVh || *(ppTh) != &pVh->Th) << endl;
        if (!pVh || *ppTh != &pVh->Th) pFES->pVh = CountPointer< FESpace >(pFES->update( ), true);
      }
      else if(tt == 3){
        pmesh3* ppTh = (pmesh3 *) vect[i]->getppTh();
        const pfes3 pFES = dynamic_cast<const pfes3>( vect[i] );
        ffassert(pFES); 
        FESpace3 * pVh =  (FESpace3*) vect[i]->getpVh();
        if (!pVh || *ppTh != &pVh->Th) pFES->pVh = CountPointer< FESpace3 >(pFES->update( ), true);
  
      }
      else if(tt == 4){
        pmeshS* ppTh = (pmeshS *) vect[i]->getppTh();
        const pfesS pFES = dynamic_cast<const pfesS>( vect[i] );
        ffassert(pFES); 
        FESpaceS * pVh =  (FESpaceS*) vect[i]->getpVh();
        if (!pVh || *ppTh != &pVh->Th) pFES->pVh = CountPointer< FESpaceS >(pFES->update( ), true);
        
      }
      else if(tt == 5){
        pmeshL* ppTh = (pmeshL *) vect[i]->getppTh();
        const pfesL pFES = dynamic_cast<const pfesL>( vect[i] );
        ffassert(pFES);
        FESpaceL * pVh =  (FESpaceL*) vect[i]->getpVh();
        if (!pVh || *ppTh != &pVh->Th) pFES->pVh = CountPointer< FESpaceL >(pFES->update( ), true);
      }
      /*else{
        cerr << "== error in the type of FESpace ===" << endl;
        ffassert(0);
      }
      */
    }
  }


  // acces to the different pgenericfes
  A &operator[](int i){ return vect[i];}

  //=== Morice :: ATTENTION : PEUT ETRE A METTRE AILLEURS ===//

  void printPointer(){
    cout <<" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% " <<endl;
    cout <<" %%%%        pointer FESpace       %%%%%%%% " <<endl;
    for(int i=0; i< N; i++){
      cout << "pointer epace i =" << vect[i]->getpVh() << endl;
    }
    cout <<" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% " <<endl;
  }

  KN<int> vectOfNbOfDF(){
    KN<int> result(N);
    for(int i=0; i< N; i++){
      int tt = typeFE[i];
      if(tt == 2){
        const FESpace * PUh =  (FESpace*) vect[i]->getpVh();
        result[i] = PUh->NbOfDF;
      }
      else if(tt == 3){
        const FESpace3 * PUh =  (FESpace3*) vect[i]->getpVh();
        result[i] =  PUh->NbOfDF;
      }
      else if(tt == 4){
        const FESpaceS * PUh =  (FESpaceS*) vect[i]->getpVh();
        result[i] = PUh->NbOfDF;
      }
      else if(tt == 5){
        const FESpaceL * PUh =  (FESpaceL*) vect[i]->getpVh();
        result[i] = PUh->NbOfDF;
      }
    }
    //cout << "vectOfNbOfDF=" << result << endl;
    return result;    
  }

  int totalNbOfDF(){
    KN<int> arrayNbOfDF = vectOfNbOfDF();
    int tNbOfDF=0;
    for(int i=0; i< N; i++){
        tNbOfDF += arrayNbOfDF[i];
    }
    //cout << "total NbOfDF=" << tNbOfDF << endl;
    return tNbOfDF;
  }

  KN<int> vectOfNbitem(){
    KN<int> result(N);
    for(int i=0; i< N; i++){
      // get the number of item of i-th component of the composite vector of FESpace
      result[i] = vect[i]->getN(); 
    }
    //cout << "vectOfNbitem=" << result << endl;
    return result;
  }

  int totalNbitem(){
    KN<int> arrayNbItem = vectOfNbitem();
    int tNbItem=0;
    for(int i=0; i< N; i++){
        tNbItem += arrayNbItem[i];
    }
    return tNbItem;
  }
};

class pfes_fes : public v_fes {
 public:
  pfes *Vh;
  int n;
  pfes_fes(pfes *Vhh, int nn, Stack s = NullStack, int n = 0, Expression *p = 0)
    : v_fes((**Vhh).N * nn, static_cast< const v_fes * >(*Vhh), s, n, p), Vh(Vhh), n(nn) {
    operator FESpace *( );
  };
  FESpace *buildupdate( ) { return new FESpace(*(FESpace *)**Vh, n); }
  FESpace *buildupdate(int &nbdfv, KN< int > &ndfv, int &nbdfe, KN< int > &ndfe) {
    InternalError(" No way to define a periodic BC in this case: tensorisation of FEspace ");
    //  return  new FESpace(***Vh,n,nbdfv,ndfv,nbdfe,ndfe);
  }
};

template< class K, class v_fes >
class FEbase;    // <<FEbase>>
template< class K, class v_fes >
class FEcomp {
 public:
  typedef typename v_fes::pfes pfes;
  typedef typename v_fes::FESpace FESpace;

  friend class FEbase< K, v_fes >;
  FEbase< K, v_fes > *base;
  int comp;
  FEcomp(FEbase< K, v_fes > *b, int c) : base(b), comp(c){};

 private:    // rule of programming
  FEcomp(const FEcomp &);
  void operator=(const FEcomp &);
};

template< class K, class v_fes >
class FEbase {
 public:
  typedef typename v_fes::pfes pfes;
  typedef typename v_fes::FESpace FESpace;

  v_fes *const *pVh;             // pointeur sur la variable stockant FESpace;
  KN< K > *xx;                   // value
  CountPointer< FESpace > Vh;    // espace courant

  KN< K > *x( ) { return xx; }

  FEbase(const pfes *ppVh) : pVh(ppVh), xx(0), Vh(0) {}

  ~FEbase( ) { delete xx; }
  void destroy( ) {    // cout << "~FEbase  destroy " << this << endl;
    delete this;
  }

  void operator=(KN< K > *y) {
    Vh = **pVh;
    throwassert((bool)Vh);
    if (xx) delete xx;
    xx = y;
    ffassert(y->N( ) == Vh->NbOfDF);
  }

  void operator=(KN_< K > &y) {
    Vh = **pVh;
    throwassert((bool)Vh);
    if (xx) {    // resize if need
      if (xx->N( ) != Vh->NbOfDF) delete xx;
      xx = 0;
    }
    if (!xx) xx = new KN< K >(Vh->NbOfDF);
    ffassert(SameShape(y, *xx));
    *xx = y;
  }

  FESpace *newVh( ) {
    throwassert(pVh);
    const pfes pp = *pVh;
    // cout << pVh << " " << *pVh << endl;
    return *pp;
  }

  operator FESpace &( ) {
    throwassert(Vh);
    return *Vh;
  }

 private:    // rule of programming
  FEbase(const FEbase &);
  void operator=(const FEbase &);
};

template< class K >
class FEbaseArrayKn {
 public:    // for eigen value
  int N;
  FEbaseArrayKn(int NN) : N(NN) {}
  virtual void set(int i, KN_< K >) = 0;
  virtual KN< K > *get(int i) const = 0;    // for P. Jolivet
  virtual void resize(int i) = 0;           // for P. Jolivet
};

template< class K, class v_fes >
class FEbaseArray : public FEbaseArrayKn< K > {
 public:
  typedef typename v_fes::pfes pfes;
  typedef typename v_fes::FESpace FESpace;

  // int N;
  FEbase< K, v_fes > **xx;
  FEbaseArray(const pfes *ppVh, int NN)
    : FEbaseArrayKn< K >(NN), xx(new FEbase< K, v_fes > *[std::max(NN, 1)]) {
    for (int i = 0; i < std::max(this->N, 1); i++) xx[i] = new FEbase< K, v_fes >(ppVh);
  }
  ~FEbaseArray( ) {
    //  cout << " ~FEbaseArray " << endl;
    for (int i = 0; i < std::max(this->N, 1); i++) xx[i]->destroy( );
    delete[] xx;
  }
  void destroy( ) {    // cout << " destroy ~FEbaseArray " << endl;
    delete this;
  }
  FEbase< K, v_fes > **operator[](int i) const {
    if (xx == 0 || i < 0 || i >= this->N) ExecError("Out of bound in FEbaseArray");
    return xx + i;
  }

  void resize(int i) {
    if (xx != 0 && i > 0 && i != this->N) {
      FEbase< K, v_fes > **yy = new FEbase< K, v_fes > *[i];
      if (i > this->N) {
        for (unsigned int j = 0; j < std::max(this->N, 1); ++j) yy[j] = xx[j];
        for (unsigned int j = std::max(this->N, 1); j < i; ++j)
          yy[j] = new FEbase< K, v_fes >(xx[0]->pVh);
      } else {
        for (unsigned int j = 0; j < i; ++j) yy[j] = xx[j];
        for (unsigned int j = i; j < this->N; ++j) xx[j]->destroy( );
      }
      FEbase< K, v_fes > **oldXx = this->xx;
      this->xx = yy;
      delete[] oldXx;
      this->N = i;
    }
  }

  void set(int i, KN_< K > v) { **(operator[](i)) = v; }

  KN< K > *get(int i) const { return (**(operator[](i))).xx; }

 private:    // rule of programming
  FEbaseArray(const FEbaseArray &);
  void operator=(const FEbaseArray &);
};

void GetPeriodic(const int d, Expression perio, int &nbcperiodic, Expression *&periodic);
int GetPeriodic(Expression bb, Expression &b);// 1dhat
int GetPeriodic(Expression bb, Expression &b, Expression &f);//2dhat
int GetPeriodic(Expression bb, Expression &b, Expression &f1, Expression &f2);//3dhat

C_F0 NewFEarray(const char *id, Block *currentblock, C_F0 &fespacetype, CC_F0 init, bool cplx,
                int dim);
C_F0 NewFEarray(ListOfId *ids, Block *currentblock, C_F0 &fespacetype, CC_F0 init, bool cplx,
                int dim);
C_F0 NewFEvariable(const char *id, Block *currentblock, C_F0 &fespacetype, CC_F0 init, bool cplx,
                   int dim);
C_F0 NewFEvariable(ListOfId *ids, Block *currentblock, C_F0 &fespacetype, CC_F0 init, bool cplx,
                   int dim);
inline C_F0 NewFEvariable(const char *id, Block *currentblock, C_F0 &fespacetype, bool cplx,
                          int dim) {
  CC_F0 init;
  init = 0;
  return NewFEvariable(id, currentblock, fespacetype, init, cplx, dim);
}

inline C_F0 NewFEvariable(ListOfId *ids, Block *currentblock, C_F0 &fespacetype, bool cplx,
                          int dim) {
  CC_F0 init;
  init = 0;
  return NewFEvariable(ids, currentblock, fespacetype, init, cplx, dim);
}

KN<size_t> dimFESpaceImage(const basicAC_F0 &args);
aType typeFESpace(const basicAC_F0 &args);

template< class K, class vv_fes, class FE = FEbase< K, vv_fes > >
class E_FEcomp : public E_F0mps {
 public:
  //  typedef FEbase<K,vv_fes> FE;
  typedef vv_fes v_fes;
  typedef pair< FE *, int > Result;
  Expression a0;
  const int comp, N;

  AnyType operator( )(Stack s) const {
    return SetAny< Result >(Result(*GetAny< FE ** >((*a0)(s)), comp));
  }
  E_FEcomp(const C_F0 &x, const int cc, int NN) : a0(x.LeftValue( )), comp(cc), N(NN) {
    if (x.left( ) != atype< FE ** >( ))
      cout << "E_FEcomp: Bug " << *x.left( ) << " != " << *atype< FE ** >( ) << "  case "
           << typeid(K).name( ) << endl;
    // CompileError("E_FEcomp: Bug ?");
    throwassert(x.left( ) == atype< FE ** >( ) && a0);
  }
  operator aType( ) const { return atype< Result >( ); }
};

// typedef double R;
typedef pair< FEbase< double, v_fes > *, int > aFEvarR;
typedef pair< FEbaseArray< double, v_fes > *, int > aFEArrayR;
typedef pair< FEbase< Complex, v_fes > *, int > aFEvarC;

template< class K >
class interpolate_f_X_1 : public OneOperator {
 public:
  //  to interpolate f o X^{-1}
  typedef FEbase< K, v_fes > *pferbase;
  typedef FEbaseArray< K, v_fes > *pferbasearray;
  typedef pair< pferbase, int > pfer;
  typedef pair< pferbasearray, int > pferarray;
  class type {};
  class CODE : public E_F0mps {
   public:
    Expression f;
    Expression x, y;

    CODE(const basicAC_F0 &args) : f(args[0].LeftValue( )) {
      const E_Array *X(dynamic_cast< const E_Array * >(args[1].LeftValue( )));
      if (!X || X->size( ) != 2) {
        CompileError("array of 2 double x,y   f([xx,yy]) = ... ");
      }
      x = to< double >((*X)[0]);
      y = to< double >((*X)[1]);
      assert(f && x && y);
    }
    AnyType operator( )(Stack) const {
      ExecError("No evaluation");
      return Nothing;
    }
    operator aType( ) const { return atype< void >( ); }

  };    //  end class CODE

  E_F0 *code(const basicAC_F0 &args) const { return new CODE(args); }

  interpolate_f_X_1( )
    : OneOperator(map_type[typeid(type).name( )], map_type[typeid(pfer).name( )],
                  map_type[typeid(E_Array).name( )]) {}
};

inline FESpace *v_fes::update( ) {
  assert(dHat == 2);
  if (!*ppTh) return 0;
  if (nbcperiodic) {
    assert(periodic);
    const Mesh &Th(**ppTh);
    KN< int > ndfv(Th.nv);
    KN< int > ndfe(Th.neb);
    int nbdfv, nbdfe;
    bool ret = buildperiodic(stack, nbdfv, ndfv, nbdfe, ndfe);
    return ret ? buildupdate(nbdfv, ndfv, nbdfe, ndfe) : buildupdate( );
  } else
    return buildupdate( );
}

inline FESpace3 *v_fes3::update( ) {
  assert(dHat == 3);
  if (!*ppTh) return 0;
  if (nbcperiodic) {
    assert(periodic);
    KN< int > ndfe;
    bool ret = buildperiodic(stack, ndfe);
    return ret ? buildupdate(ndfe) : buildupdate( );
  } else
    return buildupdate( );
}

inline FESpaceS *v_fesS::update( ) {
  assert(dHat == 2);
  if (!*ppTh) return 0;
  if (nbcperiodic) {
    assert(periodic);
    KN< int > ndfe;
    bool ret = buildperiodic(stack, ndfe);
    return ret ? buildupdate(ndfe) : buildupdate( );
  } else
    return buildupdate( );
}

inline FESpaceL *v_fesL::update( ) {
  assert(dHat == 1);
  if (!*ppTh) return 0;
  if (nbcperiodic) {
    assert(periodic);
    KN< int > ndfe;
    bool ret = buildperiodic(stack, ndfe);
    return ret ? buildupdate(ndfe) : buildupdate( );
  } else
    return buildupdate( );
}

class TabFuncArg {
 public:
  typedef double R;
  typedef Fem2D::R2 R2;
  Stack s;
  int nb;
  Expression *e;
  void eval(R *v) const {
    for (int i = 0; i < nb; i++)
      if (e[i]) {
        v[i] = GetAny< R >((*e[i])(s));
        //    cout <<" (" << i << " " << *v << ") ";
      } else
        v[i] = 0;
  }
  R2 eval_X( ) const { return R2(GetAny< R >((*e[nb - 2])(s)), GetAny< R >((*e[nb - 1])(s))); }
  void eval_2(R *v) const {
    for (int i = 0; i < nb - 2; i++)
      if (e[i])
        v[i] = GetAny< R >((*e[i])(s));
      else
        v[i] = 0;
  }
  TabFuncArg(Stack ss, int n) : s(ss), nb(n), e(new Expression[n]) {}
  void operator=(int j) {
    ffassert(j == 0);
    for (int i = 0; i < nb; i++) e[i] = 0;
  }    // resert
  Expression &operator[](int i) { return e[i]; }
  ~TabFuncArg( ) { delete[] e; }

 private:    // no copy
  TabFuncArg(const TabFuncArg &);
  void operator=(const TabFuncArg &);
};

template< class K >
class E_F_StackF0F0opt2 : public E_F0mps {
 public:
  typedef AnyType (*func)(Stack, Expression, Expression);
  func f;
  Expression a0, a1;
  E_F_StackF0F0opt2(func ff, Expression aa0, Expression aa1) : f(ff), a0(aa0), a1(aa1) {

    const E_FEcomp< K, v_fes > *e = dynamic_cast< const E_FEcomp< K, v_fes > * >(a0);
    if (e && e->N != 1) {
      cerr << " err interpolation of  no scalar FE function componant (" << e->comp << " in  0.."
           << e->N << ") <<  with scalar function \n";
      CompileError("interpolation of  no scalar FE function componant with scalar function ");
    }

    /*  //  inutil car a0 est un composante d'un vecteur ????
       // pour l'instant on a juste une erreur a l'execution
       // et non a la compilation.

         if (a0->nbitem() != a1->nbitem() )
          { // bofbof
            cerr << " err interpolation of  no scalar FE function ("<<  a0->nbitem() << " != "  <<
       a1->nbitem()  << ") <<  with scalar function \n" ; CompileError("interpolation of  no scalar
       FE function  with scalar function ");
          } */
    deque< pair< Expression, int > > ll;
    MapOfE_F0 m;
    size_t top = currentblock->OffSet(0), topbb = top;    // FH. bofbof ???
    int ret = aa1->Optimize(ll, m, top);
    a1 = new E_F0_Optimize(ll, m, ret);
    currentblock->OffSet(top - topbb);
  }
  AnyType operator( )(Stack s) const { return (*f)(s, a0, a1); }
};

template< class R >
inline ostream &ShowBound(const KN< R > &y, ostream &f) {
  f << "  -- vector function's bound  " << y.min( ) << " " << y.max( );
  return f;
}
template<>
inline ostream &ShowBound< Complex >(const KN< Complex > &y, ostream &f) {
  f << "  -- vector function's bound : (no complex Value) ";
  return f;
}

template< class K >
class Op3_K2R : public ternary_function< K, R, R, K > {
 public:
  class Op : public E_F0mps {
   public:
    Expression a, b, c;
    Op(Expression aa, Expression bb, Expression cc) : a(aa), b(bb), c(cc) {}
    AnyType operator( )(Stack s) const {
      R xx(GetAny< R >((*b)(s)));
      R yy(GetAny< R >((*c)(s)));
      MeshPoint &mp = *MeshPointStack(s), mps = mp;
      mp.set(xx, yy, 0.0);
      AnyType ret = (*a)(s);
      mp = mps;
      return ret;
    }
  };
};
template< class K >
class Op4_K2R : public quad_function< K, R, R, R, K > {
 public:
  class Op : public E_F0mps {
   public:
    Expression a, b, c, d;
    Op(Expression aa, Expression bb, Expression cc, Expression dd) : a(aa), b(bb), c(cc), d(dd) {}
    AnyType operator( )(Stack s) const {
      // cout <<"class Op4_K2R : public quad_function<K,R,R,R,K>" << endl;
      R xx(GetAny< R >((*b)(s)));
      R yy(GetAny< R >((*c)(s)));
      R zz(GetAny< R >((*d)(s)));
      MeshPoint &mp = *MeshPointStack(s), mps = mp;
      mp.set(xx, yy, zz);
      AnyType ret = (*a)(s);
      mp = mps;
      return ret;
    }
  };
};

// 3D
template< class K, class v_fes >
class E_set_fev3 : public E_F0mps {
 public:
  typedef typename v_fes::pfes pfes;
  typedef typename v_fes::FESpace FESpace;
  typedef typename FESpace::Mesh Mesh;
  typedef typename FESpace::FElement FElement;
  typedef typename Mesh::Element Element;
  typedef typename Mesh::Vertex Vertex;
  typedef typename Mesh::RdHat RdHat;
  typedef typename Mesh::Rd Rd;

  E_Array aa;
  Expression ppfe;
  bool optimize, optimizecheck;

  vector< size_t > where_in_stack_opt;
  Expression optiexp0, optiexpK;

  E_set_fev3(const E_Array *a, Expression pp);

  AnyType operator( )(Stack) const;
  operator aType( ) const { return atype< void >( ); }
};

// 2d
template< class K >
class E_set_fev : public E_F0mps {
 public:
  const int dim;
  E_Array aa;
  Expression ppfe;
  bool optimize, optimizecheck;
  vector< size_t > where_in_stack_opt;
  Expression optiexp0, optiexpK;

  E_set_fev(const E_Array *a, Expression pp, int ddim = 2);

  AnyType operator( )(Stack) const;
  AnyType Op2d(Stack) const;
  AnyType Op3d(Stack) const;
  AnyType OpS(Stack) const;

  operator aType( ) const { return atype< void >( ); }
};

bool InCircularList(const int *p, int i, int k);
template< class T >
int numeroteclink(KN_< T > &ndfv);

//  for ...   vectorial FE array ..

template< class K, class v_fes >
class E_FEcomp_get_elmnt_array : public E_F0 {
 public:
  typedef FEbaseArray< K, v_fes > *pfekbasearray;
  typedef FEbase< K, v_fes > *pfekbase;
  typedef pair< pfekbase, int > pfek;
  typedef pair< pfekbasearray, int > pfekarray;
  typedef pfek R;
  typedef pfekarray A;
  typedef long B;
  typedef FEbaseArray< K, v_fes > CFE;
  typedef E_FEcomp< K, v_fes, CFE > E_KFEArray;

  Expression a0, a1;
  const E_KFEArray *a00;
  const int comp, N;

  E_FEcomp_get_elmnt_array(Expression aa0, Expression aa1, int compp, int NN,
                           const E_KFEArray *aa00)
    : a0(aa0), a1(aa1), a00(aa00), comp(compp), N(NN) {}
  AnyType operator( )(Stack s) const {
    return SetAny< R >(get_element(GetAny< A >((*a0)(s)), GetAny< B >((*a1)(s))));
  }
  bool MeshIndependent( ) const { return a0->MeshIndependent( ) && a1->MeshIndependent( ); }    //
};

template< class K, class v_fes >
class OneOperator2_FE_get_elmnt : public OneOperator {
  // ///  Add<pferbasearray*>("[","",new
  // OneOperator2_FEcomp<pferbase*,pferbasearray*,long>(get_element)); // not used ...
  //    Add<pferarray>("[","",new OneOperator2_<pfer,pferarray,long>(get_element));
  typedef FEbase< K, v_fes > *pfekbase;
  typedef FEbaseArray< K, v_fes > *pfekbasearray;
  typedef pair< pfekbase, int > pfek;
  typedef pair< pfekbasearray, int > pfekarray;

  typedef pfek R;
  typedef pfekarray A;
  typedef long B;
  typedef E_FEcomp_get_elmnt_array< K, v_fes > CODE;
  typedef FEbaseArray< K, v_fes > CFE;
  typedef E_FEcomp< K, v_fes, CFE > E_KFEArray;
  typedef E_FEcomp< K, v_fes > E_KFEi;

  aType t0, t1;    //  return type  type de f,  f(t1, t2)
 public:
  E_F0 *code(const basicAC_F0 &args) const {
    const E_KFEArray *afe = dynamic_cast< const E_KFEArray * >(args[0].LeftValue( ));
    //  cout << " build xxx  " << afe <<  " " << *args[0].left() <<  endl;
    // afe=0;// E_KFEi n'est pas le bon type on mame le vieux code ?????
    ffassert(afe);
    return new CODE(t0->CastTo(args[0]), t1->CastTo(args[1]), afe->comp, afe->N, afe);
  }
  OneOperator2_FE_get_elmnt( )
    : OneOperator(map_type[typeid(R).name( )], map_type[typeid(A).name( )],
                  map_type[typeid(B).name( )]),
      t0(map_type[typeid(A).name( )]), t1(map_type[typeid(B).name( )]) {}
};

template< typename KK, typename vv_fes, typename CC >
struct FFset3 {
  typedef KK K;
  typedef vv_fes v_fes;
  typedef vv_fes *pfes;
  typedef FEbase< K, vv_fes > **R;
  typedef R A;
  typedef pfes *B;
  typedef CC C;
  typedef Expression MC;
  static Expression Clone(Expression e) { return e; }
  static bool Check(MC l) { return true; }
  static void f(KN< K > *x, MC a, A pp, Stack stack) {
    CC at = GetAny< C >((*a)(stack));
    *x = at;
  }
};
template< typename KK, typename vv_fes, typename CC >
struct FFset3call {
  typedef KK K;
  typedef vv_fes v_fes;
  typedef vv_fes *pfes;
  typedef FEbase< K, vv_fes > **R;
  typedef R A;
  typedef pfes *B;
  typedef CC C;
  typedef Expression MC;
  static Expression Clone(Expression e) { return e; }
  static bool Check(MC l) { return true; }
  static void f(KN< K > *x, MC a, A pp, Stack stack) {
    CC at = GetAny< C >((*a)(stack));
    at.call(*x);
  }
};

template< typename FF >
class init_FE_eqarray : public OneOperator {
 public:
  // <pferbase*,pferbase*,pfes*
  typedef typename FF::K K;
  typedef typename FF::v_fes v_fes;
  typedef v_fes *pfes;
  typedef typename v_fes::FESpace FESpace;
  typedef FEbase< K, v_fes > **R;
  typedef R A;
  typedef pfes *B;
  typedef typename FF::C C;
  typedef typename FF::MC MC;
  class CODE : public E_F0 {
   public:
    Expression a0, a1;
    MC a2;
    CODE(Expression aa0, Expression aa1, Expression aa2) : a0(aa0), a1(aa1), a2(FF::Clone(aa2)) {
      ffassert(FF::Check(a2));
    }
    AnyType operator( )(Stack s) const {

      AnyType aa0 = (*a0)(s);
      AnyType aa1 = (*a1)(s);

      A pp = GetAny< A >(aa0);
      B a = GetAny< B >(aa1);
      *pp = new FEbase< K, v_fes >(a);
      KN< K > *x = (**pp).x( );
      if (!x) {    // defined
        FESpace *Vh = (*pp)->newVh( );
        throwassert(Vh);
        **pp = x = new KN< K >(Vh->NbOfDF);
        *x = K( );
      }
      int n = x->N( );
      FF::f(x, a2, pp, s);
      int nn = x->N( );
      ffassert(n == nn);
      return aa0;
    }
    virtual size_t nbitem( ) const { return a2->nbitem( ); }
    bool MeshIndependent( ) const {
      return a0->MeshIndependent( ) && a1->MeshIndependent( ) && a2->MeshIndependent( );
    }    //
  };
  init_FE_eqarray(int ppref)
    : OneOperator(map_type[typeid(R).name( )], map_type[typeid(A).name( )],
                  map_type[typeid(B).name( )], map_type[typeid(C).name( )]) {
    pref = ppref;
  }
  E_F0 *code(const basicAC_F0 &args) const {

    return new CODE(args[0], args[1], map_type[typeid(C).name( )]->CastTo(args[2]));
  }
};

#endif
