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
// AUTHORS : Jean-Marie Mirebeau
// E-MAIL  : jean-marie.mirebeau@math.u-psud.fr

// compilation : ff-c++ FreeFemQA.cpp -I/usr/local/boost_1_47_0

// *INDENT-OFF* //
// ff-c++-LIBRARY-dep:
// ff-c++-cpp-dep: GeometryQA.cpp
// *INDENT-ON* //

#include <iostream>
#include <cfloat>
#include <cmath>
using namespace std;
#include "ff++.hpp"
using namespace Fem2D;

// #include <boost/operators.hpp>
namespace mir {
#define _FLAGGED_BOUNDARY_
#include "Geometry.hpp"
}    // namespace mir

// the main class
// details of FreeFem meshes (connectivity, etc) in file GenericMesh.hpp

class MeshGenQA : public E_F0mps {
 public:
  static basicAC_F0::name_and_type name_param[];
  static const int n_name_param = 7;
  Expression nargs[n_name_param];    // store named args

  typedef const Mesh *Result;
  Expression expTh;
  Expression expM11;
  Expression expM12;
  Expression expM22;

  MeshGenQA(const basicAC_F0 &args) {
    args.SetNameParam(n_name_param, name_param, nargs);    // named args
    expTh = to< pmesh >(args[0]);                          // a the expression to get the mesh
    expM11 = to< double >(args[1]);
    expM12 = to< double >(args[2]);
    expM22 = to< double >(args[3]);
  }

  double arg(int i, Stack stack, double a) const {
    return nargs[i] ? GetAny< double >((*nargs[i])(stack)) : a;
  }

  long arg(int i, Stack stack, long a) const {
    return nargs[i] ? GetAny< long >((*nargs[i])(stack)) : a;
  }

  bool arg(int i, Stack stack, bool a) const {
    return nargs[i] ? GetAny< bool >((*nargs[i])(stack)) : a;
  }

  KN< double > *arg(int i, Stack stack, KN< double > *a) const {
    return nargs[i] ? GetAny< KN< double > * >((*nargs[i])(stack)) : a;
  }

  ~MeshGenQA( ) {}

  static ArrayOfaType typeargs( ) {
    return ArrayOfaType(atype< pmesh >( ), atype< double >( ), atype< double >( ),
                        atype< double >( ));
    ;
  }

  static E_F0 *f(const basicAC_F0 &args) { return new MeshGenQA(args); }

  AnyType operator( )(Stack s) const;    // la vraie fonction qui fait faire le boulot
};

basicAC_F0::name_and_type MeshGenQA::name_param[MeshGenQA::n_name_param] = {
  {"noIsoRef", &typeid(bool)},
  {"finalRefine", &typeid(bool)},
  {"exportIntermediateData", &typeid(bool)},
  {"Lip", &typeid(double)},
  {"exportToMathematica", &typeid(bool)},
  {"exportMetricToMathematica", &typeid(bool)},
  {"noRef", &typeid(bool)}};
AnyType MeshGenQA::operator( )(Stack stack) const {
  const bool noIsoRef = arg(0, stack, false);
  const bool finalRefine = arg(1, stack, false);
  const bool exportIntermediateData = arg(2, stack, false);
  unsigned int flag = 0;

  if (noIsoRef) {
    flag |= mir::Triangulation::hRQA_noIsoRef;
  }

  if (finalRefine) {
    flag |= mir::Triangulation::hRQA_finalRefine;
  }

  if (exportIntermediateData) {
    flag |= mir::Triangulation::hRQA_exportIntermediateData;
  }

  const double Lip = arg(3, stack, 5.);
  const bool exportToMathematica = arg(4, stack, false);
  const bool exportMetricToMathematica = arg(5, stack, false);
  const bool noRef = arg(6, stack, false);
  const Mesh *pTh = GetAny< pmesh >((*expTh)(stack));
  ffassert(pTh);
  const Mesh &Th = *pTh;

  class FFMetric2 : public mir::Metric2 {
    const MeshGenQA &MGQA_;
    Stack stack_;

   public:
    FFMetric2(const MeshGenQA &MGQA, Stack stack, double Lip) : MGQA_(MGQA), stack_(stack) {
      lip = Lip;
    }

    const mir::sym2 operator( )(const mir::R2 &P) const {
      MeshPointStack(stack_)->set(P.x, P.y);    // needs to be done three times ?
      MeshPointStack(stack_)->set(P.x, P.y);
      MeshPointStack(stack_)->set(P.x, P.y);
      return mir::sym2(GetAny< double >((*MGQA_.expM11)(stack_)),
                       GetAny< double >((*MGQA_.expM12)(stack_)),
                       GetAny< double >((*MGQA_.expM22)(stack_)));
    }
  };

  FFMetric2 ffMetric(*this, stack, Lip);
  const mir::Metric2 &metric = ffMetric;

  mir::Triangulation triQA(Th, metric);
  if (!triQA.check( )) {
    cout << "MeshGenQA : Error while importing mesh !\n";
    return false;
  }

  if (exportToMathematica) {
    triQA.export_to_Mathematica("ThFF.txt");
  }

  if (exportMetricToMathematica) {
    triQA.export_to_Mathematica_Metric("ThFF_Metric.txt");
  }

  if (!noRef) {
    triQA.hRefineQA(1, flag);
  }

  triQA.export_to_FreeFem("triQA.msh");

  if (exportToMathematica) {
    triQA.export_to_Mathematica("TriQA.txt");
  }

  if (exportMetricToMathematica) {
    triQA.export_to_Mathematica_Metric("TriQA_Metric.txt");
  }

  // generation de la class Mesh a partir des 3 tableaux : v,t,b
  {
    Mesh *m = triQA.export_to_Mesh( );    // new Mesh(nbv+nbt,nbt*3,neb,v,t,b);
    R2 Pn, Px;
    m->BoundingBox(Pn, Px);
    m->quadtree = new Fem2D::FQuadTree(m, Pn, Px, m->nv);
    // m->decrement();
    Add2StackOfPtr2FreeRC(stack, m);
    return m;
  }
};

// Init init;
static void Load_Init( ) {
  cout << "\n  -- lood: init MeshGenQA\n";
  Global.Add("MeshGenQA", "(", new OneOperatorCode< MeshGenQA >( ));
}

LOADFUNC(Load_Init)
