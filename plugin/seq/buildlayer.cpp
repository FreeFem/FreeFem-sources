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

// NOT-USED

#include <iostream>
#include <cfloat>
#include <cmath>
#include <complex>
using namespace std;
#include "error.hpp"
#include "AFunction.hpp"
using namespace std;
#include "rgraph.hpp"
#include "RNM.hpp"
#include "fem.hpp"

#include "FESpacen.hpp"
#include "FESpace.hpp"

#include "MatriceCreuse_tpl.hpp"
#include "MeshPoint.hpp"
#include "Operator.hpp"
#include "lex.hpp"

#include "lgfem.hpp"
#include "lgmesh3.hpp"
#include "lgsolver.hpp"
#include "problem.hpp"
#include "LayerMesh.hpp"
#include "TransfoMesh_v2.hpp"

#include <set>
#include <vector>
#include <fstream>

using namespace Fem2D;

class BuildLayeMesh_Op : public E_F0mps {
 public:
  Expression eTh;
  Expression enmax, ezmin, ezmax, xx, yy, zz;
  static const int n_name_param = 9;    //
  static basicAC_F0::name_and_type name_param[];
  Expression nargs[n_name_param];
  KN_< long > arg(int i, Stack stack, KN_< long > a) const {
    return nargs[i] ? GetAny< KN_< long > >((*nargs[i])(stack)) : a;
  }

  double arg(int i, Stack stack, double a) const {
    return nargs[i] ? GetAny< double >((*nargs[i])(stack)) : a;
  }

  long arg(int i, Stack stack, long a) const {
    return nargs[i] ? GetAny< long >((*nargs[i])(stack)) : a;
  }

 public:
  BuildLayeMesh_Op(const basicAC_F0 &args, Expression tth, Expression nmaxx)
    : eTh(tth), enmax(nmaxx), ezmin(0), ezmax(0), xx(0), yy(0), zz(0) {
    cout << "construction par BuilLayeMesh_Op" << endl;
    args.SetNameParam(n_name_param, name_param, nargs);
    const E_Array *a2 = 0, *a1 = 0;
    if (nargs[0]) {
      a1 = dynamic_cast< const E_Array * >(nargs[0]);
    }

    if (nargs[1]) {
      a2 = dynamic_cast< const E_Array * >(nargs[1]);
    }

    int err = 0;
    if (a1) {
      if (a1->size( ) != 2) {
        CompileError("LayerMesh (Th,n, zbound=[zmin,zmax],) ");
      }

      ezmin = to< double >((*a1)[0]);
      ezmax = to< double >((*a1)[1]);
    }

    if (a2) {
      if (a2->size( ) != 3) {
        CompileError("LayerMesh (Th,n, transfo=[X,Y,Z],) ");
      }

      xx = to< double >((*a2)[0]);
      yy = to< double >((*a2)[1]);
      zz = to< double >((*a2)[2]);
    }
  }

  AnyType operator( )(Stack stack) const;
};

basicAC_F0::name_and_type BuildLayeMesh_Op::name_param[] = {{"zbound", &typeid(E_Array)},
                                                            {"transfo", &typeid(E_Array)},
                                                            {"coef", &typeid(double)},
                                                            {"reftet", &typeid(KN_< long >)},
                                                            {"reffacemid", &typeid(KN_< long >)},
                                                            {"reffaceup", &typeid(KN_< long >)},
                                                            {"reffacelow", &typeid(KN_< long >)},
                                                            {"facemerge", &typeid(long)},
                                                            {"ptmerge", &typeid(double)}};

class BuildLayerMesh : public OneOperator {
 public:
  BuildLayerMesh( ) : OneOperator(atype< pmesh3 >( ), atype< pmesh >( ), atype< long >( )) {}

  E_F0 *code(const basicAC_F0 &args) const {
    cout << " je suis dans code(const basicAC_F0 & args) const" << endl;
    // cout << "args: " << args << endl;
    return new BuildLayeMesh_Op(args, t[0]->CastTo(args[0]), t[1]->CastTo(args[1]));
  }
};

AnyType BuildLayeMesh_Op::operator( )(Stack stack) const {
  MeshPoint *mp(MeshPointStack(stack)), mps = *mp;
  Mesh *pTh = GetAny< Mesh * >((*eTh)(stack));
  int nlayer = (int)GetAny< long >((*enmax)(stack));

  ffassert(pTh && nlayer > 0);
  Mesh &Th = *pTh;
  Mesh *m = pTh;       // question a quoi sert *m ??
  int nbv = Th.nv;     // nombre de sommet
  int nbt = Th.nt;     // nombre de triangles
  int neb = Th.neb;    // nombre d'aretes fontiere
  cout << " " << nbv << " " << nbv << " nbe " << neb << endl;
  KN< double > zmin(nbv), zmax(nbv);
  KN< double > clayer(nbv);    // nombre de layer est nlayer*clayer

  clayer = -1;
  zmin = 0.;
  zmax = 1.;

  for (int it = 0; it < nbt; ++it) {
    for (int iv = 0; iv < 3; ++iv) {
      int i = Th(it, iv);
      if (clayer[i] < 0) {
        mp->setP(&Th, it, iv);
        // cout << "mp: fait " << endl;
        if (ezmin) {
          zmin[i] = GetAny< double >((*ezmin)(stack));
        }

        if (ezmax) {
          zmax[i] = GetAny< double >((*ezmax)(stack));
        }

        clayer[i] = Max(0., Min(1., arg(2, stack, 1.)));
      }
    }
  }

  ffassert(clayer.min( ) >= 0);

  cout << "lecture valeur des references " << endl;

  KN< long > zzempty;
  KN< long > nrtet(arg(3, stack, zzempty));
  KN< long > nrfmid(arg(4, stack, zzempty));
  KN< long > nrfup(arg(5, stack, zzempty));
  KN< long > nrfdown(arg(6, stack, zzempty));
  int point_confondus_ok(arg(7, stack, 0L));
  double precis_mesh(arg(8, stack, -1.));

  cout << nrtet.N( ) << nrfmid.N( ) << nrfup.N( ) << nrfdown.N( ) << endl;

  ffassert(nrtet.N( ) % 2 == 0);
  ffassert(nrfmid.N( ) % 2 == 0);
  ffassert(nrfup.N( ) % 2 == 0);
  ffassert(nrfdown.N( ) % 2 == 0);

  // realisation de la map par default

  map< int, int > maptet;
  map< int, int > maptrimil, maptrizmax, maptrizmin;
  map< int, int > mapemil, mapezmax, mapezmin;

  build_layer_map_tetrahedra(Th, maptet);
  build_layer_map_triangle(Th, maptrimil, maptrizmax, maptrizmin);
  build_layer_map_edge(Th, mapemil, mapezmax, mapezmin);

  // Map utilisateur
  map< int, int >::iterator imap;

  for (int ii = 0; ii < nrtet.N( ); ii += 2) {
    imap = maptet.find(nrtet[ii]);
    if (imap != maptet.end( )) {
      imap->second = nrtet[ii + 1];
    }
  }

  for (int ii = 0; ii < nrfmid.N( ); ii += 2) {
    imap = maptrimil.find(nrfmid[ii]);
    if (imap != maptrimil.end( )) {
      imap->second = nrfmid[ii + 1];
    }
  }

  for (int ii = 0; ii < nrfup.N( ); ii += 2) {
    imap = maptrizmax.find(nrfup[ii]);
    if (imap != maptrizmax.end( )) {
      imap->second = nrfup[ii + 1];
    }
  }

  for (int ii = 0; ii < nrfdown.N( ); ii += 2) {
    imap = maptrizmin.find(nrfdown[ii]);
    if (imap != maptrizmin.end( )) {
      imap->second = nrfdown[ii + 1];
    }
  }

  int nebn = 0;
  KN< int > ni(nbv);

  for (int i = 0; i < nbv; i++) {
    ni[i] = Max(0, Min(nlayer, (int)lrint(nlayer * clayer[i])));
  }

  Mesh3 *Th3 = build_layer(Th, nlayer, ni, zmin, zmax, maptet, maptrimil, maptrizmax, maptrizmin,
                           mapemil, mapezmax, mapezmin);

  if (!(xx) && !(yy) && !(zz)) {
    Th3->BuildBound( );
    Th3->BuildAdj( );
    Th3->Buildbnormalv( );
    Th3->BuildjElementConteningVertex( );
    Th3->BuildGTree( );
    Add2StackOfPtr2FreeRC(stack, m);
    *mp = mps;
    return Th3;
  } else {

    KN< double > txx(Th3->nv), tyy(Th3->nv), tzz(Th3->nv);
    KN< int > takemesh(Th3->nv);
    MeshPoint *mp3(MeshPointStack(stack));

    takemesh = 0;
    Mesh3 &rTh3 = *Th3;

    for (int it = 0; it < Th3->nt; ++it) {
      for (int iv = 0; iv < 4; ++iv) {
        int i = (*Th3)(it, iv);
        if (takemesh[i] == 0) {
          mp3->setP(Th3, it, iv);
          if (xx) {
            txx[i] = GetAny< double >((*xx)(stack));
          }

          if (yy) {
            tyy[i] = GetAny< double >((*yy)(stack));
          }

          if (zz) {
            tzz[i] = GetAny< double >((*zz)(stack));
          }

          takemesh[i] = takemesh[i] + 1;
        }
      }
    }

    int border_only = 0;
    int recollement_elem = 0, recollement_border = 1;
    Mesh3 *T_Th3 = Transfo_Mesh3(precis_mesh, rTh3, txx, tyy, tzz, border_only, recollement_elem,
                                 recollement_border, point_confondus_ok);

    T_Th3->BuildBound( );
    T_Th3->BuildAdj( );
    T_Th3->Buildbnormalv( );
    T_Th3->BuildjElementConteningVertex( );
    T_Th3->BuildGTree( );
    // T_Th3->decrement();
    Add2StackOfPtr2FreeRC(stack, T_Th3);
    *mp = mps;
    return T_Th3;
  }
}

static void Load_Init( ) {    // le constructeur qui ajoute la fonction "splitmesh3"  a freefem++
  typedef const Mesh *pmesh;
  typedef const Mesh3 *pmesh3;

  if (verbosity) {
    cout << " load: buildlayers  " << endl;
  }

  Global.Add("buildlayers", "(", new BuildLayerMesh);
}

LOADFUNC(Load_Init)
