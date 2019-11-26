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
// AUTHORS : Pierre Jolivet
// E-MAIL  : pierre.jolivet@enseeiht.fr

/* clang-format off */
//ff-c++-LIBRARY-dep: [mpi ptscotch] scotch
//ff-c++-cpp-dep:
/* clang-format on */

#ifndef _ALL_IN_ONE_
// add for missing  def of int32_t
#include <stdint.h>
#include "ff++.hpp"
#include <vector>
#include <cmath>
#endif
#ifdef WITH_PTSCOTCH
#include <mpi.h>
#include <ptscotch.h>
#else
#include <scotch.h>
#endif

template< class T, class V, class K >
class SCOTCH_Op : public E_F0mps {
 public:
  Expression partition;
  Expression Th;
  Expression lpartition;
  static const int n_name_param = 1;
  static basicAC_F0::name_and_type name_param[];
  Expression nargs[n_name_param];
  SCOTCH_Op(const basicAC_F0 &args, Expression param1, Expression param2, Expression param3)
    : partition(param1), Th(param2), lpartition(param3) {
    args.SetNameParam(n_name_param, name_param, nargs);
  }

  AnyType operator( )(Stack stack) const;
};

template< class T, class V, class K >
basicAC_F0::name_and_type SCOTCH_Op< T, V, K >::name_param[] = {{"weight", &typeid(KN< long > *)}};

template< class T, class V, class K >
class SCOTCH : public OneOperator {
 public:
  SCOTCH( ) : OneOperator(atype< K >( ), atype< KN< K > * >( ), atype< V >( ), atype< long >( )) {}

  E_F0 *code(const basicAC_F0 &args) const {
    return new SCOTCH_Op< T, V, K >(args, t[0]->CastTo(args[0]), t[1]->CastTo(args[1]),
                                    t[2]->CastTo(args[2]));
  }
};

template< class T, class V, class K >
AnyType SCOTCH_Op< T, V, K >::operator( )(Stack stack) const {
  T *pTh = GetAny< T * >((*Th)(stack));

  ffassert(pTh);
  int nt = pTh->nt;
  KN< K > *part = GetAny< KN< K > * >((*partition)(stack));
  ffassert(part);

  int nve = T::RdHat::d + 1;
  long lpart = GetAny< long >((*lpartition)(stack));
  ffassert(lpart > 1 && part->n == nt && lpart < nt);

  KN< long > *weight = nargs[0] ? GetAny< KN< long > * >((*nargs[0])(stack)) : (KN< long > *)0;

  SCOTCH_Graph GraphSCOTCH;
  SCOTCH_Strat StratSCOTCH;
  SCOTCH_graphInit(&GraphSCOTCH);
  SCOTCH_Num baseval = 0;
  SCOTCH_Num vertnbr = nt;
  SCOTCH_Num edgenbr;
  SCOTCH_Num *verttab = new SCOTCH_Num[vertnbr + 1];
  vector< SCOTCH_Num > edgevec;
  edgevec.reserve(T::Rd::d * (vertnbr - 1));
  int cptNode = 0;
  int accum = 0;
  verttab[cptNode++] = baseval;

  for (int it = 0; it < nt; ++it) {
    for (int jt = 0; jt < nve; ++jt) {
      int jtt = jt, itt = pTh->ElementAdj(it, jtt);
      if (itt != it && itt >= 0) {
        ++accum;
        edgevec.push_back(baseval + itt);
      }
    }

    verttab[cptNode++] = accum;
  }

  edgenbr = accum;
  SCOTCH_Num *edgetab = &edgevec[0];
  SCOTCH_Num *velotab;

  if (weight) {
    velotab = new SCOTCH_Num[nt];

    for (int i = 0; i < nt; ++i) {
      velotab[i] = (SCOTCH_Num)(*weight)[i];
    }
  } else {
    velotab = NULL;
  }

  SCOTCH_Num *vendtab = NULL;
  SCOTCH_Num *vlbltab = NULL;
  SCOTCH_Num *edlotab = NULL;
  SCOTCH_graphBuild(&GraphSCOTCH, baseval, vertnbr, verttab, vendtab, velotab, vlbltab, edgenbr,
                    edgetab, edlotab);
#ifdef DEBUG
  SCOTCH_graphCheck(&GraphSCOTCH);
#endif
  KN< SCOTCH_Num > epart(nt);
  SCOTCH_stratInit(&StratSCOTCH);
  SCOTCH_stratGraphMapBuild(&StratSCOTCH, SCOTCH_STRATSPEED, lpart, 0.05);
  SCOTCH_graphPart(&GraphSCOTCH, lpart, &StratSCOTCH, epart);
  SCOTCH_graphExit(&GraphSCOTCH);
  SCOTCH_stratExit(&StratSCOTCH);
  *part = epart;
  delete[] verttab;
  if (velotab) {
    delete[] velotab;
  }

  return 0L;
}

#ifndef _ALL_IN_ONE_
static void Init_Scotch( ) {
  Global.Add("scotch", "(", new SCOTCH< const Mesh, pmesh, long >);
  Global.Add("scotch", "(", new SCOTCH< const Mesh3, pmesh3, long >);
  Global.Add("scotch", "(", new SCOTCH< const MeshS, pmeshS, long >);
  Global.Add("scotch", "(", new SCOTCH< const Mesh, pmesh, double >);
  Global.Add("scotch", "(", new SCOTCH< const Mesh3, pmesh3, double >);
  Global.Add("scotch", "(", new SCOTCH< const MeshS, pmeshS, double >);
}

LOADFUNC(Init_Scotch)
#endif
