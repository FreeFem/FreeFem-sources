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
// SUMMARY  :   add interface with partionning library METSI and ParMETIS
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
// AUTHORS : Frederic Hecht, P. Jolivet
// E-MAIL  : frederic.hecht@sorbonne-universite.fr
// E-MAIL  :  P. Jolivet <pierre.jolivet@enseeiht.fr>

/* clang-format off */
//ff-c++-LIBRARY-dep: metis
//ff-c++-cpp-dep:
/* clang-format on */

#include <ff++.hpp>
#include <cmath>
typedef KNM< double > *pRnm;
typedef KN< double > *pRn;
typedef string *pstring;
extern "C" {
#include <metis.h>
}

#ifdef METIS_VER_MAJOR
// METIS_PartMeshDual(&ne, &nn, elmnts, &etype, &numflag, &nparts, &edgecut, epart, npart);
extern "C" {
real_t libmetis__ComputeElementBalance(idx_t ne, idx_t nparts, idx_t *where);
}
#else
typedef idxtype idx_t;
#endif

template< class FESPACE, int NO, typename R >
KN< R > *partmetis( KN< R > *const &part, FESPACE *const &pVh, long const &lparts) {
  ffassert(pVh);
  const FESPACE &Vh(*pVh);
   int nve = Vh[0].NbDoF( );
  const typename FESPACE::Mesh & Th = Vh.Th;
    idx_t nt = Th.nt, nv = Vh.NbOfDF;
  

  KN< idx_t > eptr(nt + 1), elmnts(nve * nt), epart(nt), npart(nv);
  if(lparts > 1) {
      for (idx_t k = 0, i = 0; k < nt; ++k) {
          eptr[k] = i;

          for (idx_t j = 0; j < nve; j++) {
              elmnts[i++] = Vh(k, j);
          }

          eptr[k + 1] = i;
      }

      idx_t numflag = 0;
      idx_t nparts = lparts;
      idx_t edgecut;
      idx_t etype = nve - 2;    // triangle or tet .  change FH fevr 2010
      idx_t ncommon = 1;
#ifdef METIS_VER_MAJOR
      if (NO == 0) {
          METIS_PartMeshNodal(&nt, &nv, eptr, (idx_t *)elmnts, 0, 0, &nparts, 0, 0, &edgecut,
                  (idx_t *)epart, (idx_t *)npart);
      } else {
          METIS_PartMeshDual(&nt, &nv, eptr, (idx_t *)elmnts, 0, 0, &ncommon, &nparts, 0, 0, &edgecut,
                  (idx_t *)epart, (idx_t *)npart);
      }

      if (verbosity) {
          printf("  --metisOA: %d-way Edge-Cut: %7d, Balance: %5.2f Nodal=0/Dual %d\n", nparts, nve,
                  libmetis__ComputeElementBalance(nt, nparts, epart), NO);
      }

#else
      if (NO == 0) {
          METIS_PartMeshNodal(&nt, &nv, elmnts, &etype, &numflag, &nparts, &edgecut, epart, npart);
      } else {
          METIS_PartMeshDual(&nt, &nv, elmnts, &etype, &numflag, &nparts, &edgecut, epart, npart);
      }

      if (verbosity) {
          printf("  --metis: %d-way Edge-Cut: %7d, Balance: %5.2f Nodal=0/Dual %d\n", nparts, nve,
                  ComputeElementBalance(nt, nparts, epart), NO);
      }

#endif
  } else epart = 0;
  part->resize(nv);
  *part = npart;
  return part;
}

template< class Mesh, int NO, typename R >
KN< R > *partmetis(Stack s, KN< R > *const &part, Mesh *const &pTh, long const &lparts) {
  ffassert(pTh);
  const Mesh &Th(*pTh);
  idx_t nt = Th.nt, nv = Th.nv;
  idx_t nve = Mesh::RdHat::d + 1;

  KN< idx_t > eptr(nt + 1), elmnts(nve * nt), epart(nt), npart(nv);
  if(lparts > 1) {
      for (idx_t k = 0, i = 0; k < nt; ++k) {
          eptr[k] = i;

          for (idx_t j = 0; j < nve; j++) {
              elmnts[i++] = Th(k, j);
          }

          eptr[k + 1] = i;
      }

      idx_t numflag = 0;
      idx_t nparts = lparts;
      idx_t edgecut;
      idx_t etype = nve - 2;    // triangle or tet .  change FH fevr 2010
      idx_t ncommon = 1;
#ifdef METIS_VER_MAJOR
      if (NO == 0) {
          METIS_PartMeshNodal(&nt, &nv, eptr, (idx_t *)elmnts, 0, 0, &nparts, 0, 0, &edgecut,
                  (idx_t *)epart, (idx_t *)npart);
      } else {
          METIS_PartMeshDual(&nt, &nv, eptr, (idx_t *)elmnts, 0, 0, &ncommon, &nparts, 0, 0, &edgecut,
                  (idx_t *)epart, (idx_t *)npart);
      }

      if (verbosity) {
          printf("  --metisOA: %d-way Edge-Cut: %7d, Balance: %5.2f Nodal=0/Dual %d\n", nparts, nve,
                  libmetis__ComputeElementBalance(nt, nparts, epart), NO);
      }

#else
      if (NO == 0) {
          METIS_PartMeshNodal(&nt, &nv, elmnts, &etype, &numflag, &nparts, &edgecut, epart, npart);
      } else {
          METIS_PartMeshDual(&nt, &nv, elmnts, &etype, &numflag, &nparts, &edgecut, epart, npart);
      }

      if (verbosity) {
          printf("  --metis: %d-way Edge-Cut: %7d, Balance: %5.2f Nodal=0/Dual %d\n", nparts, nve,
                  ComputeElementBalance(nt, nparts, epart), NO);
      }

#endif
  } else epart = 0;
  part->resize(nt);
  *part = epart;
  return part;
}

KN< long > *partmetisd(Stack s, KN< long > *const &part, Mesh *const &pTh, long const &lparts) {
  ffassert(pTh);
  const Mesh &Th(*pTh);
  idx_t nt = Th.nt, nv = Th.nv;
  idx_t nve = Mesh::Element::NbV;

  KN< idx_t > elmnts(nve * nt), epart(nt), npart(nv);
  if(lparts > 1) {
      for (idx_t k = 0, i = 0; k < nt; ++k) {
          for (idx_t j = 0; j < nve; j++) {
              elmnts[i++] = Th(k, j);
          }
      }

      idx_t numflag = 0;
      idx_t nparts = lparts;
      idx_t edgecut;
#ifdef METIS_VER_MAJOR
      printf("  %d-way Edge-Cut: %7d, Balance: %5.2f\n", nparts, nve,
              libmetis__ComputeElementBalance(nt, nparts, epart));
#else
      printf("  %d-way Edge-Cut: %7d, Balance: %5.2f\n", nparts, nve,
              ComputeElementBalance(nt, nparts, epart));
#endif
  } else epart = 0;
  part->resize(nt);
  *part = epart;
  return part;
}

template<typename pf3r,int NO>
double metisFE( pf3r const & uij,long const &npar){
    //  typedef typename v_fes::pfes pfes;
  //   typedef typename v_fes::FESpace FESpace;
   //  typedef pair<pf3rbase,int> pf3r ;
    typedef typename pf3r::first_type pf3rbase;
    typedef typename remove_pointer<typename pf3r::first_type>::type FEbase;
    typedef typename FEbase::FESpace FESpace;
    typedef typename FESpace::Mesh Mesh;
    typedef typename FEbase::pfes pfes;
 
   // typedef pf3r.pfes;
  //  typedef   v_fes3::FESpace FESpace;
    typedef double K;
    
    pf3rbase buij=uij.first;
   
    int comp =uij.second;
    cout << " composant  "<< comp << endl;
    KN<K> * pux=buij->x();
    FESpace *pVh = uij.first->newVh( );
    ffassert(pVh);
    if(pux ==0 || pux->N() != pVh->NbOfDF) {
        cout << "  FE create or recreate " << pux <<  endl;
        if(pux) delete [] pux;
        *uij.first = pux = new KN< K >(pVh->NbOfDF);
        *pux = K( );
    }
    cout << " nbdot " << pux->N() << endl;
    FESpace& Vh= *pVh;
    int nbdofK = Vh[0].NbDoF( );
    const Mesh & Th = Vh.Th;
    cout << " Th nt "<< Th.nt << " ns "<< Th.nv << endl;
    // template< class FESPACE, int NO, typename R >
    // KN< R > *partmetis( KN< R > *const &part, FESPACE *const &pVh, long const &lparts)
    partmetis<FESpace,NO,K> (pux,pVh,npar);
    return 0;
}

static void Load_Init( ) {
  if (verbosity && mpirank == 0)
#ifdef METIS_VER_MAJOR
  {
    cout << " load: init metis (v  " << METIS_VER_MAJOR << " )\n";
  }

#else
  {
    cout << " load: init metis (v 4 )\n";
  }
#endif
  Global.Add(
    "metisnodal", "(",
    new OneOperator3_< KN< long > *, KN< long > *, const Mesh *, long,
                       E_F_stackF0F0F0_< KN< long > *, KN< long > *, const Mesh *, long > >(
      &partmetis< const Mesh, 0 >));
  Global.Add(
    "metisdual", "(",
    new OneOperator3_< KN< long > *, KN< long > *, const Mesh *, long,
                       E_F_stackF0F0F0_< KN< long > *, KN< long > *, const Mesh *, long > >(
      &partmetis< const Mesh, 1 >));
  Global.Add(
    "metisnodal", "(",
    new OneOperator3_< KN< long > *, KN< long > *, const Mesh3 *, long,
                       E_F_stackF0F0F0_< KN< long > *, KN< long > *, const Mesh3 *, long > >(
      &partmetis< const Mesh3, 0 >));
  Global.Add(
    "metisdual", "(",
    new OneOperator3_< KN< long > *, KN< long > *, const Mesh3 *, long,
                       E_F_stackF0F0F0_< KN< long > *, KN< long > *, const Mesh3 *, long > >(
      &partmetis< const Mesh3, 1 >));
  Global.Add(
    "metisnodal", "(",
    new OneOperator3_< KN< long > *, KN< long > *, const MeshS *, long,
                       E_F_stackF0F0F0_< KN< long > *, KN< long > *, const MeshS *, long > >(
      &partmetis< const MeshS, 0 >));
  Global.Add(
    "metisdual", "(",
    new OneOperator3_< KN< long > *, KN< long > *, const MeshS *, long,
                       E_F_stackF0F0F0_< KN< long > *, KN< long > *, const MeshS *, long > >(
      &partmetis< const MeshS, 1 >));
  Global.Add(
    "metisnodal", "(",
    new OneOperator3_< KN< long > *, KN< long > *, const MeshL *, long,
                       E_F_stackF0F0F0_< KN< long > *, KN< long > *, const MeshL *, long > >(
      &partmetis< const MeshL, 0 >));
  Global.Add(
    "metisdual", "(",
    new OneOperator3_< KN< long > *, KN< long > *, const MeshL *, long,
                       E_F_stackF0F0F0_< KN< long > *, KN< long > *, const MeshL *, long > >(
      &partmetis< const MeshL, 1 >));

  Global.Add(
    "metisnodal", "(",
    new OneOperator3_< KN< double > *, KN< double > *, const Mesh *, long,
                       E_F_stackF0F0F0_< KN< double > *, KN< double > *, const Mesh *, long > >(
      &partmetis< const Mesh, 0 >));
  Global.Add(
    "metisdual", "(",
    new OneOperator3_< KN< double > *, KN< double > *, const Mesh *, long,
                       E_F_stackF0F0F0_< KN< double > *, KN< double > *, const Mesh *, long > >(
      &partmetis< const Mesh, 1 >));
  Global.Add(
    "metisnodal", "(",
    new OneOperator3_< KN< double > *, KN< double > *, const Mesh3 *, long,
                       E_F_stackF0F0F0_< KN< double > *, KN< double > *, const Mesh3 *, long > >(
      &partmetis< const Mesh3, 0 >));
  Global.Add(
    "metisdual", "(",
    new OneOperator3_< KN< double > *, KN< double > *, const Mesh3 *, long,
                       E_F_stackF0F0F0_< KN< double > *, KN< double > *, const Mesh3 *, long > >(
      &partmetis< const Mesh3, 1 >));
  Global.Add(
    "metisnodal", "(",
    new OneOperator3_< KN< double > *, KN< double > *, const MeshS *, long,
                       E_F_stackF0F0F0_< KN< double > *, KN< double > *, const MeshS *, long > >(
      &partmetis< const MeshS, 0 >));
  Global.Add(
    "metisdual", "(",
    new OneOperator3_< KN< double > *, KN< double > *, const MeshS *, long,
                       E_F_stackF0F0F0_< KN< double > *, KN< double > *, const MeshS *, long > >(
      &partmetis< const MeshS, 1 >));
  Global.Add(
    "metisnodal", "(",
    new OneOperator3_< KN< double > *, KN< double > *, const MeshL *, long,
                       E_F_stackF0F0F0_< KN< double > *, KN< double > *, const MeshL *, long > >(
      &partmetis< const MeshL, 0 >));
  Global.Add(
    "metisdual", "(",
    new OneOperator3_< KN< double > *, KN< double > *, const MeshL *, long,
                       E_F_stackF0F0F0_< KN< double > *, KN< double > *, const MeshL *, long > >(
      &partmetis< const MeshL, 1 >));
    Global.Add("metisnodal","(",new OneOperator2_<double,pf3r,long>(metisFE<pf3r,0>));
    Global.Add("metisdual","(",new OneOperator2_<double,pf3r,long>(metisFE<pf3r,1>));

    Global.Add("metisnodal","(",new OneOperator2_<double,pfSr,long>(metisFE<pfSr,0>));
    Global.Add("metisdual","(",new OneOperator2_<double,pfSr,long>(metisFE<pfSr,1>));

    Global.Add("metisnodal","(",new OneOperator2_<double,pfLr,long>(metisFE<pfLr,0>));
    Global.Add("metisdual","(",new OneOperator2_<double,pfLr,long>(metisFE<pfLr,1>));

    Global.Add("metisnodal","(",new OneOperator2_<double,pfer,long>(metisFE<pfer,0>));
    Global.Add("metisdual","(",new OneOperator2_<double,pfer,long>(metisFE<pfer,1>));

}

LOADFUNC(Load_Init)
