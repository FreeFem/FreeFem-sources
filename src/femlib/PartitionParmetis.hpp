/****************************************************************************/
/* This file is part of FreeFEM.                                            */
/*  ... bloc LGPLv3, cf. src/femlib/fem.hpp ...                             */
/****************************************************************************/
// SUMMARY : ParMETIS mesh partitioning via the distributed face dual graph
// LICENSE : LGPLv3
// AUTHORS : Pierre-Loïc Bacq (suivant plugin/mpi/parmetis.cpp de P. Jolivet)

#ifndef PARTITIONPARMETIS_HPP_
#define PARTITIONPARMETIS_HPP_

#include <vector>
#include <mpi.h>
#include "PartitionDual.hpp"

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wzero-as-null-pointer-constant"
#elif defined(__GNUC__) || defined(__GNUG__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wzero-as-null-pointer-constant"
#endif
#include <parmetis.h>
#if defined(__clang__)
  #pragma clang diagnostic pop
#elif defined(__GNUC__) || defined(__GNUG__)
  #pragma GCC diagnostic pop
#endif

//   nparts : nombre de parties demande
//   epart  : sortie, taille Th.nt, allouee par l'appelant ; identique sur tous
//            les rangs de comm en sortie
//   retour : 0 si succes, code METIS non nul sinon
// Prerequis : Th.nt >= taille de comm (sinon un rang recoit 0 element et
// ParMETIS_V3_PartKway echoue) ; garanti par l'appelant.
template<class Mesh>
int ffParmetisPartFaceDual(const Mesh &Th, int nparts, int *epart, MPI_Comm comm) {
  const int nt = Th.nt;

  for (int k = 0; k < nt; ++k) epart[k] = 0;
  if (nparts <= 1 || nt <= 0) return 0;

  int worker = 1, rank = 0;
  MPI_Comm_size(comm, &worker);
  MPI_Comm_rank(comm, &rank);

  // decoupage equilibre des nt elements sur les worker rangs
  std::vector<idx_t> vtxdist(worker + 1);
  for (int i = 0; i <= worker; ++i)
    vtxdist[i] = (idx_t)((long long)i * nt / worker);
  const int kdeb = (int)vtxdist[rank], kfin = (int)vtxdist[rank + 1];

  // tranche locale du graphe dual, adjncy en numerotation globale
  std::vector<idx_t> xadj, adjncy;
  ffBuildFaceDualCSR<Mesh, idx_t>(Th, kdeb, kfin, xadj, adjncy);

  idx_t wgtflag = 0, numflag = 0, ncon = 1, np = nparts, edgecut = 0;
  idx_t options[3] = {0, 0, 0};                 // options[0] == 0 : valeurs par defaut
  std::vector<real_t> tpwgts((size_t)ncon * (size_t)np, (real_t)1.0 / (real_t)np);
  std::vector<real_t> ubvec((size_t)ncon, (real_t)1.05);
  std::vector<idx_t> localpart((size_t)(kfin - kdeb), 0);

  int rc = ParMETIS_V3_PartKway(vtxdist.data(), xadj.data(),
                                adjncy.empty() ? nullptr : adjncy.data(),
                                nullptr, nullptr, &wgtflag, &numflag, &ncon, &np,
                                tpwgts.data(), ubvec.data(), options,
                                &edgecut, localpart.data(), &comm);
  if (rc != METIS_OK) return rc;

  // recollement : chaque rang possede exactement sa tranche, comptes connus
  std::vector<idx_t> allpart((size_t)nt);
  std::vector<int> counts(worker), displs(worker);
  for (int i = 0; i < worker; ++i) {
    displs[i] = (int)vtxdist[i];
    counts[i] = (int)(vtxdist[i + 1] - vtxdist[i]);
  }
  MPI_Allgatherv(localpart.data(), counts[rank], IDX_T,
                 allpart.data(), counts.data(), displs.data(), IDX_T, comm);

  for (int k = 0; k < nt; ++k) epart[k] = (int)allpart[k];
  return 0;
}

#endif  // PARTITIONPARMETIS_HPP_
