/****************************************************************************/
/* This file is part of FreeFEM.                                            */
/*  ... bloc LGPLv3, cf. src/femlib/fem.hpp ...                             */
/****************************************************************************/
// SUMMARY : Scotch mesh partitioning via the face dual graph
// LICENSE : LGPLv3
// AUTHORS : Pierre-Loïc Bacq (suivant l'algorithme original, de Pierre Jolivet dans plugin/seq/scotch.cpp)

#ifndef PARTITIONSCOTCH_HPP_
#define PARTITIONSCOTCH_HPP_

#include <stdio.h>    
#include <stdint.h>   
#include <vector>
#include "PartitionDual.hpp"
#include <scotch.h>

//   nparts : nombre de parties demande
//   epart  : sortie, taille Th.nt, allouee par l'appelant
//   weight : entree optionnelle (velotab), taille Th.nt, ou nullptr
//   retour : 0 si succes, non nul si Scotch echoue
template<class Mesh>
int ffScotchPartFaceDual(const Mesh& Th, SCOTCH_Num nparts,
                         SCOTCH_Num* epart,
                         const long* weight = nullptr) {
  const int nt = Th.nt;

  for (int k = 0; k < nt; ++k) epart[k] = 0;
  if (nparts <= 1 || nt <= 0) return 0;

  const SCOTCH_Num baseval = 0;   // le helper produit des indices bruts : baseval doit rester 0
  const SCOTCH_Num vertnbr = nt;
  std::vector<SCOTCH_Num> verttab, edgevec;
  ffBuildFaceDualCSR<Mesh, SCOTCH_Num>(Th, 0, nt, verttab, edgevec);
  const SCOTCH_Num edgenbr = (SCOTCH_Num)edgevec.size();
  SCOTCH_Num* edgetab = edgevec.empty() ? nullptr : edgevec.data();

  std::vector<SCOTCH_Num> velovec;
  SCOTCH_Num* velotab = nullptr;
  if (weight) {
    velovec.resize(nt);
    for (int i = 0; i < nt; ++i) velovec[i] = (SCOTCH_Num)weight[i];
    velotab = velovec.data();
  }

  SCOTCH_Graph GraphSCOTCH;
  SCOTCH_graphInit(&GraphSCOTCH);
  int ierr = SCOTCH_graphBuild(&GraphSCOTCH, baseval, vertnbr, verttab.data(),
                               NULL, velotab, NULL,
                               edgenbr, edgetab, NULL);
  if (ierr == 0) {
    SCOTCH_Strat StratSCOTCH;
    SCOTCH_stratInit(&StratSCOTCH);
    SCOTCH_stratGraphMapBuild(&StratSCOTCH, SCOTCH_STRATSPEED, nparts, 0.05);
    ierr = SCOTCH_graphPart(&GraphSCOTCH, nparts, &StratSCOTCH, epart);
    SCOTCH_stratExit(&StratSCOTCH);
  }
  SCOTCH_graphExit(&GraphSCOTCH);
  return ierr;
}

#endif  // PARTITIONSCOTCH_HPP_
