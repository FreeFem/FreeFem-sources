/****************************************************************************/
/* This file is part of FreeFEM.                                            */
/*  ... bloc LGPLv3, cf. src/femlib/fem.hpp ...                             */
/****************************************************************************/
// SUMMARY : Face dual graph (CSR) of a mesh element slice
// LICENSE : LGPLv3
// AUTHORS : Pierre-Loïc Bacq

#ifndef PARTITIONDUAL_HPP_
#define PARTITIONDUAL_HPP_

#include <cstddef>
#include <vector>

// Graphe dual "face" de la tranche d'elements [kdeb, kfin), au format CSR base 0.
//   xadj   : taille (kfin - kdeb + 1), offsets dans adjncy
//   adjncy : numeros GLOBAUX des elements voisins
// Les faces de bord (ElementAdj < 0) et les boucles sont ignorees.
template<class Mesh, class Idx>
void ffBuildFaceDualCSR(const Mesh &Th, int kdeb, int kfin,
                        std::vector<Idx> &xadj, std::vector<Idx> &adjncy) {
  const int nea = Th.nea;   // 4 en Tet, 3 en TriaS, 2 en SegL
  xadj.clear();
  adjncy.clear();
  xadj.reserve((size_t)(kfin - kdeb) + 1);
  adjncy.reserve((size_t)nea * (size_t)(kfin - kdeb));
  xadj.push_back(0);
  for (int k = kdeb; k < kfin; ++k) {
    for (int j = 0; j < nea; ++j) {
      int jj = j;
      int kk = Th.ElementAdj(k, jj);
      if (kk >= 0 && kk != k) adjncy.push_back((Idx)kk);
    }
    xadj.push_back((Idx)adjncy.size());
  }
}

#endif  // PARTITIONDUAL_HPP_
