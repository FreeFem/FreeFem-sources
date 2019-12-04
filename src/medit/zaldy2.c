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
/* SUMMARY : ...                                                            */
/* LICENSE : LGPLv3                                                         */
/* ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE           */
/* AUTHORS : Pascal Frey                                                    */
/* E-MAIL  : pascal.frey@sorbonne-universite.fr                             */

#ifdef __cplusplus
extern "C" {
#endif

#include "medit.h"
#include "extern.h"
#include "sproto.h"

int zaldy2(pMesh mesh) {
  int k, nbf;

  /* memory alloc. */
  mesh->sol = (pSolution)M_calloc(mesh->nbb + 1, sizeof(struct solu), "zaldy2");
  assert(mesh->sol);

  if (mesh->nfield == 1) return (1);

  if (mesh->nfield == mesh->dim)
    nbf = mesh->dim; /* vector field */
  else
    nbf = mesh->dim * (mesh->dim + 1) / 2; /* d*d matrix */

  for (k = 1; k <= mesh->nbb; k++) {
    pSolution ps;

    ps = &mesh->sol[k];
    ps->m = (float *)malloc(nbf * sizeof(float));
  }

  return (1);
}

#ifdef __cplusplus
}
#endif
