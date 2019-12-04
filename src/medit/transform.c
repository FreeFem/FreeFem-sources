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

void resetTransform(pTransform tr) {
  static float itransf[16] = {1., 0., 0., 0., 0., 1., 0., 0., 0., 0., 1., 0., 0., 0., 0., 1.};

  tr->pos[0] = tr->pos[1] = tr->pos[2] = 0.0f;
  tr->angle = 0.0f;
  tr->panx = tr->pany = 0.0f;
  tr->opanx = tr->opany = 0.0f;
  tr->mstate = 1;
  tr->manim = GL_FALSE;

  memcpy(tr->matrix, itransf, 16 * sizeof(float));
  memcpy(tr->rot, itransf, 16 * sizeof(float));
  memcpy(tr->tra, itransf, 16 * sizeof(float));
}

pTransform createTransform( ) {
  pTransform tr;

  /* default */
  if (ddebug) printf("create transformation\n");

  tr = (pTransform)M_calloc(1, sizeof(struct transform), "transform");
  assert(tr);

  /* set default values */
  resetTransform(tr);
  tr->mbutton = 0;

  return (tr);
}

#ifdef __cplusplus
}
#endif
