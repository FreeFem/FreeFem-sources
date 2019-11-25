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

#define MAX_MORPH 100

int imstep, imreverse;

int modeMorphing( ) {
  pScene scene;
  pMesh mesh1, mesh2;
  int k;
  clock_t ct;

  /* default */
  cv.nbs = 1;
  imstep = 1;
  imreverse = 1;
  if (ddebug) printf("morphing: create window\n");

  fprintf(stdout, "\n Building scene\n");
  ct = clock( );

  /* create grafix */
  scene = cv.scene[0];
  mesh1 = cv.mesh[0];
  iniopt(scene, mesh1);
  parsop(scene, mesh1);
  meshRef(scene, mesh1);
  matSort(scene);

  mesh2 = cv.mesh[1];
  parsop(scene, mesh2);
  meshRef(scene, mesh2);

  for (k = 1; k <= mesh2->np; k++) {
    pPoint ppt1, ppt2;

    ppt1 = &mesh1->point[k];
    ppt2 = &mesh2->point[k];
    ppt2->c[0] -= ppt1->c[0];
    ppt2->c[1] -= ppt1->c[1];
    ppt2->c[2] -= ppt1->c[2];
  }

  if (!createScene(scene, 0)) {
    fprintf(stderr, "  ## Unable to create scene\n");
    return (0);
  }

  ct = difftime(clock( ), ct);
  fprintf(stdout, "  Scene seconds:      %.2f\n", (double)ct / (double)CLOCKS_PER_SEC);
  return (1);
}

int morphMesh(pScene sc, pMesh mesh1) {
  pMesh mesh2;
  int k;
  static float dt = 1.0 / MAX_MORPH;

  imstep++;
  if (imstep == 0) {
    imstep = 2;
    dt = -dt;
    if (!imreverse) {
      glutIdleFunc(NULL);
      morphing = 0;
      return (0);
    }
  } else if (imstep == MAX_MORPH + 1) {
    imstep = -MAX_MORPH + 1;
    dt = -dt;
    if (!imreverse) {
      glutIdleFunc(NULL);
      morphing = 0;
      return (0);
    }
  }

  mesh2 = cv.mesh[1];

  for (k = 1; k <= mesh1->np; k++) {
    pPoint ppt1, ppt2;

    ppt1 = &mesh1->point[k];
    ppt2 = &mesh2->point[k];
    ppt1->c[0] += dt * ppt2->c[0];
    ppt1->c[1] += dt * ppt2->c[1];
    ppt1->c[2] += dt * ppt2->c[2];
  }

  doLists(sc, mesh1);
  return (1);
}

#ifdef __cplusplus
}
#endif
