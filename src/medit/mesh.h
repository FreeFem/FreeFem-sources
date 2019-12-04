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

#ifndef _MESH_H_
#define _MESH_H_

#define M_NOTAG 0
#define M_CORNER (1 << 0)
#define M_RIDGE (1 << 1)
#define M_REQUIRED (1 << 2)
#define M_TAG (1 << 3)
#define M_UNUSED (1 << 5)

#ifndef ubyte
typedef unsigned char ubyte;
#endif
#ifndef ushort
typedef unsigned short uShort;
#endif

typedef struct spoint {
  double c[3];
  int tmp;
  short ref;
  uShort mark;
  char tag, clip, flag;
} Point;
typedef Point *pPoint;
typedef struct striangle {
  int v[3], nxt;
  short ref; /* reference */
  uShort mark, cpt;
  char clip;
} Triangle;
typedef Triangle *pTriangle;
typedef struct squad {
  int v[4], nxt;
  short ref;
  char clip;
} Quad;
typedef Quad *pQuad;
typedef struct edge {
  int v[2];
  short ref;
  char tag;
} Edge;
typedef Edge *pEdge;
typedef struct stetra {
  int v[4], nxt, mark;
  short ref;
  uShort cpt;
  char clip;
} Tetra;
typedef Tetra *pTetra;
typedef struct shexa {
  int v[8], nxt, mark;
  short ref;
  uShort cpt;
  char clip;
} Hexa;
typedef Hexa *pHexa;
typedef struct extra {
  float *n, *t;
  int *nv, *nt, *nq, *tv, *te;
  int iv, it, iq, jv, je;
} Extra;
typedef Extra *pExtra;
typedef struct solu {
  float bb;
  float *m;
  int ver, dim;
} Solution;
typedef Solution *pSolution;

/* Mesh: mesh data structure */
typedef struct mesh {
  double xmin, ymin, zmin, xmax, ymax, zmax;
  double xtra, ytra, ztra;
  float bbmin, bbmax;
  int ne, nt, nq, ntet, nhex;
  int np, nc, nr, na, nre, nri;
  int nvn, ntg, dim, ver, nbb, typage, nfield;
  uShort mark;
  char name[256], typ;
  pPoint point;
  pTriangle tria;
  pQuad quad;
  pEdge edge;
  pTetra tetra;
  pHexa hexa;
  int *adja;
  int *grid;
  ubyte *voy;
  pExtra extra;
  pSolution sol;
} Mesh;
typedef Mesh *pMesh;

#endif /* _MESH_H_ */

#ifdef __cplusplus
}
#endif
