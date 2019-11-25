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
#include "libmeshb7.h"
#include "extern.h"
#include "string.h"
#include "eigenv.h"
#ifdef WIN32
#include <fcntl.h>
#endif

#define WrdSiz 4

static int debug = 0;
void getline_bin_float_vertex(int ddim, double *c, int *ref) {
  int i;
  float ff;
  size_t ret;

  for (i = 0; i < ddim; i++) {
    ret = fread((unsigned char *)(&ff), WrdSiz, 1, stdin);
    if (ret == 0) printf("fread error\n");
    c[i] = ff;
  }

  ret = fread((unsigned char *)ref, WrdSiz, 1, stdin);
  if (ret == 0) printf("fread error\n");
}

void getline_bin_double_vertex(int ddim, double *c, int *ref) {
  int i;
  size_t ret;

  for (i = 0; i < ddim; i++) {
    ret = fread((unsigned char *)(&c[i]), WrdSiz, 2, stdin);
    if (ret == 0) printf("fread error\n");
  }

  ret = fread((unsigned char *)ref, WrdSiz, 1, stdin);
  if (ret == 0) printf("fread error\n");
}

void getline_bin_elem(int ddim, int *v, int *ref) {
  int i;
  size_t ret;

  for (i = 0; i < ddim; i++) {
    ret = fread((unsigned char *)(&v[i]), WrdSiz, 1, stdin);
    if (ret == 0) printf("fread error\n");
  }

  ret = fread((unsigned char *)ref, WrdSiz, 1, stdin);
  if (ret == 0) printf("fread error\n");
}

void getline_bin_edge(int *v0, int *v1, int *ref) {
  size_t ret;

  ret = fread((unsigned char *)v0, WrdSiz, 1, stdin);
  if (ret == 0) printf("fread error\n");
  ret = fread((unsigned char *)v1, WrdSiz, 1, stdin);
  if (ret == 0) printf("fread error\n");
  ret = fread((unsigned char *)ref, WrdSiz, 1, stdin);
  if (ret == 0) printf("fread error\n");
}

void getline_bin_int_noref(int ddim, int *v) {
  int i;

  for (i = 0; i < ddim; i++) {
    size_t ret = fread((unsigned char *)(&v[i]), WrdSiz, 1, stdin);
    if (ret == 0) printf("fread error\n");
  }
}

void getline_bin_float_noref(int ddim, double *v) {
  int i;
  float ff;

  for (i = 0; i < ddim; i++) {
    size_t ret = fread((unsigned char *)(&ff), WrdSiz, 1, stdin);
    if (ret == 0) printf("fread error\n");
    if (debug) printf("value of ff %f\n", ff);

    v[i] = ff;
  }
}

void getline_bin_double_noref(int ddim, double *v) {
  int i;

  for (i = 0; i < ddim; i++) {
    size_t ret = fread((unsigned char *)(&v[i]), WrdSiz, 2, stdin);
    if (ret == 0) printf("fread error\n");
  }
}

/*   function for loadsol_popen   */

int tmptype;
void read_TypeSizeTyptab_bin(int *type, int *size, int *typtab) {
  int i;
  int tmpsize;
  size_t ret;

  ret = fread((unsigned char *)(&tmptype), WrdSiz, 1, stdin);
  if (ret == 0) printf("fread error\n");
  assert(tmptype <= GmfMaxTyp);

  tmpsize = 0;

  for (i = 0; i < tmptype; i++) {
    ret = fread((unsigned char *)(&typtab[i]), WrdSiz, 1, stdin);
    if (ret == 0) printf("fread error\n");
    tmpsize += typtab[i];
  }

  *size = tmpsize;
  *type = tmptype;
}

int loadMesh_popen_bin(pMesh mesh) {
  pPoint ppt;
  pEdge pr;
  pTriangle pt;
  pQuad pq;
  pTetra ptet;
  pHexa ph;
  double d, dn[3];
  float *n;
  int i, ia, ib, ref, is = 0, k, disc = 0, nn = 0, nt = 0, nq = 0;
  char *natureread;
  int loopdebug;
  int vatn[2];
  int tvatn[3];
  size_t ret;

  int retcode = 0;
  int cod;
  int KwdCod;
  int NulPos;

#ifdef WIN32
  _setmode(fileno(stdin), O_BINARY);
#endif

  /* parse keywords */
  mesh->np = 0;
  mesh->nt = 0;
  mesh->nq = 0;
  mesh->ntet = 0;
  mesh->nhex = 0;
  mesh->nc = 0;
  mesh->nr = 0;
  mesh->na = 0;
  mesh->nri = 0;
  mesh->nre = 0;
  mesh->nvn = 0;
  mesh->ntg = 0;
  mesh->ne = mesh->nt + mesh->nq + mesh->ntet + mesh->nhex;

  loopdebug = -1;

  // read code
  ret = fread((unsigned char *)&cod, WrdSiz, 1, stdin);
  if (ret == 0) printf("fread error\n");
  if (cod != 1) {
    printf("error in reading the binary file .meshb with popen\n");
    exit(1);
  }

  ret = fread((unsigned char *)&(mesh->ver), WrdSiz, 1, stdin);
  if (ret == 0) printf("fread error\n");
  ret = fread((unsigned char *)&KwdCod, WrdSiz, 1, stdin);
  if (ret == 0) printf("fread error\n");
  if (KwdCod != GmfDimension) {
    printf("error in reading the binary file .meshb with popen\n");
    exit(1);
  }

  ret = fread((unsigned char *)&NulPos, WrdSiz, 1, stdin);
  if (ret == 0) printf("fread error\n");
  ret = fread((unsigned char *)&(mesh->dim), WrdSiz, 1, stdin);
  if (ret == 0) printf("fread error\n");
  /*control of the dimension*/
  if ((mesh->dim != 2) && (mesh->dim != 3)) printf("the dimension is not correct");

  if (debug) printf("reading dimension %i \n", mesh->dim);

  while (!feof(stdin)) {
    loopdebug = loopdebug + 1;

    ret = fread((unsigned char *)&KwdCod, WrdSiz, 1, stdin);
    if (ret == 0) printf("fread error\n");

    /* determination of KwdCod */
    if (debug) printf("reading KwdCod %i %i\n", KwdCod, loopdebug);

    switch (KwdCod) {
      case GmfVertices:
        natureread = "Vertices";
        if (debug) printf("reading %s", natureread);

        ret = fread((unsigned char *)&NulPos, WrdSiz, 1, stdin);
        if (ret == 0) printf("fread error\n");
        ret = fread((unsigned char *)&(mesh->np), WrdSiz, 1, stdin);
        if (ret == 0) printf("fread error\n");
        if (debug) printf(": number of vertices %i\n", mesh->np);

        if (!mesh->np) {
          if (debug) fprintf(stdout, "  ## No vertex\n");

          retcode = -1;
          goto Lret;
        }

        if (ddebug) printf("allocate %d points\n", mesh->np);

        mesh->point = (pPoint)M_calloc(mesh->np + 1, sizeof(Point), "zaldy1.point");
        assert(mesh->point);

        for (k = 1; k <= mesh->np; k++) {
          ppt = &mesh->point[k];

          if (mesh->ver == GmfFloat)
            getline_bin_float_vertex(mesh->dim, ppt->c, &ref);
          else
            getline_bin_double_vertex(mesh->dim, ppt->c, &ref);

          ppt->ref = ref & 0x7fff;
          ppt->tag = M_UNUSED;
        }

        break;
      case GmfTriangles:
        natureread = "Triangles";
        if (debug) printf("reading %s", natureread);

        ret = fread((unsigned char *)&NulPos, WrdSiz, 1, stdin);
        if (ret == 0) printf("fread error\n");
        ret = fread((unsigned char *)&(mesh->nt), WrdSiz, 1, stdin);
        if (ret == 0) printf("fread error\n");
        if (debug) printf(": number of triangles %i\n", mesh->nt);

        if (ddebug) printf("allocate %d tria\n", mesh->nt);

        mesh->tria = (pTriangle)M_calloc(mesh->nt + 1, sizeof(Triangle), "zaldy1.tria");
        assert(mesh->tria);

        for (k = 1; k <= mesh->nt; k++) {
          pt = &mesh->tria[k];
          getline_bin_elem(3, pt->v, &ref);

          pt->ref = ref & 0x7fff;

          for (i = 0; i < 3; i++) {
            if (pt->v[i] < 1 || pt->v[i] > mesh->np) {
              disc++;
              pt->v[0] = 0;
              break;
            } else {
              ppt = &mesh->point[pt->v[i]];
              ppt->tag &= ~M_UNUSED;
            }
          }
        }

        break;
      case GmfQuadrilaterals:
        natureread = "Quadrilateral";
        if (debug) printf("reading %s", natureread);

        ret = fread((unsigned char *)&NulPos, WrdSiz, 1, stdin);
        if (ret == 0) printf("fread error\n");
        ret = fread((unsigned char *)&(mesh->nq), WrdSiz, 1, stdin);
        if (ret == 0) printf("fread error\n");
        if (debug) printf(": number of hexahedrons %i\n", mesh->nq);

        if (ddebug) printf("allocate %d quad\n", mesh->nq);

        if (mesh->nq) {
          mesh->quad = (pQuad)M_calloc(mesh->nq + 1, sizeof(Quad), "zaldy1.quad");
          assert(mesh->quad);
        }

        for (k = 1; k <= mesh->nq; k++) {
          pq = &mesh->quad[k];
          getline_bin_elem(4, pq->v, &ref);
          pq->ref = ref & 0x7fff;

          for (i = 0; i < 4; i++) {
            if (pq->v[i] < 1 || pq->v[i] > mesh->np) {
              disc++;
              pq->v[0] = 0;
              break;
            } else {
              ppt = &mesh->point[pq->v[i]];
              ppt->tag &= ~M_UNUSED;
            }
          }
        }

        break;
      case GmfTetrahedra:
        natureread = "Tetrahedra";
        if (debug) printf("reading %s", natureread);

        ret = fread((unsigned char *)&NulPos, WrdSiz, 1, stdin);
        if (ret == 0) printf("fread error\n");
        ret = fread((unsigned char *)&(mesh->ntet), WrdSiz, 1, stdin);
        if (ret == 0) printf("fread error\n");
        if (debug) printf(": number of tetrahedrons %i\n", mesh->ntet);

        if (mesh->ntet) {
          if (ddebug) printf("allocate %d tetra\n", mesh->ntet);

          mesh->tetra = (pTetra)M_calloc(mesh->ntet + 1, sizeof(Tetra), "zaldy1.tetra");
          assert(mesh->tetra);
        }

        for (k = 1; k <= mesh->ntet; k++) {
          ptet = &mesh->tetra[k];
          getline_bin_elem(4, ptet->v, &ref);
          ptet->ref = ref & 0x7fff;

          for (i = 0; i < 4; i++) {
            if (ptet->v[i] < 1 || ptet->v[i] > mesh->np) {
              disc++;
              ptet->v[0] = 0;
              break;
            } else {
              ppt = &mesh->point[ptet->v[i]];
              ppt->tag &= ~M_UNUSED;
            }
          }
        }

        break;
      case GmfHexahedra:
        natureread = "Hexahedra";
        if (debug) printf("reading %s", natureread);

        ret = fread((unsigned char *)&NulPos, WrdSiz, 1, stdin);
        if (ret == 0) printf("fread error\n");
        ret = fread((unsigned char *)&(mesh->nhex), WrdSiz, 1, stdin);
        if (ret == 0) printf("fread error\n");
        if (debug) printf(": number of hexahedrons %i\n", mesh->nhex);

        if (ddebug) printf("allocate %d hexa\n", mesh->nhex);

        mesh->hexa = (pHexa)M_calloc(mesh->nhex + 1, sizeof(Hexa), "zaldy1.hexa");
        assert(mesh->hexa);

        for (k = 1; k <= mesh->nhex; k++) {
          ph = &mesh->hexa[k];
          getline_bin_elem(8, ph->v, &ref);

          ph->ref = ref & 0x7fff;

          for (i = 0; i < 8; i++) {
            if (ph->v[i] < 1 || ph->v[i] > mesh->np) {
              disc++;
              ph->v[0] = 0;
              break;
            } else {
              ppt = &mesh->point[ph->v[i]];
              ppt->tag &= ~M_UNUSED;
            }
          }
        }

        break;
      case GmfCorners:
        ret = fread((unsigned char *)&NulPos, WrdSiz, 1, stdin);
        if (ret == 0) printf("fread error\n");
        ret = fread((unsigned char *)&(mesh->nc), WrdSiz, 1, stdin);
        if (ret == 0) printf("fread error\n");

        for (k = 1; k <= mesh->nc; k++) {
          ret = fread((unsigned char *)&is, WrdSiz, 1, stdin);
          if (ret == 0) printf("fread error\n");
          if (is < 1 || is > mesh->np) {
            disc++;
          } else {
            ppt = &mesh->point[is];
            ppt->tag |= M_CORNER;
            ppt->tag &= ~M_UNUSED;
          }
        }

        break;
      case GmfRequiredVertices:
        ret = fread((unsigned char *)&NulPos, WrdSiz, 1, stdin);
        if (ret == 0) printf("fread error\n");
        ret = fread((unsigned char *)&(mesh->nr), WrdSiz, 1, stdin);
        if (ret == 0) printf("fread error\n");

        for (k = 1; k <= mesh->nr; k++) {
          ret = fread((unsigned char *)&is, WrdSiz, 1, stdin);
          if (ret == 0) printf("fread error\n");
          if (is < 1 || is > mesh->np) {
            disc++;
          } else {
            ppt = &mesh->point[is];
            ppt->tag |= M_REQUIRED;
            ppt->tag &= ~M_UNUSED;
          }
        }

        break;
      case GmfEdges:
        natureread = "Edges";
        if (debug) printf("reading %s", natureread);

        ret = fread((unsigned char *)&NulPos, WrdSiz, 1, stdin);
        if (ret == 0) printf("fread error\n");
        ret = fread((unsigned char *)&(mesh->na), WrdSiz, 1, stdin);
        if (ret == 0) printf("fread error\n");
        if (debug) printf(": number of edges %i\n", mesh->na);

        if (ddebug) printf("allocate %d edges\n", mesh->na);

        mesh->edge = (pEdge)M_calloc(mesh->na + 1, sizeof(Edge), "zaldy1.edge");
        assert(mesh->edge);

        for (k = 1; k <= mesh->na; k++) {
          getline_bin_edge(&ia, &ib, &ref);
          if (ia < 1 || ia > mesh->np || ib < 1 || ib > mesh->np) {
            disc++;
          } else {
            pr = &mesh->edge[k];
            pr->v[0] = ia;
            pr->v[1] = ib;
            pr->ref = ref & 0x7fff;
            pr->tag = !pr->ref ? M_NOTAG : M_TAG;
            ppt = &mesh->point[ia];
            ppt->tag &= ~M_UNUSED;
            ppt = &mesh->point[ib];
            ppt->tag &= ~M_UNUSED;
          }
        }

        break;
      case GmfRidges:
        ret = fread((unsigned char *)&NulPos, WrdSiz, 1, stdin);
        if (ret == 0) printf("fread error\n");
        ret = fread((unsigned char *)&(mesh->nri), WrdSiz, 1, stdin);
        if (ret == 0) printf("fread error\n");

        for (k = 1; k <= mesh->nri; k++) {
          ret = fread((unsigned char *)&is, WrdSiz, 1, stdin);
          if (ret == 0) printf("fread error\n");
          if (is < 1 || is > mesh->na) {
            disc++;
          } else {
            pr = &mesh->edge[is];
            pr->tag |= M_RIDGE;
          }
        }

        break;
      case GmfRequiredEdges:
        ret = fread((unsigned char *)&NulPos, WrdSiz, 1, stdin);
        if (ret == 0) printf("fread error\n");
        ret = fread((unsigned char *)&(mesh->nre), WrdSiz, 1, stdin);
        if (ret == 0) printf("fread error\n");

        for (k = 2; k <= mesh->nre; k++) {
          ret = fread((unsigned char *)&is, WrdSiz, 1, stdin);
          if (ret == 0) printf("fread error\n");
          if (is < 1 || is > mesh->na) {
            disc++;
          } else {
            pr = &mesh->edge[is];
            pr->tag |= M_REQUIRED;
          }
        }

        break;
      case GmfNormals:
        ret = fread((unsigned char *)&NulPos, WrdSiz, 1, stdin);
        if (ret == 0) printf("fread error\n");
        ret = fread((unsigned char *)&(mesh->nvn), WrdSiz, 1, stdin);
        if (ret == 0) printf("fread error\n");

        /*  allocation  */
        if (!mesh->ntg) {
          mesh->extra = (pExtra)M_calloc(1, sizeof(Extra), "zaldy1.extra");
          assert(mesh->extra);
        }

        mesh->extra->n = (float *)M_calloc(3 * mesh->nvn + 1, sizeof(float), "inmesh");
        assert(mesh->extra->n);

        for (k = 1; k <= mesh->nvn; k++) {
          n = &mesh->extra->n[3 * (k - 1) + 1];

          if (mesh->ver == GmfFloat)
            getline_bin_float_noref(3, dn);
          else
            getline_bin_double_noref(3, dn);

          n[0] = dn[0];
          n[1] = dn[1];
          n[2] = dn[2];

          d = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
          if (d > 0.0) {
            d = 1.0 / sqrt(d);
            n[0] *= d;
            n[1] *= d;
            n[2] *= d;
          }
        }

        break;
      case GmfNormalAtVertices:
        if (!mesh->nvn) {
          printf("The field Normal need to be reading before %s", natureread);
          exit(1);
        }

        ret = fread((unsigned char *)&NulPos, WrdSiz, 1, stdin);
        if (ret == 0) printf("fread error\n");
        ret = fread((unsigned char *)&(mesh->extra->iv), WrdSiz, 1, stdin);
        if (ret == 0) printf("fread error\n");
        mesh->extra->nv = (int *)M_calloc(mesh->np + 1, sizeof(int), "inmesh");
        assert(mesh->extra->nv);

        for (k = 1; k <= mesh->extra->iv; k++) {
          getline_bin_int_noref(2, vatn);
          nn = vatn[0];
          is = vatn[1];
          if (nn < 1 || nn > mesh->np)
            disc++;
          else
            mesh->extra->nv[nn] = is;
        }

        break;
      case GmfNormalAtTriangleVertices:
        if (!mesh->nvn) {
          printf("The field Normal need to be reading before %s", natureread);
          exit(1);
        }

        ret = fread((unsigned char *)&NulPos, WrdSiz, 1, stdin);
        if (ret == 0) printf("fread error\n");
        ret = fread((unsigned char *)&(mesh->extra->it), WrdSiz, 1, stdin);
        if (ret == 0) printf("fread error\n");
        mesh->extra->nt = (int *)M_calloc(3 * mesh->nt + 1, sizeof(int), "inmesh");
        assert(mesh->extra->nt);

        for (k = 1; k <= mesh->extra->it; k++) {
          getline_bin_int_noref(3, tvatn);
          tvatn[0] = nt;
          tvatn[1] = is;
          tvatn[2] = nn;

          if (nt < 1 || nt > mesh->nt || is < 1 || is > 3 || nn < 1 || nn > mesh->nvn)
            disc++;
          else
            mesh->extra->nt[3 * (nt - 1) + is] = nn;
        }

        break;
      case GmfNormalAtQuadrilateralVertices:
        if (!mesh->nvn) {
          printf("The field Normal need to be reading before %s", natureread);
          exit(1);
        }

        ret = fread((unsigned char *)&NulPos, WrdSiz, 1, stdin);
        if (ret == 0) printf("fread error\n");
        ret = fread((unsigned char *)&(mesh->extra->iq), WrdSiz, 1, stdin);
        if (ret == 0) printf("fread error\n");
        mesh->extra->nq = (int *)M_calloc(4 * mesh->nq + 1, sizeof(int), "inmesh");
        assert(mesh->extra->nq);

        for (k = 1; k <= mesh->extra->iq; k++) {
          getline_bin_int_noref(3, tvatn);
          tvatn[0] = nq;
          tvatn[1] = is;
          tvatn[2] = nn;
          if (nq < 1 || nq > mesh->nq || is < 1 || is > 4 || nn < 1 || nn > mesh->nvn)
            disc++;
          else
            mesh->extra->nq[3 * (nq - 1) + is] = nn;
        }

        break;
      case GmfTangents:
        ret = fread((unsigned char *)&NulPos, WrdSiz, 1, stdin);
        if (ret == 0) printf("fread error\n");
        ret = fread((unsigned char *)&(mesh->ntg), WrdSiz, 1, stdin);
        if (ret == 0) printf("fread error\n");

        if (!mesh->nvn) {
          mesh->extra = (pExtra)M_calloc(1, sizeof(Extra), "zaldy1.extra");
          assert(mesh->extra);
        }

        mesh->extra->t = (float *)M_calloc(3 * mesh->ntg + 1, sizeof(float), "inmesh");
        assert(mesh->extra->t);

        for (k = 1; k <= mesh->ntg; k++) {
          n = &mesh->extra->t[3 * (k - 1) + 1];

          if (mesh->ver == GmfFloat)
            getline_bin_float_noref(3, dn);
          else
            getline_bin_double_noref(3, dn);

          n[0] = dn[0];
          n[1] = dn[1];
          n[2] = dn[2];

          d = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
          if (d > 0.0) {
            d = 1.0 / sqrt(d);
            n[0] *= d;
            n[1] *= d;
            n[2] *= d;
          }
        }

        break;
      case GmfTangentAtVertices:
        if (!mesh->ntg) {
          printf("The field Tangent need to be reading before %s", natureread);
          exit(1);
        }

        ret = fread((unsigned char *)&NulPos, WrdSiz, 1, stdin);
        if (ret == 0) printf("fread error\n");
        ret = fread((unsigned char *)&(mesh->extra->tv), WrdSiz, 1, stdin);
        if (ret == 0) printf("fread error\n");

        mesh->extra->tv = (int *)M_calloc(mesh->np + 1, sizeof(int), "inmesh");
        assert(mesh->extra->tv);

        for (k = 1; k <= mesh->extra->jv; k++) {
          getline_bin_int_noref(2, vatn);
          nn = vatn[0];
          is = vatn[1];

          if (nn < 1 || nn > mesh->np)
            disc++;
          else
            mesh->extra->tv[nn] = is;
        }

        break;
      case GmfTangentAtEdgeVertices:
        if (!mesh->ntg) {
          printf("The field Tangent need to be reading before %s", natureread);
          exit(1);
        }

        ret = fread((unsigned char *)&NulPos, WrdSiz, 1, stdin);
        if (ret == 0) printf("fread error\n");
        ret = fread((unsigned char *)&(mesh->extra->je), WrdSiz, 1, stdin);
        if (ret == 0) printf("fread error\n");

        mesh->extra->te = (int *)M_calloc(2 * mesh->na + 1, sizeof(int), "inmesh");
        assert(mesh->extra->te);

        for (k = 1; k <= mesh->extra->je; k++) {
          getline_bin_int_noref(3, tvatn);
          nt = tvatn[0];
          is = tvatn[1];
          nn = tvatn[2];

          if (nt < 1 || nt > mesh->np || is < 1 || is > 2 || nn < 1 || nn > mesh->ntg)
            disc++;
          else
            mesh->extra->te[3 * (nt - 1) + is] = is;
        }

        break;
      case GmfEnd:

        ret = fread((unsigned char *)&NulPos, WrdSiz, 1, stdin);
        if (ret == 0) printf("fread error\n");
        printf("End of mesh\n");
        break;
      default:
        printf("This data KwdCod it is not taken in this version\n");
        exit(1);
        break;
    }

    if (KwdCod == GmfEnd) break;
  }

  /* check if vertices and elements found */
  if (!mesh->np) {
    fprintf(stdout, "  ## No vertex\n");
    retcode = -1;
    goto Lret;
  }

  mesh->ne = mesh->nt + mesh->nq + mesh->ntet + mesh->nhex;

  if (disc > 0) fprintf(stdout, "  ## %d entities discarded\n", disc);

  retcode = 1;

Lret:
#ifdef WIN32
  _setmode(fileno(stdin), O_TEXT);
#endif
  return retcode;
}

/* function of lecture */

int loadScaVecTen_bin(pMesh mesh, int numsol, int dim, int ver, int nel, int type, int size,
                      int *typtab, int key) {
  pSolution sol;
  float fbuf[GmfMaxTyp];
  double m[6], lambda[3], eigv[3][3], vp[2][2];
  int k, i, iord, off;
  double ScaSol[1], VecSol[3], TenSol[9];
  int retcode = 0;

  if (ddebug) printf("numsol=%i, type=%i, size=%i\n", numsol, type, size);

  if (numsol > type) numsol = 1;

  numsol--;
  mesh->nbb = nel;
  mesh->bbmin = 1.0e20;
  mesh->bbmax = -1.0e20;
  if (!zaldy2(mesh)) retcode = 0;

  sol = mesh->sol;
  sol->dim = dim;
  sol->ver = ver;

  off = 0;

  for (i = 0; i < numsol; i++) {
    if (ddebug) printf("typtab[%i]=%i", i, typtab[i]);

    switch (typtab[i]) {
      case GmfSca:
        off++;
        break;
      case GmfVec:
        off += sol->dim;
        break;
      case GmfSymMat:
        off += sol->dim * (sol->dim + 1) / 2;
        break;
    }
  }

  if (ddebug) {
    printf("typtab[%i]=%i, off%i", i, typtab[i], off);
    printf("numsol=%i,typtab[i]=%i\n", numsol, typtab[i]);
  }

  fflush(stdout);

  switch (typtab[numsol]) {
    case GmfSca:
      mesh->nfield = 1;

      for (k = 1; k <= nel; k++) {
        // reading data must be a double !
        getline_bin_double_noref(1, ScaSol);
        mesh->sol[k].bb = ScaSol[0];

        if (mesh->sol[k].bb < mesh->bbmin) mesh->bbmin = mesh->sol[k].bb;

        if (mesh->sol[k].bb > mesh->bbmax) mesh->bbmax = mesh->sol[k].bb;
      }

      break;

    case GmfVec:
      if (ddebug) printf("   vector field\n");

      mesh->nfield = sol->dim;

      for (k = 1; k <= nel; k++) {
        // reading data must be a double !
        mesh->sol[k].bb = 0.0;
        getline_bin_double_noref(sol->dim, VecSol);

        for (i = 0; i < sol->dim; i++) fbuf[off + i] = VecSol[i];

        for (i = 0; i < sol->dim; i++) {
          mesh->sol[k].m[i] = fbuf[off + i];
          mesh->sol[k].bb += fbuf[off + i] * fbuf[off + i];
        }

        mesh->sol[k].bb = sqrt(mesh->sol[k].bb);
        if (mesh->sol[k].bb < mesh->bbmin) mesh->bbmin = mesh->sol[k].bb;

        if (mesh->sol[k].bb > mesh->bbmax) mesh->bbmax = mesh->sol[k].bb;
      }

      break;

    case GmfSymMat:
      if (ddebug) printf("   metric field\n");

      mesh->nfield = sol->dim * (sol->dim + 1) / 2;

      for (k = 1; k <= nel; k++) {
        // reading data must be double !!!
        getline_bin_double_noref(sol->dim * (sol->dim + 1) / 2, TenSol);

        for (i = 0; i < sol->dim * (sol->dim + 1) / 2; i++) {
          fbuf[off + i] = TenSol[i];
        }

        if (sol->dim == 2) {
          for (i = 0; i < 3; i++) {
            mesh->sol[k].m[i] = m[i] = fbuf[off + i];
          }

          iord = eigen2(m, lambda, vp);
          if (!iord) printf("eigen2 error\n");
          mesh->sol[k].bb = min(lambda[0], lambda[1]);
          if (mesh->sol[k].bb < mesh->bbmin) mesh->bbmin = mesh->sol[k].bb;

          if (mesh->sol[k].bb > mesh->bbmax) mesh->bbmax = mesh->sol[k].bb;
        } else {
          for (i = 0; i < 6; i++) {
            mesh->sol[k].m[i] = fbuf[off + i];
          }

          mesh->sol[k].m[2] = fbuf[off + 3];
          mesh->sol[k].m[3] = fbuf[off + 2];

          for (i = 0; i < 6; i++) {
            m[i] = mesh->sol[k].m[i];
          }

          iord = eigenv(1, m, lambda, eigv);
          if (iord) {
            mesh->sol[k].bb = lambda[0];
            mesh->sol[k].bb = max(mesh->sol[k].bb, lambda[1]);
            mesh->sol[k].bb = max(mesh->sol[k].bb, lambda[2]);
            if (mesh->sol[k].bb < mesh->bbmin) mesh->bbmin = mesh->sol[k].bb;

            if (mesh->sol[k].bb > mesh->bbmax) mesh->bbmax = mesh->sol[k].bb;
          }
        }
      }

      break;
  }

  retcode = 1;

  return retcode;
}

/*load solution (metric) */
int loadSol_popen_bin(pMesh mesh, char *filename, int numsol) {
  int inm = 0, key, nel, size, type, typtab[GmfMaxTyp], ver, dim;
  char *natureread;
  int KwdCod;
  int cod;
  int NulPos;
  int retcode = 0;
  size_t ret;

#ifdef WIN32
  _setmode(fileno(stdin), O_BINARY);
#endif
  // read code
  ret = fread((unsigned char *)&cod, WrdSiz, 1, stdin);
  if (ret == 0) printf("fread error\n");
  if (cod != 1) {
    printf("error in reading the binary file .meshb with popen\n");
    exit(1);
  }

  ret = fread((unsigned char *)&ver, WrdSiz, 1, stdin);
  if (ret == 0) printf("fread error\n");
  ret = fread((unsigned char *)&KwdCod, WrdSiz, 1, stdin);
  if (ret == 0) printf("fread error\n");
  if (KwdCod != GmfDimension) {
    printf("error in reading the binary file .meshb with popen\n");
    exit(1);
  }

  ret = fread((unsigned char *)&NulPos, WrdSiz, 1, stdin);
  if (ret == 0) printf("fread error\n");
  ret = fread((unsigned char *)&dim, WrdSiz, 1, stdin);
  if (ret == 0) printf("fread error\n");
  natureread = "Dimension";
  printf(".sol: %s %i (mesh)%i (lecture)%i \n", natureread, dim, mesh->dim, ver);
  /*control of the dimension*/
  if (dim != mesh->dim) {
    fprintf(stderr, "  %%%% Wrong dimension %d.\n", dim);
    retcode = 0;
    goto Lret;
  }

  while (!feof(stdin)) {
    ret = fread((unsigned char *)&KwdCod, WrdSiz, 1, stdin);
    if (ret == 0) printf("fread error\n");

    if (KwdCod == GmfSolAtVertices) {
      ret = fread((unsigned char *)&NulPos, WrdSiz, 1, stdin);
      if (ret == 0) printf("fread error\n");
      ret = fread((unsigned char *)&nel, WrdSiz, 1, stdin);
      if (ret == 0) printf("fread error\n");
      natureread = "SolAtVertices";
      if (debug) fprintf(stdout, "SolAtVertices : nel %i, mesh->np %i \n", nel, mesh->np);

      if (nel != mesh->np) {
        fprintf(stderr, "  %%%% Wrong number: %d Solutions discarded\n", nel - mesh->np);
        retcode = 0;
        goto Lret;
      }

      mesh->typage = 2;
      key = GmfSolAtVertices;

      /*  type,size,typetab  */
      read_TypeSizeTyptab_bin(&type, &size, typtab);
      if (debug) printf("sol: %s; type %i; size%i;\n", natureread, type, size);

      fflush(stdout);
      /* Reading solutions*/
      loadScaVecTen_bin(mesh, 1, dim, ver, nel, type, size, typtab, key);
    }

    if (mesh->dim == 2 && mesh->nt) {
      if (KwdCod == GmfSolAtTriangles) {
        natureread = "SolAtTriangles";
        ret = fread((unsigned char *)&NulPos, WrdSiz, 1, stdin);
        if (ret == 0) printf("fread error\n");
        ret = fread((unsigned char *)&nel, WrdSiz, 1, stdin);
        if (ret == 0) printf("fread error\n");
        if (debug) printf("SolAtTriangles : nel %d, mesh->nt %d \n", nel, mesh->nt);

        if (nel && nel != mesh->nt) {
          fprintf(stderr, "  %%%% Wrong number %d.\n", nel);
          retcode = 0;
          goto Lret;
        }

        mesh->typage = 1;
        key = GmfSolAtTriangles;

        /*  type,size,typetab  */
        read_TypeSizeTyptab_bin(&type, &size, typtab);
        printf("sol: %s; type %i; size%i;\n", natureread, type, size);

        /* Reading solutions*/
        loadScaVecTen_bin(mesh, 1, dim, ver, nel, type, size, typtab, key);
      }
    }

    if (mesh->dim == 2 && mesh->nq) {
      if (KwdCod == GmfSolAtQuadrilaterals) {
        natureread = "SolAtQuadrilaterals";
        ret = fread((unsigned char *)&NulPos, WrdSiz, 1, stdin);
        if (ret == 0) printf("fread error\n");
        ret = fread((unsigned char *)&nel, WrdSiz, 1, stdin);
        if (ret == 0) printf("fread error\n");
        if (debug) fprintf(stdout, "SolAtQuadrilaterals : nel %i, mesh->nq %i \n", nel, mesh->nq);

        if (nel && nel != mesh->nq) {
          fprintf(stderr, "  %%%% Wrong number %d.\n", nel);
          retcode = 0;
          goto Lret;
        }

        mesh->typage = 1;
        key = GmfSolAtQuadrilaterals;

        /*  type,size,typetab  */
        read_TypeSizeTyptab_bin(&type, &size, typtab);
        printf("sol: %s; type %i; size%i;\n", natureread, type, size);

        /* Reading solutions*/
        loadScaVecTen_bin(mesh, 1, dim, ver, nel, type, size, typtab, key);
      }
    }

    if (mesh->dim == 3 && mesh->ntet) {
      if (KwdCod == GmfSolAtTetrahedra) {
        natureread = "SolAtTetrahedra";
        ret = fread((unsigned char *)&NulPos, WrdSiz, 1, stdin);
        if (ret == 0) printf("fread error\n");
        ret = fread((unsigned char *)&nel, WrdSiz, 1, stdin);
        if (ret == 0) printf("fread error\n");
        if (debug) fprintf(stdout, "SolAtTetrahedra : nel %i, mesh->ntet %i \n", nel, mesh->ntet);

        if (nel && nel != mesh->ntet) {
          fprintf(stderr, "  %%%% Wrong number %d.\n", nel);
          retcode = 0;
          goto Lret;
        }

        mesh->typage = 1;
        key = GmfSolAtTetrahedra;

        /*  type,size,typetab  */
        read_TypeSizeTyptab_bin(&type, &size, typtab);
        printf("sol: %s; type %i; size%i;\n", natureread, type, size);

        /* Reading solutions */
        loadScaVecTen_bin(mesh, 1, dim, ver, nel, type, size, typtab, key);
      }
    }

    if (mesh->dim == 3 && mesh->nhex) {
      if (KwdCod == GmfSolAtHexahedra) {
        natureread = "SolAtHexahedra";
        ret = fread((unsigned char *)&NulPos, WrdSiz, 1, stdin);
        if (ret == 0) printf("fread error\n");
        ret = fread((unsigned char *)&nel, WrdSiz, 1, stdin);
        if (ret == 0) printf("fread error\n");
        if (debug) fprintf(stdout, "SolAtHexahedra : nel %d, mesh->nhex %d \n", nel, mesh->nhex);

        if (nel && nel != mesh->nhex) {
          fprintf(stderr, "  %%%% Wrong number %d.\n", nel);
          GmfCloseMesh(inm);
          retcode = 0;
          goto Lret;
        }

        mesh->typage = 1;
        key = GmfSolAtHexahedra;

        /*  type,size,typetab  */
        read_TypeSizeTyptab_bin(&type, &size, typtab);
        printf("sol: %s; type %i; size%i;\n", natureread, type, size);

        /* Reading solutions*/
        loadScaVecTen_bin(mesh, 1, dim, ver, nel, type, size, typtab, key);
      }
    }

    if (KwdCod == GmfEnd) {
      ret = fread((unsigned char *)&NulPos, WrdSiz, 1, stdin);
      if (ret == 0) printf("fread error\n");
      if (debug) printf("End of solution\n");

      if (ddebug) printf("Reading of mesh file is finished");

      break;
    }
  }

  retcode = 1;
Lret:
#ifdef WIN32
  _setmode(fileno(stdin), O_BINARY);
#endif
  return (retcode);
}

#ifdef __cplusplus
}
#endif
