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

#define MAX_PATH 1024
#define IMG_PER_SEC 25

/* add point to path */
int pathAdd(pScene sc, int x, int y) {
  GLdouble modelmat[16], projmat[16];
  GLdouble ptx, pty, ptz;
  GLint viewport[4];
  Trajet *path;
  pMesh mesh;
  pTransform view;
  float *p;
  int k;

  if (ddebug) fprintf(stdout, "pathAdd");

  path = &sc->path;
  mesh = cv.mesh[sc->idmesh];
  view = sc->view;
  if (path->np == MAX_PATH) return (0);

  /* memory alloc */
  if (!path->np) {
    path->pt = (float *)calloc(mesh->dim * (MAX_PATH + 1), sizeof(float));
    if (!path->pt) {
      fprintf(stderr, "   Not enough memory to store path.\n");
      return (0);
    }

    path->tg = (float *)calloc(mesh->dim * (MAX_PATH + 1), sizeof(float));
    if (!path->tg) {
      fprintf(stderr, "   Not enough memory to store tangents.\n");
      return (0);
    }
  } else if (path->np > MAX_PATH - 1) {
    return (0);
  }

  /* get 3D point */
  glGetDoublev(GL_PROJECTION_MATRIX, projmat);
  glGetIntegerv(GL_VIEWPORT, viewport);

  for (k = 0; k < 16; k++) {
    modelmat[k] = view->matrix[k];
  }

  gluUnProject(x, y, -100, modelmat, projmat, viewport, &ptx, &pty, &ptz);

  /* store point coords */
  ++path->np;
  p = &path->pt[mesh->dim * path->np];
  p[0] = ptx;
  p[1] = pty;
  p[2] = ptz;

  /* compute tangent at point-1 */
  if (path->np > 2) {
    float *p1, *a;

    p1 = &path->pt[mesh->dim * (path->np - 2)];
    a = &path->tg[mesh->dim * path->np];
    a[0] = 0.5 * (p[0] - p1[0]);
    a[1] = 0.5 * (p[1] - p1[1]);
    a[2] = 0.5 * (p[2] - p1[2]);
  }

  if (ddebug) printf(" Point %d: (%f %f %f) added\n", path->np, p[0], p[1], p[2]);

  return (1);
}

/* build list of points */
GLuint pathList(pScene sc) {
  Trajet path;
  pMesh mesh;
  GLuint dlist;
  float *p;
  int k;
  static GLfloat green[3] = {0.2, 1.0, 0.2};

  if (ddebug) fprintf(stdout, "pathList");

  path = sc->path;
  mesh = cv.mesh[sc->idmesh];

  /* build display list */
  dlist = glGenLists(1);
  glNewList(dlist, GL_COMPILE);
  if (glGetError( )) return (0);

  /* node positions */
  glPointSize(5);
  glColor3f(1.0, 0.3, 0.3);
  glBegin(GL_POINTS);

  for (k = 1; k <= path.np; k++) {
    p = &path.pt[mesh->dim * k];
    glVertex3fv(p);
  }

  glEnd( );

  /* straight path */
  glLineWidth(2.0);
  glColor3fv(green);
  glBegin(GL_LINE_STRIP);

  for (k = 1; k <= path.np; k++) {
    p = &path.pt[mesh->dim * k];
    glVertex3fv(p);
  }

  glEnd( );
  glLineWidth(1.0);

  glEndList( );
  return (dlist);
}

/* follow path */
void pathFollow(pScene sc) {
  if (!sc->path.np) return;

  if (ddebug) fprintf(stdout, "pathFollow %d\n", sc->path.np);
}

int pathLoad(char *file, pScene sc) {
  FILE *in;
  char *ptr, data[256];

  strcpy(data, file);
  ptr = (char *)strstr(data, ".mesh");
  if (ptr) *ptr = '\0';

  if (!strstr(data, ".path")) strcat(data, ".path");

  in = fopen(data, "r");
  if (!in) {
    sscanf(data, "DEFAULT.path");
    in = fopen(data, "r");
    if (!in) return (0);
  }

  if (!quiet) fprintf(stdout, "  Loading %s\n", data);

  fclose(in);
  return (1);
}

int pathSave(char *file, pScene sc) {
  FILE *out;
  pMesh mesh;
  time_t timeptr;
  int i, k;
  char *ptr, data[256];

  strcpy(data, file);
  ptr = (char *)strstr(data, ".mesh");
  if (ptr) *ptr = '\0';

  if (!strstr(data, ".path")) strcat(data, ".path");

  out = fopen(data, "w");
  if (!out) {
    fprintf(stdout, "  ## Unable to open file\n");
    ;
    return (0);
  }

  if (!quiet) fprintf(stdout, "  Writing %s\n", data);

  /* create standard path file */
  fprintf(out, "#  File created with Medit %s\n", ME_VER);
  fprintf(out, "#  Release %s\n", ME_REL);
  time(&timeptr);
  fprintf(out, "# Created: %s", ctime(&timeptr));

  mesh = cv.mesh[sc->idmesh];
  fprintf(out, "NbPoints\n%d\n", sc->path.np);

  for (k = 1; k <= sc->path.np; k++) {
    float *p;

    p = &sc->path.pt[mesh->dim * k];

    for (i = 0; i < mesh->dim; i++) {
      fprintf(out, "%f ", p[i]);
    }

    fprintf(out, "\n");
  }

  fclose(out);
  return (1);
}

#ifdef __cplusplus
}
#endif
