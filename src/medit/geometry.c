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

GLuint geomList(pScene sc, pMesh mesh) {
  GLuint list = 0;
  pMaterial pm;
  pPoint ppt, pp0, pp1;
  float n[3];
  int k, it = 0, nm;
  static float red[4] = {1.0, 0.0, 0.0, 1.0};
  static float yellow[4] = {1.0, 1.0, 0.0, 1.0};

  /* default */
  if (mesh->na + mesh->nc + mesh->np == 0) return (0);

  /* create display list */
  list = glGenLists(1);
  if (!list) return (0);

  glNewList(list, GL_COMPILE);
  if (glGetError( )) return (0);

  /* draw corners, ridges and required items */
  if (ddebug) printf("construct point list\n");

  if (mesh->ne) {
    glPointSize(sc->par.pointsize);
    glBegin(GL_POINTS);

    for (k = 1; k <= mesh->np; k++) {
      ppt = &mesh->point[k];
      if (ppt->tag & M_UNUSED && !ppt->ref) continue;

      if (ppt->tag == M_CORNER)
        glColor3fv(red);
      else if (ppt->tag == M_REQUIRED) {
        static float green[4] = {0.0, 1.0, 0.0, 1.0};
        glColor3fv(green);
      } else
        continue;

      it++;
      if (sc->par.linc == 1) glColor3fv(sc->par.edge);

      glVertex3f(ppt->c[0], ppt->c[1], ppt->c[2]);
    }

    glEnd( );
    glPointSize(1);
  } else {
    pm = &sc->material[DEFAULT_MAT];
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, pm->dif);
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, pm->amb);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, pm->spe);
    glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, pm->emi);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, &pm->shininess);
    glBegin(GL_POINTS);

    for (k = 1; k <= mesh->np; k++) {
      double dd;

      ppt = &mesh->point[k];
      n[0] = ppt->c[0] - sc->cx;
      n[1] = ppt->c[1] - sc->cy;
      n[2] = ppt->c[2] - sc->cz;
      dd = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
      if (dd > 0.0f) {
        dd = 1.0 / sqrt(dd);
        n[0] *= dd;
        n[1] *= dd;
        n[2] *= dd;
      }

      glNormal3fv(n);
      glVertex3f(ppt->c[0], ppt->c[1], ppt->c[2]);
    }

    glEnd( );
    it = mesh->np;
  }

  /* draw edges */
  if (ddebug) printf("construct edge list\n");

  glLineWidth(sc->par.linewidth);
  glBegin(GL_LINES);

  for (k = 1; k <= mesh->na; k++) {
    pEdge pr;

    pr = &mesh->edge[k];
    if (pr->v[0] > mesh->np || pr->v[1] > mesh->np) continue;

    if (pr->tag & M_RIDGE) {
      if (pr->tag & M_TAG)
        glColor3fv(yellow); /* ridge + ref in yellow */
      else
        glColor3fv(red); /* ridges in red */
    } else if (!pr->ref) {
      glColor3fv(sc->par.edge);
    } else {
      nm = matRef(sc, pr->ref);
      pm = &sc->material[nm];
      glColor3fv(pm->dif);
    }

    if (sc->par.linc == 1) glColor3fv(sc->par.edge);

    pp0 = &mesh->point[pr->v[0]];
    pp1 = &mesh->point[pr->v[1]];
    glVertex3f(pp0->c[0], pp0->c[1], pp0->c[2]);
    glVertex3f(pp1->c[0], pp1->c[1], pp1->c[2]);
    it++;
  }

  glEnd( );
  glLineWidth(1.0);
  glEndList( );

  if (it == 0) {
    glDeleteLists(list, 1);
    return (0);
  } else {
    return (list);
  }
}

#ifdef __cplusplus
}
#endif
