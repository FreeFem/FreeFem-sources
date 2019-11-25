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

void scissorScene( ) {
  pScene sc;
  pPersp p;
  pTransform view;
  int width, height;

  /* default */
  if (ddebug) printf("enable scissoring\n");

  sc = cv.scene[currentScene( )];
  view = sc->view;
  p = sc->persp;

  /* subdivide main window */
  glViewport(0, 0, sc->par.xs, sc->par.ys);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity( );
  gluOrtho2D(0, sc->par.xs, 0, sc->par.ys);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity( );
  glClearColor(sc->par.back[0], sc->par.back[1], sc->par.back[2], sc->par.back[3]);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  /* split window */
  glDisable(GL_LIGHTING);
  glColor4fv(sc->par.line);
  glBegin(GL_LINES);
  glVertex2i(sc->par.xs / 2, 0);
  glVertex2i(sc->par.xs / 2, sc->par.ys);
  glVertex2i(0, sc->par.ys / 2);
  glVertex2i(sc->par.xs, sc->par.ys / 2);
  glEnd( );
  glEnable(GL_LIGHTING);

  width = (sc->par.xs + 1) / 2;
  height = (sc->par.ys + 1) / 2;

  glDisable(GL_LIGHTING);
  glColor4fv(sc->par.line);
  output2(5, sc->par.ys / 2 + 5, "Top");
  output2(5, 5, "Front");
  output2(sc->par.xs / 2 + 5, 5, "Right");

  glEnable(GL_SCISSOR_TEST);

  /* draw top right : normal */
  glViewport(width, height, width, height);
  glScissor(width, height, width, height);
  farclip(1);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity( );
  gluLookAt(0., 0., -p->depth, 0., 0., 0., 0.0, 1.0, 0.0);
  setupView(sc);
  glMultMatrixf(view->matrix);
  glTranslatef(sc->cx, sc->cy, sc->cz);
  drawModel(sc);
  if (sc->type & S_DECO) redrawStatusBar(sc);

  /* draw top left : top view */
  glViewport(0, height, width, height);
  glScissor(0, height, width, height);
  farclip(1);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity( );
  gluLookAt(0., 0., -p->depth, 0., 0., 0., 0.0, 1.0, 0.0);
  setupView(sc);
  glTranslatef(sc->cx, sc->cy, sc->cz);
  glRotatef(-90, 0.0, 0.0, 1.0);
  glTranslatef(view->opanx, 0., 0.);
  drawModel(sc);

  /* draw bottom left : front */
  glViewport(0, 0, width, height);
  glScissor(0, 0, width, height);
  farclip(1);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity( );
  gluLookAt(0., 0., -p->depth, 0., 0., 0., 0.0, 1.0, 0.0);
  setupView(sc);
  glTranslatef(sc->cx, sc->cy, sc->cz);
  glRotatef(-90.0, 1.0, 0.0, 0.0);
  glTranslatef(view->opanx, 0., 0.);
  drawModel(sc);

  /* draw bottom right : right */
  glViewport(width, 0, width, height);
  glScissor(width, 0, width, height);
  farclip(1);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity( );
  gluLookAt(0., 0., -p->depth, 0., 0., 0., 0.0, 1.0, 0.0);
  setupView(sc);
  glRotatef(-90.0, 1.0, 0.0, 0.0);
  glRotatef(-90.0, 0.0, 0.0, 1.0);
  glTranslatef(0., view->opany, 0.);
  drawModel(sc);

  glutSwapBuffers( );
  glDisable(GL_SCISSOR_TEST);
  glViewport(0, 0, sc->par.xs, sc->par.ys);

  if (saveimg) keyFile('H', 0, 0);
}

#ifdef __cplusplus
}
#endif
