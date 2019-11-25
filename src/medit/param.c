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

#ifndef ON
#define ON 1
#define OFF 0
#endif

extern void ortho2D(pScene, ubyte);

static void parMotion(int x, int y) {
  pScene sc;

  sc = cv.scene[currentScene( )];
  glEnable(GL_COLOR_LOGIC_OP);
  glLogicOp(GL_XOR);

  glColor3ub(255, 255, 0);
  setFont("helvetica", 10);
  drwstr(10, sc->par.ys - 120, "Vector length");
  glColor3ub(0, 255, 128);
  drwstr(150, sc->par.ys - 120, "%g", 10.1);
  glFlush( );
  glDisable(GL_COLOR_LOGIC_OP);
}

static void parMouse(int button, int state, int x, int y) {
  pScene sc;

  if (button != GLUT_LEFT_BUTTON) return;

  sc = cv.scene[currentScene( )];

  if (state == GLUT_DOWN) {
    glColor3ub(0, 255, 128);
    glDrawBuffer(GL_FRONT);
    ortho2D(sc, ON);
    glutMotionFunc(parMotion);
  } else {
    glDrawBuffer(GL_BACK);
    ortho2D(sc, OFF);
    glutMotionFunc(parMotion);
  }
}

void parEdit(pScene sc) { glutMouseFunc(parMouse); }

#ifdef __cplusplus
}
#endif
