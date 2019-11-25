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

extern void drawHUD(pScene);
extern int refmat, imstep;

void initTexture(void) {
  PPMimage *imgtex;
  int typimg;
  ubyte iquiet;
  GLuint texname;

  iquiet = quiet;
  quiet = 1;
  imgtex = loadPPM("Lions.ppm", &typimg);
  quiet = iquiet;
  if (!imgtex) return;

  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
  glGenTextures(1, &texname);
  glBindTexture(GL_TEXTURE_2D, texname);

  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, imgtex->sizeX, imgtex->sizeY, 0, GL_RGB, GL_UNSIGNED_BYTE,
               imgtex->data);
  free(imgtex);
}

void backTexture(pScene sc) {
  glEnable(GL_TEXTURE_2D);
  glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  glPushAttrib(GL_ENABLE_BIT);
  glEnable(GL_TEXTURE_2D);
  glDisable(GL_LIGHTING);
  glBegin(GL_QUADS);
  glTexCoord2f(0.0f, 0.0f);
  glVertex2f(0.0, 0.0);

  glTexCoord2f(sc->par.xs, 0.0);
  glVertex2f(sc->par.xs, 0.0);

  glTexCoord2f(sc->par.xs, sc->par.ys);
  glVertex2f(sc->par.xs, sc->par.ys);

  glTexCoord2f(0.0, sc->par.ys);
  glVertex2f(0.0, sc->par.ys);
  glEnd( );
  glPopAttrib( );

  glDisable(GL_TEXTURE_2D);
}

void redrawStatusBar(pScene sc) {
  pClip clip = sc->clip;
  pMesh mesh = cv.mesh[sc->idmesh];

  if (sc->par.xs < 100) return;

  if (ddebug) fprintf(stdout, "redrawStatusBar\n");

  glDisable(GL_DEPTH_TEST);
  glDisable(GL_LIGHTING);

  glMatrixMode(GL_PROJECTION);
  glPushMatrix( );
  glLoadIdentity( );
  gluOrtho2D(1., sc->par.xs, 1., sc->par.ys);

  glMatrixMode(GL_MODELVIEW);
  glPushMatrix( );
  glLoadIdentity( );

  /* other info */
  glColor3f(1.0 - sc->par.back[0], 1.0 - sc->par.back[1], 1.0 - sc->par.back[2]);
  if (animate && !(sc->isotyp & S_PARTICLE)) {
    float frame, elpms;
    static float fps = 0.0, lastfr = 0.0;
    static int nfr = 0, pps = 0;

    nfr++;
    frame = glutGet(GLUT_ELAPSED_TIME) / 1000.0f;
    elpms = frame - lastfr;
    if (elpms > 0.999) {
      fps = nfr / elpms;
      if (mesh->nt + mesh->nq)
        pps = fps * (mesh->nt + mesh->nq);
      else if (mesh->ntet + mesh->nhex)
        pps = fps * (mesh->ntet + mesh->nhex);
      else
        pps = fps * mesh->np;

      nfr = 0;
      lastfr = frame;
    }

    output2(15, 8, "Fps: %6.2f Pps: %8d", fps, pps);
  }

  if (option == MORPHING) output2(15, 28, "%d", abs(imstep));

  if (sc->isotyp & S_STREAML && sc->par.maxtime < FLT_MAX)
    output2(15, 8, "t= %8.3f", sc->par.cumtim);
  else if (sc->isotyp & S_PARTICLE)
    output2(15, 8, "t= %8.3f", sc->par.cumtim);

  /* clip eqn */
  if (clip->active & C_ON && !(clip->active & C_HIDE)) {
    double dd;
    char buf[128], tmpbuf[128];

    sprintf(buf, "Eqn: ");
    strcpy(tmpbuf, buf);
    if (fabs(clip->eqn[0]) > EPS) sprintf(buf, "%s %+.2gx", tmpbuf, clip->eqn[0]);

    if (fabs(clip->eqn[1]) > EPS) sprintf(buf, "%s %+.2gy", tmpbuf, clip->eqn[1]);

    if (fabs(clip->eqn[2]) > EPS) sprintf(buf, "%s %+.2gz", tmpbuf, clip->eqn[2]);

    strcpy(tmpbuf, buf);
    dd = clip->eqn[3] - clip->eqn[0] * mesh->xtra - clip->eqn[1] * mesh->ytra -
         clip->eqn[2] * mesh->ztra;
    if (dd) sprintf(buf, "%s %+.2g", tmpbuf, dd);

    if (sc->par.xs > 180) output2(150, 8, "%s = 0", buf);
  }

  if ((sc->picklist) && (sc->par.xs > 390) && (!sc->isotyp) & (S_PARTICLE))
    output2(350, 8, "%15s", sc->material[refmat].name);

  if (sc->persp->pmode == PERSPECTIVE && sc->item & S_PALETTE) drawPalette(sc);

  if (sc->persp->pmode == CAMERA) drawHUD(sc);

  glPopMatrix( );

  glMatrixMode(GL_PROJECTION);
  glPopMatrix( );
  glEnable(GL_LIGHTING);
  glEnable(GL_DEPTH_TEST);
}

void mouseStatus(int button, int state, int x, int y) {
  pScene sc = cv.scene[currentScene( )];
  pTransform view = sc->view;
  ubyte axis = X_AXIS;

  /* default */
  if (ddebug) printf("control mouse %d\n", state);

  if (button == GLUT_LEFT_BUTTON) {
    if (x < 16 && x > 5)
      axis = X_AXIS;
    else if (x < 26 && x > 15)
      axis = Y_AXIS;
    else if (x < 36 && x > 25)
      axis = Z_AXIS;

    switch (axis) {
      case X_AXIS:
        view->angle = 90.0;
        break;
      case Y_AXIS:
        view->angle = 90.0;
        view->axis[0] = 1.0;
        view->axis[1] = view->axis[2] = 0.0f;
        break;
      case Z_AXIS:
        view->angle = 90.0;
        view->axis[1] = 0.0f;
        view->axis[0] = view->axis[2] = 0.0f;
        break;
    }
  }
}

#ifdef __cplusplus
}
#endif
