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

extern int refmat, reftype, refitem;
int *pilmat, ipilmat;
int refpick = -1;
int menukey = 0; /* 1: keyScene called by a menu entry (see menuKey in menus.c) */
int oldkey = 0;  /* 1: the old action of medit for this key (menu entries, see menuOld) */
GLfloat dazim = 1.0;
GLfloat delev = 1.0;
extern void mouseCamera(int button, int state, int x, int y);
extern void motionCamera(int x, int y);
extern GLboolean picking;
static void usage( ) {
  fprintf(stdout, "\n");
  fprintf(stdout, "-- Medit: OnLine Help (the keys of ffglut) --\n");
  fprintf(stdout, "Enter - next plot (FreeFEM)  p - previous plot   Esc - quit   ? - help\n");
  fprintf(stdout, "f - fill (solution in colors)   o - iso lines   n/N - fewer/more iso lines\n");
  fprintf(stdout, "m - mesh        B - border      l - lighting    L - facets / wire\n");
  fprintf(stdout, "v - palette     P - value under the mouse       j - graduated axes\n");
  fprintf(stdout, "3 - 2d elevation   h/H - lower/higher elevation   k - complex part\n");
  fprintf(stdout, "c - 2d cut line u(x) / u(y) / none   < > (or , .) move it   shift+drag: segment\n");
  fprintf(stdout, "s - iso-surfaces (3d)   V - arrows   a/A - shorter/longer arrows\n");
  fprintf(stdout, "e/E - exploded view (shrink)   T - normals   U - overlay of the windows\n");
  fprintf(stdout, "F1 clip plane  F2 edit it (mouse)  F3 freeze it on the object  F4 vclip\n");
  fprintf(stdout, "= - reset view   z/Z - zoom   +/- scale object   arrows - translate\n");
  fprintf(stdout, "w - hardcopy PPM   b - back color   G - grid   r - refresh\n");
  fprintf(stdout, "1d graph: k lines/points  G grid  F bold  x y log  z Z zoom  A axes  i auto\n");
  fprintf(stdout, "the other actions of medit: right button (menus)\n");
}

/* special keys CAMERA mode */
void specCamera(pScene sc, int key) {
  pCamera c;
  double dd, azim, elev;
  GLfloat axe[3];
  int keyact;

  c = sc->camera;
  keyact = glutGetModifiers( );

  axe[0] = c->speed[2];
  axe[1] = 0.0;
  axe[2] = -c->speed[0];
  dd = sqrt(axe[0] * axe[0] + axe[2] * axe[2]);
  if (dd != 0.0f) {
    axe[0] /= dd;
    axe[2] /= dd;
  }

  switch (key) {
    case GLUT_KEY_LEFT:
      if (keyact & GLUT_ACTIVE_SHIFT) {
        c->eye[0] += axe[0] * c->altinc;
        c->eye[2] += axe[2] * c->altinc;
        reshapeScene(sc->par.xs, sc->par.ys);
        glutPostRedisplay( );
        return;
      } else {
        azim = Azimuth(c);
        elev = Elevation(c);
        azim += dazim;
      }

      break;
    case GLUT_KEY_RIGHT:
      if (keyact & GLUT_ACTIVE_SHIFT) {
        c->eye[0] -= axe[0] * c->altinc;
        c->eye[2] -= axe[2] * c->altinc;
        reshapeScene(sc->par.xs, sc->par.ys);
        glutPostRedisplay( );
        return;
      } else {
        azim = Azimuth(c);
        elev = Elevation(c);
        azim -= dazim;
      }

      break;
    case GLUT_KEY_UP:
      if (keyact & GLUT_ACTIVE_SHIFT) {
        c->eye[1] += c->altinc;
        reshapeScene(sc->par.xs, sc->par.ys);
        glutPostRedisplay( );
        return;
      } else {
        azim = Azimuth(c);
        elev = Elevation(c);
        elev -= delev;
      }

      break;
    case GLUT_KEY_DOWN:
      if (keyact & GLUT_ACTIVE_SHIFT) {
        c->eye[1] -= c->altinc;
        reshapeScene(sc->par.xs, sc->par.ys);
        glutPostRedisplay( );
        return;
      } else {
        azim = Azimuth(c);
        elev = Elevation(c);
        elev += delev;
      }

      break;
    default:
      return;
  }

  updateCamera(sc, c, azim, elev);

  /* refresh scene */
  reshapeScene(sc->par.xs, sc->par.ys);
  glutPostRedisplay( );
}

/* change center of scene */
static void changeCenter(pScene sc, pMesh mesh) {
  pPoint p0;
  pTriangle pt;
  pQuad pq;
  pTetra ptt;
  pHexa ph;
  float cx, cy, cz;
  int i;

  if (!refitem) return;

  cx = cy = cz = 0.0;

  switch (reftype) {
    case LPoint:
      p0 = &mesh->point[refitem];
      cx = p0->c[0];
      cy = p0->c[1];
      cz = p0->c[2];
      break;
    case LTria:
      pt = &mesh->tria[refitem];

      for (i = 0; i < 3; i++) {
        p0 = &mesh->point[pt->v[i]];
        cx += 0.33 * p0->c[0];
        cy += 0.33 * p0->c[1];
        cz += 0.33 * p0->c[2];
      }

      break;
    case LQuad:
      pq = &mesh->quad[refitem];

      for (i = 0; i < 4; i++) {
        p0 = &mesh->point[pq->v[i]];
        cx += 0.25 * p0->c[0];
        cy += 0.25 * p0->c[1];
        cz += 0.25 * p0->c[2];
      }

      break;
    case LTets:
      ptt = &mesh->tetra[refitem];

      for (i = 0; i < 4; i++) {
        p0 = &mesh->point[ptt->v[i]];
        cx += 0.25 * p0->c[0];
        cy += 0.25 * p0->c[1];
        cz += 0.25 * p0->c[2];
      }

      break;
    case LHexa:
      ph = &mesh->hexa[refitem];

      for (i = 0; i < 8; i++) {
        p0 = &mesh->point[ph->v[i]];
        cx += 0.125 * p0->c[0];
        cy += 0.125 * p0->c[1];
        cz += 0.125 * p0->c[2];
      }

      break;
  }

  /* reset translation and move to center */
  sc->view->panx = sc->view->pany = 0.0;
  sc->cx = -cx;
  sc->cy = -cy;
  sc->cz = -cz;
}

/* special keys PERSPECTIVE mode */
void special(int key, int x, int y) {
  pTransform view;
  pScene sc;
  pMesh mesh;
  pClip clip;
  pCube cube;
  float pancoeff = 0.1f;
  int keyact, idw = currentScene( );
  ubyte post = TRUE;

  /* default */
  if (ddebug) printf("special key  %d\n", key);

  sc = cv.scene[idw];

  /* special mode camera */
  if (sc->persp->pmode == CAMERA) {
    specCamera(sc, key);
    return;
  }

  view = sc->view;
  mesh = cv.mesh[sc->idmesh];
  /* 1d graph: the arrows translate the graph */
  if (graph1dMesh(mesh) && graphSpecial(sc, key)) return;

  clip = sc->clip;
  cube = sc->cube;
  keyact = glutGetModifiers( );
  if (keyact & GLUT_ACTIVE_SHIFT) pancoeff = 0.01;

  switch (key) {
    case GLUT_KEY_F1: /* toggle clipping plane */
      if (mesh->dim == 3) {
        if (cube->active & C_ON)
          keyCube('C', 0, 0);
        else
          keyClip('C', 0, 0);
      }

      post = FALSE;
      break;
    case GLUT_KEY_F2: /* edit clipping plane */
      if (mesh->dim == 3) {
        if (cube->active & C_ON)
          keyCube('E', 0, 0);
        else if (clip->active & C_ON)
          keyClip('E', 0, 0);
      }

      post = FALSE;
      break;
    case GLUT_KEY_F3: /* freeze clipping plane */
      if (cube->active & C_ON)
        keyCube('F', 0, 0);
      else if (clip->active & C_ON)
        keyClip('F', 0, 0);

      post = FALSE;
      break;
    case GLUT_KEY_F4: /* toggle Vclip */
      if (mesh->dim == 3) keyClip('Z', 0, 0);

      post = FALSE;
      break;
    case GLUT_KEY_F5: /* Toggle Shrink */
      if (mesh->ntet + mesh->nhex > 0)
        keyFeature('V', 0, 0);
      else
        keyFeature('S', 0, 0);

      break;
    case GLUT_KEY_F6: /* Increase Shrink */
      keyFeature('I', 0, 0);
      break;
    case GLUT_KEY_F7: /* Decrease Shrink */
      keyFeature('i', 0, 0);
      break;

    case GLUT_KEY_LEFT: /* translate eyes or object */
      if (clip->active & C_EDIT) {
        clip->cliptr->panx -= 0.02 * sc->dmax;
        clip->cliptr->angle = 0.0;
        clip->active |= C_UPDATE + C_REDO;
      } else if (keyact & GLUT_ACTIVE_CTRL) {
        sc->par.eyesep *= 0.9;
        printf("eyesep %f\n", sc->par.eyesep);
      } else {
        view->panx -= pancoeff * sc->dmax;
      }

      break;
    case GLUT_KEY_RIGHT:
      if (clip->active & C_EDIT) {
        clip->cliptr->panx += 0.02 * sc->dmax;
        clip->cliptr->angle = 0.0;
        clip->active |= C_UPDATE + C_REDO;
      } else if (keyact & GLUT_ACTIVE_CTRL) {
        sc->par.eyesep *= 1.1;
        printf("eyesep %f\n", sc->par.eyesep);
      } else {
        view->panx += pancoeff * sc->dmax;
      }

      break;
    case GLUT_KEY_UP:
      if (clip->active & C_EDIT) {
        clip->cliptr->pany += 0.02 * sc->dmax;
        clip->cliptr->angle = 0.0;
        clip->active |= C_UPDATE + C_REDO;
      } else {
        view->pany += pancoeff * sc->dmax;
      }

      break;
    case GLUT_KEY_DOWN:
      if (clip->active & C_EDIT) {
        clip->cliptr->pany -= 0.02 * sc->dmax;
        clip->cliptr->angle = 0.0;
        clip->active |= C_UPDATE + C_REDO;
      } else {
        view->pany -= pancoeff * sc->dmax;
      }

      break;
    default:
      return;
  }

  if (post) glutPostRedisplay( );
}

/* complex solution (key k): real part, imaginary part, modulus, argument */
static void cplxKey(pScene sc, pMesh mesh) {
  extern int cplxpart;
  extern const char *cplxname[4];
  char title[512];
  double delta;
  int j;

  cplxpart = (cplxpart + 1) % 4;
  cplxSet(mesh);
  delta = mesh->bbmax - mesh->bbmin;
  for (j = 0; j < MAXISO; j++) sc->iso.val[j] = mesh->bbmin + j * delta / (MAXISO - 1);

  setupPalette(sc, mesh);
  resetLists(sc, mesh);
  doLists(sc, mesh);
  sc->glist = geomList(sc, mesh);
  if (sc->mode & (S_MAP | S_ALTITUDE)) doMapLists(sc, mesh, 1);

  if (sc->isotyp) doIsoLists(sc, mesh, 1);

  fprintf(stdout, "  complex solution: %s\n", cplxname[cplxpart]);
  snprintf(title, sizeof(title), "Medit - [%s, %s]", mesh->name, cplxname[cplxpart]);
  glutSetWindowTitle(title);
}

/* the display at the opening, as plot(u) in ffglut: for a 2d function the colored iso lines and
   the border of the domain (no mesh, no filling), for a 3d function the colors on the surface of
   the domain (no mesh), for a vector field the arrows */
void ffDefaults(pScene sc, pMesh mesh) {
  if (!mesh->nbb || graph1dMesh(mesh)) return;

  sc->isotyp &= ~(S_ISOLINE | S_VECTOR | S_ISOSURF);
  if (mesh->dim == 2) {
    sc->mode = S_FILL; /* the interior in the color of the background: no edge of the mesh */
    sc->item |= S_GEOM; /* the border of the domain */
    if (mesh->nfield == 1 && mesh->typage == 2)
      sc->isotyp |= S_ISOLINE;
    else if (mesh->nfield == mesh->dim)
      sc->isotyp |= S_VECTOR;
  } else if (mesh->ntet + mesh->nhex > 0) {
    if (mesh->nfield == 1) {
      sc->mode = S_FILL + S_MAP; /* the solution in colors on the surface, no mesh */
    } else if (mesh->nfield == mesh->dim) {
      sc->mode = HIDDEN; /* the mesh of the surface, the arrows inside */
      sc->isotyp |= S_VECTOR;
    }
  }
}

/* the keys of ffglut (same letter, same action): 1 if the key is used */
static int ffKey(pScene sc, pMesh mesh, unsigned char key) {
  extern int probeon;
  extern GLfloat altcoef;
  extern int isoncol;
  int i;

  switch (key) {
    case 'f': /* fill: the solution in colors (the iso lines stay) */
      if (!mesh->nbb) return (1);

      sc->mode |= S_FILL;
      keyMetric('m', 0, 0);
      /* the iso lines: dark on the colors, in color without them */
      for (i = 0; i < MAX_LIST; i++) {
        if (sc->ilist[i]) glDeleteLists(sc->ilist[i], 1);

        sc->ilist[i] = 0;
      }

      if (sc->isotyp) doIsoLists(sc, mesh, 0);

      return (1);
    case 'm': /* the mesh */
      sc->mode ^= S_BDRY;
      break;
    case 'B': /* the border of the domain */
      if (sc->glist) sc->item ^= S_GEOM;

      break;
    case 'o': /* iso lines */
      if (mesh->nbb) keyMetric('l', 0, 0);

      return (1);
    case 'n': /* fewer or more iso lines */
    case 'N':
      isoncol = key == 'n' ? max(1, isoncol - 1) : min(20, isoncol + 1);
      for (i = 0; i < MAX_LIST; i++) {
        if (sc->ilist[i]) glDeleteLists(sc->ilist[i], 1);

        sc->ilist[i] = 0;
      }

      if (sc->isotyp) doIsoLists(sc, mesh, 0);

      break;
    case 's': /* 3d: iso-surfaces */
      if (mesh->nbb && mesh->ntet + mesh->nhex > 0 && mesh->nfield == 1) keyMetric('s', 0, 0);

      return (1);
    case 'v': /* the values of the colors (palette) */
      if (mesh->nbb) keyMetric('p', 0, 0);

      return (1);
    case 'P': /* the value under the mouse */
      probeon = !probeon;
      break;
    case 'j': /* graduated axes */
      sc->axes = !sc->axes;
      break;
    case '3': /* 2d: elevation (3d view of the 2d function) */
      if (mesh->dim == 2 && mesh->nbb) keyMetric('k', 0, 0);

      return (1);
    case 'h': /* height of the elevation */
    case 'H':
      if (!(sc->mode & S_ALTITUDE)) return (1);

      altcoef *= key == 'H' ? 1.2 : 1.0 / 1.2;
      doMapLists(sc, mesh, 1);
      if (sc->isotyp) doIsoLists(sc, mesh, 1);

      break;
    case 'k': /* complex solution: real part, imaginary part, modulus, argument */
      if (mesh->cplx) cplxKey(sc, mesh);

      break;
    case 'l': /* lighting */
      sc->mode ^= S_COLOR;
      break;
    case 'e': /* exploded view: the elements smaller */
    case 'E':
      if (key == 'e')
        sc->shrink = sc->shrink > 0.99 ? 0.95 : max(0.1, sc->shrink - 0.05);
      else if (sc->shrink < 0.99)
        sc->shrink = sc->shrink + 0.05 > 0.96 ? 1.0 : sc->shrink + 0.05;

      doLists(sc, mesh);
      if (sc->mode & S_MAP) doMapLists(sc, mesh, 1);

      break;
    case '=': /* reset the view */
      keyView('R', 0, 0);
      return (1);
    case 'w': /* hardcopy (PPM) */
    case 'W':
      keyFile('H', 0, 0);
      return (1);
    case 'r': /* refresh */
      break;
    case 'T': /* the normals */
      if (mesh->nvn) keyItem('N', 0, 0);

      return (1);
    case 'V': /* the arrows of a vector field (no key in ffglut: the arrows are always drawn) */
      if (mesh->nbb) keyMetric('w', 0, 0);

      return (1);
    case 'L': /* faces filled or not (no key in ffglut) */
      sc->mode ^= S_FILL;
      break;
    /* the old keys of medit which are other keys of ffglut: nothing */
    case 'c': /* 2d: cut line, the 1d graph of the solution along the line (see cutKey) */
    case '<':
    case '>':
    case ',':
    case '.':
      cutKey(sc, key);
      return (1);
    case 'g':
    case 'i':
    case 'I':
    case 'q':
    case 'Q':
    case 'R':
    case 'C':
    case 't':
    case 'x':
    case 'X':
    case 'y':
    case 'Y':
    case 'F':
    case 'O':
    case '[':
    case ']':
    case '*':
    case '2':
      return (1);
    default:
      return (0);
  }

  glutPostRedisplay( );
  return (1);
}

void keyScene(unsigned char key, int x, int y) {
  pMaterial pm;
  pTetra ptt;
  pHexa ph;
  pTriangle pt;
  pQuad pq;
  pScene sc, sc1;
  pMesh mesh;
  pClip clip;
  pCube cube;
  pPersp p;
  pCamera cam;
  int ret, k, keyact, numit, idw = currentScene( );
  ubyte post = FALSE, dolist = FALSE;

  if (idw < 0) exit(0);

#ifndef WIN32
  /* server mode: n (or return) = next plot, p = previous plot */
  if (srvKey(key)) return;
#endif

  /* ESC = end medit, ? = help (q and h are other keys of ffglut) */
  if (key == 27)
    exit(0);
  else if (key == '?')
    usage( );

  /* default */
  sc = cv.scene[idw];
  mesh = cv.mesh[sc->idmesh];
  clip = sc->clip;
  cube = sc->cube;
  p = sc->persp;

  /* 1d graph: its own keys */
  if (graph1dMesh(mesh) && graphKey(sc, key)) return;

  if (key == 'U') {
    /* overlay: the data of the other scenes (same mesh) in this scene */
    sc->overlay = !sc->overlay;
    doOverlay(sc);
    glutPostRedisplay( );
    return;
  }

  /* the same keys as ffglut (not for the old actions of medit asked by a menu entry) */
  if (!oldkey && ffKey(sc, mesh, key)) return;

  /* a key sent by a menu entry (see menuKey): no modifier, glutGetModifiers is only valid
     in the keyboard and mouse callbacks */
  keyact = menukey ? 0 : glutGetModifiers( );
  if (key == ' ') {
    if (option == MORPHING) {
      morphMesh(sc, mesh);
    } else if (sc->isotyp & S_PARTICLE) {
      glutIdleFunc(0);
    } else {
      cam = sc->camera;
      cam->eye[0] += cam->spmod * cam->speed[0];
      cam->eye[1] += cam->spmod * cam->speed[1];
      cam->eye[2] += cam->spmod * cam->speed[2];
      reshapeScene(sc->par.xs, sc->par.ys);
    }

    post = TRUE;
  } else if (islower(key)) {
    switch (key) {
      case 'a': /* shorter arrows (as in ffglut; the animation is in the menu) */
        vecLength(sc, mesh, 1.0 / 1.2);
        post = TRUE;
        break;
      case 'b': /* backcolor */
        keyColor('b', 0, 0);
        break;
      case 'c':
        if (keyact & GLUT_ACTIVE_ALT) {
          keyView('C', 0, 0);
        } else {
          sc->mode ^= S_COLOR;
          post = TRUE;
        }

        break;
      case 'd':
        if (keyact & GLUT_ACTIVE_ALT) keyView('D', 0, 0);

        break;
      case 'e':
        keyColor('e', 0, 0);
        break;
      case 'f':
        sc->mode ^= S_FILL;
        post = TRUE;
        break;
      case 'g':
        keyItem('g', 0, 0);
        break;
      case 'h':
        usage( );
        break;
      case 'i':
        keyView('R', 0, 0);
        if (clip->active & C_ON) resetClip(sc, clip, mesh);

        if (cube->active & C_ON) resetCube(sc, cube, mesh);

        sc1 = sc;

        while (sc1->slave > -1) {
          sc1 = cv.scene[sc1->slave];
          glutSetWindow(sc1->idwin);
          keyScene('i', 0, 0);
        }

        glutSetWindow(sc->idwin);
        break;
      case 'j':
        sc->type ^= S_DECO;
        post = TRUE;
        break;
      case 'k':
        keyMetric('k', 0, 0);
        break;
      case 'l':
        if (keyact & GLUT_ACTIVE_ALT)
          keyView('L', 0, 0);
        else
          sc->mode ^= S_BDRY;

        post = TRUE;
        break;
      case 'm': /* toggle metric */
        if (mesh->nbb) keyMetric('m', 0, 0);

        keyMode('n', 0, 0);
        break;
      case 'o': /* iso-lines */
        if (mesh->nbb) keyMetric('l', 0, 0);

        break;
      case 'p':
        if (keyact & GLUT_ACTIVE_ALT)
          keyView('P', 0, 0);
        else
          keyMetric('p', 0, 0);

        break;

      case 'r':
        if (refmat < 0 || ipilmat == sc->par.nbmat) return;

        pilmat[++ipilmat] = refmat;
        pm = &sc->material[refmat];
        pm->flag = 1;
        updatePoints(sc, mesh, refmat);
        if (sc->picklist) glDeleteLists(sc->picklist, 1);

        sc->picklist = 0;
        refmat = -1;
        dolist = TRUE;
        post = TRUE;
        break;
      case 's':
        if (!refitem) return;

        switch (reftype) {
          case LTria:
            pt = &mesh->tria[refitem];
            pt->v[0] = 0;
            break;
          case LQuad:
            pq = &mesh->quad[refitem];
            pq->v[0] = 0;
            break;
          case LTets:
            ptt = &mesh->tetra[refitem];
            ptt->v[0] = 0;
            break;
          case LHexa:
            ph = &mesh->hexa[refitem];
            ph->v[0] = 0;
            break;
        }

        sc->picklist = 0;
        dolist = TRUE;
        post = TRUE;
        break;
      case 't': /* iso-surfaces (3d scalar solution) */
        if (mesh->nbb && mesh->ntet + mesh->nhex > 0 && mesh->nfield == 1) keyMetric('s', 0, 0);

        break;
      case 'u':
        if (keyact & GLUT_ACTIVE_ALT) keyView('U', 0, 0);

        break;
      case 'v':
        keyMetric('v', 0, 0);
        break;
      case 'w':
        keyMetric('w', 0, 0);
        break;
      case 'x':
        if (cube->active & C_EDIT) {
          cube->cma[0] += 0.1 * sc->dmax;
          post = TRUE;
          break;
        }

      case 'y':
        if (cube->active & C_EDIT) {
          cube->cma[1] += 0.1 * sc->dmax;
          post = TRUE;
          break;
        }

        if (p->pmode & CAMERA) {
          cam = sc->camera;
          cam->vecup = (cam->vecup + 1) % 3;

          switch (cam->vecup) {
            case X_AXIS:
              cam->altinc = (mesh->xmax - mesh->xmin);
              break;
            case Y_AXIS:
              cam->altinc = (mesh->ymax - mesh->ymin);
              break;
            case Z_AXIS:
              cam->altinc = (mesh->zmax - mesh->zmin);
              break;
          }

          cam->altinc *= 0.005f;
          sc->type |= S_RESET;
          reshapeScene(sc->par.xs, sc->par.ys);
        }
        /* tilt cut plane */
        else if (clip->active & C_EDIT) {
          tiltClip(sc, clip);
        } else {
          return;
        }

        post = TRUE;
        break;

      case 'z': /* zoom in */
        if (cube->active & C_EDIT) {
          cube->cma[2] += 0.1 * sc->dmax;
          post = TRUE;
          break;
        } else {
          if (p->rubber == 2)
            setPersp(sc, p, 1);
          else
            p->fovy = max(0.9 * p->fovy, 1e-05);

          farclip(1);
        }

        post = TRUE;

        /* get linked view */
        sc1 = sc;

        while (sc1->slave > -1) {
          sc1 = cv.scene[sc1->slave];
          memcpy(sc1->view, sc->view, sizeof(struct transform));
          memcpy(sc1->persp, sc->persp, sizeof(struct sperspective));
          glutSetWindow(sc1->idwin);
          reshapeScene(sc1->par.xs, sc1->par.ys);
        }

        glutSetWindow(sc->idwin);
        break;
    }
  } else if (isupper(key)) {
    switch (key) {
      case 'A': /* longer arrows (as in ffglut; the axis are in the menu) */
        vecLength(sc, mesh, 1.2);
        post = TRUE;
        break;
      case 'B': /* toggle box */
        keyItem('B', 0, 0);
        break;
      case 'C': /* toggle capping */
        keyClip('K', 0, 0);
        break;
      case 'D': /* selection cube */
        if (mesh->dim == 3) keyItem('D', 0, 0);

        break;
      case 'E':
        keyColor('E', 0, 0);
        glutSetWindow(sc->idwin);
        break;
      case 'F': /* toggle face nums */
        keyItem('F', 0, 0);
        break;
      case 'G': /* toggle grid */
        keyItem('G', 0, 0);
        break;
      case 'H': /* hardcopy PPM */
        keyFile('H', 0, 0);
        break;
      case 'I': /* follows mouse */
        sc->type ^= S_FOLLOW;
        sc1 = sc;

        while (sc1->slave > -1) {
          sc1 = cv.scene[sc1->slave];
          sc1->type = sc->type;
        }

        break;
      case 'J':
        if (mesh->dim == 2 && !(sc->mode & S_ALTITUDE)) break;

        if (p->pmode == PERSPECTIVE) {
          p->pmode = CAMERA;
          glutMouseFunc(mouseCamera);
          glutMotionFunc(motionCamera);
        } else {
          p->pmode = PERSPECTIVE;
          glutMouseFunc(mouse);
          glutMotionFunc(motion);
        }

        sc->type = S_RESET;
        reshapeScene(sc->par.xs, sc->par.ys);
        post = TRUE;
        break;
      case 'K':
        keyMetric('K', 0, 0);
        break;
      case 'L':
        keyFile('L', 0, 0);
        break;
      case 'M':
        morphing = 1 - morphing;
        post = TRUE;
        break;
      case 'N': /* draw normals */
      case 'O':
        if (mesh->nvn) keyItem(key, 0, 0);

        break;
      case 'P': /* toggle point nums */
        keyItem('P', 0, 0);
        break;
      case 'Q': /* palette (p is the previous plot with FreeFEM) */
        keyMetric('p', 0, 0);
        break;
      case 'R':
        if (ipilmat < 1) return;

        refmat = pilmat[ipilmat--];
        pm = &sc->material[refmat];
        pm->flag = 0;
        updatePoints(sc, mesh, refmat);
        dolist = TRUE;
        post = TRUE;
        break;
      case 'S': /* softcopy */
        keyAnim('S', 0, 0);
        break;
      case 'T': /* graduated axes x, y (and z) with the values */
        sc->axes = !sc->axes;
        post = TRUE;
        break;
      case 'V': /* change center of scene */
        if (!refitem) return;

        changeCenter(sc, mesh);
        reshapeScene(sc->par.xs, sc->par.ys);
        post = TRUE;
        /* linked view */
        sc1 = sc;

        while (sc1->slave > -1) {
          sc1 = cv.scene[sc1->slave];
          sc1->cx = sc->cx;
          sc1->cy = sc->cy;
          sc1->cz = sc->cz;
          memcpy(sc1->view, sc->view, sizeof(struct transform));
          memcpy(sc1->persp, sc->persp, sizeof(struct sperspective));
          glutSetWindow(sc1->idwin);
          reshapeScene(sc1->par.xs, sc1->par.ys);
        }

        glutSetWindow(sc->idwin);
        break;
      case 'W':
        keyFile('W', 0, 0);
        break;
      case 'X': /* close window */
        if (cube->active & C_EDIT) {
          cube->cma[0] -= 0.1 * sc->dmax;
          post = TRUE;
          break;
        }

        if (idw != cv.nbs) {
          deleteScene(sc);

          for (k = idw + 1; k < cv.nbs; k++) {
            cv.scene[k - 1] = cv.scene[k];
          }

          cv.scene[cv.nbs - 1] = 0;
        }

        glutHideWindow( );
        if (--cv.nbs == 0) exit(0);

        break;

      case 'Y': /* close window */
        if (cube->active & C_EDIT) {
          cube->cma[1] -= 0.1 * sc->dmax;
          post = TRUE;
          break;
        }

      case 'Z': /* zoom out */
        if (cube->active & C_EDIT) {
          cube->cma[2] -= 0.1 * sc->dmax;
          post = TRUE;
          break;
        } else {
          if (p->rubber == 2)
            setPersp(sc, p, 0);
          else
            p->fovy = min(1.1 * p->fovy, 179.0);

          farclip(1);
        }

        post = TRUE;

        /* get linked view */
        sc1 = sc;

        while (sc1->slave > -1) {
          sc1 = cv.scene[sc1->slave];
          memcpy(sc1->view, sc->view, sizeof(struct transform));
          memcpy(sc1->persp, sc->persp, sizeof(struct sperspective));
          glutSetWindow(sc1->idwin);
          reshapeScene(sc1->par.xs, sc1->par.ys);
        }

        glutSetWindow(sc->idwin);
        break;
    }
  } else {
    switch (key) {
      case '.': /* value of the solution under the mouse */
      {
        extern int probeon;

        probeon = !probeon;
        post = TRUE;
        break;
      }
      case '-':
        if (keyact & GLUT_ACTIVE_ALT) {
          keyAnim('p', 0, 0);
          break;
        }

        if (cube->active & C_EDIT) {
          cube->cma[0] *= 0.95;
          cube->cma[1] *= 0.95;
          cube->cma[2] *= 0.95;
          post = TRUE;
          break;
        }

        if (sc->persp->pmode == CAMERA) {
          cam = sc->camera;
          cam->spmod -= 1.e-04 * sc->dmax;
          post = TRUE;
          break;
        }

        if (p->depth < -20.0 * sc->dmax) break;

        p->depth -= 0.1 * sc->dmax;
        farclip(1);
        post = TRUE;

        /* get linked view */
        sc1 = sc;

        while (sc1->slave > -1) {
          sc1 = cv.scene[sc1->slave];
          memcpy(sc1->view, sc->view, sizeof(struct transform));
          memcpy(sc1->persp, sc->persp, sizeof(struct sperspective));
          glutSetWindow(sc1->idwin);
          reshapeScene(sc1->par.xs, sc1->par.ys);
        }

        glutSetWindow(sc->idwin);
        break;

      case '+':
        if (keyact & GLUT_ACTIVE_ALT) {
          keyAnim('n', 0, 0);
          break;
        }

        if (cube->active & C_EDIT) {
          cube->cma[0] *= 1.05;
          cube->cma[1] *= 1.05;
          cube->cma[2] *= 1.05;
          post = TRUE;
          break;
        }

        if (sc->persp->pmode == CAMERA) {
          cam = sc->camera;
          cam->spmod += 1.e-04 * sc->dmax;
          post = TRUE;
          break;
        }

        if (p->depth > 0.0) break;

        p->depth += 0.1 * sc->dmax;
        farclip(1);
        post = TRUE;

        /* update linked view */
        sc1 = sc;

        while (sc1->slave > -1) {
          sc1 = cv.scene[sc1->slave];
          memcpy(sc1->view, sc->view, sizeof(struct transform));
          memcpy(sc1->camera, sc->camera, sizeof(struct camera));
          memcpy(sc1->persp, sc->persp, sizeof(struct sperspective));
          glutSetWindow(sc1->idwin);
          reshapeScene(sc1->par.xs, sc1->par.ys);
        }

        glutSetWindow(sc->idwin);
        break;

      case '#': /* select entity */
        fprintf(stdout, "ENTITY NUMBER: ");
        fflush(stdout);
        while (fgetc(stdin) != EOF)
          ; /* fflush() called on input stream 'stdin' may result in undefined behaviour on
               non-linux systems */
        ret = fscanf(stdin, "%d", &numit);
        if (ret == EOF) printf("fscanf error\n");
        if (sc->picklist) glDeleteLists(sc->picklist, 1);

        if (numit > 0) sc->picklist = pickItem(mesh, sc, numit);

        post = TRUE;
        break;
      case '=':
        picking = GL_TRUE;
        if (sc->picklist) glDeleteLists(sc->picklist, 1);

        sc->picklist = pickingScene(sc, x, y, 0);
        break;

      case '!': /* clip plane */
        return;
        // if (!(clip->active & C_ON)) return;
        //
        // dd = clip->eqn[3] - clip->eqn[0] * mesh->xtra - clip->eqn[1] * mesh->ytra -
        //      clip->eqn[2] * mesh->ztra;
        // fprintf(stdout, "\nCurrent plane: %gx %+gy %+gz %+g = 0\n", clip->eqn[0], clip->eqn[1],
        //         clip->eqn[2], dd);
        // fprintf(stdout, "Plane coeffs : ");
        // fflush(stdout);
        // while (fgetc(stdin) != EOF)
        //   ; /* fflush() called on input stream 'stdin' may result in undefined behaviour on
        //        non-linux systems */
        // ret = fscanf(stdin, "%f %f %f %f", &a, &b, &c, &d);
        // if (ret == EOF) printf("fscanf error\n");
        // resetClip(sc, clip, mesh);
        // clip->eqn[0] = a;
        // clip->eqn[1] = b;
        // clip->eqn[2] = c;
        // clip->eqn[3] = d;
        //
        // fprintf(stdout, "New plane eq.: ");
        // if (clip->eqn[0]) fprintf(stdout, "%+gx", clip->eqn[0]);
        //
        // if (clip->eqn[1]) fprintf(stdout, " %+gy", clip->eqn[1]);
        //
        // if (clip->eqn[2]) fprintf(stdout, " %+gz", clip->eqn[2]);
        //
        // if (clip->eqn[3]) fprintf(stdout, " %+g", clip->eqn[3]);
        //
        // fprintf(stdout, " = 0\n");
        // clip->eqn[3] += (a * mesh->xtra + b * mesh->ytra + c * mesh->ztra);
        // post = TRUE;
        // dolist = TRUE;
        // break;

      case '@': /* add trajectoire point */
        if (p->pmode == CAMERA) pathAdd(sc, x, y);

        break;

      case '%':
        fprintf(stdout, "reference (%d): ", refpick);
        fflush(stdout);
        while (fgetc(stdin) != EOF)
          ; /* fflush() called on input stream 'stdin' may result in undefined behaviour on
               non-linux systems */
        ret = fscanf(stdin, "%d", &refpick);
        if (ret == EOF) printf("fscanf error\n");
        break;

      case '&':
        puts("ADJUST");
        parEdit(sc);
        break;
    }
  }

  if (dolist) {
    doLists(sc, mesh);
    doMapLists(sc, mesh, 1);
    doIsoLists(sc, mesh, 1);
  }

  if (post) glutPostRedisplay( );
}

#ifdef __cplusplus
}
#endif
