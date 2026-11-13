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

GLfloat altcoef = 0.0;
static int nbimg;
extern int refitem, refmat, reftype, imstep, imreverse;
extern int *pilmat, ipilmat;

#define MAX_LST 128

/* rebuild display lists */
void doLists(pScene sc, pMesh mesh) {
  int k;

  /*default */
  if (ddebug) printf("build display lists\n");

  if (!morphing) glutSetCursor(GLUT_CURSOR_WAIT);

  for (k = 0; k < MAX_LIST; k++) {
    if (sc->dlist[k]) glDeleteLists(sc->dlist[k], 1);

    if (sc->clist[k]) glDeleteLists(sc->clist[k], 1);

    sc->dlist[k] = sc->clist[k] = (GLuint)0;
  }

  /* poly lists */
  sc->dlist[LTria] = listTria(sc, mesh);
  sc->dlist[LQuad] = listQuad(sc, mesh);
  if (!mesh->nt && !mesh->nq) {
    sc->dlist[LTets] = listTetra(sc, mesh, 0);
    sc->dlist[LHexa] = listHexa(sc, mesh, 0);
  }

  if (sc->clip->active & C_VOL) sc->clip->active |= C_REDO;

  glutSetCursor(GLUT_CURSOR_INHERIT);
  checkErrors( );
}

/* build metric list */
void doMapLists(pScene sc, pMesh mesh, int reset) {
  /*default */
  if (!mesh->sol) return;

  if (ddebug) printf("build map lists\n");

  glutSetCursor(GLUT_CURSOR_WAIT);

  /* delete old lists */
  if (reset) {
    int k;

    for (k = 0; k < MAX_LIST; k++) {
      if (sc->mlist[k]) glDeleteLists(sc->mlist[k], 1);

      sc->mlist[k] = (GLuint)0;
    }
  }

  /* create map list */
  if (sc->mode & S_ALTITUDE) {
    sc->mlist[LTria] = alt2dList(sc, mesh, LTria, sc->shrink, altcoef);
    sc->mlist[LQuad] = alt2dList(sc, mesh, LQuad, sc->shrink, altcoef);
    /* the mesh on the elevation map (the list of the tets is free for a 2d mesh) */
    if (mesh->dim == 2) sc->mlist[LTets] = alt2dMesh(sc, mesh, sc->shrink, altcoef);
  } else if (sc->mode & S_MAP) {
    if (mesh->nt && !sc->mlist[LTria]) sc->mlist[LTria] = listTriaMap(sc, mesh);

    if (mesh->nq && !sc->mlist[LQuad]) sc->mlist[LQuad] = listQuadMap(sc, mesh);

    if (!mesh->nt && mesh->ntet && !sc->mlist[LTets]) sc->mlist[LTets] = listTetraMap(sc, mesh, 0);

    if (!mesh->nq && mesh->nhex && !sc->mlist[LHexa]) sc->mlist[LHexa] = listHexaMap(sc, mesh, 0);
  }

  glutSetCursor(GLUT_CURSOR_INHERIT);
  checkErrors( );
}

/* rebuild iso-values lists */
void doIsoLists(pScene sc, pMesh mesh, int reset) {
  /*default */
  if (ddebug) printf("build iso lists\n");

  glutSetCursor(GLUT_CURSOR_WAIT);

  /* delete old lists */
  if (reset) {
    int k;

    for (k = 0; k < MAX_LIST; k++) {
      if (sc->ilist[k]) glDeleteLists(sc->ilist[k], 1);

      if (sc->vlist[k]) glDeleteLists(sc->vlist[k], 1);

      if (sc->cplist) glDeleteLists(sc->cplist, 1);

      sc->ilist[k] = sc->vlist[k] = (GLuint)0;
      sc->cplist = (GLuint)0;
    }

    if (sc->stream && sc->stream->nbstl) {
      int kk;

      for (kk = 0; kk < sc->stream->nbstl; kk++) {
        if (sc->slist[kk]) glDeleteLists(sc->slist[kk], 1);

        sc->slist[kk] = (GLuint)0;
      }

      sc->slist = NULL;
      if (reset < 2) sc->stream->nbstl = 0;

      for (k = 1; k <= mesh->np; k++) {
        pPoint ppt;

        ppt = &mesh->point[k];
        ppt->flag = 0;
      }
    }
  }

  /* iso-lines */
  if (sc->isotyp & S_ISOLINE) {
    if (mesh->nt && !sc->ilist[LTria]) sc->ilist[LTria] = listTriaIso(sc, mesh);

    if (mesh->nq && !sc->ilist[LQuad]) sc->ilist[LQuad] = listQuadIso(sc, mesh);
  }

  /* iso-surfaces */
  if (sc->isotyp & S_ISOSURF)
    if (mesh->ntet && !sc->ilist[LTets]) sc->ilist[LTets] = listTetraIso(sc, mesh);

  /* vector */
  if (sc->isotyp & S_VECTOR) {
    if (mesh->dim == 2) {
      if (mesh->nt && !sc->vlist[LTria]) sc->vlist[LTria] = listTria2dVector(mesh);

      if (mesh->nq && !sc->vlist[LQuad]) sc->vlist[LQuad] = listQuad2dVector(mesh);
    } else {
      if (mesh->ntet + mesh->nhex == 0) {
        sc->vlist[LTria] = listTria3dVector(mesh);
      } else {
        /* without clipping plane: all the vectors of the mesh */
        if (!sc->vlist[LTria]) sc->vlist[LTria] = listMeshVector3d(mesh);

        if (sc->clip->active & C_ON) {
          if (mesh->ntet && !sc->vlist[LTets])
            sc->vlist[LTets] = listClipTetraVector(mesh);
          else if (mesh->nhex && !sc->vlist[LHexa])
            sc->vlist[LHexa] = listClipHexaVector(mesh);
        }
      }
    }
  }

  /* streamlines */
  if (sc->isotyp & S_STREAML) {
    if (!mesh->adja || !sc->slist) {
      sc->stream = createStream(sc, mesh);
      if (!sc->stream) sc->isotyp &= ~S_STREAML;
    }

    if (sc->stream) {
      if (reftype == LPoint) {
        streamRefPoint(sc, mesh);
      } else if (mesh->dim == 3) {
        if (reftype == LTria)
          streamRefTria(sc, mesh);
        else if (reftype == LQuad)
          streamRefQuad(sc, mesh);

        if (sc->picklist) glDeleteLists(sc->picklist, 1);

        sc->picklist = 0;
      }
    }
  }

  /* critical points */
  if (sc->isotyp & S_CRITP)
    if (mesh->dim == 2)
      if (mesh->nt && !sc->cplist) sc->cplist = listCritPoint(sc, mesh);

  glutSetCursor(GLUT_CURSOR_INHERIT);
  checkErrors( );
}

/* overlay (key U): the data of the other scenes on the same mesh are plotted in this scene:
   iso-surfaces (or iso-lines in 2d) of a scalar solution, arrows of a vector field.
   Ex: medit("u v",Th,u,[vx,vy,vz]) gives 2 windows, U in one of them shows both. */
void doOverlay(pScene sc) {
  pMesh mesh = cv.mesh[sc->idmesh];
  int j, n = 0;

  for (j = 0; j < MAX_OVER; j++) {
    if (sc->olist[j]) glDeleteLists(sc->olist[j], 1);

    sc->olist[j] = (GLuint)0;
  }

  if (!sc->overlay) return;

  for (j = 0; j < cv.nbs && n < MAX_OVER; j++) {
    pScene scj = cv.scene[j];
    pMesh mj;

    if (!scj || scj == sc) continue;

    mj = cv.mesh[scj->idmesh];
    if (!mj->nbb || mj->np != mesh->np || mj->dim != mesh->dim) continue;

    if (mj->nfield == 1) {
      if (mj->ntet) {
        sc->olist[n] = listTetraIso(scj, mj);
        sc->otyp[n] = 1;
      } else if (mj->dim == 2 && mj->nt) {
        sc->olist[n] = listTriaIso(scj, mj);
        sc->otyp[n] = 0;
      }
    } else if (mj->nfield == mj->dim) {
      /* the lists of vectors use the palette of the current scene */
      Iso iso = sc->iso;

      sc->iso = scj->iso;
      if (mj->dim == 2)
        sc->olist[n] = listTria2dVector(mj);
      else if (mj->ntet + mj->nhex)
        sc->olist[n] = listMeshVector3d(mj);

      sc->iso = iso;
      sc->otyp[n] = 0;
    }

    if (sc->olist[n]) n++;
  }

  if (!n && !graph1dMesh(mesh)) fprintf(stdout, "  ## Overlay: no other solution on the same mesh.\n");
}

void resetLists(pScene sc, pMesh mesh) {
  int kk;

  for (kk = 0; kk < MAX_OVER; kk++) {
    if (sc->olist[kk]) glDeleteLists(sc->olist[kk], 1);

    sc->olist[kk] = (GLuint)0;
  }

  for (kk = 0; kk < MAX_LIST; kk++) {
    if (sc->dlist[kk]) glDeleteLists(sc->dlist[kk], 1);

    if (sc->mlist[kk]) glDeleteLists(sc->mlist[kk], 1);

    if (sc->ilist[kk]) glDeleteLists(sc->ilist[kk], 1);

    if (sc->clist[kk]) glDeleteLists(sc->clist[kk], 1);

    if (sc->cmlist[kk]) glDeleteLists(sc->cmlist[kk], 1);

    if (sc->vlist[kk]) glDeleteLists(sc->vlist[kk], 1);

    if (sc->cplist) glDeleteLists(sc->cplist, 1);

    sc->dlist[kk] = sc->clist[kk] = sc->mlist[kk] = (GLuint)0;
    sc->ilist[kk] = sc->vlist[kk] = sc->cplist = (GLuint)0;
  }

  if (sc->glist) glDeleteLists(sc->glist, 1);

  if (sc->nlist) glDeleteLists(sc->nlist, 1);

  sc->glist = sc->nlist = (GLuint)0;
}

void keyFile(unsigned char key, int x, int y) {
  pScene sc;
  pMesh mesh;
  pTransform view;
  char *ptr, data[128], tmpdata[128];
  ;
  ubyte post = FALSE, clipon;

  /* default */
  sc = cv.scene[currentScene( )];
  mesh = cv.mesh[sc->idmesh];
  view = sc->view;

  switch (key) {
    case 'L': /* load prefs */
      if (parsop(sc, mesh)) {
        refitem = reftype = 0;
        setupPalette(sc, mesh);
        initGrafix(sc, mesh);
        doLists(sc, mesh);
        if (sc->mode & S_MAP) doMapLists(sc, mesh, 1);

        if (sc->isotyp) doIsoLists(sc, mesh, 1);

        glClearColor(sc->par.back[0], sc->par.back[1], sc->par.back[2], sc->par.back[3]);
        post = GL_TRUE;
      }

      break;
    case 'W': /* save prefs */
      saveMeditFile(mesh->name, sc);
      break;

    case 'R': /* reload mesh */
      refitem = reftype = 0;
      if (!meshUpdate(sc, mesh)) exit(1);

      meshRef(sc, mesh);
      matSort(sc);
      setupPalette(sc, mesh);
      doLists(sc, mesh);
      if (sc->mode & S_MAP) doMapLists(sc, mesh, 1);

      if (sc->isotyp) doIsoLists(sc, mesh, 1);

      glClearColor(sc->par.back[0], sc->par.back[1], sc->par.back[2], sc->par.back[3]);
      post = GL_TRUE;
      break;

    case 'S': /* save mesh */
      strcpy(tmpdata, mesh->name);
      ptr = (char *)strstr(tmpdata, ".mesh");
      if (ptr) *ptr = '\0';

      sprintf(data, "%s.d.mesh", tmpdata);
      clipon = sc->clip->active & C_ON;
      if (clipon) clipVertices(mesh, sc, sc->clip);

      saveMesh(sc, mesh, data, clipon);
      break;

    case 'B':
    case 'C':
    case 'G':
    case 'H': /* hardcopy */
    case 'T':
      strcpy(tmpdata, mesh->name);
      ptr = (char *)strstr(tmpdata, ".mesh");
      if (ptr) *ptr = '\0';

      ptr = (char *)strstr(tmpdata, ".gis");
      if (ptr) *ptr = '\0';

      if (option != SEQUENCE) {
        static int nfree = 0;

        if (key == 'H')
          nfree = filnum(tmpdata, nfree, "ppm");
        else
          nfree = filnum(tmpdata, nfree, "ps");

        if (nfree == -1) break;

        sprintf(data, "%s.%.3d", tmpdata, nfree);
      } else {
        sprintf(data, "%s.ppm", tmpdata);
      }

      if (!saveimg) glutSetCursor(GLUT_CURSOR_WAIT);

      if (key == 'B' || key == 'C' || key == 'G')
        imgTiling(sc, data, key);
      else
        imgHard(sc, data, key);

      if (animate && !(sc->persp->pmode & CAMERA) && ++nbimg > 179) {
        view->angle = 0.0;
        nbimg = 0;
        keyAnim('A', 0, 0);
      }

      if (!saveimg) glutSetCursor(GLUT_CURSOR_INHERIT);

      break;
    case 's': /* softcopy */
      sftcpy(sc, mesh);
      post = FALSE;
      break;
  }

  if (post) glutPostRedisplay( );
}

/* The actions of the menus are deferred to the next display of the window: in a menu
   callback, the OpenGL context of the window is not always the current one (GLUT on macOS),
   so the display lists built by the action (vectors, iso-surfaces, ...) were lost. */
typedef void (*MenuAction)(unsigned char, int, int);
static struct {
  int win;
  MenuAction f;
  unsigned char key;
} menupend[16];
static int nmenupend = 0;

static void menuDefer(MenuAction f, int item) {
  extern int mouseWindow(void);
  int w = glutGetWindow( ), wm = mouseWindow( ), k;

  /* the window of the menu = the window under the mouse (glutGetWindow can be another
     window when several windows are open) */
  for (k = 0; k < MAX_SCENE; k++)
    if (cv.scene[k] && cv.scene[k]->idwin == wm) break;

  if (wm > 0 && k < MAX_SCENE && wm != w) {
    fprintf(stdout, "  menu: window %d (GLUT current window %d)\n", wm, w);
    w = wm;
    glutSetWindow(w);
  }

  if (nmenupend < 16) {
    menupend[nmenupend].win = w;
    menupend[nmenupend].f = f;
    menupend[nmenupend].key = (unsigned char)item;
    nmenupend++;
  }

  glutPostRedisplay( );
}

/* run the deferred actions of the current window (at the beginning of its display): 1 if any */
int menuRun(void) {
  int w = glutGetWindow( ), i, j, ran = 0;

  for (i = 0; i < nmenupend;) {
    if (menupend[i].win == w) {
      MenuAction f = menupend[i].f;
      unsigned char key = menupend[i].key;

      for (j = i + 1; j < nmenupend; j++) menupend[j - 1] = menupend[j];

      nmenupend--;
      f(key, 0, 0);
      ran = 1;
      if (glutGetWindow( ) != w) glutSetWindow(w);
    } else {
      i++;
    }
  }

  return (ran);
}

/* keyScene from a menu: no modifier (glutGetModifiers is only valid in the input callbacks) */
static void keySceneMenu(unsigned char key, int x, int y) {
  extern int menukey;

  menukey = 1;
  keyScene(key, x, y);
  menukey = 0;
}

void menuFile(int item) { menuDefer(keyFile, item); }

void keyItem(unsigned char key, int x, int y) {
  pScene sc;
  pMesh mesh;
  pCube cube;
  ubyte post = TRUE;

  /* default */
  sc = cv.scene[currentScene( )];
  mesh = cv.mesh[sc->idmesh];
  cube = sc->cube;

  switch (key) {
    case 'A': /* draw axis */
      sc->item ^= S_AXIS;
      break;
    case 'B': /* draw bounding box */
      sc->item ^= S_BOX;
      break;
    case 'D': /*draw selection cube */
      if (cube->active & C_ON) {
        cube->active &= ~(C_ON + C_EDIT);
        dumpCube(sc, mesh, cube);
      } else {
        cube->active |= C_ON;
      }

      break;
    case 'G': /* draw grid */
      sc->item ^= S_GRID;
      break;
    case 'j': /* toggle info */
      sc->type ^= S_DECO;
      break;
    case 'P': /* toggle point nums */
      sc->item ^= S_NUMP;
      break;
    case 'F': /* toggle face nums */
      sc->item ^= S_NUMF;
      break;
    case 'T': /* graduated axes (key j) */
      sc->axes = !sc->axes;
      break;
    case 'g': /* const. items */
      if (!sc->glist) post = FALSE;

      sc->item ^= S_GEOM;
      break;
    case 'N':
      if (mesh->nvn)
        sc->type ^= S_NORMAL;
      else
        post = FALSE;

      break;
    case 'O':
      if (mesh->nvn) {
        sc->type ^= S_OPPOS;
        if (sc->nlist) {
          glDeleteLists(sc->nlist, 1);
          sc->nlist = 0;
        } else {
          post = FALSE;
        }
      } else {
        post = FALSE;
      }

      break;
  }

  if (post) glutPostRedisplay( );
}

void menuItem(int item) { menuDefer(keyItem, item); }

void keyAnim(unsigned char key, int x, int y) {
  pScene sc;
  pMesh mesh;
  pClip clip;
  ubyte post = TRUE;
  char *ptr, base[256];
  static int depart = -1;

  /* default */
  sc = cv.scene[currentScene( )];
  mesh = cv.mesh[sc->idmesh];
  clip = sc->clip;
  if (depart == -1) depart = animdep;

  switch (key) {
    case 'A': /* switch animation */
      if (animate) {
        sc->view->manim = GL_FALSE;
        clip->cliptr->manim = GL_FALSE;
        sc->view->angle = 0.0f;
        glutIdleFunc(NULL);
      } else {
        sc->view->manim = GL_TRUE;
        clip->cliptr->manim = GL_TRUE;
        if (sc->persp->pmode == CAMERA) glutIdleFunc(glutIdle);
      }

      animate = 1 - animate;
      post = FALSE;
      break;
    case 'I': /* save image */
      saveimg = 1 - saveimg;
      nbimg = 0;
      post = GL_FALSE;
      break;

    case 'M': /* morphing */
      if (morphing)
        glutIdleFunc(NULL);
      else
        glutIdleFunc(glutIdle);

      morphing = 1 - morphing;
      post = FALSE;
      break;
    case 'R': /* next morph */
      imreverse = 1 - imreverse;
      break;

    case 'S': /* start animation */
      if (ddebug) fprintf(stdout, "debut sequence %d a %d\n", animdep, animfin);

      glutSetWindow(sc->idwin);
      if (option == SEQUENCE) playAnim(sc, mesh, animdep, animfin);

      if (!saveimg) post = FALSE;

      break;

    case 'f': /* first mesh */
      /* get basename */
      ptr = (char *)strrchr(mesh->name, '.');
      if (ptr) *ptr = '\0';

      strcpy(base, mesh->name);
      resetLists(sc, mesh);
      if (!loadNextMesh(mesh, animdep, 0)) break;

      doLists(sc, mesh);
      if (sc->mode & S_MAP) doMapLists(sc, mesh, 0);

      if (sc->isotyp) doIsoLists(sc, mesh, 0);

      strcpy(mesh->name, base);
      break;

    case 'l':
      /* get basename */
      ptr = (char *)strrchr(mesh->name, '.');
      if (ptr) *ptr = '\0';

      strcpy(base, mesh->name);
      resetLists(sc, mesh);
      if (!loadNextMesh(mesh, animfin, 0)) break;

      doLists(sc, mesh);
      if (sc->mode & S_MAP) doMapLists(sc, mesh, 0);

      if (sc->isotyp) doIsoLists(sc, mesh, 0);

      strcpy(mesh->name, base);
      break;

    case 'n':
      /* get basename */
      if (++depart > animfin) depart = animdep;

      ptr = (char *)strrchr(mesh->name, '.');
      if (ptr) *ptr = '\0';

      strcpy(base, mesh->name);
      resetLists(sc, mesh);
      if (!loadNextMesh(mesh, depart, 0)) break;

      doLists(sc, mesh);
      if (sc->mode & S_MAP) doMapLists(sc, mesh, 0);

      if (sc->isotyp) doIsoLists(sc, mesh, 0);

      strcpy(mesh->name, base);
      break;

    case 'p':
      /* get basename */
      if (--depart < animdep) depart = animfin;

      ptr = (char *)strrchr(mesh->name, '.');
      if (ptr) *ptr = '\0';

      strcpy(base, mesh->name);
      resetLists(sc, mesh);
      if (!loadNextMesh(mesh, depart, 0)) break;

      doLists(sc, mesh);
      if (sc->mode & S_MAP) doMapLists(sc, mesh, 0);

      if (sc->isotyp) doIsoLists(sc, mesh, 0);

      strcpy(mesh->name, base);
      break;
  }

  if (post == TRUE) glutPostRedisplay( );
}

void menuAnim(int item) { menuDefer(keyAnim, item); }

void keyTrajet(unsigned char key, int x, int y) {
  pScene sc;
  pMesh mesh;

  sc = cv.scene[currentScene( )];
  mesh = cv.mesh[sc->idmesh];

  switch (key) {
    case 'C': /* add control point */
      pathAdd(sc, x, y);
      break;
    case 'S': /* show all points */
      sc->type ^= S_PATH;
      if (sc->type & S_PATH && !sc->path.tlist) sc->path.tlist = pathList(sc);

      break;
    case 'F': /* follow path */
      break;
    case 'L': /* load path */
      pathLoad(mesh->name, sc);
      break;
    case 'W': /* save path */
      pathSave(mesh->name, sc);
      break;
  }

  glutPostRedisplay( );
}

void menuTrajet(int item) { menuDefer(keyTrajet, item); }

void keyMode(unsigned char key, int x, int y) {
  pScene sc;
  pMesh mesh;
  ubyte post = TRUE, dolist = FALSE, material;

  sc = cv.scene[currentScene( )];
  mesh = cv.mesh[sc->idmesh];
  material = sc->mode & S_MATERIAL;

  switch (key) {
    case 'D': /* depth lines */
      if (sc->mode == DEPTH) break;

      sc->mode = DEPTH;
      break;
    case 'S': /* toggle smooth shading */
      if (sc->mode == FILL) break;

      sc->mode = FILL;
      if (material) sc->mode |= S_MATERIAL;

      break;
    case 'P': /* toggle smooth shaded polys */
      if (sc->mode == SHADED) break;

      sc->mode = SHADED;
      if (material) sc->mode |= S_MATERIAL;

      break;
    case 'H': /* hidden lines */
      if (sc->mode == HIDDEN) break;

      sc->mode = HIDDEN;
      if (material) sc->mode |= S_MATERIAL;

      break;
    case 'W': /* wireframe */
      if (sc->mode == WIRE) break;

      sc->mode = WIRE;
      if (material) sc->mode |= S_MATERIAL;

      break;
    case 'n': /* toggle normals */
      if (mesh->nvn == 0) {
        post = FALSE;
        break;
      }

      sc->type ^= S_FLAT;
      dolist = TRUE;
      break;
    default:
      post = FALSE;
      break;
  }

  if (dolist == TRUE) doLists(sc, mesh);

  if (post == TRUE) glutPostRedisplay( );
}

void menuMode(int item) { menuDefer(keyMode, item); }

void menuScene(int item) {
  if (item == 'X') menuOld(item); /* X is now a key of ffglut: the old action of medit */
  else menuKey(item);
}

void keyView(unsigned char key, int x, int y) {
  pScene sc, sc1;
  pMesh mesh;
  float dmax;
  ubyte post = FALSE;

  sc = cv.scene[currentScene( )];
  mesh = cv.mesh[sc->idmesh];

  switch (key) {
    case 'R':
      sc->type |= S_RESET;
      dmax = mesh->xmax - mesh->xmin;
      dmax = max(dmax, mesh->ymax - mesh->ymin);
      dmax = max(dmax, mesh->zmax - mesh->zmin);
      sc->cx = sc->cy = sc->cz = 0.0f;
      sc->dmax = fabs(dmax);
      if (sc->persp->pmode == PERSPECTIVE) {
        resetTransform(sc->view);
        initPersp(sc->persp, sc->dmax);
      } else if (sc->persp->pmode == CAMERA) {
        initPersp(sc->persp, sc->dmax);
        initCamera(sc, sc->camera->vecup);
        sc->persp->pmode = CAMERA;
      }

      reshapeScene(sc->par.xs, sc->par.ys);
      post = TRUE;
      break;
    case 'C':
      copyView(sc->view, sc->camera, sc->persp);
      copyClip(sc->clip);
      break;
    case 'L': /* link view */
      if (!linkView(sc)) return;

      reshapeScene(sc->par.xs, sc->par.ys);
      post = TRUE;
      break;
    case 'P':
      if (pasteView(sc->view, sc->camera, sc->persp) && pasteClip(sc->clip)) {
        reshapeScene(sc->par.xs, sc->par.ys);
        post = TRUE;
      }

      break;
    case 'U':
      unlinkView(sc);
      break;

    case 'D': /* duplicate view */
      if (cv.nbs == MAX_SCENE) break;

      if (!cv.scene[++cv.nbs]) {
        cv.scene[cv.nbs] = (pScene)M_calloc(1, sizeof(Scene), "menus.scene");
        if (!cv.scene[cv.nbs]) break;

        sc1 = cv.scene[cv.nbs];
        sc1->material = (pMaterial)calloc(2 + sc->par.nbmat, sizeof(Material));
        assert(sc1->material);
      }

      sc1 = cv.scene[cv.nbs];
      memcpy(sc1, sc, sizeof(Scene));
      memcpy(sc1->material, sc->material, sc->par.nbmat * sizeof(Material));
      memcpy(&sc1->par, &sc->par, sizeof(Param));
      if (!createScene(sc1, sc1->idmesh)) {
        fprintf(stdout, "  ## Unable to create\n");
        return;
      }

      copyView(sc->view, sc->camera, sc->persp);
      copyClip(sc->clip);
      pasteView(sc1->view, sc1->camera, sc1->persp);
      pasteClip(sc1->clip);
      break;
  }

  if (post) glutPostRedisplay( );
}

void menuView(int item) { menuDefer(keyView, item); }

void keyColor(unsigned char key, int x, int y) {
  pScene sc;
  pMesh mesh;
  pMaterial pm;
  int i, k;
  ubyte post = TRUE, dolist = FALSE;

  sc = cv.scene[currentScene( )];
  mesh = cv.mesh[sc->idmesh];

  switch (key) {
    case 'b': /* reverse backcolor */
      sc->par.back[0] = 1.0f - sc->par.back[0];
      sc->par.back[1] = 1.0f - sc->par.back[1];
      sc->par.back[2] = 1.0f - sc->par.back[2];
      glClearColor(sc->par.back[0], sc->par.back[1], sc->par.back[2], sc->par.back[3]);
      if (!sc->par.linc) {
        sc->par.line[0] = 1.0f - sc->par.line[0];
        sc->par.line[1] = 1.0f - sc->par.line[1];
        sc->par.line[2] = 1.0f - sc->par.line[2];
      }

      break;
    case 'e': /* toggle matcolors */
      sc->mode ^= S_MATERIAL;
      dolist = TRUE;
      break;
    case 'E': /* edit matcolors */
      post = FALSE;
      matEdit(sc);
      break;
    case 'r':
      if (refmat < 0 || ipilmat == sc->par.nbmat) break;

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
    case 'R': /* reset materials */

      for (k = 0; k < sc->par.nbmat; k++) {
        pm = &sc->material[k];

        for (i = LTria; i <= LHexa; i++) {
          pm->depmat[i] = abs(pm->depmat[i]);
        }
      }

      dolist = TRUE;
      break;
    default:
      post = FALSE;
  }

  if (dolist) doLists(sc, mesh);

  if (post) glutPostRedisplay( );
}

void menuColor(int item) { menuDefer(keyColor, item); }

void keyClip(unsigned char key, int x, int y) {
  pScene sc;
  pMesh mesh;
  pClip clip;
  ubyte post = TRUE;

  /* default */
  sc = cv.scene[currentScene( )];
  mesh = cv.mesh[sc->idmesh];
  clip = sc->clip;

  switch (key) {
    case 'C': /* toggle clipping */
      if (clip->active & C_ON) {
        clip->active &= ~(C_ON + C_EDIT);
      } else {
        clip->active |= C_ON;
        if (mesh->ntet + mesh->nhex) clip->active |= C_VOL;
      }

      break;
    case 'E': /* edit clip plane */
      if (!(clip->active & C_ON)) return;

      if (clip->active & C_FREEZE) clip->active ^= C_FREEZE;

      clip->active ^= C_EDIT;
      break;
    case 'F': /* freeze clip */
      if (!(clip->active & C_ON)) return;

      if (clip->active & C_EDIT) clip->active ^= C_EDIT;

      clip->active ^= C_FREEZE;
      break;
    case 'H': /* toggle draw plane */
      clip->active ^= C_HIDE;
      break;
    case 'I': /* inverse orientation */
      invertClip(sc, clip);
      clip->active |= C_REDO;
      break;
    case 'K': /* toggle capping */
      clip->active ^= C_CAP;
      clip->active |= C_REDO;
      break;
    case 'R': /* reset clip */
      if (!(clip->active & C_ON)) break;

      resetClip(sc, clip, mesh);
      break;
    case 'Z': /* toggle volclip */
      clip->active ^= C_VOL;
      if (clip->active & C_VOL) clip->active |= C_REDO;

      break;
  }

  if (post) glutPostRedisplay( );
}

void menuClip(int item) { menuDefer(keyClip, item); }

void keyCube(unsigned char key, int x, int y) {
  pScene sc;
  pMesh mesh;
  pCube cube;
  ubyte post = TRUE;

  /* default */
  sc = cv.scene[currentScene( )];
  mesh = cv.mesh[sc->idmesh];
  cube = sc->cube;

  switch (key) {
    case 'C': /* toggle clipping */
      if (cube->active & C_ON)
        cube->active &= ~(C_ON + C_EDIT);
      else
        cube->active |= C_ON;

      break;
    case 'E': /* edit clip plane */
      if (!(cube->active & C_ON)) return;

      if (cube->active & C_FREEZE) cube->active ^= C_FREEZE;

      cube->active ^= C_EDIT;
      break;
    case 'F': /* freeze cube */
      if (!(cube->active & C_ON)) return;

      if (cube->active & C_EDIT) cube->active ^= C_EDIT;

      cube->active ^= C_FREEZE;
      break;
    case 'R': /* reset cube */
      if (!(cube->active & C_ON)) break;

      resetCube(sc, cube, mesh);
      break;
  }

  if (post) glutPostRedisplay( );
}

void keyFeature(unsigned char key, int x, int y) {
  pScene sc;
  pMesh mesh;
  ubyte post = TRUE, dolist = TRUE;

  /* default */
  sc = cv.scene[currentScene( )];
  mesh = cv.mesh[sc->idmesh];

  switch (key) {
    case 'S': /* shrink mode */
    case 'V': /* volumic shrink */
      if (sc->shrink < 0.99)
        sc->shrink = 1.0f;
      else
        sc->shrink = 0.95;

      dolist = TRUE;
      break;
    case 'I': /* increase shrink value */
      if (sc->shrink > 0.99) break;

      sc->shrink -= .05;
      if (sc->shrink < 0.1) sc->shrink = 0.1;

      dolist = TRUE;
      break;
    case 'i': /* decrease shrink value */
      if (sc->shrink > 0.99) break;

      sc->shrink += .05;
      if (sc->shrink > 0.95) sc->shrink = 0.95;

      dolist = TRUE;
      break;
    case 's': /* scissor mode */
      if (mesh->dim == 2) break;

      sc->type ^= S_SCISSOR;
      if (sc->type & S_SCISSOR) {
        glutDisplayFunc(scissorScene);
      } else {
        glutDisplayFunc(redrawScene);
        reshapeScene(sc->par.xs, sc->par.ys);
      }

      break;
  }

  if (dolist == TRUE) {
    doLists(sc, mesh);
    if (sc->mode & S_MAP) doMapLists(sc, mesh, 1);
  }

  if (post == TRUE) glutPostRedisplay( );
}

void menuFeature(int item) { menuDefer(keyFeature, item); }

void menuImage(int item) { imgtype = item; }

void keyMetric(unsigned char key, int x, int y) {
  pScene sc;
  pMesh mesh;
  ubyte post = TRUE;
  int ret;

  /* default */
  sc = cv.scene[currentScene( )];
  mesh = cv.mesh[sc->idmesh];

  switch (key) {
    case 'c': /* critical points */
      if (!mesh->nbb) return;

      sc->isotyp ^= S_CRITP;
      doIsoLists(sc, mesh, 0);
      break;
    case 'f': /* flush streamlines */
      if (sc->stream->nbstl) {
        int k, kk;

        for (kk = 0; kk < sc->stream->nbstl; kk++) {
          if (sc->slist[kk]) glDeleteLists(sc->slist[kk], 1);

          sc->slist[kk] = (GLuint)0;
        }

        sc->stream->nbstl = 0;

        for (k = 1; k <= mesh->np; k++) {
          pPoint ppt;

          ppt = &mesh->point[k];
          ppt->flag = 0;
        }
      }

      break;
    case 'l': /* iso-lines */
      if (!mesh->nbb) return;

      sc->isotyp ^= S_ISOLINE;
      doIsoLists(sc, mesh, 0);
      break;
    case 's': /* iso-surfaces */
      if (!mesh->nbb) return;

      sc->isotyp ^= S_ISOSURF;
      doIsoLists(sc, mesh, 0);
      break;
    case 'u': /* field lines animation */
      if (!mesh->nbb || mesh->nfield != mesh->dim) return;

      if (refitem) {
        sc->isotyp |= S_PARTICLE;
        createParticle(sc, mesh);
      }

      glutIdleFunc(streamIdle);
      break;
    case 'd': /* displacement */
      if (!mesh->nbb || mesh->nfield != mesh->dim) return;

      sc->mode ^= S_DISPL;
      meshCoord(mesh, (sc->mode & S_DISPL) ? 1 : 0);
      meshBox(mesh, 1);
      doLists(sc, mesh);
      if (sc->mode & S_MAP) doMapLists(sc, mesh, 1);

      if (sc->isotyp) doIsoLists(sc, mesh, 1);

      break;

    case 'v': /* streamlines */
      if (!mesh->nbb || mesh->nfield != mesh->dim) return;

      if (refitem) {
        sc->isotyp |= S_STREAML;
        doIsoLists(sc, mesh, 0);
      } else {
        if (!streamIsoPoint(sc, mesh)) {
          post = FALSE;
          sc->isotyp &= ~S_STREAML;
        }
      }

      break;
    case 'w': /*vector/tensor */
      if (mesh->nfield != mesh->dim) {
        if (mesh->dim == 2 && mesh->nfield == 3) {
          if (sc->picklist) {
            glDeleteLists(sc->picklist, 1);
            sc->picklist = 0;
            break;
          } else {
            sc->picklist = drawAllEllipse(sc, mesh);
          }

          break;
        }

        return;
      }

      sc->isotyp ^= S_VECTOR;
      doIsoLists(sc, mesh, 0);
      break;

    case 'k': /* toggle elevation */
      if (mesh->dim != 2) return;

      sc->mode ^= S_ALTITUDE;
      if (altcoef == 0.0) {
        float maxd;

        float maxb;

        maxd = max(mesh->xmax - mesh->xmin, mesh->ymax - mesh->ymin);
        /* the greatest absolute value (bbmax can be null or negative) */
        maxb = max(fabs(mesh->bbmin), fabs(mesh->bbmax));
        altcoef = maxb > 0.0 ? 0.3 * maxd / maxb : 1.0;
      }

      if (!(sc->mode & S_ALTITUDE)) sc->type |= S_RESET;

      doMapLists(sc, mesh, 1);
      /* the iso-lines are on the elevation map */
      if (sc->isotyp & S_ISOLINE) doIsoLists(sc, mesh, 1);

      break;
    case 'K': /* elevation coeff */
      fprintf(stdout, "elevation coeff (%.2f): ", altcoef);
      fflush(stdout);
      while (fgetc(stdin) != EOF)
        ; /*fflush() called on input stream 'stdin' may result in undefined behaviour on non-linux
             systems*/
      ret = fscanf(stdin, "%f", &altcoef);
      if (ret == EOF) printf("fscanf error\n");
      if (altcoef == 0.0) sc->mode |= ~S_ALTITUDE;

      sc->type |= S_RESET;
      doMapLists(sc, mesh, 1);
      doIsoLists(sc, mesh, 1);
      break;

    case 'm': /* display metric */
      if (!mesh->nbb) return;

      sc->mode ^= S_MAP;
      if (sc->mode & S_MAP) sc->mode |= S_FILL;

      doMapLists(sc, mesh, 1);
      /* the iso lines: dark on the colors (as in ffglut), in color without them */
      {
        int i;

        for (i = 0; i < MAX_LIST; i++) {
          if (sc->ilist[i]) glDeleteLists(sc->ilist[i], 1);

          sc->ilist[i] = 0;
        }

        if (sc->isotyp) doIsoLists(sc, mesh, 0);
      }
      /* as in ffglut: the palette is shown only with the key v (not with the colors) */
      if (sc->mode & S_MAP && sc->clip->active & C_ON) sc->clip->active |= C_REDO;

      break;
    case 'p': /* toggle palette */
      if (!mesh->nbb) return;

      sc->item ^= S_PALETTE;
      break;

    default:
      post = FALSE;
      break;
  }

  if (post) glutPostRedisplay( );
}

void menuMetric(int item) { menuDefer(keyMetric, item); }

/* menu entry = the same action as the key */
void menuKey(int item) { menuDefer(keySceneMenu, item); }

/* menu entry = the old action of medit for this key (its letter is now a key of ffglut) */
static void keySceneOld(unsigned char key, int x, int y) {
  extern int menukey, oldkey;

  menukey = oldkey = 1;
  keyScene(key, x, y);
  menukey = oldkey = 0;
}

void menuOld(int item) { menuDefer(keySceneOld, item); }

int createMenus(pScene sc, pMesh mesh) {
  int amenu, fmenu, femenu, vmenu, mmenu, smenu;
  int clmenu, cmenu, vwmenu, gmenu = 0;
  extern char *ffserver; /* FreeFEM server mode: n and p = next and previous plot */

  /* default */
  if (ddebug) printf("create menus\n");

  smenu = 0;

  /* File management menu */
  fmenu = glutCreateMenu(menuFile);
  glutAddMenuEntry("    Load prefs", 'L');
  glutAddMenuEntry("    Save prefs", 'W');
  glutAddMenuEntry("    Update mesh", 'R');
  glutAddMenuEntry("    Save mesh", 'S');
  glutAddMenuEntry("[w] Hardcopy PPM", 'H');
  glutAddMenuEntry("    Hardcopy EPS (Color)", 'C');
  glutAddMenuEntry("    Hardcopy EPS (Grey)", 'G');
  glutAddMenuEntry("    Hardcopy EPS (B/W)", 'B');
  glutAddMenuEntry("    Softcopy EPS", 's');

  /* rendering mode selector */
  mmenu = glutCreateMenu(menuMode);
  glutAddMenuEntry("Wireframe", 'W');
  glutAddMenuEntry("Depth  lines", 'D');
  glutAddMenuEntry("Hidden lines", 'H');
  glutAddMenuEntry("Shading", 'S');
  glutAddMenuEntry("Shading+lines", 'P');
  if (mesh->nvn > 0) glutAddMenuEntry("    Smooth / flat normals", 'n');

  glutCreateMenu(menuKey);
  glutAddMenuEntry("[L] Toggle facets", 'L');
  glutAddMenuEntry("[m] Toggle mesh", 'm');
  glutAddMenuEntry("[l] Toggle lighting", 'l');
  {
    int kmenu = glutGetMenu( );

    glutSetMenu(mmenu);
    glutAddSubMenu("Facets, lines", kmenu);
  }

  /* color & material menu */
  cmenu = glutCreateMenu(menuColor);
  glutAddMenuEntry("[b] Toggle backcolor", 'b');
  glutAddMenuEntry("    Toggle matcolors", 'e');
  glutAddMenuEntry("    Edit   matcolors", 'E');
  glutAddMenuEntry("    Hide   material", 'r');
  glutAddMenuEntry("    Reset  materials", 'R');

  /* metric */
  if (mesh->nbb > 0) {
    smenu = glutCreateMenu(menuMetric);
    glutAddMenuEntry("[f] Fill (solution in colors)", 'm');
    glutAddMenuEntry("[v] Toggle palette", 'p');
    if (mesh->typage == 2) glutAddMenuEntry("[o] Toggle iso-lines", 'l');

    if (mesh->ntet + mesh->nhex > 0 && mesh->nfield == 1)
      glutAddMenuEntry("[s] Toggle iso-surfaces", 's');

    if (mesh->nfield == mesh->dim) {
      glutAddMenuEntry("[V]  Toggle vector/tensor", 'w');
      {
        int dm = glutGetMenu( ), am = glutCreateMenu(menuKey);

        glutAddMenuEntry("[a] Shorter arrows", 'a');
        glutAddMenuEntry("[A] Longer arrows", 'A');
        glutSetMenu(dm);
        glutAddSubMenu("Arrow length", am);
      }
      glutAddMenuEntry("     Toggle displacement", 'd');
      glutAddMenuEntry("     Toggle streamlines", 'v');
      glutAddMenuEntry("     Flush streamlines", 'f');
      glutAddMenuEntry("     Particle advection", 'u');
      if (mesh->dim == 2) glutAddMenuEntry("     Critical points", 'c');
    }

    if (mesh->cplx) {
      int dm = glutGetMenu( ), cm = glutCreateMenu(menuKey);

      glutAddMenuEntry("[k] Real / imaginary part / modulus / argument", 'k');
      glutSetMenu(dm);
      glutAddSubMenu("Complex solution", cm);
    }

    if (mesh->dim == 2) {
      glutAddMenuEntry("[3]  Toggle elevation", 'k');
      glutAddMenuEntry("[K]  Elevation coeff (h, H: lower, higher)", 'K');
      glutAddMenuEntry("[c]  Cut line u(x) / u(y) / none (shift+drag: segment)", 'c');
      glutAddMenuEntry("[<]  Move the cut line down / left", '<');
      glutAddMenuEntry("[>]  Move the cut line up / right", '>');
    }

    glutCreateMenu(menuKey);
    glutAddMenuEntry("[P] Value under the mouse", 'P');
    glutAddMenuEntry("[U] Overlay of the other windows", 'U');
    {
      int kmenu = glutGetMenu( );

      glutSetMenu(smenu);
      glutAddSubMenu("Probe, overlay", kmenu);
    }
  }

  /* 1d graph (solution on a mesh of edges along x) */
  if (graph1dMesh(mesh)) {
    gmenu = glutCreateMenu(menuKey);
    glutAddMenuEntry("[k] Lines / lines+points / points", 'k');
    glutAddMenuEntry("[G] Grid", 'G');
    glutAddMenuEntry("[F] Bold text", 'F');
    glutAddMenuEntry("[x] Log scale in x", 'x');
    glutAddMenuEntry("[y] Log scale in y", 'y');
    glutAddMenuEntry("[z] Zoom in", 'z');
    glutAddMenuEntry("[Z] Zoom out", 'Z');
    glutAddMenuEntry("[A] Freeze / free the axes", 'A');
    glutAddMenuEntry("[i] Automatic axes", 'i');
    glutAddMenuEntry("[b] Black / white background", 'b');
  }

  /* Show misc. items */
  vmenu = glutCreateMenu(menuItem);
  glutAddMenuEntry("    Axis (arrows)", 'A');
  glutAddMenuEntry("    Bounding box", 'B');
  glutAddMenuEntry("[j] Graduated axes", 'T');
  if (mesh->dim == 3) glutAddMenuEntry("[D] Selection cube", 'D');

  glutAddMenuEntry("[G] Grid ", 'G');
  if (sc->glist) glutAddMenuEntry("[B] Border (geometric items)", 'g');

  glutAddMenuEntry("    Toggle Info", 'j');
  glutAddMenuEntry("    Toggle Point num", 'P');
  glutAddMenuEntry("    Toggle Face num", 'F');
  if (mesh->nvn) {
    glutAddMenuEntry("[T] Toggle normals", 'N');
    glutAddMenuEntry("    Revert normals", 'O');
  }

  /* clipping menu */
  clmenu = glutCreateMenu(menuClip);
  glutAddMenuEntry("[F1] Toggle clip", 'C');
  glutAddMenuEntry("[F2] Edit clip", 'E');
  glutAddMenuEntry("[F3] Freeze clip", 'F');
  glutAddMenuEntry("     Inverse orient", 'I');
  if (mesh->ntet + mesh->nhex > 0) glutAddMenuEntry("     Toggle capping", 'K');

  glutAddMenuEntry("     Toggle plane", 'H');
  glutAddMenuEntry(" --  Reset clip", 'R');
  if (mesh->ntet + mesh->nhex > 0) {
    sc->clip->active |= C_VOL;
    glutAddMenuEntry("[F4] Toggle Vclip", 'Z');
  }

  /* feature menu */
  femenu = glutCreateMenu(menuFeature);
  if (mesh->ntet + mesh->nhex > 0)
    glutAddMenuEntry("[F5] Toggle Vshrink", 'V');
  else
    glutAddMenuEntry("[F5] Toggle shrink", 'S');

  glutAddMenuEntry("[F6] Increase shrink", 'I');
  glutAddMenuEntry("[F7] Decrease shrink", 'i');
  if (mesh->dim == 3) glutAddMenuEntry("Toggle splitview", 's');

  /* view handler menu */
  vwmenu = glutCreateMenu(menuView);
  glutAddMenuEntry("[=] Reset", 'R');
  glutAddMenuEntry("[Alt-c] Copy", 'C');
  glutAddMenuEntry("[Alt-p] Paste", 'P');
  glutAddMenuEntry("[Alt+l] Link", 'L');
  glutAddMenuEntry("[Alt+u] Unlink", 'U');
  glutAddMenuEntry("[Alt+d] Duplicate", 'D');

  glutCreateMenu(menuOld);
  glutAddMenuEntry("[z] Zoom in", 'z');
  glutAddMenuEntry("[Z] Zoom out", 'Z');
  glutAddMenuEntry("[+] Scale object up", '+');
  glutAddMenuEntry("[-] Scale object down", '-');
  glutAddMenuEntry("    Toggle interactive", 'I');
  if (mesh->dim == 3 || sc->mode & S_ALTITUDE) glutAddMenuEntry("[J] Toggle flight", 'J');

  glutAddMenuEntry("    Center on the picked item", 'V');
  glutAddMenuEntry("[#] Select entity (terminal)", '#');
  {
    int kmenu = glutGetMenu( );

    glutSetMenu(vwmenu);
    glutAddSubMenu("Zoom, center", kmenu);
  }

  /* animation menu */
  amenu = glutCreateMenu(menuAnim);
  glutAddMenuEntry("    Toggle Anim", 'A');
  glutAddMenuEntry("[S] Play animation", 'S');
  glutAddMenuEntry("Toggle ImgSave", 'I');
  if (option == SEQUENCE || option == SEQUENCE + PARTICLE) {
    glutAddMenuEntry("Play  sequence", 'S');
    glutAddMenuEntry("First mesh", 'f');
    glutAddMenuEntry("Last  mesh", 'l');
    glutAddMenuEntry("Next mesh", 'n');
    glutAddMenuEntry("Prev mesh", 'p');
  } else if (option == MORPHING) {
    glutAddMenuEntry("Start/Stop morph.", 'M');
    glutAddMenuEntry("Toggle AutoReverse", 'R');
  }

  /* trajectory menu */
  if (mesh->dim == 3 || mesh->nbb) {
    glutCreateMenu(menuTrajet);
    glutAddMenuEntry("New Ctrl point", 'C');
    glutAddMenuEntry("Toggle path", 'S');
    glutAddMenuEntry("Follow path", 'F');
    glutAddMenuEntry("Load path", 'L');
    glutAddMenuEntry("Save path", 'W');
  }

  /* main scene menu */
  glutCreateMenu(menuScene);
  glutAddSubMenu("File", fmenu);
  glutAddSubMenu("Render mode", mmenu);
  glutAddSubMenu("Colors, Materials", cmenu);
  if (mesh->nbb) glutAddSubMenu("Data", smenu);

  if (gmenu) glutAddSubMenu("Graph 1D", gmenu);

  glutAddSubMenu("Items", vmenu);
  if (mesh->dim == 3) glutAddSubMenu("Clipping", clmenu);

  glutAddSubMenu("Features", femenu);
  glutAddSubMenu("View", vwmenu);
  glutAddSubMenu("Animation", amenu);
  glutAddMenuEntry("", '\0');
  if (ffserver) {
    glutAddMenuEntry("[Enter] Next plot (or continue)", '\r');
    glutAddMenuEntry("[p] Previous plot", 'p');
  }

  glutAddMenuEntry("[?] Help (in the terminal)", '?');
  glutAddMenuEntry("    Close window", 'X');
  glutAddMenuEntry("[Esc] Quit", 27);

  return (1);
}

#ifdef __cplusplus
}
#endif
