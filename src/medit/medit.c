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
/* SUMMARY : ...  visualization tool                                        */
/* LICENSE : LGPLv3                                                         */
/* ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE           */
/* AUTHORS : Pascal Frey                                                    */
/* E-MAIL  : pascal.frey@sorbonne-universite.fr                             */

#ifdef __cplusplus
extern "C" {
#endif

#include "medit.h"
#include "compil.date"
#ifdef ppc
#include <unistd.h>
#endif

/* global variables (see extern.h) */
GLboolean hasStereo = 1;
Canvas cv;
mytime ctim[TIMEMAX];
ubyte ddebug, animate, saveimg, imgtype, infogl, fullscreen;
ubyte quiet, option, morphing, stereoMode;
int menu, amenu, fmenu, femenu, vmenu, mmenu, smenu;
int clmenu, cmenu, vwmenu, txmenu, trmenu;
int animdep, animfin;

/*  Rajout pour popen */
ubyte dpopen, dpopensol, dpopenbin;
static void excfun(int sigid) {
  fprintf(stdout, "\n Unexpected error:");
  fflush(stdout);

  switch (sigid) {
    case SIGABRT:
      fprintf(stdout, "  Abnormal stop\n");
      break;
    case SIGFPE:
      fprintf(stdout, "  Floating-point exception\n");
      break;
    case SIGILL:
      fprintf(stdout, "  Illegal instruction\n");
      break;
    case SIGSEGV:
      fprintf(stdout, "  Segmentation fault\n");
      break;
    case SIGTERM:
    case SIGINT:
      fprintf(stdout, "  Program killed\n");
      break;
  }

  exit(1);
}

static void endcod( ) {
  chrono(OFF, &ctim[0]);
  fprintf(stdout, "\n Total running seconds:  %.2f\n", gttime(ctim[0]));
  fprintf(stdout, " Thank you for using Medit.\n");
}

static void grInfo(void) {
  GLboolean b;
  GLint i;

  glutCreateWindow("Info");
  fprintf(stdout, "Graphic info:\n");
  fprintf(stdout, " GL Vendor:\t%s\n", glGetString(GL_VENDOR));
  fprintf(stdout, " GL Version:\t%s\n", glGetString(GL_VERSION));
  fprintf(stdout, " GL Renderer:\t%s\n\n", glGetString(GL_RENDERER));
  glGetBooleanv(GL_RGBA_MODE, &b);
  if (b) fprintf(stdout, "  RGBA Mode\n");

  glGetBooleanv(GL_DOUBLEBUFFER, &b);
  if (b) fprintf(stdout, "  Double Buffer\n");

  glGetBooleanv(GL_STEREO, &b);
  if (b) fprintf(stdout, "  Stereo\n");

  glGetIntegerv(GL_AUX_BUFFERS, &i);
  if (i) fprintf(stdout, "  Auxilary Buffers\t%2d\n", (int)i);

  glGetIntegerv(GL_INDEX_BITS, &i);
  if (i) fprintf(stdout, "  Index Bits\t\t%2d\n", (int)i);

  glGetIntegerv(GL_RED_BITS, &i);
  fprintf(stdout, "  RGBA Bits\t\t%2d", (int)i);
  glGetIntegerv(GL_GREEN_BITS, &i);
  fprintf(stdout, "\t%2d", (int)i);
  glGetIntegerv(GL_BLUE_BITS, &i);
  fprintf(stdout, "\t%2d", (int)i);
  glGetIntegerv(GL_ALPHA_BITS, &i);
  fprintf(stdout, "\t%2d\n", (int)i);
  glGetIntegerv(GL_ACCUM_RED_BITS, &i);
  fprintf(stdout, "  Accum RGBA Bits\t%2d", (int)i);
  glGetIntegerv(GL_ACCUM_GREEN_BITS, &i);
  fprintf(stdout, "\t%2d", (int)i);
  glGetIntegerv(GL_ACCUM_BLUE_BITS, &i);
  fprintf(stdout, "\t%2d", (int)i);
  glGetIntegerv(GL_ACCUM_ALPHA_BITS, &i);
  fprintf(stdout, "\t%2d\n", (int)i);
  glGetIntegerv(GL_DEPTH_BITS, &i);
  fprintf(stdout, "  Depth Bits\t\t%2d\n", (int)i);
  glGetIntegerv(GL_STENCIL_BITS, &i);
  fprintf(stdout, "  Stencil Bits\t\t%2d\n", (int)i);

  exit(1);
}

int medit0( ) {
  char data[128];
  int k, l;
  clock_t ct;
  char *res = "";

  /* default */
  fprintf(stdout, " Loading data file(s)\n");
  ct = clock( );

  /* enter name */
  if (!cv.nbm) {
    char *name;

    fprintf(stdout, "  File name(s) missing. Please enter : ");
    fflush(stdout);
    while (fgetc(stdin) != EOF)
      ;    // fflush() called on input stream 'stdin' may result in undefined behaviour on non-linux
           // systems
    res = fgets(data, 120, stdin);
    if (res == NULL) printf("fgets error\n");
    if (!strlen(data)) {
      fprintf(stdout, "  ## No data\n");
      return (0);
    }

    /* parse file name(s) */
    name = strtok(data, " \n");

    while (name) {
      if (!cv.mesh[cv.nbm]) {
        cv.mesh[cv.nbm] = (pMesh)M_calloc(1, sizeof(Mesh), "medit0.mesh");
        if (!cv.mesh[cv.nbm]) return (0);
      }

      strcpy(cv.mesh[cv.nbm]->name, name);
      name = strtok(NULL, " \n\0");
      if (++cv.nbm == MAX_MESH) break;
    }

    if (!cv.nbm) return (0);
  }

  if (!cv.nbm) {
    fprintf(stdout, "  Number of mesh missing:. Please enter : ");
    fflush(stdout);
    while (fgetc(stdin) != EOF)
      ;    // fflush() called on input stream 'stdin' may result in undefined behaviour on non-linux
           // systems
    res = fgets(data, 120, stdin);
    if (res == NULL) printf("fgets error\n");
    cv.nbm = atoi(data);
  }

  /* read mesh(es) */
  k = 0;

  do {
    pMesh mesh;
    int ret;

    if (!cv.mesh[k]) {
      cv.mesh[k] = M_calloc(1, sizeof(Mesh), "medit0.mesh");
      if (!cv.mesh[k]) return (0);
    }

    mesh = cv.mesh[k];
    mesh->typ = 0;
    ret = loadMesh(mesh);
    if (ret < 0) {
      mesh->typ = 1;
      ret = inmsh2(mesh);
      if (!ret) {
        mesh->typ = 2;
        ret = loadGIS(mesh);
      }
    }

    if (ret <= 0) {
      for (l = k + 1; l < cv.nbm; l++) {
        cv.mesh[l - 1] = cv.mesh[l];
      }

      cv.nbm--;
      k--;
      continue;
    }

    /* compute mesh box */
    if ((mesh->ntet && !mesh->nt) || (mesh->nhex && !mesh->nq)) meshSurf(mesh);

    meshBox(mesh, 1);
    if (!quiet) meshInfo(mesh);

    /* read metric  */
    if (!loadSol(mesh, mesh->name, 1)) bbfile(mesh);

    if (!quiet && mesh->nbb) fprintf(stdout, "    Solutions  %8d\n", mesh->nbb);
  } while (++k < cv.nbm);

  cv.nbs = cv.nbm;

  ct = difftime(clock( ), ct);
  fprintf(stdout, "  Input seconds:     %.2f\n", (double)ct / (double)CLOCKS_PER_SEC);
  ;
  return (cv.nbm);
}

int medit0_popen( ) {
  int k;
  clock_t ct;

  /* default */
  fprintf(stdout, " Loading data file(s)\n");
  ct = clock( );

  /* enter number of mesh */
  if (!cv.nbm) {
    char data[128];
    char *res;

    fprintf(stdout, "  Number of mesh missing:. Please enter : ");
    fflush(stdout);
    while (fgetc(stdin) != EOF)
      ;    // fflush() called on input stream 'stdin' may result in undefined behaviour on non-linux
           // systems
    res = fgets(data, 128, stdin);
    if (res == NULL) printf("fgets error\n");
    cv.nbm = atoi(data);
  }

  /* read mesh(es) */
  k = 0;

  do {
    pMesh mesh;

    if (!cv.mesh[k]) {
      cv.mesh[k] = M_calloc(1, sizeof(Mesh), "medit0.mesh");
      if (!cv.mesh[k]) return (0);
    }

    mesh = cv.mesh[k];
    mesh->typ = 0;

    if (dpopenbin)
      loadMesh_popen_bin(mesh);
    else
      loadMesh_popen(mesh);

    /* compute mesh box */
    if ((mesh->ntet && !mesh->nt) || (mesh->nhex && !mesh->nq)) meshSurf(mesh);

    meshBox(mesh, 1);
    if (!quiet) meshInfo(mesh);

    if (dpopensol) {
      if (dpopenbin)
        loadSol_popen_bin(mesh, mesh->name, 1);
      else
        loadSol_popen(mesh, mesh->name, 1);
    }

    if (!quiet && mesh->nbb) fprintf(stdout, "    Solutions  %8d\n", mesh->nbb);
  } while (++k < cv.nbm);

  cv.nbs = cv.nbm;

  ct = difftime(clock( ), ct);
  fprintf(stdout, "  Input seconds:     %.2f\n", (double)ct / (double)CLOCKS_PER_SEC);

  return (cv.nbm);
}

int medit1( ) {
  int k;
  clock_t ct;

  /* create grafix */
  fprintf(stdout, "\n medit1() \n");
  fprintf(stdout, "\n Building scene(s)\n");
  ct = clock( );

  for (k = 0; k < cv.nbs; k++) {
    pScene scene;
    pMesh mesh;

    if (!cv.scene[k]) {
      cv.scene[k] = (pScene)M_calloc(1, sizeof(Scene), "medit1.scene");
      if (!cv.scene[k]) return (0);
    }

    scene = cv.scene[k];
    if (!cv.mesh[k]) {
      cv.mesh[k] = (pMesh)M_calloc(1, sizeof(Mesh), "medit1.mesh");
      if (!cv.mesh[k]) return (0);
    }

    mesh = cv.mesh[k];

    fprintf(stdout, "  Creating scene %d\n", k + 1);
    parsop(scene, mesh);
    meshRef(scene, mesh);
    matSort(scene);

    if (option == ISOSURF) {
      if (!mesh->nbb) return (0);

      setupPalette(scene, mesh);
      tetraIsoPOVray(scene, mesh);
    } else if (!createScene(scene, k)) {
      fprintf(stderr, "  ## Unable to create scene\n");
      return (0);
    }
  }

  ct = difftime(clock( ), ct);
  fprintf(stdout, "  Scene seconds:     %.2f\n", (double)ct / (double)CLOCKS_PER_SEC);

  return (1);
}

int main(int argc, char *argv[]) {
  int type;

#ifdef ppc
  if (!getwd(pwd)) exit(2);

#endif

  fprintf(stdout, "  -- Medit,  Release %s (%s)\n", ME_VER, ME_REL);
  fprintf(stdout, "     %s.\n", ME_CPY);
  fprintf(stdout, "     compiled: %s.\n\n", COMPIL);

  /* trap exceptions */
  signal(SIGABRT, excfun);
  signal(SIGFPE, excfun);
  signal(SIGILL, excfun);
  signal(SIGSEGV, excfun);
  signal(SIGTERM, excfun);
  signal(SIGINT, excfun);
  atexit(endcod);

  tminit(ctim, TIMEMAX);
  chrono(ON, &ctim[0]);

  /* default values */
  option = STANDARD;
  saveimg = GL_FALSE;
  imgtype = P6;
  animate = FALSE;
  morphing = FALSE;
  fullscreen = FALSE;
  animdep = 0;
  animfin = 0;
  ddebug = FALSE;
  quiet = TRUE;    // FALSE;
  stereoMode = 0;
  cv.nbm = cv.nbs = 0;

  /* default value for popen */
  dpopen = FALSE;
  dpopenbin = FALSE;
  dpopensol = FALSE;

  /* init grafix */
  parsar(argc, argv);

  if (option == ISOSURF) {
    fprintf(stdout, "ISOSURF");
    if (!medit0( )) exit(1);

    if (!medit1( )) exit(1);

    return (0);
  }

  glutInit(&argc, argv);
#ifdef ppc
  chdir(pwd);
#endif

  chrono(ON, &ctim[0]);
  if (stereoMode == MONO)
    type = GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH;
  else
    type = GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH | GLUT_STEREO;

  glutInitDisplayMode(type);
  if (infogl) grInfo( );

  /* call animate or normal mode */
  if (dpopen == FALSE) {
    switch (option) {
      case STANDARD:

      case SCHNAUZER:
        if (!medit0( )) exit(1);

        if (!medit1( )) exit(1);

        break;
      case SEQUENCE:
        if (!animat( )) exit(1);

        break;
      case SEQUENCE + PARTICLE:

        if (!animat( )) exit(1);

        break;
      case MORPHING:

        if (!medit0( )) exit(1);

        if (!modeMorphing( )) exit(1);

        morphing = GL_FALSE;
        break;
      default:
        fprintf(stderr, "  ## Unrecognized option %d\n", option);
        exit(1);
        break;
    }
  } else {
    switch (option) {
      case STANDARD:

      case SCHNAUZER:

        if (!medit0_popen( )) exit(1);

        if (!medit1( )) exit(1);

        break;

      case SEQUENCE:

        if (!animat( )) exit(1);

        break;
      case SEQUENCE + PARTICLE:

        if (!animat( )) exit(1);

        break;
      case MORPHING:

        if (!medit0_popen( )) exit(1);

        if (!modeMorphing( )) exit(1);

        morphing = GL_FALSE;
        break;
      default:
        fprintf(stderr, "  ## Unrecognized option %d\n", option);
        exit(1);
        break;
    }
  }

  /* main grafix loop */
  fprintf(stdout, "\n Rendering scene(s)\n");
  glGetBooleanv(GL_STEREO, &hasStereo);
  glutMainLoop( );

  return (0);
}

#ifdef __cplusplus
}
#endif
