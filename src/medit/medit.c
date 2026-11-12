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

/*
 *  Server mode: one medit for all the plots of a FreeFEM script.
 *  FreeFEM launches a new process (popen) for each call of medit(...) and waits its end.
 *  Now this process is only a client: it stores the data in a temporary file, gives the file
 *  to a medit server (started the first time) through a unix socket, and stops when the server
 *  says to continue (at once without wait, otherwise when the user asks for the next plot).
 *  The server keeps the last plots, and shows them in the same window(s):
 *    n or return = next plot (or continue the script),  p = previous plot.
 *  Options:  -wait 0|1 (or the environment variable FFMEDIT_WAIT), default 1
 *            the environment variable FFMEDIT_NOSERVER = the old behaviour (a medit by plot)
 */
char *ffserver = NULL; /* -server path : name of the socket of the server */
int ffwait = -1;       /* -wait 0|1 */
extern char graphopt[1024]; /* -graph options of the 1d graph (see items.c) */

#ifndef WIN32
#include <unistd.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <sys/un.h>
#include <sys/select.h>
#include <sys/time.h>
#include <fcntl.h>
#include <errno.h>

#define SRV_MAXPLOT 30
typedef struct {
  char file[256];           /* temporary file with the data of the plot */
  char opt[1024];           /* options of the 1d graph */
  char name[MAX_MESH][256]; /* name of the meshes */
  int bin, sol, nbm;
} SrvPlot;
static SrvPlot srvplot[SRV_MAXPLOT];
static int srvnb = 0;      /* number of plots kept */
static int srvfirst = 0;   /* number of plots forgotten */
static int srvcur = -1;    /* the plot on the screen */
static int srvlisten = -1; /* socket of the server */
static int srvwait = -1;   /* socket of the client which waits to continue the script */
static char srvsock[128] = "";
static char srvtmp[256] = ""; /* temporary file of a client without server */

static int srvConnect(const char *path) {
  struct sockaddr_un addr;
  int fd = socket(AF_UNIX, SOCK_STREAM, 0);

  if (fd < 0) return (-1);

  memset(&addr, 0, sizeof(addr));
  addr.sun_family = AF_UNIX;
  strncpy(addr.sun_path, path, sizeof(addr.sun_path) - 1);
  if (connect(fd, (struct sockaddr *)&addr, sizeof(addr)) < 0) {
    close(fd);
    return (-1);
  }

  return (fd);
}

static void srvCleanup( ) {
  int i;

  for (i = 0; i < srvnb; i++) unlink(srvplot[i].file);

  if (srvsock[0]) unlink(srvsock);

  if (srvtmp[0]) unlink(srvtmp);
}

/* client: return 0 if the plot is done by the server, -1 if no server (then stdin = the data) */
static int srvClient(char *prog) {
  FILE *f;
  char path[128], buf[65536], msg[MAX_MESH * 260 + 2048];
  size_t n;
  int fd, k, w;

  /* the data in a temporary file */
  snprintf(srvtmp, sizeof(srvtmp), "/tmp/ffmedit-%d-%ld.dat", (int)getpid( ), (long)time(NULL));
  f = fopen(srvtmp, "wb");
  if (!f) {
    srvtmp[0] = '\0';
    return (-1);
  }

  while ((n = fread(buf, 1, sizeof(buf), stdin)) > 0) fwrite(buf, 1, n, f);

  fclose(f);

  /* the server of this FreeFEM run (same process group), started if necessary */
  snprintf(path, sizeof(path), "/tmp/ffmedit-%d-%d.sock", (int)getuid( ), (int)getpgrp( ));
  fd = srvConnect(path);
  if (fd < 0) {
    pid_t pid;

    unlink(path);
    pid = fork( );
    if (pid == 0) {
      if (!freopen("/dev/null", "r", stdin)) _exit(127);

      execlp(prog, prog, "-server", path, (char *)0);
      _exit(127);
    }

    for (k = 0; k < 150 && fd < 0; k++) {
      usleep(100000);
      fd = srvConnect(path);
    }
  }

  if (fd < 0) {
    fprintf(stderr, "  ## Unable to start the medit server, old behaviour.\n");
    if (!freopen(srvtmp, "rb", stdin)) exit(1);

    atexit(srvCleanup);
    return (-1);
  }

  /* message: PLOT file bin sol wait nbm name1 ... (separator: tabulation) */
  w = ffwait;
  if (w < 0) w = getenv("FFMEDIT_WAIT") ? atoi(getenv("FFMEDIT_WAIT")) : 1;

  snprintf(msg, sizeof(msg), "PLOT\t%s\t%d\t%d\t%d\t%d", srvtmp, (int)dpopenbin, (int)dpopensol, w,
           cv.nbm);

  for (k = 0; k < cv.nbm; k++) {
    strcat(msg, "\t");
    strcat(msg, cv.mesh[k]->name);
  }

  if (graphopt[0]) {
    char *c;

    for (c = graphopt; *c; c++)
      if (*c == '\t' || *c == '\n') *c = ' ';

    strcat(msg, "\tG=");
    strcat(msg, graphopt);
  }

  strcat(msg, "\n");
  srvtmp[0] = '\0'; /* the file is now the file of the server */
  if (write(fd, msg, strlen(msg)) < 0) {
    close(fd);
    return (0);
  }

  /* wait the answer of the server (or its end) */
  while (read(fd, buf, 1) == 1 && buf[0] != '\n')
    ;

  close(fd);
  return (0);
}

/* free the data of a mesh (see getmesh in animat.c), to read a new one */
static void srvFreeMesh(pMesh mesh) {
  char name[256];

  if (mesh->tria) M_free(mesh->tria);
  if (mesh->quad) M_free(mesh->quad);
  if (mesh->edge) M_free(mesh->edge);
  if (mesh->tetra) M_free(mesh->tetra);
  if (mesh->hexa) M_free(mesh->hexa);
  if (mesh->adja) M_free(mesh->adja);
  if (mesh->voy) M_free(mesh->voy);
  if (mesh->point) M_free(mesh->point);

  if (mesh->extra) {
    if (mesh->extra->iv) M_free(mesh->extra->nv);
    if (mesh->extra->it) M_free(mesh->extra->nt);
    if (mesh->extra->iq) M_free(mesh->extra->nq);
    if (mesh->extra->n) M_free(mesh->extra->n);
    M_free(mesh->extra);
  }

  if (mesh->sol && mesh->nbb) {
    /* the components m of the vectors, tensors and complex values */
    if (mesh->nfield > 1 || mesh->cplx) {
      int k;
      for (k = 1; k <= mesh->nbb; k++) free(mesh->sol[k].m);
    }

    M_free(mesh->sol);
  }

  strcpy(name, mesh->name);
  memset(mesh, 0, sizeof(Mesh));
  strcpy(mesh->name, name);
}

/* what the menus depend on (see createMenus): solution or not, its type, kind of mesh */
static int srvMenuSig(pMesh mesh) {
  return ((mesh->nbb > 0) + 2 * mesh->typage + 8 * mesh->nfield + 128 * (mesh->ntet + mesh->nhex > 0) +
          256 * (mesh->nvn > 0) + 512 * mesh->dim + 4096 * mesh->cplx);
}

static int srvsig[MAX_SCENE];
static char srvhidden[MAX_SCENE]; /* windows not used by the current plot */
extern int graph1dnb;

static void srvTitle( ) {
  char data[512];
  int k;

  for (k = 0; k < cv.nbs; k++) {
    pScene sc = cv.scene[k];

    if (!sc) continue;

    {
      extern int cplxpart;
      extern const char *cplxname[4];
      pMesh mesh = cv.mesh[sc->idmesh];

      snprintf(data, sizeof(data), "Medit - [%s%s%s] - plot %d/%d%s - n) next  p) previous", mesh->name,
               mesh->cplx ? ", " : "", mesh->cplx ? cplxname[cplxpart] : "", srvfirst + srvcur + 1,
               srvfirst + srvnb, srvwait >= 0 ? " (FreeFEM waits)" : "");
    }
    glutSetWindow(sc->idwin);
    glutSetWindowTitle(data);
    glutSetIconTitle(data);
  }
}

/* show the plot i: the meshes are read again in the same scenes (same windows, same view) */
static int srvShow(int i) {
  SrvPlot *p = &srvplot[i];
  int k, nbwin, olddim[MAX_MESH];

  if (!freopen(p->file, "rb", stdin)) return (0);

  dpopenbin = p->bin;
  dpopensol = p->sol;
  strcpy(graphopt, p->opt);

  for (k = 0; k < p->nbm; k++) {
    pMesh mesh;

    if (!cv.mesh[k]) {
      cv.mesh[k] = (pMesh)M_calloc(1, sizeof(Mesh), "srvShow.mesh");
      if (!cv.mesh[k]) return (0);
    }

    mesh = cv.mesh[k];
    olddim[k] = mesh->dim;
    srvFreeMesh(mesh);
    strcpy(mesh->name, p->name[k]);
    mesh->typ = 0;
    if (dpopenbin)
      loadMesh_popen_bin(mesh);
    else
      loadMesh_popen(mesh);

    /* compute mesh box */
    if ((mesh->ntet && !mesh->nt) || (mesh->nhex && !mesh->nq)) meshSurf(mesh);

    meshBox(mesh, 1);
    if (dpopensol) {
      if (dpopenbin)
        loadSol_popen_bin(mesh, mesh->name, 1);
      else
        loadSol_popen(mesh, mesh->name, 1);

      cplxMesh(mesh, k); /* complex solution: a scalar part of it (key c) */
    }
  }

  if (cv.nbm < p->nbm) cv.nbm = p->nbm;

  /* several solutions on a 1d mesh: all the curves in the same graph, so only one window */
  nbwin = p->nbm;
  graph1dnb = 0;
  if (p->nbm > 1) {
    for (k = 0; k < p->nbm && graph1dMesh(cv.mesh[k]); k++)
      ;

    if (k == p->nbm) {
      graph1dnb = p->nbm;
      nbwin = 1;
    }
  }

  for (k = 0; k < nbwin; k++) {
    pMesh mesh = cv.mesh[k];
    pScene sc = k < cv.nbs ? cv.scene[k] : NULL;
    int wx = -1, wy = -1, ws = 0, hs = 0;

    if (sc && olddim[k] != mesh->dim) {
      /* not the same dimension: new scene, with the size and the position of the old window
         (the user may have resized it) */
      glutSetWindow(sc->idwin);
      wx = glutGet(GLUT_WINDOW_X);
      wy = glutGet(GLUT_WINDOW_Y);
      ws = glutGet(GLUT_WINDOW_WIDTH);
      hs = glutGet(GLUT_WINDOW_HEIGHT);
      glutDestroyWindow(sc->idwin);
      deleteScene(sc);
      cv.scene[k] = sc = NULL;
    }

    if (!sc) {
      cv.scene[k] = (pScene)M_calloc(1, sizeof(Scene), "srvShow.scene");
      if (!cv.scene[k]) return (0);

      sc = cv.scene[k];
      parsop(sc, mesh);
      if (ws > 0 && hs > 0) {
        sc->par.xs = ws;
        sc->par.ys = hs;
      }

      meshRef(sc, mesh);
      matSort(sc);
      if (!createScene(sc, k)) {
        fprintf(stderr, "  ## Unable to create scene\n");
        return (0);
      }

      if (wx >= 0 && wy >= 0) glutPositionWindow(wx, wy);

      if (cv.nbs < k + 1) cv.nbs = k + 1;

      srvsig[k] = srvMenuSig(mesh);
    } else {
      /* see playAnim in animat.c */
      int j;
      double delta = mesh->bbmax - mesh->bbmin;

      glutSetWindow(sc->idwin);
      /* the menus depend on the data (ex. no menu Data for a mesh without solution): new menus */
      if (srvsig[k] != srvMenuSig(mesh)) {
        srvsig[k] = srvMenuSig(mesh);
        if (createMenus(sc, mesh)) glutAttachMenu(GLUT_RIGHT_BUTTON);

        ffDefaults(sc, mesh); /* other kind of data: the display of ffglut at the opening */
      }

      resetLists(sc, mesh);
      /* the values of the palette for the new solution (as a new medit) */
      for (j = 0; j < MAXISO; j++) sc->iso.val[j] = mesh->bbmin + j * delta / (MAXISO - 1);

      setupPalette(sc, mesh);
      meshRef(sc, mesh);
      matSort(sc);
      doLists(sc, mesh);
      sc->glist = geomList(sc, mesh);
      if (sc->mode & (S_MAP | S_ALTITUDE)) doMapLists(sc, mesh, 1);

      if (sc->isotyp) doIsoLists(sc, mesh, 1);

      glutPostRedisplay( );
    }
  }

  /* the windows not used by this plot are hidden */
  for (k = 0; k < cv.nbs; k++) {
    pScene sc = cv.scene[k];

    if (!sc) continue;

    if (k >= nbwin && !srvhidden[k]) {
      glutSetWindow(sc->idwin);
      glutHideWindow( );
      srvhidden[k] = 1;
    } else if (k < nbwin && srvhidden[k]) {
      glutSetWindow(sc->idwin);
      glutShowWindow( );
      srvhidden[k] = 0;
    }
  }

  /* overlay: to do when all the scenes are ready */
  for (k = 0; k < nbwin && k < cv.nbs; k++) {
    pScene sc = cv.scene[k];

    if (sc && sc->overlay) {
      glutSetWindow(sc->idwin);
      doOverlay(sc);
    }
  }

  srvcur = i;
  return (1);
}

static void srvGo( ) {
  if (srvwait >= 0) {
    if (write(srvwait, "GO\n", 3) < 0) fprintf(stderr, "  ## medit server: client lost\n");

    close(srvwait);
    srvwait = -1;
  }
}

/* a new plot from a client */
static int srvHandle(int fd) {
  SrvPlot *p;
  char msg[MAX_MESH * 260 + 2048], *tok;
  int n = 0, k, w;
  struct timeval tv = {5, 0};

  setsockopt(fd, SOL_SOCKET, SO_RCVTIMEO, &tv, sizeof(tv));

  while (n < (int)sizeof(msg) - 1 && read(fd, msg + n, 1) == 1 && msg[n] != '\n') n++;

  msg[n] = '\0';
  tok = strtok(msg, "\t");
  if (!tok || strcmp(tok, "PLOT")) {
    close(fd);
    return (0);
  }

  if (srvnb == SRV_MAXPLOT) {
    /* forget the oldest plot */
    unlink(srvplot[0].file);
    memmove(&srvplot[0], &srvplot[1], (SRV_MAXPLOT - 1) * sizeof(SrvPlot));
    srvnb--;
    srvfirst++;
  }

  p = &srvplot[srvnb];
  memset(p, 0, sizeof(SrvPlot));
  tok = strtok(NULL, "\t");
  if (tok) strncpy(p->file, tok, sizeof(p->file) - 1);

  tok = strtok(NULL, "\t");
  p->bin = tok ? atoi(tok) : 0;
  tok = strtok(NULL, "\t");
  p->sol = tok ? atoi(tok) : 0;
  tok = strtok(NULL, "\t");
  w = tok ? atoi(tok) : 1;
  tok = strtok(NULL, "\t");
  p->nbm = tok ? atoi(tok) : 0;
  if (p->nbm > MAX_MESH) p->nbm = MAX_MESH;

  for (k = 0; k < p->nbm; k++) {
    tok = strtok(NULL, "\t");
    strncpy(p->name[k], tok ? tok : "mesh", sizeof(p->name[k]) - 1);
  }

  tok = strtok(NULL, "\t");
  if (tok && !strncmp(tok, "G=", 2)) strncpy(p->opt, tok + 2, sizeof(p->opt) - 1);

  if (!p->file[0] || p->nbm < 1) {
    close(fd);
    return (0);
  }

  srvnb++;
  srvGo( ); /* a client which waits ? no more */
  if (!srvShow(srvnb - 1)) {
    close(fd);
    return (0);
  }

  if (w)
    srvwait = fd;
  else {
    if (write(fd, "GO\n", 3) < 0) fprintf(stderr, "  ## medit server: client lost\n");

    close(fd);
  }

  srvTitle( );
  return (1);
}

static void srvTimer(int v) {
  int fd = accept(srvlisten, NULL, NULL);

  if (fd >= 0) {
    fcntl(fd, F_SETFL, fcntl(fd, F_GETFL, 0) & ~O_NONBLOCK);
    srvHandle(fd);
  }

  glutTimerFunc(50, srvTimer, 0);
}

/* keys of the server mode */
int srvKey(unsigned char key) {
  if (!ffserver) return (0);

  if (key == '\r' || key == '\n') {
    if (srvcur < srvnb - 1) {
      if (!srvShow(srvcur + 1)) exit(1);
    } else
      srvGo( ); /* the last plot: continue the script */

    srvTitle( );
    return (1);
  } else if (key == 'p') {
    if (srvcur > 0) {
      if (!srvShow(srvcur - 1)) exit(1);

      srvTitle( );
    }

    return (1);
  }

  return (0);
}

/* start the server: socket, and wait the first plot to create the window(s) */
static int srvStart( ) {
  struct sockaddr_un addr;
  struct timeval tv = {20, 0};
  fd_set rd;
  int fd;

  strncpy(srvsock, ffserver, sizeof(srvsock) - 1);
  srvlisten = socket(AF_UNIX, SOCK_STREAM, 0);
  if (srvlisten < 0) return (0);

  memset(&addr, 0, sizeof(addr));
  addr.sun_family = AF_UNIX;
  strncpy(addr.sun_path, srvsock, sizeof(addr.sun_path) - 1);
  unlink(srvsock);
  if (bind(srvlisten, (struct sockaddr *)&addr, sizeof(addr)) < 0 || listen(srvlisten, 8) < 0) {
    fprintf(stderr, "  ## medit server: unable to create the socket %s\n", srvsock);
    srvsock[0] = '\0';
    return (0);
  }

  atexit(srvCleanup);

  FD_ZERO(&rd);
  FD_SET(srvlisten, &rd);
  if (select(srvlisten + 1, &rd, NULL, NULL, &tv) <= 0) return (0);

  fd = accept(srvlisten, NULL, NULL);
  if (fd < 0 || !srvHandle(fd)) return (0);

  fcntl(srvlisten, F_SETFL, fcntl(srvlisten, F_GETFL, 0) | O_NONBLOCK);
  glutTimerFunc(50, srvTimer, 0);
  return (1);
}
#else
int srvKey(unsigned char key) { return (0); }
#endif
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
  if (i) fprintf(stdout, "  Auxiliary Buffers\t%2d\n", (int)i);

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

      cplxMesh(mesh, k); /* complex solution: a scalar part of it (key c) */
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

#ifndef WIN32
  /* popen from FreeFEM: the plot is done by the medit server (see above) */
  if (dpopen && !ffserver && option == STANDARD && !getenv("FFMEDIT_NOSERVER") && srvClient(argv[0]) == 0)
    return (0);
#endif

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
#ifndef WIN32
  if (ffserver) {
    if (!srvStart( )) exit(1);
  } else
#endif
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
