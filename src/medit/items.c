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

/*
 *  1d graph (a la gnuplot): scalar solution on a mesh of edges along the x axis (meshL of FreeFEM)
 *  => curve (x,u(x)) with a box, tics and values, instead of a colored line in 3d.
 *  G = grid, U = overlay: the solutions of the other windows in the same graph.
 */
/* number of meshes of the current plot if they are all 1d graphs (server mode, see medit.c):
   then all the curves are in the same graph, in one window */
int graph1dnb = 0;

int graph1dMesh(pMesh mesh) {
  if (!mesh->na || !mesh->nbb || mesh->nfield != 1) return (0);

  if (mesh->nt + mesh->nq + mesh->ntet + mesh->nhex) return (0);

  if (mesh->typage != 1 && mesh->typage != 2) return (0);

  return (mesh->ymax - mesh->ymin == 0.0 && mesh->zmax - mesh->zmin == 0.0 && mesh->xmax > mesh->xmin);
}

/* options of the graph, given by FreeFEM after a | in the names of the solutions:
     medit("u w | title=... | xlabel=... | ylabel=... | legend=a,b | logx | logy | loglog",ThL,u,w);
   keys: k lines/points, x y log scale, G grid, F bold, z Z zoom, arrows translate,
         A freeze the axis (for the next plots), i reset */
char graphopt[1024] = "";
static struct {
  char title[256], xlabel[128], ylabel[128], legend[MAX_SCENE][64];
  int nleg, logx, logy;
  int fixx, fixy;             /* bound of the axis given: bb=, xrange=, yrange= */
  double xr[2], yr[2];
} gopt;
static char goptsrc[1024] = "\1";
static int gstyle = 0;             /* 0 lines, 1 lines and points, 2 points */
static int gulogx = -1, gulogy = -1; /* log scale asked with the keys (-1: no) */
static int gbold = 1, gfrozen = 0;
static double gb[4];               /* bound of the axis (log10 of the values with a log scale) */

static void graphTrim(char *dst, const char *src, int n, int size) {
  while (n > 0 && (*src == ' ' || *src == '\t')) src++, n--;

  while (n > 0 && (src[n - 1] == ' ' || src[n - 1] == '\t')) n--;

  if (n > size - 1) n = size - 1;

  memcpy(dst, src, n);
  dst[n] = '\0';
}

static void graphParse( ) {
  const char *p = graphopt;
  char f[1024];

  if (!strcmp(graphopt, goptsrc)) return;

  strcpy(goptsrc, graphopt);
  memset(&gopt, 0, sizeof(gopt));

  while (*p) {
    const char *e = strchr(p, '|');
    int n = e ? (int)(e - p) : (int)strlen(p);

    graphTrim(f, p, n, sizeof(f));
    if (!strncmp(f, "title=", 6))
      graphTrim(gopt.title, f + 6, strlen(f + 6), sizeof(gopt.title));
    else if (!strncmp(f, "xlabel=", 7))
      graphTrim(gopt.xlabel, f + 7, strlen(f + 7), sizeof(gopt.xlabel));
    else if (!strncmp(f, "ylabel=", 7))
      graphTrim(gopt.ylabel, f + 7, strlen(f + 7), sizeof(gopt.ylabel));
    else if (!strncmp(f, "legend=", 7)) {
      const char *l = f + 7;

      while (*l && gopt.nleg < MAX_SCENE) {
        const char *c = strchr(l, ',');
        int nl = c ? (int)(c - l) : (int)strlen(l);

        graphTrim(gopt.legend[gopt.nleg++], l, nl, 64);
        l += nl + (c ? 1 : 0);
      }
    } else if (!strncmp(f, "complex=", 8)) {
      ; /* the complex solutions (see cplxMesh) */
    } else if (!strncmp(f, "bb=", 3)) {
      double b[4];

      if (sscanf(f + 3, "%lf ,%lf ,%lf ,%lf", b, b + 1, b + 2, b + 3) == 4) {
        gopt.xr[0] = min(b[0], b[2]);
        gopt.xr[1] = max(b[0], b[2]);
        gopt.yr[0] = min(b[1], b[3]);
        gopt.yr[1] = max(b[1], b[3]);
        gopt.fixx = gopt.xr[1] > gopt.xr[0];
        gopt.fixy = gopt.yr[1] > gopt.yr[0];
      }
    } else if (!strncmp(f, "xrange=", 7) || !strncmp(f, "yrange=", 7)) {
      double b[2];

      if (sscanf(f + 7, "%lf :%lf", b, b + 1) == 2 && b[1] != b[0]) {
        double *r = f[0] == 'x' ? gopt.xr : gopt.yr;

        r[0] = min(b[0], b[1]);
        r[1] = max(b[0], b[1]);
        if (f[0] == 'x')
          gopt.fixx = 1;
        else
          gopt.fixy = 1;
      }
    } else if (!strcmp(f, "logx"))
      gopt.logx = 1;
    else if (!strcmp(f, "logy"))
      gopt.logy = 1;
    else if (!strcmp(f, "loglog"))
      gopt.logx = gopt.logy = 1;
    else if (f[0] && !gopt.title[0])
      graphTrim(gopt.title, f, strlen(f), sizeof(gopt.title)); /* a field without keyword = title */

    p += n + (e ? 1 : 0);
  }
}

static int graphTextWidth(const char *s) {
  int w = 0;

  for (; *s; s++) w += glutBitmapWidth(GLUT_BITMAP_HELVETICA_12, *s);

  return (w);
}

static void graphText(double x, double y, const char *s, int align) {
  const char *c;
  int b;

  x -= graphTextWidth(s) * align / 2.0;

  /* bold: the text is plotted twice */
  for (b = 0; b <= gbold; b++) {
    glRasterPos2d(floor(x) + b + 0.5, floor(y) + 0.5);

    for (c = s; *c; c++) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, *c);
  }
}

/* the tics k*step in [vmin,vmax], step = 1, 2 or 5 * 10^n; log scale: vmin, vmax = log10 */
static int graphTics(double vmin, double vmax, int nbmax, int lg, double *tics, char lab[][24],
                     double *minors, int *nminor) {
  double d, p, r, step, k, k0, k1, vm;
  int n = 0, prec, i;

  *nminor = 0;
  if (!(vmax > vmin)) return (0);

  if (nbmax < 1) nbmax = 1;

  if (lg && floor(vmax + 1e-9) - ceil(vmin - 1e-9) >= 1) {
    /* a tic by decade (or by step decades), minor tics 2..9 if step = 1 */
    int i0 = (int)ceil(vmin - 1e-9), i1 = (int)floor(vmax + 1e-9), m;
    int st = (i1 - i0 + nbmax) / nbmax;

    if (st < 1) st = 1;

    for (i = i0; i <= i1 && n < 64; i += st) {
      tics[n] = i;
      snprintf(lab[n++], 24, "%g", pow(10.0, i));
    }

    if (st == 1)
      for (i = i0 - 1; i <= i1; i++)
        for (m = 2; m < 10; m++) {
          double t = i + log10((double)m);

          if (t >= vmin && t <= vmax && *nminor < 256) minors[(*nminor)++] = t;
        }

    return (n);
  }

  d = (vmax - vmin) / nbmax;
  p = pow(10.0, floor(log10(d)));
  r = d / p;
  step = p * (r > 5 ? 10 : (r > 2 ? 5 : (r > 1 ? 2 : 1)));
  k0 = ceil(vmin / step - 1e-9);
  k1 = floor(vmax / step + 1e-9);
  vm = max(fabs(vmin), fabs(vmax));
  prec = (int)ceil(log10(vm / step)) + 1;
  if (prec < 2) prec = 2;
  if (prec > 15) prec = 15;

  for (k = k0; k <= k1 && n < 64; k++) {
    tics[n] = k * step + 0.0;
    if (lg)
      snprintf(lab[n], 24, "%.3g", pow(10.0, tics[n]));
    else
      snprintf(lab[n], 24, "%.*g", prec, tics[n]);

    n++;
  }

  return (n);
}

/* keys of the graph: return 1 if the key is used */
int graphKey(pScene sc, unsigned char key) {
  double c, h;

  switch (key) {
    case 'k':
      gstyle = (gstyle + 1) % 3;
      break;
    case 'x':
      gulogx = !(gulogx >= 0 ? gulogx : gopt.logx);
      gfrozen = 0;
      break;
    case 'y':
      gulogy = !(gulogy >= 0 ? gulogy : gopt.logy);
      gfrozen = 0;
      break;
    case 'F':
      gbold = !gbold;
      break;
    case 'A':
      gfrozen = !gfrozen;
      fprintf(stdout, "  graph: axis %s\n", gfrozen ? "frozen" : "automatic");
      break;
    case 'i':
      gfrozen = 0;
      gulogx = gulogy = -1;
      break;
    case 'z':
    case 'Z':
      h = key == 'z' ? 0.8 : 1.25;
      c = 0.5 * (gb[0] + gb[1]);
      gb[0] = c + h * (gb[0] - c);
      gb[1] = c + h * (gb[1] - c);
      c = 0.5 * (gb[2] + gb[3]);
      gb[2] = c + h * (gb[2] - c);
      gb[3] = c + h * (gb[3] - c);
      gfrozen = 1;
      break;
    default:
      return (0);
  }

  glutPostRedisplay( );
  return (1);
}

/* arrows: translation of 10 % */
int graphSpecial(pScene sc, int key) {
  double dx = 0.0, dy = 0.0;

  switch (key) {
    case GLUT_KEY_LEFT:
      dx = -0.1;
      break;
    case GLUT_KEY_RIGHT:
      dx = 0.1;
      break;
    case GLUT_KEY_DOWN:
      dy = -0.1;
      break;
    case GLUT_KEY_UP:
      dy = 0.1;
      break;
    default:
      return (0);
  }

  dx *= gb[1] - gb[0];
  dy *= gb[3] - gb[2];
  gb[0] += dx;
  gb[1] += dx;
  gb[2] += dy;
  gb[3] += dy;
  gfrozen = 1;
  glutPostRedisplay( );
  return (1);
}

#define GLOGV(v, lg) ((lg) ? log10(v) : (v))

void drawGraph1d(pScene sc) {
  static const float coldark[6][3] = {{0.35, 0.75, 1.0}, {1.0, 0.65, 0.1}, {0.3, 0.9, 0.5},
                                      {1.0, 0.4, 0.4},   {0.8, 0.6, 1.0},  {1.0, 1.0, 0.4}};
  static const float collight[6][3] = {{0.58, 0.0, 0.83}, {0.0, 0.62, 0.45}, {0.34, 0.71, 0.91},
                                       {0.9, 0.62, 0.0},  {0.0, 0.45, 0.7},  {0.9, 0.12, 0.06}};
  static double tx[64], ty[64], mx[256], my[256];
  static char lx[64][24], ly[64][24];
  pMesh mesh = cv.mesh[sc->idmesh], ms[MAX_SCENE];
  double xmin, xmax, ymin, ymax, d, ax, ay;
  int x0, y0, x1, y1, nx, ny, nmx, nmy, i, j, k, nc = 0, dark, logx, logy, pass;
  unsigned char *ypix = NULL;
  int yw = 0, yh = 0;

  graphParse( );
  logx = gulogx >= 0 ? gulogx : gopt.logx;
  logy = gulogy >= 0 ? gulogy : gopt.logy;

  /* the curves: all the meshes of the plot, or this scene (+ the other scenes with the overlay) */
  if (graph1dnb > 1) {
    for (j = 0; j < graph1dnb && nc < MAX_SCENE; j++)
      if (cv.mesh[j] && graph1dMesh(cv.mesh[j])) ms[nc++] = cv.mesh[j];
  } else
    ms[nc++] = mesh;

  if (sc->overlay && graph1dnb <= 1)
    for (j = 0; j < cv.nbs && nc < MAX_SCENE; j++) {
      pScene scj = cv.scene[j];

      if (scj && scj != sc && graph1dMesh(cv.mesh[scj->idmesh])) ms[nc++] = cv.mesh[scj->idmesh];
    }

  /* bound of the axis: 5 % of margin */
  if (!gfrozen) {
    xmin = ymin = 1e300;
    xmax = ymax = -1e300;

    for (j = 0; j < nc; j++) {
      pMesh m = ms[j];

      for (k = 1; k <= m->na; k++) {
        pEdge pe = &m->edge[k];

        if (!pe->v[0]) continue;

        for (i = 0; i < 2; i++) {
          double x = m->point[pe->v[i]].c[0] + m->xtra;
          double y = m->typage == 1 ? (k <= m->nbb ? m->sol[k].bb : 0.0) : m->sol[pe->v[i]].bb;

          if ((logx && x <= 0.0) || (logy && y <= 0.0)) continue;

          x = GLOGV(x, logx);
          y = GLOGV(y, logy);
          xmin = min(xmin, x);
          xmax = max(xmax, x);
          ymin = min(ymin, y);
          ymax = max(ymax, y);
        }
      }
    }

    if (xmin > xmax) xmin = 0.0, xmax = 1.0;

    if (ymin > ymax) ymin = 0.0, ymax = 1.0;

    if (!(xmax > xmin)) {
      d = xmin != 0.0 ? 0.05 * fabs(xmin) : 0.5;
      xmin -= d;
      xmax += d;
    }

    if (!(ymax > ymin)) {
      d = ymin != 0.0 ? 0.05 * fabs(ymin) : 0.5;
      ymin -= d;
      ymax += d;
    }

    d = 0.05 * (xmax - xmin);
    gb[0] = xmin - d;
    gb[1] = xmax + d;
    d = 0.05 * (ymax - ymin);
    gb[2] = ymin - d;
    gb[3] = ymax + d;

    /* the bound given in the options (bb=, xrange=, yrange=), without margin */
    if (gopt.fixx && (!logx || gopt.xr[0] > 0.0)) {
      gb[0] = GLOGV(gopt.xr[0], logx);
      gb[1] = GLOGV(gopt.xr[1], logx);
    }

    if (gopt.fixy && (!logy || gopt.yr[0] > 0.0)) {
      gb[2] = GLOGV(gopt.yr[0], logy);
      gb[3] = GLOGV(gopt.yr[1], logy);
    }
  }

  xmin = gb[0];
  xmax = gb[1];
  ymin = gb[2];
  ymax = gb[3];

  /* the box of the axis in pixel */
  x0 = 85 + (gopt.ylabel[0] ? 22 : 0);
  y0 = 45 + (gopt.xlabel[0] ? 22 : 0);
  x1 = max(sc->par.xs - 30, x0 + 10);
  y1 = max(sc->par.ys - 40, y0 + 10);
  ax = (x1 - x0) / (xmax - xmin);
  ay = (y1 - y0) / (ymax - ymin);
  nx = graphTics(xmin, xmax, max(2, (x1 - x0) / 80), logx, tx, lx, mx, &nmx);
  ny = graphTics(ymin, ymax, max(2, (y1 - y0) / 50), logy, ty, ly, my, &nmy);
  dark = sc->par.back[0] + sc->par.back[1] + sc->par.back[2] < 1.5;

  glMatrixMode(GL_PROJECTION);
  glPushMatrix( );
  glLoadIdentity( );
  gluOrtho2D(0., sc->par.xs, 0., sc->par.ys);
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix( );
  glLoadIdentity( );
  glDisable(GL_LIGHTING);
  glDisable(GL_DEPTH_TEST);
  glDisable(GL_CLIP_PLANE0);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

  /* the label of the y axis is vertical: the bitmap fonts can not be rotated, so the text is
     plotted on the bottom left corner (empty), read, rotated of 90 degree, and the corner cleaned */
  if (gopt.ylabel[0]) {
    unsigned char *a;
    int w = graphTextWidth(gopt.ylabel) + 3, h = 9 + 5 + 3;

    a = (unsigned char *)malloc(3 * w * h);
    ypix = (unsigned char *)malloc(3 * w * h);
    if (a && ypix) {
      glColor4fv(sc->par.line);
      graphText(1, 5, gopt.ylabel, 0);
      glPixelStorei(GL_PACK_ALIGNMENT, 1);
      glReadBuffer(GL_BACK);
      glReadPixels(0, 0, w, h, GL_RGB, GL_UNSIGNED_BYTE, a);
      glColor4fv(sc->par.back);
      glBegin(GL_QUADS);
      glVertex2d(0, 0);
      glVertex2d(w, 0);
      glVertex2d(w, h);
      glVertex2d(0, h);
      glEnd( );

      for (j = 0; j < h; j++)
        for (i = 0; i < w; i++)
          for (k = 0; k < 3; k++) ypix[3 * (i * h + (h - 1 - j)) + k] = a[3 * (j * w + i) + k];

      yw = w;
      yh = h;
    }

    free(a);
  }

  /* grid */
  if (sc->item & S_GRID) {
    glLineWidth(1.0);
    glColor3f(0.5, 0.5, 0.5);
    glEnable(GL_LINE_STIPPLE);
    glLineStipple(1, 0x3333);
    glBegin(GL_LINES);

    for (i = 0; i < nx; i++) {
      glVertex2d(x0 + (tx[i] - xmin) * ax, y0);
      glVertex2d(x0 + (tx[i] - xmin) * ax, y1);
    }

    for (i = 0; i < ny; i++) {
      glVertex2d(x0, y0 + (ty[i] - ymin) * ay);
      glVertex2d(x1, y0 + (ty[i] - ymin) * ay);
    }

    glEnd( );
    glDisable(GL_LINE_STIPPLE);
  }

  /* curves, clipped in the box: pass 0 lines, pass 1 points */
  glEnable(GL_SCISSOR_TEST);
  glScissor(x0, y0, x1 - x0, y1 - y0);
  glLineWidth(2.0);
  glPointSize(5.0);

  for (pass = 0; pass < 2; pass++) {
    if (pass == 0 && gstyle == 2) continue;

    if (pass == 1 && gstyle == 0) continue;

    for (j = 0; j < nc; j++) {
      pMesh m = ms[j];

      glColor3fv(dark ? coldark[j % 6] : collight[j % 6]);
      glBegin(pass ? GL_POINTS : GL_LINES);

      for (k = 1; k <= m->na; k++) {
        pEdge pe = &m->edge[k];
        double xa, xb, ya, yb;

        if (!pe->v[0]) continue;

        /* the points are translated at the center of the mesh (see meshBox) */
        xa = m->point[pe->v[0]].c[0] + m->xtra;
        xb = m->point[pe->v[1]].c[0] + m->xtra;
        if (m->typage == 1) {
          if (k > m->nbb) continue;

          ya = yb = m->sol[k].bb;
        } else {
          ya = m->sol[pe->v[0]].bb;
          yb = m->sol[pe->v[1]].bb;
        }

        if (pass && m->typage == 1) {
          /* solution by element: one point at the middle */
          xa = xb = 0.5 * (xa + xb);
        }

        if ((logx && (xa <= 0.0 || xb <= 0.0)) || (logy && (ya <= 0.0 || yb <= 0.0))) continue;

        glVertex2d(x0 + (GLOGV(xa, logx) - xmin) * ax, y0 + (GLOGV(ya, logy) - ymin) * ay);
        glVertex2d(x0 + (GLOGV(xb, logx) - xmin) * ax, y0 + (GLOGV(yb, logy) - ymin) * ay);
      }

      glEnd( );
    }
  }

  glPointSize(1.0);
  glDisable(GL_SCISSOR_TEST);

  /* box, tics, values */
  glLineWidth(1.0);
  glColor4fv(sc->par.line);
  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  glBegin(GL_LINE_LOOP);
  glVertex2d(x0, y0);
  glVertex2d(x1, y0);
  glVertex2d(x1, y1);
  glVertex2d(x0, y1);
  glEnd( );
  glBegin(GL_LINES);

  for (pass = 0; pass < 2; pass++) {
    /* pass 0 tics, pass 1 minor tics */
    int lt = pass ? 3 : 6, n1 = pass ? nmx : nx, n2 = pass ? nmy : ny;
    double *v1 = pass ? mx : tx, *v2 = pass ? my : ty;

    for (i = 0; i < n1; i++) {
      double x = x0 + (v1[i] - xmin) * ax;

      glVertex2d(x, y0);
      glVertex2d(x, y0 + lt);
      glVertex2d(x, y1);
      glVertex2d(x, y1 - lt);
    }

    for (i = 0; i < n2; i++) {
      double y = y0 + (v2[i] - ymin) * ay;

      glVertex2d(x0, y);
      glVertex2d(x0 + lt, y);
      glVertex2d(x1, y);
      glVertex2d(x1 - lt, y);
    }
  }

  glEnd( );

  for (i = 0; i < nx; i++) graphText(x0 + (tx[i] - xmin) * ax, y0 - 17, lx[i], 1);

  for (i = 0; i < ny; i++) graphText(x0 - 6, y0 + (ty[i] - ymin) * ay - 4, ly[i], 2);

  /* title, labels of the axis */
  if (gopt.title[0])
    graphText(0.5 * (x0 + x1), y1 + 12, gopt.title, 1);
  else if (nc == 1 && !gopt.nleg)
    graphText(0.5 * (x0 + x1), y1 + 12, mesh->name, 1);

  if (gopt.xlabel[0]) graphText(0.5 * (x0 + x1), y0 - 36, gopt.xlabel, 1);

  if (ypix && yw) {
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glRasterPos2d(6.5, floor(0.5 * (y0 + y1 - yw)) + 0.5);
    glDrawPixels(yh, yw, GL_RGB, GL_UNSIGNED_BYTE, ypix);
  }

  free(ypix);

  /* the key: name and sample of each curve */
  if (nc > 1 || gopt.nleg)
    for (j = 0; j < nc; j++) {
      double y = y1 - 19 - 17 * j;

      glColor4fv(sc->par.line);
      graphText(x1 - 58, y, j < gopt.nleg ? gopt.legend[j] : ms[j]->name, 2);
      glColor3fv(dark ? coldark[j % 6] : collight[j % 6]);
      if (gstyle != 2) {
        glLineWidth(2.0);
        glBegin(GL_LINES);
        glVertex2d(x1 - 50, y + 4);
        glVertex2d(x1 - 10, y + 4);
        glEnd( );
        glLineWidth(1.0);
      }

      if (gstyle != 0) {
        glPointSize(5.0);
        glBegin(GL_POINTS);
        glVertex2d(x1 - 30, y + 4);
        glEnd( );
        glPointSize(1.0);
      }
    }

  glMatrixMode(GL_PROJECTION);
  glPopMatrix( );
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix( );
  glEnable(GL_DEPTH_TEST);
}

void drawRulers(pScene sc) {
  pMesh mesh = cv.mesh[sc->idmesh];

  if (ddebug) printf("draw rulers\n");

  glPushMatrix( );
  glTranslatef(-mesh->xtra, -mesh->ytra, -mesh->ztra);
  glColor3fv(sc->par.line);
  glEnable(GL_LINE_STIPPLE);
  glLineStipple(1, 0x0444);
  glLineWidth(2.);
  glBegin(GL_LINES);
  glVertex3f(mesh->xmin, mesh->ymin, mesh->zmin);
  glVertex3f(mesh->xmax, mesh->ymin, mesh->zmin);
  glVertex3f(mesh->xmin, mesh->ymax, mesh->zmin);
  glVertex3f(mesh->xmax, mesh->ymax, mesh->zmin);
  glVertex3f(mesh->xmin, mesh->ymin, mesh->zmin);
  glVertex3f(mesh->xmin, mesh->ymax, mesh->zmin);
  glVertex3f(mesh->xmax, mesh->ymin, mesh->zmin);
  glVertex3f(mesh->xmax, mesh->ymax, mesh->zmin);

  glVertex3f(mesh->xmin, mesh->ymin, mesh->zmax);
  glVertex3f(mesh->xmax, mesh->ymin, mesh->zmax);
  glVertex3f(mesh->xmin, mesh->ymax, mesh->zmax);
  glVertex3f(mesh->xmax, mesh->ymax, mesh->zmax);
  glVertex3f(mesh->xmin, mesh->ymin, mesh->zmax);
  glVertex3f(mesh->xmin, mesh->ymax, mesh->zmax);
  glVertex3f(mesh->xmax, mesh->ymin, mesh->zmax);
  glVertex3f(mesh->xmax, mesh->ymax, mesh->zmax);

  glVertex3f(mesh->xmin, mesh->ymin, mesh->zmin);
  glVertex3f(mesh->xmin, mesh->ymin, mesh->zmax);
  glVertex3f(mesh->xmin, mesh->ymax, mesh->zmin);
  glVertex3f(mesh->xmin, mesh->ymax, mesh->zmax);
  glVertex3f(mesh->xmax, mesh->ymin, mesh->zmin);
  glVertex3f(mesh->xmax, mesh->ymin, mesh->zmax);
  glVertex3f(mesh->xmax, mesh->ymax, mesh->zmin);
  glVertex3f(mesh->xmax, mesh->ymax, mesh->zmax);
  glEnd( );
  glLineWidth(1.);
  glDisable(GL_LINE_STIPPLE);
  glPopMatrix( );
}

void drawAxis(pScene sc, int dim) {
  pMesh mesh;

  /* default */
  if (ddebug) printf("draw axis\n");

  mesh = cv.mesh[sc->idmesh];

  glPushMatrix( );
  glTranslatef(1.01 * (mesh->xmin - mesh->xtra), 1.01 * (mesh->ymin - mesh->ytra),
               1.01 * (mesh->zmin - mesh->ztra));
  glScalef(0.6 * sc->dmin, 0.6 * sc->dmin, 0.6 * sc->dmin);
  glLineWidth(max(2, sc->par.linewidth));
  glColor3f(1.0, 0., 0.);
  if (mesh->dim == 2) {
    glBegin(GL_LINE_STRIP);
    glVertex2f(0.0, 0.0);
    glVertex2f(1.0, 0.0);
    glVertex2f(0.95, 0.01);
    glVertex2f(0.95, -0.01);
    glVertex2f(1.0, 0.0);
    glEnd( );
    glBegin(GL_LINE_STRIP);
    glVertex2f(0.0, 0.0);
    glVertex2f(0.0, 1.0);
    glVertex2f(-0.01, 0.95);
    glVertex2f(0.01, 0.95);
    glVertex2f(0.0, 1.0);
    glEnd( );

    glColor3f(1. - sc->par.back[0], 1.0 - sc->par.back[1], 1.0 - sc->par.back[2]);
    glRasterPos3f(1.02, 0.0, 0.0);
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, 'x');
    glRasterPos3f(0.0, 1.02, 0.0);
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, 'y');
  } else {
    glBegin(GL_LINE_STRIP);
    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(1.0, 0.0, 0.0);
    glVertex3f(0.95, 0.01, 0.0);
    glVertex3f(0.95, -0.01, 0.0);
    glVertex3f(1.0, 0.0, 0.0);
    glVertex3f(0.95, 0.0, 0.01);
    glVertex3f(0.95, 0.0, -0.01);
    glVertex3f(1.0, 0.0, 0.0);
    glEnd( );
    glBegin(GL_LINE_STRIP);
    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(0.0, 1.0, 0.0);
    glVertex3f(0.0, 0.95, 0.01);
    glVertex3f(0.0, 0.95, -0.01);
    glVertex3f(0.0, 1.0, 0.0);
    glVertex3f(0.01, 0.95, 0.0);
    glVertex3f(-0.01, 0.95, 0.0);
    glVertex3f(0.0, 1.0, 0.0);
    glEnd( );
    glBegin(GL_LINE_STRIP);
    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(0.0, 0.0, 1.0);
    glVertex3f(0.01, 0.0, 0.95);
    glVertex3f(-0.01, 0.0, 0.95);
    glVertex3f(0.0, 0.0, 1.0);
    glVertex3f(0.0, 0.01, 0.95);
    glVertex3f(0.0, -0.01, 0.95);
    glVertex3f(0.0, 0.0, 1.0);
    glEnd( );

    glColor3f(1. - sc->par.back[0], 1.0 - sc->par.back[1], 1.0 - sc->par.back[2]);
    glRasterPos3f(1.02, 0.0, 0.0);
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, 'x');
    glRasterPos3f(0.0, 1.02, 0.0);
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, 'y');
    glRasterPos3f(0.0, 0.0, 1.02);
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, 'z');
  }

  glLineWidth(1.0);
  glPopMatrix( );
}

/* complex solution of FreeFEM (option complex=mask of -graph): the solution k is complex if
   mask[k] = 1, sent as the vector (real part, imaginary part); medit shows a scalar part of it:
   cplxpart = 0 real part, 1 imaginary part, 2 modulus, 3 argument (key c) */
int cplxpart = 0;
const char *cplxname[4] = {"real part", "imaginary part", "modulus", "argument"};

/* the mesh k is complex: according to the options of the plot */
void cplxMesh(pMesh mesh, int k) {
  const char *c = strstr(graphopt, "complex=");

  mesh->cplx = 0;
  if (!c || !mesh->nbb || !mesh->sol) return;

  c += 8;
  if ((int)strlen(c) > k && c[k] == '1' && mesh->nfield == mesh->dim) mesh->cplx = 1;

  cplxSet(mesh);
}

/* the scalar value bb of the solution = the part cplxpart of the complex value */
void cplxSet(pMesh mesh) {
  int k;

  if (!mesh->cplx) return;

  mesh->nfield = 1;
  mesh->bbmin = 1e30;
  mesh->bbmax = -1e30;
  for (k = 1; k <= mesh->nbb; k++) {
    pSolution ps = &mesh->sol[k];
    double re = ps->m[0], im = ps->m[1];

    ps->bb = cplxpart == 0 ? re : (cplxpart == 1 ? im : (cplxpart == 2 ? sqrt(re * re + im * im) : atan2(im, re)));
    if (ps->bb < mesh->bbmin) mesh->bbmin = ps->bb;

    if (ps->bb > mesh->bbmax) mesh->bbmax = ps->bb;
  }
}

/* graduated axes (key T): x, y (and z in 3d or for the 2d elevation) on the edges of the
   bounding box, with the tics and the values, drawn in pixel on the screen */
static GLdouble amv[16], apj[16];
static GLint avp[4];
static double abox[4 * 128]; /* boxes of the texts already drawn: no text over another one */
static int nabox = 0;

static void axisProject(const double *P, double *S) {
  GLdouble x, y, z;

  gluProject(P[0], P[1], P[2], amv, apj, avp, &x, &y, &z);
  S[0] = x - avp[0];
  S[1] = y - avp[1];
}

static int axisText(double x, double y, const char *t, int align) {
  double w = graphTextWidth(t), x0, x1, y0, y1;
  int i;

  x0 = x - w * align / 2.0 - 2;
  x1 = x0 + w + 4;
  y0 = y - 2;
  y1 = y + 9 + 2;
  for (i = 0; i < nabox; i++)
    if (x0 < abox[4 * i + 2] && abox[4 * i] < x1 && y0 < abox[4 * i + 3] && abox[4 * i + 1] < y1)
      return (0);

  if (nabox < 128) {
    abox[4 * nabox] = x0;
    abox[4 * nabox + 1] = y0;
    abox[4 * nabox + 2] = x1;
    abox[4 * nabox + 3] = y1;
    nabox++;
  }

  graphText(x, y, t, align);
  return (1);
}

/* axis from A to B (display coordinates) for the values vmin to vmax: tics and values
   on the side opposite to C (center of the box on the screen), shifted of gap pixels */
static void axisTics(double *A, double *B, double vmin, double vmax, double *C, double gap,
                     const char *name) {
  double SA[2], SB[2], d[2], o[2], M[2], P[2], Q[2], l, ext = 0.0, half, step;
  double tics[64], minors[256], ends[2];
  char lab[64][24], buf[32];
  int n, nminor, i, e, align, near;

  if (!(vmax > vmin)) return;

  axisProject(A, SA);
  axisProject(B, SB);
  d[0] = SB[0] - SA[0];
  d[1] = SB[1] - SA[1];
  l = sqrt(d[0] * d[0] + d[1] * d[1]);
  if (l < 10.0) return; /* the axis is seen end-on */

  d[0] /= l;
  d[1] /= l;
  o[0] = d[1];
  o[1] = -d[0];
  M[0] = 0.5 * (SA[0] + SB[0]);
  M[1] = 0.5 * (SA[1] + SB[1]);
  if (o[0] * (M[0] - C[0]) + o[1] * (M[1] - C[1]) < 0.0) {
    o[0] = -o[0];
    o[1] = -o[1];
  }

  for (i = 0; i < 2; i++) {
    SA[i] += o[i] * gap;
    SB[i] += o[i] * gap;
    M[i] += o[i] * gap;
  }

  glBegin(GL_LINES);
  glVertex2d(SA[0], SA[1]);
  glVertex2d(SB[0], SB[1]);
  glEnd( );
  n = graphTics(vmin, vmax, max(2, (int)(l / 70.0)), 0, tics, lab, minors, &nminor);
  align = o[0] > 0.5 ? 0 : (o[0] < -0.5 ? 2 : 1);
  step = n > 1 ? tics[1] - tics[0] : vmax - vmin;
  ends[0] = vmin;
  ends[1] = vmax;

  /* the tics, then the values at the ends of the axis if no tic is near */
  for (i = 0; i < n + 2; i++) {
    const char *t;
    double a;

    if (i < n) {
      a = (tics[i] - vmin) / (vmax - vmin);
      t = lab[i];
    } else {
      e = i - n;
      near = 0;
      for (int j = 0; j < n; j++)
        if (fabs(tics[j] - ends[e]) < 0.25 * step) near = 1;

      if (near) continue;

      a = e;
      snprintf(buf, sizeof(buf), "%.3g", ends[e]);
      t = buf;
    }

    P[0] = SA[0] + (SB[0] - SA[0]) * a;
    P[1] = SA[1] + (SB[1] - SA[1]) * a;
    if (P[0] < 0 || P[1] < 0 || P[0] > avp[2] || P[1] > avp[3]) continue;

    glBegin(GL_LINES);
    glVertex2d(P[0], P[1]);
    glVertex2d(P[0] + 5 * o[0], P[1] + 5 * o[1]);
    glEnd( );
    Q[0] = P[0] + 9 * o[0];
    Q[1] = P[1] + 9 * o[1];
    if (axisText(Q[0], Q[1] - 4.5 + o[1] * 9 * 0.7, t, align))
      ext = max(ext, fabs(o[0]) * graphTextWidth(t) + fabs(o[1]) * 1.2 * 9);
  }

  /* the name of the axis at the middle, after the values */
  half = 0.5 * (fabs(o[0]) * graphTextWidth(name) + fabs(o[1]) * 9);
  Q[0] = M[0] + o[0] * (9 + ext + 6 + half);
  Q[1] = M[1] + o[1] * (9 + ext + 6 + half);
  axisText(Q[0], Q[1] - 4.5, name, 1);
}

void drawAxesTics(pScene sc, pMesh mesh) {
  double A[3], B[3], V0[3], V1[3], C[3], SC[2], E0[3], E1[3];
  const char *name[3] = {"x", "y", "z"};
  int c, i, j, naxis, elev;

  if (graph1dMesh(mesh)) return;

  /* the box in the display coordinates (the mesh is centered on 0) and the values on the axes */
  V0[0] = mesh->xmin;
  V1[0] = mesh->xmax;
  V0[1] = mesh->ymin;
  V1[1] = mesh->ymax;
  A[0] = mesh->xmin - mesh->xtra;
  B[0] = mesh->xmax - mesh->xtra;
  A[1] = mesh->ymin - mesh->ytra;
  B[1] = mesh->ymax - mesh->ytra;
  elev = mesh->dim == 2 && (sc->mode & S_ALTITUDE) && mesh->nbb;
  if (mesh->dim == 3) {
    V0[2] = mesh->zmin;
    V1[2] = mesh->zmax;
    A[2] = mesh->zmin - mesh->ztra;
    B[2] = mesh->zmax - mesh->ztra;
    naxis = 3;
  } else if (elev) {
    /* 2d elevation: z = altcoef * u - 0.25 * ztra (see alt2dList) */
    extern GLfloat altcoef;

    V0[2] = mesh->bbmin;
    V1[2] = mesh->bbmax;
    A[2] = altcoef * mesh->bbmin - 0.25 * mesh->ztra;
    B[2] = altcoef * mesh->bbmax - 0.25 * mesh->ztra;
    naxis = 3;
  } else {
    V0[2] = V1[2] = A[2] = B[2] = 0.0;
    naxis = 2;
  }

  glGetDoublev(GL_MODELVIEW_MATRIX, amv);
  glGetDoublev(GL_PROJECTION_MATRIX, apj);
  glGetIntegerv(GL_VIEWPORT, avp);
  for (i = 0; i < 3; i++) C[i] = 0.5 * (A[i] + B[i]);

  axisProject(C, SC);

  /* the edges of the box: for x and y the lowest on the screen, for z the leftmost
     (outer edges, never behind the mesh; the values out of the window are not drawn) */
  {
    double ex0[3][3], ex1[3][3];

    for (c = 0; c < naxis; c++) {
      int c1 = (c + 1) % 3, c2 = (c + 2) % 3;
      double best = 1e100;

      for (i = 0; i < 2; i++)
        for (j = 0; j < 2; j++) {
          double Q[3], S[2], v;

          E0[c] = A[c];
          E1[c] = B[c];
          E0[c1] = E1[c1] = i ? B[c1] : A[c1];
          E0[c2] = E1[c2] = j ? B[c2] : A[c2];
          Q[0] = 0.5 * (E0[0] + E1[0]);
          Q[1] = 0.5 * (E0[1] + E1[1]);
          Q[2] = 0.5 * (E0[2] + E1[2]);
          axisProject(Q, S);
          v = c < 2 ? S[1] : S[0];
          if (v < best) {
            best = v;
            memcpy(ex0[c], E0, sizeof(E0));
            memcpy(ex1[c], E1, sizeof(E1));
          }
        }
    }

    /* drawn in pixel */
    glMatrixMode(GL_PROJECTION);
    glPushMatrix( );
    glLoadIdentity( );
    gluOrtho2D(0., avp[2], 0., avp[3]);
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix( );
    glLoadIdentity( );
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_LIGHTING);
    glColor3f(1. - sc->par.back[0], 1.0 - sc->par.back[1], 1.0 - sc->par.back[2]);
    glLineWidth(1.0);
    nabox = 0;
    for (c = 0; c < naxis; c++)
      axisTics(ex0[c], ex1[c], V0[c], V1[c], SC, naxis == 2 ? 6.0 : 0.0, name[c]);

    glMatrixMode(GL_PROJECTION);
    glPopMatrix( );
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix( );
    glEnable(GL_DEPTH_TEST);
  }
}

void drawBox(pScene sc, pMesh mesh, int mode) {
  pMaterial pm;

  glDisable(GL_LIGHTING);
  glPushMatrix( );
  glScalef(1.01 * fabs(mesh->xmax - mesh->xmin), 1.01 * fabs(mesh->ymax - mesh->ymin),
           1.01 * fabs(mesh->zmax - mesh->zmin));
  glColor3f(1.0, 0.0, 0.5);
  glutWireCube(1.0);
  glPopMatrix( );

  /* one box per sub-domain */
  if (mode) {
    int m;

    for (m = 0; m < sc->par.nbmat; m++) {
      int i;

      pm = &sc->material[m];

      for (i = 0; i < MAX_LIST; i++) {
        float cx, cy, cz;
        int k;

        k = pm->depmat[i];
        if (!k || pm->flag) continue;

        cx = 0.5 * (pm->ext[3] + pm->ext[0]);
        cy = 0.5 * (pm->ext[4] + pm->ext[1]);
        cz = 0.5 * (pm->ext[5] + pm->ext[2]);
        glPushMatrix( );
        glColor3fv(pm->dif);
        glTranslatef(cx, cy, cz);
        glScalef(pm->ext[3] - pm->ext[0], pm->ext[4] - pm->ext[1], pm->ext[5] - pm->ext[2]);
        glutWireCube(1.0);
        glPopMatrix( );
      }
    }
  }
}

void drawCube(pScene sc, pMesh mesh) {
  pTransform cubetr;
  pCube cube;
  float x1, y1, z1, z2, xd, yd;

  cube = sc->cube;
  cubetr = cube->cubetr;

  if (cube->active & C_UPDATE) updateCube(cube, mesh);

  glDisable(GL_LIGHTING);
  glPushMatrix( );
  glMultMatrixf(cubetr->matrix);

  glLineWidth(3.0);
  if (cube->active & C_EDIT)
    glColor3f(1.0, 0.0, 1.0);
  else if (cube->active & C_FREEZE)
    glColor3f(0.0, 0.6, 0.9);
  else
    glColor3f(0.0, 1.0, 0.0);

  x1 = cube->cmi[0] - mesh->xtra;
  y1 = cube->cmi[1] - mesh->ytra;
  z1 = cube->cmi[2] - mesh->ztra;
  z2 = cube->cma[2] - mesh->ztra;
  xd = cube->cma[0] - cube->cmi[0];
  yd = cube->cma[1] - cube->cmi[1];

  glBegin(GL_QUADS);
  glVertex3f(x1, y1, z1);
  glVertex3f(x1 + xd, y1, z1);
  glVertex3f(x1 + xd, y1 + yd, z1);
  glVertex3f(x1, y1 + yd, z1);

  glVertex3f(x1, y1, z2);
  glVertex3f(x1 + xd, y1, z2);
  glVertex3f(x1 + xd, y1 + yd, z2);
  glVertex3f(x1, y1 + yd, z2);
  glEnd( );
  glBegin(GL_LINES);
  glVertex3f(x1, y1, z1);
  glVertex3f(x1, y1, z2);

  glVertex3f(x1 + xd, y1, z1);
  glVertex3f(x1 + xd, y1, z2);

  glVertex3f(x1 + xd, y1 + yd, z1);
  glVertex3f(x1 + xd, y1 + yd, z2);

  glVertex3f(x1, y1 + yd, z1);
  glVertex3f(x1, y1 + yd, z2);
  glEnd( );

  glLineWidth(1.0);

  glPopMatrix( );
}

void drawGrid(pScene sc, pMesh mesh) {
  /* default */
  if (ddebug) printf("draw grid + graduation\n");

  if (!sc->grid) {
    int k;

    sc->grid = glGenLists(1);
    glNewList(sc->grid, GL_COMPILE);
    glBegin(GL_LINES);

    for (k = 0; k < 5; k++) {
      glVertex3f(k * 0.25, 0., 0.);
      glVertex3f(k * 0.25, 1., 0.);
      glVertex3f(0., k * 0.25, 0.);
      glVertex3f(1., k * 0.25, 0.);
      glVertex3f(0., k * 0.25, 0.);
      glVertex3f(0., k * 0.25, 1.);
      glVertex3f(0., 0., k * 0.25);
      glVertex3f(0., 1., k * 0.25);
      glVertex3f(k * 0.25, 0., 0.);
      glVertex3f(k * 0.25, 0., 1.);
      glVertex3f(0., 0., k * 0.25);
      glVertex3f(1., 0., k * 0.25);
    }

    glEnd( );
    glEndList( );
  }

  /* call display list */
  glPushMatrix( );
  glTranslatef(-0.3 * fabs(sc->dmax), 0., -4.7 * sc->dmax);
  glRotatef(-60., 1., 0., 0.);
  glRotatef(-120., 0., 0., 1.);
  glScalef(2.5 * sc->dmax, 2.5 * sc->dmax, 2.5 * sc->dmax);
  glDisable(GL_LIGHTING);

  glColor3f(0.4, 0.4, 0.4);
  glCallList(sc->grid);

  glColor3fv(sc->par.line);
  output3(0.0, 0.0, 0.0, "%.2f", mesh->xmin);
  output3(1.1, 0.0, 0.0, "%.2f", mesh->xmax);
  output3(0.0, 1.01, 0.0, "%.2f", mesh->ymax);
  output3(0.0, 0.0, 1.01, "%.2f", mesh->zmax);
  glEnable(GL_LIGHTING);
  glPopMatrix( );
}

void drawBase(pScene sc, pMesh mesh) {
  /* default */
  if (ddebug) printf("draw base\n");

  if (!sc->grid) {
    int k;

    sc->grid = glGenLists(1);
    glNewList(sc->grid, GL_COMPILE);
    if (glGetError( )) return;

    glColor3f(0.5, 0.5, 0.5);
    glLineWidth(2.0);
    glBegin(GL_LINES);

    for (k = 0; k < 21; k += 5) {
      glVertex3f(k * 0.05, 0., 0.);
      glVertex3f(k * 0.05, 1., 0.);
      glVertex3f(0., k * 0.05, 0.);
      glVertex3f(1., k * 0.05, 0.);
    }

    glEnd( );
    glColor3f(0.6, 0.6, 0.6);
    glLineWidth(1.0);
    glBegin(GL_LINES);

    for (k = 0; k < 21; k++) {
      if (k % 5 == 0) continue;

      glVertex3f(k * 0.05, 0., 0.);
      glVertex3f(k * 0.05, 1., 0.);
      glVertex3f(0., k * 0.05, 0.);
      glVertex3f(1., k * 0.05, 0.);
    }

    glEnd( );
    glEndList( );
  }

  glPushMatrix( );
  glTranslatef(-1.5 * sc->dmax, -1.5 * sc->dmax, -0.5 * (mesh->zmax - mesh->zmin));
  glScalef(3 * sc->dmax, 3 * sc->dmax, 3 * sc->dmax);
  glDisable(GL_LIGHTING);
  glCallList(sc->grid);
  glPopMatrix( );
}

/* draw HUD system for flight */
void drawHUD(pScene sc) {
  pCamera c;
  pMesh mesh;
  GLfloat x, y, dx, dy, alt;
  double azim, elev;
  int i, j;

  if (ddebug) fprintf(stdout, "drawHUD\n");

  glDisable(GL_DEPTH_TEST);
  glDisable(GL_LIGHTING);

  glMatrixMode(GL_PROJECTION);
  glPushMatrix( );
  glLoadIdentity( );
  gluOrtho2D(-0.5, 639.5, -0.5, 479.5);

  c = sc->camera;
  mesh = cv.mesh[sc->idmesh];
  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  glColor3f(0.0, 0.5, 0.0);
  glLineWidth(1.0f);
  glRecti(200, 160, 440, 320);

  /* altitude */
  glColor3f(1.0, 0.5, 0.0);
  output2(230, 330, "speed: %6.1f", 1000.0f * c->spmod / sc->dmax);

  switch (c->vecup) {
    case X_AXIS:
      alt = (c->eye[1] + mesh->xtra) / (mesh->xmax - mesh->xmin);
      break;
    case Y_AXIS:
      alt = (c->eye[1] + mesh->ytra) / (mesh->ymax - mesh->ymin);
      break;
    case Z_AXIS:
      alt = (c->eye[1] + mesh->ztra) / (mesh->zmax - mesh->zmin);
      break;
    default:
      alt = 0.0;
      break;
  }

  glColor3f(1.0, 0.5, 0.0);
  output2(350, 330, "alt: %9.1f", 1000.0f * alt);

  /* horiz rulers */
  output2(310, 139.0f, "azim");
  glColor3f(0.0, 1.0, 0.0);
  output2(197, 150, "-180");
  output2(257, 150, " -90");
  output2(317, 150, "0");
  output2(377, 150, " 90");
  output2(437, 150, "180");
  x = 200.0f;
  glBegin(GL_LINES);

  for (i = 1; i < 8; i++) {
    x += 240.0 / 8.0;
    glVertex2f(x, 158.0);
    glVertex2f(x, 162.0);
  }

  glEnd( );

  /* vert rulers */
  glColor3f(0.0, 1.0, 0.0);
  output2(185, 160, "-90");
  output2(185, 200, "-45");
  output2(185, 240, "0");
  output2(185, 280, "45");
  output2(185, 320, "90");
  y = 160.0f;
  glBegin(GL_LINES);

  for (i = 1; i < 8; i++) {
    y += 160.0 / 8.0;
    glVertex2f(198, y);
    glVertex2f(202, y);
  }

  glEnd( );

  /* azimuth */
  azim = Azimuth(c);
  if (azim > 0.0f)
    azim = 180.0f - azim;
  else if (azim < 0.0f)
    azim = -180.0f - azim;

  x = 2.0 / 3.0 * azim + 320.0f;
  glColor3f(1.0, 0.0, 0.0);
  glLineWidth(1.0);
  output2(x, 143.0, "%d", azim > 0 ? (int)(azim + 0.5) : (int)(azim - 0.5));
  glBegin(GL_LINES);
  glVertex2f(x, 166.0);
  glVertex2f(x, 318.0);
  glEnd( );
  y = 160.0f;
  dy = 160.0 / 8.0;
  glBegin(GL_LINES);

  for (i = 0; i < 8; i++) {
    glVertex2f(x - 4, y);
    glVertex2f(x + 4, y);

    for (j = 0; j < 5; j++) {
      glVertex2f(x - 2, y + j * dy / 5.0);
      glVertex2f(x + 2, y + j * dy / 5.0);
    }

    y += dy;
  }

  glEnd( );

  /* elevation */
  elev = Elevation(c);
  if (elev > 90.0f)
    y = 320.0f;
  else if (elev < -90.0f)
    y = 160.0f;
  else
    y = 8.0 / 9.0 * elev + 240.0;

  glColor3f(1.0, 0.0, 0.0);
  output2(175.0, y, "%5.1f", elev);
  glBegin(GL_LINES);
  glVertex2f(206.0, y);
  glVertex2f(438.0, y);
  glEnd( );
  x = 200.0f;
  dx = 240.0f / 8.0f;
  glBegin(GL_LINES);

  for (i = 1; i <= 8; i++) {
    glVertex2f(x, y - 4);
    glVertex2f(x, y + 4);

    for (j = 0; j < 5; j++) {
      glVertex2f(x + j * dx / 5.0, y - 2);
      glVertex2f(x + j * dx / 5.0, y + 2);
    }

    x += dx;
  }

  glEnd( );

  /* horizon */
  glLineWidth(2.0f);
  glColor3f(0.0, 0.0, 1.0);
  glBegin(GL_LINES);
  glVertex2f(200.0, 240.0);
  glVertex2f(440.0, 240.0);
  glEnd( );
  glLineWidth(1.0f);

  /* HUD */
  glColor3f(1.0, 0.0, 0.0);
  glBegin(GL_LINES);
  glVertex2f(310, 230);
  glVertex2f(330, 230);
  glVertex2f(310, 250);
  glVertex2f(330, 250);

  glVertex2f(310, 230);
  glVertex2f(310, 234);
  glVertex2f(330, 230);
  glVertex2f(330, 234);

  glVertex2f(310, 246);
  glVertex2f(310, 250);
  glVertex2f(330, 246);
  glVertex2f(330, 250);
  glEnd( );
  glColor3f(0.0, 1.0, 0.0);
  glLineWidth(3.0f);
  glBegin(GL_LINES);
  glVertex2f(320., 235.);
  glVertex2f(320., 245.);
  glVertex2f(315., 240.);
  glVertex2f(325., 240.);
  glEnd( );
  glLineWidth(1.0f);

  glMatrixMode(GL_PROJECTION);
  glPopMatrix( );
  glMatrixMode(GL_MODELVIEW);
  glEnable(GL_DEPTH_TEST);
}

/* ---------------------------------------------------------------------------------------
   value of the solution at the point under the mouse (key .), on the bottom left corner:
   2d: the point of the mesh (or of the elevation) under the cursor, 3d: the visible point,
   or the first translucent iso-surface crossed by the ray of the eye
   --------------------------------------------------------------------------------------- */
int probeon = 1;
static int probex = -1, probey = -1, probewin = -1;
static GLdouble pmv[16], ppj[16];
static GLint pvp[4];

/* bucket grid of the elements (triangles in 2d, tetrahedra in 3d) to locate a point */
static struct {
  pMesh mesh;
  void *elt;
  int ne, np, dim, n[3], *head, *list;
  double mn[3], h[3], key;
} pgrid = {0};

void probeMotion(int x, int y) {
  probex = x;
  probey = y;
  probewin = glutGetWindow( );
  if (probeon) glutPostRedisplay( );
}

/* the window under the mouse (the window of a right click = of a menu, see menuDefer) */
int mouseWindow(void) { return probewin; }

void probeEntry(int state) {
  if (state == GLUT_ENTERED) probewin = glutGetWindow( );

  if (state == GLUT_LEFT) {
    probex = probey = -1;
    if (probeon) glutPostRedisplay( );
  }
}

void probeMatrices(void) {
  glGetDoublev(GL_MODELVIEW_MATRIX, pmv);
  glGetDoublev(GL_PROJECTION_MATRIX, ppj);
  glGetIntegerv(GL_VIEWPORT, pvp);
}

static int probeVertex(pMesh mesh, int k, int i) {
  return mesh->dim == 3 ? mesh->tetra[k].v[i] : mesh->tria[k].v[i];
}

static void probeBuildGrid(pMesh mesh) {
  int dim = mesh->dim == 3 && mesh->ntet ? 3 : 2, ne = dim == 3 ? mesh->ntet : mesh->nt;
  int k, i, c, nc, m, nv = dim + 1, *cnt;
  double mx[3], key = 0.0;

  /* a new mesh of the server can reuse the same memory: a key with some coordinates */
  for (k = 1; k <= mesh->np; k += 1 + mesh->np / 16)
    key += mesh->point[k].c[0] + 3.1 * mesh->point[k].c[1] + 7.3 * mesh->point[k].c[2];

  if (pgrid.mesh == mesh && pgrid.elt == (dim == 3 ? (void *)mesh->tetra : (void *)mesh->tria) &&
      pgrid.ne == ne && pgrid.np == mesh->np && pgrid.key == key)
    return;

  free(pgrid.head);
  free(pgrid.list);
  pgrid.head = pgrid.list = NULL;
  pgrid.mesh = mesh;
  pgrid.elt = dim == 3 ? (void *)mesh->tetra : (void *)mesh->tria;
  pgrid.ne = ne;
  pgrid.np = mesh->np;
  pgrid.key = key;
  pgrid.dim = dim;
  if (!ne) return;

  for (c = 0; c < 3; c++) {
    pgrid.mn[c] = 1e300;
    mx[c] = -1e300;
  }

  for (k = 1; k <= mesh->np; k++)
    for (c = 0; c < dim; c++) {
      pgrid.mn[c] = min(pgrid.mn[c], mesh->point[k].c[c]);
      mx[c] = max(mx[c], mesh->point[k].c[c]);
    }

  nc = 1;
  for (c = 0; c < 3; c++) {
    pgrid.n[c] = 1;
    pgrid.h[c] = 1.0;
    if (c < dim) {
      pgrid.n[c] = max(1, (int)pow((double)ne, 1.0 / dim));
      if (pgrid.n[c] > 200) pgrid.n[c] = 200;
      pgrid.h[c] = (mx[c] - pgrid.mn[c]) / pgrid.n[c];
      if (!(pgrid.h[c] > 0)) pgrid.h[c] = 1.0;
      pgrid.mn[c] -= 1e-9 * pgrid.h[c];
    } else {
      pgrid.mn[c] = 0.0;
    }

    nc *= pgrid.n[c];
  }

  /* 2 passes: count, then fill (each element in the cells of its bounding box) */
  cnt = (int *)calloc(nc + 1, sizeof(int));
  pgrid.head = (int *)calloc(nc + 1, sizeof(int));
  for (m = 0; m < 2; m++) {
    for (k = 1; k <= ne; k++) {
      int lo[3] = {0, 0, 0}, hi[3] = {0, 0, 0}, a, b, d;

      if (!probeVertex(mesh, k, 0)) continue;

      for (c = 0; c < dim; c++) {
        double e0 = 1e300, e1 = -1e300;

        for (i = 0; i < nv; i++) {
          double v = mesh->point[probeVertex(mesh, k, i)].c[c];

          e0 = min(e0, v);
          e1 = max(e1, v);
        }

        lo[c] = max(0, min(pgrid.n[c] - 1, (int)((e0 - pgrid.mn[c]) / pgrid.h[c])));
        hi[c] = max(0, min(pgrid.n[c] - 1, (int)((e1 - pgrid.mn[c]) / pgrid.h[c])));
      }

      for (a = lo[0]; a <= hi[0]; a++)
        for (b = lo[1]; b <= hi[1]; b++)
          for (d = lo[2]; d <= hi[2]; d++) {
            int cell = a + pgrid.n[0] * (b + pgrid.n[1] * d);

            if (m == 0)
              cnt[cell + 1]++;
            else
              pgrid.list[pgrid.head[cell] + cnt[cell]++] = k;
          }
    }

    if (m == 0) {
      for (i = 0; i < nc; i++) cnt[i + 1] += cnt[i];
      memcpy(pgrid.head, cnt, (nc + 1) * sizeof(int));
      pgrid.list = (int *)malloc((cnt[nc] + 1) * sizeof(int));
      memset(cnt, 0, (nc + 1) * sizeof(int));
    }
  }

  free(cnt);
}

/* barycentric coordinates of P in the element k */
static void probeBary(pMesh mesh, int k, const double *P, double *l) {
  double *a = mesh->point[probeVertex(mesh, k, 0)].c;

  if (pgrid.dim == 2) {
    double *b = mesh->point[probeVertex(mesh, k, 1)].c, *c = mesh->point[probeVertex(mesh, k, 2)].c;
    double det = (b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1]);

    if (det == 0) det = 1e-300;
    l[1] = ((P[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (P[1] - a[1])) / det;
    l[2] = ((b[0] - a[0]) * (P[1] - a[1]) - (P[0] - a[0]) * (b[1] - a[1])) / det;
    l[0] = 1.0 - l[1] - l[2];
  } else {
    double *p[3], m[3][3], d[3], det;
    int i, j;

    for (i = 0; i < 3; i++) p[i] = mesh->point[probeVertex(mesh, k, i + 1)].c;

    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) m[j][i] = p[i][j] - a[j];

      d[i] = P[i] - a[i];
    }

    det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
          m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
          m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    if (det == 0) det = 1e-300;

    /* Cramer */
    for (i = 0; i < 3; i++) {
      double n[3][3];
      int r;

      for (r = 0; r < 3; r++)
        for (j = 0; j < 3; j++) n[r][j] = j == i ? d[r] : m[r][j];

      l[i + 1] = (n[0][0] * (n[1][1] * n[2][2] - n[1][2] * n[2][1]) -
                  n[0][1] * (n[1][0] * n[2][2] - n[1][2] * n[2][0]) +
                  n[0][2] * (n[1][0] * n[2][1] - n[1][1] * n[2][0])) /
                 det;
    }

    l[0] = 1.0 - l[1] - l[2] - l[3];
  }
}

/* the element containing P (display coordinates), 0 if none */
static int probeLocate(pMesh mesh, const double *P, double *l) {
  int c, cell = 0, i, best = 0;
  double bestmin = -1e-3, ll[4];

  probeBuildGrid(mesh);
  if (!pgrid.ne || !pgrid.head) return (0);

  for (c = pgrid.dim - 1; c >= 0; c--) {
    int ic = (int)floor((P[c] - pgrid.mn[c]) / pgrid.h[c]);

    if (ic < 0 || ic >= pgrid.n[c]) return (0);

    cell = cell * pgrid.n[c] + ic;
  }

  for (i = pgrid.head[cell]; i < pgrid.head[cell + 1]; i++) {
    int k = pgrid.list[i], j;
    double mn = 1e300;

    probeBary(mesh, k, P, ll);
    for (j = 0; j <= pgrid.dim; j++) mn = min(mn, ll[j]);

    if (mn > bestmin) {
      bestmin = mn;
      best = k;
      memcpy(l, ll, sizeof(ll));
    }
  }

  return (best);
}

/* value (and components) of the solution at P: the number of components, 0 if outside */
static int probeValue(pMesh mesh, const double *P, double *val, double *vec) {
  double l[4];
  int k = probeLocate(mesh, P, l), i, j, nc;

  if (!k) return (0);

  nc = mesh->nfield == mesh->dim && mesh->nfield > 1 ? mesh->dim : 1;
  *val = 0.0;
  for (j = 0; j < 3; j++) vec[j] = 0.0;

  if (mesh->typage == 1) { /* one value by element */
    pSolution ps = &mesh->sol[k];

    *val = ps->bb;
    if (nc > 1)
      for (j = 0; j < nc; j++) vec[j] = ps->m[j];
  } else {
    for (i = 0; i <= pgrid.dim; i++) {
      pSolution ps = &mesh->sol[probeVertex(mesh, k, i)];

      *val += l[i] * ps->bb;
      if (nc > 1)
        for (j = 0; j < nc; j++) vec[j] += l[i] * ps->m[j];
    }

    if (nc > 1) {
      *val = 0.0;
      for (j = 0; j < nc; j++) *val += vec[j] * vec[j];

      *val = sqrt(*val);
    }
  }

  return (nc);
}

void drawProbe(pScene sc) {
  pMesh mesh = cv.mesh[sc->idmesh];
  GLdouble ax, ay, az, bx, by, bz;
  GLfloat d = 1.0f;
  double P[3], val = 0.0, vec[3];
  int xp, yp, nc = 0, oniso = 0, l;
  char buf[256];

  if (!probeon || probewin != sc->idwin || probex < 0 || !mesh->nbb || graph1dMesh(mesh)) return;

  if (!(mesh->dim == 3 ? mesh->ntet : mesh->nt)) return;

  xp = probex;
  yp = pvp[3] - 1 - probey;
  if (xp < 0 || yp < 0 || xp >= pvp[2] || yp >= pvp[3]) return;

  glReadPixels(xp, yp, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &d);
  gluUnProject(xp, yp, 0.0, pmv, ppj, pvp, &ax, &ay, &az);
  gluUnProject(xp, yp, d < 1.0f ? d : 1.0, pmv, ppj, pvp, &bx, &by, &bz);

  /* 3d: the translucent iso-surfaces are not in the depth buffer: march along the ray
     of the eye up to the visible point, the first crossing of an iso-value */
  if (mesh->dim == 3 && (sc->isotyp & S_ISOSURF) && mesh->nfield == 1) {
    double A[3] = {ax, ay, az}, B[3] = {bx, by, bz}, Q[3], u0 = 0.0, u1, t0 = 0.0;
    int i, s, n = 400, in0 = 0;

    for (s = 0; s <= n && !oniso; s++) {
      double t = (double)s / n;

      for (i = 0; i < 3; i++) Q[i] = A[i] + t * (B[i] - A[i]);

      if (!probeValue(mesh, Q, &u1, vec)) {
        in0 = 0;
        continue;
      }

      if (in0)
        for (i = 0; i < MAXISO && !oniso; i++) {
          double iso = sc->iso.val[i];

          if ((u0 - iso) * (u1 - iso) <= 0.0 && u0 != u1) {
            /* dichotomy between t0 and t */
            double ta = t0, tb = t, ua = u0, um;
            int it;

            for (it = 0; it < 30; it++) {
              double tm = 0.5 * (ta + tb);
              int j;

              for (j = 0; j < 3; j++) Q[j] = A[j] + tm * (B[j] - A[j]);

              if (!probeValue(mesh, Q, &um, vec)) break;

              if ((ua - iso) * (um - iso) <= 0.0)
                tb = tm;
              else
                ta = tm, ua = um;
            }

            for (int j = 0; j < 3; j++) P[j] = A[j] + 0.5 * (ta + tb) * (B[j] - A[j]);

            val = iso;
            nc = 1;
            oniso = 1;
          }
        }

      u0 = u1;
      t0 = t;
      in0 = 1;
    }
  }

  if (!oniso && mesh->dim == 2 && !(sc->mode & S_ALTITUDE)) {
    /* flat 2d mesh (not always in the depth buffer): the ray of the eye on the plane z = 0 */
    gluUnProject(xp, yp, 1.0, pmv, ppj, pvp, &bx, &by, &bz);
    if (fabs(bz - az) < 1e-30) return;

    {
      double t = -az / (bz - az);

      P[0] = ax + t * (bx - ax);
      P[1] = ay + t * (by - ay);
      P[2] = 0.0;
    }

    nc = probeValue(mesh, P, &val, vec);
    if (!nc) return;
  } else if (!oniso) {
    if (!(d < 1.0f)) return;

    P[0] = bx;
    P[1] = by;
    P[2] = bz;
    nc = probeValue(mesh, P, &val, vec);
    if (!nc) return;
  }

  /* real coordinates (the mesh is centered on 0) */
  if (mesh->dim == 3)
    l = snprintf(buf, sizeof(buf), "x = %.4g  y = %.4g  z = %.4g", P[0] + mesh->xtra,
                 P[1] + mesh->ytra, P[2] + mesh->ztra);
  else
    l = snprintf(buf, sizeof(buf), "x = %.4g  y = %.4g", P[0] + mesh->xtra, P[1] + mesh->ytra);

  if (nc == 1 && mesh->cplx)
    snprintf(buf + l, sizeof(buf) - l, oniso ? "   :  %.6g   (%s, isosurface)" : "   :  %.6g   (%s)", val,
             cplxname[cplxpart]);
  else if (nc == 1)
    snprintf(buf + l, sizeof(buf) - l, oniso ? "   :  %.6g   (isosurface)" : "   :  %.6g", val);
  else if (nc == 2)
    snprintf(buf + l, sizeof(buf) - l, "   :  (%.4g, %.4g)   |.| = %.4g", vec[0], vec[1], val);
  else
    snprintf(buf + l, sizeof(buf) - l, "   :  (%.4g, %.4g, %.4g)   |.| = %.4g", vec[0], vec[1],
             vec[2], val);

  glMatrixMode(GL_PROJECTION);
  glPushMatrix( );
  glLoadIdentity( );
  gluOrtho2D(0., pvp[2], 0., pvp[3]);
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix( );
  glLoadIdentity( );
  glDisable(GL_DEPTH_TEST);
  glDisable(GL_LIGHTING);
  glColor3f(1. - sc->par.back[0], 1.0 - sc->par.back[1], 1.0 - sc->par.back[2]);
  graphText(6.0, (sc->type & S_DECO) ? 26.0 : 6.0, buf, 0);
  glMatrixMode(GL_PROJECTION);
  glPopMatrix( );
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix( );
  glEnable(GL_DEPTH_TEST);
}

#ifdef __cplusplus
}
#endif
