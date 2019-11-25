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

#ifndef _MEDIT_H_
#define _MEDIT_H_

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <float.h>
#include <signal.h>
#include <ctype.h>

#ifdef WIN32
#include <windows.h>
#endif
#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#endif

#include "chrono.h"
#include "memory.h"
#include "mesh.h"
#include "grafic.h"
#include "image.h"
#include "sproto.h"

#define ME_VER "3.0a"
#define ME_REL "Nov. 30, 2007"
#define ME_CPY "Copyright (c) LJLL, 1999-2007"
#define ME_STR "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"
#define DEFAULT_FILE "DEFAULT.medit"

#define MAX_MESH 32
#define MAX_SCENE 32
#define MAX_OBJ 32
#define MAX_MATERIAL 128
#define DEFAULT_MAT 0

#define DTOR 0.0174532925
#define RTOD 57.29577951308232
#define EPS 1.e-06
#define EPS2 2.e-10

#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif
#ifdef M_PI
#undef M_PI
#undef M_PI_2
#endif
#define M_PI 3.14159265358979323846   /* pi   */
#define M_PI_2 1.57079632679489661923 /* pi/2 */

#ifdef min
#undef min
#undef max
#endif
#define min(a, b) (((a) < (b)) ? (a) : (b))
#define max(a, b) (((b) > (a)) ? (b) : (a))

/* check if numbers are equal */
#define egal(x, y)      \
  ((((x) == 0.0f)       \
      ? (fabs(y) < EPS) \
      : (((y) == 0.0f) ? (fabs(x) < EPS) : (fabs((x) - (y)) / (fabs(x) + fabs(y)) < EPS2))))

/* options */
enum { STANDARD = 1, SEQUENCE, VERYBIG, MORPHING, SCHNAUZER, ISOSURF, PARTICLE };

/* structure canvas */
typedef struct canvas {
  pMesh mesh[MAX_MESH];
  pScene scene[MAX_SCENE];
  int nbm, nbs;
} Canvas;
typedef Canvas *pCanvas;

#endif

#ifdef __cplusplus
}
#endif
