#ifndef _MSHMET_H
#define _MSHMET_H

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <math.h>
#include <string.h>
#include <signal.h>
#include <ctype.h>

#include "chrono.h"
#include "memory.h"
#include "libmeshb7.h"
#include "eigenv.h"

#define MS_VER   "3.0a"
#define MS_REL   "Feb. 15, 2010"
#define MS_CPY   "Copyright (c) LJLL, 2007-10"
#define MS_STR   "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"

#define PRECI    1.0
#define EPS      1.e-6
#define EPS1     1.e-20
#define EPST    -1.e-2
#define EPSR     1.e+2
#define LONMAX   4096

#define CTE2D    2.0 / 9.0
#define CTE3D    9.0 / 32.0
  
#define MS_MAX(a,b)   ( ((a) < (b)) ? (b) : (a) )
#define MS_MIN(a,b)   ( ((a) < (b)) ? (a) : (b) )

extern char idir[5];

typedef struct {
  double         c[3];
  int            s,nv,mark;
  unsigned char  b,h;
} Point;
typedef Point * pPoint;

typedef struct {
  int     v[3];
  int     mark;
} Tria;
typedef Tria * pTria;

typedef struct {
  int     v[4];
  int     mark;
} Tetra;
typedef Tetra * pTetra;

typedef struct {
  double   delta;
  double   min[3],max[3];
  float    eps,hmin,hmax,width,ani,hgrad,map;
  int      nnu,nsol,nlis;
  char     imprim,ddebug,iso,bin,metric,ls;
} Info;

typedef struct {
  int   min,max,nxt;
} hedge;

typedef struct {
  int     size,nhmax,hnxt;
  hedge  *item;
} Hash;
typedef Hash * pHash; 

typedef struct {
  int      np,nt,ne,ver,dim;
  int     *adja,mark;
  char    *name,*mname;

  pPoint   point;
  pTria    tria;
  pTetra   tetra;
  Info     info;
} Mesh;
typedef Mesh * pMesh;

typedef struct {
  int         np,ver,dim,type,size,typtab[GmfMaxTyp];
  double     *sol,*met,*hes,*grd,*nn,umin,umax;
  char       *name,*outn,*mapname;
} Sol;
typedef Sol * pSol;


/* prototypes */
int loadMesh(pMesh ,char *);
int loadSol(pSol ,Info *,char *);
int loadMetric(pSol ,Info *,char *);
int saveMet(pSol ,Info *,char *);
int saveSol(pSol ,Info *,char *);
int outder(pMesh ,pSol );
int scaleMesh(pMesh ,pSol );
int unscaleMesh(pMesh ,pSol );

int mshme1(pMesh ,pSol );

/* function pointers */
int   boulep_3d(pMesh ,int ,int ,int *);
int   boulep_2d(pMesh ,int ,int ,int *);
int   hashel_3d(pMesh );
int   hashel_2d(pMesh );
int   gradLS_3d(pMesh ,pSol ,int ,int );
int   gradLS_2d(pMesh ,pSol ,int ,int );
int   gradLS_s(pMesh ,pSol ,int ,int );
int   hessLS_3d(pMesh ,pSol ,int ,int );
int   hessLS_2d(pMesh ,pSol ,int ,int );
int   hessLS_s(pMesh ,pSol ,int ,int );
int   avgval_3d(pMesh ,pSol ,int );
int   avgval_2d(pMesh ,pSol ,int );
int   clsval_3d(pMesh ,pSol ,int );
int   clsval_2d(pMesh ,pSol ,int );
double getSol_3d(pSol ,int ,int );
double getSol_2d(pSol ,int ,int );
int   nrmhes_3d(pMesh ,pSol ,int );
int   nrmhes_2d(pMesh ,pSol ,int );
int   redsim_3d(double *,double *,double *);
int   redsim_2d(double *,double *,double *);
int   defmet_3d(pMesh ,pSol ,int );
int   defmet_2d(pMesh ,pSol ,int );
int   defmet_s(pMesh ,pSol ,int );
int   metrLS_3d(pMesh mesh,pSol );
int   metrLS_2d(pMesh mesh,pSol );
int   lissag_3d(pMesh mesh,pSol sol,int is,int it);
int   lissag_2d(pMesh mesh,pSol sol,int is,int it);
pHash hashEdge_3d(pMesh mesh);
pHash hashEdge_2d(pMesh mesh);

extern int   (*boulep)(pMesh ,int ,int ,int *);
extern int   (*hashel)(pMesh );
extern int   (*gradLS)(pMesh ,pSol ,int ,int );
extern int   (*hessLS)(pMesh ,pSol ,int ,int );
extern int   (*avgval)(pMesh ,pSol ,int );
extern int   (*clsval)(pMesh ,pSol ,int );
extern int   (*nrmhes)(pMesh ,pSol ,int );
extern int   (*redsim)(double *,double *,double *);
extern int   (*defmet)(pMesh ,pSol ,int );
extern double (*getSol)(pSol ,int ,int );
extern int   (*metrLS)(pMesh mesh,pSol );
extern int   (*lissag)(pMesh ,pSol , int ,int );

#endif
