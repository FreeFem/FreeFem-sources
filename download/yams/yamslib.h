#ifdef __cplusplus
extern "C" {
#endif

/* Edge: Structure used to store specified mesh edges */
typedef struct yams_sedge {
  int     p1,p2;
  int     ref;
  int     tag;
} yams_Edge;
typedef yams_Edge * yams_pEdge;

#ifndef  ERR
#define  ERR   1
#define  WAR   2
#define  MSG   3
#endif

typedef struct yams_error {
  double cooerr[6];
  int    inderr[6];
  int    lerror;
  int    coderr;
} yams_Error;

#include "defines.h"

/* HashTable: hash table structure for mesh edges */
typedef struct yams_shashtab {
  int    min;           /* min(a,b)         */
  int    nxt;           /* next edge        */
  int    elt;
  int    ind;
} yams_Hashtable;
typedef yams_Hashtable * yams_pHashtable;

typedef struct yams_sstack {
  int    *t;
  int     in,out,cur;
} yams_Stack;
typedef yams_Stack * yams_pStack;

typedef struct yams_sinfo {
  double   xmin,ymin,zmin,xmax,ymax,zmax;  /* bounding box  */
  double   delta;
  double   dmin,dmax;                      /* edge lengths  */
  float    qworst;
  int      meshtype,cc,flip;
  long     nulp,nulf,nuln;                 /* not used      */
  int      qpire;
  int      nedg,nrid,ncoi,nreq,nvus;
  int      nafixe,nvrequis,ndang;
  int      manifold;
} yams_Info;

typedef struct yams_soptions {
  float    hmin,hmax;      /* desired sizes    */
  float    kmin,kmax;      /* curvature min,max*/
  float    eps,iso;        /* max. tolerance, isovalue   */
  float    alpha,gap;      /* max values allow.*/
  float    degrad;         /* max degrad. qual */
  float    ridge;          /* cosine ridge ang */
  float    geom;
  float    shock;          /* mesh gradation   */
  float    bande;          /* bandwidth        */
  float    walton;         /* angle limitation */
  float    declic;
  float    lambda,mu;      /* for smoothing    */
  int      ctrl;           /* absolute values  */
  int      minnp;
  short    iter,choix;
  unsigned char ptmult,noreff,ffem,check;
} yams_Options;

#ifndef ubyte
typedef unsigned char  ubyte;
#endif

/* Point: Structure that defines a vertex in a mesh. */
typedef struct yams_spoint {
  float   c[3];            /* coordinates           */
  float   size;            /* calculated size       */
  int     tge;             /* tangent at ridge      */
  short   color;
  int     ref;
  int     tmp;
  ubyte   tag;             /* vertex type           */
  ubyte   geom;
  ubyte   flag;
} yams_Point;
typedef yams_Point     * yams_pPoint;

/* Triangle: Structure that defines a triangle in a mesh. */
typedef struct yams_striangle {
  float   n[3];             /* face normal                */
  float   dish;             /* distance to surface        */
  float   qual;             /* triangle quality           */

  int     v[3];             /* array of vertex indices    */
  int     adj[3];           /* array of adjacent trias    */
  int     vn[3];            /* array of vertex normals    */
  int     edg[3];
  int     nxt;
  int     ref;
  short   cc;

  ubyte   voy[3];           /* array of voyeur vertices  */
  ubyte   flag1;
  ubyte   tag[3];           /* array of edge classes     */
  ubyte   flag2;
} yams_Triangle;
typedef yams_Triangle  * yams_pTriangle;

typedef struct yams_squad {
  float    qual;
  float    n[3];
  int      v[4];
  int      adj[4];
  int      ref,edg[4],vn[4];
  short    cc;
  ubyte    flag1,flag2;
  ubyte    voy[4];
  ubyte    tag[4];
} yams_Quad;
typedef yams_Quad * yams_pQuad;

typedef struct {
  int   v[4];
  int   ref;
} yams_Tetra;
typedef yams_Tetra * yams_pTetra;

typedef struct yams_geomsupp {
  float    vn[3];           /* array of vertex normals  */
  float    gap;             /* local gap value          */
  int      newnum;             /* pointer to new number    */
} yams_GeomSupp;
typedef yams_GeomSupp  * yams_pGeomSupp;

typedef struct yams_geomtge {
  float    t[3];
  int      newnum;
} yams_Geomtge;
typedef yams_Geomtge   * yams_pGeomtge;

typedef struct yams_metric {
  float k1,k2;
  float m[6];                /* anisotropic metric */
} yams_Metric;
typedef yams_Metric   * yams_pMetric;


/* SurfMesh: Structure that defines a mesh. */
typedef struct yams_smesh {
  int       dim;                /* mesh dimension (2,3)     */
  int       type;
  int       connex;             /* # connected component    */
  int       np,npfixe,npmax;    /* number of vertices       */
  int       ne,nefixe,nemax;    /* number of triangles      */
  int       nq,ntet;            /* quads, ntets             */
  int       nv,nvfixe,nvmax;    /* number of vertex normals */
  int       nafixe,nmfixe;
  int       nt,ntfixe,ntmax;    /* vertex tgtes             */
  int       mark;               /* coloring...              */
  int       ipil;

  char     *infile;
  char     *outfile;

  yams_pPoint    point;          /* array of vertices         */
  yams_pTriangle tria;           /* array of triangles        */
  yams_pTetra    tetra;
  yams_pQuad     quad;
  yams_pGeomSupp geom;           /* pointer to geometric info */
  yams_pGeomtge  tgte;           /* pointer to tge at ridge   */
  yams_pMetric   metric;         /* local metric at vertex    */

  yams_pEdge     edge;
} yams_SurfMesh;
typedef yams_SurfMesh  * yams_pSurfMesh;

#ifdef __cplusplus
namespace yams{
#endif
  int yams_main(yams_pSurfMesh sm, int intopt[23], double fopt[14], int infondang, int infocc );
  int zaldy1(int nemax,int npmax,int nvmax,int memory,yams_pSurfMesh sm,int choix);
#ifdef __cplusplus
}
#endif
#ifdef __cplusplus
}
#endif

