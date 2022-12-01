typedef struct {
  /*
  double         aire,rins;
  double         c[3];
  int            s,nv,mark;
  unsigned char  b,h;
  */
  double         c[3];
  int            s,nv,mark;
  unsigned char  b,h;
} MSHMET_Point;
typedef MSHMET_Point * MSHMET_pPoint;

typedef struct {
  int     v[3];
  int     mark;
  /*
  double  aire;
  int     v[3];
  int     mark;*/
} MSHMET_Tria;
typedef MSHMET_Tria * MSHMET_pTria;

typedef struct {
  int     v[4];
  int     mark;
} MSHMET_Tetra;
typedef MSHMET_Tetra * MSHMET_pTetra;

typedef struct {
  double   delta;
  double   min[3],max[3];
  float    eps,hmin,hmax,width,ani,hgrad,map;
  int      nnu,nsol,nlis;
  char     imprim,ddebug,iso,bin,metric,ls;

  /* double   delta;
  double   min[3],max[3];
  float    eps,hmin,hmax,width;
  int      nnu,nsol,nlis;
  char     imprim,ddebug,iso,bin,metric,ls; */
} MSHMET_Info;

typedef struct {
  /*
  int      np,nt,ne,ver,dim;
  int     *adja,mark;
  char    *name,*mname;

  MSHMET_pPoint   point;
  MSHMET_pTria    tria;
  MSHMET_pTetra   tetra;
  MSHMET_Info info;
  */

  int      np,nt,ne,ver,dim;
  int     *adja,mark;
  char    *name,*mname;

  MSHMET_pPoint   point;
  MSHMET_pTria    tria;
  MSHMET_pTetra   tetra;
  MSHMET_Info     info;
} MSHMET_Mesh;

typedef MSHMET_Mesh * MSHMET_pMesh;

typedef struct {
  int         np,ver,dim,type,size,typtab[GmfMaxTyp];
  double     *sol,*met,*hes,*grd,*nn,umin,umax;
  char       *name,*outn,*mapname;
  /*  version 2.0 
  int         np,ver,dim,type,size,typtab[GmfMaxTyp];
  double     *sol,*met,umin,umax;
  char       *name,*outn; */
} MSHMET_Sol;
typedef MSHMET_Sol * MSHMET_pSol;

typedef struct {
  double         grd[3];
  double         hes[6];
} MSHMET_Deriv;
typedef MSHMET_Deriv * MSHMET_pDeriv;

#ifdef  __cplusplus
namespace mshmet{
extern "C" {
#endif
  int MSHMET_mshmet(int intopt[7], double fopt[4], MSHMET_pMesh mesh, MSHMET_pSol sol);
#ifdef  __cplusplus
}}
#endif
