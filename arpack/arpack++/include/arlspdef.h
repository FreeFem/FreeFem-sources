/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE ARLSpDef.h.
   ALTERED version of ssp_defs.h, dsp_defs.h, csp_defs.h 
   and zsp_defs.h (from SuperLU package).
*/


#ifndef __SUPERLU_SP_DEFS /* allow multiple inclusions */
#define __SUPERLU_SP_DEFS

/*
 * File name:		sp_defs.h
 * Purpose:             Sparse matrix types and function prototypes
 * History:
 */
#include "arlnames.h"
#include "arlsupm.h"
#include "arlcomp.h"

#define USE_VENDOR_BLAS

/* No of marker arrays used in the symbolic factorization,
   each of size n */
#define NO_MARKER     3
#define NUM_TEMPV(m,w,t,b)  ( MAX(m, (t + b)*w) )

typedef enum {LUSUP, UCOL, LSUB, USUB} MemType;
typedef enum {HEAD, TAIL}              stack_end_t;
typedef enum {SYSTEM, USER}            LU_space_t;

/*
 * Global data structures used in LU factorization -
 * 
 *   nsuper: #supernodes = nsuper + 1, numbered [0, nsuper].
 *   (xsup,supno): supno[i] is the supernode no to which i belongs;
 *	xsup(s) points to the beginning of the s-th supernode.
 *	e.g.   supno 0 1 2 2 3 3 3 4 4 4 4 4   (n=12)
 *	        xsup 0 1 2 4 7 12
 *	Note: dfs will be performed on supernode rep. relative to the new 
 *	      row pivoting ordering
 *
 *   (xlsub,lsub): lsub[*] contains the compressed subscript of
 *	rectangular supernodes; xlsub[j] points to the starting
 *	location of the j-th column in lsub[*]. Note that xlsub 
 *	is indexed by column.
 *	Storage: original row subscripts
 *
 *      During the course of sparse LU factorization, we also use
 *	(xlsub,lsub) for the purpose of symmetric pruning. For each
 *	supernode {s,s+1,...,t=s+r} with first column s and last
 *	column t, the subscript set
 *		lsub[j], j=xlsub[s], .., xlsub[s+1]-1
 *	is the structure of column s (i.e. structure of this supernode).
 *	It is used for the storage of numerical values.
 *	Furthermore,
 *		lsub[j], j=xlsub[t], .., xlsub[t+1]-1
 *	is the structure of the last column t of this supernode.
 *	It is for the purpose of symmetric pruning. Therefore, the
 *	structural subscripts can be rearranged without making physical
 *	interchanges among the numerical values.
 *
 *	However, if the supernode has only one column, then we
 *	only keep one set of subscripts. For any subscript interchange
 *	performed, similar interchange must be done on the numerical
 *	values.
 *
 *	The last column structures (for pruning) will be removed
 *	after the numercial LU factorization phase.
 *
 *   (xlusup,lusup): lusup[*] contains the numerical values of the
 *	rectangular supernodes; xlusup[j] points to the starting
 *	location of the j-th column in storage vector lusup[*]
 *	Note: xlusup is indexed by column.
 *	Each rectangular supernode is stored by column-major
 *	scheme, consistent with Fortran 2-dim array storage.
 *
 *   (xusub,ucol,usub): ucol[*] stores the numerical values of
 *	U-columns outside the rectangular supernodes. The row
 *	subscript of nonzero ucol[k] is stored in usub[k].
 *	xusub[i] points to the starting location of column i in ucol.
 *	Storage: new row subscripts; that is subscripts of PA.
 */
typedef struct {
    int     *xsup;    /* supernode and column mapping */
    int     *supno;   
    int     *lsub;    /* compressed L subscripts */
    int	    *xlsub;
    float  *lusup;   /* L supernodes */
    int     *xlusup;
    float  *ucol;    /* U columns */
    int     *usub;
    int	    *xusub;
    int     nzlmax;   /* current max size of lsub */
    int     nzumax;   /*    "    "    "      ucol */
    int     nzlumax;  /*    "    "    "     lusup */
    int     n;        /* number of columns in the matrix */
    LU_space_t MemModel; /* 0 - system malloc'd; 1 - user provided */
} sGlobalLU_t;

typedef struct {
    int     *xsup;    /* supernode and column mapping */
    int     *supno;   
    int     *lsub;    /* compressed L subscripts */
    int	    *xlsub;
    double  *lusup;   /* L supernodes */
    int     *xlusup;
    double  *ucol;    /* U columns */
    int     *usub;
    int	    *xusub;
    int     nzlmax;   /* current max size of lsub */
    int     nzumax;   /*    "    "    "      ucol */
    int     nzlumax;  /*    "    "    "     lusup */
    int     n;        /* number of columns in the matrix */
    LU_space_t MemModel; /* 0 - system malloc'd; 1 - user provided */
} dGlobalLU_t;

typedef struct {
    int     *xsup;    /* supernode and column mapping */
    int     *supno;   
    int     *lsub;    /* compressed L subscripts */
    int	    *xlsub;
    lscomplex  *lusup;   /* L supernodes */
    int     *xlusup;
    lscomplex  *ucol;    /* U columns */
    int     *usub;
    int	    *xusub;
    int     nzlmax;   /* current max size of lsub */
    int     nzumax;   /*    "    "    "      ucol */
    int     nzlumax;  /*    "    "    "     lusup */
    int     n;        /* number of columns in the matrix */
    LU_space_t MemModel; /* 0 - system malloc'd; 1 - user provided */
} cGlobalLU_t;

typedef struct {
    int     *xsup;    /* supernode and column mapping */
    int     *supno;   
    int     *lsub;    /* compressed L subscripts */
    int	    *xlsub;
    ldcomplex  *lusup;   /* L supernodes */
    int     *xlusup;
    ldcomplex  *ucol;    /* U columns */
    int     *usub;
    int	    *xusub;
    int     nzlmax;   /* current max size of lsub */
    int     nzumax;   /*    "    "    "      ucol */
    int     nzlumax;  /*    "    "    "     lusup */
    int     n;        /* number of columns in the matrix */
    LU_space_t MemModel; /* 0 - system malloc'd; 1 - user provided */
} zGlobalLU_t;

typedef struct {
    int panel_size;
    int relax;
    float diag_pivot_thresh;
    float drop_tol;
} sfactor_param_t;

typedef struct {
    int panel_size;
    int relax;
    double diag_pivot_thresh;
    double drop_tol;
} dfactor_param_t;

typedef struct {
    float for_lu;
    float total_needed;
    int   expansions;
} mem_usage_t;

#ifdef __cplusplus
extern "C" {
#endif

/* Driver routines */
extern void
sgssv(SuperMatrix *, int *, int *, SuperMatrix *, SuperMatrix *, 
	SuperMatrix *, int *);
extern void
dgssv(SuperMatrix *, int *, int *, SuperMatrix *, SuperMatrix *, 
	SuperMatrix *, int *);
extern void
cgssv(SuperMatrix *, int *, int *, SuperMatrix *, SuperMatrix *, 
	SuperMatrix *, int *);
extern void
zgssv(SuperMatrix *, int *, int *, SuperMatrix *, SuperMatrix *, 
	SuperMatrix *, int *);
extern void
sgssvx(char *, char *, char *, SuperMatrix *, sfactor_param_t *,
       int *, int *, int *, char *, float *, float *,
       SuperMatrix *, SuperMatrix *, void *, int, SuperMatrix *, 
       SuperMatrix *, float *, float *, float *,
       float *, mem_usage_t *, int *);
extern void
dgssvx(char *, char *, char *, SuperMatrix *, dfactor_param_t *,
       int *, int *, int *, char *, double *, double *,
       SuperMatrix *, SuperMatrix *, void *, int, SuperMatrix *, 
       SuperMatrix *, double *, double *, double *,
       double *, mem_usage_t *, int *);
extern void
cgssvx(char *, char *, char *, SuperMatrix *, sfactor_param_t *,
       int *, int *, int *, char *, float *, float *,
       SuperMatrix *, SuperMatrix *, void *, int, SuperMatrix *, 
       SuperMatrix *, float *, float *, float *,
       float *, mem_usage_t *, int *);
extern void
zgssvx(char *, char *, char *, SuperMatrix *, dfactor_param_t *,
       int *, int *, int *, char *, double *, double *,
       SuperMatrix *, SuperMatrix *, void *, int, SuperMatrix *, 
       SuperMatrix *, double *, double *, double *,
       double *, mem_usage_t *, int *);

/* Supernodal LU factor related */
extern void
sCreate_CompCol_Matrix(SuperMatrix *, int, int, int, float *,
		       int *, int *, Stype_t, Dtype_t, Mtype_t);
extern void
dCreate_CompCol_Matrix(SuperMatrix *, int, int, int, double *,
		       int *, int *, Stype_t, Dtype_t, Mtype_t);
extern void
cCreate_CompCol_Matrix(SuperMatrix *, int, int, int, lscomplex *,
		       int *, int *, Stype_t, Dtype_t, Mtype_t);
extern void
zCreate_CompCol_Matrix(SuperMatrix *, int, int, int, ldcomplex *,
		       int *, int *, Stype_t, Dtype_t, Mtype_t);
extern void
sCopy_CompCol_Matrix(SuperMatrix *, SuperMatrix *);
extern void
dCopy_CompCol_Matrix(SuperMatrix *, SuperMatrix *);
extern void
cCopy_CompCol_Matrix(SuperMatrix *, SuperMatrix *);
extern void
zCopy_CompCol_Matrix(SuperMatrix *, SuperMatrix *);
extern void
sCreate_Dense_Matrix(SuperMatrix *, int, int, float *, int,
		     Stype_t, Dtype_t, Mtype_t);
extern void
dCreate_Dense_Matrix(SuperMatrix *, int, int, double *, int,
		     Stype_t, Dtype_t, Mtype_t);
extern void
cCreate_Dense_Matrix(SuperMatrix *, int, int, lscomplex *, int,
		     Stype_t, Dtype_t, Mtype_t);
extern void
zCreate_Dense_Matrix(SuperMatrix *, int, int, ldcomplex *, int,
		     Stype_t, Dtype_t, Mtype_t);
extern void
sCreate_SuperNode_Matrix(SuperMatrix *, int, int, int, float *, 
		         int *, int *, int *, int *, int *,
			 Stype_t, Dtype_t, Mtype_t);
extern void
dCreate_SuperNode_Matrix(SuperMatrix *, int, int, int, double *, 
		         int *, int *, int *, int *, int *,
			 Stype_t, Dtype_t, Mtype_t);
extern void
cCreate_SuperNode_Matrix(SuperMatrix *, int, int, int, lscomplex *, 
		         int *, int *, int *, int *, int *,
			 Stype_t, Dtype_t, Mtype_t);
extern void
zCreate_SuperNode_Matrix(SuperMatrix *, int, int, int, ldcomplex *, 
		         int *, int *, int *, int *, int *,
			 Stype_t, Dtype_t, Mtype_t);
extern void
sCopy_Dense_Matrix(int, int, float *, int, float *, int);
extern void
dCopy_Dense_Matrix(int, int, double *, int, double *, int);
extern void
cCopy_Dense_Matrix(int, int, lscomplex *, int, lscomplex *, int);
extern void
zCopy_Dense_Matrix(int, int, ldcomplex *, int, ldcomplex *, int);

extern void    Destroy_CompCol_Matrix(SuperMatrix *);
extern void    Destroy_SuperNode_Matrix(SuperMatrix *);
extern void    Destroy_CompCol_Permuted(SuperMatrix *);

extern void    sallocateA (int, int, float **, int **, int **);
extern void    dallocateA (int, int, double **, int **, int **);
extern void    callocateA (int, int, lscomplex **, int **, int **);
extern void    zallocateA (int, int, ldcomplex **, int **, int **);
extern void    StatInit    (int, int);
extern void    StatFree    ();
extern void    sp_preorder (char*, SuperMatrix*, int*, int*, SuperMatrix*);
extern int     sp_coletree (int *, int *, int *, int, int, int *);
extern void    sgstrf (char*, SuperMatrix*, float, float, int, int, int*,
			void *, int, int *, int *, 
                        SuperMatrix *, SuperMatrix *, int *);
extern void    dgstrf (char*, SuperMatrix*, double, double, int, int, int*,
			void *, int, int *, int *, 
                        SuperMatrix *, SuperMatrix *, int *);
extern void    cgstrf (char*, SuperMatrix*, float, float, int, int, int*,
			void *, int, int *, int *, 
                        SuperMatrix *, SuperMatrix *, int *);
extern void    zgstrf (char*, SuperMatrix*, double, double, int, int, int*,
			void *, int, int *, int *, 
                        SuperMatrix *, SuperMatrix *, int *);
extern void    relax_snode  (int, int *, int, int *, int *);
extern int     ssnode_dfs (const int, const int, const int *, const int *,
			     const int *, int *, int *, sGlobalLU_t *);
extern int     dsnode_dfs (const int, const int, const int *, const int *,
			     const int *, int *, int *, dGlobalLU_t *);
extern int     csnode_dfs (const int, const int, const int *, const int *,
			     const int *, int *, int *, cGlobalLU_t *);
extern int     zsnode_dfs (const int, const int, const int *, const int *,
			     const int *, int *, int *, zGlobalLU_t *);
extern int     ssnode_bmod (const int, const int, const int, float *,
                              float *, sGlobalLU_t *);
extern int     dsnode_bmod (const int, const int, const int, double *,
                              double *, dGlobalLU_t *);
extern int     csnode_bmod (const int, const int, const int, lscomplex *,
                              lscomplex *, cGlobalLU_t *);
extern int     zsnode_bmod (const int, const int, const int, ldcomplex *,
                              ldcomplex *, zGlobalLU_t *);
extern void    spanel_dfs (const int, const int, const int, SuperMatrix *,
			   int *, int *, float *, int *, int *, int *,
			   int *, int *, int *, int *, sGlobalLU_t *);
extern void    dpanel_dfs (const int, const int, const int, SuperMatrix *,
			   int *, int *, double *, int *, int *, int *,
			   int *, int *, int *, int *, dGlobalLU_t *);
extern void    cpanel_dfs (const int, const int, const int, SuperMatrix *,
			   int *, int *, lscomplex *, int *, int *, int *,
			   int *, int *, int *, int *, cGlobalLU_t *);
extern void    zpanel_dfs (const int, const int, const int, SuperMatrix *,
			   int *, int *, ldcomplex *, int *, int *, int *,
			   int *, int *, int *, int *, zGlobalLU_t *);
extern void    spanel_bmod (const int, const int, const int, const int,
                           float *, float *, int *, int *,
			   sGlobalLU_t *);
extern void    dpanel_bmod (const int, const int, const int, const int,
                           double *, double *, int *, int *,
			   dGlobalLU_t *);
extern void    cpanel_bmod (const int, const int, const int, const int,
                           lscomplex *, lscomplex *, int *, int *,
			   cGlobalLU_t *);
extern void    zpanel_bmod (const int, const int, const int, const int,
                           ldcomplex *, ldcomplex *, int *, int *,
			   zGlobalLU_t *);
extern int     scolumn_dfs (const int, const int, int *, int *, int *, int *,
			   int *, int *, int *, int *, int *, sGlobalLU_t *);
extern int     dcolumn_dfs (const int, const int, int *, int *, int *, int *,
			   int *, int *, int *, int *, int *, dGlobalLU_t *);
extern int     ccolumn_dfs (const int, const int, int *, int *, int *, int *,
			   int *, int *, int *, int *, int *, cGlobalLU_t *);
extern int     zcolumn_dfs (const int, const int, int *, int *, int *, int *,
			   int *, int *, int *, int *, int *, zGlobalLU_t *);
extern int     scolumn_bmod (const int, const int, float *,
			   float *, int *, int *, int, sGlobalLU_t *);
extern int     dcolumn_bmod (const int, const int, double *,
			   double *, int *, int *, int, dGlobalLU_t *);
extern int     ccolumn_bmod (const int, const int, lscomplex *,
			   lscomplex *, int *, int *, int, cGlobalLU_t *);
extern int     zcolumn_bmod (const int, const int, ldcomplex *,
			   ldcomplex *, int *, int *, int, zGlobalLU_t *);
extern int     scopy_to_ucol (int, int, int *, int *, int *,
                              float *, sGlobalLU_t *);         
extern int     dcopy_to_ucol (int, int, int *, int *, int *,
                              double *, dGlobalLU_t *);         
extern int     ccopy_to_ucol (int, int, int *, int *, int *,
                              lscomplex *, cGlobalLU_t *);         
extern int     zcopy_to_ucol (int, int, int *, int *, int *,
                              ldcomplex *, zGlobalLU_t *);         
extern int     spivotL (const int, const float, int *, int *, 
                              int *, int *, int *, sGlobalLU_t *);
extern int     dpivotL (const int, const double, int *, int *, 
                              int *, int *, int *, dGlobalLU_t *);
extern int     cpivotL (const int, const float, int *, int *, 
                              int *, int *, int *, cGlobalLU_t *);
extern int     zpivotL (const int, const double, int *, int *, 
                              int *, int *, int *, zGlobalLU_t *);
extern void    spruneL (const int, const int *, const int, const int,
			     const int *, const int *, int *, sGlobalLU_t *);
extern void    dpruneL (const int, const int *, const int, const int,
			     const int *, const int *, int *, dGlobalLU_t *);
extern void    cpruneL (const int, const int *, const int, const int,
			     const int *, const int *, int *, cGlobalLU_t *);
extern void    zpruneL (const int, const int *, const int, const int,
			     const int *, const int *, int *, zGlobalLU_t *);
extern void    resetrep_col (const int, const int *, int *);
//    extern void    countnz (const int, int *, int *, int *, sGlobalLU_t *);
//    extern void    fixupL (const int, const int *, sGlobalLU_t *);
extern int     spcoletree (int *, int *, int *, int, int, int *);
extern int     *TreePostorder (int, int *);
extern void    sreadmt (int *, int *, int *, float **, int **, int **);
extern void    dreadmt (int *, int *, int *, double **, int **, int **);
extern void    creadmt (int *, int *, int *, lscomplex **, int **, int **);
extern void    zreadmt (int *, int *, int *, ldcomplex **, int **, int **);
extern void    sGenXtrue (int, int, float *, int);
extern void    dGenXtrue (int, int, double *, int);
extern void    cGenXtrue (int, int, lscomplex *, int);
extern void    zGenXtrue (int, int, ldcomplex *, int);
extern void    sFillRHS (char *, int, float *, int, SuperMatrix *,
			SuperMatrix *);
extern void    dFillRHS (char *, int, double *, int, SuperMatrix *,
			SuperMatrix *);
extern void    cFillRHS (char *, int, lscomplex *, int, SuperMatrix *,
			SuperMatrix *);
extern void    zFillRHS (char *, int, ldcomplex *, int, SuperMatrix *,
			SuperMatrix *);
extern void    sgstrs (char *, SuperMatrix *, SuperMatrix *, int *, int *,
			SuperMatrix *, int *);
extern void    dgstrs (char *, SuperMatrix *, SuperMatrix *, int *, int *,
			SuperMatrix *, int *);
extern void    cgstrs (char *, SuperMatrix *, SuperMatrix *, int *, int *,
			SuperMatrix *, int *);
extern void    zgstrs (char *, SuperMatrix *, SuperMatrix *, int *, int *,
			SuperMatrix *, int *);


/* Driver related */

extern void    sgsequ (SuperMatrix *, float *, float *, float *,
			     float *, float *, int *);
extern void    dgsequ (SuperMatrix *, double *, double *, double *,
			     double *, double *, int *);
extern void    cgsequ (SuperMatrix *, float *, float *, float *,
			     float *, float *, int *);
extern void    zgsequ (SuperMatrix *, double *, double *, double *,
			     double *, double *, int *);
extern void    slaqgs (SuperMatrix *, float *, float *, float,
                             float, float, char *);
extern void    dlaqgs (SuperMatrix *, double *, double *, double,
                             double, double, char *);
extern void    claqgs (SuperMatrix *, float *, float *, float,
                             float, float, char *);
extern void    zlaqgs (SuperMatrix *, double *, double *, double,
                             double, double, char *);
extern void    sgscon (char *, SuperMatrix *, SuperMatrix *, 
			float, float *, int *);
extern void    dgscon (char *, SuperMatrix *, SuperMatrix *,
			double, double *, int *);
extern void    cgscon (char *, SuperMatrix *, SuperMatrix *, 
			float, float *, int *);
extern void    zgscon (char *, SuperMatrix *, SuperMatrix *,
			double, double *, int *);
extern void    sgsrfs (char *, SuperMatrix *, SuperMatrix *,
			SuperMatrix *, int *, int *, char *, float *, 
			float *, SuperMatrix *, SuperMatrix *, float *,
			float *, int *);
extern void    dgsrfs (char *, SuperMatrix *, SuperMatrix *, 
			SuperMatrix *, int *, int *, char *, double *,
			double *, SuperMatrix *, SuperMatrix *, 
			double *, double *, int *);
extern void    cgsrfs (char *, SuperMatrix *, SuperMatrix *,
			SuperMatrix *, int *, int *, char *, float *, 
			float *, SuperMatrix *, SuperMatrix *, float *,
			float *, int *);
extern void    zgsrfs (char *, SuperMatrix *, SuperMatrix *, 
			SuperMatrix *, int *, int *, char *, double *,
			double *, SuperMatrix *, SuperMatrix *, 
			double *, double *, int *);

extern int     sp_strsv (char *, char *, char *, SuperMatrix *,
			SuperMatrix *, float *, int *);
extern int     sp_dtrsv (char *, char *, char *, SuperMatrix *,
			SuperMatrix *, double *, int *);
extern int     sp_ctrsv (char *, char *, char *, SuperMatrix *,
			SuperMatrix *, lscomplex *, int *);
extern int     sp_ztrsv (char *, char *, char *, SuperMatrix *,
			SuperMatrix *, ldcomplex *, int *);
extern int     sp_sgemv (char *, float, SuperMatrix *, float *,
			int, float, float *, int);
extern int     sp_dgemv (char *, double, SuperMatrix *, double *,
			int, double, double *, int);
extern int     sp_cgemv (char *, lscomplex, SuperMatrix *, lscomplex *,
			int, lscomplex, lscomplex *, int);
extern int     sp_zgemv (char *, ldcomplex, SuperMatrix *, ldcomplex *,
			int, ldcomplex, ldcomplex *, int);

extern int     sp_sgemm (char *, char *, int, int, int, float,
			SuperMatrix *, float *, int, float, 
			float *, int);
extern int     sp_dgemm (char *, char *, int, int, int, double,
			SuperMatrix *, double *, int, double, 
			double *, int);
extern int     sp_cgemm (char *, char *, int, int, int, lscomplex,
			SuperMatrix *, lscomplex *, int, lscomplex, 
			lscomplex *, int);
extern int     sp_zgemm (char *, char *, int, int, int, ldcomplex,
			SuperMatrix *, ldcomplex *, int, ldcomplex, 
			ldcomplex *, int);

/* Memory-related */
extern void    *superlu_malloc (int);
extern void    superlu_free (void*);
extern int     sLUMemInit (char *, void *, int, int, int, int, int,
			     SuperMatrix *, SuperMatrix *,
			     sGlobalLU_t *, int **, float **);
extern int     dLUMemInit (char *, void *, int, int, int, int, int,
			     SuperMatrix *, SuperMatrix *,
			     dGlobalLU_t *, int **, double **);
extern int     cLUMemInit (char *, void *, int, int, int, int, int,
			     SuperMatrix *, SuperMatrix *,
			     cGlobalLU_t *, int **, lscomplex **);
extern int     zLUMemInit (char *, void *, int, int, int, int, int,
			     SuperMatrix *, SuperMatrix *,
			     zGlobalLU_t *, int **, ldcomplex **);
extern void    SetIWork (int, int, int, int *, int **, int **, int **,
                         int **, int **, int **, int **);
extern void    sSetRWork (int, int, float *, float **, float **);
extern void    dSetRWork (int, int, double *, double **, double **);
extern void    cSetRWork (int, int, lscomplex *, lscomplex **, lscomplex **);
extern void    zSetRWork (int, int, ldcomplex *, ldcomplex **, ldcomplex **);
extern void    sLUWorkFree (int *, float *, sGlobalLU_t *);
extern void    dLUWorkFree (int *, double *, dGlobalLU_t *);
extern void    cLUWorkFree (int *, lscomplex *, cGlobalLU_t *);
extern void    zLUWorkFree (int *, ldcomplex *, zGlobalLU_t *);
extern int     sLUMemXpand (int, int, MemType, int *, sGlobalLU_t *);
extern int     dLUMemXpand (int, int, MemType, int *, dGlobalLU_t *);
extern int     cLUMemXpand (int, int, MemType, int *, cGlobalLU_t *);
extern int     zLUMemXpand (int, int, MemType, int *, zGlobalLU_t *);

extern int     *intMalloc (int);
extern int     *intCalloc (int);
extern float  *floatMalloc(int);
extern double  *doubleMalloc(int);
extern lscomplex  *complexMalloc(int);
extern ldcomplex  *doublecomplexMalloc(int);
extern float  *floatCalloc(int);
extern double  *doubleCalloc(int);
extern lscomplex  *complexCalloc(int);
extern ldcomplex  *doublecomplexCalloc(int);
extern void    *superlu_malloc(int);
extern int     memory_usage();
extern int     sQuerySpace (SuperMatrix *, SuperMatrix *, int,
				mem_usage_t *);
extern int     dQuerySpace (SuperMatrix *, SuperMatrix *, int,
				mem_usage_t *);
extern int     cQuerySpace (SuperMatrix *, SuperMatrix *, int,
				mem_usage_t *);
extern int     zQuerySpace (SuperMatrix *, SuperMatrix *, int,
				mem_usage_t *);

/* Auxiliary routines */
extern double  SuperLU_timer_ ();
extern int     sp_ienv (int);
extern int     lsame_ (char *, char *);
extern int     xerbla_ (char *, int *);
extern void    ifill (int *, int, int);
extern void    sfill (float *, int, float);
extern void    dfill (double *, int, double);
extern void    cfill (lscomplex *, int, lscomplex);
extern void    zfill (ldcomplex *, int, ldcomplex);
//     extern void    inf_norm_error (int, SuperMatrix *, float *);
extern void    snode_profile (int, int *);
extern void    super_stats (int, int *);
extern void    PrintSumm (char *, int, int, int);
//     extern void    PrintPerf (SuperMatrix *, SuperMatrix *, mem_usage_t *,
//	         		 float, float, float *, float *, char *);

/* Routines for debugging */
//     extern void    print_lu_col(char *, int, int, int *, sGlobalLU_t *);
extern void    print_panel_seg(int, int, int, int, int *, int *);
//     extern void    check_tempv(int, float *);
extern void    check_repfnz(int, int, int, int *);

/* Reordering routine */

extern void    get_perm_c(int, SuperMatrix *, int *);  

#ifdef __cplusplus
  }
#endif

#endif /* __SUPERLU_SP_DEFS */

