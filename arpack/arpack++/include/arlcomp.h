/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE ARLComp.h.
   Unaltered copy of dcomplex.h and scomplex.h (from SuperLU package).
*/


#ifndef __SUPERLU_DCOMPLEX /* allow multiple inclusions */
#define __SUPERLU_DCOMPLEX

/* 
 * This header file is to be included in source files z*.c
 */
#ifndef DCOMPLEX_INCLUDE
#define DCOMPLEX_INCLUDE

typedef struct { double r, i; } ldcomplex;


/* Macro definitions */

/* Complex Addition c = a + b */
#define z_add(c, a, b) { (c)->r = (a)->r + (b)->r; \
			 (c)->i = (a)->i + (b)->i; }

/* Complex Subtraction c = a - b */
#define z_sub(c, a, b) { (c)->r = (a)->r - (b)->r; \
			 (c)->i = (a)->i - (b)->i; }

/* Complex-Double Multiplication */
#define zd_mult(c, a, b) { (c)->r = (a)->r * (b); \
                           (c)->i = (a)->i * (b); }

/* Complex-Complex Multiplication */
#define zz_mult(c, a, b) { \
	double cr, ci; \
    	cr = (a)->r * (b)->r - (a)->i * (b)->i; \
    	ci = (a)->i * (b)->r + (a)->r * (b)->i; \
    	(c)->r = cr; \
    	(c)->i = ci; \
    }

/* Complex equality testing */
#define z_eq(a, b)  ( (a)->r == (b)->r && (a)->i == (b)->i )


#ifdef __cplusplus
extern "C" {
#endif

/* Prototypes for functions in dcomplex.c */
void z_div(ldcomplex *, ldcomplex *, ldcomplex *);
double z_abs(ldcomplex *);     /* exact */
double z_abs1(ldcomplex *);    /* approximate */
void z_exp(ldcomplex *, ldcomplex *);
void d_cnjg(ldcomplex *r, ldcomplex *z);
double d_imag(ldcomplex *);


#ifdef __cplusplus
  }
#endif

#endif

#endif  /* __SUPERLU_DCOMPLEX */



#ifndef __SUPERLU_SCOMPLEX /* allow multiple inclusions */
#define __SUPERLU_SCOMPLEX

/* 
 * This header file is to be included in source files c*.c
 */
#ifndef SCOMPLEX_INCLUDE
#define SCOMPLEX_INCLUDE

typedef struct { float r, i; } lscomplex;


/* Macro definitions */

/* Complex Addition c = a + b */
#define c_add(c, a, b) { (c)->r = (a)->r + (b)->r; \
			 (c)->i = (a)->i + (b)->i; }

/* Complex Subtraction c = a - b */
#define c_sub(c, a, b) { (c)->r = (a)->r - (b)->r; \
			 (c)->i = (a)->i - (b)->i; }

/* Complex-Double Multiplication */
#define cs_mult(c, a, b) { (c)->r = (a)->r * (b); \
                           (c)->i = (a)->i * (b); }

/* Complex-Complex Multiplication */
#define cc_mult(c, a, b) { \
	float cr, ci; \
    	cr = (a)->r * (b)->r - (a)->i * (b)->i; \
    	ci = (a)->i * (b)->r + (a)->r * (b)->i; \
    	(c)->r = cr; \
    	(c)->i = ci; \
    }

/* Complex equality testing */
#define c_eq(a, b)  ( (a)->r == (b)->r && (a)->i == (b)->i )


#ifdef __cplusplus
extern "C" {
#endif

/* Prototypes for functions in scomplex.c */
void c_div(lscomplex *, lscomplex *, lscomplex *);
double c_abs(lscomplex *);     /* exact */
double c_abs1(lscomplex *);    /* approximate */
void c_exp(lscomplex *, lscomplex *);
void r_cnjg(lscomplex *, lscomplex *);
double r_imag(lscomplex *);


#ifdef __cplusplus
  }
#endif

#endif

#endif  /* __SUPERLU_SCOMPLEX */
