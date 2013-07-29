/* wrapper_dotblas1.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/*  a compile sans underscore   -fno-underscoring */
/*  gfortran -fno-underscoring  -O3 -c wrapper_dotblas.f */
/* ------- */
/* Double Complex */ VOID zdotc___(doublecomplex * ret_val, integer *n, 
	doublecomplex *zx, integer *incx, doublecomplex *zy, integer *incy)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static doublecomplex z__;
    extern /* Subroutine */ int cblas_zdotc_sub__(integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, doublecomplex *);

    /* Parameter adjustments */
    --zy;
    --zx;

    /* Function Body */
    i__1 = *n;
    i__2 = *incx;
    i__3 = *incy;
    cblas_zdotc_sub__(&i__1, &zx[1], &i__2, &zy[1], &i__3, &z__);
/*        print*,'cblas_zdotc_sub' */
     ret_val->r = z__.r,  ret_val->i = z__.i;
    return ;
} /* zdotc___ */

/* Double Complex */ VOID zdotu___(doublecomplex * ret_val, integer *n, 
	doublecomplex *zx, integer *incx, doublecomplex *zy, integer *incy)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static doublecomplex z__;
    extern /* Subroutine */ int cblas_zdotu_sub__(integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, doublecomplex *);

    /* Parameter adjustments */
    --zy;
    --zx;

    /* Function Body */
    i__1 = *n;
    i__2 = *incx;
    i__3 = *incy;
    cblas_zdotu_sub__(&i__1, &zx[1], &i__2, &zy[1], &i__3, &z__);
/*        print*,'cblas_zdotu_sub' */
     ret_val->r = z__.r,  ret_val->i = z__.i;
    return ;
} /* zdotu___ */

/* Complex */ VOID cdotc___(complex * ret_val, integer *n, complex *cx, 
	integer *incx, complex *cy, integer *incy)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static complex c__;
    extern /* Subroutine */ int cblas_cdotc_sub__(integer *, complex *, 
	    integer *, complex *, integer *, complex *);

    /* Parameter adjustments */
    --cy;
    --cx;

    /* Function Body */
    i__1 = *n;
    i__2 = *incx;
    i__3 = *incy;
    cblas_cdotc_sub__(&i__1, &cx[1], &i__2, &cy[1], &i__3, &c__);
/*        print*,'cblas_cdotc_sub' */
     ret_val->r = c__.r,  ret_val->i = c__.i;
    return ;
} /* cdotc___ */

/* Complex */ VOID cdotu___(complex * ret_val, integer *n, complex *cx, 
	integer *incx, complex *cy, integer *incy)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static complex c__;
    extern /* Subroutine */ int cblas_cdotu_sub__(integer *, complex *, 
	    integer *, complex *, integer *, complex *);

    /* Parameter adjustments */
    --cy;
    --cx;

    /* Function Body */
    i__1 = *n;
    i__2 = *incx;
    i__3 = *incy;
    cblas_cdotu_sub__(&i__1, &cx[1], &i__2, &cy[1], &i__3, &c__);
/*        print*,'cblas_cdotu_sub' */
     ret_val->r = c__.r,  ret_val->i = c__.i;
    return ;
} /* cdotu___ */

