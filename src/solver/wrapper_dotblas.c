typedef void doublecomplex;
typedef void complex;
typedef int integer; 
#include <vecLib/cblas.h> 

/* Double Complex */ void zdotc_(doublecomplex * ret_val, integer *n, 
	doublecomplex *zx, integer *incx, doublecomplex *zy, integer *incy)
{
    cblas_zdotc_sub(*n, zx, *incx, zy, *incx , ret_val);
} 

/* Double Complex */ void zdotu_(doublecomplex * ret_val, integer *n, 
	doublecomplex *zx, integer *incx, doublecomplex *zy, integer *incy)
{
    cblas_zdotu_sub(*n, zx, *incx, zy, *incx , ret_val);
    return ;
} /* zdotu___ */

/* Complex */ void cdotc_(complex * ret_val, integer *n, complex *cx, 
	integer *incx, complex *cy, integer *incy)
{
    cblas_cdotc_sub (*n, cx, *incx, cy, *incx , ret_val);
} /* cdotc___ */

/* Complex */ void  cdotu_(complex * ret_val, integer *n, complex *cx, 
	integer *incx, complex *cy, integer *incy)
{
  cblas_cdotu_sub(*n, cx, *incx, cy, *incx , ret_val);
} /* cdotu___ */

