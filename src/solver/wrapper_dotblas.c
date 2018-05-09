/*
 * This file is part of FreeFem++.
 *
 * FreeFem++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * FreeFem++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
 */

typedef void doublecomplex;
typedef void complex;
typedef int integer;
#include <vecLib/cblas.h>

/* Double Complex */ void zdotc_ (doublecomplex *ret_val, integer *n,
                                  doublecomplex *zx, integer *incx, doublecomplex *zy, integer *incy) {
	cblas_zdotc_sub(*n, zx, *incx, zy, *incx, ret_val);
}

/* Double Complex */ void zdotu_ (doublecomplex *ret_val, integer *n,
                                  doublecomplex *zx, integer *incx, doublecomplex *zy, integer *incy) {
	cblas_zdotu_sub(*n, zx, *incx, zy, *incx, ret_val);
	return;
}	/* zdotu___ */

/* Complex */ void cdotc_ (complex *ret_val, integer *n, complex *cx,
                           integer *incx, complex *cy, integer *incy) {
	cblas_cdotc_sub(*n, cx, *incx, cy, *incx, ret_val);
}	/* cdotc___ */

/* Complex */ void cdotu_ (complex *ret_val, integer *n, complex *cx,
                           integer *incx, complex *cy, integer *incy) {
	cblas_cdotu_sub(*n, cx, *incx, cy, *incx, ret_val);
}	/* cdotu___ */
