c This file is part of FreeFem++.
c
c FreeFem++ is free software: you can redistribute it and/or modify
c it under the terms of the GNU Lesser General Public License as published by
c the Free Software Foundation, either version 3 of the License, or
c (at your option) any later version.
c
c FreeFem++ is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU Lesser General Public License for more details.
c
c You should have received a copy of the GNU Lesser General Public License
c along with Foobar.  If not, see <http://www.gnu.org/licenses/>.

c  a compile sans underscore   -fno-underscoring
c  gfortran -fno-underscoring  -O3 -c wrapper_dotblas.f
c -------
        double complex function zdotc_(n, zx, incx, zy, incy)
        double complex zx(*), zy(*), z
        integer n, incx, incy

        call cblas_zdotc_sub((n), zx, (incx), zy, (incy), z)
c        print*,'cblas_zdotc_sub'
        zdotc_ = z
        return
        end

        double complex function zdotu_(n, zx, incx, zy, incy)
        double complex zx(*), zy(*), z
        integer n, incx, incy

        call cblas_zdotu_sub((n), zx, (incx), zy, (incy), z)
c        print*,'cblas_zdotu_sub'
        zdotu_ = z
        return
        end

        complex function cdotc_(n, cx, incx, cy, incy)
        complex cx(*), cy(*), c
        integer n, incx, incy

        call cblas_cdotc_sub((n), cx, (incx), cy, (incy), c)
c        print*,'cblas_cdotc_sub'
        cdotc_ = c
        return
        end

        complex function cdotu_(n, cx, incx, cy, incy)
        complex cx(*), cy(*), c
        integer n, incx, incy

        call cblas_cdotu_sub((n), cx, (incx), cy, (incy), c)
c        print*,'cblas_cdotu_sub'
        cdotu_ = c
        return
        end
