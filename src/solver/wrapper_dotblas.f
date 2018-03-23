c  a compile sans underscore   -fno-underscoring 
c  gfortran -fno-underscoring  -O3 -c wrapper_dotblas.f
c -------        
        double complex function zdotc_(n, zx, incx, zy, incy)
        double complex zx(*), zy(*), z
        integer n, incx, incy
        
        call cblas_zdotc_sub(%val(n), zx, %val(incx), zy, %val(incy), z)
c        print*,'cblas_zdotc_sub' 
        zdotc_ = z
        return
        end
        
        double complex function zdotu_(n, zx, incx, zy, incy)
        double complex zx(*), zy(*), z
        integer n, incx, incy
        
        call cblas_zdotu_sub(%val(n), zx, %val(incx), zy, %val(incy), z)
c        print*,'cblas_zdotu_sub' 
        zdotu_ = z
        return
        end
        
        complex function cdotc_(n, cx, incx, cy, incy)
        complex cx(*), cy(*), c
        integer n, incx, incy
        
        call cblas_cdotc_sub(%val(n), cx, %val(incx), cy, %val(incy), c)
c        print*,'cblas_cdotc_sub' 
        cdotc_ = c
        return
        end
        
        complex function cdotu_(n, cx, incx, cy, incy)
        complex cx(*), cy(*), c
        integer n, incx, incy
        
        call cblas_cdotu_sub(%val(n), cx, %val(incx), cy, %val(incy), c)
c        print*,'cblas_cdotu_sub' 
        cdotu_ = c
        return
        end
