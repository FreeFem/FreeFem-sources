SET(LIBNAME PETSC)
SET(URL http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.7.5.tar.gz)
SET(URL_MD5 bfc7a5535d5c18c6ec81ab90f3ce5074)
SET(CONFIGURE_COMMAND  cd <SOURCE_DIR> && ./configure 
                       --download-fftw
                       --download-hypre
                       --download-metis
                       --download-ml
                       --download-mpich
                       --download-mumps
                       --download-parmetis
                       --download-pastix
                       --download-ptscotch
                       --download-scalapack
                       --download-suitesparse
                       --download-superlu
                       PETSC-ARCH=arch-ff++ 
                       --prefix=<INSTALL_DIR> 
                       CXX=${CMAKE_CXX_COMPILER} 
                       CC=${CMAKE_C_COMPILER}  
                       CFLAGS=${CMAKE_C_FLAGS} 
                       CXXFLAGS=${CMAKE_CXX_FLAGS})
SET(BUILD_COMMAND cd <SOURCE_DIR> && make)
SET(INSTALL_COMMAND cd <SOURCE_DIR> && make install)

