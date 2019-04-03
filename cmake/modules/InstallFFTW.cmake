SET(LIBNAME FFTW)

SET(URL http://www.fftw.org/fftw-3.3.2.tar.gz)
SET(URL_MD5 6977ee770ed68c85698c7168ffa6e178)

SET(CONFIGURE_COMMAND  <SOURCE_DIR>/configure 
                       --disable-dependency-tracking  
                       --disable-fortran 
                       --prefix=<INSTALL_DIR> 
                       CXX=${CMAKE_CXX_COMPILER} 
                       CC=${CMAKE_C_COMPILER}  
                       CFLAGS=${CMAKE_C_FLAGS} 
                       CXXFLAGS=${CMAKE_CXX_FLAGS})

SET(LIBRARIES libfftw3.a)
