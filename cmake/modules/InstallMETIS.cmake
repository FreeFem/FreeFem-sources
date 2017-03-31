SET(LIBNAME METIS)

SET(URL http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz)
SET(URL_MD5 5465e67079419a69e0116de24fce58fe)

SET(CMAKE_ARGS -D CMAKE_C_COMPILER=${CMAKE_C_COMPILER}
               -D CMAKE_C_FLAGS=${CMAKE_C_FLAGS}
               -D CMAKE_INSTALL_PREFIX=<INSTALL_DIR>
               -D SHARED=true
               -D GKLIB_PATH=../src/GKlib)

SET(LIBRARIES libmetis.so)

