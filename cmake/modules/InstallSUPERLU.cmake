SET(LIBNAME SUPERLU)

SET(URL http://crd-legacy.lbl.gov/~xiaoye/SuperLU/superlu_5.2.1.tar.gz)
SET(URL_MD5 3a1a9bff20cb06b7d97c46d337504447) 

SET(CMAKE_ARGS -D CMAKE_C_COMPILER=${CMAKE_C_COMPILER}
               -D CMAKE_C_FLAGS=${CMAKE_C_FLAGS}
               -D CMAKE_INSTALL_PREFIX=<INSTALL_DIR>
               -D BUILD_SHARED_LIBS=true)

SET(LIBRARIES libsuperlu.a)
