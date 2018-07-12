SET(LIBNAME NLOPT)

SET(URL http://ab-initio.mit.edu/nlopt/nlopt-2.2.4.tar.gz)
SET(URL_MD5 9c60c6380a10c6d2a06895f0e8756d4f)



SET(CONFIGURE_COMMAND ../src/configure CC=${CMAKE_C_COMPILER} CFLAGS=${CMAKE_C_FLAGS} --prefix=<INSTALL_DIR>)
SET(BUILD_COMMAND make)
SET(INSTALL_COMMAND make install)

SET(LIBRARIES libnlopt.a)

