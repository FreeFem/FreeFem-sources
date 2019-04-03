SET(LIBNAME IPOPT)

SET(URL http://www.coin-or.org/download/source/Ipopt/Ipopt-3.12.4.tgz)
SET(URL_MD5 12a8ecaff8dd90025ddea6c65b49cb03)

SET(CONFIGURE_COMMAND ../src/configure --prefix=<INSTALL_DIR>)
SET(BUILD_COMMAND make)
SET(INSTALL_COMMAND make install)

SET(INCLUDE_PATHS ${FF_DOWNLOAD_DIR}/include/coin)
SET(LIBRARIES libipopt.so)
