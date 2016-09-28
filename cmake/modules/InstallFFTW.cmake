SET(LIBNAME FFTW)
IF(EXISTS ${${LIBNAME}_INSTALL_DIR})
    MESSAGE(STATUS "${LIBNAME} is declared to be installed in ${${LIBNAME}_INSTALL_DIR}")
ELSE()
    IF(ENABLE_DOWNLOAD)
        SET(URL http://www.fftw.org/fftw-3.3.2.tar.gz)
        SET(URL_MD5 6977ee770ed68c85698c7168ffa6e178)
        SET(CONFIGURE_COMMAND  <SOURCE_DIR>/configure --disable-dependency-tracking  --disable-fortran --prefix=<INSTALL_DIR> CXX=${CMAKE_CXX_COMPILER} CC=${CMAKE_CXX_COMPILER}  CFLAGS=${CMAKE_C_FLAGS} CXXFLAGS=${CMAKE_CXX_FLAGS})
    ELSE(ENABLE_DOWNLOAD)
        MESSAGE(FATAL_ERROR "${LIBNAME} not found. Please provide a value to ${LIBNAME}_DIR or set ENABLE_DOWNLOAD to True")
    ENDIF()
ENDIF()


