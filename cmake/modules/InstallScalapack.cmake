SET(LIBNAME SCALAPACK)


IF(EXISTS ${${LIBNAME}_INSTALL_DIR})
    MESSAGE(STATUS "${LIBNAME} is declared to be installed in ${${LIBNAME}_INSTALL_DIR}")
ELSE()
    IF(ENABLE_DOWNLOAD)
        SET(${LIBNAME}_URL http://www.netlib.org/scalapack/scalapack-2.0.2.tgz)
        SET(${LIBNAME}_ARGS -D CMAKE_INSTALL_PREFIX=<INSTALL_DIR>)
        SET(${LIBNAME}_INSTALL_COMMAND cd .)
    ELSE(ENABLE_DOWNLOAD)
        MESSAGE(FATAL_ERROR "${LIBNAME} not found. Please provide a value to ${LIBNAME}_DIR or set ENABLE_DOWNLOAD to True")
    ENDIF()
ENDIF()


