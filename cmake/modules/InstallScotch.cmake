SET(LIBNAME SCOTCH)


IF(EXISTS ${${LIBNAME}_INSTALL_DIR})
    MESSAGE(STATUS "${LIBNAME} is declared to be installed in ${${LIBNAME}_INSTALL_DIR}")
ELSE()
    IF(ENABLE_DOWNLOAD)
        SET(${LIBNAME}_URL https://gforge.inria.fr/frs/download.php/file/34618/scotch_6.0.4.tar.gz)
        SET(${LIBNAME}_ARGS -D CMAKE_INSTALL_PREFIX=<INSTALL_DIR>)
        SET(${LIBNAME}_INSTALL_COMMAND cd .)
    ELSE(ENABLE_DOWNLOAD)
        MESSAGE(FATAL_ERROR "${LIBNAME} not found. Please provide a value to ${LIBNAME}_DIR or set ENABLE_DOWNLOAD to True")
    ENDIF()
ENDIF()


