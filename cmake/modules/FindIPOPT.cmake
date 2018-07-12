INCLUDE(FindPackageHandleStandardArgs)
INCLUDE(PackageManagerPaths)

FIND_PATH(IPOPT_INCLUDES NAMES ipTNLP.hpp 
                               IpIpoptApplication.hpp 
                         PATHS ${PACKMAN_INCLUDE_PATHS} 
                         PATH_SUFFIXES coin)

FIND_LIBRARY(IPOPT_LIBRARIES NAMES ipopt 
                             PATHS ${PACKMAN_LIBRARIES_PATHS})

IF(IPOPT_LIBRARIES)
  SET(IPOPT_FOUND True)
  MESSAGE(STATUS "A library with IPOPT API found.")
ENDIF(IPOPT_LIBRARIES)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(IPOPT DEFAULT_MSG IPOPT_INCLUDES IPOPT_LIBRARIES)



