INCLUDE(FindPackageHandleStandardArgs)
INCLUDE(PackageManagerPaths)

FIND_PATH(MUMPS_INCLUDES NAMES smumps_c.h dmumps_c.h cmumps_c.h zmumps_c.h 
                          PATHS ${PACKMAN_INCLUDE_PATHS} 
                          PATH_SUFFIXES mumps)

FIND_LIBRARY(MUMPS_LIBRARIES NAMES mumps_common smumps dmumps cmumps zmumps 
                              PATHS ${PACKMAN_LIBRARIES_PATHS})

IF(MUMPS_INCLUDES AND MUMPS_LIBRARIES)
  SET(MUMPS_FOUND True)
  MESSAGE(STATUS "A library with SCOTCH API found.")
ENDIF(MUMPS_INCLUDES AND MUMPS_LIBRARIES)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(MUMPS DEFAULT_MSG MUMPS_INCLUDES MUMPS_LIBRARIES)



