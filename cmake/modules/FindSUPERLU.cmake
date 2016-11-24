INCLUDE(FindPackageHandleStandardArgs)

SET(SUPERLU_FOUND False)

IF(SUPERLU_INCLUDES AND SUPERLU_LIBRARIES)

  SET(SUPERLU_FOUND True)

ELSE(SUPERLU_INCLUDES AND SUPERLU_LIBRARIES)

  FIND_PATH(SUPERLU_INCLUDES NAMES scotch.h 
                                   colamd.h  
                                   slu_cdefs.h  
                                   slu_Cnames.h  
                                   slu_dcomplex.h  
                                   slu_ddefs.h  
                                   slu_scomplex.h  
                                   slu_sdefs.h  
                                   slu_util.h  
                                   slu_zdefs.h  
                                   superlu_enum_consts.h  
                                   supermatrix.h
                             PATHS /usr/include 
                             PATH_SUFFIXES superlu)

  FIND_LIBRARY(SUPERLU_LIBRARIES NAMES superlu 
                                 PATHS /usr/lib/x86_64-linux-gnu)

  FIND_PACKAGE_HANDLE_STANDARD_ARGS(SUPERLU DEFAULT_MSG SUPERLU_INCLUDES SUPERLU_LIBRARIES)

ENDIF(SUPERLU_INCLUDES AND SUPERLU_LIBRARIES)

