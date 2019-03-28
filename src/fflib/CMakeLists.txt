INCLUDE(ff_create_strversionnumber)

ff_create_strversionnumber()

FILE(GLOB FEMLIB_SRC ${CMAKE_SOURCE_DIR}/src/femlib/*.cpp)
FILE(GLOB_RECURSE FFLIB_SRC ${CMAKE_SOURCE_DIR}/src/fflib/*.cpp)

SET(LIBFF_SRC ${FEMLIB_SRC} ${FFLIB_SRC})

LIST(REMOVE_ITEM LIBFF_SRC ${CMAKE_SOURCE_DIR}/src/bamglib/Meshgibbs.cpp
                         ${CMAKE_SOURCE_DIR}/src/femlib/ConjuguedGradrientNL.cpp
                         ${CMAKE_SOURCE_DIR}/src/femlib/FESpace-v0.cpp
                         ${CMAKE_SOURCE_DIR}/src/femlib/glutdraw.cpp
                         ${CMAKE_SOURCE_DIR}/src/femlib/InvIntFunc.cpp
                         ${CMAKE_SOURCE_DIR}/src/femlib/mortar.cpp
                         ${CMAKE_SOURCE_DIR}/src/femlib/Pkorder.cpp
                         ${CMAKE_SOURCE_DIR}/src/femlib/P3korder.cpp
                         ${CMAKE_SOURCE_DIR}/src/fflib/ffapi.cpp
                         ${CMAKE_SOURCE_DIR}/src/fflib/strversionnumber.cpp)

LIST(APPEND LIBFF_SRC ${CMAKE_SOURCE_DIR}/src/Algo/lgalgo.cpp 
                    ${CMAKE_SOURCE_DIR}/src/Eigen/eigenvalue.cpp
                    ${CMAKE_SOURCE_DIR}/src/femlib/libmesh5.c
                    ${CMAKE_SOURCE_DIR}/src/Graphics/DefColor.cpp)


IF(DLOPEN_FOUND AND NOT ${CMAKE_SYSTEM_NAME} MATCHES "Windows" AND NOT MINGW)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D HAVE_DLFCN_H")
ENDIF(DLOPEN_FOUND AND NOT ${CMAKE_SYSTEM_NAME} MATCHES "Windows" AND NOT MINGW)


INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/bamglib 
                  ${CMAKE_SOURCE_DIR}/src/fflib/ 
                  ${CMAKE_SOURCE_DIR}/src/Graphics/ 
                  ${CMAKE_SOURCE_DIR}/src/lglib/
                  ${CMAKE_SOURCE_DIR}/src/femlib)

ADD_LIBRARY(libff ${LIBFF_SRC})

ADD_LIBRARY(strversionnumber STATIC strversionnumber.cpp)
TARGET_COMPILE_DEFINITIONS(strversionnumber PRIVATE VersionFreeFempp=${FREEFEM_VERSION})


IF(FF_SUITESPARSE_FOUND AND FF_CHOLMOD_FOUND)


  if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")

    TARGET_COMPILE_DEFINITIONS(libff PRIVATE HAVE_LIBUMFPACK 
                                     PRIVATE HAVE_SUITESPARSE_UMFPACK_H)

  else()

    TARGET_COMPILE_DEFINITIONS(libff PRIVATE HAVE_LIBUMFPACK 
                                     PRIVATE HAVE_UMFPACK_H)

  endif()

  if(MINGW)

    find_package(OpenMP)
    message(STATUS "OpenMP_C_LIBRARIES = ${OpenMP_C_LIBRARIES}")
    
    TARGET_LINK_LIBRARIES(libff bamglib Comdlg32 strversionnumber FREEFEM::SUITESPARSE
                                                                  FREEFEM::CHOLMOD 
                                                                  FREEFEM::AMD 
                                                                  FREEFEM::CAMD
                                                                  FREEFEM::CCOLAMD
                                                                  FREEFEM::COLAMD 
                                                                  FREEFEM::METIS
                                                                  FREEFEM::SUITESPARSECONFIG 
                                                                  ${FF_LAPACK_LIBRARIES}
                                                                  ${OpenMP_C_LIBRARIES})
  else()
    TARGET_LINK_LIBRARIES(libff bamglib dl strversionnumber FREEFEM::SUITESPARSE FREEFEM::CHOLMOD)
  endif(MINGW)

ELSE()

  if(MINGW)
    TARGET_LINK_LIBRARIES(libff bamglib Comdlg32 strversionnumber)
  else()
    TARGET_LINK_LIBRARIES(libff bamglib dl strversionnumber)
  endif(MINGW)

ENDIF(FF_SUITESPARSE_FOUND AND FF_CHOLMOD_FOUND) 


SET_TARGET_PROPERTIES(libff PROPERTIES PREFIX "")

