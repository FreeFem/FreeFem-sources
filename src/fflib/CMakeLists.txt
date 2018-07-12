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

  FIND_PACKAGE(DLOPEN)
  IF(DLOPEN_FOUND)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D HAVE_DLFCN_H")
  ENDIF(DLOPEN_FOUND)

  INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/bamglib 
                    ${CMAKE_SOURCE_DIR}/src/fflib/ 
                    ${CMAKE_SOURCE_DIR}/src/Graphics/ 
                    ${CMAKE_SOURCE_DIR}/src/lglib/
                    ${CMAKE_SOURCE_DIR}/src/femlib)
  
  ADD_LIBRARY(libff ${LIBFF_SRC})

  ADD_LIBRARY(strversionnumber SHARED strversionnumber.cpp)
  TARGET_COMPILE_DEFINITIONS(strversionnumber PRIVATE VersionFreeFempp=${FREEFEM_VERSION})

  IF(FREEFEM_SUITESPARSE_INSTALLED)

    TARGET_COMPILE_DEFINITIONS(libff PRIVATE HAVE_LIBUMFPACK 
                                     PRIVATE HAVE_SUITESPARSE_UMFPACK_H)
    TARGET_LINK_LIBRARIES(libff bamglib dl strversionnumber FREEFEM::SUITESPARSE)

  ELSE(FREEFEM_SUITESPARSE_INSTALLED)

    TARGET_LINK_LIBRARIES(libff bamglib dl strversionnumber)
 
  ENDIF(FREEFEM_SUITESPARSE_INSTALLED) 


  SET_TARGET_PROPERTIES(libff PROPERTIES PREFIX "")

