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

target_link_options(libff PUBLIC ${OPENMP_C_FLAGS} ${OPENMP_CXX_FLAGS})


IF(FREEFEM_SUITESPARSE_INSTALLED AND FREEFEM_CHOLMOD_INSTALLED)


if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")

  TARGET_COMPILE_DEFINITIONS(libff PRIVATE HAVE_LIBUMFPACK 
                                   PRIVATE HAVE_SUITESPARSE_UMFPACK_H)

else()

  TARGET_COMPILE_DEFINITIONS(libff PRIVATE HAVE_LIBUMFPACK 
                                   PRIVATE HAVE_UMFPACK_H)

endif()

  if(MINGW)
    message(STATUS "OpenMP_C_LIBRARIES = ${OpenMP_C_LIBRARIES}")
    message(STATUS "OpenMP_C_LIBRARY = ${OpenMP_C_LIBRARY}")
    message(STATUS "OpenMP_C_LIB_NAMES = ${OpenMP_C_LIB_NAMES}")
    message(STATUS "OpenMP_CXX_LIBRARIES = ${OpenMP_CXX_LIBRARIES}")
    message(STATUS "OpenMP_CXX_LIBRARY = ${OpenMP_CXX_LIBRARY}")
    message(STATUS "OpenMP_CXX_LIB_NAMES = ${OpenMP_CXX_LIB_NAMES}")

    foreach(WRONG_FILENAME ${OpenMP_C_LIBRARIES})
      string(REPLACE "C:/msys64" ""  RIGHT_FILENAME ${FILENAME})
      list(APPEND FF_OpenMP_C_LIBRARIES ${RIGHT_FILENAME})
    endforeach(WRONG_FILENAME ${OpenMP_C_LIBRARIES})

    foreach(WRONG_FILENAME ${OpenMP_CXX_LIBRARIES})
      string(REPLACE "C:/msys64" ""  RIGHT_FILENAME ${FILENAME})
      list(APPEND FF_OpenMP_CXX_LIBRARIES ${RIGHT_FILENAME})
    endforeach(WRONG_FILENAME ${OpenMP_CXX_LIBRARIES})

    message(STATUS "FF_OpenMP_C_LIBRARIES = ${FF_OpenMP_C_LIBRARIES}")
    message(STATUS "FF_OpenMP_CXX_LIBRARIES = ${FF_OpenMP_CXX_LIBRARIES}")
    
    TARGET_LINK_LIBRARIES(libff bamglib Comdlg32 strversionnumber ${FF_OpenMP_C_LIBRARIES} ${FF_OpenMP_CXX_LIBRARIES} ${LAPACK_LIBRARIES} FREEFEM::CHOLMOD FREEFEM::SUITESPARSE)
  else()
    TARGET_LINK_LIBRARIES(libff bamglib dl strversionnumber FREEFEM::CHOLMOD FREEFEM::SUITESPARSE)
  endif(MINGW)

ELSE(FREEFEM_SUITESPARSE_INSTALLED AND FREEFEM_CHOLMOD_INSTALLED)

  if(MINGW)
    TARGET_LINK_LIBRARIES(libff bamglib Comdlg32 strversionnumber)
  else()
    TARGET_LINK_LIBRARIES(libff bamglib dl strversionnumber)
  endif(MINGW)

ENDIF(FREEFEM_SUITESPARSE_INSTALLED AND FREEFEM_CHOLMOD_INSTALLED) 


SET_TARGET_PROPERTIES(libff PROPERTIES PREFIX "")

