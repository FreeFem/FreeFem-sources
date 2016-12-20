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

ADD_LIBRARY(libff ${LIBFF_SRC})

ADD_LIBRARY(strversionnumber SHARED strversionnumber.cpp)
TARGET_COMPILE_DEFINITIONS(strversionnumber PRIVATE VersionFreeFempp=${FREEFEM_VERSION})

TARGET_LINK_LIBRARIES(libff bamglib dl cblas strversionnumber)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/bamglib 
                    ${CMAKE_SOURCE_DIR}/src/fflib/ 
                    ${CMAKE_SOURCE_DIR}/src/Graphics/ 
                    ${CMAKE_SOURCE_DIR}/src/lglib/
                    ${CMAKE_SOURCE_DIR}/src/femlib)

SET_TARGET_PROPERTIES(libff PROPERTIES PREFIX "")


