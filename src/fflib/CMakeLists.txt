SET(LIBFF_SRC ${BAMGLIB_SRC} ${FEMLIB_SRC} ${FFLIB_SRC})
LIST(REMOVE_ITEM LIBFF_SRC ${CMAKE_SOURCE_DIR}/src/bamglib/Meshgibbs.cpp
                           ${CMAKE_SOURCE_DIR}/src/femlib/ConjuguedGradrientNL.cpp
                           ${CMAKE_SOURCE_DIR}/src/femlib/FESpace-v0.cpp
                           ${CMAKE_SOURCE_DIR}/src/femlib/glutdraw.cpp
                           ${CMAKE_SOURCE_DIR}/src/femlib/InvIntFunc.cpp
                           ${CMAKE_SOURCE_DIR}/src/femlib/mortar.cpp
                           ${CMAKE_SOURCE_DIR}/src/femlib/Pkorder.cpp
                           ${CMAKE_SOURCE_DIR}/src/femlib/P3korder.cpp
                           ${CMAKE_SOURCE_DIR}/src/fflib/ffapi.cpp)
LIST(APPEND LIBFF_SRC ${CMAKE_SOURCE_DIR}/src/Algo/lgalgo.cpp 
                      ${CMAKE_SOURCE_DIR}/src/Eigen/eigenvalue.cpp
                      ${CMAKE_SOURCE_DIR}/src/femlib/libmesh5.c
                      ${CMAKE_SOURCE_DIR}/src/Graphics/DefColor.cpp)

ADD_LIBRARY(libff ${LIBFF_SRC})
TARGET_LINK_LIBRARIES(libff dl cblas )
INCLUDE_DIRECTORIES(src/bamglib src/fflib/ src/Graphics/ src/lglib/)
SET_TARGET_PROPERTIES(libff PROPERTIES PREFIX "")

