CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(FreeFem++ C CXX Fortran)

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)

INCLUDE(ProcessorCount)
ProcessorCount(N)
IF(NOT N EQUAL 0)
  SET(CMAKE_BUILD_FLAGS -j${N})
ENDIF()

FILE(GLOB_RECURSE BAMG_SRC src/bamg/*.cpp)
FILE(GLOB_RECURSE BAMGLIB_SRC src/bamglib/*.cpp)
FILE(GLOB FEMLIB_SRC src/femlib/*.cpp)
FILE(GLOB_RECURSE FFLIB_SRC src/fflib/*.cpp)
FILE(GLOB_RECURSE GRAPHICS_SRC src/Graphics//*.cpp)
FILE(GLOB_RECURSE LGLIB_SRC src/lglib/*.cpp)
FILE(GLOB_RECURSE MEDIT_SRC src/medit/*.c)
FILE(GLOB_RECURSE MPI_SRC src/mpi/*.cpp)

###########
# DOWNLOAD
##########

INCLUDE(FFInstallPackage)

#FF_INSTALL_PACKAGE(MMG3D)
#FF_INSTALL_PACKAGE(Scalapack)
#FF_INSTALL_PACKAGE(Scotch)
#FF_INSTALL_PACKAGE(Tetgen)


FILE(DOWNLOAD http://crd-legacy.lbl.gov/~xiaoye/SuperLU/superlu_4.3.tar.gz 
              mydownload/superlu_4.3.tar.gz
              SHOW_PROGRESS)

#########
# LIBMESH
#########

FILE(GLOB_RECURSE LIBMESH_SRC ${CMAKE_SOURCE_DIR}/src/libMesh/*.c)
LIST(REMOVE_ITEM BAMGLIB_SRC ${CMAKE_SOURCE_DIR}/src/libMesh/test_sol.c
                             ${CMAKE_SOURCE_DIR}/src/libMesh/test_libmesh5.c)

ADD_LIBRARY(libMesh ${LIBMESH_SRC})
INCLUDE_DIRECTORIES(src/libMesh)
SET_TARGET_PROPERTIES(libMesh PROPERTIES PREFIX "")

#######
# LGLIB
#######

ADD_LIBRARY(lglib ${LGLIB_SRC})
INCLUDE_DIRECTORIES(src/bamglib src/femlib/ src/fflib/ src/Graphics/ src/lglib/)
SET_TARGET_PROPERTIES(lglib PROPERTIES PREFIX "")

#######
# LIBFF
#######

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

###########
# FREEFEM++
###########
FIND_PACKAGE(HDF5 REQUIRED)
ADD_EXECUTABLE(FreeFem++ ${CMAKE_SOURCE_DIR}/src/Graphics/sansrgraph.cpp 
                         ${CMAKE_SOURCE_DIR}/src/mpi/parallelempi-empty.cpp 
                         ${CMAKE_SOURCE_DIR}/src/fflib/ffapi.cpp)
INCLUDE_DIRECTORIES(src/lglib src/fflib src/Graphics)
TARGET_LINK_LIBRARIES(FreeFem++ arpack cblas ${HDF5_LIBRARIES} dl umfpack lglib libff)

########
# FFGLUT
########

ADD_EXECUTABLE(ffglut ${CMAKE_SOURCE_DIR}/src/Graphics/ffglut.cpp
                      ${CMAKE_SOURCE_DIR}/src/femlib/fem.cpp
                      ${CMAKE_SOURCE_DIR}/src/femlib/Mesh3dn.cpp
                      ${CMAKE_SOURCE_DIR}/src/femlib/Mesh2dn.cpp
                      ${CMAKE_SOURCE_DIR}/src/femlib/Mesh1dn.cpp
                      ${CMAKE_SOURCE_DIR}/src/femlib/GQuadTree.cpp
                      ${CMAKE_SOURCE_DIR}/src/femlib/FQuadTree.cpp
                      ${CMAKE_SOURCE_DIR}/src/femlib/Drawing.cpp
                      ${CMAKE_SOURCE_DIR}/src/femlib/mshptg.cpp
                      ${CMAKE_SOURCE_DIR}/src/Graphics/gggg.cpp
                      ${CMAKE_SOURCE_DIR}/src/fflib/ffapi.cpp
                      ${CMAKE_SOURCE_DIR}/src/Graphics/ffthreads.cpp)
FIND_PACKAGE(GLUT REQUIRED)
FIND_PACKAGE(OpenGL REQUIRED)
FIND_PACKAGE(Threads REQUIRED)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR} src/mpi src/femlib src/fflib/ src/Graphics/)
TARGET_LINK_LIBRARIES(ffglut lglib ${CMAKE_THREAD_LIBS_INIT} ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY} ${GLUT_LIBRARIES} libff)

###############
# FREEFEM++-MPI
###############

ADD_EXECUTABLE(FreeFem++-mpi ${CMAKE_SOURCE_DIR}/src/mpi/parallelempi.cpp
                             ${CMAKE_SOURCE_DIR}/src/Graphics/sansrgraph.cpp
                             ${CMAKE_SOURCE_DIR}/src/lglib/mymain.cpp
                             ${CMAKE_SOURCE_DIR}/src/lglib/lg.tab.cpp
                             ${CMAKE_SOURCE_DIR}/src/fflib/ffapi.cpp)
FIND_PACKAGE(MPI REQUIRED)
INCLUDE_DIRECTORIES(${MPI_CXX_INCLUDE_PATH})
TARGET_LINK_LIBRARIES(FreeFem++-mpi ${MPI_CXX_LIBRARIES} libff umfpack hdf5 arpack blas)

#########
# FFMEDIT
#########

ADD_EXECUTABLE(ffmedit ${MEDIT_SRC})
FILE(WRITE src/medit/compil.date "\#define COMPIL  \" jeudi 4 août 2016, 13:47:25 (UTC+0200) (with ff++ 3.47)\"")
INCLUDE_DIRECTORIES(src/libMesh src/medit)
TARGET_LINK_LIBRARIES(ffmedit libMesh m ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY} ${GLUT_LIBRARIES} )


######
# BAMG
######

LIST(REMOVE_ITEM BAMGLIB_SRC ${CMAKE_SOURCE_DIR}/src/bamglib/Meshgibbs.cpp)
ADD_EXECUTABLE(bamg src/bamg/bamg.cpp src/bamg/global.cpp ${BAMGLIB_SRC})
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR} src/bamg src/bamglib src/fflib/ src/Graphics/)
TARGET_LINK_LIBRARIES(bamg libff ${HDF5_LIBRARIES})
#${HDF5_LIBRARIES})

########
# CVMSH2
########

ADD_EXECUTABLE(cvmsh2 ${CMAKE_SOURCE_DIR}/src/bamg/cvmsh2.cpp src/bamg/global.cpp ${BAMGLIB_SRC})
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR} src/bamg src/bamglib src/fflib/ src/Graphics/ ${HDF5_INCUDE_DIRS})
TARGET_LINK_LIBRARIES(cvmsh2 libff hdf5 ${HDF5_LIBRARIES})

################
# EXAMPLE++-LOAD
################

#FILE(GLOB EXAMPLES_LOAD_SRC examples++-load/*.cpp)
#FOREACH(CPPFILE ${EXAMPLES_LOAD_SRC})
#  GET_FILENAME_COMPONENT(FILENAME ${CPPFILE} NAME_WE)
#  ADD_LIBRARY(${FILENAME} SHARED ${CPPFILE})
#  SET_TARGET_PROPERTIES(${FILENAME} PROPERTIES PREFIX "")
#ENDFOREACH(CPPFILE)

