CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(FreeFem++ C CXX Fortran)

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)

INCLUDE(ProcessorCount)
ProcessorCount(N)
IF(NOT N EQUAL 0)
  SET(CMAKE_BUILD_FLAGS -j${N})
ENDIF()

######################
# C++ COMPILER OPTIONS
######################

SET(FF_CXX_FLAGS "-DBAMG_LONG_LONG")

IF(CMAKE_BUILD_TYPE EQUAL "Debug")
    SET(FF_CXX_FLAGS "${FF_CXX_FLAGS} -DCHECK_KN -g" )
    IF(CMAKE_SYSTEM_NAME EQUAL "Darwin")
        SET(FF_CXX_FLAGS "${FF_CXX_FLAGS} -fno-inline -fexceptions" )
    ENDIF(CMAKE_SYSTEM_NAME)
ELSE()
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNCHECKPTR -O3" )
ENDIF(CMAKE_BUILD_TYPE EQUAL "Debug")

FILE(GLOB_RECURSE BAMG_SRC src/bamg/*.cpp)
FILE(GLOB_RECURSE BAMGLIB_SRC src/bamglib/*.cpp)
FILE(GLOB FEMLIB_SRC src/femlib/*.cpp)
FILE(GLOB_RECURSE FFLIB_SRC src/fflib/*.cpp)
FILE(GLOB_RECURSE GRAPHICS_SRC src/Graphics//*.cpp)
FILE(GLOB_RECURSE LGLIB_SRC src/lglib/*.cpp)
FILE(GLOB_RECURSE MEDIT_SRC src/medit/*.c)
FILE(GLOB_RECURSE MPI_SRC src/mpi/*.cpp)
FILE(GLOB_RECURSE LIBMESH_SRC ${CMAKE_SOURCE_DIR}/src/libMesh/*.c)

###########
# DOWNLOAD
##########

INCLUDE(FFInstallPackage)

#FF_INSTALL_PACKAGE(MMG3D)
#FF_INSTALL_PACKAGE(Scalapack)
#FF_INSTALL_PACKAGE(Scotch)
#FF_INSTALL_PACKAGE(Tetgen)


#FILE(DOWNLOAD http://crd-legacy.lbl.gov/~xiaoye/SuperLU/superlu_4.3.tar.gz 
#              mydownload/superlu_4.3.tar.gz
#              SHOW_PROGRESS)

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


########
# CVMSH2
########
LIST(REMOVE_ITEM BAMGLIB_SRC ${CMAKE_SOURCE_DIR}/src/bamglib/Meshgibbs.cpp)
ADD_EXECUTABLE(cvmsh2 ${CMAKE_SOURCE_DIR}/src/bamg/cvmsh2.cpp src/bamg/global.cpp ${BAMGLIB_SRC})
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR} src/bamg src/bamglib src/fflib/ src/Graphics/ ${HDF5_INCUDE_DIRS})
TARGET_LINK_LIBRARIES(cvmsh2 libff hdf5 ${HDF5_LIBRARIES})


ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(examples++-load)
