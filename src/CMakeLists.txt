###########
# FREEFEM++
###########
FIND_PACKAGE(HDF5 REQUIRED)
ADD_EXECUTABLE(FreeFem++ ${CMAKE_SOURCE_DIR}/src/Graphics/sansrgraph.cpp 
                         ${CMAKE_SOURCE_DIR}/src/mpi/parallelempi-empty.cpp 
                         ${CMAKE_SOURCE_DIR}/src/fflib/ffapi.cpp)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/lglib ${CMAKE_SOURCE_DIR}/src/fflib ${CMAKE_SOURCE_DIR}/src/Graphics)
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
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/src/mpi ${CMAKE_SOURCE_DIR}/src/femlib ${CMAKE_SOURCE_DIR}/src/fflib/ ${CMAKE_SOURCE_DIR}/src/Graphics/)
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
TARGET_LINK_LIBRARIES(FreeFem++-mpi ${MPI_CXX_LIBRARIES} libff umfpack ${HDF5_LIBRARIES} arpack blas)




ADD_SUBDIRECTORY(bamg)
ADD_SUBDIRECTORY(bamglib)
ADD_SUBDIRECTORY(fflib)
ADD_SUBDIRECTORY(Graphics)
ADD_SUBDIRECTORY(lglib)
ADD_SUBDIRECTORY(libMesh)
#ADD_SUBDIRECTORY(medit)
ADD_SUBDIRECTORY(mpi)
