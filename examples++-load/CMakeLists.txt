FILE(GLOB CPPFILES *.cpp)
LIST(REMOVE_ITEM CPPFILES ${CMAKE_SOURCE_DIR}/examples++-load/bmo_user.cpp                              # ff-cs
                          ${CMAKE_SOURCE_DIR}/examples++-load/Element_Mixte.cpp                         # ok
                          ${CMAKE_SOURCE_DIR}/examples++-load/dfft.cpp                                  # ok
                          ${CMAKE_SOURCE_DIR}/examples++-load/ff-Ipopt.cpp                              # Ipopt
                          ${CMAKE_SOURCE_DIR}/examples++-load/ff-NLopt.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/fflapack.cpp                              # ok
                          ${CMAKE_SOURCE_DIR}/examples++-load/freeyams.cpp                              # freeyams (probleme acces)
                          ${CMAKE_SOURCE_DIR}/examples++-load/ilut.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/iohdf5.cpp                                # ok
                          ${CMAKE_SOURCE_DIR}/examples++-load/lapack.cpp                                # ok
                          ${CMAKE_SOURCE_DIR}/examples++-load/lgbmo.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/metis.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/MetricPkuate.cpp                          # ne compile pas
                          ${CMAKE_SOURCE_DIR}/examples++-load/mmg3d.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/mmg3d-v4.0.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/mshmet.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/MUMPS.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/MUMPS_seq.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/netgen.cpp                                # gestion de netgen???
                          ${CMAKE_SOURCE_DIR}/examples++-load/NewSolver.cpp                             # umfpack
                          ${CMAKE_SOURCE_DIR}/examples++-load/PARDISO.cpp                               # gestion de pardiso???
                          ${CMAKE_SOURCE_DIR}/examples++-load/scotch.cpp                                # scotch
                          ${CMAKE_SOURCE_DIR}/examples++-load/SuperLu.cpp                               # superlu
                          ${CMAKE_SOURCE_DIR}/examples++-load/tetgen.cpp                                # ok
                          ${CMAKE_SOURCE_DIR}/examples++-load/UMFPACK64.cpp)                            # gestion de umfpack???

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/bamglib
                    ${CMAKE_SOURCE_DIR}/src/fflib
                    ${CMAKE_SOURCE_DIR}/src/femlib
                    ${CMAKE_SOURCE_DIR}/src/Graphics
                    ${CMAKE_SOURCE_DIR}/src/libMesh)

SET(FREEFEM_LIBRARIES libMesh libff lglib)

FOREACH(CPPFILE ${CPPFILES})
  GET_FILENAME_COMPONENT(FILENAME ${CPPFILE} NAME_WE)
  ADD_LIBRARY(${FILENAME} SHARED ${CPPFILE})
  TARGET_LINK_LIBRARIES(${FILENAME} ${FREEFEM_LIBRARIES})
ENDFOREACH(CPPFILE)

FF_INSTALL_PACKAGE(FFTW)
#FF_INSTALL_PACKAGE(FREEYAMS)         # probleme acces
#FF_INSTALL_PACKAGE(IPOPT)            # makefile
#FF_INSTALL_PACKAGE(MMG3D)            # depend de scotch
#FF_INSTALL_PACKAGE(OPENBLAS)         # probleme acces
#FF_INSTALL_PACKAGE(OPENBLAS-WIN32)   # binaire
#FF_INSTALL_PACKAGE(OPENBLAS-WIN64)   # binaire
#FF_INSTALL_PACKAGE(Scotch)           # makefile
#FF_INSTALL_PACKAGE(SUPERLU)          # makefile
FF_INSTALL_PACKAGE(Tetgen)

FIND_PACKAGE(GSL)
FIND_PACKAGE(HDF5)
FIND_PACKAGE(LAPACK)

# dfft
ADD_LIBRARY(dfft dfft.cpp)
INCLUDE_DIRECTORIES(${FF_DOWNLOAD_DIR}/fftw/install/include)


# Element_Mixte
IF(LAPACK_FOUND)
  ADD_LIBRARY(Element_Mixte SHARED Element_Mixte.cpp)
  TARGET_LINK_LIBRARIES(Element_Mixte ${FREEFEM_LIBRARIES}
                                      ${LAPACK_LIBRARIES})
ENDIF()

# fflapack
IF(LAPACK_FOUND)
	ADD_LIBRARY(fflapack SHARED fflapack.cpp)
	TARGET_LINK_LIBRARIES(fflapack ${FREEFEM_LIBRARIES}
                                 ${LAPACK_LIBRARIES})
ENDIF()
	
# lapack
IF(LAPACK_FOUND)
	ADD_LIBRARY(example-lapack SHARED lapack.cpp)
	TARGET_LINK_LIBRARIES(example-lapack ${FREEFEM_LIBRARIES}
                                       ${LAPACK_LIBRARIES})

ENDIF()

# iohdf5
IF(HDF5_FOUND)
  ADD_LIBRARY(iohdf5 SHARED iohdf5.cpp)
  TARGET_LINK_LIBRARIES(iohdf5 ${HDF5_LIBRARIES})
  INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIRS})
ENDIF()

# tetgen.cpp
ADD_LIBRARY(example_tetgen SHARED tetgen.cpp
                                  ${CMAKE_SOURCE_DIR}/src/femlib/BamgFreeFem.cpp
                                  ${CMAKE_SOURCE_DIR}/src/femlib/splitsimplex.cpp)
TARGET_LINK_LIBRARIES(example_tetgen libMesh libff lglib ${FF_DOWNLOAD_DIR}/tetgen/install/libtet.a)
ADD_DEPENDENCIES(example_tetgen tetgen)
INCLUDE_DIRECTORIES(${FF_DOWNLOAD_DIR}/tetgen/src)
