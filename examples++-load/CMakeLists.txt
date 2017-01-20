FILE(GLOB CPPFILES *.cpp)
LIST(REMOVE_ITEM CPPFILES ${CMAKE_SOURCE_DIR}/examples++-load/Element_Mixte.cpp                         # ok
                          ${CMAKE_SOURCE_DIR}/examples++-load/dfft.cpp                                  # ok
                          ${CMAKE_SOURCE_DIR}/examples++-load/ff-Ipopt.cpp                              # ok
                          ${CMAKE_SOURCE_DIR}/examples++-load/ff-NLopt.cpp                              # ok
                          ${CMAKE_SOURCE_DIR}/examples++-load/fflapack.cpp                              # ok
                          ${CMAKE_SOURCE_DIR}/examples++-load/ffrandom.cpp                              # depend de GSL
                          ${CMAKE_SOURCE_DIR}/examples++-load/freeyams.cpp                              # freeyams (probleme acces)
                          ${CMAKE_SOURCE_DIR}/examples++-load/gsl.cpp                                   # depend de gsl
                          ${CMAKE_SOURCE_DIR}/examples++-load/ilut.cpp                                  # ok (avec gmm++)
                          ${CMAKE_SOURCE_DIR}/examples++-load/iohdf5.cpp                                # ok
                          ${CMAKE_SOURCE_DIR}/examples++-load/metis.cpp                                 # ok
                          ${CMAKE_SOURCE_DIR}/examples++-load/mmg3d.cpp                                 # ok
                          ${CMAKE_SOURCE_DIR}/examples++-load/mmg3d-v4.0.cpp                            # dataff.h
                          ${CMAKE_SOURCE_DIR}/examples++-load/mshmet.cpp                                # mshmet (probleme acces)
                          ${CMAKE_SOURCE_DIR}/examples++-load/MUMPS.cpp                                 # ok
                          ${CMAKE_SOURCE_DIR}/examples++-load/MUMPS_seq.cpp                             # ok
                          ${CMAKE_SOURCE_DIR}/examples++-load/NewSolver.cpp                             # umfpack
                          ${CMAKE_SOURCE_DIR}/examples++-load/PARDISO.cpp                               # MKL
                          ${CMAKE_SOURCE_DIR}/examples++-load/scotch.cpp                                # ok
                          ${CMAKE_SOURCE_DIR}/examples++-load/SuperLu.cpp                               # ok
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

#FF_INSTALL_PACKAGE(FREEYAMS)         # probleme acces
#FF_INSTALL_PACKAGE(METIS)
#FF_INSTALL_PACKAGE(OPENBLAS)         # probleme acces
#FF_INSTALL_PACKAGE(OPENBLAS-WIN32)   # binaire
#FF_INSTALL_PACKAGE(OPENBLAS-WIN64)   # binaire
#FF_INSTALL_PACKAGE(Scotch)           # makefile
#FF_INSTALL_PACKAGE(SUPERLU)          # makefile

# etaient non commentes
#FF_INSTALL_PACKAGE(Tetgen)
#FF_INSTALL_PACKAGE(FFTW)
#FF_INSTALL_PACKAGE(MMG3D)            # depend de scotch
#FF_INSTALL_PACKAGE(IPOPT)            
#FF_INSTALL_PACKAGE(NLOPT)            # makefile



# dfft
IF(FFTW_DOWNLOADED)
  ADD_LIBRARY(dfft SHARED dfft.cpp)
  INCLUDE_DIRECTORIES(${FF_DOWNLOAD_DIR}/fftw/install/include)
ENDIF(FFTW_DOWNLOADED)

# Element_Mixte
IF(LAPACK_FOUND)
  ADD_LIBRARY(Element_Mixte SHARED Element_Mixte.cpp)
  TARGET_LINK_LIBRARIES(Element_Mixte ${FREEFEM_LIBRARIES}
                                      ${LAPACK_LIBRARIES})
ENDIF()

# ff-Ipopt.cpp
IF(IPOPT_DOWNLOADED)
  ADD_LIBRARY(ff-Ipopt SHARED ff-Ipopt.cpp)
  ADD_DEPENDENCIES(ff-Ipopt ipopt)
  INCLUDE_DIRECTORIES(${FF_DOWNLOAD_DIR}/ipopt/install/include/coin)
ENDIF(IPOPT_DOWNLOADED)

# fflapack
IF(LAPACK_FOUND)
	ADD_LIBRARY(fflapack SHARED fflapack.cpp)
	TARGET_LINK_LIBRARIES(fflapack ${FREEFEM_LIBRARIES}
                                 ${LAPACK_LIBRARIES})
ENDIF(LAPACK_FOUND)

# ffrandom
IF(GSL_FOUND)
  ADD_LIBRARY(ffrandom SHARED ffrandom.cpp)
  TARGET_LINK_LIBRARIES(ffrandom ${GSL_LIBRARIES})
  INCLUDE_DIRECTORIES(${GSL_INCLUDE_DIRS})
ENDIF(GSL_FOUND)



# ff-NLpopt.cpp
IF(NLOPT_DOWNLOADED)
  ADD_LIBRARY(ff-NLopt SHARED ff-NLopt.cpp)
  ADD_DEPENDENCIES(ff-NLopt nlopt)
  INCLUDE_DIRECTORIES(/home/cdoucet/Logiciels/modules_ff++/nlopt/install/include)
ENDIF(NLOPT_DOWNLOADED)

# gsl
IF(GSL_FOUND)
  ADD_LIBRARY(gsl SHARED gsl.cpp)
  TARGET_LINK_LIBRARIES(gsl ${GSL_LIBRARIES})
  INCLUDE_DIRECTORIES(${GSL_INCLUDE_DIRS})
ENDIF(GSL_FOUND)

# ilut.cpp
IF(GMM3D_FOUND)
  ADD_LIBRARY(ilut SHARED ilut.cpp)
ENDIF(GMM3D_FOUND)

	
# iohdf5
IF(HDF5_FOUND)
  ADD_LIBRARY(iohdf5 SHARED iohdf5.cpp)
  TARGET_LINK_LIBRARIES(iohdf5 ${HDF5_LIBRARIES})
  INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIRS})
ENDIF()

# metis.cpp
IF(METIS_FOUND)
  ADD_LIBRARY(metis SHARED metis.cpp)
ENDIF(METIS_FOUND)

# mshmet.cpp
#ADD_LIBRARY(mshmet mshmet.cpp)


# MUMPS_seq.cpp
IF(MPI_FOUND)
  ADD_LIBRARY(mumps-seq SHARED MUMPS_seq.cpp)
  INCLUDE_DIRECTORIES(${MPI_C_INCLUDE_PATH})
ENDIF(MPI_FOUND)

# MUMPS.cpp
IF(MPI_FOUND)
  ADD_LIBRARY(mumps SHARED MUMPS.cpp)
  INCLUDE_DIRECTORIES(${MPI_C_INCLUDE_PATH})
ENDIF(MPI_FOUND)

# mmg3d.cpp
IF(MMG3D_DOWNLOADED)
  ADD_LIBRARY(ff-mmg3d SHARED mmg3d.cpp)
  ADD_DEPENDENCIES(ff-mmg3d mmg3d)
  INCLUDE_DIRECTORIES(${FF_DOWNLOAD_DIR}/mmg3d/install/include)
ENDIF(MMG3D_DOWNLOADED)

# mmg3d-v4.0.cpp
#ADD_LIBRARY(ff-mmg3d-v4.0 SHARED mmg3d-v4.0.cpp)
#INCLUDE_DIRECTORIES(${FF_DOWNLOAD_DIR}/mmg3d/install/include)

# scotch.cpp
IF(SCOTCH_FOUND)
  ADD_LIBRARY(scotch SHARED scotch.cpp)
  INCLUDE_DIRECTORIES(/usr/include/scotch)
ENDIF()

#superlu.cpp
IF(SUPERLU_FOUND)
  ADD_LIBRARY(superlu SHARED SuperLu.cpp)
  INCLUDE_DIRECTORIES(${SUPERLU_INCLUDES})
ENDIF(SUPERLU_FOUND)


# tetgen.cpp
IF(TETGEN_FOUND)

  ADD_LIBRARY(example_tetgen SHARED tetgen.cpp
                             ${CMAKE_SOURCE_DIR}/src/femlib/BamgFreeFem.cpp
                             ${CMAKE_SOURCE_DIR}/src/femlib/splitsimplex.cpp)
  TARGET_LINK_LIBRARIES(example_tetgen libMesh libff lglib ${TETGEN_LIBRARIES})
  INCLUDE_DIRECTORIES(${TETGEN_INCLUDES})

ELSEIF(TETGEN_DOWNLOADED)

  ADD_LIBRARY(example_tetgen SHARED tetgen.cpp
                             ${CMAKE_SOURCE_DIR}/src/femlib/BamgFreeFem.cpp
                             ${CMAKE_SOURCE_DIR}/src/femlib/splitsimplex.cpp)
  TARGET_LINK_LIBRARIES(example_tetgen libMesh libff lglib ${FF_DOWNLOAD_DIR}/tetgen/install/libtet.a)
  ADD_DEPENDENCIES(example_tetgen tetgen)
  INCLUDE_DIRECTORIES(${FF_DOWNLOAD_DIR}/tetgen/src)

ENDIF(TETGEN_FOUND)



# UMFPACK64.cpp
#ADD_LIBRARY(umfpack64 UMFPACK64.cpp)
#TARGET_LINK_LIBRARIES(umfpack64 /usr/lib/x86_64-linux-gnu/libumfpack.a)
#INCLUDE_DIRECTORIES(/usr/include/suitesparse)
