FILE(GLOB CPPFILES *.cpp)
LIST(REMOVE_ITEM CPPFILES ${CMAKE_SOURCE_DIR}/examples++-load/Element_Mixte.cpp                         # ok
                          ${CMAKE_SOURCE_DIR}/examples++-load/dfft.cpp                                  # ok
                          ${CMAKE_SOURCE_DIR}/examples++-load/buildlayer.cpp                            # depend de LayerMesh.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/pcm2rnm.cpp                               # depend de pcm.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/ffnewuoa.cpp                              # pas compris
                          ${CMAKE_SOURCE_DIR}/examples++-load/ff-cmaes.cpp                              # depend de cmaes.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/FreeFemQA.cpp                             # depend de GeometryQA.cpp
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


ADD_LIBRARY(buildlayer buildlayer.cpp)
TARGET_LINK_LIBRARIES(buildlayer LayerMesh)

ADD_LIBRARY(pcm2rnm pcm2rnm.cpp)
TARGET_LINK_LIBRARIES(pcm2rnm pcm)

ADD_LIBRARY(ff-cmaes ff-cmaes.cpp)
TARGET_LINK_LIBRARIES(ff-cmaes cmaes)

ADD_LIBRARY(FreeFemQA FreeFemQA.cpp)
TARGET_LINK_LIBRARIES(FreeFemQA GeometryQA)

# dfft
IF(FREEFEM_FFTW_INSTALLED)
  ADD_LIBRARY(dfft SHARED dfft.cpp)
  INCLUDE_DIRECTORIES(${FREEFEM_FFTW_INCLUDE_PATHS})
ENDIF(FREEFEM_FFTW_INSTALLED)

# Element_Mixte
IF(FREEFEM_LAPACK_INSTALLED)
  ADD_LIBRARY(Element_Mixte SHARED Element_Mixte.cpp)
  TARGET_LINK_LIBRARIES(Element_Mixte ${FREEFEM_LIBRARIES}
                                      ${FREEFEM_LAPACK_LIBRARY_PATHS})
ENDIF(FREEFEM_LAPACK_INSTALLED)

# ff-Ipopt.cpp
IF(FREEFEM_IPOPT_INSTALLED)
  ADD_LIBRARY(ff-Ipopt SHARED ff-Ipopt.cpp)
  ADD_DEPENDENCIES(ff-Ipopt ipopt)
  INCLUDE_DIRECTORIES(${FREEFEM_IPOPT_INCLUDE_PATHS})
  #INCLUDE_DIRECTORIES(${FF_DOWNLOAD_DIR}/ipopt/install/include/coin)
ENDIF(FREEFEM_IPOPT_INSTALLED)

# fflapack
IF(FREEFEM_LAPACK_INSTALLED)
	ADD_LIBRARY(fflapack SHARED fflapack.cpp)
	TARGET_LINK_LIBRARIES(fflapack ${FREEFEM_LIBRARIES}
                                 ${FREEFEM_LAPACK_LIBRARIES})
ENDIF(FREEFEM_LAPACK_INSTALLED)

# ffrandom
#IF(GSL_FOUND)
#  ADD_LIBRARY(ffrandom SHARED ffrandom.cpp)
#  TARGET_LINK_LIBRARIES(ffrandom ${GSL_LIBRARIES})
#  INCLUDE_DIRECTORIES(${GSL_INCLUDE_DIRS})
#ENDIF(GSL_FOUND)



# ff-NLpopt.cpp
IF(FREEFEM_NLOPT_INSTALLED)
  ADD_LIBRARY(ff-NLopt SHARED ff-NLopt.cpp)
  #ADD_DEPENDENCIES(ff-NLopt nlopt)
  INCLUDE_DIRECTORIES(${FREEFEM_NLOPT_INCLUDE_PATHS})
  MESSAGE(STATUS "my path is ${FREEFEM_NLOPT_INCLUDE_PATHS}") 
 #INCLUDE_DIRECTORIES(/home/cdoucet/Logiciels/modules_ff++/nlopt/install/include)
ENDIF(FREEFEM_NLOPT_INSTALLED)

# gsl
#IF(GSL_FOUND)
#  ADD_LIBRARY(gsl SHARED gsl.cpp)
#  TARGET_LINK_LIBRARIES(gsl ${GSL_LIBRARIES})
#  INCLUDE_DIRECTORIES(${GSL_INCLUDE_DIRS})
#ENDIF(GSL_FOUND)

# ilut.cpp
IF(FREEFEM_GMM3D_INSTALLED)
  ADD_LIBRARY(ilut SHARED ilut.cpp)
ENDIF(FREEFEM_GMM3D_INSTALLED)

	
# iohdf5
#IF(HDF5_FOUND)
#  ADD_LIBRARY(iohdf5 SHARED iohdf5.cpp)
#  TARGET_LINK_LIBRARIES(iohdf5 libff bamglib ${HDF5_LIBRARIES})
#  INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIRS})
#ENDIF(HDF5_FOUND)

# metis.cpp
IF(FREEFEM_METIS_INSTALLED)
  ADD_LIBRARY(ff-metis SHARED metis.cpp)
  INCLUDE_DIRECTORIES(${FREEFEM_METIS_INCLUDE_PATHS})
  TARGET_LINK_LIBRARIES(ff-metis libff lglib ${FREEFEM_METIS_LIBRARIES})
ENDIF(FREEFEM_METIS_INSTALLED)

# mshmet.cpp
#ADD_LIBRARY(mshmet mshmet.cpp)


# MUMPS_seq.cpp
#IF(MPI_FOUND AND MUMPS_FOUND)
#  ADD_LIBRARY(mumps-seq SHARED MUMPS_seq.cpp)
#  INCLUDE_DIRECTORIES(${MUMPS_INCLUDES} ${MPI_C_INCLUDE_PATH})
#  TARGET_LINK_LIBRARIES(mumps-seq lglib libff ${MPI_C_LIBRARIES})
#ENDIF(MPI_FOUND AND MUMPS_FOUND)

# MUMPS.cpp
#IF(MPI_FOUND AND MUMPS_FOUND)
#  ADD_LIBRARY(mumps SHARED MUMPS.cpp)
#  INCLUDE_DIRECTORIES(${MPI_C_INCLUDE_PATH})
#  TARGET_LINK_LIBRARIES(mumps lglib libff ${MPI_C_LIBRARIES})
#ENDIF(MPI_FOUND AND MUMPS_FOUND)

# mmg3d.cpp
IF(FREEFEM_MMG3D_INSTALLED)
  ADD_LIBRARY(ff-mmg3d SHARED mmg3d.cpp)
  ADD_DEPENDENCIES(ff-mmg3d mmg3d)
  INCLUDE_DIRECTORIES(${FF_DOWNLOAD_DIR}/mmg3d/install/include)
ENDIF(FREEFEM_MMG3D_INSTALLED)

# mmg3d-v4.0.cpp
#ADD_LIBRARY(ff-mmg3d-v4.0 SHARED mmg3d-v4.0.cpp)
#INCLUDE_DIRECTORIES(${FF_DOWNLOAD_DIR}/mmg3d/install/include)

# scotch.cpp
#IF(SCOTCH_FOUND)
#  ADD_LIBRARY(scotch SHARED scotch.cpp)
#  INCLUDE_DIRECTORIES(${SCOTCH_INCLUDES})
#  TARGET_LINK_LIBRARIES(scotch libff lglib ${SCOTCH_LIBRARIES})
#ENDIF()

#superlu.cpp
#IF(SUPERLU_FOUND)
#  ADD_LIBRARY(superlu SHARED SuperLu.cpp)
#  INCLUDE_DIRECTORIES(${SUPERLU_INCLUDES})
#  MESSAGE(STATUS "slu inc = ${SUPERLU_INCLUDES}")
#  TARGET_LINK_LIBRARIES(superlu ${SUPERLU_LIBRARIES})
#ENDIF(SUPERLU_FOUND)


# tetgen.cpp
IF(FREEFEM_TETGEN_INSTALLED)

  ADD_LIBRARY(example_tetgen SHARED tetgen.cpp
                             ${CMAKE_SOURCE_DIR}/src/femlib/BamgFreeFem.cpp
                             ${CMAKE_SOURCE_DIR}/src/femlib/splitsimplex.cpp)
  TARGET_LINK_LIBRARIES(example_tetgen ${FREEFEM_LIBRARIES}
                                       ${FREEFEM_TETGEN_LIBRARY_PATHS})
  INCLUDE_DIRECTORIES(${FREEFEM_TETGEN_INCLUDE_PATHS})

#ELSEIF(TETGEN_DOWNLOADED)

#  ADD_LIBRARY(example_tetgen SHARED tetgen.cpp
#                             ${CMAKE_SOURCE_DIR}/src/femlib/BamgFreeFem.cpp
#                             ${CMAKE_SOURCE_DIR}/src/femlib/splitsimplex.cpp)
#  TARGET_LINK_LIBRARIES(example_tetgen libMesh libff lglib ${FF_DOWNLOAD_DIR}/tetgen/install/libtet.a)
#  ADD_DEPENDENCIES(example_tetgen tetgen)
#  INCLUDE_DIRECTORIES(${FF_DOWNLOAD_DIR}/tetgen/src)

ENDIF(FREEFEM_TETGEN_INSTALLED)



# UMFPACK64.cpp
#ADD_LIBRARY(umfpack64 UMFPACK64.cpp)
#TARGET_LINK_LIBRARIES(umfpack64 /usr/lib/x86_64-linux-gnu/libumfpack.a)
#INCLUDE_DIRECTORIES(/usr/include/suitesparse)
