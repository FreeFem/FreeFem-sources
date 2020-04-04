set(CMAKE_INSTALL_RPATH "${CMAKE_BINARY_DIR}/examples++-mpi ${CMAKE_BINARY_DIR}/examples++-load")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

FILE(GLOB EDP_FILES *.edp)

FILE(GLOB CPPFILES *.cpp)
LIST(REMOVE_ITEM CPPFILES ${CMAKE_SOURCE_DIR}/examples++-load/Element_Mixte.cpp                         # ok
                          ${CMAKE_SOURCE_DIR}/examples++-load/dfft.cpp                                  # ok
                          ${CMAKE_SOURCE_DIR}/examples++-load/ffnewuoa.cpp                              # pas compris
                          ${CMAKE_SOURCE_DIR}/examples++-load/ff-Ipopt.cpp                              # ok
                          ${CMAKE_SOURCE_DIR}/examples++-load/ff-NLopt.cpp                              # ok
                          ${CMAKE_SOURCE_DIR}/examples++-load/fflapack.cpp                              # ok
                          ${CMAKE_SOURCE_DIR}/examples++-load/ffrandom.cpp                              # depend de GSL
                          ${CMAKE_SOURCE_DIR}/examples++-load/freeyams.cpp                              # freeyams (probleme acces)
                          ${CMAKE_SOURCE_DIR}/examples++-load/gsl.cpp                                   # depend de gsl
                          ${CMAKE_SOURCE_DIR}/examples++-load/ilut.cpp                                  # ok (avec gmm++)
                          ${CMAKE_SOURCE_DIR}/examples++-load/iohdf5.cpp                                # ok
                          ${CMAKE_SOURCE_DIR}/examples++-load/metis.cpp                                 # ok
                          ${CMAKE_SOURCE_DIR}/examples++-load/mmg3d-v4.0.cpp                            # dataff.h
                          ${CMAKE_SOURCE_DIR}/examples++-load/mshmet.cpp                                # mshmet (probleme acces)
                          ${CMAKE_SOURCE_DIR}/examples++-load/MUMPS.cpp                                 # ok
                          ${CMAKE_SOURCE_DIR}/examples++-load/MUMPS_seq.cpp                             # ok
                          ${CMAKE_SOURCE_DIR}/examples++-load/netgen.cpp                                # probleme nglib.h introuvable
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

LIST(REMOVE_ITEM CPPFILES ${CMAKE_SOURCE_DIR}/examples++-load/Curvature.cpp)
ADD_LIBRARY(Curvature MODULE Curvature.cpp)
TARGET_LINK_LIBRARIES(Curvature libMesh)



FOREACH(CPPFILE ${CPPFILES})
  GET_FILENAME_COMPONENT(FILENAME ${CPPFILE} NAME_WE)
  ADD_LIBRARY(${FILENAME} MODULE ${CPPFILE})
  SET_TARGET_PROPERTIES(${FILENAME} PROPERTIES PREFIX "")
ENDFOREACH(CPPFILE)

# dfft
IF(FREEFEM_FFTW_INSTALLED)
  ADD_LIBRARY(dfft MODULE dfft.cpp)
  TARGET_LINK_LIBRARIES(dfft FREEFEM::FFTW)
  ADD_DEPENDENCIES(dfft FREEFEM::FFTW)
  SET_TARGET_PROPERTIES(dfft PROPERTIES PREFIX "")
ENDIF(FREEFEM_FFTW_INSTALLED)

# Element_Mixte
IF(FREEFEM_LAPACK_INSTALLED)
  ADD_LIBRARY(Element_Mixte MODULE Element_Mixte.cpp)
  TARGET_LINK_LIBRARIES(Element_Mixte ${FREEFEM_LAPACK_LIBRARY_PATHS})
  SET_TARGET_PROPERTIES(Element_Mixte PROPERTIES PREFIX "")
ELSE()
  LIST(REMOVE_ITEM EDP_FILES ${CMAKE_SOURCE_DIR}/examples++-load/lame-TD-NSS.edp)
  LIST(REMOVE_ITEM EDP_FILES ${CMAKE_SOURCE_DIR}/examples++-load/LaplaceRT1.edp)
  LIST(REMOVE_ITEM EDP_FILES ${CMAKE_SOURCE_DIR}/examples++-load/LaplaceRT2.edp)
  LIST(REMOVE_ITEM EDP_FILES ${CMAKE_SOURCE_DIR}/examples++-load/test-ElementMixte.edp)
ENDIF(FREEFEM_LAPACK_INSTALLED)

# ff-Ipopt.cpp
IF(FREEFEM_IPOPT_INSTALLED)
  ADD_LIBRARY(ff-Ipopt MODULE ff-Ipopt.cpp)
  TARGET_LINK_LIBRARIES(ff-Ipopt FREEFEM::IPOPT)
  SET_TARGET_PROPERTIES(ff-Ipopt PROPERTIES PREFIX "")
  TARGET_COMPILE_DEFINITIONS(ff-Ipopt PUBLIC "-DHAVE_CSTDDEF")
ENDIF(FREEFEM_IPOPT_INSTALLED)

# fflapack
IF(FREEFEM_LAPACK_INSTALLED)
	ADD_LIBRARY(fflapack MODULE fflapack.cpp)
	TARGET_LINK_LIBRARIES(fflapack ${FREEFEM_LAPACK_LIBRARIES})
  SET_TARGET_PROPERTIES(fflapack PROPERTIES PREFIX "")
ENDIF(FREEFEM_LAPACK_INSTALLED)

# ffrandom
IF(GSL_FOUND)
  ADD_LIBRARY(ffrandom SHARED ffrandom.cpp)
  TARGET_LINK_LIBRARIES(ffrandom ${GSL_LIBRARIES})
  INCLUDE_DIRECTORIES(${GSL_INCLUDE_DIRS})
  SET_TARGET_PROPERTIES(ffrandom PROPERTIES PREFIX "")
ENDIF(GSL_FOUND)



# ff-NLpopt.cpp
IF(FREEFEM_NLOPT_INSTALLED)
  ADD_LIBRARY(ff-NLopt MODULE ff-NLopt.cpp)
  TARGET_LINK_LIBRARIES(ff-NLopt FREEFEM::NLOPT)
  SET_TARGET_PROPERTIES(ff-NLopt PROPERTIES PREFIX "")
ENDIF(FREEFEM_NLOPT_INSTALLED)

# gsl
IF(GSL_FOUND)
  ADD_LIBRARY(gsl SHARED gsl.cpp)
  TARGET_LINK_LIBRARIES(gsl ${GSL_LIBRARIES})
  INCLUDE_DIRECTORIES(${GSL_INCLUDE_DIRS})
  SET_TARGET_PROPERTIES(gsl PROPERTIES PREFIX "")
ELSE(GSL_FOUND)
  LIST(REMOVE_ITEM EDP_FILES ${CMAKE_SOURCE_DIR}/examples++-load/gsl.edp)
ENDIF(GSL_FOUND)

# ilut.cpp
IF(FREEFEM_GMM3D_INSTALLED)
  ADD_LIBRARY(ilut MODULE ilut.cpp)
  SET_TARGET_PROPERTIES(ilut PROPERTIES PREFIX "")
ELSE()
  LIST(REMOVE_ITEM EDP_FILES ${CMAKE_SOURCE_DIR}/examples++-load/ilut.edp)
ENDIF(FREEFEM_GMM3D_INSTALLED)

	
# iohdf5
IF(HDF5_FOUND)
  ADD_LIBRARY(iohdf5 SHARED iohdf5.cpp)
  INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIRS})
  SET_TARGET_PROPERTIES(iohdf5 PROPERTIES PREFIX "")
ENDIF(HDF5_FOUND)

# metis.cpp
IF(FREEFEM_METIS_INSTALLED)
  ADD_LIBRARY(metis MODULE metis.cpp)
  TARGET_LINK_LIBRARIES(metis libff lglib FREEFEM::METIS)
  SET_TARGET_PROPERTIES(metis PROPERTIES PREFIX "")
ENDIF(FREEFEM_METIS_INSTALLED)

# mshmet.cpp
IF(FREEFEM_MSHMET_INSTALLED)
  ADD_LIBRARY(mshmet mshmet.cpp)
ENDIF()


# MUMPS_seq.cpp
IF(MPI_FOUND AND MUMPS_FOUND)
  ADD_LIBRARY(mumps-seq SHARED MUMPS_seq.cpp)
  INCLUDE_DIRECTORIES(${MUMPS_INCLUDES} ${MPI_C_INCLUDE_PATH})
  TARGET_LINK_LIBRARIES(mumps-seq lglib libff ${MPI_C_LIBRARIES})
  SET_TARGET_PROPERTIES(mumps-seq PROPERTIES PREFIX "")
ELSE()
  LIST(REMOVE_ITEM EDP_FILES ${CMAKE_SOURCE_DIR}/examples++-load/LapMUMPS_seq.edp)
ENDIF(MPI_FOUND AND MUMPS_FOUND)

# MUMPS.cpp
IF(MPI_FOUND AND MUMPS_FOUND)
  ADD_LIBRARY(mumps SHARED MUMPS.cpp)
  INCLUDE_DIRECTORIES(${MPI_C_INCLUDE_PATH})
  TARGET_LINK_LIBRARIES(mumps lglib libff ${MPI_C_LIBRARIES})
  SET_TARGET_PROPERTIES(mumps PROPERTIES PREFIX "")
ENDIF(MPI_FOUND AND MUMPS_FOUND)

# mmg3d-v4.0.cpp
IF(FREEFEM_MMG3D_INSTALLED)
  ADD_LIBRARY(mmg3d-v4.0 SHARED mmg3d-v4.0.cpp)
  INCLUDE_DIRECTORIES(${FF_DOWNLOAD_DIR}/mmg3d/install/include)
ENDIF()

# pardiso
IF(FREEFEM_PARDISO_INSTALLED)
  ADD_LIBRARY(PARDISO MODULE PARDISO.cpp)
ELSE()
  LIST(REMOVE_ITEM CPPFILES ${CMAKE_SOURCE_DIR}/examples++-load/PARDISO.cpp)
ENDIF()

# scotch.cpp
IF(SCOTCH_FOUND)
  ADD_LIBRARY(scotch SHARED scotch.cpp)
  INCLUDE_DIRECTORIES(${SCOTCH_INCLUDES})
  TARGET_LINK_LIBRARIES(scotch libff lglib ${SCOTCH_LIBRARIES})
  SET_TARGET_PROPERTIES(scotch PROPERTIES PREFIX "")
ELSE()
  LIST(REMOVE_ITEM EDP_FILES ${CMAKE_SOURCE_DIR}/examples++-load/scotch.edp)
ENDIF()

#superlu.cpp
IF(SUPERLU_FOUND)
  ADD_LIBRARY(superlu SHARED SuperLu.cpp)
  INCLUDE_DIRECTORIES(${SUPERLU_INCLUDES})
  TARGET_LINK_LIBRARIES(superlu ${SUPERLU_LIBRARIES})
  SET_TARGET_PROPERTIES(superlu PROPERTIES PREFIX "")
ELSE()
  LIST(REMOVE_ITEM EDP_FILES ${CMAKE_SOURCE_DIR}/examples++-load/lap-solvers.edp)
  LIST(REMOVE_ITEM EDP_FILES ${CMAKE_SOURCE_DIR}/examples++-load/SuperLU.edp)
ENDIF(SUPERLU_FOUND)


# tetgen.cpp
IF(FREEFEM_TETGEN_INSTALLED)
  ADD_LIBRARY(tetgen MODULE tetgen.cpp)
  TARGET_LINK_LIBRARIES(tetgen FREEFEM::TETGEN)
  SET_TARGET_PROPERTIES(tetgen PROPERTIES PREFIX "")
ENDIF(FREEFEM_TETGEN_INSTALLED)



#UMFPACK64.cpp
IF (FREEFEM_UMFPACK_INSTALLED)
  ADD_LIBRARY(umfpack64 UMFPACK64.cpp)
  TARGET_LINK_LIBRARIES(umfpack64 /usr/lib/x86_64-linux-gnu/libumfpack.a)
  SET_TARGET_PROPERTIES(umfpack64 PROPERTIES PREFIX "")
  INCLUDE_DIRECTORIES(/usr/include/suitesparse)
ELSE()
  LIST(REMOVE_ITEM EDP_FILES ${CMAKE_SOURCE_DIR}/examples++-load/LapUmfpack64.edp)
ENDIF()

#NewSolver.cpp
IF (FREEFEM_UMFPACK_INSTALLED)
  ADD_LIBRARY(NewSolver NewSolver.cpp)
  TARGET_LINK_LIBRARIES(NewSolver /usr/lib/x86_64-linux-gnu/libumfpack.a)
  SET_TARGET_PROPERTIES(NewSolver PROPERTIES PREFIX "")
  INCLUDE_DIRECTORIES(/usr/include/suitesparse)
ELSE()
  LIST(REMOVE_ITEM EDP_FILES ${CMAKE_SOURCE_DIR}/examples++-load/LapNewSolver.edp)
ENDIF()


FOREACH(EDP_FILE ${EDP_FILES})
  GET_FILENAME_COMPONENT(FILENAME ${EDP_FILE} NAME_WE)
  ADD_TEST(NAME ${FILENAME} COMMAND ${CMAKE_BINARY_DIR}/src/FreeFem++ ${EDP_FILE}
                            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/examples++-load)
ENDFOREACH(EDP_FILE)


