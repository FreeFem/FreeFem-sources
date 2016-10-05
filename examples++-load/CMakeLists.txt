FILE(GLOB CPPFILES *.cpp)
LIST(REMOVE_ITEM CPPFILES ${CMAKE_SOURCE_DIR}/examples++-load/bmo_user.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/buildlayer.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/Element_Mixte.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/dfft.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/ff-AiryBiry.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/ff-cmaes.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/ff-Ipopt.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/ff-NLopt.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/fflapack.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/ffnewuoa.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/ffrandom.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/FreeFemQA.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/freeyams.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/funcTemplate.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/gmsh.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/gsl.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/ilut.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/iohdf5.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/iovtk.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/isoline.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/isolineP1.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/lapack.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/lgbmo.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/mat_dervieux.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/mat_psi.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/medit.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/metis.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/MetricPkuate.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/mmg3d.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/mmg3d-v4.0.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/msh3.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/mshmet.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/MUMPS.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/MUMPS_seq.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/myfunction.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/myfunction2.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/netgen.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/NewSolver.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/PARDISO.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/pcm.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/pcm2rnm.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/pipe.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/ppm2rnm.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/qf11to25.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/shell.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/scotch.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/splitedges.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/splitmesh3.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/splitmesh6.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/SuperLu.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/symmetrizeCSR.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/tetgen.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/thresholdings.cpp
                          ${CMAKE_SOURCE_DIR}/examples++-load/UMFPACK64.cpp)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/bamglib
                    ${CMAKE_SOURCE_DIR}/src/fflib
                    ${CMAKE_SOURCE_DIR}/src/femlib
                    ${CMAKE_SOURCE_DIR}/src/Graphics
                    ${CMAKE_SOURCE_DIR}/src/libMesh)

FOREACH(CPPFILE ${CPPFILES})

  GET_FILENAME_COMPONENT(FILENAME ${CPPFILE} NAME_WE)
  ADD_LIBRARY(${FILENAME} SHARED ${CPPFILE})
  TARGET_LINK_LIBRARIES(${FILENAME} libMesh libff lglib)
ENDFOREACH(CPPFILE)

# Element_Mixte
FIND_PACKAGE(LAPACK)
ADD_LIBRARY(Element_Mixte SHARED Element_Mixte.cpp)
TARGET_LINK_LIBRARIES(Element_Mixte libMesh libff lglib ${LAPACK_LIBRARIES})

# lapack
ADD_LIBRARY(example-lapack SHARED lapack.cpp)
TARGET_LINK_LIBRARIES(example-lapack libMesh libff lglib ${LAPACK_LIBRARIES})


# tetgen.cpp
FF_INSTALL_PACKAGE(Tetgen)
ADD_LIBRARY(example_tetgen SHARED tetgen.cpp
                                  ${CMAKE_SOURCE_DIR}/src/femlib/BamgFreeFem.cpp
                                  ${CMAKE_SOURCE_DIR}/src/femlib/splitsimplex.cpp)
TARGET_LINK_LIBRARIES(example_tetgen libMesh libff lglib /Users/cdoucet/ff-download/tetgen/install/libtet.a)
ADD_DEPENDENCIES(example_tetgen tetgen)
INCLUDE_DIRECTORIES(${FF_DOWNLOAD_DIR}/tetgen/src)
