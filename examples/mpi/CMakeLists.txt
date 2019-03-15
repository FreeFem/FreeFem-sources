# TODO: remove this line as soon as autoconf is removed (only here for config.h and config-wrapper.h)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}) 

 
IF(MPI_FOUND)

#  FILE(GLOB CPPFILES *.cpp)
#
#  LIST(REMOVE_ITEM CPPFILES ${CMAKE_SOURCE_DIR}/examples++-mpi/complex_SuperLU_DIST_FreeFem.cpp
#                            ${CMAKE_SOURCE_DIR}/examples++-mpi/complex_pastix_FreeFem.cpp
#                            ${CMAKE_SOURCE_DIR}/examples++-mpi/dSuperLU_DIST.cpp
#                            ${CMAKE_SOURCE_DIR}/examples++-mpi/example-mumps.cpp
#                            ${CMAKE_SOURCE_DIR}/examples++-mpi/hips_FreeFem.cpp
#                            ${CMAKE_SOURCE_DIR}/examples++-mpi/hypre_FreeFem.cpp
#                            ${CMAKE_SOURCE_DIR}/examples++-mpi/interfacepastix.cpp
#                            ${CMAKE_SOURCE_DIR}/examples++-mpi/mpi-cmaes.cpp
#                            ${CMAKE_SOURCE_DIR}/examples++-mpi/MUMPS.cpp
#                            ${CMAKE_SOURCE_DIR}/examples++-mpi/MUMPS_FreeFem.cpp
#                            ${CMAKE_SOURCE_DIR}/examples++-mpi/parmetis.cpp
#                            ${CMAKE_SOURCE_DIR}/examples++-mpi/parms_FreeFem.cpp
#                            ${CMAKE_SOURCE_DIR}/examples++-mpi/real_pastix_FreeFem.cpp
#                            ${CMAKE_SOURCE_DIR}/examples++-mpi/real_SuperLU_DIST_FreeFem.cpp
#                            ${CMAKE_SOURCE_DIR}/examples++-mpi/removeDOF.cpp)
#
#  FOREACH(CPPFILE ${CPPFILES})
#    GET_FILENAME_COMPONENT(FILENAME ${CPPFILE} NAME_WE)
#    ADD_LIBRARY(${FILENAME} MODULE ${CPPFILE})
#    TARGET_LINK_LIBRARIES(${FILENAME} ${MPI_CXX_LIBRARIES})
#    INCLUDE_DIRECTORIES(${MPI_CXX_INCLUDE_PATH} 
#                        ${CMAKE_SOURCE_DIR}/src/bamglib
#                        ${CMAKE_SOURCE_DIR}/src/Graphics
#                        ${CMAKE_SOURCE_DIR}/src/femlib
#                        ${CMAKE_SOURCE_DIR}/src/fflib)
#    SET_TARGET_PROPERTIES(${FILENAME} PROPERTIES PREFIX "")
#  ENDFOREACH(CPPFILE)

#  IF(PARMETIS_FOUND)
#    ADD_LIBRARY(parmetis MODULE parmetis.cpp)
#    TARGET_LINK_LIBRARIES(parmetis ${MPI_CXX_LIBRARIES})
#    INCLUDE_DIRECTORIES(${MPI_CXX_INCLUDE_PATH} 
#                        ${CMAKE_SOURCE_DIR}/src/bamglib
#                        ${CMAKE_SOURCE_DIR}/src/Graphics
#                        ${CMAKE_SOURCE_DIR}/src/femlib
#                        ${CMAKE_SOURCE_DIR}/src/fflib)
#    SET_TARGET_PROPERTIES(parmetis PROPERTIES PREFIX "")
#  ENDIF(PARMETIS_FOUND) 

FILE(GLOB EDP_FILES *.edp)

# tests which are known to fail are removed from test list  
LIST(REMOVE_ITEM EDP_FILES ${CMAKE_CURRENT_SOURCE_DIR}/chaleur3D-hips.edp
                           ${CMAKE_CURRENT_SOURCE_DIR}/chaleur3D-mumps.edp
                           ${CMAKE_CURRENT_SOURCE_DIR}/chaleur3D-superludist.edp
                           ${CMAKE_CURRENT_SOURCE_DIR}/essai-com.edp
                           ${CMAKE_CURRENT_SOURCE_DIR}/schwarz-2d.edp
                           ${CMAKE_CURRENT_SOURCE_DIR}/schwarz-b.edp
                           ${CMAKE_CURRENT_SOURCE_DIR}/schwarz-c.edp)


IF(NOT EXISTS $ENV{FF_INCLUDEPATH}/getARGV.idp OR NOT FF_PLUGIN_MUMPS)
  LIST(REMOVE_ITEM EDP_FILES ${CMAKE_SOURCE_DIR}/examples/mpi/NSCaraCyl.edp)
ENDIF(NOT EXISTS $ENV{FF_INCLUDEPATH}/getARGV.idp OR NOT FF_PLUGIN_MUMPS)
  
IF(NOT EXISTS $ENV{FF_INCLUDEPATH}/getARGV.idp OR NOT FF_PLUGIN_METIS)
  LIST(REMOVE_ITEM EDP_FILES ${CMAKE_SOURCE_DIR}/examples/mpi/parmetis.edp)
ENDIF(NOT EXISTS $ENV{FF_INCLUDEPATH}/getARGV.idp OR NOT FF_PLUGIN_METIS)
  

IF(NOT EXISTS $ENV{FF_INCLUDEPATH}/getARGV.idp OR NOT FF_PLUGIN_MSH3)
  LIST(REMOVE_ITEM EDP_FILES ${CMAKE_SOURCE_DIR}/examples/mpi/parmetis-3d.edp)
ENDIF(NOT EXISTS $ENV{FF_INCLUDEPATH}/getARGV.idp OR NOT FF_PLUGIN_MSH3)
  


IF(NOT EXISTS $ENV{FF_INCLUDEPATH}/mortar-msh.idp)
  LIST(REMOVE_ITEM EDP_FILES ${CMAKE_SOURCE_DIR}/examples/mpi/mortar-DN-4-mpi.edp)
ENDIF(NOT EXISTS $ENV{FF_INCLUDEPATH}/mortar-msh.idp)
  



IF (NOT FF_PLUGIN_MPICG)
  LIST(REMOVE_ITEM EDP_FILES ${CMAKE_SOURCE_DIR}/examples/mpi/DDM-Schwarz-Lame-2d.edp
                             ${CMAKE_SOURCE_DIR}/examples/mpi/DDM-Schwarz-Lame-3d.edp
                             ${CMAKE_SOURCE_DIR}/examples/mpi/DDM-Schwarz-Lap-2dd.edp
                             ${CMAKE_SOURCE_DIR}/examples/mpi/DDM-Schwarz-Lap-3d.edp
                             ${CMAKE_SOURCE_DIR}/examples/mpi/DDM-Schwarz-Stokes-2d.edp
                             ${CMAKE_SOURCE_DIR}/examples/mpi/MPICGLap.edp
                             ${CMAKE_SOURCE_DIR}/examples/mpi/MPIGMRES2D.edp
                             ${CMAKE_SOURCE_DIR}/examples/mpi/MPIGMRES3D.edp)
ENDIF(NOT FF_PLUGIN_MPICG)


IF (NOT FF_PLUGIN_MUMPS)
  LIST(REMOVE_ITEM EDP_FILES ${CMAKE_SOURCE_DIR}/examples/mpi/LapMPIEigenValue.edp
                             ${CMAKE_SOURCE_DIR}/examples/mpi/testsolver_MUMPS.edp
                             ${CMAKE_SOURCE_DIR}/examples/mpi/MUMPS.edp
                             ${CMAKE_SOURCE_DIR}/examples/mpi/VG.edp
                             ${CMAKE_SOURCE_DIR}/examples/mpi/cavityNewtow-MUMPS.edp)
ENDIF(NOT FF_PLUGIN_MUMPS)


IF (NOT FF_PLUGIN_MEDIT)
  LIST(REMOVE_ITEM EDP_FILES ${CMAKE_SOURCE_DIR}/examples/mpi/beam-3d-matrix-pastix.edp
                             ${CMAKE_SOURCE_DIR}/examples/mpi/LaplaceRT-3d-matrix-mumps.edp
                             ${CMAKE_SOURCE_DIR}/examples/mpi/LaplaceRT-3d-matrix-pastix.edp
                             ${CMAKE_SOURCE_DIR}/examples/mpi/LaplaceRT-3d-matrix-superludist.edp)
ENDIF(NOT FF_PLUGIN_MEDIT)


IF (NOT FF_PLUGIN_REAL_SUPERLU_DIST_FREEFEM)
  LIST(REMOVE_ITEM EDP_FILES ${CMAKE_SOURCE_DIR}/examples/mpi/testsolver_SuperLU_DIST.edp)
ENDIF(NOT FF_PLUGIN_REAL_SUPERLU_DIST_FREEFEM)


IF (NOT FF_PLUGIN_DSUPERLU_DIST)
  LIST(REMOVE_ITEM EDP_FILES ${CMAKE_SOURCE_DIR}/examples/mpi/testsolver_dsuperlu_dist.edp)
ENDIF(NOT FF_PLUGIN_DSUPERLU_DIST)


IF (NOT FF_PLUGIN_REAL_PASTIX_FREEFEM)
  LIST(REMOVE_ITEM EDP_FILES ${CMAKE_SOURCE_DIR}/examples/mpi/testsolver_pastix.edp)
ENDIF(NOT FF_PLUGIN_REAL_PASTIX_FREEFEM)


IF (NOT FF_PLUGIN_HTOOL)
  LIST(REMOVE_ITEM EDP_FILES ${CMAKE_SOURCE_DIR}/examples/mpi/htool.edp)
ENDIF(NOT FF_PLUGIN_HTOOL)


IF (NOT FF_PLUGIN_MPI_CMAES)
  LIST(REMOVE_ITEM EDP_FILES ${CMAKE_SOURCE_DIR}/examples/mpi/cmaes-mpi-VarIneq.edp)
ENDIF(NOT FF_PLUGIN_MPI_CMAES)


IF (NOT FF_PLUGIN_MSH3)
  LIST(REMOVE_ITEM EDP_FILES ${CMAKE_SOURCE_DIR}/examples/mpi/Laplace3d-hips.edp
                             ${CMAKE_SOURCE_DIR}/examples/mpi/Laplace3d-hypre.edp
                             ${CMAKE_SOURCE_DIR}/examples/mpi/NSI3d-carac-mumps.edp
                             ${CMAKE_SOURCE_DIR}/examples/mpi/NSI3d-carac-superludist.edp
                             ${CMAKE_SOURCE_DIR}/examples/mpi/NSI3d-carac-pastix.edp
                             ${CMAKE_SOURCE_DIR}/examples/mpi/Stokes-v1-matrix-mumps.edp
                             ${CMAKE_SOURCE_DIR}/examples/mpi/Stokes-v1-matrix-pastix.edp
                             ${CMAKE_SOURCE_DIR}/examples/mpi/Stokes-v1-matrix-superludist.edp
                             ${CMAKE_SOURCE_DIR}/examples/mpi/Stokes-v2-matrix-mumps.edp
                             ${CMAKE_SOURCE_DIR}/examples/mpi/Stokes-v2-matrix-pastix.edp
                             ${CMAKE_SOURCE_DIR}/examples/mpi/Stokes-v2-matrix-superludist.edp
                             ${CMAKE_SOURCE_DIR}/examples/mpi/Stokes-v3-matrix-mumps.edp
                             ${CMAKE_SOURCE_DIR}/examples/mpi/Stokes-v3-matrix-pastix.edp
                             ${CMAKE_SOURCE_DIR}/examples/mpi/Stokes-v3-matrix-superludist.edp
                             ${CMAKE_SOURCE_DIR}/examples/mpi/beam-3d-matrix-superludist.edp
                             ${CMAKE_SOURCE_DIR}/examples/mpi/chaleur3D-hypre.edp
                             ${CMAKE_SOURCE_DIR}/examples/mpi/chaleur3D-pastix.edp
                             ${CMAKE_SOURCE_DIR}/examples/mpi/schwarz-3.edp)
ENDIF(NOT FF_PLUGIN_MSH3)

IF(NOT MPIEXEC_MAX_NUMPROCS EQUAL 4)
  LIST(REMOVE_ITEM EDP_FILES ${CMAKE_SOURCE_DIR}/examples/mpi/mortar-DN-4-mpi.edp)

  message(STATUS "${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 
  ${MPIEXEC_PREFLAGS} ${CMAKE_BINARY_DIR}/src/FreeFem++-mpi ${MPIEXEC_POSTFLAGS} -nw mortar-DN-4-mpi.edp")
  
  ADD_TEST(NAME mortar-DN-4-mpi COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 
  ${MPIEXEC_PREFLAGS} ${CMAKE_BINARY_DIR}/src/FreeFem++-mpi ${MPIEXEC_POSTFLAGS} -nw mortar-DN-4-mpi.edp
                            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/examples/mpi)

ENDIF(NOT MPIEXEC_MAX_NUMPROCS EQUAL 4)

FOREACH(EDP_FILE ${EDP_FILES})
  GET_FILENAME_COMPONENT(FILENAME ${EDP_FILE} NAME_WE)
  ADD_TEST(NAME ${FILENAME} COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} 
  ${MPIEXEC_PREFLAGS} ${CMAKE_BINARY_DIR}/src/FreeFem++-mpi ${MPIEXEC_POSTFLAGS} -nw ${EDP_FILE}
                            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/examples/mpi)
ENDFOREACH(EDP_FILE)




ENDIF(MPI_FOUND)
