# TODO: remove this line as soon as autoconf is removed (only here for config.h and config-wrapper.h)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}) 

 
IF(MPI_FOUND)

  FILE(GLOB CPPFILES *.cpp)

  LIST(REMOVE_ITEM CPPFILES ${CMAKE_SOURCE_DIR}/plugin/mpi/complex_SuperLU_DIST_FreeFem.cpp
                            ${CMAKE_SOURCE_DIR}/plugin/mpi/complex_pastix_FreeFem.cpp
                            ${CMAKE_SOURCE_DIR}/plugin/mpi/dSuperLU_DIST.cpp
                            ${CMAKE_SOURCE_DIR}plugin/mpi/example-mumps.cpp
                            ${CMAKE_SOURCE_DIR}/plugin/mpi/hips_FreeFem.cpp
                            ${CMAKE_SOURCE_DIR}/plugin/mpi/hypre_FreeFem.cpp
                            ${CMAKE_SOURCE_DIR}/plugin/mpi/interfacepastix.cpp
                            ${CMAKE_SOURCE_DIR}/plugin/mpi/mpi-cmaes.cpp
                            ${CMAKE_SOURCE_DIR}/plugin/mpi/MUMPS.cpp
                            ${CMAKE_SOURCE_DIR}/plugin/mpi/MUMPS_FreeFem.cpp
                            ${CMAKE_SOURCE_DIR}/plugin/mpi/parmetis.cpp
                            ${CMAKE_SOURCE_DIR}/plugin/mpi/parms_FreeFem.cpp
                            ${CMAKE_SOURCE_DIR}/plugin/mpi/real_pastix_FreeFem.cpp
                            ${CMAKE_SOURCE_DIR}/plugin/mpi/real_SuperLU_DIST_FreeFem.cpp
                            ${CMAKE_SOURCE_DIR}plugin/mpi/removeDOF.cpp)

  FOREACH(CPPFILE ${CPPFILES})
    GET_FILENAME_COMPONENT(FILENAME ${CPPFILE} NAME_WE)
    ADD_LIBRARY(${FILENAME} MODULE ${CPPFILE})
    TARGET_LINK_LIBRARIES(${FILENAME} ${MPI_CXX_LIBRARIES})
    INCLUDE_DIRECTORIES(${MPI_CXX_INCLUDE_PATH} 
                        ${CMAKE_SOURCE_DIR}/src/bamglib
                        ${CMAKE_SOURCE_DIR}/src/Graphics
                        ${CMAKE_SOURCE_DIR}/src/femlib
                        ${CMAKE_SOURCE_DIR}/src/fflib)
    SET_TARGET_PROPERTIES(${FILENAME} PROPERTIES PREFIX "")
  ENDFOREACH(CPPFILE)

  FIND_PACKAGE(PARMETIS)
  IF(PARMETIS_FOUND)
    ADD_LIBRARY(parmetis MODULE ${CMAKE_SOURCE_DIR}/plugin/mpi/parmetis.cpp)
    TARGET_LINK_LIBRARIES(parmetis ${MPI_CXX_LIBRARIES})
    INCLUDE_DIRECTORIES(${MPI_CXX_INCLUDE_PATH} 
                        ${CMAKE_SOURCE_DIR}/src/bamglib
                        ${CMAKE_SOURCE_DIR}/src/Graphics
                        ${CMAKE_SOURCE_DIR}/src/femlib
                        ${CMAKE_SOURCE_DIR}/src/fflib)
    SET_TARGET_PROPERTIES(parmetis PROPERTIES PREFIX "")
  ENDIF(PARMETIS_FOUND) 

ENDIF(MPI_FOUND)


