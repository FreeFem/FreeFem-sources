
IF(MPI_FOUND)

  FILE(GLOB CPPFILES *.cpp)

  LIST(REMOVE_ITEM CPPFILES ${CMAKE_SOURCE_DIR}/examples++-mpi/complex_SuperLU_DIST_FreeFem.cpp
                            ${CMAKE_SOURCE_DIR}/examples++-mpi/complex_pastix_FreeFem.cpp
                            ${CMAKE_SOURCE_DIR}/examples++-mpi/dSuperLU_DIST.cpp
                            ${CMAKE_SOURCE_DIR}/examples++-mpi/example-mumps.cpp
                            ${CMAKE_SOURCE_DIR}/examples++-mpi/hips_FreeFem.cpp
                            ${CMAKE_SOURCE_DIR}/examples++-mpi/hypre_FreeFem.cpp
                            ${CMAKE_SOURCE_DIR}/examples++-mpi/interfacepastix.cpp
                            ${CMAKE_SOURCE_DIR}/examples++-mpi/mpi-cmaes.cpp
                            ${CMAKE_SOURCE_DIR}/examples++-mpi/MUMPS.cpp
                            ${CMAKE_SOURCE_DIR}/examples++-mpi/MUMPS_FreeFem.cpp
                            ${CMAKE_SOURCE_DIR}/examples++-mpi/parms_FreeFem.cpp
                            ${CMAKE_SOURCE_DIR}/examples++-mpi/real_pastix_FreeFem.cpp
                            ${CMAKE_SOURCE_DIR}/examples++-mpi/real_SuperLU_DIST_FreeFem.cpp
                            ${CMAKE_SOURCE_DIR}/examples++-mpi/removeDOF.cpp)

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

ENDIF(MPI_FOUND)

