IF(MPI_FOUND)
  FILE(GLOB_RECURSE CPPFILES *.cpp)
  LIST(REMOVE_ITEM CPPFILES ${CMAKE_SOURCE_DIR}/src/fflib/strversionnumber.cpp)
  ADD_LIBRARY(ff-mpi ${CPPFILES})
  TARGET_LINK_LIBRARIES(ff-mpi strversionnumber)
  INCLUDE_DIRECTORIES(${MPI_C_INCLUDE_PATH}
                      ${CMAKE_SOURCE_DIR}/src/femlib
                      ${CMAKE_SOURCE_DIR}/src/fflib
                      ${CMAKE_SOURCE_DIR}/src/Graphics)
ENDIF(MPI_FOUND)
