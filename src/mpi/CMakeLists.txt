FIND_PACKAGE(MPI)
FILE(GLOB_RECURSE CPPFILES *.cpp)
ADD_LIBRARY(ff-mpi ${CPPFILES})
INCLUDE_DIRECTORIES(${MPI_C_INCLUDE_PATH}
                    ${CMAKE_SOURCE_DIR}/src/femlib
                    ${CMAKE_SOURCE_DIR}/src/fflib
                    ${CMAKE_SOURCE_DIR}/src/Graphics)
