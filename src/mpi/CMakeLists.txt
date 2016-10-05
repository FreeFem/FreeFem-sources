FILE(GLOB_RECURSE CPPFILES *.cpp)
ADD_LIBRARY(mpi ${CPPFILES})
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/femlib
                    ${CMAKE_SOURCE_DIR}/src/fflib
                    ${CMAKE_SOURCE_DIR}/src/Graphics)
