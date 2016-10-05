FILE(GLOB_RECURSE CPPFILES *.cpp)
ADD_LIBRARY(lglib ${CPPFILES})
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/bamglib
                    ${CMAKE_SOURCE_DIR}/src/femlib
                    ${CMAKE_SOURCE_DIR}/src/fflib
                    ${CMAKE_SOURCE_DIR}/src/Graphics
                    ${CMAKE_SOURCE_DIR}/src/lglib/)
SET_TARGET_PROPERTIES(lglib PROPERTIES PREFIX "")

