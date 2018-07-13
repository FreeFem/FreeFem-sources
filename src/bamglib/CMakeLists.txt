# TODO: remove this line as soon as autoconf is removed (only here for config.h and config-wrapper.h)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}) 

FILE(GLOB_RECURSE CPPFILES *.cpp)
LIST(REMOVE_ITEM CPPFILES ${CMAKE_SOURCE_DIR}/src/bamglib/Meshgibbs.cpp)
ADD_LIBRARY(bamglib ${CPPFILES})
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/fflib)
