CMAKE_MINIMUM_REQUIRED(VERSION 3.0)

PROJECT(FreeFem++ C CXX)

SET(FREEFEM_VERSION 3.5)

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)

INCLUDE(ProcessorCount)
ProcessorCount(N)
IF(NOT N EQUAL 0)
  SET(CMAKE_BUILD_FLAGS -j${N})
ENDIF()

######################
# C++ COMPILER OPTIONS
######################

INCLUDE(CompilerManager)
INCLUDE(ModuleManager)

FILE(GLOB FEMLIB_SRC src/femlib/*.cpp)
FILE(GLOB_RECURSE FFLIB_SRC src/fflib/*.cpp)

######
# TEST
######

FILE(GLOB EXAMPLE_DIRECTORIES RELATIVE ${CMAKE_SOURCE_DIR} examples*)

FIND_PACKAGE(MPI)
IF(MPI_FOUND)
  ADD_CUSTOM_TARGET(check LD_LIBRARY_PATH=${CMAKE_BINARY_DIR}/examples++-mpi:${CMAKE_BINARY_DIR}/examples++-load ${CMAKE_SOURCE_DIR}/CheckAll "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ${CMAKE_BINARY_DIR}/src/FreeFem++-mpi ")
ELSE(MPI_FOUND)
  ADD_CUSTOM_TARGET(check LD_LIBRARY_PATH=${CMAKE_BINARY_DIR}/examples++-load ${CMAKE_SOURCE_DIR}/CheckAll ${CMAKE_BINARY_DIR}/src/FreeFem++)
ENDIF(MPI_FOUND)

###########
# DOWNLOAD
##########

INCLUDE(FFInstallPackage)

IF(NOT ENABLE_DOWNLOAD)
  SET(ENABLE_DOWNLOAD False)
ENDIF(NOT ENABLE_DOWNLOAD)

ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(examples++-load)
ADD_SUBDIRECTORY(examples++-mpi)

##################
# config.h
##################

FILE(WRITE config.h "")

##################
# config-wrapper.h
##################

FILE(WRITE config-wrapper.h "")

#############################################################
# Include root for config.h and config-wrapper.h
# Only temporarily for coexistence of both autoconf and cmake
# TODO: remove as soon as autoconf is removed
#############################################################

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

######################
# strversionnumber.cpp
######################

EXECUTE_PROCESS(COMMAND "date" OUTPUT_VARIABLE VersionFreeFemDate)
STRING(STRIP ${VersionFreeFemDate} VersionFreeFemDate)
FILE(WRITE src/fflib/strversionnumber.cpp
"#include \"config-wrapper.h\"\n
#include \"strversionnumber.hpp\" // [[file:strversionnumber.hpp]] \n
#include <cstdlib> \n
using namespace std;\n
#define TOSTRING1(i) #i \n
#define TOSTRING(i) TOSTRING1(i) \n
\n
//#include <sstream>\n
#include <cstdio>\n
using namespace std;\n
\n
double VersionNumber(){\n
  return VersionFreeFempp;\n
}\n
\n
string StrVersionNumber(){\n
//  std::ostringstream buffer;\n
//  buffer.precision(8);\n
//  buffer<<VersionNumber();\n
  static char buffer[100];\n
  sprintf(buffer,\" %9f (date ' ${VersionFreeFemDate}')\",VersionNumber());\n
  return buffer; //.str()+\" (date ${VersionFreeFemDate})\" ;\n
}")

