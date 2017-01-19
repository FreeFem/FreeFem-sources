CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(FreeFem++ C CXX Fortran)

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
ADD_CUSTOM_TARGET(check ${CMAKE_SOURCE_DIR}/CheckAll)

###########
# DOWNLOAD
##########

INCLUDE(FFInstallPackage)

ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(examples++-load)
