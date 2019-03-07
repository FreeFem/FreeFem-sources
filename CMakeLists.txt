CMAKE_MINIMUM_REQUIRED(VERSION 3.0)

PROJECT(FreeFem++ C CXX)

SET(FREEFEM_VERSION 4.0)

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)

include(ff_configure_cmake)
include(ff_configure_compilers)


ff_configure_cmake()
ff_configure_compilers()


INCLUDE(ModuleManager)


FIND_PACKAGE(MPI)

###########
# DOWNLOAD
##########

INCLUDE(FFInstallPackage)

IF(NOT ENABLE_DOWNLOAD)
  SET(ENABLE_DOWNLOAD False)
ENDIF(NOT ENABLE_DOWNLOAD)

ADD_SUBDIRECTORY(plugin)
ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(examples)
