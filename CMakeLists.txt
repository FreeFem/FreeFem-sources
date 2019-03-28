cmake_minimum_required(VERSION 3.0)

project(FreeFem++ C CXX Fortran)

set(FREEFEM_VERSION 4.0)

#cmake_policy(SET CMP0060 NEW) 
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake/modules)

include(ff_configure_cmake)
include(ff_configure_compilers)
include(ff_configure_thirdparty)
include(ff_write_config_file)

if(NOT ENABLE_DOWNLOAD)
  set(ENABLE_DOWNLOAD False)
endif(NOT ENABLE_DOWNLOAD)

ff_configure_cmake()
ff_configure_compilers()
ff_configure_thirdparty()

ff_write_config_file()


###########
# DOWNLOAD
##########



# creation of "make test" command
# Note: this call MUST be performed in the main CMake script
enable_testing()

#add_subdirectory(plugin)
add_subdirectory(src)
add_subdirectory(examples)


