cmake_minimum_required(VERSION 3.0)

project(FreeFem++ C CXX Fortran)

set(FREEFEM_VERSION 4.0)

# Add the path containing cmake modules
# Note: this line cannot be put in the body of ff_configure_cmake
# because ff_configure_cmake is searched in the right repository 
# thanks to this line
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake/modules)

include(ff_configure_cmake)
include(ff_configure_compilers)
include(ff_configure_thirdparty)
include(ff_write_config_file)

ff_configure_cmake()
ff_configure_compilers()
ff_configure_thirdparty()

ff_write_config_file()

# creation of "make test" command
# Note: this call MUST be performed in the main CMake script
enable_testing()

add_subdirectory(src)
add_subdirectory(examples)


