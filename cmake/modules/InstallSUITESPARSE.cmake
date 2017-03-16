SET(LIBNAME SUITESPARSE)
SET(URL http://faculty.cse.tamu.edu/davis/SuiteSparse/SuiteSparse-4.4.4.tar.gz)
SET(URL_MD5 e0af74476935c9ff6d971df8bb6b82fc)

SET(PATCH_COMMAND cp ${CMAKE_SOURCE_DIR}/cmake/modules/suitesparse/CMakeLists.txt <SOURCE_DIR>)
