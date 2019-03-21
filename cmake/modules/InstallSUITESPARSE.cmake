SET(LIBNAME SUITESPARSE)
SET(URL http://faculty.cse.tamu.edu/davis/SuiteSparse/SuiteSparse-5.4.0.tar.gz)
SET(URL_MD5 4a6d4e74fc44c503f52996ae95cad03a)

SET(PATCH_COMMAND cp ${CMAKE_SOURCE_DIR}/cmake/modules/suitesparse/CMakeLists.txt <SOURCE_DIR>)
SET(INSTALL_COMMAND cd .)
