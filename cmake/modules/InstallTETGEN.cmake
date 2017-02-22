SET(LIBNAME TETGEN)

SET(URL http://www.tetgen.org/1.5/src/tetgen1.5.1-beta1.tar.gz)
SET(URL_MD5 3d55c197bcbfc611b7ced6f343643756)

SET(CMAKE_ARGS -D CMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} 
               -D CMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS} 
               -D CMAKE_INSTALL_PREFIX=<INSTALL_DIR>)
SET(INSTALL_COMMAND mkdir -p ../install/lib && cp libtet.a ../install/lib && mkdir -p ../install/include && cp <SOURCE_DIR>/tetgen.h ../install/include)

SET(LIBRARIES libtet.a)
