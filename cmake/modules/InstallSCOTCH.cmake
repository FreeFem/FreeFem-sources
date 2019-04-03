SET(LIBNAME SCOTCH)

SET(URL https://gforge.inria.fr/frs/download.php/file/34618/scotch_6.0.4.tar.gz)
SET(URL_MD5 d58b825eb95e1db77efe8c6ff42d329f)

SET(CMAKE_ARGS -D CMAKE_INSTALL_PREFIX=<INSTALL_DIR>
               -D CMAKE_C_COMPILER=${CMAKE_C_COMPILER})
