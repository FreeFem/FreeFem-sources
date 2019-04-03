SET(LIBNAME MMG3D)
SET(URL http://www.math.u-bordeaux1.fr/~dobrzyns/logiciels/download/mmg3d4.0.tgz)
SET(CONFIGURE_COMMAND cmake -D CMAKE_INSTALL_PREFIX=<INSTALL_DIR> 
                            -D CMAKE_C_COMPILER=${CMAKE_C_COMPILER} 
                            -D CMAKE_C_FLAGS=${CMAKE_C_FLAGS}
                            -D INCLUDE_SCOTCH=/usr/include/scotch
                            -D COMPIL_SHARED_LIBRARY=True
                            ../src/build)
SET(INSTALL_COMMAND mkdir -p <INSTALL_DIR>/include &&
                    mkdir -p <INSTALL_DIR>/lib && 
                    cp <SOURCE_DIR>/build/sources/libmmg3d.h <INSTALL_DIR>/include &&
                    cp libmmg3dlib4.0.so <INSTALL_DIR>/lib/)


