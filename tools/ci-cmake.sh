#TODO: manage LD_LIBRARY_PATH with CMake and rpath
#TODO: LD_LIBRARY_PATH is added to load msh3.so from examples++-load but it used by edp files coming from examples++-3d (e.g. sphere6.edp)
#TODO: msh3 should be removed from examples++-load and put in a common directory for all examples

mkdir -p build_cmake \
&&  cd build_cmake \
&&  cmake -D CMAKE_INSTALL_PREFIX=/builds/freefem-source-feature-cmake .. \
&&  make -j 8 \
&&  make install \
&&  export LD_LIBRARY_PATH=/builds/workspace/FreeFem-source-feature-cmake-UbuntuAll/build_cmake/examples++-load \
&&  export FF_INCLUDE_PATH=$FF_INCLUDEPATH:/builds/workspace/FreeFem-source-feature-cmake-UbuntuAll/examples++ \
&&  export FF_INCLUDE_PATH=$FF_INCLUDEPATH:/builds/workspace/FreeFem-source-feature-cmake-UbuntuAll/examples++-3d \
&&  export FF_INCLUDE_PATH=$FF_INCLUDEPATH:/builds/workspace/FreeFem-source-feature-cmake-UbuntuAll/examples++-bug \
&&  export FF_INCLUDE_PATH=$FF_INCLUDEPATH:/builds/workspace/FreeFem-source-feature-cmake-UbuntuAll/examples++-chapt3 \
&&  export FF_INCLUDE_PATH=$FF_INCLUDEPATH:/builds/workspace/FreeFem-source-feature-cmake-UbuntuAll/examples++-eigen \
&&  export FF_INCLUDE_PATH=$FF_INCLUDEPATH:/builds/workspace/FreeFem-source-feature-cmake-UbuntuAll/examples++-load \
&&  export FF_INCLUDE_PATH=$FF_INCLUDEPATH:/builds/workspace/FreeFem-source-feature-cmake-UbuntuAll/examples++-mpi \
&&  export FF_INCLUDE_PATH=$FF_INCLUDEPATH:/builds/workspace/FreeFem-source-feature-cmake-UbuntuAll/examples++-other \
&&  export FF_INCLUDE_PATH=$FF_INCLUDEPATH:/builds/workspace/FreeFem-source-feature-cmake-UbuntuAll/examples++-tutorial \
&&  make test CTEST_OUTPUT_ON_FAILURE=On 
