#TODO: manage LD_LIBRARY_PATH with CMake and rpath

mkdir -p build_cmake \
&&  cd build_cmake \
&&  cmake -D CMAKE_INSTALL_PREFIX=/builds/freefem-source-feature-cmake .. \
&&  make -j 8 \
&&  make install \
&&  export LD_LIBRARY_PATH=/builds/workspace/FreeFem-source-feature-cmake-UbuntuAll/build_cmake/example++-load \
&&  make test CTEST_OUTPUT_ON_FAILURE=On 
