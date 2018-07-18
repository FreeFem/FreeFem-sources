#TODO: manage LD_LIBRARY_PATH with CMake and rpath
#TODO: LD_LIBRARY_PATH is added to load msh3.so from examples++-load but it used by edp files coming from examples++-3d (e.g. sphere6.edp)
#TODO: msh3 should be removed from examples++-load and put in a common directory for all examples

mkdir -p build_cmake \
&&  cd build_cmake \
&&  cmake -D CMAKE_INSTALL_PREFIX=/builds/freefem-source-feature-cmake .. \
&&  make -j 8 \
&&  make install \
&&  export LD_LIBRARY_PATH=/builds/workspace/FreeFem-source-feature-cmake-UbuntuAll/build_cmake/examples++-load \
&&  make test CTEST_OUTPUT_ON_FAILURE=On 
