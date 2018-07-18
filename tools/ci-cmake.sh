mkdir -p build_cmake \
&&  cd build_cmake \
&&  cmake -D CMAKE_INSTALL_PREFIX=/builds/freefem-source-feature-cmake .. \
&&  make -j 8 \
&&  make install \
&&  make test CTEST_OUTPUT_ON_FAILURE=On 
