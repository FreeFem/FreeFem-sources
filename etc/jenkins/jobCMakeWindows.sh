#!C:\msys64\usr\bin\bash.exe --login
source shell mingw64

echo "Job CMake"

export FF_ROOT=$(pwd)
export FF_INCLUDEPATH=$FF_ROOT/idp

mkdir -p build_cmake \
  && cd build_cmake \
  && cmake -D CMAKE_CXX_COMPILER=g++ -D CMAKE_C_COMPILER=gcc -D CMAKE_INSTALL_PREFIX=/builds/workspace/freefem .. \
  && make -j 8 VERBOSE=1 \
  && make install \
  && cd $FF_ROOT/build_cmake/examples/3d \
  && make test CTEST_OUTPUT_ON_FAILURE=On \
  && cd $FF_ROOT/build_cmake/examples/eigen \
  && make test CTEST_OUTPUT_ON_FAILURE=On \
  && cd $FF_ROOT/build_cmake/examples/mpi \
  && make test CTEST_OUTPUT_ON_FAILURE=On \
  && cd $FF_ROOT/build_cmake/examples/other \
  && make test CTEST_OUTPUT_ON_FAILURE=On \
  && cd $FF_ROOT/build_cmake/examples/tutorial \
  && make test CTEST_OUTPUT_ON_FAILURE=On
