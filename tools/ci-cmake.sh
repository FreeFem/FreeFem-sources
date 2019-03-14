#TODO: manage LD_LIBRARY_PATH with CMake and rpath
#TODO: LD_LIBRARY_PATH is added to load msh3.so from examples++-load but it used by edp files coming from examples++-3d (e.g. sphere6.edp)
#TODO: msh3 should be removed from examples++-load and put in a common directory for all examples

#export FF_ROOT=/builds/workspace/FreeFem-source-feature-cmake-UbuntuAll 
#export FF_LOADPATH=$FF_ROOT/build_cmake/examples++-load 
#export LD_LIBRARY_PATH=$FF_ROOT/build_cmake/examples++-load 

export FF_ROOT=$(pwd)
export FF_INCLUDEPATH=$FF_ROOT/idp 

mkdir -p build_cmake \
&& cd build_cmake \
&& cmake -D CMAKE_CXX_COMPILER=g++ -D CMAKE_C_COMPILER=gcc -D CMAKE_INSTALL_PREFIX=/builds/freefem-source-feature-cmake .. \
&& make -j 8 VERBOSE=1 \
&& ls src \
&& ls src/fflib \
&& ldd src/FreeFem++ \
&& make install

cygcheck src/FreeFem++ 
cygcheck src/FreeFem++.exe

cd $FF_ROOT/build_cmake/examples/bug 
make test CTEST_OUTPUT_ON_FAILURE=On 
cd $FF_ROOT/build_cmake/examples/3d 
make test CTEST_OUTPUT_ON_FAILURE=On 
cd $FF_ROOT/build_cmake/examples/eigen 
make test CTEST_OUTPUT_ON_FAILURE=On 
cd $FF_ROOT/build_cmake/examples/mpi 
make test CTEST_OUTPUT_ON_FAILURE=On 
cd $FF_ROOT/build_cmake/examples/other 
make test CTEST_OUTPUT_ON_FAILURE=On 
cd $FF_ROOT/build_cmake/examples/tutorial 
make test CTEST_OUTPUT_ON_FAILURE=On 
 
 
