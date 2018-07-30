#TODO: manage LD_LIBRARY_PATH with CMake and rpath
#TODO: LD_LIBRARY_PATH is added to load msh3.so from examples++-load but it used by edp files coming from examples++-3d (e.g. sphere6.edp)
#TODO: msh3 should be removed from examples++-load and put in a common directory for all examples

mkdir -p build_cmake \
&& cd build_cmake \
&& cmake -D CMAKE_INSTALL_PREFIX=/builds/freefem-source-feature-cmake .. \
&& make -j 8 \
&& make install \
&& export FF_ROOT=/builds/workspace/FreeFem-source-feature-cmake-UbuntuAll \
&& export FF_LOADPATH=$FF_ROOT/build_cmake/examples++-load \
&& export LD_LIBRARY_PATH=$FF_ROOT/build_cmake/examples++-load \
&& export FF_INCLUDEPATH="$FF_ROOT/examples++" \
&& export FF_INCLUDEPATH="$FF_INCLUDEPATH;$FF_ROOT/examples++-3d" \
&& export FF_INCLUDEPATH="$FF_INCLUDEPATH;$FF_ROOT/examples++-bamg" \
&& export FF_INCLUDEPATH="$FF_INCLUDEPATH;$FF_ROOT/examples++-bug" \
&& export FF_INCLUDEPATH="$FF_INCLUDEPATH;$FF_ROOT/examples++-chapt3" \
&& export FF_INCLUDEPATH="$FF_INCLUDEPATH;$FF_ROOT/examples++-eigen" \
&& export FF_INCLUDEPATH="$FF_INCLUDEPATH;$FF_ROOT/examples++-load" \
&& export FF_INCLUDEPATH="$FF_INCLUDEPATH;$FF_ROOT/examples++-mpi" \
&& export FF_INCLUDEPATH="$FF_INCLUDEPATH;$FF_ROOT/examples++-other" \
&& export FF_INCLUDEPATH="$FF_INCLUDEPATH;$FF_ROOT/examples++-tutorial" \
&& make test CTEST_OUTPUT_ON_FAILURE=On 
 
