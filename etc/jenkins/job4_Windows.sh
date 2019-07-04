#!C:\msys64\usr\bin\bash.exe --login
source shell mingw64

echo "Job 4"

autoreconf -i \
&& ./configure --enable-generic --enable-optim --enable-download --enable-maintainer-mode \
        CXXFLAGS=-mtune=generic CFLAGS=-mtune=generic FFLAGS=-mtune=generic \
        --prefix=/builds/workspace/freefem \
  && ./3rdparty/getall -a \
  && cd 3rdparty/ff-petsc \
  && make petsc-slepc \
  && cd - \
  && ./reconfigure \
  && make -j2

if [ $? -eq 0 ]
then
  echo "Build process complete"
else
  echo "Build process failed"
  exit 1
fi

# check
make check

if [ $? -eq 0 ]
then
  echo "Check process complete"
else
  echo "Check process failed"
fi

# install
make install

if [ $? -eq 0 ]
then
  echo "Install process complete"
else
  echo "Install process failed"
fi
