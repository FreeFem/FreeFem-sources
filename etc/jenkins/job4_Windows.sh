#!C:\msys64\usr\bin\bash.exe --login
source shell mingw64

echo "Job 4"

autoreconf -i \
  && ./configure --with-mpi-lib=['/C/Program\ Files\ \(x86\)/Microsoft\ SDKs/MPI/Lib/x64/msmpi.lib'] \
                 --with-mpi-include=/C/Program\ Files\ \(x86\)/Microsoft\ SDKs/MPI/Include --with-fc=0 \
                 --download-metis --download-suitesparse --with-ssl=0 --with-x=0 --with-debugging=0  \
                 --download-f2cblaslapack=1 --with-mpiexec='/C/Program\ Files/Microsoft\ MPI/Bin/mpiexec'
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
