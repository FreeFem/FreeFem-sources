# Ubuntu 16.04 - VM0

sudo apt install gcc g++ gfortran m4 patch git flex bison

# Ubuntu 16.04 - VM1
sudo apt install gcc g++ gfortran m4 patch git libblas-dev liblapack-dev libsuitesparse-dev libopenmpi-dev libhdf5-dev libgsl-dev flex bison

# Ubuntu 16.04 - VM2
sudo apt install gcc g++ gfortran m4 patch git libblas-dev liblapack-dev libsuitesparse-dev libopenmpi-dev libhdf5-dev libgsl-dev flex bison wget cmake autoconf automake autotools-dev

# Ubuntu 18.04 - VM0
sudo apt install gcc g++ gfortran m4 patch git flex bison

# Ubuntu 18.04 - VM1
sudo apt install gcc g++ gfortran m4 patch git libblas-dev liblapack-dev libsuitesparse-dev libopenmpi-dev libhdf5-dev libgsl-dev flex bison

# Ubuntu 18.04 - VM2
sudo apt install gcc g++ gfortran m4 patch git libblas-dev liblapack-dev libsuitesparse-dev libopenmpi-dev libhdf5-dev libgsl-dev flex bison wget cmake autoconf automake autotools-dev

# MacOS 10.10 - VM0
TODO

# MacOS 10.10 - VM1
TODO

# MacOS 10.10 - VM2A
TODO

# MacOS 10.13 - VM0
brew install gcc m4 git flex bison

# MacOS 10.13 - VM1
brew install gcc m4 git flex bison lapack suitesparse hdf5
openmpi compilation using gcc & g++

# MacOS 10.13 - VM2
brew install gcc m4 git flex bison lapack suitesparse hdf5 cmake wget autoconf automake
openmpi compilation using gcc & g++

# Windows 7
mingw-w64-x86_64-libmicroutils
