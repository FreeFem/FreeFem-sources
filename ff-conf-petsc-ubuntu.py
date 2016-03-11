#!/usr/bin/python
if __name__ == '__main__':
  import sys
  import os
  sys.path.insert(0, os.path.abspath('config'))
  import configure
  configure_options = [
    '--CFLAGS=-O2',
    '--COPTFLAGS=-O3',
    '--CXXFLAGS=-O2 -std=c++11',
    '--CXXOPTFLAGS=-O3',
    '--FFLAGS=-O2',
    '--FOPTFLAGS=-O3',
    '--download-fftw',
    '--download-hypre',
    '--download-metis',
    '--download-ml',
    '--download-mumps',
    '--download-parmetis',
    '--download-pastix',
    '--download-ptscotch',
    '--download-scalapack',
    '--download-suitesparse',
    '--download-superlu',
    '--prefix=/usr/local/ff++/petsc',
    '--with-blas-lapack-lib=-llapack -lblas',
    'PETSC_ARCH=arch-ff++',
  ]
  configure.petsc_configure(configure_options)
