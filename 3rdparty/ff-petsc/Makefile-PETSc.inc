abs_top_builddir=@abs_top_builddir@
FF_prefix_petsc=@FF_prefix_petsc@
FF_generic_petsc=@FF_generic_petsc@
CC=@CC@
CXX=@CXX@
# FC : Fortran 90 compiler
FC=@FC@
BLASINC=@BLASINC@
FFCMAKE=@ff_cmake@
BLASLIBS=@BLASLIBS@
LAPACKLIBS=@LAPACKLIBS@
COMPILE_OPENBLAS=@COMPILE_OPENBLAS@
#  Use:
#    -DAdd_ if your Fortran compiler adds an underscore at the end
#              of symbols,
#     -DAdd__ if your Fortran compiler adds 2 underscores,
#
#     -DUPPER if your Fortran compiler uses uppercase symbols
#
#     leave empty if your Fortran compiler does not change the symbols.
#
CFLAGS=@CFLAGS@
FCFLAGS=@FCFLAGS@
CFLAGSF77=@CFLAGSF77@
MPI_INCLUDE=@MPI_INCLUDE@
MPI_INC_DIR=@MPI_INC_DIR@
MPI_LIB=@MPI_LIB@
MPI_LIBC=@MPI_LIBC@
MPI_LIBFC=@MPI_LIBFC@
MPI_LIB_DIRS=@MPI_LIB_DIRS@
MPICC=@MPICC@
MPICXX=@MPICXX@
MPIFC=@MPIFC@
MPIPROG=@MPIPROG@
MPIRUN=@MPIRUN@
MPISCRIPT=@MPISCRIPT@
prefix=@prefix@
WIN32DLLTARGET=@WIN32DLLTARGET@

