abs_top_builddir=@abs_top_builddir@
CC      = @CC@
# FC : Fortran 90 compiler
FC      = @FC@
LIBBLAS =  @BLASLIBS@
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
MPI_INCLUDE=@MPI_INCLUDE@
MPI_LIB=@MPI_LIB@


