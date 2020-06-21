############################################################################
#
#  Program:         ScaLAPACK
#
#  Module:          SLmake.inc
#
#  Purpose:         Top-level Definitions
#
#  Creation date:   February 15, 2000
#
#  Modified:        October 13, 2011
#
#  Send bug reports, comments or suggestions to scalapack@cs.utk.edu
#
############################################################################
#
#  C preprocessor definitions:  set CDEFS to one of the following:
#
#     -DNoChange (fortran subprogram names are lower case without any suffix)
#     -DUpCase   (fortran subprogram names are upper case without any suffix)
#     -DAdd_     (fortran subprogram names are lower case with "_" appended)

CDEFS         = @CFLAGSF77@ -DNO_IEEE $(USEMPI)

#
#  The fortran and C compilers, loaders, and their flags
#

FC            = @MPIFC@
CC            = @MPICC@
NOOPT         = -O0 @CNOFLAGS@

# FFCS - some return statements without value cause trouble on MacOS
# FFCS - add path to mpi.h (required for MacOS 10.8 + MacPorts OpenMPI)
# FFCS - added @CNOFLAGS@ according to upstream changes
CCFLAGS       = -O3 -Wreturn-type @CFLAGS@ -I'@MPI_INC_DIR@' @CNOFLAGS@
FCFLAGS       = -O3 @CNOFLAGS@
FCLOADER      = $(FC)
CCLOADER      = $(CC)
FCLOADFLAGS   = $(FCFLAGS)
CCLOADFLAGS   = $(CCFLAGS)

#
#  The archiver and the flag(s) to use when building archive (library)
#  Also the ranlib routine.  If your system has no ranlib, set RANLIB = echo
#

ARCH          = @AR@
ARCHFLAGS     = @ARFLAGS@
RANLIB        = @RANLIB@

#
#  The name of the ScaLAPACK library to be created
#

SCALAPACKLIB  = libscalapack.a

#
#  BLAS, LAPACK (and possibly other) libraries needed for linking test programs
#

BLASLIB       = $(FFBLASLIB)
LAPACKLIB     = $(FFLAPACKLIB)
LIBS          = $(LAPACKLIB) $(BLASLIB)

WGET = @WGET@

.NOTPARRALEL:
