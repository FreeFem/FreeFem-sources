############################################################################
#
#  Program:         SuperLU
#
#  Module:          make.inc
#
#  Purpose:         Top-level Definitions
#
#  Creation date:   October 2, 1995
#
#  Modified:	    February 4, 1997  Version 1.0
#		    November 15, 1997 Version 1.1
#		    September 1, 1999 Version 2.0
#
############################################################################
#
#
#  The name of the libraries to be created/linked to
#
TMGLIB       = $(INSTALL)/../lib/libtmglib.a
SUPERLULIB   = $(INSTALL)/../lib/libsuperlu_3.1.a
BLASLIB      =  $(BLASLIBS)

#
#  The archiver and the flag(s) to use when building archive (library)
#  If your system has no ranlib, set RANLIB = echo.
#
ARCH         = $(AR)
ARCHFLAGS    = $(ARFLAGS)


FORTRAN	     = f77
FFLAGS       = -O
LOADER       = $(CC)
LOADOPTS     = 

#
#  C preprocessor defs for compilation (-DNoChange, -DAdd_, or -DUpCase)
#
CDEFS        = -DAdd_
#
# The directory in which Matlab is installed
#
MATLAB	     = /usr/sww/pkg/matlab
