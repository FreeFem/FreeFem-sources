include $(abs_top_builddir)/download/headers-sparsesolver.inc
#############################################
#  FLAGS :: FREEFEM TO COMPILE SUPERLU_DIST
#
DSUPERLULIB     = $(abs_top_builddir)/download/lib/libsuperlu_dist_3.0.a

BLASDEF	     	= $(FFBLASDEF)
BLASLIB      	= @BLASLIBS@

#
#  The archiver and the flag(s) to use when building archive (library)
#  If your system has no ranlib, set RANLIB = echo.
#
ARCH         	= @AR@
ARCHFLAGS    	= @ARFLAGS@
RANLIB       	= @RANLIB@

############################################################################
# C compiler setup
CC           	= @CC@
# CFLAGS should be set to be the C flags that include optimization
CFLAGS          = @CFLAGS@ @MPI_INCLUDE@
# NOOPTS should be set to be the C flags that turn off any optimization
NOOPTS		=  @CNOFLAGS@ 
############################################################################
# FORTRAN compiler setup
FORTRAN         = @F77@
FFLAGS 		= @FFLAGS@
############################################################################
LOADER	 = @CXX@
LOADOPTS = 
############################################################################
#  C preprocessor defs for compilation (-DNoChange, -DAdd_, or -DUpCase)
#
#  Need follow the convention of how C calls a Fortran routine.
#
CDEFS        = @CFLAGSF77@
