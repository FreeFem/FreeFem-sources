include $(abs_top_builddir)/download/hearders_sparsesolver.inc
#################
WGET = @WGET@
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
#  Modified:
#
#  Send bug reports, comments or suggestions to scalapack@cs.utk.edu
#
############################################################################
#
SHELL         = /bin/sh
#
#  The complete path to the top level of ScaLAPACK directory, usually
#  $(HOME)/SCALAPACK
#
home          = $(abs_top_builddir)/download/scalapack/scalapack-1.8.0/
#
#  The platform identifier to suffix to the end of library names
#
PLAT          = FREEFEM
#
#  BLACS setup.  All version need the debug level (0 or 1),
#  and the directory where the BLACS libraries are
#
BLACSDBGLVL   = 0
BLACSdir      = $(home)/../../lib/blacs/
#
#  MPI setup; tailor to your system if using MPIBLACS
#  Will need to comment out these 6 lines if using PVM
#
USEMPI        = -DUsingMpiBlacs
SMPLIB        = #/usr/local/mpich-1.2.1/lib/libmpich.a
BLACSFINIT    = $(BLACSdir)/libblacsF77init_MPI-$(PLAT)-0.a
BLACSCINIT    = $(BLACSdir)/libblacsCinit_MPI-$(PLAT)-0.a
BLACSLIB      = $(BLACSdir)/libblacs_MPI-$(PLAT)-0.a
TESTINGdir    = $(home)/TESTING

#
#  PVMBLACS setup, uncomment next 6 lines if using PVM
#
#USEMPI        =
#SMPLIB        = $(PVM_ROOT)/lib/$(PLAT)/libpvm3.a
#BLACSFINIT    =
#BLACSCINIT    =
#BLACSLIB      = $(BLACSdir)/blacs_PVM-$(PLAT)-$(BLACSDBGLVL).a
#TESTINGdir    = $(HOME)/pvm3/bin/$(PLAT)

CBLACSLIB     = $(BLACSCINIT) $(BLACSLIB) $(BLACSCINIT)
FBLACSLIB     = $(BLACSFINIT) $(BLACSLIB) $(BLACSFINIT)

#
#  The directories to find the various pieces of ScaLapack
#
PBLASdir      = $(home)/PBLAS
SRCdir        = $(home)/SRC
TESTdir       = $(home)/TESTING
PBLASTSTdir   = $(TESTINGdir)
TOOLSdir      = $(home)/TOOLS
REDISTdir     = $(home)/REDIST
REDISTTSTdir  = $(TESTINGdir)
#
#  The fortran and C compilers, loaders, and their flags
#
F77           = @MPIF77@ #--prefix /Users/morice/librairie/openmpi-gcc-gfortran-4.4
CC            = @MPICC@ #--prefix /Users/morice/librairie/openmpi-gcc-gfortran-4.4
NOOPT         = 
F77FLAGS      = -O3 @FFLAGS@ $(NOOPT)
CCFLAGS       = -O3 @CFLAGS@
SRCFLAG       =
F77LOADER     = $(F77)
CCLOADER      = $(CC)
F77LOADFLAGS  = 
CCLOADFLAGS   = 
#
#  C preprocessor defs for compilation 
#  (-DNoChange, -DAdd_, -DUpCase, or -Df77IsF2C)
#  Remarque: on peut utilise -D_Add -DAdd__
CDEFS         = @CFLAGSF77@ -DNO_IEEE $(USEMPI)
#
#  The archiver and the flag(s) to use when building archive (library)
#  Also the ranlib routine.  If your system has no ranlib, set RANLIB = echo
#
ARCH          = @AR@       #ar
ARCHFLAGS     = @ARFLAGS@  #cr
RANLIB        = @RANLIB@
#
#  The name of the libraries to be created/linked to
#
SCALAPACKLIB  = $(home)/libscalapack.a
BLASLIB       = $(FFBLASLIB)
LAPACKLIB     = $(FFLAPACKLIB)
#
PBLIBS        = $(SCALAPACKLIB) $(FBLACSLIB) $(LAPACKLIB) $(BLASLIB) $(SMPLIB)
PRLIBS        = $(SCALAPACKLIB) $(CBLACSLIB) $(SMPLIB)
RLIBS         = $(SCALAPACKLIB) $(FBLACSLIB) $(CBLACSLIB) $(LAPACKLIB) $(BLASLIB) $(SMPLIB)
LIBS          = $(PBLIBS)
