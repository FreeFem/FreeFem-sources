abs_top_builddir=@abs_top_builddir@
include $(abs_top_builddir)/download/headers-sparsesolver.inc
DOWNLOADFF= $(abs_top_builddir)/download
###############################################################################
#
#  This file is part of MUMPS 4.9.2, built on Thu Nov  5 07:05:08 UTC 2009
#
################################################################################
#
#   Makefile.inc.generic
#
#   This defines some parameters dependent on your platform; you should
#   look for the approriate file in the directory ./Make.inc/ and copy it
#   into a file called Makefile.inc. For example, from the MUMPS root
#   directory, use 
#   "cp Make.inc/Makefile.inc.generic ./Makefile.inc"
#   (see the main README file for details)
#
#   If you do not find any suitable Makefile in Makefile.inc, use this file:
#   "cp Make.inc/Makefile.inc.generic ./Makefile.inc" and modify it according
#   to the comments given below. If you manage to build MUMPS on a new platform,
#   and think that this could be useful to others, you may want to send us
#   the corresponding Makefile.inc file.
#
################################################################################


########################################################################
#Begin orderings
#
# NOTE that PORD is distributed within MUMPS by default. If you would like to
# use other orderings, you need to obtain the corresponding package and modify
# the variables below accordingly.
# For example, to have Metis available within MUMPS:
#          1/ download Metis and compile it
#          2/ uncomment (suppress # in first column) lines
#             starting with LMETISDIR,  LMETIS
#          3/ add -Dmetis in line ORDERINGSF
#             ORDERINGSF  = -Dpord -Dmetis
#          4/ Compile and install MUMPS
#             make clean; make   (to clean up previous installation)
#
#          Metis/ParMetis and SCOTCH/PT-SCOTCH (ver 5.1 and later) orderings are now available for MUMPS.
#

SCOTCHDIR  = $(DOWNLOADFF)/
ISCOTCH    = -I$(SCOTCHDIR)/include/scotch
# You have to choose one among the following two lines depending on
# the type of analysis you want to perform. If you want to perform only
# sequential analysis choose the first (remember to add -Dscotch in the ORDERINGSF
# variable below); for both parallel and sequential analysis choose the second 
# line (remember to add -Dptscotch in the ORDERINGSF variable below)

#LSCOTCH    = -L$(SCOTCHDIR)/lib/scotch/ -lesmumps -lscotch -lscotcherr
LSCOTCH    = -L$(SCOTCHDIR)/lib/scotch/ -lptesmumps -lptscotch -lptscotcherr


LPORDDIR = $(topdir)/PORD/lib/
IPORD    = -I$(topdir)/PORD/include/
LPORD    = -L$(LPORDDIR) -lpord

LMETISDIR = $(DOWNLOADFF)
#IMETIS    = # Metis doesn't need include files (Fortran interface avail.)

# You have to choose one among the following two lines depending on
# the type of analysis you want to perform. If you want to perform only
# sequential analysis choose the first (remember to add -Dmetis in the ORDERINGSF
# variable below); for both parallel and sequential analysis choose the second 
# line (remember to add -Dparmetis in the ORDERINGSF variable below)

#LMETIS    = -L$(LMETISDIR) -lmetis
LMETIS    = -L$(LMETISDIR)/lib/parmetis -lparmetis -lmetis

# The following variables will be used in the compilation process.
# Please note that -Dptscotch and -Dparmetis imply -Dscotch and -Dmetis respectively.
ORDERINGSF = -Dscotch -Dmetis -Dpord -Dptscotch -Dparmetis
#ORDERINGSF  = -Dpord
ORDERINGSC  = $(ORDERINGSF)

LORDERINGS = $(LMETIS) $(LPORD) $(LSCOTCH)
IORDERINGSF = $(ISCOTCH)
IORDERINGSC = $(IMETIS) $(IPORD) $(ISCOTCH)

#End orderings
########################################################################
# DEFINE HERE SOME COMMON COMMANDS, THE COMPILER NAMES, ETC...

# PLAT : use it to add a default suffix to the generated libraries
PLAT    = # FREEFEM
# RM : remove files
RM      = /bin/rm -f
# CC : C compiler
CC      = @MPICC@
# FC : Fortran 90 compiler
FC      = @MPIFC@
# FL : Fortran linker
FL      = @MPIFC@
# AR : Archive object in a library
AR      = @AR@ @ARFLAGS@ #ar vr
# RANLIB : generate index of an archive file
#   (optionnal use "RANLIB = echo" in case of problem)
RANLIB  = @RANLIB@
#RANLIB  = echo

# SCALAP should define the SCALAPACK and  BLACS libraries.
SCALAPDIR = $(DOWNLOADFF)/lib/scalapack/
SCALAP    = $(FFBLASLIB) -L$(SCALAPDIR) -lscalapack -L$(DOWNLOADFF)/lib/blacs/ -lblacs_MPI-$(PLAT)-0 -lblacsF77init_MPI-$(PLAT)-0 -lblacsCinit_MPI-$(PLAT)-0

# INCLUDE DIRECTORY FOR MPI
INCPAR  = $(FFMPIINCLUDE)

# LIBRARIES USED BY THE PARALLEL VERSION OF MUMPS: $(SCALAP) and MPI
LIBPAR  = $(SCALAP) $(FFMPILIB)

# The parallel version is not concerned by the next two lines.
# They are related to the sequential library provided by MUMPS,
# to use instead of ScaLAPACK and MPI.
INCSEQ  = -I$(topdir)/libseq
LIBSEQ  = -L$(topdir)/libseq -lmpiseq

# DEFINE HERE YOUR BLAS LIBRARY

LIBBLAS =  $(FFBLASLIB)

# DEFINE YOUR PTHREAD LIBRARY
LIBOTHERS = @LIBSPTHREAD@ @FLIBS@ $(FFLIBOTHERSMUMPS) 

# FORTRAN/C COMPATIBILITY:
#  Use:
#    -DAdd_ if your Fortran compiler adds an underscore at the end
#              of symbols,
#     -DAdd__ if your Fortran compiler adds 2 underscores,
#
#     -DUPPER if your Fortran compiler uses uppercase symbols
#
#     leave empty if your Fortran compiler does not change the symbols.
#

CDEFS = @CFLAGSF77@

#COMPILER OPTIONS
OPTF    = -O3
OPTC    = -O3 -I.
OPTL    = -O3

# CHOOSE BETWEEN USING THE SEQUENTIAL OR THE PARALLEL VERSION.

#Sequential:
#INC = $(INCSEQ)
#LIB = $(LIBSEQ)
#LIBSEQNEEDED = libseqneeded

#Parallel:
INC = $(INCPAR)
LIB = $(LIBPAR)
LIBSEQNEEDED = 

