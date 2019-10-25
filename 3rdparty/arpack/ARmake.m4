# ARPACK ARmake.inc modified for FreeFEM
# $Id$

###########################################################################
#
#  Program:         ARPACK
#
#  Module:          ARmake.inc
#
#  Purpose:         Top-level Definitions
#
#  Creation date:   February 22, 1996
#
#  Modified:
#
#  Send bug reports, comments or suggestions to arpack@caam.rice.edu
#
############################################################################
#
# %---------------------------------%
# |  SECTION 1: PATHS AND LIBRARIES |
# %---------------------------------%
#
#
# %--------------------------------------%
# | You should change the definition of  |
# | home if ARPACK is built some place   | 
# | other than your home directory.      |
# %--------------------------------------%
#
home = FF_HOME
#
#  %--------------------------------------%
#  | The platform identifier to suffix to |
#  | the end of library names             |
#  %--------------------------------------%
#
PLAT = ff++
#
#  %------------------------------------------------------%
#  | The directories to find the various pieces of ARPACK |
#  %------------------------------------------------------%
#
BLASdir      = $(home)/BLAS
LAPACKdir    = $(home)/LAPACK
UTILdir      = $(home)/UTIL
SRCdir       = $(home)/SRC
#
#DIRS        = $(BLASdir) $(LAPACKdir) $(UTILdir) $(SRCdir)
#
# %-------------------------------------------------------------------%
# | Comment out the previous line and uncomment the following         |
# | if you already have the BLAS and LAPACK installed on your system. |
# | NOTE: ARPACK assumes the use of LAPACK version 2 codes.           |
# %-------------------------------------------------------------------%
#
DIRS         = FF_LAPACKdir  $(UTILdir) $(SRCdir)
#
# %---------------------------------------------------%
# | The name of the libraries to be created/linked to |
# %---------------------------------------------------%
#
ARPACKLIB  = FF_ARPACKLIB
LAPACKLIB  = FF_LAPACKLIB
BLASLIB = FF_BLASLIB 
#
ALIBS =  $(ARPACKLIB) $(LAPACKLIB) $(BLASLIB) 
#
# 
# %---------------------------------------------------------%
# |                  SECTION 2: COMPILERS                   |
# |                                                         |
# | The following macros specify compilers, linker/loaders, |
# | the archiver, and their options.  You need to make sure |
# | these are correct for your system.                      |
# %---------------------------------------------------------%
#
#
# %------------------------------%
# | Make our own suffixes' list. |
# %------------------------------%
#
.SUFFIXES:
.SUFFIXES:	.f	.o
#
# %------------------%
# | Default command. |
# %------------------%
#
.DEFAULT:
	@$(ECHO) "Unknown target $@, try:  make help"
#
# %-------------------------------------------%
# |  Command to build .o files from .f files. |
# %-------------------------------------------%
#
.f.o:
	@$(ECHO) Making $@ from $<
	@$(FC) -c $(FFLAGS) $<
#
# %-----------------------------------------%
# | Various compilation programs and flags. |
# | You need to make sure these are correct |
# | for your system.                        |
# %-----------------------------------------%
#
FC      = FF_FC
FFLAGS	= FF_FFLAGS
LDFLAGS = FF_LDFLAGS
SECOND_O = FF_SECOND
CD      = cd

ECHO    = echo

LN      = ln
LNFLAGS = -s

MAKE    = make

RM      = rm
RMFLAGS = -f

SHELL   = /bin/sh
#
#  %----------------------------------------------------------------%
#  | The archiver and the flag(s) to use when building an archive   |
#  | (library).  Also the ranlib routine.  If your system has no    |
#  | ranlib, set RANLIB = touch.                                    |
#  %----------------------------------------------------------------%
#
AR = FF_AR 
ARFLAGS = FF_ARFLAGS
RANLIB   = FF_RANLIB
#
# %----------------------------------%
# | This is the general help target. |
# %----------------------------------%
#
help:
	@$(ECHO) "usage: make ?"
	
.NOTPARALLEL:	
