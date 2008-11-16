# Downloading and compiling extra libraries
# -----------------------------------------

# $Id$

SUBDIRS=blas arpack umfpack fltk 
EXTRA_DIST=f2c/Makefile  f2c/fort77.sed	f2c/Makefile-MacOs tetgen/Makefile tetgen/tetgen1.4.2.patch
all-recursive: lib include pkg

lib: 
	mkdir lib
include:
	mkdir  include
pkg:
	mkdir pkg

clean-local:
	-rm -rf  lib include 
