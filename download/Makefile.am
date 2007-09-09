# Downloading and compiling extra libraries
# -----------------------------------------

# $Id$

SUBDIRS=blas arpack umfpack fltk
EXTRA_DIST=f2c/Makefile  f2c/fort77.sed  f2c/tt.f
all-recursive: lib include

lib: 
	mkdir lib
include:
	mkdir  include

clean-local:
	-rm -rf  lib include 
