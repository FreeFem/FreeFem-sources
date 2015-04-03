# Makefile using Automake + Autoconf
# ----------------------------------
# $Id$

# To create statically linked executables (see configure.ac)
CXXLD=$(STATICTOOL) $(CXX)

bin_PROGRAMS=$(GLXPROG)
EXTRA_PROGRAMS=FreeFem++-glx

# FFCS:visualization stream redirection
FreeFem___glx_SOURCES=../Graphics/xglrgraf.cpp  ../mpi/parallelempi-empty.cpp ../fflib/ffapi.cpp
FreeFem___glx_DEPENDENCIES=../fflib/libff.a ../lglib/liblg.a
LDADD=../lglib/liblg.a ../fflib/libff.a @UMFPACKLIBS@ @ARPACKLIBS@	\
	@BLASLIBS@ -lGLU -lGL @X11_LIBS@
AM_CPPFLAGS=-I$(srcdir)/../fflib
