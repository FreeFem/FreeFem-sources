# Makefile using Automake + Autoconf
# ----------------------------------
# $Id$

# To create statically linked executables (see configure.ac)
CXXLD=$(STATICTOOL) $(CXX)

bin_PROGRAMS=$(GLXPROG)
EXTRA_PROGRAMS=FreeFem++-glx
FreeFem___glx_SOURCES=../Graphics/xglrgraf.cpp
FreeFem___glx_DEPENDENCIES=../fflib/libff.a ../lglib/liblg.a
LDADD=../lglib/liblg.a ../fflib/libff.a @UMFPACKLIB@ @ARPACKLIB@	\
	@BLASLIB@ -lGLU -lGL
INCLUDES=-I../fflib
