# Makefile using Automake + Autoconf
# ----------------------------------
# $Id$

bin_PROGRAMS=$(GLXPROG)
EXTRA_PROGRAMS=FreeFem++-glx
FreeFem___glx_SOURCES=../Graphics/xglrgraf.cpp
FreeFem___glx_DEPENDENCIES=../fflib/libff.a ../lglib/liblg.a
FreeFem___glx_LDADD=../lglib/liblg.a ../fflib/libff.a
INCLUDES=-I../fflib
LIBS+=-lGL -lGLU -ldl
