# Makefile using Automake + Autoconf
# ----------------------------------
# $Id$

bin_PROGRAMS=$(AGLPROG)
EXTRA_PROGRAMS=FreeFem++-agl${FF_PROG_EXT}${EXEEXT}

FreeFem___agl_SOURCES=../Graphics/glrgraph.hpp ../Graphics/macglrgraf.cpp
FreeFem___agl_DEPENDENCIES=../fflib/libff.a ../lglib/liblg.a
LDADD=../lglib/liblg.a ../fflib/libff.a -ldl
INCLUDES=-I../fflib -I../Graphics -I../femlib

# Automake does not want to see "-framework" in LDADD
AM_LDFLAGS=-framework Carbon -framework OpenGL -framework AGL
