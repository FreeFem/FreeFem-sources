# Makefile using Automake + Autoconf
# ----------------------------------
# $Id$

bin_PROGRAMS=$(AGLPROG)
EXTRA_PROGRAMS=FreeFem++-agl

FreeFem___agl_SOURCES=../Graphics/glrgraph.hpp ../Graphics/macglrgraf.cpp
FreeFem___agl_DEPENDENCIES=../fflib/libff.a ../lglib/liblg.a
FreeFem___agl_LDADD=../lglib/liblg.a ../fflib/libff.a

INCLUDES=-I../fflib -I../Graphics -I../femlib
LIBS=-framework Carbon -framework OpenGL -framework AGL -ldl
