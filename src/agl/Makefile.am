# Makefile using Automake + Autoconf
# ----------------------------------
# $Id$

bin_PROGRAMS=$(AGLPROG)
EXTRA_PROGRAMS=FreeFem++-agl

FreeFem___mpi_SOURCES=../Graphics/macglrgraf.cpp
FreeFem___mpi_DEPENDENCIES=../fflib/libff.a ../lglib/liblg.a
FreeFem___mpi_LDADD=../lglib/liblg.a ../fflib/libff.a

INCLUDES=-I../fflib -I../Graphics -I../femlib
LIBS=-framework Carbon -framework OpenGL -framework AGL -ldl
