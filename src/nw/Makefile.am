# Makefile using Automake + Autoconf
# ----------------------------------
# $Id$

bin_PROGRAMS=FreeFem++-nw
FreeFem___nw_SOURCES=../Graphics/sansrgraph.cpp
FreeFem___nw_DEPENDENCIES=../fflib/libff.a ../lglib/liblg.a
LDADD=../lglib/liblg.a ../fflib/libff.a $(LIBDL)
INCLUDES=-I../fflib -I../Graphics
