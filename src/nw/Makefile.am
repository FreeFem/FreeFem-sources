# Makefile using Automake + Autoconf
# ----------------------------------
# $Id$

bin_PROGRAMS=FreeFem++-nw${FF_PROG_EXT}${EXEEXT}
FreeFem___nw_SOURCES=../Graphics/sansrgraph.cpp
FreeFem___nw_DEPENDENCIES=../fflib/libff.a ../lglib/liblg.a
LDADD=../lglib/liblg.a ../fflib/libff.a -ldl
INCLUDES=-I../fflib -I../Graphics
