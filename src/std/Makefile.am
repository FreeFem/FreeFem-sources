# Makefile using Automake + Autoconf
# ----------------------------------
# $Id$

bin_PROGRAMS=$(STDPROG)
EXTRA_PROGRAMS=FreeFem++${FF_PROG_EXT}${EXEEXT}
FreeFem___SOURCES=../Graphics/Xrgraph.cpp
FreeFem___DEPENDENCIES=../fflib/libff.a ../lglib/liblg.a
LDADD=../lglib/liblg.a ../fflib/libff.a $(LIBX11) -ldl
INCLUDES=-I$(INCLUDEX11)
