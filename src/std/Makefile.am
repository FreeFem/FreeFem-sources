# Makefile using Automake + Autoconf
# ----------------------------------
# $Id$

bin_PROGRAMS=$(STDPROG)
EXTRA_PROGRAMS=FreeFem++
FreeFem___SOURCES=../Graphics/Xrgraph.cpp
FreeFem___DEPENDENCIES=../fflib/libff.a ../lglib/liblg.a
LDADD=../lglib/liblg.a ../fflib/libff.a $(LIBX11) $(LIBDL)
INCLUDES=-I$(INCLUDEX11)
