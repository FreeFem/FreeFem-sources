# Makefile using Automake + Autoconf
# ----------------------------------
# $Id$

bin_PROGRAMS=$(STDPROG)
EXTRA_PROGRAMS=FreeFem++
FreeFem___SOURCES=../Graphics/Xrgraph.cpp
FreeFem___DEPENDENCIES=../fflib/libff.a ../lglib/liblg.a
FreeFem___LDADD=../lglib/liblg.a ../fflib/libff.a
INCLUDES=-I$(INCLUDEX11)
LIBS=$(LIBX11) -ldl
