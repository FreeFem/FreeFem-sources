# Makefile using Automake + Autoconf
# ----------------------------------
# $Id$

bin_PROGRAMS=$(X11PROG)
EXTRA_PROGRAMS=FreeFem++-x11
FreeFem___x11_SOURCES=../Graphics/Xrgraph.cpp
FreeFem___x11_DEPENDENCIES=../fflib/libff.a ../lglib/liblg.a
LDADD=../lglib/liblg.a ../fflib/libff.a $(LIBX11) $(LIBDL)
INCLUDES=$(INCLUDEX11)
