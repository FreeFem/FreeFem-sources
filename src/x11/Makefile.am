# Makefile using Automake + Autoconf
# ----------------------------------
# $Id$

bin_PROGRAMS=$(X11PROG)
EXTRA_PROGRAMS=FreeFem++-x11
FreeFem___x11_SOURCES=../Graphics/Xrgraph.cpp
FreeFem___x11_DEPENDENCIES=../fflib/libff.a ../lglib/liblg.a
FreeFem___x11_LDADD=../lglib/liblg.a ../fflib/libff.a
INCLUDES= -I$(INCLUDEX11)
LIBS=$(LIBX11) -ldl
