# Makefile using Automake + Autoconf
# ----------------------------------
# $Id$

bin_PROGRAMS=FreeFem++
FreeFem___SOURCES=../Graphics/Xrgraph.cpp
FreeFem___DEPENDENCIES=../fflib/libff.a ../lglib/liblg.a
FreeFem___LDADD=../lglib/liblg.a ../fflib/libff.a

INCLUDEX11 = /usr/X11R6/include
INCLUDES=-I$(INCLUDEX11)

LIBX11= -L/usr/X11R6/lib -lX11 -lm 
LIBS=$(LIBX11) -ldl
