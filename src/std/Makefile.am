# Makefile using Automake + Autoconf
# ----------------------------------
# $Id$

bin_PROGRAMS=FreeFem++
FreeFem___SOURCES=
EXTRA_FreeFem___SOURCES=../Graphics/Xrgraph.cpp ../Graphics/Pcrgraph.cpp
FreeFem___DEPENDENCIES=../fflib/libff.a ../lglib/liblg.a @STD_GRAPH_OBJ@
LDADD=@STD_GRAPH_OBJ@ ../lglib/liblg.a ../fflib/libff.a @STD_X11_LIB@ $(LIBDL)
INCLUDES=-I../fflib @STD_X11_INC@
AM_LDFLAGS=@STD_LDFLAGS@
