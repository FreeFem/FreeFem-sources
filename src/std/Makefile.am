# Makefile using Automake + Autoconf
# ----------------------------------
# $Id$

# To create statically linked executables (see configure.ac)
CXXLD=$(STATICTOOL) $(CXX)

bin_PROGRAMS=$(STDPROG)
EXTRA_PROGRAMS=FreeFem++
FreeFem___SOURCES=
EXTRA_FreeFem___SOURCES=../Graphics/Xrgraph.cpp ../Graphics/Pcrgraph.cpp
FreeFem___DEPENDENCIES=../fflib/libff.a ../lglib/liblg.a @STD_GRAPH_OBJ@
LDADD=@STD_GRAPH_OBJ@ ../lglib/liblg.a ../fflib/libff.a \
	@UMFPACKLIBS@ @ARPACKLIBS@ @BLASLIBS@ @STD_LIBS@
INCLUDES=-I../fflib
AM_LDFLAGS=@STD_LDFLAGS@
