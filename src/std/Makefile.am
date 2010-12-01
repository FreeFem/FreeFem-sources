# Makefile using Automake + Autoconf
# ----------------------------------
# $Id$

# To create statically linked executables (see configure.ac)
CXXLD=$(STATICTOOL) $(CXX)

bin_PROGRAMS=$(STDPROG)
EXTRA_PROGRAMS=FreeFem++-std
FreeFem___std_SOURCES=
EXTRA_FreeFem___std_SOURCES=../Graphics/Xrgraph.cpp ../Graphics/Pcrgraph.cpp
FreeFem___std_DEPENDENCIES=../fflib/libff.a ../lglib/liblg.a @STD_GRAPH_OBJ@ 
LDADD=@STD_GRAPH_OBJ@ ../lglib/liblg.a ../fflib/libff.a \
	@UMFPACKLIBS@ @ARPACKLIBS@ @BLASLIBS@ @STD_LIBS@ @LIBSNOCONSOLE@
INCLUDES=-I../fflib
AM_LDFLAGS=@STD_LDFLAGS@
