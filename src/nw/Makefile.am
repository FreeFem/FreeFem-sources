# Makefile using Automake + Autoconf
# ----------------------------------
# $Id$

# To create statically linked executables (see configure.ac)
CXXLD=$(STATICTOOL) $(CXX)

bin_PROGRAMS=FreeFem++-nw
FreeFem___nw_SOURCES=../Graphics/sansrgraph.cpp
FreeFem___nw_DEPENDENCIES=../fflib/libff.a ../lglib/liblg.a
LDADD=../lglib/liblg.a ../fflib/libff.a @UMFPACKLIB@ @ARPACKLIB@ @BLASLIB@
INCLUDES=-I../fflib -I../Graphics
