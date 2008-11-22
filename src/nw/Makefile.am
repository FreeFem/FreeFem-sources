# Makefile using Automake + Autoconf
# ----------------------------------
# $Id$

# To create statically linked executables (see configure.ac)
CXXLD=$(STATICTOOL) $(CXX)

bin_PROGRAMS=FreeFem++-nw @FFGLUTPROG@
EXTRA_PROGRAMS=ffglut
ffglut_SOURCES=../Graphics/ffglut.cpp ../Graphics/gggg.cpp ../Graphics/ffthreads.cpp \
../Graphics/ffthreads.hpp 

ffglut_DEPENDENCIES=../fflib/fem.o ../fflib/FQuadTree.o ../fflib/Drawing.o ../fflib/mshptg.o  
ffglut_LDADD= ../fflib/fem.o ../fflib/FQuadTree.o ../fflib/Drawing.o ../fflib/mshptg.o  @LIBSGLUT@ @LIBSPTHREAD@
FreeFem___nw_SOURCES=../Graphics/sansrgraph.cpp

FreeFem___nw_DEPENDENCIES=../fflib/libff.a ../lglib/liblg.a
LDADD=../lglib/liblg.a ../fflib/libff.a @UMFPACKLIBS@ @ARPACKLIBS@ @BLASLIBS@
INCLUDES=-I../fflib -I../Graphics -I../.. -I../femlib 
