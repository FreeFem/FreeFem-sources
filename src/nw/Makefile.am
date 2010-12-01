# Makefile using Automake + Autoconf
# ----------------------------------
# $Id$

# To create statically linked executables (see configure.ac)
CXXLD=$(STATICTOOL) $(CXX)
bin_PROGRAMS=FreeFem++ FreeFem++-nw @FFGLUTPROG@
EXTRA_PROGRAMS=ffglut
ffglut_SOURCES=../Graphics/ffglut.cpp ../Graphics/gggg.cpp ../Graphics/ffthreads.cpp \
../Graphics/ffthreads.hpp \
../femlib/fem.cpp ../femlib/Mesh3dn.cpp ../femlib/Mesh2dn.cpp  ../femlib/Mesh1dn.cpp ../femlib/GQuadTree.cpp  ../femlib/FQuadTree.cpp \
../femlib/Drawing.cpp ../femlib/mshptg.cpp  

ffglut_DEPENDENCIES= ../libMesh/libMesh.a
ffglut_LDADD= ../libMesh/libMesh.a @LIBSGLUT@ @LIBSPTHREAD@ 
FreeFem___nw_SOURCES=../Graphics/sansrgraph.cpp 
FreeFem___SOURCES=../Graphics/sansrgraph.cpp 

FreeFem___nw_DEPENDENCIES=../fflib/libff.a ../lglib/liblg.a
FreeFem___DEPENDENCIES=../fflib/libff.a ../lglib/liblg.a

LDADD=../lglib/liblg.a ../fflib/libff.a @UMFPACKLIBS@ @ARPACKLIBS@ @BLASLIBS@
INCLUDES=-I../fflib -I../Graphics -I../.. -I../femlib 
