# Makefile using Automake + Autoconf
# ----------------------------------
# $Id$
CXXLD=$(STATICTOOL) $(CXX)

bin_PROGRAMS=$(BAMGPROG)
EXTRA_PROGRAMS=bamg  cvmsh2 drawbdmesh
bamg_SOURCES=bamg.cpp 
bamg_DEPENDENCIES=../fflib/libff.a 

cvmsh2_SOURCES=cvmsh2.cpp
cvmsh2_DEPENDENCIES=../fflib/libff.a 

EXTRA_drawbdmesh_SOURCES=../Graphics/Xrgraph.cpp ../Graphics/Pcrgraph.cpp
drawbdmesh_DEPENDENCIES=../fflib/libff.a  @STD_GRAPH_OBJ@
drawbdmesh_SOURCES=drawbdmesh.cpp 
drawbdmesh_LDADD = @STD_GRAPH_OBJ@ $(LDADD)  @STD_LIBS@ @LIBSNOCONSOLE@

LDADD= ../fflib/libff.a
INCLUDES=-I$(srcdir)/../fflib -I$(srcdir)/../bamglib -I$(srcdir)/../Graphics


