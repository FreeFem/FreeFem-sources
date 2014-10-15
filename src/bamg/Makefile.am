# Makefile using Automake + Autoconf
# ----------------------------------
# $Id$
CXXLD=$(STATICTOOL) $(CXX)

bin_PROGRAMS=$(BAMGPROG)
EXTRA_PROGRAMS=bamg  cvmsh2 #drawbdmesh
bamg_SOURCES=bamg.cpp global.cpp

# ALH - 7/5/14 - dependencies directory changed to follow the new Automake policy
# [[file:../../configure.ac::subdir-objects]]

bamg_DEPENDENCIES=../bamglib/Mesh2.o ../bamglib/MeshDraw.o ../bamglib/MeshGeom.o ../bamglib/MeshQuad.o		\
../bamglib/MeshRead.o ../bamglib/MeshWrite.o ../bamglib/Meshio.o ../bamglib/Metric.o ../bamglib/QuadTree.o	\
../bamglib/R2.o ../bamglib/SetOfE4.o

cvmsh2_SOURCES=cvmsh2.cpp global.cpp
cvmsh2_DEPENDENCIES=../bamglib/Mesh2.o ../bamglib/MeshDraw.o ../bamglib/MeshGeom.o ../bamglib/MeshQuad.o	\
../bamglib/MeshRead.o ../bamglib/MeshWrite.o ../bamglib/Meshio.o ../bamglib/Metric.o ../bamglib/QuadTree.o	\
../bamglib/R2.o ../bamglib/SetOfE4.o

#EXTRA_drawbdmesh_SOURCES= # ../Graphics/Xrgraph.cpp ../Graphics/Pcrgraph.cpp 
#drawbdmesh_DEPENDENCIES=../bamglib/libff.a  @STD_GRAPH_OBJ@
#drawbdmesh_SOURCES= drawbdmesh.cpp 
#drawbdmesh_LDADD = @STD_GRAPH_OBJ@ $(LDADD)  @STD_LIBS@ @LIBSNOCONSOLE@

LDADD= ../bamglib/Mesh2.o ../bamglib/MeshDraw.o ../bamglib/MeshGeom.o ../bamglib/MeshQuad.o ../bamglib/MeshRead.o	\
../bamglib/MeshWrite.o ../bamglib/Meshio.o ../bamglib/Metric.o ../bamglib/QuadTree.o ../bamglib/R2.o			\
../bamglib/SetOfE4.o
AM_CPPFLAGS=-I$(srcdir)/../fflib -I$(srcdir)/../bamglib -I$(srcdir)/../Graphics


