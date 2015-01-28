# Makefile using Automake + Autoconf
# ----------------------------------
# $Id$
CXXLD=$(STATICTOOL) $(CXX)

bin_PROGRAMS=$(BAMGPROG)
EXTRA_PROGRAMS=bamg  cvmsh2 #drawbdmesh
bamg_SOURCES=bamg.cpp global.cpp
bamg_DEPENDENCIES=../fflib/Mesh2.o ../fflib/MeshDraw.o ../fflib/MeshGeom.o ../fflib/MeshQuad.o ../fflib/MeshRead.o ../fflib/MeshWrite.o  ../fflib/Meshio.o ../fflib/Metric.o ../fflib/QuadTree.o ../fflib/R2.o ../fflib/SetOfE4.o ../fflib/write_hdf5.o	../fflib/write_xdmf.o

cvmsh2_SOURCES=cvmsh2.cpp global.cpp
cvmsh2_DEPENDENCIES=../fflib/Mesh2.o ../fflib/MeshDraw.o ../fflib/MeshGeom.o ../fflib/MeshQuad.o ../fflib/MeshRead.o ../fflib/MeshWrite.o  ../fflib/Meshio.o ../fflib/Metric.o ../fflib/QuadTree.o ../fflib/R2.o ../fflib/SetOfE4.o ../fflib/write_hdf5.o	../fflib/write_xdmf.o

#EXTRA_drawbdmesh_SOURCES= # ../Graphics/Xrgraph.cpp ../Graphics/Pcrgraph.cpp 
#drawbdmesh_DEPENDENCIES=../fflib/libff.a  @STD_GRAPH_OBJ@
#drawbdmesh_SOURCES= drawbdmesh.cpp 
#drawbdmesh_LDADD = @STD_GRAPH_OBJ@ $(LDADD)  @STD_LIBS@ @LIBSNOCONSOLE@

LDADD= ../fflib/Mesh2.o ../fflib/MeshDraw.o ../fflib/MeshGeom.o ../fflib/MeshQuad.o ../fflib/MeshRead.o ../fflib/MeshWrite.o  ../fflib/Meshio.o ../fflib/Metric.o ../fflib/QuadTree.o ../fflib/R2.o ../fflib/SetOfE4.o ../fflib/write_hdf5.o	../fflib/write_xdmf.o
AM_CPPFLAGS=-I$(srcdir)/../fflib -I$(srcdir)/../bamglib -I$(srcdir)/../Graphics $(HDF5_CPPFLAGS)


