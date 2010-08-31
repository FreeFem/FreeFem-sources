# Makefile using Automake + Autoconf
# ----------------------------------
# $Id$

noinst_LIBRARIES=libff.a

libff_a_SOURCES2=  UMFPack_Solver.cpp \
AFunction.cpp AFunction2.cpp  \
 array_long.cpp          array_real.cpp        array_complex.cpp \
lex.cpp lgmesh.cpp	lgmesh3.cpp 	\
CodeAlloc.cpp lgmat.cpp global.cpp 					\
../femlib/Drawing.cpp ../femlib/gibbs.cpp				\
../femlib/CheckPtr.cpp ../femlib/fem.cpp				\
../femlib/QuadratureFormular.cpp ../femlib/FESpace.cpp			\
../femlib/Element_RT.cpp ../femlib/mshptg.cpp ../femlib/FQuadTree.cpp	\
../bamglib/QuadTree.cpp ../bamglib/R2.cpp ../bamglib/Meshio.cpp		\
../bamglib/Mesh2.cpp ../bamglib/Metric.cpp ../femlib/BamgFreeFem.cpp	\
../bamglib/MeshDraw.cpp ../bamglib/MeshGeom.cpp				\
../bamglib/MeshQuad.cpp ../bamglib/SetOfE4.cpp ../bamglib/MeshRead.cpp	\
../bamglib/MeshWrite.cpp problem.cpp mt19937ar.cpp \
../Graphics/DefColor.cpp \
 InitFunct.cpp ../Algo/lgalgo.cpp			\
../femlib/Element_P2h.cpp load.cpp lgfem.cpp AFunction.hpp AnyType.hpp	\
error.hpp ../femlib/gmres.hpp InitFunct.hpp lex.hpp lgfem.hpp lgmesh3.hpp \
lgsolver.hpp	\
Operator.hpp problem.hpp Serialize.hpp showverb.hpp String.hpp		\
throwassert.hpp versionnumber.hpp  CodeAlloc.hpp \
array_init.hpp  array_tlp.hpp array_resize.hpp \
strversionnumber.hpp ffstack.hpp AddNewFE.h \
environment.cpp environment.hpp string_def.cpp \
../femlib/FESpacen.cpp \
../femlib/P012_1d.cpp \
../femlib/P012_2d.cpp \
../femlib/P012_3d.cpp \
../femlib/Mesh1dn.cpp \
../femlib/Mesh2dn.cpp \
../femlib/Mesh3dn.cpp \
../femlib/GQuadTree.cpp \
../femlib/libmesh5.c \
glumesh2D.cpp \
PlotStream.hpp \
endian.hpp \
ff++.hpp \
../Eigen/arpackff.hpp \
../femlib/splitsimplex.cpp AFunction_ext.hpp

libff_a_SOURCES=$(libff_a_SOURCES2)  strversionnumber.cpp 

# eigenvalue.cpp is optional (see configure.ac)
EXTRA_libff_a_SOURCES=../Eigen/eigenvalue.cpp
libff_a_LIBADD=@EIGENOBJ@
libff_a_DEPENDENCIES=@EIGENOBJ@

INCLUDES=-I../lglib -I../bamglib -I../Graphics -I../femlib
# -I$(top_srcdir)/arpack/arpack++/include

# Build date changes at every change
EXTRA_DIST=strversionnumber.m4
BUILT_SOURCES=strversionnumber.cpp
strversionnumber.cpp: $(libff_a_SOURCES2) ../../configure
	m4 -DVersionFreeFemDate="`date`" strversionnumber.m4 > $@
FORCE:
