# Makefile using Automake + Autoconf
# ----------------------------------
# $Id$

noinst_LIBRARIES=libff.a

libff_a_SOURCES=AFunction.cpp AFunction2.cpp lex.cpp lgmesh.cpp \
lgmat.cpp ../femlib/Drawing.cpp	\
../femlib/gibbs.cpp ../femlib/CheckPtr.cpp ../femlib/fem.cpp		\
../femlib/QuadratureFormular.cpp ../femlib/FESpace.cpp			\
../femlib/Element_RT.cpp ../femlib/mshptg.cpp ../femlib/FQuadTree.cpp	\
../bamglib/QuadTree.cpp ../bamglib/R2.cpp ../bamglib/Meshio.cpp		\
../bamglib/Mesh2.cpp ../bamglib/Metric.cpp ../femlib/BamgFreeFem.cpp	\
../bamglib/MeshDraw.cpp ../bamglib/MeshGeom.cpp				\
../bamglib/MeshQuad.cpp ../bamglib/SetOfE4.cpp ../bamglib/MeshRead.cpp	\
../bamglib/MeshWrite.cpp problem.cpp strversionnumber.cpp		\
InitFunct.cpp ../Algo/lgalgo.cpp ../femlib/Element_P2h.cpp load.cpp	\
lgfem.cpp AFunction.hpp AnyType.hpp error.hpp gmres.hpp InitFunct.hpp	\
lex.hpp lgfem.hpp lgsolver.hpp Operator.hpp problem.hpp Serialize.hpp	\
showverb.hpp String.hpp throwassert.hpp versionnumber.hpp


# eigenvalue.cpp is optional (see configure.ac)
EXTRA_libff_a_SOURCES=../Eigen/eigenvalue.cpp
libff_a_LIBADD=@EIGENOBJ@
libff_a_DEPENDENCIES=@EIGENOBJ@

INCLUDES=-I../lglib -I../bamglib -I../Graphics -I../femlib	\
-I$(top_srcdir)/arpack/arpack++/include

# Build date changes at every build
BUILT_SOURCES=strversionnumber.cpp
strversionnumber.cpp: FORCE
	m4 -DVersionFreeFemDate="`date`" strversionnumber.m4 > $@
FORCE:
