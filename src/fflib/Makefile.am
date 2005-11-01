# Makefile using Automake + Autoconf
# ----------------------------------
# $Id$

noinst_LIBRARIES=libff.a

libff_a_SOURCES2=AFunction.cpp AFunction2.cpp  \
 array_long.cpp          array_real.cpp        array_complex.cpp \
lex.cpp lgmesh.cpp		\
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
error.hpp gmres.hpp InitFunct.hpp lex.hpp lgfem.hpp lgsolver.hpp	\
Operator.hpp problem.hpp Serialize.hpp showverb.hpp String.hpp		\
throwassert.hpp versionnumber.hpp  CodeAlloc.hpp \
array_init.hpp  array_tlp.hpp \
strversionnumber.hpp

libff_a_SOURCES=$(libff_a_SOURCES2)  strversionnumber.cpp 

# eigenvalue.cpp is optional (see configure.ac)
EXTRA_libff_a_SOURCES=../Eigen/eigenvalue.cpp
libff_a_LIBADD=@EIGENOBJ@
libff_a_DEPENDENCIES=@EIGENOBJ@

INCLUDES=-I../lglib -I../bamglib -I../Graphics -I../femlib	\
-I$(top_srcdir)/arpack/arpack++/include

# Build date changes at every change
EXTRA_DIST=strversionnumber.m4
BUILT_SOURCES=strversionnumber.cpp
strversionnumber.cpp: $(libff_a_SOURCES2) ../../configure
	m4 -DVersionFreeFemDate="`date`" strversionnumber.m4 > $@
FORCE:
