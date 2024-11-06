############################################################################
# This file is part of FreeFEM.                                            #
#                                                                          #
# FreeFEM is free software: you can redistribute it and/or modify          #
# it under the terms of the GNU Lesser General Public License as           #
# published by the Free Software Foundation, either version 3 of           #
# the License, or (at your option) any later version.                      #
#                                                                          #
# FreeFEM is distributed in the hope that it will be useful,               #
# but WITHOUT ANY WARRANTY; without even the implied warranty of           #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
# GNU Lesser General Public License for more details.                      #
#                                                                          #
# You should have received a copy of the GNU Lesser General Public License #
# along with FreeFEM. If not, see <http://www.gnu.org/licenses/>.          #
############################################################################
# SUMMARY : Makefile for FreeFem++
# LICENSE : LGPLv3
# ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
# AUTHORS : ...
# E-MAIL  : ...

noinst_LIBRARIES=libff.a

# FFCS - 28/11/11 - ffapi.cpp cannot be part of libff.a because it
# needs to be compiled with different options depending on the
# executable it is included in (eg with/without MPI)

libff_a_SOURCES2=  UMFPack_Solver.cpp \
	AFunction.cpp AFunction2.cpp \
	array_long.cpp array_real.cpp array_complex.cpp \
	lex.cpp lgmesh.cpp lgmesh3.cpp \
	CodeAlloc.cpp lgmat.cpp global.cpp \
	../femlib/Drawing.cpp ../femlib/gibbs.cpp \
	../femlib/CheckPtr.cpp ../femlib/fem.cpp \
	../femlib/QuadratureFormular.cpp ../femlib/FESpace.cpp \
	../femlib/Element_RT.cpp ../femlib/mshptg.cpp ../femlib/FQuadTree.cpp \
	../bamglib/QuadTree.cpp ../bamglib/R2.cpp ../bamglib/Meshio.cpp \
	../bamglib/Mesh2.cpp ../bamglib/Metric.cpp ../femlib/BamgFreeFem.cpp \
	../bamglib/MeshDraw.cpp ../bamglib/MeshGeom.cpp \
	../bamglib/MeshQuad.cpp ../bamglib/SetOfE4.cpp ../bamglib/MeshRead.cpp \
	../bamglib/write_hdf5.cpp ../bamglib/write_hdf5.hpp ../bamglib/write_xdmf.cpp \
	../bamglib/MeshWrite.cpp problem.cpp mt19937ar.cpp \
	../Graphics/DefColor.cpp \
	InitFunct.cpp ../Algo/lgalgo.cpp \
	../femlib/Element_P2h.cpp load.cpp lgfem.cpp AFunction.hpp AnyType.hpp \
	error.hpp ../femlib/gmres.hpp InitFunct.hpp lex.hpp lgfem.hpp lgmesh3.hpp \
	lgmesh.hpp lgsolver.hpp \
	Operator.hpp problem.hpp Serialize.hpp showverb.hpp String.hpp \
	throwassert.hpp versionnumber.hpp  CodeAlloc.hpp \
	array_init.hpp  array_tlp.hpp array_resize.hpp \
	strversionnumber.hpp ffstack.hpp AddNewFE.h \
	environment.cpp environment.hpp string_def.cpp \
	../femlib/FESpacen.cpp \
	../femlib/P012_1d.cpp \
	../femlib/P012_2d.cpp \
	../femlib/P012_3d.cpp \
	../femlib/P012_3dSurf.cpp \
	../femlib/P012_3dCurve.cpp \
	../femlib/Mesh1dn.cpp \
	../femlib/Mesh2dn.cpp \
	../femlib/Mesh3dn.cpp \
	../femlib/MeshSn.cpp \
	../femlib/MeshLn.cpp \
	../femlib/GQuadTree.cpp \
	../femlib/libmeshb7.c \
	glumesh2D.cpp \
	PlotStream.hpp \
	endian.hpp \
	ff++.hpp \
	../Eigen/arpackff.hpp \
	../femlib/splitsimplex.cpp AFunction_ext.hpp \
	ffapi.hpp P1IsoValue.cpp P1IsoValue.hpp \
	../femlib/HashMatrix.hpp \
	../femlib/HashMatrix.cpp \
	../femlib/VirtualMatrix.hpp \
	../femlib/SparseLinearSolver.hpp \
	../femlib/VirtualSolverCG.hpp \
	../femlib/VirtualSolver.hpp \
	../femlib/MatriceCreuse_tpl.hpp \
	../femlib/MatriceElementaire.hpp \
	../femlib/CG.hpp \
	../femlib/CG.cpp \
	../femlib/SkyLineSolver.hpp \
	../femlib/SparseLinearSolver.cpp \
	../femlib/SolverSkyLine.hpp \
	../femlib/VirtualSolverSparseSuite.hpp \
	../femlib/VirtualSolverSkyLine.hpp \
	msh3.cpp msh3.hpp renumb.hpp

libff_a_SOURCES=$(libff_a_SOURCES2)  strversionnumber.cpp

# eigenvalue.cpp is optional (see [[file:../../configure.ac::EIGENOBJ]])
if ENABLE_EIGEN
EXTRA_libff_a_SOURCES=../Eigen/eigenvalue.cpp
libff_a_LIBADD=../Eigen/@EIGENOBJ@
libff_a_DEPENDENCIES=../Eigen/@EIGENOBJ@
endif

AM_CPPFLAGS=-I$(srcdir)/../lglib -I$(srcdir)/../bamglib -I$(srcdir)/../Graphics -I$(srcdir)/../femlib $(HDF5_CPPFLAGS) $(UMFPACK_CPPFLAGS)
if ENABLE_EIGEN
AM_CPPFLAGS+=-I$(srcdir)/../Eigen
endif
# -I$(top_srcdir)/arpack/arpack++/include

# Build date changes at every change
EXTRA_DIST=strversionnumber.m4
BUILT_SOURCES=strversionnumber.cpp
strversionnumber.cpp: $(libff_a_SOURCES2) ../../Makefile 
	m4 -DVersionFreeFemDate="`date`" -DGitVersion="`git describe --tags 2>/dev/null|| echo 'no git'`" strversionnumber.m4 > $@
FORCE:
