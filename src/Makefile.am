# Makefile using Automake + Autoconf
# ----------------------------------
# $Id$

# Debug versions
noinst_PROGRAMS=FreeFem++-g FreeFem++-nw-g FreeFem++-x11-g	\
FreeFem++-glx-g

# Optimized versions
if OPTIM
bin_PROGRAMS=FreeFem++ FreeFem++-nw FreeFem++-x11 FreeFem++-glx
endif

# Parallel versions
if WITHMPI
noinst_PROGRAMS+=FreeFem++-mpi-g
if OPTIM
bin_PROGRAMS+=FreeFem++-mpi
endif
endif

BASESOURCES=lg.tab.cpp versionnumber.hpp AFunction.cpp lex.cpp		\
lgmesh.cpp femlib/Drawing.cpp femlib/gibbs.cpp femlib/CheckPtr.cpp	\
femlib/fem.cpp femlib/QuadratureFormular.cpp femlib/FESpace.cpp		\
femlib/Element_RT.cpp femlib/mshptg.cpp femlib/FQuadTree.cpp		\
bamglib/QuadTree.cpp bamglib/R2.cpp bamglib/Meshio.cpp			\
bamglib/Mesh2.cpp bamglib/Metric.cpp femlib/BamgFreeFem.cpp		\
bamglib/MeshDraw.cpp bamglib/MeshGeom.cpp bamglib/MeshQuad.cpp		\
bamglib/SetOfE4.cpp bamglib/MeshRead.cpp bamglib/MeshWrite.cpp		\
problem.cpp strversionnumber.cpp InitFunct.cpp Algo/lgalgo.cpp		\
femlib/Element_P2h.cpp load.cpp lgfem.cpp

FreeFem___SOURCES=$(BASESOURCES) Graphics/Xrgraph.cpp
FreeFem___g_SOURCES=$(BASESOURCES) Graphics/Xrgraph.cpp
FreeFem___nw_SOURCES=$(BASESOURCES) Graphics/sansrgraph.cpp
FreeFem___nw_g_SOURCES=$(BASESOURCES) Graphics/sansrgraph.cpp
FreeFem___x11_SOURCES=$(BASESOURCES) Graphics/Xrgraph.cpp
FreeFem___x11_g_SOURCES=$(BASESOURCES) Graphics/Xrgraph.cpp
FreeFem___glx_SOURCES=$(BASESOURCES) Graphics/xglrgraf.cpp
FreeFem___glx_g_SOURCES=$(BASESOURCES) Graphics/xglrgraf.cpp
FreeFem___mpi_SOURCES=$(BASESOURCES) Graphics/sansrgraph.cpp	\
mpi/parallelempi.cpp
FreeFem___mpi_g_SOURCES=$(BASESOURCES) Graphics/sansrgraph.cpp	\
mpi/parallelempi.cpp

BUILT_SOURCES=versionnumber.hpp
versionnumber.hpp:
	(echo "#define VersionFreeFempp " '$(VERSION)'; \
	echo "#define VersionFreeFemDate " '"'`date`'"') >$@

# Includes
# --------

INCLUDEX11 = /usr/X11R6/include
INCLUDES= -I$(INCLUDEX11) -Ifemlib -Ibamglib \
	-IGraphics -Impi -IAlgo $(IUMFPACK)

# ------------------------------
# compilation flags of FreeFem++ 
# ------------------------------
#  -DDRAWING  bamg active some drawing facility in bamg (inquire mesh)
#  -DLONG_LONG active the use of  long long in bamg
#  -DDEBUG active the  checking in bamg  
#  -DNCHECKPTR remove check pointer facility 
#  -DCHECK_KN active subscripting of some array  
#  -DWITHCHECK of the dynamic type ckecking of the langague (very slow)
#  -DEIGENVALUE  to compile the eigen value part
#  ------------------------------------------------
# FLAGS for optimisated version 
OFFFLAGS= -DDRAWING -DLONG_LONG  -DNCHECKPTR 
#  flags for debug version  
GFFFLAGS= -DDRAWING -DLONG_LONG -DCHECK_KN #-DWITHCHECK

OOPTFLAGS= -O  # compilation with optimisation
GOPTFLAGS= -g  # compilation with debugging 

FreeFem___CXXFLAGS=$(OOPTFLAGS) $(OFFFLAGS)
FreeFem___g_CXXFLAGS=$(GOPTFLAGS) $(GFFFLAGS)
FreeFem___nw_CXXFLAGS=$(OOPTFLAGS) $(OFFFLAGS)
FreeFem___nw_g_CXXFLAGS=$(GOPTFLAGS) $(GFFFLAGS)
FreeFem___x11_CXXFLAGS=$(OOPTFLAGS) $(OFFFLAGS)
FreeFem___x11_g_CXXFLAGS=$(GOPTFLAGS) $(GFFFLAGS)
FreeFem___glx_CXXFLAGS=$(OOPTFLAGS) $(OFFFLAGS)
FreeFem___glx_g_CXXFLAGS=$(GOPTFLAGS) $(GFFFLAGS)
FreeFem___mpi_CXXFLAGS=$(OOPTFLAGS) $(OFFFLAGS) $(MPIFLAGS) -DPARALLELE
FreeFem___mpi_g_CXXFLAGS=$(GOPTFLAGS) $(GFFFLAGS) $(MPIFLAGS) -DPARALLELE

# Libraries
# ---------

LIBX11= -L/usr/X11R6/lib -lX11 -lm 
LIBS=$(LIBX11) -ldl

FreeFem___glx_LDADD=-lGL -lGLU
FreeFem___glx_g_LDADD=-lGL -lGLU
