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
# SUMMARY : ...
# LICENSE : LGPLv3
# ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
# AUTHORS : ...
# E-MAIL  : ...

EXTRA_DIST = launchff++.cpp
LDADD2 = $(ARPACKLIBS) $(UMFPACKLIBS) $(BLASLIBS) $(G2CLIB)
# remark the libff.dll was splited in 6 because ld trap
# to much routine => memory fault in ld
# but now I see thant the couple allocation/delete
# must be make in same .dll file under windows
# is impossible to manage with the template generation
# so one .dll (this work with gcc version 3.4.4 (cygming special) (gdc 0.12, using dmd 0.125))
# change the entry points mymain to mainff
# to be compatibile with bamg.
#

FF_OBJ0 = CodeAlloc.o AFunction2.o CheckPtr.o lex.o global.o environment.o load.o
FF_OBJ_AC = array_complex.o
FF_OBJ_AR = array_real.o
FF_OBJ_AL = array_long.o
FF_OBJ1 = AFunction.o InitFunct.o strversionnumber.o mt19937ar.o string_def.o
FF_OBJ2 = BamgFreeFem.o Drawing.o Element_P2h.o Element_RT.o FESpace.o FQuadTree.o InitFunct.o Mesh2.o MeshDraw.o MeshGeom.o MeshQuad.o MeshRead.o MeshWrite.o Meshio.o Metric.o QuadTree.o QuadratureFormular.o R2.o SetOfE4.o @EIGENOBJ@ fem.o gibbs.o  lgalgo.o lgfem.o lgmat.o lgmesh.o  mshptg.o problem.o DefColor.o lgmesh3.o  P012_2d.o P012_3d.o P012_3dSurf.o P012_3dCurve.o FESpacen.o HashMatrix.o CG.o SparseLinearSolver.o chrono.o libmeshb7.o memory.o eigenv.o Mesh1dn.o Mesh2dn.o Mesh3dn.o MeshSn.o MeshSL.o GQuadTree.o glumesh2D.o UMFPack_Solver.o splitsimplex.o P1IsoValue.o
FF_OBJIDE = FreeFem___cs-cs.o FreeFem___cs-commands.o FreeFem___cs-draw.o FreeFem___cs-strversionnumber.o FreeFem___cs-socket.o FreeFem___cs-threads.o FreeFem___cs-hl_yacc.o FreeFem___cs-hl_lex.o FreeFem___cs-highlight.o FreeFem___cs-spawn.o FreeFem___cs-editor.o FreeFem___cs-server.o    windres.o

FFD_OBJ0 = ../fflib/CodeAlloc.o ../fflib/AFunction2.o ../femlib/CheckPtr.o ../fflib/lex.o ../fflib/global.o ../fflib/environment.o ../fflib/load.o
FFD_OBJ_AC = $(patsubst %.o,../fflib/%.o,$(FF_OBJ_AC))
FFD_OBJ_AL = $(patsubst %.o,../fflib/%.o,$(FF_OBJ_AL))
FFD_OBJ_AR = $(patsubst %.o,../fflib/%.o,$(FF_OBJ_AR))
FFD_OBJ1 = $(patsubst %.o,../fflib/%.o,$(FF_OBJ1))
if ENABLE_EIGEN
FFD_OBJ2 = ../femlib/BamgFreeFem.o ../femlib/Drawing.o ../femlib/Element_P2h.o ../femlib/Element_RT.o ../femlib/FESpace.o ../femlib/FQuadTree.o ../fflib/InitFunct.o ../bamglib/Mesh2.o ../bamglib/MeshDraw.o ../bamglib/MeshGeom.o ../bamglib/MeshQuad.o ../bamglib/MeshRead.o ../bamglib/MeshWrite.o ../bamglib/Meshio.o ../bamglib/Metric.o ../bamglib/QuadTree.o ../femlib/QuadratureFormular.o ../bamglib/R2.o ../bamglib/SetOfE4.o ../Eigen/eigenvalue.o ../femlib/fem.o ../femlib/gibbs.o ../Algo/lgalgo.o ../fflib/lgfem.o ../fflib/lgmat.o ../fflib/lgmesh.o ../femlib/mshptg.o ../fflib/problem.o ../Graphics/DefColor.o ../fflib/lgmesh3.o ../femlib/P012_2d.o ../femlib/P012_3d.o ../femlib/P012_3dSurf.o ../femlib/FESpacen.o ../femlib/HashMatrix.o ../femlib/CG.o ../femlib/SparseLinearSolver.o ../libMesh/chrono.o ../libMesh/libmeshb7.o ../libMesh/memory.o ../libMesh/eigenv.o ../femlib/Mesh1dn.o ../femlib/Mesh2dn.o ../femlib/Mesh3dn.o ../femlib/MeshSn.o ../femlib/GQuadTree.o ../fflib/glumesh2D.o ../fflib/UMFPack_Solver.o ../femlib/splitsimplex.o ../fflib/P1IsoValue.o
else
FFD_OBJ2 = ../femlib/BamgFreeFem.o ../femlib/Drawing.o ../femlib/Element_P2h.o ../femlib/Element_RT.o ../femlib/FESpace.o ../femlib/FQuadTree.o ../fflib/InitFunct.o ../bamglib/Mesh2.o ../bamglib/MeshDraw.o ../bamglib/MeshGeom.o ../bamglib/MeshQuad.o ../bamglib/MeshRead.o ../bamglib/MeshWrite.o ../bamglib/Meshio.o ../bamglib/Metric.o ../bamglib/QuadTree.o ../femlib/QuadratureFormular.o ../bamglib/R2.o ../bamglib/SetOfE4.o  ../femlib/fem.o ../femlib/gibbs.o ../Algo/lgalgo.o ../fflib/lgfem.o ../fflib/lgmat.o ../fflib/lgmesh.o ../femlib/mshptg.o ../fflib/problem.o ../Graphics/DefColor.o ../fflib/lgmesh3.o ../femlib/P012_2d.o ../femlib/P012_3d.o ../femlib/P012_3dSurf.o ../femlib/FESpacen.o ../femlib/HashMatrix.o ../femlib/CG.o ../femlib/SparseLinearSolver.o ../libMesh/chrono.o ../libMesh/libmeshb7.o ../libMesh/memory.o ../libMesh/eigenv.o ../femlib/Mesh1dn.o ../femlib/Mesh2dn.o ../femlib/Mesh3dn.o ../femlib/MeshSn.o ../femlib/GQuadTree.o ../fflib/glumesh2D.o ../fflib/UMFPack_Solver.o ../femlib/splitsimplex.o ../fflib/P1IsoValue.o

endif

FFD_OBJIDE = $(patsubst %.o,../ide/%.o,$(FF_OBJIDE))
FFD_OBJLMSH = $(patsubst %.o,../ide/%.o,$(FF_OBJLMSH))

LIBS_FF = ../fflib/ffapi.o libff.dll
LIBS_FFmpi = ../mpi/ffapi.o libff.dll

all-local: $(WIN32DLLTARGET)
	echo done

# FFCS does not need most Windows executables since it compiles
# executables in ffcs/src and ffcs/pack. But still keep one executable
# to make sure that the whole process works and that libff.dll is
# correctly built.

if ENABLE_FFCS
win32-dll-target: FreeFem++.exe ff-c++
	-for i in $(LDADD2) ; do case $$i in *.dll) echo cp -p $$i .; cp -p $$i . ;; esac ; done
	./copy-dll.sh
	echo "on a fini"
else
win32-dll-target: FreeFem++.exe FreeFem++-nw.exe bamg.exe cvmsh2.exe launchff++.exe ff-c++ $(MPIPROG)
	-for i in $(LDADD2) ; do case $$i in *.dll) echo cp -p $$i .; cp -p $$i . ;; esac ; done
	./copy-dll.sh
	echo "on a fini"
endif

libff.dll: $(FFD_OBJ0) $(FFD_OBJ_AR) $(FFD_OBJ_AC) $(FFD_OBJ_AL) $(FFD_OBJ1) $(FFD_OBJ2) ../lglib/lg.tab.o
	$(CXX) $(GCCNOCYGWIN) -shared -Wl,--enable-auto-import $^ -o $@ $(LDADD2)

parallelempi-empty.o:../mpi/parallelempi-empty.cpp
	$(CXX) $(CXXFLAGS) $< -c -o $@

FreeFem++-nw.exe: ../Graphics/sansrgraph.o parallelempi-empty.o ../lglib/mymain.o $(LIBS_FF)
	$(CXX) $(GCCNOCYGWIN) libff.dll $^ -o $@ -lcomdlg32

# ALH - 28/7/14 - the new automake default is to place the objects in
# the same directory as the corresponding source file

FreeFem++.exe: ../Graphics/sansrgraph.o parallelempi-empty.o ../lglib/mymain.o $(LIBS_FF)
	$(CXX) $(GCCNOCYGWIN) libff.dll $^ -o $@ -lcomdlg32

## FH.. Strange : pb with MPI finaly with ffapi-mpi.dll ( J'ai compris)
## dur dur ...
## on mais init ffpt dans mymain.cpp
FreeFem++-mpi.exe:	../Graphics/sansrgraph.o ../mpi/parallelempi.o ../lglib/mymain.o $(LIBS_FFmpi)
	$(MPICXX) $(GCCNOCYGWIN) libff.dll $^ -o $@ -lcomdlg32 $(MPI_LIB)
bamg.exe: libff.dll $(LIBS_FF)
	$(CXX) $(GCCNOCYGWIN) ../bamg/bamg.o -Wl,--enable-auto-import $^ -o $@
cvmsh2.exe: $(LIBS_FF)
	$(CXX) $(GCCNOCYGWIN) -Wl,--enable-auto-import ../bamg/cvmsh2.o $^ -o $@
drawbdmesh.exe:../std/Pcrgraph.o ../bamg/drawbdmesh.o $(LIBS_FF)
	$(CXX) $(GCCNOCYGWIN) -mwindows -Wl,--enable-auto-import $^ -o $@
launchff++.exe:launchff++.cpp
	$(CXX) $(GCCNOCYGWIN) launchff++.cpp -o launchff++.exe -lcomdlg32

ff-c++:../../plugin/seq/load.link.in
	../../config.status --file=ff-c++:$<
	chmod a+x ff-c++
install-exec-local:: install-exec-local--$(WIN32DLLTARGET)

install-exec-local--win32-dll-target::
	$(INSTALL_SCRIPT) launchff++.exe FreeFem++.exe FreeFem++-nw.exe bamg.exe cvmsh2.exe ff-c++ $(MPIPROG) $(DESTDIR)${bindir}
	$(INSTALL_SCRIPT) *.dll $(DESTDIR)${bindir}
install-exec-local--::

# FFCS: automake does not clean the exes automatically because they are not built automatically.
clean-local::
	-rm *.exe
	-rm *.dll
