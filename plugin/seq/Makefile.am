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

all-local: @LOAD_COMPILE@

LIST_COMPILE=myfunction.$(DYLIB_SUFFIX) BernardiRaugel.$(DYLIB_SUFFIX) \
	Morley.$(DYLIB_SUFFIX) funcTemplate.$(DYLIB_SUFFIX) addNewType.$(DYLIB_SUFFIX) \
	qf11to25.$(DYLIB_SUFFIX) Element_P3.$(DYLIB_SUFFIX) Element_P4.$(DYLIB_SUFFIX) \
	Element_P3dc.$(DYLIB_SUFFIX) Element_P4dc.$(DYLIB_SUFFIX) \
	Element_PkEdge.$(DYLIB_SUFFIX) msh3.$(DYLIB_SUFFIX) \
	splitmesh3.$(DYLIB_SUFFIX) splitmesh6.$(DYLIB_SUFFIX) \
	splitmesh4.$(DYLIB_SUFFIX) splitmesh12.$(DYLIB_SUFFIX) \
	ffrandom.$(DYLIB_SUFFIX) \
	medit.$(DYLIB_SUFFIX) \
	mat_dervieux.$(DYLIB_SUFFIX) lgbmo.$(DYLIB_SUFFIX) mat_psi.$(DYLIB_SUFFIX) \
	ppm2rnm.$(DYLIB_SUFFIX) DxWriter.$(DYLIB_SUFFIX) \
	pcm2rnm.$(DYLIB_SUFFIX) $(DYLIB_OTHER_COMPILE) \
	iovtk.$(DYLIB_SUFFIX) gmsh.$(DYLIB_SUFFIX) MetricKuate.$(DYLIB_SUFFIX) \
	Element_P1dc1.$(DYLIB_SUFFIX) BinaryIO.$(DYLIB_SUFFIX) \
	isoline.$(DYLIB_SUFFIX) thresholdings.$(DYLIB_SUFFIX) \
	VTK_writer_3d.$(DYLIB_SUFFIX) VTK_writer.$(DYLIB_SUFFIX) \
	splitedges.$(DYLIB_SUFFIX) Element_Mixte.$(DYLIB_SUFFIX) Element_Mixte3d.$(DYLIB_SUFFIX) \
	myfunction2.$(DYLIB_SUFFIX) \
	MetricPk.$(DYLIB_SUFFIX) FreeFemQA.$(DYLIB_SUFFIX) shell.$(DYLIB_SUFFIX) \
	@TOOL_DYLIB_pipe@ symmetrizeCSR.$(DYLIB_SUFFIX) BEC.$(DYLIB_SUFFIX) \
	bfstream.$(DYLIB_SUFFIX) ff-AiryBiry.$(DYLIB_SUFFIX) \
	exactpartition.$(DYLIB_SUFFIX) ClosePoints.$(DYLIB_SUFFIX) \
	Element_HCT.$(DYLIB_SUFFIX) Curvature.$(DYLIB_SUFFIX) Element_P1bl.$(DYLIB_SUFFIX) \
	Element_QF.$(DYLIB_SUFFIX) distance.$(DYLIB_SUFFIX) Element_P1ncdc.$(DYLIB_SUFFIX) \
	Element_P2bulle3.$(DYLIB_SUFFIX) \
	Element_P2pnc.$(DYLIB_SUFFIX) SaveHB.$(DYLIB_SUFFIX) IncompleteCholesky.$(DYLIB_SUFFIX) \
	Schur-Complement.$(DYLIB_SUFFIX) biofunc.$(DYLIB_SUFFIX) mat_edgeP1.$(DYLIB_SUFFIX) \
	geophysics.$(DYLIB_SUFFIX) CircumCenter.$(DYLIB_SUFFIX)

# FFCS - some libraries are skipped because the corresponding tool is deactivated.

LIST_COMPILE_PKG=@TOOL_DYLIB_tetgen@ @TOOL_DYLIB_superlu@ dfft.$(DYLIB_SUFFIX) \
	@TOOL_DYLIB_umfpack@ NewSolver.$(DYLIB_SUFFIX) @TOOL_DYLIB_lapack@ \
	ffnewuoa.$(DYLIB_SUFFIX) @TOOL_DYLIB_gmm@ @TOOL_DYLIB_yams@ @TOOL_DYLIB_mmg3d@ @TOOL_DYLIB_mmg@ @TOOL_DYLIB_mshmet@ \
	@TOOL_DYLIB_gsl@ @TOOL_DYLIB_mumps_seq@ @TOOL_DYLIB_ipopt@ @TOOL_DYLIB_nlopt@ ff-cmaes.$(DYLIB_SUFFIX) \
	@TOOL_DYLIB_scotch@ @TOOL_DYLIB_metis@ @TOOL_DYLIB_mkl@ @TOOL_DYLIB_MMAP@ @TOOL_DYLIB_iohdf5@
	

EXTRA_PROGRAMS= ffmaster
bin_PROGRAMS= @BIN_ffmaster@

# FFCS parallel make: some targets need to be built first, sequentially

load_compile: ff-c++ WHERE_LIBRARY-download include.done 
#
# max_load
	-@rm "Missing-plugins-$(DYLIB_SUFFIX).log" 2>/dev/null || true
	@$(MAKE) -s $(AM_MAKEFLAGS) $(LIST_COMPILE) $(LIST_COMPILE_PKG)
	@if test -f Missing-plugins-$(DYLIB_SUFFIX).log ; then cat Missing-plugins-$(DYLIB_SUFFIX).log; fi; exit 0
	@echo " Warning missing plugin:" `for i in $(LIST_COMPILE) $(LIST_COMPILE_PKG); do if test ! -s $i ; then j=1; echo "$i," ;fi; done`
	@echo " finish build list $(DYLIB_SUFFIX)"

.cpp.$(DYLIB_SUFFIX): ff-c++ Makefile
	 ./ff-c++ -auto $<

#to force the recompile if change ... FH
$(LIST_COMPILE): ../seq/ff-c++ ../seq/WHERE_LIBRARY-download

# FFCS - 26/10/11 - Unpacking include.tar.gz is very often buggy under Cygwin (softlinks are randomly replaced with
# empty files without any access right). So just replace the whole thing with a plain copy.

allheaders=../../src/fflib/AddNewFE.h ../../src/fflib/AFunction_ext.hpp ../../src/fflib/AFunction.hpp ../../src/fflib/AnyType.hpp	\
	../../src/fflib/array_init.hpp ../../src/fflib/array_resize.hpp ../../src/fflib/array_tlp.hpp				\
	../../src/femlib/assertion.hpp ../../src/femlib/BamgFreeFem.hpp ../../src/Algo/BFGS.hpp ../../src/Algo/BrentLS.hpp		\
	../../src/femlib/CGNL.hpp ../../src/fflib/CodeAlloc.hpp ../../config.h 	\
	../../src/Algo/CubicLS.hpp ../../src/Algo/defs.hpp ../../src/femlib/DOperator.hpp ../../src/libMesh/eigenv.h		\
	../../src/fflib/endian.hpp ../../src/fflib/environment.hpp ../../src/fflib/error.hpp ../../src/femlib/fem3.hpp		\
	../../src/femlib/fem.hpp ../../src/femlib/FESpace.hpp ../../src/femlib/FESpacen.hpp ../../src/fflib/ff++.hpp		\
	../../src/fflib/ffstack.hpp ../../src/femlib/FQuadTree.hpp ../../src/femlib/GenericMesh.hpp				\
	../../src/Graphics/getprog-unix.hpp ../../src/Graphics/glrgraph.hpp ../../src/femlib/gmres.hpp				\
	../../src/femlib/GQuadTree.hpp ../../src/femlib/HashTable.hpp ../../src/femlib/HeapSort.hpp ../../src/fflib/InitFunct.hpp	\
	../../src/fflib/ffapi.hpp ../../src/femlib/Label.hpp ../../src/fflib/lex.hpp ../../src/fflib/lgfem.hpp			\
	../../src/fflib/lgmesh3.hpp ../../src/fflib/lgsolver.hpp ../../src/lglib/lg.tab.hpp ../../src/femlib/libmesh5.h		\
	../../src/Algo/LineSearch.hpp ../../src/femlib/MatriceCreuse.hpp ../../src/femlib/MatriceCreuse_tpl.hpp			\
	../../src/femlib/Mesh1dn.hpp ../../src/femlib/Mesh2dn.hpp ../../src/bamglib/Mesh2.h ../../src/femlib/Mesh3dn.hpp		\
	../../src/femlib/MeshSn.hpp ../../src/femlib/MeshLn.hpp   												\
	../../src/bamglib/Meshio.h ../../src/femlib/MeshPoint.hpp ../../src/bamglib/meshtype.h ../../src/bamglib/Metric.h		\
	../../src/Graphics/mode_open.hpp ../../src/Algo/NewtonRaphson.hpp ../../src/Algo/NRJ.hpp ../../src/fflib/Operator.hpp	\
	../../src/Algo/Optima.hpp ../../src/Algo/Param.hpp ../../src/femlib/PkLagrange.hpp ../../src/fflib/PlotStream.hpp		\
	../../src/fflib/problem.hpp ../../src/femlib/QuadratureFormular.hpp ../../src/bamglib/QuadTree.h ../../src/femlib/R1.hpp	\
	../../src/bamglib/R2.h ../../src/femlib/R2.hpp ../../src/femlib/R3.hpp ../../src/femlib/RefCounter.hpp			\
	../../src/Graphics/rgraph.hpp ../../src/femlib/RNM.hpp ../../src/femlib/RNM_opc.hpp ../../src/femlib/RNM_op.hpp		\
	../../src/femlib/RNM_tpl.hpp ../../src/Algo/RosenBrock.hpp ../../src/fflib/Serialize.hpp ../../src/bamglib/SetOfE4.h	\
	../../src/fflib/showverb.hpp ../../src/femlib/splitsimplex.hpp ../../src/fflib/String.hpp				\
	../../src/fflib/strversionnumber.hpp ../../src/fflib/throwassert.hpp ../../src/femlib/ufunction.hpp			\
	../../src/fflib/versionnumber.hpp ../../src/bamglib/write_hdf5.hpp	../../src/bamglib/write_xdmf.hpp                   \
	../../src/fflib/P1IsoValue.hpp \
	../../src/femlib/SkyLineSolver.hpp \
	../../src/femlib/SolverSkyLine.hpp \
	../../src/femlib/SparseLinearSolver.hpp \
	../../src/femlib/VirtualMatrix.hpp \
	../../src/femlib/VirtualSolver.hpp \
	../../src/femlib/CG.hpp \
	../../src/femlib/VirtualSolverCG.hpp \
	../../src/femlib/VirtualSolverSkyLine.hpp \
	../../src/femlib/VirtualSolverSparseSuite.hpp \
	../../src/femlib/HashMatrix.hpp \
	../../src/femlib/MatriceElementaire.hpp \
	../../src/femlib/DataFindBoundary.hpp

if ENABLE_FFCS

include.done: $(allheaders)
	-rm -rf include
	mkdir -p include
	cp $^ include
	touch $@

else
include.done: $(allheaders)
	-rm -rf include
	mkdir -p  include
	for i in $(allheaders); do  ln -s ../$$i include/. ; done
	touch $@
endif

clean-local::
	-rm -r include
	-rm include.done


# FFCS - 27/2/13 - remove dependency on Makefile to avoid recompiling everything everytime something is changed in the
# configuration step

$(LIST_COMPILE):ff-c++ include.done

clean-local::
	-rm *.o *.$(DYLIB_SUFFIX) load.link WHERE_LIBRARY-download ff-get-dep ff-c++ ff-pkg-download \
	  $(LIST_COMPILE) $(LIST_COMPILE_PKG)
	-rm -rf include include.done
ff-c++:load.link.in load.link WHERE_LIBRARY-download ff-get-dep ../../config.status
	../../config.status  --file=ff-c++:load.link.in
	chmod a+x ff-c++
load.link:load.link.in ../../config.status
	../../config.status  --file=load.link:load.link.in
	chmod a+x load.link
ff-pkg-download:ff-pkg-download.in ../../config.status
	../../config.status  --file=$@:$@.in
	chmod a+x $@
	cp $@ ../../3rdparty/bin
ff-get-dep:ff-get-dep.in ../../config.status
	../../config.status  --file=$@:$@.in
	chmod a+x $@
	-if [ -d ../../3rdparty/bin ] ;then cp $@ ../../3rdparty/bin; fi
WHERE_LIBRARY-download:ff-pkg-download FORCE  WHERE_LIBRARY-config WHERE_LIBRARY
	@./ff-pkg-download >$@-new
	@diff $@-new $@  || (mv $@-new $@ ;touch WHERE_LIBRARY)
WHERE_LIBRARY:
	touch $@
FORCE: ;

install-exec-local:: load_compile WHERE_LIBRARY
	$(mkinstalldirs) -m 755 $(DESTDIR)$(ff_prefix_dir)/lib
	test -n "$(MPIPROG)" && $(mkinstalldirs) -m 755 $(DESTDIR)$(ff_prefix_dir)/lib/mpi || true
	$(mkinstalldirs) -m 755 $(DESTDIR)$(ff_prefix_dir)/include
	$(mkinstalldirs) -m 755 $(DESTDIR)$(ff_prefix_dir)/etc
	$(INSTALL_DATA)  renumb.hpp  clapack.h ppmimg.h  bmo.hpp msh3.hpp pcm.hpp  include/*   $(DESTDIR)$(ff_prefix_dir)/include
	$(INSTALL)  -m 555 $(LIST_COMPILE) $(DESTDIR)$(ff_prefix_dir)/lib
	$(INSTALL)  -m 555 ff-get-dep.awk  WHERE_LIBRARY-config WHERE_LIBRARY  $(DESTDIR)$(ff_prefix_dir)/lib
	for i in  WHERE_LIBRARY-config WHERE_LIBRARY WHERE_LIBRARY-download; do \
	  sed <$$i >$(DESTDIR)$(ff_prefix_dir)/lib/$$i 's#$(abs_top_builddir)/3rdparty#$(ff_prefix_dir)#' ;\
	 done
	echo loadpath += \"./\"  >$(DESTDIR)$(ff_prefix_dir)/etc/freefem++.pref
	echo loadpath += \"$(ff_prefix_dir)/lib\"  >>$(DESTDIR)$(ff_prefix_dir)/etc/freefem++.pref
	echo includepath += \"$(ff_prefix_dir)/idp\"  >>$(DESTDIR)$(ff_prefix_dir)/etc/freefem++.pref
	$(INSTALL_SCRIPT) ff-c++  $(DESTDIR)${bindir}
	$(INSTALL_SCRIPT) ff-pkg-download  $(DESTDIR)${bindir}
	$(INSTALL_SCRIPT) ff-get-dep   $(DESTDIR)${bindir}
##	$(DESTDIR)${bindir}/ff-pkg-download $(DESTDIR) >$(DESTDIR)$(ff_prefix_dir)/lib/WHERE_LIBRARY-download
	-for i in $(LIST_COMPILE_PKG); do \
	 if [ -f $$i ] ; then 	$(INSTALL)  -m 555 $$i $(DESTDIR)$(ff_prefix_dir)/lib; fi; done
	- if [ -n "$(MPIPROG)" ] ; then \
		for i in $(LIST_COMPILE_PKG) $(LIST_COMPILE); do \
		 if [ -f ../mpi/$$i ] ; then \
	    	    $(INSTALL)  -m 555 ../mpi/$$i $(DESTDIR)$(ff_prefix_dir)/lib/mpi; \
	        else echo missing install mpi version ???  ../mpi/$$i ;\
	         fi;  \
	        done ; \
	fi

ffmaster$(EXEEXT): libff-mmap-semaphore.o ffmaster.c libff-mmap-semaphore.h
	$(CC) $(CFLAGS) ffmaster.c libff-mmap-semaphore.o -o ffmaster $(LIBSPTHREAD)

clean-local::
	-rm *.fg *.eps *~ *.ps  *.mesh *.mesh.gmsh ListOfAllocPtr-8.bin ffglut*.ppm Th.o.meshb	Th3.d.meshb	dep.sol.meshbgh \
ipopt.out	mmg.out xxxx ThFF.txt		TriQA.txt		cavidadTT.txt		mass.txt		xy.txt \
ThFF_Metric.txt		TriQA_Metric.txt	datanc.txt		u100m8.txt \
TTh.sol		Th.o.sol	Th3.sol		dep.o.sol	dep.sol		gradient.sol	hessien.sol pippo.data pippo.dx testsavemedit.solb med.gp bidule.vtk	disque.vtk g-iso
	-rm $(LIST_COMPILE_PKG)
	-cd ../mpi; rm $(LIST_COMPILE) $(LIST_COMPILE_PKG)

# FFCS: add tags for files that are not specified to automake by a SOURCE directive (this tags file is automatically
# taken into account by automake when it sees it).

TAGS:
	etags *.?pp
