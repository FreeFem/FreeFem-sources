
LDADD2 = $(ARPACKLIB) $(UMFPACKLIB)  $(BLASLIB) $(G2CLIB)
# remark the libff.dll  is split in 6 because  ld trap
#  to much routine => memory fault 
FF_OBJ0 = CodeAlloc.o AFunction2.o  CheckPtr.o lex.o  global.o 
FF_OBJ_AC = array_complex.o
FF_OBJ_AR = array_real.o
FF_OBJ_AL = array_long.o
FF_OBJ1 = AFunction.o   InitFunct.o strversionnumber.o  mt19937ar.o
FF_OBJ2 = BamgFreeFem.o  Drawing.o Element_P2h.o Element_RT.o FESpace.o FQuadTree.o InitFunct.o Mesh2.o MeshDraw.o MeshGeom.o MeshQuad.o MeshRead.o MeshWrite.o Meshio.o Metric.o QuadTree.o QuadratureFormular.o R2.o SetOfE4.o eigenvalue.o fem.o gibbs.o  lgalgo.o lgfem.o lgmat.o lgmesh.o load.o mshptg.o problem.o DefColor.o environment.o
FF_OBJIDE = FreeFem___cs-cs.o FreeFem___cs-commands.o FreeFem___cs-draw.o FreeFem___cs-strversionnumber.o FreeFem___cs-socket.o FreeFem___cs-threads.o FreeFem___cs-hl_yacc.o FreeFem___cs-hl_lex.o FreeFem___cs-highlight.o FreeFem___cs-spawn.o FreeFem___cs-editor.o FreeFem___cs-server.o    windres.o
FFD_OBJ0 = $(patsubst %.o,../fflib/%.o,$(FF_OBJ0))
FFD_OBJ_AC = $(patsubst %.o,../fflib/%.o,$(FF_OBJ_AC))
FFD_OBJ_AL = $(patsubst %.o,../fflib/%.o,$(FF_OBJ_AL))
FFD_OBJ_AR = $(patsubst %.o,../fflib/%.o,$(FF_OBJ_AR))
FFD_OBJ1 = $(patsubst %.o,../fflib/%.o,$(FF_OBJ1))
FFD_OBJ2 = $(patsubst %.o,../fflib/%.o,$(FF_OBJ2))
FFD_OBJIDE = $(patsubst %.o,../ide/%.o,$(FF_OBJIDE))





all-local: $(WIN32DLLTARGET)
	-for i in $(LDADD2) ; do case $$i in  *.dll) echo cp -p $$i .; cp -p $$i . ;;	  esac ; done
	echo done
win32-dll-target: FreeFem++.exe FreeFem++-nw.exe  FreeFem++-cs.exe
	echo "on a fini"

libff0.dll: $(FFD_OBJ0)
	g++ -mno-cygwin -shared $^  -o $@
libffAR.dll: $(FFD_OBJ_AR) libff0.dll
	g++ -mno-cygwin -shared $^  -o $@
libffAC.dll: $(FFD_OBJ_AC) libff0.dll
	g++ -mno-cygwin -shared $^  -o $@
libffAL.dll: $(FFD_OBJ_AL) libff0.dll
	g++ -mno-cygwin -shared $^  -o $@
libff1.dll: $(FFD_OBJ1) libff0.dll libffAR.dll libffAC.dll libffAL.dll
	g++ -mno-cygwin -shared $^  -o $@
libff.dll: $(FFD_OBJ2)  ../lglib/liblg.a libff0.dll libffAR.dll libffAC.dll libffAL.dll  libff1.dll
	g++  -shared  -mno-cygwin -mwindows -mthreads  $^  $(LDADD2) -o $@ 
FreeFem++.exe: ../std/Pcrgraph.o libff0.dll libffAR.dll libffAC.dll libffAL.dll  libff1.dll  libff.dll 
	g++   -mno-cygwin -mwindows -mthreads  $^  -o $@
FreeFem++-nw.exe: ../nw/sansrgraph.o  libff0.dll libffAR.dll libffAC.dll libffAL.dll  libff1.dll libff.dll
	g++   -mno-cygwin libff.dll $^ -o $@ 
FreeFem++-cs.exe:$(FFD_OBJIDE) libff0.dll libffAR.dll libffAC.dll libffAL.dll  libff1.dll libff.dll
	g++  $^  `$(FLTK_CONFIG) --ldflags` -mno-cygwin    -o $@
