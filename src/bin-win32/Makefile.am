LDADD2=  $(ARPACKLIB) $(UMFPACKLIB)  $(BLASLIB) $(G2CLIB)
# remark the libff.dll  is split in 6 because  ld trap
#  to much routine => memory fault 
FF_OBJ0=CodeAlloc.o AFunction2.o  CheckPtr.o lex.o  global.o 
FF_OBJ_AC=array_complex.cpp
FF_OBJ_AR=array_real.cpp
FF_OBJ_AL=array_long.cpp
FF_OBJ1=AFunction.o  libff0.dll  InitFunct.o strversionnumber.o  mt19937ar.o
FF_OBJ2=BamgFreeFem.o  Drawing.o Element_P2h.o Element_RT.o FESpace.o FQuadTree.o InitFunct.o Mesh2.o MeshDraw.o MeshGeom.o MeshQuad.o MeshRead.o MeshWrite.o Meshio.o Metric.o QuadTree.o QuadratureFormular.o R2.o SetOfE4.o eigenvalue.o fem.o gibbs.o  lgalgo.o lgfem.o lgmat.o lgmesh.o load.o mshptg.o problem.o DefColor.o

FF_OBJIDE=FreeFem___cs-cs.o FreeFem___cs-commands.o FreeFem___cs-draw.o FreeFem___cs-strversionnumber.o FreeFem___cs-socket.o FreeFem___cs-threads.o FreeFem___cs-hl_yacc.o FreeFem___cs-hl_lex.o FreeFem___cs-highlight.o FreeFem___cs-spawn.o FreeFem___cs-editor.o FreeFem___cs-server.o    windres.o

FFO_OBJ0=$(patsubst %.o,../fflib/%.o,$(FF_OBJ0))
FFO_OBJ_AC=$(patsubst %.o,../fflib/%.o,$(FF_OBJ_AC))
FFO_OBJ_AL=$(patsubst %.o,../fflib/%.o,$(FF_OBJ_AL))
FFO_OBJ_AR=$(patsubst %.o,../fflib/%.o,$(FF_OBJ_AR))
FFO_OBJ1=$(patsubst %.o,../fflib/%.o,$(FF_OBJ1))
FFO_OBJ2=$(patsubst %.o,../fflib/%.o,$(FF_OBJ2))
FFO_OBJIDES=$(patsubst %.o,../ide/%.o,$(FF_OBJIDE))
# the prefixe of lib freefem++
FF_LIB=libff 

all-local: $(WIN32DLLTARGET)
	-for i in $(LDADD2) ; do case $$i in  *.dll) echo cp -p $$i .; cp -p $$i . ;;	  esac ; done
	echo done
win32-dll-target:$(FF_LIB)0.dll $(FF_LIB)1.dll $(FF_LIB)2.dll FreeFem++.exe FreeFem++-nw.exe  FreeFem++-cs.exe
	echo "on a fini"
$(FF_LIB)0.dll: $(FFO_OBJ0)
	g++ -mno-cygwin -shared $^  -o $@
$(FF_LIB)AR.dll: $(FFO_OBJ_AR) $(FF_LIB)0.dll
	g++ -mno-cygwin -shared $^  -o $
$(FF_LIB)AC.dll: $(FFO_OBJ_AC) $(FF_LIB)0.dll
	g++ -mno-cygwin -shared $^  -o $@
$(FF_LIB)AL.dll: $(FFO_OBJ_AL) $(FF_LIB)0.dll
	g++ -mno-cygwin -shared $^  -o $@
$(FF_LIB)1.dll: $(FFO_OBJ1) $(FF_LIB)0.dll $(FF_LIB)AR.dll $(FF_LIB)AC.dll $(FF_LIB)AL.dll
	g++ -mno-cygwin -shared $^  -o $@
$(FF_LIB).dll: $(FFO_OBJ2)  ../lglib/liblg.a  $(FF_LIB)1.dll
	g++  -shared  -mno-cygwin -mwindows -mthreads  $(FFO_OBJ2)  $^  $(LDADD2) -o $@ 
FreeFem++.exe: ../std/Pcrgraph.o  $(FF_LIB).dll 
	g++   -mno-cygwin -mwindows -mthreads  $^  -o $@
FreeFem++-nw.exe: ../nw/sansrgraph.o   $(FF_LIB).dll
	g++   -mno-cygwin $(FF_LIB).dll $^ -o $@ 
FreeFem++-cs.exe:$(FFO_OBJIDE) $(FFO_DLL)
	g++  $^  `$(FLTK_CONFIG) --ldflags` -mno-cygwin    -o $@


