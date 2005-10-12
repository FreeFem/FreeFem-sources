LDADD2=  $(ARPACKLIB) $(UMFPACKLIB)  $(BLASLIB) $(G2CLIB)
OBJ0=CodeAlloc.o AFunction2.o AFunction3.o CheckPtr.o lex.o  global.o 
OBJ1=AFunction.o  libff0.dll  InitFunct.o strversionnumber.o  mt19937ar.o
OBJ2=BamgFreeFem.o  Drawing.o Element_P2h.o Element_RT.o FESpace.o FQuadTree.o InitFunct.o Mesh2.o MeshDraw.o MeshGeom.o MeshQuad.o MeshRead.o MeshWrite.o Meshio.o Metric.o QuadTree.o QuadratureFormular.o R2.o SetOfE4.o eigenvalue.o fem.o gibbs.o  lgalgo.o lgfem.o lgmat.o lgmesh.o load.o mshptg.o problem.o DefColor.o

OBJIDE=FreeFem___cs-cs.o FreeFem___cs-commands.o FreeFem___cs-draw.o FreeFem___cs-strversionnumber.o FreeFem___cs-socket.o FreeFem___cs-threads.o FreeFem___cs-hl_yacc.o FreeFem___cs-hl_lex.o FreeFem___cs-highlight.o FreeFem___cs-spawn.o FreeFem___cs-editor.o FreeFem___cs-server.o    windres.o

OBJ00=$(patsubst %.o,../fflib/%.o,$(OBJ0))
OBJ11=$(patsubst %.o,../fflib/%.o,$(OBJ1))
OBJ22=$(patsubst %.o,../fflib/%.o,$(OBJ2))
OBJIDE1=$(patsubst %.o,../ide/%.o,$(OBJIDE))


all-local: $(WIN32DLLTARGET)
	echo done
win32-dll-target:libff0.dll libff1.dll libff2.dll FreeFem++.exe FreeFem++-nw.exe  FreeFem++-cs.exe
	echo "on a fini"
libff0.dll: $(OBJ00)
	g++ -mno-cygwin -shared $(OBJ00)  -o libff0.dll	
libff1.dll: $(OBJ11) libff0.dll
	g++ -mno-cygwin -shared  $(OBJ11) libff0.dll  -o libff1.dll
libff2.dll: $(OBJ22)  ../lglib/liblg.a libff0.dll libff1.dll
	g++  -shared  -mno-cygwin -mwindows -mthreads  $(OBJ22)  ../lglib/liblg.a libff0.dll libff1.dll $(LDADD2) -o libff2.dll 
FreeFem++.exe: ../std/Pcrgraph.o  libff0.dll libff1.dll libff2.dll 
	g++  ../std/Pcrgraph.o -mno-cygwin -mwindows -mthreads  libff0.dll libff1.dll libff2.dll  -o FreeFem++.exe 
FreeFem++-nw.exe: ../nw/sansrgraph.o  libff0.dll libff1.dll libff2.dll
	g++  ../nw/sansrgraph.o -mno-cygwin  libff0.dll libff1.dll libff2.dll -o FreeFem++-nw.exe   
FreeFem++-cs.exe:$(OBJIDE1) libff0.dll libff1.dll libff2.dll
	g++  $(OBJIDE1) `$(FLTK_CONFIG) --ldflags` -mno-cygwin  libff0.dll libff1.dll libff2.dll -o FreeFem++-cs.exe 


