DOTUVECLIB = wrapper_dotblas.o
F77 = g95
ARCH = ar
RANLIB = ranlib
ARCHFLAGS = cr
LDOTUVECLIB = libwrapperdotblas.a
all: dotu

#APPBLAS = pilaenv.o

dotu: $(DOTUVECLIB)
	$(ARCH) $(ARCHFLAGS) $(LDOTUVECLIB) $(DOTUVECLIB)
	$(RANLIB) $(LDOTUVECLIB)

clean :
	rm -f *.o

.f.o : ; $(F77) -fno-underscoring -O3 -c $*.f