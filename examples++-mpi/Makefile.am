# $Id$

TESTS=regtests.sh
TESTS_ENVIRONMENT=PROGLOC=$(top_srcdir)/src EXEEXT=$(EXEEXT) MPIPROG=$(MPIPROG)
EXTRA_DIST=*.edp *.idp  regtests.sh  \
MPICG.cpp \
MUMPS_FreeFem.cpp \
MUMPS_all.cpp \
complex_SuperLU_DIST_FreeFem.cpp \
complex_pastix_FreeFem.cpp \
dSuperLU_DIST.cpp \
ffsuperludistoption-1.hpp \
ffsuperludistoption.hpp \
generaldefs.h \
hips_FreeFem.cpp \
hypre_FreeFem.cpp \
interfacepastix.cpp \
parms_FreeFem.cpp \
real_SuperLU_DIST_FreeFem.cpp \
real_pastix_FreeFem.cpp 




LIST_COMPILE=MUMPS_FreeFem.$(DYLIB_SUFFIX) hips_FreeFem.$(DYLIB_SUFFIX) \
MUMPS_all.$(DYLIB_SUFFIX) interfacepastix.$(DYLIB_SUFFIX) \
complex_SuperLU_DIST_FreeFem.$(DYLIB_SUFFIX) real_SuperLU_DIST_FreeFem.$(DYLIB_SUFFIX) \
complex_pastix_FreeFem.$(DYLIB_SUFFIX) real_pastix_FreeFem.$(DYLIB_SUFFIX) \
dSuperLU_DIST.$(DYLIB_SUFFIX) MPICG.$(DYLIB_SUFFIX) parms_FreeFem.$(DYLIB_SUFFIX)

all-local:$(LIST_COMPILE)
	echo "finish compile load mpi solver !"

.cpp.$(DYLIB_SUFFIX): ff-c++
	 ../examples++-load//ff-c++ -auto  $< 
