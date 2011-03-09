# Downloading and compiling extra libraries
# -----------------------------------------
# $Id: Makefile,v 1.2 2010/03/31 08:50:56 hecht Exp $
include makefile.inc

all-local: hips

# Downloading and compiling HIPS
# ------------------------------
# http://hips.gforge.inria.fr/release/hips-1.2b-rc4.tar.gz
# Hips information

DIRPKG=../pkg
SRCDIR=hips-$(hips_VERSION)
PACKAGE=$(DIRPKG)/hips-$(hips_VERSION).tar.gz
SERVER=http://hips.gforge.inria.fr/release/
INSTALL=../..
hips_VERSION=1.2b-rc4

# rappel ::  SRC= $(wildcard *.c) ==> permet de chopper tout les points *.c d'un repertoire
# OBJ =  $(SRC:.c=.o)
# ==> implique *.o dans les repertoire pour recuperer

SRCPhidalCommon = $(SRCDIR)/SRC/PHIDAL/COMMON
ObjPhidalCommon = $(SRCPhidalCommon)/*.o

SRCPhidalOrdering = $(SRCDIR)/SRC/PHIDAL/ORDERING
ObjPhidalOrdering = $(SRCPhidalOrdering)/*.o 

SRCPhidalSequential = $(SRCDIR)/SRC/PHIDAL/SEQUENTIAL
ObjPhidalSequential  = $(SRCPhidalSequential)/*.o 

SRCDBMATRIX = $(SRCDIR)/SRC/BLOCK/DBMATRIX/
ObjBlockDBMATRIX   = $(SRCDIR)/SRC/BLOCK/DBMATRIX/amalgamate.o \
	$(SRCDIR)/SRC/BLOCK/DBMATRIX/DBMatrix2PhidalMatrix.o \
	$(SRCDIR)/SRC/BLOCK/DBMATRIX/DBMatrix_Build.o \
	$(SRCDIR)/SRC/BLOCK/DBMATRIX/DBMatrix_Copy.o \
	$(SRCDIR)/SRC/BLOCK/DBMATRIX/DBMatrix_FACTO2.o \
	$(SRCDIR)/SRC/BLOCK/DBMATRIX/DBMatrix_FACTO.o \
	$(SRCDIR)/SRC/BLOCK/DBMATRIX/DBMatrix_FACTOu.o \
	$(SRCDIR)/SRC/BLOCK/DBMATRIX/DBMatrix_func.o \
	$(SRCDIR)/SRC/BLOCK/DBMATRIX/DBMatrix_GEMM.o \
	$(SRCDIR)/SRC/BLOCK/DBMATRIX/DBMatrix_GEMM_N.o \
	$(SRCDIR)/SRC/BLOCK/DBMATRIX/DBMatrix_GEMMu.o \
	$(SRCDIR)/SRC/BLOCK/DBMATRIX/DBMATRIX_MatVec2.o \
	$(SRCDIR)/SRC/BLOCK/DBMATRIX/DBMATRIX_MatVec.o \
	$(SRCDIR)/SRC/BLOCK/DBMATRIX/DBMATRIX_MLILUPrec.o \
	$(SRCDIR)/SRC/BLOCK/DBMATRIX/DBMATRIX_MLIutil.o \
	$(SRCDIR)/SRC/BLOCK/DBMATRIX/DBMATRIX_PrecSolve.o \
	$(SRCDIR)/SRC/BLOCK/DBMATRIX/DBMATRIX_SEQ.o \
	$(SRCDIR)/SRC/BLOCK/DBMATRIX/DBMatrix_Setup.o \
	$(SRCDIR)/SRC/BLOCK/DBMATRIX/DBMatrix_SetupHID.o \
  	$(SRCDIR)/SRC/BLOCK/DBMATRIX/DBMatrix_solve2.o \
	$(SRCDIR)/SRC/BLOCK/DBMATRIX/DBMatrix_solve.o \
	$(SRCDIR)/SRC/BLOCK/DBMATRIX/DBMatrix_TRSM.o \
	$(SRCDIR)/SRC/BLOCK/DBMATRIX/DBPrec_FACTO_TRSM_GEMM_LL.o \
	$(SRCDIR)/SRC/BLOCK/DBMATRIX/DBPrec_FACTO_TRSM_GEMM_RL.o \
	$(SRCDIR)/SRC/BLOCK/DBMATRIX/DBPrec_FACTO_TRSM_GEMM_RL_DROP.o \
	$(SRCDIR)/SRC/BLOCK/DBMATRIX/DBPrec_func.o \
	$(SRCDIR)/SRC/BLOCK/DBMATRIX/DBPrec_pic.o \
	$(SRCDIR)/SRC/BLOCK/DBMATRIX/DBPrec_SchurProd.o \
	$(SRCDIR)/SRC/BLOCK/DBMATRIX/db_struct.o \
	$(SRCDIR)/SRC/BLOCK/DBMATRIX/dumpDBMatrix.o \
	$(SRCDIR)/SRC/BLOCK/DBMATRIX/HIPS_SymbolMatrix.o \
	$(SRCDIR)/SRC/BLOCK/DBMATRIX/PHIDAL_SymbolMatrix.o \
	$(SRCDIR)/SRC/BLOCK/DBMATRIX/print.o \
	$(SRCDIR)/SRC/BLOCK/DBMATRIX/SF_Direct.o \
	$(SRCDIR)/SRC/BLOCK/DBMATRIX/size.o \
	$(SRCDIR)/SRC/BLOCK/DBMATRIX/DBMATRIX_MLICCPrec.o \
	$(SRCDIR)/SRC/BLOCK/DBMATRIX/DBMatrix2PhidalMatrix_drop.o \
	$(SRCDIR)/SRC/BLOCK/DBMATRIX/HIPS_Fgmresd_PH_DB.o \
	$(SRCDIR)/SRC/BLOCK/DBMATRIX/DBMatrix_SetupHID.o \
	$(SRCDIR)/SRC/BLOCK/DBMATRIX/DBMATRIX_MLICCPrec.o \
	$(SRCDIR)/SRC/BLOCK/DBMATRIX/DBMATRIX_MLILUPrec.o \
	$(SRCDIR)/SRC/BLOCK/DBMATRIX/DBPrec_FACTO_TRSM_GEMM_LL.o \
	$(SRCDIR)/SRC/BLOCK/DBMATRIX/DBPrec_FACTO_TRSM_GEMM_LLu.o \
	$(SRCDIR)/SRC/BLOCK/DBMATRIX/DBMatrix_FACTO.o \
	$(SRCDIR)/SRC/BLOCK/DBMATRIX/DBMatrix_FACTO2u.o \
	$(SRCDIR)/SRC/BLOCK/DBMATRIX/DBPrec_FACTO_TRSM_GEMM_RL.o \
	$(SRCDIR)/SRC/BLOCK/DBMATRIX/DBPrec_FACTO_TRSM_GEMM_RLu.o \
	$(SRCDIR)/SRC/BLOCK/DBMATRIX/DBPrec_FACTO_TRSM_GEMM_RL_DROP.o \
	$(SRCDIR)/SRC/BLOCK/DBMATRIX/DBPrec_FACTO_TRSM_GEMM_RL_DROPu.o \
	$(SRCDIR)/SRC/BLOCK/DBMATRIX/HIPS_Fgmresd_DB_PH.o \
	$(SRCDIR)/SRC/BLOCK/DBMATRIX/HIPS_Fgmresd_PH_DB.o \
	$(SRCDIR)/SRC/BLOCK/DBMATRIX/HIPS_Fgmresd_DB.o \
	$(SRCDIR)/SRC/BLOCK/DBMATRIX/HIPS_Fgmresd_DB_DB.o \
	$(SRCDIR)/SRC/BLOCK/DBMATRIX/HIPS_PCG_DB_PH.o \
	$(SRCDIR)/SRC/BLOCK/DBMATRIX/HIPS_PCG_PH_DB.o \
	$(SRCDIR)/SRC/BLOCK/DBMATRIX/HIPS_PCG_DB.o \
	$(SRCDIR)/SRC/BLOCK/DBMATRIX/HIPS_PCG_DB_DB.o \
	$(SRCDIR)/SRC/BLOCK/DBMATRIX/DBMatrix2PhidalMatrix_drop.o



SRCBlockSOLVMATRIX = $(SRCDIR)/SRC/BLOCK/SOLVMATRIX
ObjBlockSOLVMATRIX = $(SRCBlockSOLVMATRIX)/*.o

SRCBlockSOLVMATRIX2 = $(SRCDIR)/SRC/BLOCK/SOLVMATRIX2
ObjBlockSOLVMATRIX2 = $(SRCBlockSOLVMATRIX2)/*.o

#/SRC/IO;
ObjIO = $(SRCDIR)/SRC/IO/*.o

#/SRC/PHIDAL/PARALLEL
SRCPhidalParallel = $(SRCDIR)/SRC/PHIDAL/PARALLEL
ObjPhidalParallel  = $(SRCPhidalParallel)/*.o 

#SRC/BLOCK/PARALLEL
ObjBlockParallel = $(SRCDIR)/SRC/BLOCK/PARALLEL/DBDistrMatrix_func.o $(SRCDIR)/SRC/BLOCK/PARALLEL/DBDistrMatrix_Build.o $(SRCDIR)/SRC/BLOCK/PARALLEL/DBDistrMatrix_Setup.o $(SRCDIR)/SRC/BLOCK/PARALLEL/DBDistrMatrix_Copy.o $(SRCDIR)/SRC/BLOCK/PARALLEL/DBDISTRMATRIX_SEQ.o $(SRCDIR)/SRC/BLOCK/PARALLEL/DBDISTRMATRIX_MLICCPrec.o $(SRCDIR)/SRC/BLOCK/PARALLEL/DBDISTRMATRIX_MLILUPrec.o $(SRCDIR)/SRC/BLOCK/PARALLEL/db_parallel_struct.o $(SRCDIR)/SRC/BLOCK/PARALLEL/DBcomm.o $(SRCDIR)/SRC/BLOCK/PARALLEL/DBDistrMatrix2PhidalDistrMatrix.o $(SRCDIR)/SRC/BLOCK/PARALLEL/DBDISTRMATRIX_PrecSolve.o $(SRCDIR)/SRC/BLOCK/PARALLEL/DBDistrMatrix_solve.o $(SRCDIR)/SRC/BLOCK/PARALLEL/DBDistrMatrix_ICCT.o $(SRCDIR)/SRC/BLOCK/PARALLEL/CellDBDistr.o $(SRCDIR)/SRC/BLOCK/PARALLEL/VS_RowDiv.o $(SRCDIR)/SRC/BLOCK/PARALLEL/HIPS_DistrFgmresd_PH_DB.o $(SRCDIR)/SRC/BLOCK/PARALLEL/HIPS_DistrFgmresd_DB_PH.o $(SRCDIR)/SRC/BLOCK/PARALLEL/HIPS_DistrFgmresd_DB.o $(SRCDIR)/SRC/BLOCK/PARALLEL/HIPS_DistrFgmresd_DB_DB.o $(SRCDIR)/SRC/BLOCK/PARALLEL/DBDistrPrec_SchurProd.o $(SRCDIR)/SRC/BLOCK/PARALLEL/DBDistrMatrix_SetupHID.o $(SRCDIR)/SRC/BLOCK/PARALLEL/HIPS_DistrPCG_PH_DB.o $(SRCDIR)/SRC/BLOCK/PARALLEL/HIPS_DistrPCG_DB_PH.o $(SRCDIR)/SRC/BLOCK/PARALLEL/DBDistrPrec_func.o $(SRCDIR)/SRC/BLOCK/PARALLEL/DBDistrMatrix_FACTO.o $(SRCDIR)/SRC/BLOCK/PARALLEL/DBDistrMatrix_FACTOu.o \
	$(SRCDIR)/SRC/BLOCK/PARALLEL/DBDistrPrec_FACTO_TRSM_GEMM_LL.o \
	$(SRCDIR)/SRC/BLOCK/PARALLEL/DBDistrPrec_FACTO_TRSM_GEMM_RL.o \
	$(SRCDIR)/SRC/BLOCK/PARALLEL/DBDistrPrec_FACTO_TRSM_GEMM_RL_DROP.o \
	$(SRCDIR)/SRC/BLOCK/PARALLEL/DBDistrPrec_FACTO_TRSM_GEMM_LLu.o \
	$(SRCDIR)/SRC/BLOCK/PARALLEL/DBDistrPrec_FACTO_TRSM_GEMM_RLu.o \
	$(SRCDIR)/SRC/BLOCK/PARALLEL/DBDistrPrec_FACTO_TRSM_GEMM_RL_DROPu.o \
	$(SRCDIR)/SRC/BLOCK/PARALLEL/DBDistrPrec_func2.o $(SRCDIR)/SRC/BLOCK/PARALLEL/DBDistrMatrix_solve_sub.o $(SRCDIR)/SRC/BLOCK/PARALLEL/DBcomm2.o $(SRCDIR)/SRC/BLOCK/PARALLEL/DBDISTRMATRIX_MatVec2.o


#SRC/HIPS_WRAPPER
ObjHIPS_WRAPPER= $(SRCDIR)/SRC/HIPS_WRAPPER/*.o

#SRC/BLAS
ObjBLAS= $(SRCDIR)/SRC/BLAS/zsymv.o $(SRCDIR)/SRC/BLAS/csymv.o

#SRC/SPKIT pas besoin
ObjSPKIT= $(SRCDIR)/SRC/SPKIT/skitfc.o $(SRCDIR)/SRC/SPKIT/skitf.o


# definition of obj for librairies
OBJLIBHIPSSEQUENTIAL = $(ObjPhidalCommon) $(ObjPhidalOrdering) $(ObjPhidalSequential) $(ObjBlockDBMATRIX) $(ObjBlockSOLVMATRIX) $(ObjBlockSOLVMATRIX2) $(ObjHIPS_WRAPPER) $(ObjBLAS)
OBJLIBHIPS = $(ObjPhidalCommon) $(ObjPhidalOrdering) $(ObjPhidalSequential) $(ObjBlockDBMATRIX) $(ObjBlockSOLVMATRIX) $(ObjBlockSOLVMATRIX2) $(ObjPhidalParallel) $(ObjBlockParallel) $(ObjHIPS_WRAPPER) $(ObjBLAS)

OBJLIBIO = $(ObjIO)

ifndef FC
OBJLIBSPKIT = 
else
OBJLIBSPKIT = $(ObjSPKIT)
endif

hips: FAIRE  

FAIT:  
	$(MAKE) $(SRCDIR)/$(INSTALL)  WHERE
	touch FAIT

FAIRE: FAIT 

makefile.inc: ../Makefile   ../../config.status makefile-hips.inc Makefile     
	../../config.status  --file="makefile.inc:makefile-hips.inc"


install:$(SRCDIR)/$(INSTALL)

$(SRCDIR)/$(INSTALL): $(SRCDIR)/tag-ff++
	cp makefile.inc $(SRCDIR)
	cd $(SRCDIR); $(MAKE)
	-rm $(SRCDIR)/LIB/*.a
	$(AR) $(ARFLAGS) $(SRCDIR)/LIB/libhips.a $(OBJLIBHIPS)
	$(RANLIB) $(SRCDIR)/LIB/libhips.a
	$(AR) $(ARFLAGS) $(SRCDIR)/LIB/libhipssequential.a $(OBJLIBHIPSSEQUENTIAL)
	$(RANLIB) $(SRCDIR)/LIB/libhipssequential.a
	$(AR) $(ARFLAGS) $(SRCDIR)/LIB/libio.a $(OBJLIBIO) 
	$(RANLIB) $(SRCDIR)/LIB/libio.a
	if [ -n "$(OBJLIBSPKIT)" ] ; then \
	$(AR) $(ARFLAGS) $(SRCDIR)/LIB/libspkit.a $(OBJLIBSPKIT) ;\
	$(RANLIB) $(SRCDIR)/LIB/libspkit.a;\
	fi
	mkdir -p $(SRCDIR)/$(INSTALL)/include/hips
	cp $(SRCDIR)/LIB/*.h $(SRCDIR)/$(INSTALL)/include/hips/
	mkdir -p $(SRCDIR)/$(INSTALL)/lib/hips
	cp $(SRCDIR)/LIB/*.a $(SRCDIR)/$(INSTALL)/lib/hips
WHERE:
	echo hips LD '-L@DIR@/lib/hips -lio -lhips'  >$(SRCDIR)/$(INSTALL)/lib/WHERE.hips
	echo hips INCLUDE '-I@DIR@/include/hips'>> $(SRCDIR)/$(INSTALL)/lib/WHERE.hips
	echo hipssequential LD '-L@DIR@/lib/hips -lio -lhipssequential'  >>$(SRCDIR)/$(INSTALL)/lib/WHERE.hips
	echo hipssequential INCLUDE '-I@DIR@/include/hips'>> $(SRCDIR)/$(INSTALL)/lib/WHERE.hips


$(SRCDIR)/tag-ff++: $(PACKAGE)
	tar xvzf $(PACKAGE)
	cd $(SRCDIR); patch -p1 <../hips-1.2b-rc4.patch
#	cd $(SRCDIR)/SRC/INCLUDE/; patch -p1 < ../../../hips-1.2b-rc4-1.patch;
#	cd $(SRCDIR)/SRC/IO/; patch -p1 < ../../../hips-1.2b-rc4-2.patch;
#	cd $(SRCDIR); for d in SRC/INCLUDE SRC/IO SRC/SPKIT TESTS/PARALLEL; do \
#	  cp $$d/makefile $$d/makefile.orig; \
#	  sed 's/\$$(CC).*-E/$$(CPP) $$(CFLAGS) $$(CC_OPT)/' >$$d/makefile <$$d/makefile.orig;\
#	done
#	cp  SRC_SPKIT_makefile $(SRCDIR)/SRC/SPKIT/makefile
	touch $(SRCDIR)/tag-ff++

$(PACKAGE):
	-mkdir $(DIRPKG)
	cd $(DIRPKG);$(WGET)   $(SERVER)/`basename $(PACKAGE)`

clean-local:
	-cd $(SRCDIR) && $(MAKE) clean  -C $(SRCDIR)

clean: clean-local
	rm makefile.inc
	-rm $(SRCDIR)/$(INSTALL)/WHERE/hips
	-rm $(SRCDIR)/$(INSTALL)/lib/hips/*
	-rm -rf $(SRCDIR)/$(INSTALL)/lib/hips/
	-rm $(SRCDIR)/$(INSTALL)/include/hips/*
	-rm -rf $(SRCDIR)/$(INSTALL)/include/hips/
	-rm -rf $(SRCDIR)
	rm FAIT
.PHONY:$(SRCDIR)/$(INSTALL)