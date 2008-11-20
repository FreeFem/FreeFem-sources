# Downloading and compiling extra libraries
# -----------------------------------------
# $Id$

all-local: superlu

# Downloading and compiling Tetgen
# ------------------------------

# SUPERLU information
#http://crd.lbl.gov/~xiaoye/SuperLU/superlu_3.1.tar.gz
DIRPKG=../pkg
SRCDIR=SuperLU_3.1
PACKAGE=superlu_3.1.tar.gz
PACKAGE_PATH=$(DIRPKG)/$(PACKAGE)
SERVER=http://crd.lbl.gov/~xiaoye/SuperLU/
INSTALL=../..
SUPERLU_VERSION=3.1
SUPERLULIB=../$(INSTALL)/lib/libsuperlu_3.1.a
include  ff-FLAGS

### LAPACK 
LAAUX 	= lsame.o xerbla.o
SLASRC 	= slacon.o 
DLASRC	= dlacon.o
CLASRC	= clacon.o scsum1.o icmax1.o
ZLASRC	= zlacon.o dzsum1.o izmax1.o
SCLAUX 	= slamch.o
DZLAUX 	= dlamch.o

### SuperLU 
ALLAUX 	= superlu_timer.o util.o memory.o get_perm_c.o mmd.o \
	  sp_coletree.o sp_preorder.o sp_ienv.o relax_snode.o \
	  heap_relax_snode.o colamd.o

SLUSRC = \
	sgssv.o sgssvx.o \
	ssp_blas2.o ssp_blas3.o sgscon.o  \
	slangs.o sgsequ.o slaqgs.o spivotgrowth.o \
	sgsrfs.o sgstrf.o sgstrs.o scopy_to_ucol.o \
	ssnode_dfs.o ssnode_bmod.o \
	spanel_dfs.o spanel_bmod.o sreadhb.o \
	scolumn_dfs.o scolumn_bmod.o spivotL.o spruneL.o \
	smemory.o sutil.o smyblas2.o

DLUSRC = \
	dgssv.o dgssvx.o \
	dsp_blas2.o dsp_blas3.o dgscon.o \
	dlangs.o dgsequ.o dlaqgs.o dpivotgrowth.o  \
	dgsrfs.o dgstrf.o dgstrs.o dcopy_to_ucol.o \
	dsnode_dfs.o dsnode_bmod.o \
	dpanel_dfs.o dpanel_bmod.o dreadhb.o \
	dcolumn_dfs.o dcolumn_bmod.o dpivotL.o dpruneL.o \
	dmemory.o dutil.o dmyblas2.o ## dgstrsL.o dgstrsU.o

CLUSRC = \
	scomplex.o cgssv.o cgssvx.o csp_blas2.o csp_blas3.o cgscon.o \
	clangs.o cgsequ.o claqgs.o cpivotgrowth.o  \
	cgsrfs.o cgstrf.o cgstrs.o ccopy_to_ucol.o \
	csnode_dfs.o csnode_bmod.o \
	cpanel_dfs.o cpanel_bmod.o creadhb.o \
	ccolumn_dfs.o ccolumn_bmod.o cpivotL.o cpruneL.o \
	cmemory.o cutil.o cmyblas2.o

ZLUSRC = \
	dcomplex.o zgssv.o zgssvx.o zsp_blas2.o zsp_blas3.o zgscon.o \
	zlangs.o zgsequ.o zlaqgs.o zpivotgrowth.o  \
	zgsrfs.o zgstrf.o zgstrs.o zcopy_to_ucol.o \
	zsnode_dfs.o zsnode_bmod.o \
	zpanel_dfs.o zpanel_bmod.o zreadhb.o \
	zcolumn_dfs.o zcolumn_bmod.o zpivotL.o zpruneL.o \
	zmemory.o zutil.o zmyblas2.o


OBJ_SUPERLU = $(ALLAUX) $(LAAUX) \
 $(SLUSRC)  $(SLASRC) $(SCLAUX) \
 $(DLUSRC)  $(DLASRC) $(DZLAUX) \
 $(CLUSRC)  $(CLASRC) \
 $(ZLUSRC)  $(ZLASRC) 

superlu: $(SRCDIR)/SRC/$(SUPERLULIB) 

$(SRCDIR)/SRC/$(SUPERLULIB): $(SRCDIR)
	cd $(SRCDIR)/SRC; $(MAKE) -f ../../ff-FLAGS  $(OBJ_SUPERLU)
	cd $(SRCDIR)/SRC; $(AR) $(ARFLAGS) $(SUPERLULIB)  $(OBJ_SUPERLU)
	cd $(SRCDIR)/SRC; $(RANLIB) $(SUPERLULIB) 
	cp $(SRCDIR)/SRC/s*h  $(SRCDIR)/$(INSTALL)/include


$(SRCDIR)/$(INSTALL): 
	mkdir $(SRCDIR)/$(INSTALL)


$(SRCDIR): $(PACKAGE_PATH)
	tar xvzf $(PACKAGE_PATH)
	touch $(SRCDIR)

$(PACKAGE_PATH):
	-mkdir $(DIRPKG);
	cd $(DIRPKG);$(WGET)  $(SERVER)$(PACKAGE)

clean-local:
	-rm -rf superlu*
clean:
	-rm $(SRCDIR)/SRC/$(SUPERLULIB)
	-rm $(SRCDIR)/$(INSTALL)/include/slu*h 
	-rm $(SRCDIR)/$(INSTALL)/include/supermatrix.h
	-rm -r $(SRCDIR) ff-FLAGS


ff-FLAGS: ../Makefile  Makefile 
	grep '^CXX *=' ../Makefile >>$@
	grep '^CC *=' ../Makefile >>$@
	grep '^CFLAGS *=' ../Makefile >>$@
	grep '^BLAS[A-Z ]*='  ../Makefile|grep =  >>$@
	grep '^AR *='  ../Makefile|grep =  >>$@
	grep '^ARFLAGS *='  ../Makefile|grep =  >>$@
	grep '^RANLIB *='    ../Makefile|grep =  >>$@
	grep '^WGET *='    ../Makefile|grep =  >>$@

