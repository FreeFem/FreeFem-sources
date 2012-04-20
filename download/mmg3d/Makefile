# Downloading and compiling extra libraries
# -----------------------------------------

# $Id$
all-local: mmg3d

include ff-flags

# Downloading and compiling mmg3d
# -------------------------------
# 
DIRPKG= ../pkg
SRCDIR= ./mmg3d4# lib#-$(mmg3d_VERSION)
PACKAGE=$(DIRPKG)/mmg3d4.0.tgz #mmg3dlib.tar.gz #-$(mmg3d_VERSION).tar.gz
INSTALL=..
mmg3d_VERSION=
# mmg3d pas sur internet
LIBMMG3D=$(INSTALL)/lib/libmmg3d-v4.a
OPT=4

OBJS= analar.o	chkmsh.o	hash.o		memory.o	optcte.o	outqua.o	simu44.o	swap44.o	zaldy.o \
analarcutting.o	chrono.o	heap.o		mmg3d1.o	optlap.o	pattern.o	simu56.o	swap56.o \
baryct.o	colpoi.o	inout.o		mmg3d4.o	optlen.o	quality.o	simu68.o	swap68.o \
boulep.o	coquil.o	length.o	mmg3d9.o	optlentet.o	queue.o		simu710.o	swap710.o \
bucket.o	cutelt.o	librnbg.o	movevertex.o	optra4.o	ratio.o		solmap.o	swapar.o \
		delaunay.o	locate.o	optbdry.o	opttet.o	scalem.o	spledg.o	swaptet.o \
cenrad.o	eigenv.o	matrix.o	optcoq.o	opttyp.o	simu23.o	swap23.o	typelt.o 

OBJSNOP =  cendel.o swapar.o

mmg3d: build-4


build-4: 
	-$(MAKE) tag-tar
	if [ -d mmg3d4  ] ; then  $(MAKE)  $(LIBMMG3D) install-4 WHERE ; else echo "Sorry no source " ; fi
$(LIBMMG3D):
	cd mmg3d4/build/sources/; $(MAKE) CC='$(CC)' CFLAGS='$(CNOFLAGS)'  $(OBJSNOP) 
	cd mmg3d4/build/sources/; $(MAKE) CC='$(CC)' CFLAGS='$(CFLAGS)'  $(OBJS) 
	cd mmg3d4/build/sources/mmg3dmain; $(CC) -c $(CNOFLAGS)   mmg3d.c -I..
	$(AR) $(ARFLAGS) $(LIBMMG3D) mmg3d4/build/sources/*.o mmg3d4/build/sources/mmg3dmain/mmg3d.o
	-$(CC) $(CNOFLAGS) mmg3d4/build/sources/mmg3dmain/mmg3d.o   mmg3d4/build/sources/*.o -o ../bin/mmg3d $(STD_LIBS)
	touch build-4 FAIT-4
install-4:
	-mkdir ../include/mmg3d-v4/		
	cp mmg3d4/build/sources/*.h ../include/mmg3d-v4/
	touch FAIT-4

mmg3d-4:$(PACKAGE)

$(PACKAGE):
	-mkdir $(DIRPKG);
	cd $(DIRPKG);$(WGET) http://www.math.u-bordeaux1.fr/~dobj/logiciels/download/mmg3d4.0.tgz 
#http://www.math.u-bordeaux1.fr/~dobj/logiciels/download.php?file=`basename $(PACKAGE)`


install:install-4 WHERE

WHERE: 
	@-if [ -f FAIT-4 ] ; then \
	make install-4;  \
	echo mmg3d-v4  LD -L@DIR@/lib -lmmg3d-v4  >../lib/WHERE.mmg3d ;\
	echo mmg3d-v4 INCLUDE -I@DIR@/include/mmg3d-v4>> ../lib/WHERE.mmg3d ;\
	echo build WHERE ./lib/WHERE.mmg3d ;\
	fi


FAIRE: FAIT-4 install-4

tag-tar: $(PACKAGE)
	if tar tfz $(PACKAGE)  mmg3d4/build/sources/mmg3d.c 2>/dev/null  1>/dev/null  ;\
	 then  rm $(PACKAGE);   $(MAKe) $(PACKAGE); fi	
	tar xvzf $(PACKAGE)
	touch mmg3d4/build/sources/dataff.h 	
	cd mmg3d4;patch -p1 <../patch-mmg3dv4.diff
	cat </dev/null >mmg3d4/build/sources/mmg3dConfig.h
	touch tag-tar
#	cp makefile-mmg3d.inc $(SRCDIR)/makefile



clean: 
	-rm ff-flags
	-rm -rf $(SRCDIR)/$(INSTALL)/lib/libmmg3d*.a  ../lib/libmmg3d-v4.a
	-rm -rf $(SRCDIR)/$(INSTALL)/include/mmg3d-v4
	-rm -r $(SRCDIR)
	-rm FAIT* mmg* flags-* build-4 tag-tar

ff-flags: ../Makefile Makefile
	grep 'abs_top_builddir *=' ../Makefile >> ff-flags
	grep 'CC *=' ../Makefile >> ff-flags
	grep 'CFLAGS *=' ../Makefile >> ff-flags
	grep 'LDFLAGS *=' ../Makefile >> ff-flags
	grep 'AR *=' ../Makefile >> ff-flags
	grep 'ARFLAGS *=' ../Makefile >> ff-flags
	grep 'RANLIB *=' ../Makefile >> ff-flags
	grep 'WGET *=' ../Makefile >> ff-flags
	grep 'STD_LIBS *=' ../Makefile >> ff-flags

.PHONY:$(SRCDIR)/$(INSTALL)