# Downloading and compiling extra libraries
# -----------------------------------------

all-local: mmg3d

include ff-flags

# Downloading and compiling mmg3d
# -------------------------------
#
DIRPKG= ../pkg
SRCDIR= ./mmg3d4
PACKAGE=$(DIRPKG)/mmg3d4.0.tgz
INSTALL=..
mmg3d_VERSION=
# mmg3d pas sur internet
LIBMMG3D=$(INSTALL)/lib/libmmg3d-v4.a
OPT=4
# size of the PKG file ( this file change See Cecile.) FFCS - 19/2/13 - the test for the file size (`stat -f "%z" file`)
# is not portable from MacOS, so just leave it out in FFCS (SIZEPKG is still useful to force a remake when the package
# changes)
SIZEPKG=158547

OBJS= analar.o	chkmsh.o	hash.o		memory.o	optcte.o	outqua.o	simu44.o	swap44.o	zaldy.o \
analarcutting.o	chrono.o	heap.o		mmg3d1.o	optlap.o	pattern.o	simu56.o	swap56.o \
baryct.o	colpoi.o	inout.o		mmg3d4.o	optlen.o	quality.o	simu68.o	swap68.o \
boulep.o	coquil.o	length.o	mmg3d9.o	optlentet.o	queue.o		simu710.o	swap710.o \
bucket.o	cutelt.o	librnbg.o	movevertex.o	optra4.o	ratio.o		solmap.o	swapar.o \
		delaunay.o	locate.o	optbdry.o	opttet.o	scalem.o	spledg.o	swaptet.o \
cenrad.o	eigenv.o	matrix.o	optcoq.o	opttyp.o	simu23.o	swap23.o	typelt.o

OBJSNOP =  cendel.o swapar.o

# FFCS - simplify Makefile structure
mmg3d:WHERE.done

mmg3d4/FAIT-4.done:tag-tar-$(SIZEPKG)
	echo "#define COMPIL " '"' `date` '"' > ./mmg3d4/build/sources/compil.date
	cd mmg3d4/build/sources/; $(MAKE) CC='$(CC)' CFLAGS='$(CNOFLAGS) -g'   $(OBJSNOP)
	cd mmg3d4/build/sources/; $(MAKE) CC='$(CC)' CFLAGS='$(CFLAGS) -g'   $(OBJS)
	cd mmg3d4/build/sources/mmg3dmain; $(CC) -c $(CFLAGS)   mmg3d.c -I..
	$(AR) $(ARFLAGS) $(LIBMMG3D) mmg3d4/build/sources/*.o mmg3d4/build/sources/mmg3dmain/mmg3d.o
#
#	FFCS - ranlib required on Windows 7 64 bits
#
	$(RANLIB) $(LIBMMG3D)
	-$(CC) $(CNOFLAGS) mmg3d4/build/sources/mmg3dmain/mmg3d.o   mmg3d4/build/sources/*.o -o ../bin/mmg3d $(STD_LIBS)
	touch mmg3d4/FAIT-4.done

# FFCS - simplify Makefile structure
install-4.done:mmg3d4/FAIT-4.done
	-mkdir ../include/mmg3d-v4/
	cp mmg3d4/build/sources/*.h ../include/mmg3d-v4/
	touch $@
clean::
	-rm install-4.done

mmg3d-4:$(PACKAGE)

install:install-4.done WHERE.done

# FFCS - keep simplest makefile structure for automatic recompilations
WHERE.done: install-4.done
	echo mmg3d-v4  LD -L@DIR@/lib -lmmg3d-v4  >../lib/WHERE.mmg3d ;
	echo mmg3d-v4 INCLUDE -I@DIR@/include/mmg3d-v4>> ../lib/WHERE.mmg3d ;
	echo build WHERE ./lib/WHERE.mmg3d ;
	touch $@
clean::
	-rm WHERE.done

FAIRE: mmg3d4/FAIT-4.done install-4.done

# FFCS - keep it simple
tag-tar-$(SIZEPKG): $(PACKAGE) patch-mmg3dv4.diff
	-rm -rf  mmg3d4
	tar xvzf $(PACKAGE)
	-rm mmg3d4/build/libexamples/main.
	touch mmg3d4/build/sources/dataff.h
#
#	ALH - clean-up all CR/LF to make patching more successful
#
	../../bin/build/cleancrlf mmg3d4
#
	cd mmg3d4;patch -p1 <../patch-mmg3dv4.diff
	cat </dev/null >mmg3d4/build/sources/mmg3dConfig.h
	touch tag-tar-$(SIZEPKG)
#	cp makefile-mmg3d.inc $(SRCDIR)/makefile

$(PACKAGE):
	../getall -o MMG3D -a

clean::
	-rm ff-flags
#
#	FFCS - make sure that all directories are cleaned. Thisis especially important under Windows because there is no
#	compilation dependencies control there (see
#	[[file:c:/cygwin/home/alh/ffcs/dist/configure.ac::dependency_tracking]])
#
	-rm -r mmg3d4
	-rm FAIT* mmg* flags-* tag-tar*
#FH 	-rm $(PACKAGE)

ff-flags: ../Makefile Makefile
	grep 'abs_top_builddir *=' ../Makefile > ff-flags
	grep 'CC *=' ../Makefile >> ff-flags
	grep 'CFLAGS *=' ../Makefile >> ff-flags
	grep 'LDFLAGS *=' ../Makefile >> ff-flags
	grep 'AR *=' ../Makefile >> ff-flags
	grep 'ARFLAGS *=' ../Makefile >> ff-flags
	grep 'RANLIB *=' ../Makefile >> ff-flags
	grep 'WGET *=' ../Makefile >> ff-flags
	grep 'STD_LIBS *=' ../Makefile >> ff-flags
	grep 'CNOFLAGS *=' ../Makefile >> ff-flags
.PHONY:$(SRCDIR)/$(INSTALL)
