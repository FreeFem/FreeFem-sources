# Downloading and compiling extra libraries
# -----------------------------------------

include SLmake.inc

all-local: scalapack

# FFCS: parallel compilation crashes on Win32 (same archive updated from 2 different parallel makes)
#MAKEFLAGS=-j 1
.NOTPARALLEL:
# Downloading and compiling scalapack
# ------------------------------
# http://www.netlib.org/scalapack/
# Hips information
DIRPKG=../pkg
SRCDIR=scalapack-2.0.2
PACKAGE1=$(DIRPKG)/scalapack-2.0.2.tgz
SERVER=http://www.netlib.org/scalapack/
INSTALL=../..

# FFCS: shorten argument length to avoid "sh: ../ar: Argument list too long" under Cygwin
DIR1 = TOOLS/LAPACK/*.o
DIR2 = TOOLS/*.o
DIR3 = PBLAS/SRC/PBBLAS/*.o
DIR4 = PBLAS/SRC/PTZBLAS/*.o
DIR5 = PBLAS/SRC/PTOOLS/*.o
DIR6 = PBLAS/SRC/*.o
DIR7 = REDIST/SRC/*.o
DIR8 = SRC/*.o
DIR9 = BLACS/SRC/*.o BLACS/SRC/*.oo

scalapack: FAIRE

# FFCS: shorten argument length to avoid "sh: ../ar: Argument list too long" under Cygwin
$(SRCDIR)/FAIT-202:$(SRCDIR)/tag-tar
	cp SLmake.inc $(SRCDIR)
	cd $(SRCDIR);$(MAKE) MAKEFLAGS=-j1 lib
	rm $(SRCDIR)/libscalapack.a


	cd $(SRCDIR) && $(ARCH) $(ARCHFLAGS) libscalapack.a $(DIR1) $(DIR2) $(DIR3) $(DIR4) $(DIR5) $(DIR6) $(DIR7) $(DIR8) $(DIR9)
	$(RANLIB) $(SRCDIR)/libscalapack.a
	touch $(SRCDIR)/FAIT-202
install:
	mkdir -p ../include
	cp $(SRCDIR)/SRC/*.h ../include
	mkdir -p ../lib
	cp $(SRCDIR)/libscalapack.a ../lib

WHERE:
	if [ -f $(SRCDIR)/FAIT-202 ] ; then \
	make install;  \
	echo scalapack LD -L@DIR@/lib  -lscalapack   >../lib/WHERE.scalapack ;\
	echo scalapack INCLUDE -I@DIR@/include  >> ../lib/WHERE.scalapack ;\
	fi

FAIRE: $(SRCDIR)/FAIT-202
	$(MAKE) WHERE
	touch FAIRE 

SLmake.inc: ../../config.status	Makefile SLmake-scalapack.inc
	../../config.status  --file="SLmake.inc:SLmake-scalapack.inc"

$(SRCDIR)/tag-tar: $(PACKAGE1)
	gunzip -c $(PACKAGE1) | tar xvf -
	#  fto pach 
	patch -p0 <scalapack.patch
	touch $(SRCDIR)/tag-tar

$(PACKAGE1):
	../getall -o ScaLAPACK -a

clean:
#	FFCS - need to clean completely even in case of error
	-rm SLmake.inc FAIRE FAIT
#	FFCS - make sure that all directories are cleaned. Thisis especially important under Windows because there is no
#	compilation dependencies control there (see
#	[[file:c:/cygwin/home/alh/ffcs/dist/configure.ac::dependency_tracking]])
	-rm -rf scalapack-*
	-rm config.log
#	-rm $(PACKAGE1)

.PHONY:$(SRCDIR)/$(INSTALL)
.NOTPARALLEL: