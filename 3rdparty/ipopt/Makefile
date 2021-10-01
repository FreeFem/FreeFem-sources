# Downloading and compiling extra libraries
# -----------------------------------------

include Makefile.inc

all-local: Ipopt

MAKEFLAGS += -j1

# Downloading and compiling mumps
# ------------------------------

DIR=$(DOWNLOADFF)/ipopt
DIRPKG=../pkg
SRCDIR=Ipopt-$(VERSION)
PACKAGE=$(DIRPKG)/Ipopt-$(VERSION).tgz
INSTALL=../..
VERSION=3.12.4
# 3.10.2
URL=http://www.coin-or.org/download/source/Ipopt
FHSL=#$(DIRPKG)/ddeps.f $(DIRPKG)/ma27ad.f  $(DIRPKG)/mc19d.f 
WHERE_MUMPSSEQ=$(wildcard ../lib/WHERE.mumpsseq)

ifeq ($(WHERE_MUMPSSEQ),)
LIBMUMPS=-ldmumps_seq -lzmumps_seq  -lmumps_common_seq  -lpord_seq -lmpiseq_seq 
INCMUMPS=/usr/include/mumps_seq
else
LIBMUMPS=-L$(DOWNLOADFF)/lib  -ldmumpsFREEFEM-SEQ -lzmumpsFREEFEM-SEQ -lmumps_commonFREEFEM-SEQ -lpordFREEFEM-SEQ -lmpiseqFREEFEM-SEQ
INCMUMPS=$(DOWNLOADFF)/include/mumps_seq
endif

#/Ipopt-3.10.2.tgz
Ipopt: $(SRCDIR)/FAIRE
#  --enable-static  --disable-shared  

# ALH - 6/1/14 - We need to specify FLIBS explicitely because Ipopt
# fails to guess FLIBS correctly on Windows+Cygwin+Mingw64, see [[file:Makefile.inc.in::FLIBS]]

$(SRCDIR)/FAIT: $(SRCDIR)/tag-tar
#	FFCS - disable dependency tracking like in FFCS itself for MinGW compilation (problem with backslashes, see
#	[[file:../../../../configure.ac::enable_dependency_tracking]])
	cd $(SRCDIR) ; \
	 ./configure --disable-dependency-tracking \
		   --disable-shared --enable-static \
                  --with-mumps='$(LIBMUMPS)' \
	          --without-hsl \
	          --with-mumps-incdir='$(INCMUMPS)' \
	  CXX='$(CXX)' CXXFLAGS='$(CXXFLAGS) -I$(INCMUMPS)' \
	  CC='$(CC)' CFLAGS='$(CFLAGS) -I$(INCMUMPS)' \
	  F77='$(FC)' FFLAGS='$(FCFLAGS)' \
	  FLIBS='$(FLIBS)' \
	  CXXCPP=''  CPP='' MPICC='' MPICXX='' MPIFC='' \
	 --with-blas-lib='$(LIBBLAS)' --with-lapack='$(LIBLAPACK)' --prefix='$(DOWNLOADFF)'
	 -rm -rf ../include/coin
	touch $(SRCDIR)/FAIT	

# FFCS - avoid remaking install every time
install.done: Makefile $(SRCDIR)/FAIT
	$(MAKE)  -C  $(SRCDIR) install	
	touch $@
clean-local::
	-rm *.done

# FFCS - install and WHERE need to be sequential
WHERE.done: install.done Makefile
	echo ipopt LD -L@DIR@/lib   -lipopt   >$(SRCDIR)/$(INSTALL)/lib/WHERE.Ipopt;
	echo ipopt INCLUDE -I@DIR@/include  >> $(SRCDIR)/$(INSTALL)/lib/WHERE.Ipopt ;
	touch $@

Makefile.inc:
	../../config.status  --file="Makefile.inc:Makefile.inc.in"

# FFCS - install and WHERE need to be sequential
$(SRCDIR)/FAIRE: install.done WHERE.done
	touch $@


$(SRCDIR)/$(INSTALL): $(SRCDIR)/tag-tar

$(SRCDIR)/tag-tar:$(PACKAGE) $(FHSL) 
	tar xzf $(PACKAGE)
#	patch -p0 <patch-IpBlas
	touch $(SRCDIR)/tag-tar

$(PACKAGE):
	../getall -o Ipopt -a

clean-local::
	rm -rf  $(SRCDIR)  *~ Makefile.inc

clean: clean-local
	-rm -rf ../include/coin
	-rm ../lib/libipopt* ../lib/liblcoinhsl* 
#FH	-rm $(PACKAGE)

.PHONY:$(SRCDIR)/$(INSTALL)
#.NOTPARALLEL:$(SRCDIR)/$(INSTALL)
