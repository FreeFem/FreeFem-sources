# Downloading and compiling extra libraries
# -----------------------------------------

# $Id$
include Makefile.inc

all-local: Ipopt

# Downloading and compiling mumps
# ------------------------------

DIR=$(DOWNLOADFF)/ipopt
DIRPKG=../pkg
SRCDIR=Ipopt-$(VERSION)
PACKAGE=$(DIRPKG)/Ipopt-$(VERSION).tgz
INSTALL=../..
VERSION=3.10.2
# 3.10.2
URL=http://www.coin-or.org/download/source/Ipopt
FHSL=#$(DIRPKG)/ddeps.f $(DIRPKG)/ma27ad.f  $(DIRPKG)/mc19d.f 
LIBMUMPS=-L$(DOWNLOADFF)/lib  -ldmumpsFREEFEM-SEQ -lzmumpsFREEFEM-SEQ -lmumps_commonFREEFEM-SEQ -lpordFREEFEM-SEQ -lmpiseqFREEFEM-SEQ
INCMUMPS=$(DOWNLOADFF)/include/libseq 
#/Ipopt-3.10.2.tgz
Ipopt: $(SRCDIR)/FAIRE

$(SRCDIR)/FAIT: $(SRCDIR)/tag-tar
	cd $(SRCDIR) ; \
	 ./configure  --enable-static  --disable-shared  \
                  --with-mumps='$(LIBMUMPS)' \
	          --without-hsl \
	          --with-mumps-incdir='$(INCMUMPS)' \
	  CXX='$(CXX)' CXXFLAGS='$(CXXFLAGS)' \
	  F77='$(FC)' FCFLAGS='$(FCFLAGS)' \
	 --with-blas='$(LIBBLAS)' --with-lapack='$(LIBLAPACK)' --prefix='$(DOWNLOADFF)'
	touch $(SRCDIR)/FAIT	
install: $(SRCDIR)/FAIT
	$(MAKE)  -C  $(SRCDIR) install	

WHERE: 
	if [ -f $(SRCDIR)/FAIT  ] ;then \
	$(MAKE)  -C  $(SRCDIR) install ; \
	echo Ipopt LD -L@DIR@/lib   -lipopt   >$(SRCDIR)/$(INSTALL)/lib/WHERE.Ipopt;\
	echo Ipopt INCLUDE -I@DIR@/include/coin  >> $(SRCDIR)/$(INSTALL)/lib/WHERE.Ipopt ;\
	fi

Makefile.inc:
	../../config.status  --file="Makefile.inc:Makefile.inc.in"

$(SRCDIR)/FAIRE:
	$(MAKE) install WHERE 
	touch $(SRCDIR)/FAIRE


$(SRCDIR)/$(INSTALL): $(SRCDIR)/tag-tar

$(SRCDIR)/tag-tar:$(PACKAGE) $(FHSL) 
	tar xvzf $(PACKAGE)
	#cp $(FHSL) $(SRCDIR)/ThirdParty/HSL/
	touch $(SRCDIR)/tag-tar

$(PACKAGE):
	cd `dirname $@`; $(WGET)   $(URL)/`basename $(PACKAGE)`
clean-local:
	rm -rf  $(SRCDIR)  *~ Makefile.inc

clean: clean-local
	-rm -rf ../include/coin
	-rm ../lib/libipopt* ../lib/liblcoinhsl* 


.PHONY:$(SRCDIR)/$(INSTALL)