# Downloading and compiling extra libraries
# -----------------------------------------

include Makefile.in

all-local: parms

# Downloading and compiling parms
# -------------------------------
# http://www-users.cs.umn.edu/~saad/software/pARMS/pARMS_2.2.php
# Parms information
DIRPKG=../pkg
SRCDIR=pARMS_$(parms_VERSION)
PACKAGE=$(DIRPKG)/pARMS_$(parms_VERSION).tar.gz
SERVER=http://www-users.cs.umn.edu/~saad/software/pARMS/pARMS_2.2.php
INSTALL=../..
parms_VERSION=2.2

parms: FAIRE

# FFCS: more dependencies for parallel builds
$(SRCDIR)/FAIT:$(SRCDIR)
	$(MAKE) $(SRCDIR)/$(INSTALL)  WHERE
	touch $(SRCDIR)/FAIT
install:$(SRCDIR)/$(INSTALL)

FAIRE: $(SRCDIR)/FAIT

Makefile.in: ../../config.status	Makefile makefile-parms.in
	../../config.status  --file="Makefile.in:makefile-parms.in"
	case `uname` in *CYGWIN*)  cp Makefile.in Makefile.inn; sed "s/COPTIONS =/COPTIONS = -D__VC__/" <Makefile.inn > Makefile.in; rm Makefile.inn ;; esac 

$(SRCDIR)/$(INSTALL): $(SRCDIR)
	cp Makefile.in $(SRCDIR)/makefile.in
	cd $(SRCDIR); $(MAKE) ./LIB/libparms$(DBG).a
	mkdir -p ../include/parms
	cp $(SRCDIR)/INCLUDE/*.h ../include/parms/
	cp $(SRCDIR)/LIB/*.a ../lib/

$(SRCDIR): $(PACKAGE)
	tar xvzf $(PACKAGE)
	touch $(SRCDIR)

$(PACKAGE):
	../getall -o pARMS -a

clean:
	-make cleanall -C $(SRCDIR)/src
	-rm Makefile.in ../lib/WHERE.parms
	-rm ../lib/libparms.*
	-rm -rf ../include/parms
	-rm -rf $(SRCDIR)
#
#	FFCS - make sure that all directories are cleaned. This is especially important under Windows because there is no
#	compilation dependencies control there (see
#	[[file:c:/cygwin/home/alh/ffcs/dist/configure.ac::dependency_tracking]])
#
	-rm -rf pARMS_*
	-rm config.log
#	-rm $(PACKAGE)

WHERE:
	echo parms LD -L@DIR@/lib -lparms$(DBG)  >../lib/WHERE.parms
	echo parms INCLUDE -I@DIR@/include/parms >> ../lib/WHERE.parms

.PHONY:$(SRCDIR)/$(INSTALL)

