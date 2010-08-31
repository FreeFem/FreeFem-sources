# Downloading and compiling extra libraries
# -----------------------------------------

# $Id$
include config.in

all-local: pastix

# Downloading and compiling pastix
# http://gforge.inria.fr/frs/download.php/21873/
# ------------------------------
DIRPKG   = ../pkg
SRCDIR   = pastix_release_$(VERSION)
PACKAGE1 = $(DIRPKG)/pastix_release_$(VERSION).tar.bz2
SERVER   = http://gforge.inria.fr/frs/download.php/21873/
INSTALL  = ../..
VERSION  = 2200

VERSIONTYPE=$(VERSIONINT)$(VERSIONPRC)$(VERSIONFLT)
VERSIONTYPECOMPLEX=$(VERSIONINT)_complex
VERSIONNAME=$(VERSIONBIT)$(VERSIONMPI)$(VERSIONSMP)$(VERSIONBUB)$(VERSIONTYPE)$(VERSIONORD)_$(HOSTARCH)
PASTIX_DIR = $(abs_top_builddir)/download/pastix/$(SRCDIR)/install

ifneq (,$(findstring libtool,$(ARPROG)))
	OBJLIBTOOL = $(SRCDIR)/src/sopalin/obj/$(HOSTARCH)/*.o $(SRCDIR)/src/blend/obj/$(HOSTARCH)/*.o $(SRCDIR)/src/fax/obj/$(HOSTARCH)/*.o $(SRCDIR)/src/kass/obj/$(HOSTARCH)/*.o $(SRCDIR)/src/order/obj/$(HOSTARCH)/*.o $(SRCDIR)/src/symbol/obj/$(HOSTARCH)/*.o $(SRCDIR)/src/common/obj/$(HOSTARCH)/*.o
else
	OBJLIBTOOL = 
endif

pastix: FAIRE


WHERE: 
	if [ -f $(SRCDIR)/FAIT  ] ;then \
	echo "double_pastix LD '-L@DIR@/lib/pastix/double   -lpastix'"  >$(SRCDIR)/$(INSTALL)/lib/WHERE.pastix ;\
	echo "double_pastix INCLUDE '-I@DIR@/include/pastix/double'">> $(SRCDIR)/$(INSTALL)/lib/WHERE.patix ;\
	echo "complex_pastix LD '-L@DIR@/lib/pastix/complex   -lpastix'"  >>$(SRCDIR)/$(INSTALL)/lib/WHERE.pastix ;\
	echo "complex_pastix INCLUDE '-I@DIR@/include/pastix/complex'">> $(SRCDIR)/$(INSTALL)/lib/WHERE.patix ;\
	fi

FAIRE: $(SRCDIR)/FAIT  WHERE

config.in: ../Makefile ../../config.status	Makefile
	../../config.status  --file="config.in:config-pastix-real.in"
	../../config.status  --file="config-complex.in:config-pastix-complex.in"


$(SRCDIR)/FAIT: $(SRCDIR) 
	cp config.in $(SRCDIR)/src/
	cd $(SRCDIR)/src; $(MAKE) expor
	#ifndef (,$(OBJLIBTOOL)) # redefinition of "cd $(SRCDIR)/src; $(MAKE) install" for libtool
	rm -f $(SRCDIR)/install/pastix*.h $(SRCDIR)/install/pastix*.in $(SRCDIR)/install/murge*.inc
	$(abs_top_builddir)/download/pastix/$(SRCDIR)/src/utils/bin/${HOSTARCH}/genheader $(PASTIX_DIR)/pastix$(VERSIONTYPE).h $(PASTIX_DIR)/pastix_fortran$(VERSIONTYPE).h \
	$(PASTIX_DIR)/murge$(VERSIONTYPE).inc $(SRCDIR)/src/murge/include/murge.h $(SRCDIR)/src/murge/scripts/genfort.pl
	cat $(SRCDIR)/src/../bin/$(HOSTARCH)/pastix.h         >> $(PASTIX_DIR)/pastix$(VERSIONTYPE).h
	cat $(SRCDIR)/src/common/src/api.h                    >> $(PASTIX_DIR)/pastix$(VERSIONTYPE).h
	cat $(SRCDIR)/src/../bin/$(HOSTARCH)/pastix_fortran.inc >> $(PASTIX_DIR)/murge$(VERSIONTYPE).inc
	ln -sf $(PASTIX_DIR)/pastix$(VERSIONTYPE).h         	$(PASTIX_DIR)/pastix.h
	ln -sf $(PASTIX_DIR)/pastix_fortran$(VERSIONTYPE).h 	$(PASTIX_DIR)/pastix_fortran.h
	ln -sf $(PASTIX_DIR)/murge$(VERSIONTYPE).inc	 	$(PASTIX_DIR)/murge.inc
	cp  $(SRCDIR)/src/sopalin/src/csc_utils.h                     $(PASTIX_DIR)/csc_utils$(VERSIONTYPE).h
	ln -fs $(PASTIX_DIR)/csc_utils$(VERSIONTYPE).h $(PASTIX_DIR)/csc_utils.h
	cp  $(SRCDIR)/src/sopalin/src/cscd_utils.h                    $(PASTIX_DIR)/cscd_utils${VERSIONTYPE}.h
	ln -fs $(PASTIX_DIR)/cscd_utils$(VERSIONTYPE).h $(PASTIX_DIR)/cscd_utils.h
	cp  $(SRCDIR)/src/common/src/nompi.h                          $(PASTIX_DIR)/pastix_nompi.h
	cp  $(SRCDIR)/src/murge/include/murge.h 			 $(PASTIX_DIR)/murge.h
	$(ARPROG) $(ARFLAGS) $(SRCDIR)/install/libpastix$(VERSIONNAME).a $(OBJLIBTOOL)
	ranlib $(SRCDIR)/install/libpastix$(VERSIONNAME).a
	ln -sf $(PASTIX_DIR)/libpastix$(VERSIONNAME).a $(PASTIX_DIR)/libpastix.a
	cp $(PASTIX_DIR)/../bin/$(HOSTARCH)/libpastix_murge.a $(PASTIX_DIR)/libpastix_murge$(VERSIONNAME).a
	ln -sf $(PASTIX_DIR)/libpastix_murge$(VERSIONNAME).a $(PASTIX_DIR)/libpastix_murge.a
	#else
	#cd $(SRCDIR)/src; $(MAKE) install
	#endif 
	mkdir -p $(SRCDIR)/$(INSTALL)/include/pastix/double
	cp $(SRCDIR)/install/*.h $(SRCDIR)/$(INSTALL)/include/pastix/double/
	mkdir -p $(SRCDIR)/$(INSTALL)/lib/pastix/double
	cp $(SRCDIR)/install/*.a $(SRCDIR)/$(INSTALL)/lib/pastix/double/
	cd $(SRCDIR)/src/; $(MAKE) clean
	#complex version
	cp config-complex.in $(SRCDIR)/src/config.in	
	cd $(SRCDIR)/src/;$(MAKE) expor
#ifdef ($(OBJLIBTOOL))
	rm -f $(SRCDIR)/install/pastix*.h $(SRCDIR)/install/pastix*.in $(SRCDIR)/install/murge*.inc
	$(abs_top_builddir)/download/pastix/$(SRCDIR)/src/utils/bin/${HOSTARCH}/genheader $(PASTIX_DIR)/pastix$(VERSIONTYPECOMPLEX).h $(PASTIX_DIR)/pastix_fortran$(VERSIONTYPECOMPLEX).h \
	$(PASTIX_DIR)/murge$(VERSIONTYPECOMPLEX).inc $(SRDIR)/src/murge/include/murge.h $(SRDIR)/src/murge/scripts/genfort.pl
	cat $(SRCDIR)/src/../bin/$(HOSTARCH)/pastix.h         >> $(PASTIX_DIR)/pastix$(VERSIONTYPECOMPLEX).h
	cat $(SRCDIR)/src/common/src/api.h                    >> $(PASTIX_DIR)/pastix$(VERSIONTYPECOMPLEX).h
	cat $(SRCDIR)/src/../bin/$(HOSTARCH)/pastix_fortran.inc >> $(PASTIX_DIR)/murge$(VERSIONTYPECOMPLEX).inc
	ln -sf $(PASTIX_DIR)/pastix$(VERSIONTYPECOMPLEX).h         	$(PASTIX_DIR)/pastix.h
	ln -sf $(PASTIX_DIR)/pastix_fortran$(VERSIONTYPECOMPLEX).h 	$(PASTIX_DIR)/pastix_fortran.h
	ln -sf $(PASTIX_DIR)/murge$(VERSIONTYPECOMPLEX).inc	 	$(PASTIX_DIR)/murge.inc
	cp  $(SRCDIR)/src/sopalin/src/csc_utils.h                     $(PASTIX_DIR)/csc_utils$(VERSIONTYPECOMPLEX).h
	ln -fs $(PASTIX_DIR)/csc_utils$(VERSIONTYPECOMPLEX).h $(PASTIX_DIR)/csc_utils.h
	cp  $(SRCDIR)/src/sopalin/src/cscd_utils.h                    $(PASTIX_DIR)/cscd_utils${VERSIONTYPECOMPLEX}.h
	ln -fs $(PASTIX_DIR)/cscd_utils$(VERSIONTYPECOMPLEX).h $(PASTIX_DIR)/cscd_utils.h
	cp  $(SRCDIR)/src/common/src/nompi.h                          $(PASTIX_DIR)/pastix_nompi.h
	cp  $(SRCDIR)/src/murge/include/murge.h 			 $(PASTIX_DIR)/murge.h
	$(ARPROG) $(ARFLAGS) $(SRCDIR)/install/libpastix$(VERSIONNAME).a $(OBJLIBTOOL)
	ranlib $(SRCDIR)/install/libpastix$(VERSIONNAME).a
	ln -sf $(PASTIX_DIR)/libpastix$(VERSIONNAME).a $(PASTIX_DIR)/libpastix.a
	cp $(PASTIX_DIR)/../bin/$(HOSTARCH)/libpastix_murge.a $(PASTIX_DIR)/libpastix_murge$(VERSIONNAME).a
	ln -sf $(PASTIX_DIR)/libpastix_murge$(VERSIONNAME).a $(PASTIX_DIR)/libpastix_murge.a
#else
	#cd $(SRCDIR)/src/;$(MAKE) install
#endif
	mkdir -p $(SRCDIR)/$(INSTALL)/include/pastix/complex
	cp $(SRCDIR)/install/*.h $(SRCDIR)/$(INSTALL)/include/pastix/complex/
	mkdir -p $(SRCDIR)/$(INSTALL)/lib/pastix/complex
	cp $(SRCDIR)/install/*.a $(SRCDIR)/$(INSTALL)/lib/pastix/complex/
	patch $(SRCDIR)/$(INSTALL)/include/pastix/complex/pastix_long_complex.h   <patch-pastix_long_complex.h
	patch $(SRCDIR)/$(INSTALL)/include/pastix/complex/pastix.h <patch-pastix_long_complex.h
	touch $(SRCDIR)/FAIT



$(SRCDIR): $(PACKAGE1)
	tar xvfj $(PACKAGE1)
	#cd $(SRCDIR)/src/; patch -p1 < ../../pastix_release_2200-gen.patch
	cd $(SRCDIR)/src/blend/src/; patch -p2 < ../../../../pastix_release_2200-blend.patch
	touch $(SRCDIR)

$(PACKAGE1):
	-mkdir $(DIRPKG);
	cd $(DIRPKG);$(WGET)   $(SERVER)/`basename $(PACKAGE1)`

clean-local:
	-cd $(SRCDIR)/src/; $(MAKE) clean

clean: clean-local
	-rm config.in
	-rm $(SRCDIR)/FAIT:
	-rm -rf $(SRCDIR)/$(INSTALL)/include/pastix/double
	-rm -rf $(SRCDIR)/$(INSTALL)/lib/pastix/double
	-rm -rf $(SRCDIR)/$(INSTALL)/include/pastix/complex
	-rm -rf $(SRCDIR)/$(INSTALL)/lib/pastix/complex
	-rm -rf $(SRCDIR)

.PHONY:$(SRCDIR)/$(INSTALL)