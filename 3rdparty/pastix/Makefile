# ======================================================================
# Laboratoire Jacques-Louis Lions
# Université Pierre et Marie Curie-Paris6, UMR 7598, Paris, F-75005 France
# http://www.ljll.math.upmc.fr/lehyaric
# ======================================================================
# This file is part of Freefem++
# 
# Freefem++ is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of
# the License, or (at your option) any later version.
# 
# Freefem++ is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with Freefem++; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
# ======================================================================
# headeralh default=0 freefem make multipleauthors start=19/03/10 upmc

# Downloading and compiling extra libraries
# -----------------------------------------

-include config.in

all-local: pastix

# Downloading and compiling pastix
# ------------------------------
DIRPKG   = ../pkg
SRCDIR   = pastix_release_$(VERSION)
PACKAGE1 = $(DIRPKG)/pastix_release_$(VERSION).tar.bz2
SERVER   = https://gforge.inria.fr/frs/download.php/21873/
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


# FFCS - parallel make rules
WHERE: install
	echo double_pastix LD -L@DIR@/lib/pastix/double   -lpastix   >../lib/WHERE.pastix ;
	echo double_pastix INCLUDE -I@DIR@/include/pastix/double >> ../lib/WHERE.pastix ;
	echo complex_pastix LD -L@DIR@/lib/pastix/complex   -lpastix   >>../lib/WHERE.pastix ;
	echo complex_pastix INCLUDE -I@DIR@/include/pastix/complex >> ../lib/WHERE.pastix ;

FAIRE:$(SRCDIR)/FAIT 
	$(MAKE) WHERE
	touch FAIRE

config.in: ../Makefile ../../config.status	Makefile
	../../config.status  --file="config.in:config-pastix-real.in"
	../../config.status  --file="config-complex.in:config-pastix-complex.in"

install:$(SRCDIR)/FAIT

$(SRCDIR)/FAIT:  $(SRCDIR)/FAITC2
	touch $(SRCDIR)/FAIT

## change F. Hecht
## WARNING target $(SRCDIR)/FAITR2 $(SRCDIR)/FAITC2  not in //
#  => $(SRCDIR)/FAITC2 depend on $(SRCDIR)/FAITR2 

$(SRCDIR)/FAITR1: $(SRCDIR)/tag-tar
	cp config.in $(SRCDIR)/src/
#
#	FFCS - 16/1/13 - this crashes in parallel
#
	cd $(SRCDIR)/src; $(MAKE) -j1 expor
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
	#cd $(SRCDIR)/src; $(MAKE) -j1 install
	#endif 
	touch $(SRCDIR)/FAITR1
$(SRCDIR)/FAITR2: $(SRCDIR)/FAITR1
	mkdir -p ../include/pastix/double
	cp $(SRCDIR)/install/*.h ../include/pastix/double/
	mkdir -p ../lib/pastix/double
	cp $(SRCDIR)/install/*.a ../lib/pastix/double/
	cd $(SRCDIR)/src/ &&  $(MAKE) clean 
	touch $(SRCDIR)/FAITR2;

$(SRCDIR)/FAITC1:$(SRCDIR)/tag-tar
#
# complex version
#
	cp config-complex.in $(SRCDIR)/src/config.in	
#
#	FFCS - 16/1/13 - this crashes in parallel
#
	cd $(SRCDIR)/src/;$(MAKE) -j1 expor
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
	#cd $(SRCDIR)/src/ && $(MAKE) install -j1 -C $(SRCDIR)/src
#endif
	touch $(SRCDIR)/FAITC1

$(SRCDIR)/FAITC2: $(SRCDIR)/FAITR2 $(SRCDIR)/FAITC1
	mkdir -p ../include/pastix/complex
	cp $(SRCDIR)/install/*.h ../include/pastix/complex/
	mkdir -p ../lib/pastix/complex
	cp $(SRCDIR)/install/*.a ../lib/pastix/complex/
	(cd ../include/pastix/complex/; patch pastix_int_complex.h)   <patch-pastix_long_complex.h
	(cd ../include/pastix/complex/; patch pastix.h) <patch-pastix_long_complex.h
	touch $(SRCDIR)/FAITC2



$(SRCDIR) $(SRCDIR)/tag-tar: $(PACKAGE1)
	tar xvfj $(PACKAGE1)
#
# 31/7/13 - ALH - patch from Fred to avoid picking up a locally installed version of Umfpack when FF compiles its own on
# Slackware64-14.0. For more information see FFCS comments from Fred on
# [[http://www.ljll.math.upmc.fr//~lehyaric/ffcs/news/20130528-1139.php]] (dated 30/7/13)
#
	patch -p0 < all_macros.diff
#
	cd $(SRCDIR)/src/blend/src/; patch -p2 < ../../../../pastix_release_2200-blend.patch
	touch $(SRCDIR) touch $(SRCDIR)/tag-tar 

$(PACKAGE1):
	../getall -o PaStiX -a

# FFCS: avoid loops when SRCDIR does not exist
clean-local:
	-cd $(SRCDIR)/src &&  $(MAKE) clean 

clean: clean-local
	-rm config.in FAIRE 
	-rm $(SRCDIR)/FAIT
	-rm -rf ../include/pastix/double
	-rm -rf ../lib/pastix/double
	-rm -rf ../include/pastix/complex
	-rm -rf ../lib/pastix/complex
	-rm -rf $(SRCDIR)
#
# FFCS - also force reconfiguration
#
	-rm config.in config-complex.in config.log
	-rm config.log
#	-rm $(PACKAGE1)

# Local Variables:
# mode:makefile
# ispell-local-dictionary:"british"
# coding:utf-8
# End:
