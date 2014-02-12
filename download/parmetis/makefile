# Downloading and compiling extra libraries
# -----------------------------------------

include Makefile.in

all-local: parmetis

# Downloading and compiling Tetgen
# ------------------------------
# http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/parmetis-4.0.tar.gz
# Parmetis information
DIRPKG=../pkg
SRCDIR=ParMetis-$(parmetis_VERSION)
PACKAGE=$(DIRPKG)/ParMetis-$(parmetis_VERSION).tar.gz
SERVER=http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/OLD

INSTALL=../..
parmetis_VERSION=3.1.1

parmetis: WHERE.done

compile.done:$(SRCDIR)/tags Makefile.in
	-mkdir -p "../lib/parmetis"
	-mkdir -p "../include/parmetis"
	cp Makefile.in $(SRCDIR)
# 	FFCS - we need to know when errors occur
	cd $(SRCDIR)/METISLib && $(MAKE) 
	cd $(SRCDIR)/ParMETISLib && $(MAKE) 
	touch $@

# FFCS - simpler makefile
install.done:compile.done
#	cp $(SRCDIR)/libparmetis.a $(SRCDIR)/$(INSTALL)/lib
#	cp $(SRCDIR)/libmetis.a $(SRCDIR)/$(INSTALL)/lib
	cp $(SRCDIR)/parmetis.h ../include/parmetis
	cp $(SRCDIR)/METISLib/*.h ../include/parmetis
	sed  's;../parmetis.h;parmetis.h;' <$(SRCDIR)/METISLib/metis.h >../include/parmetis/metis.h
#
#	FFCS - 23/5/12 - cannot keep name libmetis.a because it is identical to the library created by
#	[[file:../metis/Makefile]] and library path mechanisms at link time pick one for the other on MinGW.
	mv $(SRCDIR)/$(INSTALL)/lib/parmetis/libmetis.a $(SRCDIR)/$(INSTALL)/lib/parmetis/libmetispar.a
#
	-cd $(SRCDIR)/Programs && make 
	touch $@

# FFCS - simpler makefile
WHERE.done:install.done
	echo parmetis LD -L@DIR@/lib/parmetis -lparmetis -lmetispar  >$(SRCDIR)/$(INSTALL)/lib/WHERE.parmetis ;
	echo parmetis INCLUDE -I@DIR@/include/parmetis   >> $(SRCDIR)/$(INSTALL)/lib/WHERE.parmetis ;
	touch $@

Makefile.in: ../../config.status Makefile-parmetis.in
	../../config.status  --file="Makefile.in:Makefile-parmetis.in"



$(SRCDIR)/tags: $(PACKAGE)
	tar xvzf $(PACKAGE)
	patch -p0  <parmetis-3.1.1.patch
	mv  $(SRCDIR)/Programs/Makefile $(SRCDIR)/Programs/Makefile-orig
	sed -e 's;BINDIR = .*$$;BINDIR = ../$(INSTALL)/bin;' \
            -e 's;../libparmetis.a;../$(INSTALL)/lib/parmetis/libparmetis.a;' \
            -e 's;../libmetis.a;../$(INSTALL)/lib/parmetis/libmetis.a;' \
            -e 's;-L[.][.];-L../$(INSTALL)/lib/parmetis;' \
          <$(SRCDIR)/Programs/Makefile-orig \
          >$(SRCDIR)/Programs/Makefile
	mv $(SRCDIR)/ParMETISLib/Makefile $(SRCDIR)/ParMETISLib/Makefile-orig
	sed   -e 's;../libparmetis.a;../$(INSTALL)/lib/parmetis/libparmetis.a;' \
	      -e 's;../libmetis.a;../$(INSTALL)/lib/parmetis/libmetis.a;' \
          <$(SRCDIR)/ParMETISLib/Makefile-orig \
          >$(SRCDIR)/ParMETISLib/Makefile
	mv $(SRCDIR)/METISLib/Makefile $(SRCDIR)/METISLib/Makefile-orig
	sed   -e 's;../libparmetis.a;../$(INSTALL)/lib/parmetis/libparmetis.a;' \
	      -e 's;../libmetis.a;../$(INSTALL)/lib/parmetis/libmetis.a;' \
          <$(SRCDIR)/METISLib/Makefile-orig \
          >$(SRCDIR)/METISLib/Makefile
	touch $(SRCDIR)/tags

$(PACKAGE):
	../getall -o ParMETIS -a

# FFCS - 23/5/12 - corrected bug in removing .a from ../lib
clean:
	-cd $(SRCDIR) && $(MAKE) realclean  -C  $(SRCDIR)
	-rm Makefile.in $(SRCDIR)/tags
	-rm -rf ../lib/parmetis
	-rm -rf ../include/parmetis
	-rm -rf ../include/metis
	-rm -rf $(SRCDIR)
	-rm *.done
